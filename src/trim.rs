// This file contains functions related to trimming paths to remove start-end overlap.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

#![allow(clippy::needless_range_loop)]

use colored::Colorize;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::{HashMap, VecDeque};
use std::fmt;
use std::path::{Path, PathBuf};

use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::metrics::TrimmedClusterMetrics;
use crate::misc::{check_if_dir_exists, check_if_file_exists, format_float, quit_with_error,
                  median_isize, mad_isize, reverse_path};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


const GAP: i32 = 0;
const NONE: usize = usize::MAX;


pub fn trim(cluster_dir: PathBuf, min_identity: f64, max_unitigs: usize, mad: f64, threads: usize) {
    let untrimmed_gfa = cluster_dir.join("1_untrimmed.gfa");
    let trimmed_gfa = cluster_dir.join("2_trimmed.gfa");
    let trimmed_yaml = cluster_dir.join("2_trimmed.yaml");
    check_settings(&cluster_dir, &untrimmed_gfa, min_identity, mad, threads);
    starting_message();
    print_settings(&cluster_dir, min_identity, max_unitigs, mad, threads);
    let (mut graph, sequences) = load_graph(&untrimmed_gfa);
    let unitig_lengths: HashMap<_, _> = graph.unitigs.iter().map(|rc| {let u = rc.borrow(); (u.number as i32, u.length())}).collect();
    let start_end_results = trim_start_end_overlap(&graph, &sequences, &unitig_lengths, min_identity, max_unitigs);
    let hairpin_results = trim_harpin_overlap(&graph, &sequences, &unitig_lengths, min_identity, max_unitigs);
    let sequences = choose_trim_type(start_end_results, hairpin_results, &mut graph, &sequences);
    let sequences = exclude_outliers_in_length(&mut graph, &sequences, mad);
    clean_up_graph(&mut graph, &sequences);
    graph.save_gfa(&trimmed_gfa, &sequences, false).unwrap();
    save_metrics(&trimmed_yaml, &sequences);
    finished_message(&trimmed_gfa);
}


fn check_settings(cluster_dir: &Path, untrimmed_gfa: &Path, min_identity: f64, mad: f64,
                  threads: usize) {
    check_if_dir_exists(cluster_dir);
    check_if_file_exists(untrimmed_gfa);
    if !(0.0..=1.0).contains(&min_identity) {
        quit_with_error("--min_identity must be between 0.0 and 1 (inclusive)");
    }
    if threads < 1   { quit_with_error("--threads cannot be less than 1"); }
    if threads > 100 { quit_with_error("--threads cannot be greater than 100"); }
    if mad < 0.0     { quit_with_error("--mad cannot be less than 0"); }
    ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
}


fn starting_message() {
    section_header("Starting autocycler trim");
    explanation("This command takes a single-cluster unitig graph (made by autocycler cluster) and \
                 trims any overlaps. It looks for both start-end overlaps (can occur with circular \
                 sequences) and hairpin overlaps (can occur with linear sequences).");
}


fn finished_message(trimmed_gfa: &Path) {
    section_header("Finished!");
    explanation("You can now run autocycler resolve on this cluster. If you want to manually \
                 inspect the trimming, you can run autocycler dotplot on the sequences both before \
                 and after trimming.");
    eprintln!("Unitig graph of trimmed sequences: {}", trimmed_gfa.display());
    eprintln!();
}


fn print_settings(cluster_dir: &Path, min_identity: f64, max_unitigs: usize, mad: f64,
                  threads: usize) {
    eprintln!("Settings:");
    eprintln!("  --cluster_dir {}", cluster_dir.display());
    eprintln!("  --min_identity {}", format_float(min_identity));
    eprintln!("  --max_unitigs {}", max_unitigs);
    eprintln!("  --mad {}", format_float(mad));
    eprintln!("  --threads {}", threads);
    eprintln!();
    if max_unitigs == 0 {
        eprintln!("Since --max_unitigs was set to 0, trimming is disabled.");
        eprintln!();
    }
}


fn load_graph(gfa: &Path) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading graph");
    explanation("The unitig graph is now loaded into memory.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(gfa);
    unitig_graph.print_basic_graph_info();
    (unitig_graph, sequences)
}


fn trim_start_end_overlap(graph: &UnitigGraph, sequences: &Vec<Sequence>, weights: &HashMap<i32, u32>,
                          min_identity: f64, max_unitigs: usize) -> Vec<Option<(Vec<i32>, u32)>> {
    if max_unitigs == 0 {
        return vec![None; sequences.len()];
    }
    section_header("Trim start-end overlaps");
    explanation("Paths for circular replicons may contain start-end overlaps. These overlaps \
                 are searched for and trimmed if found.");
    let paths: Vec<_> = sequences.iter().map(|seq| graph.get_unitig_path_for_sequence_i32(seq)).collect();
    let results: Vec<_> = sequences.par_iter().zip(paths.par_iter()).map(|(seq, path)| {  // parallel for loop with rayon
        let trimmed_path = trim_path_start_end(path, weights, min_identity, max_unitigs);
        if let Some(trimmed_path) = trimmed_path {
            let trimmed_length: u32 = trimmed_path.iter().map(|&u| weights[&u.abs()]).sum();
            (Some((trimmed_path, trimmed_length)), format!("{}: {}", seq, format!("trimmed to {} bp", trimmed_length).red()))
        } else {
            (None, format!("{}: {}", seq, "not trimmed".green()))
        }
    }).collect();
    for (_, message) in &results {
        eprintln!("{}", message);
    }
    eprintln!();
    results.into_iter().map(|(result, _)| result).collect()
}


fn trim_harpin_overlap(graph: &UnitigGraph, sequences: &Vec<Sequence>, weights: &HashMap<i32, u32>,
                       min_identity: f64, max_unitigs: usize) -> Vec<Option<(Vec<i32>, u32)>> {
    if max_unitigs == 0 {
        return vec![None; sequences.len()];
    }
    section_header("Trim hairpin overlaps");
    explanation("Paths for linear replicons may contain hairpin overlaps at the start and/or end \
                 of the contig. These overlaps are searched for and trimmed if found.");
    let paths: Vec<_> = sequences.iter().map(|seq| graph.get_unitig_path_for_sequence_i32(seq)).collect();
    let results: Vec<_> = sequences.par_iter().zip(paths.par_iter()).map(|(seq, path)| {  // parallel for loop with rayon
        let mut trimmed_start = false;
        let mut trimmed_end = false;

        let path_2;
        let trimmed_path_start = trim_path_hairpin_start(path, weights, min_identity, max_unitigs);
        if trimmed_path_start.is_some() {
            trimmed_start = true;
            path_2 = trimmed_path_start.unwrap();
        } else {
            path_2 = path.clone();
        }

        let path_3;
        let trimmed_path_end = trim_path_hairpin_end(&path_2, weights, min_identity, max_unitigs);
        if trimmed_path_end.is_some() {
            trimmed_end = true;
            path_3 = trimmed_path_end.unwrap();
        } else {
            path_3 = path_2;
        }

        if !trimmed_start && !trimmed_end {
            (None, format!("{}: {}", seq, "not trimmed".green()))
        } else {
            let message;
            let trimmed_length: u32 = path_3.iter().map(|&u| weights[&u.abs()]).sum();
            if trimmed_start && trimmed_end {
                message = format!("{}: {}", seq, format!("trimmed from start and end to {} bp", trimmed_length).red());
            } else if trimmed_start {
                message = format!("{}: {}", seq, format!("trimmed from start to {} bp", trimmed_length).red());
            } else {
                message = format!("{}: {}", seq, format!("trimmed from end to {} bp", trimmed_length).red());
            }
            (Some((path_3, trimmed_length)), message)
        }
    }).collect();
    for (_, message) in &results {
        eprintln!("{}", message);
    }
    eprintln!();
    results.into_iter().map(|(result, _)| result).collect()
}


fn choose_trim_type(start_end_results: Vec<Option<(Vec<i32>, u32)>>, hairpin_results: Vec<Option<(Vec<i32>, u32)>>,
                    graph: &mut UnitigGraph, sequences: &[Sequence]) -> Vec<Sequence> {
    let start_end_count = start_end_results.iter().filter(|x| x.is_some()).count();
    let hairpin_count = hairpin_results.iter().filter(|x| x.is_some()).count();
    if start_end_count == 0 && hairpin_count == 0 {
        return sequences.to_owned();
    }

    let mut trimmed_sequences = vec![];
    let results;
    if start_end_count >= hairpin_count {
        results = start_end_results;
        if hairpin_count > 0 {
            eprintln!("Start-end trimming was more successful than hairpin trimming. Discarding \
                       hairpin trimming.\n");
        }
    } else {  // hairpin_count > start_end_count
        results = hairpin_results;
        if start_end_count > 0 {
            eprintln!("Hairpin trimming was more successful than start-end trimming. Discarding \
                      start-end trimming.\n");
        }
    }
    assert!(sequences.len() == results.len());

    for (seq, result) in sequences.iter().zip(results.iter()) {
        if result.is_none() {
            trimmed_sequences.push(seq.clone());
        } else {
            graph.remove_sequence_from_graph(seq.id);
            let (path, trimmed_length) = result.as_ref().unwrap();
            let trimmed_sequence = graph.create_sequence_and_positions(seq.id, *trimmed_length, seq.filename.clone(),
                                                                        seq.contig_header.clone(), seq.cluster, path_to_tuples(path));
            trimmed_sequences.push(trimmed_sequence);
        }
    }
    trimmed_sequences
}


fn exclude_outliers_in_length(graph: &mut UnitigGraph, sequences: &Vec<Sequence>,
                              mad_threshold: f64) -> Vec<Sequence> {
    if mad_threshold == 0.0 {
        return sequences.clone();
    }
    section_header("Exclude outliers");
    explanation("Sequences which vary too much in their length are now excluded from the cluster.");
    let lengths: Vec<_> = sequences.iter().map(|s| s.length as isize).collect();
    let median = median_isize(&lengths);
    let median_absolute_deviation = mad_isize(&lengths);
    let min_length = (median as f64 - (median_absolute_deviation as f64 * mad_threshold)).round() as usize;
    let max_length = (median as f64 + (median_absolute_deviation as f64 * mad_threshold)).round() as usize;
    eprintln!("Median sequence length:    {} bp", median);
    eprintln!("Median absolute deviation: {} bp", median_absolute_deviation);
    eprintln!("Allowed length range:      {}-{} bp", min_length, max_length);
    eprintln!();
    let mut new_sequences = vec![];
    for seq in sequences {
        if min_length <= seq.length && seq.length <= max_length {
            new_sequences.push(seq.clone());
            eprintln!("{}: {}", seq, "kept".green());
        } else {
            eprintln!("{} {}", format!("{}:", seq).dimmed(), "excluded".red());
            graph.remove_sequence_from_graph(seq.id);
        }
    }
    eprintln!();
    new_sequences
}


fn clean_up_graph(graph: &mut UnitigGraph, sequences: &Vec<Sequence>) {
    section_header("Clean graph");
    explanation("The unitig graph is now cleaned up based on any trimming and/or exclusion that \
                 has occurred above.");
    graph.recalculate_depths();
    graph.remove_zero_depth_unitigs();
    merge_linear_paths(graph, sequences);
    graph.print_basic_graph_info();
    graph.renumber_unitigs();
}


fn save_metrics(trimmed_yaml: &Path, sequences: &[Sequence]) {
    let seq_lengths = sequences.iter().map(|s| s.length).collect();
    let metrics = TrimmedClusterMetrics::new(seq_lengths);
    metrics.save_to_yaml(trimmed_yaml);
}


fn path_to_tuples(path: &[i32]) -> Vec<(u32, bool)> {
    // Paths can be represented either as a vector of tuples or as vector of signed integers:
    //   [(1, true), (2, false), (3, true)]
    //   [1, -2, 3]
    // This function converts the latter to the former.
    path.iter().map(|p| if *p > 0 {(*p as u32, true)} else {(-*p as u32, false)}).collect()
}


fn trim_path_start_end(path: &[i32], weights: &HashMap<i32, u32>, min_identity: f64,
                       max_unitigs: usize) -> Option<Vec<i32>> {
    let alignment = overlap_alignment(path, path, weights, min_identity, max_unitigs, true);
    if alignment.is_empty() { return None; }
    let midpoint_index = find_midpoint(&alignment, weights);
    let start = alignment[midpoint_index].a_index;
    let end = alignment[midpoint_index].b_index;
    Some(path[start..end].to_vec())
}


fn trim_path_hairpin_end(path: &[i32], weights: &HashMap<i32, u32>, min_identity: f64,
                         max_unitigs: usize) -> Option<Vec<i32>> {
    let rev_path = reverse_path(path);
    let mut alignment = overlap_alignment(&rev_path, path, weights, min_identity, max_unitigs, false);
    if alignment.is_empty() { return None; }
    let mut end = 0;
    while !alignment.is_empty() {
        trim_gaps_a_front(&mut alignment);
        trim_gaps_b_back(&mut alignment);
        if alignment.is_empty() { break; }
        let back = alignment.pop_back().unwrap();
        assert!(back.b_unitig == -alignment.front().unwrap().a_unitig);
        if back.a_unitig != GAP {
            end = back.b_index;
        }
        alignment.pop_front();
    }
    Some(path[..end].to_vec())
}


fn trim_path_hairpin_start(path: &[i32], weights: &HashMap<i32, u32>, min_identity: f64,
                           max_unitigs: usize) -> Option<Vec<i32>> {
    let rev_path = reverse_path(path);
    let trimmed_reverse_path = trim_path_hairpin_end(&rev_path, weights, min_identity, max_unitigs);
    trimmed_reverse_path.as_ref()?;
    Some(reverse_path(&trimmed_reverse_path.unwrap()))
}


#[derive(PartialEq, Eq)]
struct AlignmentPiece {
    a_unitig: i32,
    a_index: usize,
    b_unitig: i32,
    b_index: usize,
}

impl fmt::Display for AlignmentPiece {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let a_unitig = if self.a_unitig == 0 { "GAP".to_string() } else { self.a_unitig.to_string() };
        let b_unitig = if self.b_unitig == 0 { "GAP".to_string() } else { self.b_unitig.to_string() };
        let a_index = if self.a_index == usize::MAX { "NONE".to_string() } else { self.a_index.to_string() };
        let b_index = if self.b_index == usize::MAX { "NONE".to_string() } else { self.b_index.to_string() };
        write!(f, "{},{},{},{}", a_unitig, a_index, b_unitig, b_index)
    }
}

impl fmt::Debug for AlignmentPiece {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { fmt::Display::fmt(self, f) }
}


fn trim_gaps_a_front(alignment: &mut VecDeque<AlignmentPiece>) {
    while !alignment.is_empty() && alignment.front().unwrap().a_unitig == GAP {
        alignment.pop_front();
    }
}


fn trim_gaps_b_back(alignment: &mut VecDeque<AlignmentPiece>) {
    while !alignment.is_empty() && alignment.back().unwrap().b_unitig == GAP {
        alignment.pop_back();
    }
}


fn overlap_alignment(path_a: &[i32], path_b: &[i32], weights: &HashMap<i32, u32>, min_identity: f64,
                     max_unitigs: usize, skip_diagonal: bool) -> VecDeque<AlignmentPiece> {
    // This function performs a dynamic-programming alignment to find overlaps between two paths,
    // with some special logic for Autocycler trim:
    // * Scores are weighted by the length of the unitig (matches are positive scores, mismatches
    //   and indels are negative scores).
    // * Looks only for overlap alignments that extend from the right edge of the matrix to the top
    //   edge of the matrix.
    // * The matrix is limited in size (by max_unitigs) to prevent bad scaling with huge paths.
    // * Only sufficiently high identity alignments are returned (min_identity).
    // This function can be used on a path vs itself (same strand), in which case path_a and path_b
    // will be the same. Or it can be used on a path vs its opposite strand, in which case path_b
    // will be the reverse complement of path_a.
    assert!(path_a.len() == path_b.len());
    let n = path_a.len();
    let k = max_unitigs.min(n);
    let mut scoring_matrix = vec![vec![-f64::INFINITY; k + 1]; k + 1];

    // Initialize the top and left edges to zero.
    for i in 0..=k {
        scoring_matrix[i][0] = 0.0;
        scoring_matrix[0][i] = 0.0;
    }

    // Fill the scoring matrix.
    for i in 1..=k {
        for j in 1..=k {
            let global_i = i - 1;
            let global_j = n - k + j - 1;
            if skip_diagonal && global_i == global_j {continue;}  // skipping the diagonal avoids whole-vs-whole alignment
            let weight_i = *weights.get(&path_a[global_i].abs()).unwrap() as f64;
            let weight_j = *weights.get(&path_b[global_j].abs()).unwrap() as f64;
            let match_score = scoring_matrix[i - 1][j - 1] + if path_a[global_i] == path_b[global_j] {
                weight_i
            } else {
                -(weight_i + weight_j) / 2.0
            };
            let delete_score = scoring_matrix[i - 1][j] - weight_i;
            let insert_score = scoring_matrix[i][j - 1] - weight_j;
            scoring_matrix[i][j] = match_score.max(delete_score).max(insert_score);
        }
    }

    // Find the maximum score from the right edge.
    let mut max_score = f64::NEG_INFINITY;
    let mut max_i = 0;
    let mut max_j = 0;
    for i in 1..=k {
        if scoring_matrix[i][k] > max_score {
            max_score = scoring_matrix[i][k];
            max_i = i;
            max_j = k;
        }
    }

    // A negative score indicates a very poor alignment.
    if max_score <= 0.0 {
        return VecDeque::new();
    }

    // Traceback to get the alignment and record indices
    let mut alignment_a = Vec::new();
    let mut alignment_b = Vec::new();
    let mut alignment_a_indices = Vec::new();
    let mut alignment_b_indices = Vec::new();
    let mut i = max_i;
    let mut j = max_j;
    while i > 0 && j > 0 {
        let global_i = i - 1;
        let global_j = n - k + j - 1;
        if path_a[global_i] == path_b[global_j] {
            alignment_a.push(path_a[global_i]);
            alignment_b.push(path_b[global_j]);
            alignment_a_indices.push(global_i);
            alignment_b_indices.push(global_j);
            i -= 1;
            j -= 1;
        } else if scoring_matrix[i - 1][j] >= scoring_matrix[i][j - 1] {
            alignment_a.push(path_a[global_i]);
            alignment_b.push(GAP);
            alignment_a_indices.push(global_i);
            alignment_b_indices.push(NONE);
            i -= 1;
        } else {
            alignment_a.push(GAP);
            alignment_b.push(path_b[global_j]);
            alignment_a_indices.push(NONE);
            alignment_b_indices.push(global_j);
            j -= 1;
        }
    }

    // Ensure that the traceback hit the top edge and not the left edge.
    if i > 0 {
        return VecDeque::new();
    }

    alignment_a.reverse();
    alignment_b.reverse();
    alignment_a_indices.reverse();
    alignment_b_indices.reverse();

    let alignment_a_length: u32 = alignment_a.iter().filter_map(|&u| if u != GAP { Some(weights[&u.abs()]) } else { None }).sum();
    let alignment_b_length: u32 = alignment_b.iter().filter_map(|&u| if u != GAP { Some(weights[&u.abs()]) } else { None }).sum();
    let mean_length = (alignment_a_length as f64 + alignment_b_length as f64) / 2.0;
    let total_matches = alignment_a.iter().zip(alignment_b.iter()).fold(0, |acc, (a, b)| { if a == b { acc + weights[&a.abs()] } else { acc }});
    let alignment_identity = total_matches as f64 / mean_length;
    if alignment_identity < min_identity {
        return VecDeque::new();
    }

    alignment_a.iter().zip(alignment_a_indices.iter()).zip(alignment_b.iter()).zip(alignment_b_indices.iter())
        .map(|(((a_u, a_i), b_u), b_i)| AlignmentPiece {a_unitig: *a_u, a_index: *a_i, b_unitig: *b_u, b_index: *b_i}).collect()
}


fn find_midpoint(alignment: &VecDeque<AlignmentPiece>, weights: &HashMap<i32, u32>) -> usize {
    let total_weight = alignment.iter().map(|p| {
        let mut weight = 0;
        if p.a_unitig != GAP { weight += weights[&p.a_unitig.abs()]; }
        if p.b_unitig != GAP { weight += weights[&p.b_unitig.abs()]; }
        weight
    }).sum::<u32>();
    let mut cumulative_weight = 0;
    let mut best_index = 0;
    let mut best_closeness = 1.0;

    for (i, p) in alignment.iter().enumerate() {
        if p.a_unitig != GAP {
            cumulative_weight += weights[&p.a_unitig.abs()];
        }
        if p.b_unitig != GAP {
            cumulative_weight += weights[&p.b_unitig.abs()];
        }
        let closeness_to_midpoint = (0.5 - (cumulative_weight as f64 / total_weight as f64)).abs();
        if p.a_unitig == p.b_unitig && closeness_to_midpoint < best_closeness {
            best_index = i;
            best_closeness = closeness_to_midpoint;
        }
    }
    best_index
}


#[cfg(test)]
mod tests {
    use maplit::hashmap;
    use super::*;
    use crate::misc::strand;

    #[test]
    fn test_path_to_tuples() {
        let path = vec![1, 2, -3];
        assert_eq!(path_to_tuples(&path), vec![(1, strand::FORWARD), (2, strand::FORWARD), (3, strand::REVERSE)]);

        let path = vec![-4, 5, -6];
        assert_eq!(path_to_tuples(&path), vec![(4, strand::REVERSE), (5, strand::FORWARD), (6, strand::REVERSE)]);

        assert_eq!(path_to_tuples(&[]), vec![]);
    }

    #[test]
    fn test_overlap_alignment() {
        // No alignment
        let path = vec![1, -2, 3, -4, 5];
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10};
        let alignment = overlap_alignment(&path, &path, &weights, 0.9, 100, true);
        assert!(alignment.is_empty());

        // Exact overlap of two unitigs
        let path = vec![1, -2, 3, -4, 5, 1, -2];
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10};
        let alignment = overlap_alignment(&path, &path, &weights, 0.9, 100, true);
        let expected_alignment = VecDeque::from(vec![
            AlignmentPiece { a_unitig: 1,  a_index: 0, b_unitig: 1,  b_index: 5 },
            AlignmentPiece { a_unitig: -2, a_index: 1, b_unitig: -2, b_index: 6 },
        ]);
        assert_eq!(alignment, expected_alignment);

        // Same as above but with max_unitigs of 4.
        let path = vec![1, -2, 3, -4, 5, 1, -2];
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10};
        let alignment = overlap_alignment(&path, &path, &weights, 0.9, 4, true);
        let expected_alignment = VecDeque::from(vec![
            AlignmentPiece { a_unitig: 1,  a_index: 0, b_unitig: 1,  b_index: 5 },
            AlignmentPiece { a_unitig: -2, a_index: 1, b_unitig: -2, b_index: 6 },
        ]);
        assert_eq!(alignment, expected_alignment);

        // Same as above but with max_unitigs of 2 (will just manage to catch the overlap)
        let path = vec![1, -2, 3, -4, 5, 1, -2];
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10};
        let alignment = overlap_alignment(&path, &path, &weights, 0.9, 2, true);
        let expected_alignment = VecDeque::from(vec![
            AlignmentPiece { a_unitig: 1,  a_index: 0, b_unitig: 1,  b_index: 5 },
            AlignmentPiece { a_unitig: -2, a_index: 1, b_unitig: -2, b_index: 6 },
        ]);
        assert_eq!(alignment, expected_alignment);

        // Same as above but with max_unitigs of 1 resulting in no alignment
        let path = vec![1, -2, 3, -4, 5, 1, -2];
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10};
        let alignment = overlap_alignment(&path, &path, &weights, 0.9, 1, true);
        assert!(alignment.is_empty());

        // Inexact overlap of three unitigs
        let path = vec![1, -2, 3, -4, 5, 1, 6, 3];
        let weights = hashmap!{1 => 30, 2 => 1, 3 => 10, 4 => 10, 5 => 10, 6 => 1};
        let alignment = overlap_alignment(&path, &path, &weights, 0.9, 100, true);
        let expected_alignment = VecDeque::from(vec![
            AlignmentPiece { a_unitig: 1,   a_index: 0,    b_unitig: 1,   b_index: 5 },
            AlignmentPiece { a_unitig: GAP, a_index: NONE, b_unitig: 6,   b_index: 6 },
            AlignmentPiece { a_unitig: -2,  a_index: 1,    b_unitig: GAP, b_index: NONE },
            AlignmentPiece { a_unitig: 3,   a_index: 2,    b_unitig: 3,   b_index: 7 },
        ]);
        assert_eq!(alignment, expected_alignment);

        // Same as above but with a high min_alignment_identity resulting in no alignment
        let path = vec![1, -2, 3, -4, 5, 1, 6, 3];
        let weights = hashmap!{1 => 30, 2 => 1, 3 => 10, 4 => 10, 5 => 10, 6 => 1};
        let alignment = overlap_alignment(&path, &path, &weights, 0.99, 100, true);
        assert!(alignment.is_empty());

        // Same as above but with max_unitigs of 2 resulting in no alignment
        let path = vec![1, -2, 3, -4, 5, 1, 6, 3];
        let weights = hashmap!{1 => 30, 2 => 1, 3 => 10, 4 => 10, 5 => 10, 6 => 1};
        let alignment = overlap_alignment(&path, &path, &weights, 0.9, 2, true);
        assert!(alignment.is_empty());
    }

    #[test]
    fn test_trim_path_start_end_1() {
        let weights = hashmap!{653 => 541, 728 => 413, 757 => 366, 977 => 185, 1010 => 170, 1058 => 153, 1105 => 138, 1133 => 133, 1492 => 79, 1552 => 74, 1637 => 68, 1667 => 65, 1913 => 51, 1943 => 50, 1949 => 50, 1952 => 50, 1967 => 50, 1982 => 50, 1993 => 50, 2012 => 49, 2018 => 48, 2065 => 45, 2070 => 45, 2110 => 42, 2148 => 39, 2276 => 32, 2289 => 32, 2499 => 25, 2640 => 21, 2826 => 15, 2937 => 11, 3148 => 6, 3208 => 5, 3456 => 2, 3578 => 2, 4216 => 1, 4238 => 1, 4575 => 1, 4875 => 1, 4876 => 1, 5191 => 1};

        let path = vec![-653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());

        let path = vec![-1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());

        let path = vec![-728, 2640, -4216, 2640, -4216, 2640, 757, -1993, 977, 1492, 1913, -2018, 3456, -4876, 653, -1943, -1058, -2070, 2148, -2012, 2065, -3578, -2110, 3148, -2937, 2276, 5191, -1952, -2499, -1105, 1133, 1949, 1667, -2826, 1010, -1637, -1982, -4875, 2289, 4575, 1552, 4238, -1967, -728, 2640, -4216, 2640, -4216, 2640, 757, -1993, 977, 1492, 1913, -2018, 3456, -4876, 653, -1943, -1058, -2070, 2148, -2012, 2065, -3578, -2110, 3148, -2937, 2276, 5191, -1952, -2499, -1105, 1133, 1949, 1667, -2826, 1010, -1637, -1982, -4875, 2289, 4575, 1552, 4238];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![653, -1943, -1058, -2070, 2148, -2012, 2065, -3578, -2110, 3148, -2937, 2276, 5191, -1952, -2499, -1105, 1133, 1949, 1667, -2826, 1010, -1637, -1982, -4875, 2289, 4575, 1552, 4238, -1967, -728, 2640, -4216, 2640, -4216, 2640, 757, -1993, 977, 1492, 1913, -2018, 3456, -4876]);

        let path = vec![-977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, -3208, 2018, -1913];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010]);
    }

    #[test]
    fn test_trim_path_start_end_2() {
        let weights = hashmap!{608 => 662, 645 => 560, 693 => 461, 697 => 453, 710 => 442, 722 => 424, 884 => 228, 892 => 224, 992 => 181, 1028 => 163, 1098 => 140, 1110 => 137, 1125 => 133, 1141 => 131, 1147 => 129, 1177 => 123, 1230 => 111, 1241 => 110, 1242 => 110, 1257 => 109, 1319 => 99, 1322 => 99, 1328 => 98, 1405 => 89, 1427 => 86, 1431 => 86, 1456 => 83, 1457 => 83, 1607 => 71, 1633 => 68, 1646 => 67, 1722 => 62, 1731 => 61, 1751 => 60, 1770 => 59, 1782 => 58, 1792 => 57, 1801 => 57, 1854 => 54, 1871 => 53, 1890 => 52, 1911 => 51, 1959 => 50, 1992 => 50, 2021 => 48, 2041 => 47, 2069 => 45, 2075 => 45, 2128 => 40, 2129 => 40, 2146 => 39, 2152 => 39, 2167 => 38, 2186 => 36, 2204 => 35, 2265 => 32, 2269 => 32, 2296 => 31, 2297 => 31, 2319 => 30, 2384 => 27, 2392 => 27, 2415 => 27, 2453 => 26, 2472 => 26, 2544 => 25, 2572 => 24, 2653 => 21, 2681 => 20, 2723 => 19, 2729 => 19, 2803 => 16, 2817 => 16, 2867 => 14, 2980 => 10, 3088 => 7, 3214 => 5, 3286 => 4, 3291 => 4, 3576 => 2, 3725 => 2, 4233 => 1, 4236 => 1, 4237 => 1, 4241 => 1, 4242 => 1, 4243 => 1, 4244 => 1, 4697 => 1, 5192 => 1};

        let path = vec![-2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());

        let path = vec![-992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, -2803, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -3725, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -2472];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![-2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, -2803, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453]);

        let path = vec![-2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, 1782];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177]);
    }

    #[test]
    fn test_trim_path_start_end_3() {
        let weights = hashmap!{475 => 1316, 535 => 986, 650 => 543, 709 => 443, 734 => 404, 742 => 381, 760 => 363, 767 => 344, 820 => 275, 832 => 265, 878 => 231, 938 => 203, 941 => 202, 1186 => 121, 1190 => 120, 1265 => 108, 1347 => 96, 1368 => 94, 1384 => 92, 1434 => 85, 1473 => 81, 1576 => 73, 1680 => 64, 1920 => 51, 1947 => 50, 1958 => 50, 1970 => 50, 1989 => 50, 2009 => 49, 2082 => 44, 2128 => 40, 2191 => 36, 2233 => 34, 2269 => 32, 2292 => 31, 2301 => 31, 2342 => 29, 2396 => 27, 2431 => 26, 2535 => 25, 2545 => 25, 2582 => 24, 2585 => 24, 2590 => 23, 2593 => 23, 2595 => 23, 2596 => 23, 2602 => 23, 2626 => 22, 2638 => 22, 2641 => 21, 2665 => 21, 2796 => 16, 2814 => 16, 2835 => 15, 2908 => 12, 3075 => 7, 3122 => 6, 3219 => 5, 3223 => 5, 3312 => 4, 3328 => 3, 3351 => 3, 3379 => 3, 3405 => 3, 3483 => 2, 3564 => 2, 3722 => 2, 3808 => 2, 3928 => 1, 3929 => 1, 3930 => 1, 3931 => 1, 3932 => 1, 3933 => 1, 3935 => 1, 3974 => 1, 4217 => 1, 4218 => 1, 4333 => 1, 4334 => 1, 4335 => 1, 4336 => 1, 4337 => 1, 4338 => 1, 4363 => 1, 4559 => 1, 4560 => 1, 4574 => 1, 4637 => 1, 4638 => 1, 4640 => 1, 4958 => 1, 4959 => 1, 4983 => 1, 5177 => 1};

        let path = vec![2396, 760, -4560, 1473, -4574, 1190, 709, 938, -1970, -475, -2269, -1265, 2128, 832, 1958, -1384, -2191, -820, 3808, -1680, -2641, -4217, -2641, 535, -1186, 650, 3328, 2582, 4336, -2292, -3122, -4334, 3933, -2602, -4959, -2593, -3974, 2545, -4363, 767, -3075, -1920, 2596, -3930, 3928, -742, 4640, 1368, 4983, 941, 1947, 878, 4637, 2431, 2665, 3932, 2835, 3929, -4333, 2638, -4958, 1434, -2082, 3931, -3312, 2908, 2233, -2301, 3223, 2814, 3379, -4638, -2796, 4335, 3935, -4337, -4338, 2585, -3351, 1576, -1989, -1347, -4559, 734];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());

        let path = vec![-3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191]);

        let path = vec![-1190, -3722, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, -3722, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186];
        let trimmed_path = trim_path_start_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![-1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878]);
    }

    #[test]
    fn test_trim_path_hairpin_end_1() {
        // Tests some exact hairpin overlaps.
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10,
                               6 => 10, 7 => 10, 8 => 10, 9 => 10, 10 => 10};

        let path = vec![1, 2, 3, 4, 5];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());

        let path = vec![1, 2, 3, 4, 5, -5];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![1, 2, 3, 4, 5, -5, -4];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![1, 2, 3, 4, 5, -5, -4, -3, -2, -1];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![7, 8, 9, 10, -10, -9, -8];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![7, 8, 9, 10]);

        let path = vec![7, 8, 9, 10, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());
    }

    #[test]
    fn test_trim_path_hairpin_end_2() {
        // Tests some inexact hairpin overlaps.
        let weights = hashmap!{1 => 100, 2 => 100, 3 => 10, 4 => 100, 5 => 100,
                               6 => 1, 7 => 1, 8 => 1, 9 => 1, 10 => 1};

        let path = vec![1, 2, 3, 6, 4, 5, -5, 7, -4, -3, -2, -1];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 6, 4, 5]);

        let path = vec![1, 2, 3, 6, 4, 7, 5, -5, 8, 9, 10, -4, -3, -2, -1];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 6, 4, 7, 5]);

        let path = vec![1, 2, 3, 6, 7, 4, 8, 9, 5, -5, -4, -3, -2, 10, -1];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 6, 7, 4, 8, 9, 5]);

        let path = vec![1, 2, 3, 4, 6, -4, -3];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 6]);

        let path = vec![1, 2, 3, 4, -4, -3, 6, 7, 8, 9, 10];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4]);

        let path = vec![6, 5, 4, 3, 2, 1, -1, -2, -3, 9];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![6, 5, 4, 3, 2, 1]);
    }

    #[test]
    fn test_trim_path_hairpin_end_3() {
        // When a start overlap is so big that it could be processed as a low-identity end overlap,
        // it should be not be trimmed.
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10,
                               6 => 10, 7 => 10, 8 => 10, 9 => 10, 10 => 10};
        let path = vec![-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8];
        let trimmed_path = trim_path_hairpin_end(&path, &weights, 0.2, 1000);
        assert!(trimmed_path.is_none());
    }

    #[test]
    fn test_trim_path_hairpin_start_1() {
        // Tests some exact hairpin overlaps.
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10,
                               6 => 10, 7 => 10, 8 => 10, 9 => 10, 10 => 10};

        let path = vec![1, 2, 3, 4, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());

        let path = vec![-1, 1, 2, 3, 4, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![-2, -1, 1, 2, 3, 4, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![-5, -4, -3, -2, -1, 1, 2, 3, 4, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert!(trimmed_path.is_none());
    }

    #[test]
    fn test_trim_path_hairpin_start_2() {
        // Tests some inexact hairpin overlaps.
        let weights = hashmap!{1 => 100, 2 => 100, 3 => 10, 4 => 100, 5 => 100,
                               6 => 1, 7 => 1, 8 => 1, 9 => 1, 10 => 1};

        let path = vec![-5, 7, -4, -3, -2, -1, 1, 2, 3, 6, 4, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 6, 4, 5]);

        let path = vec![-5, 8, 9, 10, -4, -3, -2, -1, 1, 2, 3, 6, 4, 7, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 6, 4, 7, 5]);

        let path = vec![-5, -4, -3, -2, 10, -1, 1, 2, 3, 6, 7, 4, 8, 9, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 6, 7, 4, 8, 9, 5]);

        let path = vec![-2, -1, 6, 1, 2, 3, 4];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![6, 1, 2, 3, 4]);

        let path = vec![6, 7, 8, 9, 10, -2, -1, 1, 2, 3, 4];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4]);

        let path = vec![-9, 3, 2, 1, -1, -2, -3, -4, -5, -6];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![-1, -2, -3, -4, -5, -6]);
    }

    #[test]
    fn test_trim_path_hairpin_start_3() {
        // When an end overlap is so big that it could be processed as a low-identity start overlap,
        // it should be not be trimmed.
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10,
                               6 => 10, 7 => 10, 8 => 10, 9 => 10, 10 => 10};
        let path = vec![-8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.2, 1000);
        assert!(trimmed_path.is_none());
    }

    #[test]
    fn test_trim_path_hairpin_start_and_end() {
        // Tests some exact hairpin overlaps on both ends of a path.
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10,
                               6 => 10, 7 => 10, 8 => 10, 9 => 10, 10 => 10};

        let path = vec![-1, 1, 2, 3, 4, 5, -5];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        let trimmed_path = trim_path_hairpin_end(&trimmed_path.unwrap(), &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![-2, -1, 1, 2, 3, 4, 5, -5, -4];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        let trimmed_path = trim_path_hairpin_end(&trimmed_path.unwrap(), &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![-3, -2, -1, 1, 2, 3, 4, 5, -5, -4, -3];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        let trimmed_path = trim_path_hairpin_end(&trimmed_path.unwrap(), &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![-4, -3, -2, -1, 1, 2, 3, 4, 5, -5, -4, -3, -2];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        let trimmed_path = trim_path_hairpin_end(&trimmed_path.unwrap(), &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);

        let path = vec![-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, -5, -4, -3, -2, -1];
        let trimmed_path = trim_path_hairpin_start(&path, &weights, 0.95, 1000);
        let trimmed_path = trim_path_hairpin_end(&trimmed_path.unwrap(), &weights, 0.95, 1000);
        assert_eq!(trimmed_path.unwrap(), vec![1, 2, 3, 4, 5]);
    }
}
