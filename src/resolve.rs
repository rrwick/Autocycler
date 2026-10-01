// This file contains the code for the autocycler resolve subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use colored::Colorize;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists, reverse_path, load_file_lines,
                  sign_at_end, sign_at_end_vec};
use crate::sequence::Sequence;
use crate::unitig::{Unitig, UnitigType};
use crate::unitig_graph::UnitigGraph;


pub fn resolve(cluster_dir: PathBuf, verbose: bool) {
    let trimmed_gfa = cluster_dir.join("2_trimmed.gfa");
    let bridged_gfa = cluster_dir.join("3_bridged.gfa");
    let merged_gfa = cluster_dir.join("4_merged.gfa");
    let final_gfa = cluster_dir.join("5_final.gfa");

    check_settings(&cluster_dir, &trimmed_gfa);
    starting_message();
    print_settings(&cluster_dir, verbose);

    let gfa_lines = load_file_lines(&trimmed_gfa);
    let (mut unitig_graph, sequences) = load_graph(&gfa_lines, true, None);

    let anchors = find_anchor_unitigs(&mut unitig_graph, &sequences);
    let mut bridges = create_bridges(&unitig_graph, &sequences, &anchors, verbose);
    let bridge_count = bridges.len();
    let bridge_depth = sequences.len() as f64;
    determine_ambiguity(&mut bridges);
    print_bridges(&bridges, verbose);

    apply_unique_message();
    apply_bridges(&mut unitig_graph, &bridges, bridge_depth);
    unitig_graph.save_gfa(&bridged_gfa, &[], false).unwrap();
    merge_after_bridging(&mut unitig_graph);
    unitig_graph.save_gfa(&merged_gfa, &[], false).unwrap();

    let cull_count = cull_ambiguity(&mut bridges, verbose);
    if cull_count > 0 {
        (unitig_graph, _) = load_graph(&gfa_lines, false, Some(&anchors));
        apply_final_message();
        apply_bridges(&mut unitig_graph, &bridges, bridge_depth);
        merge_after_bridging(&mut unitig_graph);
    } else if bridge_count > 0 {
        eprintln!("All bridges were unique, no culling necessary.\n");
    }

    unitig_graph.save_gfa(&final_gfa, &[], true).unwrap();
    finished_message(&final_gfa);
}


fn check_settings(cluster_dir: &Path, trimmed_gfa: &Path) {
    check_if_dir_exists(cluster_dir);
    check_if_file_exists(trimmed_gfa);
}


fn starting_message() {
    section_header("Starting autocycler resolve");
    explanation("This command resolves repeats in the unitig graph.");
}


fn apply_unique_message() {
    section_header("Applying unique bridges");
    explanation("All unique bridges (those that do not conflict with other bridges) are now \
                 applied to the graph, with linear paths merged to create consentigs.");
}

fn apply_final_message() {
    section_header("Applying final bridges");
    explanation("Now that conflicting bridges have been removed, bridges are applied one more \
                 time to create the final graph.");
}


fn finished_message(final_gfa: &Path) {
    section_header("Finished!");
    eprintln!("Final consensus graph: {}", final_gfa.display());
    eprintln!();
}


fn print_settings(cluster_dir: &Path, verbose: bool) {
    eprintln!("Settings:");
    eprintln!("  --cluster_dir {}", cluster_dir.display());
    if verbose {
        eprintln!("  --verbose");
    }
    eprintln!();
}


fn load_graph(gfa_lines: &[String], print_info: bool,
              anchors: Option<&[u32]>) -> (UnitigGraph, Vec<Sequence>) {
    if print_info {
        section_header("Loading graph");
        explanation("The unitig graph is now loaded into memory.");
    }
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_lines(gfa_lines);
    if let Some(anchors) = anchors {
        for num in anchors {
            unitig_graph.unitig_index.get(num).unwrap()
                        .borrow_mut().unitig_type = UnitigType::Anchor;
        }
    }
    if print_info {
        unitig_graph.print_basic_graph_info();
    }
    (unitig_graph, sequences)
}


fn find_anchor_unitigs(graph: &mut UnitigGraph, sequences: &[Sequence]) -> Vec<u32> {
    section_header("Finding anchor unitigs");
    explanation("Anchor unitigs are those that occur once and only once in each sequence. They \
                 will definitely be present in the final sequence and will serve as the connection \
                 points for bridges.");
    let mut all_seq_ids: Vec<_> = sequences.iter().map(|s| s.id).collect();
    all_seq_ids.sort();
    let mut anchor_ids = Vec::new();
    for unitig_rc in &graph.unitigs {
        let mut unitig = unitig_rc.borrow_mut();
        let mut forward_seq_ids: Vec<_> = unitig.forward_positions.iter().map(|p| p.seq_id()).collect();
        forward_seq_ids.sort();
        if forward_seq_ids == all_seq_ids {
            unitig.unitig_type = UnitigType::Anchor;
            anchor_ids.push(unitig.number);
        }
    }

    // TODO: add additional logic to better handle linear replicons?
    //       Specifically, I think it could be good to extend the set of anchor unitigs by adding
    //       unitigs which are definitely in the same path as an anchor. E.g. if unitig A in an
    //       anchor, and all of unitig A's paths lead to unitig B and vice versa, then unitig B
    //       should also be an anchor. Do this iteratively until no more anchors are added.

    // TODO: I should also allow users to manually specify anchor unitigs via their ID.

    eprintln!("{} anchor unitig{} found", anchor_ids.len(), match anchor_ids.len() { 1 => "", _ => "s" });
    eprintln!();
    anchor_ids
}


fn create_bridges(graph: &UnitigGraph, sequences: &[Sequence], anchors: &[u32], verbose: bool)
        -> Vec<Bridge> {
    section_header("Building bridges");
    explanation("Bridges connect one anchor unitig to the next.");
    let anchor_set: HashSet<u32> = anchors.iter().copied().collect();
    let mut sequence_paths = Vec::new();
    for sequence in sequences {
        let weight = sequence.consensus_weight();
        if verbose { eprintln!("{sequence} consensus weight = {weight}"); }
        if weight > 0 {
            let path = graph.get_unitig_path_for_sequence_i32(sequence);
            sequence_paths.extend(std::iter::repeat_n(path, weight));
        }
    }
    if verbose { eprintln!(); }

    let anchor_to_anchor_paths = get_anchor_to_anchor_paths(&sequence_paths, &anchor_set);
    let grouped_paths = group_paths_by_start_end(anchor_to_anchor_paths);
    let unitig_lengths: HashMap<_, _> = graph.unitigs.iter().map(|rc| {
        let unitig = rc.borrow();
        (unitig.number as i32, unitig.length())
    }).collect();
    let mut bridges: Vec<_> = grouped_paths.into_iter()
        .map(|((start, end), paths)| Bridge::new(start, end, paths, &unitig_lengths)).collect();
    bridges.sort();
    bridges
}


fn determine_ambiguity(bridges: &mut [Bridge]) {
    // Counting starts in both orientations also detects conflicts at ends.
    let mut start_count = HashMap::new();
    for bridge in bridges.iter() {
        *start_count.entry(bridge.start).or_insert(0) += 1;
        *start_count.entry(bridge.rev_start()).or_insert(0) += 1;
    }
    for bridge in bridges {
        bridge.conflicting = start_count[&bridge.start] > 1 || start_count[&bridge.rev_start()] > 1;
    }
}


fn apply_bridges(graph: &mut UnitigGraph, bridges: &[Bridge], bridge_depth: f64) {
    graph.clear_positions();
    for bridge in bridges.iter().filter(|b| !b.conflicting) {
        graph.delete_outgoing_links(bridge.start);
        graph.delete_incoming_links(bridge.end);

        if bridge.best_path.is_empty() {
            graph.create_link(bridge.start, bridge.end);
        } else {
            let bridge_seq = graph.get_sequence_from_path_signed(&bridge.best_path);
            let bridge_num = graph.max_unitig_number() + 1;
            let bridge_unitig = Unitig::bridge(bridge_num, bridge_seq, bridge_depth);
            let bridge_unitig_rc = Rc::new(RefCell::new(bridge_unitig));
            graph.unitigs.push(bridge_unitig_rc.clone());
            graph.unitig_index.insert(bridge_num, bridge_unitig_rc);
            reduce_depths(graph, bridge);
            graph.create_link(bridge.start, bridge_num as i32);
            graph.create_link(bridge_num as i32, bridge.end)
        }
    }
    delete_unitigs_not_connected_to_anchor(graph);
    graph.remove_zero_depth_unitigs();

    // TODO: add logic for removing non-anchor tips to handle open ends?
}


fn merge_after_bridging(graph: &mut UnitigGraph) {
    merge_linear_paths(graph, &[]);
    graph.print_basic_graph_info();
    graph.renumber_unitigs();
}


fn reduce_depths(graph: &mut UnitigGraph, bridge: &Bridge) {
    for signed_num in bridge.all_paths.iter().flatten() {
        let mut unitig = graph.unitig_index[&signed_num.unsigned_abs()].borrow_mut();
        unitig.reduce_depth_by_one();
    }
}


fn delete_unitigs_not_connected_to_anchor(graph: &mut UnitigGraph) {
    let to_delete: HashSet<u32> = graph.connected_components().into_iter()
        .filter(|component| component.iter().all(|num|
            graph.unitig_index[num].borrow().unitig_type != UnitigType::Anchor))
        .flatten().collect();
    graph.remove_unitigs_by_number(to_delete);
}


fn cull_ambiguity(bridges: &mut Vec<Bridge>, verbose: bool) -> usize {
    if !bridges.iter().any(|b| b.conflicting) {
        return 0;
    }
    section_header("Culling conflicting bridges");
    explanation("The least-supported conflicting bridges are now culled until no bridges \
                 conflict.");
    let mut cull_count = 0;
    if verbose {
        eprintln!("Culled bridges:");
    }
    while let Some(to_cull) = bridges.iter().filter(|b| b.conflicting)
        .min_by_key(|b| (b.depth(), *b)) {
        if verbose {
            eprintln!("  {to_cull}");
        }
        let index = bridges.iter().position(|b| (b.start, b.end) == (to_cull.start, to_cull.end)).unwrap();
        bridges.remove(index);
        cull_count += 1;
        determine_ambiguity(bridges);
    }
    if verbose { eprintln!(); }
    eprintln!("{} conflicting bridge{} culled", cull_count, match cull_count { 1 => "", _ => "s" });
    eprintln!();
    cull_count
}


fn print_bridges(bridges: &[Bridge], verbose: bool) {
    let unique_count = bridges.iter().filter(|b| !b.conflicting).count();
    let conflicting_count = bridges.len() - unique_count;
    if verbose {
        if unique_count > 0 {
            eprintln!("Unique bridges:");
            for bridge in bridges.iter().filter(|b| !b.conflicting) {
                eprintln!("  {bridge}");
            }
        }
        if conflicting_count > 0 {
            eprintln!("\nConflicting bridges:");
            for bridge in bridges.iter().filter(|b| b.conflicting) {
                eprintln!("  {bridge}");
            }
        }
    } else {
        eprintln!("     Unique bridges: {unique_count}");
        eprintln!("Conflicting bridges: {conflicting_count}");
    }
    eprintln!();
}


fn get_anchor_to_anchor_paths(sequence_paths: &[Vec<i32>], anchor_set: &HashSet<u32>)
        -> Vec<Vec<i32>> {
    let mut anchor_to_anchor_paths = Vec::new();
    for path in sequence_paths {
        let mut last_anchor_i: Option<usize> = None;
        for (i, &value) in path.iter().enumerate() {
            if anchor_set.contains(&value.unsigned_abs()) {
                if let Some(start) = last_anchor_i {
                    let a_to_a_forward = &path[start..=i];
                    let a_to_a_reverse = reverse_path(a_to_a_forward);
                    if a_to_a_forward > &a_to_a_reverse {
                        anchor_to_anchor_paths.push(a_to_a_forward.to_vec());
                    } else {
                        anchor_to_anchor_paths.push(a_to_a_reverse);
                    }
                }
                last_anchor_i = Some(i);
            }
        }
    }
    anchor_to_anchor_paths
}


fn group_paths_by_start_end(anchor_to_anchor_paths: Vec<Vec<i32>>)
        -> HashMap<(i32, i32), Vec<Vec<i32>>> {
    let mut grouped_paths: HashMap<(i32, i32), Vec<Vec<i32>>> = HashMap::new();
    for path in anchor_to_anchor_paths {
        if let (Some(&start), Some(&end)) = (path.first(), path.last()) {
            grouped_paths.entry((start, end)).or_default().push(path);
        }
    }
    grouped_paths
}


fn choose_best_path(paths: &[Vec<i32>], unitig_lengths: &HashMap<i32, u32>) -> Vec<i32> {
    let mut best_path: &[i32] = &[];
    let mut best_total = u32::MAX;
    for path in paths {
        let total = paths.iter().filter(|other| *other != path)
            .map(|other| global_alignment_distance(path, other, unitig_lengths)).sum();
        if total < best_total || (total == best_total && path.as_slice() < best_path) {
            best_total = total;
            best_path = path;
        }
    }
    best_path.to_vec()
}


fn global_alignment_distance(path_a: &[i32], path_b: &[i32], weights: &HashMap<i32, u32>) -> u32 {
    // Unitig lengths weight gaps and mismatches in this two-row dynamic-programming alignment.
    let weighted_b: Vec<_> = path_b.iter().map(|&num| (num, weights[&num.abs()])).collect();
    let mut prev = vec![0u32; path_b.len() + 1];
    let mut curr = vec![0u32; path_b.len() + 1];

    for (j, &(_, weight_b)) in weighted_b.iter().enumerate() {
        prev[j + 1] = prev[j] + weight_b;
    }

    for &unitig_a in path_a {
        let weight_a = weights[&unitig_a.abs()];
        curr[0] = prev[0] + weight_a;
        for (j, &(unitig_b, weight_b)) in weighted_b.iter().enumerate() {
            let substitution = if unitig_a == unitig_b { 0 } else { weight_a.max(weight_b) };
            curr[j + 1] = (prev[j] + substitution)
                .min(prev[j + 1] + weight_a).min(curr[j] + weight_b);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[path_b.len()]
}


pub struct Bridge {
    start: i32,
    end: i32,
    all_paths: Vec<Vec<i32>>,
    best_path: Vec<i32>,
    conflicting: bool,
}

impl Bridge {
    fn new(start: i32, end: i32, mut all_paths: Vec<Vec<i32>>, unitig_lengths: &HashMap<i32, u32>)
            -> Self {
        // Remove the start and end unitigs from the paths.
        for path in &mut all_paths {
            path.remove(0);
            path.pop();
        }

        Bridge {
            start,
            end,
            best_path: choose_best_path(&all_paths, unitig_lengths),
            all_paths,
            conflicting: false,
        }
    }

    fn rev_start(&self) -> i32 {
        -self.end
    }

    fn depth(&self) -> usize {
        self.all_paths.len()
    }
}

impl fmt::Display for Bridge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.best_path.is_empty() {
            write!(f, "{} → {} ({}×)", sign_at_end(self.start),
                   sign_at_end(self.end), self.depth())
        } else {
            write!(f, "{} → {} → {} ({}×)", sign_at_end(self.start),
                   sign_at_end_vec(&self.best_path).dimmed(), sign_at_end(self.end), self.depth())
        }
    }
}

impl PartialEq for Bridge {
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start &&
        self.end == other.end &&
        self.best_path == other.best_path
    }
}

impl Eq for Bridge {}

impl PartialOrd for Bridge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Bridge {
    fn cmp(&self, other: &Self) -> Ordering {
        self.start.abs().cmp(&other.start.abs())
            .then_with(|| other.start.cmp(&self.start))
            .then_with(|| self.end.abs().cmp(&other.end.abs()))
            .then_with(|| other.end.cmp(&self.end))
            .then_with(|| self.best_path.cmp(&other.best_path))
    }
}


#[cfg(test)]
mod tests {
    use maplit::hashmap;
    use super::*;

    #[test]
    fn test_get_anchor_to_anchor_paths() {
        let sequence_paths = vec![vec![1, -10, 4, 6, -5, -2, -9, 3, 8, -7],
                                  vec![-2, -9, 12, 8, -7, 1, -10, 4, 6, -5],
                                  vec![7, -8, -3, 9, 2, 11, -6, -4, 10, -1]];
        let anchor_set: HashSet<u32> = HashSet::from([1, 2, 6, 8]);
        let anchor_to_anchor_paths = get_anchor_to_anchor_paths(&sequence_paths, &anchor_set);
        assert_eq!(anchor_to_anchor_paths,
                   vec![vec![1, -10, 4, 6], vec![6, -5, -2], vec![-2, -9, 3, 8],
                        vec![-2, -9, 12, 8], vec![8, -7, 1], vec![1, -10, 4, 6],
                        vec![-2, -9, 3, 8], vec![6, -11, -2], vec![1, -10, 4, 6]]);
    }

    #[test]
    fn test_group_paths_by_start_end() {
        let anchor_to_anchor_paths = vec![vec![1, -10, 4, 6], vec![6, -5, -2], vec![-2, -9, 3, 8],
                                          vec![-2, -9, 12, 8], vec![8, -7, 1], vec![1, -10, 4, 6],
                                          vec![-2, -9, 3, 8], vec![6, -11, -2], vec![1, -10, 4, 6]];
        let grouped_paths = group_paths_by_start_end(anchor_to_anchor_paths);

        assert_eq!(grouped_paths,
                   hashmap!{(1, 6) => vec![vec![1, -10, 4, 6], vec![1, -10, 4, 6], vec![1, -10, 4, 6]],
                            (6, -2) => vec![vec![6, -5, -2], vec![6, -11, -2]],
                            (-2, 8) => vec![vec![-2, -9, 3, 8], vec![-2, -9, 12, 8], vec![-2, -9, 3, 8]],
                            (8, 1) => vec![vec![8, -7, 1]]});
    }

    #[test]
    fn test_bridge_unitig_nums() {
        let unitig_lengths = hashmap!{1 => 10, 12 => 10, 23 => 10, 8 => 10, 41 => 10, 2 => 10,
                                      17 => 10, 123 => 10};
        let paths = vec![vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, 17, 123, 41, 2]];
        let bridge = Bridge::new(1, 2, paths, &unitig_lengths);
        assert_eq!(bridge.rev_start(), -2);
        assert_eq!(bridge.depth(), 4);
    }

    #[test]
    fn test_determine_ambiguity_1() {
        let unitig_lengths = hashmap!{1 => 10, 2 => 10, 4 => 10, 5 => 10,
                                      6 => 10, 11 => 10, 12 => 10};
        let bridge_a = Bridge::new(1, -2, vec![vec![1, 12, 2]], &unitig_lengths);
        let bridge_b = Bridge::new(-2, 5, vec![vec![-2, 6, 5]], &unitig_lengths);
        let bridge_c = Bridge::new(4, -5, vec![vec![4, -5]], &unitig_lengths);
        let bridge_d = Bridge::new(-4, 6, vec![vec![-4, 12, 6]], &unitig_lengths);
        let bridge_e = Bridge::new(-1, -6, vec![vec![-1, 11, -6]], &unitig_lengths);
        let mut bridges = vec![bridge_a, bridge_b, bridge_c, bridge_d, bridge_e];
        determine_ambiguity(&mut bridges);
        assert!(!bridges[0].conflicting);
        assert!(!bridges[1].conflicting);
        assert!(!bridges[2].conflicting);
        assert!(!bridges[3].conflicting);
        assert!(!bridges[4].conflicting);
    }

    #[test]
    fn test_determine_ambiguity_2() {
        let unitig_lengths = hashmap!{1 => 10, 2 => 10, 4 => 10, 5 => 10, 6 => 10, 7 => 10,
                                      8 => 10, 9 => 10, 11 => 10, 12 => 10, 13 => 10, 14 => 10};
        let bridge_a = Bridge::new(1, -2, vec![vec![1, 12, 2]], &unitig_lengths);
        let bridge_b = Bridge::new(-2, 5, vec![vec![-2, 6, 5]], &unitig_lengths);
        let bridge_c = Bridge::new(4, -5, vec![vec![4, -5]], &unitig_lengths);
        let bridge_d = Bridge::new(-4, 6, vec![vec![-4, 12, 6]], &unitig_lengths);
        let bridge_e = Bridge::new(-1, -6, vec![vec![-1, 11, -6]], &unitig_lengths);
        let bridge_f = Bridge::new(-4, 7, vec![vec![-4, 13, 7]], &unitig_lengths);
        let bridge_g = Bridge::new(1, 8, vec![vec![1, 14, 8]], &unitig_lengths);
        let bridge_h = Bridge::new(4, -8, vec![vec![4, 9, -8]], &unitig_lengths);
        let mut bridges = vec![bridge_a, bridge_b, bridge_c, bridge_d,
                               bridge_e, bridge_f, bridge_g, bridge_h];
        determine_ambiguity(&mut bridges);
        assert!(bridges[0].conflicting);
        assert!(!bridges[1].conflicting);
        assert!(bridges[2].conflicting);
        assert!(bridges[3].conflicting);
        assert!(!bridges[4].conflicting);
        assert!(bridges[5].conflicting);
        assert!(bridges[6].conflicting);
        assert!(bridges[7].conflicting);
    }

    #[test]
    fn test_cull_ambiguity() {
        let lengths = HashMap::new();
        let mut bridges = vec![
            Bridge::new(2, 3, vec![vec![2, 3]; 2], &lengths),
            Bridge::new(1, 3, vec![vec![1, 3]; 2], &lengths),
            Bridge::new(1, 4, vec![vec![1, 4]], &lengths),
            Bridge::new(-2, 5, vec![vec![-2, 5]], &lengths),
        ];
        determine_ambiguity(&mut bridges);
        assert_eq!(cull_ambiguity(&mut bridges, false), 2);
        assert_eq!(bridges.iter().map(|b| (b.start, b.end)).collect::<Vec<_>>(),
                   vec![(2, 3), (-2, 5)]);
        assert!(bridges.iter().all(|b| !b.conflicting));
        assert_eq!(cull_ambiguity(&mut bridges, true), 0);
        assert_eq!(cull_ambiguity(&mut Vec::new(), false), 0);
    }

    #[test]
    fn test_hairpin_and_circular_bridge_conflicts() {
        let lengths = HashMap::new();
        let mut bridges = vec![Bridge::new(1, -1, vec![vec![1, -1]], &lengths),
                               Bridge::new(2, 2, vec![vec![2, 2]], &lengths)];
        determine_ambiguity(&mut bridges);
        assert!(bridges[0].conflicting);
        assert!(!bridges[1].conflicting);
        assert_eq!(cull_ambiguity(&mut bridges, true), 1);
        assert_eq!((bridges[0].start, bridges[0].end), (2, 2));
    }

    #[test]
    fn test_best_path_with_empty_alternative() {
        let bridge = Bridge::new(1, 2, vec![vec![1, 3, 2], vec![1, 2]], &hashmap!{3 => 10});
        assert!(bridge.best_path.is_empty());
        assert_eq!(bridge.all_paths, vec![vec![3], vec![]]);
        assert_eq!(bridge.depth(), 2);
    }

    #[test]
    fn test_best_path_1() {
        // Easy case: 3 vs 1
        let unitig_lengths = hashmap!{1 => 10, 12 => 10, 23 => 10, 8 => 10, 41 => 10, 2 => 10,
                                      17 => 10, 123 => 10};

        let paths = vec![vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, 17, 123, 41, 2]];
        let bridge = Bridge::new(1, 2, paths, &unitig_lengths);
        assert_eq!(bridge.best_path, vec![12, -23, -8, 41]);
    }

    #[test]
    fn test_best_path_2() {
        // 2 vs 2 tie, go with the lexographically smaller path.
        let unitig_lengths = hashmap!{1 => 10, 12 => 10, 23 => 10, 8 => 10, 41 => 10, 2 => 10,
                                      17 => 10, 123 => 10};

        let paths = vec![vec![1, 12, 17, 123, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, 17, 123, 41, 2]];
        let bridge = Bridge::new(1, 2, paths, &unitig_lengths);
        assert_eq!(bridge.best_path, vec![12, -23, -8, 41]);
    }

    #[test]
    fn test_best_path_3() {
        // Tricky case: the most common path is the not the best. Even though 1,13,12 appears twice,
        // the best path should be 1,2,3,4,5,6,7,8,9,10,11,12 since it minimises the total distance
        // to all paths.
        let unitig_lengths = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10, 6 => 10,
                                      7 => 10, 8 => 10, 9 => 10, 10 => 10, 11 => 10, 12 => 10,
                                      13 => 10, 14 => 10, 15 => 10, 16 => 10, 17 => 10, 18 => 10,
                                      19 => 10, 20 => 10, 21 => 10};

        let paths = vec![vec![1, 2, 3, 4, 5, 6, 7, 8, 20, 10, 11, 12],
                         vec![1, 13, 12],
                         vec![1, 2, 3, 4, 16, 6, 7, 8, 9, 10, 11, 12],
                         vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 21, 11, 12],
                         vec![1, 2, 3, 4, 5, 17, 7, 8, 9, 10, 11, 12],
                         vec![1, 13, 12],
                         vec![1, 2, 3, 4, 5, 6, 18, 8, 9, 10, 11, 12],
                         vec![1, 2, 14, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                         vec![1, 2, 3, 15, 5, 6, 7, 8, 9, 10, 11, 12],
                         vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                         vec![1, 2, 3, 4, 5, 6, 7, 19, 9, 10, 11, 12]];
        let bridge = Bridge::new(1, 2, paths, &unitig_lengths);
        assert_eq!(bridge.best_path, vec![2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
    }

    #[test]
    fn test_global_alignment_distance_1() {
        // Two identical paths have a distance of 0.
        let unitig_lengths = hashmap!{1 => 10, 2 => 1, 3 => 2, 4 => 3, 5 => 4, 6 => 10};
        let a = vec![1, 2, 3, 4, 5, 6];
        let b = vec![1, 2, 3, 4, 5, 6];
        assert_eq!(global_alignment_distance(a.as_slice(), b.as_slice(), &unitig_lengths), 0);
        let a = vec![];
        let b = vec![];
        assert_eq!(global_alignment_distance(a.as_slice(), b.as_slice(), &unitig_lengths), 0);
    }

    #[test]
    fn test_global_alignment_distance_2() {
        // Gaps create distance equal to their unitig lengths.
        let unitig_lengths = hashmap!{1 => 10, 2 => 1, 3 => 2, 4 => 3, 5 => 4, 6 => 10};
        let a = vec![1, 2, 3, 4, 5, 6];
        let b = vec![1, 2, 3, 4, 6];
        assert_eq!(global_alignment_distance(a.as_slice(), b.as_slice(), &unitig_lengths), 4);
        assert_eq!(global_alignment_distance(b.as_slice(), a.as_slice(), &unitig_lengths), 4);
        let a = vec![1, 2, 4, 5, 6];
        let b = vec![1, 2, 3, 4, 5, 6];
        assert_eq!(global_alignment_distance(a.as_slice(), b.as_slice(), &unitig_lengths), 2);
        assert_eq!(global_alignment_distance(b.as_slice(), a.as_slice(), &unitig_lengths), 2);
        let a = vec![1, 3, 4, 5, 6];
        let b = vec![1, 2, 3, 5, 6];
        assert_eq!(global_alignment_distance(a.as_slice(), b.as_slice(), &unitig_lengths), 4);
        assert_eq!(global_alignment_distance(b.as_slice(), a.as_slice(), &unitig_lengths), 4);
        let a = vec![1, 2, 3, 4, 5, 6];
        let b = vec![];
        assert_eq!(global_alignment_distance(a.as_slice(), b.as_slice(), &unitig_lengths), 30);
        assert_eq!(global_alignment_distance(b.as_slice(), a.as_slice(), &unitig_lengths), 30);
    }

    #[test]
    fn test_global_alignment_distance_3() {
        // Mismatches create distance equal to the longer unitig length.
        let unitig_lengths = hashmap!{1 => 10, 2 => 1, 3 => 2, 4 => 3, 5 => 4, 6 => 10};
        let a = vec![1, 2, 3, 5, 6];
        let b = vec![1, 2, 4, 5, 6];
        assert_eq!(global_alignment_distance(a.as_slice(), b.as_slice(), &unitig_lengths), 3);
        assert_eq!(global_alignment_distance(b.as_slice(), a.as_slice(), &unitig_lengths), 3);
    }
}
