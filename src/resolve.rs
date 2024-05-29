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

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::cluster::load_graph;
use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists, quit_with_error};
use crate::sequence::Sequence;
use crate::trim_path_overlap::trim_path_overlap;
use crate::unitig_graph::UnitigGraph;


pub fn resolve(out_dir: PathBuf) {
    let gfa = out_dir.join("01_unitig_graph.gfa");
    let clusters = out_dir.join("05_clusters.tsv");
    let clean_gfa = out_dir.join("06_cleaned_graph.gfa");
    check_settings(&out_dir, &gfa, &clusters);
    starting_message();
    print_settings(&out_dir);
    let (mut unitig_graph, mut sequences) = load_graph(&gfa);
    load_clusters(&clusters, &mut sequences);
    let sequences = remove_excluded_contigs_from_graph(&mut unitig_graph, &sequences);
    let sequences = trim_path_overlap(&mut unitig_graph, &sequences, 0.95);  // TODO: make overlap alignment identity an argument
    // TODO: trim over-long circular contigs (unless suppressed with an argument)
    unitig_graph.save_gfa(&clean_gfa, &sequences).unwrap();
    // TODO: find initial single-copy contigs
    // TODO: expand set of single-copy contigs based on graph structure
    // TODO: find the ordering of single-copy contigs
    // TODO: resolve repeats by duplicating unitigs
}


fn check_settings(out_dir: &PathBuf, gfa: &PathBuf, clusters: &PathBuf) {
    check_if_dir_exists(&out_dir);
    check_if_file_exists(&gfa);
    check_if_file_exists(&clusters);
}


fn starting_message() {
    section_header("Starting autocycler resolve");
    explanation("This command resolves repeats in the unitig graph.");
}


fn print_settings(out_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!();
}


fn load_clusters(clusters: &PathBuf, sequences: &mut Vec<Sequence>) {
    let file = File::open(clusters).unwrap();
    let reader = BufReader::new(file);
    let mut header = true;
    let mut cluster_nums: HashMap<(String, String), i32> = HashMap::new();
    for line_result in reader.lines() {
        let line = line_result.unwrap();
        let parts: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if header {
            check_cluster_header(&parts, clusters);
            header = false; continue;
        }
        let (assembly, contig_name, cluster) = get_cluster_line_parts(&parts, clusters);
        cluster_nums.insert((assembly, contig_name), cluster);
    }
    for s in sequences {
        if let Some(&cluster_num) = cluster_nums.get(&(s.filename.clone(), s.contig_name())) {
            s.cluster = cluster_num;
        } else {
            quit_with_error(&format!("missing cluster for {}", s.filename));
        }
    }
}


fn check_cluster_header(parts: &Vec<&str>, clusters: &PathBuf) {
    if parts != &vec!["assembly", "contig_name", "length", "cluster"] {
        quit_with_error(&format!("cluster file ({}) does not contain expected header", clusters.display()));
    }
}


fn get_cluster_line_parts(parts: &Vec<&str>, clusters: &PathBuf) -> (String, String, i32) {
    if parts.len() != 4 {
        quit_with_error(&format!("cluster file ({}) contains an unexpected number of columns", clusters.display()));
    }
    let assembly = parts[0].to_string();
    let contig_name = parts[1].to_string();
    let cluster_str = parts[3];
    if cluster_str == "none" {
        return (assembly, contig_name, -1);
    }
    if let Ok(num) = cluster_str.parse::<i32>() {
        if num > 0 {
            return (assembly, contig_name, num);
        }
    }
    quit_with_error(&format!("cluster file ({}) contains an invalid cluster number: {}",
                    clusters.display(), cluster_str));
    unreachable!()
}


pub fn remove_excluded_contigs_from_graph(graph: &mut UnitigGraph, sequences: &Vec<Sequence>) -> Vec<Sequence> {
    section_header("Cleaning graph");
    explanation("Excluded contigs (those which could not be clustered) are now removed from the \
                 unitig graph.");
    let seqs_to_remove: Vec<_> = sequences.iter().filter(|s| s.cluster == -1).collect();
    if seqs_to_remove.len() == 0 {
        eprintln!("No contigs needed to be removed");
    } else {
        eprintln!("Removed contigs:");
        for s in seqs_to_remove {
            eprintln!("  {}", s);
            remove_contig_from_graph(graph, s.id);
        }
    }
    graph.remove_zero_depth_unitigs();
    let sequences = sequences.iter().filter(|s| s.cluster != -1).cloned().collect();
    eprintln!();
    merge_linear_paths(graph, &sequences);
    graph.print_basic_graph_info();
    graph.renumber_unitigs();
    sequences
}


pub fn remove_contig_from_graph(graph: &mut UnitigGraph, seq_id: u16) {
    // Removes all Positions from the Unitigs which have the given sequence ID. This reduces
    // depths of affected Unitigs, and Unitigs which have their depths reduced to zero are
    // removed from the graph.
    for u in &graph.unitigs {
        u.borrow_mut().remove_sequence(seq_id);
    }
}
