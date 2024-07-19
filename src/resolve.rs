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

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists, reverse_path};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn resolve(cluster_dir: PathBuf) {
    let trimmed_gfa = cluster_dir.join("3_trimmed.gfa");
    let resolved_gfa = cluster_dir.join("5_resolved.gfa");
    check_settings(&cluster_dir, &trimmed_gfa);
    starting_message();
    print_settings(&cluster_dir);
    let (mut unitig_graph, sequences) = load_graph(&trimmed_gfa);
    let anchors = find_anchor_unitigs(&mut unitig_graph, &sequences);
    let bridges = create_bridges(&unitig_graph, &sequences, &anchors);

    for bridge in bridges {  // TEMP
        eprintln!("{} -> {}: {:?}", bridge.start, bridge.end, bridge.all_paths);  // TEMP
    }  // TEMP

    // TODO: save a graph with unambiguous bridges applied, keeping variants in. This graph can
    //       still contain the input sequence paths.
    
    // TODO: save a graph with unambiguous bridges applied, with each bridge reduced to its best
    //       path. This will allow for a lot of graph simplification (linear path merging), making
    //       any remaining ambiguities nice and obvious. This graph can't contain the paths.

    // TODO: gather up ambiguous bridges and cull them (in order of increasing support) until there
    //       is no more ambiguity.

    // TODO: save a graph with all bridges applied, with each bridge reduced to its best
    //       path. This will hopefully be a completed genome.

    unitig_graph.save_gfa(&resolved_gfa, &sequences).unwrap();
}


fn check_settings(cluster_dir: &PathBuf, trimmed_gfa: &PathBuf) {
    check_if_dir_exists(&cluster_dir);
    check_if_file_exists(&trimmed_gfa);
}


fn starting_message() {
    section_header("Starting autocycler resolve");
    explanation("This command resolves repeats in the unitig graph.");
}


fn print_settings(cluster_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --cluster_dir {}", cluster_dir.display());
    eprintln!();
}


fn load_graph(gfa: &PathBuf) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading graph");
    explanation("The unitig graph is now loaded into memory.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(&gfa);
    unitig_graph.print_basic_graph_info();
    (unitig_graph, sequences)
}


fn find_anchor_unitigs(graph: &mut UnitigGraph, sequences: &Vec<Sequence>) -> Vec<u32> {
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

        // TODO: I can probably remove this later for performance - it's just a sanity check.
        let mut reverse_seq_ids: Vec<_> = unitig.reverse_positions.iter().map(|p| p.seq_id()).collect();
        reverse_seq_ids.sort();
        assert!(forward_seq_ids == reverse_seq_ids);

        if forward_seq_ids == all_seq_ids {
            unitig.anchor = true;
            anchor_ids.push(unitig.number);
        }

        // TODO: dead-ends (with or without a hairpin connection) should also become anchors.
    }
    eprintln!("{} anchor unitig{} found", anchor_ids.len(), match anchor_ids.len() { 1 => "", _ => "s" });
    eprintln!();
    anchor_ids
}


fn create_bridges(graph: &UnitigGraph, sequences: &Vec<Sequence>, anchors: &Vec<u32>) -> Vec<Bridge> {
    section_header("Building bridges");
    explanation("Bridges connect one anchor unitig to the next.");
    let anchor_set: HashSet<u32> = anchors.iter().cloned().collect();
    let sequence_paths: Vec<_> = sequences.iter().map(|s| graph.get_unitig_path_for_sequence_i32(s)).collect();
    let anchor_to_anchor_paths = get_anchor_to_anchor_paths(&sequence_paths, &anchor_set);
    let grouped_paths = group_paths_by_start_end(anchor_to_anchor_paths);
    let mut bridges = Vec::new();
    for ((start, end), paths) in grouped_paths {
        bridges.push(Bridge::new(start, end, paths));
    }
    bridges
}


fn get_anchor_to_anchor_paths(sequence_paths: &Vec<Vec<i32>>, anchor_set: &HashSet<u32>) -> Vec<Vec<i32>> {
    let mut anchor_to_anchor_paths = Vec::new();
    for path in sequence_paths {
        let mut last_anchor_i: Option<usize> = None;
        for (i, &value) in path.iter().enumerate() {
            if anchor_set.contains(&(value.abs() as u32)) {
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


fn group_paths_by_start_end(anchor_to_anchor_paths: Vec<Vec<i32>>) -> HashMap<(i32, i32), Vec<Vec<i32>>> {
    let mut grouped_paths: HashMap<(i32, i32), Vec<Vec<i32>>> = HashMap::new();
    for path in anchor_to_anchor_paths {
        if let (Some(&start), Some(&end)) = (path.first(), path.last()) {
            grouped_paths.entry((start, end)).or_insert_with(Vec::new).push(path);
        }
    }
    grouped_paths
}


pub struct Bridge {
    start: i32,
    end: i32,
    all_paths: Vec<Vec<i32>>,
    best_path: Vec<i32>,
}

impl Bridge {
    pub fn new(start: i32, end: i32, all_paths: Vec<Vec<i32>>) -> Self {

        // TODO: trim start/end off all_paths (should be same as start and end)

        // TODO: set a best path

        Bridge {
            start,
            end,
            all_paths,
            best_path: vec![],  // TEMP
        }
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
}
