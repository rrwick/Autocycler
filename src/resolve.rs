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
    let mut bridges = create_bridges(&unitig_graph, &sequences, &anchors);
    let bridge_depth = sequences.len() as f64;
    determine_ambiguity(&mut bridges);
    print_bridges(&bridges, verbose);

    apply_unique_message();
    apply_bridges(&mut unitig_graph, &bridges, bridge_depth);
    unitig_graph.save_gfa(&bridged_gfa, &vec![], false).unwrap();
    merge_after_bridging(&mut unitig_graph);
    unitig_graph.save_gfa(&merged_gfa, &vec![], false).unwrap();

    let cull_count = cull_ambiguity(&mut bridges, verbose);
    if cull_count > 0 {
        (unitig_graph, _) = load_graph(&gfa_lines, false, Some(&anchors));
        apply_final_message();
        apply_bridges(&mut unitig_graph, &bridges, bridge_depth);
        merge_after_bridging(&mut unitig_graph);
    } else {
        eprintln!("All bridges were unique, no culling necessary.\n");
    }
    unitig_graph.save_gfa(&final_gfa, &vec![], true).unwrap();
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


fn load_graph(gfa_lines: &Vec<String>, print_info: bool,
              anchors: Option<&Vec<u32>>) -> (UnitigGraph, Vec<Sequence>) {
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


fn create_bridges(graph: &UnitigGraph, sequences: &[Sequence], anchors: &[u32]) -> Vec<Bridge> {
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
    bridges.sort();
    bridges
}


fn determine_ambiguity(bridges: &mut [Bridge]) -> usize {
    // This function classifies each Bridge as conflicting or not. A Bridge is conflicting if it
    // shares its start or end unitig with another Bridge. The return value is the number of
    // conflicting Bridges.
    let mut start_count = HashMap::new();
    for bridge in bridges.iter() {
        *start_count.entry(bridge.start).or_insert(0) += 1;
        *start_count.entry(bridge.rev_start()).or_insert(0) += 1;
    }
    let mut end_count = HashMap::new();
    for bridge in bridges.iter() {
        *end_count.entry(bridge.end).or_insert(0) += 1;
        *end_count.entry(bridge.rev_end()).or_insert(0) += 1;
    }
    let mut ambi_count = 0;
    let ambi_starts: HashSet<_> = start_count.iter().filter(|&(_, &c)| c > 1).map(|(&n, _)| n).collect();
    let ambi_ends: HashSet<_> = end_count.iter().filter(|&(_, &c)| c > 1).map(|(&n, _)| n).collect();
    for bridge in bridges.iter_mut() {
        if ambi_starts.contains(&bridge.start) || ambi_starts.contains(&bridge.rev_start()) ||
           ambi_ends.contains(&bridge.end) || ambi_ends.contains(&bridge.rev_end()) {
            bridge.conflicting = true;
            ambi_count += 1;
        } else {
            bridge.conflicting = false;
        }
    }
    ambi_count
}


fn apply_bridges(graph: &mut UnitigGraph, bridges: &Vec<Bridge>, bridge_depth: f64) {
    // This function applies bridges to the graph. If keep_all_paths is true, then unitigs in all
    // bridge paths are kept. If keep_all_paths is false, then only the best path is kept.
    graph.clear_positions();
    for bridge in bridges {
        if bridge.conflicting {
            continue;
        }
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

    // TODO: add logic for removing non-anchor tips to handle blunt ends?
}


fn merge_after_bridging(graph: &mut UnitigGraph) {
    merge_linear_paths(graph, &vec![]);
    graph.print_basic_graph_info();
    graph.renumber_unitigs();
}


fn reduce_depths(graph: &mut UnitigGraph, bridge: &Bridge) {
    // This function is run after a bridge has been applied. It reduces the depth of unitigs in the
    // bridge.
    for path in &bridge.all_paths {
        for signed_num in path {
            let mut unitig = graph.unitig_index.get(&signed_num.unsigned_abs()).unwrap().borrow_mut();
            unitig.reduce_depth_by_one();
        }
    }
}


fn delete_unitigs_not_connected_to_anchor(graph: &mut UnitigGraph) {
    let to_delete: HashSet<u32> = graph.connected_components().into_iter()
        .filter_map(|component| {
            if component.iter().all(|&num| graph.unitig_index.get(&num).unwrap().borrow()
                               .unitig_type != UnitigType::Anchor) { Some(component) }
            else { None } })
        .flat_map(|component| component.into_iter())
        .collect();
    graph.remove_unitigs_by_number(to_delete);
}


fn cull_ambiguity(bridges: &mut Vec<Bridge>, verbose: bool) -> usize {
    let mut ambi_bridges: Vec<_> = bridges.iter().filter(|b| b.conflicting).collect();
    if ambi_bridges.is_empty() {
        return 0;
    }
    section_header("Culling conflicting bridges");
    explanation("The least-supported conflicting bridges are now culled until no bridges \
                 conflict.");
    ambi_bridges.sort_by(|a, b| a.depth().cmp(&b.depth()).then(a.cmp(b)));
    let mut cull_count = 0;
    if verbose {
        eprintln!("Culled bridges:");
    }
    while !ambi_bridges.is_empty() {
        let to_cull = ambi_bridges[0];
        if verbose {
            eprintln!("  {}", to_cull);
        }
        bridges.remove(bridges.iter().position(|b| b.start == to_cull.start && b.end == to_cull.end).unwrap());
        cull_count += 1;
        determine_ambiguity(bridges);
        ambi_bridges = bridges.iter().filter(|b| b.conflicting).collect();
        ambi_bridges.sort_by(|a, b| a.depth().cmp(&b.depth()).then(a.cmp(b)));
    }
    if verbose { eprintln!(); }
    eprintln!("{} conflicting bridge{} culled", cull_count, match cull_count { 1 => "", _ => "s" });
    eprintln!();
    cull_count
}


fn print_bridges(bridges: &Vec<Bridge>, verbose: bool) {
    let unique_count = bridges.iter().filter(|b| !b.conflicting).count();
    let conflicting_count = bridges.iter().filter(|b| b.conflicting).count();
    if verbose {
        if unique_count > 0 {
            eprintln!("Unique bridges:");
            for b in bridges {
                if !b.conflicting {
                    eprintln!("  {}", b);
                }
            }
        }
        if conflicting_count > 0 {
            eprintln!("\nConflicting bridges:");
            for b in bridges {
                if b.conflicting {
                    eprintln!("  {}", b);
                }
            }
        }
    } else {
        eprintln!("     Unique bridges: {}", unique_count);
        eprintln!("Conflicting bridges: {}", conflicting_count);
    }
    eprintln!();
}


fn get_anchor_to_anchor_paths(sequence_paths: &Vec<Vec<i32>>, anchor_set: &HashSet<u32>) -> Vec<Vec<i32>> {
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


fn group_paths_by_start_end(anchor_to_anchor_paths: Vec<Vec<i32>>) -> HashMap<(i32, i32), Vec<Vec<i32>>> {
    let mut grouped_paths: HashMap<(i32, i32), Vec<Vec<i32>>> = HashMap::new();
    for path in anchor_to_anchor_paths {
        if let (Some(&start), Some(&end)) = (path.first(), path.last()) {
            grouped_paths.entry((start, end)).or_default().push(path);
        }
    }
    grouped_paths
}


fn compare_paths(path_a: &Vec<i32>, path_b: &Vec<i32>) -> bool {
    // Returns true if path_a is 'better' than path_b. Used to break ties when two paths are equally
    // common in a bridge.

    // TODO: explore different tie-breaking strategies (just using lexographic order at the moment).
    path_a < path_b
}


pub struct Bridge {
    start: i32,
    end: i32,
    all_paths: Vec<Vec<i32>>,
    best_path: Vec<i32>,
    conflicting: bool,
}

impl Bridge {
    fn new(start: i32, end: i32, all_paths: Vec<Vec<i32>>) -> Self {

        // Remove the start and end unitigs from the paths.
        let mut trimmed_paths = all_paths;
        for path in &mut trimmed_paths {
            path.remove(0);
            path.pop();
        }

        // Set the best path to the most common path, using compare_paths to break ties.
        let mut path_counts = HashMap::new();
        for path in &trimmed_paths { *path_counts.entry(path).or_insert(0) += 1; }
        let mut best_path = Vec::new();
        let mut max_count = 0;
        for (path, count) in path_counts {
            if count > max_count || (count == max_count && compare_paths(path, &best_path)) {
                best_path = path.clone();
                max_count = count;
            }
        }

        Bridge {
            start,
            end,
            all_paths: trimmed_paths,
            best_path,
            conflicting: false,
        }
    }

    fn rev_start(&self) -> i32 {
        -self.end
    }

    fn rev_end(&self) -> i32 {
        -self.start
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
        let paths = vec![vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, -23, -8, 41, 2],
                         vec![1, 12, 17, 123, 41, 2]];
        let bridge = Bridge::new(1, 2, paths);
        assert_eq!(bridge.rev_start(), -2);
        assert_eq!(bridge.rev_end(), -1);
        assert_eq!(bridge.depth(), 4);
    }

    #[test]
    fn test_determine_ambiguity_1() {
        let bridge_a = Bridge::new(1, -2, vec![vec![1, 12, 2]]);
        let bridge_b = Bridge::new(-2, 5, vec![vec![-2, 6, 5]]);
        let bridge_c = Bridge::new(4, -5, vec![vec![4, -5]]);
        let bridge_d = Bridge::new(-4, 6, vec![vec![-4, 12, 6]]);
        let bridge_e = Bridge::new(-1, -6, vec![vec![-1, 11, -6]]);
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
        let bridge_a = Bridge::new(1, -2, vec![vec![1, 12, 2]]);
        let bridge_b = Bridge::new(-2, 5, vec![vec![-2, 6, 5]]);
        let bridge_c = Bridge::new(4, -5, vec![vec![4, -5]]);
        let bridge_d = Bridge::new(-4, 6, vec![vec![-4, 12, 6]]);
        let bridge_e = Bridge::new(-1, -6, vec![vec![-1, 11, -6]]);
        let bridge_f = Bridge::new(-4, 7, vec![vec![-4, 13, 7]]);
        let bridge_g = Bridge::new(1, 8, vec![vec![1, 14, 8]]);
        let bridge_h = Bridge::new(4, -8, vec![vec![4, 9, -8]]);
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
}
