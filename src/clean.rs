// This file contains the code for the autocycler combine subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::misc::{check_if_file_exists, quit_with_error};
use crate::unitig_graph::UnitigGraph;


pub fn clean(in_gfa: PathBuf, out_gfa: PathBuf, remove: Option<String>, duplicate: Option<String>) {
    check_settings(&in_gfa);
    starting_message();
    let remove = parse_tig_numbers(remove);
    let duplicate = parse_tig_numbers(duplicate);
    print_settings(&in_gfa, &out_gfa, &remove, &duplicate);
    let mut graph = load_graph(&in_gfa);
    check_tig_numbers_are_valid(&in_gfa, &graph, &remove);
    check_tig_numbers_are_valid(&in_gfa, &graph, &duplicate);
    if !remove.is_empty() {
        remove_tigs(&mut graph, &remove);
    }
    if !duplicate.is_empty() {
        duplicate_tigs(&mut graph, &duplicate);
    }
    merge_graph(&mut graph);
    graph.save_gfa(&out_gfa, &vec![], true).unwrap();
    finished_message(&out_gfa);
}


fn check_settings(in_gfa: &Path) {
    check_if_file_exists(in_gfa);
}


fn starting_message() {
    section_header("Starting autocycler clean");
    explanation("This command removes user-specified tigs from a combined Autocycler graph and \
                 then merges all linear paths to produce a clean output graph.");
}


fn print_settings(in_gfa: &Path, out_gfa: &Path, remove: &[u32], duplicate: &[u32]) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_gfa {}", out_gfa.display());
    if !remove.is_empty() {
        eprintln!("  --remove {}", remove.iter().map(|c| c.to_string())
                                         .collect::<Vec<String>>() .join(","));
    }
    if !duplicate.is_empty() {
        eprintln!("  --remove {}", duplicate.iter().map(|c| c.to_string())
                                            .collect::<Vec<String>>() .join(","));
    }
    eprintln!();
}


fn finished_message(out_gfa: &Path) {
    section_header("Finished!");
    eprintln!("Cleaned graph: {}", out_gfa.display());
    eprintln!();
}


fn remove_tigs(graph: &mut UnitigGraph, remove: &[u32]) {
    section_header("Removing sequences");
    explanation("The user-specified tigs are now removed from the graph.");
    let remove_set: HashSet<u32> = HashSet::from_iter(remove.iter().cloned());
    graph.remove_unitigs_by_number(remove_set);
    graph.print_basic_graph_info();
}


fn duplicate_tigs(graph: &mut UnitigGraph, duplicate: &[u32]) {
    section_header("Duplicating sequences");
    explanation("The user-specified tigs are now duplicated in the graph.");
    for tig_num in duplicate {
        graph.duplicate_unitig_by_number(tig_num);
    }
    graph.print_basic_graph_info();
}


fn merge_graph(graph: &mut UnitigGraph) {
    section_header("Merging linear paths");
    explanation("Linear paths in the graph are now merged.");
    merge_linear_paths(graph, &vec![]);
    graph.print_basic_graph_info();
    graph.renumber_unitigs();
}


fn load_graph(gfa: &Path) -> UnitigGraph {
    section_header("Loading graph");
    explanation("The unitig graph is now loaded into memory.");
    let (graph, _) = UnitigGraph::from_gfa_file(gfa);
    graph.print_basic_graph_info();
    graph
}


fn check_tig_numbers_are_valid(in_gfa: &Path, graph: &UnitigGraph, tig_numbers: &[u32]) {
    let mut all_tig_numbers = HashSet::new();
    for unitig in &graph.unitigs {
        all_tig_numbers.insert(unitig.borrow().number);
    }
    for tig in tig_numbers {
        if !all_tig_numbers.contains(tig) {
            quit_with_error(&format!("{} does not contain tig {}", in_gfa.display(), tig));
        }
    }
}


fn parse_tig_numbers(tig_num_str: Option<String>) -> Vec<u32> {
    if tig_num_str.is_none() {
        return Vec::new();
    }
    let tig_num_str = tig_num_str.unwrap().replace(' ', "");
    let mut tig_numbers: Vec<_> = tig_num_str.split(',')
            .map(|s| s.parse::<u32>().unwrap_or_else(|_| quit_with_error(
                &format!("failed to parse '{}' as a node number", s)))).collect();
    tig_numbers.sort();
    tig_numbers
}
