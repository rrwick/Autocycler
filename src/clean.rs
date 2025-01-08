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


pub fn clean(in_gfa: PathBuf, out_gfa: PathBuf, remove: Option<String>) {
    check_settings(&in_gfa);
    starting_message();
    let remove_tigs = parse_remove_tigs(remove);
    print_settings(&in_gfa, &out_gfa, &remove_tigs);
    let mut graph = load_graph(&in_gfa);
    clean_graph(&mut graph, &remove_tigs, &in_gfa);
    merge_graph(&mut graph);
    graph.save_gfa(&out_gfa, &vec![]).unwrap();
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


fn print_settings(in_gfa: &Path, out_gfa: &Path, remove_tigs: &[u32]) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_gfa {}", out_gfa.display());
    if !remove_tigs.is_empty() {
        eprintln!("  --remove {}", remove_tigs.iter().map(|c| c.to_string())
                                              .collect::<Vec<String>>() .join(","));
    }
    eprintln!();
}


fn finished_message(out_gfa: &Path) {
    section_header("Finished!");
    eprintln!("Cleaned graph: {}", out_gfa.display());
    eprintln!();
}


fn clean_graph(graph: &mut UnitigGraph, remove_tigs: &[u32], in_gfa: &Path) {
    section_header("Cleaning graph");
    explanation("The user-specified tigs are now removed from the graph.");
    check_tig_nums_are_valid(in_gfa, graph, remove_tigs);
    let remove_set: HashSet<u32> = HashSet::from_iter(remove_tigs.iter().cloned());
    graph.remove_unitigs_by_number(remove_set);
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


fn check_tig_nums_are_valid(in_gfa: &Path, graph: &UnitigGraph, remove_tigs: &[u32]) {
    let mut all_tig_numbers = HashSet::new();
    for unitig in &graph.unitigs {
        all_tig_numbers.insert(unitig.borrow().number);
    }
    for tig in remove_tigs {
        if !all_tig_numbers.contains(tig) {
            quit_with_error(&format!("{} does not contain tig {}", in_gfa.display(), tig));
        }
    }
}


fn parse_remove_tigs(remove: Option<String>) -> Vec<u32> {
    if remove.is_none() {
        return Vec::new();
    }
    let remove = remove.unwrap().replace(' ', "");
    let mut remove_tigs: Vec<_> = remove.split(',')
            .map(|s| s.parse::<u32>().unwrap_or_else(|_| quit_with_error(
                &format!("failed to parse '{}' as a node number", s)))).collect();
    remove_tigs.sort();
    remove_tigs
}
