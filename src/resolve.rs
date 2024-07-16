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

use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists};
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

    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO

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
    }
    eprintln!("{} anchor unitig{} found", anchor_ids.len(), match anchor_ids.len() { 1 => "", _ => "s" });
    eprintln!();
    anchor_ids
}
