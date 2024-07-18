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

use std::cell::RefCell;
use std::collections::HashSet;
use std::path::PathBuf;
use std::rc::Rc;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists, strand};
use crate::position::Position;
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;
use crate::unitig::{Unitig, UnitigStrand};


pub fn resolve(cluster_dir: PathBuf) {
    let trimmed_gfa = cluster_dir.join("3_trimmed.gfa");
    let resolved_gfa = cluster_dir.join("5_resolved.gfa");
    check_settings(&cluster_dir, &trimmed_gfa);
    starting_message();
    print_settings(&cluster_dir);
    let (mut unitig_graph, sequences) = load_graph(&trimmed_gfa);
    let anchors = find_anchor_unitigs(&mut unitig_graph, &sequences);
    let bridges = create_bridges(&unitig_graph, &anchors);
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

        // TODO: dead-ends (with or without a hairpin connection) should also become anchors.
    }
    eprintln!("{} anchor unitig{} found", anchor_ids.len(), match anchor_ids.len() { 1 => "", _ => "s" });
    eprintln!();
    anchor_ids
}


fn create_bridges(graph: &UnitigGraph, anchors: &Vec<u32>) -> Vec<Bridge> {
    section_header("Building bridges");
    explanation("Bridges connect one anchor unitig to the next.");
    let mut bridges = Vec::new();
    let anchor_set: HashSet<u32> = anchors.iter().cloned().collect();
    for a in anchors {
        let unitig_rc = graph.unitig_index.get(a).unwrap();
        let unitig = unitig_rc.borrow();
        for strand in [strand::FORWARD, strand::REVERSE] {
            let unitig_strand = UnitigStrand::new(unitig_rc, strand);
            eprintln!("\n{}", unitig_strand);  // TEMP
            let positions = match strand {
                strand::FORWARD => &unitig.forward_positions,
                _ => &unitig.reverse_positions
            };
            for p in positions {
                let path = get_path_to_next_anchor(graph, unitig_rc, strand, p, &anchor_set);
                eprintln!("    {:?}", path);  // TEMP
            }


        }
    }
    bridges
}


fn get_path_to_next_anchor(graph: &UnitigGraph, unitig_rc: &Rc<RefCell<Unitig>>, strand: bool,
                           pos: &Position, anchor_set: &HashSet<u32>) -> Vec<i32> {
    let seq_id = pos.seq_id();
    let seq_strand = pos.strand();
    let mut u = UnitigStrand::new(unitig_rc, strand);
    let mut p = pos.pos;
    let mut path = vec![u.signed_num()];
    while let Some((next_u, next_p)) = graph.get_next_unitig(seq_id, seq_strand, &u.unitig, u.strand, p) {
        u = next_u; p = next_p;
        path.push(u.signed_num());
        if anchor_set.contains(&u.number()) {  // the path reached an anchor
            return path;
        }
    }
    vec![]  // the path ran out before reaching an anchor
}


pub struct Bridge {
    start: i32,
    end: i32,
    all_paths: Vec<Vec<i32>>,
    best_path: Vec<i32>,
}
