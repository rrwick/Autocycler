// This file contains functions related to manipulating a UnitigGraph in order to simplify its
// structure.

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
use std::rc::Rc;

use crate::sequence::Sequence;
use crate::unitig::Unitig;
use crate::unitig_graph::UnitigGraph;


pub fn simplify_structure(graph: &mut UnitigGraph, seqs: &Vec<Sequence>) {
    // This function simplifies the graph structure by expanding repeats.
    //
    // For example, it will turn this:
    //    ACTACTCAACT
    //               \
    //                ATCGACTACGCTACG
    //               /
    //    GACTACGAACT
    //
    // Into this:
    //    ACTACTC
    //           \
    //            AACTATCGACTACGCTACG
    //           /
    //    GACTACG
    //
    let (fixed_starts, fixed_ends) = get_fixed_unitig_starts_and_ends(graph, seqs);
    for unitig_rc in &graph.unitigs {
        let unitig_number = unitig_rc.borrow().number;
        let inputs = get_exclusive_inputs(graph, &unitig_rc);
        if inputs.len() >= 2 && !fixed_starts.contains(&unitig_number) {
            // TODO: check inputs for fixed_starts/fixed_ends.
            shift_sequence_1(&inputs, &unitig_rc);
        }
        let outputs = get_exclusive_outputs(graph, &unitig_rc);
        if outputs.len() >= 2 && !fixed_ends.contains(&unitig_number) {
            // TODO: check puts for fixed_starts/fixed_ends.
            shift_sequence_2(&unitig_rc, &outputs);
        }
    }
    graph.renumber_unitigs();
}

fn shift_sequence_1(sources: &Vec<(Rc<RefCell<Unitig>>, bool)>, destination_rc: &Rc<RefCell<Unitig>>) {
    // This function:
    // * removes any common sequence from the ends of the source unitigs
    // * adds that common sequence to the start of the destination unitig
    eprintln!("\n"); // TEMP
    for (source_rc, strand) in sources {
        let mut source = source_rc.borrow_mut();
        eprintln!("SOURCE: {}", source); // TEMP
    }
    let mut destination = destination_rc.borrow_mut();
    eprintln!("DESTINATION: {}", destination); // TEMP
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    eprintln!("\n"); // TEMP
}

fn shift_sequence_2(destination_rc: &Rc<RefCell<Unitig>>, sources: &Vec<(Rc<RefCell<Unitig>>, bool)>) {
    // This function:
    // * removes any common sequence from the starts of the source unitigs
    // * adds that common sequence to the end of the destination unitig
    eprintln!("\n"); // TEMP
    for (source_rc, strand) in sources {
        let mut source = source_rc.borrow_mut();
        eprintln!("SOURCE: {}", source); // TEMP
    }
    let mut destination = destination_rc.borrow_mut();
    eprintln!("DESTINATION: {}", destination); // TEMP
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    eprintln!("\n"); // TEMP
}

fn get_fixed_unitig_starts_and_ends(graph: &UnitigGraph,
                                    sequences: &Vec<Sequence>) -> (HashSet<u32>, HashSet<u32>) {
    // Returns two sets of unitig IDs: all unitigs where the start can't be changed and all
    // unitigs where the end can't be changed. All results are in terms of the unitig's forward
    // strand.
    let mut fixed_starts = HashSet::new();
    let mut fixed_ends = HashSet::new();
    for seq in sequences {
        let unitig_path = graph.get_unitig_path_for_sequence(seq);
        if unitig_path.len() == 0 {
            continue
        }
        let (first_unitig, first_strand) = unitig_path[0];
        if first_strand {
            fixed_starts.insert(first_unitig);
        } else {
            fixed_ends.insert(first_unitig);
        }
        let (last_unitig, last_strand) = unitig_path.last().unwrap();
        if *last_strand {
            fixed_ends.insert(*last_unitig);
        } else {
            fixed_starts.insert(*last_unitig);
        }
    }
    (fixed_starts, fixed_ends)
}

fn get_exclusive_inputs(graph: &UnitigGraph,
                        unitig_rc: &Rc<RefCell<Unitig>>) -> Vec<(Rc<RefCell<Unitig>>, bool)> {
    // This function returns a vector of unitigs which exclusively input to the given unitig.
    // Exclusive input means the unitig leads only to the given unitig. If any of the given
    // unitig's inputs are not exclusive inputs, then this function returns an empty vector.
    let mut inputs = Vec::new();
    let unitig = unitig_rc.borrow();
    for (prev_unitig_rc, prev_strand) in &unitig.forward_prev {
        let prev_unitig = &prev_unitig_rc.borrow();
        let prev_next_unitigs = if *prev_strand { &prev_unitig.forward_next } else { &prev_unitig.reverse_next };
        if prev_next_unitigs.len() != 1 {
            return Vec::new();
        }
        let (prev_next_unitig_rc, prev_next_strand) = &prev_next_unitigs[0];
        if *prev_next_strand && prev_next_unitig_rc.borrow().number == unitig.number {
            inputs.push((Rc::clone(&prev_unitig_rc), *prev_strand));
        } else {
            return Vec::new();
        }
    }
    inputs
}

fn get_exclusive_outputs(graph: &UnitigGraph,
                         unitig_rc: &Rc<RefCell<Unitig>>) -> Vec<(Rc<RefCell<Unitig>>, bool)> {
    // This function returns a vector of unitigs which exclusively output from the given unitig.
    // Exclusive output means the given unitig leads only to the unitig. If any of the given
    // unitig's outputs are not exclusive outputs, then this function returns an empty vector.
    let mut outputs = Vec::new();
    let unitig = unitig_rc.borrow();
    for (next_unitig_rc, next_strand) in &unitig.forward_next {
        let next_unitig = &next_unitig_rc.borrow();
        let next_prev_unitigs = if *next_strand { &next_unitig.forward_prev } else { &next_unitig.reverse_prev };
        if next_prev_unitigs.len() != 1 {
            return Vec::new();
        }
        let (next_prev_unitig_rc, next_prev_strand) = &next_prev_unitigs[0];
        if *next_prev_strand && next_prev_unitig_rc.borrow().number == unitig.number {
            outputs.push((Rc::clone(&next_unitig_rc), *next_strand));
        } else {
            return Vec::new();
        }
    }
    outputs
}
