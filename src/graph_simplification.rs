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

use crate::misc::{reverse_complement, strand};
use crate::position::Position;
use crate::sequence::Sequence;
use crate::unitig::{Unitig, UnitigStrand, UnitigType};
use crate::unitig_graph::UnitigGraph;


pub fn simplify_structure(graph: &mut UnitigGraph, seqs: &Vec<Sequence>) {
    while expand_repeats(graph, seqs) > 0 {}

    // TODO: sometimes the simplified graph ends up with a little redundant dead-end contig. This
    //       occurs because graph simplification won't allow contigs to be shortened to 0-bp. So
    //       some additional graph-simplification logic to clean these up could be nice.
    //
    //       GACTACG - T
    //                  \
    //                   ATCGACTACGCTACG
    //                  /
    //                 T

    graph.renumber_unitigs();
}


fn expand_repeats(graph: &mut UnitigGraph, seqs: &Vec<Sequence>) -> usize {
    // This function simplifies the graph structure by expanding repeats.
    //
    // For example, it will turn this:
    //    ACTACTCAACT                 GCTACGACTAC
    //               \               /
    //                ATCGACTACGCTACG
    //               /               \
    //    GACTACGAACT                 GCTATTGTACC
    //
    // Into this:
    //    ACTACTC                         CGACTAC
    //           \                       /
    //            AACTATCGACTACGCTACGGCTA
    //           /                       \
    //    GACTACG                         TTGTACC
    //
    // To avoid messing with input sequence paths, this function will not shift sequences at the
    // start/ends of such paths. The return value is the total amount of sequence shifted.
    let (fixed_starts, fixed_ends) = get_fixed_unitig_starts_and_ends(graph, seqs);
    let mut total_shifted_seq = 0;
    for unitig_rc in &graph.unitigs {
        let unitig_number = unitig_rc.borrow().number;
        let inputs = get_exclusive_inputs(unitig_rc);
        if inputs.len() >= 2 && !fixed_starts.contains(&unitig_number) {
            let can_shift = inputs.iter().all(|input| {
                !(input.strand && fixed_ends.contains(&input.number()) ||
                !input.strand && fixed_starts.contains(&input.number()))});
            if can_shift {
                total_shifted_seq += shift_sequence_1(&inputs, unitig_rc);
            }
        }
        let outputs = get_exclusive_outputs(unitig_rc);
        if outputs.len() >= 2 && !fixed_ends.contains(&unitig_number) {
            let can_shift = outputs.iter().all(|output| {
                !(output.strand && fixed_starts.contains(&output.number()) ||
                !output.strand && fixed_ends.contains(&output.number()))});
            if can_shift {
                total_shifted_seq += shift_sequence_2(unitig_rc, &outputs);
            }
        }
    }
    total_shifted_seq
}


fn shift_sequence_1(sources: &Vec<UnitigStrand>, destination_rc: &Rc<RefCell<Unitig>>) -> usize {
    // This function:
    // * removes any common sequence from the ends of the source unitigs
    // * adds that common sequence to the start of the destination unitig
    //
    // This function also guards against a couple of potential complications with sequence paths
    // (which could result in a path having more than one starting unitig):
    // * It won't let unitigs get down to a length of zero. This requires some extra logic for when
    //   both strands of one unitig appear in the sources (and will therefore have sequence
    //   removed from both ends).
    // * It won't add sequence to the destination unitig causing any of its positions to reach the
    //   start of a path.
    //
    // The return value is the amount of sequence shifted.
    let mut common_seq = get_common_end_seq(sources);
    avoid_zero_len_unitigs(&mut common_seq, sources, true);
    avoid_start_of_path(&mut common_seq, destination_rc, true);
    let shifted_amount = common_seq.len();
    if shifted_amount == 0 {
        return 0;
    }
    for source in sources {
        if source.strand {
            source.unitig.borrow_mut().remove_seq_from_end(shifted_amount);
        } else {
            source.unitig.borrow_mut().remove_seq_from_start(shifted_amount);
        }
    }
    destination_rc.borrow_mut().add_seq_to_start(common_seq);
    shifted_amount
}


fn shift_sequence_2(destination_rc: &Rc<RefCell<Unitig>>, sources: &[UnitigStrand]) -> usize {
    // This function does the same thing as shift_sequence_1, but for the other side of a unitig:
    // * removes any common sequence from the starts of the source unitigs
    // * adds that common sequence to the end of the destination unitig
    let mut common_seq = get_common_start_seq(sources);
    avoid_zero_len_unitigs(&mut common_seq, sources, false);
    avoid_start_of_path(&mut common_seq, destination_rc, false);
    let shifted_amount = common_seq.len();
    if shifted_amount == 0 {
        return 0;
    }
    for source in sources {
        if source.strand {
            source.unitig.borrow_mut().remove_seq_from_start(shifted_amount);
        } else {
            source.unitig.borrow_mut().remove_seq_from_end(shifted_amount);
        }
    }
    destination_rc.borrow_mut().add_seq_to_end(common_seq);
    shifted_amount
}


fn avoid_zero_len_unitigs(common_seq: &mut Vec<u8>, sources: &[UnitigStrand], trim_from_start: bool) {
    // This function takes some common sequence (sequence that will be shifted from some unitigs
    // onto another) and trims it down to ensure that none of the source unitigs will end up with
    // a length of zero.
    if common_seq.is_empty() {
        return;
    }
    let dup = if check_for_duplicates(sources) { 2 } else { 1 };
    let min_source_len = sources.iter().map(|source| source.length()).min().unwrap();
    while min_source_len <= (common_seq.len() as u32) * dup {
        if trim_from_start {
            common_seq.remove(0);
        } else {
            common_seq.pop();
        }
    }
}


fn avoid_start_of_path(common_seq: &mut Vec<u8>, dest_rc: &Rc<RefCell<Unitig>>,
                       trim_from_start: bool) {
    // This function takes some common sequence (sequence that will be shifted from some unitigs
    // onto another) and trims it down to ensure that the destination unitig's positions will not
    // end up at the start of a path, as this can cause problems.
    if common_seq.is_empty() {
        return;
    }
    if trim_from_start {
        while dest_rc.borrow().forward_positions.iter().any(|p| p.pos <= common_seq.len() as u32) {
            common_seq.remove(0);
        }
    } else {
        while dest_rc.borrow().reverse_positions.iter().any(|p| p.pos <= common_seq.len() as u32) {
            common_seq.pop();
        }
    }
}


fn check_for_duplicates(unitigs: &[UnitigStrand]) -> bool {
    // Returns true if any two unitigs in the vector have the same number.
    unitigs.iter().map(|u| u.number()).collect::<HashSet<_>>().len() != unitigs.len()
}


fn get_fixed_unitig_starts_and_ends(graph: &UnitigGraph,
                                    sequences: &Vec<Sequence>) -> (HashSet<u32>, HashSet<u32>) {
    // Returns two sets of unitig IDs: all unitigs where the start can't be changed and all
    // unitigs where the end can't be changed. All results are in terms of the unitig's forward
    // strand.
    let mut fixed_starts = HashSet::new();
    let mut fixed_ends = HashSet::new();

    // The starts/end of sequence paths are fixed.
    for seq in sequences {
        let unitig_path = graph.get_unitig_path_for_sequence(seq);
        if unitig_path.is_empty() { continue; }
        let (first_unitig, first_strand) = unitig_path[0];
        if first_strand { fixed_starts.insert(first_unitig); }
                   else { fixed_ends.insert(first_unitig); }
        let (last_unitig, last_strand) = unitig_path.last().unwrap();
        if *last_strand { fixed_ends.insert(*last_unitig); }
                   else { fixed_starts.insert(*last_unitig); }
    }

    let fixed_starts_copy = fixed_starts.clone();
    let fixed_ends_copy = fixed_ends.clone();

    // Any unitig which is upstream of a fixed start has a fixed end.
    for u in &fixed_starts_copy {
        for upstream in &graph.unitig_index.get(u).unwrap().borrow().forward_prev {
            if upstream.strand { fixed_ends.insert(upstream.number()); }
                          else { fixed_starts.insert(upstream.number()); }
        }
    }

    // Any unitig which is downstream of a fixed end has a fixed start.
    for u in &fixed_ends_copy {
        for downstream in &graph.unitig_index.get(u).unwrap().borrow().forward_next {
            if downstream.strand { fixed_starts.insert(downstream.number()); }
                            else { fixed_ends.insert(downstream.number()); }
        }
    }

    (fixed_starts, fixed_ends)
}


fn get_exclusive_inputs(unitig_rc: &Rc<RefCell<Unitig>>) -> Vec<UnitigStrand> {
    // This function returns a vector of unitigs which exclusively input to the given unitig.
    // Exclusive input means the unitig leads only to the given unitig. If any of the given
    // unitig's inputs are not exclusive inputs, then this function returns an empty vector.
    let mut inputs = Vec::new();
    let unitig = unitig_rc.borrow();
    for prev in &unitig.forward_prev {
        let prev_unitig = &prev.unitig.borrow();
        let prev_next_unitigs = if prev.strand { &prev_unitig.forward_next } else { &prev_unitig.reverse_next };
        if prev_next_unitigs.len() != 1 {
            return Vec::new();
        }
        let prev_next = &prev_next_unitigs[0];
        if prev_next.strand && prev_next.number() == unitig.number {
            inputs.push(UnitigStrand::new(&prev.unitig, prev.strand));
        } else {
            return Vec::new();
        }
    }
    if inputs.iter().any(|input| { input.number() == unitig.number }) {
        return Vec::new();
    }
    inputs
}


fn get_exclusive_outputs(unitig_rc: &Rc<RefCell<Unitig>>) -> Vec<UnitigStrand> {
    // This function returns a vector of unitigs which exclusively output from the given unitig.
    // Exclusive output means the given unitig leads only to the unitig. If any of the given
    // unitig's outputs are not exclusive outputs, then this function returns an empty vector.
    let mut outputs = Vec::new();
    let unitig = unitig_rc.borrow();
    for next in &unitig.forward_next {
        let next_unitig = &next.unitig.borrow();
        let next_prev_unitigs = if next.strand { &next_unitig.forward_prev } else { &next_unitig.reverse_prev };
        if next_prev_unitigs.len() != 1 {
            return Vec::new();
        }
        let next_prev = &next_prev_unitigs[0];
        if next_prev.strand && next_prev.number() == unitig.number {
            outputs.push(UnitigStrand::new(&next.unitig, next.strand));
        } else {
            return Vec::new();
        }
    }
    if outputs.iter().any(|output| { output.number() == unitig.number }) {
        return Vec::new();
    }
    outputs
}


fn get_common_start_seq(unitigs: &[UnitigStrand]) -> Vec<u8> {
    // This function returns the common sequence at the start of all given unitigs.
    let seqs: Vec<_> = unitigs.iter().map(|u| u.get_seq()).collect();
    if seqs.is_empty() { return Vec::new(); }
    let mut prefix = seqs[0].clone();
    for seq in seqs.iter() {
        while !seq.starts_with(&prefix) {
            prefix.pop();
            if prefix.is_empty() { return Vec::new(); }
        }
    }
    prefix
}


fn get_common_end_seq(unitigs: &[UnitigStrand]) -> Vec<u8> {
    // This function returns the common sequence at the end of all given unitigs.
    let seqs: Vec<Vec<u8>> = unitigs.iter().map(|u| u.get_seq())
        .map(|mut seq| { seq.reverse(); seq }).collect();
    if seqs.is_empty() { return Vec::new(); }
    let mut suffix = seqs[0].clone();
    for seq in seqs.iter() {
        while !seq.starts_with(&suffix) {
            suffix.pop();
            if suffix.is_empty() { return Vec::new(); }
        }
    }
    suffix.reverse();
    suffix
}


pub fn merge_linear_paths(graph: &mut UnitigGraph, seqs: &Vec<Sequence>) {
    // This function looks for linear paths in the graph (where one Unitig leads only to another
    // and vice versa) and merges them together when possible.
    //
    // For example, it will turn this:
    //    ACTACTCAACT - ATCGACTACGCTACG
    //
    // Into this:
    //    ACTACTCAACTATCGACTACGCTACG
    //
    // To avoid messing with input sequence paths, this function will not merge sequences at the
    // start/ends of such paths. If no sequences are provided, then all possible linear paths will
    // be merged.
    let (mut fixed_starts, fixed_ends) = get_fixed_unitig_starts_and_ends(graph, seqs);
    fix_circular_loops(graph, &mut fixed_starts);
    let mut already_used = HashSet::new();
    let mut merge_paths = Vec::new();
    for unitig_rc in &graph.unitigs {
        let unitig_number = unitig_rc.borrow().number;
        for unitig_strand in [strand::FORWARD, strand::REVERSE] {

            // Find unitigs which can potentially start a mergeable path.
            if already_used.contains(&unitig_number) { continue; }
            if has_single_exclusive_input(unitig_rc, unitig_strand) && can_merge_start(unitig_number, unitig_strand, &fixed_starts, &fixed_ends) { continue; }
            let mut current_path = vec![UnitigStrand::new(unitig_rc, unitig_strand)];
            already_used.insert(unitig_number);

            // Extend the path as far as possible.
            loop {
                let unitig = current_path.last().unwrap();
                if cannot_merge_end(unitig.number(), unitig.strand, &fixed_starts, &fixed_ends) { break; }
                let mut outputs = if unitig.strand { get_exclusive_outputs(&unitig.unitig) } else { get_exclusive_inputs(&unitig.unitig) };
                if outputs.len() != 1 { break; }
                let output = &mut outputs[0];
                if !unitig.strand { output.strand = !output.strand; }
                let output_number = output.number();
                if already_used.contains(&output_number) { break; }
                if cannot_merge_start(output_number, output.strand, &fixed_starts, &fixed_ends) { break; }
                current_path.push(UnitigStrand::new(&output.unitig, output.strand));
                already_used.insert(output_number);
            }

            if current_path.len() > 1 {
                merge_paths.push(current_path.clone());
            }
        }
    }

    let mut new_unitig_number: u32 = graph.max_unitig_number();
    for path in merge_paths {
        new_unitig_number += 1;
        merge_path(graph, &path, new_unitig_number);
    }
    graph.delete_dangling_links();
    graph.build_unitig_index();
    graph.check_links();
}


fn fix_circular_loops(graph: &UnitigGraph, fixed_starts: &mut HashSet<u32>) {
    // This function looks for components of the graph which form a simple circular loop, and if
    // any are found, adds their lowest numbered unitig to fixed_starts. This will allow for the
    // component to be merged into a single unitig.
    let components = graph.connected_components();
    for component in components {
        if graph.component_is_circular_loop(&component) {
            fixed_starts.insert(component[0]);
        }
    }
}


fn cannot_merge_start(unitig_number: u32, unitig_strand: bool, fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> bool {
    // Checks whether a given unitig (by number and strand) can have its start merged with other unitigs.
    (unitig_strand && fixed_starts.contains(&unitig_number)) || (!unitig_strand && fixed_ends.contains(&unitig_number))
}


fn can_merge_start(unitig_number: u32, unitig_strand: bool, fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> bool {
    !cannot_merge_start(unitig_number, unitig_strand, fixed_starts, fixed_ends)
}


fn cannot_merge_end(unitig_number: u32, unitig_strand: bool, fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> bool {
    // Checks whether a given unitig (by number and strand) can have its end merged with other unitigs.
    (unitig_strand && fixed_ends.contains(&unitig_number)) || (!unitig_strand && fixed_starts.contains(&unitig_number))
}


fn has_single_exclusive_input(unitig_rc: &Rc<RefCell<Unitig>>, unitig_strand: bool) -> bool {
    let inputs = if unitig_strand {get_exclusive_inputs(unitig_rc)} else {get_exclusive_outputs(unitig_rc)};
    inputs.len() == 1
}


fn merge_path(graph: &mut UnitigGraph, path: &Vec<UnitigStrand>, new_unitig_number: u32) {
    let merged_seq = merge_unitig_seqs(path);
    let first = &path[0];
    let last = path.last().unwrap();
    let forward_positions = if first.strand {first.unitig.borrow().forward_positions.clone()} else {first.unitig.borrow().reverse_positions.clone()};
    let reverse_positions = if last.strand {last.unitig.borrow().reverse_positions.clone()} else {last.unitig.borrow().forward_positions.clone()};

    // Check to see if the path has any self links, so we can make those after the merge if needed.
    let end_to_start_link = graph.link_exists(last.number(), last.strand, first.number(), first.strand);
    let start_flip_link = graph.link_exists(first.number(), !first.strand, first.number(), first.strand);
    let end_flip_link = graph.link_exists(last.number(), last.strand, last.number(), !last.strand);

    // For the new unitig, we take links (forward_prev, reverse_next, forward_next, reverse_prev)
    // from the first/last unitigs in the path.
    let forward_prev = if first.strand {first.unitig.borrow().forward_prev.clone()} else {first.unitig.borrow().reverse_prev.clone()};
    let reverse_next = if first.strand {first.unitig.borrow().reverse_next.clone()} else {first.unitig.borrow().forward_next.clone()};
    let forward_next = if last.strand {last.unitig.borrow().forward_next.clone()} else {last.unitig.borrow().reverse_next.clone()};
    let reverse_prev = if last.strand {last.unitig.borrow().reverse_prev.clone()} else {last.unitig.borrow().forward_prev.clone()};

    let mut unitig = Unitig {
        number: new_unitig_number,
        reverse_seq: reverse_complement(&merged_seq),
        forward_seq: merged_seq,
        depth: get_merge_path_depth(path, &forward_positions),
        forward_positions, reverse_positions,
        forward_next, forward_prev, reverse_next, reverse_prev,
        ..Default::default()
    };

    if path.iter().any(|p| p.is_anchor() || p.is_consentig()) {
        unitig.unitig_type = UnitigType::Consentig;
    }

    let unitig_rc = Rc::new(RefCell::new(unitig));
    graph.unitigs.push(unitig_rc.clone());

    // Create links to the new unitig from its neighbours.
    for u in &unitig_rc.borrow().forward_next {
        if u.strand {u.unitig.borrow_mut().forward_prev.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));}
              else {u.unitig.borrow_mut().reverse_prev.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));}
    }
    for u in &unitig_rc.borrow().forward_prev {
        if u.strand {u.unitig.borrow_mut().forward_next.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));}
              else {u.unitig.borrow_mut().reverse_next.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));}
    }
    for u in &unitig_rc.borrow().reverse_next {
        if u.strand {u.unitig.borrow_mut().forward_prev.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));}
              else {u.unitig.borrow_mut().reverse_prev.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));}
    }
    for u in &unitig_rc.borrow().reverse_prev {
        if u.strand {u.unitig.borrow_mut().forward_next.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));}
              else {u.unitig.borrow_mut().reverse_next.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));}
    }

    // Create any needed links from the new unitig to itself.
    if end_to_start_link {
        let mut u = unitig_rc.borrow_mut();
        u.forward_next.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));
        u.forward_prev.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));
        u.reverse_next.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));
        u.reverse_prev.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));
    }
    if start_flip_link {
        let mut u = unitig_rc.borrow_mut();
        u.reverse_next.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));
        u.forward_prev.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));
    }
    if end_flip_link {
        let mut u = unitig_rc.borrow_mut();
        u.forward_next.push(UnitigStrand::new(&unitig_rc, strand::REVERSE));
        u.reverse_prev.push(UnitigStrand::new(&unitig_rc, strand::FORWARD));
    }

    let path_numbers: HashSet<_> = path.iter().map(|u| u.number()).collect();
    graph.unitigs.retain(|u| !path_numbers.contains(&u.borrow().number));
}


fn merge_unitig_seqs(path: &Vec<UnitigStrand>) -> Vec<u8> {
    // Given a path of unitigs (with their strand), this function returns their merged sequence. It
    // assumes no overlap and it does not check that the given unitigs are actually linked to each
    // other.
    let total_length: usize = path.iter().map(|u| u.length()).sum::<u32>().try_into().unwrap();
    let mut merged_seq = Vec::with_capacity(total_length);
    for u in path {
        merged_seq.extend(u.get_seq());
    }
    merged_seq
}


fn get_merge_path_depth(path: &Vec<UnitigStrand>, forward_positions: &[Position]) -> f64 {
    // If the unitigs have position information, use that to determine depth.
    if !forward_positions.is_empty() {
        return forward_positions.len() as f64;
    }

    // If the path contains an anchor unitig, set the merged depth to the anchor's depth.
    for u in path {
        if u.is_anchor() {
            return u.depth();
        }
    }

    // Otherwise, give a weighted mean depth of the unitigs in the path.
    weighted_mean_depth(path)
}


fn weighted_mean_depth(path: &Vec<UnitigStrand>) -> f64 {
    let total_length = path.iter().map(|u| u.length()).sum::<u32>() as f64;
    let mut depth_sum = 0.0;
    for u in path {
        depth_sum += u.depth() * u.length() as f64;
    }
    depth_sum / total_length
}


#[cfg(test)]
mod tests {
    use crate::test_gfa::*;
    use super::*;

    fn unitig_vec_to_str(mut unitigs: Vec<UnitigStrand>) -> String {
        // Converts a vector of Unitigs to a string (makes my tests easier to write).
        unitigs.sort_by(|a, b| {a.number().cmp(&b.number()).then_with(|| a.strand.cmp(&b.strand))});
        unitigs.iter().map(|u| {format!("{}{}", u.number(), if u.strand {'+'} else {'-'})}).collect::<Vec<String>>().join(",")
    }

    #[test]
    fn test_get_common_start_seq() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&get_common_start_seq(&unitigs)).unwrap(), "AC");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::REVERSE)];
        assert_eq!(std::str::from_utf8(&get_common_start_seq(&unitigs)).unwrap(), "A");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::REVERSE), UnitigStrand::new(&c, strand::REVERSE)];
        assert_eq!(std::str::from_utf8(&get_common_start_seq(&unitigs)).unwrap(), "");
    }

    #[test]
    fn test_get_common_end_seq() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&get_common_end_seq(&unitigs)).unwrap(), "");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::REVERSE), UnitigStrand::new(&b, strand::REVERSE), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&get_common_end_seq(&unitigs)).unwrap(), "T");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::REVERSE), UnitigStrand::new(&b, strand::REVERSE), UnitigStrand::new(&c, strand::REVERSE)];
        assert_eq!(std::str::from_utf8(&get_common_end_seq(&unitigs)).unwrap(), "GT");
    }

    #[test]
    fn test_get_exclusive_inputs_and_outputs() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());

        let unitig_1 = &graph.unitigs[0];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_1)), "2+,3-");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_1)), "");

        let unitig_2 = &graph.unitigs[1];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_2)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_2)), "");

        let unitig_3 = &graph.unitigs[2];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_3)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_3)), "");

        let unitig_4 = &graph.unitigs[3];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_4)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_4)), "7-,8+");

        let unitig_5 = &graph.unitigs[4];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_5)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_5)), "");

        let unitig_6 = &graph.unitigs[5];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_6)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_6)), "");

        let unitig_7 = &graph.unitigs[6];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_7)), "9-,9+");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_7)), "");

        let unitig_8 = &graph.unitigs[7];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_8)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_8)), "10-");

        let unitig_9 = &graph.unitigs[8];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_9)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_9)), "");

        let unitig_10 = &graph.unitigs[9];
        assert_eq!(unitig_vec_to_str(get_exclusive_inputs(unitig_10)), "");
        assert_eq!(unitig_vec_to_str(get_exclusive_outputs(unitig_10)), "8-");
    }

    #[test]
    fn test_simplify_structure_1() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        let sequences: Vec<Sequence> = vec![];

        assert_eq!(std::str::from_utf8(&graph.unitigs[0].borrow().forward_seq).unwrap(), "TTCGCTGCGCTCGCTTCGCTTT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[1].borrow().forward_seq).unwrap(), "TGCCGTCGTCGCTGTGCA");
        assert_eq!(std::str::from_utf8(&graph.unitigs[2].borrow().forward_seq).unwrap(), "TGCCTGAATCGCCTA");
        assert_eq!(std::str::from_utf8(&graph.unitigs[3].borrow().forward_seq).unwrap(), "GCTCGGCTCG");
        assert_eq!(std::str::from_utf8(&graph.unitigs[4].borrow().forward_seq).unwrap(), "CGAACCAT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[5].borrow().forward_seq).unwrap(), "TACTTGT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[6].borrow().forward_seq).unwrap(), "GCCTT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[7].borrow().forward_seq).unwrap(), "ATCT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[8].borrow().forward_seq).unwrap(), "GC");
        assert_eq!(std::str::from_utf8(&graph.unitigs[9].borrow().forward_seq).unwrap(), "T");

        simplify_structure(&mut graph, &sequences);

        assert_eq!(std::str::from_utf8(&graph.unitigs[0].borrow().forward_seq).unwrap(), "GCATTCGCTGCGCTCGCTTCGCTTT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[1].borrow().forward_seq).unwrap(), "TGCCGTCGTCGCTGT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[2].borrow().forward_seq).unwrap(), "CTGAATCGCCTA");
        assert_eq!(std::str::from_utf8(&graph.unitigs[3].borrow().forward_seq).unwrap(), "GCTCGGCTCGA");
        assert_eq!(std::str::from_utf8(&graph.unitigs[4].borrow().forward_seq).unwrap(), "CGAACCAT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[5].borrow().forward_seq).unwrap(), "TACTTGT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[6].borrow().forward_seq).unwrap(), "GCCT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[7].borrow().forward_seq).unwrap(), "TCT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[8].borrow().forward_seq).unwrap(), "GC");
        assert_eq!(std::str::from_utf8(&graph.unitigs[9].borrow().forward_seq).unwrap(), "T");
    }

    #[test]
    fn test_simplify_structure_2() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        let sequences: Vec<Sequence> = vec![];

        assert_eq!(std::str::from_utf8(&graph.unitigs[0].borrow().forward_seq).unwrap(), "ACCGCTGCGCTCGCTTCGCTCT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[1].borrow().forward_seq).unwrap(), "ATGAT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[2].borrow().forward_seq).unwrap(), "GCGC");

        simplify_structure(&mut graph, &sequences);

        assert_eq!(std::str::from_utf8(&graph.unitigs[0].borrow().forward_seq).unwrap(), "CACCGCTGCGCTCGCTTCGCTCTAT");
        assert_eq!(std::str::from_utf8(&graph.unitigs[1].borrow().forward_seq).unwrap(), "CG"); // formerly unitig 3
        assert_eq!(std::str::from_utf8(&graph.unitigs[2].borrow().forward_seq).unwrap(), "G");  // formerly unitig 2
    }

    #[test]
    fn test_check_for_duplicates() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));

        let unitigs_1 = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::FORWARD)];
        assert!(!check_for_duplicates(&unitigs_1));
        
        let unitigs_2 = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&a, strand::REVERSE)];
        assert!(check_for_duplicates(&unitigs_2));
    }

    #[test]
    fn test_merge_unitig_seqs() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let path = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&merge_unitig_seqs(&path)).unwrap(), "ACGATCAGCACTATCAGCACTACGACT");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let path = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::REVERSE), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&merge_unitig_seqs(&path)).unwrap(), "ACGATCAGCGCTGATAGTACTACGACT");
    }

    #[test]
    fn test_can_merge() {
        let (graph, seqs) = UnitigGraph::from_gfa_lines(&get_test_gfa_14());
        let (mut fixed_starts, fixed_ends) = get_fixed_unitig_starts_and_ends(&graph, &seqs);
        fix_circular_loops(&graph, &mut fixed_starts);
        assert_eq!(fixed_starts, HashSet::from([5, 8, 12, 19, 22]));
        assert_eq!(fixed_ends, HashSet::from([8, 17, 19, 22, 37]));

        assert!(cannot_merge_start(5, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(8, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(8, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(12, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(17, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(19, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(19, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(22, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(22, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_start(37, strand::REVERSE, &fixed_starts, &fixed_ends));

        assert!(cannot_merge_end(5, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(8, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(8, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(12, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(17, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(19, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(19, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(22, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(22, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(cannot_merge_end(37, strand::FORWARD, &fixed_starts, &fixed_ends));

        assert!(!cannot_merge_start(12, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(!cannot_merge_start(21, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(!cannot_merge_start(21, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(!cannot_merge_start(37, strand::FORWARD, &fixed_starts, &fixed_ends));

        assert!(!cannot_merge_end(12, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(!cannot_merge_end(21, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(!cannot_merge_end(21, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(!cannot_merge_end(37, strand::REVERSE, &fixed_starts, &fixed_ends));
    }

    #[test]
    fn test_merge_linear_paths_1() {
        let (mut graph, seqs) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
        assert_eq!(graph.unitigs.len(), 7);
        merge_linear_paths(&mut graph, &seqs);
        assert_eq!(graph.unitigs.len(), 3);
        assert_eq!(std::str::from_utf8(&graph.unitig_index.get(&8).unwrap().borrow().forward_seq).unwrap(),
                   "TTCGCTGCGCTCGCTTCGCTTTTGCACAGCGACGACGGCATGCCTGAATCGCCTA");
        assert_eq!(std::str::from_utf8(&graph.unitig_index.get(&9).unwrap().borrow().forward_seq).unwrap(),
                    "GCTCGGCTCGATGGTTCG");
        assert_eq!(std::str::from_utf8(&graph.unitig_index.get(&10).unwrap().borrow().forward_seq).unwrap(),
                    "TACTTGTAAGGC");
        let mut links = graph.get_links_for_gfa(0);
        let mut expected_links = vec![("8".to_string(), "+".to_string(), "9".to_string(), "+".to_string()),
                                      ("9".to_string(), "-".to_string(), "8".to_string(), "-".to_string()),
                                      ("9".to_string(), "+".to_string(), "9".to_string(), "-".to_string()),
                                      ("8".to_string(), "+".to_string(), "10".to_string(), "+".to_string()),
                                      ("10".to_string(), "-".to_string(), "8".to_string(), "-".to_string()),
                                      ("10".to_string(), "+".to_string(), "10".to_string(), "+".to_string()),
                                      ("10".to_string(), "-".to_string(), "10".to_string(), "-".to_string())];
        links.sort(); expected_links.sort();
        assert_eq!(links, expected_links);
    }

    #[test]
    fn test_merge_linear_paths_2() {
        let (mut graph, seqs) = UnitigGraph::from_gfa_lines(&get_test_gfa_4());
        assert_eq!(graph.unitigs.len(), 5);
        merge_linear_paths(&mut graph, &seqs);
        assert_eq!(graph.unitigs.len(), 2);
        assert_eq!(std::str::from_utf8(&graph.unitig_index.get(&6).unwrap().borrow().forward_seq).unwrap(),
                   "ACGACTACGAGCACGAGTCGTCGTCGTAACTGACT");
        assert_eq!(std::str::from_utf8(&graph.unitig_index.get(&7).unwrap().borrow().forward_seq).unwrap(),
                   "GCTCGGTG");
        let mut links = graph.get_links_for_gfa(0);
        let mut expected_links = vec![("6".to_string(), "+".to_string(), "6".to_string(), "+".to_string()),
                                      ("6".to_string(), "-".to_string(), "6".to_string(), "-".to_string()),
                                      ("7".to_string(), "+".to_string(), "7".to_string(), "+".to_string()),
                                      ("7".to_string(), "-".to_string(), "7".to_string(), "-".to_string())];
        links.sort(); expected_links.sort();
        assert_eq!(links, expected_links);
    }

    #[test]
    fn test_merge_linear_paths_3() {
        let (mut graph, seqs) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        assert_eq!(graph.unitigs.len(), 6);
        merge_linear_paths(&mut graph, &seqs);
        assert_eq!(graph.unitigs.len(), 5);
        assert_eq!(std::str::from_utf8(&graph.unitig_index.get(&7).unwrap().borrow().forward_seq).unwrap(),
                   "AAATGCGACTGTG");
    }

    #[test]
    fn test_merge_linear_paths_4() {
        let (mut graph, seqs) = UnitigGraph::from_gfa_lines(&get_test_gfa_14());
        assert_eq!(graph.unitigs.len(), 13);
        merge_linear_paths(&mut graph, &seqs);
        assert_eq!(graph.unitigs.len(), 11);
    }
}
