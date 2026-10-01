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


pub fn simplify_structure(graph: &mut UnitigGraph, seqs: &[Sequence]) {
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


fn expand_repeats(graph: &mut UnitigGraph, seqs: &[Sequence]) -> usize {
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
            let can_shift = inputs.iter().all(|input|
                !end_is_fixed(input.number(), input.strand, &fixed_starts, &fixed_ends));
            if can_shift {
                total_shifted_seq += shift_common_sequence(&inputs, unitig_rc, true);
            }
        }
        let outputs = get_exclusive_outputs(unitig_rc);
        if outputs.len() >= 2 && !fixed_ends.contains(&unitig_number) {
            let can_shift = outputs.iter().all(|output|
                !start_is_fixed(output.number(), output.strand, &fixed_starts, &fixed_ends));
            if can_shift {
                total_shifted_seq += shift_common_sequence(&outputs, unitig_rc, false);
            }
        }
    }
    total_shifted_seq
}


fn shift_common_sequence(sources: &[UnitigStrand], destination: &Rc<RefCell<Unitig>>,
                         to_start: bool) -> usize {
    let mut common_seq = get_common_seq(sources, !to_start);
    avoid_zero_len_unitigs(&mut common_seq, sources, to_start);
    avoid_start_of_path(&mut common_seq, destination, to_start);
    let shifted_amount = common_seq.len();
    if shifted_amount == 0 { return 0; }
    for source in sources {
        let source_rc = source.unitig();
        let mut unitig = source_rc.borrow_mut();
        if source.strand == to_start {
            unitig.remove_seq_from_end(shifted_amount);
        } else {
            unitig.remove_seq_from_start(shifted_amount);
        }
    }
    if to_start { destination.borrow_mut().add_seq_to_start(common_seq); }
           else { destination.borrow_mut().add_seq_to_end(common_seq); }
    shifted_amount
}


fn avoid_zero_len_unitigs(common_seq: &mut Vec<u8>, sources: &[UnitigStrand], trim_from_start: bool) {
    if common_seq.is_empty() {
        return;
    }
    // If both strands occur, sequence is removed from both ends of that unitig.
    let removals = if check_for_duplicates(sources) { 2 } else { 1 };
    let min_source_len = sources.iter().map(|source| source.length()).min().unwrap();
    while min_source_len <= (common_seq.len() as u32) * removals {
        if trim_from_start {
            common_seq.remove(0);
        } else {
            common_seq.pop();
        }
    }
}


fn avoid_start_of_path(common_seq: &mut Vec<u8>, dest_rc: &Rc<RefCell<Unitig>>,
                       trim_from_start: bool) {
    // Reaching position zero would introduce another starting unitig for the path.
    if common_seq.is_empty() {
        return;
    }
    let destination = dest_rc.borrow();
    let positions = if trim_from_start { &destination.forward_positions }
                                 else { &destination.reverse_positions };
    while positions.iter().any(|p| p.pos <= common_seq.len() as u32) {
        if trim_from_start { common_seq.remove(0); }
                      else { common_seq.pop(); }
    }
}


fn check_for_duplicates(unitigs: &[UnitigStrand]) -> bool {
    let mut seen = HashSet::new();
    unitigs.iter().any(|u| !seen.insert(u.number()))
}


fn get_fixed_unitig_starts_and_ends(graph: &UnitigGraph,
                                    sequences: &[Sequence]) -> (HashSet<u32>, HashSet<u32>) {
    // Boundaries refer to each unitig's forward strand.
    let mut fixed_starts = HashSet::new();
    let mut fixed_ends = HashSet::new();

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
    get_exclusive_neighbours(unitig_rc, false)
}


fn get_exclusive_outputs(unitig_rc: &Rc<RefCell<Unitig>>) -> Vec<UnitigStrand> {
    get_exclusive_neighbours(unitig_rc, true)
}


fn get_exclusive_neighbours(unitig_rc: &Rc<RefCell<Unitig>>, outgoing: bool) -> Vec<UnitigStrand> {
    // Every neighbour must link back exclusively to this unitig; self links are excluded.
    let unitig = unitig_rc.borrow();
    let neighbours = if outgoing { &unitig.forward_next } else { &unitig.forward_prev };
    for neighbour in neighbours {
        let Some(neighbour_rc) = neighbour.unitig.upgrade() else { return Vec::new(); };
        let neighbour_unitig = neighbour_rc.borrow();
        let links = match (outgoing, neighbour.strand) {
            (true, true) => &neighbour_unitig.forward_prev,
            (true, false) => &neighbour_unitig.reverse_prev,
            (false, true) => &neighbour_unitig.forward_next,
            (false, false) => &neighbour_unitig.reverse_next,
        };
        if links.len() != 1 || !links[0].strand || links[0].number() != unitig.number {
            return Vec::new();
        }
    }
    if neighbours.iter().any(|n| n.number() == unitig.number) { return Vec::new(); }
    neighbours.clone()
}


fn get_common_seq(unitigs: &[UnitigStrand], from_start: bool) -> Vec<u8> {
    let Some(first) = unitigs.first() else { return Vec::new(); };
    let first_rc = first.unitig();
    let first_unitig = first_rc.borrow();
    let mut common = first_unitig.get_seq(first.strand);
    for unitig in &unitigs[1..] {
        let unitig_rc = unitig.unitig();
        let borrowed = unitig_rc.borrow();
        let seq = borrowed.get_seq(unitig.strand);
        let matching = if from_start {
            common.iter().zip(seq).take_while(|(a, b)| a == b).count()
        } else {
            common.iter().rev().zip(seq.iter().rev()).take_while(|(a, b)| a == b).count()
        };
        common = if from_start { &common[..matching] } else { &common[common.len()-matching..] };
        if common.is_empty() { break; }
    }
    common.to_vec()
}


pub fn merge_linear_paths(graph: &mut UnitigGraph, seqs: &[Sequence]) {
    // Merge unbranching paths without crossing input sequence boundaries.
    let (mut fixed_starts, fixed_ends) = get_fixed_unitig_starts_and_ends(graph, seqs);
    fix_circular_loops(graph, &mut fixed_starts);
    let mut already_used = HashSet::new();
    let mut merge_paths = Vec::new();
    for unitig_rc in &graph.unitigs {
        let unitig_number = unitig_rc.borrow().number;
        for unitig_strand in [strand::FORWARD, strand::REVERSE] {
            if already_used.contains(&unitig_number) { continue; }
            if has_single_exclusive_input(unitig_rc, unitig_strand)
                && !start_is_fixed(unitig_number, unitig_strand, &fixed_starts, &fixed_ends) { continue; }
            let current_path = extend_merge_path(UnitigStrand::new(unitig_rc, unitig_strand),
                                                 &mut already_used, &fixed_starts, &fixed_ends);

            if current_path.len() > 1 {
                merge_paths.push(current_path);
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


fn extend_merge_path(start: UnitigStrand, already_used: &mut HashSet<u32>,
                     fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> Vec<UnitigStrand> {
    already_used.insert(start.number());
    let mut path = vec![start];
    loop {
        let unitig = path.last().unwrap();
        if end_is_fixed(unitig.number(), unitig.strand, fixed_starts, fixed_ends) { break; }
        let mut outputs = if unitig.strand { get_exclusive_outputs(&unitig.unitig()) }
                                     else { get_exclusive_inputs(&unitig.unitig()) };
        if outputs.len() != 1 { break; }
        let mut output = outputs.pop().unwrap();
        if !unitig.strand { output.strand = !output.strand; }
        let output_number = output.number();
        if already_used.contains(&output_number) { break; }
        if start_is_fixed(output_number, output.strand, fixed_starts, fixed_ends) { break; }
        path.push(output);
        already_used.insert(output_number);
    }
    path
}


fn fix_circular_loops(graph: &UnitigGraph, fixed_starts: &mut HashSet<u32>) {
    // Break each simple loop at its lowest numbered unitig so it can become one unitig.
    for component in graph.connected_components() {
        if graph.component_is_circular_loop(&component) {
            fixed_starts.insert(component[0]);
        }
    }
}


fn start_is_fixed(unitig_number: u32, unitig_strand: bool,
                  fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> bool {
    let fixed = if unitig_strand { fixed_starts } else { fixed_ends };
    fixed.contains(&unitig_number)
}


fn end_is_fixed(unitig_number: u32, unitig_strand: bool,
                fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> bool {
    start_is_fixed(unitig_number, !unitig_strand, fixed_starts, fixed_ends)
}


fn has_single_exclusive_input(unitig_rc: &Rc<RefCell<Unitig>>, unitig_strand: bool) -> bool {
    let inputs = if unitig_strand { get_exclusive_inputs(unitig_rc) }
                            else { get_exclusive_outputs(unitig_rc) };
    inputs.len() == 1
}


fn merge_path(graph: &mut UnitigGraph, path: &[UnitigStrand], new_unitig_number: u32) {
    let merged_seq = merge_unitig_seqs(path);
    let first = &path[0];
    let last = path.last().unwrap();
    let forward_positions = if first.strand {first.unitig().borrow().forward_positions.clone()} else {first.unitig().borrow().reverse_positions.clone()};
    let reverse_positions = if last.strand {last.unitig().borrow().reverse_positions.clone()} else {last.unitig().borrow().forward_positions.clone()};

    let end_to_start_link = graph.link_exists(last.number(), last.strand, first.number(), first.strand);
    let start_flip_link = graph.link_exists(first.number(), !first.strand, first.number(), first.strand);
    let end_flip_link = graph.link_exists(last.number(), last.strand, last.number(), !last.strand);

    let forward_prev = if first.strand {first.unitig().borrow().forward_prev.clone()} else {first.unitig().borrow().reverse_prev.clone()};
    let reverse_next = if first.strand {first.unitig().borrow().reverse_next.clone()} else {first.unitig().borrow().forward_next.clone()};
    let forward_next = if last.strand {last.unitig().borrow().forward_next.clone()} else {last.unitig().borrow().reverse_next.clone()};
    let reverse_prev = if last.strand {last.unitig().borrow().reverse_prev.clone()} else {last.unitig().borrow().forward_prev.clone()};

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

    add_reciprocal_links(&unitig_rc);
    restore_self_links(&unitig_rc, end_to_start_link, start_flip_link, end_flip_link);

    let path_numbers: HashSet<_> = path.iter().map(|u| u.number()).collect();
    graph.unitigs.retain(|u| !path_numbers.contains(&u.borrow().number));
}


fn add_reciprocal_links(unitig_rc: &Rc<RefCell<Unitig>>) {
    let unitig = unitig_rc.borrow();
    for (unitig_strand, next, prev) in [
        (strand::FORWARD, &unitig.forward_next, &unitig.forward_prev),
        (strand::REVERSE, &unitig.reverse_next, &unitig.reverse_prev),
    ] {
        for neighbour in next {
            let neighbour_rc = neighbour.unitig();
            let mut neighbour_unitig = neighbour_rc.borrow_mut();
            let links = if neighbour.strand { &mut neighbour_unitig.forward_prev }
                                      else { &mut neighbour_unitig.reverse_prev };
            links.push(UnitigStrand::new(unitig_rc, unitig_strand));
        }
        for neighbour in prev {
            let neighbour_rc = neighbour.unitig();
            let mut neighbour_unitig = neighbour_rc.borrow_mut();
            let links = if neighbour.strand { &mut neighbour_unitig.forward_next }
                                      else { &mut neighbour_unitig.reverse_next };
            links.push(UnitigStrand::new(unitig_rc, unitig_strand));
        }
    }
}


fn restore_self_links(unitig_rc: &Rc<RefCell<Unitig>>, end_to_start_link: bool,
                      start_flip_link: bool, end_flip_link: bool) {
    if end_to_start_link {
        let mut u = unitig_rc.borrow_mut();
        u.forward_next.push(UnitigStrand::new(unitig_rc, strand::FORWARD));
        u.forward_prev.push(UnitigStrand::new(unitig_rc, strand::FORWARD));
        u.reverse_next.push(UnitigStrand::new(unitig_rc, strand::REVERSE));
        u.reverse_prev.push(UnitigStrand::new(unitig_rc, strand::REVERSE));
    }
    if start_flip_link {
        let mut u = unitig_rc.borrow_mut();
        u.reverse_next.push(UnitigStrand::new(unitig_rc, strand::FORWARD));
        u.forward_prev.push(UnitigStrand::new(unitig_rc, strand::REVERSE));
    }
    if end_flip_link {
        let mut u = unitig_rc.borrow_mut();
        u.forward_next.push(UnitigStrand::new(unitig_rc, strand::REVERSE));
        u.reverse_prev.push(UnitigStrand::new(unitig_rc, strand::FORWARD));
    }
}


fn merge_unitig_seqs(path: &[UnitigStrand]) -> Vec<u8> {
    let total_length: usize = path.iter().map(|u| u.length()).sum::<u32>().try_into().unwrap();
    let mut merged_seq = Vec::with_capacity(total_length);
    for u in path {
        let unitig_rc = u.unitig();
        let unitig = unitig_rc.borrow();
        merged_seq.extend_from_slice(unitig.get_seq(u.strand));
    }
    merged_seq
}


fn get_merge_path_depth(path: &[UnitigStrand], forward_positions: &[Position]) -> f64 {
    if !forward_positions.is_empty() {
        return forward_positions.len() as f64;
    }

    for u in path {
        if u.is_anchor() {
            return u.depth();
        }
    }

    weighted_mean_depth(path)
}


fn weighted_mean_depth(path: &[UnitigStrand]) -> f64 {
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
        unitigs.sort_by_key(|u| (u.number(), u.strand));
        unitigs.iter().map(ToString::to_string).collect::<Vec<_>>().join(",")
    }

    #[test]
    fn test_common_sequence_boundaries() {
        for (seqs, prefix, suffix) in [
            (vec![], "", ""), (vec![""], "", ""), (vec!["ACGT"], "ACGT", "ACGT"),
            (vec!["ACGT", "ACGT"], "ACGT", "ACGT"),
            (vec!["ACGT", "AC"], "AC", ""), (vec!["GT", "ACGT"], "", "GT"),
            (vec!["ACGT", ""], "", ""), (vec!["", "ACGT"], "", ""),
        ] {
            let unitigs: Vec<_> = seqs.iter().enumerate().map(|(i, seq)| {
                Rc::new(RefCell::new(Unitig::from_segment_line(&format!("S\t{i}\t{seq}\tDP:f:1"))))
            }).collect();
            let strands: Vec<_> = unitigs.iter().map(|u| UnitigStrand::new(u, true)).collect();
            assert_eq!(get_common_seq(&strands, true), prefix.as_bytes());
            assert_eq!(get_common_seq(&strands, false), suffix.as_bytes());
        }
    }

    #[test]
    fn test_shift_sequence_preserves_paths_and_nonempty_sources() {
        for to_start in [true, false] {
            for (source_seq, positions, expected_shift) in [
                ("AT", vec![], 0), ("AATT", vec![], 1), ("AAATTT", vec![], 2),
                ("AAATTT", vec![1], 0), ("AAATTT", vec![5, 2], 1),
            ] {
                let source = Rc::new(RefCell::new(Unitig::from_segment_line(
                    &format!("S\t1\t{source_seq}\tDP:f:1"))));
                let sources = vec![UnitigStrand::new(&source, true), UnitigStrand::new(&source, false)];
                let destination = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tCG\tDP:f:1")));
                let path_positions = positions.iter().map(|&p| Position::new(1, true, p)).collect();
                if to_start { destination.borrow_mut().forward_positions = path_positions; }
                       else { destination.borrow_mut().reverse_positions = path_positions; }
                let shifted = shift_common_sequence(&sources, &destination, to_start);
                assert_eq!(shifted, expected_shift);
                assert_eq!(source.borrow().forward_seq, source_seq.as_bytes()[shifted..source_seq.len()-shifted]);
                let expected_dest = if to_start { format!("{}CG", "T".repeat(shifted)) }
                                          else { format!("CG{}", "A".repeat(shifted)) };
                let destination = destination.borrow();
                assert_eq!(destination.forward_seq, expected_dest.as_bytes());
                assert_eq!(destination.reverse_seq, reverse_complement(expected_dest.as_bytes()));
                let shifted_positions = if to_start { &destination.forward_positions }
                                               else { &destination.reverse_positions };
                assert_eq!(shifted_positions.iter().map(|p| p.pos as usize).collect::<Vec<_>>(),
                           positions.iter().map(|p| p - shifted).collect::<Vec<_>>());
            }
        }
    }

    #[test]
    fn test_get_common_start_seq() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&get_common_seq(&unitigs, true)).unwrap(), "AC");

        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::REVERSE)];
        assert_eq!(std::str::from_utf8(&get_common_seq(&unitigs, true)).unwrap(), "A");

        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::REVERSE), UnitigStrand::new(&c, strand::REVERSE)];
        assert_eq!(std::str::from_utf8(&get_common_seq(&unitigs, true)).unwrap(), "");
    }

    #[test]
    fn test_get_common_end_seq() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![UnitigStrand::new(&a, strand::FORWARD), UnitigStrand::new(&b, strand::FORWARD), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&get_common_seq(&unitigs, false)).unwrap(), "");

        let unitigs = vec![UnitigStrand::new(&a, strand::REVERSE), UnitigStrand::new(&b, strand::REVERSE), UnitigStrand::new(&c, strand::FORWARD)];
        assert_eq!(std::str::from_utf8(&get_common_seq(&unitigs, false)).unwrap(), "T");

        let unitigs = vec![UnitigStrand::new(&a, strand::REVERSE), UnitigStrand::new(&b, strand::REVERSE), UnitigStrand::new(&c, strand::REVERSE)];
        assert_eq!(std::str::from_utf8(&get_common_seq(&unitigs, false)).unwrap(), "GT");
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

        assert!(start_is_fixed(5, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(8, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(8, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(12, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(17, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(19, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(19, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(22, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(22, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(start_is_fixed(37, strand::REVERSE, &fixed_starts, &fixed_ends));

        assert!(end_is_fixed(5, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(8, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(8, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(12, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(17, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(19, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(19, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(22, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(22, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(end_is_fixed(37, strand::FORWARD, &fixed_starts, &fixed_ends));

        assert!(!start_is_fixed(12, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(!start_is_fixed(21, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(!start_is_fixed(21, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(!start_is_fixed(37, strand::FORWARD, &fixed_starts, &fixed_ends));

        assert!(!end_is_fixed(12, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(!end_is_fixed(21, strand::FORWARD, &fixed_starts, &fixed_ends));
        assert!(!end_is_fixed(21, strand::REVERSE, &fixed_starts, &fixed_ends));
        assert!(!end_is_fixed(37, strand::REVERSE, &fixed_starts, &fixed_ends));
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
        let mut expected_links = vec![(8, "+", 9, "+"),
                                      (9, "-", 8, "-"),
                                      (9, "+", 9, "-"),
                                      (8, "+", 10, "+"),
                                      (10, "-", 8, "-"),
                                      (10, "+", 10, "+"),
                                      (10, "-", 10, "-")];
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
        let mut expected_links = vec![(6, "+", 6, "+"),
                                      (6, "-", 6, "-"),
                                      (7, "+", 7, "+"),
                                      (7, "-", 7, "-")];
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
