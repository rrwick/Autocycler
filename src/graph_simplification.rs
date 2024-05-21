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
    while expand_repeats(graph, seqs) > 0 {}
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
    let half_k = graph.k_size / 2;
    let mut total_shifted_seq = 0;
    for unitig_rc in &graph.unitigs {
        let unitig_number = unitig_rc.borrow().number;
        let inputs = get_exclusive_inputs(&unitig_rc);
        if inputs.len() >= 2 && !fixed_starts.contains(&unitig_number) {
            let can_shift = inputs.iter().all(|(input_rc, input_strand)| {
                !(*input_strand && fixed_ends.contains(&input_rc.borrow().number) ||
                  !*input_strand && fixed_starts.contains(&input_rc.borrow().number))});
            if can_shift {
                total_shifted_seq += shift_sequence_1(&inputs, &unitig_rc, half_k);
            }
        }
        let outputs = get_exclusive_outputs(&unitig_rc);
        if outputs.len() >= 2 && !fixed_ends.contains(&unitig_number) {
            let can_shift = outputs.iter().all(|(output_rc, output_strand)| {
                !(*output_strand && fixed_starts.contains(&output_rc.borrow().number) ||
                  !*output_strand && fixed_ends.contains(&output_rc.borrow().number))});
            if can_shift {
                total_shifted_seq += shift_sequence_2(&unitig_rc, &outputs, half_k);
            }
        }
    }
    total_shifted_seq
}


fn shift_sequence_1(sources: &Vec<(Rc<RefCell<Unitig>>, bool)>,
                    destination_rc: &Rc<RefCell<Unitig>>, half_k: u32) -> usize {
    // This function:
    // * removes any common sequence from the ends of the source unitigs
    // * adds that common sequence to the start of the destination unitig
    //
    // This function also guards against a couple of potential complications with sequence paths
    // (which could result in a path having more than one starting unitig):
    // * It won't let unitigs get down to a length of zero. This requires some extra logic for when
    //   both strands of one unitig appears in the sources (and will therefore be have sequence
    //   removed from both ends).
    // * It won't add sequence to the destination unitig causing any of its positions to reach the
    //   start of a path.
    //
    // The return value is the amount of sequence shifted.
    let mut common_seq = get_common_end_seq(sources);
    avoid_zero_len_unitigs(&mut common_seq, &sources, true);
    avoid_start_of_path(&mut common_seq, &destination_rc, half_k, true);
    let shifted_amount = common_seq.len();
    if shifted_amount == 0 {
        return 0;
    }
    for (source_rc, strand) in sources {
        if *strand {
            source_rc.borrow_mut().remove_seq_from_end(shifted_amount);
        } else {
            source_rc.borrow_mut().remove_seq_from_start(shifted_amount);
        }
    }
    destination_rc.borrow_mut().add_seq_to_start(common_seq);
    shifted_amount
}


fn shift_sequence_2(destination_rc: &Rc<RefCell<Unitig>>,
                    sources: &Vec<(Rc<RefCell<Unitig>>, bool)>, half_k: u32) -> usize {
    // This function does the same thing as shift_sequence_1, but for the other side of a unitig:
    // * removes any common sequence from the starts of the source unitigs
    // * adds that common sequence to the end of the destination unitig
    let mut common_seq = get_common_start_seq(sources);
    avoid_zero_len_unitigs(&mut common_seq, &sources, false);
    avoid_start_of_path(&mut common_seq, &destination_rc, half_k, false);
    let shifted_amount = common_seq.len();
    if shifted_amount == 0 {
        return 0;
    }
    for (source_rc, strand) in sources {
        if *strand {
            source_rc.borrow_mut().remove_seq_from_start(shifted_amount);
        } else {
            source_rc.borrow_mut().remove_seq_from_end(shifted_amount);
        }
    }
    destination_rc.borrow_mut().add_seq_to_end(common_seq);
    shifted_amount
}


fn avoid_zero_len_unitigs(common_seq: &mut Vec<u8>, sources: &Vec<(Rc<RefCell<Unitig>>, bool)>,
                          trim_from_start: bool) {
    // This function takes some common sequence (sequence that will be shifted from some unitigs
    // onto another) and trims it down to ensure that none of the source unitigs will end up with
    // a length of zero.
    if common_seq.len() == 0 {
        return;
    }
    let dup = if check_for_duplicates(&sources) { 2 } else { 1 };
    let min_source_len = sources.iter().map(|(source_rc, _)| source_rc.borrow().length()).min().unwrap();
    while min_source_len <= (common_seq.len() as u32) * dup {
        if trim_from_start {
            common_seq.remove(0);
        } else {
            common_seq.pop();
        }
    }
}


fn avoid_start_of_path(common_seq: &mut Vec<u8>, destination_rc: &Rc<RefCell<Unitig>>, half_k: u32,
                       trim_from_start: bool) {
    // This function takes some common sequence (sequence that will be shifted from some unitigs
    // onto another) and trims it down to ensure that the destination unitig's positions will not
    // end up at the start of a path, as this can cause problems.
    if common_seq.len() == 0 {
        return;
    }
    if trim_from_start {
        while let Some(_) = destination_rc.borrow().forward_positions.iter().find(|p| p.pos <= common_seq.len() as u32 + half_k) {
            common_seq.remove(0);
        }
    } else {
        while let Some(_) = destination_rc.borrow().reverse_positions.iter().find(|p| p.pos <= common_seq.len() as u32 + half_k) {
            common_seq.pop();
        }
    }
}


fn check_for_duplicates(unitigs: &Vec<(Rc<RefCell<Unitig>>, bool)>) -> bool {
    // Returns true if any two unitigs in the vector have the same number.
    unitigs.iter().map(|(u, _)| u.borrow().number).collect::<HashSet<_>>().len() != unitigs.len()
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


fn get_exclusive_inputs(unitig_rc: &Rc<RefCell<Unitig>>) -> Vec<(Rc<RefCell<Unitig>>, bool)> {
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
    if inputs.iter().any(|(input_rc, _)| { input_rc.borrow().number == unitig.number }) {
        return Vec::new();
    }
    inputs
}


fn get_exclusive_outputs(unitig_rc: &Rc<RefCell<Unitig>>) -> Vec<(Rc<RefCell<Unitig>>, bool)> {
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
    if outputs.iter().any(|(output_rc, _)| { output_rc.borrow().number == unitig.number }) {
        return Vec::new();
    }
    outputs
}


fn get_common_start_seq(unitigs: &Vec<(Rc<RefCell<Unitig>>, bool)>) -> Vec<u8> {
    // This function returns the common sequence at the start of all given unitigs.
    let seqs: Vec<_> = unitigs.iter().map(|(u, strand)| u.borrow().get_seq(*strand, 0, 0)).collect();
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


fn get_common_end_seq(unitigs: &Vec<(Rc<RefCell<Unitig>>, bool)>) -> Vec<u8> {
    // This function returns the common sequence at the end of all given unitigs.
    let seqs: Vec<Vec<u8>> = unitigs.iter().map(|(u, strand)| u.borrow().get_seq(*strand, 0, 0))
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
    // start/ends of such paths.
    let (fixed_starts, fixed_ends) = get_fixed_unitig_starts_and_ends(graph, seqs);
    let mut already_used = HashSet::new();
    let mut merge_paths = Vec::new();
    for unitig_rc in &graph.unitigs {
        let unitig_number = unitig_rc.borrow().number;
        for unitig_strand in [true, false] {

            // Find unitigs which can potentially start a mergeable path.
            if already_used.contains(&unitig_number) {continue;}
            if has_single_exclusive_input(&unitig_rc, unitig_strand) {continue;}
            let mut current_path = vec![(Rc::clone(unitig_rc), unitig_strand)];
            already_used.insert(unitig_number);

            // Extend the path as far as possible.
            loop {
                let &(ref unitig_rc, unitig_strand) = current_path.last().unwrap();
                if cannot_merge_end(unitig_rc.borrow().number, unitig_strand, &fixed_starts, &fixed_ends) {break;}
                let outputs = if unitig_strand {get_exclusive_outputs(&unitig_rc)} else {get_exclusive_inputs(&unitig_rc)};
                if outputs.len() != 1 {break;}
                let (output_rc, mut output_strand) = &outputs[0];
                if !unitig_strand {output_strand = !output_strand;}
                let output_number = output_rc.borrow().number;
                if already_used.contains(&output_number) {break;}
                if cannot_merge_start(output_number, output_strand, &fixed_starts, &fixed_ends) {break;}
                current_path.push((Rc::clone(output_rc), output_strand));
                already_used.insert(output_number);
            }

            if current_path.len() > 1 {
                merge_paths.push(current_path.clone());
            }
        }
    }

    let mut new_unitig_number: u32 = graph.unitigs.iter().map(|u| u.borrow().number).max().unwrap();
    for path in merge_paths {
        new_unitig_number += 1;
        merge_path(graph, &path, new_unitig_number);
    }
    graph.delete_dangling_links();
}


fn cannot_merge_start(unitig_number: u32, unitig_strand: bool, fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> bool {
    (unitig_strand && fixed_starts.contains(&unitig_number)) || (!unitig_strand && fixed_ends.contains(&unitig_number))
}


fn cannot_merge_end(unitig_number: u32, unitig_strand: bool, fixed_starts: &HashSet<u32>, fixed_ends: &HashSet<u32>) -> bool {
    (unitig_strand && fixed_ends.contains(&unitig_number)) || (!unitig_strand && fixed_starts.contains(&unitig_number))
}


fn has_single_exclusive_input(unitig_rc: &Rc<RefCell<Unitig>>, unitig_strand: bool) -> bool {
    let inputs = if unitig_strand {get_exclusive_inputs(&unitig_rc)} else {get_exclusive_outputs(&unitig_rc)};
    inputs.len() == 1
}


fn print_path(path: &Vec<(Rc<RefCell<Unitig>>, bool)>) {
    let path_str: Vec<String> = path.iter().map(|(unitig_rc, unitig_strand)| format!("{}{}", unitig_rc.borrow().number, if *unitig_strand { "+" } else { "-" })).collect();
    eprintln!("{}", path_str.join(","));
}


fn merge_path(graph: &mut UnitigGraph, path: &Vec<(Rc<RefCell<Unitig>>, bool)>, new_unitig_number: u32) {
    print_path(&path);  // TEMP

    let merged_seq = merge_unitig_seqs(path);
    let (first_unitig, first_strand) = &path[0];
    let (last_unitig, last_strand) = path.last().unwrap();
    let first_positions = if *first_strand {first_unitig.borrow().forward_positions.clone()} else {first_unitig.borrow().reverse_positions.clone()};
    let last_positions = if *last_strand {last_unitig.borrow().forward_positions.clone()} else {last_unitig.borrow().reverse_positions.clone()};

    let forward_prev = if *first_strand {first_unitig.borrow().forward_prev.clone()} else {first_unitig.borrow().reverse_prev.clone()};
    let reverse_next = if *first_strand {first_unitig.borrow().reverse_next.clone()} else {first_unitig.borrow().forward_next.clone()};
    let forward_next = if *last_strand {last_unitig.borrow().forward_next.clone()} else {last_unitig.borrow().reverse_next.clone()};
    let reverse_prev = if *last_strand {last_unitig.borrow().reverse_prev.clone()} else {last_unitig.borrow().forward_prev.clone()};

    eprintln!("{}\n", new_unitig_number); // TEMP
    let unitig = Unitig::manual(new_unitig_number, merged_seq, first_positions, last_positions,
                                forward_next, forward_prev, reverse_next, reverse_prev);

    // TODO: go through neighbouring unitigs (forward_prev, reverse_next, forward_next, reverse_prev)
    //       and create new links to the new Unitig.

    graph.unitigs.push(Rc::new(RefCell::new(unitig)));

    let path_numbers: HashSet<_> = path.iter().map(|(u, _)| u.borrow().number).collect();
    graph.unitigs.retain(|u| !path_numbers.contains(&u.borrow().number));
}


fn merge_unitig_seqs(path: &Vec<(Rc<RefCell<Unitig>>, bool)>) -> Vec<u8> {
    let total_length: usize = path.iter().map(|(unitig_rc, _)| unitig_rc.borrow().length()).sum::<u32>().try_into().unwrap();
    let mut merged_seq = Vec::with_capacity(total_length);
    for (unitig_rc, unitig_strand) in path {
        if *unitig_strand {
            merged_seq.extend(unitig_rc.borrow().forward_seq.iter());
        } else {
            merged_seq.extend(unitig_rc.borrow().reverse_seq.iter());
        }
    }
    merged_seq
}


#[cfg(test)]
mod tests {
    use std::io::Write;
    use std::fs::File;
    use std::path::PathBuf;
    use tempfile::tempdir;

    use super::*;

    fn make_test_file(file_path: &PathBuf, contents: &str) {
        let mut file = File::create(&file_path).unwrap();
        write!(file, "{}", contents).unwrap();
    }

    fn unitig_vec_to_str(mut unitigs: Vec<(Rc<RefCell<Unitig>>, bool)>) -> String {
        // Converts a vector of Unitigs to a string (makes my tests easier to write).
        unitigs.sort_by(|a, b| {
            let num_a = a.0.borrow().number;
            let num_b = b.0.borrow().number;
            num_a.cmp(&num_b).then_with(|| a.1.cmp(&b.1))
        });
        unitigs.iter()
            .map(|(unitig, strand)| {
                let num = unitig.borrow().number;
                let sign = if *strand { '+' } else { '-' };
                format!("{}{}", num, sign)
            }).collect::<Vec<String>>().join(",")
    }

    fn get_test_gfa_1() -> String {
        "S\t1\tTTCGCTGCGCTCGCTTCGCTTT\tDP:f:1\n\
        S\t2\tTGCCGTCGTCGCTGTGCA\tDP:f:1\n\
        S\t3\tTGCCTGAATCGCCTA\tDP:f:1\n\
        S\t4\tGCTCGGCTCG\tDP:f:1\n\
        S\t5\tCGAACCAT\tDP:f:1\n\
        S\t6\tTACTTGT\tDP:f:1\n\
        S\t7\tGCCTT\tDP:f:1\n\
        S\t8\tATCT\tDP:f:1\n\
        S\t9\tGC\tDP:f:1\n\
        S\t10\tT\tDP:f:1\n\
        L\t1\t+\t4\t+\t0M\n\
        L\t4\t-\t1\t-\t0M\n\
        L\t1\t+\t5\t-\t0M\n\
        L\t5\t+\t1\t-\t0M\n\
        L\t2\t+\t1\t+\t0M\n\
        L\t1\t-\t2\t-\t0M\n\
        L\t3\t-\t1\t+\t0M\n\
        L\t1\t-\t3\t+\t0M\n\
        L\t4\t+\t7\t-\t0M\n\
        L\t7\t+\t4\t-\t0M\n\
        L\t4\t+\t8\t+\t0M\n\
        L\t8\t-\t4\t-\t0M\n\
        L\t6\t-\t5\t-\t0M\n\
        L\t5\t+\t6\t+\t0M\n\
        L\t6\t+\t6\t-\t0M\n\
        L\t7\t-\t9\t+\t0M\n\
        L\t9\t-\t7\t+\t0M\n\
        L\t8\t+\t10\t-\t0M\n\
        L\t10\t+\t8\t-\t0M\n\
        L\t9\t+\t7\t+\t0M\n\
        L\t7\t-\t9\t-\t0M\n".to_string()
    }

    fn get_test_gfa_2() -> String {
        "S\t1\tACCGCTGCGCTCGCTTCGCTCT\tDP:f:1\n\
        S\t2\tATGAT\tDP:f:1\n\
        S\t3\tGCGC\tDP:f:1\n\
        L\t1\t+\t2\t+\t0M\n\
        L\t2\t-\t1\t-\t0M\n\
        L\t1\t+\t2\t-\t0M\n\
        L\t2\t+\t1\t-\t0M\n\
        L\t1\t-\t3\t+\t0M\n\
        L\t3\t-\t1\t+\t0M\n\
        L\t1\t-\t3\t-\t0M\n\
        L\t3\t+\t1\t+\t0M\n".to_string()
    }

    #[test]
    fn test_get_common_start_seq() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![(a, true), (b, true), (c, true)];
        assert_eq!(std::str::from_utf8(&get_common_start_seq(&unitigs)).unwrap(), "AC");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![(a, true), (b, true), (c, false)];
        assert_eq!(std::str::from_utf8(&get_common_start_seq(&unitigs)).unwrap(), "A");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![(a, true), (b, false), (c, false)];
        assert_eq!(std::str::from_utf8(&get_common_start_seq(&unitigs)).unwrap(), "");
    }

    #[test]
    fn test_get_common_end_seq() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![(a, true), (b, true), (c, true)];
        assert_eq!(std::str::from_utf8(&get_common_end_seq(&unitigs)).unwrap(), "");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![(a, false), (b, false), (c, true)];
        assert_eq!(std::str::from_utf8(&get_common_end_seq(&unitigs)).unwrap(), "T");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let unitigs = vec![(a, false), (b, false), (c, false)];
        assert_eq!(std::str::from_utf8(&get_common_end_seq(&unitigs)).unwrap(), "GT");
    }

    #[test]
    fn test_get_exclusive_inputs_and_outputs() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_1());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);

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
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_1());
        let (mut graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);
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
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_2());
        let (mut graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);
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

        let unitigs_1 = vec![(a.clone(), true), (b.clone(), true), (c.clone(), true)];
        assert!(!check_for_duplicates(&unitigs_1));
        
        let unitigs_2 = vec![(a.clone(), true), (b.clone(), true), (a.clone(), false)];
        assert!(check_for_duplicates(&unitigs_2));
    }

    #[test]
    fn test_merge_unitig_seqs() {
        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let path = vec![(a, true), (b, true), (c, true)];
        assert_eq!(std::str::from_utf8(&merge_unitig_seqs(&path)).unwrap(), "ACGATCAGCACTATCAGCACTACGACT");

        let a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tACGATCAGC\tDP:f:1")));
        let b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tACTATCAGC\tDP:f:1")));
        let c = Rc::new(RefCell::new(Unitig::from_segment_line("S\t3\tACTACGACT\tDP:f:1")));
        let path = vec![(a, true), (b, false), (c, true)];
        assert_eq!(std::str::from_utf8(&merge_unitig_seqs(&path)).unwrap(), "ACGATCAGCGCTGATAGTACTACGACT");
    }
}
