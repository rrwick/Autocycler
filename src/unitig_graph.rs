// This file defines the UnitigGraph struct for building a compacted unitig graph from a KmerGraph.

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
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write, BufRead, BufReader};
use std::path::PathBuf;
use std::rc::Rc;

use crate::kmer_graph::KmerGraph;
use crate::position::Position;
use crate::sequence::Sequence;
use crate::unitig::{Unitig, UnitigStrand};
use crate::misc::{reverse_complement, quit_with_error, strand, load_file_lines};


pub struct UnitigGraph {
    pub unitigs: Vec<Rc<RefCell<Unitig>>>,
    pub k_size: u32,
    pub unitig_index: HashMap<u32, Rc<RefCell<Unitig>>>,
}

impl UnitigGraph {
    pub fn from_kmer_graph(k_graph: &KmerGraph) -> Self {
        let mut u_graph = UnitigGraph {
            unitigs: Vec::new(),
            k_size: k_graph.k_size,
            unitig_index: HashMap::new(),
        };
        u_graph.build_unitigs_from_kmer_graph(k_graph);
        u_graph.simplify_seqs();
        u_graph.create_links();
        u_graph.trim_overlaps();
        u_graph.renumber_unitigs();
        u_graph.check_links();
        u_graph
    }

    pub fn from_gfa_file(gfa_filename: &PathBuf) -> (Self, Vec<Sequence>) {
        let gfa_lines = load_file_lines(gfa_filename);
        Self::from_gfa_lines(&gfa_lines)
    }

    pub fn from_gfa_lines(gfa_lines: &Vec<String>) -> (Self, Vec<Sequence>) {
        let mut u_graph = UnitigGraph {
            unitigs: Vec::new(),
            k_size: 0,
            unitig_index: HashMap::new(),
        };
        let mut link_lines: Vec<&str> = Vec::new();
        let mut path_lines: Vec<&str> = Vec::new();
        for line in gfa_lines {
            let parts: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
            match parts.get(0) {
                Some(&"H") => u_graph.read_gfa_header_line(&parts),
                Some(&"S") => u_graph.unitigs.push(Rc::new(RefCell::new(Unitig::from_segment_line(&line)))),
                Some(&"L") => link_lines.push(line),
                Some(&"P") => path_lines.push(line),
                _ => {}
            }
        }
        u_graph.build_unitig_index();
        u_graph.build_links_from_gfa(&link_lines);
        let sequences = u_graph.build_paths_from_gfa(&path_lines);
        u_graph.check_links();
        (u_graph, sequences)
    }

    pub fn build_unitig_index(&mut self) {
        self.unitig_index = self.unitigs.iter().map(|u| {(u.borrow().number, Rc::clone(u))}).collect();
    }

    fn read_gfa_header_line(&mut self, parts: &Vec<&str>) {
        for &p in parts {
            if p.starts_with("KM:i:") {
                if let Ok(k) = p[5..].parse::<u32>() {
                    self.k_size = k;
                    return;
                }
            }
        }
        quit_with_error("could not find a valid k-mer tag (e.g. KM:i:51) in the GFA header line.\n\
                         Are you sure this is an Autocycler-generated GFA file?");
    }

    fn build_links_from_gfa(&mut self, link_lines: &[&str]) {
        for line in link_lines {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 6 || parts[5] != "0M" {
                quit_with_error("non-zero overlap found on the GFA link line.\n\
                                 Are you sure this is an Autocycler-generated GFA file?");
            }
            let seg_1: u32 = parts[1].parse().expect("Error parsing segment 1 as integer");
            let seg_2: u32 = parts[3].parse().expect("Error parsing segment 2 as integer");
            let strand_1 = parts[2] == "+";
            let strand_2 = parts[4] == "+";
            if let Some(unitig_1) = self.unitig_index.get(&seg_1) {
                if let Some(unitig_2) = self.unitig_index.get(&seg_2) {
                    if strand_1 {unitig_1.borrow_mut().forward_next.push(UnitigStrand::new(unitig_2, strand_2));
                         } else {unitig_1.borrow_mut().reverse_next.push(UnitigStrand::new(unitig_2, strand_2));}
                    if strand_2 {unitig_2.borrow_mut().forward_prev.push(UnitigStrand::new(unitig_1, strand_1));
                         } else {unitig_2.borrow_mut().reverse_prev.push(UnitigStrand::new(unitig_1, strand_1));}
                } else {
                    quit_with_error(&format!("link refers to nonexistent unitig: {}", seg_2));
                }
            } else {
                quit_with_error(&format!("link refers to nonexistent unitig: {}", seg_1));
            }
        }
    }

    fn build_paths_from_gfa(&mut self, path_lines: &[&str]) -> Vec<Sequence> {
        let mut sequences = Vec::new();
        for line in path_lines {
            let parts: Vec<&str> = line.split('\t').collect();
            let seq_id: u16 = parts[1].parse().expect("Error parsing sequence ID as integer");
            let mut length = None;
            let mut filename = None;
            let mut header = None;
            let mut cluster = 0;
            let mut extend = true;
            for p in &parts[2..] {
                if p.starts_with("LN:i:") {
                    length = Some(p[5..].parse::<u32>().expect("Error parsing length"));
                } else if p.starts_with("FN:Z:") {
                    filename = Some(p[5..].to_string());
                } else if p.starts_with("HD:Z:") {
                    header = Some(p[5..].to_string());
                } else if p.starts_with("CL:i:") {
                    cluster = p[5..].parse::<u16>().expect("Error parsing cluster");
                } else if *p == "EX:Z:false" {
                    extend = false;
                }
            }
            if length.is_none() || filename.is_none() || header.is_none() {
                quit_with_error("missing required tag in GFA path line.");
            }
            let length = length.unwrap();
            let filename = filename.unwrap();
            let header = header.unwrap();
            let path = parse_unitig_path(parts[2]);
            let sequence = self.create_sequence_and_positions(seq_id, length, filename, header, cluster, extend, path);
            sequences.push(sequence);
        }
        sequences
    }

    pub fn create_sequence_and_positions(&mut self, seq_id: u16, length: u32,
                                         filename: String, header: String, cluster: u16, extend: bool,
                                         forward_path: Vec<(u32, bool)>) -> Sequence {
        let reverse_path = reverse_path(&forward_path);
        self.add_positions_from_path(&forward_path, strand::FORWARD, seq_id, length, extend);
        self.add_positions_from_path(&reverse_path, strand::REVERSE, seq_id, length, extend);
        Sequence::new_without_seq(seq_id, filename, header, length as usize, cluster, extend)
    }

    fn add_positions_from_path(&mut self, path: &[(u32, bool)], path_strand: bool, seq_id: u16,
                               length: u32, extend: bool) {
        let half_k = self.k_size / 2;
        let mut pos = half_k;
        for (unitig_num, unitig_strand) in path {
            if let Some(unitig) = self.unitig_index.get(unitig_num) {
                let mut u = unitig.borrow_mut();
                let positions = if *unitig_strand {&mut u.forward_positions} 
                                             else {&mut u.reverse_positions};
                positions.push(Position::new(seq_id, path_strand, pos as usize));
                pos += u.length();
                if u.dead_end_start(*unitig_strand) {pos -= half_k;}
                if u.dead_end_end(*unitig_strand) {pos -= half_k;}
            } else {
                quit_with_error(&format!("unitig {} not found in unitig index", unitig_num));
            }
        }
        if extend {
            assert!(pos + half_k == length, "Position calculation mismatch");
        } else {
            assert!(pos - half_k == length, "Position calculation mismatch");
        }
    }

    fn build_unitigs_from_kmer_graph(&mut self, k_graph: &KmerGraph) {
        let mut seen: HashSet<&[u8]> = HashSet::new();
        let mut unitig_number = 0;
        let half_k = self.k_size as usize / 2;
        for forward_kmer in k_graph.iterate_kmers() {
            if seen.contains(forward_kmer.seq()) {
                continue;
            }
            let reverse_kmer = k_graph.reverse(forward_kmer);
            unitig_number += 1;
            let mut unitig = Unitig::from_kmers(unitig_number, forward_kmer, &reverse_kmer);
            seen.insert(forward_kmer.seq());
            seen.insert(reverse_kmer.seq());

            // Extend unitig forward
            let mut for_k = forward_kmer;
            let mut rev_k = reverse_kmer;
            loop {
                if rev_k.first_position(half_k) { break; }
                let next_kmers = k_graph.next_kmers(for_k.seq());
                if next_kmers.len() != 1 { break; }
                for_k = &next_kmers[0];
                if seen.contains(for_k.seq()) { break; }
                let prev_kmers = k_graph.prev_kmers(for_k.seq());
                if prev_kmers.len() != 1 { break; }
                rev_k = k_graph.reverse(for_k);
                if for_k.first_position(half_k) { break; }
                unitig.add_kmer_to_end(for_k, rev_k);
                seen.insert(for_k.seq());
                seen.insert(rev_k.seq());
            }

            // Extend unitig backward
            let mut for_k = forward_kmer;
            let mut rev_k;
            loop {
                if for_k.first_position(half_k) { break; }
                let prev_kmers = k_graph.prev_kmers(for_k.seq());
                if prev_kmers.len() != 1 { break; }
                for_k = &prev_kmers[0];
                if seen.contains(for_k.seq()) { break; }
                let next_kmers = k_graph.next_kmers(for_k.seq());
                if next_kmers.len() != 1 { break; }
                rev_k = k_graph.reverse(for_k);
                if rev_k.first_position(half_k) { break; }
                unitig.add_kmer_to_start(for_k, rev_k);
                seen.insert(for_k.seq());
                seen.insert(rev_k.seq());
            }
            self.unitigs.push(Rc::new(RefCell::new(unitig)));
        }
    }

    fn simplify_seqs(&mut self) {
        for unitig in &self.unitigs {
            unitig.borrow_mut().simplify_seqs();
        }
    }

    fn create_links(&mut self) {
        let piece_len = self.k_size as usize - 1;

        // Index unitigs by their k-1 starting sequences.
        let mut forward_starts = HashMap::new();
        let mut reverse_starts = HashMap::new();
        for (i, unitig) in self.unitigs.iter().enumerate() {
            let forward_key = unitig.borrow().forward_seq[..piece_len].to_vec();
            let reverse_key = unitig.borrow().reverse_seq[..piece_len].to_vec();
            forward_starts.entry(forward_key).or_insert_with(Vec::new).push(i);
            reverse_starts.entry(reverse_key).or_insert_with(Vec::new).push(i);
        }

        // Use the indices to find connections between unitigs.
        for i in 0..self.unitigs.len() {
            let unitig_a = Rc::clone(&self.unitigs[i]);
            let ending_forward_seq = unitig_a.borrow().forward_seq[unitig_a.borrow().forward_seq.len() - piece_len..].to_vec();
            let ending_reverse_seq = unitig_a.borrow().reverse_seq[unitig_a.borrow().reverse_seq.len() - piece_len..].to_vec();

            if let Some(next_idxs) = forward_starts.get(&ending_forward_seq) {
                for &j in next_idxs {
                    let unitig_b = Rc::clone(&self.unitigs[j]);

                    // unitig_a+ -> unitig_b+
                    unitig_a.borrow_mut().forward_next.push(UnitigStrand::new(&unitig_b, strand::FORWARD));
                    unitig_b.borrow_mut().forward_prev.push(UnitigStrand::new(&unitig_a, strand::FORWARD));

                    // unitig_b- -> unitig_a-
                    unitig_b.borrow_mut().reverse_next.push(UnitigStrand::new(&unitig_a, strand::REVERSE));
                    unitig_a.borrow_mut().reverse_prev.push(UnitigStrand::new(&unitig_b, strand::REVERSE));
                }
            }

            if let Some(next_idxs) = reverse_starts.get(&ending_forward_seq) {
                for &j in next_idxs {
                    let unitig_b = Rc::clone(&self.unitigs[j]);

                    // unitig_a+ -> unitig_b-
                    unitig_a.borrow_mut().forward_next.push(UnitigStrand::new(&unitig_b, strand::REVERSE));
                    unitig_b.borrow_mut().reverse_prev.push(UnitigStrand::new(&unitig_a, strand::FORWARD));
                }
            }

            if let Some(next_idxs) = forward_starts.get(&ending_reverse_seq) {
                for &j in next_idxs {
                    let unitig_b = Rc::clone(&self.unitigs[j]);

                    // unitig_a- -> unitig_b+
                    unitig_a.borrow_mut().reverse_next.push(UnitigStrand::new(&unitig_b, strand::FORWARD));
                    unitig_b.borrow_mut().forward_prev.push(UnitigStrand::new(&unitig_a, strand::REVERSE));
                }
            }
        }
    }

    pub fn trim_overlaps(&mut self) {
        // This method trims all overlaps from the graph. Dead-ends will not be trimmed, as that
        // would cause irreversible loss of sequence.
        for unitig in &self.unitigs {
            unitig.borrow_mut().trim_overlaps(self.k_size as usize);
        }
    }

    pub fn restore_overlaps(&mut self) {
        // This method puts all overlaps back onto the graph. When deleting unitigs, it is necessary
        // to first restore overlaps with this method, then delete the unitigs, then trim overlaps
        // again. This is because deleting unitigs can create new dead-ends which changes the
        // trimming.
        let overlap = self.k_size as usize / 2;
        let new_seqs: HashMap<_, _> = self.unitigs.iter().map(|rc| {let u = rc.borrow(); (u.number, u.get_extended_seq(strand::FORWARD, overlap, overlap))}).collect();
        for rc in &self.unitigs {
            let mut u = rc.borrow_mut();
            u.forward_seq = new_seqs.get(&u.number).unwrap().clone();
            u.reverse_seq = reverse_complement(&u.forward_seq);
            u.trimmed = false;
        }
    }

    pub fn renumber_unitigs(&mut self) {
        // This method sorts and renumbers Unitigs by: length (decreasing), sequence (lexicographic)
        // and depth (decreasing).
        self.unitigs.sort_by(|a_rc, b_rc| {
            let a = a_rc.borrow();
            let b = b_rc.borrow();
            let length_cmp = a.length().cmp(&b.length()).reverse();
            if length_cmp != std::cmp::Ordering::Equal {
                return length_cmp;
            }
            let seq_cmp = a.forward_seq.cmp(&b.forward_seq);
            if seq_cmp != std::cmp::Ordering::Equal {
                return seq_cmp;
            }
            a.depth.partial_cmp(&b.depth).unwrap_or(std::cmp::Ordering::Equal).reverse()
        });
        for (new_number, unitig) in self.unitigs.iter().enumerate() {
            unitig.borrow_mut().number = (new_number + 1) as u32;
        }
        self.build_unitig_index();
    }

    pub fn save_gfa(&self, gfa_filename: &PathBuf, sequences: &Vec<Sequence>) -> io::Result<()> {
        let mut file = File::create(gfa_filename)?;
        writeln!(file, "H\tVN:Z:1.0\tKM:i:{}", self.k_size)?;
        for unitig in &self.unitigs {
            writeln!(file, "{}", unitig.borrow().gfa_segment_line())?;
        }
        for (a, a_strand, b, b_strand) in self.get_links_for_gfa() {
            writeln!(file, "L\t{}\t{}\t{}\t{}\t0M", a, a_strand, b, b_strand)?;
        }
        for s in sequences {
            writeln!(file, "{}", self.get_gfa_path_line(&s))?;
        }
        Ok(())
    }

    pub fn get_links_for_gfa(&self) -> Vec<(String, String, String, String)> {
        let mut links = Vec::new();
        for a_rc in &self.unitigs {
            let a = a_rc.borrow();
            for b in &a.forward_next {
                links.push((a.number.to_string(), "+".to_string(), b.number().to_string(),
                            (if b.strand {"+"} else {"-"}).to_string()));
            }
            for b in &a.reverse_next {
                links.push((a.number.to_string(), "-".to_string(), b.number().to_string(),
                            (if b.strand {"+"} else {"-"}).to_string()));
            }
        }
        links
    }

    fn get_gfa_path_line(&self, seq: &Sequence) -> String {
        let unitig_path = self.get_unitig_path_for_sequence(seq);
        let path_str: Vec<String> = unitig_path.iter()
            .map(|(num, strand)| format!("{}{}", num, if *strand { "+" } else { "-" })).collect();
        let path_str = path_str.join(",");
        let cluster_tag = if seq.cluster > 0 {format!("\tCL:i:{}", seq.cluster)} else {"".to_string()};
        let extend_tag = if !seq.extend {"\tEX:Z:false"} else {""};
        format!("P\t{}\t{}\t*\tLN:i:{}\tFN:Z:{}\tHD:Z:{}{}{}",
                seq.id, path_str, seq.length, seq.filename, seq.contig_header, cluster_tag, extend_tag)
    }

    pub fn reconstruct_original_sequences(&self, seqs: &Vec<Sequence>) -> HashMap<String, Vec<(String, String)>> {
        let mut original_seqs: HashMap<String, Vec<(String, String)>> = HashMap::new();
        for seq in seqs {
            let (filename, header, sequence) = self.reconstruct_original_sequence(&seq);
            original_seqs.entry(filename).or_insert_with(Vec::new).push((header, sequence));
        }
        original_seqs
    }

    fn reconstruct_original_sequence(&self, seq: &Sequence) -> (String, String, String) {
        eprintln!("  {}: {} ({} bp)", seq.filename, seq.contig_name(), seq.length);
        let path = self.get_unitig_path_for_sequence(&seq);
        let sequence = self.get_sequence_from_path(&path, seq.extend);
        assert_eq!(sequence.len(), seq.length, "reconstructed sequence does not have expected length");
        (seq.filename.clone(), seq.contig_header.clone(), sequence)
    }

    fn get_sequence_from_path(&self, path: &Vec<(u32, bool)>, extend: bool) -> String {
        // Given a path (vector of unitig IDs and strands), this function returns the sequence
        // traced by that path. It also requires a unitig index so it can quickly look up unitigs
        // by their number.
        let half_k = (self.k_size / 2) as usize;
        let mut sequence = Vec::new();
        for (i, (unitig_num, strand)) in path.iter().enumerate() {
            let unitig = self.unitig_index.get(unitig_num).unwrap();
            let upstream = if !extend { 0 } else if i == 0 { half_k } else { 0 };
            let downstream = if !extend { 0 } else if i == path.len() - 1 { half_k } else { 0 };
            sequence.push(String::from_utf8(unitig.borrow().get_extended_seq(*strand, upstream, downstream)).unwrap());
        }
        sequence.into_iter().collect()
    }

    fn find_starting_unitig(&self, seq_id: u16) -> UnitigStrand {
        // For a given sequence ID, this function returns the Unitig and strand where that sequence
        // begins.
        let half_k = self.k_size / 2;
        let mut starting_unitigs = Vec::new();
        for unitig in &self.unitigs {
            for p in &unitig.borrow().forward_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == half_k {
                    starting_unitigs.push(UnitigStrand::new(unitig, strand::FORWARD));
                }
            }
            for p in &unitig.borrow().reverse_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == half_k {
                    starting_unitigs.push(UnitigStrand::new(unitig, strand::REVERSE));
                }
            }
        }
        assert_eq!(starting_unitigs.len(), 1);
        starting_unitigs[0].clone()
    }

    fn get_next_unitig(&self, seq_id: u16, unitig_rc: &Rc<RefCell<Unitig>>, strand: bool,
                       pos: u32) -> Option<(UnitigStrand, u32)> {
        // For a given unitig that's part of a sequence's path, this function will return the next
        // unitig in that sequence's path.
        let unitig = unitig_rc.borrow();
        let next_pos = pos + unitig.untrimmed_length(self.k_size as usize) - self.k_size as u32 + 1;
        let next_unitigs = if strand { &unitig.forward_next } else { &unitig.reverse_next };
        for next in next_unitigs {
            let u = next.unitig.borrow();
            let positions = if next.strand { &u.forward_positions } else { &u.reverse_positions};
            for p in positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == next_pos {
                    return Some((UnitigStrand::new(&next.unitig, next.strand), next_pos));
                }
            }
        }
        None
    }

    pub fn get_unitig_path_for_sequence(&self, seq: &Sequence) -> Vec<(u32, bool)> {
        let half_k = self.k_size / 2;
        let mut unitig_path = Vec::new();
        let mut u = self.find_starting_unitig(seq.id);
        let mut pos = half_k;
        loop {
            unitig_path.push((u.number(), u.strand));
            match self.get_next_unitig(seq.id, &u.unitig, u.strand, pos) {
                None => break,
                Some((next, next_pos)) => {
                    (u, pos) = (next, next_pos);
                }
            }
        }
        unitig_path
    }

    pub fn get_total_length(&self) -> u32 {
        self.unitigs.iter().map(|u| u.borrow().length()).sum()
    }

    pub fn get_link_count(&self) -> u32 {
        let mut link_count = 0;
        for unitig in &self.unitigs {
            link_count += unitig.borrow().forward_next.len();
            link_count += unitig.borrow().reverse_next.len();
        }
        link_count.try_into().unwrap()
    }

    pub fn print_basic_graph_info(&self) {
        eprintln!("{} unitigs, {} links", self.unitigs.len(), self.get_link_count());
        eprintln!("total length: {} bp", self.get_total_length());
        eprintln!();
    }

    pub fn delete_dangling_links(&mut self) {
        // This method deletes any links to no-longer-existing unitigs. It should be run after any
        // code which deletes Unitigs from the graph.
        let unitig_numbers: HashSet<u32> = self.unitigs.iter().map(|u| u.borrow().number).collect();
        for unitig_rc in &self.unitigs {
            let unitig = unitig_rc.borrow();
            let forward_next_to_remove = unitig.forward_next.iter().enumerate().filter_map(|(index, u)| {if !unitig_numbers.contains(&u.number()) {Some(index)} else {None}}).collect::<Vec<_>>();
            let forward_prev_to_remove = unitig.forward_prev.iter().enumerate().filter_map(|(index, u)| {if !unitig_numbers.contains(&u.number()) {Some(index)} else {None}}).collect::<Vec<_>>();
            let reverse_next_to_remove = unitig.reverse_next.iter().enumerate().filter_map(|(index, u)| {if !unitig_numbers.contains(&u.number()) {Some(index)} else {None}}).collect::<Vec<_>>();
            let reverse_prev_to_remove = unitig.reverse_prev.iter().enumerate().filter_map(|(index, u)| {if !unitig_numbers.contains(&u.number()) {Some(index)} else {None}}).collect::<Vec<_>>();
            drop(unitig);
            let mut unitig = unitig_rc.borrow_mut();
            for index in forward_next_to_remove.into_iter().rev() { unitig.forward_next.remove(index); }
            for index in forward_prev_to_remove.into_iter().rev() { unitig.forward_prev.remove(index); }
            for index in reverse_next_to_remove.into_iter().rev() { unitig.reverse_next.remove(index); }
            for index in reverse_prev_to_remove.into_iter().rev() { unitig.reverse_prev.remove(index); }
        }
    }

    pub fn remove_sequence_from_graph(&mut self, seq_id: u16) {
        // Removes all Positions from the Unitigs which have the given sequence ID. This reduces
        // depths of affected Unitigs, and can result in zero-depth unitigs, so it may be necessary
        // to run remove_zero_depth_unitigs after this.
        for u in &self.unitigs {
            u.borrow_mut().remove_sequence(seq_id);
        }
    }

    pub fn recalculate_depths(&mut self) {
        // Sets each unitig's depth based on its Positions. Useful after adding/removing paths.
        for u in &self.unitigs {
            u.borrow_mut().recalculate_depth();
        }
    }

    pub fn remove_zero_depth_unitigs(&mut self) {
        // Removes zero-depth unitigs from the graph. Doing so can create new dead-ends, so this
        // function first un-trims the contigs (adds overlap back on) and then re-trims after the
        // unitigs have been deleted.
        self.restore_overlaps();
        self.unitigs.retain(|u| u.borrow().depth > 0.0);
        self.delete_dangling_links();
        self.trim_overlaps();
        self.build_unitig_index();
    }

    pub fn link_exists(&self, a_num: u32, a_strand: bool, b_num: u32, b_strand: bool) -> bool {
        // Checks if the given link exists (looks for it in forward_next/reverse_next).
        if let Some(unitig_a) = self.unitig_index.get(&a_num) {
            let unitig_a = unitig_a.borrow();
            let next_links = if a_strand {&unitig_a.forward_next} else {&unitig_a.reverse_next};
            for next in next_links {
                if next.number() == b_num && next.strand == b_strand {
                    return true;
                }
            }
        }
        false
    }

    pub fn link_exists_prev(&self, a_num: u32, a_strand: bool, b_num: u32, b_strand: bool) -> bool {
        // This is like the link_exists method, but it checks in the opposite direction (looks for
        // it in forward_prev/reverse_prev).
        if let Some(unitig_b) = self.unitig_index.get(&b_num) {
            let unitig_b = unitig_b.borrow();
            let prev_links = if b_strand {&unitig_b.forward_prev} else {&unitig_b.reverse_prev};
            for prev in prev_links {
                if prev.number() == a_num && prev.strand == a_strand {
                    return true;
                }
            }
        }
        false
    }

    pub fn check_links(&self) {
        // Makes sure that all of the graph's links are valid:
        // * Each link should have a corresponding link on the opposite strand.
        // * Each next link should be matched with a prev link.
        // * All linked Unitigs should be in the unitig_index.
        // If any of the above aren't true, this method will panic.
        for a_rc in &self.unitigs {
            let a = a_rc.borrow();
            for b in &a.forward_next {
                let a_strand = strand::FORWARD;
                if !self.link_exists(a.number, a_strand, b.number(), b.strand) {panic!("missing next link");}
                if !self.link_exists_prev(a.number, a_strand, b.number(), b.strand) {panic!("missing prev link");}
                if !self.link_exists(b.number(), !b.strand, a.number, !a_strand) {panic!("missing next link");}
                if !self.link_exists_prev(b.number(), !b.strand, a.number, !a_strand) {panic!("missing prev link");}
                if !self.unitig_index.contains_key(&b.number()) {panic!("unitig missing from index");}
            }
            for b in &a.reverse_next {
                let a_strand = strand::REVERSE;
                if !self.link_exists(a.number, a_strand, b.number(), b.strand) {panic!("missing next link");}
                if !self.link_exists_prev(a.number, a_strand, b.number(), b.strand) {panic!("missing prev link");}
                if !self.link_exists(b.number(), !b.strand, a.number, !a_strand) {panic!("missing next link");}
                if !self.link_exists_prev(b.number(), !b.strand, a.number, !a_strand) {panic!("missing prev link");}
                if !self.unitig_index.contains_key(&b.number()) {panic!("unitig missing from index");}
            }
            for b in &a.forward_prev {
                let a_strand = strand::FORWARD;
                if !self.link_exists(b.number(), b.strand, a.number, a_strand) {panic!("missing next link");}
                if !self.link_exists_prev(b.number(), b.strand, a.number, a_strand) {panic!("missing prev link");}
                if !self.link_exists(a.number, !a_strand, b.number(), !b.strand) {panic!("missing next link");}
                if !self.link_exists_prev(a.number, !a_strand, b.number(), !b.strand) {panic!("missing prev link");}
                if !self.unitig_index.contains_key(&b.number()) {panic!("unitig missing from index");}
            }
            for b in &a.reverse_prev {
                let a_strand = strand::REVERSE;
                if !self.link_exists(b.number(), b.strand, a.number, a_strand) {panic!("missing next link");}
                if !self.link_exists_prev(b.number(), b.strand, a.number, a_strand) {panic!("missing prev link");}
                if !self.link_exists(a.number, !a_strand, b.number(), !b.strand) {panic!("missing next link");}
                if !self.link_exists_prev(a.number, !a_strand, b.number(), !b.strand) {panic!("missing prev link");}
                if !self.unitig_index.contains_key(&b.number()) {panic!("unitig missing from index");}
            }
        }
    }
}


fn parse_unitig_path(path_str: &str) -> Vec<(u32, bool)> {
    path_str.split(',')
        .map(|u| {
            let strand = if u.ends_with('+') { strand::FORWARD } else if u.ends_with('-') { strand::REVERSE }
                         else { panic!("Invalid path strand") };
            let num = u[..u.len() - 1].parse::<u32>().expect("Error parsing unitig number");
            (num, strand)
        }).collect()
}


fn reverse_path(path: &[(u32, bool)]) -> Vec<(u32, bool)> {
    path.iter().rev().map(|&(num, strand)| (num, !strand)).collect()
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

    fn get_test_gfa_1() -> String {
        "H\tVN:Z:1.0\tKM:i:9\n\
        S\t1\tTTCGCTGCGCTCGCTTCGCTTT\tDP:f:1\n\
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
        "H\tVN:Z:1.0\tKM:i:9\n\
        S\t1\tACCGCTGCGCTCGCTTCGCTCT\tDP:f:1\n\
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

    fn get_test_gfa_3() -> String {
        "H\tVN:Z:1.0\tKM:i:9\n\
        S\t1\tTTCGCTGCGCTCGCTTCGCTTT\tDP:f:1\n\
        S\t2\tTGCCGTCGTCGCTGTGCA\tDP:f:1\n\
        S\t3\tTGCCTGAATCGCCTA\tDP:f:1\n\
        S\t4\tGCTCGGCTCG\tDP:f:1\n\
        S\t5\tCGAACCAT\tDP:f:1\n\
        S\t6\tTACTTGT\tDP:f:1\n\
        S\t7\tGCCTT\tDP:f:1\n\
        L\t1\t+\t2\t-\t0M\n\
        L\t2\t+\t1\t-\t0M\n\
        L\t2\t-\t3\t+\t0M\n\
        L\t3\t-\t2\t+\t0M\n\
        L\t3\t+\t4\t+\t0M\n\
        L\t4\t-\t3\t-\t0M\n\
        L\t4\t+\t5\t-\t0M\n\
        L\t5\t+\t4\t-\t0M\n\
        L\t5\t-\t5\t+\t0M\n\
        L\t3\t+\t6\t+\t0M\n\
        L\t6\t-\t3\t-\t0M\n\
        L\t6\t+\t7\t-\t0M\n\
        L\t7\t+\t6\t-\t0M\n\
        L\t7\t-\t6\t+\t0M\n\
        L\t6\t-\t7\t+\t0M\n".to_string()
    }

    #[test]
    fn test_graph_stats() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");

        make_test_file(&gfa_filename, &get_test_gfa_1());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);
        graph.check_links();
        assert_eq!(graph.k_size, 9);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.get_total_length(), 92);
        assert_eq!(graph.get_link_count(), 21);

        make_test_file(&gfa_filename, &get_test_gfa_2());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);
        graph.check_links();
        assert_eq!(graph.k_size, 9);
        assert_eq!(graph.unitigs.len(), 3);
        assert_eq!(graph.get_total_length(), 31);
        assert_eq!(graph.get_link_count(), 8);

        make_test_file(&gfa_filename, &get_test_gfa_3());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);
        graph.check_links();
        assert_eq!(graph.k_size, 9);
        assert_eq!(graph.unitigs.len(), 7);
        assert_eq!(graph.get_total_length(), 85);
        assert_eq!(graph.get_link_count(), 15);
    }

    #[test]
    fn test_parse_unitig_path() {
        assert_eq!(parse_unitig_path("2+,1-"), vec![(2, strand::FORWARD), (1, strand::REVERSE)]);
        assert_eq!(parse_unitig_path("3+,8-,4-"), vec![(3, strand::FORWARD), (8, strand::REVERSE), (4, strand::REVERSE)]);
    }

    #[test]
    fn test_reverse_path() {
        assert_eq!(reverse_path(&vec![(1, strand::FORWARD), (2, strand::REVERSE)]),
                                 vec![(2, strand::FORWARD), (1, strand::REVERSE)]);
        assert_eq!(reverse_path(&vec![(4, strand::FORWARD), (8, strand::FORWARD), (3, strand::REVERSE)]),
                                 vec![(3, strand::FORWARD), (8, strand::REVERSE), (4, strand::REVERSE)]);
    }

    #[test]
    fn test_link_exists_1() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_1());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);

        assert!(graph.link_exists(1, strand::FORWARD, 4, strand::FORWARD));
        assert!(graph.link_exists(4, strand::REVERSE, 1, strand::REVERSE));
        assert!(graph.link_exists(1, strand::FORWARD, 5, strand::REVERSE));
        assert!(graph.link_exists(5, strand::FORWARD, 1, strand::REVERSE));
        assert!(graph.link_exists(2, strand::FORWARD, 1, strand::FORWARD));
        assert!(graph.link_exists(1, strand::REVERSE, 2, strand::REVERSE));
        assert!(graph.link_exists(3, strand::REVERSE, 1, strand::FORWARD));
        assert!(graph.link_exists(1, strand::REVERSE, 3, strand::FORWARD));
        assert!(graph.link_exists(4, strand::FORWARD, 7, strand::REVERSE));
        assert!(graph.link_exists(7, strand::FORWARD, 4, strand::REVERSE));
        assert!(graph.link_exists(4, strand::FORWARD, 8, strand::FORWARD));
        assert!(graph.link_exists(8, strand::REVERSE, 4, strand::REVERSE));
        assert!(graph.link_exists(6, strand::REVERSE, 5, strand::REVERSE));
        assert!(graph.link_exists(5, strand::FORWARD, 6, strand::FORWARD));
        assert!(graph.link_exists(6, strand::FORWARD, 6, strand::REVERSE));
        assert!(graph.link_exists(7, strand::REVERSE, 9, strand::FORWARD));
        assert!(graph.link_exists(9, strand::REVERSE, 7, strand::FORWARD));
        assert!(graph.link_exists(8, strand::FORWARD, 10, strand::REVERSE));
        assert!(graph.link_exists(10, strand::FORWARD, 8, strand::REVERSE));
        assert!(graph.link_exists(9, strand::FORWARD, 7, strand::FORWARD));
        assert!(graph.link_exists(7, strand::REVERSE, 9, strand::REVERSE));

        assert!(!graph.link_exists(5, strand::REVERSE, 5, strand::FORWARD));
        assert!(!graph.link_exists(7, strand::FORWARD, 9, strand::FORWARD));
        assert!(!graph.link_exists(123, strand::FORWARD, 456, strand::FORWARD));
    }

    #[test]
    fn test_link_exists_2() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_2());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);

        assert!(graph.link_exists(1, strand::FORWARD, 2, strand::FORWARD));
        assert!(graph.link_exists(2, strand::REVERSE, 1, strand::REVERSE));
        assert!(graph.link_exists(1, strand::FORWARD, 2, strand::REVERSE));
        assert!(graph.link_exists(2, strand::FORWARD, 1, strand::REVERSE));
        assert!(graph.link_exists(1, strand::REVERSE, 3, strand::FORWARD));
        assert!(graph.link_exists(3, strand::REVERSE, 1, strand::FORWARD));
        assert!(graph.link_exists(1, strand::REVERSE, 3, strand::REVERSE));
        assert!(graph.link_exists(3, strand::FORWARD, 1, strand::FORWARD));

        assert!(!graph.link_exists(2, strand::FORWARD, 1, strand::FORWARD));
        assert!(!graph.link_exists(2, strand::FORWARD, 2, strand::REVERSE));
        assert!(!graph.link_exists(2, strand::REVERSE, 3, strand::REVERSE));
        assert!(!graph.link_exists(4, strand::FORWARD, 5, strand::FORWARD));
        assert!(!graph.link_exists(6, strand::REVERSE, 7, strand::REVERSE));
    }

    #[test]
    fn test_link_exists_3() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_3());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);

        assert!(graph.link_exists(1, strand::FORWARD, 2, strand::REVERSE));
        assert!(graph.link_exists(2, strand::FORWARD, 1, strand::REVERSE));
        assert!(graph.link_exists(2, strand::REVERSE, 3, strand::FORWARD));
        assert!(graph.link_exists(3, strand::REVERSE, 2, strand::FORWARD));
        assert!(graph.link_exists(3, strand::FORWARD, 4, strand::FORWARD));
        assert!(graph.link_exists(4, strand::REVERSE, 3, strand::REVERSE));
        assert!(graph.link_exists(4, strand::FORWARD, 5, strand::REVERSE));
        assert!(graph.link_exists(5, strand::FORWARD, 4, strand::REVERSE));
        assert!(graph.link_exists(5, strand::REVERSE, 5, strand::FORWARD));
        assert!(graph.link_exists(3, strand::FORWARD, 6, strand::FORWARD));
        assert!(graph.link_exists(6, strand::REVERSE, 3, strand::REVERSE));
        assert!(graph.link_exists(6, strand::FORWARD, 7, strand::REVERSE));
        assert!(graph.link_exists(7, strand::FORWARD, 6, strand::REVERSE));
        assert!(graph.link_exists(7, strand::REVERSE, 6, strand::FORWARD));
        assert!(graph.link_exists(6, strand::REVERSE, 7, strand::FORWARD));

        assert!(!graph.link_exists(1, strand::FORWARD, 3, strand::FORWARD));
        assert!(!graph.link_exists(5, strand::FORWARD, 5, strand::REVERSE));
        assert!(!graph.link_exists(7, strand::REVERSE, 4, strand::REVERSE));
        assert!(!graph.link_exists(8, strand::FORWARD, 9, strand::FORWARD));
    }

    #[test]
    fn test_link_exists_prev_1() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_1());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);

        assert!(graph.link_exists_prev(1, strand::FORWARD, 4, strand::FORWARD));
        assert!(graph.link_exists_prev(4, strand::REVERSE, 1, strand::REVERSE));
        assert!(graph.link_exists_prev(1, strand::FORWARD, 5, strand::REVERSE));
        assert!(graph.link_exists_prev(5, strand::FORWARD, 1, strand::REVERSE));
        assert!(graph.link_exists_prev(2, strand::FORWARD, 1, strand::FORWARD));
        assert!(graph.link_exists_prev(1, strand::REVERSE, 2, strand::REVERSE));
        assert!(graph.link_exists_prev(3, strand::REVERSE, 1, strand::FORWARD));
        assert!(graph.link_exists_prev(1, strand::REVERSE, 3, strand::FORWARD));
        assert!(graph.link_exists_prev(4, strand::FORWARD, 7, strand::REVERSE));
        assert!(graph.link_exists_prev(7, strand::FORWARD, 4, strand::REVERSE));
        assert!(graph.link_exists_prev(4, strand::FORWARD, 8, strand::FORWARD));
        assert!(graph.link_exists_prev(8, strand::REVERSE, 4, strand::REVERSE));
        assert!(graph.link_exists_prev(6, strand::REVERSE, 5, strand::REVERSE));
        assert!(graph.link_exists_prev(5, strand::FORWARD, 6, strand::FORWARD));
        assert!(graph.link_exists_prev(6, strand::FORWARD, 6, strand::REVERSE));
        assert!(graph.link_exists_prev(7, strand::REVERSE, 9, strand::FORWARD));
        assert!(graph.link_exists_prev(9, strand::REVERSE, 7, strand::FORWARD));
        assert!(graph.link_exists_prev(8, strand::FORWARD, 10, strand::REVERSE));
        assert!(graph.link_exists_prev(10, strand::FORWARD, 8, strand::REVERSE));
        assert!(graph.link_exists_prev(9, strand::FORWARD, 7, strand::FORWARD));
        assert!(graph.link_exists_prev(7, strand::REVERSE, 9, strand::REVERSE));

        assert!(!graph.link_exists_prev(5, strand::REVERSE, 5, strand::FORWARD));
        assert!(!graph.link_exists_prev(7, strand::FORWARD, 9, strand::FORWARD));
        assert!(!graph.link_exists_prev(123, strand::FORWARD, 456, strand::FORWARD));
    }

    #[test]
    fn test_link_exists_prev_2() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_2());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);

        assert!(graph.link_exists_prev(1, strand::FORWARD, 2, strand::FORWARD));
        assert!(graph.link_exists_prev(2, strand::REVERSE, 1, strand::REVERSE));
        assert!(graph.link_exists_prev(1, strand::FORWARD, 2, strand::REVERSE));
        assert!(graph.link_exists_prev(2, strand::FORWARD, 1, strand::REVERSE));
        assert!(graph.link_exists_prev(1, strand::REVERSE, 3, strand::FORWARD));
        assert!(graph.link_exists_prev(3, strand::REVERSE, 1, strand::FORWARD));
        assert!(graph.link_exists_prev(1, strand::REVERSE, 3, strand::REVERSE));
        assert!(graph.link_exists_prev(3, strand::FORWARD, 1, strand::FORWARD));

        assert!(!graph.link_exists_prev(2, strand::FORWARD, 1, strand::FORWARD));
        assert!(!graph.link_exists_prev(2, strand::FORWARD, 2, strand::REVERSE));
        assert!(!graph.link_exists_prev(2, strand::REVERSE, 3, strand::REVERSE));
        assert!(!graph.link_exists_prev(4, strand::FORWARD, 5, strand::FORWARD));
        assert!(!graph.link_exists_prev(6, strand::REVERSE, 7, strand::REVERSE));
    }

    #[test]
    fn test_link_exists_prev_3() {
        let temp_dir = tempdir().unwrap();
        let gfa_filename = temp_dir.path().join("graph.gfa");
        make_test_file(&gfa_filename, &get_test_gfa_3());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa_filename);

        assert!(graph.link_exists_prev(1, strand::FORWARD, 2, strand::REVERSE));
        assert!(graph.link_exists_prev(2, strand::FORWARD, 1, strand::REVERSE));
        assert!(graph.link_exists_prev(2, strand::REVERSE, 3, strand::FORWARD));
        assert!(graph.link_exists_prev(3, strand::REVERSE, 2, strand::FORWARD));
        assert!(graph.link_exists_prev(3, strand::FORWARD, 4, strand::FORWARD));
        assert!(graph.link_exists_prev(4, strand::REVERSE, 3, strand::REVERSE));
        assert!(graph.link_exists_prev(4, strand::FORWARD, 5, strand::REVERSE));
        assert!(graph.link_exists_prev(5, strand::FORWARD, 4, strand::REVERSE));
        assert!(graph.link_exists_prev(5, strand::REVERSE, 5, strand::FORWARD));
        assert!(graph.link_exists_prev(3, strand::FORWARD, 6, strand::FORWARD));
        assert!(graph.link_exists_prev(6, strand::REVERSE, 3, strand::REVERSE));
        assert!(graph.link_exists_prev(6, strand::FORWARD, 7, strand::REVERSE));
        assert!(graph.link_exists_prev(7, strand::FORWARD, 6, strand::REVERSE));
        assert!(graph.link_exists_prev(7, strand::REVERSE, 6, strand::FORWARD));
        assert!(graph.link_exists_prev(6, strand::REVERSE, 7, strand::FORWARD));

        assert!(!graph.link_exists_prev(1, strand::FORWARD, 3, strand::FORWARD));
        assert!(!graph.link_exists_prev(5, strand::FORWARD, 5, strand::REVERSE));
        assert!(!graph.link_exists_prev(7, strand::REVERSE, 4, strand::REVERSE));
        assert!(!graph.link_exists_prev(8, strand::FORWARD, 9, strand::FORWARD));
    }
}
