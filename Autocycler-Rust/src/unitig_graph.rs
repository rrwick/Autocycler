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
use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write, BufRead, BufReader};
use std::path::PathBuf;
use std::rc::Rc;

use crate::kmer_graph::KmerGraph;
use crate::position::Position;
use crate::sequence::Sequence;
use crate::unitig::Unitig;
use crate::misc::quit_with_error;


pub struct UnitigGraph {
    pub unitigs: Vec<Rc<RefCell<Unitig>>>,
    pub k_size: u32,
    pub link_count: usize,
}

impl UnitigGraph {
    pub fn from_kmer_graph(k_graph: &KmerGraph) -> Self {
        let mut u_graph = UnitigGraph {
            unitigs: Vec::new(),
            k_size: k_graph.k_size,
            link_count: 0,
        };
        u_graph.build_unitigs_from_kmer_graph(k_graph);
        u_graph.simplify_seqs();
        u_graph.create_links();
        u_graph.trim_overlaps();
        u_graph.renumber_unitigs();
        u_graph
    }

    pub fn from_gfa_file(gfa_filename: &PathBuf) -> Self {
        let mut u_graph = UnitigGraph {
            unitigs: Vec::new(),
            k_size: 0,
            link_count: 0,
        };
        let file = File::open(gfa_filename).unwrap();
        let reader = BufReader::new(file);
        let mut link_lines: Vec<String> = Vec::new();
        let mut path_lines: Vec<String> = Vec::new();
        for line_result in reader.lines() {
            let line = line_result.unwrap();
            let parts: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
            match parts.get(0) {
                Some(&"H") => u_graph.read_gfa_header_line(&parts),
                Some(&"S") => u_graph.unitigs.push(Rc::new(RefCell::new(
                                                   Unitig::from_segment_line(&line)))),
                Some(&"L") => link_lines.push(line),
                Some(&"P") => path_lines.push(line),
                _ => {}
            }
        }
        let unitig_index: HashMap<u32, Rc<RefCell<Unitig>>> =
            u_graph.unitigs.iter().map(|u| {(u.borrow().number, Rc::clone(u))}).collect();
        u_graph.build_links_from_gfa(&link_lines, &unitig_index);
        u_graph.build_paths_from_gfa(&path_lines, &unitig_index);
        u_graph
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

    fn build_links_from_gfa(&mut self, link_lines: &Vec<String>,
                            unitig_index: &HashMap<u32, Rc<RefCell<Unitig>>>) {
        self.link_count = 0;
        for line in link_lines {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 6 || parts[5] != "0M" {
                quit_with_error("non-zero overlap found on the GFA link line.\n\
                                 Are you sure this is an Autocycler-generated GFA file?");
                return;
            }
            let seg_1: u32 = parts[1].parse().expect("Error parsing segment 1 as integer");
            let seg_2: u32 = parts[3].parse().expect("Error parsing segment 2 as integer");
            let strand_1 = parts[2] == "+";
            let strand_2 = parts[4] == "+";
            if let Some(unitig_1) = unitig_index.get(&seg_1) {
                if let Some(unitig_2) = unitig_index.get(&seg_2) {
                    if strand_1 {
                        unitig_1.borrow_mut().forward_next.push((Rc::clone(unitig_2), strand_2));
                    } else {
                        unitig_1.borrow_mut().reverse_next.push((Rc::clone(unitig_2), strand_2));
                    }
                    if strand_2 {
                        unitig_2.borrow_mut().forward_prev.push((Rc::clone(unitig_1), strand_1));
                    } else {
                        unitig_2.borrow_mut().reverse_prev.push((Rc::clone(unitig_1), strand_1));
                    }
                    self.link_count += 1;
                } else {
                    quit_with_error(&format!("link refers to nonexistent unitig: {}", seg_2));
                }
            } else {
                quit_with_error(&format!("link refers to nonexistent unitig: {}", seg_1));
            }
        }
    }

    fn build_paths_from_gfa(&mut self, path_lines: &Vec<String>,
                            unitig_index: &HashMap<u32, Rc<RefCell<Unitig>>>) {
        for line in path_lines {
            let parts: Vec<&str> = line.split('\t').collect();
            let seq_id: u16 = parts[1].parse().expect("Error parsing sequence ID as integer");
            let mut length = None;
            let mut filename = None;
            let mut header = None;
            for p in &parts[2..] {
                if p.starts_with("LN:i:") {
                    length = Some(p[5..].parse::<u32>().expect("Error parsing length"));
                } else if p.starts_with("FN:Z:") {
                    filename = Some(p[5..].to_string());
                } else if p.starts_with("HD:Z:") {
                    header = Some(p[5..].to_string());
                }
            }
            if length.is_none() || filename.is_none() || header.is_none() {
                quit_with_error("missing required tag in GFA path line.");
            }
            let length = length.unwrap();
            let forward_path = Self::parse_unitig_path(parts[2]);
            let reverse_path = Self::reverse_path(&forward_path);
    
            self.add_positions_from_path(&forward_path, true, seq_id, unitig_index, length);
            self.add_positions_from_path(&reverse_path, false, seq_id, unitig_index, length);
        }
    }

    fn add_positions_from_path(&mut self, path: &[(u32, bool)], path_strand: bool, seq_id: u16,
                               unitig_index: &HashMap<u32, Rc<RefCell<Unitig>>>, length: u32) {
        let half_k = self.k_size / 2;
        let mut pos = half_k;

        for (unitig_num, unitig_strand) in path {
            if let Some(unitig) = unitig_index.get(unitig_num) {
                let mut u = unitig.borrow_mut();
                let positions = if *unitig_strand {
                    &mut u.forward_positions
                } else {
                    &mut u.reverse_positions
                };
                positions.push(Position::new(seq_id, path_strand, pos as usize));
                pos += u.length();
                if u.dead_end_start(*unitig_strand) {
                    pos -= half_k;
                }
                if u.dead_end_end(*unitig_strand) {
                    pos -= half_k;
                }
            } else {
                quit_with_error(&format!("unitig {} not found in unitig index", unitig_num));
            }
        }
    
        assert!(pos + half_k == length, "Position calculation mismatch");

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
            let starting_kmer = forward_kmer;

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
            let mut rev_k = reverse_kmer;
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
        self.link_count = 0;
        for i in 0..self.unitigs.len() {
            let unitig_a = Rc::clone(&self.unitigs[i]);
            let ending_forward_seq = unitig_a.borrow().forward_seq[unitig_a.borrow().forward_seq.len() - piece_len..].to_vec();
            let ending_reverse_seq = unitig_a.borrow().reverse_seq[unitig_a.borrow().reverse_seq.len() - piece_len..].to_vec();

            if let Some(next_idxs) = forward_starts.get(&ending_forward_seq) {
                for &j in next_idxs {
                    let unitig_b = Rc::clone(&self.unitigs[j]);

                    // unitig_a+ -> unitig_b+
                    unitig_a.borrow_mut().forward_next.push((Rc::clone(&unitig_b), true));
                    unitig_b.borrow_mut().forward_prev.push((Rc::clone(&unitig_a), true));
                    self.link_count += 1;

                    // unitig_b- -> unitig_a-
                    unitig_b.borrow_mut().reverse_next.push((Rc::clone(&unitig_a), false));
                    unitig_a.borrow_mut().reverse_prev.push((Rc::clone(&unitig_b), false));
                    self.link_count += 1;
                }
            }

            if let Some(next_idxs) = reverse_starts.get(&ending_forward_seq) {
                for &j in next_idxs {
                    let unitig_b = Rc::clone(&self.unitigs[j]);

                    // unitig_a+ -> unitig_b-
                    unitig_a.borrow_mut().forward_next.push((Rc::clone(&unitig_b), false));
                    unitig_b.borrow_mut().reverse_prev.push((Rc::clone(&unitig_a), true));
                    self.link_count += 1;
                }
            }

            if let Some(next_idxs) = forward_starts.get(&ending_reverse_seq) {
                for &j in next_idxs {
                    let unitig_b = Rc::clone(&self.unitigs[j]);

                    // unitig_a- -> unitig_b+
                    unitig_a.borrow_mut().reverse_next.push((Rc::clone(&unitig_b), true));
                    unitig_b.borrow_mut().forward_prev.push((Rc::clone(&unitig_a), false));
                    self.link_count += 1;
                }
            }
        }
    }

    fn trim_overlaps(&mut self) {
        for unitig in &self.unitigs {
            unitig.borrow_mut().trim_overlaps(self.k_size as usize);
        }
    }

    fn renumber_unitigs(&mut self) {
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

    fn get_links_for_gfa(&self) -> Vec<(String, String, String, String)> {
        let mut links = Vec::new();
        for a_rc in &self.unitigs {
            let a = a_rc.borrow();
            for (b_rc, b_strand) in &a.forward_next {
                let b = b_rc.borrow();
                links.push((a.number.to_string(), "+".to_string(), b.number.to_string(),
                            (if *b_strand {"+"} else {"-"}).to_string()));
            }

            for (b_rc, b_strand) in &a.reverse_next {
                let b = b_rc.borrow();
                links.push((a.number.to_string(), "-".to_string(), b.number.to_string(),
                            (if *b_strand {"+"} else {"-"}).to_string()));
            }
        }
        links
    }

    fn get_gfa_path_line(&self, seq: &Sequence) -> String {
        let unitig_path = self.get_unitig_path_for_sequence(seq);
        let path_str: Vec<String> = unitig_path.iter()
            .map(|(num, strand)| format!("{}{}", num, if *strand { "+" } else { "-" })).collect();
        let path_str = path_str.join(",");
        format!("P\t{}\t{}\t*\tLN:i:{}\tFN:Z:{}\tHD:Z:{}",
                seq.id, path_str, seq.length, seq.filename, seq.contig_header)
    }

    pub fn reconstruct_original_sequences(&self) -> Vec<Sequence> {
        let mut original_seqs = Vec::new();
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
        original_seqs
    }

    fn reconstruct_original_sequence(&self) -> Sequence {
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
        Sequence::new(1, "ACGACTGACATCAGCACTGA".to_string(),
                      "assembly.fasta".to_string(), "contig_1".to_string(), 20)  // TEMP
    }

    fn find_starting_unitig(&self, seq_id: u16) -> (Rc<RefCell<Unitig>>, bool) {
        // For a given sequence ID, this function returns the Unitig and strand where that sequence
        // begins.
        let half_k = self.k_size / 2;
        let mut starting_unitigs = Vec::new();
        for unitig in &self.unitigs {
            for p in &unitig.borrow().forward_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == half_k {
                    starting_unitigs.push((Rc::clone(unitig), true));
                }
            }
            for p in &unitig.borrow().reverse_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == half_k {
                    starting_unitigs.push((Rc::clone(unitig), false));
                }
            }
        }
        assert_eq!(starting_unitigs.len(), 1);
        starting_unitigs[0].clone()
    }

    fn get_next_unitig(&self, seq_id: u16, unitig_rc: &Rc<RefCell<Unitig>>, strand: bool,
                       pos: u32) -> Option<(Rc<RefCell<Unitig>>, bool, u32)> {
        let unitig = unitig_rc.borrow();
        let next_pos = pos + unitig.untrimmed_length(self.k_size as usize) - self.k_size as u32 + 1;
        let next_unitigs = if strand { &unitig.forward_next } else { &unitig.reverse_next };
        for (next_unitig_rc, next_strand) in next_unitigs {
            let u = next_unitig_rc.borrow();
            let positions = if *next_strand { &u.forward_positions } else { &u.reverse_positions};
            for p in positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == next_pos {
                    return Some((Rc::clone(next_unitig_rc), *next_strand, next_pos));
                }
            }
        }
        None
    }

    fn get_unitig_path_for_sequence(&self, seq: &Sequence) -> Vec<(u32, bool)> {
        let half_k = self.k_size / 2;
        let mut unitig_path = Vec::new();
        let (mut unitig, mut strand) = self.find_starting_unitig(seq.id);
        let mut pos = half_k;
        loop {
            unitig_path.push((unitig.borrow().number, strand));
            match self.get_next_unitig(seq.id, &unitig, strand, pos) {
                None => break,
                Some((next_unitig, next_strand, next_pos)) => {
                    (unitig, strand, pos) = (next_unitig, next_strand, next_pos);
                }
            }
        }
        unitig_path
    }

    fn parse_unitig_path(path_str: &str) -> Vec<(u32, bool)> {
        path_str.split(',')
            .map(|u| {
                let strand = if u.ends_with('+') { true } else if u.ends_with('-') { false }
                             else { panic!("Invalid path strand") };
                let num = u[..u.len() - 1].parse::<u32>().expect("Error parsing unitig number");
                (num, strand)
            }).collect()
    }

    fn reverse_path(path: &[(u32, bool)]) -> Vec<(u32, bool)> {
        path.iter().rev().map(|&(num, strand)| (num, !strand)).collect()
    }
}
