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
use std::io::{self, Write};
use std::path::Path;
use std::rc::Rc;

use crate::kmer_graph::KmerGraph;
use crate::position::Position;
use crate::sequence::Sequence;
use crate::unitig::{Unitig, UnitigStrand};
use crate::misc::{quit_with_error, strand, load_file_lines, find_replace_i32_tuple};


#[derive(Default)]
pub struct UnitigGraph {
    pub unitigs: Vec<Rc<RefCell<Unitig>>>,
    pub k_size: u32,
    pub unitig_index: HashMap<u32, Rc<RefCell<Unitig>>>,
}

impl UnitigGraph {
    pub fn from_kmer_graph(k_graph: &KmerGraph) -> Self {
        let mut u_graph = UnitigGraph {
            k_size: k_graph.k_size,
            ..Default::default()
        };
        u_graph.build_unitigs_from_kmer_graph(k_graph);
        u_graph.simplify_seqs();
        u_graph.create_links();
        u_graph.trim_overlaps();
        u_graph.renumber_unitigs();
        u_graph.check_links();
        u_graph
    }

    pub fn from_gfa_file(gfa_filename: &Path) -> (Self, Vec<Sequence>) {
        let gfa_lines = load_file_lines(gfa_filename);
        Self::from_gfa_lines(&gfa_lines)
    }

    pub fn from_gfa_lines(gfa_lines: &Vec<String>) -> (Self, Vec<Sequence>) {
        let mut u_graph = UnitigGraph::default();
        let mut link_lines: Vec<&str> = Vec::new();
        let mut path_lines: Vec<&str> = Vec::new();
        for line in gfa_lines {
            let parts: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
            match parts.first() {
                Some(&"H") => u_graph.read_gfa_header_line(&parts),
                Some(&"S") => u_graph.unitigs.push(Rc::new(RefCell::new(Unitig::from_segment_line(line)))),
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
            if let Some(tag_val) = p.strip_prefix("KM:i:") {
                if let Ok(k) = tag_val.parse::<u32>() {
                    self.k_size = k;
                    return;
                }
            }
        }
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
            for p in &parts[2..] {
                if let Some(tag_val) = p.strip_prefix("LN:i:") {
                    length = Some(tag_val.parse::<u32>().expect("Error parsing length"));
                } else if let Some(tag_val) = p.strip_prefix("FN:Z:") {
                    filename = Some(tag_val.to_string());
                } else if let Some(tag_val) = p.strip_prefix("HD:Z:") {
                    header = Some(tag_val.to_string());
                } else if let Some(tag_val) = p.strip_prefix("CL:i:") {
                    cluster = tag_val.parse::<u16>().expect("Error parsing cluster");
                }
            }
            if length.is_none() || filename.is_none() || header.is_none() {
                quit_with_error("missing required tag in GFA path line.");
            }
            let length = length.unwrap();
            let filename = filename.unwrap();
            let header = header.unwrap();
            let path = parse_unitig_path(parts[2]);
            let sequence = self.create_sequence_and_positions(seq_id, length, filename, header,
                                                              cluster, path);
            sequences.push(sequence);
        }
        sequences
    }

    pub fn create_sequence_and_positions(&mut self, seq_id: u16, length: u32,
                                         filename: String, header: String, cluster: u16,
                                         forward_path: Vec<(u32, bool)>) -> Sequence {
        let reverse_path = reverse_path(&forward_path);
        self.add_positions_from_path(&forward_path, strand::FORWARD, seq_id, length);
        self.add_positions_from_path(&reverse_path, strand::REVERSE, seq_id, length);
        Sequence::new_without_seq(seq_id, filename, header, length as usize, cluster)
    }

    fn add_positions_from_path(&mut self, path: &[(u32, bool)], path_strand: bool, seq_id: u16, length: u32) {
        let mut pos = 0;
        for (unitig_num, unitig_strand) in path {
            if let Some(unitig) = self.unitig_index.get(unitig_num) {
                let mut u = unitig.borrow_mut();
                let positions = if *unitig_strand {&mut u.forward_positions} 
                                             else {&mut u.reverse_positions};
                positions.push(Position::new(seq_id, path_strand, pos as usize));
                pos += u.length();
            } else {
                quit_with_error(&format!("unitig {} not found in unitig index", unitig_num));
            }
        }
        assert!(pos == length, "Position calculation mismatch");
    }

    fn build_unitigs_from_kmer_graph(&mut self, k_graph: &KmerGraph) {
        let mut seen: HashSet<&[u8]> = HashSet::new();
        let mut unitig_number = 0;
        for forward_kmer in k_graph.iterate_kmers() {
            if seen.contains(forward_kmer.seq()) {
                continue;
            }
            let reverse_kmer = k_graph.reverse(forward_kmer);
            unitig_number += 1;
            let mut unitig = Unitig::from_kmers(unitig_number, forward_kmer, reverse_kmer);
            seen.insert(forward_kmer.seq());
            seen.insert(reverse_kmer.seq());

            // Extend unitig forward
            let mut for_k = forward_kmer;
            let mut rev_k = reverse_kmer;
            loop {
                if rev_k.first_position() { break; }
                let next_kmers = k_graph.next_kmers(for_k.seq());
                if next_kmers.len() != 1 { break; }
                for_k = next_kmers[0];
                if seen.contains(for_k.seq()) { break; }
                let prev_kmers = k_graph.prev_kmers(for_k.seq());
                if prev_kmers.len() != 1 { break; }
                rev_k = k_graph.reverse(for_k);
                if for_k.first_position() { break; }
                unitig.add_kmer_to_end(for_k, rev_k);
                seen.insert(for_k.seq());
                seen.insert(rev_k.seq());
            }

            // Extend unitig backward
            let mut for_k = forward_kmer;
            let mut rev_k;
            loop {
                if for_k.first_position() { break; }
                let prev_kmers = k_graph.prev_kmers(for_k.seq());
                if prev_kmers.len() != 1 { break; }
                for_k = prev_kmers[0];
                if seen.contains(for_k.seq()) { break; }
                let next_kmers = k_graph.next_kmers(for_k.seq());
                if next_kmers.len() != 1 { break; }
                rev_k = k_graph.reverse(for_k);
                if rev_k.first_position() { break; }
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
        for unitig in &self.unitigs {
            unitig.borrow_mut().trim_overlaps(self.k_size as usize);
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

    pub fn save_gfa(&self, gfa_filename: &Path, sequences: &Vec<Sequence>,
                    use_other_colour: bool) -> io::Result<()> {
        let mut file = File::create(gfa_filename)?;
        writeln!(file, "H\tVN:Z:1.0\tKM:i:{}", self.k_size)?;
        for unitig in &self.unitigs {
            writeln!(file, "{}", unitig.borrow().gfa_segment_line(use_other_colour))?;
        }
        for (a, a_strand, b, b_strand) in self.get_links_for_gfa(0) {
            writeln!(file, "L\t{}\t{}\t{}\t{}\t0M", a, a_strand, b, b_strand)?;
        }
        for s in sequences {
            writeln!(file, "{}", self.get_gfa_path_line(s))?;
        }
        Ok(())
    }

    pub fn get_links_for_gfa(&self, offset: u32) -> Vec<(String, String, String, String)> {
        let mut links = Vec::new();
        for a_rc in &self.unitigs {
            let a = a_rc.borrow();
            let a_num = a.number + offset;
            for b in &a.forward_next {
                let b_num = b.number() + offset;
                links.push((a_num.to_string(), "+".to_string(), b_num.to_string(),
                            (if b.strand {"+"} else {"-"}).to_string()));
            }
            for b in &a.reverse_next {
                let b_num = b.number() + offset;
                links.push((a_num.to_string(), "-".to_string(), b_num.to_string(),
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
        format!("P\t{}\t{}\t*\tLN:i:{}\tFN:Z:{}\tHD:Z:{}{}",
                seq.id, path_str, seq.length, seq.filename, seq.contig_header, cluster_tag)
    }

    pub fn reconstruct_original_sequences(&self, seqs: &Vec<Sequence>)
            -> HashMap<String, Vec<(String, String)>> {
        let mut original_seqs = HashMap::new();
        for seq in seqs {
            let (filename, header, sequence) = self.reconstruct_original_sequence(seq);
            original_seqs.entry(filename).or_insert_with(Vec::new).push((header, sequence));
        }
        original_seqs
    }

    pub fn reconstruct_original_sequences_u8(&self, seqs: &Vec<Sequence>)
            -> Vec<((String, String), Vec<u8>)> {
        let mut original_seqs = Vec::new();
        for seq in seqs {
            let (filename, _header, sequence) = self.reconstruct_original_sequence(seq);
            original_seqs.push(((filename, seq.contig_name()), sequence.as_bytes().to_owned()));
        }
        original_seqs.sort();
        original_seqs
    }

    fn reconstruct_original_sequence(&self, seq: &Sequence) -> (String, String, String) {
        let path = self.get_unitig_path_for_sequence(seq);
        let sequence = self.get_sequence_from_path(&path);
        assert_eq!(sequence.len(), seq.length, "reconstructed sequence does not have expected length");
        (seq.filename.clone(), seq.contig_header.clone(), sequence)
    }

    fn get_sequence_from_path(&self, path: &[(u32, bool)]) -> String {
        // Given a path (vector of unitig IDs and strands), this function returns the sequence
        // traced by that path. It also requires a unitig index so it can quickly look up unitigs
        // by their number.
        let mut sequence = Vec::new();
        for (unitig_num, strand) in path.iter() {
            let unitig = self.unitig_index.get(unitig_num).unwrap();
            sequence.push(String::from_utf8(unitig.borrow().get_seq(*strand)).unwrap());
        }
        sequence.into_iter().collect()
    }

    pub fn get_sequence_from_path_signed(&self, path: &[i32]) -> Vec<u8> {
        let path: Vec<_> = path.iter().map(|&x| (x.unsigned_abs(), x >= 0)).collect();
        self.get_sequence_from_path(&path).as_bytes().to_owned()
    }

    fn find_starting_unitig(&self, seq_id: u16) -> UnitigStrand {
        // For a given sequence ID, this function returns the Unitig and strand where that sequence
        // begins.
        let mut starting_unitigs = Vec::new();
        for unitig in &self.unitigs {
            for p in &unitig.borrow().forward_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == 0 {
                    starting_unitigs.push(UnitigStrand::new(unitig, strand::FORWARD));
                }
            }
            for p in &unitig.borrow().reverse_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == 0 {
                    starting_unitigs.push(UnitigStrand::new(unitig, strand::REVERSE));
                }
            }
        }
        assert_eq!(starting_unitigs.len(), 1);
        starting_unitigs[0].clone()
    }

    pub fn get_next_unitig(&self, seq_id: u16, seq_strand: bool, unitig_rc: &Rc<RefCell<Unitig>>,
                           strand: bool, pos: u32) -> Option<(UnitigStrand, u32)> {
        // For a given unitig that's part of a sequence's path, this function will return the next
        // unitig in that sequence's path.
        let unitig = unitig_rc.borrow();
        let next_pos = pos + unitig.length();
        let next_unitigs = if strand { &unitig.forward_next } else { &unitig.reverse_next };
        for next in next_unitigs {
            let u = next.unitig.borrow();
            let positions = if next.strand { &u.forward_positions } else { &u.reverse_positions};
            for p in positions {
                if p.seq_id() == seq_id && p.strand() == seq_strand && p.pos == next_pos {
                    return Some((UnitigStrand::new(&next.unitig, next.strand), next_pos));
                }
            }
        }
        None
    }

    pub fn get_unitig_path_for_sequence(&self, seq: &Sequence) -> Vec<(u32, bool)> {
        let mut unitig_path = Vec::new();
        let mut u = self.find_starting_unitig(seq.id);
        let mut pos = 0;
        loop {
            unitig_path.push((u.number(), u.strand));
            match self.get_next_unitig(seq.id, strand::FORWARD, &u.unitig, u.strand, pos) {
                None => break,
                Some((next, next_pos)) => {
                    (u, pos) = (next, next_pos);
                }
            }
        }
        unitig_path
    }

    pub fn get_unitig_path_for_sequence_i32(&self, seq: &Sequence) -> Vec<i32> {
        // Same as the above function, but instead of giving unitig IDs and strands as a (u32, bool)
        // tuple, it gives them as i32 (negative numbers for reverse strand).
        let unitig_path = self.get_unitig_path_for_sequence(seq);
        unitig_path.iter().map(|(u, s)| if *s { *u as i32 } else { -(*u as i32)}).collect()
    }

    pub fn total_length(&self) -> u64 {
        self.unitigs.iter().map(|u| u.borrow().length() as u64).sum()
    }

    pub fn link_count(&self) -> (usize, usize) {
        // Returns the number of links in the graph in two ways:
        // * All links (both directions, as Bandage would show in double mode)
        // * Single-direction links (no redundancy, as Bandage would show in single mode)
        // Usually the former is twice the latter, but not when there are hairpin links (which don't
        // have a reverse).
        let mut all_links = HashSet::new();
        let mut one_way_links = HashSet::new();
        for a_rc in &self.unitigs {
            let a = a_rc.borrow();
            let a_num = a_rc.borrow().number as i32;
            for b in &a.forward_next {
                let b_num = b.signed_number();
                let link = (a_num, b_num);
                let rev_link = (-b_num, -a_num);
                all_links.insert(link);
                all_links.insert(rev_link);
                one_way_links.insert(if link > rev_link { link } else { rev_link });
            }
            for b in &a.reverse_next {
                let b_num = b.signed_number();
                let link = (-a_num, b_num);
                let rev_link = (-b_num, a_num);
                all_links.insert(link);
                all_links.insert(rev_link);
                one_way_links.insert(if link > rev_link { link } else { rev_link });
            }
        }
        (all_links.len(), one_way_links.len())
    }

    pub fn print_basic_graph_info(&self) {
        let link_count = self.link_count().1;
        eprintln!("{} unitig{}, {} link{}",
                  self.unitigs.len(), match self.unitigs.len() { 1 => "", _ => "s" },
                  link_count, match link_count { 1 => "", _ => "s" });
        eprintln!("total length: {} bp", self.total_length());
        eprintln!();
    }

    pub fn topology(&self) -> String {
        // Returns one of the following:
        // * circular: the graph contains one unitig with a simple circularising link
        // * linear_blunt_blunt: the graph contains one unitig with no link (both ends are blunt)
        // * linear_hairpin_hairpin: the graph contains one unitig with hairpin links on both ends
        // * linear_blunt_hairpin: the graph contains one unitig with a hairpin link on one end
        // * fragmented: the graph contains multiple unitigs
        // * empty: the graph contains no unitigs
        // * other: none of the above (e.g. a circularising link and a hairpin link, should be rare)
        if self.unitigs.is_empty() { return "empty".to_string(); }
        if self.unitigs.len() > 1 { return "fragmented".to_string(); }
        let u = self.unitigs[0].borrow();  // the only unitig in the graph
        if self.link_count().0 == 0 { return "linear_blunt_blunt".to_string(); }
        if u.is_isolated_and_circular() { return "circular".to_string(); }
        if u.hairpin_start() && u.hairpin_end() { return "linear_hairpin_hairpin".to_string(); }
        if u.hairpin_start() && u.blunt_end() { return "linear_blunt_hairpin".to_string(); }
        if u.blunt_start() && u.hairpin_end() { return "linear_blunt_hairpin".to_string(); }
        "other".to_string()
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
        self.unitigs.retain(|u| u.borrow().depth > 0.0);
        self.delete_dangling_links();
        self.build_unitig_index();
    }

    pub fn remove_unitigs_by_number(&mut self, to_remove: HashSet<u32>) {
        self.unitigs.retain(|u| !to_remove.contains(&u.borrow().number));
        self.delete_dangling_links();
        self.build_unitig_index();
    }

    pub fn duplicate_unitig_by_number(&mut self, unitig_num: &u32) {
        // This method duplicates the specified unitig. It is a requirement that the unitig has
        // exactly two non-self links, because each copy of the unitig will keep one of these links.
        self.check_if_unitig_can_be_duplicated(unitig_num);

        // Create two copies of the target unitig.
        let u = self.unitig_index.get(unitig_num).unwrap().clone();
        let target = u.borrow();
        let target_num = target.number;
        let mut copy_a = target.clone();
        let mut copy_b = target.clone();
        copy_a.depth /= 2.0;
        copy_b.depth /= 2.0;
        let a_num = self.max_unitig_number() + 1;
        let b_num = a_num + 1;
        copy_a.number = a_num;
        copy_b.number = b_num;
        copy_a.clear_all_links();
        copy_b.clear_all_links();
        self.unitigs.push(Rc::new(RefCell::new(copy_a)));
        self.unitigs.push(Rc::new(RefCell::new(copy_b)));

        // Remove the original unitig.
        self.remove_unitigs_by_number(std::iter::once(*unitig_num).collect());
        self.delete_dangling_links();
        self.build_unitig_index();

        // Add self-links (loops and hairpins) to the copies
        for link in &target.forward_next {
            if link.number() == *unitig_num {
                self.create_link(a_num as i32, a_num as i32 * if link.strand { 1 } else { -1 });
                self.create_link(b_num as i32, b_num as i32 * if link.strand { 1 } else { -1 });
            }
        }
        for link in &target.reverse_next {
            if link.number() == *unitig_num {
                self.create_link(-(a_num as i32), a_num as i32 * if link.strand { 1 } else { -1 });
                self.create_link(-(b_num as i32), b_num as i32 * if link.strand { 1 } else { -1 });
            }
        }

        // Distribute the two non-self links to the copies.
        let mut non_self_links = Vec::new();
        for link in &target.forward_next {
            if link.number() != target_num {
                non_self_links.push((target_num as i32, link.signed_number()));
            }
        }
        for link in &target.reverse_next {
            if link.number() != target_num {
                non_self_links.push((-(target_num as i32), link.signed_number()));
            }
        }
        assert!(non_self_links.len() == 2);
        let new_link_a = find_replace_i32_tuple(non_self_links[0], target_num as i32, a_num as i32);
        let new_link_b = find_replace_i32_tuple(non_self_links[1], target_num as i32, b_num as i32);
        self.create_link(new_link_a.0, new_link_a.1);
        self.create_link(new_link_b.0, new_link_b.1);
        self.check_links();
    }

    fn check_if_unitig_can_be_duplicated(&self, unitig_num: &u32) {
        let unitig = self.unitig_index.get(unitig_num).unwrap().borrow();
        let unitig_num = unitig.number;
        let mut links = Vec::new();
        for link in &unitig.forward_next {
            if unitig_num != link.number() { links.push(link); }
        }
        for link in &unitig.reverse_next {
            if unitig_num != link.number() { links.push(link); }
        }
        if links.len() != 2 {
            quit_with_error(&format!("unitig {} does not contain exactly two non-self links",
                                     unitig_num));
        }
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

    pub fn delete_outgoing_links(&mut self, signed_num: i32) {
        let strand = if signed_num > 0 { strand::FORWARD } else { strand::REVERSE };
        let unitig_num = signed_num.unsigned_abs();
        let next_numbers: Vec<i32> = {
            let unitig = self.unitig_index.get(&unitig_num).unwrap().borrow();
            let next_unitigs = if strand { &unitig.forward_next } else { &unitig.reverse_next }; 
            next_unitigs.iter().map(|u| u.signed_number()).collect()
        };
        for next_num in next_numbers {
            self.delete_link(signed_num, next_num);
        }
    }

    pub fn delete_incoming_links(&mut self, signed_num: i32) {
        let strand = if signed_num > 0 { strand::FORWARD } else { strand::REVERSE };
        let unitig_num = signed_num.unsigned_abs();
        let prev_numbers: Vec<i32> = {
            let unitig = self.unitig_index.get(&unitig_num).unwrap().borrow();
            let prev_unitigs = if strand { &unitig.forward_prev } else { &unitig.reverse_prev }; 
            prev_unitigs.iter().map(|u| u.signed_number()).collect()
        };
        for prev_num in prev_numbers {
            self.delete_link(prev_num, signed_num);
        }
    }

    pub fn delete_link(&mut self, start_num: i32, end_num: i32) {
        self.delete_link_one_way(start_num, end_num);
        self.delete_link_one_way(-end_num, -start_num);
    }

    fn delete_link_one_way(&mut self, start_num: i32, end_num: i32) {
        let start_strand = if start_num > 0 { strand::FORWARD } else { strand::REVERSE };
        let end_strand = if end_num > 0 { strand::FORWARD } else { strand::REVERSE };
        let start_num = start_num.unsigned_abs();
        let end_num = end_num.unsigned_abs();
        let start_rc = self.unitig_index.get(&start_num).unwrap();
        let end_rc = self.unitig_index.get(&end_num).unwrap();

        // Collect the indices to remove for start unitig
        let start_indices: Vec<usize> = {
            let start = start_rc.borrow();
            let next_unitigs = if start_strand { &start.forward_next } else { &start.reverse_next };
            next_unitigs.iter().enumerate().filter_map(|(i, connection)| { if connection.unitig.borrow().number == end_num && connection.strand == end_strand { Some(i) } else { None } }).collect()
        };

        // Remove the elements from start unitig
        {
            let mut start = start_rc.borrow_mut();
            let next_unitigs = if start_strand { &mut start.forward_next } else { &mut start.reverse_next };
            for &i in start_indices.iter().rev() {
                next_unitigs.remove(i);
            }
        }

        // Collect the indices to remove for end unitig
        let end_indices: Vec<usize> = {
            let end = end_rc.borrow();
            let prev_unitigs = if end_strand { &end.forward_prev } else { &end.reverse_prev };
            prev_unitigs.iter().enumerate().filter_map(|(i, connection)| { if connection.unitig.borrow().number == start_num && connection.strand == start_strand { Some(i) } else { None } }).collect()
        };

        // Remove the elements from end unitig
        {
            let mut end = end_rc.borrow_mut();
            let prev_unitigs = if end_strand { &mut end.forward_prev } else { &mut end.reverse_prev };
            for &i in end_indices.iter().rev() {
                prev_unitigs.remove(i);
            }
        }
    }

    pub fn create_link(&mut self, start_num: i32, end_num: i32) {
        self.create_link_one_way(start_num, end_num);
        if start_num != -end_num {
            self.create_link_one_way(-end_num, -start_num);
        }
    }

    fn create_link_one_way(&mut self, start_num: i32, end_num: i32) {
        let start_strand = if start_num > 0 { strand::FORWARD } else { strand::REVERSE };
        let end_strand = if end_num > 0 { strand::FORWARD } else { strand::REVERSE };
        let start_num = start_num.unsigned_abs();
        let end_num = end_num.unsigned_abs();
        let start_rc = self.unitig_index.get(&start_num).unwrap();
        let end_rc = self.unitig_index.get(&end_num).unwrap();
        {
            let mut start = start_rc.borrow_mut();
            let connection = UnitigStrand { unitig: Rc::clone(end_rc), strand: end_strand };
            let next_unitigs = if start_strand { &mut start.forward_next } else { &mut start.reverse_next };
            next_unitigs.push(connection);
        }
        {
            let mut end = end_rc.borrow_mut();
            let reverse_connection = UnitigStrand { unitig: Rc::clone(start_rc), strand: start_strand };
            let prev_unitigs = if end_strand { &mut end.forward_prev } else { &mut end.reverse_prev };
            prev_unitigs.push(reverse_connection);
        }
    }

    pub fn clear_positions(&mut self) {
        for u in &self.unitigs {
            u.borrow_mut().clear_positions();
        }
    }

    pub fn max_unitig_number(&self) -> u32 {
        self.unitigs.iter().map(|u| u.borrow().number).max().unwrap_or(0)
    }

    pub fn connected_components(&self) -> Vec<Vec<u32>> {
        let mut visited = HashSet::new();
        let mut components = Vec::new();
        for unitig in &self.unitigs {
            let unitig_num = unitig.borrow().number;
            if !visited.contains(&unitig_num) {
                let mut component = Vec::new();
                self.dfs(unitig_num, &mut visited, &mut component);
                component.sort();
                components.push(component);
            }
        }
        components.sort();
        components
    }

    fn dfs(&self, unitig_num: u32, visited: &mut HashSet<u32>, component: &mut Vec<u32>) {
        let mut stack = vec![unitig_num];
        while let Some(current) = stack.pop() {
            if visited.insert(current) {
                component.push(current);
                for neighbor in self.connected_unitigs(current) {
                    if !visited.contains(&neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }
    }

    fn connected_unitigs(&self, unitig_num: u32) -> HashSet<u32> {
        // Given a unitig (by number), this function returns the unitigs (by number) it is directly
        // connected to.
        let mut connections = HashSet::new();
        if let Some(unitig_rc) = self.unitig_index.get(&unitig_num) {
            let unitig = unitig_rc.borrow();
            for c in &unitig.forward_next { connections.insert(c.number()); }
            for c in &unitig.forward_prev { connections.insert(c.number()); }
            for c in &unitig.reverse_next { connections.insert(c.number()); }
            for c in &unitig.reverse_prev { connections.insert(c.number()); }
        }
        connections
    }

    pub fn component_is_circular_loop(&self, component: &[u32]) -> bool {
        // Given a connected component of the graph, this function returns whether or not it forms
        // a simple circular loop.
        if component.is_empty() { return false; }
        let first = component[0];
        let mut num = first;
        let mut strand = strand::FORWARD;
        let mut visited = HashSet::new();
        while num != first || visited.is_empty() {
            if !visited.insert(num) { return false; }
            let unitig = self.unitig_index.get(&num).unwrap().borrow();
            if unitig.forward_next.len() != 1 || unitig.forward_prev.len() != 1 ||
               unitig.reverse_next.len() != 1 || unitig.reverse_prev.len() != 1 { return false; }
            let next = if strand { &unitig.forward_next[0] } else { &unitig.reverse_next[0] };
            num = next.number();
            strand = next.strand;
        }
        visited.len() == component.len()
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
    use crate::test_gfa::*;
    use super::*;

    #[test]
    fn test_graph_stats() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        graph.check_links();
        assert_eq!(graph.k_size, 9);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (21, 11));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        graph.check_links();
        assert_eq!(graph.k_size, 9);
        assert_eq!(graph.unitigs.len(), 3);
        assert_eq!(graph.total_length(), 31);
        assert_eq!(graph.link_count(), (8, 4));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
        graph.check_links();
        assert_eq!(graph.k_size, 9);
        assert_eq!(graph.unitigs.len(), 7);
        assert_eq!(graph.total_length(), 85);
        assert_eq!(graph.link_count(), (15, 8));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_4());
        graph.check_links();
        assert_eq!(graph.k_size, 3);
        assert_eq!(graph.unitigs.len(), 5);
        assert_eq!(graph.total_length(), 43);
        assert_eq!(graph.link_count(), (10, 5));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        graph.check_links();
        assert_eq!(graph.k_size, 3);
        assert_eq!(graph.unitigs.len(), 6);
        assert_eq!(graph.total_length(), 60);
        assert_eq!(graph.link_count(), (8, 4));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_6());
        graph.check_links();
        assert_eq!(graph.k_size, 3);
        assert_eq!(graph.unitigs.len(), 2);
        assert_eq!(graph.total_length(), 34);
        assert_eq!(graph.link_count(), (2, 1));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_7());
        graph.check_links();
        assert_eq!(graph.k_size, 3);
        assert_eq!(graph.unitigs.len(), 2);
        assert_eq!(graph.total_length(), 34);
        assert_eq!(graph.link_count(), (2, 1));
    }

    #[test]
    fn test_parse_unitig_path() {
        assert_eq!(parse_unitig_path("2+,1-"), vec![(2, strand::FORWARD), (1, strand::REVERSE)]);
        assert_eq!(parse_unitig_path("3+,8-,4-"), vec![(3, strand::FORWARD), (8, strand::REVERSE), (4, strand::REVERSE)]);
    }

    #[test]
    fn test_reverse_path() {
        assert_eq!(reverse_path(&[(1, strand::FORWARD), (2, strand::REVERSE)]),
                             vec![(2, strand::FORWARD), (1, strand::REVERSE)]);
        assert_eq!(reverse_path(&[(4, strand::FORWARD), (8, strand::FORWARD), (3, strand::REVERSE)]),
                             vec![(3, strand::FORWARD), (8, strand::REVERSE), (4, strand::REVERSE)]);
    }

    #[test]
    fn test_link_exists_1() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
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
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
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
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
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
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
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
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
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
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
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

    #[test]
    fn test_max_unitig_number() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        assert_eq!(graph.max_unitig_number(), 10);

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        assert_eq!(graph.max_unitig_number(), 3);

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
        assert_eq!(graph.max_unitig_number(), 7);
    }

    #[test]
    fn test_delete_link_and_create_link() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());

        graph.delete_link(-3, 1);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (19, 10));

        graph.delete_link(6, -6);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (18, 9));

        graph.delete_link(5, 6);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (16, 8));

        graph.delete_link(-1, 7);  // link doesn't exist, should do nothing
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (16, 8));

        graph.create_link(5, 6);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (18, 9));

        graph.create_link(6, -6);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (19, 10));

        graph.create_link(-3, 1);
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (21, 11));
    }

    #[test]
    fn test_get_sequence_from_path() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());

        assert_eq!(graph.get_sequence_from_path(&[(10, true), (8, false), (4, false), (1, false), (3, true)]),
                   "TAGATCGAGCCGAGCAAAGCGAAGCGAGCGCAGCGAATGCCTGAATCGCCTA".to_string());
        assert_eq!(graph.get_sequence_from_path(&[(5, true), (6, true), (6, false), (5, false)]),
                   "CGAACCATTACTTGTACAAGTAATGGTTCG".to_string());
        assert_eq!(graph.get_sequence_from_path(&[(3, false), (1, true), (4, true), (7, false), (9, false), (7, true), (4, false), (1, false), (2, false)]),
                   "TAGGCGATTCAGGCATTCGCTGCGCTCGCTTCGCTTTGCTCGGCTCGAAGGCGCGCCTTCGAGCCGAGCAAAGCGAAGCGAGCGCAGCGAATGCACAGCGACGACGGCA".to_string());

        assert_eq!(graph.get_sequence_from_path_signed(&[10, -8, -4, -1, 3]),
                   "TAGATCGAGCCGAGCAAAGCGAAGCGAGCGCAGCGAATGCCTGAATCGCCTA".as_bytes());
        assert_eq!(graph.get_sequence_from_path_signed(&[5, 6, -6, -5]),
                   "CGAACCATTACTTGTACAAGTAATGGTTCG".as_bytes());
        assert_eq!(graph.get_sequence_from_path_signed(&[-3, 1, 4, -7, -9, 7, -4, -1, -2]),
                   "TAGGCGATTCAGGCATTCGCTGCGCTCGCTTCGCTTTGCTCGGCTCGAAGGCGCGCCTTCGAGCCGAGCAAAGCGAAGCGAGCGCAGCGAATGCACAGCGACGACGGCA".as_bytes());
    }

    #[test]
    fn test_connected_components() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]);

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3]]);

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3, 4, 5, 6, 7]]);

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_4());
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3], vec![4, 5]]);

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        assert_eq!(graph.connected_components(), vec![vec![1, 5], vec![2], vec![3, 6], vec![4]]);
    }

    #[test]
    fn test_component_is_circular_loop() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        assert!(!graph.component_is_circular_loop(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        assert!(!graph.component_is_circular_loop(&[1, 2, 3]));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
        assert!(!graph.component_is_circular_loop(&[1, 2, 3, 4, 5, 6, 7]));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_4());
        assert!(graph.component_is_circular_loop(&[1, 2, 3]));
        assert!(graph.component_is_circular_loop(&[3, 2, 1]));
        assert!(graph.component_is_circular_loop(&[2, 3, 1]));
        assert!(graph.component_is_circular_loop(&[4, 5]));
        assert!(graph.component_is_circular_loop(&[5, 4]));

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        assert!(!graph.component_is_circular_loop(&[1, 5]));
        assert!(!graph.component_is_circular_loop(&[2]));
        assert!(!graph.component_is_circular_loop(&[3, 6]));
        assert!(graph.component_is_circular_loop(&[4]));
        assert!(!graph.component_is_circular_loop(&[]));
    }

    #[test]
    fn test_delete_link_break_into_components() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_6());
        assert_eq!(graph.connected_components(), vec![vec![1, 2]]);
        graph.delete_link(1, -2);
        assert_eq!(graph.connected_components(), vec![vec![1], vec![2]]);

        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_7());
        assert_eq!(graph.connected_components(), vec![vec![1, 2]]);
        graph.delete_link(-1, 2);
        assert_eq!(graph.connected_components(), vec![vec![1], vec![2]]);

        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]);
        graph.delete_link(4, 8);
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3, 4, 5, 6, 7, 9], vec![8, 10]]);
        graph.delete_link(-3, 1);
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 4, 5, 6, 7, 9], vec![3], vec![8, 10]]);
    }

    #[test]
    fn test_remove_unitigs_by_number() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]);
        graph.remove_unitigs_by_number(HashSet::from([4, 5]));
        assert_eq!(graph.connected_components(), vec![vec![1, 2, 3], vec![6], vec![7, 9], vec![8, 10]]);
        graph.remove_unitigs_by_number(HashSet::from([1]));
        assert_eq!(graph.connected_components(), vec![vec![2], vec![3], vec![6], vec![7, 9], vec![8, 10]]);
    }

    #[test]
    fn test_is_isolated_and_circular() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        for unitig in &graph.unitigs {
            assert!(!unitig.borrow().is_isolated_and_circular());
        }

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        for unitig in &graph.unitigs {
            assert!(!unitig.borrow().is_isolated_and_circular());
        }

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
        for unitig in &graph.unitigs {
            assert!(!unitig.borrow().is_isolated_and_circular());
        }

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_4());
        for unitig in &graph.unitigs {
            assert!(!unitig.borrow().is_isolated_and_circular());
        }

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        assert!(!graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_circular());
        assert!(!graph.unitig_index.get(&2).unwrap().borrow().is_isolated_and_circular());
        assert!(!graph.unitig_index.get(&3).unwrap().borrow().is_isolated_and_circular());
        assert!(graph.unitig_index.get(&4).unwrap().borrow().is_isolated_and_circular());
        assert!(!graph.unitig_index.get(&5).unwrap().borrow().is_isolated_and_circular());
        assert!(!graph.unitig_index.get(&6).unwrap().borrow().is_isolated_and_circular());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_6());
        for unitig in &graph.unitigs {
            assert!(!unitig.borrow().is_isolated_and_circular());
        }

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_7());
        for unitig in &graph.unitigs {
            assert!(!unitig.borrow().is_isolated_and_circular());
        }
    }

    #[test]
    fn test_topology() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        assert_eq!(graph.topology(), "fragmented".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        assert_eq!(graph.topology(), "fragmented".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_3());
        assert_eq!(graph.topology(), "fragmented".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_4());
        assert_eq!(graph.topology(), "fragmented".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        assert_eq!(graph.topology(), "fragmented".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_6());
        assert_eq!(graph.topology(), "fragmented".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_7());
        assert_eq!(graph.topology(), "fragmented".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());
        assert_eq!(graph.topology(), "circular".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        assert_eq!(graph.topology(), "linear_blunt_blunt".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_10());
        assert_eq!(graph.topology(), "linear_hairpin_hairpin".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_11());
        assert_eq!(graph.topology(), "linear_blunt_hairpin".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_12());
        assert_eq!(graph.topology(), "linear_blunt_hairpin".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_13());
        assert_eq!(graph.topology(), "other".to_string());
    }

    #[test]
    fn test_duplicate_unitig_by_number() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_4());
        assert_eq!(graph.unitigs.len(), 5);
        assert_eq!(graph.total_length(), 43);
        assert_eq!(graph.link_count(), (10, 5));
        graph.duplicate_unitig_by_number(&5);
        assert_eq!(graph.unitigs.len(), 6);
        assert_eq!(graph.total_length(), 46);
        assert_eq!(graph.link_count(), (10, 5));

        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        assert_eq!(graph.unitigs.len(), 6);
        assert_eq!(graph.total_length(), 60);
        assert_eq!(graph.link_count(), (8, 4));
        graph.duplicate_unitig_by_number(&1);
        assert_eq!(graph.unitigs.len(), 7);
        assert_eq!(graph.total_length(), 79);
        assert_eq!(graph.link_count(), (8, 4));
    }
}
