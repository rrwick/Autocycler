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
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use std::rc::Rc;

use crate::kmer_graph::KmerGraph;
use crate::log::{section_header, explanation};
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

    pub fn load_gfa_with_summary(gfa_filename: &Path) -> (Self, Vec<Sequence>) {
        section_header("Loading graph");
        explanation("The unitig graph is now loaded into memory.");
        let (graph, sequences) = Self::from_gfa_file(gfa_filename);
        graph.print_basic_graph_info();
        (graph, sequences)
    }

    pub fn from_gfa_lines(gfa_lines: &[String]) -> (Self, Vec<Sequence>) {
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
        self.unitig_index = self.unitigs.iter().map(|u| (u.borrow().number, Rc::clone(u))).collect();
    }

    fn read_gfa_header_line(&mut self, parts: &[&str]) {
        if let Some(k) = parts.iter().find_map(|p| p.strip_prefix("KM:i:")?.parse().ok()) {
            self.k_size = k;
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
            let get_unitig = |number| self.unitig_index.get(&number).unwrap_or_else(|| {
                quit_with_error(&format!("link refers to nonexistent unitig: {number}"));
            });
            add_link(get_unitig(seg_1), strand_1, get_unitig(seg_2), strand_2);
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
            let (Some(length), Some(filename), Some(header)) = (length, filename, header) else {
                quit_with_error("missing required tag in GFA path line.");
            };
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
                quit_with_error(&format!("unitig {unitig_num} not found in unitig index"));
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

            // Walking forward on the reverse strand extends the start of the unitig.
            for (mut kmer, mut reverse, strand) in [(forward_kmer, reverse_kmer, strand::FORWARD),
                                                   (reverse_kmer, forward_kmer, strand::REVERSE)] {
                loop {
                    if reverse.first_position() { break; }
                    let next_kmers = k_graph.next_kmers(kmer.seq());
                    if next_kmers.len() != 1 { break; }
                    kmer = next_kmers[0];
                    if seen.contains(kmer.seq()) { break; }
                    if k_graph.prev_kmers(kmer.seq()).len() != 1 { break; }
                    reverse = k_graph.reverse(kmer);
                    if kmer.first_position() { break; }
                    if strand {
                        unitig.add_kmer_to_end(kmer, reverse);
                    } else {
                        unitig.add_kmer_to_start(reverse, kmer);
                    }
                    seen.insert(kmer.seq());
                    seen.insert(reverse.seq());
                }
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
            let unitig = unitig.borrow();
            let forward_key = unitig.forward_seq[..piece_len].to_vec();
            let reverse_key = unitig.reverse_seq[..piece_len].to_vec();
            forward_starts.entry(forward_key).or_insert_with(Vec::new).push(i);
            reverse_starts.entry(reverse_key).or_insert_with(Vec::new).push(i);
        }

        // Use the indices to find connections between unitigs.
        for unitig_a in &self.unitigs {
            let (ending_forward_seq, ending_reverse_seq) = {
                let unitig = unitig_a.borrow();
                (unitig.forward_seq[unitig.forward_seq.len() - piece_len..].to_vec(),
                 unitig.reverse_seq[unitig.reverse_seq.len() - piece_len..].to_vec())
            };

            if let Some(next_idxs) = forward_starts.get(&ending_forward_seq) {
                for &j in next_idxs {
                    let unitig_b = &self.unitigs[j];
                    add_link(unitig_a, strand::FORWARD, unitig_b, strand::FORWARD);
                    add_link(unitig_b, strand::REVERSE, unitig_a, strand::REVERSE);
                }
            }

            if let Some(next_idxs) = reverse_starts.get(&ending_forward_seq) {
                for &j in next_idxs {
                    add_link(unitig_a, strand::FORWARD, &self.unitigs[j], strand::REVERSE);
                }
            }

            if let Some(next_idxs) = forward_starts.get(&ending_reverse_seq) {
                for &j in next_idxs {
                    add_link(unitig_a, strand::REVERSE, &self.unitigs[j], strand::FORWARD);
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
        self.unitigs.sort_by(|a_rc, b_rc| {
            let a = a_rc.borrow();
            let b = b_rc.borrow();
            b.length().cmp(&a.length())
                .then_with(|| a.forward_seq.cmp(&b.forward_seq))
                .then_with(|| b.depth.partial_cmp(&a.depth).unwrap_or(std::cmp::Ordering::Equal))
        });
        for (new_number, unitig) in self.unitigs.iter().enumerate() {
            unitig.borrow_mut().number = (new_number + 1) as u32;
        }
        self.build_unitig_index();
    }

    pub fn save_gfa(&self, gfa_filename: &Path, sequences: &[Sequence],
                    use_other_colour: bool) -> io::Result<()> {
        let mut file = File::create(gfa_filename)?;
        writeln!(file, "H\tVN:Z:1.0\tKM:i:{}", self.k_size)?;
        for unitig in &self.unitigs {
            writeln!(file, "{}", unitig.borrow().gfa_segment_line(use_other_colour))?;
        }
        for (a, a_strand, b, b_strand) in self.get_links_for_gfa(0) {
            writeln!(file, "L\t{a}\t{a_strand}\t{b}\t{b_strand}\t0M")?;
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
            for (a_strand, next) in [("+", &a.forward_next), ("-", &a.reverse_next)] {
                for b in next {
                    let b_num = b.number() + offset;
                    links.push((a_num.to_string(), a_strand.to_string(), b_num.to_string(),
                                (if b.strand {"+"} else {"-"}).to_string()));
                }
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

    pub fn reconstruct_original_sequences(&self, seqs: &[Sequence])
            -> BTreeMap<String, Vec<(String, String)>> {
        let mut original_seqs = BTreeMap::new();
        for seq in seqs {
            let sequence = self.reconstruct_sequence(seq);
            original_seqs.entry(seq.filename.clone()).or_insert_with(Vec::new)
                .push((seq.contig_header.clone(), sequence));
        }
        original_seqs
    }

    pub fn reconstruct_original_sequences_u8(&self, seqs: &[Sequence])
            -> Vec<((String, String), Vec<u8>)> {
        let mut original_seqs: Vec<_> = seqs.iter().map(|seq| {
            let sequence = self.reconstruct_sequence(seq);
            ((seq.filename.clone(), seq.contig_name()), sequence.into_bytes())
        }).collect();
        original_seqs.sort();
        original_seqs
    }

    fn reconstruct_sequence(&self, seq: &Sequence) -> String {
        let path = self.get_unitig_path_for_sequence(seq);
        let sequence = self.get_sequence_from_path(&path);
        assert_eq!(sequence.len(), seq.length, "reconstructed sequence does not have expected length");
        sequence
    }

    fn get_sequence_from_path(&self, path: &[(u32, bool)]) -> String {
        let mut sequence = String::new();
        for (unitig_num, strand) in path {
            let unitig = self.unitig_index.get(unitig_num).unwrap().borrow();
            sequence.push_str(std::str::from_utf8(unitig.get_seq(*strand)).unwrap());
        }
        sequence
    }

    pub fn get_sequence_from_path_signed(&self, path: &[i32]) -> Vec<u8> {
        let path: Vec<_> = path.iter().map(|&x| (x.unsigned_abs(), x >= 0)).collect();
        self.get_sequence_from_path(&path).into_bytes()
    }

    fn find_starting_unitig(&self, seq_id: u16) -> UnitigStrand {
        let mut starting_unitigs = Vec::new();
        for unitig_rc in &self.unitigs {
            let unitig = unitig_rc.borrow();
            for (strand, positions) in [(strand::FORWARD, &unitig.forward_positions),
                                        (strand::REVERSE, &unitig.reverse_positions)] {
                for p in positions {
                    if p.seq_id() == seq_id && p.strand() && p.pos == 0 {
                        starting_unitigs.push(UnitigStrand::new(unitig_rc, strand));
                    }
                }
            }
        }
        assert_eq!(starting_unitigs.len(), 1);
        starting_unitigs[0].clone()
    }

    pub fn get_next_unitig(&self, seq_id: u16, seq_strand: bool, unitig_rc: &Rc<RefCell<Unitig>>,
                           strand: bool, pos: u32) -> Option<(UnitigStrand, u32)> {
        let unitig = unitig_rc.borrow();
        let next_pos = pos + unitig.length();
        let next_edges = if strand { &unitig.forward_next } else { &unitig.reverse_next };
        for next in next_edges {
            let next_rc = next.unitig.upgrade()?;
            let next_u  = next_rc.borrow();
            let positions = if next.strand { &next_u.forward_positions }
                                      else { &next_u.reverse_positions };
            if positions.iter().any(|p|
                p.seq_id() == seq_id && p.strand() == seq_strand && p.pos == next_pos){
                return Some((next.clone(), next_pos));
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
            let Some(current_rc) = u.unitig.upgrade() else { break };
            let Some(next) = self.get_next_unitig(seq.id, strand::FORWARD, &current_rc, u.strand, pos)
                else { break };
            (u, pos) = next;
        }
        unitig_path
    }

    pub fn get_unitig_path_for_sequence_i32(&self, seq: &Sequence) -> Vec<i32> {
        let unitig_path = self.get_unitig_path_for_sequence(seq);
        unitig_path.iter().map(|(u, s)| if *s { *u as i32 } else { -(*u as i32)}).collect()
    }

    pub fn total_length(&self) -> u64 {
        self.unitigs.iter().map(|u| u.borrow().length() as u64).sum()
    }

    pub fn link_count(&self) -> (usize, usize) {
        // Counts both orientations and unique links; hairpins have only one orientation.
        let mut all_links = HashSet::new();
        let mut one_way_links = HashSet::new();
        for a_rc in &self.unitigs {
            let a = a_rc.borrow();
            let a_num = a.number as i32;
            for (direction, next) in [(1, &a.forward_next), (-1, &a.reverse_next)] {
                for b in next {
                    let link = (a_num * direction, b.signed_number());
                    let rev_link = (-link.1, -link.0);
                    all_links.insert(link);
                    all_links.insert(rev_link);
                    one_way_links.insert(link.max(rev_link));
                }
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

    pub fn print_basic_graph_info_with_topology(&self) {
        let link_count = self.link_count().1;
        eprintln!("{} unitig{}, {} link{} ({})",
                  self.unitigs.len(), match self.unitigs.len() { 1 => "", _ => "s" },
                  link_count, match link_count { 1 => "", _ => "s" }, self.topology());
        eprintln!("total length: {} bp", self.total_length());
        eprintln!();
    }

    pub fn topology(&self) -> String {
        if self.unitigs.is_empty() { return "empty".to_string(); }
        if self.unitigs.len() > 1 { return "fragmented".to_string(); }
        let u = self.unitigs[0].borrow();
        if self.link_count().0 == 0 { return "linear-open-open".to_string(); }
        if u.is_isolated_and_circular() { return "circular".to_string(); }
        if u.hairpin_start() && u.hairpin_end() { return "linear-hairpin-hairpin".to_string(); }
        if u.hairpin_start() && u.open_end() { return "linear-open-hairpin".to_string(); }
        if u.open_start() && u.hairpin_end() { return "linear-open-hairpin".to_string(); }
        "other".to_string()
    }

    pub fn delete_dangling_links(&mut self) {
        // Run after removing unitigs, before rebuilding the index that keeps them alive.
        let unitig_numbers: HashSet<u32> = self.unitigs.iter().map(|u| u.borrow().number).collect();
        for unitig_rc in &self.unitigs {
            let unitig = unitig_rc.borrow();
            let keep = [&unitig.forward_next, &unitig.forward_prev,
                        &unitig.reverse_next, &unitig.reverse_prev].map(|links| {
                links.iter().map(|link| unitig_numbers.contains(&link.number())).collect()
            });
            drop(unitig);
            let unitig = &mut *unitig_rc.borrow_mut();
            for (links, keep) in [&mut unitig.forward_next, &mut unitig.forward_prev,
                                 &mut unitig.reverse_next, &mut unitig.reverse_prev].into_iter().zip(keep) {
                retain_links(links, keep);
            }
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
        // Each copy keeps one of the two non-self links, plus all loops and hairpins.
        self.check_if_unitig_can_be_duplicated(unitig_num);

        let u = self.unitig_index.get(unitig_num).unwrap().clone();
        let target = u.borrow();
        let target_num = target.number;
        let a_num = self.max_unitig_number() + 1;
        let b_num = a_num + 1;
        for number in [a_num, b_num] {
            let mut copy = target.clone();
            copy.number = number;
            copy.depth /= 2.0;
            copy.clear_all_links();
            self.unitigs.push(Rc::new(RefCell::new(copy)));
        }

        self.remove_unitigs_by_number(std::iter::once(*unitig_num).collect());

        let mut non_self_links = Vec::new();
        for (direction, next) in [(1, &target.forward_next), (-1, &target.reverse_next)] {
            for link in next {
                if link.number() == *unitig_num {
                    for number in [a_num, b_num] {
                        // Both orientations are already present in the target's links.
                        self.create_link_one_way(number as i32 * direction,
                                                 number as i32 * if link.strand { 1 } else { -1 });
                    }
                } else {
                    non_self_links.push((target_num as i32 * direction, link.signed_number()));
                }
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
        let link_count = unitig.forward_next.iter().chain(&unitig.reverse_next)
            .filter(|link| link.number() != unitig_num).count();
        if link_count != 2 {
            quit_with_error(&format!("unitig {unitig_num} does not contain exactly two non-self links"));
        }
    }

    pub fn remove_low_depth_unitigs(&mut self, min_depth: f64) {
        // Work backwards to preferentially keep longer unitigs in a sorted graph.
        for idx in (0..self.unitigs.len()).rev() {
            let Some(unitig_rc) = self.unitigs.get(idx) else { continue };
            let unitig = unitig_rc.borrow();
            if unitig.depth > min_depth || removal_creates_dead_end(&unitig) { continue; }
            let number = unitig.number;
            drop(unitig);
            self.unitigs.retain(|u| u.borrow().number != number);
            self.delete_dangling_links();
            self.build_unitig_index();
        }
    }

    pub fn link_exists(&self, a_num: u32, a_strand: bool, b_num: u32, b_strand: bool) -> bool {
        let Some(unitig_a) = self.unitig_index.get(&a_num) else { return false };
        let unitig_a = unitig_a.borrow();
        let next_links = if a_strand { &unitig_a.forward_next } else { &unitig_a.reverse_next };
        next_links.iter().any(|next| next.number() == b_num && next.strand == b_strand)
    }

    pub fn link_exists_prev(&self, a_num: u32, a_strand: bool, b_num: u32, b_strand: bool) -> bool {
        let Some(unitig_b) = self.unitig_index.get(&b_num) else { return false };
        let unitig_b = unitig_b.borrow();
        let prev_links = if b_strand { &unitig_b.forward_prev } else { &unitig_b.reverse_prev };
        prev_links.iter().any(|prev| prev.number() == a_num && prev.strand == a_strand)
    }

    pub fn check_links(&self) {
        for a_rc in &self.unitigs {
            let a = a_rc.borrow();
            for (a_strand, next) in [(strand::FORWARD, &a.forward_next),
                                     (strand::REVERSE, &a.reverse_next)] {
                for b in next {
                    self.check_link(a.number, a_strand, b.number(), b.strand);
                    assert!(self.unitig_index.contains_key(&b.number()), "unitig missing from index");
                }
            }
            for (a_strand, prev) in [(strand::FORWARD, &a.forward_prev),
                                     (strand::REVERSE, &a.reverse_prev)] {
                for b in prev {
                    self.check_link(b.number(), b.strand, a.number, a_strand);
                    assert!(self.unitig_index.contains_key(&b.number()), "unitig missing from index");
                }
            }
        }
    }

    fn check_link(&self, a_num: u32, a_strand: bool, b_num: u32, b_strand: bool) {
        // Both orientations must have matching next/prev entries.
        assert!(self.link_exists(a_num, a_strand, b_num, b_strand), "missing next link");
        assert!(self.link_exists_prev(a_num, a_strand, b_num, b_strand), "missing prev link");
        assert!(self.link_exists(b_num, !b_strand, a_num, !a_strand), "missing next link");
        assert!(self.link_exists_prev(b_num, !b_strand, a_num, !a_strand), "missing prev link");
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
        let start_strand = start_num > 0;
        let end_strand = end_num > 0;
        let start_num = start_num.unsigned_abs();
        let end_num = end_num.unsigned_abs();
        let start_rc = self.unitig_index.get(&start_num).unwrap();
        let end_rc = self.unitig_index.get(&end_num).unwrap();

        let keep_next = {
            let start = start_rc.borrow();
            let next_unitigs = if start_strand { &start.forward_next } else { &start.reverse_next };
            next_unitigs.iter().map(|link| link.number() != end_num || link.strand != end_strand).collect()
        };
        {
            let mut start = start_rc.borrow_mut();
            let next_unitigs = if start_strand { &mut start.forward_next } else { &mut start.reverse_next };
            retain_links(next_unitigs, keep_next);
        }
        let keep_prev = {
            let end = end_rc.borrow();
            let prev_unitigs = if end_strand { &end.forward_prev } else { &end.reverse_prev };
            prev_unitigs.iter().map(|link| link.number() != start_num || link.strand != start_strand).collect()
        };
        let mut end = end_rc.borrow_mut();
        let prev_unitigs = if end_strand { &mut end.forward_prev } else { &mut end.reverse_prev };
        retain_links(prev_unitigs, keep_prev);
    }

    pub fn create_link(&mut self, start_num: i32, end_num: i32) {
        self.create_link_one_way(start_num, end_num);
        if start_num != -end_num {
            self.create_link_one_way(-end_num, -start_num);
        }
    }

    fn create_link_one_way(&mut self, start_num: i32, end_num: i32) {
        let start_rc = self.unitig_index.get(&start_num.unsigned_abs()).unwrap();
        let end_rc = self.unitig_index.get(&end_num.unsigned_abs()).unwrap();
        add_link(start_rc, start_num > 0, end_rc, end_num > 0);
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
        let Some(unitig_rc) = self.unitig_index.get(&unitig_num) else { return HashSet::new() };
        let unitig = unitig_rc.borrow();
        [&unitig.forward_next, &unitig.forward_prev, &unitig.reverse_next, &unitig.reverse_prev]
            .into_iter().flatten().map(UnitigStrand::number).collect()
    }

    pub fn component_is_circular_loop(&self, component: &[u32]) -> bool {
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


fn retain_links(links: &mut Vec<UnitigStrand>, keep: Vec<bool>) {
    // Decide which links to keep before borrowing mutably: a link can point to its own unitig.
    let mut keep = keep.into_iter();
    links.retain(|_| keep.next().unwrap());
}


fn add_link(start_rc: &Rc<RefCell<Unitig>>, start_strand: bool,
            end_rc: &Rc<RefCell<Unitig>>, end_strand: bool) {
    // Borrow each end separately because links can connect a unitig to itself.
    {
        let mut start = start_rc.borrow_mut();
        let next = if start_strand { &mut start.forward_next } else { &mut start.reverse_next };
        next.push(UnitigStrand::new(end_rc, end_strand));
    }
    let mut end = end_rc.borrow_mut();
    let prev = if end_strand { &mut end.forward_prev } else { &mut end.reverse_prev };
    prev.push(UnitigStrand::new(start_rc, start_strand));
}


fn removal_creates_dead_end(unitig: &Unitig) -> bool {
    for (outgoing, links) in [(true, &unitig.forward_next), (false, &unitig.forward_prev)] {
        for link in links {
            let neighbor_rc = link.unitig();
            let neighbor = neighbor_rc.borrow();
            if neighbor.number == unitig.number { continue; }
            let alternatives = match (outgoing, link.strand) {
                (true, true) => &neighbor.forward_prev,
                (true, false) => &neighbor.reverse_prev,
                (false, true) => &neighbor.forward_next,
                (false, false) => &neighbor.reverse_next,
            };
            if !alternatives.iter().any(|link| link.number() != unitig.number) { return true; }
        }
    }
    false
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
    use crate::graph_simplification::merge_linear_paths;
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
    fn test_self_link_order_and_duplicates() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        for (start, end) in [(1, 1), (1, -1), (-1, 1), (1, -1)] {
            graph.create_link(start, end);
        }
        graph.check_links();
        let expected: Vec<_> = [("+", "+"), ("+", "-"), ("+", "-"), ("-", "-"), ("-", "+")]
            .into_iter().map(|(a, b)| ("11".to_string(), a.to_string(), "11".to_string(), b.to_string()))
            .collect();
        assert_eq!(graph.get_links_for_gfa(10), expected);
        assert_eq!(graph.link_count(), (4, 3));
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("graph.gfa");
        graph.save_gfa(&path, &[], false).unwrap();
        let (reloaded, _) = UnitigGraph::from_gfa_file(&path);
        assert_eq!(reloaded.get_links_for_gfa(10), expected);
        graph.delete_link(1, -1);
        graph.check_links();
        assert_eq!(graph.link_count(), (3, 2));
    }

    #[test]
    fn test_dangling_link_removal_preserves_order_and_duplicates() {
        for remove_by_number in [true, false] {
            let lines = (1..=4).map(|n| format!("S\t{n}\tACGT\tDP:f:{}", n % 2)).collect::<Vec<_>>();
            let (mut graph, _) = UnitigGraph::from_gfa_lines(&lines);
            for start in [-4, -3, -2, -1, 1, 2, 3, 4] {
                for end in [1, -2, 3, 4, -1, 2, -3, -4, 3] {
                    graph.create_link(start, end);
                }
            }
            let expected = graph.get_links_for_gfa(0).into_iter()
                .filter(|(a, _, b, _)| (a == "1" || a == "3") && (b == "1" || b == "3"))
                .collect::<Vec<_>>();
            if remove_by_number {
                graph.remove_unitigs_by_number(HashSet::from([2, 4]));
            } else {
                graph.remove_zero_depth_unitigs();
            }
            graph.check_links();
            assert_eq!(graph.get_links_for_gfa(0), expected);
        }
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
    fn test_is_isolated_and_linear() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        assert!(!graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_linear());
        assert!(graph.unitig_index.get(&2).unwrap().borrow().is_isolated_and_linear());
        assert!(!graph.unitig_index.get(&3).unwrap().borrow().is_isolated_and_linear());
        assert!(!graph.unitig_index.get(&4).unwrap().borrow().is_isolated_and_linear());
        assert!(!graph.unitig_index.get(&5).unwrap().borrow().is_isolated_and_linear());
        assert!(!graph.unitig_index.get(&6).unwrap().borrow().is_isolated_and_linear());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());
        assert!(!graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_linear());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        assert!(graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_linear());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_10());
        assert!(graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_linear());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_11());
        assert!(graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_linear());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_12());
        assert!(graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_linear());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_13());
        assert!(!graph.unitig_index.get(&1).unwrap().borrow().is_isolated_and_linear());
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
        assert_eq!(graph.topology(), "linear-open-open".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_10());
        assert_eq!(graph.topology(), "linear-hairpin-hairpin".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_11());
        assert_eq!(graph.topology(), "linear-open-hairpin".to_string());

        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_12());
        assert_eq!(graph.topology(), "linear-open-hairpin".to_string());

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

    #[test]
    fn test_duplicate_unitig_with_self_links() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        for number in [2, 3] {
            graph.unitigs.push(Rc::new(RefCell::new(Unitig::bridge(number, b"A".to_vec(), 1.0))));
        }
        graph.build_unitig_index();
        for (start, end) in [(2, 1), (1, -3), (1, 1), (1, -1), (-1, 1)] {
            graph.create_link(start, end);
        }
        graph.duplicate_unitig_by_number(&1);
        assert!(!graph.unitig_index.contains_key(&1));
        for (number, next, reverse_next) in [(4, vec![4, -4, -3], vec![-4, 4]),
                                             (5, vec![5, -5], vec![-5, 5, -2])] {
            let unitig = graph.unitig_index[&number].borrow();
            assert_eq!(unitig.depth, 0.5);
            assert_eq!(unitig.forward_next.iter().map(|u| u.signed_number()).collect::<Vec<_>>(), next);
            assert_eq!(unitig.reverse_next.iter().map(|u| u.signed_number()).collect::<Vec<_>>(), reverse_next);
        }
        graph.remove_unitigs_by_number(HashSet::from([2, 3]));
        for number in [4, 5] {
            graph.delete_link(number, -number);
            graph.delete_link(-number, number);
            assert!(graph.unitig_index[&(number as u32)].borrow().is_isolated_and_circular());
            assert!(graph.component_is_circular_loop(&[number as u32]));
        }
        graph.check_links();
    }

    #[test]
    fn test_remove_low_depth_unitigs_keeps_alternative_path() {
        for directions in 0..16 {
            let signed: Vec<i32> = (1..=4)
                .map(|n| if directions & (1 << (n - 1)) == 0 { n } else { -n }).collect();
            let mut graph = UnitigGraph::default();
            for number in 1..=4 {
                let depth = if number == 1 || number == 4 { 2.0 } else { 1.0 };
                graph.unitigs.push(Rc::new(RefCell::new(Unitig::bridge(number, b"A".to_vec(), depth))));
            }
            graph.build_unitig_index();
            for (start, end) in [(0, 1), (1, 3), (0, 2), (2, 3)] {
                graph.create_link(signed[start], signed[end]);
            }
            graph.remove_low_depth_unitigs(1.0);
            graph.check_links();
            assert_eq!(graph.connected_components(), vec![vec![1, 2, 4]]);
            assert_eq!(graph.link_count(), (4, 2));
        }
    }

    #[test]
    fn test_remove_low_depth_unitigs() {
        let (mut graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        assert_eq!(graph.unitigs.len(), 10);
        assert_eq!(graph.total_length(), 92);
        assert_eq!(graph.link_count(), (21, 11));

        graph.remove_low_depth_unitigs(1.0);
        assert_eq!(graph.unitigs.len(), 8);
        assert_eq!(graph.total_length(), 70);
        assert_eq!(graph.link_count(), (16, 8));

        merge_linear_paths(&mut graph, &[]);
        assert_eq!(graph.unitigs.len(), 6);
        assert_eq!(graph.total_length(), 70);
        assert_eq!(graph.link_count(), (12, 6));

        graph.remove_low_depth_unitigs(1.0);
        assert_eq!(graph.unitigs.len(), 5);
        assert_eq!(graph.total_length(), 65);
        assert_eq!(graph.link_count(), (10, 5));
    }
}
