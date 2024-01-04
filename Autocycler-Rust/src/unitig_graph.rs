// This file defines the UnitigGraph struct for building a compacted unitig graph from a KmerGraph.

// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write, BufRead, BufReader};
use std::path::PathBuf;

use crate::kmer_graph::KmerGraph;
use crate::position::Position;
use crate::sequence::Sequence;
use crate::unitig::Unitig;


pub struct UnitigGraph {
    pub unitigs: Vec<Unitig>,
    k_size: u32,
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
                Some(&"H") => u_graph.read_gfa_header_line(&line),
                Some(&"S") => u_graph.unitigs.push(Unitig::from_segment_line(&line)),
                Some(&"L") => link_lines.push(line),
                Some(&"P") => path_lines.push(line),
                _ => {}
            }
        }
        u_graph.build_links_from_gfa(&link_lines);
        u_graph.build_paths_from_gfa(&path_lines);
        u_graph
    }

    fn read_gfa_header_line(&mut self, line: &str) {
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
    }

    fn build_links_from_gfa(&mut self, link_lines: &Vec<String>) {
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
    }

    fn build_paths_from_gfa(&mut self, path_lines: &Vec<String>) {
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
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
            self.unitigs.push(unitig);
        }
    }

    fn simplify_seqs(&mut self) {
        for unitig in &mut self.unitigs {
            unitig.simplify_seqs();
        }
    }

    fn create_links(&mut self) {
        let piece_len = self.k_size as usize - 1;

        // Index unitigs by their k-1 starting sequences.
        let mut forward_starts = HashMap::new();
        let mut reverse_starts = HashMap::new();
        for (i, unitig) in self.unitigs.iter().enumerate() {
            let forward_key = unitig.forward_seq[..piece_len].to_vec();
            let reverse_key = unitig.reverse_seq[..piece_len].to_vec();
            forward_starts.entry(forward_key).or_insert_with(Vec::new).push(i);
            reverse_starts.entry(reverse_key).or_insert_with(Vec::new).push(i);
        }

        // Use the indices to find connections between unitigs.
        self.link_count = 0;
        for i in 0..self.unitigs.len() {
            unsafe {
                let unitig_a = self.unitigs.get_unchecked_mut(i) as *mut Unitig;
                let ending_forward_seq = self.unitigs[i].forward_seq[self.unitigs[i].len() - piece_len..].to_vec();
                let ending_reverse_seq = self.unitigs[i].reverse_seq[self.unitigs[i].len() - piece_len..].to_vec();

                if let Some(next_idxs) = forward_starts.get(&ending_forward_seq) {
                    for &j in next_idxs {
                        let unitig_b = self.unitigs.get_unchecked_mut(j) as *mut Unitig;

                        // unitig_a+ -> unitig_b+
                        (*unitig_a).forward_next.push((unitig_b, true));
                        (*unitig_b).forward_prev.push((unitig_a, true));
                        self.link_count += 1;

                        // unitig_b- -> unitig_a-
                        (*unitig_b).reverse_next.push((unitig_a, false));
                        (*unitig_a).reverse_prev.push((unitig_b, false));
                        self.link_count += 1;
                    }
                }

                if let Some(next_idxs) = reverse_starts.get(&ending_forward_seq) {
                    for &j in next_idxs {
                        let unitig_b = self.unitigs.get_unchecked_mut(j) as *mut Unitig;

                        // unitig_a+ -> unitig_b-
                        (*unitig_a).forward_next.push((unitig_b, false));
                        (*unitig_b).reverse_prev.push((unitig_a, true));
                        self.link_count += 1;
                    }
                }

                if let Some(next_idxs) = forward_starts.get(&ending_reverse_seq) {
                    for &j in next_idxs {
                        let unitig_b = self.unitigs.get_unchecked_mut(j) as *mut Unitig;

                        // unitig_a- -> unitig_b+
                        (*unitig_a).reverse_next.push((unitig_b, true));
                        (*unitig_b).forward_prev.push((unitig_a, false));
                        self.link_count += 1;
                    }
                }
            }
        }
    }

    fn trim_overlaps(&mut self) {
        for unitig in &mut self.unitigs {
            unitig.trim_overlaps(self.k_size as usize);
        }
    }

    fn renumber_unitigs(&mut self) {
        // This method reorders Unitigs by: length (decreasing), sequence (lexicographic) and
        // depth (decreasing). Importantly, it does not sort the self.unitigs vector, as that would
        // invalidate the raw pointers between Unitig objects.
        let mut unitig_data: Vec<(f64, Vec<u8>, f64, usize)> = self.unitigs.iter().enumerate()
            .map(|(i, unitig)| {
                let length_inverse = 1.0 / (unitig.length() as f64);
                let seq_clone = unitig.forward_seq.clone();
                let depth_inverse = 1.0 / unitig.depth;
                (length_inverse, seq_clone, depth_inverse, i)
            }).collect();
        unitig_data.sort_by(|(a_length, a_seq, a_depth, _), (b_length, b_seq, b_depth, _)| {
            (a_length, a_seq, a_depth).partial_cmp(&(b_length, b_seq, b_depth)).unwrap_or(std::cmp::Ordering::Equal)
        });
        for (new_number, (_, _, _, i)) in unitig_data.into_iter().enumerate() {
            self.unitigs[i].number = (new_number + 1) as u32;
        }
    }

    pub fn iterate_unitigs(&self) -> impl Iterator<Item = &Unitig> {
        // This method allows for iterating over the Unitigs in their number order, despite the
        // fact that the self.unitigs vector is not sorted in number order.
        let mut unitig_refs: Vec<&Unitig> = self.unitigs.iter().collect();
        unitig_refs.sort_by_key(|u| u.number);
        unitig_refs.into_iter()
    }

    pub fn save_gfa(&self, gfa_filename: &PathBuf, sequences: &Vec<Sequence>) -> io::Result<()> {
        let mut file = File::create(gfa_filename)?;
        writeln!(file, "H\tVN:Z:1.0\tKM:i:{}", self.k_size)?;
        for unitig in self.iterate_unitigs() {
            writeln!(file, "{}", unitig.gfa_segment_line())?;
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
        for a in self.iterate_unitigs() {
            unsafe {
                for &(b_ptr, b_strand) in &a.forward_next {
                    if let Some(b) = b_ptr.as_ref() {
                        links.push((a.number.to_string(), "+".to_string(), b.number.to_string(),
                                   (if b_strand {"+"} else {"-"}).to_string()));
                    }
                }
                for &(b_ptr, b_strand) in &a.reverse_next {
                    if let Some(b) = b_ptr.as_ref() {
                        links.push((a.number.to_string(), "-".to_string(), b.number.to_string(),
                                   (if b_strand {"+"} else {"-"}).to_string()));
                    }
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

    fn find_starting_unitig(&self, seq_id: u16) -> (&Unitig, bool) {
        // For a given sequence ID, this function returns the Unitig and strand where that sequence
        // begins.
        let half_k = self.k_size / 2;
        let mut starting_unitigs = Vec::new();
        for unitig in &self.unitigs {
            for p in &unitig.forward_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == half_k {
                    starting_unitigs.push((unitig, true));
                }
            }
            for p in &unitig.reverse_positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == half_k {
                    starting_unitigs.push((unitig, false));
                }
            }
        }
        assert_eq!(starting_unitigs.len(), 1);
        starting_unitigs[0]
    }

    fn get_next_unitig(&self, seq_id: u16, unitig: &Unitig, strand: bool, pos: u32) -> Option<(&Unitig, bool, u32)> {
        let next_pos = pos + unitig.untrimmed_length(self.k_size as usize) - self.k_size as u32 + 1;
        let next_unitigs = if strand { &unitig.forward_next } else { &unitig.reverse_next };
        for &(next_unitig, next_strand) in next_unitigs {
            let u = unsafe{next_unitig.as_ref()}.unwrap();
            let positions = if next_strand { &u.forward_positions } else { &u.reverse_positions};
            for p in positions {
                if p.seq_id() == seq_id && p.strand() && p.pos == next_pos {
                    return Some((&u, next_strand, next_pos));
                }
            }
        }
        None
    }

    fn get_unitig_path_for_sequence(&self, seq: &Sequence) -> Vec<(u32, bool)> {
        let total_length = seq.length as u32;
        let half_k = self.k_size / 2;
        let mut unitig_path = Vec::new();
        let (mut unitig, mut strand) = self.find_starting_unitig(seq.id);
        let mut pos = half_k;
        loop {
            unitig_path.push((unitig.number, strand));
            match self.get_next_unitig(seq.id, unitig, strand, pos) {
                None => break,
                Some((next_unitig, next_strand, next_pos)) => {
                    (unitig, strand, pos) = (next_unitig, next_strand, next_pos);
                }
            }
        }
        unitig_path
    }
}
