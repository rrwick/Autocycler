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
use std::fs::File;
use std::io::{self, Write, BufRead, BufReader};
use std::path::PathBuf;

use crate::kmer_graph::KmerGraph;
use crate::sequence::Sequence;
use crate::unitig::Unitig;


pub struct UnitigGraph {
    pub unitigs: Vec<Unitig>,
    k_size: u32,
}

impl UnitigGraph {
    pub fn from_kmer_graph(k_graph: &KmerGraph) -> Self {
        let mut u_graph = UnitigGraph {
            unitigs: Vec::new(),
            k_size: k_graph.k_size,
        };
        u_graph.build_unitigs_from_kmer_graph(k_graph);
        u_graph.create_links();
        u_graph.connect_positions();
        u_graph.trim_overlaps();
        u_graph.renumber_unitigs();
        u_graph
    }

    pub fn from_gfa_file(gfa_filename: &PathBuf) -> Self {
        let mut u_graph = UnitigGraph {
            unitigs: Vec::new(),
            k_size: 0,
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
        u_graph.connect_positions();
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
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
    }

    fn create_links(&mut self) {
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
    }

    fn connect_positions(&mut self) {
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
    }

    fn trim_overlaps(&mut self) {
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
    }

    fn renumber_unitigs(&mut self) {
        self.unitigs.sort_by(|a, b| {
            let a_key = (1.0 / (a.length() as f64), &a.forward_seq, 1.0 / a.depth);
            let b_key = (1.0 / (b.length() as f64), &b.forward_seq, 1.0 / b.depth);
            a_key.partial_cmp(&b_key).unwrap_or(std::cmp::Ordering::Equal)
        });
        for (number, unitig) in self.unitigs.iter_mut().enumerate() {
            unitig.number = (number + 1) as u32;
        }
    }

    pub fn save_gfa(&self, gfa_filename: &PathBuf, sequences: &Vec<Sequence>) -> io::Result<()> {
        let mut file = File::create(gfa_filename)?;
        writeln!(file, "H\tVN:Z:1.0\tKM:i:{}", self.k_size)?;
        for unitig in &self.unitigs {
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
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
        links
    }

    fn get_gfa_path_line(&self, seq: &Sequence) -> String {
        let path_str = String::new();  // TEMP
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
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
}
