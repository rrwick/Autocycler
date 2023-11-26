// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use regex::Regex;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;

use crate::position::Position;
use crate::misc::{load_fasta, reverse_complement};


pub struct Kmer {
    seq: String,
    positions: Vec<Position>,
}

impl Kmer {
    pub fn new(seq: String) -> Kmer {
        Kmer {
            seq,
            positions: Vec::new(),
        }
    }

    pub fn add_position(&mut self, seq_id: u32, strand: i32, pos: u32) {
        let position = Position::new(seq_id, strand, pos, None, None, None);
        self.positions.push(position);
    }

    pub fn count(&self) -> usize {
        self.positions.len()
    }

    pub fn first_position(&self, half_k: u32) -> bool {
        self.positions.iter().any(|p| p.pos == half_k)
    }
}

impl<'a> fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let positions = self.positions.iter()
                                      .map(|p| p.to_string())
                                      .collect::<Vec<String>>()
                                      .join(",");
        write!(f, "{}:{}", self.seq, positions)
    }
}



pub struct KmerGraph {
    k_size: u32,
    kmers: HashMap<String, Kmer>,
    id_to_contig_info: HashMap<u32, (String, String, usize)>,
}

impl KmerGraph {
    pub fn new(k_size: u32) -> KmerGraph {
        KmerGraph {
            k_size,
            kmers: HashMap::new(),
            id_to_contig_info: HashMap::new(),
        }
    }

    pub fn add_assemblies(&mut self, assemblies: Vec<PathBuf>) -> io::Result<()> {
        let mut seq_id = 0u32;
        let whitespace_re = Regex::new(r"\s+").unwrap();

        for assembly in assemblies {
            println!("\nAdding {:?} to graph:", assembly);
            let fasta = load_fasta(&assembly);
            for (name, info, seq) in &fasta {
                if seq.len() < self.k_size as usize {
                    continue;
                }
                seq_id += 1;
                println!("  {}: {} ({} bp)", seq_id, name, seq.len());

                self.add_sequence(&seq, seq_id);

                let contig_header = name.to_string() + " " + &info;
                let contig_header = whitespace_re.replace_all(&contig_header, " ");
                let filename = assembly.file_name().unwrap().to_string_lossy().into_owned();
                self.id_to_contig_info.insert(seq_id, (filename, contig_header.into_owned(), seq.len()));
            }
        }

        println!("\nGraph contains {} k-mers", self.kmers.len());
        Ok(())
    }

    fn add_sequence(&mut self, seq: &str, seq_id: u32) {
        let k_size = self.k_size as usize;
        let half_k = (self.k_size / 2) as usize;
        let rev_seq = reverse_complement(seq);
        for forward_pos in 0..seq.len() - k_size + 1 {
            let reverse_pos = seq.len() - forward_pos - k_size;

            let forward_seq = &seq[forward_pos..forward_pos + k_size];
            let reverse_seq = &rev_seq[reverse_pos..reverse_pos + k_size];

            self.kmers.entry(forward_seq.to_string()).or_insert_with(|| Kmer::new(forward_seq.to_string()));
            self.kmers.entry(reverse_seq.to_string()).or_insert_with(|| Kmer::new(reverse_seq.to_string()));

            let forward_centre_pos = forward_pos + half_k;
            let reverse_centre_pos = reverse_pos + half_k;

            self.kmers.get_mut(forward_seq).unwrap().add_position(seq_id, 1, forward_centre_pos as u32);
            self.kmers.get_mut(reverse_seq).unwrap().add_position(seq_id, -1, reverse_centre_pos as u32);
        }
    }
}
