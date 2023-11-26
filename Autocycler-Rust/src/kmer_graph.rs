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
use crate::misc::load_fasta;


pub struct Kmer<'a> {
    seq: &'a str,
    positions: Vec<Position>,
}

impl<'a> Kmer<'a> {
    pub fn new(seq: &'a str) -> Kmer<'a> {
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

impl<'a> fmt::Display for Kmer<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let positions = self.positions.iter()
                                      .map(|p| p.to_string())
                                      .collect::<Vec<String>>()
                                      .join(",");
        write!(f, "{}:{}", self.seq, positions)
    }
}



pub struct KmerGraph<'a> {
    k_size: u32,
    kmers: HashMap<&'a str, Kmer<'a>>,
    id_to_contig_info: HashMap<u32, (String, String, u32)>,
}

impl<'a> KmerGraph<'a> {
    pub fn new(k_size: u32) -> KmerGraph<'a> {
        KmerGraph {
            k_size,
            kmers: HashMap::new(),
            id_to_contig_info: HashMap::new(),
        }
    }

    pub fn add_assemblies(&mut self, assemblies: Vec<PathBuf>) -> io::Result<()> {
        let mut seq_id = 0u32;
        let whitespace_re = Regex::new(r"\s+").unwrap();

        for assembly_filename in assemblies {
            println!("\nAdding {:?} to graph:", assembly_filename);
            let fasta = load_fasta(&assembly_filename);
            for (name, seq) in &fasta {
                if seq.len() < self.k_size as usize {
                    continue;
                }
                seq_id += 1;
                println!("  {}: {} ({} bp)...", seq_id, name, seq.len());

                // TODO
                // TODO
                // TODO
                // TODO
                // TODO
                // TODO
                // TODO
            }
        }

        println!("\nGraph contains {} k-mers", self.kmers.len());
        Ok(())
    }
}
