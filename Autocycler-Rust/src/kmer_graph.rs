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
use crate::sequence::Sequence;

// TODO: this file currently uses String a lot instead of &str. This keeps things simple (no
// lifetimes needed) but increases the time and memory to run. It would be great if I could use
// &str more for efficiency. I might be able to do this by saving the loaded sequences into
// the KmerGraph object, so I can continue to reference into them as long as that object lives.

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
    pub k_size: u32,
    pub kmers: HashMap<String, Kmer>,
}

impl KmerGraph {
    pub fn new(k_size: u32) -> KmerGraph {
        KmerGraph {
            k_size,
            kmers: HashMap::new(),
        }
    }

    pub fn add_sequence(&mut self, seq: &Sequence) {
        let k_size = self.k_size as usize;
        let half_k = (self.k_size / 2) as usize;
        let forward_seq = &seq.seq;
        let reverse_seq = reverse_complement(&forward_seq);
        for forward_pos in 0..seq.length - k_size + 1 {
            let reverse_pos = seq.length - forward_pos - k_size;

            let forward_k = &forward_seq[forward_pos..forward_pos + k_size];
            let reverse_k = &reverse_seq[reverse_pos..reverse_pos + k_size];

            self.kmers.entry(forward_k.to_string()).or_insert_with(|| Kmer::new(forward_k.to_string()));
            self.kmers.entry(reverse_k.to_string()).or_insert_with(|| Kmer::new(reverse_k.to_string()));

            let forward_centre_pos = forward_pos + half_k;
            let reverse_centre_pos = reverse_pos + half_k;

            self.kmers.get_mut(forward_k).unwrap().add_position(seq.id, 1, forward_centre_pos as u32);
            self.kmers.get_mut(reverse_k).unwrap().add_position(seq.id, -1, reverse_centre_pos as u32);
        }
    }
}
