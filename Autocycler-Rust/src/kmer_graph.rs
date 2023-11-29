// This file defines structs for building a k-mer De Bruijn graph from the input assemblies.

// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use fxhash::FxHashMap;  // a bit faster than Rust's built-in HashMap
use std::collections::hash_map::Entry;
use std::fmt;
use std::slice::from_raw_parts;

use crate::position::KmerPos;
use crate::sequence::Sequence;


pub struct Kmer {
    // Instead of storing a slice of the sequence, Kmer objects store a raw pointer to sequence.
    // This requires unsafe code to access the k-mer sequence, but it avoids a bunch of tricky
    // lifetimes.
    pointer: *const u8,
    length: usize,
    pub positions: Vec<KmerPos>,
}

impl Kmer {
    pub fn new(pointer: *const u8, length: usize, assembly_count: usize) -> Kmer {
        Kmer {
            pointer,
            length,
            positions: Vec::with_capacity(assembly_count), // most k-mers occur once per assembly
        }
    }

    pub fn seq(&self) -> &[u8] {
        let kseq = unsafe{ from_raw_parts(self.pointer, self.length) };
        kseq
    }

    pub fn add_position(&mut self, seq_id: u16, strand: bool, pos: usize) {
        let position = KmerPos::new(seq_id, strand, pos);
        self.positions.push(position);
    }

    pub fn count(&self) -> usize {
        self.positions.len()
    }

    pub fn first_position(&self, half_k: usize) -> bool {
        self.positions.iter().any(|p| p.pos as usize == half_k)
    }
}

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let seq = std::str::from_utf8(self.seq()).unwrap();
        let positions = self.positions.iter()
                                      .map(|p| p.to_string())
                                      .collect::<Vec<String>>()
                                      .join(",");
        write!(f, "{}:{}", seq, positions)
    }
}


pub struct KmerGraph<'a> {
    pub k_size: u32,
    pub kmers: FxHashMap<&'a [u8], Kmer>,
}

impl<'a> KmerGraph<'a> {
    pub fn new(k_size: u32) -> KmerGraph<'a> {
        KmerGraph {
            k_size,
            kmers: FxHashMap::default(),
        }
    }

    pub fn add_sequences(&mut self, seqs: &'a Vec<Sequence>, assembly_count: usize) {
        for seq in seqs {
            self.add_sequence(seq, assembly_count)
        }
    }

    pub fn add_sequence(&mut self, seq: &'a Sequence, assembly_count: usize) {
        let k_size = self.k_size as usize;
        let half_k = (self.k_size / 2) as usize;

        let forward_raw = seq.forward_seq.as_ptr();
        let reverse_raw = seq.reverse_seq.as_ptr();

        for forward_start in 0..seq.length - k_size + 1 {
            let reverse_start = seq.length - forward_start - k_size;
            let forward_end = forward_start + k_size;
            let reverse_end = reverse_start + k_size;

            let forward_k = &seq.forward_seq[forward_start..forward_end];
            let reverse_k = &seq.reverse_seq[reverse_start..reverse_end];

            match self.kmers.entry(forward_k) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().add_position(seq.id, true, forward_start + half_k);
                },
                Entry::Vacant(entry) => {
                    let mut kmer = unsafe { Kmer::new(forward_raw.add(forward_start), k_size, assembly_count) };
                    kmer.add_position(seq.id, true, forward_start + half_k);
                    entry.insert(kmer);
                }
            }

            match self.kmers.entry(reverse_k) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().add_position(seq.id, false, reverse_start + half_k);
                },
                Entry::Vacant(entry) => {
                    let mut kmer = unsafe { Kmer::new(reverse_raw.add(reverse_start), k_size, assembly_count) };
                    kmer.add_position(seq.id, false, reverse_start + half_k);
                    entry.insert(kmer);
                }
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer() {
        let seq = String::from("ACGACTGACATCAGCACTGA").into_bytes();
        let raw = seq.as_ptr();
        let mut k = Kmer::new(raw, 4, 2);
        k.add_position(1, true, 123);
        k.add_position(2, false, 456);
        assert_eq!(format!("{}", k), "ACGA:1+123,2-456");
    }

    #[test]
    fn test_kmer_graph() {
        let mut kmer_graph = KmerGraph::new(4);
        let seq = Sequence::new(1, "ACGACTGACATCAGCACTGA".to_string(),
                                "assembly.fasta".to_string(), "contig_1".to_string(), 20);
        kmer_graph.add_sequence(&seq, 1);
        assert_eq!(kmer_graph.kmers.len(), 28);
    }
}
