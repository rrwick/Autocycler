// This file defines structs for building a k-mer De Bruijn graph from the input assemblies.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
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

use crate::misc::{reverse_complement, strand};
use crate::position::Position;
use crate::sequence::Sequence;

pub static ALPHABET: [u8; 5] = [b'.', b'A', b'C', b'G', b'T'];


pub struct Kmer {
    // Kmer objects store a raw pointer to sequence. This is faster and uses less memory than
    // storing a copy, and it avoids a bunch of tricky lifetimes which would be needed to store a
    // slice. However, it requires unsafe code to access the k-mer sequence.
    pointer: *const u8,
    length: usize,
    pub positions: Vec<Position>,
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
        unsafe{ from_raw_parts(self.pointer, self.length) }
    }

    pub fn add_position(&mut self, seq_id: u16, strand: bool, pos: usize) {
        self.positions.push(Position::new(seq_id, strand, pos));
    }

    pub fn depth(&self) -> usize {
        // Returns how many times this k-mer appears in the input sequences.
        self.positions.len()
    }

    pub fn first_position(&self) -> bool {
        // Returns true if any of this k-mer's positions are at the start of an input sequence.
        self.positions.iter().any(|p| p.pos as usize == 0)
    }
}

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let seq = std::str::from_utf8(self.seq()).unwrap();
        let positions = self.positions.iter().map(|p| p.to_string())
                                      .collect::<Vec<String>>().join(",");
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
        // Adds a sequence to the KmerGraph. For each k-mer in the sequence, a Kmer object and its
        // reverse complement are created (if necessary), and then the position of that k-mer in
        // the sequence is added to the Kmer object.
        let k_size = self.k_size as usize;
        let half_k = (self.k_size / 2) as usize;
        let two_half_k = half_k + half_k;

        let forward_raw = seq.forward_seq.as_ptr();
        let reverse_raw = seq.reverse_seq.as_ptr();

        for forward_start in 0..seq.length {
            let forward_end = forward_start + k_size;
            let reverse_start = seq.length + two_half_k - forward_end;
            let reverse_end = reverse_start + k_size;
            let forward_k = &seq.forward_seq[forward_start..forward_end];
            let reverse_k = &seq.reverse_seq[reverse_start..reverse_end];

            match self.kmers.entry(forward_k) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().add_position(seq.id, strand::FORWARD, forward_start);
                },
                Entry::Vacant(entry) => {
                    let mut kmer = unsafe { Kmer::new(forward_raw.add(forward_start), k_size,
                                                      assembly_count) };
                    kmer.add_position(seq.id, strand::FORWARD, forward_start);
                    entry.insert(kmer);
                }
            }

            match self.kmers.entry(reverse_k) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().add_position(seq.id, strand::REVERSE, reverse_start);
                },
                Entry::Vacant(entry) => {
                    let mut kmer = unsafe { Kmer::new(reverse_raw.add(reverse_start), k_size,
                                                      assembly_count) };
                    kmer.add_position(seq.id, strand::REVERSE, reverse_start);
                    entry.insert(kmer);
                }
            }
        }
    }

    pub fn next_kmers(&self, kmer: &[u8]) -> Vec<&Kmer> {
        // Given an input k-mer, this function returns all k-mers in the graph which overlap by k-1
        // bases on the right side. For example, ACGACT -> CGACTA, CGACTG.
        let mut next_kmers = Vec::new();
        let mut next_kmer = kmer[1..].to_vec();
        next_kmer.push(b'N');
        for &base in &ALPHABET {
            *next_kmer.last_mut().unwrap() = base;
            if let Some(k) = self.kmers.get(next_kmer.as_slice()) {
                next_kmers.push(k);
            }
        }
        debug_assert!(next_kmers.len() <= 4);
        next_kmers
    }

    pub fn prev_kmers(&self, kmer: &[u8]) -> Vec<&Kmer> {
        // Given an input k-mer, this function returns all k-mers in the graph which overlap by k-1
        // bases on the left side. For example, ACGACT -> AACGAC, GACGAC.
        let mut prev_kmers = Vec::new();
        let mut prev_kmer = vec![b'N'];
        prev_kmer.extend_from_slice(&kmer[..kmer.len() - 1]);
        for &base in &ALPHABET {
            *prev_kmer.first_mut().unwrap() = base;
            if let Some(k) = self.kmers.get(prev_kmer.as_slice()) {
                prev_kmers.push(k);
            }
        }
        debug_assert!(prev_kmers.len() <= 4);
        prev_kmers
    }

    pub fn iterate_kmers(&self) -> impl Iterator<Item = &Kmer> {
        // Iterates through the Kmer objects in alphabetical order.
        let mut sorted_keys: Vec<&&[u8]> = self.kmers.keys().collect();
        sorted_keys.sort_unstable();
        sorted_keys.into_iter().map(move |&k| self.kmers.get(k).unwrap())
    }

    pub fn reverse(&self, kmer: &Kmer) -> &Kmer {
        // Given a Kmer object, this function returns the reverse-complement Kmer object. Since all
        // k-mers are added on both strands, it can be assumed that the reverse-complement Kmer
        // object exists.
        let reverse_seq: &[u8] = &reverse_complement(kmer.seq());
        self.kmers.get(reverse_seq).unwrap()
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
        k.add_position(1, strand::FORWARD, 123);
        k.add_position(2, strand::REVERSE, 456);
        assert_eq!(format!("{}", k), "ACGA:1+123,2-456");
    }

    #[test]
    fn test_kmer_graph() {
        let k_size = 5; let half_k = k_size / 2;
        let mut kmer_graph = KmerGraph::new(k_size);
        let seq = Sequence::new_with_seq(1, "ACGACTGACATCAGCACTGC".to_string(),
                                         "assembly.fasta".to_string(), "contig_1".to_string(), 20, half_k);
        kmer_graph.add_sequence(&seq, 1);
        // Graph contains these 40 5-mers:
        // ..ACG ..GCA .ACGA .GCAG ACATC ACGAC ACTGA ACTGC AGCAC AGTCG
        // AGTGC ATCAG ATGTC CACTG CAGCA CAGTC CAGTG CATCA CGACT CGT..
        // CTGAC CTGAT CTGC. GACAT GACTG GATGT GCACT GCAGT GCTGA GTCAG
        // GTCGT GTGCT TCAGC TCAGT TCGT. TGACA TGATG TGC.. TGCTG TGTCA
        assert_eq!(kmer_graph.kmers.len(), 40);
    }

    #[test]
    fn test_next_kmers() {
        let k_size = 5; let half_k = k_size / 2;
        let mut kmer_graph = KmerGraph::new(k_size);
        let seq = Sequence::new_with_seq(1, "ACGACTGACATCAGCACTGC".to_string(),
                                         "assembly.fasta".to_string(), "contig_1".to_string(), 20, half_k);
        kmer_graph.add_sequence(&seq, 1);

        let next = kmer_graph.next_kmers(b"ACATC");
        assert_eq!(next.len(), 1);
        assert_eq!(next[0].seq(), b"CATCA".as_slice());

        let next = kmer_graph.next_kmers(b"CACTG");
        assert_eq!(next.len(), 2);
        assert_eq!(next[0].seq(), b"ACTGA".as_slice());
        assert_eq!(next[1].seq(), b"ACTGC".as_slice());

        let next = kmer_graph.next_kmers(b"ACTGA");
        assert_eq!(next.len(), 2);
        assert_eq!(next[0].seq(), b"CTGAC".as_slice());
        assert_eq!(next[1].seq(), b"CTGAT".as_slice());

        let next = kmer_graph.next_kmers(b"AAAAA");
        assert_eq!(next.len(), 0);
    }

    #[test]
    fn test_prev_kmers() {
        let k_size = 5; let half_k = k_size / 2;
        let mut kmer_graph = KmerGraph::new(k_size);
        let seq = Sequence::new_with_seq(1, "ACGACTGACATCAGCACTGC".to_string(),
                                         "assembly.fasta".to_string(), "contig_1".to_string(), 20, half_k);
        kmer_graph.add_sequence(&seq, 1);

        let prev = kmer_graph.prev_kmers(b"CATCA");
        assert_eq!(prev.len(), 1);
        assert_eq!(prev[0].seq(), b"ACATC".as_slice());

        let prev = kmer_graph.prev_kmers(b"CTGAC");
        assert_eq!(prev.len(), 2);
        assert_eq!(prev[0].seq(), b"ACTGA".as_slice());
        assert_eq!(prev[1].seq(), b"GCTGA".as_slice());

        let prev = kmer_graph.prev_kmers(b"ACTGC");
        assert_eq!(prev.len(), 2);
        assert_eq!(prev[0].seq(), b"CACTG".as_slice());
        assert_eq!(prev[1].seq(), b"GACTG".as_slice());

        let prev = kmer_graph.prev_kmers(b"AAAAA");
        assert_eq!(prev.len(), 0);
    }

    #[test]
    fn test_iterate_kmers() {
        let k_size = 5; let half_k = k_size / 2;
        let mut kmer_graph = KmerGraph::new(k_size);
        let seq = Sequence::new_with_seq(1, "ACGACTGACATCAGCACTGC".to_string(),
                                         "assembly.fasta".to_string(), "contig_1".to_string(), 20, half_k);
        kmer_graph.add_sequence(&seq, 1);
        let expected_kmers = vec![
            "..ACG", "..GCA", ".ACGA", ".GCAG", "ACATC", "ACGAC", "ACTGA", "ACTGC", "AGCAC", "AGTCG",
            "AGTGC", "ATCAG", "ATGTC", "CACTG", "CAGCA", "CAGTC", "CAGTG", "CATCA", "CGACT", "CGT..",
            "CTGAC", "CTGAT", "CTGC.", "GACAT", "GACTG", "GATGT", "GCACT", "GCAGT", "GCTGA", "GTCAG",
            "GTCGT", "GTGCT", "TCAGC", "TCAGT", "TCGT.", "TGACA", "TGATG", "TGC..", "TGCTG", "TGTCA"
        ];
        let expected_kmers: Vec<&[u8]> = expected_kmers.iter().map(|s| s.as_bytes()).collect();
        let actual_kmers: Vec<&[u8]> = kmer_graph.iterate_kmers().map(|kmer| kmer.seq()).collect();
        assert_eq!(expected_kmers, actual_kmers);
    }
}
