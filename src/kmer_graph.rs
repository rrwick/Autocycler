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

use fxhash::FxHashMap;
use std::fmt;
use std::slice::from_raw_parts;

use crate::misc::{reverse_complement, strand};
use crate::position::Position;
use crate::sequence::Sequence;

pub static ALPHABET: [u8; 5] = *b".ACGT";


pub struct Kmer {
    // Borrows sequence storage without copying. The source sequence must outlive this k-mer
    // and must not be modified while the pointer is in use.
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
        unsafe { from_raw_parts(self.pointer, self.length) }
    }

    pub fn add_position(&mut self, seq_id: u16, strand: bool, pos: usize) {
        self.positions.push(Position::new(seq_id, strand, pos));
    }

    pub fn depth(&self) -> usize {
        self.positions.len()
    }

    pub fn first_position(&self) -> bool {
        // Returns true if any of this k-mer's positions are at the start of an input sequence.
        self.positions.iter().any(|p| p.pos == 0)
    }
}

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let seq = std::str::from_utf8(self.seq()).unwrap();
        let positions = self.positions.iter().map(|p| p.to_string())
                                      .collect::<Vec<String>>().join(",");
        write!(f, "{seq}:{positions}")
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

    pub fn add_sequences(&mut self, seqs: &'a [Sequence], assembly_count: usize) {
        for seq in seqs {
            self.add_sequence(seq, assembly_count)
        }
    }

    pub fn add_sequence(&mut self, seq: &'a Sequence, assembly_count: usize) {
        let k_size = self.k_size as usize;
        let padding = 2 * (k_size / 2);
        for forward_start in 0..seq.length {
            let reverse_start = seq.length + padding - (forward_start + k_size);
            for (bases, strand, start) in [(&seq.forward_seq, strand::FORWARD, forward_start),
                                          (&seq.reverse_seq, strand::REVERSE, reverse_start)] {
                let kmer = &bases[start..start + k_size];
                self.kmers.entry(kmer)
                    .or_insert_with(|| Kmer::new(kmer.as_ptr(), k_size, assembly_count))
                    .add_position(seq.id, strand, start);
            }
        }
    }

    pub fn next_kmers(&self, kmer: &[u8]) -> Vec<&Kmer> {
        // Given an input k-mer, this function returns all k-mers in the graph which overlap by k-1
        // bases on the right side. For example, ACGACT -> CGACTA, CGACTG.
        let mut next_kmer = kmer.to_vec();
        next_kmer.rotate_left(1);
        self.matching_kmers(next_kmer, kmer.len() - 1)
    }

    pub fn prev_kmers(&self, kmer: &[u8]) -> Vec<&Kmer> {
        // Given an input k-mer, this function returns all k-mers in the graph which overlap by k-1
        // bases on the left side. For example, ACGACT -> AACGAC, GACGAC.
        let mut prev_kmer = kmer.to_vec();
        prev_kmer.rotate_right(1);
        self.matching_kmers(prev_kmer, 0)
    }

    fn matching_kmers(&self, mut sequence: Vec<u8>, variable_base: usize) -> Vec<&Kmer> {
        let mut kmers = Vec::new();
        for base in ALPHABET {
            sequence[variable_base] = base;
            if let Some(kmer) = self.kmers.get(sequence.as_slice()) {
                kmers.push(kmer);
            }
        }
        debug_assert!(kmers.len() <= 4);
        kmers
    }

    pub fn iterate_kmers(&self) -> impl Iterator<Item = &Kmer> {
        // Iterates through the Kmer objects in alphabetical order.
        let mut kmers: Vec<_> = self.kmers.values().collect();
        kmers.sort_unstable_by_key(|kmer| kmer.seq());
        kmers.into_iter()
    }

    pub fn reverse(&self, kmer: &Kmer) -> &Kmer {
        // Every k-mer is added on both strands, so its reverse complement must exist.
        self.kmers.get(reverse_complement(kmer.seq()).as_slice()).unwrap()
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
        assert_eq!(format!("{k}"), "ACGA:1+123,2-456");
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
    fn test_kmer_position_order() {
        let sequences: Vec<_> = (1..=2).map(|id| Sequence::new_with_seq(
            id, "ATATAT".to_string(), "assembly.fasta".to_string(),
            format!("contig_{id}"), 6, 1)).collect();
        let mut graph = KmerGraph::new(3);
        graph.add_sequences(&sequences, 1);
        let kmer = &graph.kmers[b"ATA".as_slice()];
        assert_eq!(kmer.to_string(), "ATA:1+1,1-3,1+3,1-1,2+1,2-3,2+3,2-1");
        assert_eq!(kmer.depth(), 8);
        assert!(!kmer.first_position());
        assert_eq!(graph.reverse(kmer).seq(), b"TAT");
        assert!(graph.kmers[b".AT".as_slice()].first_position());
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

        assert_eq!(kmer_graph.next_kmers(b"CTGC.")[0].seq(), b"TGC..");
        assert!(kmer_graph.next_kmers(b"TGC..").is_empty());
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

        assert_eq!(kmer_graph.prev_kmers(b".ACGA")[0].seq(), b"..ACG");
        assert!(kmer_graph.prev_kmers(b"..ACG").is_empty());
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
