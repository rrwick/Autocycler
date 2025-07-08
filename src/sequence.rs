// This file defines a struct for storing the input assembly sequences and relevant information.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::fmt;

use crate::misc::{quit_with_error, reverse_complement, up_to_first_space, after_first_space};


#[derive(Clone)]
pub struct Sequence {
    pub id: u16,
    pub forward_seq: Vec<u8>,
    pub reverse_seq: Vec<u8>,
    pub filename: String,
    pub contig_header: String,
    pub length: usize,
    pub cluster: u16,
}

impl Sequence {
    pub fn new_with_seq(id: usize, seq: String, filename: String, contig_header: String,
                        length: usize, half_k: u32) -> Sequence {
        // This constructor creates a Sequence object with the actual sequence stored. This is used
        // when creating a k-mer graph from Sequences, because the actual sequence is needed to get
        // the k-mers.
        // It pads the sequence with a half-k number of dots so the entire sequence will be present
        // when overlaps are trimmed off. Dots are used because they will be regex wildcards for
        // substituting this padding with real sequence (as much as possible).
        let mut forward_seq = seq.into_bytes();
        if !forward_seq.iter().all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
            quit_with_error(&format!("{} contains non-ACGT characters", filename));
        }

        let padding = vec![b'.'; half_k as usize];
        forward_seq.splice(0..0, padding.iter().cloned());
        forward_seq.extend(padding.iter().cloned());

        let reverse_seq = reverse_complement(&forward_seq);

        Sequence {
            id: id as u16,
            forward_seq,
            reverse_seq,
            filename,
            contig_header,
            length,
            cluster: 0,
        }
    }

    pub fn new_without_seq(id: u16, filename: String, contig_header: String, length: usize,
                           cluster: u16) -> Sequence {
        // This constructor creates a Sequence object without storing the sequence. This is used at
        // later stages in Autocycler where the sequence is stored in the UnitigGraph and so doesn't
        // need to be stored here as well.
        Sequence {
            id,
            forward_seq: vec![],
            reverse_seq: vec![],
            filename,
            contig_header,
            length,
            cluster,
        }
    }

    pub fn contig_name(&self) -> String {
        up_to_first_space(&self.contig_header)
    }

    pub fn contig_description(&self) -> String {
        after_first_space(&self.contig_header)
    }

    pub fn string_for_newick(&self) -> String {
        format!("{}__{}__{}__{}_bp", self.id, self.filename, self.contig_name(), self.length)
    }

    pub fn is_trusted(&self) -> bool {
        self.contig_header.to_lowercase().contains("autocycler_trusted")
    }

    pub fn is_ignored(&self) -> bool {
        self.contig_header.to_lowercase().contains("autocycler_ignore")
    }

    pub fn cluster_weight(&self) -> usize {
        self.contig_header.to_lowercase().split_whitespace()
            .find_map(|token| { token.strip_prefix("autocycler_cluster_weight=")
                                     .and_then(|s| s.parse::<usize>().ok()) })
            .unwrap_or(1)
    }

    pub fn consensus_weight(&self) -> usize {
        self.contig_header.to_lowercase().split_whitespace()
            .find_map(|token| { token.strip_prefix("autocycler_consensus_weight=")
                                     .and_then(|s| s.parse::<usize>().ok()) })
            .unwrap_or(1)
    }
}


impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut extras = Vec::new();
        if self.is_trusted() {
            extras.push("trusted".to_string());
        }
        if self.is_ignored() {
            extras.push("ignored".to_string());
        }
        if self.cluster_weight() != 1 {
            extras.push(format!("cluster weight = {}", self.cluster_weight()));
        }
        if self.consensus_weight() != 1 {
            extras.push(format!("consensus weight = {}", self.consensus_weight()));
        }
        if extras.is_empty() {
            write!(f, "{} {} ({} bp)", self.filename, self.contig_name(), self.length)
        } else {
            write!(f, "{} {} ({} bp) [{}]", self.filename, self.contig_name(), self.length,
                   extras.join(", "))
        }
    }
}


impl fmt::Debug for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { fmt::Display::fmt(self, f) }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_trusted() {
        let mut s = Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(),
                                           "c123".to_string(), 1, 1);
        assert!(!s.is_trusted());

        s.contig_header = "c123 other stuff".to_string();
        assert!(!s.is_trusted());

        s.contig_header = "c123 Autocycler_trusted".to_string();
        assert!(s.is_trusted());

        s.contig_header = "c123 other stuff autocycler_trusted".to_string();
        assert!(s.is_trusted());

        s.contig_header = "c123 AUTOCYCLER_TRUSTED other stuff".to_string();
        assert!(s.is_trusted());
    }

    #[test]
    fn test_cluster_weight() {
        let mut s = Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(),
                                           "c123".to_string(), 1, 1);
        assert_eq!(s.cluster_weight(), 1);

        s.contig_header = "c123 other stuff".to_string();
        assert_eq!(s.cluster_weight(), 1);

        s.contig_header = "c123 Autocycler_cluster_weight=1".to_string();
        assert_eq!(s.cluster_weight(), 1);

        s.contig_header = "c123 other stuff Autocycler_cluster_weight=2 other stuff".to_string();
        assert_eq!(s.cluster_weight(), 2);

        s.contig_header = "c123 AUTOCYCLER_CLUSTER_WEIGHT=5".to_string();
        assert_eq!(s.cluster_weight(), 5);

        s.contig_header = "c123 Autocycler_cluster_weight=0".to_string();
        assert_eq!(s.cluster_weight(), 0);

        s.contig_header = "c123 autocycler_cluster_weight=1234".to_string();
        assert_eq!(s.cluster_weight(), 1234);

        // Non-integer values result in 1
        s.contig_header = "c123 Autocycler_cluster_weight=0.1".to_string();
        assert_eq!(s.cluster_weight(), 1);

        // Non-numeric values result in 1
        s.contig_header = "c123 Autocycler_cluster_weight=abc".to_string();
        assert_eq!(s.cluster_weight(), 1);
    }

    #[test]
    fn test_consensus_weight() {
        let mut s = Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(),
                                           "c123".to_string(), 1, 1);
        assert_eq!(s.consensus_weight(), 1);

        s.contig_header = "c123 other stuff".to_string();
        assert_eq!(s.consensus_weight(), 1);

        s.contig_header = "c123 Autocycler_consensus_weight=1".to_string();
        assert_eq!(s.consensus_weight(), 1);

        s.contig_header = "c123 AUTOCYCLER_CONSENSUS_WEIGHT=2".to_string();
        assert_eq!(s.consensus_weight(), 2);

        s.contig_header = "c123 other stuff Autocycler_consensus_weight=0 other stuff".to_string();
        assert_eq!(s.consensus_weight(), 0);

        s.contig_header = "c123 autocycler_consensus_weight=1234".to_string();
        assert_eq!(s.consensus_weight(), 1234);

        // Non-integer values result in 1
        s.contig_header = "c123 Autocycler_consensus_weight=23.456".to_string();
        assert_eq!(s.consensus_weight(), 1);

        // Negative values result in 1
        s.contig_header = "c123 Autocycler_consensus_weight=-1".to_string();
        assert_eq!(s.consensus_weight(), 1);

        // Non-numeric values result in 1
        s.contig_header = "c123 Autocycler_consensus_weight=abc".to_string();
        assert_eq!(s.consensus_weight(), 1);
    }

    #[test]
    fn test_sequence_display() {
        let mut s = Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(),
                                           "c123".to_string(), 1, 1);
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp)");

        s.contig_header = "c123 other stuff".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp)");

        s.contig_header = "c123 Autocycler_trusted".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp) [trusted]");

        s.contig_header = "c123 Autocycler_ignore".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp) [ignored]");

        s.contig_header = "c123 Autocycler_cluster_weight=2".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp) [cluster weight = 2]");

        s.contig_header = "c123 Autocycler_consensus_weight=2".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp) [consensus weight = 2]");

        s.contig_header = "c123 Autocycler_trusted Autocycler_cluster_weight=2".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp) [trusted, cluster weight = 2]");

        s.contig_header = "c123 Autocycler_trusted Autocycler_consensus_weight=2".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp) [trusted, consensus weight = 2]");

        s.contig_header = "c123 Autocycler_consensus_weight=0 \
                           Autocycler_cluster_weight=0".to_string();
        assert_eq!(s.to_string(), "assembly_1.fasta c123 (1 bp) [cluster weight = 0, \
                                   consensus weight = 0]");
    }
}
