// This file contains the code for setting consentig depths using reads. It is used by the
// autocycler combine subcommand when the user supplies reads with the --reads option.

// Copyright 2026 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use fxhash::FxHashMap;
use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::unitig_graph::UnitigGraph;


pub fn set_read_depths(graphs: &[&UnitigGraph], reads: &[PathBuf], k_size: u32) -> DepthInfo {
    // Sets the depth of each tig in the given graphs using the given reads. The graphs are all
    // handled together (not one at a time) because repeats must be found across the entire
    // consensus assembly, not just within a single cluster.
    section_header("Setting depths from reads");
    explanation("The reads are now used to set the depth of each sequence in the consensus \
                 assembly, replacing the input-assembly-count depths from the cluster graphs.");

    let kmers = build_kmer_table(graphs, k_size);
    let unique_count = kmers.values().filter(|&&count| count == 1).count();
    eprintln!("Tallying {}-mers in the consensus assembly:", k_size);
    eprintln!("  {} distinct {}-mers", kmers.len(), k_size);
    eprintln!("  {} unique {}-mers", unique_count, k_size);
    eprintln!();

    // TODO: count the reads' k-mers.
    // TODO: set each tig's depth from its k-mer counts.

    let _ = reads;
    for graph in graphs {
        for unitig in &graph.unitigs {
            unitig.borrow_mut().depth = 1.0;  // TEMP
        }
    }
    DepthInfo { mean_depth: 1.0 }
}


fn build_kmer_table(graphs: &[&UnitigGraph], k_size: u32) -> FxHashMap<u64, u32> {
    // Builds a table of all of the consensus assembly's k-mers, where the value is how many times
    // that k-mer occurs in the assembly. The assembly's total length is an upper bound on the
    // number of distinct k-mers, so reserving that much space saves the table from repeatedly
    // rehashing itself as it grows.
    let capacity: usize = graphs.iter().map(|g| g.total_length() as usize).sum();
    let mut kmers = FxHashMap::with_capacity_and_hasher(capacity, Default::default());
    for graph in graphs {
        for tig in &graph.unitigs {
            let tig = tig.borrow();
            add_seq_kmers(&tig.forward_seq, k_size, tig.is_isolated_and_circular(), &mut kmers);
        }
    }
    kmers
}


fn add_seq_kmers(seq: &[u8], k_size: u32, circular: bool, kmers: &mut FxHashMap<u64, u32>) {
    // Adds each of a sequence's k-mers to the table. K-mers are stored in canonical form (the
    // lesser of the k-mer and its reverse complement) so a tig's depth doesn't depend on which
    // strand it happens to be stored in. If the sequence is circular, the k-mers which span its
    // start-end junction are included too.
    let k = k_size as usize;
    debug_assert!(k <= 31);
    if seq.len() < k { return; }
    let mask = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);
    let wrap = if circular { &seq[..k-1] } else { &seq[..0] };
    let (mut forward, mut reverse, mut valid) = (0u64, 0u64, 0usize);
    for &base in seq.iter().chain(wrap.iter()) {
        match base_to_bits(base) {
            Some(bits) => {
                forward = ((forward << 2) | bits) & mask;
                reverse = (reverse >> 2) | ((3 - bits) << shift);
                valid += 1;
            },
            // Anything which isn't an unambiguous base (e.g. an N) breaks the run of k-mers.
            None => { forward = 0; reverse = 0; valid = 0; continue; },
        }
        if valid >= k {
            *kmers.entry(forward.min(reverse)).or_insert(0) += 1;
        }
    }
}


fn base_to_bits(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}


#[derive(Debug, Default)]
pub struct DepthInfo {
    // Assembly-wide values which will go into the log and the metrics file. More will be added
    // here later, e.g. the implied read accuracy and how well the reads are explained by the
    // consensus assembly.
    pub mean_depth: f64, // length-weighted mean depth of the consensus assembly
}


#[cfg(test)]
mod tests {
    use crate::misc::reverse_complement;
    use crate::test_gfa::*;
    use super::*;

    fn kmer(seq: &str) -> u64 {
        // Independently encodes a k-mer into the canonical value used by add_seq_kmers.
        fn encode(seq: &[u8]) -> u64 {
            seq.iter().fold(0, |total, &b| (total << 2) | base_to_bits(b).unwrap())
        }
        encode(seq.as_bytes()).min(encode(&reverse_complement(seq.as_bytes())))
    }

    fn get_kmers(seq: &str, k_size: u32, circular: bool) -> FxHashMap<u64, u32> {
        let mut kmers = FxHashMap::default();
        add_seq_kmers(seq.as_bytes(), k_size, circular, &mut kmers);
        kmers
    }

    fn total_count(kmers: &FxHashMap<u64, u32>) -> u32 {
        kmers.values().sum()
    }

    #[test]
    fn test_add_seq_kmers_canonical() {
        // ACGT contains ACG and CGT, which are reverse complements of each other, so they share a
        // single canonical k-mer.
        let kmers = get_kmers("ACGT", 3, false);
        assert_eq!(kmers.len(), 1);
        assert_eq!(kmers[&kmer("ACG")], 2);
        assert_eq!(kmers[&kmer("CGT")], 2);
    }

    #[test]
    fn test_add_seq_kmers_both_strands() {
        // A sequence and its reverse complement give the same table.
        let seq = "AGCATCGACATCGACTACG";
        let rev = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
        assert_eq!(get_kmers(seq, 5, false), get_kmers(&rev, 5, false));
    }

    #[test]
    fn test_add_seq_kmers_repeat() {
        // AAA occurs at both the start and the end, so it is the only repeated k-mer.
        let kmers = get_kmers("AAACCCAAA", 3, false);
        assert_eq!(kmers.len(), 6);
        assert_eq!(total_count(&kmers), 7);
        assert_eq!(kmers[&kmer("AAA")], 2);
        assert_eq!(kmers[&kmer("CCC")], 1);
    }

    #[test]
    fn test_add_seq_kmers_circular() {
        // A linear sequence has length-k+1 k-mers, while a circular one has length k-mers, because
        // of the k-mers spanning its start-end junction.
        let seq = "AGCATCGACATCGACTACG"; // 19 bp
        assert_eq!(total_count(&get_kmers(seq, 5, false)), 15);
        assert_eq!(total_count(&get_kmers(seq, 5, true)), 19);
        assert_eq!(get_kmers(seq, 5, true)[&kmer("TACGA")], 1);  // spans the junction
    }

    #[test]
    fn test_add_seq_kmers_too_short() {
        assert!(get_kmers("ACGT", 5, false).is_empty());
        assert!(get_kmers("ACGT", 5, true).is_empty());
        assert_eq!(get_kmers("ACGTA", 5, false).len(), 1);
    }

    #[test]
    fn test_add_seq_kmers_ambiguous() {
        // The N breaks the run of k-mers, so only the ACGT on either side of it contributes.
        let kmers = get_kmers("ACGTNACGT", 3, false);
        assert_eq!(kmers.len(), 1);
        assert_eq!(kmers[&kmer("ACG")], 4);
    }

    #[test]
    fn test_build_kmer_table_one_graph() {
        // This graph has one circular 19 bp tig which contains three repeated 5-mers.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());
        let kmers = build_kmer_table(&[&graph], 5);
        assert_eq!(total_count(&kmers), 19);
        assert_eq!(kmers.len(), 16);
        assert_eq!(kmers[&kmer("ATCGA")], 2);
        assert_eq!(kmers[&kmer("TACGA")], 1);
    }

    #[test]
    fn test_build_kmer_table_two_graphs() {
        // These two graphs contain the same sequence, so its k-mers are repeats even though
        // neither graph repeats them on its own. The only k-mers left with a count of one are the
        // four which span the circular tig's junction, as they are not in the linear tig.
        let (graph_1, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());  // circular
        let (graph_2, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());  // linear
        let kmers = build_kmer_table(&[&graph_1, &graph_2], 5);
        assert_eq!(total_count(&kmers), 19 + 15);
        assert_eq!(kmers[&kmer("AGCAT")], 2);
        assert_eq!(kmers.values().filter(|&&count| count == 1).count(), 4);
        assert_eq!(kmers[&kmer("TACGA")], 1);
    }
}
