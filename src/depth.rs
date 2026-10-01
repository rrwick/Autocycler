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

use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::atomic::{AtomicU32, Ordering::Relaxed};

use crate::log::{section_header, explanation};
use crate::misc::{quit_with_error, read_batches, reverse_complement};
use crate::unitig::{Unitig, UnitigStrand};
use crate::unitig_graph::UnitigGraph;


// At 80% accuracy and k=19, only about 1.4% of a read's k-mers survive. A foreign read should
// match nothing, so this low threshold still excludes it.
const MIN_READ_HIT_RATE: f64 = 0.005;

const READ_BATCH_BASES: usize = 4 * 1024 * 1024;

// Bounds work in pathological graphs with long chains of branching, sub-k-mer-sized tigs.
const MAX_WALK_STEPS: usize = 10000;


pub fn set_read_depths(graphs: &[&UnitigGraph], reads: &[PathBuf], k_size: u32) {
    // Handle graphs together to recognise repeats across the whole consensus assembly.
    section_header("Setting depths from reads");
    explanation("The reads are now used to set the depth of each sequence in the consensus \
                 assembly, replacing the input-assembly-count depths from the cluster graphs.");

    let mut kmers = build_kmer_table(graphs, k_size);
    let repeats = find_repeats_and_reset_counts(&mut kmers);
    eprintln!("Tallying {}-mers in the consensus assembly:", k_size);
    eprintln!("  {} distinct {}-mers", kmers.len(), k_size);
    eprintln!("  {} unique {}-mers", kmers.len() - repeats.len(), k_size);
    eprintln!();

    let totals = count_read_kmers(reads, k_size, &kmers);
    eprintln!("Tallying {}-mers in the reads:", k_size);
    eprintln!("  {} reads ({} bp) used", totals.reads, totals.span_bases);
    eprintln!("  {} reads excluded", totals.rejected_reads);
    if totals.span_kmers > 0 {
        eprintln!("  {:.1}% of read {}-mers found in the assembly",
                  100.0 * totals.hits as f64 / totals.span_kmers as f64, k_size);
    }
    eprintln!();
    check_read_totals(&totals, reads);

    // Bases per hit corrects for sequencing errors and for reads having fewer k-mers than bases.
    let scale = if totals.hits > 0 { totals.span_bases as f64 / totals.hits as f64 } else { 0.0 };
    set_tig_depths(graphs, k_size, &kmers, &repeats, scale);
}


fn build_kmer_table(graphs: &[&UnitigGraph], k_size: u32) -> FxHashMap<u64, AtomicU32> {
    // Count occurrences in the assembly, including k-mers spanning links.
    let capacity: usize = graphs.iter().map(|g| g.total_length() as usize).sum();
    let mut kmers = FxHashMap::with_capacity_and_hasher(capacity, Default::default());
    for tig in graphs.iter().flat_map(|graph| &graph.unitigs) {
        add_seq_kmers(&tig.borrow().forward_seq, k_size, &mut kmers);
        for kmer in junction_kmers(tig, k_size) {
            *kmers.entry(kmer).or_default().get_mut() += 1;
        }
    }
    kmers
}


fn junction_kmers(tig: &Rc<RefCell<Unitig>>, k_size: u32) -> Vec<u64> {
    // Each k-mer crossing a link is spelled in both directions. Keep only the canonical
    // direction, but count separate walks even if they spell the same k-mer: those are repeats.
    let k = k_size as usize;
    let tig = tig.borrow();
    let mut kmers = Vec::new();
    for (seq, next) in [(&tig.forward_seq, &tig.forward_next),
                       (&tig.reverse_seq, &tig.reverse_next)] {
        let mut steps = MAX_WALK_STEPS;
        let walks = extensions(next, k - 1, &mut steps);
        for start in seq.len().saturating_sub(k - 1)..seq.len() {
            let needed = k - (seq.len() - start);
            for walk in &walks {
                if needed > walk.len() { continue; }
                if let Some((forward, reverse)) = encode_kmer(&seq[start..], &walk[..needed]) {
                    if forward < reverse { kmers.push(forward); }
                }
            }
        }
    }
    kmers
}


fn encode_kmer(prefix: &[u8], suffix: &[u8]) -> Option<(u64, u64)> {
    // Encode both strands across a junction without allocating a concatenated sequence.
    let k = prefix.len() + suffix.len();
    let mask = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);
    let (mut forward, mut reverse) = (0u64, 0u64);
    for &base in prefix.iter().chain(suffix) {
        let bits = base_to_bits(base)?;
        forward = ((forward << 2) | bits) & mask;
        reverse = (reverse >> 2) | ((3 - bits) << shift);
    }
    Some((forward, reverse))
}


fn add_seq_kmers(seq: &[u8], k_size: u32, kmers: &mut FxHashMap<u64, AtomicU32>) {
    each_kmer(seq, k_size, |_, kmer| {
        *kmers.entry(kmer).or_default().get_mut() += 1;
    });
}


fn each_kmer(seq: &[u8], k_size: u32, mut f: impl FnMut(usize, u64)) {
    // Visit each unambiguous k-mer's start position and canonical (strand-independent) value.
    let k = k_size as usize;
    debug_assert!(k <= 31);
    if seq.len() < k { return; }
    let mask = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);
    let (mut forward, mut reverse, mut valid) = (0u64, 0u64, 0usize);
    for (i, &base) in seq.iter().enumerate() {
        match base_to_bits(base) {
            Some(bits) => {
                forward = ((forward << 2) | bits) & mask;
                reverse = (reverse >> 2) | ((3 - bits) << shift);
                valid += 1;
            },
            None => { forward = 0; reverse = 0; valid = 0; continue; },
        }
        if valid >= k { f(i + 1 - k, forward.min(reverse)); }
    }
}


fn context_kmers(tig: &Rc<RefCell<Unitig>>, k_size: u32) -> BTreeMap<i32, Vec<u64>> {
    // K-mers overlapping a tig's ends, grouped by start offset relative to the tig. A read
    // passing through the tig contains exactly one of the alternatives at each offset.
    let tig = tig.borrow();
    let (k, n) = (k_size as usize, tig.forward_seq.len());
    let mut steps = MAX_WALK_STEPS;
    let mut left = extensions(&tig.reverse_next, k - 1, &mut steps);
    let mut right = extensions(&tig.forward_next, k - 1, &mut steps);

    // A walk from the tig's start spells the reverse complement of the sequence preceding it.
    for seq in &mut left { *seq = reverse_complement(seq); }

    // Tig ends with no links give a single empty context, so that the other end still gets used.
    if left.is_empty()  { left.push(Vec::new()); }
    if right.is_empty() { right.push(Vec::new()); }

    let mut kmers = BTreeMap::new();
    if n >= k {
        // Handle each end separately: no k-mer reaches across both, and the middle needs no scan.
        for seq in &left {
            add_context_kmers(&[&seq[..], &tig.forward_seq[..k-1]].concat(),
                              k_size, -(seq.len() as i32), n, &mut kmers);
        }
        for seq in &right {
            add_context_kmers(&[&tig.forward_seq[n-k+1..], &seq[..]].concat(),
                              k_size, (n - k + 1) as i32, n, &mut kmers);
        }
    } else {
        // Short tigs can need sequence from both sides, so try every context combination.
        for l in &left {
            for r in &right {
                add_context_kmers(&[&l[..], &tig.forward_seq[..], &r[..]].concat(),
                                  k_size, -(l.len() as i32), n, &mut kmers);
            }
        }
    }
    kmers
}


fn add_context_kmers(context: &[u8], k_size: u32, first_offset: i32, tig_length: usize,
                     kmers: &mut BTreeMap<i32, Vec<u64>>) {
    // Skip k-mers wholly inside the tig and deduplicate alternatives at each offset.
    let last_inside = tig_length as i32 - k_size as i32;
    each_kmer(context, k_size, |i, kmer| {
        let offset = first_offset + i as i32;
        if offset < 0 || offset > last_inside {
            let variants = kmers.entry(offset).or_default();
            if !variants.contains(&kmer) { variants.push(kmer); }
        }
    });
}


fn extensions(next: &[UnitigStrand], length: usize, steps: &mut usize) -> Vec<Vec<u8>> {
    // Preserve separate walks even when they spell the same sequence. Stop at the requested
    // length, a dead end or the shared step limit, which can leave a walk shorter than requested.
    let mut seqs = Vec::new();
    for strand in next {
        if *steps == 0 { break; }
        *steps -= 1;
        let unitig = strand.unitig();
        let unitig = unitig.borrow();
        let (seq, next) = if strand.strand { (&unitig.forward_seq, &unitig.forward_next) }
                                    else { (&unitig.reverse_seq, &unitig.reverse_next) };
        let seq = &seq[..length.min(seq.len())];
        if seq.len() == length {
            seqs.push(seq.to_vec());
            continue;
        }
        let further = extensions(next, length - seq.len(), steps);
        if further.is_empty() {
            seqs.push(seq.to_vec());
        } else {
            for suffix in further {
                seqs.push([seq, &suffix].concat());
            }
        }
    }
    seqs
}


fn set_tig_depths(graphs: &[&UnitigGraph], k_size: u32, kmers: &FxHashMap<u64, AtomicU32>,
                  repeats: &FxHashSet<u64>, scale: f64) {
    for tig in graphs.iter().flat_map(|graph| &graph.unitigs) {
        // Gather counts before borrowing mutably: a graph walk can return to this tig.
        let counts = tig_kmer_counts(tig, k_size, kmers, repeats);
        tig.borrow_mut().read_depth = clipped_mean(&counts).map(|mean| mean * scale);
    }
}


fn tig_kmer_counts(tig: &Rc<RefCell<Unitig>>, k_size: u32, kmers: &FxHashMap<u64, AtomicU32>,
                   repeats: &FxHashSet<u64>) -> Vec<u32> {
    // Exclude positions with repeated k-mers, whose counts include reads from elsewhere.
    let mut counts = Vec::new();
    each_kmer(&tig.borrow().forward_seq, k_size, |_, kmer| {
        if !repeats.contains(&kmer) {
            if let Some(count) = kmers.get(&kmer) { counts.push(count.load(Relaxed)); }
        }
    });
    // Sum alternatives at each offset: a read takes exactly one path through the tig.
    for variants in context_kmers(tig, k_size).into_values() {
        if variants.iter().any(|kmer| repeats.contains(kmer)) { continue; }
        counts.push(variants.iter().filter_map(|kmer| kmers.get(kmer))
                            .map(|count| count.load(Relaxed)).sum());
    }
    counts
}


fn clipped_mean(counts: &[u32]) -> Option<f64> {
    // Cap high outliers from collapsed repeats without discarding low counts. Six Poisson
    // standard deviations allow low-depth scatter; twice the mean allows high-depth variation.
    if counts.is_empty() { return None; }
    let count = counts.len() as f64;
    let mean = counts.iter().map(|&c| c as f64).sum::<f64>() / count;
    let limit = (mean + 6.0 * mean.sqrt()).max(2.0 * mean);
    Some(counts.iter().map(|&c| (c as f64).min(limit)).sum::<f64>() / count)
}


fn count_read_kmers(reads: &[PathBuf], k_size: u32,
                    kmers: &FxHashMap<u64, AtomicU32>) -> ReadTotals {
    read_batches(reads, READ_BATCH_BASES).map(|batch| batch.par_iter()
        .map_init(Vec::new, |hits, r| count_one_read(&r.seq, k_size, kmers, hits))
        .reduce(ReadTotals::default, |a, b| a + b))
        .fold(ReadTotals::default(), |a, b| a + b)
}


fn check_read_totals(totals: &ReadTotals, reads: &[PathBuf]) {
    if totals.reads == 0 {
        let read_str = reads.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ");
        quit_with_error(&format!("no reads were found in {} which match the consensus assembly",
                                 read_str));
    }
}


fn count_one_read<'a>(seq: &[u8], k_size: u32, kmers: &'a FxHashMap<u64, AtomicU32>,
                      hits: &mut Vec<&'a AtomicU32>) -> ReadTotals {
    // Exclude reads with too few hits before updating counts. Measure accepted reads from first
    // hit to last so non-assembly sequence (e.g. adapters) doesn't inflate depths.
    hits.clear();
    let (mut first_hit, mut last_hit, mut valid_kmers) = (0, 0, 0u64);
    each_kmer(seq, k_size, |i, kmer| {
        valid_kmers += 1;
        if let Some(count) = kmers.get(&kmer) {
            if hits.is_empty() { first_hit = i; }
            last_hit = i;
            hits.push(count);
        }
    });
    if valid_kmers == 0 { return ReadTotals::default(); }
    if (hits.len() as f64) < MIN_READ_HIT_RATE * valid_kmers as f64 {
        return ReadTotals { rejected_reads: 1, ..Default::default() };
    }
    for count in hits.iter() { count.fetch_add(1, Relaxed); }
    let span_kmers = (last_hit - first_hit + 1) as u64;
    ReadTotals { reads: 1, rejected_reads: 0, span_kmers,
                 span_bases: span_kmers + k_size as u64 - 1, hits: hits.len() as u64 }
}


fn find_repeats_and_reset_counts(kmers: &mut FxHashMap<u64, AtomicU32>) -> FxHashSet<u64> {
    let mut repeats = FxHashSet::default();
    for (kmer, count) in kmers.iter_mut() {
        let count = count.get_mut();
        if *count > 1 { repeats.insert(*kmer); }
        *count = 0;
    }
    repeats
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
struct ReadTotals {
    reads: u64,          // reads which look like they came from the consensus assembly
    rejected_reads: u64, // reads which don't
    span_bases: u64,     // bases from each included read's first k-mer hit to its last
    span_kmers: u64,     // k-mers from each included read's first k-mer hit to its last
    hits: u64,           // read k-mers found in the consensus assembly
}

impl std::ops::Add for ReadTotals {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        ReadTotals { reads:          self.reads          + other.reads,
                     rejected_reads: self.rejected_reads + other.rejected_reads,
                     span_bases:     self.span_bases     + other.span_bases,
                     span_kmers:     self.span_kmers     + other.span_kmers,
                     hits:           self.hits           + other.hits }
    }
}


#[cfg(test)]
mod tests {
    use tempfile::tempdir;
    use crate::misc::reverse_complement;
    use crate::test_gfa::*;
    use crate::tests::assert_almost_eq;
    use super::*;

    fn kmer(seq: &str) -> u64 {
        // Independently encodes a k-mer into the canonical value used by add_seq_kmers.
        fn encode(seq: &[u8]) -> u64 {
            seq.iter().fold(0, |total, &b| (total << 2) | base_to_bits(b).unwrap())
        }
        encode(seq.as_bytes()).min(encode(&reverse_complement(seq.as_bytes())))
    }

    fn to_plain(kmers: FxHashMap<u64, AtomicU32>) -> FxHashMap<u64, u32> {
        // Drops the atomics, so tests can compare counts directly.
        kmers.into_iter().map(|(kmer, count)| (kmer, count.into_inner())).collect()
    }

    fn get_kmers(seq: &str, k_size: u32) -> FxHashMap<u64, u32> {
        let mut kmers = FxHashMap::default();
        add_seq_kmers(seq.as_bytes(), k_size, &mut kmers);
        to_plain(kmers)
    }

    fn total_count(kmers: &FxHashMap<u64, u32>) -> u32 {
        kmers.values().sum()
    }

    #[test]
    fn test_add_seq_kmers_canonical() {
        // ACGT contains ACG and CGT, which are reverse complements of each other, so they share a
        // single canonical k-mer.
        let kmers = get_kmers("ACGT", 3);
        assert_eq!(kmers.len(), 1);
        assert_eq!(kmers[&kmer("ACG")], 2);
        assert_eq!(kmers[&kmer("CGT")], 2);
    }

    #[test]
    fn test_add_seq_kmers_both_strands() {
        // A sequence and its reverse complement give the same table.
        let seq = "AGCATCGACATCGACTACG";
        let rev = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
        assert_eq!(get_kmers(seq, 5), get_kmers(&rev, 5));
    }

    #[test]
    fn test_add_seq_kmers_repeat() {
        // AAA occurs at both the start and the end, so it is the only repeated k-mer.
        let kmers = get_kmers("AAACCCAAA", 3);
        assert_eq!(kmers.len(), 6);
        assert_eq!(total_count(&kmers), 7);
        assert_eq!(kmers[&kmer("AAA")], 2);
        assert_eq!(kmers[&kmer("CCC")], 1);
    }

    #[test]
    fn test_junction_kmers_circular() {
        // A circular tig links to itself, so its junction k-mers are the four which span its
        // start-end junction. Each is spelled in both directions but is counted only once.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());  // 19 bp circular tig
        let kmers = junction_kmers(&graph.unitig_index[&1], 5);
        assert_eq!(kmers.len(), 4);
        assert!(kmers.contains(&kmer("TACGA")));
    }

    #[test]
    fn test_junction_kmers_isolated() {
        // A tig with no links has no k-mers running off its ends.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        assert!(junction_kmers(&graph.unitig_index[&1], 5).is_empty());
    }

    #[test]
    fn test_build_kmer_table_shared_branch_kmers() {
        // Tigs 3 and 4 are alternative paths differing only by two extra bases on tig 3, so the
        // k-mers running off their shared tail into tig 2 are spelled by both paths. They occur
        // twice in the assembly and have to be counted twice, or each branch would be credited
        // with the other branch's reads.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_16());
        let kmers = to_plain(build_kmer_table(&[&graph], 5));
        assert_eq!(kmers[&kmer("ATGCC")], 2);  // spans into tig 2 on both paths
        assert_eq!(kmers[&kmer("ACGTT")], 1);  // spans into tig 3 only
        assert_eq!(kmers[&kmer("ACGTG")], 1);  // spans into tig 4 only
    }

    #[test]
    fn test_add_seq_kmers_too_short() {
        assert!(get_kmers("ACGT", 5).is_empty());
        assert_eq!(get_kmers("ACGTA", 5).len(), 1);
    }

    #[test]
    fn test_add_seq_kmers_ambiguous() {
        // The N breaks the run of k-mers, so only the ACGT on either side of it contributes.
        let kmers = get_kmers("ACGTNACGT", 3);
        assert_eq!(kmers.len(), 1);
        assert_eq!(kmers[&kmer("ACG")], 4);
    }

    #[test]
    fn test_build_kmer_table_one_graph() {
        // This graph has one circular 19 bp tig which contains three repeated 5-mers.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());
        let kmers = to_plain(build_kmer_table(&[&graph], 5));
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
        let kmers = to_plain(build_kmer_table(&[&graph_1, &graph_2], 5));
        assert_eq!(total_count(&kmers), 19 + 15);
        assert_eq!(kmers[&kmer("AGCAT")], 2);
        assert_eq!(kmers.values().filter(|&&count| count == 1).count(), 4);
        assert_eq!(kmers[&kmer("TACGA")], 1);
    }

    fn ready_kmer_table() -> FxHashMap<u64, AtomicU32> {
        // A table built from one circular 19 bp tig, with its counts zeroed and ready for reads.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());
        let mut kmers = build_kmer_table(&[&graph], 5);
        find_repeats_and_reset_counts(&mut kmers);
        kmers
    }

    fn count_read(seq: &str, kmers: &FxHashMap<u64, AtomicU32>) -> ReadTotals {
        count_one_read(seq.as_bytes(), 5, kmers, &mut Vec::new())
    }

    #[test]
    fn test_count_one_read_included() {
        // This read is the assembly's sequence, so all 15 of its 5-mers are hits.
        let kmers = ready_kmer_table();
        let totals = count_read("AGCATCGACATCGACTACG", &kmers);
        assert_eq!(totals.reads, 1);
        assert_eq!(totals.rejected_reads, 0);
        assert_eq!(totals.hits, 15);
        assert_eq!(totals.span_bases, 19);
        assert_eq!(totals.span_kmers, 15);
    }

    #[test]
    fn test_count_one_read_excluded() {
        // This read has nothing in common with the assembly, so it is excluded and its k-mers are
        // left uncounted.
        let kmers = ready_kmer_table();
        let totals = count_read("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", &kmers);
        assert_eq!(totals.reads, 0);
        assert_eq!(totals.rejected_reads, 1);
        assert_eq!(totals.hits, 0);
        assert_eq!(totals.span_bases, 0);
        assert!(kmers.values().all(|count| count.load(Relaxed) == 0));
    }

    #[test]
    fn test_count_one_read_too_short() {
        // A read shorter than the k-mer size is neither used nor counted as excluded.
        let kmers = ready_kmer_table();
        let totals = count_read("AGCA", &kmers);
        assert_eq!(totals.reads, 0);
        assert_eq!(totals.rejected_reads, 0);
    }

    #[test]
    fn test_count_one_read_adapter() {
        // This read has 10 bp of non-assembly sequence on its front, like an untrimmed adapter.
        // It is used, but only the part of it from the first k-mer hit to the last counts towards
        // the read length, so the adapter cannot inflate depths.
        let kmers = ready_kmer_table();
        let totals = count_read("TTTTTTTTTTAGCATCGACATCGACTACG", &kmers);
        assert_eq!(totals.reads, 1);
        assert_eq!(totals.hits, 15);
        assert_eq!(totals.span_bases, 19);
        assert_eq!(totals.span_kmers, 15);
    }

    #[test]
    fn test_count_one_read_commits_counts() {
        // ATCGA occurs twice in the read, TACGA not at all (it only spans the assembly's junction).
        let kmers = ready_kmer_table();
        count_read("AGCATCGACATCGACTACG", &kmers);
        assert_eq!(kmers[&kmer("ATCGA")].load(Relaxed), 2);
        assert_eq!(kmers[&kmer("AGCAT")].load(Relaxed), 1);
        assert_eq!(kmers[&kmer("TACGA")].load(Relaxed), 0);
    }

    #[test]
    fn test_count_read_kmers() {
        // Two of these reads are from the assembly and one isn't.
        let temp_dir = tempdir().unwrap();
        let fastq = temp_dir.path().join("reads.fastq");
        std::fs::write(&fastq, "@read1\nAGCATCGACATCGACTACG\n+\nIIIIIIIIIIIIIIIIIII\n\
                                @read2\nAGCATCGACATCGACTACG\n+\nIIIIIIIIIIIIIIIIIII\n\
                                @read3\nTTTTTTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIIIIIIII\n").unwrap();
        let kmers = ready_kmer_table();
        let totals = count_read_kmers(&[fastq], 5, &kmers);
        assert_eq!(totals.reads, 2);
        assert_eq!(totals.rejected_reads, 1);
        assert_eq!(totals.hits, 30);
        assert_eq!(totals.span_bases, 38);
        assert_eq!(kmers[&kmer("ATCGA")].load(Relaxed), 4);
    }

    #[test]
    fn test_clipped_mean() {
        assert_eq!(clipped_mean(&[]), None);
        assert_eq!(clipped_mean(&[7]), Some(7.0));
        assert_eq!(clipped_mean(&[10, 10, 10, 10]), Some(10.0));
        assert_eq!(clipped_mean(&[0, 0, 0, 0]), Some(0.0));
        assert_eq!(clipped_mean(&[8, 9, 10, 11, 12]), Some(10.0));
    }

    #[test]
    fn test_clipped_mean_clips_outlier() {
        // One wildly high count (e.g. from sequence repeated in the reads but not the assembly)
        // is pulled down to the limit, so it can't drag the depth up with it.
        let mut counts = vec![10; 100];
        assert_eq!(clipped_mean(&counts), Some(10.0));
        counts.push(10000);
        let mean: f64 = 11000.0 / 101.0;               // the mean before clipping
        let limit = 2.0 * mean;                        // higher than six standard deviations here
        assert_almost_eq(clipped_mean(&counts).unwrap(), (1000.0 + limit) / 101.0, 1e-9);
        assert!(clipped_mean(&counts).unwrap() < 13.0);
    }

    #[test]
    fn test_clipped_mean_limit_at_high_counts() {
        // At high counts, six standard deviations is a narrow margin (here 5% of the mean), so
        // the limit of twice the mean applies instead and ordinary variation survives untouched.
        let mut counts = vec![10000; 100];
        counts.extend([5000, 15000]);  // well beyond six sigma, but ordinary variation in depth
        let mean = counts.iter().sum::<u32>() as f64 / counts.len() as f64;
        assert_almost_eq(clipped_mean(&counts).unwrap(), mean, 1e-9);
    }

    #[test]
    fn test_clipped_mean_keeps_ordinary_variation() {
        // Poisson-like scatter around a depth of 10 is well within six sigma, so nothing is
        // clipped and the mean is unchanged.
        let counts = vec![4, 6, 7, 8, 9, 10, 10, 11, 12, 13, 14, 16];
        let mean = counts.iter().sum::<u32>() as f64 / counts.len() as f64;
        assert_almost_eq(clipped_mean(&counts).unwrap(), mean, 1e-9);
    }

    #[test]
    fn test_set_tig_depths() {
        // Three reads of the tig's sequence, so each of its k-mers is seen three times.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());  // linear 19 bp tig
        let mut kmers = build_kmer_table(&[&graph], 5);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        for _ in 0..3 { count_read("AGCATCGACATCGACTACG", &kmers); }
        set_tig_depths(&[&graph], 5, &kmers, &repeats, 1.0);
        assert_eq!(graph.unitigs[0].borrow().read_depth, Some(3.0));
    }

    #[test]
    fn test_set_tig_depths_scaled() {
        // The scale factor converts mean k-mer counts into read depths.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        let mut kmers = build_kmer_table(&[&graph], 5);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        count_read("AGCATCGACATCGACTACG", &kmers);
        set_tig_depths(&[&graph], 5, &kmers, &repeats, 2.5);
        assert_eq!(graph.unitigs[0].borrow().read_depth, Some(2.5));
    }

    #[test]
    fn test_set_tig_depths_no_kmers() {
        // This graph is a single 19 bp tig with no links, so with a k-mer size of 25 there is no
        // way to spell a k-mer anywhere in it and it cannot be given a depth.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        let mut kmers = build_kmer_table(&[&graph], 25);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        set_tig_depths(&[&graph], 25, &kmers, &repeats, 1.0);
        assert_eq!(graph.unitigs[0].borrow().read_depth, None);
    }

    #[test]
    fn test_set_tig_depths_short_tigs() {
        // Tigs 3 and 4 are single bases, too short to hold a k-mer of their own, so their depths
        // come from k-mers which run out through their neighbours. This read follows the path
        // through tig 3, so tig 3 gets its depth and tig 4 is shown to have no support.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_15());
        let mut kmers = build_kmer_table(&[&graph], 5);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        count_read("ATTGTAGGTACCGATCGATCGT", &kmers);
        set_tig_depths(&[&graph], 5, &kmers, &repeats, 1.0);
        assert_eq!(graph.unitig_index[&3].borrow().read_depth, Some(1.0));
        assert_eq!(graph.unitig_index[&4].borrow().read_depth, Some(0.0));
    }

    #[test]
    fn test_set_tig_depths_sums_variants() {
        // Tig 2 is preceded by either tig 3 or tig 4, so the k-mers reaching back off its start
        // have two alternatives. Two reads take one path and three take the other, and a read
        // through tig 2 contains exactly one of the alternatives, so its depth is the total of
        // the two, not the average.
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_15());
        let mut kmers = build_kmer_table(&[&graph], 5);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        for _ in 0..2 { count_read("ATTGTAGGTACCGATCGATCGT", &kmers); }
        for _ in 0..3 { count_read("CTTGTAGGTACCGATCGATCGT", &kmers); }
        set_tig_depths(&[&graph], 5, &kmers, &repeats, 1.0);
        assert_eq!(graph.unitig_index[&2].borrow().read_depth, Some(5.0));
        assert_eq!(graph.unitig_index[&3].borrow().read_depth, Some(2.0));
        assert_eq!(graph.unitig_index[&4].borrow().read_depth, Some(3.0));
    }

    fn contexts(gfa: &[String], number: u32, k_size: u32) -> Vec<(i32, Vec<u64>)> {
        let (graph, _) = UnitigGraph::from_gfa_lines(gfa);
        context_kmers(&graph.unitig_index[&number], k_size).into_iter().collect()
    }

    fn offsets(kmers: &[(i32, Vec<u64>)]) -> Vec<i32> {
        kmers.iter().map(|(offset, _)| *offset).collect()
    }

    #[test]
    fn test_extensions_step_budget() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_15());
        let tig = graph.unitig_index[&1].borrow();
        for (budget, expected) in [(0, vec![]), (1, vec!["CAA"]),
                                   (2, vec!["CAAT"]), (3, vec!["CAAT", "CAAG"])] {
            let mut steps = budget;
            let walks = extensions(&tig.forward_next, 5, &mut steps);
            let expected: Vec<_> = expected.iter().map(|s| s.as_bytes().to_vec()).collect();
            assert_eq!(walks, expected);
            assert_eq!(steps, 0);
        }
    }

    #[test]
    fn test_extensions_preserves_duplicate_walks() {
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_16());
        let tig = graph.unitig_index[&2].borrow();
        let mut steps = 3;
        assert_eq!(extensions(&tig.reverse_next, 5, &mut steps), vec![b"CATGC", b"CATGC"]);
        assert_eq!(steps, 1);
    }

    #[test]
    fn test_context_kmers_isolated() {
        // A tig with no links has no neighbouring sequence, so it gets no context k-mers.
        assert!(contexts(&get_test_gfa_9(), 1, 5).is_empty());
    }

    #[test]
    fn test_context_kmers_circular() {
        // A circular tig links to itself, so its context k-mers are the ones which span its
        // start-end junction. They are found from both ends, hence the two ranges of offsets.
        let kmers = contexts(&get_test_gfa_8(), 1, 5);  // 19 bp circular tig
        assert_eq!(offsets(&kmers), vec![-4, -3, -2, -1, 15, 16, 17, 18]);
        assert!(kmers.iter().all(|(_, variants)| variants.len() == 1));
        // TACGA spans the junction: TACG ends the tig and A starts it.
        assert!(kmers.iter().any(|(_, variants)| variants[0] == kmer("TACGA")));
    }

    #[test]
    fn test_context_kmers_hairpin() {
        // A hairpin link joins a tig to its own reverse complement, so the sequence before the
        // tig is the reverse complement of the sequence after it.
        let seq = "AGCATCGACATCGACTACG";
        let kmers = contexts(&get_test_gfa_10(), 1, 5);  // hairpin links on both ends
        assert_eq!(offsets(&kmers), vec![-4, -3, -2, -1, 15, 16, 17, 18]);
        let rev_start = String::from_utf8(reverse_complement(&seq.as_bytes()[..4])).unwrap();
        assert_eq!(kmers[0].1, vec![kmer(&format!("{}{}", rev_start, &seq[..1]))]);
    }

    #[test]
    fn test_context_kmers_dead_end() {
        // Tigs 3 and 4 are single-base dead ends, so each has exactly one k-mer: the one which
        // starts at its base and runs out through tig 2 into tig 1. That single k-mer is what
        // tells the two of them apart.
        let tig_3 = contexts(&get_test_gfa_15(), 3, 5);
        let tig_4 = contexts(&get_test_gfa_15(), 4, 5);
        assert_eq!(tig_3, vec![(0, vec![kmer("ATTGT")])]);
        assert_eq!(tig_4, vec![(0, vec![kmer("CTTGT")])]);
        assert_ne!(tig_3, tig_4);
    }

    #[test]
    fn test_context_kmers_short_tig() {
        // Tig 2 is 3 bp, shorter than the k-mer size, so all of its k-mers need context, and
        // those which reach back into tigs 3 and 4 have a variant for each.
        let kmers = contexts(&get_test_gfa_15(), 2, 5);
        assert_eq!(offsets(&kmers), vec![-1, 0, 1, 2]);
        assert_eq!(kmers[0].1, vec![kmer("ATTGT"), kmer("CTTGT")]);
        assert_eq!(kmers[1].1, vec![kmer("TTGTA")]);
        assert_eq!(kmers[2].1, vec![kmer("TGTAG")]);
        assert_eq!(kmers[3].1, vec![kmer("GTAGG")]);
    }

    #[test]
    fn test_context_kmers_long_tig() {
        // Tig 1 is longer than the k-mer size, so only the k-mers at its two ends need context:
        // its hairpin end on the left, and the fork into tigs 3 and 4 on the right.
        let kmers = contexts(&get_test_gfa_15(), 1, 5);  // 18 bp tig
        assert_eq!(offsets(&kmers), vec![-4, -3, -2, -1, 14, 15, 16, 17]);
        assert_eq!(kmers[0].1, vec![kmer("TCGTA")]);  // reverse complement of the tig's own start
        assert_eq!(kmers[4].1, vec![kmer("CCTAC")]);
        assert_eq!(kmers[7].1, vec![kmer("ACAAT"), kmer("ACAAG")]);  // the fork
    }

    #[test]
    fn test_find_repeat_kmers() {
        // AAA occurs twice in this sequence, so it is the only repeat k-mer.
        let mut kmers = FxHashMap::default();
        add_seq_kmers(b"AAACCCAAA", 3, &mut kmers);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        assert_eq!(repeats.len(), 1);
        assert!(repeats.contains(&kmer("AAA")));
        assert!(!repeats.contains(&kmer("CCC")));
    }

    #[test]
    fn test_find_repeat_kmers_zeroes_counts() {
        // Repeats stay in the table (they are needed to test whether a read k-mer is from the
        // assembly) and every count is zeroed, ready to tally read k-mers.
        let mut kmers = FxHashMap::default();
        add_seq_kmers(b"AAACCCAAA", 3, &mut kmers);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        assert_eq!(kmers.len(), 6);
        assert!(kmers.contains_key(&kmer("AAA")));
        assert!(kmers.values().all(|count| count.load(Relaxed) == 0));
        assert!(repeats.iter().all(|kmer| kmers.contains_key(kmer)));
    }

    #[test]
    fn test_find_repeat_kmers_two_graphs() {
        // The same two graphs as above: 16 distinct k-mers, of which only the four spanning the
        // circular tig's junction occur once.
        let (graph_1, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());
        let (graph_2, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        let mut kmers = build_kmer_table(&[&graph_1, &graph_2], 5);
        let repeats = find_repeats_and_reset_counts(&mut kmers);
        assert_eq!(kmers.len(), 16);
        assert_eq!(repeats.len(), 12);
        assert!(repeats.contains(&kmer("AGCAT")));
        assert!(!repeats.contains(&kmer("TACGA")));
    }
}
