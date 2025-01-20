// This file contains the code for the autocycler compress subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use regex::bytes::Regex;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::str;
use std::time::Instant;

use crate::log::{section_header, explanation};
use crate::graph_simplification::simplify_structure;
use crate::kmer_graph::KmerGraph;
use crate::misc::{check_if_dir_exists, check_if_dir_is_not_dir, create_dir, find_all_assemblies,
                  load_fasta, format_duration, spinner, quit_with_error, reverse_complement};
use crate::metrics::{InputAssemblyMetrics, InputAssemblyDetails, InputContigDetails};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;



pub fn compress(assemblies_dir: PathBuf, autocycler_dir: PathBuf, k_size: u32, max_contigs: u32,
                threads: usize) {
    let start_time = Instant::now();
    check_settings(&assemblies_dir, &autocycler_dir, k_size, threads);
    starting_message();
    print_settings(&assemblies_dir, &autocycler_dir, k_size, threads);
    create_dir(&autocycler_dir);
    let mut metrics = InputAssemblyMetrics::default();
    let (sequences, assembly_count) =
        load_sequences(&assemblies_dir, k_size, &mut metrics, max_contigs);
    let kmer_graph = build_kmer_graph(k_size, assembly_count, &sequences);
    let mut unitig_graph = build_unitig_graph(kmer_graph);
    simplify_unitig_graph(&mut unitig_graph, &sequences);
    let out_gfa = autocycler_dir.join("input_assemblies.gfa");
    let out_yaml = autocycler_dir.join("input_assemblies.yaml");
    unitig_graph.save_gfa(&out_gfa, &sequences, false).unwrap();
    save_metrics(&mut metrics, assembly_count, &sequences, &unitig_graph, &out_yaml);
    finished_message(start_time, out_gfa, out_yaml);
}


fn check_settings(assemblies_dir: &Path, autocycler_dir: &Path, k_size: u32, threads: usize) {
    check_if_dir_exists(assemblies_dir);
    check_if_dir_is_not_dir(autocycler_dir);
    if k_size < 11   { quit_with_error("--kmer cannot be less than 11"); }
    if k_size > 501  { quit_with_error("--kmer cannot be greater than 501"); }
    if threads < 1   { quit_with_error("--threads cannot be less than 1"); }
    if threads > 100 { quit_with_error("--threads cannot be greater than 100"); }
    ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
}


fn starting_message() {
    section_header("Starting autocycler compress");
    explanation("This command finds all assemblies in the given input directory and compresses \
                 them into a compacted De Bruijn graph. This graph can then be used to recover \
                 the assemblies (with autocycler decompress) or generate a consensus assembly \
                 (with autocycler resolve).");
}


fn print_settings(assemblies_dir: &Path, autocycler_dir: &Path, k_size: u32, threads: usize) {
    eprintln!("Settings:");
    eprintln!("  --assemblies_dir {}", assemblies_dir.display());
    eprintln!("  --autocycler_dir {}", autocycler_dir.display());
    eprintln!("  --kmer {}", k_size);
    eprintln!("  --threads {}", threads);
    eprintln!();
}


fn check_sequence_count(sequences: &[Sequence], assembly_count: usize, max_contigs: u32) {
    let sequence_count = sequences.len() as f64;
    if sequence_count == 0.0 {
        quit_with_error("no sequences found in input assemblies")
    }
    let mean_seqs_per_assembly = sequence_count / (assembly_count as f64);
    if mean_seqs_per_assembly > max_contigs as f64 {
        let e = format!("the mean number of contigs per input assembly ({:.1}) exceeds the allowed \
                         threshold ({}). Are your input assemblies fragmented or contaminated?",
                        mean_seqs_per_assembly, max_contigs);
        quit_with_error(&e);
    }
}


pub fn load_sequences(assemblies_dir: &Path, k_size: u32, metrics: &mut InputAssemblyMetrics,
                      max_contigs: u32) -> (Vec<Sequence>, usize) {
    section_header("Loading input assemblies");
    explanation("Input assemblies are now loaded and each contig is given a unique ID.");
    let assemblies = find_all_assemblies(assemblies_dir);
    let half_k = k_size / 2;
    let mut seq_id: usize = 0;
    let mut sequences = Vec::new();
    for assembly in &assemblies {
        let mut assembly_details = InputAssemblyDetails::new(assembly);
        for (name, header, seq) in load_fasta(assembly) {
            let seq_len = seq.len();
            if seq_len < k_size as usize { continue; }
            seq_id += 1;
            eprintln!(" {:>3}: {} {} ({} bp)", seq_id, assembly.display(), name, seq_len);
            if seq_id > 32767 {
                quit_with_error("no more than 32767 input sequences are allowed");
            }
            let contig_header = header.split_whitespace().collect::<Vec<&str>>().join(" ");
            let filename = assembly.file_name().unwrap().to_string_lossy().into_owned();
            let seq = Sequence::new_with_seq(seq_id, seq, filename, contig_header, seq_len, half_k);
            assembly_details.contigs.push(InputContigDetails::new(&seq));
            sequences.push(seq);
        }
        metrics.input_assembly_details.push(assembly_details);
    }
    eprintln!();
    check_sequence_count(&sequences, assemblies.len(), max_contigs);
    let pb = spinner("repairing sequence ends...");
    sequence_end_repair(&mut sequences, k_size);
    pb.finish_and_clear();
    print_sequence_info(seq_id, assemblies.len());
    (sequences, assemblies.len())
}


fn print_sequence_info(sequence_count: usize, assembly_count: usize) {
    eprintln!("{} sequence{} loaded from {} assembl{}",
              sequence_count, match sequence_count { 1 => "", _ => "s" },
              assembly_count, match assembly_count { 1 => "y", _ => "ies" });
    eprintln!();
}


fn build_kmer_graph(k_size: u32, assembly_count: usize, sequences: &Vec<Sequence>) -> KmerGraph {
    section_header("Building k-mer De Bruijn graph");
    explanation("K-mers in the input sequences are now hashed to make a De Bruijn graph.");
    let mut kmer_graph = KmerGraph::new(k_size);
    let pb = spinner("adding k-mers to graph...");
    kmer_graph.add_sequences(sequences, assembly_count);
    pb.finish_and_clear();
    eprintln!("Graph contains {} k-mers", kmer_graph.kmers.len());
    eprintln!();
    kmer_graph
}


fn build_unitig_graph(kmer_graph: KmerGraph) -> UnitigGraph {
    section_header("Building compacted unitig graph");
    explanation("All non-branching paths are now collapsed to form a compacted De Bruijn graph, \
                 a.k.a. a unitig graph.");
    let pb = spinner("building graph...");
    let unitig_graph = UnitigGraph::from_kmer_graph(&kmer_graph);
    pb.finish_and_clear();
    unitig_graph.print_basic_graph_info();
    unitig_graph
}


fn simplify_unitig_graph(unitig_graph: &mut UnitigGraph, sequences: &Vec<Sequence>) {
    section_header("Simplifying unitig graph");
    explanation("The graph structure is now simplified by moving sequence into repeat unitigs \
                 when possible.");
    let pb = spinner("simplifying graph...");
    simplify_structure(unitig_graph, sequences);
    pb.finish_and_clear();
    unitig_graph.print_basic_graph_info();
}


fn save_metrics(metrics: &mut InputAssemblyMetrics, assembly_count: usize,
                sequences: &[Sequence], graph: &UnitigGraph, out_yaml: &Path) {
    metrics.input_assemblies_count = assembly_count as u32;
    metrics.input_assemblies_total_contigs = sequences.len() as u32;
    metrics.input_assemblies_total_length = sequences.iter().map(|s| s.length as u64).sum();
    metrics.compressed_unitig_count = graph.unitigs.len() as u32;
    metrics.compressed_unitig_total_length = graph.total_length();
    metrics.save_to_yaml(out_yaml);
}


fn finished_message(start_time: Instant, out_gfa: PathBuf, out_yaml: PathBuf) {
    section_header("Finished!");
    explanation("You can now run autocycler cluster to group contigs based on their similarity.");
    eprintln!("Compressed unitig graph: {}", out_gfa.display());
    eprintln!("Input assembly stats:    {}", out_yaml.display());
    eprintln!("Time to run: {}", format_duration(start_time.elapsed()));
    eprintln!();
}


fn sequence_end_repair(sequences: &mut Vec<Sequence>, k_size: u32) {
    // Since each sequence ends with a half-k string of dots, these will create a dead-end tip for
    // the sequence's start and end in the graph. To prevent this, this function looks for matching
    // sequences to replace the dots in other sequences, and if found, replaces the dots. Since the
    // half-k ends will be trimmed off during overlap trimming, it doesn't matter if the replacing
    // sequences are 'wrong'.
    let overlap_size = (k_size - 1) as usize;
    let all_seqs: Vec<_> = sequences.iter().flat_map(|s| vec![s.forward_seq.clone(), s.reverse_seq.clone()]).collect();
    sequences.par_iter_mut().for_each(|seq| {  // parallel for loop with rayon
        let start = &seq.forward_seq[..overlap_size];
        let start_re = Regex::new(str::from_utf8(start).unwrap()).unwrap();
        let end = &seq.forward_seq[seq.forward_seq.len() - overlap_size..];
        let end_re = Regex::new(str::from_utf8(end).unwrap()).unwrap();

        let mut all_matches = Vec::new();
        for s in &all_seqs {
            for m in start_re.find_iter(s) {
                all_matches.push(m.as_bytes().to_vec());
            }
        }
        let best_match = find_best_match(all_matches);
        seq.forward_seq.splice(..overlap_size, best_match.iter().cloned());

        let mut all_matches = Vec::new();
        for s in &all_seqs {
            for m in end_re.find_iter(s) {
                all_matches.push(m.as_bytes().to_vec());
            }
        }
        let best_match = find_best_match(all_matches);
        seq.forward_seq.splice(seq.forward_seq.len() - overlap_size.., best_match.iter().cloned());

        seq.reverse_seq = reverse_complement(&seq.forward_seq);
    });
}


fn find_best_match(matches: Vec<Vec<u8>>) -> Vec<u8> {
    // This function takes all of the regex matches and returns the best as defined by:
    // 1. fewest dots
    // 2. most occurrences
    // 3. first alphabetically
    let mut match_counts = HashMap::new();
    for m in &matches {
        let entry = match_counts.entry(m.clone()).or_insert((0, 0));
        entry.0 += 1;
        entry.1 = m.iter().filter(|&&c| c == b'.').count();
    }
    matches.into_iter().min_by(|a, b| {
        // Compare by number of `.` characters (fewer is better)
        let dot_count_a = match_counts.get(a).unwrap().1;
        let dot_count_b = match_counts.get(b).unwrap().1;
        match dot_count_a.cmp(&dot_count_b) {
            std::cmp::Ordering::Equal => {
                // Compare by frequency (higher is better)
                let freq_a = match_counts.get(a).unwrap().0;
                let freq_b = match_counts.get(b).unwrap().0;
                match freq_a.cmp(&freq_b).reverse() {
                    std::cmp::Ordering::Equal => {
                        // Compare alphabetically
                        a.cmp(b)
                    }
                    other => other,
                }
            }
            other => other,
        }
    }).expect("There should be at least one match")
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;
    use tempfile::tempdir;

    use crate::tests::make_test_file;

    #[test]
    fn test_find_best_match_1() {
        let all_matches = vec![b"...ACGT".to_vec()];
        assert_eq!(find_best_match(all_matches), b"...ACGT");

        let all_matches = vec![b"...ACGT".to_vec(), b"..GACGT".to_vec()];
        assert_eq!(find_best_match(all_matches), b"..GACGT");

        let all_matches = vec![b"..GACGT".to_vec(), b"...ACGT".to_vec()];
        assert_eq!(find_best_match(all_matches), b"..GACGT");

        let all_matches = vec![b"...GAAA".to_vec(), b"...CAAA".to_vec(), b"...TAAA".to_vec()];
        assert_eq!(find_best_match(all_matches), b"...CAAA");

        let all_matches = vec![b"...ACGT".to_vec(),
                               b"..GACGT".to_vec(),
                               b"..CACGT".to_vec(),
                               b"..GACGT".to_vec(),
                               b"..CACGT".to_vec()];
        assert_eq!(find_best_match(all_matches), b"..CACGT");

        let all_matches = vec![b"...ACGT".to_vec(),
                               b"..GACGT".to_vec(),
                               b"..GACGT".to_vec(),
                               b".AGACGT".to_vec(),
                               b".CGACGT".to_vec()];
        assert_eq!(find_best_match(all_matches), b".AGACGT");

        let all_matches = vec![b"...ACGT".to_vec(),
                               b".CGACGT".to_vec(),
                               b"..GACGT".to_vec(),
                               b".AGACGT".to_vec(),
                               b".CGACGT".to_vec()];
        assert_eq!(find_best_match(all_matches), b".CGACGT");
    }

    #[test]
    fn test_find_best_match_2() {
        let all_matches = vec![b"ACGT...".to_vec()];
        assert_eq!(find_best_match(all_matches), b"ACGT...");

        let all_matches = vec![b"ACGT...".to_vec(), b"ACGTT..".to_vec()];
        assert_eq!(find_best_match(all_matches), b"ACGTT..");

        let all_matches = vec![b"..GACGT".to_vec(), b"...ACGT".to_vec()];
        assert_eq!(find_best_match(all_matches), b"..GACGT");

        let all_matches = vec![b"GAAA...".to_vec(), b"CAAA...".to_vec(), b"TAAA...".to_vec()];
        assert_eq!(find_best_match(all_matches), b"CAAA...");

        let all_matches = vec![b"CACG...".to_vec(),
                               b"GACGT..".to_vec(),
                               b"CACGT..".to_vec(),
                               b"GACGT..".to_vec(),
                               b"CACGT..".to_vec()];
        assert_eq!(find_best_match(all_matches), b"CACGT..");

        let all_matches = vec![b"AGAC...".to_vec(),
                               b"AGACG..".to_vec(),
                               b"AGACG..".to_vec(),
                               b"AGACGT.".to_vec(),
                               b"CGACGT.".to_vec()];
        assert_eq!(find_best_match(all_matches), b"AGACGT.");
    }

    #[test]
    fn test_load_sequences_1() {
        let assembly_dir = tempdir().unwrap();
        make_test_file(&assembly_dir.path().join("a.fasta"), ">a1\nACGT\n");
        make_test_file(&assembly_dir.path().join("b.fasta"), ">b1\nACGT\n>b2\nACGT\n");
        make_test_file(&assembly_dir.path().join("c.fasta"), ">c1\nACGT\n>c2\nACGT\n>c3\nACGT\n");
        let mut metrics = InputAssemblyMetrics::default();
        let (sequences, count) = load_sequences(&assembly_dir.into_path(), 3, &mut metrics, 25);
        assert_eq!(sequences.len(), 6);
        assert_eq!(count, 3);
    }

    #[test]
    fn test_load_sequences_2() {
        // In this test, c.fasta has a duplicate sequence name which causes an error.
        let assembly_dir = tempdir().unwrap();
        make_test_file(&assembly_dir.path().join("a.fasta"), ">a1\nACGT\n");
        make_test_file(&assembly_dir.path().join("b.fasta"), ">b1\nACGT\n>b2\nACGT\n");
        make_test_file(&assembly_dir.path().join("c.fasta"), ">c1\nACGT\n>c1\nACGT\n>c3\nACGT\n");
        assert!(panic::catch_unwind(|| {
            let mut metrics = InputAssemblyMetrics::default();
            load_sequences(&assembly_dir.into_path(), 3, &mut metrics, 25);
        }).is_err());
    }
}
