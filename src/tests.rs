// This file contains some high-level tests for Autocycler and functions common to other tests.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use rand::{rngs::StdRng, SeedableRng};
use rand::seq::SliceRandom;
use std::fs::{File, read_to_string};
use std::io::{Read, Write};
use std::path::Path;
use tempfile::tempdir;

use crate::compress::load_sequences;
use crate::decompress::save_original_seqs_to_dir;
use crate::graph_simplification::simplify_structure;
use crate::kmer_graph::KmerGraph;
use crate::metrics::InputAssemblyMetrics;
use crate::unitig_graph::UnitigGraph;


pub fn assert_almost_eq(a: f64, b: f64, epsilon: f64) {
    assert!((a - b).abs() < epsilon,
            "Numbers are not within {:?} of each other: {} vs {}", epsilon, a, b);
}


pub fn make_test_file(file_path: &Path, contents: &str) {
    let mut file = File::create(file_path).unwrap();
    write!(file, "{}", contents).unwrap();
}


pub fn make_gzipped_test_file(file_path: &Path, contents: &str) {
    let mut file = File::create(file_path).unwrap();
    let mut e = GzEncoder::new(Vec::new(), Compression::default());
    e.write_all(contents.as_bytes()).unwrap();
    let _ = file.write_all(&e.finish().unwrap());
}


fn random_seq(length: usize, seed: u64) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    let mut rng = StdRng::seed_from_u64(seed);
    (0..length).map(|_| *bases.choose(&mut rng).unwrap()).collect()
}


fn assert_same_content(a: &Path, b: &Path) {
    assert_eq!(read_to_string(a).unwrap(), read_to_string(b).unwrap());
}


fn assert_same_content_gzipped(a: &Path, b: &Path) {
    let mut gz_a = MultiGzDecoder::new(File::open(a).unwrap());
    let mut gz_b = MultiGzDecoder::new(File::open(b).unwrap());
    let mut content_a = String::new();
    let mut content_b = String::new();
    gz_a.read_to_string(&mut content_a).unwrap();
    gz_b.read_to_string(&mut content_b).unwrap();
    assert_eq!(content_a, content_b);
}


fn test_high_level(seq_a: &str, seq_b: &str, seq_c: &str, seq_d: &str, seq_e: &str,
                    k_size: u32) {
    let assembly_dir = tempdir().unwrap();
    let graph_dir = tempdir().unwrap();
    let reconstructed_dir = tempdir().unwrap();

    // Save the sequences to the assembly directory.
    let original_a = assembly_dir.path().join("a.fasta");
    let original_b = assembly_dir.path().join("b.fna");
    let original_c = assembly_dir.path().join("c.fa");
    let original_d = assembly_dir.path().join("d.fasta.gz");
    let original_e = assembly_dir.path().join("e.fna.gz");
    let original_f = assembly_dir.path().join("e.xyz");  // bad file extension, not included
    make_test_file(&original_a, seq_a);
    make_test_file(&original_b, seq_b);
    make_test_file(&original_c, seq_c);
    make_gzipped_test_file(&original_d, seq_d);
    make_gzipped_test_file(&original_e, seq_e);
    make_test_file(&original_f, seq_a);

    // Build a k-mer graph from the sequences.
    let mut metrics = InputAssemblyMetrics::default();
    let (sequences, assembly_count) = load_sequences(&assembly_dir.path().to_path_buf(),
                                                     k_size, &mut metrics, 25);
    assert_eq!(assembly_count, 5);
    let mut kmer_graph = KmerGraph::new(k_size);
    kmer_graph.add_sequences(&sequences, assembly_count);

    // Build a unitig graph and save it to file.
    let mut unitig_graph = UnitigGraph::from_kmer_graph(&kmer_graph);
    simplify_structure(&mut unitig_graph, &sequences);
    let gfa_1 = graph_dir.path().join("graph_1.gfa");
    unitig_graph.save_gfa(&gfa_1, &sequences, false).unwrap();

    // Load the unitig graph from file, save it back to file and ensure the files are the same.
    let gfa_2 = graph_dir.path().join("graph_2.gfa");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(&gfa_1);
    unitig_graph.save_gfa(&gfa_2, &sequences, false).unwrap();
    assert_same_content(&gfa_1, &gfa_2);

    // Reconstruct the sequences from the unitig graph.
    save_original_seqs_to_dir(reconstructed_dir.path(), &unitig_graph, &sequences);
    let reconstructed_a = reconstructed_dir.path().join("a.fasta");
    let reconstructed_b = reconstructed_dir.path().join("b.fna");
    let reconstructed_c = reconstructed_dir.path().join("c.fa");
    let reconstructed_d = reconstructed_dir.path().join("d.fasta.gz");
    let reconstructed_e = reconstructed_dir.path().join("e.fna.gz");

    // Make sure original sequences match reconstruction.
    assert_same_content(&original_a, &reconstructed_a);
    assert_same_content(&original_b, &reconstructed_b);
    assert_same_content(&original_c, &reconstructed_c);
    assert_same_content_gzipped(&original_d, &reconstructed_d);
    assert_same_content_gzipped(&original_e, &reconstructed_e);
}


#[test]
fn test_fixed_seqs() {
    let seq_a = ">a\nCTTATGAGCAGTCCTTAACGTAGCGGTGTGTGGCTTTGAGAA\
                     GTTAGCGGTGGCGAGCTACATCCTGGCTCCAAT\n".to_string();
    let seq_b = ">b\nACCGTTACGTTAAGGACTGCTCATAAGATTGGAGCCAGGATG\
                     TAGCTCGCCACGGCTAACTTCTCAAAGCGGCAC\n".to_string();
    let seq_c = ">c\nCATCCTGGCTCCAATCTTATGAGCAGTCCTTAACGTAACGGT\
                     GTGTGGCTTTGAGAAGTTAGCCGTGGCGAGATA\n".to_string();
    let seq_d = ">d\nGGACTGCTCATAAGATTGGAGCCAGGATGTAGCTCGCCACGG\
                     CTAACTTCTCAAAGCCACACACCGTTACGTTAA\n".to_string();
    let seq_e = ">e\nTTGAGAAGTTAGCCGTGGCGAGCTACATCCTGGCTCCAATCT\
                     TATGAGCAGTCCTTAACGTAACGGTGTGTGGCC\n".to_string();
    test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 1);
    test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 5);
    test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 9);
    test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 13);
    test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 51);
}


#[test]
fn test_random_seqs() {
    for length in [10, 20, 50, 100] {
        for seed in [0, 5, 10, 15, 20] {
            eprintln!("{}", seed);
            let seq_a = format!(">a\n{}\n", random_seq(length, seed));
            let seq_b = format!(">b\n{}\n", random_seq(length, seed+1));
            let seq_c = format!(">c\n{}\n", random_seq(length, seed+2));
            let seq_d = format!(">d\n{}\n", random_seq(length, seed+3));
            let seq_e = format!(">e\n{}\n", random_seq(length, seed+4));
            test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 3);
            test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 5);
            test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 7);
            test_high_level(&seq_a, &seq_b, &seq_c, &seq_d, &seq_e, 9);
        }
    }
}


#[test]
fn test_whitespace() {
    // This test checks that each instance of whitespace in contig headers is turned into a
    // single space. This is because tabs in contig headers will screw up the GFA file.
    let temp_dir = tempdir().unwrap();
    let temp_file = temp_dir.path().join("assembly.fasta");
    let fasta = ">name abc  def\tghi\nCTTATGAGCAGTCCTTAACGTAGCGGT\n".to_string();
    make_test_file(&temp_file, &fasta);
    let k_size = 11;
    let mut metrics = InputAssemblyMetrics::default();
    let (sequences, assembly_count) = load_sequences(&temp_dir.path().to_path_buf(), k_size,
                                                     &mut metrics, 25);
    assert_eq!(assembly_count, 1);
    let sequence = sequences.first().unwrap();
    assert_eq!(sequence.filename, "assembly.fasta");
    assert_eq!(sequence.contig_name(), "name");
    assert_eq!(sequence.contig_header, "name abc def ghi");
    assert_eq!(sequence.forward_seq,
               String::from(".....CTTATGAGCAGTCCTTAACGTAGCGGT.....").into_bytes());
}
