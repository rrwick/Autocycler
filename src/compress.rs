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

use std::path::PathBuf;
use std::time::Instant;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_dir_is_not_dir, create_dir, find_all_assemblies,
                  load_fasta, format_duration, spinner, quit_with_error};
use crate::kmer_graph::KmerGraph;
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;
use crate::graph_simplification::simplify_structure;


pub fn compress(assemblies_dir: PathBuf, autocycler_dir: PathBuf, k_size: u32) {
    let start_time = Instant::now();
    check_settings(&assemblies_dir, &autocycler_dir, k_size);
    starting_message();
    print_settings(&assemblies_dir, &autocycler_dir, k_size);
    create_dir(&autocycler_dir);
    let (sequences, assembly_count) = load_sequences(&assemblies_dir, k_size);
    let kmer_graph = build_kmer_graph(k_size, assembly_count, &sequences);
    let mut unitig_graph = build_unitig_graph(kmer_graph);
    simplify_unitig_graph(&mut unitig_graph, &sequences);
    let out_gfa = autocycler_dir.join("1_input_assemblies.gfa");
    unitig_graph.save_gfa(&out_gfa, &sequences).unwrap();
    finished_message(start_time, out_gfa);
}


fn check_settings(assemblies_dir: &PathBuf, autocycler_dir: &PathBuf, k_size: u32) {
    check_if_dir_exists(&assemblies_dir);
    check_if_dir_is_not_dir(&autocycler_dir);
    if k_size < 11 {
        quit_with_error("--kmer cannot be less than 11");
    }
    if k_size > 501 {
        quit_with_error("--kmer cannot be greater than 501");
    }
}


fn starting_message() {
    section_header("Starting autocycler compress");
    explanation("This command finds all assemblies in the given input directory and compresses \
                 them into a compacted De Bruijn graph. This graph can then be used to recover \
                 the assemblies (with autocycler decompress) or generate a consensus assembly \
                 (with autocycler resolve).");
}


fn print_settings(assemblies_dir: &PathBuf, autocycler_dir: &PathBuf, k_size: u32) {
    eprintln!("Settings:");
    eprintln!("  --assemblies_dir {}", assemblies_dir.display());
    eprintln!("  --autocycler_dir {}", autocycler_dir.display());
    eprintln!("  --kmer {}", k_size);
    eprintln!();
}


pub fn load_sequences(assemblies_dir: &PathBuf, k_size: u32) -> (Vec<Sequence>, usize) {
    section_header("Loading input assemblies");
    explanation("Input assemblies are now loaded and each contig is given a unique ID.");
    let assemblies = find_all_assemblies(assemblies_dir);
    let mut seq_id = 0usize;
    let mut sequences = Vec::new();
    for assembly in &assemblies {
        for (name, header, seq) in load_fasta(&assembly) {
            let seq_len = seq.len();
            if seq_len < k_size as usize {
                continue;
            }
            seq_id += 1;
            eprintln!(" {:>3}: {} {} ({} bp)", seq_id, assembly.display(), name, seq_len);
            if seq_id > 32767 {
                quit_with_error("no more than 32767 input sequences are allowed");
            }
            let contig_header = header.split_whitespace().collect::<Vec<&str>>().join(" ");
            let filename = assembly.file_name().unwrap().to_string_lossy().into_owned();
            sequences.push(Sequence::new_with_seq(seq_id as u16, seq, filename, contig_header, seq_len));
        }
    }
    // TODO: I should make sure that all sequences have a unique string (assembly filename
    // followed by contig name), because any duplicates could cause problems later.
    eprintln!();
    (sequences, assemblies.len())
}


fn build_kmer_graph(k_size: u32, assembly_count: usize, sequences: &Vec<Sequence>) -> KmerGraph {
    section_header("Building k-mer De Bruijn graph");
    explanation("K-mers in the input sequences are now hashed to make a De Bruijn graph.");
    let mut kmer_graph = KmerGraph::new(k_size);
    let pb = spinner("adding k-mers to graph...");
    kmer_graph.add_sequences(&sequences, assembly_count);
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
    simplify_structure(unitig_graph, &sequences);
    pb.finish_and_clear();
    unitig_graph.print_basic_graph_info();
}


fn finished_message(start_time: Instant, out_gfa: PathBuf) {
    section_header("Finished!");
    explanation("You can now run autocycler cluster to group contigs based on their similarity.");
    eprintln!("Final unitig graph: {}", out_gfa.display());
    eprintln!("Time to run: {}", format_duration(start_time.elapsed()));
    eprintln!();
}
