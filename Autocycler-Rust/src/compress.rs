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
use crate::misc::{find_all_assemblies, load_fasta, format_duration, quit_with_error};
use crate::kmer_graph::KmerGraph;
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;
use crate::graph_simplification::simplify_structure;


pub fn compress(in_dir: PathBuf, out_gfa: PathBuf, k_size: u32) {
    let start_time = Instant::now();
    starting_message(&in_dir, &out_gfa, k_size);
    let (sequences, assembly_count) = load_sequences(&in_dir, k_size);
    let kmer_graph = build_kmer_graph(k_size, assembly_count, &sequences);
    let mut unitig_graph = build_unitig_graph(kmer_graph);
    simplify_unitig_graph(&mut unitig_graph, &sequences);
    unitig_graph.save_gfa(&out_gfa, &sequences).unwrap();
    finished_message(start_time, out_gfa);
}


fn starting_message(in_dir: &PathBuf, out_gfa: &PathBuf, k_size: u32) {
    section_header("Starting autocycler compress");
    explanation("This command will find all assemblies in the given input directory and compress \
                 them into a compacted De Bruijn graph. This graph can then be used to recover \
                 the assemblies (with autocycler decompress) or generate a consensus assembly \
                 (with autocycler resolve).");
    eprintln!("Settings:");
    eprintln!("  --in_dir {}", in_dir.display());
    eprintln!("  --out_gfa {}", out_gfa.display());
    eprintln!("  --kmer {}", k_size);
    eprintln!();
}


pub fn load_sequences(in_dir: &PathBuf, k_size: u32) -> (Vec<Sequence>, usize) {
    section_header("Loading input assemblies");
    explanation("All contigs in the input assemblies are now loaded and given a unique integer \
                 ID.");
    eprintln!("Loading sequences:");
    let assemblies = find_all_assemblies(in_dir);
    let mut seq_id = 0usize;
    let mut sequences = Vec::new();
    for assembly in &assemblies {
        for (name, header, seq) in load_fasta(&assembly) {
            let seq_len = seq.len();
            if seq_len < k_size as usize {
                continue;
            }
            seq_id += 1;
            eprintln!("  {:>2}: {} {} ({} bp)", seq_id, assembly.display(), name, seq_len);
            if seq_id > 32767 {
                quit_with_error("no more than 32767 input sequences are allowed");
            }
            let contig_header = header.split_whitespace().collect::<Vec<&str>>().join(" ");
            let filename = assembly.file_name().unwrap().to_string_lossy().into_owned();
            sequences.push(Sequence::new(seq_id as u16, seq, filename, contig_header, seq_len));
        }
    }
    // TODO: I should make sure that all sequences have a unique string (assembly filename
    // followed by contig name), because any duplicates could cause problems later.
    eprintln!();
    (sequences, assemblies.len())
}


fn build_kmer_graph(k_size: u32, assembly_count: usize, sequences: &Vec<Sequence>) -> KmerGraph {
    section_header("Building k-mer De Bruijn graph");
    explanation("All k-mers in the input sequences are now hashed to make a De Bruijn graph.");
    let mut kmer_graph = KmerGraph::new(k_size);
    eprintln!("Adding k-mers to graph...");
    kmer_graph.add_sequences(&sequences, assembly_count);
    eprintln!("Graph contains {} k-mers", kmer_graph.kmers.len());
    eprintln!();
    kmer_graph
}


fn build_unitig_graph(kmer_graph: KmerGraph) -> UnitigGraph {
    section_header("Building compacted unitig graph");
    explanation("All non-branching paths are now collapsed to form a compacted De Bruijn graph, \
                 a.k.a. a unitig graph.");
    let unitig_graph = UnitigGraph::from_kmer_graph(&kmer_graph);
    eprintln!("{} unitigs", unitig_graph.unitigs.len());
    eprintln!("{} links", unitig_graph.get_link_count());
    eprintln!("total length: {} bp", unitig_graph.get_total_length());
    eprintln!();
    unitig_graph
}


fn simplify_unitig_graph(unitig_graph: &mut UnitigGraph, sequences: &Vec<Sequence>) {
    section_header("Simplifying unitig graph");
    explanation("Then graph structure is now simplified by moving sequence into repeat unitigs \
                 when possible.");
    simplify_structure(unitig_graph, &sequences);
    eprintln!("{} unitigs", unitig_graph.unitigs.len());
    eprintln!("{} links", unitig_graph.get_link_count());
    eprintln!("total length: {} bp", unitig_graph.get_total_length());
    eprintln!();
}


fn finished_message(start_time: Instant, out_gfa: PathBuf) {
    section_header("Finished!");
    explanation("You can now run autocycler cluster to group contigs based on their similarity.");
    eprintln!("Final unitig graph: {}", out_gfa.display());
    eprintln!();
    eprintln!("Time to run: {}", format_duration(start_time.elapsed()));
    eprintln!();
}
