// This file contains the code for the autocycler compress subcommand.

// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use regex::Regex;
use std::path::PathBuf;
use std::time::Instant;

use crate::log::{section_header, explanation};
use crate::misc::{find_all_assemblies, load_fasta, format_duration};
use crate::kmer_graph::KmerGraph;
use crate::sequence::Sequence;


pub fn compress(in_dir: PathBuf, out_gfa: PathBuf, k_size: u32) {
    let start_time = Instant::now();
    starting_message(&in_dir, &out_gfa, k_size);
    let (sequences, assembly_count) = load_sequences(&in_dir, k_size);
    let kmer_graph = build_kmer_graph(k_size, assembly_count, &sequences);
    let unitig_graph = build_unitig_graph(kmer_graph, &out_gfa);
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


fn load_sequences(in_dir: &PathBuf, k_size: u32) -> (Vec<Sequence>, usize) {
    section_header("Loading input assemblies");
    explanation("All contigs in the input assemblies are now loaded and given a unique integer \
                 ID.");
    eprintln!("Loading sequences:");
    let assemblies = find_all_assemblies(in_dir);
    let mut seq_id = 0u16;
    let whitespace_re = Regex::new(r"\s+").unwrap();
    let mut sequences = Vec::new();
    for assembly in &assemblies {
        for (name, info, seq) in load_fasta(&assembly) {
            let seq_len = seq.len();
            if seq_len < k_size as usize {
                continue;
            }
            seq_id += 1;
            eprintln!("  {:>2}: {} {} ({} bp)", seq_id, assembly.display(), name, seq_len);
            let contig_header = name.to_string() + " " + &info;
            let contig_header = whitespace_re.replace_all(&contig_header, " ").to_string();
            let filename = assembly.file_name().unwrap().to_string_lossy().into_owned();
            sequences.push(Sequence::new(seq_id, seq, filename, contig_header, seq_len));
        }
    }
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


fn build_unitig_graph(kmer_graph: KmerGraph, out_gfa: &PathBuf) {
    section_header("Building compacted unitig graph");
    explanation("All non-branching paths are now collapsed to form a compacted De Bruijn graph, \
                 a.k.a. a unitig graph.");
    // let unitig_graph = UnitigGraph::new(&kmer_graph);
    // unitig_graph.save_gfa(&out_gfa);
    eprintln!();
}


fn finished_message(start_time: Instant, out_gfa: PathBuf) {
    section_header("Finished!");
    explanation("You can now run autocycler resolve to generate a consensus assembly from the \
                 unitig graph.");
    eprintln!("Final unitig graph: {}", out_gfa.display());
    eprintln!();
    eprintln!("Time to run: {}", format_duration(start_time.elapsed()));
    eprintln!();
}
