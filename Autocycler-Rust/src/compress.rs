// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use crate::log::{section_header, explanation};
use crate::misc::{find_all_assemblies, load_fasta};
use crate::kmer_graph::KmerGraph;
use crate::sequence::Sequence;

use regex::Regex;
use std::path::PathBuf;


pub fn compress(in_dir: PathBuf, out_gfa: PathBuf, k_size: u32) {
    section_header("Starting autocycler compress");
    explanation("This command will find all assemblies in the given input directory and compress \
                 them into a compacted De Bruijn graph. This graph can then be used to recover \
                 the assemblies (with autocycler decompress) or generate a consensus assembly \
                 (with autocycler resolve).");
    print_settings(&in_dir, &out_gfa, k_size);

    let (sequences, assembly_count) = load_sequences(&in_dir, k_size);
    let mut kmer_graph = KmerGraph::new(k_size);
    eprintln!("\nAdding k-mers to graph...");
    kmer_graph.add_sequences(&sequences, assembly_count);
    eprintln!("Graph contains {} k-mers", kmer_graph.kmers.len());

    // let unitig_graph = UnitigGraph::new(&kmer_graph);
    // unitig_graph.save_gfa(&out_gfa);
}


fn load_sequences(in_dir: &PathBuf, k_size: u32) -> (Vec<Sequence>, usize) {
    let mut seq_id = 0u16;
    let whitespace_re = Regex::new(r"\s+").unwrap();
    let mut sequences = Vec::new();
    let assemblies = find_all_assemblies(in_dir);
    eprintln!("\nLoading sequences:");
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
    (sequences, assemblies.len())
}

fn print_settings(in_dir: &PathBuf, out_gfa: &PathBuf, k_size: u32) {
    eprintln!("Settings:");
    eprintln!("  --in_dir {}", in_dir.display());
    eprintln!("  --out_gfa {}", out_gfa.display());
    eprintln!("  --kmer {}", k_size);
}
