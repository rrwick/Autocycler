// This file contains the code for the autocycler decompress subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::io::{BufWriter, Write};
use std::fs::File;
use std::path::PathBuf;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_is_not_dir, check_if_file_exists, create_dir};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn decompress(in_gfa: PathBuf, out_dir: PathBuf) {
    check_settings(&in_gfa, &out_dir);
    starting_message();
    print_settings(&in_gfa, &out_dir);
    create_dir(&out_dir);
    let (unitig_graph, sequences) = load_graph(&in_gfa);
    save_original_seqs(&out_dir, unitig_graph, sequences);
}


fn check_settings(in_gfa: &PathBuf, out_dir: &PathBuf) {
    check_if_file_exists(&in_gfa);
    check_if_dir_is_not_dir(&out_dir);
}


fn starting_message() {
    section_header("Starting autocycler decompress");
    explanation("This command will take a unitig graph (made by autocycler compress), reconstruct \
                 the assemblies used to build that graph and save them in the specified \
                 directory.");
}


fn print_settings(in_gfa: &PathBuf, out_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_dir {}", out_dir.display());
}


fn load_graph(in_gfa: &PathBuf) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading graph");
    explanation("The compressed sequence graph is now loaded into memory.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(&in_gfa);
    unitig_graph.print_basic_graph_info();
    (unitig_graph, sequences)
}


pub fn save_original_seqs(out_dir: &PathBuf, unitig_graph: UnitigGraph, sequences: Vec<Sequence>) {
    section_header("Reconstructing assemblies from unitig graph");
    explanation("Each contig is reconstructed by tracing its path through the unitig graph");
    let original_seqs = unitig_graph.reconstruct_original_sequences(&sequences);
    for (filename, headers_seqs) in original_seqs {
        let file_path = out_dir.join(filename);
        let file = File::create(&file_path).unwrap();
        if file_path.extension().and_then(|s| s.to_str()) == Some("gz") {
            let writer = GzEncoder::new(file, Compression::default());
            let mut buf_writer = BufWriter::new(writer);
            for (header, seq) in headers_seqs {
                writeln!(buf_writer, ">{}", header).unwrap();
                writeln!(buf_writer, "{}", seq).unwrap();
            }
        } else {
            let mut buf_writer = BufWriter::new(file);
            for (header, seq) in headers_seqs {
                writeln!(buf_writer, ">{}", header).unwrap();
                writeln!(buf_writer, "{}", seq).unwrap();
            }
        }
    }
    eprintln!();
}
