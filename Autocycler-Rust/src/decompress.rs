// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::fs;
use std::fs::File;
use std::path::PathBuf;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_file_exists, quit_with_error};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn decompress(in_gfa: PathBuf, out_dir: PathBuf) {
    section_header("Starting autocycler decompress");
    explanation("This command will take a compacted De Bruijn graph (made by autocycler \
                 compress), reconstruct the assemblies used to build that graph and save them \
                 in the specified directory.");
    print_settings(&in_gfa, &out_dir);
    create_output_dir(&out_dir);
    let (unitig_graph, sequences) = load_graph(&in_gfa);
    unitig_graph.save_gfa(&out_dir.join("temp_test.gfa"), &sequences);  // TEMP - this graph should be identical to the input graph
    let original_seqs = unitig_graph.reconstruct_original_sequences(&sequences);
    save_original_seqs(&out_dir, original_seqs);
}


fn print_settings(in_gfa: &PathBuf, out_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_dir {}", out_dir.display());
    check_if_file_exists(&in_gfa);
}


fn create_output_dir(out_dir: &PathBuf) {
    match fs::create_dir_all(&out_dir) {
        Ok(_) => {},
        Err(e) => quit_with_error(&format!("failed to create output directory\n{:?}", e)),
    }
}


fn load_graph(in_gfa: &PathBuf) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading graph");
    explanation("The compressed sequence graph is now loaded into memory.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(&in_gfa);
    eprintln!("{} unitigs", unitig_graph.unitigs.len());
    eprintln!("{} links", unitig_graph.link_count);
    (unitig_graph, sequences)
}


fn save_original_seqs(out_dir: &PathBuf, original_seqs: HashMap<String, Vec<(String, String)>>) {
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
}