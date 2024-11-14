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
use std::path::{Path, PathBuf};
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_is_not_dir, check_if_file_exists, create_dir, quit_with_error,
                  up_to_first_space};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn decompress(in_gfa: PathBuf, out_dir: Option<PathBuf>, out_file: Option<PathBuf>) {
    check_settings(&in_gfa, &out_dir, &out_file);
    starting_message();
    print_settings(&in_gfa, &out_dir, &out_file);
    let (unitig_graph, sequences) = load_graph(&in_gfa);
    if let Some(out_dir) = out_dir {
        create_dir(&out_dir);
        save_original_seqs_to_dir(&out_dir, &unitig_graph, &sequences);
    }
    if let Some(out_file) = out_file {
        save_original_seqs_to_file(&out_file, &unitig_graph, &sequences);
    }
}


fn check_settings(in_gfa: &Path, out_dir: &Option<PathBuf>, out_file: &Option<PathBuf>) {
    check_if_file_exists(in_gfa);
    if out_dir.is_none() && out_file.is_none() {
        quit_with_error("either --out_dir or --out_file is required")
    }
    if let Some(out_dir) = out_dir {
        check_if_dir_is_not_dir(out_dir);
    }
}


fn starting_message() {
    section_header("Starting autocycler decompress");
    explanation("This command will take a unitig graph (made by autocycler compress), reconstruct \
                 the assemblies used to build that graph and save them in the specified \
                 directory and/or file.");
}


fn print_settings(in_gfa: &Path, out_dir: &Option<PathBuf>, out_file: &Option<PathBuf>) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    if let Some(out_dir) = out_dir {
        eprintln!("  --out_dir {}", out_dir.display());
    }
    if let Some(out_file) = out_file {
        eprintln!("  --out_file {}", out_file.display());
    }
    eprintln!();
}


fn load_graph(gfa: &Path) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading graph");
    explanation("The unitig graph is now loaded into memory.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(gfa);
    unitig_graph.print_basic_graph_info();
    (unitig_graph, sequences)
}


pub fn save_original_seqs_to_dir(out_dir: &Path, unitig_graph: &UnitigGraph,
                                 sequences: &Vec<Sequence>) {
    section_header("Reconstructing assemblies from unitig graph");
    explanation("Each contig is reconstructed by tracing its path through the unitig graph, with \
                 the results saved to a directory.");
    let original_seqs = unitig_graph.reconstruct_original_sequences(sequences);
    let mut filenames: Vec<&String> = original_seqs.keys().collect();
    filenames.sort();
    for filename in filenames {
        let headers_seqs = &original_seqs[filename];
        let file_path = out_dir.join(filename.clone());
        let file = File::create(&file_path).unwrap();
        eprintln!("{}:", file_path.display());
        if file_path.extension().and_then(|s| s.to_str()) == Some("gz") {
            let buf_writer = BufWriter::new(GzEncoder::new(file, Compression::default()));
            write_sequences(buf_writer, headers_seqs);
        } else {
            let buf_writer = BufWriter::new(file);
            write_sequences(buf_writer, headers_seqs);
        }
        eprintln!();
    }
}


fn write_sequences<W: Write>(mut writer: BufWriter<W>, headers_seqs: &Vec<(String, String)>) {
    for (header, seq) in headers_seqs {
        eprintln!("  {} ({} bp)", up_to_first_space(header), seq.len());
        writeln!(writer, ">{}", header).unwrap();
        writeln!(writer, "{}", seq).unwrap();
    }
}


fn save_original_seqs_to_file(out_file: &Path, unitig_graph: &UnitigGraph,
                              sequences: &Vec<Sequence>) {
    section_header("Reconstructing assemblies from unitig graph");
    explanation("Each contig is reconstructed by tracing its path through the unitig graph, with \
                 the results saved to a file.");
    eprintln!("{}:", out_file.display());
    let original_seqs = unitig_graph.reconstruct_original_sequences(sequences);
    let file = File::create(out_file).unwrap();
    let mut buf_writer = BufWriter::new(file);
    let mut filenames: Vec<&String> = original_seqs.keys().collect();
    filenames.sort();
    for filename in filenames {
        let headers_seqs = &original_seqs[filename];
        let clean_filename = filename.replace(" ", "_");
        for (header, seq) in headers_seqs {
            eprintln!("  {}__{} ({} bp)", filename, up_to_first_space(header), seq.len());
            writeln!(buf_writer, ">{}__{}", clean_filename, header).unwrap();
            writeln!(buf_writer, "{}", seq).unwrap();
        }
    }
    eprintln!();
}
