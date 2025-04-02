// This file contains the code for the autocycler gfa2fasta subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::log::{section_header, explanation};
use crate::misc::check_if_file_exists;
use crate::unitig_graph::UnitigGraph;


pub fn gfa2fasta(in_gfa: PathBuf, out_fasta: PathBuf) {
    check_if_file_exists(&in_gfa);
    starting_message();
    print_settings(&in_gfa, &out_fasta);
    let graph = load_graph(&in_gfa);
    save_graph_to_fasta(&graph, &out_fasta);
}


fn starting_message() {
    section_header("Starting autocycler gfa2fasta");
    explanation("This command loads an Autocycler graph and saves it as a FASTA file with \
                 topological information in the sequence headers.");
}


fn print_settings(in_gfa: &Path, out_fasta: &Path) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_fasta {}", out_fasta.display());
    eprintln!();
}


fn load_graph(gfa: &Path) -> UnitigGraph {
    section_header("Loading graph");
    explanation("The unitig graph is now loaded into memory.");
    let (graph, _) = UnitigGraph::from_gfa_file(gfa);
    graph.print_basic_graph_info();
    graph
}


fn save_graph_to_fasta(graph: &UnitigGraph, out_fasta: &Path) {
    section_header("Saving to FASTA");
    explanation("The unitig graph is now saved to a FASTA file.");
    let mut fasta_file = File::create(out_fasta).unwrap();
    let (mut circ_count, mut non_circ_count) = (0, 0);
    for unitig in &graph.unitigs {
        let unitig = unitig.borrow();
        let seq = String::from_utf8_lossy(&unitig.forward_seq);
        let circ = if unitig.is_isolated_and_circular() {
            circ_count += 1;
            " circular=true".to_string()
        } else {
            non_circ_count += 1;
            "".to_string()
        };
        writeln!(fasta_file, ">{} length={}{}", unitig.number, unitig.length(), circ).unwrap();
        writeln!(fasta_file, "{}", seq).unwrap();
    }
    eprintln!("{} circular sequence{}", circ_count, if circ_count == 1 { "" } else { "s" });
    eprintln!("{} non-circular sequence{}", non_circ_count,
              if non_circ_count == 1 { "" } else { "s" });
    eprintln!();
}
