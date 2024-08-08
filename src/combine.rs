// This file contains the code for the autocycler combine subcommand.

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
use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_file_exists, create_dir};
use crate::unitig_graph::UnitigGraph;


pub fn combine(in_gfas: Vec<PathBuf>, out_prefix: PathBuf) {
    let combined_gfa = out_prefix.with_extension("gfa");
    let combined_fasta = out_prefix.with_extension("fasta");

    check_settings(&in_gfas);
    if let Some(parent) = combined_gfa.parent() {
        create_dir(&parent.to_path_buf());
    }
    starting_message();
    print_settings(&in_gfas, &out_prefix);

    // TODO: add an optional argument for reads, which will add depth values to the combined
    //       assembly. Find unique k-mers in the combined assembly and then count the occurrences
    //       of those k-mers in the reads.

    combine_clusters(&in_gfas, &combined_gfa, &combined_fasta);
    finished_message(&combined_gfa, &combined_fasta);
}


fn check_settings(in_gfas: &Vec<PathBuf>) {
    for gfa in in_gfas {
        check_if_file_exists(&gfa);
    }
}


fn starting_message() {
    section_header("Starting autocycler combine");
    explanation("This command combines different clusters into a single assembly file.");
}


fn print_settings(in_gfas: &Vec<PathBuf>, out_prefix: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --in_gfas");
    for gfa in in_gfas {
        eprintln!("    {}", gfa.display());
    }
    eprintln!("  --out_prefix {}", out_prefix.display());
    eprintln!();
}


fn finished_message(combined_gfa: &PathBuf, combined_fasta: &PathBuf) {
    section_header("Finished!");
    eprintln!("Combined graph: {}", combined_gfa.display());
    eprintln!("Combined fasta: {}", combined_fasta.display());
    eprintln!();
}


fn combine_clusters(in_gfas: &Vec<PathBuf>, combined_gfa: &PathBuf, combined_fasta: &PathBuf) {
    section_header("Combining clusters");
    explanation("This command combines different clusters into a single assembly file.");

    let mut gfa_file = File::create(combined_gfa).unwrap();
    let mut fasta_file = File::create(combined_fasta).unwrap();
    writeln!(gfa_file, "H\tVN:Z:1.0").unwrap();

    let mut offset = 0;
    for gfa in in_gfas {
        eprintln!("{}", gfa.display());
        let (graph, _) = UnitigGraph::from_gfa_file(&gfa);
        graph.print_basic_graph_info();
        for unitig in &graph.unitigs {
            let unitig = unitig.borrow();
            let unitig_number = unitig.number + offset;
            let unitig_seq = String::from_utf8_lossy(&unitig.forward_seq);
            let circ = if unitig.is_isolated_and_circular() { "\tcircular=true".to_string() }
                                                       else { "".to_string() };
            writeln!(gfa_file, "S\t{}\t{}", unitig_number, unitig_seq).unwrap();
            writeln!(fasta_file, ">{} length={}{}", unitig_number, unitig.length(), circ).unwrap();
            writeln!(fasta_file, "{}", unitig_seq).unwrap();
        }
        for (a, a_strand, b, b_strand) in &graph.get_links_for_gfa(offset) {
            writeln!(gfa_file, "L\t{}\t{}\t{}\t{}\t0M", a, a_strand, b, b_strand).unwrap();
        }
        offset += graph.max_unitig_number();
    }
}
