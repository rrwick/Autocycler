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
    let (graph, _) = UnitigGraph::load_gfa_with_summary(&in_gfa);
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


fn save_graph_to_fasta(graph: &UnitigGraph, out_fasta: &Path) {
    section_header("Saving to FASTA");
    explanation("The unitig graph is now saved to a FASTA file.");
    let mut fasta_file = File::create(out_fasta).unwrap();
    let (mut circ_count, mut linear_count, mut other_count) = (0, 0, 0);
    for unitig in &graph.unitigs {
        let unitig = unitig.borrow();
        let seq = String::from_utf8_lossy(&unitig.forward_seq);
        if seq.is_empty() { continue; }
        let topology = if unitig.is_isolated_and_circular() {
            circ_count += 1;
            " circular=true topology=circular"
        } else if unitig.is_isolated_and_linear() {
            linear_count += 1;
            " circular=false topology=linear"
        } else {
            other_count += 1;
            ""
        };
        writeln!(fasta_file, ">{} length={} depth={:.1}{}\n{seq}", unitig.number, unitig.length(),
                 unitig.depth, topology).unwrap();
    }
    for (count, topology) in [(circ_count, "circular"), (linear_count, "linear"), (other_count, "other")] {
        eprintln!("{count} {topology} sequence{}", if count == 1 { "" } else { "s" });
    }
    eprintln!();
}


#[cfg(test)]
mod tests {
    use tempfile::tempdir;
    use crate::test_gfa::*;
    use super::*;

    fn assert_fasta(gfa_lines: &[String], expected: &str) {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(gfa_lines);
        save_graph_to_fasta(&graph, &fasta_file);
        assert_eq!(std::fs::read_to_string(fasta_file).unwrap(), expected);
    }

    #[test]
    fn test_gfa2fasta_empty_sequences() {
        assert_fasta(&[], "");
        let empty = "S\t1\t\tDP:f:1".to_string();
        assert_fasta(std::slice::from_ref(&empty), "");
        assert_fasta(&[empty, "S\t10\tACGT\tDP:f:1.25".to_string(),
                      "S\t2\tGATTAC\tDP:f:2".to_string()],
                     ">10 length=4 depth=1.2 circular=false topology=linear\nACGT\n\
                      >2 length=6 depth=2.0 circular=false topology=linear\nGATTAC\n");
    }

    #[test]
    fn test_gfa2fasta_1() {
        assert_fasta(&get_test_gfa_1(),
                     ">1 length=22 depth=5.0\nTTCGCTGCGCTCGCTTCGCTTT\n\
                      >2 length=18 depth=4.0\nTGCCGTCGTCGCTGTGCA\n\
                      >3 length=15 depth=1.0\nTGCCTGAATCGCCTA\n\
                      >4 length=10 depth=4.0\nGCTCGGCTCG\n\
                      >5 length=8 depth=2.0\nCGAACCAT\n\
                      >6 length=7 depth=1.0\nTACTTGT\n\
                      >7 length=5 depth=2.0\nGCCTT\n\
                      >8 length=4 depth=1.0\nATCT\n\
                      >9 length=2 depth=1.0\nGC\n\
                      >10 length=1 depth=1.0\nT\n");
    }


    #[test]
    fn test_gfa2fasta_2() {
        assert_fasta(&get_test_gfa_2(),
                     ">1 length=22 depth=1.0\nACCGCTGCGCTCGCTTCGCTCT\n\
                      >2 length=5 depth=1.0\nATGAT\n\
                      >3 length=4 depth=1.0\nGCGC\n");
    }


    #[test]
    fn test_gfa2fasta_5() {
        assert_fasta(&get_test_gfa_5(),
                     ">1 length=19 depth=1.0\nAGCATCGACATCGACTACG\n\
                      >2 length=15 depth=1.0 circular=false topology=linear\nAGCATCAGCATCAGC\n\
                      >3 length=9 depth=1.0\nGTCGCATTT\n\
                      >4 length=7 depth=1.0 circular=true topology=circular\nTCGCGAA\n\
                      >5 length=6 depth=1.0\nTTAAAC\n\
                      >6 length=4 depth=1.0\nCACA\n");
    }


    #[test]
    fn test_gfa2fasta_8() {
        assert_fasta(&get_test_gfa_8(),
                     ">1 length=19 depth=1.0 circular=true topology=circular\nAGCATCGACATCGACTACG\n");
    }

    #[test]
    fn test_gfa2fasta_9() {
        assert_fasta(&get_test_gfa_9(),
                     ">1 length=19 depth=1.0 circular=false topology=linear\nAGCATCGACATCGACTACG\n");
    }

    #[test]
    fn test_gfa2fasta_10() {
        assert_fasta(&get_test_gfa_10(),
                     ">1 length=19 depth=1.0 circular=false topology=linear\nAGCATCGACATCGACTACG\n");
    }

    #[test]
    fn test_gfa2fasta_13() {
        assert_fasta(&get_test_gfa_13(),
                     ">1 length=19 depth=1.0\nAGCATCGACATCGACTACG\n");
    }
}
