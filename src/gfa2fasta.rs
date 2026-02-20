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
    let (mut circ_count, mut linear_count, mut other_count) = (0, 0, 0);
    for unitig in &graph.unitigs {
        let unitig = unitig.borrow();
        let seq = String::from_utf8_lossy(&unitig.forward_seq);
        if seq.is_empty() { continue; }
        let topology = if unitig.is_isolated_and_circular() {
            circ_count += 1;
            " circular=true topology=circular".to_string()
        } else if unitig.is_isolated_and_linear() {
            linear_count += 1;
            " circular=false topology=linear".to_string()
        } else {
            other_count += 1;
            "".to_string()
        };
        writeln!(fasta_file, ">{} length={}{}", unitig.number, unitig.length(), topology).unwrap();
        writeln!(fasta_file, "{seq}").unwrap();
    }
    eprintln!("{} circular sequence{}", circ_count, if circ_count == 1 { "" } else { "s" });
    eprintln!("{} linear sequence{}", linear_count, if linear_count == 1 { "" } else { "s" });
    eprintln!("{} other sequence{}", other_count, if other_count == 1 { "" } else { "s" });
    eprintln!();
}


#[cfg(test)]
mod tests {
    use tempfile::tempdir;
    use crate::test_gfa::*;
    use super::*;

    #[test]
    fn test_gfa2fasta_1() {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_1());
        save_graph_to_fasta(&graph, &fasta_file);
        let contents = std::fs::read_to_string(&fasta_file).unwrap();
        assert_eq!(contents, ">1 length=22\nTTCGCTGCGCTCGCTTCGCTTT\n\
                              >2 length=18\nTGCCGTCGTCGCTGTGCA\n\
                              >3 length=15\nTGCCTGAATCGCCTA\n\
                              >4 length=10\nGCTCGGCTCG\n\
                              >5 length=8\nCGAACCAT\n\
                              >6 length=7\nTACTTGT\n\
                              >7 length=5\nGCCTT\n\
                              >8 length=4\nATCT\n\
                              >9 length=2\nGC\n\
                              >10 length=1\nT\n");
    }


    #[test]
    fn test_gfa2fasta_2() {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_2());
        save_graph_to_fasta(&graph, &fasta_file);
        let contents = std::fs::read_to_string(&fasta_file).unwrap();
        assert_eq!(contents, ">1 length=22\nACCGCTGCGCTCGCTTCGCTCT\n\
                              >2 length=5\nATGAT\n\
                              >3 length=4\nGCGC\n");
    }


    #[test]
    fn test_gfa2fasta_5() {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_5());
        save_graph_to_fasta(&graph, &fasta_file);
        let contents = std::fs::read_to_string(&fasta_file).unwrap();
        assert_eq!(contents, ">1 length=19\nAGCATCGACATCGACTACG\n\
                              >2 length=15 circular=false topology=linear\nAGCATCAGCATCAGC\n\
                              >3 length=9\nGTCGCATTT\n\
                              >4 length=7 circular=true topology=circular\nTCGCGAA\n\
                              >5 length=6\nTTAAAC\n\
                              >6 length=4\nCACA\n");
    }


    #[test]
    fn test_gfa2fasta_8() {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_8());
        save_graph_to_fasta(&graph, &fasta_file);
        let contents = std::fs::read_to_string(&fasta_file).unwrap();
        assert_eq!(contents, ">1 length=19 circular=true topology=circular\nAGCATCGACATCGACTACG\n");
    }

    #[test]
    fn test_gfa2fasta_9() {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_9());
        save_graph_to_fasta(&graph, &fasta_file);
        let contents = std::fs::read_to_string(&fasta_file).unwrap();
        assert_eq!(contents, ">1 length=19 circular=false topology=linear\nAGCATCGACATCGACTACG\n");
    }

    #[test]
    fn test_gfa2fasta_10() {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_10());
        save_graph_to_fasta(&graph, &fasta_file);
        let contents = std::fs::read_to_string(&fasta_file).unwrap();
        assert_eq!(contents, ">1 length=19 circular=false topology=linear\nAGCATCGACATCGACTACG\n");
    }

    #[test]
    fn test_gfa2fasta_13() {
        let temp_dir = tempdir().unwrap();
        let fasta_file = temp_dir.path().join("temp.fasta");
        let (graph, _) = UnitigGraph::from_gfa_lines(&get_test_gfa_13());
        save_graph_to_fasta(&graph, &fasta_file);
        let contents = std::fs::read_to_string(&fasta_file).unwrap();
        assert_eq!(contents, ">1 length=19\nAGCATCGACATCGACTACG\n");
    }
}
