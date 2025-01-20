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

use colored::Colorize;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::log::{section_header, explanation};
use crate::metrics::{CombineMetrics, ResolvedClusterDetails};
use crate::misc::{check_if_file_exists, create_dir};
use crate::unitig_graph::UnitigGraph;


pub fn combine(autocycler_dir: PathBuf, in_gfas: Vec<PathBuf>) {
    let combined_gfa = autocycler_dir.join("consensus_assembly.gfa");
    let combined_fasta = autocycler_dir.join("consensus_assembly.fasta");
    let combined_yaml = autocycler_dir.join("consensus_assembly.yaml");

    check_settings(&in_gfas);
    if let Some(parent) = combined_gfa.parent() {
        create_dir(parent);
    }
    starting_message();
    print_settings(&autocycler_dir, &in_gfas);

    // TODO: add an optional argument for reads. When used, the reads will be used to set unitig
    //       depths (instead of input assembly count). Find unique k-mers in the combined assembly
    //       and then count the occurrences of those k-mers in the reads.

    // TODO: sort the input GFA files in order of decreasing size. They may already be in this
    //       order (from the clustering step), but not necessarily (e.g. due to small plasmid
    //       duplication).

    let mut metrics = CombineMetrics::default();
    combine_clusters(&in_gfas, &combined_gfa, &combined_fasta, &mut metrics);
    metrics.save_to_yaml(&combined_yaml);
    finished_message(&combined_gfa, &combined_fasta, &metrics);
}


fn check_settings(in_gfas: &Vec<PathBuf>) {
    for gfa in in_gfas {
        check_if_file_exists(gfa);
    }
}


fn starting_message() {
    section_header("Starting autocycler combine");
    explanation("This command combines different clusters into a single assembly file.");
}


fn print_settings(autocycler_dir: &Path, in_gfas: &[PathBuf]) {
    eprintln!("Settings:");
    eprintln!("  --autocycler_dir {}", autocycler_dir.display());
    eprintln!("  --in_gfas {}", in_gfas[0].display());
    for gfa in &in_gfas[1..] {
        eprintln!("            {}", gfa.display());
    }
    eprintln!();
}


fn finished_message(combined_gfa: &Path, combined_fasta: &Path, metrics: &CombineMetrics) {
    section_header("Finished!");
    eprintln!("Combined graph: {}", combined_gfa.display());
    eprintln!("Combined fasta: {}", combined_fasta.display());
    eprintln!();
    if metrics.consensus_assembly_fully_resolved {
        eprintln!("{}", "Consensus assembly is fully resolved 😄".green().bold());
    } else {
        eprintln!("{}", "One or more clusters failed to fully resolve 😟".red().bold());
    }
    eprintln!();
}


fn combine_clusters(in_gfas: &Vec<PathBuf>, combined_gfa: &Path, combined_fasta: &Path,
                    metrics: &mut CombineMetrics) {
    section_header("Combining clusters");
    explanation("This command combines different clusters into a single assembly file.");
    let mut gfa_file = File::create(combined_gfa).unwrap();
    let mut fasta_file = File::create(combined_fasta).unwrap();
    writeln!(gfa_file, "H\tVN:Z:1.0").unwrap();
    metrics.consensus_assembly_fully_resolved = true;
    let mut offset = 0;
    for gfa in in_gfas {
        eprintln!("{}", gfa.display());
        let (graph, _) = UnitigGraph::from_gfa_file(gfa);
        graph.print_basic_graph_info();
        for unitig in &graph.unitigs {
            let unitig = unitig.borrow();
            let unitig_num = unitig.number + offset;
            let unitig_seq = String::from_utf8_lossy(&unitig.forward_seq);
            let circ = if unitig.is_isolated_and_circular() { " circular=true".to_string() }
                                                       else { "".to_string() };
            let depth_tag = format!("\tDP:f:{:.2}", unitig.depth);
            let mut colour_tag = unitig.colour_tag(true);
            if colour_tag.is_empty() {
                colour_tag = "\tCL:z:orangered".to_string();
            }
            writeln!(gfa_file, "S\t{}\t{}{}{}",
                     unitig_num, unitig_seq, depth_tag, colour_tag).unwrap();
            writeln!(fasta_file, ">{} length={}{}", unitig_num, unitig.length(), circ).unwrap();
            writeln!(fasta_file, "{}", unitig_seq).unwrap();
        }
        for (a, a_strand, b, b_strand) in &graph.get_links_for_gfa(offset) {
            writeln!(gfa_file, "L\t{}\t{}\t{}\t{}\t0M", a, a_strand, b, b_strand).unwrap();
        }
        offset += graph.max_unitig_number();
        let component_length = graph.total_length();
        let unitig_count = graph.unitigs.len() as u32;
        metrics.consensus_assembly_bases += component_length;
        metrics.consensus_assembly_unitigs += unitig_count;
        let cluster_metrics = ResolvedClusterDetails { length: component_length,
                                                       unitigs: unitig_count,
                                                       topology: graph.topology() };
        metrics.consensus_assembly_clusters.push(cluster_metrics);
        if unitig_count > 1 { metrics.consensus_assembly_fully_resolved = false; }
    }
}
