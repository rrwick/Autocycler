// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::misc::quit_with_error;
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn cluster(in_gfa: PathBuf, out_dir: PathBuf) {
    section_header("Starting autocycler cluster");
    explanation("This command will take a compacted De Bruijn graph (made by autocycler \
                 compress) and cluster the contigs based on their similarity.");
    print_settings(&in_gfa, &out_dir);
    create_output_dir(&out_dir);
    let (unitig_graph, sequences) = load_graph(&in_gfa);
    let asymmetrical_distances = pairwise_contig_distances(unitig_graph, &sequences);
    save_distance_matrix(&asymmetrical_distances, &sequences, &out_dir,
                         "distances_asymmetrical.phylip");
    let symmetrical_distances = make_symmetrical_distances(&asymmetrical_distances, &sequences);
    save_distance_matrix(&symmetrical_distances, &sequences, &out_dir,
                         "distances_symmetrical.phylip");
    // TODO: cluster using distance
}


fn print_settings(in_gfa: &PathBuf, out_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_dir {}", out_dir.display());
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
    eprintln!("{} links", unitig_graph.get_link_count());
    eprintln!("{} contigs", sequences.len());
    (unitig_graph, sequences)
}


fn pairwise_contig_distances(unitig_graph: UnitigGraph, sequences: &Vec<Sequence>) -> HashMap<(u16, u16), f64> {
    let unitig_lengths: HashMap<u32, u32> = unitig_graph.unitigs.iter()
        .map(|rc| {let u = rc.borrow(); (u.number, u.length())}).collect();
    let sequence_unitigs: HashMap<u16, HashSet<u32>> = sequences.iter()
        .map(|s| (s.id, unitig_graph.get_unitig_path_for_sequence(s).iter().map(|(number, _)| *number).collect::<HashSet<u32>>()))
        .collect();
    let mut distances: HashMap<(u16, u16), f64> = HashMap::new();
    for seq_a in sequences {
        let a = sequence_unitigs.get(&seq_a.id).unwrap();
        let a_len = total_unitig_length(&a, &unitig_lengths) as f64;
        for seq_b in sequences {
            let b = sequence_unitigs.get(&seq_b.id).unwrap();
            let ab: HashSet<u32> = a.intersection(&b).cloned().collect();
            let ab_len = total_unitig_length(&ab, &unitig_lengths) as f64;
            let distance = 1.0 - (ab_len / a_len);
            distances.insert((seq_a.id, seq_b.id), distance);
        }
    }
    distances
}


fn total_unitig_length(unitigs: &HashSet<u32>, unitig_lengths: &HashMap<u32, u32>) -> u32 {
    unitigs.iter().map(|u| unitig_lengths.get(u).unwrap()).sum()
}


fn save_distance_matrix(distances: &HashMap<(u16, u16), f64>, sequences: &Vec<Sequence>,
                        out_dir: &PathBuf, filename: &str) {
    let file_path = out_dir.join(filename);
    let mut f = File::create(&file_path).unwrap();
    write!(f, "{}\n", sequences.len()).unwrap();
    for seq_a in sequences {
        write!(f, "{}", seq_a).unwrap();
        for seq_b in sequences {
            write!(f, "\t{:.8}", distances.get(&(seq_a.id, seq_b.id)).unwrap()).unwrap();
        }
        write!(f, "\n").unwrap();
    }
}


fn make_symmetrical_distances(asymmetrical_distances: &HashMap<(u16, u16), f64>, sequences: &Vec<Sequence>) -> HashMap<(u16, u16), f64> {
    // This function takes in an asymmetrical distance matrix (where A vs B is not necessarily
    // equal to B vs A) and makes it symmetric by setting each distance (both orders) to the
    // minimum distance of the two orders.
    let mut symmetrical_distances: HashMap<(u16, u16), f64> = HashMap::new();
    for seq_a in sequences {
        for seq_b in sequences {
            let a_vs_b = asymmetrical_distances.get(&(seq_a.id, seq_b.id)).unwrap();
            let b_vs_a = asymmetrical_distances.get(&(seq_b.id, seq_a.id)).unwrap();
            // let distance = f64::min(*a_vs_b, *b_vs_a);
            let distance = f64::max(*a_vs_b, *b_vs_a);
            // let distance = (a_vs_b + b_vs_a) / 2.0;
            symmetrical_distances.insert((seq_a.id, seq_b.id), distance);
        }
    }
    symmetrical_distances
}
