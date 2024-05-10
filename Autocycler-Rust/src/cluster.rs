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
use crate::misc::{format_float, quit_with_error, usize_division_rounded};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn cluster(in_gfa: PathBuf, out_dir: PathBuf, eps: f64, minpts: Option<usize>) {
    check_settings(eps, &minpts);
    section_header("Starting autocycler cluster");
    explanation("This command will take a compacted De Bruijn graph (made by autocycler \
                 compress) and cluster the contigs based on their similarity.");
    print_settings(&in_gfa, &out_dir, eps, &minpts);
    create_output_dir(&out_dir);
    let (unitig_graph, mut sequences) = load_graph(&in_gfa);
    let asymmetrical_distances = pairwise_contig_distances(unitig_graph, &sequences);
    save_distance_matrix(&asymmetrical_distances, &sequences, &out_dir,
                         "distances_asymmetrical.phylip");
    let symmetrical_distances = make_symmetrical_distances(&asymmetrical_distances, &sequences);
    save_distance_matrix(&symmetrical_distances, &sequences, &out_dir,
                         "distances_symmetrical.phylip");
    dbscan_cluster(&symmetrical_distances, &mut sequences, eps, minpts);
    // TODO: exclude any clusters that aren't spread over enough assemblies
    print_clusters(&sequences);
}


fn check_settings(eps: f64, minpts: &Option<usize>) {
    if eps <= 0.0 || eps >= 1.0 {
        quit_with_error("--eps must be between 0 and 1 (exclusive)");
    }
    if minpts.is_some() {
        if minpts.unwrap() < 2 {
            quit_with_error("--minpts must be 2 or greater");
        }
    }
}


fn print_settings(in_gfa: &PathBuf, out_dir: &PathBuf, eps: f64, minpts: &Option<usize>) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!("  --eps {}", format_float(eps));
    if minpts.is_none() {
        eprintln!("  --minpts (automatically set)");
    } else {
        eprintln!("  --minpts {}", minpts.unwrap());
    }
    eprintln!();
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
    let assembly_count = sequences.iter().map(|s| &s.filename).collect::<HashSet<_>>().len();
    eprintln!("{} unitigs", unitig_graph.unitigs.len());
    eprintln!("{} links", unitig_graph.get_link_count());
    eprintln!("{} contigs from {} assemblies", sequences.len(), assembly_count);
    eprintln!();
    (unitig_graph, sequences)
}


fn pairwise_contig_distances(unitig_graph: UnitigGraph, sequences: &Vec<Sequence>) -> HashMap<(u16, u16), f64> {
    section_header("Calculating pairwise distances");
    explanation("Every pairwise distance between contigs is calculated based on the similarity of \
                 their paths through the compressed sequence graph.");
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
            let distance = (a_vs_b + b_vs_a) / 2.0;
            symmetrical_distances.insert((seq_a.id, seq_b.id), distance);
        }
    }
    symmetrical_distances
}


fn dbscan_cluster(distances: &HashMap<(u16, u16), f64>, sequences: &mut Vec<Sequence>,
                  eps: f64, min_pts_option: Option<usize>) {
                    section_header("Clustering sequences");
                    explanation("Contigs are clustered using the DBSCAN* algorithm. This is a \
                                 variant of the DBSCAN algorithm where all border points are \
                                 treated as noise. Contigs will be either put into a cluster \
                                 (defined by sufficient density) or be classified as noise (no \
                                 cluster).");
    // Based on https://en.wikipedia.org/wiki/DBSCAN#Algorithm, but I'm using DBSCAN* instead of
    // DBSCAN, so border points become noise. This serves to make the algorithm deterministic (not
    // sensitive to the order of sequences).
    eprintln!("ε value: {}", format_float(eps));
    let min_pts = set_min_pts(min_pts_option, sequences);
    let ids: Vec<u16> = sequences.iter().map(|s| s.id).collect();
    let mut index: HashMap<u16, &mut Sequence> = sequences.iter_mut().map(|s| (s.id, s)).collect();
    let mut c = 0;
    for p in &ids {
        let p_seq: &mut Sequence = index.get_mut(p).unwrap();
        if p_seq.cluster != 0 { continue; }
        let n = range_query(&ids, distances, *p, eps);
        if n.len() < min_pts {
            p_seq.cluster = -1;
            continue;
        }
        c += 1;
        p_seq.cluster = c;
        let mut seed_set: Vec<u16> = n.into_iter().filter(|&id| id != *p).collect();
        while let Some(q) = seed_set.pop() {
            let q_seq = index.get_mut(&q).unwrap();
            if q_seq.cluster != 0 { continue; }
            q_seq.cluster = c;
            let n_prime = range_query(&ids, distances, q, eps);
            if n_prime.len() >= min_pts {
                seed_set.extend(n_prime);
            }
        }
    }
}


fn range_query(ids: &Vec<u16>, distances: &HashMap<(u16, u16), f64>, q: u16, eps: f64) -> Vec<u16> {
    let mut neighbours = Vec::new();
    for p in ids {
        if distances.get(&(q, *p)).unwrap() <= &eps {
            neighbours.push(*p);
        }
    }
    neighbours
}


fn set_min_pts(min_pts_option: Option<usize>, sequences: &mut Vec<Sequence>) -> usize {
    // This function automatically sets the minPts parameter, if the user didn't supply one.
    if min_pts_option.is_some() {
        eprintln!("minPts value: {} (user-suppled)\n", min_pts_option.unwrap());
        return min_pts_option.unwrap();
    }
    let assembly_count = sequences.iter().map(|s| &s.filename).collect::<HashSet<_>>().len();
    let mut auto_min_pts = usize_division_rounded(assembly_count, 4);
    if auto_min_pts < 3 {
        auto_min_pts = 3;
    }
    eprintln!("minPts value: {} (automatically set)\n", auto_min_pts);
    auto_min_pts
}


fn print_clusters(sequences: &Vec<Sequence>) {
    let max_cluster = sequences.iter().map(|s| s.cluster).max().unwrap();
    for c in 1..max_cluster+1 {
        eprintln!("Cluster {}:", c);
        for s in sequences {
            if s.cluster == c {
                eprintln!("  {}", s);
            }
        }
        eprintln!();
    }
    let noise_sequences: Vec<_> = sequences.iter().filter(|s| s.cluster == -1).collect();
    if noise_sequences.len() == 0 {
        eprintln!("No noise")
    } else {
        eprintln!("Noise:");
        for s in noise_sequences {
            eprintln!("  {}", s);
        }
    }
    eprintln!("\n");
}
