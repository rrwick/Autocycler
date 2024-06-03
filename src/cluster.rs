// This file contains the code for the autocycler cluster subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use maud::{html, Markup, DOCTYPE, PreEscaped};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists, format_float, median_f64, median_usize,
                  round_float, quit_with_error, usize_division_rounded};
use crate::sequence::Sequence;
use crate::trim_path_overlap::trim_path_overlap;
use crate::unitig_graph::UnitigGraph;


pub fn cluster(out_dir: PathBuf, min_overlap_id: f64, max_overlap_unitigs: u32,
               eps: Option<f64>, minpts: Option<usize>, max_len_var: f64) {
    let gfa = out_dir.join("01_unitig_graph.gfa");
    let clusters_html = out_dir.join("04_clusters.html");
    let clusters_tsv = out_dir.join("05_clusters.tsv");
    let clean_gfa = out_dir.join("06_cleaned_graph.gfa");
    check_settings(&out_dir, &gfa, min_overlap_id, &eps, &minpts, max_len_var);
    starting_message();
    print_settings(&out_dir, min_overlap_id, max_overlap_unitigs, &eps, &minpts, max_len_var);
    let (mut unitig_graph, sequences) = load_graph(&gfa);
    let mut sequences = trim_path_overlap(&mut unitig_graph, &sequences, min_overlap_id, max_overlap_unitigs);
    // TODO: add in some way for a user to manually define clusters.
    let distances = pairwise_contig_distances(&unitig_graph, &sequences, &out_dir);
    let minpts = dbscan_cluster(&distances, &mut sequences, eps, minpts);
    select_clusters(&mut sequences, minpts);
    exclude_sequences_by_length(&mut sequences, max_len_var);
    reorder_clusters(&mut sequences);
    print_clusters(&sequences);
    save_html_table(&sequences, &clusters_html);
    save_tsv_table(&sequences, &clusters_tsv);
    let sequences = remove_excluded_contigs_from_graph(&mut unitig_graph, &sequences);
    unitig_graph.save_gfa(&clean_gfa, &sequences).unwrap();
    finished_message(&clusters_tsv, &clusters_html, &clean_gfa);
}


fn check_settings(out_dir: &PathBuf, gfa: &PathBuf, min_overlap_id: f64,
                  eps: &Option<f64>, minpts: &Option<usize>, max_len_var: f64) {
    check_if_dir_exists(&out_dir);
    check_if_file_exists(&gfa);
    if min_overlap_id < 0.0 || min_overlap_id > 1.0 {
        quit_with_error("--min_overlap_id must be between 0 and 1 (inclusive)");
    }
    if max_len_var < 0.0 {
        quit_with_error("--max_len_var must be a positive value");
    }
    if eps.is_some() && (eps.unwrap() <= 0.0 || eps.unwrap() >= 1.0) {
        quit_with_error("--eps must be between 0 and 1 (exclusive)");
    }
    if minpts.is_some() && minpts.unwrap() < 2 {
        quit_with_error("--minpts must be 2 or greater");
    }
}


fn starting_message() {
    section_header("Starting autocycler cluster");
    explanation("This command takes a compacted De Bruijn graph (made by autocycler compress), \
                 trims any start-end overlap from the sequences and then clusters the sequences \
                 based on their similarity. Ideally, each cluster will then contain sequences \
                 which can be combined into a consensus.");
}


fn print_settings(out_dir: &PathBuf, min_overlap_id: f64, max_overlap_unitigs: u32,
                  eps: &Option<f64>, minpts: &Option<usize>, max_len_var: f64) {
    eprintln!("Settings:");
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!("  --min_overlap_id {}", format_float(min_overlap_id));
    eprintln!("  --max_overlap_unitigs {}", max_overlap_unitigs);
    if eps.is_none() {
        eprintln!("  --eps (automatically set)");
    } else {
        eprintln!("  --eps {}", format_float(eps.unwrap()));
    }
    if minpts.is_none() {
        eprintln!("  --minpts (automatically set)");
    } else {
        eprintln!("  --minpts {}", minpts.unwrap());
    }
    eprintln!("  --max_len_var {}", format_float(max_len_var));
    eprintln!();
}


pub fn load_graph(gfa: &PathBuf) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading graph");
    explanation("The compressed sequence graph is now loaded into memory.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(&gfa);
    unitig_graph.print_basic_graph_info();
    (unitig_graph, sequences)
}


fn pairwise_contig_distances(unitig_graph: &UnitigGraph, sequences: &Vec<Sequence>, out_dir: &PathBuf) -> HashMap<(u16, u16), f64> {
    section_header("Pairwise distances");
    explanation("Every pairwise distance between contigs is calculated based on the similarity of \
                 their paths through the graph.");
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
    eprintln!("{} sequences, {} total pairwise distances", sequences.len(), distances.len());
    eprintln!();
    eprintln!("Saving distance matrix files:");
    save_distance_matrix(&distances, &sequences, &out_dir, "02_pairwise_distances_asymmetrical.phylip");
    let distances = make_symmetrical_distances(&distances, &sequences);
    save_distance_matrix(&distances, &sequences, &out_dir, "03_pairwise_distances_symmetrical.phylip");
    eprintln!();
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
    eprintln!("  {}", file_path.display());
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
                  eps_option: Option<f64>, min_pts_option: Option<usize>) -> usize {
    section_header("Clustering sequences");
    explanation("Contigs are clustered using the DBSCAN* algorithm. Contigs will be either put \
                 into a cluster (defined by sufficient density) or be classified as noise (no \
                 cluster).");
    // Based on https://en.wikipedia.org/wiki/DBSCAN#Algorithm, but I'm using DBSCAN* instead of
    // DBSCAN, so border points become noise. This serves to make the algorithm deterministic (not
    // sensitive to the order of sequences).
    let eps = set_eps(eps_option, distances);
    let min_pts = set_min_pts(min_pts_option, sequences);
    eprintln!();
    let ids: Vec<u16> = sequences.iter().map(|s| s.id).collect();
    let mut index: HashMap<u16, &mut Sequence> = sequences.iter_mut().map(|s| (s.id, s)).collect();
    let mut c = 0;
    for p in &ids {
        let p_seq: &mut Sequence = index.get_mut(p).unwrap();
        if p_seq.cluster != 0 { continue; }
        let n = range_query(&ids, distances, *p, eps);
        if n.len() < min_pts {
            p_seq.cluster = -1;
            p_seq.cluster_fail_reason = "dbscan noise".to_string();
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
    min_pts
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


fn set_eps(eps_option: Option<f64>, distances: &HashMap<(u16, u16), f64>) -> f64 {
    // This function automatically sets the ε parameter, if the user didn't supply one.
    // The auto-set value be 2x the median close distance, where close distances are all distances
    // less than 0.2. The auto-set value can't be less than 0.01 or more than 0.1.
    if eps_option.is_some() {
        eprintln!("ε parameter: {} (user-suppled)", eps_option.unwrap());
        return eps_option.unwrap();
    }
    let close_distances: Vec<_> = distances.values().cloned().filter(|&val| val > 0.0 && val < 0.2).collect();
    let mut auto_eps = 0.05;
    if close_distances.len() > 0 {
        auto_eps = 2.0 * round_float(median_f64(&close_distances), 3);
    }
    if auto_eps < 0.01 {
        auto_eps = 0.01
    }
    if auto_eps > 0.1 {
        auto_eps = 0.1
    }
    eprintln!("ε parameter: {} (automatically set)", format_float(auto_eps));
    auto_eps
}


fn set_min_pts(min_pts_option: Option<usize>, sequences: &Vec<Sequence>) -> usize {
    // This function automatically sets the minPts parameter, if the user didn't supply one.
    // The auto-set value will be one-quarter of the assembly count (rounded) but no less than 3.
    if min_pts_option.is_some() {
        eprintln!("minPts parameter: {} (user-suppled)", min_pts_option.unwrap());
        return min_pts_option.unwrap();
    }
    let assembly_count = get_assembly_count(&sequences);
    let mut auto_min_pts = usize_division_rounded(assembly_count, 4);
    if auto_min_pts < 3 {
        auto_min_pts = 3;
    }
    eprintln!("minPts parameter: {} (automatically set)", auto_min_pts);
    auto_min_pts
}


fn select_clusters(sequences: &mut Vec<Sequence>, minpts: usize) {
    // Clusters that appear in too few assemblies are excluded.
    for c in 1..get_max_cluster(sequences)+1 {
        let assemblies: HashSet<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.filename.clone()).collect();
        if assemblies.len() < minpts {
            for s in &mut *sequences {
                if s.cluster == c {
                    s.cluster = -1;
                    s.cluster_fail_reason = "cluster in too few assemblies".to_string();
                }
            }
        }
    }
}

fn exclude_sequences_by_length(sequences: &mut Vec<Sequence>, max_len_var: f64) {
    // This function will exclude sequences to ensure that within-cluster sequence lengths don't
    // have too much variation.

    for c in 1..get_max_cluster(sequences)+1 {
        let lengths: Vec<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.length).collect();
        let cluster_median = median_usize(&lengths);
        let (cluster_min, cluster_max) = cluster_min_max_lengths(cluster_median, max_len_var);
        for s in &mut *sequences {
            if s.cluster == c && s.length < cluster_min {
                s.cluster = -1;
                s.cluster_fail_reason = "too short".to_string();
            }
            if s.cluster == c && s.length > cluster_max {
                s.cluster = -1;
                s.cluster_fail_reason = "too long".to_string();
            }
        }
    }
}


fn cluster_min_max_lengths(median_length: usize, max_len_var: f64) -> (usize, usize) {
    let cluster_min = (median_length as f64 * (1.0 - max_len_var)).round() as usize;
    let cluster_max = (median_length as f64 * (1.0 + max_len_var)).round() as usize;
    (cluster_min, cluster_max)
}


fn reorder_clusters(sequences: &mut Vec<Sequence>) {
    // Reorder clusters based on their median sequence length (large to small).
    let mut cluster_lengths = HashMap::new();
    for c in 1..get_max_cluster(sequences)+1 {
        let lengths: Vec<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.length).collect();
        cluster_lengths.insert(c, median_usize(&lengths));
    }
    let mut sorted_cluster_lengths: Vec<_> = cluster_lengths.iter().collect();
    sorted_cluster_lengths.sort_by(|a, b| {match b.1.cmp(a.1) {std::cmp::Ordering::Equal => a.0.cmp(b.0), other => other,}});
    let mut old_to_new = HashMap::new();
    let mut new_c: i32 = 0;
    for (old_c, _length) in sorted_cluster_lengths {
        new_c += 1;
        old_to_new.insert(old_c, new_c);
    }
    for s in sequences {
        if s.cluster < 1 {continue;}
        s.cluster = *old_to_new.get(&s.cluster).unwrap();
    }
}


fn print_clusters(sequences: &Vec<Sequence>) {
    for c in 1..get_max_cluster(sequences)+1 {
        eprintln!("Cluster {}:", c);
        for s in sequences {
            assert!(s.cluster != 0);
            if s.cluster == c {
                eprintln!("  {}", s);
            }
        }
        eprintln!();
    }
    let noise_sequences: Vec<_> = sequences.iter().filter(|s| s.cluster == -1).collect();
    if noise_sequences.len() == 0 {
        eprintln!("No contigs excluded")
    } else {
        eprintln!("Excluded:");
        for s in noise_sequences {
            eprintln!("  {}", s);
        }
    }
    eprintln!();
}


fn save_html_table(sequences: &Vec<Sequence>, file_path: &PathBuf) {
    let markup = create_html(sequences);
    std::fs::write(file_path.clone(), markup.into_string()).unwrap();
}


fn save_tsv_table(sequences: &Vec<Sequence>, file_path: &PathBuf) {
    let mut f = File::create(&file_path).unwrap();
    write!(f, "assembly\tcontig_name\tlength\tcluster\tfail_reason\n").unwrap();
    for c in 1..get_max_cluster(sequences)+1 {
        for s in sequences {
            if s.cluster == c {
                write!(f, "{}\t{}\t{}\t{}\t{}\n", s.filename, s.contig_name(), s.length, c, s.cluster_fail_reason).unwrap();
            }
        }
    }
    for s in sequences {
        if s.cluster == -1 {
            write!(f, "{}\t{}\t{}\tnone\t{}\n", s.filename, s.contig_name(), s.length, s.cluster_fail_reason).unwrap();
        }
    }
}


pub fn remove_excluded_contigs_from_graph(graph: &mut UnitigGraph, sequences: &Vec<Sequence>) -> Vec<Sequence> {
    section_header("Cleaning graph");
    explanation("Excluded contigs (those which could not be clustered) are now removed from the \
                 unitig graph.");
    let seqs_to_remove: Vec<_> = sequences.iter().filter(|s| s.cluster == -1).collect();
    for s in seqs_to_remove {
        graph.remove_sequence_from_graph(s.id);
    }
    graph.remove_zero_depth_unitigs();
    let sequences = sequences.iter().filter(|s| s.cluster != -1).cloned().collect();
    eprintln!();
    merge_linear_paths(graph, &sequences);
    graph.print_basic_graph_info();
    graph.renumber_unitigs();
    sequences
}


fn get_assembly_count(sequences: &Vec<Sequence>) -> usize {
    sequences.iter().map(|s| &s.filename).collect::<HashSet<_>>().len()
}


fn get_max_cluster(sequences: &Vec<Sequence>) -> i32 {
    sequences.iter().map(|s| s.cluster).max().unwrap()
}


fn get_assemblies(sequences: &[Sequence]) -> Vec<String> {
    // Returns the assembly filenames, without duplicates, in the same order as they are first seen
    // in Sequence objects.
    let mut seen = HashSet::new();
    let mut assemblies = Vec::new();
    for s in sequences.iter() {
        let a = s.filename.clone();
        if seen.insert(a.clone()) {
            assemblies.push(a);
        }
    }
    assemblies
}


fn create_html(sequences: &Vec<Sequence>) -> Markup {
    let headers: Vec<String> = std::iter::once("Assembly".to_string())
        .chain((1..=get_max_cluster(&sequences)).map(|c| format!("cluster {}", c)))
        .chain(std::iter::once("excluded".to_string())).collect();
    let mut rows = vec![];
    for a in get_assemblies(sequences) {
        let mut row = vec![a.clone()];
        for c in 1..get_max_cluster(sequences)+1 {
            let cell_seqs: Vec<String> = sequences.iter().filter(|s| s.filename == a && s.cluster == c)
                .map(|s| format!("{} ({} bp)", s.contig_name(), s.length)).collect();
            row.push(cell_seqs.join("<br>"));
        }
        let cell_seqs: Vec<String> = sequences.iter().filter(|s| s.filename == a && s.cluster == -1)
            .map(|s| format!("{} ({} bp)", s.contig_name(), s.length)).collect();
        row.push(cell_seqs.join("<br>"));
        rows.push(row);
    }

    html! {
        (DOCTYPE)
        html {
            head {
                meta charset="utf-8";
                title { "Autocycler clustering" }
                link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css";
                style { "td { vertical-align: middle !important; }" }
            }
            body {
                h2 { "Autocycler clustering" }
                table class="table" {
                    thead class="thead-dark" { tr { @for h in headers { th { (h) } } } }
                    tbody { @for row in rows { tr { @for d in row { td { (PreEscaped(d)) } } } } }
                }
            }
        }
    }
}


fn finished_message(cluster_tsv: &PathBuf, cluster_html: &PathBuf, clean_gfa: &PathBuf) {
    section_header("Finished!");
    explanation("You can now run autocycler resolve to resolve repeats.");
    eprintln!("Cluster table (HTML): {}", cluster_html.display());
    eprintln!("Cluster table (TSV):  {}", cluster_tsv.display());
    eprintln!("Cleaned unitig graph: {}", clean_gfa.display());
    eprintln!();
}


#[cfg(test)]
mod tests {
    use super::*;

    fn assert_almost_eq(a: f64, b: f64, epsilon: f64) {
        assert!((a - b).abs() < epsilon, "Numbers are not within {:?} of each other: {} vs {}", epsilon, a, b);
    }

    #[test]
    fn test_get_assembly_count() {
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 1)];
        assert_eq!(get_assembly_count(&sequences), 2);
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_3.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_4.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(5, "A".to_string(), "assembly_5.fasta".to_string(), "contig_1".to_string(), 1)];
        assert_eq!(get_assembly_count(&sequences), 5);
    }

    #[test]
    fn test_get_max_cluster() {
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                                 Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 1)];
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3; sequences[3].cluster = 1; sequences[4].cluster = 2;
        assert_eq!(get_max_cluster(&sequences), 3);
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3; sequences[3].cluster = 5; sequences[4].cluster = 4;
        assert_eq!(get_max_cluster(&sequences), 5);
    }

    #[test]
    fn test_get_assemblies() {
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 1)];
        assert_eq!(get_assemblies(&sequences), vec!["assembly_1.fasta".to_string(), "assembly_2.fasta".to_string()]);
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_4.fasta".to_string(), "contig_3".to_string(), 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_4.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(5, "A".to_string(), "assembly_3.fasta".to_string(), "contig_2".to_string(), 1)];
        assert_eq!(get_assemblies(&sequences), vec!["assembly_2.fasta".to_string(), "assembly_1.fasta".to_string(), "assembly_4.fasta".to_string(), "assembly_3.fasta".to_string()]);
    }

    #[test]
    fn test_select_clusters() {
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                                 Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 1),
                                 Sequence::new_with_seq(6, "A".to_string(), "assembly_2.fasta".to_string(), "contig_3".to_string(), 1),
                                 Sequence::new_with_seq(7, "A".to_string(), "assembly_3.fasta".to_string(), "contig_1".to_string(), 1),
                                 Sequence::new_with_seq(8, "A".to_string(), "assembly_3.fasta".to_string(), "contig_2".to_string(), 1)];
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3;
        sequences[3].cluster = 1; sequences[4].cluster = 2; sequences[5].cluster = 3;
        sequences[6].cluster = 1; sequences[7].cluster = 2;
        select_clusters(&mut sequences, 1);
        assert_eq!(sequences[0].cluster, 1); assert_eq!(sequences[1].cluster, 2); assert_eq!(sequences[2].cluster, 3);
        assert_eq!(sequences[3].cluster, 1); assert_eq!(sequences[4].cluster, 2); assert_eq!(sequences[5].cluster, 3);
        assert_eq!(sequences[6].cluster, 1); assert_eq!(sequences[7].cluster, 2);
        
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3;
        sequences[3].cluster = 1; sequences[4].cluster = 2; sequences[5].cluster = 3;
        sequences[6].cluster = 1; sequences[7].cluster = 2;
        select_clusters(&mut sequences, 2);
        assert_eq!(sequences[0].cluster, 1); assert_eq!(sequences[1].cluster, 2); assert_eq!(sequences[2].cluster, 3);
        assert_eq!(sequences[3].cluster, 1); assert_eq!(sequences[4].cluster, 2); assert_eq!(sequences[5].cluster, 3);
        assert_eq!(sequences[6].cluster, 1); assert_eq!(sequences[7].cluster, 2);
        
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3;
        sequences[3].cluster = 1; sequences[4].cluster = 2; sequences[5].cluster = 3;
        sequences[6].cluster = 1; sequences[7].cluster = 2;
        select_clusters(&mut sequences, 3);
        assert_eq!(sequences[0].cluster, 1); assert_eq!(sequences[1].cluster, 2); assert_eq!(sequences[2].cluster, -1);
        assert_eq!(sequences[3].cluster, 1); assert_eq!(sequences[4].cluster, 2); assert_eq!(sequences[5].cluster, -1);
        assert_eq!(sequences[6].cluster, 1); assert_eq!(sequences[7].cluster, 2);
        
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3;
        sequences[3].cluster = 1; sequences[4].cluster = 2; sequences[5].cluster = 3;
        sequences[6].cluster = 4; sequences[7].cluster = 4;
        select_clusters(&mut sequences, 1);
        assert_eq!(sequences[0].cluster, 1); assert_eq!(sequences[1].cluster, 2); assert_eq!(sequences[2].cluster, 3);
        assert_eq!(sequences[3].cluster, 1); assert_eq!(sequences[4].cluster, 2); assert_eq!(sequences[5].cluster, 3);
        assert_eq!(sequences[6].cluster, 4); assert_eq!(sequences[7].cluster, 4);
        
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3;
        sequences[3].cluster = 1; sequences[4].cluster = 2; sequences[5].cluster = 3;
        sequences[6].cluster = 4; sequences[7].cluster = 4;
        select_clusters(&mut sequences, 2);
        assert_eq!(sequences[0].cluster, 1); assert_eq!(sequences[1].cluster, 2); assert_eq!(sequences[2].cluster, 3);
        assert_eq!(sequences[3].cluster, 1); assert_eq!(sequences[4].cluster, 2); assert_eq!(sequences[5].cluster, 3);
        assert_eq!(sequences[6].cluster, -1); assert_eq!(sequences[7].cluster, -1);
    }

    #[test]
    fn test_reorder_clusters() {
        let mut sequences = vec![Sequence::new_with_seq(1, "CGCGA".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 5),
                                 Sequence::new_with_seq(2, "T".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1),
                                 Sequence::new_with_seq(3, "AACGACTACG".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 10),
                                 Sequence::new_with_seq(4, "CGCGA".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 5),
                                 Sequence::new_with_seq(5, "T".to_string(), "assembly_2.fasta".to_string(), "contig_3".to_string(), 1),
                                 Sequence::new_with_seq(6, "AACGACTACG".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 10)];
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3;
        sequences[3].cluster = 1; sequences[4].cluster = 2; sequences[5].cluster = 3;
        reorder_clusters(&mut sequences);
        assert_eq!(sequences[0].cluster, 2); assert_eq!(sequences[1].cluster, 3); assert_eq!(sequences[2].cluster, 1);
        assert_eq!(sequences[3].cluster, 2); assert_eq!(sequences[4].cluster, 3); assert_eq!(sequences[5].cluster, 1);
        reorder_clusters(&mut sequences);
        assert_eq!(sequences[0].cluster, 2); assert_eq!(sequences[1].cluster, 3); assert_eq!(sequences[2].cluster, 1);
        assert_eq!(sequences[3].cluster, 2); assert_eq!(sequences[4].cluster, 3); assert_eq!(sequences[5].cluster, 1);
    }

    #[test]
    fn test_set_eps() {
        let distances = HashMap::from_iter(vec![]);
        assert_almost_eq(set_eps(Some(0.123), &distances), 0.123, 1e-8);
        assert_almost_eq(set_eps(Some(0.321), &distances), 0.321, 1e-8);

        let distances = HashMap::from_iter(vec![((1, 1), 0.00),((1, 2), 0.01),((1, 3), 0.90),((1, 4), 0.02),
                                                ((2, 1), 0.01),((2, 2), 0.00),((2, 3), 0.80),((2, 4), 0.01),
                                                ((3, 1), 0.90),((3, 2), 0.80),((3, 3), 0.00),((3, 4), 0.90),
                                                ((4, 1), 0.02),((4, 2), 0.01),((4, 3), 0.90),((4, 4), 0.00)]);
        assert_almost_eq(set_eps(None, &distances), 0.02, 1e-8);

        let distances = HashMap::from_iter(vec![((1, 1), 0.00),((1, 2), 0.01),((1, 3), 0.90),((1, 4), 0.01),
                                                ((2, 1), 0.02),((2, 2), 0.00),((2, 3), 0.80),((2, 4), 0.01),
                                                ((3, 1), 0.90),((3, 2), 0.80),((3, 3), 0.00),((3, 4), 0.90),
                                                ((4, 1), 0.02),((4, 2), 0.02),((4, 3), 0.90),((4, 4), 0.00)]);
        assert_almost_eq(set_eps(None, &distances), 0.03, 1e-8);

        let distances = HashMap::from_iter(vec![((1, 1), 0.00),((1, 2), 0.02),((1, 3), 0.90),((1, 4), 0.02),
                                                ((2, 1), 0.02),((2, 2), 0.00),((2, 3), 0.80),((2, 4), 0.01),
                                                ((3, 1), 0.90),((3, 2), 0.80),((3, 3), 0.00),((3, 4), 0.90),
                                                ((4, 1), 0.02),((4, 2), 0.01),((4, 3), 0.90),((4, 4), 0.00)]);
        assert_almost_eq(set_eps(None, &distances), 0.04, 1e-8);

    }

    #[test]
    fn test_set_minpts() {
        let sequences = vec![];
        assert_eq!(set_min_pts(Some(2), &sequences), 2);
        assert_eq!(set_min_pts(Some(3), &sequences), 3);
        assert_eq!(set_min_pts(Some(4), &sequences), 4);

        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_3.fasta".to_string(), "contig_1".to_string(), 1)];
        assert_eq!(set_min_pts(None, &sequences), 3);

        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_3.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(5, "A".to_string(), "assembly_4.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_4.fasta".to_string(), "contig_2".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_4.fasta".to_string(), "contig_3".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_5.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_6.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_7.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_8.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_9.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_10.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_11.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_12.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_13.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_14.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_15.fasta".to_string(), "contig_1".to_string(), 1),
                             Sequence::new_with_seq(6, "A".to_string(), "assembly_16.fasta".to_string(), "contig_1".to_string(), 1)];
        assert_eq!(set_min_pts(None, &sequences), 4);
    }

    #[test]
    fn test_cluster_min_max_lengths() {
        let (cluster_min, cluster_max) = cluster_min_max_lengths(1000, 0.1);
        assert_eq!(cluster_min, 900); assert_eq!(cluster_max, 1100);

        let (cluster_min, cluster_max) = cluster_min_max_lengths(10000, 0.1);
        assert_eq!(cluster_min, 9000); assert_eq!(cluster_max, 11000);

        let (cluster_min, cluster_max) = cluster_min_max_lengths(1000, 0.01);
        assert_eq!(cluster_min, 990); assert_eq!(cluster_max, 1010);

        let (cluster_min, cluster_max) = cluster_min_max_lengths(1000, 0.5);
        assert_eq!(cluster_min, 500); assert_eq!(cluster_max, 1500);

        let (cluster_min, cluster_max) = cluster_min_max_lengths(1000, 0.9);
        assert_eq!(cluster_min, 100); assert_eq!(cluster_max, 1900);
    }
}
