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

use colored::Colorize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists, format_float, median_usize,
                  quit_with_error, usize_division_rounded, create_dir, delete_dir_if_exists};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn cluster(out_dir: PathBuf, cutoff: f64, min_assemblies_option: Option<usize>) {
    let gfa = out_dir.join("01_unitig_graph.gfa");
    let clustering_dir = out_dir.join("02_clustering");
    delete_dir_if_exists(&clustering_dir);
    create_dir(&clustering_dir);
    let pairwise_phylip = clustering_dir.join("pairwise_distances.phylip");
    let clustering_newick = clustering_dir.join("clustering.newick");
    // let clustering_png = clustering_dir.join("clustering.png");
    check_settings(&out_dir, &gfa, cutoff, &min_assemblies_option);
    starting_message();
    let (unitig_graph, mut sequences) = load_graph(&gfa, true);
    let assembly_count = get_assembly_count(&sequences);
    let min_assemblies = set_min_assemblies(min_assemblies_option, assembly_count);
    print_settings(&out_dir, cutoff, min_assemblies, min_assemblies_option);
    let asymmetrical_distances = pairwise_contig_distances(&unitig_graph, &sequences, &pairwise_phylip);
    let symmetrical_distances = make_symmetrical_distances(&asymmetrical_distances, &sequences);
    let mut tree = upgma_clustering(&symmetrical_distances, &mut sequences);
    normalise_tree(&mut tree);
    save_tree_to_newick(&tree, &sequences, &clustering_newick);
    define_clusters_from_tree(&tree, &mut sequences, cutoff);
    reorder_clusters(&mut sequences);
    let cluster_qc_results = cluster_qc(&sequences, &asymmetrical_distances, cutoff, min_assemblies);
    save_clusters(&sequences, &cluster_qc_results, &clustering_dir, &gfa);

    // TODO: save a PNG (or maybe SVG/PDF?) drawing of the tree to file

    finished_message(&pairwise_phylip, &clustering_newick);
}


fn check_settings(out_dir: &PathBuf, gfa: &PathBuf, cutoff: f64, min_assemblies: &Option<usize>) {
    check_if_dir_exists(&out_dir);
    check_if_file_exists(&gfa);
    if cutoff <= 0.0 || cutoff >= 1.0 {
        quit_with_error("--min_overlap_id must be between 0 and 1 (exclusive)");
    }
    if min_assemblies.is_some() && min_assemblies.unwrap() < 1 {
        quit_with_error("--min_assemblies must be 1 or greater");
    }
}


fn starting_message() {
    section_header("Starting autocycler cluster");
    explanation("This command takes a compacted De Bruijn graph (made by autocycler compress), \
                 trims any start-end overlap from the sequences and then clusters the sequences \
                 based on their similarity. Ideally, each cluster will then contain sequences \
                 which can be combined into a consensus.");
}


fn print_settings(out_dir: &PathBuf, cutoff: f64, min_assemblies: usize,
                  min_assemblies_option: Option<usize>) {
    eprintln!("Settings:");
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!("  --cutoff {}", format_float(cutoff));
    if min_assemblies_option.is_none() {
        eprintln!("  --min_assemblies {} (automatically set)", min_assemblies);
    } else {
        eprintln!("  --min_assemblies {}", min_assemblies);
    }
    eprintln!();
}


pub fn load_graph(gfa: &PathBuf, print_info: bool) -> (UnitigGraph, Vec<Sequence>) {
    if print_info {
        section_header("Loading graph");
        explanation("The compressed sequence graph is now loaded into memory.");
    }
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(&gfa);
    if print_info {
        unitig_graph.print_basic_graph_info();
    }
    (unitig_graph, sequences)
}


fn pairwise_contig_distances(unitig_graph: &UnitigGraph, sequences: &Vec<Sequence>, file_path: &PathBuf) -> HashMap<(u16, u16), f64> {
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
    save_distance_matrix(&distances, &sequences, file_path);
    distances
}


fn total_unitig_length(unitigs: &HashSet<u32>, unitig_lengths: &HashMap<u32, u32>) -> u32 {
    unitigs.iter().map(|u| unitig_lengths.get(u).unwrap()).sum()
}


fn save_distance_matrix(distances: &HashMap<(u16, u16), f64>, sequences: &Vec<Sequence>, file_path: &PathBuf) {
    eprintln!("Saving distance matrix:");
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
    eprintln!();
}


fn make_symmetrical_distances(asymmetrical_distances: &HashMap<(u16, u16), f64>,
                              sequences: &Vec<Sequence>) -> HashMap<(u16, u16), f64> {
    // This function takes in an asymmetrical distance matrix (where A vs B is not necessarily
    // equal to B vs A) and makes it symmetric by setting each distance (both orders) to the
    // maximum distance of the two orders.
    let mut symmetrical_distances: HashMap<(u16, u16), f64> = HashMap::new();
    for seq_a in sequences {
        for seq_b in sequences {
            let a_vs_b = asymmetrical_distances.get(&(seq_a.id, seq_b.id)).unwrap();
            let b_vs_a = asymmetrical_distances.get(&(seq_b.id, seq_a.id)).unwrap();
            let distance = a_vs_b.max(*b_vs_a);
            symmetrical_distances.insert((seq_a.id, seq_b.id), distance);
        }
    }
    symmetrical_distances
}


#[derive(Debug)]
struct TreeNode {
    id: u16,
    left: Option<Box<TreeNode>>,
    right: Option<Box<TreeNode>>,
    distance: f64,  // distance from this node to the tree tips
}


fn save_tree_to_newick(root: &TreeNode, sequences: &Vec<Sequence>, file_path: &PathBuf) {
    // Saves the tree to a NEWICK file. If necessary, it will add an additional node to specify the
    // length of the root, in order to ensure that root-to-tip distances are 0.5.
    eprintln!("Saving clustering tree:");
    let index: HashMap<u16, &Sequence> = sequences.iter().map(|s| (s.id, s)).collect();
    let newick_string = tree_to_newick(root, &index);
    let mut file = File::create(file_path).unwrap();
    if root.distance < 0.5 {
        let root_length = 0.5 - root.distance;
        write!(file, "({}:{});\n", newick_string, root_length).unwrap();
    } else {
        write!(file, "{};\n", newick_string).unwrap();
    }
    eprintln!("  {}", file_path.display());
    eprintln!();
}


fn tree_to_newick(node: &TreeNode, index: &HashMap<u16, &Sequence>) -> String {
    match (&node.left, &node.right) {
        (Some(left), Some(right)) => {
            let left_str = tree_to_newick(left, &index);
            let right_str = tree_to_newick(right, &index);
            format!("({}:{},{}:{})", left_str, node.distance - left.distance,
                                     right_str, node.distance - right.distance)
        }
        _ => format!("{}", index.get(&node.id).unwrap().string_for_newick()),
    }
}


fn upgma_clustering(distances: &HashMap<(u16, u16), f64>, sequences: &mut Vec<Sequence>) -> TreeNode {
    section_header("Clustering sequences");
    explanation("Contigs are organise into a tree using UPGMA. Then clusters are defined from the \
                 tree using the distance cutoff.");
    let mut clusters: HashMap<u16, HashSet<u16>> = HashMap::new();
    let mut cluster_distances: HashMap<(u16, u16), f64> = distances.clone();
    let mut nodes: HashMap<u16, TreeNode> = HashMap::new();


    // Initialise each sequence as its own cluster and create initial nodes.
    for seq in sequences.iter() {
        clusters.insert(seq.id, HashSet::from([seq.id]));
        nodes.insert(seq.id, TreeNode {
            id: seq.id,
            left: None,
            right: None,
            distance: 0.0,
        });
    }

    while clusters.len() > 1 {
        // Find the closest pair of clusters.
        let (a, b, a_b_distance) = get_closest_pair(&cluster_distances);
        let cluster_a = clusters.remove(&a).unwrap();
        let cluster_b = clusters.remove(&b).unwrap();
        let new_id = a.min(b);
        let mut new_cluster = HashSet::new();
        new_cluster.extend(cluster_a.iter());
        new_cluster.extend(cluster_b.iter());
        clusters.insert(new_id, new_cluster.clone());

        // Create a new tree node for the merged cluster.
        let new_node = TreeNode {
            id: new_id,
            left: Some(Box::new(nodes.remove(&a).unwrap())),
            right: Some(Box::new(nodes.remove(&b).unwrap())),
            distance: a_b_distance / 2.0,
        };
        nodes.insert(new_id, new_node);

        // Update distances between the new cluster and remaining clusters.
        let mut new_distances = HashMap::new();
        for (&(a, b), &dist) in cluster_distances.iter() {
            if clusters.contains_key(&a) && clusters.contains_key(&b) {
                new_distances.insert((a, b), dist);
            }
        }
        for &other_id in clusters.keys() {
            if other_id != new_id {
                let mut avg_dist = 0.0;
                let mut count = 0;
                for &id1 in new_cluster.iter() {
                    for &id2 in clusters[&other_id].iter() {
                        let dist = *distances.get(&(id1, id2)).unwrap_or(&distances[&(id2, id1)]);
                        avg_dist += dist;
                        count += 1;
                    }
                }
                avg_dist /= count as f64;
                new_distances.insert((new_id, other_id), avg_dist);
                new_distances.insert((other_id, new_id), avg_dist);
            }
        }
        cluster_distances = new_distances;
    }

    nodes.into_iter().next().unwrap().1  // return the root of the tree
}


fn get_closest_pair(distances: &HashMap<(u16, u16), f64>) -> (u16, u16, f64) {
    let mut min_distance = f64::INFINITY;
    let mut closest_pair = (0, 0);

    let mut unique_keys: Vec<u16> = distances.keys().flat_map(|&(a, b)| vec![a, b]).collect();
    unique_keys.sort_unstable();
    unique_keys.dedup();

    for i in 0..unique_keys.len() {
        let a = unique_keys[i];
        for j in i + 1..unique_keys.len() {
            let b = unique_keys[j];
            if let Some(&dist) = distances.get(&(a, b)).or_else(|| distances.get(&(b, a))) {
                if dist < min_distance {
                    min_distance = dist;
                    closest_pair = (a, b);
                }
            }
        }
    }
    (closest_pair.0, closest_pair.1, min_distance)
}


fn normalise_tree(root: &mut TreeNode) {
    if root.distance > 0.5 {
        scale_node_distance(root, 0.5 / root.distance);
    }
}


fn scale_node_distance(node: &mut TreeNode, scaling_factor: f64) {
    node.distance *= scaling_factor;
    if let Some(left) = &mut node.left {
        scale_node_distance(left, scaling_factor);
    }
    if let Some(right) = &mut node.right {
        scale_node_distance(right, scaling_factor);
    }
}


fn define_clusters_from_tree(tree: &TreeNode, sequences: &mut Vec<Sequence>, cutoff: f64) {
    let mut seq_id_to_cluster = HashMap::new();
    let mut current_cluster = 0;
    define_clusters_from_node(tree, cutoff / 2.0, &mut seq_id_to_cluster, &mut current_cluster);
    for s in sequences.iter_mut() {
        s.cluster = *seq_id_to_cluster.get(&s.id).unwrap();
    }
    eprintln!("{} clusters", current_cluster);
    eprintln!();
}


fn define_clusters_from_node(node: &TreeNode, cutoff: f64,
                             seq_id_to_cluster: &mut HashMap<u16, u16>, current_cluster: &mut u16) {
    if node.left.is_none() || node.distance < cutoff {  // node is either a tip or a tight group
        *current_cluster += 1;
        define_clusters(node, seq_id_to_cluster, *current_cluster);
    } else {  // node is a loose group, so continue recursively
        define_clusters_from_node(node.left.as_ref().unwrap(), cutoff, seq_id_to_cluster, current_cluster);
        define_clusters_from_node(node.right.as_ref().unwrap(), cutoff, seq_id_to_cluster, current_cluster);
    }
}


fn define_clusters(node: &TreeNode, seq_id_to_cluster: &mut HashMap<u16, u16>, cluster_number: u16) {
    if node.left.is_none() {  // is a tip
        seq_id_to_cluster.insert(node.id, cluster_number);
    } else {
        define_clusters(node.left.as_ref().unwrap(), seq_id_to_cluster, cluster_number);
        define_clusters(node.right.as_ref().unwrap(), seq_id_to_cluster, cluster_number);
    }
}


fn set_min_assemblies(min_assemblies_option: Option<usize>, assembly_count: usize) -> usize {
    // This function automatically sets the --min_assemblies parameter, if the user didn't
    // explicitly supply one. The auto-set value will be one-quarter of the assembly count (rounded)
    // but no less than 2.
    if min_assemblies_option.is_some() {
        return min_assemblies_option.unwrap();
    }
    let mut min_assemblies = usize_division_rounded(assembly_count, 4);
    if min_assemblies < 2 {
        min_assemblies = 2;
    }
    min_assemblies
}


fn cluster_qc(sequences: &Vec<Sequence>, distances: &HashMap<(u16, u16), f64>, cutoff: f64,
              min_assemblies: usize) -> HashMap<u16, Vec<String>> {
    // This function chooses which clusters pass and which clusters fail. It returns the result as
    // a vector of failure reasons, with an empty vector meaning the cluster passed QC.
    section_header("Cluster QC");
    explanation("Clusters are rejected if they occur in too few assemblies (set by \
                 --min_assemblies) or if they appear to be contained within a larger better \
                 cluster.");
    for s in sequences {
        assert!(s.cluster != 0);
    }
    let max_cluster = get_max_cluster(&sequences) + 1;
    let mut failure_reasons: HashMap<u16, Vec<String>> = (1..=max_cluster).map(|c| (c, Vec::new())).collect();

    for c in 1..max_cluster {
        let assemblies: HashSet<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.filename.clone()).collect();
        if assemblies.len() < min_assemblies {
            failure_reasons.get_mut(&c).unwrap().push("present in too few assemblies".to_string());
        }
    }

    for c in 1..max_cluster {
        let container = cluster_is_contained_in_another(c, sequences, distances, cutoff, &failure_reasons);
        if container > 0 {
            failure_reasons.get_mut(&c).unwrap().push(format!("contained within cluster {}", container));
        }
    }

    failure_reasons
}


fn cluster_is_contained_in_another(cluster_num: u16, sequences: &Vec<Sequence>,
                                   distances: &HashMap<(u16, u16), f64>, cutoff: f64,
                                   failure_reasons: &HashMap<u16, Vec<String>>) -> u16 {
    // Checks whether this cluster is contained within another cluster that has so-far passed QC.
    // If so, it returns the id of the containing cluster. If not, it returns 0.
    let passed_clusters: Vec<u16> = failure_reasons.iter().filter(|(_, v)| v.is_empty()).map(|(&k, _)| k).collect();
    for seq_a in sequences.iter().filter(|s| s.cluster == cluster_num) {
        for seq_b in sequences.iter().filter(|s| s.cluster != cluster_num && passed_clusters.contains(&s.cluster)) {
            let distance_a_b = distances.get(&(seq_a.id, seq_b.id)).unwrap();
            let distance_b_a = distances.get(&(seq_b.id, seq_a.id)).unwrap();
            if *distance_a_b < *distance_b_a && *distance_a_b < cutoff {
                return seq_b.cluster;
            }
        }
    }
    0
}


fn save_clusters(sequences: &Vec<Sequence>, qc_results: &HashMap<u16, Vec<String>>,
                 clustering_dir: &PathBuf, gfa_path: &PathBuf) {
    let pass_dir = clustering_dir.join("qc_pass");
    let fail_dir = clustering_dir.join("qc_fail");
    let all_seq_ids: Vec<_> = sequences.iter().map(|s| s.id).collect();
    save_qc_pass_clusters(sequences, &all_seq_ids, qc_results, gfa_path, &pass_dir);
    save_qc_fail_clusters(sequences, &all_seq_ids, qc_results, gfa_path, &fail_dir);
}


fn save_qc_pass_clusters(sequences: &Vec<Sequence>, all_seq_ids: &Vec<u16>, qc_results: &HashMap<u16, Vec<String>>,
                         gfa_path: &PathBuf, pass_dir: &PathBuf) {
    for c in 1..get_max_cluster(sequences)+1 {
        let failure_reasons = qc_results.get(&c).unwrap();
        if failure_reasons.len() == 0 {
            eprintln!("Cluster {:03}:", c);
            for s in sequences.iter().filter(|s| s.cluster == c) {
                eprintln!("  {}", s);
            }
            eprintln!("{}", "  passed QC".green());
            eprintln!();
            let cluster_dir = pass_dir.join(format!("cluster_{:03}", c));
            create_dir(&cluster_dir);
            save_cluster_gfa(&sequences, all_seq_ids, c, gfa_path, cluster_dir.join("1_untrimmed.gfa"));
        }
    }
}


fn save_qc_fail_clusters(sequences: &Vec<Sequence>, all_seq_ids: &Vec<u16>, qc_results: &HashMap<u16, Vec<String>>,
                         gfa_path: &PathBuf, fail_dir: &PathBuf) {
    for c in 1..get_max_cluster(sequences)+1 {
        let failure_reasons = qc_results.get(&c).unwrap();
        if failure_reasons.len() > 0 {
            eprintln!("Cluster {:03}:", c);
            for s in sequences.iter().filter(|s| s.cluster == c) {
                eprintln!("  {}", s.to_string().dimmed());
            }
            for f in failure_reasons {
                eprintln!("  {}", format!("failed QC: {}", f).red());
            }
            let cluster_dir = fail_dir.join(format!("cluster_{:03}", c));
            create_dir(&cluster_dir);
            save_cluster_gfa(&sequences, all_seq_ids, c, gfa_path, cluster_dir.join("1_untrimmed.gfa"));
            eprintln!();
        }
    }
}


fn save_cluster_gfa(sequences: &Vec<Sequence>, all_seq_ids: &Vec<u16>, cluster_num: u16,
                    in_gfa: &PathBuf, out_gfa: PathBuf) {
    let cluster_seqs: Vec<Sequence> = sequences.iter().filter(|s| s.cluster == cluster_num).cloned().collect();
    let (mut cluster_graph, _) = load_graph(&in_gfa, false);
    let seq_ids_to_remove:Vec<_> = sequences.iter().filter(|s| s.cluster != cluster_num).map(|s| s.id).collect();
    for id in seq_ids_to_remove {
        cluster_graph.remove_sequence_from_graph(id);
    }
    cluster_graph.remove_zero_depth_unitigs();
    merge_linear_paths(&mut cluster_graph, &cluster_seqs);
    cluster_graph.save_gfa(&out_gfa, &cluster_seqs).unwrap();
}


fn reorder_clusters(sequences: &mut Vec<Sequence>) {
    // Reorder clusters based on their median sequence length (large to small).
    // TODO: perhaps sort not based on sequence length but on UNIQUE sequence length (the sum of 
    //       the lengths of all unitigs). This will prevent doubled plasmids from being longer than
    //       they should be.
    let mut cluster_lengths = HashMap::new();
    for c in 1..get_max_cluster(sequences)+1 {
        let lengths: Vec<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.length).collect();
        cluster_lengths.insert(c, median_usize(&lengths));
    }
    let mut sorted_cluster_lengths: Vec<_> = cluster_lengths.iter().collect();
    sorted_cluster_lengths.sort_by(|a, b| {match b.1.cmp(a.1) {std::cmp::Ordering::Equal => a.0.cmp(b.0), other => other,}});
    let mut old_to_new = HashMap::new();
    let mut new_c: u16 = 0;
    for (old_c, _length) in sorted_cluster_lengths {
        new_c += 1;
        old_to_new.insert(old_c, new_c);
    }
    for s in sequences {
        if s.cluster < 1 {continue;}
        s.cluster = *old_to_new.get(&s.cluster).unwrap();
    }
}


fn get_assembly_count(sequences: &Vec<Sequence>) -> usize {
    sequences.iter().map(|s| &s.filename).collect::<HashSet<_>>().len()
}


fn get_max_cluster(sequences: &Vec<Sequence>) -> u16 {
    sequences.iter().map(|s| s.cluster).max().unwrap()
}


fn finished_message(pairwise_phylip: &PathBuf, clustering_newick: &PathBuf) {
    section_header("Finished!");
    explanation("You can now run autocycler trim on each cluster.");
    eprintln!("Pairwise distances: {}", pairwise_phylip.display());
    eprintln!("Clustering tree:    {}", clustering_newick.display());
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
    fn test_upgma_clustering_1() {
        // Uses example from https://en.wikipedia.org/wiki/UPGMA
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "a".to_string(), "a".to_string(), 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "b".to_string(), "b".to_string(), 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "c".to_string(), "c".to_string(), 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "d".to_string(), "d".to_string(), 1),
                                 Sequence::new_with_seq(5, "A".to_string(), "e".to_string(), "e".to_string(), 1)];
        let distances = HashMap::from_iter(vec![((1, 1), 00.0), ((1, 2), 17.0), ((1, 3), 21.0), ((1, 4), 31.0), ((1, 5), 23.0),
                                                ((2, 1), 17.0), ((2, 2), 00.0), ((2, 3), 30.0), ((2, 4), 34.0), ((2, 5), 21.0),
                                                ((3, 1), 21.0), ((3, 2), 30.0), ((3, 3), 00.0), ((3, 4), 28.0), ((3, 5), 39.0),
                                                ((4, 1), 31.0), ((4, 2), 34.0), ((4, 3), 28.0), ((4, 4), 00.0), ((4, 5), 43.0),
                                                ((5, 1), 23.0), ((5, 2), 21.0), ((5, 3), 39.0), ((5, 4), 43.0), ((5, 5), 00.0)]);
        let root = upgma_clustering(&distances, &mut sequences);
        assert_almost_eq(root.distance, 16.5, 1e-8);

        let index: HashMap<u16, &Sequence> = sequences.iter().map(|s| (s.id, s)).collect();
        let newick_string = tree_to_newick(&root, &index);
        assert_eq!(newick_string, "(((a__a__1_bp:8.5,b__b__1_bp:8.5):2.5,e__e__1_bp:11):5.5,(c__c__1_bp:14,d__d__1_bp:14):2.5)");
    }

    #[test]
    fn test_upgma_clustering_2() {
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "a".to_string(), "a".to_string(), 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "b".to_string(), "b".to_string(), 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "c".to_string(), "c".to_string(), 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "d".to_string(), "d".to_string(), 1)];
        let distances = HashMap::from_iter(vec![((1, 1), 0.0), ((1, 2), 0.1), ((1, 3), 0.5), ((1, 4), 0.5),
                                                ((2, 1), 0.1), ((2, 2), 0.0), ((2, 3), 0.5), ((2, 4), 0.5),
                                                ((3, 1), 0.5), ((3, 2), 0.5), ((3, 3), 0.0), ((3, 4), 0.2),
                                                ((4, 1), 0.5), ((4, 2), 0.5), ((4, 3), 0.2), ((4, 4), 0.0)]);
        let mut root = upgma_clustering(&distances, &mut sequences);
        normalise_tree(&mut root);
        assert_almost_eq(root.distance, 0.25, 1e-8);

        let index: HashMap<u16, &Sequence> = sequences.iter().map(|s| (s.id, s)).collect();
        let newick_string = tree_to_newick(&root, &index);
        assert_eq!(newick_string, "((a__a__1_bp:0.05,b__b__1_bp:0.05):0.2,(c__c__1_bp:0.1,d__d__1_bp:0.1):0.15)");
    }

    // #[test]
    // fn test_reorder_clusters() {
    //     let mut sequences = vec![Sequence::new_with_seq(1, "CGCGA".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 5),
    //                              Sequence::new_with_seq(2, "T".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1),
    //                              Sequence::new_with_seq(3, "AACGACTACG".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 10),
    //                              Sequence::new_with_seq(4, "CGCGA".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 5),
    //                              Sequence::new_with_seq(5, "T".to_string(), "assembly_2.fasta".to_string(), "contig_3".to_string(), 1),
    //                              Sequence::new_with_seq(6, "AACGACTACG".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 10)];
    //     sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3;
    //     sequences[3].cluster = 1; sequences[4].cluster = 2; sequences[5].cluster = 3;
    //     reorder_clusters(&mut sequences);
    //     assert_eq!(sequences[0].cluster, 2); assert_eq!(sequences[1].cluster, 3); assert_eq!(sequences[2].cluster, 1);
    //     assert_eq!(sequences[3].cluster, 2); assert_eq!(sequences[4].cluster, 3); assert_eq!(sequences[5].cluster, 1);
    //     reorder_clusters(&mut sequences);
    //     assert_eq!(sequences[0].cluster, 2); assert_eq!(sequences[1].cluster, 3); assert_eq!(sequences[2].cluster, 1);
    //     assert_eq!(sequences[3].cluster, 2); assert_eq!(sequences[4].cluster, 3); assert_eq!(sequences[5].cluster, 1);
    // }

    #[test]
    fn test_set_minpts() {
        assert_eq!(set_min_assemblies(Some(2), 10), 2);
        assert_eq!(set_min_assemblies(Some(11), 10), 11);
        assert_eq!(set_min_assemblies(Some(321), 10), 321);
        assert_eq!(set_min_assemblies(None, 4), 2);
        assert_eq!(set_min_assemblies(None, 8), 2);
        assert_eq!(set_min_assemblies(None, 12), 3);
        assert_eq!(set_min_assemblies(None, 13), 3);
        assert_eq!(set_min_assemblies(None, 15), 4);
        assert_eq!(set_min_assemblies(None, 16), 4);
    }
}
