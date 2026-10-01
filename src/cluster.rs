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
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::metrics::{ClusteringMetrics, UntrimmedClusterMetrics};
use crate::misc::{check_if_dir_exists, check_if_file_exists, format_float, median,
                  quit_with_error, usize_division_rounded, create_dir, delete_dir_if_exists,
                  load_file_lines, parse_node_numbers};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn cluster(autocycler_dir: PathBuf, cutoff: f64, min_assemblies_option: Option<usize>,
               max_contigs: u32, manual_clusters: Option<String>) {
    let gfa = autocycler_dir.join("input_assemblies.gfa");
    let clustering_dir = autocycler_dir.join("clustering");
    let pairwise_phylip = clustering_dir.join("pairwise_distances.phylip");
    let clustering_newick = clustering_dir.join("clustering.newick");
    let clustering_tsv = clustering_dir.join("clustering.tsv");
    let clustering_yaml = clustering_dir.join("clustering.yaml");
    check_settings(&autocycler_dir, &gfa, cutoff, &min_assemblies_option);
    delete_dir_if_exists(&clustering_dir);
    create_dir(&clustering_dir);
    starting_message();
    let gfa_lines = load_file_lines(&gfa);
    let (graph, mut sequences) = UnitigGraph::from_gfa_lines(&gfa_lines);
    let min_assemblies = set_min_assemblies(min_assemblies_option, &sequences);
    let manual_clusters = parse_node_numbers(manual_clusters);
    print_settings(&autocycler_dir, cutoff, min_assemblies, min_assemblies_option, max_contigs,
                   &manual_clusters);
    check_sequence_count(&sequences, max_contigs);
    let asymmetrical_distances = pairwise_contig_distances(&graph, &sequences, &pairwise_phylip);
    let symmetrical_distances = make_symmetrical_distances(&asymmetrical_distances, &sequences);
    let mut tree = upgma(&symmetrical_distances, &sequences);
    normalise_tree(&mut tree);
    save_tree_to_newick(&tree, &sequences, &clustering_newick);
    let qc_results = generate_clusters(&tree, &mut sequences, &asymmetrical_distances, cutoff,
                                       min_assemblies, &manual_clusters);
    save_clusters(&sequences, &qc_results, &clustering_dir, &gfa_lines);
    save_data_to_tsv(&sequences, &qc_results, &clustering_tsv);
    let metrics = clustering_metrics(&sequences, &qc_results);
    metrics.save_to_yaml(&clustering_yaml);

    finished_message(&pairwise_phylip, &clustering_newick, &clustering_tsv);
}


fn check_settings(autocycler_dir: &Path, gfa: &Path, cutoff: f64, min_assemblies: &Option<usize>) {
    check_if_dir_exists(autocycler_dir);
    check_if_file_exists(gfa);
    if cutoff <= 0.0 || cutoff >= 1.0 {
        quit_with_error("--cutoff must be between 0 and 1 (exclusive)");
    }
    if *min_assemblies == Some(0) {
        quit_with_error("--min_assemblies must be 1 or greater");
    }
}


fn starting_message() {
    section_header("Starting autocycler cluster");
    explanation("This command takes a unitig graph (made by autocycler compress) and clusters the \
                 sequences based on their similarity. Ideally, each cluster will then contain \
                 sequences which can be combined into a consensus.");
}


fn finished_message(pairwise_phylip: &Path, clustering_newick: &Path, clustering_tsv: &Path) {
    section_header("Finished!");
    explanation("You can now run autocycler trim on each cluster. If you want to manually \
                 inspect the clustering, you can view the following files.");
    eprintln!("Pairwise distances:         {}", pairwise_phylip.display());
    eprintln!("Clustering tree (Newick):   {}", clustering_newick.display());
    eprintln!("Clustering tree (metadata): {}", clustering_tsv.display());
    eprintln!();
}


fn print_settings(autocycler_dir: &Path, cutoff: f64, min_assemblies: usize,
                  min_assemblies_option: Option<usize>, max_contigs: u32, manual_clusters: &[u16]) {
    eprintln!("Settings:");
    eprintln!("  --autocycler_dir {}", autocycler_dir.display());
    eprintln!("  --cutoff {}", format_float(cutoff));
    if min_assemblies_option.is_none() {
        eprintln!("  --min_assemblies {min_assemblies} (automatically set)");
    } else {
        eprintln!("  --min_assemblies {min_assemblies}");
    }
    eprintln!("  --max_contigs {max_contigs}");
    if !manual_clusters.is_empty() {
        eprintln!("  --manual {}", manual_clusters.iter().map(|c| c.to_string())
                                                  .collect::<Vec<String>>().join(","));
    }
    eprintln!();
}


fn check_sequence_count(sequences: &[Sequence], max_contigs: u32) {
    let assembly_count = get_assembly_count(sequences) as f64;
    let sequence_count = sequences.len() as f64;
    if sequence_count == 0.0 {
        quit_with_error("no sequences found in input_assemblies.gfa")
    }
    let mean_seqs_per_assembly = sequence_count / assembly_count;
    if mean_seqs_per_assembly > max_contigs as f64 {
        let e = format!("the mean number of contigs per input assembly ({mean_seqs_per_assembly:.1}) exceeds the allowed \
                         threshold ({max_contigs}). Are your input assemblies fragmented or contaminated?");
        quit_with_error(&e);
    }
}


fn pairwise_contig_distances(graph: &UnitigGraph, sequences: &[Sequence], file_path: &Path)
        -> HashMap<(u16, u16), f64> {
    section_header("Pairwise distances");
    explanation("Every pairwise distance between contigs is calculated based on the similarity of \
                 their paths through the graph.");
    let unitig_lengths: HashMap<u32, u32> = graph.unitigs.iter()
        .map(|rc| {let u = rc.borrow(); (u.number, u.length())}).collect();
    let sequence_unitigs: HashMap<u16, HashSet<u32>> = sequences.iter()
        .map(|s| (s.id, graph.get_unitig_path_for_sequence(s).iter()
        .map(|(number, _)| *number).collect::<HashSet<u32>>())).collect();
    let mut distances: HashMap<(u16, u16), f64> = HashMap::new();
    for seq_a in sequences {
        let a = sequence_unitigs.get(&seq_a.id).unwrap();
        let a_len = f64::from(a.iter().map(|u| unitig_lengths[u]).sum::<u32>());
        for seq_b in sequences {
            let b = sequence_unitigs.get(&seq_b.id).unwrap();
            let ab_len: f64 = a.intersection(b).map(|u| f64::from(unitig_lengths[u])).sum();
            let distance = 1.0 - (ab_len / a_len);
            distances.insert((seq_a.id, seq_b.id), distance);
        }
    }
    eprintln!("{} sequences, {} total pairwise distances", sequences.len(), distances.len());
    eprintln!();
    save_distance_matrix(&distances, sequences, file_path);
    distances
}


fn save_distance_matrix(distances: &HashMap<(u16, u16), f64>, sequences: &[Sequence],
                        file_path: &Path) {
    eprintln!("Saving distance matrix:");
    let mut f = BufWriter::new(File::create(file_path).unwrap());
    writeln!(f, "{}", sequences.len()).unwrap();
    for seq_a in sequences {
        write!(f, "{seq_a}").unwrap();
        for seq_b in sequences {
            write!(f, "\t{:.8}", distances.get(&(seq_a.id, seq_b.id)).unwrap()).unwrap();
        }
        writeln!(f).unwrap();
    }
    f.flush().unwrap();
    eprintln!("  {}", file_path.display());
    eprintln!();
}


fn make_symmetrical_distances(asymmetrical_distances: &HashMap<(u16, u16), f64>,
                              sequences: &[Sequence]) -> HashMap<(u16, u16), f64> {
    // Use the larger distance for both directions.
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


#[derive(Debug, Default)]
struct TreeNode {
    id: u16,
    left: Option<Box<TreeNode>>,
    right: Option<Box<TreeNode>>,
    distance: f64,  // distance from this node to the tree tips
}

impl TreeNode {
    fn is_tip(&self) -> bool {
        self.left.is_none()
    }

    fn max_pairwise_distance(&self, node_num: u16) -> f64 {
        self.find_node(node_num).map_or(-1.0, |node| node.distance * 2.0)
    }

    fn automatic_clustering(&self, cutoff: f64) -> Vec<u16> {
        self.manual_clustering(cutoff, &[])
    }

    fn manual_clustering(&self, cutoff: f64, manual_clusters: &[u16]) -> Vec<u16> {
        for &id in manual_clusters {
            if self.find_node(id).is_none() {
                quit_with_error(&format!("clustering tree does not contain a node with id {id}"));
            }
        }
        if !manual_clusters.is_empty() {
            self.check_consistency(manual_clusters);
        }
        let mut clusters = Vec::new();
        self.collect_clusters(cutoff / 2.0, manual_clusters, &mut clusters);
        clusters.sort();
        clusters
    }

    fn collect_clusters(&self, cutoff: f64, manual_clusters: &[u16], clusters: &mut Vec<u16>) {
        if manual_clusters.contains(&self.id) ||
                (self.distance <= cutoff && !self.contains_manual_cluster(manual_clusters)) {
            clusters.push(self.id);
        } else if !self.is_tip() {
            self.left.as_ref().unwrap().collect_clusters(cutoff, manual_clusters, clusters);
            self.right.as_ref().unwrap().collect_clusters(cutoff, manual_clusters, clusters);
        }
    }

    fn contains_manual_cluster(&self, manual_clusters: &[u16]) -> bool {
        manual_clusters.contains(&self.id) || (!self.is_tip() &&
            (self.left.as_ref().unwrap().contains_manual_cluster(manual_clusters) ||
             self.right.as_ref().unwrap().contains_manual_cluster(manual_clusters)))
    }

    fn check_consistency(&self, manual_clusters: &[u16]) {
        // Ensures that no manual cluster is contained within another manual cluster.
        if !self.is_tip() {
            if manual_clusters.contains(&self.id) &&
                    (self.left.as_ref().unwrap().contains_manual_cluster(manual_clusters) ||
                     self.right.as_ref().unwrap().contains_manual_cluster(manual_clusters)) {
                quit_with_error("manual clusters cannot be nested");
            }
            self.left.as_ref().unwrap().check_consistency(manual_clusters);
            self.right.as_ref().unwrap().check_consistency(manual_clusters);
        }
    }

    fn get_tips(&self, node_num: u16) -> Vec<u16> {
        let mut tips = Vec::new();
        if let Some(node) = self.find_node(node_num) {
            node.collect_tips(&mut tips);
        }
        tips
    }

    fn collect_tips(&self, tips: &mut Vec<u16>) {
        if self.is_tip() {
            tips.push(self.id);
        } else {
            self.left.as_ref().unwrap().collect_tips(tips);
            self.right.as_ref().unwrap().collect_tips(tips);
        }
    }

    fn check_complete_coverage(&self, clusters: &[u16]) {
        // Every tip must belong to exactly one cluster.
        let all_tips: HashSet<u16> = self.get_tips(self.id).into_iter().collect();
        let mut covered_tips = HashSet::new();
        for &c in clusters {
            for tip in self.get_tips(c) {
                if !covered_tips.insert(tip) { panic!("overlap detected"); }
            }
        }
        if covered_tips != all_tips { panic!("incomplete coverage");}
    }

    fn split_clusters(&self, clusters: &[u16]) -> Vec<Vec<u16>> {
        // Each alternative splits one cluster into its two children.
        self.check_complete_coverage(clusters);
        let mut result = Vec::new();

        for &cluster in clusters {
            let node = self.find_node(cluster).unwrap();
            if !node.is_tip() {
                let mut new_cluster: Vec<_> = clusters.iter().copied()
                    .filter(|&other| other != cluster).collect();
                new_cluster.push(node.left.as_ref().unwrap().id);
                new_cluster.push(node.right.as_ref().unwrap().id);
                new_cluster.sort();
                result.push(new_cluster);
            }
        }
        result.sort();
        result
    }

    fn find_node(&self, node_num: u16) -> Option<&TreeNode> {
        if self.id == node_num { return Some(self); }
        if self.is_tip() { return None; }
        self.left.as_ref().unwrap().find_node(node_num)
            .or_else(|| self.right.as_ref().unwrap().find_node(node_num))
    }
}


fn save_tree_to_newick(root: &TreeNode, sequences: &[Sequence], file_path: &Path) {
    // Saves the tree to a NEWICK file. If necessary, it will add an additional node to specify the
    // length of the root, in order to ensure that root-to-tip distances are 0.5.
    eprintln!("Saving clustering tree:");
    let index: HashMap<u16, &Sequence> = sequences.iter().map(|s| (s.id, s)).collect();
    let newick_string = tree_to_newick(root, &index);
    let mut file = File::create(file_path).unwrap();
    if root.distance < 0.5 {
        let root_length = 0.5 - root.distance;
        writeln!(file, "({newick_string}:{root_length});").unwrap();
    } else {
        writeln!(file, "{newick_string};").unwrap();
    }
    eprintln!("  {}", file_path.display());
    eprintln!();
}


fn tree_to_newick(node: &TreeNode, index: &HashMap<u16, &Sequence>) -> String {
    match (&node.left, &node.right) {
        (Some(left), Some(right)) => {
            let left_str = tree_to_newick(left, index);
            let right_str = tree_to_newick(right, index);
            format!("({}:{},{}:{}){}", left_str, node.distance - left.distance,
                                       right_str, node.distance - right.distance,
                                       node.id)
        }
        _ => index[&node.id].string_for_newick(),
    }
}


fn upgma(distances: &HashMap<(u16, u16), f64>, sequences: &[Sequence]) -> TreeNode {
    section_header("Clustering sequences");
    explanation("Contigs are organised into a tree using UPGMA. Then clusters are defined from the \
                 tree using the distance cutoff.");
    let mut clusters: HashMap<u16, HashSet<u16>> = HashMap::new();
    let mut cluster_distances: HashMap<(u16, u16), f64> = distances.clone();
    let mut nodes: HashMap<u16, TreeNode> = HashMap::new();
    let mut internal_node_num: u16 = sequences.iter().map(|s| s.id).max().unwrap();

    for seq in sequences {
        clusters.insert(seq.id, HashSet::from([seq.id]));
        nodes.insert(seq.id, TreeNode { id: seq.id, ..Default::default() });
    }

    while clusters.len() > 1 {
        let (a, b, a_b_distance) = get_closest_pair(&cluster_distances);
        let cluster_a = clusters.remove(&a).unwrap();
        let cluster_b = clusters.remove(&b).unwrap();
        let new_id = a.min(b);
        let mut new_cluster = HashSet::new();
        new_cluster.extend(cluster_a.iter());
        new_cluster.extend(cluster_b.iter());
        clusters.insert(new_id, new_cluster);

        internal_node_num += 1;
        let new_node = TreeNode {
            id: internal_node_num,
            left: Some(Box::new(nodes.remove(&a).unwrap())),
            right: Some(Box::new(nodes.remove(&b).unwrap())),
            distance: a_b_distance / 2.0,
        };
        nodes.insert(new_id, new_node);

        cluster_distances.retain(|(a, b), _| clusters.contains_key(a) && clusters.contains_key(b));
        for &other_id in clusters.keys() {
            if other_id != new_id {
                let distance = mean_cluster_distance(&clusters[&new_id], &clusters[&other_id], distances);
                cluster_distances.insert((new_id, other_id), distance);
                cluster_distances.insert((other_id, new_id), distance);
            }
        }
    }

    nodes.into_values().next().unwrap()
}


fn mean_cluster_distance(cluster_a: &HashSet<u16>, cluster_b: &HashSet<u16>,
                         distances: &HashMap<(u16, u16), f64>) -> f64 {
    let mut total = 0.0;
    let mut count = 0;
    for &a in cluster_a {
        for &b in cluster_b {
            total += distances.get(&(a, b)).unwrap_or(&distances[&(b, a)]);
            count += 1;
        }
    }
    total / count as f64
}


fn get_closest_pair(distances: &HashMap<(u16, u16), f64>) -> (u16, u16, f64) {
    let mut min_distance = f64::INFINITY;
    let mut closest_pair = (0, 0);

    let mut unique_keys: Vec<u16> = distances.keys().flat_map(|&(a, b)| [a, b]).collect();
    unique_keys.sort_unstable();
    unique_keys.dedup();

    for (i, &a) in unique_keys.iter().enumerate() {
        for &b in unique_keys.iter().skip(i + 1) {
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
    if let Some(left)  = &mut node.left  { scale_node_distance(left,  scaling_factor); }
    if let Some(right) = &mut node.right { scale_node_distance(right, scaling_factor); }
}


fn generate_clusters(tree: &TreeNode, sequences: &mut [Sequence],
                     distances: &HashMap<(u16, u16), f64>, cutoff: f64, min_assemblies: usize,
                     manual_clusters: &[u16]) -> HashMap<u16, ClusterQC> {
    let clusters = if manual_clusters.is_empty() {
        let auto_clusters = tree.automatic_clustering(cutoff);
        refine_auto_clusters(tree, sequences, distances, &auto_clusters, cutoff, min_assemblies)
    } else {
        tree.manual_clustering(cutoff, manual_clusters)
    };
    tree.check_complete_coverage(&clusters);
    qc_clusters(tree, sequences, distances, &clusters, manual_clusters, cutoff, min_assemblies)
}


fn qc_clusters(tree: &TreeNode, sequences: &mut [Sequence], distances: &HashMap<(u16, u16), f64>,
               cluster_nodes: &[u16], manual_clusters: &[u16], cutoff: f64,
               min_assemblies: usize) -> HashMap<u16, ClusterQC> {
    let mut qc_results = initialise_cluster_qc(tree, sequences, cluster_nodes, manual_clusters);
    if !manual_clusters.is_empty() { return qc_results; }

    // Finish the assembly-count checks before looking for passing clusters that contain others.
    let max_cluster = get_max_cluster(sequences);
    for c in 1..=max_cluster {
        if cluster_assembly_count(sequences, c) < min_assemblies && !cluster_is_trusted(sequences, c) {
            qc_results.get_mut(&c).unwrap().failure_reasons.push("present in too few assemblies".to_string());
        }
    }
    for c in 1..=max_cluster {
        let container = cluster_is_contained_in_another(c, sequences, distances, cutoff, &qc_results);
        if container > 0 && !cluster_is_trusted(sequences, c) {
            qc_results.get_mut(&c).unwrap().failure_reasons.push(format!("contained within cluster {container}"));
        }
    }
    qc_results
}


fn initialise_cluster_qc(tree: &TreeNode, sequences: &mut [Sequence], cluster_nodes: &[u16],
                          manual_clusters: &[u16]) -> HashMap<u16, ClusterQC> {
    let mut current_cluster = 0;
    let mut qc_results = HashMap::new();
    for n in cluster_nodes {
        if let Some(node) = tree.find_node(*n) {
            current_cluster += 1;
            assign_cluster_to_node(node, sequences, current_cluster);
            let mut qc = ClusterQC::new(tree.max_pairwise_distance(*n));
            if !manual_clusters.is_empty() && !manual_clusters.contains(n) {
                qc.failure_reasons.push("not included in manual clusters".to_string());
            }
            qc_results.insert(current_cluster, qc);
        } else {
            quit_with_error(&format!("clustering tree does not contain a node with id {n}"));
        }
    }

    let old_to_new = reorder_clusters(sequences);
    qc_results.into_iter().filter_map(|(old, qc)| old_to_new.get(&old).map(|&new| (new, qc))).collect()
}


fn cluster_assembly_count(sequences: &[Sequence], c: u16) -> usize {
    // Each assembly contributes the largest weight among its contigs in this cluster.
    let mut assembly_weights: HashMap<&str, usize> = HashMap::new();
    for seq in sequences.iter().filter(|s| s.cluster == c) {
        let weight = seq.cluster_weight();
        assembly_weights.entry(&seq.filename)
            .and_modify(|existing| *existing = (*existing).max(weight)).or_insert(weight);
    }
    assembly_weights.values().sum()
}


fn cluster_is_trusted(sequences: &[Sequence], c: u16) -> bool {
    sequences.iter().any(|s| s.cluster == c && s.is_trusted())
}


fn score_clustering(tree: &TreeNode, sequences: &mut [Sequence],
                    distances: &HashMap<(u16, u16), f64>, clusters: &[u16], cutoff: f64,
                    min_assemblies: usize) -> f64 {
    let qc_results = qc_clusters(tree, sequences, distances, clusters, &[], cutoff,
                                 min_assemblies);
    clustering_metrics(sequences, &qc_results).overall_clustering_score
}


fn refine_auto_clusters(tree: &TreeNode, sequences: &mut [Sequence],
                        distances: &HashMap<(u16, u16), f64>, clusters: &[u16], cutoff: f64,
                        min_assemblies: usize) -> Vec<u16> {
    // Keep splitting clusters while the score improves.
    let mut best_clusters = clusters.to_vec();
    let mut best_score = score_clustering(tree, sequences, distances, &best_clusters, cutoff,
                                          min_assemblies);
    let mut improved = true;
    while improved {
        improved = false;
        for alt_clusters in tree.split_clusters(&best_clusters) {
            let alt_score = score_clustering(tree, sequences, distances, &alt_clusters, cutoff,
                                             min_assemblies);
            if alt_score > best_score + 1e-12 {  // add a tolerance to avoid floating point issues
                best_clusters = alt_clusters;
                best_score = alt_score;
                improved = true;
            }
        }
    }
    best_clusters
}


fn assign_cluster_to_node(node: &TreeNode, sequences: &mut [Sequence], cluster: u16) {
    for s in sequences.iter_mut() {
        if s.id == node.id {
            s.cluster = cluster;
        }
    }
    if let Some(ref left)  = node.left  { assign_cluster_to_node(left,  sequences, cluster); }
    if let Some(ref right) = node.right { assign_cluster_to_node(right, sequences, cluster); }
}


fn set_min_assemblies(min_assemblies_option: Option<usize>, sequences: &[Sequence]) -> usize {
    // Default to one-quarter of the assemblies, with a minimum of two (one for a single assembly).
    if let Some(min_assemblies) = min_assemblies_option {
        return min_assemblies;
    }
    let assembly_count = get_assembly_count(sequences);
    if assembly_count == 1 {
        return 1;
    }
    usize_division_rounded(assembly_count, 4).max(2)
}


#[derive(Default)]
struct ClusterQC {
    pub failure_reasons: Vec<String>,
    pub cluster_dist: f64,
}

impl ClusterQC {
    pub fn new(cluster_dist: f64) -> Self {
        ClusterQC {
            failure_reasons: Vec::new(),
            cluster_dist,
        }
    }
    pub fn pass(&self) -> bool { self.failure_reasons.is_empty() }
}


fn cluster_is_contained_in_another(cluster_num: u16, sequences: &[Sequence],
                                   distances: &HashMap<(u16, u16), f64>, cutoff: f64,
                                   qc_results: &HashMap<u16, ClusterQC>) -> u16 {
    // Checks whether this cluster is contained within another cluster that has so-far passed QC.
    // If so, it returns the id of the containing cluster. If not, it returns 0.
    // A cluster counts as contained if the majority of the pairwise comparisons to another cluster
    // are asymmetrical and below the cutoff.
    let passed_clusters = qc_results.iter().filter(|(_, q)| q.pass()).map(|(&k, _)| k);
    for passed_cluster in passed_clusters {
        if passed_cluster == cluster_num {
            continue;
        }
        let mut contain_count = 0;
        let mut total_count = 0;
        for seq_a in sequences.iter().filter(|s| s.cluster == cluster_num) {
            for seq_b in sequences.iter().filter(|s| s.cluster == passed_cluster) {
                total_count += 1;
                let distance_a_b = distances.get(&(seq_a.id, seq_b.id)).unwrap();
                let distance_b_a = distances.get(&(seq_b.id, seq_a.id)).unwrap();
                if *distance_a_b < *distance_b_a && *distance_a_b < cutoff {
                    contain_count += 1;
                }
            }
        }
        let contained_fraction = contain_count as f64 / total_count as f64;
        if contained_fraction > 0.5 {
            return passed_cluster;
        }
    }
    0
}


fn save_clusters(sequences: &[Sequence], qc_results: &HashMap<u16, ClusterQC>,
                 clustering_dir: &Path, gfa_lines: &[String]) {
    for (passing, dirname) in [(true, "qc_pass"), (false, "qc_fail")] {
        let qc_dir = clustering_dir.join(dirname);
        for c in 1..=get_max_cluster(sequences) {
            let qc = &qc_results[&c];
            if qc.pass() == passing {
                save_cluster(sequences, c, qc, gfa_lines, &qc_dir.join(format!("cluster_{c:03}")));
            }
        }
    }
}


fn save_cluster(sequences: &[Sequence], cluster: u16, qc: &ClusterQC, gfa_lines: &[String],
                 cluster_dir: &Path) {
    eprintln!("Cluster {cluster:03}:");
    let mut seq_lengths = Vec::new();
    for seq in sequences.iter().filter(|s| s.cluster == cluster) {
        if qc.pass() {
            eprintln!("  {seq}");
        } else {
            eprintln!("  {}", seq.to_string().dimmed());
        }
        seq_lengths.push(seq.length);
    }
    if seq_lengths.len() > 1 {
        let message = format!("cluster distance: {}", format_float(qc.cluster_dist));
        if qc.pass() {
            eprintln!("  {message}");
        } else {
            eprintln!("  {}", message.dimmed());
        }
    }
    if qc.pass() {
        eprintln!("{}", "  passed QC".green());
    } else {
        for reason in &qc.failure_reasons {
            eprintln!("  {}", format!("failed QC: {reason}").red());
        }
    }
    create_dir(cluster_dir);
    save_cluster_gfa(sequences, cluster, gfa_lines, cluster_dir.join("1_untrimmed.gfa"));
    let metrics = UntrimmedClusterMetrics::new(seq_lengths, qc.cluster_dist);
    metrics.save_to_yaml(&cluster_dir.join("1_untrimmed.yaml"));
    eprintln!();
}


fn save_cluster_gfa(sequences: &[Sequence], cluster_num: u16, gfa_lines: &[String],
                    out_gfa: PathBuf) {
    let cluster_seqs: Vec<Sequence> = sequences.iter().filter(|s| s.cluster == cluster_num)
                                               .cloned().collect();
    let seq_ids_to_remove:Vec<_> = sequences.iter().filter(|s| s.cluster != cluster_num)
                                            .map(|s| s.id).collect();
    let filtered_gfa_lines: Vec<String> = filter_gfa_lines(gfa_lines, &seq_ids_to_remove);
    let (mut cluster_graph, _) = UnitigGraph::from_gfa_lines(&filtered_gfa_lines);
    cluster_graph.recalculate_depths();
    cluster_graph.remove_zero_depth_unitigs();
    merge_linear_paths(&mut cluster_graph, &cluster_seqs);
    cluster_graph.save_gfa(&out_gfa, &cluster_seqs, false).unwrap();
}


fn filter_gfa_lines(gfa_lines: &[String], paths_to_remove: &[u16]) -> Vec<String> {
    let paths_to_remove: HashSet<_> = paths_to_remove.iter().copied().collect();
    gfa_lines.iter().filter(|line| {
        let Some(rest) = line.strip_prefix("P\t") else { return true; };
        let path_name = rest.split('\t').next().unwrap_or_default();
        path_name.parse::<u16>().map_or(true, |id| !paths_to_remove.contains(&id))
    }).cloned().collect()
}


fn save_data_to_tsv(sequences: &[Sequence], qc_results: &HashMap<u16, ClusterQC>,
                    file_path: &Path) {
    let mut file = File::create(file_path).unwrap();
    writeln!(file, "node_name\tpassing_clusters\tall_clusters\tsequence_id\t\
                    file_name\tcontig_name\tlength\ttrusted\tcluster_weight\t\
                    consensus_weight").unwrap();
    for seq in sequences {
        assert!(seq.cluster != 0);
        let qc = qc_results.get(&seq.cluster).unwrap();
        let pass_cluster = if qc.pass() { seq.cluster.to_string() }
                                   else { "none".to_string() };
        writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                 seq.string_for_newick(), pass_cluster, seq.cluster, seq.id, seq.filename,
                 seq.contig_name(), seq.length, seq.is_trusted(), seq.cluster_weight(),
                 seq.consensus_weight()).unwrap();
    }
}


fn clustering_metrics(sequences: &[Sequence], qc_results: &HashMap<u16, ClusterQC>)
        -> ClusteringMetrics {
    let mut metrics = ClusteringMetrics::default();
    let mut cluster_filenames = HashMap::new();
    for seq in sequences {
        let qc = qc_results.get(&seq.cluster).unwrap();
        cluster_filenames.entry(seq.cluster).or_insert_with(Vec::new).push(seq.filename.clone());
        if qc.pass() {
            metrics.pass_contig_count += 1;
        } else {
            metrics.fail_contig_count += 1;
        }
    }
    let mut pass_cluster_stats = Vec::new();
    for c in 1..=get_max_cluster(sequences) {
        let qc = qc_results.get(&c).unwrap();
        if qc.pass() {
            metrics.pass_cluster_count += 1;
            let cluster_size = cluster_filenames.get(&c).map_or(0, Vec::len);
            pass_cluster_stats.push((qc.cluster_dist, cluster_size));
        } else {
            metrics.fail_cluster_count += 1;
        }
    }
    metrics.calculate_fractions();
    metrics.calculate_scores(cluster_filenames, pass_cluster_stats);
    metrics
}


fn reorder_clusters(sequences: &mut [Sequence]) -> HashMap<u16, u16> {
    // Reorder clusters based on their median sequence length (large to small). Returns the mapping
    // of old cluster numbers to new cluster numbers.
    let mut cluster_lengths: Vec<_> = (1..=get_max_cluster(sequences)).map(|c| {
        let lengths: Vec<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.length).collect();
        (c, median(&lengths))
    }).collect();
    cluster_lengths.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let old_to_new: HashMap<_, _> = cluster_lengths.iter().enumerate()
        .map(|(i, &(old, _))| (old, (i + 1) as u16)).collect();
    for seq in sequences.iter_mut().filter(|s| s.cluster != 0) {
        seq.cluster = old_to_new[&seq.cluster];
    }
    old_to_new
}


fn get_assembly_count(sequences: &[Sequence]) -> usize {
    sequences.iter().map(|s| &s.filename).collect::<HashSet<_>>().len()
}


fn get_max_cluster(sequences: &[Sequence]) -> u16 {
    sequences.iter().map(|s| s.cluster).max().unwrap()
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;
    use crate::tests::assert_almost_eq;

    #[test]
    fn test_get_assembly_count() {
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1, 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1, 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1, 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1, 1),
                             Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 1, 1)];
        assert_eq!(get_assembly_count(&sequences), 2);
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1, 1),
                             Sequence::new_with_seq(2, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1, 1),
                             Sequence::new_with_seq(3, "A".to_string(), "assembly_3.fasta".to_string(), "contig_1".to_string(), 1, 1),
                             Sequence::new_with_seq(4, "A".to_string(), "assembly_4.fasta".to_string(), "contig_1".to_string(), 1, 1),
                             Sequence::new_with_seq(5, "A".to_string(), "assembly_5.fasta".to_string(), "contig_1".to_string(), 1, 1)];
        assert_eq!(get_assembly_count(&sequences), 5);
    }

    #[test]
    fn test_get_max_cluster() {
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1, 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1, 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 1, 1)];
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3; sequences[3].cluster = 1; sequences[4].cluster = 2;
        assert_eq!(get_max_cluster(&sequences), 3);
        sequences[0].cluster = 1; sequences[1].cluster = 2; sequences[2].cluster = 3; sequences[3].cluster = 5; sequences[4].cluster = 4;
        assert_eq!(get_max_cluster(&sequences), 5);
    }

    #[test]
    fn test_get_closest_pair_ties_and_reverse_distances() {
        let distances = HashMap::from([
            ((1, 1), 0.0), ((3, 1), 0.2), ((2, 1), 0.2), ((2, 3), 0.2),
        ]);
        assert_eq!(get_closest_pair(&distances), (1, 2, 0.2));
        let mut distances = distances;
        distances.insert((1, 2), 0.4);
        assert_eq!(get_closest_pair(&distances), (1, 3, 0.2));
    }

    #[test]
    fn test_qc_clusters_automatic_and_manual() {
        let tree = test_tree_1();
        let mut sequences = vec![
            Sequence::new_without_seq(1, "a".into(), "a Autocycler_cluster_weight=2".into(), 100, 0),
            Sequence::new_without_seq(2, "b".into(), "b Autocycler_trusted".into(), 100, 0),
            Sequence::new_without_seq(3, "a".into(), "c".into(), 40, 0),
            Sequence::new_without_seq(4, "b".into(), "d".into(), 10, 0),
            Sequence::new_without_seq(5, "c".into(), "e".into(), 10, 0),
        ];
        let mut distances: HashMap<_, _> = (1..=5).flat_map(|a| (1..=5)
            .map(move |b| ((a, b), if a == b { 0.0 } else { 1.0 }))).collect();
        distances.insert((3, 1), 0.0);
        let clusters = vec![1, 2, 3, 6];
        let qc = qc_clusters(&tree, &mut sequences, &distances, &clusters, &[], 0.2, 2);
        assert_eq!(sequences.iter().map(|s| s.cluster).collect::<Vec<_>>(), vec![1, 2, 3, 4, 4]);
        assert!(qc[&1].pass() && qc[&2].pass() && qc[&4].pass());
        assert_eq!(qc[&3].failure_reasons,
                   vec!["present in too few assemblies", "contained within cluster 1"]);
        assert_eq!(qc[&4].cluster_dist, 0.2);

        let qc = qc_clusters(&tree, &mut sequences, &distances, &clusters, &[1, 3, 6], 0.2, 2);
        assert!(qc[&1].pass() && qc[&3].pass() && qc[&4].pass());
        assert_eq!(qc[&2].failure_reasons, vec!["not included in manual clusters"]);
    }

    #[test]
    fn test_upgma_1() {
        // Uses example from https://en.wikipedia.org/wiki/UPGMA
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "a".to_string(), "a".to_string(), 1, 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "b".to_string(), "b".to_string(), 1, 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "c".to_string(), "c".to_string(), 1, 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "d".to_string(), "d".to_string(), 1, 1),
                                 Sequence::new_with_seq(5, "A".to_string(), "e".to_string(), "e".to_string(), 1, 1)];
        let distances = HashMap::from_iter(vec![((1, 1), 00.0), ((1, 2), 17.0), ((1, 3), 21.0), ((1, 4), 31.0), ((1, 5), 23.0),
                                                ((2, 1), 17.0), ((2, 2), 00.0), ((2, 3), 30.0), ((2, 4), 34.0), ((2, 5), 21.0),
                                                ((3, 1), 21.0), ((3, 2), 30.0), ((3, 3), 00.0), ((3, 4), 28.0), ((3, 5), 39.0),
                                                ((4, 1), 31.0), ((4, 2), 34.0), ((4, 3), 28.0), ((4, 4), 00.0), ((4, 5), 43.0),
                                                ((5, 1), 23.0), ((5, 2), 21.0), ((5, 3), 39.0), ((5, 4), 43.0), ((5, 5), 00.0)]);
        let mut root = upgma(&distances, &sequences);
        assert_almost_eq(root.distance, 16.5, 1e-8);

        let index: HashMap<u16, &Sequence> = sequences.iter().map(|s| (s.id, s)).collect();
        let newick_string = tree_to_newick(&root, &index);
        assert_eq!(newick_string, "(((1__a__a__1_bp:8.5,2__b__b__1_bp:8.5)6:2.5,5__e__e__1_bp:11)7:5.5,(3__c__c__1_bp:14,4__d__d__1_bp:14)8:2.5)9");

        normalise_tree(&mut root);
        assert_almost_eq(root.distance, 0.5, 1e-8);
    }

    #[test]
    fn test_upgma_2() {
        let sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "a".to_string(), "a".to_string(), 1, 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "b".to_string(), "b".to_string(), 1, 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "c".to_string(), "c".to_string(), 1, 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "d".to_string(), "d".to_string(), 1, 1)];
        let distances = HashMap::from_iter(vec![((1, 1), 0.0), ((1, 2), 0.1), ((1, 3), 0.5), ((1, 4), 0.5),
                                                ((2, 1), 0.1), ((2, 2), 0.0), ((2, 3), 0.5), ((2, 4), 0.5),
                                                ((3, 1), 0.5), ((3, 2), 0.5), ((3, 3), 0.0), ((3, 4), 0.2),
                                                ((4, 1), 0.5), ((4, 2), 0.5), ((4, 3), 0.2), ((4, 4), 0.0)]);
        let mut root = upgma(&distances, &sequences);
        normalise_tree(&mut root);
        assert_almost_eq(root.distance, 0.25, 1e-8);

        let index: HashMap<u16, &Sequence> = sequences.iter().map(|s| (s.id, s)).collect();
        let newick_string = tree_to_newick(&root, &index);
        assert_eq!(newick_string, "((1__a__a__1_bp:0.05,2__b__b__1_bp:0.05)5:0.2,(3__c__c__1_bp:0.1,4__d__d__1_bp:0.1)6:0.15)7");
    }

    #[test]
    fn test_reorder_clusters() {
        let mut sequences = vec![Sequence::new_with_seq(1, "CGCGA".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 5, 1),
                                 Sequence::new_with_seq(2, "T".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1, 1),
                                 Sequence::new_with_seq(3, "AACGACTACG".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 10, 1),
                                 Sequence::new_with_seq(4, "CGCGA".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 5, 1),
                                 Sequence::new_with_seq(5, "T".to_string(), "assembly_2.fasta".to_string(), "contig_3".to_string(), 1, 1),
                                 Sequence::new_with_seq(6, "AACGACTACG".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 10, 1)];
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
    fn test_set_minpts() {
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "assembly_3.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "assembly_4.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(5, "A".to_string(), "assembly_5.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(6, "A".to_string(), "assembly_6.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(7, "A".to_string(), "assembly_7.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(8, "A".to_string(), "assembly_8.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(9, "A".to_string(), "assembly_9.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(10, "A".to_string(), "assembly_10.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(11, "A".to_string(), "assembly_10.fasta".to_string(), "contig_2".to_string(), 1, 1),
                                 Sequence::new_with_seq(12, "A".to_string(), "assembly_10.fasta".to_string(), "contig_3".to_string(), 1, 1),
                                 Sequence::new_with_seq(10, "A".to_string(), "assembly_11.fasta".to_string(), "contig_1".to_string(), 1, 1),
                                 Sequence::new_with_seq(11, "A".to_string(), "assembly_11.fasta".to_string(), "contig_2".to_string(), 1, 1),
                                 Sequence::new_with_seq(12, "A".to_string(), "assembly_12.fasta".to_string(), "contig_1".to_string(), 1, 1)];

        assert_eq!(set_min_assemblies(Some(2), &sequences), 2);
        assert_eq!(set_min_assemblies(Some(11), &sequences), 11);
        assert_eq!(set_min_assemblies(Some(321), &sequences), 321);
        assert_eq!(set_min_assemblies(None, &sequences), 3);  // 12 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 3);  // 11 assemblies
        sequences.truncate(9);
        assert_eq!(set_min_assemblies(None, &sequences), 2);  // 9 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 2);  // 8 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 2);  // 7 assemblies
        sequences.truncate(5);
        assert_eq!(set_min_assemblies(None, &sequences), 2);  // 5 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 2);  // 4 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 2);  // 3 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 2);  // 2 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 1);  // 1 assembly
    }

    fn test_tree_1() -> TreeNode {
        // Generates this tree: (1:0.5,(2:0.3,(3:0.2,(4:0.1,5:0.1):0.1):0.1):0.2);
        let n1 = TreeNode { id: 1, ..Default::default() };
        let n2 = TreeNode { id: 2, ..Default::default() };
        let n3 = TreeNode { id: 3, ..Default::default() };
        let n4 = TreeNode { id: 4, ..Default::default() };
        let n5 = TreeNode { id: 5, ..Default::default() };
        let n6 = TreeNode { id: 6, left: Some(Box::new(n4)), right: Some(Box::new(n5)), distance: 0.1 };
        let n7 = TreeNode { id: 7, left: Some(Box::new(n3)), right: Some(Box::new(n6)), distance: 0.2 };
        let n8 = TreeNode { id: 8, left: Some(Box::new(n2)), right: Some(Box::new(n7)), distance: 0.3 };
        TreeNode { id: 9, left: Some(Box::new(n1)), right: Some(Box::new(n8)), distance: 0.5 }
    }

    fn test_tree_2() -> TreeNode {
        // Generates this tree: (1:0.5,((2:0.1,3:0.1):0.2,(4:0.2,(5:0.1,6:0.1):0.1):0.1):0.2);
        let n1 = TreeNode { id: 1, ..Default::default() };
        let n2 = TreeNode { id: 2, ..Default::default() };
        let n3 = TreeNode { id: 3, ..Default::default() };
        let n4 = TreeNode { id: 4, ..Default::default() };
        let n5 = TreeNode { id: 5, ..Default::default() };
        let n6 = TreeNode { id: 6, ..Default::default() };
        let n7 = TreeNode { id: 7, left: Some(Box::new(n2)), right: Some(Box::new(n3)), distance: 0.1 };
        let n8 = TreeNode { id: 8, left: Some(Box::new(n5)), right: Some(Box::new(n6)), distance: 0.1 };
        let n9 = TreeNode { id: 9, left: Some(Box::new(n4)), right: Some(Box::new(n8)), distance: 0.2 };
        let n10 = TreeNode { id: 10, left: Some(Box::new(n7)), right: Some(Box::new(n9)), distance: 0.3 };
        TreeNode { id: 11, left: Some(Box::new(n1)), right: Some(Box::new(n10)), distance: 0.5 }
    }

    #[test]
    fn test_automatic_clustering() {
        let tree = test_tree_1();
        assert_eq!(tree.automatic_clustering(0.8), vec![1, 8]);
        assert_eq!(tree.automatic_clustering(0.5), vec![1, 2, 7]);
        assert_eq!(tree.automatic_clustering(0.3), vec![1, 2, 3, 6]);
        assert_eq!(tree.automatic_clustering(0.1), vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_contains_manual_cluster() {
        let tree = test_tree_1();
        assert!(!tree.contains_manual_cluster(&[]));
        for n in 1..=9   { assert!( tree.contains_manual_cluster(&[n])); }
        for n in 10..=19 { assert!(!tree.contains_manual_cluster(&[n])); }
    }

    #[test]
    fn test_manual_clustering() {
        let tree = test_tree_1();
        assert_eq!(tree.manual_clustering(0.5, &[]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &[1]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &[1, 2]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &[1, 2, 7]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &[3]), vec![1, 2, 3, 6]);
        assert_eq!(tree.manual_clustering(0.5, &[4]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.5, &[5]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.5, &[1, 2, 3, 4, 5]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.8, &[]), vec![1, 8]);
        assert_eq!(tree.manual_clustering(0.8, &[1]), vec![1, 8]);
        assert_eq!(tree.manual_clustering(0.8, &[2]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.8, &[3]), vec![1, 2, 3, 6]);
        assert_eq!(tree.manual_clustering(0.8, &[4]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.8, &[5]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.8, &[6]), vec![1, 2, 3, 6]);
        assert_eq!(tree.manual_clustering(0.8, &[7]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.8, &[8]), vec![1, 8]);
    }

    #[test]
    fn test_manual_clustering_missing_nodes() {
        let tree = test_tree_1();
        for id in [0, 10, u16::MAX] {
            for clusters in [vec![id], vec![1, id]] {
                let error = panic::catch_unwind(|| tree.manual_clustering(0.5, &clusters)).unwrap_err();
                assert_eq!(error.downcast_ref::<String>().unwrap(),
                           &format!("clustering tree does not contain a node with id {id}"));
            }
        }
    }

    #[test]
    fn test_check_consistency() {
        let tree = test_tree_1();
        tree.check_consistency(&[1, 2, 3, 4, 5]);
        tree.check_consistency(&[1, 2, 3, 6]);
        tree.check_consistency(&[1, 2, 7]);
        tree.check_consistency(&[1, 8]);
        tree.check_consistency(&[9]);
        assert!(panic::catch_unwind(|| {
            tree.check_consistency(&[5, 6]);
        }).is_err());
        assert!(panic::catch_unwind(|| {
            tree.check_consistency(&[6, 8]);
        }).is_err());
        assert!(panic::catch_unwind(|| {
            tree.check_consistency(&[1, 9]);
        }).is_err());
    }

    #[test]
    fn test_max_pairwise_distance() {
        let tree = test_tree_1();
        assert_almost_eq(tree.max_pairwise_distance(1), 0.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(2), 0.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(3), 0.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(4), 0.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(5), 0.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(6), 0.2, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(7), 0.4, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(8), 0.6, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(9), 1.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(10), -1.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(11), -1.0, 1e-8);
        assert_almost_eq(tree.max_pairwise_distance(12), -1.0, 1e-8);
    }

    #[test]
    fn test_get_tips() {
        let tree = test_tree_1();
        assert_eq!(tree.get_tips(1), vec![1]);
        assert_eq!(tree.get_tips(2), vec![2]);
        assert_eq!(tree.get_tips(3), vec![3]);
        assert_eq!(tree.get_tips(4), vec![4]);
        assert_eq!(tree.get_tips(5), vec![5]);
        assert_eq!(tree.get_tips(6), vec![4, 5]);
        assert_eq!(tree.get_tips(7), vec![3, 4, 5]);
        assert_eq!(tree.get_tips(8), vec![2, 3, 4, 5]);
        assert_eq!(tree.get_tips(9), vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_check_complete_coverage() {
        let tree = test_tree_1();
        tree.check_complete_coverage(&[1, 2, 3, 4, 5]);
        tree.check_complete_coverage(&[1, 2, 3, 6]);
        tree.check_complete_coverage(&[1, 2, 7]);
        tree.check_complete_coverage(&[1, 8]);
        tree.check_complete_coverage(&[9]);
        assert!(panic::catch_unwind(|| {
            tree.check_complete_coverage(&[1, 2, 3, 4, 5, 6]);
        }).is_err());
        assert!(panic::catch_unwind(|| {
            tree.check_complete_coverage(&[1, 2, 3, 4]);
        }).is_err());
        assert!(panic::catch_unwind(|| {
            tree.check_complete_coverage(&[1, 6, 7]);
        }).is_err());
    }

    #[test]
    fn test_split_clusters() {
        let tree = test_tree_1();
        assert_eq!(tree.split_clusters(&[1, 2, 3, 6]), vec![vec![1, 2, 3, 4, 5]]);
        assert_eq!(tree.split_clusters(&[1, 2, 7]), vec![vec![1, 2, 3, 6]]);
        assert_eq!(tree.split_clusters(&[1, 8]), vec![vec![1, 2, 7]]);
        assert_eq!(tree.split_clusters(&[9]), vec![vec![1, 8]]);

        let tree = test_tree_2();
        assert_eq!(tree.split_clusters(&[1, 4, 5, 6, 7]), vec![vec![1, 2, 3, 4, 5, 6]]);
        assert_eq!(tree.split_clusters(&[1, 2, 3, 4, 8]), vec![vec![1, 2, 3, 4, 5, 6]]);
        assert_eq!(tree.split_clusters(&[1, 4, 7, 8]), vec![vec![1, 2, 3, 4, 8], vec![1, 4, 5, 6, 7]]);
    }

    #[test]
    fn test_find_node() {
        let tree = test_tree_1();
        assert_eq!(tree.find_node(1).unwrap().id, 1);
        assert_eq!(tree.find_node(2).unwrap().id, 2);
        assert_eq!(tree.find_node(3).unwrap().id, 3);
        assert_eq!(tree.find_node(4).unwrap().id, 4);
        assert_eq!(tree.find_node(5).unwrap().id, 5);
        assert_eq!(tree.find_node(6).unwrap().id, 6);
        assert_eq!(tree.find_node(7).unwrap().id, 7);
        assert_eq!(tree.find_node(8).unwrap().id, 8);
        assert_eq!(tree.find_node(9).unwrap().id, 9);
        assert!(tree.find_node(10).is_none());
        assert!(tree.find_node(11).is_none());
        assert!(tree.find_node(12).is_none());
    }

    #[test]
    fn test_cluster_assembly_count_1() {
        // Simple case without any weights
        let mut seq_1 = Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1".to_string(), 1, 1); seq_1.cluster = 1;
        let mut seq_2 = Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1, 1); seq_2.cluster = 2;
        let mut seq_3 = Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1, 1); seq_3.cluster = 3;
        let mut seq_4 = Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1, 1); seq_4.cluster = 1;
        let mut seq_5 = Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2".to_string(), 1, 1); seq_5.cluster = 3;
        let sequences = vec![seq_1, seq_2, seq_3, seq_4, seq_5];
        assert_eq!(cluster_assembly_count(&sequences, 1), 2);
        assert_eq!(cluster_assembly_count(&sequences, 2), 1);
        assert_eq!(cluster_assembly_count(&sequences, 3), 2);
    }

    #[test]
    fn test_cluster_assembly_count_2() {
        // Various weights
        let mut seq_1 = Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1 Autocycler_cluster_weight=3 other stuff".to_string(), 1, 1); seq_1.cluster = 1;
        let mut seq_2 = Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2 other stuff autocycler_cluster_weight=6".to_string(), 1, 1); seq_2.cluster = 2;
        let mut seq_3 = Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3".to_string(), 1, 1); seq_3.cluster = 3;
        let mut seq_4 = Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1, 1); seq_4.cluster = 1;
        let mut seq_5 = Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2 AuToCyCleR_cluster_weight=0".to_string(), 1, 1); seq_5.cluster = 3;
        let sequences = vec![seq_1, seq_2, seq_3, seq_4, seq_5];
        assert_eq!(cluster_assembly_count(&sequences, 1), 4);
        assert_eq!(cluster_assembly_count(&sequences, 2), 6);
        assert_eq!(cluster_assembly_count(&sequences, 3), 1);
    }

    #[test]
    fn test_cluster_assembly_count_3() {
        // Multiple different weights for the same assembly
        let mut seq_1 = Sequence::new_with_seq(1, "A".to_string(), "assembly_1.fasta".to_string(), "contig_1 Autocycler_cluster_weight=3".to_string(), 1, 1); seq_1.cluster = 1;
        let mut seq_2 = Sequence::new_with_seq(2, "A".to_string(), "assembly_1.fasta".to_string(), "contig_2".to_string(), 1, 1); seq_2.cluster = 1;
        let mut seq_3 = Sequence::new_with_seq(3, "A".to_string(), "assembly_1.fasta".to_string(), "contig_3 other stuff Autocycler_cluster_weight=2".to_string(), 1, 1); seq_3.cluster = 1;
        let mut seq_4 = Sequence::new_with_seq(4, "A".to_string(), "assembly_2.fasta".to_string(), "contig_1".to_string(), 1, 1); seq_4.cluster = 2;
        let mut seq_5 = Sequence::new_with_seq(5, "A".to_string(), "assembly_2.fasta".to_string(), "contig_2 Autocycler_cluster_weight=0 other stuff".to_string(), 1, 1); seq_5.cluster = 2;
        let sequences = vec![seq_1, seq_2, seq_3, seq_4, seq_5];
        assert_eq!(cluster_assembly_count(&sequences, 1), 3);
        assert_eq!(cluster_assembly_count(&sequences, 2), 1);
    }
}
