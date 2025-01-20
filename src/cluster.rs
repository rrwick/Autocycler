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
use std::path::{Path, PathBuf};

use crate::graph_simplification::merge_linear_paths;
use crate::log::{section_header, explanation};
use crate::metrics::{ClusteringMetrics, UntrimmedClusterMetrics};
use crate::misc::{check_if_dir_exists, check_if_file_exists, format_float, median_usize,
                  quit_with_error, usize_division_rounded, create_dir, delete_dir_if_exists,
                  load_file_lines};
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
    let manual_clusters = parse_manual_clusters(manual_clusters);
    print_settings(&autocycler_dir, cutoff, min_assemblies, min_assemblies_option, max_contigs,
                   &manual_clusters);
    check_sequence_count(&sequences, max_contigs);
    let asymmetrical_distances = pairwise_contig_distances(&graph, &sequences, &pairwise_phylip);
    let symmetrical_distances = make_symmetrical_distances(&asymmetrical_distances, &sequences);
    let mut tree = upgma(&symmetrical_distances, &mut sequences);
    normalise_tree(&mut tree);
    save_tree_to_newick(&tree, &sequences, &clustering_newick);
    let qc_results = generate_clusters(&tree, &mut sequences, &asymmetrical_distances, cutoff,
                                       min_assemblies, &manual_clusters);
    save_clusters(&sequences, &qc_results, &clustering_dir, &gfa_lines);
    save_data_to_tsv(&sequences, &qc_results, &clustering_tsv);
    let metrics = clustering_metrics(&sequences, &qc_results);
    metrics.save_to_yaml(&clustering_yaml);

    // TODO: create a PDF of the tree with clusters? printpdf?

    finished_message(&pairwise_phylip, &clustering_newick, &clustering_tsv);
}


fn check_settings(autocycler_dir: &Path, gfa: &Path, cutoff: f64, min_assemblies: &Option<usize>) {
    check_if_dir_exists(autocycler_dir);
    check_if_file_exists(gfa);
    if cutoff <= 0.0 || cutoff >= 1.0 {
        quit_with_error("--cutoff must be between 0 and 1 (exclusive)");
    }
    if min_assemblies.is_some() && min_assemblies.unwrap() < 1 {
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
        eprintln!("  --min_assemblies {} (automatically set)", min_assemblies);
    } else {
        eprintln!("  --min_assemblies {}", min_assemblies);
    }
    eprintln!("  --max_contigs {}", max_contigs);
    if !manual_clusters.is_empty() {
        eprintln!("  --manual {}", manual_clusters.iter().map(|c| c.to_string())
                                                  .collect::<Vec<String>>() .join(","));
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
        let e = format!("the mean number of contigs per input assembly ({:.1}) exceeds the allowed \
                         threshold ({}). Are your input assemblies fragmented or contaminated?",
                        mean_seqs_per_assembly, max_contigs);
        quit_with_error(&e);
    }
}


fn pairwise_contig_distances(graph: &UnitigGraph, sequences: &Vec<Sequence>, file_path: &Path)
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
        let a_len = total_unitig_length(a, &unitig_lengths) as f64;
        for seq_b in sequences {
            let b = sequence_unitigs.get(&seq_b.id).unwrap();
            let ab: HashSet<u32> = a.intersection(b).cloned().collect();
            let ab_len = total_unitig_length(&ab, &unitig_lengths) as f64;
            let distance = 1.0 - (ab_len / a_len);
            distances.insert((seq_a.id, seq_b.id), distance);
        }
    }
    eprintln!("{} sequences, {} total pairwise distances", sequences.len(), distances.len());
    eprintln!();
    save_distance_matrix(&distances, sequences, file_path);
    distances
}


fn total_unitig_length(unitigs: &HashSet<u32>, unitig_lengths: &HashMap<u32, u32>) -> u32 {
    unitigs.iter().map(|u| unitig_lengths.get(u).unwrap()).sum()
}


fn save_distance_matrix(distances: &HashMap<(u16, u16), f64>, sequences: &Vec<Sequence>,
                        file_path: &Path) {
    eprintln!("Saving distance matrix:");
    let mut f = File::create(file_path).unwrap();
    writeln!(f, "{}", sequences.len()).unwrap();
    for seq_a in sequences {
        write!(f, "{}", seq_a).unwrap();
        for seq_b in sequences {
            write!(f, "\t{:.8}", distances.get(&(seq_a.id, seq_b.id)).unwrap()).unwrap();
        }
        writeln!(f).unwrap();
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
        // This method is run on the root of the tree, and it returns the distance for the given
        // node. The distance is doubled because this function returns the max pairwise distance
        // within the clade, not the node-to-tip-distance.
        if self.id == node_num { return self.distance * 2.0; }
        if self.is_tip() { return -1.0; }  // -1 means the node_num wasn't found
        let left_dist = self.left.as_ref().unwrap().max_pairwise_distance(node_num);
        let right_dist = self.right.as_ref().unwrap().max_pairwise_distance(node_num);
        left_dist.max(right_dist)
    }

    fn automatic_clustering(&self, cutoff: f64) -> Vec<u16> {
        // This method is run on the root of the tree, and it divides the tree into clusters using
        // only the cutoff distance. 
        let mut clusters = Vec::new();
        self.collect_clusters(cutoff / 2.0, &[], &mut clusters);
        clusters.sort();
        clusters
    }

    fn manual_clustering(&self, cutoff: f64, manual_clusters: &[u16]) -> Vec<u16> {
        // This method is run on the root of the tree, and it divides the tree into clusters. Any
        // manual clusters specified by the user are included in the final clusters, and the rest
        // of the tree is clustered by distance.
        let mut clusters = Vec::new();
        self.check_consistency(manual_clusters);
        self.collect_clusters(cutoff / 2.0, manual_clusters, &mut clusters);
        clusters.sort();
        clusters
    }

    fn collect_clusters(&self, cutoff: f64, manual_clusters: &[u16], clusters: &mut Vec<u16>) {
        if manual_clusters.contains(&self.id) ||
                (self.distance <= cutoff && !self.has_manual_child(manual_clusters)) {
            clusters.push(self.id);
        } else if !self.is_tip() {
            self.left.as_ref().unwrap().collect_clusters(cutoff, manual_clusters, clusters);
            self.right.as_ref().unwrap().collect_clusters(cutoff, manual_clusters, clusters);
        }
    }

    fn has_manual_child(&self, manual_clusters: &[u16]) -> bool {
        // Returns true if this node is in the user-specified manual clusters or any of its
        // children are.
        if manual_clusters.contains(&self.id) { return true; }
        if !self.is_tip() {
            if self.left.as_ref().unwrap().has_manual_child(manual_clusters) { return true; }
            if self.right.as_ref().unwrap().has_manual_child(manual_clusters) { return true; }
        }
        false
    }

    fn check_consistency(&self, manual_clusters: &[u16]) {
        // Ensures that no manual cluster is contained within another manual cluster.
        if !self.is_tip() {
            if manual_clusters.contains(&self.id) &&
                    (self.left.as_ref().unwrap().has_manual_child(manual_clusters) ||
                     self.right.as_ref().unwrap().has_manual_child(manual_clusters)) {
                quit_with_error("manual clusters cannot be nested");
            }
            self.left.as_ref().unwrap().check_consistency(manual_clusters);
            self.right.as_ref().unwrap().check_consistency(manual_clusters);
        }
    }

    fn get_tips(&self, node_num: u16) -> Vec<u16> {
        // This method is run on the root of the tree, and it returns the node IDs of the tips under
        // the given node.
        if self.id == node_num {
            let mut tips = Vec::new();
            self.collect_tips(&mut tips);
            return tips;
        }
        if self.is_tip() { return vec![]; }
        let mut left_tips =  self.left.as_ref().unwrap().get_tips(node_num);
        let right_tips =  self.right.as_ref().unwrap().get_tips(node_num);
        left_tips.extend(right_tips);
        left_tips
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
        // This method is run on the root of the tree, and it ensures that the given clusters nicely
        // cover the entire tree, i.e. every tip is in one and only one cluster.
        let all_tips: HashSet<u16> = self.get_tips(self.id).into_iter().collect();
        let mut covered_tips = HashSet::new();
        for &c in clusters {
            let cluster_tips = self.get_tips(c);
            for tip in cluster_tips {
                if !covered_tips.insert(tip) { panic!("overlap detected"); }
            }
        }
        if covered_tips != all_tips { panic!("incomplete coverage");}
    }

    fn split_clusters(&self, clusters: &[u16]) -> Vec<Vec<u16>> {
        // This method is run on the root of the tree, and given a clustering, it returns all
        // possible clusterings where a splitable cluster has been split.
        self.check_complete_coverage(clusters);
        let mut result = Vec::new();

        for &cluster in clusters {
            let node = self.find_node(cluster).unwrap();
            if !node.is_tip() {
                // Create a new clustering with this cluster split into its children
                let mut new_cluster = Vec::new();
                for &other_cluster in clusters {
                    if other_cluster != cluster {
                        new_cluster.push(other_cluster);
                    }
                }
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
        // Given a number, this method will search for a node with the matching number, and return
        // a reference to that node if found.
        if self.id == node_num { return Some(self); }
        if self.is_tip() { return None; }
        let left_result = self.left.as_ref().unwrap().find_node(node_num);
        if let Some(left_node) = left_result { return Some(left_node); }
        let right_result = self.right.as_ref().unwrap().find_node(node_num);
        if let Some(right_node) = right_result { return Some(right_node); }
        None
    }
}


fn find_node_by_id(node: &TreeNode, id: u16) -> Option<&TreeNode> {
    if node.id == id { return Some(node); }
    if let Some(ref left) = node.left {
        if let Some(found) = find_node_by_id(left, id) { return Some(found); }
    }
    if let Some(ref right) = node.right {
        if let Some(found) = find_node_by_id(right, id) { return Some(found); }
    }
    None
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
        writeln!(file, "({}:{});", newick_string, root_length).unwrap();
    } else {
        writeln!(file, "{};", newick_string).unwrap();
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
        _ => index.get(&node.id).unwrap().string_for_newick().to_string(),
    }
}


fn upgma(distances: &HashMap<(u16, u16), f64>, sequences: &mut Vec<Sequence>) -> TreeNode {
    section_header("Clustering sequences");
    explanation("Contigs are organise into a tree using UPGMA. Then clusters are defined from the \
                 tree using the distance cutoff.");
    let mut clusters: HashMap<u16, HashSet<u16>> = HashMap::new();
    let mut cluster_distances: HashMap<(u16, u16), f64> = distances.clone();
    let mut nodes: HashMap<u16, TreeNode> = HashMap::new();
    let mut internal_node_num: u16 = sequences.iter().map(|s| (s.id)).max().unwrap();

    // Initialise each sequence as its own cluster and create initial nodes.
    for seq in sequences {
        clusters.insert(seq.id, HashSet::from([seq.id]));
        nodes.insert(seq.id, TreeNode { id: seq.id, ..Default::default() });
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
        internal_node_num += 1;
        let new_node = TreeNode {
            id: internal_node_num,
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


fn generate_clusters(tree: &TreeNode, sequences: &mut Vec<Sequence>,
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


fn qc_clusters(tree: &TreeNode, sequences: &mut Vec<Sequence>, distances: &HashMap<(u16, u16), f64>,
               cluster_nodes: &Vec<u16>, manual_clusters: &[u16], cutoff: f64,
               min_assemblies: usize) -> HashMap<u16, ClusterQC> {
    // Given a set of node numbers for the tree which define clusters, this function returns the
    // QC-results HashMap.

    // TODO: I could add some additional logic input assembly preferences. For example, I could use
    //       Plassembler as one of the input assemblies, and any cluster which contains a
    //       Plassembler contig should not fail QC due to low contig count. Could do this via
    //       string matching on assembly names, e.g. 'plassembler'.

    // Create the ClusterQC object for each cluster. If using manual clustering, the clusters will
    // fail if they aren't included in the user-supplied clusters. If using automatic clustering,
    // all clusters will initially pass.
    let mut current_cluster = 0;
    let mut qc_results = HashMap::new();
    for n in cluster_nodes {
        if let Some(node) = find_node_by_id(tree, *n) {
            current_cluster += 1;
            assign_cluster_to_node(node, sequences, current_cluster);
            let mut qc = ClusterQC::new(tree.max_pairwise_distance(*n));
            if !manual_clusters.is_empty() && !manual_clusters.contains(n) {
                qc.failure_reasons.push("not included in manual clusters".to_string());
            }
            qc_results.insert(current_cluster, qc);
        } else {
            quit_with_error(&format!("clustering tree does not contain a node with id {}", n));
        }
    }

    // Reorder clusters from large to small.
    let old_to_new = reorder_clusters(sequences);
    let mut reordered_qc_results = HashMap::new();
    for (old_key, qc) in qc_results {
        if let Some(&new_key) = old_to_new.get(&old_key) {
            reordered_qc_results.insert(new_key, qc);
        }
    }
    qc_results = reordered_qc_results;

    // If using automatic clustering, clusters are now failed for being contained in other clusters
    // or appearing in too few input assemblies.
    if manual_clusters.is_empty() {
        let max_cluster = get_max_cluster(sequences);
        for c in 1..=max_cluster {
            let assemblies: HashSet<_> = sequences.iter().filter(|s| s.cluster == c)
                                                  .map(|s| s.filename.clone()).collect();
            if assemblies.len() < min_assemblies {
                let fail_reason = "present in too few assemblies".to_string();
                qc_results.get_mut(&c).unwrap().failure_reasons.push(fail_reason);
            }
        }
        for c in 1..=max_cluster {
            let container = cluster_is_contained_in_another(c, sequences, distances, cutoff,
                                                            &qc_results);
            if container > 0 {
                let fail_reason = format!("contained within cluster {}", container);
                qc_results.get_mut(&c).unwrap().failure_reasons.push(fail_reason);
            }
        }
    }
    qc_results
}


fn score_clustering(tree: &TreeNode, sequences: &mut Vec<Sequence>,
                    distances: &HashMap<(u16, u16), f64>, clusters: &Vec<u16>, cutoff: f64,
                    min_assemblies: usize) -> f64 {
    // Given a set of node numbers for the tree which define clusters, this function returns the
    // overall score for that clustering (higher is better).
    let qc_results = qc_clusters(tree, sequences, distances, clusters, &[], cutoff,
                                 min_assemblies);
    let metrics = clustering_metrics(sequences, &qc_results);
    metrics.overall_clustering_score
}


fn refine_auto_clusters(tree: &TreeNode, sequences: &mut Vec<Sequence>,
                        distances: &HashMap<(u16, u16), f64>, clusters: &[u16], cutoff: f64,
                        min_assemblies: usize) -> Vec<u16> {
    // Given a set of node numbers for the tree which define clusters, this function tries to
    // improve the clustering by splitting each cluster and checking if the score gets better,
    // repeating until no improvements can be made.
    let mut best_clusters = clusters.to_vec();
    let mut best_score = score_clustering(tree, sequences, distances, &best_clusters, cutoff,
                                          min_assemblies);
    let mut improved = true;
    while improved {
        improved = false;
        for alt_clusters in tree.split_clusters(&best_clusters) {
            let alt_score = score_clustering(tree, sequences, distances, &alt_clusters, cutoff,
                                             min_assemblies);
            if alt_score > best_score {
                best_clusters = alt_clusters;
                best_score = alt_score;
                improved = true;
            }
        }
    }
    best_clusters
}


fn assign_cluster_to_node(node: &TreeNode, sequences: &mut Vec<Sequence>, cluster: u16) {
    // This function assigns all sequences under a given node to the given cluster.
    for s in sequences.iter_mut() {
        if s.id == node.id {
            s.cluster = cluster;
        }
    }
    if let Some(ref left)  = node.left  { assign_cluster_to_node(left,  sequences, cluster); }
    if let Some(ref right) = node.right { assign_cluster_to_node(right, sequences, cluster); }
}


fn set_min_assemblies(min_assemblies_option: Option<usize>, sequences: &[Sequence]) -> usize {
    // This function automatically sets the --min_assemblies parameter, if the user didn't
    // explicitly supply one. The auto-set value will be one-quarter of the assembly count (rounded)
    // but no less than 2, unless there is only one input assembly, in which case it will be 1.
    if let Some(min_assemblies) = min_assemblies_option {
        return min_assemblies;
    }
    let assembly_count = get_assembly_count(sequences);
    if assembly_count == 1 {
        return 1;
    }
    let mut min_assemblies = usize_division_rounded(assembly_count, 4);
    if min_assemblies < 2 {
        min_assemblies = 2;
    }
    min_assemblies
}


fn parse_manual_clusters(manual_clusters: Option<String>) -> Vec<u16> {
    if manual_clusters.is_none() {
        return Vec::new();
    }
    let manual_clusters = manual_clusters.unwrap().replace(' ', "");
    let mut clusters: Vec<_> = manual_clusters.split(',')
            .map(|s| s.parse::<u16>().unwrap_or_else(|_| quit_with_error(
                &format!("failed to parse '{}' as a node number", s)))).collect();
    clusters.sort();
    clusters
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
    pub fn fail(&self) -> bool { !self.failure_reasons.is_empty() }
}


fn cluster_is_contained_in_another(cluster_num: u16, sequences: &[Sequence],
                                   distances: &HashMap<(u16, u16), f64>, cutoff: f64,
                                   qc_results: &HashMap<u16, ClusterQC>) -> u16 {
    // Checks whether this cluster is contained within another cluster that has so-far passed QC.
    // If so, it returns the id of the containing cluster. If not, it returns 0.
    // A cluster counts as contained if the majority of the pairwise comparisons to another cluster
    // are asymmetrical and below the cutoff.
    let passed_clusters: Vec<u16> = qc_results.iter().filter(|(_, q)| q.pass())
                                              .map(|(&k, _)| k).collect();
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
                 clustering_dir: &Path, gfa_lines: &Vec<String>) {
    let pass_dir = clustering_dir.join("qc_pass");
    let fail_dir = clustering_dir.join("qc_fail");
    save_qc_pass_clusters(sequences, qc_results, gfa_lines, &pass_dir);
    save_qc_fail_clusters(sequences, qc_results, gfa_lines, &fail_dir);
}


fn save_qc_pass_clusters(sequences: &[Sequence], qc_results: &HashMap<u16, ClusterQC>,
                         gfa_lines: &Vec<String>, pass_dir: &Path) {
    for c in 1..=get_max_cluster(sequences) {
        let qc = qc_results.get(&c).unwrap();
        if qc.pass() {
            eprintln!("Cluster {:03}:", c);
            let mut seq_count = 0;
            let mut seq_lengths = Vec::new();
            for s in sequences.iter().filter(|s| s.cluster == c) {
                eprintln!("  {}", s);
                seq_count += 1;
                seq_lengths.push(s.length);
            }
            if seq_count > 1 {
                eprintln!("  cluster distance: {}", format_float(qc.cluster_dist));
            }
            eprintln!("{}", "  passed QC".green());
            let cluster_dir = pass_dir.join(format!("cluster_{:03}", c));
            create_dir(&cluster_dir);
            save_cluster_gfa(sequences, c, gfa_lines, cluster_dir.join("1_untrimmed.gfa"));
            save_untrimmed_cluster_metrics(seq_lengths, qc.cluster_dist,
                                           cluster_dir.join("1_untrimmed.yaml"));
            eprintln!();
        }
    }
}


fn save_qc_fail_clusters(sequences: &[Sequence], qc_results: &HashMap<u16, ClusterQC>,
                         gfa_lines: &Vec<String>, fail_dir: &Path) {
    for c in 1..=get_max_cluster(sequences) {
        let qc = qc_results.get(&c).unwrap();
        if qc.fail() {
            eprintln!("Cluster {:03}:", c);
            let mut seq_count = 0;
            let mut seq_lengths = Vec::new();
            for s in sequences.iter().filter(|s| s.cluster == c) {
                eprintln!("  {}", s.to_string().dimmed());
                seq_count += 1;
                seq_lengths.push(s.length);
            }
            if seq_count > 1 {
                let s = format!("cluster distance: {}", format_float(qc.cluster_dist));
                eprintln!("  {}", s.to_string().dimmed());
            }
            for f in &qc.failure_reasons {
                eprintln!("  {}", format!("failed QC: {}", f).red());
            }
            let cluster_dir = fail_dir.join(format!("cluster_{:03}", c));
            create_dir(&cluster_dir);
            save_cluster_gfa(sequences, c, gfa_lines, cluster_dir.join("1_untrimmed.gfa"));
            save_untrimmed_cluster_metrics(seq_lengths, qc.cluster_dist,
                                           cluster_dir.join("1_untrimmed.yaml"));
            eprintln!();
        }
    }
}


fn save_cluster_gfa(sequences: &[Sequence], cluster_num: u16, gfa_lines: &Vec<String>,
                    out_gfa: PathBuf) {
    let cluster_seqs: Vec<Sequence> = sequences.iter().filter(|s| s.cluster == cluster_num)
                                               .cloned().collect();
    let (mut cluster_graph, _) = UnitigGraph::from_gfa_lines(gfa_lines);
    let seq_ids_to_remove:Vec<_> = sequences.iter().filter(|s| s.cluster != cluster_num)
                                            .map(|s| s.id).collect();
    for id in seq_ids_to_remove {
        cluster_graph.remove_sequence_from_graph(id);
    }
    cluster_graph.remove_zero_depth_unitigs();
    merge_linear_paths(&mut cluster_graph, &cluster_seqs);
    cluster_graph.save_gfa(&out_gfa, &cluster_seqs, false).unwrap();
}


fn save_untrimmed_cluster_metrics(seq_lengths: Vec<usize>, cluster_dist: f64,
                                  cluster_yaml: PathBuf) {
    let metrics = UntrimmedClusterMetrics::new(seq_lengths, cluster_dist);
    metrics.save_to_yaml(&cluster_yaml);
}


fn save_data_to_tsv(sequences: &Vec<Sequence>, qc_results: &HashMap<u16, ClusterQC>,
                    file_path: &Path) {
    let mut file = File::create(file_path).unwrap();
    writeln!(file, "node_name\tpassing_clusters\tall_clusters\tsequence_id\t\
                  file_name\tcontig_name\tlength").unwrap();
    for seq in sequences {
        assert!(seq.cluster != 0);
        let qc = qc_results.get(&seq.cluster).unwrap();
        let all_cluster = format!("{}", seq.cluster);
        let pass_cluster = if qc.pass() {
            "none".to_string()
        } else {
            format!("{}", seq.cluster)
        };
        writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}", seq.string_for_newick(),
               pass_cluster, all_cluster, seq.id, seq.filename, seq.contig_name(),
               seq.length).unwrap();
    }
}


fn clustering_metrics(sequences: &Vec<Sequence>, qc_results: &HashMap<u16, ClusterQC>)
        -> ClusteringMetrics {
    let mut metrics = ClusteringMetrics::default();
    let mut pass_cluster_distances = Vec::new();
    for c in 1..=get_max_cluster(sequences) {
        let qc = qc_results.get(&c).unwrap();
        if qc.pass() {
            metrics.pass_cluster_count += 1;
            pass_cluster_distances.push(qc.cluster_dist);
        } else {
            metrics.fail_cluster_count += 1;
        }
    }
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
    metrics.calculate_fractions();
    metrics.calculate_scores(cluster_filenames, pass_cluster_distances);
    metrics
}


fn reorder_clusters(sequences: &mut Vec<Sequence>) -> HashMap<u16, u16>{
    // Reorder clusters based on their median sequence length (large to small). Returns the mapping
    // of old cluster numbers to new cluster numbers.
    let mut cluster_lengths = HashMap::new();
    for c in 1..=get_max_cluster(sequences) {
        let lengths: Vec<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.length).collect();
        cluster_lengths.insert(c, median_usize(&lengths));
    }
    let mut sorted_cluster_lengths: Vec<_> = cluster_lengths.iter().collect();
    sorted_cluster_lengths.sort_by(|a, b| {match b.1.cmp(a.1) {std::cmp::Ordering::Equal => a.0.cmp(b.0), other => other,}});
    let mut old_to_new = HashMap::new();
    let mut new_c: u16 = 0;
    for (old_c, _length) in sorted_cluster_lengths {
        new_c += 1;
        old_to_new.insert(*old_c, new_c);
    }
    for s in sequences {
        if s.cluster < 1 {continue;}
        s.cluster = *old_to_new.get(&s.cluster).unwrap();
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
    fn test_upgma_1() {
        // Uses example from https://en.wikipedia.org/wiki/UPGMA
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "a".to_string(), "a".to_string(), 1, 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "b".to_string(), "b".to_string(), 1, 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "c".to_string(), "c".to_string(), 1, 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "d".to_string(), "d".to_string(), 1, 1),
                                 Sequence::new_with_seq(5, "A".to_string(), "e".to_string(), "e".to_string(), 1, 1)];
        let distances = HashMap::from_iter(vec![((1, 1), 00.0), ((1, 2), 17.0), ((1, 3), 21.0), ((1, 4), 31.0), ((1, 5), 23.0),
                                                ((2, 1), 17.0), ((2, 2), 00.0), ((2, 3), 30.0), ((2, 4), 34.0), ((2, 5), 21.0),
                                                ((3, 1), 21.0), ((3, 2), 30.0), ((3, 3), 00.0), ((3, 4), 28.0), ((3, 5), 39.0),
                                                ((4, 1), 31.0), ((4, 2), 34.0), ((4, 3), 28.0), ((4, 4), 00.0), ((4, 5), 43.0),
                                                ((5, 1), 23.0), ((5, 2), 21.0), ((5, 3), 39.0), ((5, 4), 43.0), ((5, 5), 00.0)]);
        let mut root = upgma(&distances, &mut sequences);
        assert_almost_eq(root.distance, 16.5, 1e-8);

        let index: HashMap<u16, &Sequence> = sequences.iter().map(|s| (s.id, s)).collect();
        let newick_string = tree_to_newick(&root, &index);
        assert_eq!(newick_string, "(((1__a__a__1_bp:8.5,2__b__b__1_bp:8.5)6:2.5,5__e__e__1_bp:11)7:5.5,(3__c__c__1_bp:14,4__d__d__1_bp:14)8:2.5)9");

        normalise_tree(&mut root);
        assert_almost_eq(root.distance, 0.5, 1e-8);
    }

    #[test]
    fn test_upgma_2() {
        let mut sequences = vec![Sequence::new_with_seq(1, "A".to_string(), "a".to_string(), "a".to_string(), 1, 1),
                                 Sequence::new_with_seq(2, "A".to_string(), "b".to_string(), "b".to_string(), 1, 1),
                                 Sequence::new_with_seq(3, "A".to_string(), "c".to_string(), "c".to_string(), 1, 1),
                                 Sequence::new_with_seq(4, "A".to_string(), "d".to_string(), "d".to_string(), 1, 1)];
        let distances = HashMap::from_iter(vec![((1, 1), 0.0), ((1, 2), 0.1), ((1, 3), 0.5), ((1, 4), 0.5),
                                                ((2, 1), 0.1), ((2, 2), 0.0), ((2, 3), 0.5), ((2, 4), 0.5),
                                                ((3, 1), 0.5), ((3, 2), 0.5), ((3, 3), 0.0), ((3, 4), 0.2),
                                                ((4, 1), 0.5), ((4, 2), 0.5), ((4, 3), 0.2), ((4, 4), 0.0)]);
        let mut root = upgma(&distances, &mut sequences);
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
    fn test_has_manual_child() {
        let tree = test_tree_1();
        assert!(!tree.has_manual_child(&[]));
        for n in 1..=9   { assert!( tree.has_manual_child(&[n])); }
        for n in 10..=19 { assert!(!tree.has_manual_child(&[n])); }
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
    fn test_parse_manual_clusters() {
        assert_eq!(parse_manual_clusters(Some("1,2,3".to_string())), vec![1, 2, 3]);
        assert_eq!(parse_manual_clusters(None), Vec::<u16>::new());
        assert!(panic::catch_unwind(|| {
            parse_manual_clusters(Some("x,y,z".to_string()));
        }).is_err());
        assert!(panic::catch_unwind(|| {
            parse_manual_clusters(Some("^&%^*".to_string()));
        }).is_err());
    }
}
