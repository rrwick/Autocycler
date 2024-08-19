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
use crate::metrics::{ClusteringMetrics, UntrimmedClusterMetrics, save_yaml};
use crate::misc::{check_if_dir_exists, check_if_file_exists, format_float, median_usize,
                  quit_with_error, usize_division_rounded, create_dir, delete_dir_if_exists,
                  load_file_lines};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn cluster(autocycler_dir: PathBuf, cutoff: f64, min_assemblies_option: Option<usize>,
               manual_clusters: Option<String>) {
    let gfa = autocycler_dir.join("input_assemblies.gfa");
    let clustering_dir = autocycler_dir.join("clustering");
    delete_dir_if_exists(&clustering_dir);
    create_dir(&clustering_dir);
    let pairwise_phylip = clustering_dir.join("pairwise_distances.phylip");
    let clustering_newick = clustering_dir.join("clustering.newick");
    let clustering_tsv = clustering_dir.join("clustering.tsv");
    let clustering_yaml = clustering_dir.join("clustering.yaml");
    check_settings(&autocycler_dir, &gfa, cutoff, &min_assemblies_option);
    starting_message();
    let gfa_lines = load_file_lines(&gfa);
    let (unitig_graph, mut sequences) = load_graph(&gfa_lines, true);
    let min_assemblies = set_min_assemblies(min_assemblies_option, &sequences);
    let manual_clusters = parse_manual_clusters(manual_clusters);
    print_settings(&autocycler_dir, cutoff, min_assemblies, min_assemblies_option,
                   &manual_clusters);
    let asymmetrical_distances = pairwise_contig_distances(&unitig_graph, &sequences,
                                                           &pairwise_phylip);
    let symmetrical_distances = make_symmetrical_distances(&asymmetrical_distances, &sequences);
    let mut tree = upgma(&symmetrical_distances, &mut sequences);
    normalise_tree(&mut tree);
    save_tree_to_newick(&tree, &sequences, &clustering_newick);
    let qc_results = generate_clusters(&tree, &mut sequences, &asymmetrical_distances, cutoff,
                                       &manual_clusters);




    // define_clusters_from_node_numbers(&tree, &mut sequences, clusters);





    // // let qc_results = if manual_clusters.is_empty() {
    // //     define_clusters_automatically(&tree, &mut sequences, &asymmetrical_distances, min_assemblies, cutoff)
    // // } else {
    // //     define_clusters_manually(&tree, &mut sequences, manual_clusters)
    // // };





    // save_clusters(&sequences, &qc_results, &clustering_dir, &gfa_lines);
    // save_data_to_tsv(&sequences, &qc_results, &clustering_tsv);
    // save_clustering_metrics(&clustering_yaml, &sequences, &qc_results);

    // // TODO: create a PDF of the tree with clusters? printpdf?

    // finished_message(&pairwise_phylip, &clustering_newick, &clustering_tsv);
}


fn check_settings(autocycler_dir: &PathBuf, gfa: &PathBuf, cutoff: f64, min_assemblies: &Option<usize>) {
    check_if_dir_exists(&autocycler_dir);
    check_if_file_exists(&gfa);
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


fn finished_message(pairwise_phylip: &PathBuf, clustering_newick: &PathBuf, clustering_tsv: &PathBuf) {
    section_header("Finished!");
    explanation("You can now run autocycler trim on each cluster. If you want to manually \
                 inspect the clustering, you can view the following files.");
    eprintln!("Pairwise distances:         {}", pairwise_phylip.display());
    eprintln!("Clustering tree (Newick):   {}", clustering_newick.display());
    eprintln!("Clustering tree (metadata): {}", clustering_tsv.display());
    eprintln!();
}


fn print_settings(autocycler_dir: &PathBuf, cutoff: f64, min_assemblies: usize,
                  min_assemblies_option: Option<usize>, manual_clusters: &Vec<u16>) {
    eprintln!("Settings:");
    eprintln!("  --autocycler_dir {}", autocycler_dir.display());
    eprintln!("  --cutoff {}", format_float(cutoff));
    if min_assemblies_option.is_none() {
        eprintln!("  --min_assemblies {} (automatically set)", min_assemblies);
    } else {
        eprintln!("  --min_assemblies {}", min_assemblies);
    }
    if !manual_clusters.is_empty() {
        eprintln!("  --manual {}", manual_clusters.iter().map(|c| c.to_string())
                                                  .collect::<Vec<String>>() .join(","));
    }
    eprintln!();
}


fn load_graph(gfa_lines: &Vec<String>, print_info: bool) -> (UnitigGraph, Vec<Sequence>) {
    if print_info {
        section_header("Loading graph");
        explanation("The unitig graph is now loaded into memory.");
    }
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_lines(&gfa_lines);
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
        self.collect_clusters(cutoff / 2.0, &vec![], &mut clusters);
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
        if manual_clusters.contains(&self.id) {
            clusters.push(self.id);
        } else if self.distance <= cutoff && !self.has_manual_child(manual_clusters) {
            clusters.push(self.id);
        } else {
            if !self.is_tip() {
                self.left.as_ref().unwrap().collect_clusters(cutoff, manual_clusters, clusters);
                self.right.as_ref().unwrap().collect_clusters(cutoff, manual_clusters, clusters);
            }
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
            format!("({}:{},{}:{}){}", left_str, node.distance - left.distance,
                                       right_str, node.distance - right.distance,
                                       node.id)
        }
        _ => format!("{}", index.get(&node.id).unwrap().string_for_newick()),
    }
}


fn upgma(distances: &HashMap<(u16, u16), f64>, sequences: &mut Vec<Sequence>)
        -> TreeNode {
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
    if let Some(left)  = &mut node.left  { scale_node_distance(left,  scaling_factor); }
    if let Some(right) = &mut node.right { scale_node_distance(right, scaling_factor); }
}









fn generate_clusters(tree: &TreeNode, sequences: &mut Vec<Sequence>,
                     distances: &HashMap<(u16, u16), f64>, cutoff: f64, manual_clusters: &Vec<u16>)
        -> HashMap<u16, ClusterQC> {
    let clusters = if manual_clusters.is_empty() {
        let auto_clusters = tree.automatic_clustering(cutoff);
        refine_auto_clusters(tree, sequences, distances, &auto_clusters)
    } else {
        tree.manual_clustering(cutoff, &manual_clusters)
    };

    // TODO: sanity check that the clusters cover all tips in the tree.

    let mut qc_results = HashMap::new();
    for c in &clusters {
        let dist = tree.max_pairwise_distance(*c);

        // If using manual clustering, then the only reason clusters fail QC is not being included in
        // the user-supplied clusters.
        if manual_clusters.is_empty() {
            if manual_clusters.contains(&c) {
                qc_results.insert(*c, ClusterQC::new_pass(dist));
            } else {
                qc_results.insert(*c, ClusterQC::new_manual_fail(dist));
            }
        }

        // If using automatic clustering, then clusters can fail QC for being contained in other
        // clusters or appearing in too few input assemblies.
        else {
            // TODO
            // TODO
            // TODO
            // TODO
        }
    }


    eprintln!("{:?}", clusters);  // TEMP



    // TODO
    // TODO
    // TODO
    // TODO
    // TODO

    qc_results
}










fn score_clustering(tree: &TreeNode, sequences: &mut Vec<Sequence>,
                    distances: &HashMap<(u16, u16), f64>, clusters: &Vec<u16>) -> f64 {
    // Given a set of node numbers for the tree which define clusters, this function returns the
    // overall score for that clustering (higher is better).


    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO

    0.0  // TEMP
}


fn refine_auto_clusters(tree: &TreeNode, sequences: &mut Vec<Sequence>,
                        distances: &HashMap<(u16, u16), f64>, clusters: &Vec<u16>) -> Vec<u16> {
    // Given a set of node numbers for the tree which define clusters, this function tries to
    // improve the clustering. It does this by trying to split each cluster and seeing if the score
    // gets better, repeating until no improvements can be made.

    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO

    clusters.clone()  // TEMP
}





















fn define_clusters_from_node_numbers(tree: &TreeNode, sequences: &mut Vec<Sequence>,
                                     cluster_nodes: Vec<u16>) -> HashMap<u16, ClusterQC> {
    let mut current_cluster = 0;
    let mut qc_results = HashMap::new();
    for n in cluster_nodes {
        if let Some(node) = find_node_by_id(tree, n) {
            current_cluster += 1;
            assign_cluster_to_node(node, sequences, current_cluster);
            qc_results.insert(current_cluster, ClusterQC::new_pass(node.distance * 2.0));
        } else {
            quit_with_error(&format!("clustering tree does not contain a node with id {}", n));
        }
    }
    for s in sequences.iter_mut() {
        if s.cluster == 0 {  // wasn't assigned a cluster by the above loop
            current_cluster += 1;
            s.cluster = current_cluster;
            let mut qc = ClusterQC::new();
            qc.failure_reasons.push("not included in manual clusters".to_string());
            qc_results.insert(current_cluster, qc);
        }
    }
    let old_to_new = reorder_clusters(sequences);
    let mut reordered_qc_results = HashMap::new();
    for (old_key, qc) in qc_results {
        if let Some(&new_key) = old_to_new.get(&old_key) {
            reordered_qc_results.insert(new_key, qc);
        }
    }
    reordered_qc_results
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


fn define_clusters_automatically(tree: &TreeNode, sequences: &mut Vec<Sequence>,
                                 distances: &HashMap<(u16, u16), f64>, min_assemblies: usize,
                                 cutoff: f64) -> HashMap<u16, ClusterQC> {
    // TODO: It might be nice to rethink this function and cluster more intelligently. Currently,
    //       I'm just using a fixed distance threshold for the entire tree to define clusters. But
    //       one tree might contain clusters with a mix of different ideal thresholds. So instead,
    //       I could start with a high threshold and then explore restricting clusters by moving
    //       towards the tips. A good result would balance both cluster distance (smaller is better)
    //       and the number of assemblies in a cluster (more is better). For example, a more
    //       restrictive cluster that excludes an outlier might be worth it, because while it would
    //       reduce the assembly count by 1 (slightly bad) it could make the cluster much tighter
    //       (good).
    let mut seq_id_to_cluster = HashMap::new();
    let mut current_cluster = 0;
    let mut cluster_distances = HashMap::new();
    define_clusters_from_node(tree, cutoff / 2.0, &mut seq_id_to_cluster, &mut current_cluster,
                              &mut cluster_distances);
    for s in sequences.iter_mut() {
        s.cluster = *seq_id_to_cluster.get(&s.id).unwrap();
    }
    eprintln!("{} cluster{}", current_cluster, match current_cluster { 1 => "", _ => "s" });
    eprintln!();
    let old_to_new = reorder_clusters(sequences);
    let mut reordered_cluster_distances = HashMap::new();
    for (old_key, cluster_distance) in cluster_distances {
        if let Some(&new_key) = old_to_new.get(&old_key) {
            reordered_cluster_distances.insert(new_key, cluster_distance);
        }
    }
    cluster_qc(&sequences, &distances, &reordered_cluster_distances, cutoff, min_assemblies)
}


fn define_clusters_from_node(node: &TreeNode, cutoff: f64,
                             seq_id_to_cluster: &mut HashMap<u16, u16>, current_cluster: &mut u16,
                             cluster_distances: &mut HashMap<u16, f64>) {
    if node.is_tip() || node.distance < cutoff {
        *current_cluster += 1;
        define_clusters(node, seq_id_to_cluster, *current_cluster);
        cluster_distances.insert(*current_cluster, node.distance * 2.0);
    } else {  // node is a loose group, so continue recursively
        define_clusters_from_node(node.left.as_ref().unwrap(), cutoff, seq_id_to_cluster,
                                  current_cluster, cluster_distances);
        define_clusters_from_node(node.right.as_ref().unwrap(), cutoff, seq_id_to_cluster,
                                  current_cluster, cluster_distances);
    }
}


fn define_clusters(node: &TreeNode, seq_id_to_cluster: &mut HashMap<u16, u16>, cluster_number: u16) {
    if node.is_tip() {
        seq_id_to_cluster.insert(node.id, cluster_number);
    } else {
        define_clusters(node.left.as_ref().unwrap(), seq_id_to_cluster, cluster_number);
        define_clusters(node.right.as_ref().unwrap(), seq_id_to_cluster, cluster_number);
    }
}


fn set_min_assemblies(min_assemblies_option: Option<usize>, sequences: &Vec<Sequence>) -> usize {
    // This function automatically sets the --min_assemblies parameter, if the user didn't
    // explicitly supply one. The auto-set value will be one-quarter of the assembly count (rounded)
    // but no less than 1.
    if min_assemblies_option.is_some() {
        return min_assemblies_option.unwrap();
    }
    let assembly_count = get_assembly_count(&sequences);
    let mut min_assemblies = usize_division_rounded(assembly_count, 4);
    if min_assemblies < 1 {
        min_assemblies = 1;
    }
    min_assemblies
}


fn parse_manual_clusters(manual_clusters: Option<String>) -> Vec<u16> {
    if manual_clusters.is_none() {
        return Vec::new();
    }
    let mut clusters: Vec<_> = manual_clusters.unwrap().split(',')
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
    pub fn new() -> Self { Self::default() }
    pub fn new_pass(cluster_dist: f64) -> Self {
        ClusterQC {
            failure_reasons: Vec::new(),
            cluster_dist,
        }
    }
    pub fn new_manual_fail(cluster_dist: f64) -> Self {
        ClusterQC {
            failure_reasons: vec!["not included in manual clusters".to_string()],
            cluster_dist,
        }
    }
    pub fn pass(&self) -> bool { self.failure_reasons.is_empty() }
    pub fn fail(&self) -> bool { !self.failure_reasons.is_empty() }
}


fn cluster_qc(sequences: &Vec<Sequence>, distances: &HashMap<(u16, u16), f64>,
              cluster_distances: &HashMap<u16, f64>, cutoff: f64, min_assemblies: usize)
        -> HashMap<u16, ClusterQC> {
    // This function chooses which clusters pass and which clusters fail. It returns the result as
    // a vector of failure reasons, with an empty vector meaning the cluster passed QC.
    section_header("Cluster QC");
    explanation("Clusters are rejected if they occur in too few assemblies (set by \
                 --min_assemblies) or if they appear to be contained within a larger better \
                 cluster.");
    for s in sequences { assert!(s.cluster != 0); }
    let max_cluster = get_max_cluster(&sequences);
    let mut qc_results: HashMap<_, _> = (1..=max_cluster).map(|c| (c, ClusterQC::new())).collect();
    for c in 1..=max_cluster {
        let assemblies: HashSet<_> = sequences.iter().filter(|s| s.cluster == c).map(|s| s.filename.clone()).collect();
        if assemblies.len() < min_assemblies {
            qc_results.get_mut(&c).unwrap().failure_reasons.push("present in too few assemblies".to_string());
        }
    }
    for c in 1..=max_cluster {
        let container = cluster_is_contained_in_another(c, sequences, distances, cutoff, &qc_results);
        if container > 0 {
            qc_results.get_mut(&c).unwrap().failure_reasons.push(format!("contained within cluster {}", container));
        }
    }
    for c in 1..=max_cluster {
        qc_results.get_mut(&c).unwrap().cluster_dist = *cluster_distances.get(&c).unwrap();
    }
    qc_results
}


fn cluster_is_contained_in_another(cluster_num: u16, sequences: &Vec<Sequence>,
                                   distances: &HashMap<(u16, u16), f64>, cutoff: f64,
                                   qc_results: &HashMap<u16, ClusterQC>) -> u16 {
    // Checks whether this cluster is contained within another cluster that has so-far passed QC.
    // If so, it returns the id of the containing cluster. If not, it returns 0.
    // A cluster counts as contained if the majority of the pairwise comparisons to another cluster
    // are asymmetrical and below the cutoff.
    let passed_clusters: Vec<u16> = qc_results.iter().filter(|(_, q)| q.pass()).map(|(&k, _)| k).collect();
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


fn save_clusters(sequences: &Vec<Sequence>, qc_results: &HashMap<u16, ClusterQC>,
                 clustering_dir: &PathBuf, gfa_lines: &Vec<String>) {
    let pass_dir = clustering_dir.join("qc_pass");
    let fail_dir = clustering_dir.join("qc_fail");
    save_qc_pass_clusters(sequences, qc_results, gfa_lines, &pass_dir);
    save_qc_fail_clusters(sequences, qc_results, gfa_lines, &fail_dir);
}


fn save_qc_pass_clusters(sequences: &Vec<Sequence>, qc_results: &HashMap<u16, ClusterQC>,
                         gfa_lines: &Vec<String>, pass_dir: &PathBuf) {
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
            save_cluster_gfa(&sequences, c, gfa_lines, cluster_dir.join("1_untrimmed.gfa"));
            save_untrimmed_cluster_metrics(seq_lengths, qc.cluster_dist,
                                           cluster_dir.join("1_untrimmed.yaml"));
            eprintln!();
        }
    }
}


fn save_qc_fail_clusters(sequences: &Vec<Sequence>, qc_results: &HashMap<u16, ClusterQC>,
                         gfa_lines: &Vec<String>, fail_dir: &PathBuf) {
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
            save_cluster_gfa(&sequences, c, gfa_lines, cluster_dir.join("1_untrimmed.gfa"));
            save_untrimmed_cluster_metrics(seq_lengths, qc.cluster_dist,
                                           cluster_dir.join("1_untrimmed.yaml"));
            eprintln!();
        }
    }
}


fn save_cluster_gfa(sequences: &Vec<Sequence>, cluster_num: u16, gfa_lines: &Vec<String>, out_gfa: PathBuf) {
    let cluster_seqs: Vec<Sequence> = sequences.iter().filter(|s| s.cluster == cluster_num).cloned().collect();
    let (mut cluster_graph, _) = load_graph(&gfa_lines, false);
    let seq_ids_to_remove:Vec<_> = sequences.iter().filter(|s| s.cluster != cluster_num).map(|s| s.id).collect();
    for id in seq_ids_to_remove {
        cluster_graph.remove_sequence_from_graph(id);
    }
    cluster_graph.remove_zero_depth_unitigs();
    merge_linear_paths(&mut cluster_graph, &cluster_seqs, None);
    cluster_graph.save_gfa(&out_gfa, &cluster_seqs).unwrap();
}


fn save_untrimmed_cluster_metrics(seq_lengths: Vec<usize>, cluster_dist: f64,
                                  cluster_yaml: PathBuf) {
    let metrics = UntrimmedClusterMetrics::new(seq_lengths, cluster_dist);
    save_yaml(&cluster_yaml, &metrics).unwrap();
}


fn save_data_to_tsv(sequences: &Vec<Sequence>, qc_results: &HashMap<u16, ClusterQC>,
                    file_path: &PathBuf) {
    let mut file = File::create(file_path).unwrap();
    write!(file, "node_name\tpassing_clusters\tall_clusters\tsequence_id\tfile_name\tcontig_name\tlength\n").unwrap();
    for seq in sequences {
        assert!(seq.cluster != 0);
        let qc = qc_results.get(&seq.cluster).unwrap();
        let all_cluster = format!("{}", seq.cluster);
        let pass_cluster = if qc.pass() {
            "none".to_string()
        } else {
            format!("{}", seq.cluster)
        };
        write!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", seq.string_for_newick(),
               pass_cluster, all_cluster, seq.id, seq.filename, seq.contig_name(), seq.length).unwrap();
    }
}


fn save_clustering_metrics(clustering_yaml: &PathBuf, sequences: &Vec<Sequence>,
                           qc_results: &HashMap<u16, ClusterQC>) {
    let mut metrics = ClusteringMetrics::new();
    for c in 1..=get_max_cluster(sequences) {
        let qc = qc_results.get(&c).unwrap();
        if qc.pass() { metrics.pass_cluster_count += 1; }
                else { metrics.fail_cluster_count += 1; }
    }
    let mut cluster_filenames = HashMap::new();
    for seq in sequences {
        let qc = qc_results.get(&seq.cluster).unwrap();
        cluster_filenames.entry(seq.cluster).or_insert_with(Vec::new).push(seq.filename.clone());
        if qc.pass() { metrics.pass_contig_count += 1; }
                else { metrics.fail_contig_count += 1; }
    }
    metrics.calculate_fractions();
    metrics.calculate_balance(cluster_filenames);
    save_yaml(&clustering_yaml, &metrics).unwrap();
}


fn reorder_clusters(sequences: &mut Vec<Sequence>) -> HashMap<u16, u16>{
    // Reorder clusters based on their median sequence length (large to small). Returns the mapping
    // of old cluster numbers to new cluster numbers.
    // TODO: perhaps sort not based on sequence length but on UNIQUE sequence length (the sum of 
    //       the lengths of all unitigs). This will prevent doubled plasmids from being longer than
    //       they should be.
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


fn get_assembly_count(sequences: &Vec<Sequence>) -> usize {
    sequences.iter().map(|s| &s.filename).collect::<HashSet<_>>().len()
}


fn get_max_cluster(sequences: &Vec<Sequence>) -> u16 {
    sequences.iter().map(|s| s.cluster).max().unwrap()
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;

    fn assert_almost_eq(a: f64, b: f64, epsilon: f64) {
        assert!((a - b).abs() < epsilon,
                "Numbers are not within {:?} of each other: {} vs {}", epsilon, a, b);
    }

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
        assert_eq!(set_min_assemblies(None, &sequences), 1);  // 5 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 1);  // 4 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 1);  // 3 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 1);  // 2 assemblies
        sequences.pop();
        assert_eq!(set_min_assemblies(None, &sequences), 1);  // 1 assembly
    }

    fn test_tree() -> TreeNode {
        // Generates this tree: (1:0.5,(2:0.3,(3:0.2,(4:0.1,5:0.1):0.1):0.1):0.2);
        let n1 = TreeNode { id: 1, left: None, right: None, distance: 0.0 };
        let n2 = TreeNode { id: 2, left: None, right: None, distance: 0.0 };
        let n3 = TreeNode { id: 3, left: None, right: None, distance: 0.0 };
        let n4 = TreeNode { id: 4, left: None, right: None, distance: 0.0 };
        let n5 = TreeNode { id: 5, left: None, right: None, distance: 0.0 };
        let n6 = TreeNode { id: 6, left: Some(Box::new(n4)), right: Some(Box::new(n5)), distance: 0.1 };
        let n7 = TreeNode { id: 7, left: Some(Box::new(n3)), right: Some(Box::new(n6)), distance: 0.2 };
        let n8 = TreeNode { id: 8, left: Some(Box::new(n2)), right: Some(Box::new(n7)), distance: 0.3 };
        let n9 = TreeNode { id: 9, left: Some(Box::new(n1)), right: Some(Box::new(n8)), distance: 0.5 };
        n9
    }

    #[test]
    fn test_automatic_clustering() {
        let tree = test_tree();
        assert_eq!(tree.automatic_clustering(0.8), vec![1, 8]);
        assert_eq!(tree.automatic_clustering(0.5), vec![1, 2, 7]);
        assert_eq!(tree.automatic_clustering(0.3), vec![1, 2, 3, 6]);
        assert_eq!(tree.automatic_clustering(0.1), vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_has_manual_child() {
        let tree = test_tree();
        assert!(!tree.has_manual_child(&vec![]));
        for n in 1..=9   { assert!( tree.has_manual_child(&vec![n])); }
        for n in 10..=19 { assert!(!tree.has_manual_child(&vec![n])); }
    }

    #[test]
    fn test_manual_clustering() {
        let tree = test_tree();
        assert_eq!(tree.manual_clustering(0.5, &vec![]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &vec![1]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &vec![1, 2]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &vec![1, 2, 7]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.5, &vec![3]), vec![1, 2, 3, 6]);
        assert_eq!(tree.manual_clustering(0.5, &vec![4]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.5, &vec![5]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.5, &vec![1, 2, 3, 4, 5]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.8, &vec![]), vec![1, 8]);
        assert_eq!(tree.manual_clustering(0.8, &vec![1]), vec![1, 8]);
        assert_eq!(tree.manual_clustering(0.8, &vec![2]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.8, &vec![3]), vec![1, 2, 3, 6]);
        assert_eq!(tree.manual_clustering(0.8, &vec![4]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.8, &vec![5]), vec![1, 2, 3, 4, 5]);
        assert_eq!(tree.manual_clustering(0.8, &vec![6]), vec![1, 2, 3, 6]);
        assert_eq!(tree.manual_clustering(0.8, &vec![7]), vec![1, 2, 7]);
        assert_eq!(tree.manual_clustering(0.8, &vec![8]), vec![1, 8]);
    }

    #[test]
    fn test_check_consistency() {
        let tree = test_tree();
        tree.check_consistency(&vec![1, 2, 3, 4, 5]);
        tree.check_consistency(&vec![1, 2, 3, 6]);
        tree.check_consistency(&vec![1, 2, 7]);
        tree.check_consistency(&vec![1, 8]);
        tree.check_consistency(&vec![9]);
        assert!(panic::catch_unwind(|| {
            tree.check_consistency(&vec![5, 6]);
        }).is_err());
        assert!(panic::catch_unwind(|| {
            tree.check_consistency(&vec![6, 8]);
        }).is_err());
        assert!(panic::catch_unwind(|| {
            tree.check_consistency(&vec![1, 9]);
        }).is_err());
    }

    #[test]
    fn test_max_pairwise_distance() {
        let tree = test_tree();
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
    fn test_parse_manual_clusters() {
        assert_eq!(parse_manual_clusters(Some("1,2,3".to_string())), vec![1, 2, 3]);
        assert_eq!(parse_manual_clusters(None), vec![]);
        assert!(panic::catch_unwind(|| {
            parse_manual_clusters(Some("x,y,z".to_string()));
        }).is_err());
        assert!(panic::catch_unwind(|| {
            parse_manual_clusters(Some("^&%^*".to_string()));
        }).is_err());
    }
}
