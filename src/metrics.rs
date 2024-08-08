// This file contains the code for reading and writing Autocycler's YAML files of metrics.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::PathBuf;


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct InputAssemblyMetrics {
    pub assembly_count: u32,
    pub total_contigs: u32,
    pub total_length: u64,
    pub compressed_unitig_count: u32,
    pub compressed_unitig_total_length: u64,
}
impl InputAssemblyMetrics {
    pub fn new() -> Self { Self::default() }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ClusteringMetrics {
    pub qc_pass_cluster_count: u32,
    pub qc_fail_cluster_count: u32,
    pub qc_pass_contig_count: u32,
    pub qc_fail_contig_count: u32,
    pub qc_pass_contig_fraction: f64,
    pub qc_fail_contig_fraction: f64,
    pub cluster_balance: f64,
}
impl ClusteringMetrics {
    pub fn new() -> Self { Self::default() }

    pub fn calculate_fractions(&mut self) {
        let total_contigs = self.qc_pass_contig_count + self.qc_fail_contig_count;
        if total_contigs > 0 {
            self.qc_pass_contig_fraction = self.qc_pass_contig_count  as f64 / total_contigs as f64;
            self.qc_fail_contig_fraction = self.qc_fail_contig_count  as f64 / total_contigs as f64;
        }
    }

    pub fn calculate_balance(&mut self, all_filenames: HashSet<String>,
                              pass_cluster_filenames: HashMap<u16, Vec<String>>) {
        // The cluster balance score ranges from 0 to 1, with 0 being worst and 1 being best. It
        // quantifies how well each assembly is represented in each cluster, with the ideal case
        // being that each assembly contributes one contig to each cluster.
        // 1. For each cluster, count the occurrences of each filename.
        // 2. Calculate the proportion of unique filenames present in each cluster relative to all
        //    unique filenames.
        // 3. Penalize clusters for repeated filenames by computing a penalty ratio.
        // 4. Compute a score for each cluster based on the proportion of unique filename and
        //    penalty ratio.
        // 5. The overall balance metric is the average score of all clusters.
        let mut cluster_scores = Vec::new();
        let total_unique_labels = all_filenames.len() as f64;
        for cluster in pass_cluster_filenames.values() {
            let mut label_counter = HashMap::new();
            for item in cluster { *label_counter.entry(item).or_insert(0) += 1; }
            let unique_labels_in_cluster = label_counter.len() as f64;
            let proportion_unique = unique_labels_in_cluster / total_unique_labels;
            let mut penalty = 0.0;
            for &count in label_counter.values() {
                if count > 1 { penalty += (count - 1) as f64; }
            }
            let penalty_ratio = penalty / cluster.len() as f64;
            let cluster_score = proportion_unique * (1.0 - penalty_ratio);
            cluster_scores.push(cluster_score);
        }
        self.cluster_balance = cluster_scores.iter().sum::<f64>() / cluster_scores.len() as f64;
    }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ClusterMetrics {
    pub cluster_tightness: f64,
    pub cluster_mad: f64,
    pub conflicting_bridges: u32,
    pub resolved_unitigs: u32,
}
impl ClusterMetrics { pub fn new() -> Self { Self::default() } }


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ConsensusMetrics {
    pub total_length: u32,
    pub component_count: u32,
    pub sequence_count: u32,
    pub component_lengths: Vec<u32>,
    pub sequence_lengths: Vec<u32>,
    pub dead_ends: u32,
    pub overall_score: f64,
}
impl ConsensusMetrics { pub fn new() -> Self { Self::default() } }


pub fn save_yaml<T: Serialize>(yaml_filename: &PathBuf, data: T) -> io::Result<()> {
    let yaml_string = serde_yaml::to_string(&data).unwrap();
    let mut file = File::create(yaml_filename)?;
    file.write_all(yaml_string.as_bytes())?;
    Ok(())
}


#[cfg(test)]
mod tests {
    use maplit::hashmap;
    use super::*;

    fn assert_almost_eq(a: f64, b: f64, epsilon: f64) {
        assert!((a - b).abs() < epsilon,
                "Numbers are not within {:?} of each other: {} vs {}", epsilon, a, b);
    }

    #[test]
    fn test_calculate_balance() {
        let mut metrics = ClusteringMetrics::new();
        let pass_cluster_filenames = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                              2 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                              3 => vec!["a".to_string(), "b".to_string(), "c".to_string()]};
        let all_filenames: HashSet<String> = pass_cluster_filenames.values()
            .flat_map(|cluster| cluster.iter().cloned()).collect();
        metrics.calculate_balance(all_filenames, pass_cluster_filenames);
        assert_almost_eq(metrics.cluster_balance, 1.0, 1e-8);

        let mut metrics = ClusteringMetrics::new();
        let pass_cluster_filenames = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                              2 => vec!["a".to_string(), "b".to_string(), "c".to_string(), "a".to_string()],
                                              3 => vec!["a".to_string(), "b".to_string()]};
        let all_filenames: HashSet<String> = pass_cluster_filenames.values()
            .flat_map(|cluster| cluster.iter().cloned()).collect();
        metrics.calculate_balance(all_filenames, pass_cluster_filenames);
        assert_almost_eq(metrics.cluster_balance, 0.8055555555555555, 1e-8);

        let mut metrics = ClusteringMetrics::new();
        let pass_cluster_filenames = hashmap!{1 => vec!["a".to_string(), "b".to_string()],
                                              2 => vec!["c".to_string(), "d".to_string()],
                                              3 => vec!["e".to_string(), "f".to_string()],
                                              4 => vec!["g".to_string(), "h".to_string()]};
        let all_filenames: HashSet<String> = pass_cluster_filenames.values()
            .flat_map(|cluster| cluster.iter().cloned()).collect();
        metrics.calculate_balance(all_filenames, pass_cluster_filenames);
        assert_almost_eq(metrics.cluster_balance, 0.25, 1e-8);
    }
}