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
use std::io::Write;
use std::path::Path;

use crate::misc::{median_usize, mad_usize};
use crate::sequence::Sequence;


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct SubsampleMetrics {
    pub input_read_count: usize,
    pub input_read_bases: u64,
    pub input_read_n50: u64,

    // TODO: add input_read_min_length andinput_read_max_length

    pub output_reads: Vec<ReadSetDetails>,
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ReadSetDetails {
    pub count: usize,
    pub bases: u64,
    pub n50: u64,
}

impl ReadSetDetails {
    pub fn new(sorted_read_lengths: &Vec<u64>) -> Self {
        let bases: u64 = sorted_read_lengths.iter().sum();
        let n50_target_bases = bases / 2;
        let mut running_total = 0;
        let mut n50 = 0;
        for read_length in sorted_read_lengths {
            running_total += read_length;
            if running_total >= n50_target_bases {
                n50 = *read_length;
                break;
            }
        }
        ReadSetDetails {
            count: sorted_read_lengths.len(),
            bases,
            n50,
        }
    }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct InputAssemblyMetrics {
    pub input_assemblies_count: u32,
    pub input_assemblies_total_contigs: u32,
    pub input_assemblies_total_length: u64,
    pub compressed_unitig_count: u32,
    pub compressed_unitig_total_length: u64,
    pub input_assembly_details: Vec<InputAssemblyDetails>,
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct InputAssemblyDetails {
    pub filename: String,
    pub contigs: Vec<InputContigDetails>
}

impl InputAssemblyDetails {
    pub fn new(filename: &Path) -> Self {
        InputAssemblyDetails {
            filename: filename.to_string_lossy().to_string(),
            contigs: Vec::new(),
        }
    }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct InputContigDetails {
    pub name: String,
    pub description: String,
    pub length: u64,
}

impl InputContigDetails {
    pub fn new(seq: &Sequence) -> Self {
        InputContigDetails {
            name: seq.contig_name(),
            description: seq.contig_description(),
            length: seq.length as u64,
        }
    }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ClusteringMetrics {
    pub pass_cluster_count: u32,
    pub fail_cluster_count: u32,
    pub pass_contig_count: u32,
    pub fail_contig_count: u32,
    pub pass_contig_fraction: f64,
    pub fail_contig_fraction: f64,
    pub cluster_balance_score: f64,
    pub cluster_tightness_score: f64,
    pub overall_clustering_score: f64,
}

impl ClusteringMetrics {
    pub fn calculate_fractions(&mut self) {
        let total_contigs = self.pass_contig_count + self.fail_contig_count;
        if total_contigs > 0 {
            self.pass_contig_fraction = self.pass_contig_count  as f64 / total_contigs as f64;
            self.fail_contig_fraction = self.fail_contig_count  as f64 / total_contigs as f64;
        }
    }

    pub fn calculate_scores(&mut self, cluster_filenames: HashMap<u16, Vec<String>>,
                            pass_cluster_distances: Vec<f64>) {
        self.calculate_balance(cluster_filenames);
        self.calculate_tightness(pass_cluster_distances);
        self.overall_clustering_score = (self.cluster_balance_score +
                                         self.cluster_tightness_score) / 2.0;
    }

    pub fn calculate_balance(&mut self, cluster_filenames: HashMap<u16, Vec<String>>) {
        // Calculates the balance score for clustering, indicating how evenly filenames are
        // distributed.
        // * For each cluster:
        //   * Count the occurrences of each filename.
        //   * Score each filename based on the count: 0 => 0.0, 1 => 1.0, 2+ => 0.0
        //     (a count of 1 is the best number, i.e. the cluster contains one contig from the file)
        //   * Average (mean) the filename scores to get the cluster score.
        // * The overall balance score is a weighted mean of the cluster scores, weighted by the
        //   number of sequences in the cluster.
        fn count_score(count: u32) -> f64 {
            if count == 0 { 0.0 } else if count == 1 { 1.0 } else { 0.0 }
        }
        let all_filenames: HashSet<String> = cluster_filenames.values()
            .flat_map(|cluster| cluster.iter().cloned()).collect();
        let mut cluster_scores = Vec::new();
        let mut total_weight = 0.0;
        for cluster in cluster_filenames.values() {
            let mut counter = HashMap::new();
            for filename in cluster { *counter.entry(filename).or_insert(0) += 1; }
            let scores: Vec<f64> = all_filenames.iter()
                .map(|filename| count_score(*counter.get(filename).unwrap_or(&0))).collect();
            let cluster_score: f64 = scores.iter().sum::<f64>() / all_filenames.len() as f64;
            cluster_scores.push((cluster_score, cluster.len() as f64));
            total_weight += cluster.len() as f64;
        }
        self.cluster_balance_score = cluster_scores.iter().map(|(score, weight)| score * weight)
            .sum::<f64>() / total_weight;
    }

    pub fn calculate_tightness(&mut self, pass_cluster_distances: Vec<f64>) {
        // Calculates the tightness score for clustering, indicating how tight the QC-pass clusters
        // are. Each QC-pass cluster gets a score using this formula: 
        // https://www.desmos.com/calculator/jfplrrjeyu
        // And the overall tightness score is the mean of these scores.
        if pass_cluster_distances.is_empty() {
            self.cluster_tightness_score = 0.0;
        } else {
            let sum_scores: f64 = pass_cluster_distances.iter()
                .map(|d| 1.0 - d.sqrt()).sum();
            self.cluster_tightness_score = sum_scores / pass_cluster_distances.len() as f64;
        }
    }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct UntrimmedClusterMetrics {
    pub untrimmed_cluster_size: u32,
    pub untrimmed_cluster_lengths: Vec<usize>,
    pub untrimmed_cluster_median: u32,
    pub untrimmed_cluster_mad: u32,
    pub untrimmed_cluster_distance: f64,
}

impl UntrimmedClusterMetrics {
    pub fn new(sequence_lengths: Vec<usize>, untrimmed_cluster_distance: f64) -> Self {
        UntrimmedClusterMetrics {
            untrimmed_cluster_size: sequence_lengths.len() as u32,
            untrimmed_cluster_median: median_usize(&sequence_lengths) as u32,
            untrimmed_cluster_mad: mad_usize(&sequence_lengths) as u32,
            untrimmed_cluster_lengths: sequence_lengths,
            untrimmed_cluster_distance,
        }
    }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct TrimmedClusterMetrics {
    pub trimmed_cluster_size: u32,
    pub trimmed_cluster_lengths: Vec<usize>,
    pub trimmed_cluster_median: u32,
    pub trimmed_cluster_mad: u32,
}

impl TrimmedClusterMetrics {
    pub fn new(sequence_lengths: Vec<usize>) -> Self {
        TrimmedClusterMetrics {
            trimmed_cluster_size: sequence_lengths.len() as u32,
            trimmed_cluster_median: median_usize(&sequence_lengths) as u32,
            trimmed_cluster_mad: mad_usize(&sequence_lengths) as u32,
            trimmed_cluster_lengths: sequence_lengths,
        }
    }
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct CombineMetrics {
    pub consensus_assembly_bases: u64,
    pub consensus_assembly_unitigs: u32,
    pub consensus_assembly_fully_resolved: bool,
    pub consensus_assembly_clusters: Vec<ResolvedClusterDetails>,
}


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ResolvedClusterDetails {
    pub length: u64,
    pub unitigs: u32,
    pub topology: String,
}


// This macro adds some common methods to the metric structs allowing them to be used with
// Autocycler table.
macro_rules! impl_metrics_helpers {
    ($struct_name:ty) => {
        impl $struct_name {
            pub fn save_to_yaml(&self, filename: &Path) {
                let yaml_string = serde_yaml::to_string(&self).unwrap();
                let mut file = File::create(filename).unwrap();
                file.write_all(yaml_string.as_bytes()).unwrap();
            }

            pub fn get_field_names() -> Vec<String> {
                let mut field_names: Vec<String> = match serde_json::to_value(Self::default())
                    .expect("serialisation failed").as_object()
                {
                    Some(map) => map.keys().cloned().collect(),
                    None => Vec::new(),
                };
                field_names.sort();
                field_names
            }
        }
    };
}
impl_metrics_helpers!(SubsampleMetrics);
impl_metrics_helpers!(InputAssemblyMetrics);
impl_metrics_helpers!(ClusteringMetrics);
impl_metrics_helpers!(UntrimmedClusterMetrics);
impl_metrics_helpers!(TrimmedClusterMetrics);
impl_metrics_helpers!(CombineMetrics);


#[cfg(test)]
mod tests {
    use maplit::hashmap;
    use super::*;
    use crate::tests::assert_almost_eq;

    #[test]
    fn test_calculate_balance() {
        // These test cases are ordered by decreasing balance score.
        let filenames_1 = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                   2 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                   3 => vec!["a".to_string(), "b".to_string(), "c".to_string()]};

        let filenames_2 = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                   2 => vec!["a".to_string(), "b".to_string(), "c".to_string(), "a".to_string()],
                                   3 => vec!["a".to_string(), "b".to_string(), "c".to_string()]};

        let filenames_3 = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                   2 => vec!["a".to_string(), "b".to_string(), "c".to_string(), "a".to_string()],
                                   3 => vec!["a".to_string(), "b".to_string()]};

        let filenames_4 = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                   2 => vec!["a".to_string(), "b".to_string(), "c".to_string(), "a".to_string()],
                                   3 => vec!["a".to_string()]};

        let filenames_5 = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                   2 => vec!["a".to_string(), "b".to_string(), "c".to_string(), "a".to_string()],
                                   3 => vec!["a".to_string(), "a".to_string()]};

        let filenames_6 = hashmap!{1 => vec!["a".to_string(), "b".to_string(), "c".to_string()],
                                   2 => vec!["d".to_string(), "e".to_string()],
                                   3 => vec!["f".to_string()]};

        let mut metrics_1 = ClusteringMetrics::default();
        let mut metrics_2 = ClusteringMetrics::default();
        let mut metrics_3 = ClusteringMetrics::default();
        let mut metrics_4 = ClusteringMetrics::default();
        let mut metrics_5 = ClusteringMetrics::default();
        let mut metrics_6 = ClusteringMetrics::default();

        metrics_1.calculate_balance(filenames_1);
        metrics_2.calculate_balance(filenames_2);
        metrics_3.calculate_balance(filenames_3);
        metrics_4.calculate_balance(filenames_4);
        metrics_5.calculate_balance(filenames_5);
        metrics_6.calculate_balance(filenames_6);

        assert_almost_eq(metrics_1.cluster_balance_score, 1.0, 1e-8);
        assert!(metrics_2.cluster_balance_score < metrics_1.cluster_balance_score);
        assert!(metrics_3.cluster_balance_score < metrics_2.cluster_balance_score);
        assert!(metrics_4.cluster_balance_score < metrics_3.cluster_balance_score);
        assert!(metrics_5.cluster_balance_score < metrics_4.cluster_balance_score);
        assert!(metrics_6.cluster_balance_score < metrics_5.cluster_balance_score);
    }

    #[test]
    fn test_get_field_names() {
        assert_eq!(SubsampleMetrics::get_field_names(),
                   vec!["input_read_bases", "input_read_count", "input_read_n50", "output_reads"]);

        assert_eq!(InputAssemblyMetrics::get_field_names(),
                   vec!["compressed_unitig_count", "compressed_unitig_total_length",
                        "input_assemblies_count", "input_assemblies_total_contigs",
                        "input_assemblies_total_length", "input_assembly_details"]);

        assert_eq!(ClusteringMetrics::get_field_names(),
                   vec!["cluster_balance_score", "cluster_tightness_score", "fail_cluster_count",
                        "fail_contig_count", "fail_contig_fraction", "overall_clustering_score",
                        "pass_cluster_count", "pass_contig_count", "pass_contig_fraction"]);

        assert_eq!(UntrimmedClusterMetrics::get_field_names(),
                   vec!["untrimmed_cluster_distance", "untrimmed_cluster_lengths",
                        "untrimmed_cluster_mad", "untrimmed_cluster_median",
                        "untrimmed_cluster_size"]);

        assert_eq!(TrimmedClusterMetrics::get_field_names(),
                   vec!["trimmed_cluster_lengths", "trimmed_cluster_mad", "trimmed_cluster_median",
                        "trimmed_cluster_size"]);

        assert_eq!(CombineMetrics::get_field_names(),
                   vec!["consensus_assembly_bases", "consensus_assembly_clusters",
                        "consensus_assembly_fully_resolved", "consensus_assembly_unitigs"]);
    }
}
