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
impl InputAssemblyMetrics { pub fn new() -> Self { Self::default() } }


#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ClusteringMetrics {
    pub qc_pass_cluster_count: u32,
    pub qc_pass_contig_count: u32,
    pub qc_fail_cluster_count: u32,
    pub qc_fail_contig_count: u32,
    pub clustered_fraction: f64,
    pub cluster_evenness: f64,
}
impl ClusteringMetrics { pub fn new() -> Self { Self::default() } }


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
