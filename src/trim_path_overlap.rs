// This file contains functions related to trimming paths to remove start-end overlap.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::collections::HashSet;

use crate::log::{section_header, explanation};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


pub fn trim_path_overlap(graph: &mut UnitigGraph, sequences: &Vec<Sequence>,
                         min_alignment_identity: f64) -> Vec<Sequence> {
    section_header("Trim path overlaps");
    explanation("Paths for circular replicons may contain start-end overlap. These overlaps are \
                 searched for and trimmed if found.");

    let mut clusters: Vec<_> = sequences.iter().map(|s| s.cluster).collect::<HashSet<_>>().into_iter().collect();
    clusters.sort();
    for c in clusters {
        eprintln!("Cluster {}", c); // TEMP
        let cluster_sequences: Vec<_> = sequences.iter().filter(|s| s.cluster == c).cloned().collect();
        let cluster_lengths: Vec<_> = cluster_sequences.iter().map(|s| s.length).collect();
        eprintln!("  untrimmed lengths: {:?}", cluster_lengths); // TEMP
        for s in cluster_sequences {
            let sequence_path = graph.get_unitig_path_for_sequence(&s);
            eprintln!("  {:?}", sequence_path); // TEMP
        }
        eprintln!();
    }

    let sequences = sequences.iter().filter(|s| s.cluster != -1).cloned().collect();  // TEMP

    eprintln!();
    sequences
}
