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

use std::collections::{HashMap, HashSet};

use crate::log::{section_header, explanation};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


const GAP: i32 = 0;
const NONE: usize = usize::MAX;


pub fn trim_path_overlap(graph: &mut UnitigGraph, sequences: &Vec<Sequence>, min_identity: f64) -> Vec<Sequence> {
    section_header("Trim path overlaps");
    explanation("Paths for circular replicons may contain start-end overlap. These overlaps are \
                 searched for and trimmed if found.");
    let weights: HashMap<_, _> = graph.unitigs.iter().map(|rc| {let u = rc.borrow(); (u.number as i32, u.length())}).collect();
    let mut clusters: Vec<_> = sequences.iter().map(|s| s.cluster).collect::<HashSet<_>>().into_iter().collect();
    clusters.sort();
    let mut trimmed_sequences = vec![];
    for c in clusters {
        eprintln!("Cluster {}", c);
        let cluster_sequences: Vec<_> = sequences.iter().filter(|s| s.cluster == c).cloned().collect();
        for s in cluster_sequences {
            let path = graph.get_unitig_path_for_sequence(&s);
            let path = path_to_signed_numbers(&path);

            let trimmed_path = trim_path(&path, &weights, min_identity);
            if trimmed_path.is_some() {
                let trimmed_path = trimmed_path.unwrap();
                let trimmed_length: u32 = trimmed_path.iter().filter_map(|&u| Some(weights[&u.abs()])).sum();
                eprintln!("  {} - trimmed to {} bp", s, trimmed_length);

                // TODO: trim the sequence in the unitig graph
                //       first remove the sequence, then add it back in with its smaller path
            } else {
                eprintln!("  {} - not trimmed", s);
                trimmed_sequences.push(s.clone());
            }
        }
        eprintln!();
    }
    trimmed_sequences
}


fn path_to_signed_numbers(path: &Vec<(u32, bool)>) -> Vec<i32> {
    path.iter().map(|p| if p.1 {p.0 as i32} else {-(p.0 as i32)}).collect()
}


fn trim_path(path: &Vec<i32>, weights: &HashMap<i32, u32>, min_identity: f64) -> Option<Vec<i32>> {
    let (alignment_a, alignment_b, alignment_a_indices, alignment_b_indices) = overlap_alignment(&path, &weights, min_identity);
    if alignment_a.is_empty() {
        return None;
    }
    let midpoint_index = find_midpoint(&alignment_a, &alignment_b, weights);
    let start = alignment_a_indices[midpoint_index];
    let end = alignment_b_indices[midpoint_index];
    Some(path[start..end].to_vec())
}


fn overlap_alignment(path: &Vec<i32>, weights: &HashMap<i32, u32>, min_identity: f64) -> (Vec<i32>, Vec<i32>, Vec<usize>, Vec<usize>) {
    let n = path.len();
    let mut scoring_matrix = vec![vec![-f64::INFINITY; n + 1]; n + 1];

    // Initialize the top and left edges to zero.
    for i in 0..=n {
        scoring_matrix[i][0] = 0.0;
        scoring_matrix[0][i] = 0.0;
    }

    // Fill the scoring matrix.
    for i in 1..=n {
        for j in 1..=n {
            if i == j {continue;} // skip diagonal
            let weight_i = *weights.get(&path[i - 1].abs()).unwrap() as f64;
            let weight_j = *weights.get(&path[j - 1].abs()).unwrap() as f64;
            let match_score = scoring_matrix[i - 1][j - 1] + if path[i - 1] == path[j - 1] {
                weight_i
            } else {
                -(weight_i + weight_j) / 2.0
            };
            let delete_score = scoring_matrix[i - 1][j] - weight_i;
            let insert_score = scoring_matrix[i][j - 1] - weight_j;
            scoring_matrix[i][j] = match_score.max(delete_score).max(insert_score);
        }
    }

    // Find the maximum score from the right and bottom edges.
    let mut max_score = f64::NEG_INFINITY;
    let mut max_i = 0;
    let mut max_j = 0;
    for i in 1..=n {
        if scoring_matrix[i][n] > max_score {
            max_score = scoring_matrix[i][n];
            max_i = i;
            max_j = n;
        }
        if scoring_matrix[n][i] > max_score {
            max_score = scoring_matrix[n][i];
            max_i = n;
            max_j = i;
        }
    }

    // A negative score indicates a very poor alignment.
    if max_score <= 0.0 {
        return (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    }

    // Traceback to get the alignment and record indices
    let mut alignment_a = Vec::new();
    let mut alignment_b = Vec::new();
    let mut alignment_a_indices = Vec::new();
    let mut alignment_b_indices = Vec::new();
    let mut i = max_i;
    let mut j = max_j;
    while i > 0 && j > 0 {
        if path[i - 1] == path[j - 1] {
            alignment_a.push(path[i - 1]);
            alignment_b.push(path[j - 1]);
            alignment_a_indices.push(i - 1);
            alignment_b_indices.push(j - 1);
            i -= 1;
            j -= 1;
        } else if scoring_matrix[i - 1][j] >= scoring_matrix[i][j - 1] {
            alignment_a.push(path[i - 1]);
            alignment_b.push(GAP);
            alignment_a_indices.push(i - 1);
            alignment_b_indices.push(NONE);
            i -= 1;
        } else {
            alignment_a.push(GAP);
            alignment_b.push(path[j - 1]);
            alignment_a_indices.push(NONE);
            alignment_b_indices.push(j - 1);
            j -= 1;
        }
    }

    alignment_a.reverse();
    alignment_b.reverse();
    alignment_a_indices.reverse();
    alignment_b_indices.reverse();

    let alignment_a_length: u32 = alignment_a.iter().filter_map(|&u| if u != GAP { Some(weights[&u.abs()]) } else { None }).sum();
    let alignment_b_length: u32 = alignment_b.iter().filter_map(|&u| if u != GAP { Some(weights[&u.abs()]) } else { None }).sum();
    let mean_length = (alignment_a_length as f64 + alignment_b_length as f64) / 2.0;
    let alignment_identity = max_score / mean_length;

    if alignment_identity < min_identity {
        return (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    }

    (alignment_a, alignment_b, alignment_a_indices, alignment_b_indices)
}


fn find_midpoint(alignment_a: &Vec<i32>, alignment_b: &Vec<i32>, weights: &HashMap<i32, u32>) -> usize {
    assert!(alignment_a.len() == alignment_b.len());

    let total_weight = alignment_a.iter().filter_map(|&u| if u != GAP { Some(weights[&u.abs()]) } else { None }).sum::<u32>()
                     + alignment_b.iter().filter_map(|&u| if u != GAP { Some(weights[&u.abs()]) } else { None }).sum::<u32>();
    let mut cumulative_weight = 0;
    let mut best_index = 0;
    let mut best_closeness = 1.0;

    for (i, &unitig_a) in alignment_a.iter().enumerate() {
        let unitig_b = alignment_b[i];
        if unitig_a != GAP {
            cumulative_weight += weights[&unitig_a.abs()];
        }
        if unitig_b != GAP {
            cumulative_weight += weights[&unitig_b.abs()];
        }
        let closeness_to_midpoint = (0.5 - (cumulative_weight as f64 / total_weight as f64)).abs();
        if unitig_a == unitig_b && closeness_to_midpoint < best_closeness {
            best_index = i;
            best_closeness = closeness_to_midpoint;
        }
    }
    best_index
}


#[cfg(test)]
mod tests {
    use maplit::hashmap;
    use super::*;
    use crate::misc::strand;

    #[test]
    fn test_path_to_signed_numbers() {
        let path = vec![(1, strand::FORWARD), (2, strand::FORWARD), (3, strand::REVERSE)];
        assert_eq!(path_to_signed_numbers(&path), vec![1, 2, -3]);

        let path = vec![(4, strand::REVERSE), (5, strand::FORWARD), (6, strand::REVERSE)];
        assert_eq!(path_to_signed_numbers(&path), vec![-4, 5, -6]);

        assert_eq!(path_to_signed_numbers(&vec![]), vec![]);
    }

    #[test]
    fn test_overlap_alignment() {
        // No alignment
        let path = vec![1, -2, 3, -4, 5];
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10};
        let (alignment_a, alignment_b, alignment_a_indices, alignment_b_indices) = overlap_alignment(&path, &weights, 0.9);
        assert_eq!(alignment_a, vec![]);
        assert_eq!(alignment_b, vec![]);
        assert_eq!(alignment_a_indices, vec![]);
        assert_eq!(alignment_b_indices, vec![]);

        // Exact overlap of two unitigs
        let path = vec![1, -2, 3, -4, 5, 1, -2];
        let weights = hashmap!{1 => 10, 2 => 10, 3 => 10, 4 => 10, 5 => 10};
        let (alignment_a, alignment_b, alignment_a_indices, alignment_b_indices) = overlap_alignment(&path, &weights, 0.9);
        assert_eq!(alignment_a, vec![1, -2]);
        assert_eq!(alignment_b, vec![1, -2]);
        assert_eq!(alignment_a_indices, vec![0, 1]);
        assert_eq!(alignment_b_indices, vec![5, 6]);

        // Inexact overlap of three unitigs
        let path = vec![1, -2, 3, -4, 5, 1, 6, 3];
        let weights = hashmap!{1 => 30, 2 => 1, 3 => 10, 4 => 10, 5 => 10, 6 => 1};
        let (alignment_a, alignment_b, alignment_a_indices, alignment_b_indices) = overlap_alignment(&path, &weights, 0.9);
        assert_eq!(alignment_a, vec![1, GAP,  -2, 3]);
        assert_eq!(alignment_b, vec![1,   6, GAP, 3]);
        assert_eq!(alignment_a_indices, vec![0, NONE,    1, 2]);
        assert_eq!(alignment_b_indices, vec![5,    6, NONE, 7]);

        // Same as above but with a high min_alignment_identity resulting in no alignment
        let path = vec![1, -2, 3, -4, 5, 1, 6, 3];
        let weights = hashmap!{1 => 30, 2 => 1, 3 => 10, 4 => 10, 5 => 10, 6 => 1};
        let (alignment_a, alignment_b, alignment_a_indices, alignment_b_indices) = overlap_alignment(&path, &weights, 0.99);
        assert_eq!(alignment_a, vec![]);
        assert_eq!(alignment_b, vec![]);
        assert_eq!(alignment_a_indices, vec![]);
        assert_eq!(alignment_b_indices, vec![]);
    }

    #[test]
    fn test_trim_path_1() {
        let weights = hashmap!{653 => 541, 728 => 413, 757 => 366, 977 => 185, 1010 => 170, 1058 => 153, 1105 => 138, 1133 => 133, 1492 => 79, 1552 => 74, 1637 => 68, 1667 => 65, 1913 => 51, 1943 => 50, 1949 => 50, 1952 => 50, 1967 => 50, 1982 => 50, 1993 => 50, 2012 => 49, 2018 => 48, 2065 => 45, 2070 => 45, 2110 => 42, 2148 => 39, 2276 => 32, 2289 => 32, 2499 => 25, 2640 => 21, 2826 => 15, 2937 => 11, 3148 => 6, 3208 => 5, 3456 => 2, 3578 => 2, 4216 => 1, 4238 => 1, 4575 => 1, 4875 => 1, 4876 => 1, 5191 => 1};

        let path = vec![-653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert!(trimmed_path.is_none());

        let path = vec![-1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert!(trimmed_path.is_none());

        let path = vec![-728, 2640, -4216, 2640, -4216, 2640, 757, -1993, 977, 1492, 1913, -2018, 3456, -4876, 653, -1943, -1058, -2070, 2148, -2012, 2065, -3578, -2110, 3148, -2937, 2276, 5191, -1952, -2499, -1105, 1133, 1949, 1667, -2826, 1010, -1637, -1982, -4875, 2289, 4575, 1552, 4238, -1967, -728, 2640, -4216, 2640, -4216, 2640, 757, -1993, 977, 1492, 1913, -2018, 3456, -4876, 653, -1943, -1058, -2070, 2148, -2012, 2065, -3578, -2110, 3148, -2937, 2276, 5191, -1952, -2499, -1105, 1133, 1949, 1667, -2826, 1010, -1637, -1982, -4875, 2289, 4575, 1552, 4238];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert_eq!(trimmed_path.unwrap(), vec![653, -1943, -1058, -2070, 2148, -2012, 2065, -3578, -2110, 3148, -2937, 2276, 5191, -1952, -2499, -1105, 1133, 1949, 1667, -2826, 1010, -1637, -1982, -4875, 2289, 4575, 1552, 4238, -1967, -728, 2640, -4216, 2640, -4216, 2640, 757, -1993, 977, 1492, 1913, -2018, 3456, -4876]);

        let path = vec![-977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010, 2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, -3208, 2018, -1913];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert_eq!(trimmed_path.unwrap(), vec![2826, -1667, -1949, -1133, 1105, 2499, 1952, -5191, -2276, 2937, -3148, 2110, 3578, -2065, 2012, -2148, 2070, 1058, 1943, -653, 4876, -3456, 2018, -1913, -1492, -977, 1993, -757, -2640, 4216, -2640, 4216, -2640, 728, 1967, -4238, -1552, -4575, -2289, 4875, 1982, 1637, -1010]);
    }

    #[test]
    fn test_trim_path_2() {
        let weights = hashmap!{608 => 662, 645 => 560, 693 => 461, 697 => 453, 710 => 442, 722 => 424, 884 => 228, 892 => 224, 992 => 181, 1028 => 163, 1098 => 140, 1110 => 137, 1125 => 133, 1141 => 131, 1147 => 129, 1177 => 123, 1230 => 111, 1241 => 110, 1242 => 110, 1257 => 109, 1319 => 99, 1322 => 99, 1328 => 98, 1405 => 89, 1427 => 86, 1431 => 86, 1456 => 83, 1457 => 83, 1607 => 71, 1633 => 68, 1646 => 67, 1722 => 62, 1731 => 61, 1751 => 60, 1770 => 59, 1782 => 58, 1792 => 57, 1801 => 57, 1854 => 54, 1871 => 53, 1890 => 52, 1911 => 51, 1959 => 50, 1992 => 50, 2021 => 48, 2041 => 47, 2069 => 45, 2075 => 45, 2128 => 40, 2129 => 40, 2146 => 39, 2152 => 39, 2167 => 38, 2186 => 36, 2204 => 35, 2265 => 32, 2269 => 32, 2296 => 31, 2297 => 31, 2319 => 30, 2384 => 27, 2392 => 27, 2415 => 27, 2453 => 26, 2472 => 26, 2544 => 25, 2572 => 24, 2653 => 21, 2681 => 20, 2723 => 19, 2729 => 19, 2803 => 16, 2817 => 16, 2867 => 14, 2980 => 10, 3088 => 7, 3214 => 5, 3286 => 4, 3291 => 4, 3576 => 2, 3725 => 2, 4233 => 1, 4236 => 1, 4237 => 1, 4241 => 1, 4242 => 1, 4243 => 1, 4244 => 1, 4697 => 1, 5192 => 1};

        let path = vec![-2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert!(trimmed_path.is_none());

        let path = vec![-992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, -2803, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -3725, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -2472];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert_eq!(trimmed_path.unwrap(), vec![-2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, -2803, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453]);

        let path = vec![-2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177, 3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, 1782];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert_eq!(trimmed_path.unwrap(), vec![3291, 608, -1792, -4236, -1731, -1427, 2297, -1319, 4241, -1230, 3088, -2319, 1457, -2392, 645, -3214, -2075, -884, 2653, 4237, 2653, 1405, 4244, -2867, -2186, -1242, 892, -2167, 2544, 1328, 2069, 693, -2384, 1098, 2204, 722, -1722, 2453, -2146, -1241, -1801, 1646, -1633, -1257, -1607, -1028, -1751, -5192, -2129, 2980, -2152, 4242, -2729, -1456, -1992, -2041, -2681, 1770, 1854, -2128, -1890, -4243, 1871, 2269, 1322, -2723, -710, -3286, 1959, 1147, 4233, -4697, 2817, 2265, -2415, -1431, -2572, 1110, 1141, -2021, -1125, -697, 1911, 2296, -992, -3576, -1177]);
    }

    #[test]
    fn test_trim_path_3() {
        let weights = hashmap!{475 => 1316, 535 => 986, 650 => 543, 709 => 443, 734 => 404, 742 => 381, 760 => 363, 767 => 344, 820 => 275, 832 => 265, 878 => 231, 938 => 203, 941 => 202, 1186 => 121, 1190 => 120, 1265 => 108, 1347 => 96, 1368 => 94, 1384 => 92, 1434 => 85, 1473 => 81, 1576 => 73, 1680 => 64, 1920 => 51, 1947 => 50, 1958 => 50, 1970 => 50, 1989 => 50, 2009 => 49, 2082 => 44, 2128 => 40, 2191 => 36, 2233 => 34, 2269 => 32, 2292 => 31, 2301 => 31, 2342 => 29, 2396 => 27, 2431 => 26, 2535 => 25, 2545 => 25, 2582 => 24, 2585 => 24, 2590 => 23, 2593 => 23, 2595 => 23, 2596 => 23, 2602 => 23, 2626 => 22, 2638 => 22, 2641 => 21, 2665 => 21, 2796 => 16, 2814 => 16, 2835 => 15, 2908 => 12, 3075 => 7, 3122 => 6, 3219 => 5, 3223 => 5, 3312 => 4, 3328 => 3, 3351 => 3, 3379 => 3, 3405 => 3, 3483 => 2, 3564 => 2, 3722 => 2, 3808 => 2, 3928 => 1, 3929 => 1, 3930 => 1, 3931 => 1, 3932 => 1, 3933 => 1, 3935 => 1, 3974 => 1, 4217 => 1, 4218 => 1, 4333 => 1, 4334 => 1, 4335 => 1, 4336 => 1, 4337 => 1, 4338 => 1, 4363 => 1, 4559 => 1, 4560 => 1, 4574 => 1, 4637 => 1, 4638 => 1, 4640 => 1, 4958 => 1, 4959 => 1, 4983 => 1, 5177 => 1};

        let path = vec![2396, 760, -4560, 1473, -4574, 1190, 709, 938, -1970, -475, -2269, -1265, 2128, 832, 1958, -1384, -2191, -820, 3808, -1680, -2641, -4217, -2641, 535, -1186, 650, 3328, 2582, 4336, -2292, -3122, -4334, 3933, -2602, -4959, -2593, -3974, 2545, -4363, 767, -3075, -1920, 2596, -3930, 3928, -742, 4640, 1368, 4983, 941, 1947, 878, 4637, 2431, 2665, 3932, 2835, 3929, -4333, 2638, -4958, 1434, -2082, 3931, -3312, 2908, 2233, -2301, 3223, 2814, 3379, -4638, -2796, 4335, 3935, -4337, -4338, 2585, -3351, 1576, -1989, -1347, -4559, 734];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert!(trimmed_path.is_none());

        let path = vec![-3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert_eq!(trimmed_path.unwrap(), vec![1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191]);

        let path = vec![-1190, -3722, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, -3722, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878, -1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186];
        let trimmed_path = trim_path(&path, &weights, 0.95);
        assert_eq!(trimmed_path.unwrap(), vec![-1947, -941, -4983, -1368, -4640, 742, -3928, 3930, -2596, 1920, 3075, -767, 4363, -2545, 3974, 2593, 4959, 2602, -3933, 4334, 3122, 2292, -4336, -2582, -3328, -650, 1186, -535, 2641, 4217, 2641, 1680, -3808, 820, 2191, 1384, -1958, -832, -2128, 1265, 2269, 475, 1970, -938, -709, -1190, 4574, -1473, 4560, -760, -2396, -2590, -4218, -2535, -5177, -734, 4559, 1347, 1989, -1576, 3351, -2585, 4338, 4337, -3935, -4335, 2796, 4638, -3379, -2814, -3223, 2301, -2233, -2908, 3312, -3931, 2082, -1434, 4958, -2638, 4333, -3929, -2835, -3932, -2665, -2431, -4637, -878]);
    }
}
