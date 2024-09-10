// This file contains some GFA files for Autocycler's unit tests.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

#[cfg(test)]
pub fn get_test_gfa_1() -> Vec<String> {
    // One connected component
    // Branching structure
    // Includes a hairpin link
    vec!["H\tVN:Z:1.0\tKM:i:9",
         "S\t1\tTTCGCTGCGCTCGCTTCGCTTT\tDP:f:1",
         "S\t2\tTGCCGTCGTCGCTGTGCA\tDP:f:1",
         "S\t3\tTGCCTGAATCGCCTA\tDP:f:1",
         "S\t4\tGCTCGGCTCG\tDP:f:1",
         "S\t5\tCGAACCAT\tDP:f:1",
         "S\t6\tTACTTGT\tDP:f:1",
         "S\t7\tGCCTT\tDP:f:1",
         "S\t8\tATCT\tDP:f:1",
         "S\t9\tGC\tDP:f:1",
         "S\t10\tT\tDP:f:1",
         "L\t1\t+\t4\t+\t0M",
         "L\t4\t-\t1\t-\t0M",
         "L\t1\t+\t5\t-\t0M",
         "L\t5\t+\t1\t-\t0M",
         "L\t2\t+\t1\t+\t0M",
         "L\t1\t-\t2\t-\t0M",
         "L\t3\t-\t1\t+\t0M",
         "L\t1\t-\t3\t+\t0M",
         "L\t4\t+\t7\t-\t0M",
         "L\t7\t+\t4\t-\t0M",
         "L\t4\t+\t8\t+\t0M",
         "L\t8\t-\t4\t-\t0M",
         "L\t6\t-\t5\t-\t0M",
         "L\t5\t+\t6\t+\t0M",
         "L\t6\t+\t6\t-\t0M",
         "L\t7\t-\t9\t+\t0M",
         "L\t9\t-\t7\t+\t0M",
         "L\t8\t+\t10\t-\t0M",
         "L\t10\t+\t8\t-\t0M",
         "L\t9\t+\t7\t+\t0M",
         "L\t7\t-\t9\t-\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_2() -> Vec<String> {
    // One connected component
    // Linear structure with loop-back unitigs on the ends
    vec!["H\tVN:Z:1.0\tKM:i:9",
         "S\t1\tACCGCTGCGCTCGCTTCGCTCT\tDP:f:1",
         "S\t2\tATGAT\tDP:f:1",
         "S\t3\tGCGC\tDP:f:1",
         "L\t1\t+\t2\t+\t0M",
         "L\t2\t-\t1\t-\t0M",
         "L\t1\t+\t2\t-\t0M",
         "L\t2\t+\t1\t-\t0M",
         "L\t1\t-\t3\t+\t0M",
         "L\t3\t-\t1\t+\t0M",
         "L\t1\t-\t3\t-\t0M",
         "L\t3\t+\t1\t+\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_3() -> Vec<String> {
    // One connected component
    // Includes a loop
    // Includes a hairpin link
    vec!["H\tVN:Z:1.0\tKM:i:9",
         "S\t1\tTTCGCTGCGCTCGCTTCGCTTT\tDP:f:1",
         "S\t2\tTGCCGTCGTCGCTGTGCA\tDP:f:1",
         "S\t3\tTGCCTGAATCGCCTA\tDP:f:1",
         "S\t4\tGCTCGGCTCG\tDP:f:1",
         "S\t5\tCGAACCAT\tDP:f:1",
         "S\t6\tTACTTGT\tDP:f:1",
         "S\t7\tGCCTT\tDP:f:1",
         "L\t1\t+\t2\t-\t0M",
         "L\t2\t+\t1\t-\t0M",
         "L\t2\t-\t3\t+\t0M",
         "L\t3\t-\t2\t+\t0M",
         "L\t3\t+\t4\t+\t0M",
         "L\t4\t-\t3\t-\t0M",
         "L\t4\t+\t5\t-\t0M",
         "L\t5\t+\t4\t-\t0M",
         "L\t5\t-\t5\t+\t0M",
         "L\t3\t+\t6\t+\t0M",
         "L\t6\t-\t3\t-\t0M",
         "L\t6\t+\t7\t-\t0M",
         "L\t7\t+\t6\t-\t0M",
         "L\t7\t-\t6\t+\t0M",
         "L\t6\t-\t7\t+\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_4() -> Vec<String> {
    // Two connected components
    // Both form simple loops of multiple unitigs
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tACGACTACGAGCACG\tDP:f:1",
         "S\t2\tTACGACGACGACT\tDP:f:1",
         "S\t3\tACTGACT\tDP:f:1",
         "S\t4\tGCTCG\tDP:f:1",
         "S\t5\tCAC\tDP:f:1",
         "L\t1\t+\t2\t-\t0M",
         "L\t2\t+\t1\t-\t0M",
         "L\t2\t-\t3\t+\t0M",
         "L\t3\t-\t2\t+\t0M",
         "L\t3\t+\t1\t+\t0M",
         "L\t1\t-\t3\t-\t0M",
         "L\t4\t+\t5\t-\t0M",
         "L\t5\t+\t4\t-\t0M",
         "L\t5\t-\t4\t+\t0M",
         "L\t4\t-\t5\t+\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_5() -> Vec<String> {
    // Four connected components
    // One simple loop of one unitig
    // One isolated linear unitig
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "S\t2\tAGCATCAGCATCAGC\tDP:f:1",
         "S\t3\tGTCGCATTT\tDP:f:1",
         "S\t4\tTCGCGAA\tDP:f:1",
         "S\t5\tTTAAAC\tDP:f:1",
         "S\t6\tCACA\tDP:f:1",
         "L\t1\t+\t5\t+\t0M",
         "L\t5\t-\t1\t-\t0M",
         "L\t1\t+\t5\t-\t0M",
         "L\t5\t+\t1\t-\t0M",
         "L\t3\t-\t6\t-\t0M",
         "L\t6\t+\t3\t+\t0M",
         "L\t4\t+\t4\t+\t0M",
         "L\t4\t-\t4\t-\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_6() -> Vec<String> {
    // One connected component
    // Two unitigs on opposite strands
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "S\t2\tAGCATCAGCATCAGC\tDP:f:1",
         "L\t1\t+\t2\t-\t0M",
         "L\t2\t+\t1\t-\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_7() -> Vec<String> {
    // One connected component
    // Two unitigs on opposite strands
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "S\t2\tAGCATCAGCATCAGC\tDP:f:1",
         "L\t1\t-\t2\t+\t0M",
         "L\t2\t-\t1\t+\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_8() -> Vec<String> {
    // One connected component
    // Single circular unitig
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "L\t1\t+\t1\t+\t0M",
         "L\t1\t-\t1\t-\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_9() -> Vec<String> {
    // One connected component
    // Single linear unitig, no links
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_10() -> Vec<String> {
    // One connected component
    // Single linear unitig, hairpin on both ends
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "L\t1\t+\t1\t-\t0M",
         "L\t1\t-\t1\t+\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_11() -> Vec<String> {
    // One connected component
    // Single linear unitig, hairpin on one end
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "L\t1\t+\t1\t-\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_12() -> Vec<String> {
    // One connected component
    // Single linear unitig, hairpin on one end
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "L\t1\t-\t1\t+\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_13() -> Vec<String> {
    // One connected component
    // Single linear unitig, circular connection and a hairpin (odd case)
    vec!["H\tVN:Z:1.0\tKM:i:3",
         "S\t1\tAGCATCGACATCGACTACG\tDP:f:1",
         "L\t1\t+\t1\t+\t0M",
         "L\t1\t-\t1\t-\t0M",
         "L\t1\t-\t1\t+\t0M"].into_iter().map(String::from).collect()
}


#[cfg(test)]
pub fn get_test_gfa_14() -> Vec<String> {
    // Cluster 2 from the toy example before merging linear paths.
    vec!["H\tVN:Z:1.0\tKM:i:13",
         "S\t5\tTGCTCAAAGCCTCGTATTGAG\tDP:f:4.00",
         "S\t8\tGCAGTTCAATCCAATAA\tDP:f:4.00",
         "S\t12\tCATTCGTAACTTGCA\tDP:f:3.00",
         "S\t17\tCCAACGTGTACT\tDP:f:4.00",
         "S\t18\tGGAGTTAGCTTC\tDP:f:4.00",
         "S\t19\tAAGTAGGCG\tDP:f:4.00",
         "S\t21\tGTTTAG\tDP:f:3.00",
         "S\t22\tATACC\tDP:f:3.00",
         "S\t27\tAT\tDP:f:1.00",
         "S\t34\tG\tDP:f:3.00",
         "S\t36\tT\tDP:f:3.00",
         "S\t37\tT\tDP:f:3.00",
         "S\t38\tT\tDP:f:1.00",
         "L\t5\t+\t34\t-\t0M",
         "L\t5\t+\t38\t-\t0M",
         "L\t5\t-\t12\t+\t0M",
         "L\t8\t+\t22\t+\t0M",
         "L\t8\t-\t19\t-\t0M",
         "L\t12\t+\t21\t-\t0M",
         "L\t12\t-\t5\t+\t0M",
         "L\t17\t+\t22\t-\t0M",
         "L\t17\t-\t27\t+\t0M",
         "L\t17\t-\t36\t+\t0M",
         "L\t18\t+\t27\t-\t0M",
         "L\t18\t+\t36\t-\t0M",
         "L\t18\t-\t34\t+\t0M",
         "L\t18\t-\t38\t+\t0M",
         "L\t19\t+\t8\t+\t0M",
         "L\t19\t-\t37\t-\t0M",
         "L\t21\t+\t12\t-\t0M",
         "L\t21\t-\t37\t+\t0M",
         "L\t22\t+\t17\t-\t0M",
         "L\t22\t-\t8\t-\t0M",
         "L\t27\t+\t18\t-\t0M",
         "L\t27\t-\t17\t+\t0M",
         "L\t34\t+\t5\t-\t0M",
         "L\t34\t-\t18\t+\t0M",
         "L\t36\t+\t18\t-\t0M",
         "L\t36\t-\t17\t+\t0M",
         "L\t37\t+\t19\t+\t0M",
         "L\t37\t-\t21\t+\t0M",
         "L\t38\t+\t5\t-\t0M",
         "L\t38\t-\t18\t+\t0M",
         "P\t2\t8+,22+,17-,27+,18-,34+,5-,12+,21-,37+,19+\t*\tLN:i:101\tFN:Z:a.fasta\tHD:Z:a_2\tCL:i:2",
         "P\t4\t5+,38-,18+,36-,17+,22-,8-,19-,37-,21+,12-,5+,34-,18+,36-,17+,22-,8-,19-\t*\tLN:i:178\tFN:Z:b.fasta\tHD:Z:b_2\tCL:i:2",
         "P\t7\t17-,36+,18-,34+,5-,12+,21-,37+,19+,8+\t*\tLN:i:95\tFN:Z:d.fasta\tHD:Z:d_2\tCL:i:2"].into_iter().map(String::from).collect()
}
