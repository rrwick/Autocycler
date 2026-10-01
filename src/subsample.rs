// This file contains the code for the autocycler subsample subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use rand::{rngs::StdRng, SeedableRng};
use rand::seq::SliceRandom;
use seq_io::fastq::Record;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::log::{section_header, explanation};
use crate::metrics::{ReadSetDetails, SubsampleMetrics};
use crate::misc::{check_if_dir_is_not_dir, check_if_file_exists, create_dir, format_float,
                  quit_with_error, read_iter, spinner, unique_read_iter};


pub fn subsample(reads: Vec<PathBuf>, out_dir: PathBuf, genome_size_str: String,
                 subset_count: usize, min_read_depth: f64, seed: u64) {
    let subsample_yaml = out_dir.join("subsample.yaml");
    let genome_size = parse_genome_size(&genome_size_str);
    check_settings(&reads, &out_dir, genome_size, subset_count, min_read_depth);
    create_dir(&out_dir);
    starting_message();
    print_settings(&reads, &out_dir, genome_size, subset_count, min_read_depth, seed);
    let mut metrics = SubsampleMetrics::default();
    let (input_count, input_bases) = input_read_stats(&reads, &mut metrics);
    let reads_per_subset = calculate_subsets(input_count, input_bases, genome_size, min_read_depth);
    save_subsets(&reads, subset_count, input_count, reads_per_subset, &out_dir, seed,
                 &mut metrics);
    metrics.save_to_yaml(&subsample_yaml);
    finished_message();
}


fn check_settings(reads: &[PathBuf], out_dir: &Path, genome_size: u64, subset_count: usize,
                  min_read_depth: f64) {
    for read_file in reads { check_if_file_exists(read_file); }
    check_if_dir_is_not_dir(out_dir);
    if genome_size < 1       { quit_with_error("--genome_size must be at least 1"); }
    if subset_count < 2      { quit_with_error("--count must be at least 2"); }
    if min_read_depth <= 0.0 { quit_with_error("--min_read_depth must be greater than 0"); }
}


fn starting_message() {
    section_header("Starting autocycler subsample");
    explanation("This command subsamples a long-read set into subsets that are maximally \
                 independent from each other.");
}


fn print_settings(reads: &[PathBuf], out_dir: &Path, genome_size: u64, subset_count: usize,
                  min_read_depth: f64, seed: u64) {
    eprintln!("Settings:");
    eprintln!("  --reads {}", reads[0].display());
    for read_file in &reads[1..] {
        eprintln!("          {}", read_file.display());
    }
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!("  --genome_size {genome_size}");
    eprintln!("  --count {subset_count}");
    eprintln!("  --min_read_depth {}", format_float(min_read_depth));
    eprintln!("  --seed {seed}");
    eprintln!();
}


pub fn parse_genome_size(genome_size_str: &str) -> u64 {
    let genome_size_str = genome_size_str.trim().to_lowercase();
    if let Ok(size) = genome_size_str.parse::<f64>() {
        return size.round() as u64;
    }
    let multiplier = match genome_size_str.chars().last() {
        Some('k') => 1_000.0,
        Some('m') => 1_000_000.0,
        Some('g') => 1_000_000_000.0,
        _ => { quit_with_error("cannot interpret genome size"); }
    };
    let number_part = &genome_size_str[..genome_size_str.len() - 1];
    if let Ok(size) = number_part.parse::<f64>() {
        return (size * multiplier).round() as u64;
    }
    quit_with_error("cannot interpret genome size");
}


fn input_read_stats(reads: &[PathBuf], metrics: &mut SubsampleMetrics) -> (usize, u64) {
    let mut read_lengths: Vec<u64> = unique_read_iter(reads).map(|r| r.seq.len() as u64).collect();
    read_lengths.sort_unstable();
    let details = ReadSetDetails::new(&read_lengths);
    metrics.input_read_count = details.count;
    metrics.input_read_bases = details.bases;
    metrics.input_read_n50 = details.n50;
    eprintln!("Input reads:");
    eprintln!("  Read count: {}", details.count);
    eprintln!("  Read bases: {}", details.bases);
    eprintln!("  Read N50 length: {} bp", details.n50);
    eprintln!();
    (details.count, details.bases)
}


fn calculate_subsets(read_count: usize, read_bases: u64, genome_size: u64, min_depth: f64)
        -> usize {
    section_header("Calculating subset size");
    explanation("Autocycler will now calculate the number of reads to put in each subset.");
    let total_depth = read_bases as f64 / genome_size as f64;
    let mean_read_length = (read_bases as f64 / read_count as f64).round() as u64;
    eprintln!("Total read depth: {total_depth:.1}×");
    eprintln!("Mean read length: {mean_read_length} bp");
    eprintln!();
    if total_depth < min_depth {
        quit_with_error("input reads are too shallow to subset");
    }
    eprintln!("Calculating subset sizes:");
    eprintln!("  subset_depth = {} * log_2(4 * total_depth / {}) / 2",
              format_float(min_depth), format_float(min_depth));
    let subset_depth = min_depth * (4.0 * total_depth / min_depth).log2() / 2.0;
    eprintln!("               = {subset_depth:.1}x");
    let subset_ratio = subset_depth / total_depth;
    let reads_per_subset = (subset_ratio * read_count as f64).round() as usize;
    eprintln!("  reads per subset: {reads_per_subset}");
    eprintln!();
    reads_per_subset
}


fn save_subsets(reads: &[PathBuf], subset_count: usize, input_count: usize,
                reads_per_subset: usize, out_dir: &Path, seed: u64,
                metrics: &mut SubsampleMetrics) {
    section_header("Subsetting reads");
    explanation("The reads are now shuffled and grouped into subset files.");
    let mut rng = StdRng::seed_from_u64(seed);
    let mut read_order: Vec<usize> = (0..input_count).collect();
    read_order.shuffle(&mut rng);
    let mut subset_indices = Vec::new();
    let mut subset_files = Vec::new();
    for i in 0..subset_count {
        eprintln!("subset {}:", i+1);
        subset_indices.push(subsample_indices(subset_count, reads_per_subset, &read_order, i));
        let subset_filename = out_dir.join(format!("sample_{:02}.fastq", i + 1));
        eprintln!("  {}", subset_filename.display());
        let subset_file = BufWriter::new(File::create(subset_filename)
            .expect("Failed to create subset file"));
        subset_files.push(subset_file);
        eprintln!();
    }
    let sample_read_lengths = write_subsampled_reads(reads, &subset_indices, &mut subset_files);
    metrics.output_reads.extend(sample_read_lengths.iter().map(ReadSetDetails::new));
}


fn subsample_indices(subset_count: usize, reads_per_subset: usize, read_order: &[usize], i: usize)
        -> HashSet<usize> {
    let input_count = read_order.len();
    let mut subsample_indices = HashSet::new();
    let start = ((i * input_count) as f64 / subset_count as f64).round() as usize;
    let end = start + reads_per_subset;
    if end > input_count {
        let wrapped_end = end - input_count;
        eprintln!("  reads {}-{} and 1-{}", start + 1, input_count, wrapped_end);
        subsample_indices.extend(&read_order[..wrapped_end]);
    } else {
        eprintln!("  reads {}-{}", start + 1, end);
    }
    subsample_indices.extend(&read_order[start..end.min(input_count)]);
    assert_eq!(subsample_indices.len(), reads_per_subset);
    subsample_indices
}


fn write_subsampled_reads(reads: &[PathBuf], subset_indices: &[HashSet<usize>],
                          subset_files: &mut [BufWriter<File>]) -> Vec<Vec<u64>> {
    let mut sample_read_lengths: Vec<Vec<u64>> = vec![Vec::new(); subset_indices.len()];
    let pb = spinner("writing subsampled reads to files...");
    for (read_i, record) in read_iter(reads).enumerate() {
        for (subset_i, (indices, file)) in subset_indices.iter().zip(subset_files.iter_mut()).enumerate() {
            if indices.contains(&read_i) {
                record.write(file).unwrap();
                sample_read_lengths[subset_i].push(record.seq.len() as u64);
            }
        }
    }
    for file in subset_files.iter_mut() { file.flush().unwrap(); }
    for lengths in &mut sample_read_lengths { lengths.sort_unstable(); }
    pb.finish_and_clear();
    sample_read_lengths
}


fn finished_message() {
    section_header("Finished!");
    explanation("You can now assemble each of the subsampled read sets to produce a set of \
                 assemblies for input into Autocycler compress.")
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;

    #[test]
    fn test_parse_genome_size() {
        assert_eq!(parse_genome_size("100"), 100);
        assert_eq!(parse_genome_size("5000"), 5000);
        assert_eq!(parse_genome_size("5000.1"), 5000);
        assert_eq!(parse_genome_size("5000.9"), 5001);
        assert_eq!(parse_genome_size(" 435 "), 435);
        assert_eq!(parse_genome_size("1234567890"), 1234567890);
        assert_eq!(parse_genome_size("12.0k"), 12000);
        assert_eq!(parse_genome_size("47K"), 47000);
        assert_eq!(parse_genome_size("2m"), 2000000);
        assert_eq!(parse_genome_size("13.1M"), 13100000);
        assert_eq!(parse_genome_size("3g"), 3000000000);
        assert_eq!(parse_genome_size("1.23456G"), 1234560000);
        assert!(panic::catch_unwind(|| {
            parse_genome_size("abcd");
        }).is_err());
        assert!(panic::catch_unwind(|| {
            parse_genome_size("12q");
        }).is_err());
        assert!(panic::catch_unwind(|| {
            parse_genome_size("m123");
        }).is_err());
        assert!(panic::catch_unwind(|| {
            parse_genome_size("15kg");
        }).is_err());
    }

    #[test]
    fn test_subsample_indices() {
        let read_order = vec![4, 2, 3, 1, 0, 5];

        assert_eq!(subsample_indices(6, 2, &read_order, 0), HashSet::from([4, 2]));
        assert_eq!(subsample_indices(6, 2, &read_order, 1), HashSet::from([2, 3]));
        assert_eq!(subsample_indices(6, 2, &read_order, 2), HashSet::from([3, 1]));
        assert_eq!(subsample_indices(6, 2, &read_order, 3), HashSet::from([1, 0]));
        assert_eq!(subsample_indices(6, 2, &read_order, 4), HashSet::from([0, 5]));
        assert_eq!(subsample_indices(6, 2, &read_order, 5), HashSet::from([5, 4]));

        assert_eq!(subsample_indices(3, 5, &read_order, 0), HashSet::from([4, 2, 3, 1, 0]));
        assert_eq!(subsample_indices(3, 5, &read_order, 1), HashSet::from([3, 1, 0, 5, 4]));
        assert_eq!(subsample_indices(3, 5, &read_order, 2), HashSet::from([0, 5, 4, 2, 3]));

        assert_eq!(subsample_indices(2, 5, &read_order, 0), HashSet::from([4, 2, 3, 1, 0]));
        assert_eq!(subsample_indices(2, 5, &read_order, 1), HashSet::from([1, 0, 5, 4, 2]));

        assert_eq!(subsample_indices(4, 2, &read_order, 1), HashSet::from([3, 1]));
        assert_eq!(subsample_indices(4, 2, &read_order, 3), HashSet::from([5, 4]));
        assert_eq!(subsample_indices(12, 1, &read_order, 11), HashSet::from([4]));
        assert_eq!(subsample_indices(4, 6, &read_order, 1), HashSet::from([0, 1, 2, 3, 4, 5]));
        assert_eq!(subsample_indices(3, 0, &read_order, 1), HashSet::new());
        assert_eq!(subsample_indices(2, 0, &[], 0), HashSet::new());
    }
}
