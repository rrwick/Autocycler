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

#![allow(clippy::needless_range_loop)]

use rand::{rngs::StdRng, SeedableRng};
use rand::seq::SliceRandom;
use seq_io::fastq::Record;
use std::collections::HashSet;
use std::fs::File;
use std::path::{Path, PathBuf};

use crate::log::{section_header, explanation};
use crate::metrics::{ReadSetDetails, SubsampleMetrics};
use crate::misc::{check_if_dir_is_not_dir, check_if_file_exists, create_dir, fastq_reader,
                  format_float, quit_with_error, spinner};


pub fn subsample(fastq_file: PathBuf, out_dir: PathBuf, genome_size_str: String,
                 subset_count: usize, min_read_depth: f64, seed: u64) {
    let subsample_yaml = out_dir.join("subsample.yaml");
    let genome_size = parse_genome_size(&genome_size_str);
    check_settings(&fastq_file, &out_dir, genome_size, subset_count, min_read_depth);
    create_dir(&out_dir);
    starting_message();
    print_settings(&fastq_file, &out_dir, genome_size, subset_count, min_read_depth, seed);

    // TODO: add automatic genome size estimation

    let mut metrics = SubsampleMetrics::default();
    let (input_count, input_bases) = input_fastq_stats(&fastq_file, &mut metrics);
    let reads_per_subset = calculate_subsets(input_count, input_bases, genome_size, min_read_depth);
    save_subsets(&fastq_file, subset_count, input_count, reads_per_subset, &out_dir, seed,
                 &mut metrics);
    metrics.save_to_yaml(&subsample_yaml);
    finished_message();
}


fn check_settings(fastq_file: &Path, out_dir: &Path, genome_size: u64, subset_count: usize,
                  min_read_depth: f64) {
    check_if_file_exists(fastq_file);
    check_if_dir_is_not_dir(out_dir);
    if genome_size < 1 {       quit_with_error("--genome_size must be at least 1"); }
    if subset_count < 1 {      quit_with_error("--count must be at least 2"); }
    if min_read_depth <= 0.0 { quit_with_error("--min_read_depth must be greater than 0"); }
}


fn starting_message() {
    section_header("Starting autocycler subsample");
    explanation("This command subsamples a long-read set into subsets that are maximally \
                 independent from each other.");
}


fn print_settings(fastq_file: &Path, out_dir: &Path, genome_size: u64, subset_count: usize,
                  min_read_depth: f64, seed: u64) {
    eprintln!("Settings:");
    eprintln!("  --reads {}", fastq_file.display());
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!("  --genome_size {}", genome_size);
    eprintln!("  --count {}", subset_count);
    eprintln!("  --min_read_depth {}", format_float(min_read_depth));
    eprintln!("  --seed {}", seed);
    eprintln!();
}


fn parse_genome_size(genome_size_str: &str) -> u64 {
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


fn input_fastq_stats(fastq_file: &Path, metrics: &mut SubsampleMetrics) -> (usize, u64) {
    let mut read_lengths: Vec<u64> = fastq_reader(fastq_file).records()
        .map(|record| record.expect("Error reading FASTQ file").seq().len() as u64).collect();
    read_lengths.sort_unstable();
    let details = ReadSetDetails::new(&read_lengths);
    metrics.input_read_count = details.count;
    metrics.input_read_bases = details.bases;
    metrics.input_read_n50 = details.n50;
    eprintln!("Input FASTQ:");
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
    eprintln!("Total read depth: {:.1}×", total_depth);
    eprintln!("Mean read length: {} bp", mean_read_length);
    eprintln!();
    if total_depth < min_depth {
        quit_with_error("input reads are too shallow to subset");
    }
    eprintln!("Calculating subset sizes:");
    eprintln!("  subset_depth = {} * log_2(4 * total_depth / {}) / 2",
              format_float(min_depth), format_float(min_depth));
    let subset_depth = min_depth * (4.0 * total_depth / min_depth).log2() / 2.0;
    eprintln!("               = {:.1}x", subset_depth);
    let subset_ratio = subset_depth / total_depth;
    let reads_per_subset = (subset_ratio * read_count as f64).round() as usize;
    eprintln!("  reads per subset: {}", reads_per_subset);
    eprintln!();
    reads_per_subset
}


fn save_subsets(input_fastq: &Path, subset_count: usize, input_count: usize,
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
        let subset_file = File::create(subset_filename).expect("Failed to create subset file");
        subset_files.push(subset_file);
        eprintln!();
    }
    let sample_read_lengths = write_subsampled_reads(input_fastq, subset_count, &subset_indices,
                                                     &mut subset_files);
    for i in 0..subset_count {
        metrics.output_reads.push(ReadSetDetails::new(&sample_read_lengths[i]));
    }
}


fn subsample_indices(subset_count: usize, reads_per_subset: usize, read_order: &[usize], i: usize)
        -> HashSet<usize> {
    // For a given subsample (index i), this function returns a HashSet of the read indices which
    // will go in that subsample.
    let input_count = read_order.len();
    let mut subsample_indices = HashSet::new();
    let start_1 = ((i * input_count) as f64 / subset_count as f64).round() as usize;
    let mut end_1 = start_1 + reads_per_subset;
    if end_1 > input_count {
        let start_2 = 0;
        let end_2 = end_1 - input_count;
        end_1 = input_count;
        eprintln!("  reads {}-{} and {}-{}", start_1 + 1, end_1, start_2 + 1, end_2);
        for j in start_2..end_2 {
            subsample_indices.insert(read_order[j]);
        }
    } else {
        eprintln!("  reads {}-{}", start_1 + 1, end_1);
    }
    for j in start_1..end_1 {
        subsample_indices.insert(read_order[j]);
    }
    assert_eq!(subsample_indices.len(), reads_per_subset);
    subsample_indices
}


fn write_subsampled_reads(input_fastq: &Path, subset_count: usize,
                          subset_indices: &[HashSet<usize>], subset_files: &mut [File])
        -> Vec<Vec<u64>> {
    // This function loops through the input reads, and saves each read to the appropriate output
    // file. It also gathers up and returns the sorted read lengths for each subsampled read set.
    let mut sample_read_lengths: Vec<Vec<u64>> = vec![Vec::new(); subset_count];
    let pb = spinner("writing subsampled reads to files...");
    let mut read_i = 0;
    let mut reader = fastq_reader(input_fastq);
    while let Some(record) = reader.next() {
        let record = record.expect("Error reading FASTQ file");
        for subset_i in 0..subset_count {
            if subset_indices[subset_i].contains(&read_i) {
                record.write(&subset_files[subset_i]).unwrap();
                sample_read_lengths[subset_i].push(record.seq().len() as u64);
            }
        }
        read_i += 1;
    }
    for i in 0..subset_count {
        sample_read_lengths[i].sort_unstable();
    }
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
    }
}
