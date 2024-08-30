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
use seq_io::fastq::Record;
use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::metrics::{ReadSetMetrics, SubsampleMetrics};
use crate::misc::{check_if_dir_is_not_dir, create_dir, fastq_reader, format_float, quit_with_error};


pub fn subsample(fastq_file: PathBuf, out_dir: PathBuf, genome_size: u64, subset_count: u32,
                 min_read_depth: f64, seed: u64) {
    let subsample_yaml = out_dir.join("subsample.yaml");
    check_settings(&out_dir, genome_size, subset_count, min_read_depth);
    create_dir(&out_dir);
    starting_message();
    print_settings(&fastq_file, &out_dir, genome_size, subset_count, min_read_depth, seed);
    let mut metrics = SubsampleMetrics::new();
    let (input_count, input_bases) = input_fastq_stats(&fastq_file, &mut metrics);
    let reads_per_subset = calculate_subsets(input_count, input_bases, genome_size, min_read_depth);
    save_subsets(&fastq_file, subset_count, reads_per_subset, &out_dir, seed, &mut metrics);
    metrics.save_to_yaml(&subsample_yaml);
    finished_message(&out_dir);
}


fn check_settings(out_dir: &PathBuf, genome_size: u64, subset_count: u32, min_read_depth: f64) {
    check_if_dir_is_not_dir(out_dir);
    if genome_size < 1 {
        quit_with_error("--genome_size must be at least 1");
    }
    if subset_count < 1 {
        quit_with_error("--count must be at least 2");
    }
    if min_read_depth <= 0.0 {
        quit_with_error("--min_read_depth must be greater than 0");
    }
}


fn starting_message() {
    section_header("Starting autocycler subsample");
    explanation("This command subsamples a long-read set into subsets that are maximally \
                 independent from each other.");
}


fn print_settings(fastq_file: &PathBuf, out_dir: &PathBuf, genome_size: u64, subset_count: u32,
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


fn input_fastq_stats(fastq_file: &PathBuf, metrics: &mut SubsampleMetrics) -> (u64, u64) {
    let mut read_lengths: Vec<u64> = fastq_reader(fastq_file).records()
        .map(|record| record.expect("Error reading FASTQ file").seq().len() as u64).collect();
    read_lengths.sort_unstable();
    let total_bases = read_lengths.iter().sum();
    let n50_target_bases = total_bases / 2;
    let mut running_total = 0;
    let mut n50 = 0;
    for read_length in &read_lengths {
        running_total += read_length;
        if running_total >= n50_target_bases {
            n50 = *read_length;
            break;
        }
    }
    let total_count = read_lengths.len() as u64;
    eprintln!("Input FASTQ:");
    eprintln!("  Read count: {}", total_count);
    eprintln!("  Read bases: {}", total_bases);
    eprintln!("  Read N50 length: {} bp", n50);
    eprintln!();
    metrics.input_reads = ReadSetMetrics { count: total_count,
                                           bases: total_bases,
                                           n50: n50 };
    (total_count, total_bases)
}


fn calculate_subsets(read_count: u64, read_bases: u64, genome_size: u64, min_depth: f64) -> u64 {
    section_header("Calculating subset size");
    explanation("Autocycler will now calculate the number of reads to put in each subset.");

    let total_depth = read_bases as f64 / genome_size as f64;
    let mean_read_length = (read_bases as f64 / read_count as f64).round() as u64;
    eprintln!("Total read depth: {:.1}×", total_depth);
    eprintln!("Mean read length: {} bp", mean_read_length);
    eprintln!();
    if total_depth < min_depth {
        quit_with_error("Error: input reads are too shallow to subset");
    }
    eprintln!("Calculating subset sizes:");
    eprintln!("  subset_depth = {} * log_2(4 * total_depth / {}) / 2",
              format_float(min_depth), format_float(min_depth));
    let subset_depth = min_depth * (4.0 * total_depth / min_depth).log2() / 2.0;
    eprintln!("               = {:.1}x", subset_depth);
    let subset_ratio = subset_depth / total_depth;
    let reads_per_subset = (subset_ratio * read_count as f64).round() as u64;
    eprintln!("  reads per subset: {}", reads_per_subset);
    eprintln!();
    reads_per_subset
}


fn save_subsets(fastq_file: &PathBuf, subset_count: u32, reads_per_subset: u64, out_dir: &PathBuf,
                seed: u64, metrics: &mut SubsampleMetrics) {
    let mut rng = StdRng::seed_from_u64(seed);

    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
}


fn finished_message(out_dir: &PathBuf) {
    section_header("Finished!");
    eprintln!("Directory with read subsets: {}", out_dir.display());
    eprintln!();
}
