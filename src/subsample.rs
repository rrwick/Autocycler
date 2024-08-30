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
use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::metrics::SubsampleMetrics;
use crate::misc::{check_if_dir_is_not_dir, create_dir, format_float, quit_with_error};



pub fn subsample(fastq_file: PathBuf, out_dir: PathBuf, genome_size: u64, count: u32,
                 min_read_depth: f64, seed: u64) {
    let subsample_yaml = out_dir.join("subsample.yaml");
    check_settings(&out_dir, genome_size, count, min_read_depth);
    create_dir(&out_dir);
    starting_message();
    print_settings(&fastq_file, &out_dir, genome_size, count, min_read_depth, seed);
    let (read_count, read_bases) = input_fastq_stats(&fastq_file);
    let mut rng = StdRng::seed_from_u64(seed);

    // TODO
    // TODO
    // TODO
    // TODO
    // TODO

    let mut metrics = SubsampleMetrics::new();
    metrics.save_to_yaml(&subsample_yaml);
    finished_message(&out_dir);
}


fn check_settings(out_dir: &PathBuf, genome_size: u64, count: u32, min_read_depth: f64) {
    check_if_dir_is_not_dir(out_dir);
    if genome_size < 1 {
        quit_with_error("--genome_size must be at least 1");
    }
    if count < 1 {
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


fn print_settings(fastq_file: &PathBuf, out_dir: &PathBuf, genome_size: u64, count: u32,
                  min_read_depth: f64, seed: u64) {
    eprintln!("Settings:");
    eprintln!("  --reads {}", fastq_file.display());
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!("  --genome_size {}", genome_size);
    eprintln!("  --count {}", count);
    eprintln!("  --min_read_depth {}", format_float(min_read_depth));
    eprintln!("  --seed {}", seed);
    eprintln!();
}


fn input_fastq_stats(fastq_file: &PathBuf) -> (u64, u64) {

    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO

    (0, 0)  // TEMP
}


fn finished_message(out_dir: &PathBuf) {
    section_header("Finished!");
    eprintln!("Directory with read subsets: {}", out_dir.display());
    eprintln!();
}
