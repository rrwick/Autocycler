// This file contains the code for the autocycler correct subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use crate::log::{section_header, explanation};
use crate::misc::{check_if_dir_exists, check_if_file_exists};

use std::path::PathBuf;


pub fn correct(out_dir: PathBuf) {
    let gfa = out_dir.join("06_repeats_resolved.gfa");
    check_settings(&out_dir, &gfa);
    starting_message();
    print_settings(&out_dir);
    // TODO: load graph
    // TODO: load clustering
    // TODO: remove excluded contigs from the graph
    // TODO: find initial single-copy contigs
    // TODO: expand set of single-copy contigs based on graph structure
    // TODO: find the ordering of single-copy contigs
    // TODO: resolve repeats by duplicating unitigs
}


fn check_settings(out_dir: &PathBuf, gfa: &PathBuf) {
    check_if_dir_exists(&out_dir);
    check_if_file_exists(&gfa);
}


fn starting_message() {
    section_header("Starting autocycler correct");
    explanation("This command corrects errors in the unitig graph.");
}


fn print_settings(out_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --out_dir {}", out_dir.display());
    eprintln!();
}