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

use std::path::PathBuf;


pub fn resolve(in_gfa: PathBuf, out_dir: PathBuf) {
    section_header("Starting autocycler resolve");
    explanation("This command will take a compacted De Bruijn graph (made by autocycler \
                 compress) and simplify it into a consensus assembly. It does this by resolving \
                 repeats and removing errors from the graph.");
    print_settings(&in_gfa, &out_dir);
    // TODO: load graph
    // TODO: load clustering
    // TODO: remove excluded contigs from the graph
    // TODO: find initial single-copy contigs
    // TODO: expand set of single-copy contigs based on graph structure
    // TODO: find the ordering of single-copy contigs
    // TODO: resolve repeats by duplicating unitigs
}


fn print_settings(in_gfa: &PathBuf, out_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_dir {}", out_dir.display());
}