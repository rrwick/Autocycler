// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::misc::check_if_file_exists;
use crate::unitig_graph::UnitigGraph;


pub fn decompress(in_gfa: PathBuf, out_dir: PathBuf) {
    section_header("Starting autocycler decompress");
    explanation("This command will take a compacted De Bruijn graph (made by autocycler \
                 compress), reconstruct the assemblies used to build that graph and save them \
                 in the specified directory.");
    print_settings(&in_gfa, &out_dir);
    let unitig_graph = load_graph(&in_gfa);
}


fn print_settings(in_gfa: &PathBuf, out_dir: &PathBuf) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --out_dir {}", out_dir.display());
    check_if_file_exists(&in_gfa);
}


fn load_graph(in_gfa: &PathBuf) -> UnitigGraph{
    section_header("Loading graph");
    explanation("The compressed sequence graph is now loaded into memory.");
    let unitig_graph = UnitigGraph::from_gfa_file(&in_gfa);
    eprintln!("{} unitigs", unitig_graph.unitigs.len());
    eprintln!("{} links", unitig_graph.link_count);
    unitig_graph
}