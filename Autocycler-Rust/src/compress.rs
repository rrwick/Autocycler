// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use crate::log::{section_header, explanation};
use crate::misc::{find_all_assemblies};

use std::path::PathBuf;


pub fn compress(in_dir: PathBuf, out_gfa: PathBuf, kmer: u32) {
    section_header("Starting autocycler compress");
    explanation("This command will find all assemblies in the given input directory and compress \
                 them into a compacted De Bruijn graph. This graph can then be used to recover \
                 the assemblies (with autocycler decompress) or generate a consensus assembly \
                 (with autocycler resolve).");
    print_settings(&in_dir, &out_gfa, &kmer);

    let assemblies = find_all_assemblies(&in_dir);
}


fn print_settings(in_dir: &PathBuf, out_gfa: &PathBuf, kmer: &u32) {
    eprintln!("Settings:");
    eprintln!("  --in_dir {}", in_dir.display());
    eprintln!("  --out_gfa {}", out_gfa.display());
    eprintln!("  --kmer {}", kmer);
}
