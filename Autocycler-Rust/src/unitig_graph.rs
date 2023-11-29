// This file defines the UnitigGraph struct for building a compacted unitig graph from a KmerGraph.

// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

use crate::kmer_graph::KmerGraph;
use crate::sequence::Sequence;
use crate::unitig::Unitig;


pub struct UnitigGraph {
    pub unitigs: Vec<Unitig>,
    k_size: u32,
}

impl UnitigGraph {
    pub fn from_kmer_graph(kmer_graph: &KmerGraph) -> Self {
        UnitigGraph {
            unitigs: Vec::new(),
            k_size: kmer_graph.k_size,
        }
    }

    pub fn from_gfa_file(gfa_filename: &PathBuf) -> Self {
        UnitigGraph {
            unitigs: Vec::new(),
            k_size: 0,
        }
    }

    pub fn save_gfa(&self, gfa_filename: &PathBuf, sequences: &Vec<Sequence>) -> io::Result<()> {
        let mut file = File::create(gfa_filename)?;
        writeln!(file, "H\tVN:Z:1.0\tKM:i:{}", self.k_size)?;
        for unitig in &self.unitigs {
            writeln!(file, "{}", unitig.gfa_segment_line())?;
        }
        for (a, a_strand, b, b_strand) in self.get_links_for_gfa() {
            writeln!(file, "L\t{}\t{}\t{}\t{}\t0M", a, a_strand, b, b_strand)?;
        }
        for s in sequences {
            writeln!(file, "{}", self.get_gfa_path_line(&s))?;
        }

        Ok(())
    }

    fn get_links_for_gfa(&self) -> Vec<(String, String, String, String)> {
        let mut links = Vec::new();
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
        links
    }

    fn get_gfa_path_line(&self, seq: &Sequence) -> String {
        let path_str = String::new();  // TEMP
        // TODO
        // TODO
        // TODO
        // TODO
        // TODO
        format!("P\t{}\t{}\t*\tLN:i:{}\tFN:Z:{}\tHD:Z:{}",
                seq.id, path_str, seq.length, seq.filename, seq.contig_header)
    }
}
