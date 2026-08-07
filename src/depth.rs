// This file contains the code for setting consentig depths using reads. It is used by the
// autocycler combine subcommand when the user supplies reads with the --reads option.

// Copyright 2026 Ryan Wick (rrwick@gmail.com)
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
use crate::unitig_graph::UnitigGraph;


pub fn set_read_depths(graphs: &[&UnitigGraph], reads: &[PathBuf], k_size: u32) -> DepthInfo {
    // Sets the depth of each tig in the given graphs using the given reads. The graphs are all
    // handled together (not one at a time) because repeats must be found across the entire
    // consensus assembly, not just within a single cluster.
    section_header("Setting depths from reads");
    explanation("The reads are now used to set the depth of each sequence in the consensus \
                 assembly, replacing the input-assembly-count depths from the cluster graphs.");

    // TODO: build a table of the assembly's k-mers, masking those which occur more than once.
    // TODO: count the reads' k-mers.
    // TODO: set each unitig's depth from its k-mer counts.

    let _ = (reads, k_size);
    for graph in graphs {
        for unitig in &graph.unitigs {
            unitig.borrow_mut().depth = 1.0;
        }
    }
    DepthInfo { mean_depth: 1.0 }
}


#[derive(Debug, Default)]
pub struct DepthInfo {
    // Assembly-wide values which will go into the log and the metrics file. More will be added
    // here later, e.g. the implied read accuracy and how well the reads are explained by the
    // consensus assembly.
    pub mean_depth: f64, // length-weighted mean depth of the consensus assembly
}
