// Position objects store the sequence, strand and position for contigs in the input assemblies.
// They are used:
// * In KmerGraph objects, where each Kmer object has one or more KmerPos objects.
// * In UnitigGraph objects, where each Unitig has one or more starting and ending UnitigPos
//   objects on both strands. These form a doubly linked list, tracing the input contig through
//   the UnitigGraph.

// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::fmt;


pub struct KmerPos {
    seq_id: u16,
    strand: bool, // true for forward strand, false for reverse strand
    pub pos: u32, // 0-based indexing
}

impl KmerPos {
    pub fn new(seq_id: u16, strand: bool, pos: usize) -> KmerPos {
        KmerPos {
            seq_id,
            strand,
            pos: pos as u32,
        }
    }
}

impl fmt::Display for KmerPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.seq_id, if self.strand { "+" } else { "-" }, self.pos)
    }
}


pub struct UnitigPos {
    seq_id: u32,
    strand: i32, // 1 for forward strand, -1 for reverse strand
    pub pos: usize, // 0-based indexing
    prev: Option<Box<UnitigPos>>,
    next: Option<Box<UnitigPos>>,
    unitig: Option<i32>,  // TODO: change to Unitig
    unitig_strand: Option<i32>, // 1 for forward, -1 for reverse
    unitig_start_end: Option<u32>, // 0 for start, 1 for end
}

impl UnitigPos {
    pub fn new(seq_id: u32, strand: i32, pos: usize, unitig: Option<i32>,
           unitig_strand: Option<i32>, unitig_start_end: Option<u32>) -> UnitigPos {
        UnitigPos {
            seq_id,
            strand,
            pos,
            prev: None,
            next: None,
            unitig,
            unitig_strand,
            unitig_start_end,
        }
    }

    fn on_unitig_start(&self) -> bool {
        self.unitig_start_end == Some(0)
    }

    fn on_unitig_end(&self) -> bool {
        self.unitig_start_end == Some(1)
    }
}

impl fmt::Display for UnitigPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.seq_id, if self.strand == 1 { "+" } else { "-" }, self.pos)
    }
}
