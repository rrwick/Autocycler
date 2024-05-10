// This file defines a struct for storing the input assembly sequences and relevant information.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::fmt;

use crate::misc::{quit_with_error, reverse_complement};


pub struct Sequence {
    pub id: u16,
    pub forward_seq: Vec<u8>,
    pub reverse_seq: Vec<u8>,
    pub filename: String,
    pub contig_header: String,
    pub length: usize,
    pub cluster: i32,
}

impl Sequence {
    pub fn new(id: u16, seq: String, filename: String, contig_header: String, length: usize) -> Sequence {

        let forward_seq = seq.into_bytes();
        if !forward_seq.iter().all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
            quit_with_error(&format!("{} contains non-ACGT characters", filename));
        }
        let reverse_seq = reverse_complement(&forward_seq);

        Sequence {
            id,
            forward_seq: forward_seq,
            reverse_seq: reverse_seq,
            filename,
            contig_header,
            length,
            cluster: 0,
        }
    }

    pub fn contig_name(&self) -> String {
        self.contig_header.split_whitespace().next().unwrap_or("").to_string()
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} ({}bp)", self.filename, self.contig_name(), self.length)
    }
}

impl fmt::Debug for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { fmt::Display::fmt(self, f) }
}
