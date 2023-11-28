// This file defines a struct for storing the input assembly sequences and relevant information.

// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use crate::misc::reverse_complement;


pub struct Sequence {
    pub id: u16,
    pub forward_seq: Vec<u8>,
    pub reverse_seq: Vec<u8>,
    pub filename: String,
    pub contig_header: String,
    pub length: usize,
}

impl Sequence {
    pub fn new(id: u16, seq: String, filename: String, contig_header: String, length: usize) -> Sequence {
        Sequence {
            id,
            reverse_seq: reverse_complement(&seq).into_bytes(),
            forward_seq: seq.into_bytes(),
            filename,
            contig_header,
            length,
        }
    }
}
