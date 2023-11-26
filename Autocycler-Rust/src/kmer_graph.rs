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

use crate::position::Position;

pub struct Kmer<'a> {
    seq: &'a str,
    positions: Vec<Position>,
}

impl<'a> Kmer<'a> {
    pub fn new(seq: &'a str) -> Kmer<'a> {
        Kmer {
            seq,
            positions: Vec::new(),
        }
    }

    pub fn add_position(&mut self, seq_id: u32, strand: i32, pos: u32) {
        let position = Position::new(seq_id, strand, pos, None, None, None);
        self.positions.push(position);
    }

    pub fn count(&self) -> usize {
        self.positions.len()
    }

    pub fn first_position(&self, half_k: u32) -> bool {
        self.positions.iter().any(|p| p.pos == half_k)
    }
}

impl<'a> fmt::Display for Kmer<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let positions = self.positions.iter()
                                      .map(|p| p.to_string())
                                      .collect::<Vec<String>>()
                                      .join(",");
        write!(f, "{}:{}", self.seq, positions)
    }
}
