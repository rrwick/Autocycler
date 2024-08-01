// This file defines structs which store the sequence, strand and position for contigs in the input
// assemblies.

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


#[derive(Clone, Eq, PartialEq, Hash)]
pub struct Position {
    pub pos: u32,
    seq_id_and_strand: u16,  // seq_id (15 bits) and strand (1 bit) are packed into a u16
}

impl Position {
    const STRAND_BIT_MASK: u16 = 0b1000_0000_0000_0000; // highest bit stores the strand

    pub fn new(seq_id: u16, strand: bool, pos: usize) -> Position {
        let mut seq_id_and_strand = seq_id;
        if strand {
            seq_id_and_strand |= Position::STRAND_BIT_MASK;  // set the strand bit
        }
        Position {
            pos: pos as u32,
            seq_id_and_strand,
        }
    }

    pub fn seq_id(&self) -> u16 {
        self.seq_id_and_strand & !Position::STRAND_BIT_MASK  // mask out the strand bit
    }

    pub fn strand(&self) -> bool {
        // true for forward strand, false for reverse strand
        (self.seq_id_and_strand & Position::STRAND_BIT_MASK) != 0
    }
}

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.seq_id(), if self.strand() { "+" } else { "-" }, self.pos)
    }
}

impl fmt::Debug for Position {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { fmt::Display::fmt(self, f) }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::misc::strand;

    #[test]
    fn test_position() {
        let p1 = Position::new(1, strand::FORWARD, 123);
        let p2 = Position::new(2, strand::REVERSE, 456);
        let p3 = Position::new(32767, strand::FORWARD, 4294967295);  // max values for sed_id and pos
        assert_eq!(format!("{}", p1), "1+123");
        assert_eq!(format!("{}", p2), "2-456");
        assert_eq!(format!("{}", p3), "32767+4294967295");
    }
}
