// This file defines structs which store the sequence, strand and position for contigs in the input
// assemblies:
// * KmerPos objects:
//   * Hold just sequence, strand and position.
//   * Are used for every k-mer in the input assemblies (so there a lot of them).
//   * For this reason, I keep them as small as possible to save memory.
// * UnitigPos objects:
//   * In addition to sequence, strand and position, they also store a reference to their Unitig
//     (including the strand and whether they are on the start or end of the Unitig) and links to 
//     the next and previous UnitigPos objects (forming a doubly-linked list).
//   * Are used only at the start/end of each Unitig (so there are fewer of them).
//   * For this reason, memory efficiency matters less.

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
use std::ptr;

use crate::unitig::Unitig;


pub struct KmerPos {
    pub pos: u32,
    seq_id_and_strand: u16, // seq_id (15 bits) and strand (1 bit) are packed into a u16
}

impl KmerPos {
    const STRAND_BIT_MASK: u16 = 0b1000_0000_0000_0000; // highest bit stores the strand

    pub fn new(seq_id: u16, strand: bool, pos: usize) -> KmerPos {
        let mut seq_id_and_strand = seq_id;
        if strand {
            seq_id_and_strand |= KmerPos::STRAND_BIT_MASK; // Set the strand bit
        }
        KmerPos {
            pos: pos as u32,
            seq_id_and_strand,
        }
    }

    pub fn seq_id(&self) -> u16 {
        self.seq_id_and_strand & !KmerPos::STRAND_BIT_MASK // Mask out the strand bit
    }

    pub fn strand(&self) -> bool {
        // true for forward strand, false for reverse strand
        (self.seq_id_and_strand & KmerPos::STRAND_BIT_MASK) != 0
    }
}

impl fmt::Display for KmerPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.seq_id(), if self.strand() { "+" } else { "-" }, self.pos)
    }
}


pub struct UnitigPos {
    pub seq_id: u16,
    pub strand: bool, // true for forward strand, false for reverse strand
    pub pos: u32,
    pub prev: *mut UnitigPos,
    pub next: *mut UnitigPos,
    unitig: *mut Unitig,
    unitig_strand: bool, // true for forward strand, false for reverse strand
    unitig_start_end: bool, // true for start, false for end
}

impl UnitigPos {
    pub fn new(kmer_pos: &KmerPos, unitig: *mut Unitig, unitig_strand: bool,
               unitig_start_end: bool) -> UnitigPos {
        UnitigPos {
            seq_id: kmer_pos.seq_id(),
            strand: kmer_pos.strand(),
            pos: kmer_pos.pos,
            prev: ptr::null_mut(),
            next: ptr::null_mut(),
            unitig,
            unitig_strand,
            unitig_start_end,
        }
    }

    fn on_unitig_start(&self) -> bool {
        self.unitig_start_end
    }

    fn on_unitig_end(&self) -> bool {
        !self.unitig_start_end
    }
}

impl fmt::Display for UnitigPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.seq_id, if self.strand { "+" } else { "-" }, self.pos)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmerpos() {
        let p1 = KmerPos::new(1, true, 123);
        let p2 = KmerPos::new(2, false, 456);
        let p3 = KmerPos::new(32767, true, 4294967295);  // max values for sed_id and pos
        assert_eq!(format!("{}", p1), "1+123");
        assert_eq!(format!("{}", p2), "2-456");
        assert_eq!(format!("{}", p3), "32767+4294967295");
    }

    #[test]
    fn test_unitigpos() {
        let k1 = KmerPos::new(1, true, 123);
        let k2 = KmerPos::new(2, false, 456);
        let k3 = KmerPos::new(32767, true, 4294967295);
        let p1 = UnitigPos::new(&k1, ptr::null_mut(), true, true);
        let p2 = UnitigPos::new(&k2, ptr::null_mut(), true, false);
        let p3 = UnitigPos::new(&k3, ptr::null_mut(), false, true);
        assert_eq!(format!("{}", p1), "1+123");
        assert_eq!(format!("{}", p2), "2-456");
        assert_eq!(format!("{}", p3), "32767+4294967295");
    }
}
