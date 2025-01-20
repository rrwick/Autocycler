// This file defines the Unitig struct.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::cell::RefCell;
use std::collections::VecDeque;
use std::fmt;
use std::rc::Rc;

use crate::kmer_graph::Kmer;
use crate::misc::{quit_with_error, reverse_complement, strand};
use crate::position::Position;


static ANCHOR_COLOUR: &str = "forestgreen";
static BRIDGE_COLOUR: &str = "pink";
static CONSENTIG_COLOUR: &str = "steelblue";
static OTHER_COLOUR: &str = "orangered";


#[derive(Clone, Default)]
pub struct Unitig {
    pub number: u32,
    pub forward_kmers: VecDeque<*const Kmer>,
    pub reverse_kmers: VecDeque<*const Kmer>,
    pub forward_seq: Vec<u8>,
    pub reverse_seq: Vec<u8>,
    pub depth: f64,
    pub unitig_type: UnitigType,  // anchor, bridge, consentig or other
    pub forward_positions: Vec<Position>,
    pub reverse_positions: Vec<Position>,
    pub forward_next: Vec<UnitigStrand>,
    pub forward_prev: Vec<UnitigStrand>,
    pub reverse_next: Vec<UnitigStrand>,
    pub reverse_prev: Vec<UnitigStrand>,
}

impl Unitig {
    pub fn from_kmers(number: u32, forward_kmer: &Kmer, reverse_kmer: &Kmer) -> Self {
        // This constructor is for Unitig objects built from k-mers. Happens in multiple stages:
        // 1. Initialised with a starting k-mer (forward and reverse).
        // 2. K-mers are then added with add_kmer_to_end and add_kmer_to_start methods.
        // 3. The simplify_seqs method combines the k-mers into forward and reverse sequences.
        // 4. The trim_overlaps method removes overlapping sequences from both ends.
        Unitig {
            number,
            forward_kmers: VecDeque::from(vec![forward_kmer as *const Kmer]),
            reverse_kmers: VecDeque::from(vec![reverse_kmer as *const Kmer]),
            ..Default::default()
        }
    }

    pub fn from_segment_line(segment_line: &str) -> Self {
        let parts: Vec<&str> = segment_line.split('\t').collect();
        if parts.len() < 3 {
            quit_with_error("Segment line does not have enough parts.");
        }
        let number = parts[1].parse::<u32>().unwrap_or_else(|_| {
            quit_with_error("Unable to parse unitig number.");
        });
        let forward_seq = parts[2].as_bytes().to_owned();
        let reverse_seq = reverse_complement(&forward_seq);
        let depth = parts.iter()
            .find(|&p| p.starts_with("DP:f:")).and_then(|p| p[5..].parse::<f64>().ok())
            .unwrap_or_else(|| {
                quit_with_error("Could not find a depth tag (e.g. DP:f:10.00) in the GFA segment \
                                 line.\nAre you sure this is an Autocycler-generated GFA file?");
            });
        let unitig_type = if parts.iter().any(|p| *p == format!("CL:z:{}", CONSENTIG_COLOUR)) {
            UnitigType::Consentig
        } else if parts.iter().any(|p| *p == format!("CL:z:{}", ANCHOR_COLOUR)) {
            UnitigType::Anchor
        } else if parts.iter().any(|p| *p == format!("CL:z:{}", BRIDGE_COLOUR)) {
            UnitigType::Bridge
        } else {
            UnitigType::Other
        };
        Unitig {
            number, forward_seq, reverse_seq, depth, unitig_type,
            ..Default::default()
        }
    }

    pub fn bridge(number: u32, forward_seq: Vec<u8>, depth: f64) -> Self {
        // This constructor is for manually building a Unitig object when creating bridges.
        let reverse_seq = reverse_complement(&forward_seq);
        Unitig {
            number, forward_seq, reverse_seq, depth, unitig_type: UnitigType::Bridge,
            ..Default::default()
        }
    }

    pub fn add_kmer_to_end(&mut self, forward_kmer: &Kmer, reverse_kmer: &Kmer) {
        self.forward_kmers.push_back(forward_kmer);
        self.reverse_kmers.push_front(reverse_kmer);
    }

    pub fn add_kmer_to_start(&mut self, forward_kmer: &Kmer, reverse_kmer: &Kmer) {
        self.forward_kmers.push_front(forward_kmer);
        self.reverse_kmers.push_back(reverse_kmer);
    }

    pub fn simplify_seqs(&mut self) {
        self.combine_kmers_into_sequences();
        self.set_positions();
        self.set_average_depth();
        self.forward_kmers.clear();
        self.reverse_kmers.clear();
    }

    fn combine_kmers_into_sequences(&mut self) {
        if let Some(first_kmer) = self.forward_kmers.front() {
            self.forward_seq = unsafe{&**first_kmer}.seq().to_vec();
            self.forward_kmers.iter().skip(1).for_each(|kmer| {
                self.forward_seq.push(*unsafe{&**kmer}.seq().last().unwrap());
            });
        }
        if let Some(first_kmer) = self.reverse_kmers.front() {
            self.reverse_seq = unsafe{&**first_kmer}.seq().to_vec();
            self.reverse_kmers.iter().skip(1).for_each(|kmer| {
                self.reverse_seq.push(*unsafe{&**kmer}.seq().last().unwrap());
            });
        }
    }

    fn set_positions(&mut self) {
        // Sets the Unitig's positions on each strand to be the same as the positions of the first
        // Kmer on that strand.
        if let Some(&kmer_ptr) = self.forward_kmers.front() {
            let kmer = unsafe { &*kmer_ptr };
            self.forward_positions.extend_from_slice(&kmer.positions);
        }
        if let Some(&kmer_ptr) = self.reverse_kmers.front() {
            let kmer = unsafe { &*kmer_ptr };
            self.reverse_positions.extend_from_slice(&kmer.positions);
        }
    }

    fn set_average_depth(&mut self) {
        let for_depths: Vec<f64> = self.forward_kmers.iter().map(|k| unsafe{&**k}.depth() as f64).collect();
        let rev_depths: Vec<f64> = self.reverse_kmers.iter().map(|k| unsafe{&**k}.depth() as f64).collect();
        let forward_avg = for_depths.iter().sum::<f64>() / for_depths.len() as f64;
        let reverse_avg = rev_depths.iter().sum::<f64>() / rev_depths.len() as f64;
        assert_eq!(forward_avg, reverse_avg);
        self.depth = forward_avg;
    }

    pub fn trim_overlaps(&mut self, k_size: usize) {
        let overlap = k_size / 2;
        assert!(self.forward_seq.len() >= k_size);
        self.forward_seq = self.forward_seq[overlap..].to_vec();
        self.reverse_seq = self.reverse_seq[..self.reverse_seq.len() - overlap].to_vec();
        self.forward_seq = self.forward_seq[..self.forward_seq.len() - overlap].to_vec();
        self.reverse_seq = self.reverse_seq[overlap..].to_vec();
        assert!(!self.forward_seq.is_empty());
    }

    pub fn gfa_segment_line(&self, use_other_colour: bool) -> String {
        let seq_str = String::from_utf8_lossy(&self.forward_seq);
        format!("S\t{}\t{}\tDP:f:{:.2}{}", self.number, seq_str, self.depth,
                                           self.colour_tag(use_other_colour))
    }

    pub fn colour_tag(&self, use_other_colour: bool) -> String {
        match self.unitig_type {
            UnitigType::Consentig => format!("\tCL:z:{}", CONSENTIG_COLOUR),
            UnitigType::Anchor => format!("\tCL:z:{}", ANCHOR_COLOUR),
            UnitigType::Bridge => format!("\tCL:z:{}", BRIDGE_COLOUR),
            UnitigType::Other => { if use_other_colour { format!("\tCL:z:{}", OTHER_COLOUR) }
                                                  else { String::new() } }
        }
    }

    pub fn length(&self) -> u32 {
        self.forward_seq.len() as u32
    }

    pub fn get_seq(&self, strand: bool) -> Vec<u8> {
        // This function returns the unitig's sequence on the given strand.
        if strand {
            self.forward_seq.clone()
        } else {
            self.reverse_seq.clone()
        }
    }

    pub fn blunt_start(&self) -> bool {
        self.reverse_next.is_empty()
    }

    pub fn blunt_end(&self) -> bool {
        self.forward_next.is_empty()
    }

    pub fn hairpin_start(&self) -> bool {
        self.reverse_next.len() == 1 &&
        self.reverse_next[0].strand == strand::FORWARD &&
        self.reverse_next[0].unitig.borrow().number == self.number
    }

    pub fn hairpin_end(&self) -> bool {
        self.forward_next.len() == 1 &&
        self.forward_next[0].strand == strand::REVERSE &&
        self.forward_next[0].unitig.borrow().number == self.number
    }

    pub fn remove_seq_from_start(&mut self, amount: usize) {
        for p in &mut self.forward_positions {
            p.pos += amount as u32
        }
        assert!(amount <= self.forward_seq.len());
        self.forward_seq.drain(0..amount);
        self.reverse_seq.truncate(self.reverse_seq.len() - amount);
    }

    pub fn remove_seq_from_end(&mut self, amount: usize) {
        for p in &mut self.reverse_positions {
            p.pos += amount as u32
        }
        assert!(amount <= self.forward_seq.len());
        self.forward_seq.truncate(self.reverse_seq.len() - amount);
        self.reverse_seq.drain(0..amount);
    }

    pub fn add_seq_to_start(&mut self, seq: Vec<u8>) {
        for p in &mut self.forward_positions {
            p.pos -= seq.len() as u32;
        }
        self.forward_seq.splice(0..0, seq.iter().cloned());
        self.reverse_seq = reverse_complement(&self.forward_seq);
    }

    pub fn add_seq_to_end(&mut self, mut seq: Vec<u8>) {
        for p in &mut self.reverse_positions {
            p.pos -= seq.len() as u32;
        }
        self.forward_seq.append(&mut seq);
        self.reverse_seq = reverse_complement(&self.forward_seq);
    }

    pub fn remove_sequence(&mut self, id: u16) {
        // Removes all Positions from the Unitig which have the given sequence ID. This can reduce
        // the Unitig's depth.
        self.forward_positions.retain(|p| p.seq_id() != id);
        self.reverse_positions.retain(|p| p.seq_id() != id);
        assert_eq!(self.forward_positions.len(), self.reverse_positions.len());
        self.recalculate_depth();
    }

    pub fn recalculate_depth(&mut self) {
        self.depth = self.forward_positions.len() as f64;
    }

    pub fn clear_positions(&mut self) {
        self.forward_positions.clear();
        self.reverse_positions.clear();
    }

    pub fn reduce_depth_by_one(&mut self) {
        self.depth -= 1.0;
        if self.depth < 0.0 {
            self.depth = 0.0;
        }
    }

    pub fn is_isolated_and_circular(&self) -> bool {
        // Returns whether or not this unitig has a circularising link and no other links.
        if self.forward_next.len() != 1 || self.forward_prev.len() != 1 {
            return false;
        }
        let next = &self.forward_next[0];
        let prev = &self.forward_prev[0];
        next.number() == self.number && next.strand && prev.number() == self.number && prev.strand
    }

    pub fn clear_all_links(&mut self) {
        self.forward_next.clear();
        self.forward_prev.clear();
        self.reverse_next.clear();
        self.reverse_prev.clear();
    }
}

impl fmt::Display for Unitig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let display_seq = if self.forward_seq.len() < 15 {
            String::from_utf8_lossy(&self.forward_seq).to_string()
        } else {
            format!("{}...{}",
                    String::from_utf8_lossy(&self.forward_seq[..6]),
                    String::from_utf8_lossy(&self.forward_seq[self.forward_seq.len() - 6..]))
        };
        write!(f, "unitig {}: {}, {} bp, {:.2}x", self.number, display_seq, self.forward_seq.len(), self.depth)
    }
}

impl fmt::Debug for Unitig {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { fmt::Display::fmt(self, f) }
}


// Since Unitigs are often dealt with in a strand-specific manner, this struct bundles up a Unitig
// (via a reference-counted reference) along with its strand.
#[derive(Clone)]
pub struct UnitigStrand {
    pub unitig: Rc<RefCell<Unitig>>,
    pub strand: bool,
}

impl UnitigStrand {
    pub fn new(unitig: &Rc<RefCell<Unitig>>, strand: bool) -> Self {
        UnitigStrand {
            unitig: Rc::clone(unitig),
            strand,
        }
    }

    pub fn number(&self) -> u32 {
        self.unitig.borrow().number
    }

    pub fn signed_number(&self) -> i32 {
        if self.strand {
            self.unitig.borrow().number as i32
        } else {
            -(self.unitig.borrow().number as i32)
        }
    }

    pub fn length(&self) -> u32 {
        self.unitig.borrow().length()
    }

    pub fn depth(&self) -> f64 {
        self.unitig.borrow().depth
    }

    pub fn get_seq(&self) -> Vec<u8> {
        self.unitig.borrow().get_seq(self.strand)
    }

    pub fn is_anchor(&self) -> bool {
        self.unitig.borrow().unitig_type == UnitigType::Anchor
    }

    pub fn is_consentig(&self) -> bool {
        self.unitig.borrow().unitig_type == UnitigType::Consentig
    }
}

impl fmt::Display for UnitigStrand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.unitig.borrow().number, if self.strand { "+" } else { "-" })
    }
}

impl fmt::Debug for UnitigStrand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { fmt::Display::fmt(self, f) }
}


#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UnitigType {
    Anchor,
    Bridge,
    Consentig,
    #[default]
    Other,
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::Sequence;
    use crate::misc::strand;

    #[test]
    fn test_from_segment_line() {
        let u1 = Unitig::from_segment_line("S\t123\tACGATCGACTACGT\tDP:f:4.56");
        assert_eq!(format!("{}", u1), "unitig 123: ACGATCGACTACGT, 14 bp, 4.56x");

        let u1 = Unitig::from_segment_line("S\t321\tATCGACTACGACTACGACATCG\tDP:f:6.54");
        assert_eq!(format!("{}", u1), "unitig 321: ATCGAC...ACATCG, 22 bp, 6.54x");
    }

    #[test]
    fn test_from_kmers() {
        let k_size = 5; let half_k = k_size / 2;
        let seq = Sequence::new_with_seq(1, "ACGCATAGCACTAGCTACGA".to_string(),
                                         "assembly.fasta".to_string(), "contig_1".to_string(), 20, half_k);
        let forward_raw = seq.forward_seq.as_ptr();
        let reverse_raw = seq.reverse_seq.as_ptr();

        let forward_k1 = Kmer::new(unsafe{forward_raw.add(4)}, 5, 1);
        let reverse_k1 = Kmer::new(unsafe{reverse_raw.add(15)}, 5, 1);

        let forward_k2 = Kmer::new(unsafe{forward_raw.add(5)}, 5, 1);
        let reverse_k2 = Kmer::new(unsafe{reverse_raw.add(14)}, 5, 1);

        let forward_k3 = Kmer::new(unsafe{forward_raw.add(6)}, 5, 1);
        let reverse_k3 = Kmer::new(unsafe{reverse_raw.add(13)}, 5, 1);

        let mut u = Unitig::from_kmers(123, &forward_k2, &reverse_k2);
        u.add_kmer_to_start(&forward_k1, &reverse_k1);
        u.add_kmer_to_end(&forward_k3, &reverse_k3);
        u.simplify_seqs();

        assert_eq!(u.length(), 7_u32);
        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "GCATAGC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GCTATGC");

        u.trim_overlaps(k_size as usize);
        assert_eq!(u.length(), 3_u32);
        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "ATA");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "TAT");
    }

    #[test]
    fn test_get_seq() {
        let unitig_a = Rc::new(RefCell::new(Unitig::from_segment_line("S\t1\tGCTGAAGGGC\tDP:f:1")));
        let unitig_b = Rc::new(RefCell::new(Unitig::from_segment_line("S\t2\tCGCGTTCGAC\tDP:f:1")));
        unitig_a.borrow_mut().forward_next.push(UnitigStrand::new(&unitig_b, strand::FORWARD));
        unitig_b.borrow_mut().forward_prev.push(UnitigStrand::new(&unitig_a, strand::FORWARD));
        unitig_b.borrow_mut().reverse_next.push(UnitigStrand::new(&unitig_a, strand::REVERSE));
        unitig_a.borrow_mut().reverse_prev.push(UnitigStrand::new(&unitig_b, strand::REVERSE));

        assert_eq!(std::str::from_utf8(&unitig_a.borrow().get_seq(strand::FORWARD)).unwrap(), "GCTGAAGGGC");
        assert_eq!(std::str::from_utf8(&unitig_a.borrow().get_seq(strand::REVERSE)).unwrap(), "GCCCTTCAGC");
        assert_eq!(std::str::from_utf8(&unitig_b.borrow().get_seq(strand::FORWARD)).unwrap(), "CGCGTTCGAC");
        assert_eq!(std::str::from_utf8(&unitig_b.borrow().get_seq(strand::REVERSE)).unwrap(), "GTCGAACGCG");
    }

    #[test]
    fn test_remove_seq_from_start() {
        let mut u = Unitig::from_segment_line("S\t1\tGCTGAAGGGC\tDP:f:1");
        u.forward_positions.push(Position::new(1, strand::FORWARD, 100));
        u.reverse_positions.push(Position::new(2, strand::REVERSE, 890));
        u.forward_positions.push(Position::new(2, strand::REVERSE, 200));
        u.reverse_positions.push(Position::new(2, strand::FORWARD, 790));

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "GCTGAAGGGC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GCCCTTCAGC");
        assert_eq!(u.forward_positions[0].pos, 100);
        assert_eq!(u.reverse_positions[0].pos, 890);
        assert_eq!(u.forward_positions[1].pos, 200);
        assert_eq!(u.reverse_positions[1].pos, 790);

        u.remove_seq_from_start(2);

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "TGAAGGGC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GCCCTTCA");
        assert_eq!(u.forward_positions[0].pos, 102);
        assert_eq!(u.reverse_positions[0].pos, 890);
        assert_eq!(u.forward_positions[1].pos, 202);
        assert_eq!(u.reverse_positions[1].pos, 790);
    }

    #[test]
    fn test_remove_seq_from_end() {
        let mut u = Unitig::from_segment_line("S\t1\tGCTGAAGGGC\tDP:f:1");
        u.forward_positions.push(Position::new(1, strand::FORWARD, 100));
        u.reverse_positions.push(Position::new(2, strand::REVERSE, 890));
        u.forward_positions.push(Position::new(2, strand::REVERSE, 200));
        u.reverse_positions.push(Position::new(2, strand::FORWARD, 790));

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "GCTGAAGGGC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GCCCTTCAGC");
        assert_eq!(u.forward_positions[0].pos, 100);
        assert_eq!(u.reverse_positions[0].pos, 890);
        assert_eq!(u.forward_positions[1].pos, 200);
        assert_eq!(u.reverse_positions[1].pos, 790);

        u.remove_seq_from_end(2);

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "GCTGAAGG");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "CCTTCAGC");
        assert_eq!(u.forward_positions[0].pos, 100);
        assert_eq!(u.reverse_positions[0].pos, 892);
        assert_eq!(u.forward_positions[1].pos, 200);
        assert_eq!(u.reverse_positions[1].pos, 792);
    }

    #[test]
    fn test_add_seq_to_start() {
        let mut u = Unitig::from_segment_line("S\t1\tGCTGAAGGGC\tDP:f:1");
        u.forward_positions.push(Position::new(1, strand::FORWARD, 100));
        u.reverse_positions.push(Position::new(2, strand::REVERSE, 890));
        u.forward_positions.push(Position::new(2, strand::REVERSE, 200));
        u.reverse_positions.push(Position::new(2, strand::FORWARD, 790));

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "GCTGAAGGGC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GCCCTTCAGC");
        assert_eq!(u.forward_positions[0].pos, 100);
        assert_eq!(u.reverse_positions[0].pos, 890);
        assert_eq!(u.forward_positions[1].pos, 200);
        assert_eq!(u.reverse_positions[1].pos, 790);

        u.add_seq_to_start("AC".into());

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "ACGCTGAAGGGC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GCCCTTCAGCGT");
        assert_eq!(u.forward_positions[0].pos, 98);
        assert_eq!(u.reverse_positions[0].pos, 890);
        assert_eq!(u.forward_positions[1].pos, 198);
        assert_eq!(u.reverse_positions[1].pos, 790);
    }

    #[test]
    fn test_add_seq_to_end() {
        let mut u = Unitig::from_segment_line("S\t1\tGCTGAAGGGC\tDP:f:1");
        u.forward_positions.push(Position::new(1, strand::FORWARD, 100));
        u.reverse_positions.push(Position::new(2, strand::REVERSE, 890));
        u.forward_positions.push(Position::new(2, strand::REVERSE, 200));
        u.reverse_positions.push(Position::new(2, strand::FORWARD, 790));

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "GCTGAAGGGC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GCCCTTCAGC");
        assert_eq!(u.forward_positions[0].pos, 100);
        assert_eq!(u.reverse_positions[0].pos, 890);
        assert_eq!(u.forward_positions[1].pos, 200);
        assert_eq!(u.reverse_positions[1].pos, 790);

        u.add_seq_to_end("AC".into());

        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "GCTGAAGGGCAC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GTGCCCTTCAGC");
        assert_eq!(u.forward_positions[0].pos, 100);
        assert_eq!(u.reverse_positions[0].pos, 888);
        assert_eq!(u.forward_positions[1].pos, 200);
        assert_eq!(u.reverse_positions[1].pos, 788);
    }
}
