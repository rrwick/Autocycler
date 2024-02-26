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
use crate::misc::{reverse_complement_u8, quit_with_error};
use crate::position::Position;


pub struct Unitig {
    pub number: u32,
    forward_kmers: VecDeque<*const Kmer>,
    reverse_kmers: VecDeque<*const Kmer>,
    pub forward_seq: Vec<u8>,
    pub reverse_seq: Vec<u8>,
    pub depth: f64,
    pub forward_positions: Vec<Position>,
    pub reverse_positions: Vec<Position>,
    pub forward_next: Vec<(Rc<RefCell<Unitig>>, bool)>,
    pub forward_prev: Vec<(Rc<RefCell<Unitig>>, bool)>,
    pub reverse_next: Vec<(Rc<RefCell<Unitig>>, bool)>,
    pub reverse_prev: Vec<(Rc<RefCell<Unitig>>, bool)>,
    trimmed: bool,
}

impl Unitig {
    pub fn from_kmers(number: u32, forward_kmer: &Kmer, reverse_kmer: &Kmer) -> Self {
        // This constructor is for Unitig objects built from k-mers. Happens in multiple stages:
        // 1. Initialised with a starting k-mer (forward and reverse).
        // 2. K-mers are then added with add_kmer_to_end and add_kmer_to_start methods.
        // 3. The simplify_seqs method combines the k-mers into forward and reverse sequences.
        // 4. The trim_overlaps method removes overlapping sequences from both ends.
        Unitig {
            number: number,
            forward_kmers: VecDeque::from(vec![forward_kmer as *const Kmer]),
            reverse_kmers: VecDeque::from(vec![reverse_kmer as *const Kmer]),
            forward_seq: Vec::new(),
            reverse_seq: Vec::new(),
            depth: 0.0,
            forward_positions: Vec::new(),
            reverse_positions: Vec::new(),
            forward_next: Vec::new(),
            forward_prev: Vec::new(),
            reverse_next: Vec::new(),
            reverse_prev: Vec::new(),
            trimmed: false,
        }
    }

    pub fn from_segment_line(segment_line: &str) -> Self {
        let parts: Vec<&str> = segment_line.split('\t').collect();
        if parts.len() < 3 {
            quit_with_error("Segment line does not have enough parts.");
        }
        let number = parts[1].parse::<u32>().unwrap_or_else(|_| {
            quit_with_error("Unable to parse unitig number.");
            std::process::exit(1)
        });
        let forward_seq = parts[2].as_bytes().to_owned();
        let reverse_seq = reverse_complement_u8(&forward_seq);
        let depth = parts.iter()
            .find(|&p| p.starts_with("DP:f:"))
            .and_then(|p| p[5..].parse::<f64>().ok())
            .unwrap_or_else(|| {
                quit_with_error("Could not find a depth tag (e.g. DP:f:10.00) in the GFA segment \
                                 line.\nAre you sure this is an Autocycler-generated GFA file?");
                std::process::exit(1)
            });
        Unitig {
            number: number,
            forward_kmers: VecDeque::new(),
            reverse_kmers: VecDeque::new(),
            forward_seq,
            reverse_seq,
            depth,
            forward_positions: Vec::new(),
            reverse_positions: Vec::new(),
            forward_next: Vec::new(),
            forward_prev: Vec::new(),
            reverse_next: Vec::new(),
            reverse_prev: Vec::new(),
            trimmed: true,
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
        let trim_start = !self.forward_prev.is_empty();
        let trim_end = !self.forward_next.is_empty();
        if trim_start {
            self.forward_seq = self.forward_seq[overlap..].to_vec();
            self.reverse_seq = self.reverse_seq[..self.reverse_seq.len() - overlap].to_vec();
        }
        if trim_end {
            self.forward_seq = self.forward_seq[..self.forward_seq.len() - overlap].to_vec();
            self.reverse_seq = self.reverse_seq[overlap..].to_vec();
        }
        assert!(self.forward_seq.len() >= 1);
        self.trimmed = true;
    }

    pub fn gfa_segment_line(&self) -> String {
        let seq_str = String::from_utf8_lossy(&self.forward_seq);
        format!("S\t{}\t{}\tDP:f:{:.2}", self.number, seq_str, self.depth)
    }

    pub fn dead_end_start(&self, strand: bool) -> bool {
        match strand {
            true => self.forward_prev.is_empty(),
            false => self.reverse_prev.is_empty(),
        }
    }

    pub fn dead_end_end(&self, strand: bool) -> bool {
        match strand {
            true => self.forward_next.is_empty(),
            false => self.reverse_next.is_empty(),
        }
    }

    pub fn length(&self) -> u32 {
        self.forward_seq.len() as u32
    }

    pub fn untrimmed_length(&self, k_size: usize) -> u32 {
        let half_k = k_size / 2;
        let mut untrimmed_length = self.forward_seq.len();
        if !self.dead_end_start(true) {
            untrimmed_length += half_k;
        }
        if !self.dead_end_end(true) {
            untrimmed_length += half_k;
        }
        untrimmed_length as u32
    }

    pub fn get_seq(&self, strand: bool, upstream: usize, downstream: usize) -> Vec<u8> {
        // This function returns the unitig's sequence on the given strand. It can also add on a
        // bit of upstream or downstream sequence, if available. Note that this only works up to
        // the overlap size, because this is the amount of upstream/downstream sequence that can be
        // reliably found, regardless of path.
        let mut seq = Vec::new();
        if upstream > 0 {
            let upstream_seq = self.get_upstream_seq(strand, upstream);
            seq.extend(&upstream_seq);
        }
        match strand {
            true => seq.extend(&self.forward_seq),
            false => seq.extend(&self.reverse_seq),
        };
        if downstream > 0 {
            let downstream_seq = self.get_downstream_seq(strand, downstream);
            seq.extend(&downstream_seq);
        }
        seq
    }

    fn get_upstream_seq(&self, strand: bool, amount: usize) -> Vec<u8> {
        let mut upstream_seq = Vec::new();
        let mut current_unitig = self as *const Unitig;
        let mut current_strand = strand;
        unsafe {
            loop {
                let prev_unitigs = match current_strand {
                    true => &(*current_unitig).forward_prev,
                    false => &(*current_unitig).reverse_prev,
                };
                if prev_unitigs.len() == 0 { break; }
                let (prev_unitig, prev_strand) = prev_unitigs.first().unwrap();
                let mut seq = prev_unitig.borrow().get_seq(*prev_strand, 0, 0);
                seq.extend(upstream_seq);
                upstream_seq = seq;
                if upstream_seq.len() >= amount { break; }
                current_unitig = &*prev_unitig.borrow() as *const Unitig;
                current_strand = *prev_strand;
            }
        }
        let start = upstream_seq.len().saturating_sub(amount);
        upstream_seq[start..].to_vec()
    }

    fn get_downstream_seq(&self, strand: bool, amount: usize) -> Vec<u8> {
        let mut downstream_seq = Vec::new();
        let mut current_unitig = self as *const Unitig;
        let mut current_strand = strand;
        unsafe {
            loop {
                let next_unitigs = match current_strand {
                    true => &(*current_unitig).forward_next,
                    false => &(*current_unitig).reverse_next,
                };
                if next_unitigs.len() == 0 { break; }
                let (next_unitig, next_strand) = next_unitigs.first().unwrap();
                let seq = next_unitig.borrow().get_seq(*next_strand, 0, 0);
                downstream_seq.extend(seq);
                if downstream_seq.len() >= amount { break; }
                current_unitig = &*next_unitig.borrow() as *const Unitig;
                current_strand = *next_strand;
            }
        }
        let end = if amount > downstream_seq.len() {
            downstream_seq.len()
        } else {
            amount
        };
        downstream_seq[..end].to_vec()
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::Sequence;

    #[test]
    fn test_from_segment_line() {
        let u1 = Unitig::from_segment_line("S\t123\tACGATCGACTACGT\tDP:f:4.56");
        assert_eq!(format!("{}", u1), "unitig 123: ACGATCGACTACGT, 14 bp, 4.56x");

        let u1 = Unitig::from_segment_line("S\t321\tATCGACTACGACTACGACATCG\tDP:f:6.54");
        assert_eq!(format!("{}", u1), "unitig 321: ATCGAC...ACATCG, 22 bp, 6.54x");
    }

    #[test]
    fn test_from_kmers() {
        let seq = Sequence::new(1, "ACGCATAGCACTAGCTACGA".to_string(),
                                "assembly.fasta".to_string(), "contig_1".to_string(), 20);
        let forward_raw = seq.forward_seq.as_ptr();
        let reverse_raw = seq.reverse_seq.as_ptr();

        let forward_k1 = Kmer::new(unsafe{forward_raw.add(4)}, 5, 1);
        let reverse_k1 = Kmer::new(unsafe{reverse_raw.add(11)}, 5, 1);

        let forward_k2 = Kmer::new(unsafe{forward_raw.add(5)}, 5, 1);
        let reverse_k2 = Kmer::new(unsafe{reverse_raw.add(10)}, 5, 1);

        let forward_k3 = Kmer::new(unsafe{forward_raw.add(6)}, 5, 1);
        let reverse_k3 = Kmer::new(unsafe{reverse_raw.add(9)}, 5, 1);

        let mut u = Unitig::from_kmers(123, &forward_k2, &reverse_k2);
        u.add_kmer_to_start(&forward_k1, &reverse_k1);
        u.add_kmer_to_end(&forward_k3, &reverse_k3);
        u.simplify_seqs();

        assert_eq!(u.length(), 7 as u32);
        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "ATAGCAC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GTGCTAT");

        u.trim_overlaps(5);  // trimming does nothing because the unitig has dead-ends
        assert_eq!(u.length(), 7 as u32);
        assert_eq!(std::str::from_utf8(&u.forward_seq).unwrap(), "ATAGCAC");
        assert_eq!(std::str::from_utf8(&u.reverse_seq).unwrap(), "GTGCTAT");

        assert_eq!(u.untrimmed_length(5), 7 as u32);
    }
}
