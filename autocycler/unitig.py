"""
Copyright 2023 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Autocycler

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import collections
import copy
import statistics

from .misc import reverse_complement


class Unitig(object):
    """
    Unitig objects are built in multiple stages:
    1. Initialised with a starting k-mer (forward and reverse).
    2. K-mers are added with add_kmer_to_end and add_kmer_to_start methods.
    3. The simplify_seqs method combines the k-mers into forward and reverse sequences.
    4. The trim_overlaps method removes overlapping sequences from both ends.
    """
    def __init__(self, number, forward_kmer, reverse_kmer):
        self.number = number

        # These k-mer deques are used when building the unique (will be set to None afterward):
        self.forward_kmers = collections.deque([forward_kmer])
        self.reverse_kmers = collections.deque([reverse_kmer])

        # After building is complete, only the unitig sequences are stored:
        self.forward_seq, self.reverse_seq = '', ''
        self.depth = 0.0

        # Position objects (connecting back to the original sequences) are stored for the start and
        # end of both strands:
        self.forward_start_positions, self.forward_end_positions = [], []
        self.reverse_start_positions, self.reverse_end_positions = [], []

        # Pointers to the preceding and following Unitig objects, along with their strand:
        self.forward_next, self.forward_prev = [], []
        self.reverse_next, self.reverse_prev = [], []

    def __repr__(self):
        if len(self.forward_seq) < 15:
            seq = self.forward_seq
        else:
            seq = self.forward_seq[:6] + '...' + self.forward_seq[-6:]
        return f'unitig {self.number}: {seq}, {len(self.forward_seq)} bp, {self.depth:.2f}x'

    def length(self):
        """
        This method will only work after the unitig has been built (simplify_seqs has been run).
        """
        return len(self.forward_seq)

    def get_seq(self, strand, upstream=0, downstream=0):
        """
        This function returns the unitig's sequence on the given strand. It can also add on a bit
        of upstream or downstream sequence, if available. Note that this only works up to the
        overlap size, because this is the amount of upstream/downstream sequence that can be
        reliably found, regardless of path.
        """
        if strand == 1:
            seq = self.forward_seq
        elif strand == -1:
            seq = self.reverse_seq
        else:
            assert False
        if upstream:
            seq = self.get_upstream_seq(strand, upstream) + seq
        if downstream:
            seq += self.get_downstream_seq(strand, downstream)
        return seq

    def get_upstream_seq(self, strand, amount):
        upstream_seq = ''
        unitig = self
        while True:
            prev_unitigs = unitig.forward_prev if strand == 1 else unitig.reverse_prev
            if not prev_unitigs:
                break
            unitig, strand = prev_unitigs[0]
            upstream_seq = unitig.get_seq(strand) + upstream_seq
            if len(upstream_seq) >= amount:
                break
        return upstream_seq[-amount:]

    def get_downstream_seq(self, strand, amount):
        downstream_seq = ''
        unitig = self
        while True:
            next_unitigs = unitig.forward_next if strand == 1 else unitig.reverse_next
            if not next_unitigs:
                break
            unitig, strand = next_unitigs[0]
            downstream_seq += unitig.get_seq(strand)
            if len(downstream_seq) >= amount:
                break
        return downstream_seq[:amount]

    def add_kmer_to_end(self, forward_kmer, reverse_kmer):
        self.forward_kmers.append(forward_kmer)
        self.reverse_kmers.appendleft(reverse_kmer)

    def add_kmer_to_start(self, forward_kmer, reverse_kmer):
        self.forward_kmers.appendleft(forward_kmer)
        self.reverse_kmers.append(reverse_kmer)

    def simplify_seqs(self):
        self.combine_kmers_into_sequences()
        self.set_start_end_positions()
        self.set_average_depth()
        self.forward_kmers, self.reverse_kmers = None, None

    def combine_kmers_into_sequences(self):
        forward_seq, reverse_seq = [], []
        for k in self.forward_kmers:
            if not forward_seq:
                forward_seq.append(k.seq)
            else:
                forward_seq.append(k.seq[-1])
        for k in self.reverse_kmers:
            if not reverse_seq:
                reverse_seq.append(k.seq)
            else:
                reverse_seq.append(k.seq[-1])
        self.forward_seq = ''.join(forward_seq)
        self.reverse_seq = ''.join(reverse_seq)
        assert reverse_complement(self.forward_seq) == self.reverse_seq

    def set_start_end_positions(self):
        """
        Sets this unitig's start and end for both strands. The end positions get a 1 added to make
        them exclusive-end Pythonic ranges.
        """
        self.forward_start_positions = copy.deepcopy(self.forward_kmers[0].positions)
        self.forward_end_positions = copy.deepcopy(self.forward_kmers[-1].positions)
        self.reverse_start_positions = copy.deepcopy(self.reverse_kmers[0].positions)
        self.reverse_end_positions = copy.deepcopy(self.reverse_kmers[-1].positions)
        for p in self.forward_start_positions:
            p.unitig = self
            p.unitig_strand = 1
            p.unitig_start_end = 0
        for p in self.forward_end_positions:
            p.unitig = self
            p.unitig_strand = 1
            p.unitig_start_end = 1
            p.pos += 1
        for p in self.reverse_start_positions:
            p.unitig = self
            p.unitig_strand = -1
            p.unitig_start_end = 0
        for p in self.reverse_end_positions:
            p.unitig = self
            p.unitig_strand = -1
            p.unitig_start_end = 1
            p.pos += 1

    def set_average_depth(self):
        forward_depths, reverse_depths = [], []
        for k in self.forward_kmers:
            forward_depths.append(k.count())
        for k in self.reverse_kmers:
            reverse_depths.append(k.count())
        forward_depth = statistics.mean(forward_depths)
        reverse_depth = statistics.mean(reverse_depths)
        assert forward_depth == reverse_depth
        self.depth = forward_depth

    def trim_overlaps(self, k_size):
        overlap = k_size // 2
        assert len(self.forward_seq) >= k_size
        trim_start = len(self.forward_prev) > 0
        trim_end = len(self.forward_next) > 0
        if trim_start:
            self.forward_seq = self.forward_seq[overlap:]
            self.reverse_seq = self.reverse_seq[:-overlap]
        if trim_end:
            self.forward_seq = self.forward_seq[:-overlap]
            self.reverse_seq = self.reverse_seq[overlap:]
        assert reverse_complement(self.forward_seq) == self.reverse_seq
        assert len(self.forward_seq) >= 1

    def gfa_segment_line(self):
        return f'S\t{self.number}\t{self.forward_seq}\t' \
               f'DP:f:{self.depth:.2f}\n'

    def connect_positions(self, k_size):
        """
        Connects the start and end positions for this unitig on both strands, when they line up.
        There should always be a 1-to-1 relationship between starting and ending positions, i.e.
        every starting position will connect to an ending position.
        """
        adjusted_length = self.length() - k_size + 1
        for start in self.forward_start_positions:
            assert start.prev is None and start.next is None
            matches = [end for end in self.forward_end_positions
                       if start.seq_id == end.seq_id and start.strand == end.strand and \
                       start.pos + adjusted_length == end.pos]
            assert len(matches) == 1
            end = matches[0]
            start.next = end
            end.prev = start
        for start in self.reverse_start_positions:
            assert start.prev is None and start.next is None
            matches = [end for end in self.reverse_end_positions
                       if start.seq_id == end.seq_id and start.strand == end.strand and \
                       start.pos + adjusted_length == end.pos]
            assert len(matches) == 1
            end = matches[0]
            start.next = end
            end.prev = start
