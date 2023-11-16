#!/usr/bin/env python3
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
        self.forward_kmers = collections.deque([forward_kmer])  # only used when building
        self.reverse_kmers = collections.deque([reverse_kmer])  # only used when building
        self.forward_seq, self.reverse_seq = '', ''
        self.forward_start_positions, self.forward_end_positions = [], []
        self.reverse_start_positions, self.reverse_end_positions = [], []
        self.depth = 0.0

    def __repr__(self):
        if len(self.forward_seq) < 15:
            seq = self.forward_seq
        else:
            seq = self.forward_seq[:6] + '...' + self.forward_seq[-6:]
        return f'unitig {self.number}: {seq}, {len(self.forward_seq)} bp, {self.depth:.2f}x'

    def length(self):
        assert len(self.forward_seq) == len(self.reverse_seq)
        return len(self.forward_seq)

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
        self.forward_start_positions = self.forward_kmers[0].positions
        self.forward_end_positions = self.forward_kmers[-1].positions
        self.reverse_start_positions = self.reverse_kmers[0].positions
        self.reverse_end_positions = self.reverse_kmers[-1].positions
        for p in self.forward_end_positions:
            p.pos += 1
        for p in self.reverse_end_positions:
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
        self.forward_seq = self.forward_seq[overlap:-overlap]
        self.reverse_seq = self.reverse_seq[overlap:-overlap]
        assert len(self.forward_seq) >= 1

    def gfa_segment_line(self):
        forward_start_positions_str = ','.join(str(p) for p in self.forward_start_positions)
        forward_end_positions_str = ','.join(str(p) for p in self.forward_end_positions)
        reverse_start_positions_str = ','.join(str(p) for p in self.reverse_start_positions)
        reverse_end_positions_str = ','.join(str(p) for p in self.reverse_end_positions)

        return f'S\t{self.number}\t{self.forward_seq}\t' \
               f'DP:f:{self.depth:.2f}\t' \
               f'FS:z:{forward_start_positions_str}\t' \
               f'FE:z:{forward_end_positions_str}\t' \
               f'RS:z:{reverse_start_positions_str}\t' \
               f'RE:z:{reverse_end_positions_str}\n'


class UnitigGraph(object):
    """
    This class builds a unitig graph from a k-mer graph, where all nonbranching paths are merged
    into unitigs. This simplifies things and saves memory.
    """
    def __init__(self, kmer_graph):
        self.unitigs = []
        self.links = []
        self.k_size = kmer_graph.k_size
        print('\nBuilding unitig graph from k-mer graph:')
        self.build_unitigs_from_kmer_graph(kmer_graph)
        self.renumber_unitigs()
        self.create_links()
        self.trim_overlaps()
        self.connect_positions()

    def save_gfa(self, gfa_filename, id_to_contig_info):
        with open(gfa_filename, 'wt') as f:
            f.write('H\tVN:Z:1.0\n')
            for id, contig in id_to_contig_info.items():
                assembly, contig_header, seq_len = contig
                f.write(f'# {id}: {seq_len} bp, {assembly}, {contig_header}\n')
            for unitig in self.unitigs:
                f.write(unitig.gfa_segment_line())
            for a_num, a_strand, b_num, b_strand in self.links:
                f.write(f'L\t{a_num}\t{a_strand}\t{b_num}\t{b_strand}\t0M\n')

    def build_unitigs_from_kmer_graph(self, kmer_graph):
        seen = set()
        unitig_number = 0
        for forward_kmer in kmer_graph.iterate_kmers():
            if forward_kmer in seen:
                continue
            reverse_kmer = kmer_graph.reverse(forward_kmer)

            # Initialise unitig with k-mer
            unitig_number += 1
            unitig = Unitig(unitig_number, forward_kmer, reverse_kmer)
            seen.add(forward_kmer)
            seen.add(reverse_kmer)
            starting_kmer = forward_kmer

            # Extend unitig forward
            while True:
                next_kmers = kmer_graph.next_kmers(forward_kmer)
                if len(next_kmers) != 1:
                    break
                forward_kmer = next_kmers[0]
                if forward_kmer in seen:
                    break
                prev_kmers = kmer_graph.prev_kmers(forward_kmer)
                if len(prev_kmers) != 1:
                    break
                reverse_kmer = kmer_graph.reverse(forward_kmer)
                unitig.add_kmer_to_end(forward_kmer, reverse_kmer)
                seen.add(forward_kmer)
                seen.add(reverse_kmer)

            # Extend unitig backward
            forward_kmer = starting_kmer
            while True:
                prev_kmers = kmer_graph.prev_kmers(forward_kmer)
                if len(prev_kmers) != 1:
                    break
                forward_kmer = prev_kmers[0]
                if forward_kmer in seen:
                    break
                next_kmers = kmer_graph.next_kmers(forward_kmer)
                if len(next_kmers) != 1:
                    break
                reverse_kmer = kmer_graph.reverse(forward_kmer)
                unitig.add_kmer_to_start(forward_kmer, reverse_kmer)
                seen.add(forward_kmer)
                seen.add(reverse_kmer)

            unitig.simplify_seqs()
            self.unitigs.append(unitig)
        print(f'  {len(self.unitigs)} unitigs')

    def renumber_unitigs(self):
        self.unitigs = sorted(self.unitigs, key=lambda u: u.length(), reverse=True)
        unitig_number = 0
        for unitig in self.unitigs:
            unitig_number += 1
            unitig.number = unitig_number

    def create_links(self):
        self.links = []
        piece_len = self.k_size - 1

        # Index unitigs by their k-1 starting and ending sequences
        starting_forward = collections.defaultdict(list)
        starting_reverse = collections.defaultdict(list)
        for unitig in self.unitigs:
            starting_forward[unitig.forward_seq[:piece_len]].append(unitig.number)
            starting_reverse[unitig.reverse_seq[:piece_len]].append(unitig.number)

        for unitig in self.unitigs:
            ending_forward_seq = unitig.forward_seq[-piece_len:]
            for next_unitig in starting_forward[ending_forward_seq]:
                self.links.append((unitig.number, '+', next_unitig, '+'))
                self.links.append((next_unitig, '-', unitig.number, '-'))
            for next_unitig in starting_reverse[ending_forward_seq]:
                self.links.append((unitig.number, '+', next_unitig, '-'))
            ending_reverse_seq = unitig.reverse_seq[-piece_len:]
            for next_unitig in starting_forward[ending_reverse_seq]:
                self.links.append((unitig.number, '-', next_unitig, '+'))

        self.links = sorted(self.links, key=lambda x: (sorted([x[0], x[2]]), x[0]))
        print(f'  {len(self.links)} links')

    def trim_overlaps(self):
        for unitig in self.unitigs:
            unitig.trim_overlaps(self.k_size)
    
    def connect_positions(self):
        pass
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
