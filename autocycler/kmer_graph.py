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

from .misc import reverse_complement, reverse_complement_position, iterate_fasta


class Position(object):
    """
    Position objects store the sequence, strand and position for contigs in the input assemblies.
    They are used:
    * In KmerGraph objects, where each Kmer object has a single Position.
    * In UnitigGraph objects, where each Unitig has one or more starting and ending Position
      objects on both strands. These form a doubly linked list, tracing the input contig through
      the UnitigGraph.
    """
    def __init__(self, seq_id, strand, pos):
        self.seq_id = seq_id
        self.strand = strand  # 1 for forward strand, -1 for reverse strand
        self.pos = pos  # 0-based indexing

        # Pointers to the preceding and following Position objects:
        self.prev_kmer_position = None
        self.next_kmer_position = None

    def __repr__(self):
        return f'{self.seq_id}{"+" if self.strand == 1 else "-"}{self.pos}'


class Kmer(object):
    def __init__(self, seq):
        self.seq = seq
        self.positions = []

    def __repr__(self):
        positions = ','.join([str(p) for p in sorted(self.positions)])
        return f'{positions}'

    def add_position(self, seq_id, strand, pos):
        self.positions.append(Position(seq_id, strand, pos))

    def count(self):
        return len(self.positions)


class KmerGraph(object):
    """
    This class stores all k-mers in the sequences that were added to it (a De Bruijn graph), and
    for each k-mer it stores the positions in the input sequences where that k-mer occurred.
    """
    def __init__(self, k_size):
        self.k_size = k_size
        self.kmers = {}

    def add_sequence(self, seq, seq_id):
        half_k = self.k_size // 2  # actually k_size/2 - 0.5, because k_size is odd
        for forward_pos in range(len(seq) - self.k_size + 1):
            reverse_pos = len(seq) - forward_pos - self.k_size

            forward_seq = seq[forward_pos:forward_pos+self.k_size]
            reverse_seq = reverse_complement(forward_seq)

            if forward_seq not in self.kmers:
                self.kmers[forward_seq] = Kmer(forward_seq)
                self.kmers[reverse_seq] = Kmer(reverse_seq)

            forward_centre_pos = forward_pos + half_k
            reverse_centre_pos = reverse_pos + half_k
            assert reverse_centre_pos == reverse_complement_position(forward_centre_pos, len(seq))

            self.kmers[forward_seq].add_position(seq_id, 1, forward_centre_pos)
            self.kmers[reverse_seq].add_position(seq_id, -1, reverse_centre_pos)

    def add_assemblies(self, assemblies):
        id_to_contig_info = {}
        seq_id = 0
        for assembly in assemblies:
            print(f'\nAdding {assembly} to graph:')
            for name, info, seq in iterate_fasta(assembly):
                seq_id += 1
                print(f'  {seq_id}: {name} ({len(seq)} bp)...', flush=True, end='')
                self.add_sequence(seq, seq_id)
                contig_header = name + ' ' + info if info else name
                id_to_contig_info[seq_id] = (assembly, contig_header, len(seq))
                print(' done')
        print(f'\nGraph contains {len(self.kmers)} k-mers')
        return id_to_contig_info

    def next_kmers(self, kmer):
        """
        Returns a list of the k-mers in the graph which follow the given one. The list will have a
        length from 0 to 4.
        """
        a = kmer.seq[1:] + 'A'
        c = kmer.seq[1:] + 'C'
        g = kmer.seq[1:] + 'G'
        t = kmer.seq[1:] + 'T'
        next_list = []
        if a in self.kmers:
            next_list.append(self.kmers[a])
        if c in self.kmers:
            next_list.append(self.kmers[c])
        if g in self.kmers:
            next_list.append(self.kmers[g])
        if t in self.kmers:
            next_list.append(self.kmers[t])
        assert 0 <= len(next_list) <= 4
        return next_list

    def prev_kmers(self, kmer):
        """
        Returns a list of the k-mers in the graph which precede the given one. The list will have a
        length from 0 to 4.
        """
        a = 'A' + kmer.seq[:-1]
        c = 'C' + kmer.seq[:-1]
        g = 'G' + kmer.seq[:-1]
        t = 'T' + kmer.seq[:-1]
        prev_list = []
        if a in self.kmers:
            prev_list.append(self.kmers[a])
        if c in self.kmers:
            prev_list.append(self.kmers[c])
        if g in self.kmers:
            prev_list.append(self.kmers[g])
        if t in self.kmers:
            prev_list.append(self.kmers[t])
        assert 0 <= len(prev_list) <= 4
        return prev_list

    def iterate_kmers(self):
        """
        This generator yields all the k-mers in the graph in lexographical order.
        """
        for kmer in sorted(self.kmers):
            yield self.kmers[kmer]

    def reverse(self, kmer):
        """
        Takes a Kmer object and returns its reverse complement Kmer object.
        """
        return self.kmers[reverse_complement(kmer.seq)]
