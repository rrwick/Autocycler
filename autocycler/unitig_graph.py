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

from .unitig import Unitig


class UnitigGraph(object):
    """
    This class builds a unitig graph from a k-mer graph, where all nonbranching paths are merged
    into unitigs. This simplifies things and saves memory.
    """
    def __init__(self, k_graph):
        self.unitigs = []
        self.k_size = k_graph.k_size

        self.contig_ids = sorted(k_graph.id_to_contig_info.keys())
        self.contig_ids_to_assembly = {i: k_graph.id_to_contig_info[i][0]for i in self.contig_ids}
        self.contig_ids_to_header = {i: k_graph.id_to_contig_info[i][1]for i in self.contig_ids}
        self.contig_ids_to_seq_len = {i: k_graph.id_to_contig_info[i][2]for i in self.contig_ids}

        print('\nBuilding unitig graph from k-mer graph:')
        self.build_unitigs_from_kmer_graph(k_graph)
        self.renumber_unitigs()
        self.create_links()
        self.trim_overlaps()
        self.connect_positions()

    def create_from_kmer_graph(self):
        pass
        # TODO
        # TODO
        # TODO

    def create_from_gfa_file(self):
        pass
        # TODO
        # TODO
        # TODO

    def save_gfa(self, gfa_filename):
        with open(gfa_filename, 'wt') as f:
            f.write('H\tVN:Z:1.0\n')
            f.write(f'# Unitig graph built using {self.k_size}-mers from these sequences:\n')
            for i in self.contig_ids:
                f.write(f'# {i}: {self.contig_ids_to_seq_len[i]} bp, '
                        f'{self.contig_ids_to_assembly[i]}, {self.contig_ids_to_header[i]}\n')
            links = []
            for unitig in self.unitigs:
                f.write(unitig.gfa_segment_line())
                for next_unitig, next_strand in unitig.forward_next:
                    links.append((unitig.number, '+',
                                  next_unitig.number, '+' if next_strand == 1 else '-'))
                for next_unitig, next_strand in unitig.reverse_next:
                    links.append((unitig.number, '-',
                                  next_unitig.number, '+' if next_strand == 1 else '-'))
            for a, a_strand, b, b_strand in links:
                f.write(f'L\t{a}\t{a_strand}\t{b}\t{b_strand}\t0M\n')

    def build_unitigs_from_kmer_graph(self, kmer_graph):
        seen = set()
        unitig_number = 0
        half_k = self.k_size // 2
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
                if forward_kmer.first_position(half_k) or reverse_kmer.first_position(half_k):
                    break
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
                if forward_kmer.first_position(half_k) or reverse_kmer.first_position(half_k):
                    break
                unitig.add_kmer_to_start(forward_kmer, reverse_kmer)
                seen.add(forward_kmer)
                seen.add(reverse_kmer)

            unitig.simplify_seqs()
            self.unitigs.append(unitig)
        print(f'  {len(self.unitigs)} unitigs')

    def renumber_unitigs(self):
        overlap_size = self.k_size // 2
        self.unitigs = sorted(self.unitigs,
                              key=lambda u: (1.0/u.length(), u.forward_seq[overlap_size:]))
        unitig_number = 0
        for unitig in self.unitigs:
            unitig_number += 1
            unitig.number = unitig_number

    def create_links(self):
        piece_len = self.k_size - 1

        # Index unitigs by their k-1 starting and ending sequences.
        starting_forward = collections.defaultdict(list)
        starting_reverse = collections.defaultdict(list)
        for unitig in self.unitigs:
            starting_forward[unitig.forward_seq[:piece_len]].append(unitig)
            starting_reverse[unitig.reverse_seq[:piece_len]].append(unitig)

        # Use the indices to find connections between unitigs.
        link_count = 0
        for unitig in self.unitigs:
            ending_forward_seq = unitig.forward_seq[-piece_len:]
            for next_unitig in starting_forward[ending_forward_seq]:
                # unitig+ -> next_unitig+
                unitig.forward_next.append((next_unitig, 1))
                next_unitig.forward_prev.append((unitig, 1))
                link_count += 1
                # next_unitig- -> unitig-
                next_unitig.reverse_next.append((unitig, -1))
                unitig.reverse_prev.append((next_unitig, -1))
                link_count += 1
            for next_unitig in starting_reverse[ending_forward_seq]:
                # unitig+ -> next_unitig-
                unitig.forward_next.append((next_unitig, -1))
                next_unitig.reverse_prev.append((unitig, 1))
                link_count += 1
            ending_reverse_seq = unitig.reverse_seq[-piece_len:]
            for next_unitig in starting_forward[ending_reverse_seq]:
                # unitig- -> next_unitig+
                unitig.reverse_next.append((next_unitig, 1))
                next_unitig.forward_prev.append((unitig, -1))
                link_count += 1
        print(f'  {link_count} links')

    def trim_overlaps(self):
        for unitig in self.unitigs:
            unitig.trim_overlaps(self.k_size)

    def connect_positions(self):
        """
        This method connects up all Position objects in the UnitigGraph. This happens:
        * Within unitigs, where a starting position is connected to an end position (carried out by
          the Unitig connect_positions method).
        * Between unitigs, where an ending position of one unitig is connected to a starting
          position of a linked unitig (carried out by this method).
        """
        for unitig in self.unitigs:
            unitig.connect_positions()

        links = []
        for unitig in self.unitigs:
            for next_unitig, next_strand in unitig.forward_next:
                links.append((unitig, 1, next_unitig, next_strand))
            for next_unitig, next_strand in unitig.reverse_next:
                links.append((unitig, -1, next_unitig, next_strand))

        for a, a_strand, b, b_strand in links:
            if a_strand == 1 and b_strand == 1:
                first, second = a.forward_end_positions, b.forward_start_positions
            elif a_strand == 1 and b_strand == -1:
                first, second = a.forward_end_positions, b.reverse_start_positions
            elif a_strand == -1 and b_strand == 1:
                first, second = a.reverse_end_positions, b.forward_start_positions
            elif a_strand == -1 and b_strand == -1:
                first, second = a.reverse_end_positions, b.reverse_start_positions
                pass
            else:
                assert False
            for a in first:
                matches = [b for b in second
                        if a.seq_id == b.seq_id and a.strand == b.strand and a.pos == b.pos]
                assert len(matches) <= 1
                if len(matches) == 1:
                    match = matches[0]
                    if a.next is None:
                        a.next = match
                    else:
                        assert a.next == match
                    if match.prev is None:
                        match.prev = a
                    else:
                        assert match.prev_kmer == a

    def reconstruct_original_sequences(self):
        """
        Returns a dictionary of the original sequences used to build the unitig graph.
        """
        print(f'\nReconstructing input assemblies from unitig graph')
        return {i: self.reconstruct_original_sequence(i) for i in self.contig_ids}

    def reconstruct_original_sequence(self, seq_id):
        """
        Returns the original sequence for the given sequence ID. Works by finding the first
        forward-strand position for that sequence and then following the position-to-position
        links through the unitig graph, building sequence as it goes.
        """
        total_length = self.contig_ids_to_seq_len[seq_id]
        print(f'  {seq_id}: {total_length} bp')
        sequence = []
        p = self.find_first_position(seq_id)
        assert p.on_unitig_end()
        starting_seq = p.unitig.get_seq(p.unitig_strand, upstream=self.k_size//2)
        assert p.pos <= len(starting_seq)
        sequence.append(starting_seq[-p.pos:])
        while p.next is not None:
            p = p.next
            assert p.on_unitig_start()
            if p.next is not None:
                p = p.next
                assert p.on_unitig_end()
                sequence.append(p.unitig.get_seq(p.unitig_strand))
            else:  # reached the last position
                remaining_length = total_length - sum(len(s) for s in sequence)
                final_seq = p.unitig.get_seq(p.unitig_strand, downstream=self.k_size//2)
                assert remaining_length <= len(final_seq)
                sequence.append(final_seq[:remaining_length])
        sequence = ''.join(sequence)
        assert len(sequence) == total_length
        return ''.join(sequence)

    def find_first_position(self, i):
        """
        This method looks through all of the unitig-end positions (both strands), looking for the
        first instance of the positive strand for the given sequence ID.
        """
        best = None
        for unitig in self.unitigs:
            for p in unitig.forward_end_positions:
                if p.seq_id == i and p.strand == 1 and (best is None or p.pos < best.pos):
                    best = p
            for p in unitig.reverse_end_positions:
                if p.seq_id == i and p.strand == 1 and (best is None or p.pos < best.pos):
                    best = p
        return best
