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

class Position(object):
    """
    Position objects store the sequence, strand and position for contigs in the input assemblies.
    They are used:
    * In KmerGraph objects, where each Kmer object has one or more Position objects.
    * In UnitigGraph objects, where each Unitig has one or more starting and ending Position
      objects on both strands. These form a doubly linked list, tracing the input contig through
      the UnitigGraph.
    """
    def __init__(self, seq_id, strand, pos, unitig=None, unitig_strand=None, unitig_start_end=None):
        self.seq_id = seq_id
        self.strand = strand  # 1 for forward strand, -1 for reverse strand
        self.pos = pos  # 0-based indexing

        # Pointers to the preceding and following Position objects:
        self.prev = None
        self.next = None

        # Pointers to the associated unitig:
        self.unitig = unitig
        self.unitig_strand = unitig_strand  # 1 for forward, -1 for reverse
        self.unitig_start_end = unitig_start_end  # 0 for start, 1 for end

    def __repr__(self):
        return f'{self.seq_id}{"+" if self.strand == 1 else "-"}{self.pos}'

    def on_unitig_start(self):
        return self.unitig_start_end == 0

    def on_unitig_end(self):
        return self.unitig_start_end == 1
