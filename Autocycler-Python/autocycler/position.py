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

from typing import Optional


class Position(object):
    """
    Position objects store the sequence, strand and position for contigs in the input assemblies.
    They are used:
    * In KmerGraph objects, where each Kmer object has one or more Position objects.
    * In UnitigGraph objects, where each Unitig has one or more Position objects on both strands.
    """
    def __init__(self, seq_id: int, strand: int, pos: int):
        self.seq_id = seq_id
        self.strand = strand  # 1 for forward strand, -1 for reverse strand
        self.pos = pos  # 0-based indexing

    def __repr__(self):
        return f'{self.seq_id}{"+" if self.strand == 1 else "-"}{self.pos}'
