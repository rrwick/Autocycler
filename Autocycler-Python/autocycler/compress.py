"""
Copyright 2024 Ryan Wick (rrwick@gmail.com)
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

from .kmer_graph import KmerGraph
from .misc import find_all_assemblies
from .unitig_graph import UnitigGraph

def compress(args):
    assemblies = find_all_assemblies(args.in_dir)

    kmer_graph = KmerGraph(args.kmer)
    kmer_graph.add_assemblies(assemblies)

    unitig_graph = UnitigGraph(kmer_graph)
    unitig_graph.save_gfa(args.out_gfa)
