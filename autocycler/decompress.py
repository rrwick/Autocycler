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

from .unitig_graph import UnitigGraph

def decompress(args):
    args.out_dir.mkdir(exist_ok=True)
    unitig_graph = UnitigGraph(args.in_gfa)
    unitig_graph.save_gfa(args.out_dir / 'temp_test.gfa')  # TEMP
    seqs = unitig_graph.reconstruct_original_sequences()
    # TODO: split output into separate files
    with open(args.out_dir / 'reconstructed.fasta', 'wt') as f:
        for name, seq in seqs.items():
            f.write(f'>{name}\n{seq}\n')
