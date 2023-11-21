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

import argparse
import pathlib
import sys

from .help_formatter import MyParser, MyHelpFormatter
from .kmer_graph import KmerGraph
from .misc import find_all_assemblies, get_default_thread_count
from .unitig_graph import UnitigGraph

__version__ = '0.0.0'


def get_arguments(args):
    parser = MyParser(description='Autocycler', add_help=False,
                      formatter_class=MyHelpFormatter)

    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-i', '--in_dir', type=pathlib.Path, required=True,
                               help='Directory containing input assemblies')
    required_args.add_argument('-o', '--out_dir', type=pathlib.Path, required=True,
                               help='Directory where consensus assembly will be constructed')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('-k', '--kmer', type=int, default=91,
                              help='K-mer size for De Bruijn graph (default: DEFAULT)')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of CPU threads (default: DEFAULT)')
    setting_args.add_argument('--verbose', action='store_true',
                              help='Display more output information')

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='Autocycler v' + __version__,
                            help="Show program's version number and exit")

    args = parser.parse_args(args)
    return args


def main(args=None):
    args = get_arguments(args)
    check_args(args)
    args.out_dir.mkdir(exist_ok=True)
    assemblies = find_all_assemblies(args.in_dir)

    kmer_graph = KmerGraph(args.kmer)
    kmer_graph.add_assemblies(assemblies)

    unitig_graph = UnitigGraph(kmer_graph)
    unitig_graph.save_gfa(args.out_dir / f'001_unitig_graph.gfa')

    seqs = unitig_graph.reconstruct_original_sequences()  # TEMP
    with open(args.out_dir / 'reconstructed.fasta', 'wt') as f:  # TEMP
        for name, seq in seqs.items():  # TEMP
            f.write(f'>{name}\n{seq}\n')  # TEMP


def check_args(args):
    if args.kmer < 11:
        sys.exit('Error: --kmer must be 11 or greater')
    if args.kmer % 2 == 0:
        sys.exit('Error: --kmer must be odd')


if __name__ == '__main__':
    main()
