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

from .compress import compress
from .decompress import decompress
from .help_formatter import MyParser, MyHelpFormatter
from .misc import get_default_thread_count
from .resolve import resolve

__version__ = '0.0.0'

def main():
    args = parse_args(sys.argv[1:])

    if args.subparser_name == 'compress':
        check_compress_args(args)
        compress(args)
    elif args.subparser_name == 'decompress':
        check_decompress_args(args)
        decompress(args)
    elif args.subparser_name == 'resolve':
        check_resolve_args(args)
        resolve(args)


def parse_args(args):
    description = 'Autocycler: a tool for producing consensus bacterial genome assemblies'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    compress_subparser(subparsers)
    decompress_subparser(subparsers)
    resolve_subparser(subparsers)

    longest_choice_name = max(len(c) for c in subparsers.choices)
    subparsers.help = 'R|'
    for choice, choice_parser in subparsers.choices.items():
        padding = ' ' * (longest_choice_name - len(choice))
        subparsers.help += choice + ': ' + padding
        d = choice_parser.description
        subparsers.help += d[0].lower() + d[1:]  # don't capitalise the first letter
        subparsers.help += '\n'

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='Autocycler v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def compress_subparser(subparsers):
    group = subparsers.add_parser('compress',
                                  description='compress input assemblies into a De Bruijn graph',
                                  formatter_class=MyHelpFormatter, add_help=False)
    
    required_args = group.add_argument_group('Required')
    required_args.add_argument('-i', '--in_dir', type=pathlib.Path, required=True,
                               help='Directory containing input assemblies')
    required_args.add_argument('-o', '--out_gfa', type=pathlib.Path, required=True,
                               help='Filename of output assembly graph')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('-k', '--kmer', type=int, default=91,
                              help='K-mer size for De Bruijn graph (default: DEFAULT)')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of CPU threads (default: DEFAULT)')
    setting_args.add_argument('--verbose', action='store_true',
                              help='Display more output information')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Autocycler v' + __version__,
                            help="Show program's version number and exit")


def decompress_subparser(subparsers):
    group = subparsers.add_parser('decompress',
                                  description='decompress De Bruijn graph back into assemblies',
                                  formatter_class=MyHelpFormatter, add_help=False)
    
    required_args = group.add_argument_group('Required')
    required_args.add_argument('-i', '--in_gfa', type=pathlib.Path, required=True,
                               help='Autocycler GFA file')
    required_args.add_argument('-o', '--out_dir', type=pathlib.Path, required=True,
                               help='Directory where decompressed assemblies will be saved')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of CPU threads (default: DEFAULT)')
    setting_args.add_argument('--verbose', action='store_true',
                              help='Display more output information')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Autocycler v' + __version__,
                            help="Show program's version number and exit")


def resolve_subparser(subparsers):
    group = subparsers.add_parser('resolve',
                                  description='resolve De Bruijn graph to a consensus assembly',
                                  formatter_class=MyHelpFormatter, add_help=False)
    
    required_args = group.add_argument_group('Required')
    required_args.add_argument('-i', '--in_gfa', type=pathlib.Path, required=True,
                               help='Autocycler GFA file')
    required_args.add_argument('-o', '--out_dir', type=pathlib.Path, required=True,
                               help='Directory where resolved assembly graphs will be saved')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of CPU threads (default: DEFAULT)')
    setting_args.add_argument('--verbose', action='store_true',
                              help='Display more output information')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Autocycler v' + __version__,
                            help="Show program's version number and exit")


def check_compress_args(args):
    if args.kmer < 11:
        sys.exit('Error: --kmer must be 11 or greater')
    if args.kmer % 2 == 0:
        sys.exit('Error: --kmer must be odd')


def check_decompress_args(args):
    pass


def check_resolve_args(args):
    pass


if __name__ == '__main__':
    main()
