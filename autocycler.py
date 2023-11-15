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

import argparse
import collections
import gzip
import multiprocessing
import os
import pathlib
import shutil
import statistics
import subprocess
import sys

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
    unitig_graph.save_gfa(args.out_dir / f'001_de_bruijn_graph.gfa')


def check_args(args):
    if args.kmer < 11:
        sys.exit('Error: --kmer must be 11 or greater')
    if args.kmer % 2 == 0:
        sys.exit('Error: --kmer must be odd')


def find_all_assemblies(in_dir):
    print(f'\nLooking for assembly files in {in_dir}...', flush=True, end='')
    all_assemblies = [str(x) for x in sorted(pathlib.Path(in_dir).glob('**/*'))
                      if x.is_file()]
    all_assemblies = [x for x in all_assemblies if
                      x.endswith('.fasta') or x.endswith('.fasta.gz') or
                      x.endswith('.fna') or x.endswith('.fna.gz') or
                      x.endswith('.fa') or x.endswith('.fa.gz')]
    # TODO: also look for GFA-format assemblies
    plural = 'assembly' if len(all_assemblies) == 1 else 'assemblies'
    print(f' found {len(all_assemblies)} {plural}')
    if len(all_assemblies) == 0:
        sys.exit(f'Error: no assemblies found in {in_dir}')
    return sorted(all_assemblies)


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
        seq_id = 0
        for assembly in assemblies:
            print(f'\nAdding {assembly} to graph:')
            for name, _, seq in iterate_fasta(assembly):
                seq_id += 1
                print(f'  {seq_id}: {name} ({len(seq)} bp)...', flush=True, end='')
                self.add_sequence(seq, seq_id)
                print(' done')
        print(f'\nGraph contains {len(self.kmers)} k-mers')

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

    def previous_kmers(self, kmer):
        """
        Returns a list of the k-mers in the graph which precede the given one. The list will have a
        length from 0 to 4.
        """
        a = 'A' + kmer.seq[:-1]
        c = 'C' + kmer.seq[:-1]
        g = 'G' + kmer.seq[:-1]
        t = 'T' + kmer.seq[:-1]
        previous_list = []
        if a in self.kmers:
            previous_list.append(self.kmers[a])
        if c in self.kmers:
            previous_list.append(self.kmers[c])
        if g in self.kmers:
            previous_list.append(self.kmers[g])
        if t in self.kmers:
            previous_list.append(self.kmers[t])
        assert 0 <= len(previous_list) <= 4
        return previous_list

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


class Kmer(object):
    def __init__(self, seq):
        self.seq = seq
        self.positions = []

    def __repr__(self):
        positions = ','.join([str(p) for p in sorted(self.positions)])
        return f'{positions}'

    def add_position(self, seq_id, strand, pos):
        self.positions.append(KmerPosition(seq_id, strand, pos))

    def count(self):
        return len(self.positions)


class KmerPosition(object):
    def __init__(self, seq_id, strand, pos):
        self.seq_id = seq_id
        self.strand = strand  # 1 for forward strand, -1 for reverse strand
        self.pos = pos  # position of the k-mer's centre base (0-based indexing)

    def __repr__(self):
        return f'{self.seq_id}{"+" if self.strand == 1 else "-"}{self.pos}'


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
        self.create_links()
        self.trim_overlaps()

    def save_gfa(self, gfa_filename):
        with open(gfa_filename, 'wt') as f:
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
                previous_kmers = kmer_graph.previous_kmers(forward_kmer)
                if len(previous_kmers) != 1:
                    break
                reverse_kmer = kmer_graph.reverse(forward_kmer)
                unitig.add_kmer_to_end(forward_kmer, reverse_kmer)
                seen.add(forward_kmer)
                seen.add(reverse_kmer)

            # Extend unitig backward
            forward_kmer = starting_kmer
            while True:
                previous_kmers = kmer_graph.previous_kmers(forward_kmer)
                if len(previous_kmers) != 1:
                    break
                forward_kmer = previous_kmers[0]
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
            forward_depths.append(len(k.positions))
        for k in self.reverse_kmers:
            reverse_depths.append(len(k.positions))
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


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def reverse_complement_position(pos, seq_len):
    """
    Returns the position of a base on the reverse complement strand using 0-based indexing.
    """
    return seq_len - pos - 1


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def iterate_fasta(filename):
    """
    Takes a FASTA file as input and yields the contents as (name, info, seq) tuples.
    """
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    name_parts = name.split(maxsplit=1)
                    contig_name = name_parts[0]
                    info = '' if len(name_parts) == 1 else name_parts[1]
                    yield contig_name, info, ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            name_parts = name.split(maxsplit=1)
            contig_name = name_parts[0]
            info = '' if len(name_parts) == 1 else name_parts[1]
            yield contig_name, info, ''.join(sequence)


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
DIM = '\033[2m'


class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if the script is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """
    def error(self, message):
        if len(sys.argv) == 1:  # if no arguments were given.
            self.print_help(file=sys.stderr)
            sys.exit(1)
        else:
            super().error(message)


class MyHelpFormatter(argparse.HelpFormatter):

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        self.colours = get_colours_from_tput()
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if 'default: DEFAULT' in help_text:
                help_text = help_text.replace('default: DEFAULT', f'default: {action.default}')
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = BOLD + heading + END_FORMATTING
        super().start_section(heading)

    def _format_action(self, action):
        """
        Override this method to make help descriptions dim.
        """
        help_position = min(self._action_max_length + 2, self._max_help_position)
        help_width = self._width - help_position
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)
        if not action.help:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = 0
        elif len(action_header) <= action_width:
            tup = self._current_indent, '', action_width, action_header
            action_header = '%*s%-*s  ' % tup
            indent_first = 0
        else:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = help_position
        parts = [action_header]
        if action.help:
            help_text = self._expand_help(action)
            help_lines = self._split_lines(help_text, help_width)
            first_line = help_lines[0]
            if self.colours > 8:
                first_line = DIM + first_line + END_FORMATTING
            parts.append('%*s%s\n' % (indent_first, '', first_line))
            for line in help_lines[1:]:
                if self.colours > 8:
                    line = DIM + line + END_FORMATTING
                parts.append('%*s%s\n' % (help_position, '', line))
        elif not action_header.endswith('\n'):
            parts.append('\n')
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))
        return self._join_parts(parts)


def get_colours_from_tput():
    try:
        return int(subprocess.check_output(['tput', 'colors']).decode().strip())
    except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
        return 1


def get_default_thread_count():
    return min(multiprocessing.cpu_count(), 16)


if __name__ == '__main__':
    main()


# Pytest tests
# ------------

def test_reverse_complement():
    assert reverse_complement('ACGACTACG') == 'CGTAGTCGT'
    assert reverse_complement('AXX???XXG') == 'CNN???NNT'


def test_reverse_complement_position():
    assert reverse_complement_position(0, 10) == 9
    assert reverse_complement_position(1, 10) == 8
    assert reverse_complement_position(8, 10) == 1
    assert reverse_complement_position(9, 10) == 0
    assert reverse_complement_position(2, 5) == 2
