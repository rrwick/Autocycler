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

import gzip
import multiprocessing
import pathlib
import sys


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


def get_default_thread_count():
    return min(multiprocessing.cpu_count(), 16)
