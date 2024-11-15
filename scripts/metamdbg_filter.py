#!/usr/bin/env python3
"""
This script removes low-depth contigs from a metaMDBG assembly. The user can supply a threshold
with --min_depth, or else the default threshold will be 1/4 the depth of the longest contig.

Copyright 2024 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Autocycler

This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Autocycler is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Autocycler.
If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import gzip
import re
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='metaMDBG contig depth filter')
    parser.add_argument('input', type=str,
                        help='Filename of metaMDBG assembly in FASTA format (can be gzipped)')
    parser.add_argument('-d', '--min_depth', type=float,
                        help='Exclude contigs with a depth less than this value '
                             '(default: 1/4 the depth of the longest contig)')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    assembly = load_fasta(args.input)
    if len(assembly) == 0:
        sys.exit('Error: no contigs in assembly')
    headers = [h for h, _ in assembly]
    if args.min_depth is None:
        args.min_depth = auto_depth_threshold(headers)
    for header, seq in assembly:
        depth = length_and_depth(header)[1]
        if depth >= args.min_depth:
            print(f'>{header}')
            print(f'{seq}')


def auto_depth_threshold(headers):
    # If the user didn't supply a depth threshold, use 25% of the depth of the longest contig.
    lengths_and_depths = [length_and_depth(h) for h in headers]
    _, max_depth = max(lengths_and_depths, key=lambda x: x[0])
    return max_depth / 4


def length_and_depth(header):
    match = re.search(r'length=(\d+).*coverage=(\d+)', header)
    if match:
        return int(match.group(1)), float(match.group(2))
    else:
        sys.exit('Error: no length or coverage values found in header')


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    https://stackoverflow.com/questions/13044562
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
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_fasta(fasta_filename):
    fasta_seqs = []
    with get_open_func(fasta_filename)(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name, ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            fasta_seqs.append((name, ''.join(sequence)))
    return fasta_seqs


if __name__ == '__main__':
    main()


# Unit tests for Pytest
# =====================

def test_length_and_depth():
    assert length_and_depth('ctg0 length=5113 coverage=2 circular=no') == (5113, 2.0)
    assert length_and_depth('ctg1 length=8684 coverage=381 circular=yes') == (8684, 381.0)
    assert length_and_depth('ctg2 length=94195 coverage=55 circular=yes') == (94195, 55.0)
    assert length_and_depth('ctg3 length=143270 coverage=65 circular=yes') == (143270, 65.0)


def test_auto_depth_threshold():
    assert auto_depth_threshold(['ctg0 length=5113 coverage=2 circular=no',
                                 'ctg1 length=8684 coverage=381 circular=yes',
                                 'ctg2 length=94195 coverage=55 circular=yes',
                                 'ctg3 length=143270 coverage=65 circular=yes']) == 65.0 / 4.0
