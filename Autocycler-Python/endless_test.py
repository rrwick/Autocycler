#!/usr/bin/env python3
"""
This script run endless random tests on Autocycler's graph construction, graph loading from file
and sequence reconstruction.

Run it from Autocycler's root directory like this:
./endless_test.py

Copyright 2023 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Autocycler

This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Autocycler is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Autocycler.
If not, see <https://www.gnu.org/licenses/>.
"""

from test.test_reconstruction import full_test_protocol

import random


def main():
    i = 0
    while True:
        i += 1
        print(f'\n\n\n\n\n\n\n\nTEST {i}')
        print('--------------------------------------------------------------------------------')
        seed = random.randint(0, 1000000)
        kmer = 2 * random.randint(8, 50) + 1
        seq_count =int(round(random.uniform(1, 5) ** 2))
        seq_length = int(round(random.uniform(10, 300) ** 2))
        repeat_count = int(round(random.uniform(0, 10) ** 2))
        repeat_length = random.randint(1, seq_length//10)
        repeat_multiplicity = int(round(random.uniform(0, 10) ** 2))
        sub_chance = random.uniform(0.0, 0.4) ** 3
        indel_chance = random.uniform(0.0, 0.4) ** 3
        gap_chance = random.uniform(0.0, 1.0) ** 3
        max_gap_size = int(round(random.uniform(0, 10) ** 3))
        overlap_chance = random.uniform(0.0, 1.0) ** 3
        max_overlap_size = int(round(random.uniform(0, 10) ** 3))
        junk_chance = random.uniform(0.0, 1.0) ** 3
        max_junk_size = int(round(random.uniform(0, 10) ** 3))
        print(f'seed:                {seed}')
        print(f'kmer:                {kmer}')
        print(f'seq_count:           {seq_count}')
        print(f'seq_length:          {seq_length}')
        print(f'repeat_count:        {repeat_count}')
        print(f'repeat_length:       {repeat_length}')
        print(f'repeat_multiplicity: {repeat_multiplicity}')
        print(f'sub_chance:          {sub_chance}')
        print(f'indel_chance:        {indel_chance}')
        print(f'gap_chance:          {gap_chance}')
        print(f'max_gap_size:        {max_gap_size}')
        print(f'overlap_chance:      {overlap_chance}')
        print(f'max_overlap_size:    {max_overlap_size}')
        print(f'junk_chance:         {junk_chance}')
        print(f'max_junk_size:       {max_junk_size}')
        full_test_protocol(seed=seed, kmer=kmer, seq_count=seq_count, seq_length=seq_length,
                           repeat_count=repeat_count, repeat_length=repeat_length,
                           repeat_multiplicity=repeat_multiplicity, sub_chance=sub_chance,
                           indel_chance=indel_chance, gap_chance=gap_chance,
                           max_gap_size=max_gap_size, overlap_chance=overlap_chance,
                           max_overlap_size=max_overlap_size, junk_chance=junk_chance,
                           max_junk_size=max_junk_size)
        print('\nTEST PASSED!')


if __name__ == '__main__':
    main()
