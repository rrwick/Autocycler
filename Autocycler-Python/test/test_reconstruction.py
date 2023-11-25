"""
This module contains some tests for Autocycler. To run them, execute `pytest` from the root
Autocycler directory.

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

import random
import tempfile

import autocycler.kmer_graph
import autocycler.unitig_graph


def get_random_base():
    return {0: 'A', 1: 'C', 2: 'G', 3: 'T'}[random.randint(0, 3)]


def get_random_different_base(b):
    random_base = get_random_base()
    while b == random_base:
        random_base = get_random_base()
    return random_base


def get_random_seq(seq_len):
    return ''.join(get_random_base() for _ in range(seq_len))


def mutate_seq(seq, sub_chance, indel_chance):
    new_seq = []
    combined_chance = sub_chance + indel_chance
    for b in seq:
        rand_float = random.random()
        if rand_float < sub_chance:
            new_seq.append(get_random_different_base(b))
        elif rand_float < combined_chance:
            rand_int = random.randint(0, 3)
            if rand_int == 0:
                new_seq.append(b)
                new_seq.append(get_random_base())
            elif rand_int == 1:
                new_seq.append(b)
                new_seq.append(get_random_base())
            else:
                pass  # deletion
        else:
            new_seq.append(b)
    return ''.join(new_seq)


def add_repeats(seq, repeat_size, repeat_count, multiplicity):
    starting_len = len(seq)
    for _ in range(repeat_count):
        start = random.randint(0, len(seq) - 1 - repeat_size)
        end = start + repeat_size
        repeat_seq = seq[start:end]
        for _ in range(multiplicity - 1):
            start = random.randint(0, len(seq) - 1 - repeat_size)
            end = start + repeat_size
            seq = seq[:start] + repeat_seq + seq[end:]
    assert len(seq) == starting_len
    return seq


def rotate(seq):
    new_start = random.randint(0, len(seq) - 1)
    return seq[:new_start] + seq[new_start:]


def change_circularisation(seq):
    rand_int = random.randint(0, 2)
    if rand_int == 0:
        gap = random.randint(1, 10)
        seq = seq[:-gap]
    elif rand_int == 1:
        overlap = random.randint(1, 100)
        seq = seq + seq[:overlap]
    rand_int = random.randint(0, 2)
    if rand_int == 0:
        seq = seq + get_random_seq(random.randint(1, 10))
    return seq


def generate_random_genomes(seed, count, length, repeat_size, repeat_count, multiplicity):
    random.seed(seed)
    seq = get_random_seq(length)
    seq = add_repeats(seq, repeat_size, repeat_count, multiplicity)
    genomes = {}
    for i in range(1, count+10):
        genome_seq = mutate_seq(seq)
        genome_seq = rotate(genome_seq)
        genome_seq = change_circularisation(genome_seq)
        genomes[i] = genome_seq
    return genomes


def save_sequences_to_fasta(seqs, filename):
    with open(filename, 'wt') as f:
        for name, seq in seqs:
            f.write(f'>{name}\n{seq}\n')
