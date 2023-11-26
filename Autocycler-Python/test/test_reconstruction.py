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

import filecmp
import pathlib
import random
import tempfile

import autocycler.kmer_graph
import autocycler.misc
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
    for b in seq:
        rand_float = random.random()
        if rand_float < sub_chance:
            new_seq.append(get_random_different_base(b))
        elif rand_float < sub_chance + indel_chance:
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
    seq = seq[new_start:] + seq[:new_start]
    if random.randint(0, 1) == 0:
        seq = autocycler.misc.reverse_complement(seq)
    return seq


def change_circularisation(seq, gap_chance, max_gap_size, overlap_chance, max_overlap_size):
    rand_float = random.random()
    if rand_float < gap_chance and max_gap_size > 0:
        max_gap_size = min(len(seq) - 1, max_gap_size)
        gap = random.randint(1, max_gap_size)
        return seq[:-gap]
    elif rand_float < gap_chance + overlap_chance and max_overlap_size > 0:
        overlap = random.randint(1, max_overlap_size)
        return seq + seq[:overlap]
    else:
        return seq


def add_junk_to_ends(seq, junk_chance, max_junk_size):
    if random.random() < junk_chance and max_junk_size > 0:
        seq = seq + get_random_seq(random.randint(1, max_junk_size))
    if random.random() < junk_chance and max_junk_size > 0:
        seq = get_random_seq(random.randint(1, max_junk_size)) + seq
    return seq


def generate_random_genomes(seed, seq_count, seq_length,
                            repeat_count, repeat_length, repeat_multiplicity,
                            sub_chance, indel_chance,
                            gap_chance, max_gap_size, overlap_chance, max_overlap_size,
                            junk_chance, max_junk_size):
    random.seed(seed)
    seq = get_random_seq(seq_length)
    seq = add_repeats(seq, repeat_length, repeat_count, repeat_multiplicity)
    genomes = {}
    for i in range(1, seq_count+1):
        genome_seq = mutate_seq(seq, sub_chance, indel_chance)
        genome_seq = rotate(genome_seq)
        genome_seq = change_circularisation(genome_seq, gap_chance, max_gap_size,
                                            overlap_chance, max_overlap_size)
        genome_seq = add_junk_to_ends(genome_seq, junk_chance, max_junk_size)
        genomes[str(i)] = genome_seq
    return genomes


def save_sequences_to_fasta(seqs, tmp_dir):
    filenames = []
    for name, seq in seqs.items():
        filename = tmp_dir / (name + '.fasta')
        with open(filename, 'wt') as f:
            f.write(f'>{name}\n{seq}\n')
        filenames.append(filename)
    return filenames


def full_test_protocol(seed = 0, kmer = 51,
                       seq_count = 5, seq_length = 10000,
                       repeat_count = 10, repeat_length = 100, repeat_multiplicity = 10,
                       sub_chance = 0.001, indel_chance = 0.001,
                       gap_chance = 0.25, max_gap_size = 100,
                       overlap_chance = 0.25, max_overlap_size = 100,
                       junk_chance = 0.1, max_junk_size = 100):
    """
    This function:
    * Generates random sequences
    * Builds a kmer-graph and then a unitig graph
    * Saves the unitig graph to file
    * Reconstructs the sequences from the unitig graph and checks they match the originals
    * Loads the unitig graph from file
    * Reconstructs the sequences from the loaded unitig graph and checks they match the originals
    * Saves the loaded unitig graph to file and checks that it's identical to the first graph file
    """
    seqs = generate_random_genomes(seed, seq_count, seq_length, repeat_count, repeat_length,
                                   repeat_multiplicity, sub_chance, indel_chance, gap_chance,
                                   max_gap_size, overlap_chance, max_overlap_size, junk_chance,
                                   max_junk_size)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        filenames = save_sequences_to_fasta(seqs, tmp_dir)
        kmer_graph = autocycler.kmer_graph.KmerGraph(kmer)
        kmer_graph.add_assemblies(filenames)
        unitig_graph = autocycler.unitig_graph.UnitigGraph(kmer_graph)
        gfa_filename = tmp_dir / 'graph.gfa'
        unitig_graph.save_gfa(gfa_filename)

        original_seqs = {name: seq for name, seq in seqs.items() if len(seq) >= kmer}
        reconstructed_seqs = unitig_graph.reconstruct_original_sequences()
        reconstructed_seqs = {x[0][0]: x[0][1] for _, x in reconstructed_seqs.items()}
        assert original_seqs == reconstructed_seqs

        loaded_graph = autocycler.unitig_graph.UnitigGraph(gfa_filename)
        reconstructed_seqs = loaded_graph.reconstruct_original_sequences()
        reconstructed_seqs = {x[0][0]: x[0][1] for _, x in reconstructed_seqs.items()}
        assert original_seqs == reconstructed_seqs

        gfa_filename_2 = tmp_dir / 'graph2.gfa'
        loaded_graph.save_gfa(gfa_filename_2)
        assert filecmp.cmp(gfa_filename, gfa_filename_2)


def test_k11():
    full_test_protocol(seed=0, kmer=11)


def test_k31():
    full_test_protocol(seed=1, kmer=31)


def test_k51():
    full_test_protocol(seed=2, kmer=51)


def test_k91():
    full_test_protocol(seed=3, kmer=91)


def test_k301():
    full_test_protocol(seed=4, kmer=301)


def test_high_repeats():
    full_test_protocol(seed=5, repeat_count = 200, repeat_multiplicity=100)


def test_high_divergence():
    full_test_protocol(seed=6, sub_chance = 0.1, indel_chance = 0.1)


def test_high_seq_count():
    full_test_protocol(seed=7, seq_count=5, repeat_count=50)


def test_short_seqs():
    full_test_protocol(seed=8, seq_length = 100, repeat_length = 10,
                       gap_chance = 1.0, max_gap_size = 99)
