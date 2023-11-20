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

import pytest

from autocycler.kmer_graph import Kmer
import autocycler.unitig


def build_unitig_1():
    # 5mers, CAACTTAG
    u = autocycler.unitig.Unitig(1, Kmer('AACTT'), Kmer('AAGTT'))
    u.add_kmer_to_end(Kmer('ACTTA'), Kmer('TAAGT'))
    u.add_kmer_to_end(Kmer('CTTAG'), Kmer('CTAAG'))
    u.add_kmer_to_start(Kmer('CAACT'), Kmer('AGTTG'))
    u.simplify_seqs()
    return u


def build_unitig_2():
    # 5mers, TTAGACTCA
    u = autocycler.unitig.Unitig(1, Kmer('AGACT'), Kmer('AGTCT'))
    u.add_kmer_to_end(Kmer('GACTC'), Kmer('GAGTC'))
    u.add_kmer_to_end(Kmer('ACTCA'), Kmer('TGAGT'))
    u.add_kmer_to_start(Kmer('TAGAC'), Kmer('GTCTA'))
    u.add_kmer_to_start(Kmer('TTAGA'), Kmer('TCTAA'))
    u.simplify_seqs()
    return u


def build_two_linked_unitig_same_strand():
    # 5mers, CAACTTAG, TTAGACTCA
    u_1 = build_unitig_1()
    u_2 = build_unitig_2()
    u_1.forward_next.append((u_2, 1))
    u_2.forward_prev.append((u_1, 1))
    u_2.reverse_next.append((u_1, -1))
    u_1.reverse_prev.append((u_2, -1))
    return u_1, u_2


def test_get_seq():
    u = build_unitig_1()
    assert u.get_seq(1) == 'CAACTTAG'
    assert u.get_seq(-1) == 'CTAAGTTG'
    u = build_unitig_2()
    assert u.get_seq(1) == 'TTAGACTCA'
    assert u.get_seq(-1) == 'TGAGTCTAA'
    with pytest.raises(AssertionError):
        u.get_seq(0)


def test_trim_overlaps_1():
    u = build_unitig_1()
    u.trim_overlaps(5)  # does nothing because this unitig has dead ends
    assert u.get_seq(1) == 'CAACTTAG'
    assert u.get_seq(-1) == 'CTAAGTTG'


def test_trim_overlaps_2():
    u_1, u_2 = build_two_linked_unitig_same_strand()
    u_1.trim_overlaps(5)
    u_2.trim_overlaps(5)
    assert u_1.get_seq(1) == 'CAACTT'
    assert u_2.get_seq(1) == 'AGACTCA'
    assert u_1.get_seq(-1) == 'AAGTTG'
    assert u_2.get_seq(-1) == 'TGAGTCT'


def test_length():
    u = build_unitig_1()
    assert u.length() == 8
    u = build_unitig_2()
    assert u.length() == 9


def test_depth():
    u = build_unitig_1()
    assert u.depth == 0.0
