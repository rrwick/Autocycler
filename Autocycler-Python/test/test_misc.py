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

import gzip
import pytest

import autocycler.misc


def test_get_compression_type_1():
    assert autocycler.misc.get_compression_type('test/test_misc/test.txt') == 'plain'


def test_get_compression_type_2():
    assert autocycler.misc.get_compression_type('test/test_misc/test.gz') == 'gz'


def test_get_compression_type_3():
    with pytest.raises(SystemExit) as e:
        autocycler.misc.get_compression_type('test/test_misc/test.bz2')
    assert 'cannot use bzip2' in str(e.value)


def test_get_compression_type_4():
    with pytest.raises(SystemExit) as e:
        autocycler.misc.get_compression_type('test/test_misc/test.zip')
    assert 'cannot use zip' in str(e.value)


def test_get_open_func_1():
    assert autocycler.misc.get_open_func('test/test_misc/test.txt') == open


def test_get_open_func_2():
    assert autocycler.misc.get_open_func('test/test_misc/test.gz') == gzip.open


def test_reverse_complement():
    assert autocycler.misc.reverse_complement('ACGACTACG') == 'CGTAGTCGT'
    assert autocycler.misc.reverse_complement('AXX???XXG') == 'CNN???NNT'


def test_reverse_complement_position():
    assert autocycler.misc.reverse_complement_position(0, 10) == 9
    assert autocycler.misc.reverse_complement_position(1, 10) == 8
    assert autocycler.misc.reverse_complement_position(8, 10) == 1
    assert autocycler.misc.reverse_complement_position(9, 10) == 0
    assert autocycler.misc.reverse_complement_position(2, 5) == 2


def test_iterate_fasta_1():
    fasta = list(autocycler.misc.iterate_fasta('test/test_misc/test.fasta'))
    assert fasta == [('A', 'info', 'TTGCCTGTAGTCGGGACCCC'),
                     ('B', 'stuff and more stuff', 'ATTCTCAGAATGGCGTAGTA'),
                     ('C', '', 'TACGCAGCTACG')]


def test_iterate_fasta_2():
    fasta = list(autocycler.misc.iterate_fasta('test/test_misc/test.fasta.gz'))
    assert fasta == [('A', 'info', 'TTGCCTGTAGTCGGGACCCC'),
                     ('B', 'stuff and more stuff', 'ATTCTCAGAATGGCGTAGTA'),
                     ('C', '', 'TACGCAGCTACG')]


def test_get_default_thread_count():
    assert 1 <= autocycler.misc.get_default_thread_count() <= 16
