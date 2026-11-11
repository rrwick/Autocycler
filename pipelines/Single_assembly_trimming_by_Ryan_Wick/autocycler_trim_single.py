#!/usr/bin/env python3
"""
Trim overlaps from each contig in a single assembly using Autocycler

Usage:
    autocycler_trim_single.py [options] <assembly.fasta[.gz]> <out_dir>
    autocycler_trim_single.py --help

Copyright 2026 Ryan Wick (rrwick@gmail.com)
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
import gzip
import logging
import os
import shlex
import shutil
import subprocess
import sys
import zlib
from pathlib import Path


logger = logging.getLogger('autocycler_trim')
KMER_SIZE = 51


def get_arguments(args):
    parser = argparse.ArgumentParser(description='Trim a single assembly with Autocycler',
                                     add_help=False,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('assembly', type=str, help='Input assembly FASTA (can be gzipped)')
    required_args.add_argument('out_dir', type=str,
                               help='Output directory (must not already exist)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--min_identity', type=float, default=0.75,
                              help='Minimum unitig alignment identity for trimming')

    output_args = parser.add_argument_group('Output')
    output_args.add_argument('--keep', action='store_true',
                             help='Keep intermediate files')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    return parser.parse_args(args)


def main(args=None):
    args = get_arguments(args)
    check_args(args)
    check_dependencies()
    out_dir = create_output_dir(args.out_dir)
    work_dir = out_dir / 'working_files'
    try:
        configure_logging(out_dir)
        log_file = out_dir / 'autocycler.log'
        print(file=sys.stderr)
        logger.info(f'Input assembly: {Path(args.assembly).resolve()}')
        logger.info(f'Autocycler log: {log_file}')
        records = load_fasta(args.assembly)
        logger.info(f'Input: {len(records):,} contigs; '
                    f'{sum(len(seq) for _, seq in records):,} bp')
        run_command(['autocycler', '--version'], log_file)

        work_dir.mkdir()
        partial_fasta = work_dir / 'assembly.fasta'
        trimmed_count, copied_count, output_bases = 0, 0, 0
        with partial_fasta.open('w', encoding='utf-8') as output:
            for number, (header, sequence) in enumerate(records, 1):
                name = header.split()[0]
                print(file=sys.stderr)
                logger.info(f'[{number}/{len(records)}] {name}: {len(sequence):,} bp')
                reason = passthrough_reason(sequence)
                if reason:
                    logger.warning(f'  Warning: copied unchanged ({reason})')
                    result_header, result_seq = header, sequence
                    copied_count += 1
                else:
                    contig_dir = work_dir / f'contig_{number:06d}'
                    result_header, result_seq = trim_contig(header, sequence, contig_dir,
                                                            args, log_file)
                    removed = len(sequence) - len(result_seq)
                    if removed:
                        logger.info(f'  Trimmed to {len(result_seq):,} bp '
                                    f'({removed:,} bp removed)')
                        trimmed_count += 1
                    else:
                        logger.info('  Unchanged')
                    if not args.keep:
                        shutil.rmtree(contig_dir)
                write_fasta(output, result_header, result_seq)
                output_bases += len(result_seq)

        final_fasta = out_dir / 'assembly.fasta'
        partial_fasta.replace(final_fasta)
        print(file=sys.stderr)
        logger.info(f'Output: {len(records):,} contigs; {output_bases:,} bp; '
                    f'{trimmed_count:,} trimmed; {copied_count:,} copied without processing')
        logger.info(f'Assembly: {final_fasta}')
        print(file=sys.stderr)
    except KeyboardInterrupt:
        quit_with_error('interrupted', exit_code=130)
    except (OSError, EOFError, UnicodeError, ValueError, RuntimeError, zlib.error) as error:
        quit_with_error(str(error))
    finally:
        if not args.keep:
            clean_up(work_dir)


def check_args(args):
    if not 0.0 <= args.min_identity <= 1.0:
        quit_with_error('--min_identity must be between 0 and 1 (inclusive)')
    if not Path(args.assembly).is_file():
        quit_with_error(f"input assembly does not exist or is not a file: '{args.assembly}'")
    out_dir = Path(args.out_dir)
    if out_dir.exists() or out_dir.is_symlink():
        quit_with_error(f"output directory already exists: '{out_dir}'")


def check_dependencies():
    if shutil.which('autocycler') is None:
        quit_with_error('autocycler was not found in PATH')


def create_output_dir(out_dir):
    out_dir = Path(out_dir).resolve()
    try:
        out_dir.mkdir(parents=True)
    except OSError as error:
        quit_with_error(f"could not create output directory '{out_dir}': {error}")
    return out_dir


def configure_logging(out_dir):
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(message)s'))
    file_handler = logging.FileHandler(out_dir / 'autocycler.log', encoding='utf-8')
    file_handler.setFormatter(logging.Formatter('%(asctime)s  %(message)s',
                                                datefmt='%Y-%m-%d %H:%M:%S'))
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)


def get_open_func(filename):
    with open(filename, 'rb') as input_file:
        start = input_file.read(4)
    if start.startswith(b'\x1f\x8b'):
        return gzip.open
    if start.startswith(b'BZh') or start.startswith(b'PK\x03\x04'):
        raise ValueError('unsupported compression format: use plain FASTA or gzip')
    return open


def load_fasta(filename):
    records = []
    header, sequence = None, []
    with get_open_func(filename)(filename, 'rt', encoding='utf-8') as fasta:
        for line_number, line in enumerate(fasta, 1):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(sequence)))
                header, sequence = line[1:], []
            elif header is None:
                raise ValueError(f"malformed FASTA '{filename}': sequence before the first "
                                 f'header on line {line_number}')
            else:
                sequence.append(line.strip())
    if header is not None:
        records.append((header, ''.join(sequence)))
    if not records:
        raise ValueError(f"FASTA contains no contigs: '{filename}'")
    names = set()
    for header, sequence in records:
        if not header.strip():
            raise ValueError(f"FASTA contains an unnamed contig: '{filename}'")
        name = header.split()[0]
        if not sequence:
            raise ValueError(f"FASTA contig '{name}' has an empty sequence: '{filename}'")
        if name in names:
            raise ValueError(f"duplicate FASTA contig name '{name}': '{filename}'")
        names.add(name)
    return records


def passthrough_reason(sequence):
    reasons = []
    if len(sequence) < KMER_SIZE:
        reasons.append(f'shorter than {KMER_SIZE} bp')
    if set(sequence.upper()) - set('ACGT'):
        reasons.append('contains non-ACGT bases')
    return '; '.join(reasons)


def write_fasta(output, header, sequence):
    output.write(f'>{header}\n{sequence}\n')


def trim_contig(header, sequence, contig_dir, args, log_file):
    assemblies_dir = contig_dir / 'assemblies'
    assemblies_dir.mkdir(parents=True)
    with (assemblies_dir / 'input.fasta').open('w', encoding='utf-8') as fasta:
        write_fasta(fasta, header, sequence)

    autocycler_dir = contig_dir / 'autocycler'
    run_command(['autocycler', 'compress', '-i', assemblies_dir, '-a', autocycler_dir,
                 '--kmer', KMER_SIZE, '--threads', 1], log_file)
    run_command(['autocycler', 'cluster', '-a', autocycler_dir], log_file)
    clusters = sorted(path for path in (autocycler_dir / 'clustering').glob('qc_*/cluster_*')
                      if path.is_dir())
    if len(clusters) != 1:
        raise RuntimeError(f"expected one cluster for '{header.split()[0]}', found "
                           f"{len(clusters)} (see '{log_file}')")
    cluster_dir = clusters[0]
    run_command(['autocycler', 'trim', '-c', cluster_dir,
                 '--min_identity', args.min_identity,
                 '--mad', 0, '--threads', 1], log_file)
    trimmed_dir = contig_dir / 'trimmed_sequences'
    run_command(['autocycler', 'decompress', '-i', cluster_dir / '2_trimmed.gfa',
                 '-o', trimmed_dir], log_file)
    records = load_fasta(trimmed_dir / 'input.fasta')
    if len(records) != 1 or len(records[0][1]) > len(sequence):
        raise RuntimeError(f"unexpected trimming output for '{header.split()[0]}' "
                           f"(see '{log_file}')")
    return header, records[0][1]


def run_command(command, log_file):
    command = [str(part) for part in command]
    environment = dict(os.environ, NO_COLOR='1')
    with log_file.open('a', encoding='utf-8') as log:
        log.write(f'\n$ {shlex.join(command)}\n')
        log.flush()
        result = subprocess.run(command, stdin=subprocess.DEVNULL, stdout=log,
                                stderr=subprocess.STDOUT, env=environment)
    if result.returncode != 0:
        raise RuntimeError(f"Autocycler {command[1]} failed with exit code {result.returncode} "
                           f"(see '{log_file}')")


def clean_up(work_dir):
    if work_dir.exists():
        try:
            shutil.rmtree(work_dir)
        except OSError as error:
            logger.warning(f"Warning: could not remove working directory '{work_dir}': {error}")


def quit_with_error(message, exit_code=1):
    if logger.handlers:
        logger.error(f'Error: {message}')
    else:
        print(f'Error: {message}', file=sys.stderr)
    sys.exit(exit_code)


if __name__ == '__main__':
    main()
