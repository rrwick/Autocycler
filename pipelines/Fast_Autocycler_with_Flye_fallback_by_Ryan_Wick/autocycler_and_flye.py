#!/usr/bin/env python3
"""
This pipeline builds both Flye and Autocycler assemblies. It uses Autocycler when it was
successful, otherwise falling back to Flye.

Usage:
    autocycler_and_flye.py [options] <reads.fastq[.gz]> <out_dir> 
    autocycler_and_flye.py --help

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
import csv
import gzip
import logging
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


logger = logging.getLogger('autocycler_and_flye')

AUTOCYCLER_ASSEMBLERS = ('plassembler', 'raven', 'myloasm', 'miniasm', 'flye', 'metamdbg')

DEPENDENCIES = ('autocycler', 'flye', 'lrge', 'metaMDBG', 'miniasm', 'minimap2', 'minipolish',
                'myloasm', 'nice', 'parallel', 'plassembler', 'racon', 'rasusa', 'raven')


def get_arguments(args):
    parser = argparse.ArgumentParser(description='Fast Autocycler with Flye fallback',
                                     add_help=False,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('reads', type=str,
                               help='Read FASTQ file (can be gzipped)')
    required_args.add_argument('out_dir', type=str,
                               help='Directory for working files and final assembly (will be created)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--read-type',
                              choices=('ont_r9', 'ont_r10', 'pacbio_clr', 'pacbio_hifi'),
                              default='ont_r10',
                              help='Type of long reads')
    setting_args.add_argument('--genome_size', type=int,
                              help='Genome size in bp (skips LRGE when supplied)')
    setting_args.add_argument('--min-size-ratio', type=float, default=0.75,
                              help='Reject assemblies smaller than this multiple of genome size')
    setting_args.add_argument('--max-size-ratio', type=float, default=1.25,
                              help='Reject assemblies larger than this multiple of genome size')
    setting_args.add_argument('--seed', type=int, default=0,
                              help='Random seed for reproducible read subsampling')
    setting_args.add_argument('--subset_count', type=int, default=2,
                              help='Number of Autocycler read subsets')

    resources_args = parser.add_argument_group('Resources')
    resources_args.add_argument('--threads', type=int, default=16,
                                help='Maximum number of CPU threads')
    resources_args.add_argument('--jobs', type=int, default=4,
                                help='Number of simultaneous Autocycler input assembly jobs')
    resources_args.add_argument('--max-job-time', type=str, default='4h',
                                help='Maximum runtime for each Autocycler input assembly job')

    output_args = parser.add_argument_group('Output')
    output_args.add_argument('--clean', type=int, choices=range(4), default=2,
                             help='Cleanup level from 0 (keep everything) to 3 '
                                  '(keep only final assembly and logs)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')

    args = parser.parse_args(args)
    return args


def main(args=None):
    args = get_arguments(args)
    check_args(args)
    check_dependencies(skip_lrge=args.genome_size is not None)
    out_dir = create_output_dir(args.out_dir)
    configure_logging(out_dir)
    input_read_stats = get_fastq_stats(args.reads)
    logger.info(format_read_stats('Input reads', input_read_stats))
    if args.genome_size is None:
        genome_size = estimate_genome_size_with_lrge(args.reads, out_dir, args.threads, args.seed,
                                                     args.read_type)
    else:
        genome_size = args.genome_size
        logger.info(f'Genome size: {genome_size:,} bp (user-supplied)')
    logger.info(f'Estimated read depth: {input_read_stats[1] / genome_size:.1f}×')
    try:
        rasusa_reads = subsample_with_rasusa(args.reads, out_dir, genome_size,
                                             input_read_stats[1], args.seed)
        flye_fasta = flye_assembly(rasusa_reads, out_dir, args.threads, args.read_type)
        flye_successful = check_flye_assembly(flye_fasta, genome_size,
                                              args.min_size_ratio, args.max_size_ratio)
        autocycler_fasta = autocycler_assembly(args.reads, flye_fasta, out_dir, genome_size,
                                               args.threads, args.jobs, args.subset_count,
                                               args.max_job_time, args.seed, args.read_type)
        create_autocycler_metrics(args.reads, out_dir)
        autocycler_successful = check_autocycler_assembly(autocycler_fasta, out_dir, genome_size,
                                                          args.min_size_ratio, args.max_size_ratio)
        create_final_assembly(out_dir, flye_fasta, flye_successful, autocycler_fasta,
                              autocycler_successful)

    finally:
        clean_up(out_dir, args.clean)


def check_args(args):
    if args.threads < 1:
        quit_with_error('--threads must be at least 1')
    if args.threads > 100:
        quit_with_error('--threads cannot be greater than 100')
    if args.jobs < 1:
        quit_with_error('--jobs must be at least 1')
    if args.jobs > args.threads:
        quit_with_error('--jobs cannot be greater than --threads')
    if args.subset_count < 2:
        quit_with_error('--subset_count must be at least 2')
    if args.seed < 0:
        quit_with_error('--seed must be at least 0')
    if args.genome_size is not None and args.genome_size < 1:
        quit_with_error('--genome_size must be at least 1')
    if args.min_size_ratio <= 0:
        quit_with_error('--min-size-ratio must be greater than 0')
    if args.max_size_ratio < args.min_size_ratio:
        quit_with_error('--max-size-ratio must be greater than or equal to --min-size-ratio')
    if Path(args.out_dir).exists():
        quit_with_error(f"output directory already exists: '{args.out_dir}'")


def check_dependencies(skip_lrge=False):
    dependencies = (dependency for dependency in DEPENDENCIES
                    if not skip_lrge or dependency != 'lrge')
    missing = [dependency for dependency in dependencies if shutil.which(dependency) is None]
    if missing:
        missing_list = '\n'.join(f'  {dependency}' for dependency in missing)
        quit_with_error(f'the following required dependencies were not found in PATH:\n'
                        f'{missing_list}')


def create_output_dir(out_dir):
    out_dir = Path(out_dir).resolve()
    try:
        out_dir.mkdir(parents=True)
    except OSError as error:
        quit_with_error(f"could not create output directory '{out_dir}': {error}")
    return out_dir


def configure_logging(out_dir):
    (out_dir / 'logs').mkdir()
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(message)s'))
    log_file = out_dir / 'assembly.log'
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setFormatter(logging.Formatter('%(asctime)s  %(message)s',
                                                datefmt='%Y-%m-%d %H:%M:%S'))
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)


def estimate_genome_size_with_lrge(reads, out_dir, threads, seed, read_type):
    logger.info('Genome size estimation started')
    platform = 'pb' if read_type in ('pacbio_clr', 'pacbio_hifi') else 'ont'
    command = ['lrge', '--threads', threads, '--seed', seed, '--platform', platform, reads]
    lrge_log = out_dir / 'logs' / 'lrge.log'
    result = run_command(command, capture_stdout=True, check=False, output_log=lrge_log)
    if result.returncode != 0:
        quit_with_error(f"LRGE failed (see '{lrge_log}')")
    try:
        genome_size = int(result.stdout.strip())
    except (AttributeError, ValueError):
        quit_with_error(f"could not parse LRGE genome size (see '{lrge_log}')")
    if genome_size < 1:
        quit_with_error(f'LRGE returned an invalid genome size: {genome_size}')
    logger.info(f'Estimated genome size: {genome_size:,} bp')
    return genome_size


def subsample_with_rasusa(reads, out_dir, genome_size, read_bases, seed):
    read_depth = read_bases / genome_size
    if read_bases <= genome_size * 100:
        logger.info(f'Rasusa skipped: input read depth {read_depth:.1f}× is at or below 100×')
        return Path(reads)
    logger.info('Rasusa subsampling started')
    rasusa_reads = out_dir / 'rasusa_reads.fastq.gz'
    command = ['rasusa', 'reads', '--coverage', 100, '--genome-size', genome_size,
               '--seed', seed, '--output', rasusa_reads, reads]
    rasusa_log = out_dir / 'logs' / 'rasusa.log'
    result = run_command(command, check=False, output_log=rasusa_log)
    if result.returncode != 0:
        quit_with_error(f"Rasusa failed (see '{rasusa_log}')")
    if not rasusa_reads.is_file() or rasusa_reads.stat().st_size == 0:
        quit_with_error(f"Rasusa did not create a non-empty output file (see '{rasusa_log}')")
    subsampled_read_stats = get_fastq_stats(rasusa_reads)
    logger.info(format_read_stats('Rasusa subsampling successful', subsampled_read_stats))
    return rasusa_reads


def flye_assembly(reads, out_dir, threads, read_type):
    logger.info('Flye assembly started')
    flye_read_flags = {'ont_r9': '--nano-raw', 'ont_r10': '--nano-hq',
                       'pacbio_clr': '--pacbio-raw', 'pacbio_hifi': '--pacbio-hifi'}
    try:
        read_flag = flye_read_flags[read_type]
    except KeyError:
        quit_with_error(f"unsupported read type for Flye: '{read_type}'")
    flye_out_dir = out_dir / 'flye'
    flye_log = out_dir / 'logs' / 'flye.log'
    command = ['flye', read_flag, reads, '--threads', threads, '--out-dir', flye_out_dir]
    result = run_command(command, check=False, output_log=flye_log)
    copy_flye_log(flye_out_dir / 'flye.log', flye_log)
    if result.returncode != 0:
        quit_with_error(f"Flye assembly failed (see '{flye_log}')")
    flye_fasta = flye_out_dir / 'assembly.fasta'
    if not flye_fasta.is_file() or flye_fasta.stat().st_size == 0:
        quit_with_error(f"Flye did not create a non-empty assembly (see '{flye_log}')")
    return flye_fasta


def check_flye_assembly(flye_fasta, genome_size, min_size_ratio, max_size_ratio):
    assembly_size, size_ratio = get_assembly_size_and_ratio(flye_fasta, genome_size)
    size_is_acceptable = min_size_ratio <= size_ratio <= max_size_ratio
    message = f'{assembly_size:,} bp ({size_ratio:.2f}× estimated genome size)'
    if size_is_acceptable:
        logger.info(f'Flye assembly successful: {message}')
    else:
        logger.warning(f'Warning: Flye assembly unsuccessful: {message}')
    return size_is_acceptable


def copy_flye_log(source, destination):
    try:
        if not source.is_file() or source.stat().st_size == 0:
            return
        shutil.copy2(source, destination)
    except OSError as error:
        logger.warning(f"Warning: could not copy Flye log '{source}': {error}")


def autocycler_assembly(reads, flye_fasta, out_dir, genome_size, threads, jobs, subset_count,
                        max_job_time, seed, read_type):
    logger.info('Autocycler assembly started')
    assembly_job_count = len(AUTOCYCLER_ASSEMBLERS) * subset_count
    parallel_jobs = min(jobs, assembly_job_count)
    threads_per_job = threads // parallel_jobs
    autocycler_dir = out_dir / 'autocycler'
    subsampled_reads_dir = autocycler_dir / 'subsampled_reads'
    assemblies_dir = autocycler_dir / 'assemblies'
    autocycler_log = out_dir / 'logs' / 'autocycler.log'
    try:
        autocycler_log.touch()
    except OSError as error:
        logger.warning(f"Warning: could not create Autocycler log '{autocycler_log}': {error}")
        return None

    subsample_command = [
        'autocycler', 'subsample',
        '--reads', reads,
        '--out_dir', subsampled_reads_dir,
        '--genome_size', genome_size,
        '--count', subset_count,
        '--seed', seed,
    ]
    if run_autocycler_stage('subsample', subsample_command, autocycler_log) is None:
        return None

    try:
        assemblies_dir.mkdir()
        flye_100x_fasta = assemblies_dir / 'flye_100x.fasta'
        shutil.copy2(flye_fasta, flye_100x_fasta)
    except OSError as error:
        logger.warning(f'Warning: could not add the 100× Flye assembly to Autocycler: {error}')
        return None

    jobs_file = assemblies_dir / 'jobs.txt'
    try:
        write_autocycler_jobs(jobs_file, subsampled_reads_dir, assemblies_dir, genome_size,
                              subset_count, threads_per_job, read_type)
    except OSError as error:
        logger.warning(f"Warning: could not write Autocycler job file '{jobs_file}': {error}")
        return None
    joblog = assemblies_dir / 'joblog.tsv'
    parallel_command = [
        'nice', '-n', 19,
        'parallel',
        '--jobs', parallel_jobs,
        '--joblog', joblog,
        '--results', assemblies_dir / 'logs',
        '--timeout', max_job_time,
        '::::', jobs_file,
    ]
    parallel_result = run_command(parallel_command, check=False)
    failed_jobs = warn_about_failed_assembly_jobs(joblog, assembly_job_count)
    if parallel_result.returncode != 0 and failed_jobs == 0:
        logger.warning('Warning: one or more Autocycler input assembly jobs failed or timed '
                       f'out (GNU Parallel exit code {parallel_result.returncode})')

    weight_autocycler_input_assemblies(assemblies_dir)
    input_assemblies = [fasta for fasta in assemblies_dir.glob('*.fasta')
                        if fasta.is_file() and fasta.stat().st_size > 0]
    expected_input_assemblies = assembly_job_count + 1
    logger.info(f'Autocycler input assemblies available: {len(input_assemblies)} of '
                f'{expected_input_assemblies}')

    compress_command = [
        'autocycler', 'compress',
        '--assemblies_dir', assemblies_dir,
        '--autocycler_dir', autocycler_dir,
        '--threads', threads,
    ]
    if run_autocycler_stage('compress', compress_command, autocycler_log) is None:
        return None

    cluster_command = [
        'autocycler', 'cluster',
        '--autocycler_dir', autocycler_dir,
    ]
    if run_autocycler_stage('cluster', cluster_command, autocycler_log) is None:
        return None

    cluster_dirs = sorted(path for path in
                          (autocycler_dir / 'clustering' / 'qc_pass').glob('cluster_*')
                          if path.is_dir())
    if not cluster_dirs:
        logger.warning('Warning: Autocycler produced no QC-pass clusters')
        return None

    for cluster_dir in cluster_dirs:
        trim_command = [
            'autocycler', 'trim',
            '--cluster_dir', cluster_dir,
            '--threads', threads,
        ]
        if run_autocycler_stage(f'trim ({cluster_dir.name})', trim_command,
                                autocycler_log) is None:
            return None

        resolve_command = [
            'autocycler', 'resolve',
            '--cluster_dir', cluster_dir,
        ]
        if run_autocycler_stage(f'resolve ({cluster_dir.name})', resolve_command,
                                autocycler_log) is None:
            return None

    final_gfas = [cluster_dir / '5_final.gfa' for cluster_dir in cluster_dirs]
    missing_gfas = [gfa for gfa in final_gfas if not gfa.is_file() or gfa.stat().st_size == 0]
    if missing_gfas:
        logger.warning('Warning: Autocycler did not produce all expected final cluster graphs')
        return None

    combine_command = [
        'autocycler', 'combine',
        '--autocycler_dir', autocycler_dir,
        '--in_gfas', *final_gfas,
    ]
    if run_autocycler_stage('combine', combine_command, autocycler_log) is None:
        return None

    autocycler_fasta = autocycler_dir / 'consensus_assembly.fasta'
    if not autocycler_fasta.is_file() or autocycler_fasta.stat().st_size == 0:
        logger.warning('Warning: Autocycler did not create a non-empty consensus assembly')
        return None

    return autocycler_fasta


def check_autocycler_assembly(autocycler_fasta, out_dir, genome_size, min_size_ratio,
                              max_size_ratio):
    if autocycler_fasta is None:
        logger.warning('Warning: Autocycler assembly unsuccessful: no assembly produced')
        return False
    assembly_size, size_ratio = get_assembly_size_and_ratio(autocycler_fasta, genome_size)
    size_is_acceptable = min_size_ratio <= size_ratio <= max_size_ratio
    fully_resolved = autocycler_completed_successfully(out_dir / 'logs' / 'autocycler.log')
    message = f'{assembly_size:,} bp ({size_ratio:.2f}× estimated genome size)'
    resolution = 'fully resolved' if fully_resolved else 'not fully resolved'
    if size_is_acceptable and fully_resolved:
        logger.info(f'Autocycler assembly successful: {message}; {resolution}')
    else:
        logger.warning(f'Warning: Autocycler assembly unsuccessful: {message}; {resolution}')
    return size_is_acceptable and fully_resolved


def autocycler_completed_successfully(autocycler_log):
    autocycler_success_message = 'Consensus assembly is fully resolved 😄'
    try:
        return autocycler_success_message in autocycler_log.read_text(encoding='utf-8')
    except (OSError, UnicodeError) as error:
        logger.warning(f"Warning: could not read Autocycler log '{autocycler_log}': {error}")
        return False


def create_autocycler_metrics(reads, out_dir):
    header_result = run_command(['autocycler', 'table'], capture_stdout=True, check=False)
    row_result = run_command([
        'autocycler', 'table',
        '--autocycler_dir', out_dir / 'autocycler',
        '--name', reads,
    ], capture_stdout=True, check=False)

    if header_result.returncode != 0 or row_result.returncode != 0:
        logger.warning('Warning: could not create Autocycler metrics table')
        return False

    header = (header_result.stdout or '').rstrip('\r\n')
    row = (row_result.stdout or '').rstrip('\r\n')
    if not header or not row:
        logger.warning('Warning: could not create Autocycler metrics table: empty table output')
        return False

    metrics_tsv = out_dir / 'autocycler' / 'metrics.tsv'
    try:
        metrics_tsv.write_text(f'{header}\n{row}\n', encoding='utf-8')
    except OSError as error:
        logger.warning(f"Warning: could not write Autocycler metrics table '{metrics_tsv}': "
                       f'{error}')
        return False

    return True


def create_final_assembly(out_dir, flye_fasta, flye_successful, autocycler_fasta,
                          autocycler_successful):
    if autocycler_successful:
        assembler = 'Autocycler'
        source_fasta = autocycler_fasta
        source_gfa = out_dir / 'autocycler' / 'consensus_assembly.gfa'
    elif flye_successful:
        assembler = 'Flye'
        source_fasta = flye_fasta
        source_gfa = out_dir / 'flye' / 'assembly_graph.gfa'
    else:
        quit_with_error('neither Autocycler nor Flye produced a successful assembly')

    source_files = ((source_fasta, 'FASTA'), (source_gfa, 'GFA'))
    for source, file_type in source_files:
        if source is None or not source.is_file() or source.stat().st_size == 0:
            quit_with_error(f"{assembler} {file_type} is missing or empty: '{source}'")

    try:
        shutil.copy2(source_fasta, out_dir / 'assembly.fasta')
        shutil.copy2(source_gfa, out_dir / 'assembly.gfa')
    except OSError as error:
        quit_with_error(f'could not copy final {assembler} assembly: {error}')

    logger.info(f'Final assembly: {assembler}')


def write_autocycler_jobs(jobs_file, subsampled_reads_dir, assemblies_dir, genome_size,
                          subset_count, threads_per_job, read_type):
    with jobs_file.open('w', encoding='utf-8') as jobs:
        for assembler in AUTOCYCLER_ASSEMBLERS:
            for sample_number in range(1, subset_count + 1):
                sample = subsampled_reads_dir / f'sample_{sample_number:02d}.fastq'
                out_prefix = assemblies_dir / f'{assembler}_{sample_number:02d}'
                command = [
                    'autocycler', 'helper', assembler,
                    '--reads', sample,
                    '--out_prefix', out_prefix,
                    '--threads', threads_per_job,
                    '--genome_size', genome_size,
                    '--read_type', read_type,
                    '--min_depth_rel', 0.1,
                ]
                jobs.write(f'{shlex.join(str(part) for part in command)}\n')


def warn_about_failed_assembly_jobs(joblog, assembly_job_count):
    if not joblog.is_file():
        return 0

    failed_jobs = []
    try:
        with joblog.open(encoding='utf-8', newline='') as joblog_file:
            for row in csv.DictReader(joblog_file, delimiter='\t'):
                if row.get('Exitval') != '0' or row.get('Signal') != '0':
                    failed_jobs.append(row)
    except (OSError, csv.Error) as error:
        logger.warning(f"Warning: could not read GNU Parallel job log '{joblog}': {error}")
        return 0

    if failed_jobs:
        logger.warning(f'Warning: {len(failed_jobs)} of {assembly_job_count} Autocycler input '
                       'assembly jobs failed or timed out (see autocycler/assemblies/joblog.tsv and '
                       'autocycler/assemblies/logs/)')
    return len(failed_jobs)


def weight_autocycler_input_assemblies(assemblies_dir):
    for fasta in assemblies_dir.glob('plassembler*.fasta'):
        add_fasta_header_attribute(fasta, 'Autocycler_cluster_weight=3',
                                   required_text='circular=True')
    for fasta in assemblies_dir.glob('flye*.fasta'):
        add_fasta_header_attribute(fasta, 'Autocycler_consensus_weight=2')


def add_fasta_header_attribute(fasta, attribute, required_text=None):
    try:
        text = fasta.read_text(encoding='utf-8')
        if not text:
            return
        lines = text.splitlines()
        changed = False
        for i, line in enumerate(lines):
            if (line.startswith('>') and attribute not in line
                    and (required_text is None or required_text in line)):
                lines[i] = f'{line} {attribute}'
                changed = True
        if changed:
            fasta.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    except (OSError, UnicodeError) as error:
        logger.warning(f"Warning: could not update FASTA headers in '{fasta}': {error}")


def run_autocycler_stage(stage, command, autocycler_log):
    result = run_command(command, check=False, output_log=autocycler_log, append_log=True)
    if result.returncode != 0:
        logger.warning(f'Warning: Autocycler {stage} failed with exit code '
                       f'{result.returncode}; the pipeline will fall back to Flye')
        return None
    return result


def get_fastq_stats(reads):
    reads = Path(reads)
    read_lengths = []
    try:
        if reads.suffix == '.gz':
            fastq = gzip.open(reads, 'rt', encoding='utf-8')
        else:
            fastq = reads.open(encoding='utf-8')
        with fastq:
            while True:
                header = fastq.readline()
                if not header:
                    break
                record_number = len(read_lengths) + 1
                if not header.startswith('@'):
                    quit_with_error(f"malformed FASTQ '{reads}': record {record_number} does "
                                    "not start with '@'")
                sequence_length = 0
                while True:
                    line = fastq.readline()
                    if not line:
                        quit_with_error(f"truncated FASTQ '{reads}' in record {record_number}")
                    if line.startswith('+'):
                        break
                    sequence_length += len(line.rstrip('\r\n'))
                if sequence_length == 0:
                    quit_with_error(f"malformed FASTQ '{reads}': record {record_number} has "
                                    'an empty sequence')
                quality_length = 0
                while quality_length < sequence_length:
                    line = fastq.readline()
                    if not line:
                        quit_with_error(f"truncated FASTQ '{reads}' in record {record_number}")
                    quality_length += len(line.rstrip('\r\n'))
                if quality_length != sequence_length:
                    quit_with_error(f"malformed FASTQ '{reads}': sequence and quality lengths "
                                    f'differ in record {record_number}')
                read_lengths.append(sequence_length)
    except (OSError, EOFError, UnicodeError) as error:
        quit_with_error(f"could not read FASTQ '{reads}': {error}")

    if not read_lengths:
        quit_with_error(f"FASTQ contains no reads: '{reads}'")

    total_bases = sum(read_lengths)
    halfway = (total_bases + 1) // 2
    cumulative_bases = 0
    read_n50 = 0
    for read_length in sorted(read_lengths, reverse=True):
        cumulative_bases += read_length
        if cumulative_bases >= halfway:
            read_n50 = read_length
            break

    return len(read_lengths), total_bases, read_n50


def format_read_stats(label, read_stats):
    read_count, total_bases, read_n50 = read_stats
    return f'{label}: {read_count:,} reads; {total_bases:,} bp; N50 {read_n50:,} bp'


def get_assembly_size_and_ratio(fasta, genome_size):
    assembly_size = get_fasta_size(fasta)
    return assembly_size, assembly_size / genome_size


def get_fasta_size(fasta):
    try:
        with fasta.open(encoding='utf-8') as fasta_file:
            size = sum(len(line.strip()) for line in fasta_file if not line.startswith('>'))
    except (OSError, UnicodeError) as error:
        quit_with_error(f"could not read assembly FASTA '{fasta}': {error}")
    if size < 1:
        quit_with_error(f"assembly FASTA contains no sequence: '{fasta}'")
    return size


def clean_up(out_dir, clean_level):
    if clean_level >= 1:
        clean_up_reads(out_dir)
    if clean_level >= 2:
        clean_up_working_directories(out_dir)
    if clean_level >= 3:
        clean_up_assembler_directories(out_dir)


def clean_up_reads(out_dir):
    rasusa_reads = out_dir / 'rasusa_reads.fastq.gz'
    remove_path(rasusa_reads)
    subsampled_reads_dir = out_dir / 'autocycler' / 'subsampled_reads'
    for pattern in ('*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz'):
        for reads in subsampled_reads_dir.glob(pattern):
            remove_path(reads)


def clean_up_working_directories(out_dir):
    flye_dir = out_dir / 'flye'
    for directory in ['00-assembly', '10-consensus', '20-repeat', '30-contigger', '40-polishing']:
        remove_path(flye_dir / directory)
    autocycler_dir = out_dir / 'autocycler'
    remove_path(autocycler_dir / 'subsampled_reads')
    remove_path(autocycler_dir / 'assemblies')
    remove_path(autocycler_dir / 'clustering')


def clean_up_assembler_directories(out_dir):
    remove_path(out_dir / 'autocycler')
    remove_path(out_dir / 'flye')


def remove_path(path):
    try:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        elif path.exists() or path.is_symlink():
            path.unlink()
    except OSError as error:
        logger.warning(f"Warning: could not remove '{path}': {error}")


def run_command(command, capture_stdout=False, check=True, output_log=None, append_log=False):
    command = [str(part) for part in command]

    try:
        if capture_stdout:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True, encoding='utf-8', errors='replace')
            if output_log is not None:
                write_command_output(output_log, result, append_log)
        else:
            if output_log is None:
                result = subprocess.run(command, stdout=subprocess.DEVNULL,
                                        stderr=subprocess.STDOUT)
            else:
                mode = 'a' if append_log else 'w'
                with Path(output_log).open(mode, encoding='utf-8') as log_file:
                    result = subprocess.run(command, stdout=log_file, stderr=subprocess.STDOUT,
                                            text=True, encoding='utf-8', errors='replace')
    except OSError as error:
        if check:
            quit_with_error(f"could not run command '{command[0]}': {error}")
        if output_log is not None:
            try:
                mode = 'a' if append_log else 'w'
                with Path(output_log).open(mode, encoding='utf-8') as log_file:
                    log_file.write(f"Could not run command '{command[0]}': {error}\n")
            except OSError:
                pass
        logger.warning(f"Warning: could not run command '{command[0]}': {error}")
        result = subprocess.CompletedProcess(command, 127, stderr=str(error))

    if check and result.returncode != 0:
        quit_with_error(f"command failed with exit code {result.returncode}: "
                        f"'{shlex.join(command)}'")
    return result


def write_command_output(output_log, result, append):
    mode = 'a' if append else 'w'
    with Path(output_log).open(mode, encoding='utf-8') as log_file:
        log_file.write(result.stderr or '')
        log_file.write(result.stdout or '')


def quit_with_error(message):
    if logger.handlers:
        logger.error(f'Error: {message}')
    else:
        print(f'Error: {message}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
    main()
