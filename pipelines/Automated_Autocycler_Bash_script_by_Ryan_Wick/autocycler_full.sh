#!/usr/bin/env bash

# This script is a wrapper for running a fully-automated Autocycler assembly.

# Usage:
#   autocycler_full.sh <read_fastq> <threads> <jobs>

# Copyright 2025 Ryan Wick (rrwick@gmail.com)
# Licensed under the GNU General Public License v3.
# See https://www.gnu.org/licenses/gpl-3.0.html.

# Ensure script exits on error.
set -e

# Get arguments.
reads=$1    # input reads FASTQ
threads=$2  # threads per job
jobs=$3     # number of simultaneous jobs

# Validate input parameters.
if [[ -z "$reads" || -z "$threads" || -z "$jobs" ]]; then
    >&2 echo "Usage: $0 <read_fastq> <threads> <jobs>"
    exit 1
fi
if [[ ! -f "$reads" ]]; then
    >&2 echo "Error: Input file '$reads' does not exist."
    exit 1
fi

genome_size=$(genome_size_raven.sh "$reads" "$threads")

# Step 1: subsample the long-read set into multiple files
autocycler subsample --reads "$reads" --out_dir subsampled_reads --genome_size "$genome_size" 2>> autocycler.stderr

# Step 2: assemble each subsampled file
mkdir -p assemblies
rm -f assemblies/jobs.txt
for assembler in raven miniasm flye metamdbg necat nextdenovo canu; do
    for i in 01 02 03 04; do
        echo "nice -n 19 $assembler.sh subsampled_reads/sample_$i.fastq assemblies/${assembler}_$i $threads $genome_size" >> assemblies/jobs.txt
    done
done
set +e
parallel --jobs "$jobs" --joblog assemblies/joblog.txt --results assemblies/logs < assemblies/jobs.txt
set -e
find assemblies/ -maxdepth 1 -type f -name "*.fasta" -empty -delete

# Optional step: remove the subsampled reads to save space
rm subsampled_reads/*.fastq

# Step 3: compress the input assemblies into a unitig graph
autocycler compress -i assemblies -a autocycler_out 2>> autocycler.stderr

# Step 4: cluster the input contigs into putative genomic sequences
autocycler cluster -a autocycler_out 2>> autocycler.stderr

# Steps 5 and 6: trim and resolve each QC-pass cluster
for c in autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c" 2>> autocycler.stderr
    autocycler resolve -c "$c" 2>> autocycler.stderr
done

# Step 7: combine resolved clusters into a final assembly
autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa 2>> autocycler.stderr
