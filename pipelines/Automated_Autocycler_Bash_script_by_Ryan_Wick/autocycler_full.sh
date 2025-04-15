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
    echo "Usage: $0 <read_fastq> <threads> <jobs>" 1>&2
    exit 1
fi
if [[ ! -f "$reads" ]]; then
    echo "Error: Input file '$reads' does not exist." 1>&2
    exit 1
fi
if (( threads > 128 )); then threads=128; fi  # Flye won't work with more than 128 threads

genome_size=$(genome_size_raven.sh "$reads" "$threads")

# Step 1: subsample the long-read set into multiple files
autocycler subsample --reads "$reads" --out_dir subsampled_reads --genome_size "$genome_size" 2>> autocycler.stderr

# Step 2: assemble each subsampled file
mkdir -p assemblies
rm -f assemblies/jobs.txt
for assembler in raven miniasm flye metamdbg necat nextdenovo plassembler canu; do
    for i in 01 02 03 04; do
        echo "$assembler.sh subsampled_reads/sample_$i.fastq assemblies/${assembler}_$i $threads $genome_size" >> assemblies/jobs.txt
    done
done
set +e
nice -n 19 parallel --jobs "$jobs" --joblog assemblies/joblog.txt --results assemblies/logs < assemblies/jobs.txt
set -e
find assemblies/ -maxdepth 1 -type f -name "*.fasta" -empty -delete

# Give circular contigs from Plassembler extra clustering weight
for f in assemblies/plassembler*.fasta; do
    sed -i 's/circular=True/circular=True Autocycler_cluster_weight=2/' "$f"
done

# Give contigs from Canu and Flye extra consensus weight
for f in assemblies/canu*.fasta assemblies/flye*.fasta; do
    sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
done

# Remove the subsampled reads to save space
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
