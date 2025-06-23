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
reads=$1                 # input reads FASTQ
threads=$2               # threads per job
jobs=$3                  # number of simultaneous jobs
read_type=${4:-ont_r10}  # read type (default = ont_r10)

# Validate input parameters.
if [[ -z "$reads" || -z "$threads" || -z "$jobs" ]]; then
    echo "Usage: $0 <read_fastq> <threads> <jobs> [read_type]" 1>&2
    exit 1
fi
if [[ ! -f "$reads" ]]; then
    echo "Error: Input file '$reads' does not exist." 1>&2
    exit 1
fi
if (( threads > 128 )); then threads=128; fi  # Flye won't work with more than 128 threads
case $read_type in
    ont_r9|ont_r10|pacbio_clr|pacbio_hifi) ;;
    *) echo "Error: read_type must be ont_r9, ont_r10, pacbio_clr or pacbio_hifi" 1>&2; exit 1 ;;
esac

genome_size=$(autocycler helper genome_size --reads "$reads" --threads "$threads")

# Step 1: subsample the long-read set into multiple files
autocycler subsample --reads "$reads" --out_dir subsampled_reads --genome_size "$genome_size" 2>> autocycler.stderr

# Step 2: assemble each subsampled file
mkdir -p assemblies
rm -f assemblies/jobs.txt
for assembler in raven miniasm flye metamdbg necat nextdenovo plassembler canu; do
    for i in 01 02 03 04; do
        echo "autocycler helper $assembler --reads subsampled_reads/sample_$i.fastq --out_prefix assemblies/${assembler}_$i --threads $threads --genome_size $genome_size --read_type $read_type" --min_depth_rel 0.1 >> assemblies/jobs.txt
    done
done
set +e
nice -n 19 parallel --jobs "$jobs" --joblog assemblies/joblog.tsv --results assemblies/logs < assemblies/jobs.txt
set -e

# Give circular contigs from Plassembler extra clustering weight
shopt -s nullglob
for f in assemblies/plassembler*.fasta; do
    sed -i 's/circular=True/circular=True Autocycler_cluster_weight=2/' "$f"
done

# Give contigs from Canu and Flye extra consensus weight
for f in assemblies/canu*.fasta assemblies/flye*.fasta; do
    sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
done
shopt -u nullglob

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
