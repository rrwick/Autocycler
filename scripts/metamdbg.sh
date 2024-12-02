#!/usr/bin/env bash

# This script is a wrapper for running metaMDBG in a single command.

# Usage:
#   metamdbg.sh <read_fastq> <assembly_prefix> <threads>

# Requirements:
#   metaMDBG: https://github.com/GaetanBenoitDev/metaMDBG

# Copyright 2024 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Autocycler

# This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
# is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details. You should have received a copy of the GNU General Public
# License along with Autocycler. If not, see <https://www.gnu.org/licenses/>.


# Ensure script exits on error.
set -e

# The minimum contig depth is hard-coded here.
min_cov=10

# Get arguments.
reads=$1        # input reads FASTQ
assembly=$2     # output assembly prefix (not including file extension)
threads=$3      # thread count

# Validate input parameters.
if [[ -z "$reads" || -z "$assembly" || -z "$threads" ]]; then
    >&2 echo "Usage: $0 <read_fastq> <assembly_prefix> <threads>"
    exit 1
fi

# Check that the reads file exists.
if [[ ! -f "$reads" ]]; then
    >&2 echo "Error: $reads does not exist"
    exit 1
fi

# Ensure the requirements are met.
for cmd in metaMDBG; do
    if ! command -v "$cmd" &> /dev/null; then
        >&2 echo "Error: $cmd not found in PATH"
        exit 1
    fi
done

# Ensure the output prefix will work.
if ! touch "$assembly".fasta &> /dev/null; then
    >&2 echo "Error: cannot write to this location: $assembly"
    exit 1
fi

# Create a temporary directory which is deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

# Run metaMDBG.
metaMDBG asm --out-dir "$temp_dir" --in-ont "$reads" --threads "$threads"

# Check if metaMDBG ran successfully.
if [[ ! -s "$temp_dir"/contigs.fasta.gz ]]; then
    >&2 echo "Error: metaMDBG assembly failed."
    exit 1
fi

# Remove low-depth contigs.
metamdbg_filter.py "$temp_dir"/contigs.fasta.gz > "$assembly".fasta
