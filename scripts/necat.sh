#!/usr/bin/env bash

# This script is a wrapper for running NECAT in a single command.

# Usage:
#   necat.sh <read_fastq> <assembly_prefix> <threads> <genome_size>

# Requirements:
#   NECAT: https://github.com/xiaochuanle/NECAT

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

# Get arguments.
reads=$1        # input reads FASTQ
assembly=$2     # output assembly prefix (not including file extension)
threads=$3      # thread count
genome_size=$4  # estimated genome size

# Validate input parameters.
if [[ -z "$reads" || -z "$assembly" || -z "$threads" || -z "$genome_size" ]]; then
    >&2 echo "Usage: $0 <read_fastq> <assembly_prefix> <threads> <genome_size>"
    exit 1
fi
assembly_abs=$(realpath "$assembly")

# Check that the reads file exists.
if [[ ! -f "$reads" ]]; then
    >&2 echo "Error: $reads does not exist"
    exit 1
fi
reads_abs=$(realpath "$reads")

# Ensure the requirements are met (either necat.pl or necat is okay).
found_necat=false
for cmd in necat.pl necat; do
    if command -v "$cmd" &> /dev/null; then
        necat_executable="$cmd"
        found_necat=true
        break
    fi
done
if [[ "$found_necat" = false ]]; then
    >&2 echo "Error: neither necat.pl nor necat found in PATH"
    exit 1
fi

# Ensure the output prefix will work.
if ! touch "$assembly".fasta &> /dev/null; then
    >&2 echo "Error: cannot write to this location: $assembly"
    exit 1
fi

# Create a temporary directory which is deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    popd  # return to original directory
    rm -rf "$temp_dir"
}
trap cleanup EXIT

# Work from the temp directory.
pushd "$temp_dir"

# Create config and read list files.
"$necat_executable" config config.txt
sed -i "s/PROJECT=/PROJECT=necat/" config.txt
sed -i "s/ONT_READ_LIST=/ONT_READ_LIST=read_list.txt/" config.txt
sed -i "s/GENOME_SIZE=/GENOME_SIZE=$genome_size/" config.txt
sed -i "s/THREADS=4/THREADS=$threads/" config.txt
echo "$reads_abs" > read_list.txt

# Run NECAT.
"$necat_executable" bridge config.txt

# Check if NECAT ran successfully.
if [[ ! -s necat/6-bridge_contigs/polished_contigs.fasta ]]; then
    >&2 echo "Error: NECAT assembly failed."
    exit 1
fi

# Copy output file.
cp "$temp_dir"/necat/6-bridge_contigs/polished_contigs.fasta "$assembly_abs".fasta
