#!/usr/bin/env bash

# This script is a wrapper for running miniasm and Minipolish in a single command.

# Usage:
#   miniasm.sh <read_fastq> <assembly_prefix> <threads>

# Requirements:
#   miniasm: https://github.com/lh3/miniasm
#   Minipolish: https://github.com/rrwick/Minipolish
#   minimap2: https://github.com/lh3/minimap2
#   Racon: https://github.com/lbcb-sci/racon
#   any2fasta: https://github.com/tseemann/any2fasta

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
for cmd in miniasm minipolish minimap2 racon any2fasta; do
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

# Find read overlaps with minimap2.
minimap2 -x ava-ont -t "$threads" "$reads" "$reads" > "$temp_dir"/overlap.paf

# Run miniasm to make an unpolished assembly.
miniasm -f "$reads" "$temp_dir"/overlap.paf > "$temp_dir"/unpolished.gfa

# Check if miniasm ran successfully.
if [[ ! -s "$temp_dir"/unpolished.gfa ]]; then
    >&2 echo "Error: miniasm assembly failed."
    exit 1
fi

# Polish the assembly with Minipolish, outputting the result to stdout.
minipolish --threads "$threads" "$reads" "$temp_dir"/unpolished.gfa > "$assembly".gfa

# Check if Minipolish ran successfully.
if [[ ! -s "$assembly".gfa ]]; then
    >&2 echo "Error: Minipolish assembly failed."
    exit 1
fi

# Convert the GFA assembly to FASTA format.
any2fasta "$assembly".gfa > "$assembly".fasta

# Check if any2fasta ran successfully.
if [[ ! -s "$assembly".fasta ]]; then
    >&2 echo "Error: any2fasta assembly failed."
    exit 1
fi
