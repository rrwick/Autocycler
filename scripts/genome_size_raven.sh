#!/usr/bin/env bash

# This script is a wrapper for getting a genome size estimate in a single command using Raven.

# Usage:
#   raven.sh <read_fastq> <threads>

# Requirements:
#   Raven: https://github.com/lbcb-sci/raven
#   seqtk: https://github.com/lh3/seqtk

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
threads=$2      # thread count

# Validate input parameters.
if [[ -z "$reads" || -z "$threads" ]]; then
    >&2 echo "Usage: $0 <read_fastq> <threads>"
    exit 1
fi

# Check that the reads file exists.
if [[ ! -f "$reads" ]]; then
    >&2 echo "Error: $reads does not exist"
    exit 1
fi

# Ensure the requirements are met.
for cmd in raven seqtk cut; do
    if ! command -v "$cmd" &> /dev/null; then
        >&2 echo "Error: $cmd not found in PATH"
        exit 1
    fi
done

# Create a temporary directory which is deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

# Run Raven.
raven --threads "$threads" --disable-checkpoints "$reads" > "$temp_dir"/raven.fasta

# Check if Raven ran successfully.
if [[ ! -s "$temp_dir"/raven.fasta ]]; then
    >&2 echo "Error: Raven assembly failed."
    exit 1
fi

# Print genome size.
seqtk size "$temp_dir"/raven.fasta | cut -f2
