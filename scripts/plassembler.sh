#!/usr/bin/env bash

# This script is a wrapper for running Plassembler in a single command.

# Usage:
#   plassembler.sh <read_fastq> <assembly_prefix> <threads>

# Requirements:
#   Plassembler: https://github.com/gbouras13/plassembler
#   Plassembler database (in $CONDA_PREFIX/plassembler_db or location set with PLASSEMBLER_DB)
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
for cmd in plassembler; do
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

# Find the Plassembler database.
if [[ -n "$PLASSEMBLER_DB" && -d "$PLASSEMBLER_DB" ]]; then
    database="$PLASSEMBLER_DB"
elif [[ -n "$CONDA_PREFIX" && -d "$CONDA_PREFIX/plassembler_db" ]]; then
    database="$CONDA_PREFIX/plassembler_db"
else
    >&2 echo "Error: No Plassembler database found."
    >&2 echo "       Please set PLASSEMBLER_DB or ensure \$CONDA_PREFIX/plassembler_db exists."
    exit 1
fi

# Create a temporary directory which will be deleted when the script exits.
temp_dir=$(mktemp -d)
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT

# Ensure the reads are gzipped.
if file "$reads" | grep -q 'gzip compressed data'; then
    gzipped_reads="$reads"
else
    gzip -c "$reads" > "$temp_dir/reads.fastq.gz"
    gzipped_reads="$temp_dir/reads.fastq.gz"
fi

# Run Plassembler.
plassembler long -d "$database" -l "$gzipped_reads" -o "$temp_dir"/out -t "$threads" --skip_qc

# Check if Plassembler ran successfully.
if [[ ! -s "$temp_dir"/out/plassembler_plasmids.fasta ]]; then
    >&2 echo "Error: Plassembler assembly failed."
    exit 1
fi

# Copy the GFA file as is.
cp "$temp_dir"/out/plassembler_plasmids.gfa "$assembly".gfa

# For the FASTA file, rotate any circular sequences by a random amount.
seqtk seq "$temp_dir"/out/plassembler_plasmids.fasta | paste - - | awk 'BEGIN {FS="\t"; srand(); }
{
  if ($1 ~ /circular=True/) {
    r = int(rand() * (length($2) - 1)) + 1;
    $2 = substr($2, r+1) substr($2, 1, r);
  }
  print $1; print $2;
}' > "$assembly".fasta
