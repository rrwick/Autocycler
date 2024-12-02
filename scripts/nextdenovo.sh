#!/usr/bin/env bash

# This script is a wrapper for running NextDenovo and NextPolish in a single command.

# Usage:
#   nextdenovo.sh <read_fastq> <assembly_prefix> <threads> <genome_size>

# Requirements:
#   NextDenovo: https://github.com/Nextomics/NextDenovo
#   NextPolish: https://github.com/Nextomics/NextPolish

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

# Ensure the requirements are met.
for cmd in nextDenovo nextPolish; do
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
    popd  # return to original directory
    rm -rf "$temp_dir"
}
trap cleanup EXIT

# Work from the temp directory.
pushd "$temp_dir"

# Create the NextDenovo config and read list files.
cat <<EOF > nextdenovo_run.cfg
[General]
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 1
input_type = raw
read_type = ont
input_fofn = input.fofn
workdir = nextdenovo

[correct_option]
read_cutoff = 1k
genome_size = $genome_size
sort_options = -m 20g -t $threads
minimap2_options_raw = -t $threads
pa_correction = 1
correction_options = -p $threads

[assemble_option]
minimap2_options_cns = -t $threads
nextgraph_options = -a 1
EOF
echo "$reads_abs" > input.fofn

# Run NextDenovo.
nextDenovo nextdenovo_run.cfg

# Check if NextDenovo ran successfully.
if [[ ! -s nextdenovo/03.ctg_graph/nd.asm.fasta ]]; then
    >&2 echo "Error: NextDenovo assembly failed."
    exit 1
fi

# Create the NextPolish config and read list files.
cat <<EOF > nextpolish_run.cfg
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
deltmp = yes
rerun = 3
parallel_jobs = 1
multithread_jobs = $threads
genome = nextdenovo/03.ctg_graph/nd.asm.fasta
genome_size = auto
workdir = nextpolish
polish_options = -p $threads

[lgs_option]
lgs_fofn = lgs.fofn
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-ont -t $threads
EOF
echo "$reads_abs" > lgs.fofn

# Run NextPolish.
nextPolish nextpolish_run.cfg

# Check if NextPolish ran successfully.
if [[ ! -s nextpolish/genome.nextpolish.fasta ]]; then
    >&2 echo "Error: NextPolish failed."
    exit 1
fi

# Copy output file.
cp "$temp_dir"/nextpolish/genome.nextpolish.fasta "$assembly_abs".fasta
