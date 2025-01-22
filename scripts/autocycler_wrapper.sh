#!/usr/bin/env bash

# Command line wrapper for autocycler.
# Usage:
#   autocycler_wrapper.sh <reads> <output_folder> <read_partitions> <threads> 
#
# Example usage:
#   autocycler_wrapper.sh ont.fastq.gz output_autocycler 4 16
#
# The example  runs  autocycler on the reads 'ont.fastq.gz' and
# outputs everything into the folder 'output_autocycler'.
# It divides the reads into 4 partitions and uses 16 threads.

reads=$1  # Your read set goes here
output=$2  # Name you output folder here
subset=$3  # Number of read subset partitions (default is 4. Auto-prefixes 0 for any number below 10)
threads=$4  # set as appropriate for your system

# Add 0 prefix for read sub-sets if lower than 10
if [ "${#1}" -lt "2" ]; then
    subset="0${subset}"
fi

# Input parameter print
echo Running autocycler with following parameters\:
echo Reads\: $reads
echo Output folder\: $output
echo Read partitions\: $subset
echo Threads\: $threads
mkdir -p $output

# Estimate genome size: Use lrge if it is installed, otherwise use the slower bundled Raven estimator
echo -e "\nEstimating genome size:"
if ! command -v lrge 2>&1 >/dev/null
then
    genome_size=$(genome_size_raven.sh "$reads" "$threads")
else
    genome_size=$(lrge -t "$threads" "$reads")
fi

# Step 1: subsample the long-read set into multiple files
autocycler subsample --count $subset  --reads "$reads" --out_dir ${output}/subsampled_reads --genome_size "$genome_size"

# Step 2: assemble each subsampled file
mkdir -p $output/assemblies
for assembler in canu flye miniasm necat nextdenovo raven; do
    for i in `eval echo {01..$subset}`; do
        "$assembler".sh ${output}/subsampled_reads/sample_"$i".fastq ${output}/assemblies/"$assembler"_"$i" "$threads" "$genome_size"
    done
done

# Optional step: remove the subsampled reads to save space
rm ${output}/subsampled_reads/*.fastq

# Step 3: compress the input assemblies into a unitig graph
autocycler compress -i ${output}/assemblies -a ${output}/autocycler_out

# Step 4: cluster the input contigs into putative genomic sequences
autocycler cluster -a ${output}/autocycler_out

# Steps 5 and 6: trim and resolve each QC-pass cluster
for c in ${output}/autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c"
    autocycler resolve -c "$c"
done

# Step 7: combine resolved clusters into a final assembly
echo Running Autocycler combine
autocycler combine -a ${output}/autocycler_out -i ${output}/autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

