# Automated Autocycler Bash script (by Ryan Wick)

This is a simple Bash script that automates running a full Autocycler assembly workflow. The script is designed to be minimalistic, without many frills and assumes that the input reads are ready for assembly.

I wrote this script in January 2025 for Autocycler v0.2.1. I'll do my best to update it if future Autocycler versions break compatibility, but I can't make any promises.



## Key Features

* Does not perform any quality control on the input reads before starting the assembly.
* Uses Raven to [estimate the genome size](https://github.com/rrwick/Autocycler/wiki/Genome-size-estimation).
* Hard-codes the use of four read subsets for the assembly process.
* Uses [GNU Parallel](https://github.com/rrwick/Autocycler/wiki/Parallelising-input-assemblies#gnu-parallel) to run multiple assembly jobs at once.
* Runs assemblies with `nice -n 19` to give them lower priority with the operating system.
* Uses seven different assemblers (in this order): Raven, miniasm, Flye, MetaMDBG, NECAT, NextDenovo and Canu. This order was chosen to put the faster assemblers first and slower assemblers last, in case you want to check a run in progress to see if the assemblies look okay.
* 4 read subsets × 7 assemblers = 28 total input assemblies (assume all are successful).



## Usage

The script takes the following three arguments:
1. **Read filename**: Path to the input FASTQ file (can be gzipped).
2. **Thread count**: Number of threads per assembly.
3. **Job count**: Number of simultaneous assemblies to run. Note: Each assembly will use up to the given thread count, so the maximum total threads in use will be threads × jobs.


**Example command:** `./autocycler_full.sh reads.fastq.gz 16 4`



## Output

When run, the script will create the following directories and files in the working directory:

* **`subsampled_reads/`**: Directory containing the read subsets for assembly. The actual FASTQ files will be deleted after the assemblies are complete (to save disk space), but the directory and its YAML file will remain.
* **`assemblies/`**: Directory containing the input assemblies for Autocycler. The logs for each assembler can be found in the `assemblies/logs` directory.
* **`autocycler_out/`**: Autocycler output directory which will include the final combined assembly as `consensus_assembly.gfa` and `consensus_assembly.fasta`.
* **`autocycler.stderr`**: File containing all `stderr` output from Autocycler across all steps.
