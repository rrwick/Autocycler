# Automated Autocycler Slurm Bash script (by Michael Hall @mbhall88)

This script is heavily inspired (copied) from Ryan's ([`full_autocycler.sh`](../Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh)).

It extends/changes the functionality in the following ways:
- It allows optional read filtering up front with `filtlong`
- It estimates genome size with [`lrge`](https://github.com/mbhall88/lrge). This reduces the memory and time requirements of this step drastically.
- Each assembly job is submitted to the Slurm job scheduler using [`ssubmit`](https://github.com/mbhall88/ssubmit/). In the original script, each assembly job was run using `parallel` in batches based on the `jobs` value provided to the script. This script therefore allows running all assemblies concurrently (depending on job queues on your cluster). The script will wait for the jobs to complete, exiting if any of the jobs failed. This ultimately makes the full autocycler process faster - I am impatient...
- There is a resume-after-assembly functionality. Sometimes I would find an assembly job might fail due to memory (or other) issues and then it is a bit annoying to manually run the remaining steps of the original script. This script has an option `-R` that allows you to manually complete that assembly job outside of this script and then rerun the autocycler process, skipping the filtering, subsampling, and assembly parts.

## Dependencies

All requirements are listed in the [`environment.yaml`](environment.yaml) file. You can create a conda environment with all the required software using:

```bash
conda env create -f environment.yaml
```

## Usage

The only required argument is the input FASTQ file containing long reads for assembly. All other options are optional.

**Basic usage example:**
```bash
./autocycler_slurm.sh -o output_directory -t 16 -l 1000 reads.fastq
```

```
Usage: autocycler_slurm.sh <read_fastq> [options]
This script performs long-read assembly using Autocycler on a Slurm cluster. 
In particular, it submits each assembly job to the cluster scheduler, thus 
enabling parallel execution of multiple assemblers and subsamples.

Required:
  <read_fastq>               Input FASTQ file for long-read assembly

General options:
  -o, --outdir <path>         Output directory [default: current directory]
  -t, --threads <int>         Threads per assembly job [default: 8]
  -r, --read-type <type>      Read type: ont_r9 | ont_r10 | pacbio_clr | pacbio_hifi [default: ont_r10]
  -k, --keep-intermediate     Keep subsampled and filtered FASTQ files after assembly
  -w, --overwrite             Delete output directory if it exists (use with caution!)
  -T, --max-time <duration>   Slurm time limit per job [default: 8h]
  -M, --max-mem <size>        Slurm memory per job [default: 32g]
  -R, --resume-after-assembly Resume pipeline after assembly step (skip read filtering, subsampling, and assembly)

Read filtering with filtlong:
  -l, --min-length <int>      Filter out reads shorter than this length
  -b, --target-bases <int>    Keep top reads until total base count is reached
  -p, --keep-percent <float>  Keep best X%% of reads (e.g. 90)

Help:
  -h, --help                  Show this help message and exit
```

## Output

The script creates the following outputs in the working directory:

* **`subsampled_reads/`**: Contains the read subsets used for assembly. The FASTQ files are deleted after assembly to save space, but the directory (and its [YAML file](https://github.com/rrwick/Autocycler/wiki/Metrics#read-subsampling-metrics)) remains. You can use the `-k` option to keep these files.
* **`filtered_reads/`**: Contains the filtered reads if read filtering was performed. You can use the `-k` option to keep these files, otherwise they are deleted after assembly.
* **`assemblies/`**:
  * Input assemblies for Autocycler in FASTA format
  * GFA and log files for assemblers that produce them
  * `slurm_logs/`: full stdout/stderr for each assembly slurm job
* **`autocycler_out/`**: Output directory for Autocycler. Final results are saved as `consensus_assembly.gfa` and `consensus_assembly.fasta`.
* **`autocycler.stderr`**: stderr output from all Autocycler steps.
