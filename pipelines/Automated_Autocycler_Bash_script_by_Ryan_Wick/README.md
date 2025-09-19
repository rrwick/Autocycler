# Automated Autocycler Bash script (by Ryan Wick)

This `autocycler_full.sh` Bash script runs a complete Autocycler assembly workflow from start to finish. It's minimalistic by design, with no frills or read quality control – just a straightforward way to go from reads to consensus assembly.

I wrote the first version in January 2025 (for Autocycler v0.2.1) and updated it in June 2025 for v0.5.0. I'll try to keep it in sync with future Autocycler versions.



## Notes and key features

* No quality control is performed on the input reads – they must be ready for assembly.
* Genome size is estimated using Raven ([details](https://github.com/rrwick/Autocycler/wiki/Genome-size-estimation)).
* Uses four read subsets to generate input assemblies.
* Uses nine assemblers in this order: Raven, Myloasm, miniasm, Flye, metaMDBG, NECAT, NextDenovo, Plassembler and Canu. Faster assemblers come first so you can preview results early in a run.
* With 4 read subsets × 9 assemblers, the script generates 36 input assemblies.
* Runs multiple assemblies in parallel using [GNU Parallel](https://github.com/rrwick/Autocycler/wiki/Parallelising-input-assemblies#gnu-parallel).
* Assemblies are run via the [Autocycler helper](https://github.com/rrwick/Autocycler/wiki/Autocycler-helper) command.
* Assemblies are run with `nice -n 19` to reduce their impact on other system processes.
* Assemblies are run using GNU Parallel's `--timeout` option, which kills any job that takes too long to complete. The timeout is hard-coded to a generous 8 hours by default, but you can change the `max_time` variable near the top of the script to adjust this.
* [Plassembler](https://github.com/gbouras13/plassembler) is included to help recover small plasmids that other long-read assemblers may miss.
  * Plassembler requires a reference database. Its helper will look for it via the `PLASSEMBLER_DB` environment variable, or else in a `plassembler_db` directory inside the active conda environment.
  * Circular Plassembler contigs are given extra clustering weight which helps small plasmids (which may only be assembled by Plassembler) to be included in the final assembly.
* Canu and Flye contigs are given extra consensus weight, since those two assemblers often produce the most accurate sequence.
* The `--min_depth_rel 0.1` filter is used to exclude contigs with <10% of the read depth of the longest contig. This helps to exclude low-level contamination.
* Uses GNU `sed` syntax. On macOS (BSD `sed`), the `sed` commands may need adjustment.



## Dependencies

This script assumes the following are available in your `$PATH`:
* `autocycler`
* [GNU Parallel](https://www.gnu.org/software/parallel): `parallel`
* Long-read assemblers and supporting tools: `canu`, `flye`, `metaMDBG`, `miniasm`, `minimap2`, `minipolish`, `myloasm`, `necat`, `nextDenovo`, `nextPolish`, `plassembler`, `racon`, `raven`

Installing all of these into a single conda environment is usually possible – see [Conda environment file (by Ryan Wick)](https://github.com/rrwick/Autocycler/tree/main/pipelines/Conda_environment_file_by_Ryan_Wick).

If you use separate conda environments for each assembler, modify Step 2 of the script to activate/deactivate as needed. For example:
```bash
for i in 01 02 03 04; do
    conda activate canu
    autocycler helper canu --reads subsampled_reads/sample_"$i".fastq --out_prefix assemblies/canu_"$i" --threads "$threads" --genome_size "$genome_size"
    conda deactivate

    conda activate flye
    autocycler helper flye --reads subsampled_reads/sample_"$i".fastq --out_prefix assemblies/flye_"$i" --threads "$threads" --genome_size "$genome_size"
    conda deactivate

    # and so on...
done
```



## Usage

The script takes the following three required arguments:
1. **Reads**: Path to the input FASTQ file (can be gzipped).
2. **Threads**: Number of threads to use per assembly.
3. **Jobs**: Number of assemblies to run in parallel.

And one optional argument:
4. **Read type**: either `ont_r9`, `ont_r10`, `pacbio_clr` or `pacbio_hifi`. If not specified, then `ont_r10` is used.

Note: each assembly will use up to the specified thread count, so the maximum total threads in use is threads × jobs. A job count of 4 works well and minimises mixing assemblers in the same batch. This prevents slow assemblers (e.g. Canu) from dragging out the end of the run when only one job remains active.

**Example command:**
```bash
autocycler_full.sh reads.fastq.gz 16 4
```



## Output

The script creates the following outputs in the working directory:

* **`subsampled_reads/`**: Contains the read subsets used for assembly. The FASTQ files are deleted after assembly to save space, but the directory (and its [YAML file](https://github.com/rrwick/Autocycler/wiki/Metrics#read-subsampling-metrics)) remains.
* **`assemblies/`**:
  * Input assemblies for Autocycler in FASTA format
  * GFA and log files for assemblers that produce them
  * `joblog.tsv`: job metadata (e.g. commands, timing, exit codes)
  * `logs/`: full stdout/stderr for each job
* **`autocycler_out/`**: Output directory for Autocycler. Final results are saved as `consensus_assembly.gfa` and `consensus_assembly.fasta`.
* **`autocycler.stderr`**: stderr output from all Autocycler steps.
