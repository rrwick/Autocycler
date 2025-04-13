# Automated Autocycler Bash script (by Ryan Wick)

This `autocycler_full.sh` Bash script runs a complete Autocycler assembly workflow from start to finish. It's minimalistic by design, with no frills or read quality control – just a straightforward way to go from reads to consensus assembly.

I wrote the first version in January 2025 (for Autocycler v0.2.1) and updated it in April 2025 for v0.3.1. I’ll try to keep it in sync with future Autocycler versions, but no guarantees.



## Notes and key features

* No quality control is performed on the input reads – they must be ready for assembly.
* Genome size is estimated using Raven ([details](https://github.com/rrwick/Autocycler/wiki/Genome-size-estimation)).
* Uses four hard-coded read subsets for the input assemblies.
* Runs multiple assemblies in parallel using [GNU Parallel](https://github.com/rrwick/Autocycler/wiki/Parallelising-input-assemblies#gnu-parallel).
* Assemblies are run with `nice -n 19` to reduce their impact on other system processes.
* Uses eight assemblers in this order: Raven, miniasm, Flye, MetaMDBG, NECAT, NextDenovo, Plassembler and Canu. Faster assemblers come first so you can preview results early in a run.
* Assemblies are launched via the [Autocycler helper scripts](https://github.com/rrwick/Autocycler/wiki/Generating-input-assemblies#assembly-helper-scripts).
* With 4 read subsets × 8 assemblers, the script generates 32 input assemblies.
* [Plassembler](github.com/gbouras13/plassembler) is included to help recover small plasmids that other long-read assemblers may miss.
* Plassembler requires a reference database. Its helper script will look for it via the `PLASSEMBLER_DB` environment variable, or else in a `plassembler_db` directory inside the active conda environment.
* Circular Plassembler contigs are given extra clustering weight to help small plasmids (which may only be assembled by Plassembler) be included in the final assembly.
* Canu and Flye contigs are given extra consensus weight, since those two assemblers often produce the most accurate sequence.
* This script uses GNU sed syntax. On macOS (BSD sed), the `sed` commands may need adjustment.



## Dependencies

This script assumes the following are available in your `$PATH`:
* `autocycler`
* [Autocycler's helper scripts](https://github.com/rrwick/Autocycler/tree/main/scripts): `raven.sh`, `flye.sh`, etc.
* [GNU Parallel](https://www.gnu.org/software/parallel): `parallel`
* Long-read assemblers and supporting tools: `any2fasta`, `canu`, `flye`, `metaMDBG`, `miniasm`, `minimap2`, `minipolish`, `necat`, `nextDenovo`, `nextPolish`, `plassembler`, `racon`, `raven`, `seqtk`

The last point is often the trickiest, especially if you want everything installed in a single conda environment. It can be done, but may require some fiddling – see [Autocycler's installation instructions](https://github.com/rrwick/Autocycler/wiki/Software-requirements-and-installation) for guidance.

If you prefer to use separate conda environments (or can't get all dependencies into one), you'll need to modify the `Step 2: assemble each subsampled file` section of the script to activate/deactivate the appropriate conda environment for each assembler. For example:
```bash
for i in 01 02 03 04; do
    conda activate canu
    canu.sh subsampled_reads/sample_"$i".fastq assemblies/canu_"$i" $threads $genome_size
    conda deactivate

    conda activate flye
    flye.sh subsampled_reads/sample_"$i".fastq assemblies/flye_"$i" $threads $genome_size
    conda deactivate

    # and so on...
done
```



## Usage

The script takes the following three arguments:
1. **Read filename**: Path to the input FASTQ file (can be gzipped).
2. **Thread count**: Number of threads to use per assembly.
3. **Job count**: Number of assemblies to run in parallel.

Note: each assembly will use up to the specified thread count, so the maximum total threads in use is `threads × jobs`. A job count of 4 works well, matching the four read subsets. This keeps each batch of parallel jobs within the same assembler, avoiding situations where one slow job delays everything else.

**Example command:**
```bash
autocycler_full.sh reads.fastq.gz 16 4
```



## Output

The script creates the following outputs in the working directory:

* **`subsampled_reads/`**: Contains the read subsets used for assembly. The FASTQ files are deleted after assembly to save space, but the directory (and its [YAML file](https://github.com/rrwick/Autocycler/wiki/Metrics#read-subsampling-metrics)) remains.
* **`assemblies/`**: Contains the input assemblies for Autocycler. Logs from each assembler are saved in `assemblies/logs/`.
* **`autocycler_out/`**: Output directory for Autocycler. Final results are saved as `consensus_assembly.gfa` and `consensus_assembly.fasta`.
* **`autocycler.stderr`**: Contains `stderr` output from all Autocycler steps.
