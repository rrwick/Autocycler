# Fast Autocycler with Flye fallback (by Ryan Wick)

This `autocycler_and_flye.py` Python script builds both Flye and Autocycler assemblies. If the Autocycler assembly was completely successful, it uses that as the final result. Otherwise, it returns the Flye assembly as the final result.

This script requires Autocycler v0.7.0 or later.



## Notes and key features

* No quality control is performed on the input reads – they must be ready for assembly.
* Basic input read statistics (count, total bases and N50) are reported in `assembly.log`.
* Genome size is estimated using [LRGE](https://github.com/mbhall88/lrge), unless supplied with  `--genome_size`.
* [Rasusa](https://github.com/mbhall88/rasusa) downsamples the reads to 100× depth for the initial Flye assembly. This is because Flye can crash with excessively deep read sets. If the input depth is already ≤100×, this step is skipped.
* Autocycler's read subsets are made from the original full read set, not the Rasusa-downsampled reads.
* The fast Autocycler workflow uses six assemblers: Flye, metaMDBG, miniasm, Myloasm, Plassembler and Raven.
* Two read subsets are used by default, producing up to 12 Autocycler input assemblies. The initial Flye assembly is also included, giving up to 13 input assemblies in total when all jobs succeed.
* Input assemblies are generated with the [Autocycler helper](https://github.com/rrwick/Autocycler/wiki/Autocycler-helper) and run in parallel using [GNU Parallel](https://github.com/rrwick/Autocycler/wiki/Parallelising-input-assemblies#gnu-parallel).
* Input assembly jobs are run with `nice -n 19` and a default timeout of 4 hours.
* Circular Plassembler contigs receive extra clustering weight, helping small plasmids which may only be assembled by Plassembler.
* Flye contigs receive extra consensus weight.
* The `--min_depth_rel 0.1` filter excludes contigs with <10% of the read depth of the longest contig, helping to remove low-level contamination.
* The Autocycler workflow follows the usual subsample, compress, cluster, trim, resolve and combine steps. See [Fully automated assembly](https://github.com/rrwick/Autocycler/wiki/Fully-automated-assembly) for an overview.
* An [`autocycler table`](https://github.com/rrwick/Autocycler/wiki/Autocycler-table) summary is saved as `autocycler/metrics.tsv`.



## Assembly selection

The Flye and Autocycler assemblies are assessed against the estimated (or user-supplied) genome size. By default, an acceptable assembly must be between 0.75× and 1.25× the genome size. These limits can be changed with `--min-size-ratio` and `--max-size-ratio`.

The final assembly is selected as follows:

1. Use Autocycler if its size is acceptable and `autocycler combine` reported `Consensus assembly is fully resolved`.
2. Otherwise, use Flye if its size is acceptable.
3. If neither assembly is acceptable, exit with an error.

The selected assembly is copied to `assembly.fasta` and `assembly.gfa` in the output directory.



## Dependencies

This script uses only the Python standard library, but assumes the following command-line tools are available in your `$PATH`:

* `autocycler`
* [LRGE](https://github.com/mbhall88/lrge): `lrge` (not required when using `--genome_size`)
* [Rasusa](https://github.com/mbhall88/rasusa): `rasusa`
* [GNU Parallel](https://www.gnu.org/software/parallel): `parallel`
* `nice`
* Long-read assemblers and supporting tools: `flye`, `metaMDBG`, `miniasm`, `minimap2`,
  `minipolish`, `myloasm`, `plassembler`, `racon`, `raven`

Installing these into a single conda environment is usually possible – see [Conda environment file (by Ryan Wick)](https://github.com/rrwick/Autocycler/tree/main/pipelines/Conda_environment_file_by_Ryan_Wick).



## Usage

The script takes two positional arguments:

1. **Reads**: path to the input FASTQ file (can be gzipped).
2. **Output directory**: directory for working files and final outputs. This directory must not already exist.

**Example command:**
```bash
./autocycler_and_flye.py reads.fastq.gz out_dir
```

**Full usage:**
```
usage: autocycler_and_flye.py [--read-type {ont_r9,ont_r10,pacbio_clr,pacbio_hifi}] [--genome_size GENOME_SIZE]
                              [--min-size-ratio MIN_SIZE_RATIO] [--max-size-ratio MAX_SIZE_RATIO] [--seed SEED]
                              [--subset_count SUBSET_COUNT] [--threads THREADS] [--jobs JOBS]
                              [--max-job-time MAX_JOB_TIME] [--clean {0,1,2,3}] [-h]
                              reads out_dir

Fast Autocycler with Flye fallback

Positional arguments:
  reads                 Read FASTQ file (can be gzipped)
  out_dir               Directory for working files and final assembly (will be created)

Settings:
  --read-type {ont_r9,ont_r10,pacbio_clr,pacbio_hifi}
                        Type of long reads (default: ont_r10)
  --genome_size GENOME_SIZE
                        Genome size in bp (skips LRGE when supplied) (default: None)
  --min-size-ratio MIN_SIZE_RATIO
                        Reject assemblies smaller than this multiple of genome size (default: 0.75)
  --max-size-ratio MAX_SIZE_RATIO
                        Reject assemblies larger than this multiple of genome size (default: 1.25)
  --seed SEED           Random seed for reproducible read subsampling (default: 0)
  --subset_count SUBSET_COUNT
                        Number of Autocycler read subsets (default: 2)

Resources:
  --threads THREADS     Maximum number of CPU threads (default: 16)
  --jobs JOBS           Number of simultaneous Autocycler input assembly jobs (default: 4)
  --max-job-time MAX_JOB_TIME
                        Maximum runtime for each Autocycler input assembly job (default: 4h)

Output:
  --clean {0,1,2,3}     Cleanup level from 0 (keep everything) to 3 (keep only final assembly and logs)
                        (default: 2)

Help:
  -h, --help            Show this help message and exit
```



## Output

The most important outputs are:

* **`assembly.fasta`**: final selected assembly in FASTA format.
* **`assembly.gfa`**: final selected assembly in GFA format.
* **`assembly.log`**: concise, timestamped pipeline summary.
* **`logs/`**: detailed tool logs:
  * `lrge.log` when LRGE was run
  * `rasusa.log` when Rasusa was run
  * `flye.log`
  * `autocycler.log`

Depending on the cleanup level, the output directory can also contain:

* **`rasusa_reads.fastq.gz`**: reads downsampled to 100× for Flye (only created when needed and only
  retained with `--clean 0`).
* **`flye/`**: Flye assembly and working files.
* **`autocycler/`**:
  * `consensus_assembly.fasta`, `consensus_assembly.gfa` and
    `consensus_assembly.yaml`
  * `input_assemblies.gfa` and `input_assemblies.yaml`
  * `metrics.tsv`: one-row summary of key Autocycler metrics, with a header
  * `assemblies/`: Autocycler input assemblies, GNU Parallel job metadata and per-assembly logs
  * `clustering/`: Autocycler cluster, trim and resolve outputs
  * `subsampled_reads/`: Autocycler read subsets and subsampling metrics

With the default `--clean 2`, the final assemblies and logs are retained, along with the top-level Flye and Autocycler results, but their large working directories are removed.
