# Single-assembly trimming (by Ryan Wick)

This script came from a discussion I had with [Torsten Seemann](https://github.com/tseemann). We were chatting about how [Flye](https://github.com/mikolmogorov/Flye) assemblies often contain duplicated small plasmids. E.g. a genome might contain a 3 kbp plasmid but Flye will produce a 6 kbp contig with two entire copies. This is a well known issue with Flye and other long-read assemblers ([issue](https://github.com/mikolmogorov/Flye/issues/706), [issue](https://github.com/mikolmogorov/Flye/issues/769), [repo](https://github.com/ncezid-narst/messer), [paper](https://doi.org/10.1099/mgen.0.001024), [paper](https://www.nature.com/articles/s41467-021-21656-9), [paper](https://link.springer.com/article/10.1186/s12864-023-09729-z)). I mentioned to Torsten that Autocycler can usually fix these plasmid duplications, thanks to its [`autocycler trim`](https://github.com/rrwick/Autocycler/wiki/Autocycler-trim) step. He asked if there was any way to perform this on just a single assembly (not a full Autocycler pipeline).

It turns out that the answer is yes! But it took a bit of scripting.

The `autocycler_trim_single.py` script in this directory will:
- Take an assembly FASTA as input.
- Run Autocycler trimming (`compress` → `cluster` → `trim` → `decompress`) individually on each sequence in that assembly. Doing this one sequence at a time ensures that sequences are not lost during the clustering step.
- Save the result as a FASTA in an output directory, using the same sequence headers as the input assembly.

Quick usage:
```bash
python3 autocycler_trim_single.py assembly.fasta out_dir
```

Full usage:
```
usage: autocycler_trim_single.py [--min_identity MIN_IDENTITY] [--keep] [-h] assembly out_dir

Trim a single assembly with Autocycler

Positional arguments:
  assembly              Input assembly FASTA (can be gzipped)
  out_dir               Output directory (must not already exist)

Settings:
  --min_identity MIN_IDENTITY
                        Minimum unitig alignment identity for trimming (default: 0.75)

Output:
  --keep                Keep intermediate files (default: False)

Help:
  -h, --help            Show this help message and exit
```

It's pretty quick: on my laptop, it usually takes less than 10 seconds to run.

Requirements:
- Python 3.8+
- Autocycler (tested with Autocycler v0.8.0)

Options:
- `--min_identity`: Minimum unitig alignment identity for trimming (default: 0.75). Note that this is unitig-level identity (not sequence-level identity), so 0.75 represents a reasonably strong alignment. Increase to 1.0 if you want to limit trimming to exact overlaps.
- `--keep`: Keep intermediate files.

See the [Autocycler trim docs](https://github.com/rrwick/Autocycler/wiki/Autocycler-trim) for more details.

After successfully running, the output directory will contain:
- `assembly.fasta`: the assembly with any detected overlaps removed.
- `autocycler.log`: all Autocycler output.
- `working_files`: a directory with the intermediate files (if `--keep` was used).

Caveats:
- Genuine duplications/repeats can be trimmed. E.g. if you provide a sequence that looks like this: `REPEAT→MIDDLE→REPEAT`, the output will contain just one copy of the repeat.
- Short and/or inexact overlaps can be missed.
- Excessively short sequences cannot be trimmed. Any sequence below 51 bp (the default Autocycler _k_-mer size) will be copied to the output without any trimming.
- Sequences with non-canonical characters (anything other than `A`, `C`, `G` or `T`) cannot be trimmed and will be copied to the output unchanged.
- It will check for both start-end (circular) overlaps and hairpin overlaps.
- Trimming can remove sequence from both ends, so the starting position can change.
- Sequence case (upper vs lower) and FASTA wrapping are not preserved.
- Any sequences with `Autocycler_ignore` in their header will cause the script to fail.
- Header metadata may become stale after running this. For example, if the header contains `length=8424`, that will be copied to the output FASTA, even if the sequence was trimmed to a shorter length.
