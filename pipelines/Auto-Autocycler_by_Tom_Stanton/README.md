# autoautocycler

`autoautocycler` is a simple shell wrapper around [Autocycler](https://github.com/rrwick/Autocycler.git) for running all steps on multiple assemblies with a single command.

## Installation
No installation required! Just make sure [Autocycler](https://github.com/rrwick/Autocycler.git) is installed and make the script excecutable with 
`chmod +x autoautocycler.sh`.

## Usage
Running `./autoautocycler.sh` with `-h`/`--help` will print all options:

```sh
Usage: autoautocycler <reads> [<reads> ...] [options]
  <reads>            Long reads to assemble in fastq(.gz) format
  -o, --out          Output directory (default: .)
  -t, --threads      #threads to use (default: all available)
  -c, --count        #subsampled read sets to output (default: 4)
  -k, --kmer         K-mer size for De Bruijn graph (default: 51)
  -s, --size         Genome size (default: AUTO)
  -a, --assemblers   Assemblers to use (default: all available)
                     Possible assemblers: canu flye lja metamdbg miniasm necat nextdenovo raven redbean
                     Note: this argument MUST BE WRAPPED in quotes
```

E.g. for making assemblies with flye, miniasm and raven, just run:
```sh
./autoautocycler.sh -o autocycler -c 2 -s 5500000 -a "flye miniasm raven" reads/*.fastq.gz
```
Here I've set the subsampling count to 2 and set the genome size to 5500000 to skip the automatic genome size calculation with raven.

## Bug reporting
Bugs can be reported using the [Autocycler issues page](https://github.com/rrwick/Autocycler/issues) on GitHub,
just make sure you tag me with @tomdstanton.
