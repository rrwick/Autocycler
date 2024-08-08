# Autocycler




## Usage

```bash
# Compress input sequences into a unitig graph:
autocycler compress -i assemblies -o autocycler

# Cluster input contigs:
autocycler cluster -o autocycler

# For each QC-pass cluster:
for c in autocycler/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c"
    autocycler resolve -c "$c"
    autocycler correct -c "$c"
done

# Reconstruct input assemblies from unitig graph:
autocycler decompress -i autocycler/input_assemblies.gfa -o reconstructed
```


## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
