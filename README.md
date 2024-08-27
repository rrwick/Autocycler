# Autocycler




## Usage

```bash
# Compress input sequences into a unitig graph:
autocycler compress -i assemblies -a autocycler

# Cluster input contigs:
autocycler cluster -a autocycler

# For each QC-pass cluster:
for c in autocycler/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c"
    autocycler resolve -c "$c"
done

# Combine clusters into a final assembly:
autocycler combine -i autocycler/clustering/qc_pass/cluster_*/5_final.gfa -o autocycler/final
```


## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
