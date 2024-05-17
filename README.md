# Autocycler




## Usage

```bash
# Create a consensus assembly:
autocycler compress -i assemblies -o autocycler
autocycler cluster -o autocycler
autocycler resolve -o autocycler
autocycler correct -o autocycler

# Reconstruct input assemblies from unitig graph:
autocycler decompress -i autocycler/01_unitig_graph.gfa -o reconstructed
```


## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
