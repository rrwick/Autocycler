# Conda environment file (by Ryan Wick)

This directory contains an `environment.yml` file for installing Autocycler along with a collection of long-read assemblers and related tools, all into a single conda environment.



## Usage

To create the environment, run:

```bash
conda env create --file environment.yml --name autocycler
conda activate autocycler
plassembler download -d "$CONDA_PREFIX"/plassembler_db
```

All tools will be installed into the same environment (`autocycler`).



## Notes

* If you encounter dependency conflicts on your platform, try commenting out one or more lines in the `environment.yml` file (e.g. for tools you don't plan to use).
* [Plassembler](https://github.com/gbouras13/plassembler) requires a database. The above `plassembler download` command saves it in the conda environment directory.
* You can use `mamba` instead of `conda` in the installation command if you have [mamba](https://mamba.readthedocs.io/en/latest) installed. It resolves dependencies more quickly.
