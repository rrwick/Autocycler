// This is the main file of Autocycler and where execution starts. It mainly handles the CLI and
// then calls into other files to run whichever subcommand the user chose.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::path::PathBuf;
use clap::{Parser, Subcommand, crate_version};

mod clean;
mod cluster;
mod combine;
mod compress;
mod decompress;
mod dotplot;
mod graph_simplification;
mod kmer_graph;
mod log;
mod metrics;
mod misc;
mod position;
mod resolve;
mod sequence;
mod subsample;
mod table;
mod test_gfa;
mod trim;
mod unitig;
mod unitig_graph;

#[cfg(test)]
mod tests;

#[derive(Parser)]
#[clap(name = "Autocycler",
       version = concat!("v", crate_version!()),
       about = "a tool for generating consensus bacterial genome assemblies\n\
                Documenation: https://github.com/rrwick/Autocycler/wiki",
       before_help = concat!(r#"                _                        _           "#, "\n",
                             r#"     /\        | |                      | |          "#, "\n",
                             r#"    /  \  _   _| |_ ___   ___ _   _  ___| | ___ _ __ "#, "\n",
                             r#"   / /\ \| | | | __/ _ \ / __| | | |/ __| |/ _ \ '__|"#, "\n",
                             r#"  / ____ \ |_| | || (_) | (__| |_| | (__| |  __/ |   "#, "\n",
                             r#" /_/    \_\__,_|\__\___/ \___|\__, |\___|_|\___|_|   "#, "\n",
                             r#"                               __/ |                 "#, "\n",
                             r#"                              |___/                  "#))]
#[command(author, version, long_about = None, disable_help_subcommand = true,
          propagate_version = true)]
#[clap(subcommand_required = true)]
#[clap(arg_required_else_help = true)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {

    /// manual manipulation of the final consensus assebly graph
    Clean {
        /// Autocycler GFA file (required)
        #[clap(short = 'i', long = "in_gfa", required = true)]
        in_gfa: PathBuf,

        /// Output GFA file (required)
        #[clap(short = 'o', long = "out_gfa", required = true)]
        out_gfa: PathBuf,

        /// Tig numbers to remove from the input graph
        #[clap(short = 'r', long = "remove")]
        remove: Option<String>,

        /// Tig numbers to duplication in the input graph
        #[clap(short = 'd', long = "duplicate")]
        duplicate: Option<String>,
    },

    /// cluster contigs in the unitig graph based on similarity
    Cluster {
        /// Autocycler directory containing input_assemblies.gfa file (required)
        #[clap(short = 'a', long = "autocycler_dir", required = true)]
        autocycler_dir: PathBuf,

        /// cutoff distance threshold for hierarchical clustering
        #[clap(long = "cutoff", default_value = "0.2")]
        cutoff: f64,

        /// exclude clusters with fewer than this many assemblies
        #[clap(long = "min_assemblies", hide_default_value = true,
               help = "exclude clusters with fewer than this many assemblies [default: automatic]")]
        min_assemblies: Option<usize>,

        /// refuse to run if mean contigs per assembly exceeds this value
        #[clap(long = "max_contigs", default_value = "25")]
        max_contigs: u32,

        /// manually define clusters using tree node numbers
        #[clap(long = "manual", hide_default_value = true,
               help = "manually define clusters using tree node numbers [default: automatic]")]
        manual: Option<String>,
    },

    /// combine Autocycler GFAs into one assembly
    Combine {
        /// Autocycler directory (required)
        #[clap(short = 'a', long = "autocycler_dir", required = true)]
        autocycler_dir: PathBuf,

        /// Autocycler cluster GFA files (one or more required)
        #[clap(short = 'i', long = "in_gfas", required = true, num_args = 1..)]
        in_gfas: Vec<PathBuf>,
    },

    /// compress input contigs into a unitig graph
    Compress {
        /// Directory containing input assemblies (required)
        #[clap(short = 'i', long = "assemblies_dir", required = true)]
        assemblies_dir: PathBuf,

        /// Autocycler directory to be created (required)
        #[clap(short = 'a', long = "autocycler_dir", required = true)]
        autocycler_dir: PathBuf,

        /// K-mer size for De Bruijn graph
        #[clap(long = "kmer", default_value = "51")]
        kmer: u32,

        /// refuse to run if mean contigs per assembly exceeds this value
        #[clap(long = "max_contigs", default_value = "25")]
        max_contigs: u32,

        /// Number of CPU threads
        #[clap(short = 't', long = "threads", default_value = "8")]
        threads: usize,
    },

    /// decompress contigs from a unitig graph
    Decompress {
        /// Autocycler GFA file (required)
        #[clap(short = 'i', long = "in_gfa", required = true)]
        in_gfa: PathBuf,

        /// Directory where decompressed sequences will be saved (either -o or -f is required)
        #[clap(short = 'o', long = "out_dir")]
        out_dir: Option<PathBuf>,

        /// FASTA file where decompressed sequences will be saved (either -o or -f is required)
        #[clap(short = 'f', long = "out_file")]
        out_file: Option<PathBuf>,
    },

    /// generate an all-vs-all dotplot from a unitig graph
    Dotplot {
        /// Input Autocycler GFA file, FASTA file or directory (required)
        #[clap(short = 'i', long = "input", required = true)]
        input: PathBuf,

        /// File path where dotplot PNG will be saved (required)
        #[clap(short = 'o', long = "out_png", required = true)]
        out_png: PathBuf,

        /// Size (in pixels) of dotplot image
        #[clap(long = "res", default_value = "2000")]
        res: u32,

        /// K-mer size to use in dotplot
        #[clap(long = "kmer", default_value = "32")]
        kmer: u32,
    },

    /// resolve repeats in the the unitig graph
    Resolve {
        /// Autocycler directory (required)
        #[clap(short = 'c', long = "cluster_dir", required = true)]
        cluster_dir: PathBuf,

        /// Enable verbose output
        #[clap(long = "verbose")]
        verbose: bool,
    },

    /// subsample a long-read set
    Subsample {
        /// Input long reads in FASTQ format (required)
        #[clap(short = 'r', long = "reads", required = true)]
        reads: PathBuf,

        /// Output directory (required)
        #[clap(short = 'o', long = "out_dir")]
        out_dir: PathBuf,

        /// Estimated genome size (required)
        #[clap(short = 'g', long = "genome_size")]
        genome_size: String,

        /// Number of subsampled read sets to output
        #[clap(short = 'c', long = "count", default_value = "4")]
        count: usize,

        /// Minimum allowed read depth
        #[clap(short = 'd', long = "min_read_depth", default_value = "25.0")]
        min_read_depth: f64,

        /// Seed for random number generator
        #[clap(short = 's', long = "seed", default_value = "0")]
        seed: u64,
    },

    /// create TSV line from YAML files
    Table {
        /// Autocycler directory (if absent, a header line will be output)
        #[clap(short = 'a', long = "autocycler_dir")]
        autocycler_dir: Option<PathBuf>,

        /// Sample name
        #[clap(short = 'n', long = "name", default_value = "", hide_default_value = true,
               help = "Sample name [default: blank]")]
        name: String,

        /// Comma-delimited list of YAML fields to include
        #[clap(short = 'f', long = "fields",
               default_value = "input_read_count, input_read_bases, input_read_n50, \
                                pass_cluster_count, fail_cluster_count, overall_clustering_score, \
                                untrimmed_cluster_size, untrimmed_cluster_distance, \
                                trimmed_cluster_size, trimmed_cluster_median, trimmed_cluster_mad, \
                                consensus_assembly_bases, consensus_assembly_unitigs, \
                                consensus_assembly_fully_resolved")]
        fields: String,

        /// Significant figures to use for floating point numbers
        #[clap(short = 's', long = "sigfigs", default_value = "3")]
        sigfigs: usize,
    },

    /// trim contigs in a cluster
    Trim {
        /// Autocycler cluster directory containing 1_untrimmed.gfa file (required)
        #[clap(short = 'c', long = "cluster_dir", required = true)]
        cluster_dir: PathBuf,

        /// Minimum alignment identity for trimming alignments
        #[clap(long = "min_identity", default_value = "0.75")]
        min_identity: f64,

        /// Maximum unitigs used for overlap alignment, set to 0 to disable trimming
        #[clap(long = "max_unitigs", default_value = "5000")]
        max_unitigs: usize,

        /// Allowed variability in cluster length, measured in median absolute deviations, set to 0 to disable exclusion of length outliers
        #[clap(long = "mad", default_value = "5.0")]
        mad: f64,

        /// Number of CPU threads
        #[clap(short = 't', long = "threads", default_value = "8")]
        threads: usize,
    },
}


fn main() {
    let cli = Cli::parse();

    match cli.command {
        Some(Commands::Clean { in_gfa, out_gfa, remove, duplicate }) => {
            clean::clean(in_gfa, out_gfa, remove, duplicate);
        },
        Some(Commands::Cluster { autocycler_dir, cutoff, min_assemblies, max_contigs, manual }) => {
            cluster::cluster(autocycler_dir, cutoff, min_assemblies, max_contigs, manual);
        },
        Some(Commands::Combine { autocycler_dir, in_gfas }) => {
            combine::combine(autocycler_dir, in_gfas);
        },
        Some(Commands::Compress { assemblies_dir, autocycler_dir, kmer, max_contigs, threads }) => {
            compress::compress(assemblies_dir, autocycler_dir, kmer, max_contigs, threads);
        },
        Some(Commands::Decompress { in_gfa, out_dir, out_file }) => {
            decompress::decompress(in_gfa, out_dir, out_file);
        },
        Some(Commands::Dotplot { input, out_png, res, kmer }) => {
            dotplot::dotplot(input, out_png, res, kmer);
        },
        Some(Commands::Resolve { cluster_dir, verbose }) => {
            resolve::resolve(cluster_dir, verbose);
        },
        Some(Commands::Subsample { reads, out_dir, genome_size, count, min_read_depth, seed }) => {
            subsample::subsample(reads, out_dir, genome_size, count, min_read_depth, seed);
        },
        Some(Commands::Table { autocycler_dir, name, fields, sigfigs }) => {
            table::table(autocycler_dir, name, fields, sigfigs);
        },
        Some(Commands::Trim { cluster_dir, min_identity, max_unitigs, mad, threads }) => {
            trim::trim(cluster_dir, min_identity, max_unitigs, mad, threads);
        },
        None => {}
    }
}
