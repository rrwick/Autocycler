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
mod tests;
mod test_gfa;
mod trim;
mod unitig;
mod unitig_graph;

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

        /// manually define clusters using tree node numbers
        #[clap(long = "manual", hide_default_value = true,
               help = "manually define clusters using tree node numbers [default: automatic]")]
        manual: Option<String>,
    },

    /// combine Autocycler GFAs into one assembly
    Combine {
        /// Autocycler GFA files (one or more required)
        #[clap(short = 'i', long = "in_gfas", required = true, num_args = 1..)]
        in_gfas: Vec<PathBuf>,

        /// Output prefix (required)
        #[clap(short = 'o', long = "out_prefix", required = true)]
        out_prefix: PathBuf,
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

    /// generate an all-vs-all dot plot from a unitig graph
    Dotplot {
        /// Autocycler GFA file (required)
        #[clap(short = 'i', long = "in_gfa", required = true)]
        in_gfa: PathBuf,

        /// File path where dot plot PNG will be saved (required)
        #[clap(short = 'o', long = "out_png", required = true)]
        out_png: PathBuf,

        /// Size (in pixels) of dot plot image
        #[clap(long = "res", default_value = "2000")]
        res: u32,

        /// K-mer size to use in dot plot
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

    // TODO: add subset command

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
        Some(Commands::Compress { assemblies_dir, autocycler_dir, kmer, threads }) => {
            compress::compress(assemblies_dir, autocycler_dir, kmer, threads);
        },
        Some(Commands::Decompress { in_gfa, out_dir, out_file }) => {
            decompress::decompress(in_gfa, out_dir, out_file);
        },
        Some(Commands::Dotplot { in_gfa, out_png, res, kmer }) => {
            dotplot::dotplot(in_gfa, out_png, res, kmer);
        },
        Some(Commands::Cluster { autocycler_dir, cutoff, min_assemblies, manual }) => {
            cluster::cluster(autocycler_dir, cutoff, min_assemblies, manual);
        },
        Some(Commands::Trim { cluster_dir, min_identity, max_unitigs, mad, threads }) => {
            trim::trim(cluster_dir, min_identity, max_unitigs, mad, threads);
        },
        Some(Commands::Resolve { cluster_dir, verbose }) => {
            resolve::resolve(cluster_dir, verbose);
        },
        Some(Commands::Combine { in_gfas, out_prefix }) => {
            combine::combine(in_gfas, out_prefix);
        },
        None => {}
    }
}
