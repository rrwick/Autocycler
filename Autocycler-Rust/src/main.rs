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
mod compress;
mod decompress;
mod kmer_graph;
mod log;
mod misc;
mod position;
mod resolve;
mod sequence;
mod tests;
mod unitig;
mod unitig_graph;
mod graph_simplification;


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
    /// compress input assemblies into a De Bruijn graph
    Compress {
        /// Directory containing input assemblies (required)
        #[clap(short = 'i', long = "in_dir", required = true)]
        in_dir: PathBuf,

        /// Filename of output assembly graph (required)
        #[clap(short = 'o', long = "out_gfa", required = true)]
        out_gfa: PathBuf,

        /// K-mer size for De Bruijn graph
        #[clap(short = 'k', long = "kmer", default_value = "51")]
        kmer: u32,
    },

    /// decompress the De Bruijn graph back into assemblies
    Decompress {
        /// Autocycler GFA file (required)
        #[clap(short = 'i', long = "in_gfa", required = true)]
        in_gfa: PathBuf,

        /// Directory where decompressed assemblies will be saved (required)
        #[clap(short = 'o', long = "out_dir", required = true)]
        out_dir: PathBuf,
    },

    /// cluster contigs in the De Bruijn graph based on similarity
    Cluster {
        /// Autocycler GFA file (required)
        #[clap(short = 'i', long = "in_gfa", required = true)]
        in_gfa: PathBuf,

        /// Directory where clustering results will be saved (required)
        #[clap(short = 'o', long = "out_dir", required = true)]
        out_dir: PathBuf,

        /// ε parameter for DBSCAN* clustering
        #[clap(short = 's', long = "eps", hide_default_value = true,
        help = "ε parameter for DBSCAN* clustering [default: automatic]")]
        eps: Option<f64>,

        /// minPts parameter for DBSCAN* clustering
        #[clap(short = 'm', long = "minpts", hide_default_value = true,
               help = "minPts parameter for DBSCAN* clustering [default: automatic]")]
        minpts: Option<usize>,
    },

    /// resolve the De Bruijn graph to a consensus assembly
    Resolve {
        /// Autocycler GFA file (required)
        #[clap(short = 'i', long = "in_gfa", required = true)]
        in_gfa: PathBuf,

        /// Directory where resolved assembly graphs will be saved (required)
        #[clap(short = 'o', long = "out_dir", required = true)]
        out_dir: PathBuf,
    },
}


fn main() {
    let cli = Cli::parse();

    match cli.command {
        Some(Commands::Compress { in_dir, out_gfa, kmer }) => {
            compress::compress(in_dir, out_gfa, kmer);
        },
        Some(Commands::Decompress { in_gfa, out_dir }) => {
            decompress::decompress(in_gfa, out_dir);
        },
        Some(Commands::Cluster { in_gfa, out_dir, eps, minpts }) => {
            cluster::cluster(in_gfa, out_dir, eps, minpts);
        },
        Some(Commands::Resolve { in_gfa, out_dir }) => {
            resolve::resolve(in_gfa, out_dir);
        },
        None => {}
    }
}
