// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

mod compress;
mod decompress;
mod kmer_graph;
mod log;
mod misc;
mod position;
mod sequence;
mod resolve;

use std::path::PathBuf;
use std::collections::HashMap;
use std::time::Instant;
use std::fs::File;
use std::io::prelude::*;
use clap::{Parser, ArgGroup, Subcommand, crate_version};
use num_format::{Locale, ToFormattedString};


#[derive(Parser)]
#[clap(name = "Autocycler",
       version = concat!("v", crate_version!()),
       about = "consensus bacterial genome assemblies\ngithub.com/rrwick/Autocycler",
       before_help = concat!(r#"                _                        _           "#, "\n",
                             r#"     /\        | |                      | |          "#, "\n",
                             r#"    /  \  _   _| |_ ___   ___ _   _  ___| | ___ _ __ "#, "\n",
                             r#"   / /\ \| | | | __/ _ \ / __| | | |/ __| |/ _ \ '__|"#, "\n",
                             r#"  / ____ \ |_| | || (_) | (__| |_| | (__| |  __/ |   "#, "\n",
                             r#" /_/    \_\__,_|\__\___/ \___|\__, |\___|_|\___|_|   "#, "\n",
                             r#"                               __/ |                 "#, "\n",
                             r#"                              |___/                  "#))]
#[command(author, version, about, long_about = None, disable_help_subcommand = true,
          propagate_version = true)]
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

    /// decompress De Bruijn graph back into assemblies
    Decompress {
        /// Autocycler GFA file (required)
        #[clap(short = 'i', long = "in_gfa", required = true)]
        in_gfa: PathBuf,

        /// Directory where decompressed assemblies will be saved (required)
        #[clap(short = 'o', long = "out_dir", required = true)]
        out_dir: PathBuf,
    },

    /// cresolve De Bruijn graph to a consensus assembly
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
        Some(Commands::Resolve { in_gfa, out_dir }) => {
            resolve::resolve(in_gfa, out_dir);
        },
        None => {}
    }
}
