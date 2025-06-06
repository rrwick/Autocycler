// This file contains the code for the autocycler helper subcommand.

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
use clap::ValueEnum;

use crate::misc::{check_if_file_exists, quit_with_error};
use crate::subsample::parse_genome_size;


#[derive(ValueEnum, Clone, Debug)]
pub enum Task {
    Genomesize,  // calculate genome size using a Raven assembly
    Canu,        // assemble using Canu and clean results
    Flye,        // assemble using Flye
    Lja,         // assemble using LJA
    Metamdbg,    // assemble using metaMDBG
    Miniasm,     // assemble using miniasm and Minipolish
    Myloasm,     // assemble using Myloasm
    Necat,       // assemble using NECAT
    Nextdenovo,  // assemble using NextDenovo and NextPolish
    Plassembler, // assemble using Plassembler
    Raven,       // assemble using Raven
    Redbean,     // assemble using Redbean (aka wtdbg2)
}


pub fn helper(task: Task, reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: Option<PathBuf>, args: String) {
    check_if_file_exists(&reads);
    match task {
        Task::Genomesize => { genome_size_raven(reads, threads, dir, args); }
        Task::Canu => { canu(reads, out_prefix, genome_size, threads, dir, args); }
        Task::Flye => { flye(reads, out_prefix, threads, dir, args); }
        Task::Lja => { lja(reads, out_prefix, threads, dir, args); }
        Task::Metamdbg => { metamdbg(reads, out_prefix, threads, dir, args); }
        Task::Miniasm => { miniasm(reads, out_prefix, threads, dir, args); }
        Task::Myloasm => { myloasm(reads, out_prefix, threads, dir, args); }
        Task::Necat => { necat(reads, out_prefix, genome_size, threads, dir, args); }
        Task::Nextdenovo => { nextdenovo(reads, out_prefix, genome_size, threads, dir, args); }
        Task::Plassembler => { plassembler(reads, out_prefix, threads, dir, args); }
        Task::Raven => { raven(reads, out_prefix, threads, dir, args); }
        Task::Redbean => { redbean(reads, out_prefix, genome_size, threads, dir, args); }
    }
}


fn genome_size_raven(reads: PathBuf, threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn canu(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
        threads: usize, dir: Option<PathBuf>, args: String) {
    let genome_size = get_genome_size(genome_size, "Canu");
    // TODO
}


fn flye(reads: PathBuf, out_prefix: Option<PathBuf>,
        threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn lja(reads: PathBuf, out_prefix: Option<PathBuf>,
       threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn metamdbg(reads: PathBuf, out_prefix: Option<PathBuf>,
            threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn miniasm(reads: PathBuf, out_prefix: Option<PathBuf>,
           threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn myloasm(reads: PathBuf, out_prefix: Option<PathBuf>,
           threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn necat(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
         threads: usize, dir: Option<PathBuf>, args: String) {
    let genome_size = get_genome_size(genome_size, "NECAT");
    // TODO
}


fn nextdenovo(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: Option<PathBuf>, args: String) {
    let genome_size = get_genome_size(genome_size, "NextDenovo");
    // TODO
}


fn plassembler(reads: PathBuf, out_prefix: Option<PathBuf>,
               threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn raven(reads: PathBuf, out_prefix: Option<PathBuf>,
         threads: usize, dir: Option<PathBuf>, args: String) {
    // TODO
}


fn redbean(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
           threads: usize, dir: Option<PathBuf>, args: String) {
    let genome_size = get_genome_size(genome_size, "Redbean");
    // TODO
}


fn get_genome_size(genome_size: Option<String>, assembler_name: &str) -> u64 {
    if genome_size.is_none() {
        quit_with_error(&format!("assembly with {assembler_name} requires --genomesize"));
    }
    parse_genome_size(&genome_size.unwrap())
}
