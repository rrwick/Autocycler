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

use clap::ValueEnum;
use std::fs::{OpenOptions, remove_file, create_dir_all};
use std::io::ErrorKind;
use std::path::PathBuf;
use which::which;
use tempfile::{tempdir, NamedTempFile};

use crate::misc::{check_if_file_exists, quit_with_error};
use crate::subsample::parse_genome_size;


pub fn helper(task: Task, reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: Option<PathBuf>, read_type: ReadType, args: String) {
    check_if_file_exists(&reads);
    let dir = get_working_dir(dir);
    match task {
        Task::Genomesize => {
            genome_size_raven(reads, threads, dir, args);
        }
        Task::Canu => {
            canu(reads, out_prefix, genome_size, threads, dir, read_type, args);
        }
        Task::Flye => {
            flye(reads, out_prefix, threads, dir, read_type, args);
        }
        Task::Lja => {
            lja(reads, out_prefix, threads, dir, args);
        }
        Task::Metamdbg => {
            metamdbg(reads, out_prefix, threads, dir, read_type, args);
        }
        Task::Miniasm => {
            miniasm(reads, out_prefix, threads, dir, read_type, args);
        }
        Task::Myloasm => {
            myloasm(reads, out_prefix, threads, dir, read_type, args);
        }
        Task::Necat => {
            necat(reads, out_prefix, genome_size, threads, dir, read_type, args);
        }
        Task::Nextdenovo => {
            nextdenovo(reads, out_prefix, genome_size, threads, dir, read_type, args);
        }
        Task::Plassembler => {
            plassembler(reads, out_prefix, threads, dir, read_type, args);
        }
        Task::Raven => {
            raven(reads, out_prefix, threads, dir, args);
        }
        Task::Redbean => {
            redbean(reads, out_prefix, genome_size, threads, dir, read_type, args);
        }
    }
}


#[derive(ValueEnum, Clone, Debug)]
pub enum Task {
    Genomesize,   // calculate genome size using a Raven assembly
    Canu,         // assemble using Canu and clean results
    Flye,         // assemble using Flye
    Lja,          // assemble using LJA
    Metamdbg,     // assemble using metaMDBG
    Miniasm,      // assemble using miniasm and Minipolish
    Myloasm,      // assemble using Myloasm
    Necat,        // assemble using NECAT
    Nextdenovo,   // assemble using NextDenovo and NextPolish
    Plassembler,  // assemble using Plassembler
    Raven,        // assemble using Raven
    Redbean,      // assemble using Redbean (aka wtdbg2)
}

#[derive(ValueEnum, Clone, Debug)]
#[value(rename_all = "snake_case")]
pub enum ReadType {
    OntR9,       // older ONT reads, e.g. R9.4.1
    OntR10,      // newer ONT reads, e.g. R10.4.1
    PacbioClr,   // older PacBio CLR reads
    PacbioHifi,  // newer PacBio HiFi reads
}


fn genome_size_raven(reads: PathBuf, threads: usize, dir: PathBuf, args: String) {
    // TODO
}


fn canu(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
        threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "Canu");
    check_requirements(&["canu"]);
    // TODO
}


fn flye(reads: PathBuf, out_prefix: Option<PathBuf>,
        threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["flye"]);
    // TODO
}


fn lja(reads: PathBuf, out_prefix: Option<PathBuf>,
       threads: usize, dir: PathBuf, args: String) {
    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["lja"]);
    // TODO
}


fn metamdbg(reads: PathBuf, out_prefix: Option<PathBuf>,
            threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["metaMDBG"]);
    // TODO
}


fn miniasm(reads: PathBuf, out_prefix: Option<PathBuf>,
           threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["miniasm", "minipolish", "minimap2", "racon", "any2fasta"]);
    // TODO
}


fn myloasm(reads: PathBuf, out_prefix: Option<PathBuf>,
           threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    // TODO
}


fn necat(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
         threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "NECAT");
    let necat = find_necat();
    // TODO
}


fn nextdenovo(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "NextDenovo");
    check_requirements(&["nextDenovo", "nextPolish"]);
    // TODO
}


fn plassembler(reads: PathBuf, out_prefix: Option<PathBuf>,
               threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["plassembler"]);
    // TODO
}


fn raven(reads: PathBuf, out_prefix: Option<PathBuf>,
         threads: usize, dir: PathBuf, args: String) {
    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["raven"]);
    // TODO
}


fn redbean(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
           threads: usize, dir: PathBuf, read_type: ReadType, args: String) {
    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "Redbean");
    check_requirements(&["wtdbg2", "wtpoa-cns"]);
    // TODO
}


fn check_prefix(out_prefix: Option<PathBuf>) -> PathBuf {
    let prefix = out_prefix.unwrap_or_else(|| {
        quit_with_error("assembly helper commands require --out_prefix")
    });

    let mut fasta = prefix.clone();
    fasta.set_extension("fasta");

    let writable = match OpenOptions::new().write(true).open(&fasta) {
        Ok(_) => true,
        Err(e) if e.kind() == ErrorKind::NotFound => match OpenOptions::new()
            .write(true).create_new(true).open(&fasta)
        {
            Ok(_) => { let _ = remove_file(&fasta); true }
            Err(_) => false,
        },
        _ => false,
    };
    if !writable {
        quit_with_error(&format!("cannot write to this location: {}", fasta.display()));
    }

    prefix
}


fn get_genome_size(genome_size: Option<String>, assembler_name: &str) -> u64 {
    if genome_size.is_none() {
        quit_with_error(&format!("assembly with {assembler_name} requires --genome_size"));
    }
    parse_genome_size(&genome_size.unwrap())
}


fn check_requirements(reqs: &[&str]) {
    for cmd in reqs {
        if which(cmd).is_err() {
            quit_with_error(&format!("required program '{cmd}' not found in $PATH"));
        }
    }
}


fn find_necat() -> PathBuf {
    ["necat", "necat.pl"]
        .into_iter()
        .find_map(|cmd| which(cmd).ok())
        .unwrap_or_else(|| quit_with_error("required program 'necat' (or 'necat.pl') not found in \
                                            $PATH"))
}


fn get_working_dir(dir: Option<PathBuf>) -> PathBuf {
    let dir = dir.unwrap_or_else(|| {
        tempdir().unwrap().into_path()
    });
    create_dir_all(&dir).unwrap_or_else(|e| {
        quit_with_error(&format!("cannot create directory {}: {e}", dir.display()))
    });
    if !dir.is_dir() {
        quit_with_error(&format!("{} exists but is not a directory", dir.display()));
    }
    if NamedTempFile::new_in(&dir).is_err() {
        quit_with_error(&format!("cannot write inside directory {}", dir.display()));
    }
    dir
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;

    #[test]
    fn test_check_prefix_1() {
        // prefix.fasta doesn't exist - should work
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("prefix");
        check_prefix(Some(prefix.clone()));
    }

    #[test]
    fn test_check_prefix_2() {
    // no prefix provided - should fail
        assert!(panic::catch_unwind(|| {
            check_prefix(None);
        }).is_err());
    }

    #[test]
    fn test_check_prefix_3() {
        // can't create nested directories - should fail
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("nonexistent/prefix");
        assert!(panic::catch_unwind(|| {
            check_prefix(Some(prefix.clone()));
        }).is_err());
    }
}
