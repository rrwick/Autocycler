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
use ctrlc::set_handler;
use std::fs::{File, OpenOptions, copy, remove_file, create_dir_all, remove_dir_all};
use std::io::ErrorKind;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Once;
use which::which;
use tempfile::{tempdir, NamedTempFile, TempDir};

use crate::log::bold;
use crate::misc::{check_if_file_exists, quit_with_error, total_fasta_length};
use crate::subsample::parse_genome_size;


pub fn helper(task: Task, reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: Option<PathBuf>, read_type: ReadType, extra_args: Vec<String>) {
    check_if_file_exists(&reads);
    let (dir, _guard) = get_working_dir(dir);

    match task {
        Task::Genomesize => {
            genome_size_raven(reads, threads, dir, extra_args);
        }
        Task::Canu => {
            canu(reads, out_prefix, genome_size, threads, dir, read_type, extra_args);
        }
        Task::Flye => {
            flye(reads, out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Lja => {
            lja(reads, out_prefix, threads, dir, extra_args);
        }
        Task::Metamdbg => {
            metamdbg(reads, out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Miniasm => {
            miniasm(reads, out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Myloasm => {
            myloasm(reads, out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Necat => {
            necat(reads, out_prefix, genome_size, threads, dir, read_type, extra_args);
        }
        Task::Nextdenovo => {
            nextdenovo(reads, out_prefix, genome_size, threads, dir, read_type, extra_args);
        }
        Task::Plassembler => {
            plassembler(reads, out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Raven => {
            raven(reads, out_prefix, threads, extra_args);
        }
        Task::Redbean => {
            redbean(reads, out_prefix, genome_size, threads, dir, read_type, extra_args);
        }
    }
}


fn canu(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
        threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/marbl/canu
    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "Canu");
    check_requirements(&["canu"]);
    // TODO
}


fn flye(reads: PathBuf, out_prefix: Option<PathBuf>,
        threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/mikolmogorov/Flye

    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["flye"]);
    let fasta = out_prefix.with_extension("fasta");
    let gfa = out_prefix.with_extension("gfa");
    let log = out_prefix.with_extension("log");

    let input_flag = match read_type {
        ReadType::OntR9      => "--nano-raw",
        ReadType::OntR10     => "--nano-hq",
        ReadType::PacbioClr  => "--pacbio-raw",
        ReadType::PacbioHifi => "--pacbio-hifi",
    };

    // Build the Flye command
    let mut cmd = Command::new("flye");
    cmd.arg(input_flag).arg(&reads)
       .arg("--threads").arg(threads.to_string())
       .arg("--out-dir").arg(&dir);
    for token in extra_args { cmd.arg(token); }

    // Redirect Flye's stdout and stderr to the terminal
    cmd.stdin(Stdio::null());
    cmd.stdout(Stdio::inherit());
    cmd.stderr(Stdio::inherit());

    run_command(&mut cmd, "flye");

    // Copy the output files
    check_fasta(&dir.join("assembly.fasta"));
    copy_or_die(&dir.join("assembly.fasta"), &fasta);
    copy_or_die(&dir.join("assembly_graph.gfa"), &gfa);
    copy_or_die(&dir.join("flye.log"), &log);
}


fn lja(reads: PathBuf, out_prefix: Option<PathBuf>,
       threads: usize, dir: PathBuf, extra_args: Vec<String>) {
    // https://github.com/AntonBankevich/LJA

    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["lja"]);
    let fasta = out_prefix.with_extension("fasta");
    let gfa = out_prefix.with_extension("gfa");
    let log = out_prefix.with_extension("log");

    // Build the LJA command
    let mut cmd = Command::new("lja");
    cmd.arg("--output-dir").arg(&dir)
       .arg("--reads").arg(&reads)
       .arg("--threads").arg(threads.to_string());
    for token in extra_args { cmd.arg(token); }

    // Redirect LJA's stdout and stderr to the terminal
    cmd.stdin(Stdio::null());
    cmd.stdout(Stdio::inherit());
    cmd.stderr(Stdio::inherit());

    run_command(&mut cmd, "lja");

    check_fasta(&dir.join("assembly.fasta"));

    // Copy the output files
    copy_or_die(&dir.join("assembly.fasta"), &fasta);
    copy_or_die(&dir.join("mdbg.gfa"), &gfa);
    copy_or_die(&dir.join("dbg.log"), &log);
}


fn metamdbg(reads: PathBuf, out_prefix: Option<PathBuf>,
            threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/GaetanBenoitDev/metaMDBG

    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["metaMDBG"]);
    // TODO
}


fn miniasm(reads: PathBuf, out_prefix: Option<PathBuf>,
           threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/lh3/miniasm https://github.com/rrwick/Minipolish

    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["miniasm", "minipolish", "minimap2", "racon", "any2fasta"]);
    // TODO
}


fn myloasm(reads: PathBuf, out_prefix: Option<PathBuf>,
           threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/bluenote-1577/myloasm

    let out_prefix = check_prefix(out_prefix);
    // TODO
}


fn necat(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
         threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/xiaochuanle/NECAT

    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "NECAT");
    let necat = find_necat();
    // TODO
}


fn nextdenovo(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/Nextomics/NextDenovo https://github.com/Nextomics/NextPolish

    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "NextDenovo");
    check_requirements(&["nextDenovo", "nextPolish"]);
    // TODO
}


fn plassembler(reads: PathBuf, out_prefix: Option<PathBuf>,
               threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/gbouras13/plassembler

    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["plassembler"]);
    // TODO
}


fn raven(reads: PathBuf, out_prefix: Option<PathBuf>, threads: usize, extra_args: Vec<String>) {
    // https://github.com/lbcb-sci/raven

    let out_prefix = check_prefix(out_prefix);
    check_requirements(&["raven"]);
    let fasta = out_prefix.with_extension("fasta");
    let gfa = out_prefix.with_extension("gfa");

    // Build the Raven command
    let mut cmd = Command::new("raven");
    cmd.arg("--threads").arg(threads.to_string())
       .arg("--disable-checkpoints")
       .arg("--graphical-fragment-assembly").arg(&gfa)
       .arg(&reads);
    for token in extra_args { cmd.arg(token); }

    // Redirect Raven's stdout to the FASTA file, stderr to the terminal
    let out_fasta = File::create(&fasta).unwrap_or_else(|e| {
        quit_with_error(&format!("cannot create {}: {e}", fasta.display()))
    });
    cmd.stdin(Stdio::null());
    cmd.stdout(Stdio::from(out_fasta));
    cmd.stderr(Stdio::inherit());

    run_command(&mut cmd, "raven");

    check_fasta(&fasta);
}


fn genome_size_raven(reads: PathBuf, threads: usize, dir: PathBuf, extra_args: Vec<String>) {
    check_requirements(&["raven"]);
    let fasta = dir.join("assembly.fasta");

    // Build the Raven command
    let mut cmd = Command::new("raven");
    cmd.arg("--threads").arg(threads.to_string())
       .arg("--disable-checkpoints")
       .arg(&reads);
    for token in extra_args { cmd.arg(token); }

    // Redirect Raven's stdout to the FASTA file, stderr to the terminal
    let out_fasta = File::create(&fasta).unwrap_or_else(|e| {
        quit_with_error(&format!("cannot create {}: {e}", fasta.display()))
    });
    cmd.stdin(Stdio::null());
    cmd.stdout(Stdio::from(out_fasta));
    cmd.stderr(Stdio::inherit());

    run_command(&mut cmd, "raven");

    // Print the genome size to stdout
    check_fasta(&fasta);
    println!("{}", total_fasta_length(&fasta));
}


fn redbean(reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
           threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/ruanjue/wtdbg2

    let out_prefix = check_prefix(out_prefix);
    let genome_size = get_genome_size(genome_size, "Redbean");
    check_requirements(&["wtdbg2", "wtpoa-cns"]);
    let fasta = out_prefix.with_extension("fasta");

    let preset = match read_type {
        ReadType::OntR9      => "preset2",
        ReadType::OntR10     => "preset2",
        ReadType::PacbioClr  => "preset1",
        ReadType::PacbioHifi => "preset4",
    };

    // Build the wtdbg2 command
    let mut cmd = Command::new("wtdbg2");
    cmd.arg("-x").arg(preset)
       .arg("-g").arg(genome_size.to_string())
       .arg("-i").arg(&reads)
       .arg("-t").arg(threads.to_string())
       .arg("-f")
       .arg("-o").arg(&dir.join("dbg"));
    for token in extra_args { cmd.arg(token); }

    // Redirect wtdbg2's stdout and stderr to the terminal
    cmd.stdin(Stdio::null());
    cmd.stdout(Stdio::inherit());
    cmd.stderr(Stdio::inherit());

    run_command(&mut cmd, "wtdbg2");

    // Build the wtpoa-cns command
    let mut cmd = Command::new("wtpoa-cns");
    cmd.arg("-t").arg(threads.to_string())
       .arg("-i").arg(&dir.join("dbg.ctg.lay.gz"))
       .arg("-f")
       .arg("-o").arg(&dir.join("assembly.fasta"));

    // Redirect wtpoa-cns's stdout and stderr to the terminal
    cmd.stdin(Stdio::null());
    cmd.stdout(Stdio::inherit());
    cmd.stderr(Stdio::inherit());

    run_command(&mut cmd, "wtpoa-cns");

    // Copy the output files
    check_fasta(&dir.join("assembly.fasta"));
    copy_or_die(&dir.join("assembly.fasta"), &fasta);
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


fn check_prefix(out_prefix: Option<PathBuf>) -> PathBuf {
    // The output prefix will be used to create the output FASTA file (prefix.fasta) and possibly
    // other files (prefix.gfa, prefix.log). It can contain a directory path, e.g. path/to/prefix.
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
    // Some assemblers require a genome size to be specified, others do not.
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
    // Either 'necat' or 'necat.pl' are acceptable commands for running NECAT.
    ["necat", "necat.pl"]
        .into_iter()
        .find_map(|cmd| which(cmd).ok())
        .unwrap_or_else(|| quit_with_error("required program 'necat' (or 'necat.pl') not found in \
                                            $PATH"))
}


fn get_working_dir(dir: Option<PathBuf>) -> (PathBuf, Option<TempDir>) {
    // If the user provided a directory, use it. Otherwise, create a temporary directory that will
    // be cleaned up when Autocycler helper exits.
    let (path, guard) = match dir {
        Some(p) => (p, None),
        None => {
            let temp_dir = tempdir().expect("cannot create temp dir");
            let p = temp_dir.path().to_path_buf();
            sigint_cleanup(&p);
            (p, Some(temp_dir))
        }
    };
    create_dir_all(&path).unwrap_or_else(|e| {
        quit_with_error(&format!("cannot create directory {}: {e}", path.display()))
    });
    if !path.is_dir() {
        quit_with_error(&format!("{} exists but is not a directory", path.display()));
    }
    if NamedTempFile::new_in(&path).is_err() {
        quit_with_error(&format!("cannot write inside directory {}", path.display()));
    }
    (path, guard)
}


fn sigint_cleanup(dir: &Path) {
    // Ensures that the temporary directory is removed when the user presses Ctrl-C.
    static ONCE: Once = Once::new();
    let dir = dir.to_path_buf();
    ONCE.call_once(|| {
        set_handler(move || {
            let _ = remove_dir_all(&dir);
            std::process::exit(130);
        }).expect("failed to set Ctrl-C handler");
    });
}


fn check_fasta(fasta: &Path) {
    let ok = std::fs::metadata(&fasta).map(|m| m.len() > 0).unwrap_or(false);
    if !ok {
        let _ = remove_file(&fasta);
        quit_with_error("assembly failed or produced empty output");
    }
}


fn copy_or_die(src: &Path, dest: &Path) {
    copy(src, dest).unwrap_or_else(|e| {
        quit_with_error(&format!("failed to copy {} → {}: {e}", src.display(), dest.display()))
    });
}


fn print_command(cmd: &Command) {
    let mut parts = Vec::new();
    parts.push(cmd.get_program().to_string_lossy().into_owned());
    for arg in cmd.get_args() {
        let s = arg.to_string_lossy();
        if s.contains(' ') { parts.push(format!("\"{s}\"")); }
                      else { parts.push(s.into_owned()); }
    }
    eprintln!();
    bold(&format!("{}", parts.join(" ")));
    eprintln!();
}


fn run_command(cmd: &mut Command, name: &str) {
    // Runs a command and checks if it was successful.
    print_command(&cmd);
    let status = cmd.status().unwrap_or_else(|e| {
        quit_with_error(&format!("failed to launch {name}: {e}"))
    });
    if !status.success() {
        quit_with_error(&format!("{name} exited with status {status}"));
    }
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
