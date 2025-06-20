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
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use regex::Regex;
use std::collections::HashMap;
use std::env;
use std::fs::{File, OpenOptions, copy, remove_file, create_dir_all, remove_dir_all, read_dir,
              metadata, rename};
use std::io::{BufRead, BufReader, BufWriter, ErrorKind, Write, copy as io_copy};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Once;
use std::time::SystemTime;
use which::which;
use tempfile::{tempdir, NamedTempFile, TempDir};

use crate::log::{bold, underline};
use crate::misc::{check_if_file_exists, quit_with_error, total_fasta_length, load_fasta,
                  is_file_empty, is_fasta_empty};
use crate::subsample::parse_genome_size;


#[allow(clippy::too_many_arguments)]
pub fn helper(task: Task, reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: Option<PathBuf>, read_type: ReadType,
              min_depth_absolute: Option<f64>, min_depth_relative: Option<f64>,
              extra_args: Vec<String>) {
    check_if_file_exists(&reads);
    let (dir, _guard) = get_working_dir(dir);

    if task == Task::Genomesize {
        genome_size_raven(reads, threads, dir, extra_args);
        return;
    }

    let out_prefix = check_prefix(out_prefix);
    match task {
        Task::Genomesize => unreachable!(),
        Task::Canu => {
            canu(reads, &out_prefix, genome_size, threads, dir, read_type, extra_args);
        }
        Task::Flye => {
            flye(reads, &out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Lja => {
            lja(reads, &out_prefix, threads, dir, extra_args);
        }
        Task::Metamdbg => {
            metamdbg(reads, &out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Miniasm => {
            miniasm(reads, &out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Myloasm => {
            myloasm(reads, &out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Necat => {
            necat(reads, &out_prefix, genome_size, threads, dir, extra_args);
        }
        Task::Nextdenovo => {
            nextdenovo(reads, &out_prefix, genome_size, threads, dir, read_type, extra_args);
        }
        Task::Plassembler => {
            plassembler(reads, &out_prefix, threads, dir, read_type, extra_args);
        }
        Task::Raven => {
            raven(reads, &out_prefix, threads, extra_args);
        }
        Task::Redbean => {
            redbean(reads, &out_prefix, genome_size, threads, dir, read_type, extra_args);
        }
    }

    depth_filter(&out_prefix, &min_depth_absolute, &min_depth_relative);
    delete_fasta_if_empty(&out_prefix);
}


fn canu(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>, threads: usize,
        dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/marbl/canu

    let genome_size = get_genome_size(genome_size, "Canu");
    check_requirements(&["canu"]);

    let input_flag = match read_type {
        ReadType::OntR9      => "-nanopore",
        ReadType::OntR10     => "-nanopore",
        ReadType::PacbioClr  => "-pacbio",
        ReadType::PacbioHifi => "-pacbio-hifi",
    };

    let mut cmd = Command::new("canu");
    cmd.arg("-p").arg("canu")
       .arg("-d").arg(&dir)
       .arg("-fast")
       .arg(format!("genomeSize={genome_size}"))
       .arg("useGrid=false")
       .arg(format!("maxThreads={threads}"))
       .arg(input_flag).arg(&reads);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_canu_fasta(&dir.join("canu.contigs.fasta"), &dir.join("canu.contigs.layout.tigInfo"),
                    &out_prefix.with_extension("fasta"));
    copy_output_file(&dir.join("canu.report"), &out_prefix.with_extension("log"));
}


fn flye(reads: PathBuf, out_prefix: &Path, threads: usize, dir: PathBuf, read_type: ReadType,
        extra_args: Vec<String>) {
    // https://github.com/mikolmogorov/Flye

    check_requirements(&["flye"]);

    let input_flag = match read_type {
        ReadType::OntR9      => "--nano-raw",
        ReadType::OntR10     => "--nano-hq",
        ReadType::PacbioClr  => "--pacbio-raw",
        ReadType::PacbioHifi => "--pacbio-hifi",
    };

    let mut cmd = Command::new("flye");
    cmd.arg(input_flag).arg(&reads)
       .arg("--threads").arg(threads.to_string())
       .arg("--out-dir").arg(&dir);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_flye_fasta(&dir.join("assembly.fasta"), &dir.join("assembly_info.txt"),
                    &out_prefix.with_extension("fasta"));
    copy_output_file(&dir.join("assembly_graph.gfa"), &out_prefix.with_extension("gfa"));
    copy_output_file(&dir.join("flye.log"), &out_prefix.with_extension("log"));
}


fn lja(reads: PathBuf, out_prefix: &Path, threads: usize, dir: PathBuf, extra_args: Vec<String>) {
    // https://github.com/AntonBankevich/LJA

    check_requirements(&["lja"]);

    let mut cmd = Command::new("lja");
    cmd.arg("--output-dir").arg(&dir)
       .arg("--reads").arg(&reads)
       .arg("--threads").arg(threads.to_string());
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_fasta(&dir.join("assembly.fasta"), &out_prefix.with_extension("fasta"));
    copy_output_file(&dir.join("mdbg.gfa"), &out_prefix.with_extension("gfa"));
    copy_output_file(&dir.join("dbg.log"), &out_prefix.with_extension("log"));
}


fn metamdbg(reads: PathBuf, out_prefix: &Path, threads: usize, dir: PathBuf, read_type: ReadType,
            extra_args: Vec<String>) {
    // https://github.com/GaetanBenoitDev/metaMDBG

    check_requirements(&["metaMDBG"]);

    let input_flag = match read_type {
        ReadType::OntR9      => "--in-ont",
        ReadType::OntR10     => "--in-ont",
        ReadType::PacbioClr  => "--in-ont",
        ReadType::PacbioHifi => "--in-hifi",
    };

    let mut cmd = Command::new("metaMDBG");
    cmd.arg("asm")
       .arg("--out-dir").arg(&dir)
       .arg(input_flag).arg(&reads)
       .arg("--threads").arg(threads.to_string());
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_fasta(&dir.join("contigs.fasta.gz"), &out_prefix.with_extension("fasta"));
    copy_output_file(&dir.join("metaMDBG.log"), &out_prefix.with_extension("log"));
}


fn miniasm(reads: PathBuf, out_prefix: &Path, threads: usize, dir: PathBuf, read_type: ReadType,
           extra_args: Vec<String>) {
    // https://github.com/lh3/miniasm
    // https://github.com/rrwick/Minipolish

    check_requirements(&["miniasm", "minipolish", "minimap2", "racon"]);

    let ava_arg = match read_type {
        ReadType::OntR9      => "ava-ont",
        ReadType::OntR10     => "-k19 -Xw7 -e0 -m100",
        ReadType::PacbioClr  => "ava-pb",
        ReadType::PacbioHifi => "-k23 -Xw11 -e0 -m100",
    };
    let map_preset = match read_type {
        ReadType::OntR9      => "map-ont",
        ReadType::OntR10     => "lr:hq",
        ReadType::PacbioClr  => "map-pb",
        ReadType::PacbioHifi => "map-hifi",
    };

    let mut cmd = Command::new("minimap2");
    cmd.arg("-t").arg(threads.to_string());
    if ava_arg.starts_with('-') { cmd.args(ava_arg.split_whitespace()); }
                           else { cmd.arg("-x").arg(ava_arg); }
    cmd.arg(&reads).arg(&reads);
    redirect_stderr_and_stdout(&mut cmd, Some(&dir.join("overlap.paf")));
    run_command(&mut cmd);

    let mut cmd = Command::new("miniasm");
    cmd.arg("-f").arg(&reads)
       .arg(dir.join("overlap.paf"));
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, Some(&dir.join("unpolished.gfa")));
    run_command(&mut cmd);

    let mut cmd = Command::new("minipolish");
    cmd.arg("--threads").arg(threads.to_string())
       .arg("--minimap2-preset").arg(map_preset)
       .arg(&reads)
       .arg(dir.join("unpolished.gfa"));
    redirect_stderr_and_stdout(&mut cmd, Some(&out_prefix.with_extension("gfa")));
    run_command(&mut cmd);

    minipolish_gfa_to_fasta(&out_prefix.with_extension("gfa"), &out_prefix.with_extension("fasta"));
}


fn myloasm(reads: PathBuf, out_prefix: &Path, threads: usize, dir: PathBuf, read_type: ReadType,
           extra_args: Vec<String>) {
    // https://github.com/bluenote-1577/myloasm

    check_requirements(&["myloasm"]);

    let mut cmd = Command::new("myloasm");
    cmd.arg("--output-dir").arg(&dir)
       .arg(&reads)
       .arg("--threads").arg(threads.to_string());
    if read_type == ReadType::PacbioHifi { cmd.arg("--hifi"); }
    else if read_type == ReadType::OntR10 { cmd.arg("--nano-r10"); }
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_fasta(&dir.join("assembly_primary.fa"), &out_prefix.with_extension("fasta"));
    replace_underscores_with_spaces(&out_prefix.with_extension("fasta"));
    copy_output_file(&dir.join("final_contig_graph.gfa"), &out_prefix.with_extension("gfa"));
    copy_output_file(&find_log_file(&dir, "myloasm"), &out_prefix.with_extension("log"));
}


fn necat(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>, threads: usize,
         dir: PathBuf, extra_args: Vec<String>) {
    // https://github.com/xiaochuanle/NECAT

    let genome_size = get_genome_size(genome_size, "NECAT");
    make_necat_files(&reads, &dir, genome_size, threads);

    let mut cmd = Command::new(find_necat());
    cmd.arg("bridge").arg("config.txt");
    for token in extra_args { cmd.arg(token); }
    cmd.current_dir(&dir);
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_fasta(&dir.join("necat/6-bridge_contigs/polished_contigs.fasta"),
               &out_prefix.with_extension("fasta"));
}


fn nextdenovo(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>, threads: usize,
              dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/Nextomics/NextDenovo
    // https://github.com/Nextomics/NextPolish

    let genome_size = get_genome_size(genome_size, "NextDenovo");
    check_requirements(&["nextDenovo", "nextPolish"]);
    make_nextdenovo_files(&dir, &reads, genome_size, threads, read_type);

    let mut cmd = Command::new("nextDenovo");
    cmd.arg("nextdenovo_run.cfg");
    for token in extra_args { cmd.arg(token); }
    cmd.current_dir(&dir);
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    let mut cmd = Command::new("nextPolish");
    cmd.arg("nextpolish_run.cfg");
    cmd.current_dir(&dir);
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_fasta(&dir.join("nextpolish/genome.nextpolish.fasta"),
               &out_prefix.with_extension("fasta"));
    combine_nextdenovo_logs(&dir, &out_prefix.with_extension("log"));
}


fn plassembler(reads: PathBuf, out_prefix: &Path, threads: usize, dir: PathBuf, read_type: ReadType,
               extra_args: Vec<String>) {
    // https://github.com/gbouras13/plassembler

    check_requirements(&["plassembler", "chopper", "dnaapler", "fastp", "mash", "minimap2",
                         "raven", "samtools", "unicycler"]);
    let db = find_plassembler_db();

    let mut cmd = Command::new("plassembler");
    cmd.arg("long")
       .arg("-d").arg(&db)
       .arg("-l").arg(&reads)
       .arg("-o").arg(&dir)
       .arg("-t").arg(threads.to_string())
       .arg("--force")
       .arg("--skip_qc");
    if read_type == ReadType::OntR9      { cmd.arg("--raw_flag"); }
    if read_type == ReadType::PacbioClr  { cmd.arg("--pacbio_model").arg("pacbio-raw"); }
    if read_type == ReadType::PacbioHifi { cmd.arg("--pacbio_model").arg("pacbio-hifi"); }
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_output_file(&dir.join("plassembler_plasmids.gfa"), &out_prefix.with_extension("gfa"));
    rotate_plassembler_contigs(&dir.join("plassembler_plasmids.fasta"),
                               &out_prefix.with_extension("fasta"));
    copy_output_file(&find_log_file(&dir, "plassembler"), &out_prefix.with_extension("log"));
}


fn raven(reads: PathBuf, out_prefix: &Path, threads: usize, extra_args: Vec<String>) {
    // https://github.com/lbcb-sci/raven

    check_requirements(&["raven"]);

    let mut cmd = Command::new("raven");
    cmd.arg("--threads").arg(threads.to_string())
       .arg("--disable-checkpoints")
       .arg("--graphical-fragment-assembly").arg(out_prefix.with_extension("gfa"))
       .arg(&reads);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, Some(&out_prefix.with_extension("fasta")));
    run_command(&mut cmd);
}


fn genome_size_raven(reads: PathBuf, threads: usize, dir: PathBuf, extra_args: Vec<String>) {
    check_requirements(&["raven"]);

    let mut cmd = Command::new("raven");
    cmd.arg("--threads").arg(threads.to_string())
       .arg("--disable-checkpoints")
       .arg(&reads);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, Some(&dir.join("assembly.fasta")));
    run_command(&mut cmd);

    if is_fasta_empty(&dir.join("assembly.fasta")) {
        quit_with_error("Raven assembly failed");
    }
    println!("{}", total_fasta_length(&dir.join("assembly.fasta")));
}


fn redbean(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>, threads: usize,
           dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/ruanjue/wtdbg2

    let genome_size = get_genome_size(genome_size, "Redbean");
    check_requirements(&["wtdbg2", "wtpoa-cns"]);

    let preset = match read_type {
        ReadType::OntR9      => "preset2",
        ReadType::OntR10     => "preset2",
        ReadType::PacbioClr  => "preset1",
        ReadType::PacbioHifi => "preset4",
    };

    let mut cmd = Command::new("wtdbg2");
    cmd.arg("-x").arg(preset)
       .arg("-g").arg(genome_size.to_string())
       .arg("-i").arg(&reads)
       .arg("-t").arg(threads.to_string())
       .arg("-f")
       .arg("-o").arg(dir.join("dbg"));
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    let mut cmd = Command::new("wtpoa-cns");
    cmd.arg("-t").arg(threads.to_string())
       .arg("-i").arg(dir.join("dbg.ctg.lay.gz"))
       .arg("-f")
       .arg("-o").arg(dir.join("assembly.fasta"));
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    copy_fasta(&dir.join("assembly.fasta"), &out_prefix.with_extension("fasta"));
}


#[derive(ValueEnum, Clone, Debug, PartialEq)]
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

#[derive(ValueEnum, Clone, Debug, PartialEq)]
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

    // Make parent directories if they don't exist.
    if let Some(parent) = prefix.parent() {
        if let Err(e) = std::fs::create_dir_all(parent) {
            quit_with_error(&format!("cannot create directory {}: {e}", parent.display()));
        }
    }

    // Ensure that we can create/overwrite prefix.fasta
    let fasta = prefix.clone().with_extension("fasta");
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
    if !writable { quit_with_error(&format!("cannot write to {}", fasta.display())); }

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


fn copy_output_file(src: &Path, dest: &Path) {
    if !src.exists() || is_file_empty(src) { let _ = remove_file(dest); return; }
    copy(src, dest).unwrap_or_else(|e| {
        quit_with_error(&format!("failed to copy {} → {}: {e}", src.display(), dest.display()))
    });
}


fn copy_fasta(src: &Path, dest: &Path) {
    // Can handle gzipped and uncompressed FASTA files - destination will always be uncompressed.
    // Multi-line-per-seq FASTA files are converted to single-line-per-seq.
    if !src.exists() || is_fasta_empty(src) { let _ = remove_file(dest); return; }
    let mut writer = BufWriter::new(File::create(dest).unwrap());
    for (_, header, seq) in load_fasta(src) {
        writeln!(writer, ">{header}\n{seq}").unwrap();
    }
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
    bold(&parts.join(" "));
    eprintln!();
}


fn run_command(cmd: &mut Command) {
    let name = cmd.get_program().to_string_lossy().into_owned();
    print_command(cmd);
    let status = cmd.status().unwrap_or_else(|e| {
        quit_with_error(&format!("failed to launch {name}: {e}"))
    });
    if !status.success() {
        eprintln!("{name} failed with status {status}");
    }
}


fn redirect_stderr_and_stdout(cmd: &mut Command, stdout_file: Option<&Path>) {
    // Redirects the command's:
    // * stderr to the terminal
    // * stdout to a file if stdout_file is provided, otherwise to the terminal
    cmd.stdin(Stdio::null());
    if let Some(file) = stdout_file {
        let out_file = File::create(file).unwrap_or_else(|e| {
            quit_with_error(&format!("cannot create {}: {e}", file.display()))
        });
        cmd.stdout(Stdio::from(out_file));
    } else {
        cmd.stdout(Stdio::inherit());
    }
    cmd.stderr(Stdio::inherit());
}


fn find_log_file(dir: &Path, prefix: &str) -> PathBuf {
    read_dir(dir).unwrap().filter_map(|e| e.ok().map(|e| e.path()))
      .find(|p| p.file_name().and_then(|s| s.to_str())
                 .is_some_and(|n| n.starts_with(prefix) && n.ends_with(".log")))
      .unwrap_or_else(|| quit_with_error(&format!("{prefix} log file not found")))
}


fn minipolish_gfa_to_fasta(gfa: &Path, fasta: &Path) {
    if !gfa.exists() || is_file_empty(gfa) { return; }
    let reader = BufReader::new(File::open(gfa).unwrap());
    let mut writer = BufWriter::new(File::create(fasta).unwrap());
    for line in reader.lines().map_while(Result::ok) {
        if !line.starts_with('S') { continue; }
        let mut cols = line.split('\t');
        cols.next();
        let (name, seq) = (cols.next().unwrap_or(""), cols.next().unwrap_or(""));
        let depth = cols.find_map(|field| field.strip_prefix("dp:f:"));
        let mut header = format!(">{name}");
        if name.ends_with('c') { header.push_str(" circular=true"); }
        if let Some(d) = depth { header.push_str(&format!(" depth={d}")); }
        writeln!(writer, "{header}\n{seq}").unwrap();
    }
}


fn copy_flye_fasta(src: &Path, assembly_info: &Path, dest: &Path) {
    if !src.exists() || is_fasta_empty(src) { return; }
    let mut writer = BufWriter::new(File::create(dest).unwrap());
    let info = load_flye_assembly_info(assembly_info);
    for (name, _, seq) in load_fasta(src) {
        let mut header = name.to_string();
        if let Some((is_circ, depth)) = info.get(&name) {
            if *is_circ { header.push_str(" circular=true"); }
            header.push_str(&format!(" depth={depth}"));
        }
        writeln!(writer, ">{header}\n{seq}").unwrap();
    }
}


fn load_flye_assembly_info(assembly_info: &Path) -> HashMap<String, (bool, String)> {
    // Loads Flye's assembly_info.txt file, returns a map of contig names to circularity and depth.
    let mut info: HashMap<String, (bool, String)> = HashMap::new();
    for line in BufReader::new(File::open(assembly_info).unwrap()).lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 { continue; }
        let name  = cols[0].to_string();
        let depth = cols[2].to_string();
        let circ  = cols[3] == "Y";
        info.insert(name, (circ, depth));
    }
    info
}


fn copy_canu_fasta(src: &Path, assembly_info: &Path, dest: &Path) {
    // Copies the Canu assembly FASTA file, with the following modifications:
    // * Adds the depth from the assembly_info file to the header.
    // * Excludes repeat/bubble contigs.
    // * Trims overlap from circular contigs.
    if !src.exists() || is_fasta_empty(src) { return; }
    let mut writer = BufWriter::new(File::create(dest).unwrap());
    let depth = load_canu_assembly_depth(assembly_info);
    for (name, header, seq) in load_fasta(src) {
        if header.contains("suggestRepeat=yes") || header.contains("suggestBubble=yes") {
            continue;
        }
        let (mut header, seq) = trim_canu_contig(header, seq);
        if let Some(d) = depth.get(&name) {
            header.push_str(&format!(" depth={d}"));
        }
        writeln!(writer, ">{header}\n{seq}").unwrap();
    }
}


fn load_canu_assembly_depth(assembly_info: &Path) -> HashMap<String, String> {
    // Loads Canu's *.contigs.layout.tigInfo file, returns a map of contig names to depth.
    let mut info: HashMap<String, String> = HashMap::new();
    for line in BufReader::new(File::open(assembly_info).unwrap()).lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 { continue; }
        let id: u64 = match cols[0].parse() { Ok(n) => n, Err(_) => continue };
        let name = format!("tig{:08}", id);
        let depth = cols[2].to_string();
        info.insert(name, depth);
    }
    info
}


fn trim_canu_contig(mut header: String, mut seq: String) -> (String, String) {
    if !header.contains("suggestCircular=yes") {
        return (header, seq);
    }
    let re_trim = Regex::new(r"trim=(\d+)-(\d+)").unwrap();
    let re_len = Regex::new(r"len=\d+").unwrap();
    if let Some(cap) = re_trim.captures(&header) {
        let start: usize = cap[1].parse().unwrap_or(0);
        let end: usize = cap[2].parse().unwrap_or(seq.len());
        if start < end && end <= seq.len() {
            seq = seq[start..end].to_owned();
            header = re_trim.replace(&header, format!("trim=0-{}", seq.len())).into_owned();
            header = re_len.replace(&header, format!("len={}", seq.len())).into_owned();
        }
    }
    (header, seq)
}


fn make_necat_files(reads: &Path, dir: &Path, genome_size: u64, threads: usize) {
    let mut r = BufWriter::new(File::create(dir.join("read_list.txt")).unwrap());
    writeln!(r, "{}", reads.canonicalize().unwrap_or_else(|_| reads.to_path_buf())
                           .display()).unwrap();

    let mut w = BufWriter::new(File::create(dir.join("config.txt")).unwrap());
    writeln!(w, "PROJECT=necat").unwrap();
    writeln!(w, "ONT_READ_LIST=read_list.txt").unwrap();
    writeln!(w, "GENOME_SIZE={}", genome_size).unwrap();
    writeln!(w, "THREADS={threads}").unwrap();
    writeln!(w, "MIN_READ_LENGTH=3000").unwrap();
    writeln!(w, "PREP_OUTPUT_COVERAGE=40").unwrap();
    writeln!(w, "OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000").unwrap();
    writeln!(w, "OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000").unwrap();
    writeln!(w, "CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0").unwrap();
    writeln!(w, "CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0").unwrap();
    writeln!(w, "TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400").unwrap();
    writeln!(w, "ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400").unwrap();
    writeln!(w, "NUM_ITER=2").unwrap();
    writeln!(w, "CNS_OUTPUT_COVERAGE=30").unwrap();
    writeln!(w, "CLEANUP=1").unwrap();
    writeln!(w, "USE_GRID=false").unwrap();
    writeln!(w, "GRID_NODE=0").unwrap();
    writeln!(w, "GRID_OPTIONS=").unwrap();
    writeln!(w, "SMALL_MEMORY=0").unwrap();
    writeln!(w, "FSA_OL_FILTER_OPTIONS=").unwrap();
    writeln!(w, "FSA_ASSEMBLE_OPTIONS=").unwrap();
    writeln!(w, "FSA_CTG_BRIDGE_OPTIONS=").unwrap();
    writeln!(w, "POLISH_CONTIGS=true").unwrap();
}


fn make_nextdenovo_files(dir: &Path, reads: &Path, genome_size: u64, threads: usize,
                         read_type: ReadType) {
    let (lgs_or_hifi, ont_clr_or_hifi, map_preset) = match read_type {
            ReadType::OntR9      => ("lgs",  "ont",  "map-ont"),
            ReadType::OntR10     => ("lgs",  "ont",  "map-ont"),  // lr:hq breaks NextPolish
            ReadType::PacbioClr  => ("lgs",  "clr",  "map-pb"),
            ReadType::PacbioHifi => ("hifi", "hifi", "map-hifi"),
    };

    let mut r = BufWriter::new(File::create(dir.join("input.fofn")).unwrap());
    writeln!(r, "{}", reads.canonicalize().unwrap_or_else(|_| reads.to_path_buf())
                           .display()).unwrap();

    let mut c1 = BufWriter::new(File::create(dir.join("nextdenovo_run.cfg")).unwrap());
    writeln!(c1, "[General]").unwrap();
    writeln!(c1, "job_type = local\njob_prefix = nextDenovo\ntask = all").unwrap();
    writeln!(c1, "rewrite = yes\ndeltmp = yes\nparallel_jobs = 1\ninput_type = raw").unwrap();
    writeln!(c1, "read_type = {ont_clr_or_hifi}").unwrap();
    writeln!(c1, "input_fofn = input.fofn\nworkdir = nextdenovo").unwrap();
    writeln!(c1).unwrap();
    writeln!(c1, "[correct_option]").unwrap();
    writeln!(c1, "read_cutoff = 1k").unwrap();
    writeln!(c1, "genome_size = {genome_size}").unwrap();
    writeln!(c1, "sort_options = -m 20g -t {threads}").unwrap();
    writeln!(c1, "minimap2_options_raw = -t {threads}").unwrap();
    writeln!(c1, "pa_correction = 1").unwrap();
    writeln!(c1, "correction_options = -p {threads}").unwrap();
    writeln!(c1).unwrap();
    writeln!(c1, "[assemble_option]").unwrap();
    writeln!(c1, "minimap2_options_cns = -t {threads}").unwrap();
    writeln!(c1, "nextgraph_options = -a 1").unwrap();

    let mut c2 = BufWriter::new(File::create(dir.join("nextpolish_run.cfg")).unwrap());
    writeln!(c2, "[General]").unwrap();
    writeln!(c2, "job_type = local\njob_prefix = nextPolish\ntask = best").unwrap();
    writeln!(c2, "rewrite = yes\ndeltmp = yes\nrerun = 3\nparallel_jobs = 1").unwrap();
    writeln!(c2, "multithread_jobs = {threads}").unwrap();
    writeln!(c2, "genome = nextdenovo/03.ctg_graph/nd.asm.fasta").unwrap();
    writeln!(c2, "genome_size = auto\nworkdir = nextpolish").unwrap();
    writeln!(c2, "polish_options = -p {threads}").unwrap();
    writeln!(c2).unwrap();
    writeln!(c2, "[{lgs_or_hifi}_option]").unwrap();
    writeln!(c2, "{lgs_or_hifi}_fofn = input.fofn").unwrap();
    writeln!(c2, "{lgs_or_hifi}_options = -min_read_len 1k -max_depth 100").unwrap();
    writeln!(c2, "{lgs_or_hifi}_minimap2_options = -x {map_preset} -t {threads}").unwrap();
}


fn combine_nextdenovo_logs(dir: &Path, dest: &Path) {
    // Combines the two pid*.log.info files from NextDenovo and NextPolish into a single file.
    let mut logs: Vec<_> = read_dir(dir).unwrap().filter_map(|e| {
        let p = e.ok()?.path();
        let is_log = p.file_name().and_then(|s| s.to_str())
            .map(|name| name.starts_with("pid") && name.ends_with(".log.info")).unwrap_or(false);
        if is_log { Some(p) } else { None }
    }).collect();
    if logs.is_empty() { return; }
    logs.sort_by_key(|p| {
        metadata(p).and_then(|m| m.modified()).unwrap_or(SystemTime::UNIX_EPOCH)
    });
    let mut out = BufWriter::new(File::create(dest).unwrap());
    for log in logs {
        let mut rdr = BufReader::new(File::open(&log).unwrap());
        io_copy(&mut rdr, &mut out).unwrap();
    }
}


fn find_plassembler_db() -> PathBuf {
    if let Ok(p) = env::var("PLASSEMBLER_DB") {
        let pb = PathBuf::from(&p);
        if pb.is_dir() { return pb; }
    }
    if let Ok(p) = env::var("CONDA_PREFIX") {
        let pb = Path::new(&p).join("plassembler_db");
        if pb.is_dir() { return pb; }
    }
    quit_with_error("No Plassembler database found. \
                     Set PLASSEMBLER_DB or ensure $CONDA_PREFIX/plassembler_db exists.")
}


fn rotate_plassembler_contigs(src: &Path, dest: &Path) {
    if !src.exists() || is_fasta_empty(src) { return; }
    let mut w = BufWriter::new(File::create(dest).unwrap());
    let mut rng = StdRng::seed_from_u64(0);
    for (_name, header, seq) in load_fasta(src) {
        if header.to_ascii_lowercase().contains("circular=true") && seq.len() > 1 {
            let r = rng.gen_range(1..seq.len());
            let rotated = seq[r..].to_owned() + &seq[..r];
            writeln!(w, ">{header}\n{rotated}").unwrap();
        } else {
            writeln!(w, ">{header}\n{seq}").unwrap();
        }
    }
}


fn replace_underscores_with_spaces(filename: &Path) {
    if !filename.exists() || is_file_empty(filename) { return; }
    let in_file = BufReader::new(File::open(filename).unwrap());
    let tmp_path = filename.with_extension("tmp");
    let mut out_file = BufWriter::new(File::create(&tmp_path).unwrap());
    for line in in_file.lines().map_while(Result::ok) {
        writeln!(out_file, "{}", line.replace('_', " ")).unwrap();
    }
    rename(tmp_path, filename).unwrap();
}


fn depth_filter(out_prefix: &Path, min_depth_absolute: &Option<f64>,
                min_depth_relative: &Option<f64>) {
    // Filters the final FASTA file by depth, overwriting the original file. If any contig does
    // not have depth, the function does nothing.
    if min_depth_absolute.is_none() && min_depth_relative.is_none() { return; }
    let fasta = out_prefix.with_extension("fasta");
    if !fasta.exists() || is_fasta_empty(&fasta) { return; }

    let mut records = Vec::new();
    let (mut longest_len, mut longest_depth) = (0usize, 0.0);
    for (name, header, seq) in load_fasta(&fasta) {
        let depth = match depth_from_header(&header) { Some(d) => d, None => return };
        let len = seq.len();
        if len > longest_len { longest_len = len; longest_depth = depth; }
        records.push((name, header, seq, depth));
    }

    let mut threshold = min_depth_absolute.unwrap_or(0.0);
    if let Some(r) = min_depth_relative { threshold = threshold.max(r * longest_depth); }
    eprintln!();
    underline("Autocycler helper depth filter");
    eprintln!("threshold = {:.3}", threshold);

    let kept: Vec<_> = records.into_iter().filter_map(|(name, header, seq, depth)| {
        let pass = depth >= threshold;
        eprintln!("{name}: depth={:.3}, {}", depth, if pass { "PASS" } else { "FAIL" });
        if pass { Some((header, seq)) } else { None }
    }).collect();

    if kept.is_empty() { let _ = remove_file(&fasta); return; }
    let mut w = BufWriter::new(File::create(&fasta).unwrap());
    for (header, seq) in kept { writeln!(w, ">{header}\n{seq}").unwrap(); }
}


fn depth_from_header(header: &str) -> Option<f64> {
    fn parse_num(s: &str) -> Option<f64> {
        s.split(['-', '_', ' ']).next()?.parse().ok()
    }
    if let Some(i) = header.find("depth=")    { return parse_num(&header[i + 6..]); }
    if let Some(i) = header.find("depth-")    { return parse_num(&header[i + 6..]); }
    if let Some(i) = header.find("coverage=") { return parse_num(&header[i + 9..]); }
    None
}


fn delete_fasta_if_empty(out_prefix: &Path) {
    let fasta = out_prefix.with_extension("fasta");
    if fasta.exists() && is_fasta_empty(&fasta) {
        let _ = std::fs::remove_file(&fasta);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;
    use crate::tests::{make_test_file, make_gzipped_test_file};

    #[test]
    fn test_check_prefix() {
        // prefix.fasta doesn't exist - should work
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("prefix");
        check_prefix(Some(prefix.clone()));

        // nested directories should work
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("abc/def/prefix");
        check_prefix(Some(prefix.clone()));

        // no prefix provided - should fail
        assert!(panic::catch_unwind(|| {
            check_prefix(None);
        }).is_err());
    }

    #[test]
    fn test_depth_from_header() {
        assert_eq!(depth_from_header(">contig depth=10.5"), Some(10.5));
        assert_eq!(depth_from_header(">contig circular=true depth=5.0"), Some(5.0));
        assert_eq!(depth_from_header(">contig"), None);
        assert_eq!(depth_from_header(">a_len-12_circular-no_depth-37-37-37_mult-2.00"), Some(37.0));
        assert_eq!(depth_from_header(">b_len-9_circular-yes_depth-25-24-23_mult-1.00"), Some(25.0));
        assert_eq!(depth_from_header(">a len-12 circular-no depth-37-37-37 mult-2.00"), Some(37.0));
        assert_eq!(depth_from_header(">b len-9 circular-yes depth-25-24-23 mult-1.00"), Some(25.0));
        assert_eq!(depth_from_header(">ctg15 length=123 coverage=49.70 circular=yes"), Some(49.7));
    }

    #[test]
    fn test_depth_filter() {
        let dir = tempdir().unwrap();
        let out_prefix = dir.path().join("test");
        let fasta = out_prefix.with_extension("fasta");

        make_test_file(&fasta, ">a depth=20\nACGT\n\
                                >b depth=120\nCGA\n\
                                >c depth=200\nACAGACTACGACTACGACGACGATCAGCGACATCGACGT\n\
                                >d depth=100\nCGATCGACTACC\n");

        depth_filter(&out_prefix, &None, &None);
        assert_eq!(load_fasta(&fasta).len(), 4);

        depth_filter(&out_prefix, &None, &Some(0.09));
        assert_eq!(load_fasta(&fasta).len(), 4);

        depth_filter(&out_prefix, &None, &Some(0.11));
        assert_eq!(load_fasta(&fasta).len(), 3);

        depth_filter(&out_prefix, &Some(99.0), &None);
        assert_eq!(load_fasta(&fasta).len(), 3);

        depth_filter(&out_prefix, &Some(101.0), &None);
        assert_eq!(load_fasta(&fasta).len(), 2);

        depth_filter(&out_prefix, &None, &Some(0.61));
        assert_eq!(load_fasta(&fasta).len(), 1);

        depth_filter(&out_prefix, &Some(201.0), &None);
        assert!(panic::catch_unwind(|| {
            load_fasta(&fasta).len();
        }).is_err());
    }

    #[test]
    fn test_trim_canu_contig_1() {
        let header = ">tig00000001 len=60 reads=50 class=contig suggestRepeat=no \
                      suggestBubble=no suggestCircular=no trim=0-60".to_string();
        let seq = "AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC".to_string();
        let (new_header, new_seq) = trim_canu_contig(header.clone(), seq.clone());
        assert_eq!(new_header, header);
        assert_eq!(new_seq, seq);
    }

    #[test]
    fn test_trim_canu_contig_2() {
        let header = ">tig00000001 len=60 reads=50 class=contig suggestRepeat=no \
                      suggestBubble=no suggestCircular=yes trim=0-50".to_string();
        let seq = "AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC".to_string();
        let (new_header, new_seq) = trim_canu_contig(header, seq);
        assert_eq!(new_header, ">tig00000001 len=50 reads=50 class=contig suggestRepeat=no \
                                suggestBubble=no suggestCircular=yes trim=0-50");
        assert_eq!(new_seq, "AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCA");
    }

    #[test]
    fn test_trim_canu_contig_3() {
        let header = ">tig00000001 len=60 reads=50 class=contig suggestRepeat=no \
                      suggestBubble=no suggestCircular=yes trim=10-60".to_string();
        let seq = "AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC".to_string();
        let (new_header, new_seq) = trim_canu_contig(header, seq);
        assert_eq!(new_header, ">tig00000001 len=50 reads=50 class=contig suggestRepeat=no \
                                suggestBubble=no suggestCircular=yes trim=0-50");
        assert_eq!(new_seq, "CTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC");
    }

    #[test]
    fn test_trim_canu_contig_4() {
        let header = ">tig00000001 len=60 reads=50 class=contig suggestRepeat=no \
                      suggestBubble=no suggestCircular=yes trim=10-50".to_string();
        let seq = "AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC".to_string();
        let (new_header, new_seq) = trim_canu_contig(header, seq);
        assert_eq!(new_header, ">tig00000001 len=40 reads=50 class=contig suggestRepeat=no \
                                suggestBubble=no suggestCircular=yes trim=0-40");
        assert_eq!(new_seq, "CTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCA");
    }

    #[test]
    fn test_rotate_plassembler_contigs_1() {
        // non-circular contigs should not be modified
        let dir = tempdir().unwrap();
        let in_fasta = dir.path().join("input.fasta");
        let out_fasta = dir.path().join("output.fasta");
        make_test_file(&in_fasta, ">a\nACGATCGCT\n\
                                   >b\nCGATCGACTAC\n");
        rotate_plassembler_contigs(&in_fasta, &out_fasta);
        let in_seqs: Vec<String> = load_fasta(&in_fasta).into_iter().map(|(_, _, s)| s).collect();
        let out_seqs: Vec<String> = load_fasta(&out_fasta).into_iter().map(|(_, _, s)| s).collect();
        assert_eq!(in_seqs, out_seqs);
    }

    #[test]
    fn test_rotate_plassembler_contigs_2() {
        // circular contigs should be modified
        let dir = tempdir().unwrap();
        let in_fasta = dir.path().join("input.fasta");
        let out_fasta = dir.path().join("output.fasta");
        make_test_file(&in_fasta, ">a circular=True\nACGATCGCT\n\
                                   >b circular=True\nCGATCGACTAC\n");
        rotate_plassembler_contigs(&in_fasta, &out_fasta);
        let in_seqs: Vec<String> = load_fasta(&in_fasta).into_iter().map(|(_, _, s)| s).collect();
        let out_seqs: Vec<String> = load_fasta(&out_fasta).into_iter().map(|(_, _, s)| s).collect();
        assert_ne!(in_seqs, out_seqs);
    }

    #[test]
    fn test_replace_underscores_with_spaces() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("test.fasta");
        make_test_file(&filename, ">a_len-12_circular-no_depth-37-37-37_mult-2.00\nACGATCGCT\n\
                                   >b_len-9_circular-yes_depth-25-24-23_mult-1.00\nCGATCGACTAC\n");
        replace_underscores_with_spaces(&filename);
        let expected = ">a len-12 circular-no depth-37-37-37 mult-2.00\nACGATCGCT\n\
                        >b len-9 circular-yes depth-25-24-23 mult-1.00\nCGATCGACTAC\n";
        assert_eq!(std::fs::read_to_string(&filename).unwrap(), expected);
    }

    #[test]
    fn test_copy_fasta() {
        let dir = tempdir().unwrap();
        let in_fasta = dir.path().join("in.fasta");
        let out_fasta = dir.path().join("out.fasta");

        // Empty FASTA files are not copied
        make_test_file(&in_fasta, "");
        copy_fasta(&in_fasta, &out_fasta);
        assert!(!out_fasta.exists());

        // Multi-line-per-seq FASTA file is converted to single-line-per-seq FASTA
        make_test_file(&in_fasta, ">a\nACGA\nTCGC\nT\n>b\nCGAT\nCGAC\nTAC\n");
        let expected = ">a\nACGATCGCT\n>b\nCGATCGACTAC\n";
        copy_fasta(&in_fasta, &out_fasta);
        assert_eq!(std::fs::read_to_string(&out_fasta).unwrap(), expected);

        // Gzipped FASTA file is converted to uncompressed FASTA
        let in_fasta = dir.path().join("in.fasta.gz");
        make_gzipped_test_file(&in_fasta, ">a\nACGATCGCT\n>b\nCGATCGACTAC\n");
        let expected = ">a\nACGATCGCT\n>b\nCGATCGACTAC\n";
        copy_fasta(&in_fasta, &out_fasta);
        assert_eq!(std::fs::read_to_string(&out_fasta).unwrap(), expected);
    }
}
