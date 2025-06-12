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
use std::collections::HashMap;
use std::fs::{File, OpenOptions, copy, remove_file, create_dir_all, remove_dir_all, read_dir};
use std::io::{BufRead, BufReader, BufWriter, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Once;
use which::which;
use tempfile::{tempdir, NamedTempFile, TempDir};

use crate::log::bold;
use crate::misc::{check_if_file_exists, quit_with_error, total_fasta_length, load_fasta};
use crate::subsample::parse_genome_size;


pub fn helper(task: Task, reads: PathBuf, out_prefix: Option<PathBuf>, genome_size: Option<String>,
              threads: usize, dir: Option<PathBuf>, read_type: ReadType,
              min_depth_absolute: Option<f64>, min_depth_relative: Option<f64>,
              extra_args: Vec<String>) {
    check_if_file_exists(&reads);
    let (dir, _guard) = get_working_dir(dir);
    let out_prefix = check_prefix(out_prefix);

    match task {
        Task::Genomesize => {
            genome_size_raven(reads, threads, dir, extra_args);
        }
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
            necat(reads, &out_prefix, genome_size, threads, dir, read_type, extra_args);
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

    if min_depth_absolute.is_some() || min_depth_relative.is_some() {
        depth_filter(&out_prefix, &min_depth_absolute, &min_depth_relative);
    }
}


fn canu(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>,
        threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/marbl/canu

    let genome_size = get_genome_size(genome_size, "Canu");
    check_requirements(&["canu"]);
    // TODO
}


fn flye(reads: PathBuf, out_prefix: &Path,
        threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/mikolmogorov/Flye

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

    let mut cmd = Command::new("flye");
    cmd.arg(input_flag).arg(&reads)
       .arg("--threads").arg(threads.to_string())
       .arg("--out-dir").arg(&dir);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    check_fasta(&dir.join("assembly.fasta"));
    copy_flye_fasta(&dir.join("assembly.fasta"), &dir.join("assembly_info.txt"), &fasta);
    copy_output_file(&dir.join("assembly_graph.gfa"), &gfa);
    copy_output_file(&dir.join("flye.log"), &log);
}


fn lja(reads: PathBuf, out_prefix: &Path, threads: usize, dir: PathBuf, extra_args: Vec<String>) {
    // https://github.com/AntonBankevich/LJA

    check_requirements(&["lja"]);
    let fasta = out_prefix.with_extension("fasta");
    let gfa = out_prefix.with_extension("gfa");
    let log = out_prefix.with_extension("log");

    let mut cmd = Command::new("lja");
    cmd.arg("--output-dir").arg(&dir)
       .arg("--reads").arg(&reads)
       .arg("--threads").arg(threads.to_string());
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    check_fasta(&dir.join("assembly.fasta"));
    copy_output_file(&dir.join("assembly.fasta"), &fasta);
    copy_output_file(&dir.join("mdbg.gfa"), &gfa);
    copy_output_file(&dir.join("dbg.log"), &log);
}


fn metamdbg(reads: PathBuf, out_prefix: &Path,
            threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/GaetanBenoitDev/metaMDBG

    check_requirements(&["metaMDBG"]);
    // TODO
}


fn miniasm(reads: PathBuf, out_prefix: &Path,
           threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/lh3/miniasm
    // https://github.com/rrwick/Minipolish

    check_requirements(&["miniasm", "minipolish", "minimap2", "racon"]);
    let paf = dir.join("overlap.paf");
    let unpolished = dir.join("unpolished.gfa");
    let fasta = out_prefix.with_extension("fasta");
    let gfa = out_prefix.with_extension("gfa");

    let ava_preset = match read_type {
        ReadType::OntR9      => "ava-ont",
        ReadType::OntR10     => "ava-pb",
        ReadType::PacbioClr  => "ava-pb",
        ReadType::PacbioHifi => "ava-pb",
    };

    let mut cmd = Command::new("minimap2");
    cmd.arg("-x").arg(ava_preset)
       .arg("-t").arg(threads.to_string())
       .arg(&reads).arg(&reads);
    redirect_stderr_and_stdout(&mut cmd, Some(&paf));
    run_command(&mut cmd);

    let mut cmd = Command::new("miniasm");
    cmd.arg("-f").arg(&reads)
       .arg(&paf);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, Some(&unpolished));
    run_command(&mut cmd);

    let mut cmd = Command::new("minipolish");
    cmd.arg("--threads").arg(threads.to_string())
       .arg(&reads)
       .arg(&unpolished);
    redirect_stderr_and_stdout(&mut cmd, Some(&gfa));
    run_command(&mut cmd);

    minpolish_gfa_to_fasta(&gfa, &fasta);
    check_fasta(&fasta);
}


fn myloasm(reads: PathBuf, out_prefix: &Path,
           threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/bluenote-1577/myloasm

    check_requirements(&["myloasm"]);
    let fasta = out_prefix.with_extension("fasta");
    let gfa = out_prefix.with_extension("gfa");
    let log = out_prefix.with_extension("log");

    let mut cmd = Command::new("myloasm");
    cmd.arg("--output-dir").arg(&dir)
       .arg(&reads)
       .arg("--threads").arg(threads.to_string());
    if read_type == ReadType::PacbioHifi { cmd.arg("--hifi"); }
    else if read_type == ReadType::OntR10 { cmd.arg("--nano-r10"); }
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    check_fasta(&dir.join("assembly_primary.fa"));
    copy_output_file(&dir.join("assembly_primary.fa"), &fasta);
    copy_output_file(&dir.join("final_contig_graph.gfa"), &gfa);
    copy_output_file(&find_myloasm_log(&dir), &log);
}


fn necat(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>,
         threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/xiaochuanle/NECAT

    let genome_size = get_genome_size(genome_size, "NECAT");
    let necat = find_necat();
    // TODO
}


fn nextdenovo(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>,
              threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/Nextomics/NextDenovo
    // https://github.com/Nextomics/NextPolish

    let genome_size = get_genome_size(genome_size, "NextDenovo");
    check_requirements(&["nextDenovo", "nextPolish"]);
    // TODO
}


fn plassembler(reads: PathBuf, out_prefix: &Path,
               threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/gbouras13/plassembler

    check_requirements(&["plassembler"]);
    // TODO
}


fn raven(reads: PathBuf, out_prefix: &Path, threads: usize, extra_args: Vec<String>) {
    // https://github.com/lbcb-sci/raven

    check_requirements(&["raven"]);
    let fasta = out_prefix.with_extension("fasta");
    let gfa = out_prefix.with_extension("gfa");

    let mut cmd = Command::new("raven");
    cmd.arg("--threads").arg(threads.to_string())
       .arg("--disable-checkpoints")
       .arg("--graphical-fragment-assembly").arg(&gfa)
       .arg(&reads);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, Some(&fasta));
    run_command(&mut cmd);

    check_fasta(&fasta);
}


fn genome_size_raven(reads: PathBuf, threads: usize, dir: PathBuf, extra_args: Vec<String>) {
    check_requirements(&["raven"]);
    let fasta = dir.join("assembly.fasta");

    let mut cmd = Command::new("raven");
    cmd.arg("--threads").arg(threads.to_string())
       .arg("--disable-checkpoints")
       .arg(&reads);
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, Some(&fasta));
    run_command(&mut cmd);

    check_fasta(&fasta);
    println!("{}", total_fasta_length(&fasta));
}


fn redbean(reads: PathBuf, out_prefix: &Path, genome_size: Option<String>,
           threads: usize, dir: PathBuf, read_type: ReadType, extra_args: Vec<String>) {
    // https://github.com/ruanjue/wtdbg2

    let genome_size = get_genome_size(genome_size, "Redbean");
    check_requirements(&["wtdbg2", "wtpoa-cns"]);
    let fasta = out_prefix.with_extension("fasta");

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
       .arg("-o").arg(&dir.join("dbg"));
    for token in extra_args { cmd.arg(token); }
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    let mut cmd = Command::new("wtpoa-cns");
    cmd.arg("-t").arg(threads.to_string())
       .arg("-i").arg(&dir.join("dbg.ctg.lay.gz"))
       .arg("-f")
       .arg("-o").arg(&dir.join("assembly.fasta"));
    redirect_stderr_and_stdout(&mut cmd, None);
    run_command(&mut cmd);

    check_fasta(&dir.join("assembly.fasta"));
    copy_output_file(&dir.join("assembly.fasta"), &fasta);
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


fn copy_output_file(src: &Path, dest: &Path) {
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


fn run_command(cmd: &mut Command) {
    let name = cmd.get_program().to_string_lossy().into_owned();
    print_command(&cmd);
    let status = cmd.status().unwrap_or_else(|e| {
        quit_with_error(&format!("failed to launch {name}: {e}"))
    });
    if !status.success() {
        quit_with_error(&format!("{name} exited with status {status}"));
    }
}


fn redirect_stderr_and_stdout(cmd: &mut Command, stdout_file: Option<&Path>) {
    // Redirects the command's stdout and stderr to the terminal.
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


fn find_myloasm_log(dir: &Path) -> PathBuf {
    read_dir(dir).unwrap().filter_map(|e| e.ok().map(|e| e.path()))
        .find(|p| {
            p.file_name().and_then(|s| s.to_str())
             .map_or(false, |name| { name.starts_with("myloasm_") && name.ends_with(".log") })
        })
        .unwrap_or_else(|| quit_with_error("myloasm log file not found"))
}


fn minpolish_gfa_to_fasta(gfa: &Path, fasta: &Path) {
    let reader = BufReader::new(File::open(gfa).unwrap());
    let mut writer = BufWriter::new(File::create(fasta).unwrap());
    for line in reader.lines() {
        let line = line.unwrap();
        if !line.starts_with('S') { continue; }
        let mut cols = line.split('\t');
        cols.next();
        let name = cols.next().unwrap_or("");
        let seq  = cols.next().unwrap_or("");
        let depth = cols.find_map(|field| { if field.starts_with("dp:f:") { Some(&field[5..]) }
                                                                     else { None } });
        let mut header = format!(">{name}");
        if name.ends_with('c') { header.push_str(" circular=true"); }
        if let Some(d) = depth { header.push_str(&format!(" depth={d}")); }
        writeln!(writer, "{header}\n{seq}").unwrap();
    }
}


fn copy_flye_fasta(flye_fasta: &Path, assembly_info: &Path, out_fasta: &Path) {
    let mut writer = BufWriter::new(File::create(out_fasta).unwrap());
    let info = load_flye_assembly_info(assembly_info);
    for (name, _, seq) in load_fasta(flye_fasta) {
        let mut header = format!(">{name}");
        if let Some((is_circ, depth)) = info.get(&name) {
            if *is_circ { header.push_str(" circular=true"); }
            header.push_str(&format!(" depth={depth}"));
        }
        writeln!(writer, "{header}\n{seq}").unwrap();
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


fn depth_filter(out_prefix: &Path, min_depth_absolute: &Option<f64>,
                min_depth_relative: &Option<f64>) {
    // Filters the final FASTA file by depth, overwriting the original file. If no depths are
    // available, does nothing.
    if min_depth_absolute.is_none() && min_depth_relative.is_none() { return; }
    let fasta = out_prefix.with_extension("fasta");
    if !fasta.exists() { return; }

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
    eprintln!("\nAutocycler helper depth filter threshold = {:.3}", threshold);

    let kept: Vec<_> = records.into_iter().filter_map(|(name, header, seq, depth)| {
        let pass = depth >= threshold;
        eprintln!("  {name}: depth={:.3}, {}", depth, if pass { "PASS" } else { "FAIL" });
        if pass { Some((header, seq)) } else { None }
    }).collect();

    if kept.is_empty() { let _ = remove_file(&fasta); return; }
    let mut w = BufWriter::new(File::create(&fasta).unwrap());
    for (header, seq) in kept { writeln!(w, ">{header}\n{seq}").unwrap(); }
}


fn depth_from_header(header: &str) -> Option<f64> {
    fn parse_num(s: &str) -> Option<f64> {
        s.split(|c: char| c == '-' || c == '_' || c == ' ')
            .next().and_then(|n| n.parse::<f64>().ok())
    }
    if let Some(i) = header.find("depth=")    { return parse_num(&header[i + 6..]); }
    if let Some(i) = header.find("depth-")    { return parse_num(&header[i + 6..]); }
    if let Some(i) = header.find("coverage=") { return parse_num(&header[i + 9..]); }
    None
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;

    use crate::tests::make_test_file;

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

    #[test]
    fn test_depth_from_header() {
        assert_eq!(depth_from_header(">contig depth=10.5"), Some(10.5));
        assert_eq!(depth_from_header(">contig circular=true depth=5.0"), Some(5.0));
        assert_eq!(depth_from_header(">contig"), None);
        assert_eq!(depth_from_header(">a_len-12_circular-no_depth-37-37-37_mult-2.00"), Some(37.0));
        assert_eq!(depth_from_header(">b_len-9_circular-yes_depth-25-24-23_mult-1.00"), Some(25.0));
        assert_eq!(depth_from_header(">ctg15 length=123 coverage=49.70 circular=yes"), Some(49.7));
    }


    #[test]
    fn test_depth_filter() {
        let dir = tempdir().unwrap();
        let out_prefix = dir.path().join("test");
        let fasta = out_prefix.with_extension("fasta");

        make_test_file(&fasta, ">a depth=20\n\
                                ACGT\n\
                                >b depth=120\n\
                                CGA\n\
                                >c depth=200\n\
                                ACAGACTACGACTACGACGACGATCAGCGACATCGACGT\n\
                                >d depth=100\n\
                                CGATCGACTACC\n");

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
}
