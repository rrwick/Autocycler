// This file contains miscellaneous functions used by various parts of Autocycler.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use indicatif::{ProgressBar, ProgressStyle};
use flate2::read::MultiGzDecoder;
use noodles_bam as bam;
use seq_io::fastq::{OwnedRecord, Reader, Record};
use std::collections::HashSet;
use std::fs;
use std::fs::{File, read_dir, create_dir_all, remove_dir_all};
use std::io;
use std::io::{prelude::*, BufReader, BufWriter, Read};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, Once};
use std::time::Duration;
use tempfile::{tempdir, TempDir};


pub mod strand {
    pub const FORWARD: bool = true;
    pub const REVERSE: bool = false;
}


pub fn create_dir(dir_path: &Path) {
    create_dir_all(dir_path).unwrap_or_else(|e| {
        quit_with_error(&format!("failed to create directory {}\n{}", dir_path.display(), e))
    });
}


pub fn delete_dir_if_exists(dir_path: &Path) {
    if dir_path.exists() && dir_path.is_dir() {
        remove_dir_all(dir_path).unwrap_or_else(|e| {
            quit_with_error(&format!("failed to delete directory {}\n{}", dir_path.display(), e))
        });
    }
}


pub fn load_file_lines(filename: &Path) -> Vec<String> {
    let file = File::open(filename).unwrap_or_else(|e| {
        quit_with_error(&format!("failed to open file {}\n{}", filename.display(), e));
    });
    let reader = BufReader::new(file);
    reader.lines().map(|line_result| {
        line_result.unwrap_or_else(|e| {
            quit_with_error(&format!("failed to read line\n{e}"));
        })
    }).collect()
}


pub fn find_all_assemblies(in_dir: &Path) -> Vec<PathBuf> {
    let paths = read_dir(in_dir).unwrap_or_else(|e| {
        quit_with_error(&format!("unable to read directory {}\n{}", in_dir.display(), e))
    });
    let mut all_assemblies: Vec<_> = paths.map(|entry| entry.unwrap().path())
        .filter(|path| is_assembly_file(path)).collect();
    all_assemblies.sort_unstable();
    if all_assemblies.is_empty() {
        quit_with_error(&format!("no assemblies found in {}", in_dir.display()));
    }
    all_assemblies
}


fn is_assembly_file(path: &Path) -> bool {
    let extension = path.extension().unwrap_or_default();
    let stem = path.file_stem().and_then(|stem| stem.to_str()).unwrap_or_default();
    path.is_file() &&
        (extension == "fasta" || extension == "fna" || extension == "fa" ||
         (extension == "gz" &&
          (stem.ends_with(".fasta") || stem.ends_with(".fna") || stem.ends_with(".fa"))))
}


pub fn check_if_file_exists(path: &Path) {
    if !path.exists() {
        quit_with_error(&format!("file does not exist: {}", path.display()));
    }
    if !path.is_file() {
        quit_with_error(&format!("{} is not a file", path.display()));
    }
}


pub fn check_if_dir_exists(path: &Path) {
    if !path.exists() {
        quit_with_error(&format!("directory does not exist: {}", path.display()));
    }
    if !path.is_dir() {
        quit_with_error(&format!("{} is not a directory", path.display()));
    }
}


pub fn check_if_dir_is_not_dir(dir: &Path) {
    if dir.exists() && !dir.is_dir() {
        quit_with_error(&format!("{} exists but is not a directory", dir.display()));
    }
}


#[cfg(not(test))]
pub fn quit_with_error(text: &str) -> ! {
    remove_temp_dirs();
    eprintln!();
    eprintln!("Error: {text}");
    std::process::exit(1);
}
#[cfg(test)]
pub fn quit_with_error(text: &str) -> ! {
    // Panicking lets tests catch errors without exiting the test process.
    panic!("{}", text);
}


pub fn load_fasta(filename: &Path) -> Vec<(String, String, String)> {
    if is_file_empty(filename) {
        quit_with_error(&format!("{} is an empty file", filename.display()));
    }
    let fasta_seqs = load_fasta_allow_empty(filename);
    validate_fasta(&fasta_seqs, filename);
    fasta_seqs
}


pub fn load_fasta_allow_empty(filename: &Path) -> Vec<(String, String, String)> {
    // Skips validation of empty sequences and duplicate names as well as empty files.
    text_file_reader(filename).and_then(|reader| parse_fasta(reader, filename))
        .unwrap_or_else(|e| quit_with_error(&format!("unable to load {}\n{}", filename.display(), e)))
}


fn validate_fasta(fasta_seqs: &[(String, String, String)], filename: &Path) {
    if fasta_seqs.is_empty() {
        quit_with_error(&format!("{} contains no sequences", filename.display()));
    }
    for (name, _, sequence) in fasta_seqs {
        if name.is_empty() {
            quit_with_error(&format!("{} has an unnamed sequence", filename.display()));
        }
        if sequence.is_empty() {
            quit_with_error(&format!("{} has an empty sequence", filename.display()));
        }
    }
    let mut names = HashSet::new();
    for (name, _, _) in fasta_seqs {
        if !names.insert(name) {
            quit_with_error(&format!("{} has a duplicate name: {}", filename.display(), name));
        }
    }
}


fn fastq_reader(fastq_file: &Path) -> Reader<Box<dyn Read>> {
    // Returns a reader for a FASTQ file that works on both unzipped and gzipped files.
    let file = File::open(fastq_file).expect("Error opening file");
    let reader: Box<dyn Read> = if is_file_gzipped(fastq_file) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Reader::new(reader)
}


pub fn single_fastq(files: Vec<PathBuf>) -> (PathBuf, Option<TempDir>) {
    // Returns the reads as a single FASTQ file. A single FASTQ input is returned as is, but BAM or
    // multiple input files are first converted to an uncompressed FASTQ in a temporary directory.
    if files.len() == 1 && !is_file_bam(&files[0]) {
        return (files.into_iter().next().unwrap(), None);
    }
    eprintln!("Converting input reads to a temporary FASTQ file...");
    let dir = temp_dir();
    let fastq = dir.path().join("reads.fastq");
    File::create(&fastq).and_then(|file| {
        let mut writer = BufWriter::new(file);
        for read in unique_read_iter(&files) { read.write(&mut writer)?; }
        writer.flush()
    }).unwrap_or_else(|e| quit_with_error(&format!("unable to write {}: {e}", fastq.display())));
    (fastq, Some(dir))
}


pub fn read_iter(files: &[PathBuf]) -> impl Iterator<Item = OwnedRecord> + '_ {
    // Iterates over the reads in one or more files, each of which can be FASTQ (optionally
    // gzipped) or unaligned BAM.
    files.iter().flat_map(|file| -> Box<dyn Iterator<Item = OwnedRecord>> {
        if is_file_bam(file) { Box::new(bam_reads(file)) } else { Box::new(fastq_reads(file)) }
    })
}


pub fn unique_read_iter(files: &[PathBuf]) -> impl Iterator<Item = OwnedRecord> + '_ {
    // Same as read_iter, but quits with an error if any read name occurs more than once.
    let mut names = HashSet::new();
    read_iter(files).inspect(move |read| {
        if !names.insert(read.id_bytes().to_vec()) {
            quit_with_error(&format!("duplicate read name: {}",
                                     String::from_utf8_lossy(read.id_bytes())));
        }
    })
}


pub fn read_batches(files: &[PathBuf], batch_bases: usize)
        -> impl Iterator<Item = Vec<OwnedRecord>> + '_ {
    // Same as read_iter, but groups the reads into batches of about batch_bases, which is useful
    // for processing reads in parallel.
    let mut reads = read_iter(files);
    std::iter::from_fn(move || {
        let (mut batch, mut bases) = (Vec::new(), 0);
        while bases < batch_bases {
            let Some(read) = reads.next() else { break };
            bases += read.seq.len();
            batch.push(read);
        }
        (!batch.is_empty()).then_some(batch)
    })
}


fn fastq_reads(file: &Path) -> impl Iterator<Item = OwnedRecord> {
    let file = file.to_path_buf();
    fastq_reader(&file).into_records().map(move |r| r.unwrap_or_else(|e| quit_with_error(
        &format!("unable to read {}: {e}\nAre you sure this is a FASTQ file?", file.display()))))
}


fn bam_reads(file: &Path) -> impl Iterator<Item = OwnedRecord> {
    let file = file.to_path_buf();
    let mut reader = bam::io::reader::Builder.build_from_path(&file)
        .and_then(|mut reader| reader.read_header().map(|_| reader))
        .unwrap_or_else(|e| quit_with_error(&format!("unable to read {}: {e}", file.display())));
    let mut record = bam::Record::default();
    std::iter::from_fn(move || match reader.read_record(&mut record) {
        Ok(0) => None,
        Ok(_) => Some(bam_record_to_fastq(&record, &file)),
        Err(e) => quit_with_error(&format!("unable to read {}: {e}", file.display())),
    })
}


fn bam_record_to_fastq(record: &bam::Record, file: &Path) -> OwnedRecord {
    let name = record.name().unwrap_or_else(|| {
        quit_with_error(&format!("{} contains a read with no name", file.display()))
    }).to_vec();
    let name_str = String::from_utf8_lossy(&name);
    if !record.flags().is_unmapped() {
        quit_with_error(&format!("{} contains aligned reads (e.g. {name_str}), but only unaligned \
                                  BAM is supported", file.display()));
    }
    if record.flags().is_segmented() {
        quit_with_error(&format!("{} contains paired reads (e.g. {name_str}), which are not \
                                  supported in BAM format", file.display()));
    }
    let seq: Vec<u8> = record.sequence().iter().collect();
    let qual: Vec<u8> = record.quality_scores().iter().map(|q| q.saturating_add(33)).collect();
    if qual.len() != seq.len() {
        quit_with_error(&format!("read {name_str} in {} has no quality scores", file.display()));
    }
    OwnedRecord { head: name, seq, qual }
}


fn is_file_bam(filename: &Path) -> bool {
    // BAM files are BGZF-compressed (a type of gzip) and start with a magic string.
    if !is_file_gzipped(filename) { return false; }
    let file = File::open(filename).unwrap_or_else(|e| {
        quit_with_error(&format!("unable to open {}: {e}", filename.display()))
    });
    let mut magic = [0u8; 4];
    MultiGzDecoder::new(file).read_exact(&mut magic).is_ok() && &magic == b"BAM\x01"
}


pub fn decompress_if_gzipped(filename: &Path) -> Option<(PathBuf, TempDir)> {
    // The returned temporary directory is automatically deleted when dropped.
    if !is_file_gzipped(filename) { return None; }
    let file = File::open(filename).unwrap_or_else(|e| {
        quit_with_error(&format!("unable to open {}: {e}", filename.display()))
    });
    let temp_dir = temp_dir();
    let uncompressed_name = if filename.extension().unwrap_or_default() == "gz" {
        filename.file_stem().unwrap_or_default()
    } else {
        filename.file_name().unwrap_or_default()
    };
    let temp_path = temp_dir.path().join(uncompressed_name);
    let mut temp_file = File::create(&temp_path).unwrap_or_else(|e| {
        quit_with_error(&format!("unable to create temporary file: {e}"))
    });
    io::copy(&mut MultiGzDecoder::new(file), &mut temp_file).unwrap_or_else(|e| {
        quit_with_error(&format!("unable to decompress {}: {e}", filename.display()))
    });
    Some((temp_path, temp_dir))
}


static TEMP_DIRS: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());

pub fn temp_dir() -> TempDir {
    // Creates a temporary directory which is deleted when dropped, but also when Autocycler quits
    // with an error or is interrupted with Ctrl-C (neither of which run drop).
    static CTRL_C: Once = Once::new();
    CTRL_C.call_once(|| ctrlc::set_handler(|| { remove_temp_dirs(); std::process::exit(130); })
                     .expect("failed to set Ctrl-C handler"));
    let dir = tempdir().unwrap_or_else(|e| {
        quit_with_error(&format!("unable to create temporary directory: {e}"))
    });
    TEMP_DIRS.lock().unwrap().push(dir.path().to_path_buf());
    dir
}


fn remove_temp_dirs() {
    if let Ok(dirs) = TEMP_DIRS.lock() {
        for dir in dirs.iter() { let _ = remove_dir_all(dir); }
    }
}


pub fn is_file_empty(filename: &Path) -> bool {
    fs::metadata(filename).is_ok_and(|metadata| metadata.len() == 0)
}


pub fn total_fasta_length(filename: &Path) -> usize {
    if !filename.exists() { return 0; }
    let fasta_seqs = load_fasta_allow_empty(filename);
    fasta_seqs.iter().map(|(_, _, seq)| seq.len()).sum()
}


pub fn is_fasta_empty(filename: &Path) -> bool {
    total_fasta_length(filename) == 0
}


fn is_file_gzipped(filename: &Path) -> bool {
    let file = File::open(filename).unwrap_or_else(|e| {
        quit_with_error(&format!("unable to open {}: {e}", filename.display()))
    });
    let mut buf = [0u8; 2];
    let n = BufReader::new(file).read(&mut buf)
        .unwrap_or_else(|e| {
            quit_with_error(&format!("error reading {}: {e}", filename.display()))
        });
    n == 2 && buf == [0x1f, 0x8b]
}


fn text_file_reader(filename: &Path) -> io::Result<BufReader<Box<dyn Read>>> {
    let gzipped = is_file_gzipped(filename);
    let file = File::open(filename)?;
    let reader: Box<dyn Read> = if gzipped { Box::new(MultiGzDecoder::new(file)) }
                                    else { Box::new(file) };
    Ok(BufReader::new(reader))
}


fn parse_fasta(reader: impl BufRead, filename: &Path) -> io::Result<Vec<(String, String, String)>> {
    let mut fasta_seqs = Vec::new();
    let mut name = String::new();
    let mut header = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let text = line?;
        if text.is_empty() { continue; }
        if let Some(text) = text.strip_prefix('>') {
            if !name.is_empty() {
                sequence.make_ascii_uppercase();
                fasta_seqs.push((name, header, sequence));
                sequence = String::new();
            }
            header = text.to_string();
            name = text.split_whitespace().next().unwrap_or_else(|| {
                quit_with_error(&format!("{} is not correctly formatted", filename.display()))
            }).to_string();
        } else {
            if name.is_empty() {
                quit_with_error(&format!("{} is not correctly formatted", filename.display()));
            }
            sequence.push_str(&text);
        }
    }
    if !name.is_empty() {
        sequence.make_ascii_uppercase();
        fasta_seqs.push((name, header, sequence));
    }
    Ok(fasta_seqs)
}


fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'.' => b'.',
        _ => b'N'
    }
}


pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&base| complement_base(base)).collect()
}


pub fn format_duration(duration: Duration) -> String {
    let microseconds = duration.subsec_micros();
    let seconds = duration.as_secs() % 60;
    let minutes = duration.as_secs() / 60 % 60;
    let hours = duration.as_secs() / 3600;
    format!("{hours}:{minutes:02}:{seconds:02}.{microseconds:06}")
}


pub fn usize_division_rounded(dividend: usize, divisor: usize) -> usize {
    // Divides an integer by another integer, giving the result rounded to the nearest integer.
    if divisor == 0 {
        panic!("Attempt to divide by zero");
    }
    (dividend + divisor / 2) / divisor
}


pub fn format_float(num: f64) -> String {
    // Formats a float with up to six decimal places but then drops trailing zeros.
    let mut formatted = format!("{num:.6}");
    if !formatted.contains('.') { return formatted }
    while formatted.ends_with('0') { formatted.pop(); }
    if formatted.ends_with('.') { formatted.pop(); }
    formatted
}


pub fn format_float_sigfigs(value: f64, sigfigs: usize) -> String {
    // Formats a float with the specified significant figures.
    if value == 0.0 {
        return format!("{:.*}", sigfigs - 1, 0.0);
    }
    let decimals = sigfigs as i32 - value.abs().log10().floor() as i32 - 1;
    let factor = 10f64.powi(decimals);
    let rounded_value = (value * factor).round() / factor;
    if decimals > 0 {
        format!("{:.*}", decimals as usize, rounded_value)
    } else {
        format!("{rounded_value}")
    }
}


pub fn median_usize(values: &[usize]) -> usize {
    if values.is_empty() { return 0; }
    let mut sorted_values = values.to_vec();
    sorted_values.sort_unstable();
    let len = sorted_values.len();
    if len.is_multiple_of(2) { (sorted_values[len / 2 - 1] + sorted_values[len / 2]) / 2 }
                        else { sorted_values[len / 2] }
}


pub fn median_isize(values: &[isize]) -> isize {
    if values.is_empty() { return 0; }
    let mut sorted_values = values.to_vec();
    sorted_values.sort_unstable();
    let len = sorted_values.len();
    if len.is_multiple_of(2) { (sorted_values[len / 2 - 1] + sorted_values[len / 2]) / 2 }
                        else { sorted_values[len / 2] }
}


pub fn mad_usize(values: &[usize]) -> usize {
    if values.is_empty() { return 0; }
    let median = median_usize(values);
    let absolute_deviations: Vec<_> = values.iter()
        .map(|v| (*v as isize - median as isize).abs()).collect();
    median_isize(&absolute_deviations) as usize
}


pub fn mad_isize(values: &[isize]) -> isize {
    if values.is_empty() { return 0; }
    let median = median_isize(values);
    let absolute_deviations: Vec<_> = values.iter().map(|v| (*v - median).abs()).collect();
    median_isize(&absolute_deviations)
}


pub fn spinner(message: &str) -> ProgressBar {
    if cfg!(test) {
        ProgressBar::hidden() // don't show a spinner during unit tests
    } else {
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(Duration::from_millis(100));
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠚", "⠞", "⠖", "⠦", "⠴", "⠲", "⠳", "⠓"])  // dots3 from github.com/sindresorhus/cli-spinners
                .template("{spinner} {msg}").unwrap(),
        );
        pb.set_message(message.to_string());
        pb
    }
}


pub fn reverse_path(path: &[i32]) -> Vec<i32> {
    path.iter().rev().map(|&num| -num).collect()
}


pub fn parse_node_numbers<T: std::str::FromStr + Ord>(numbers: Option<String>) -> Vec<T> {
    let Some(numbers) = numbers else { return Vec::new(); };
    let mut numbers: Vec<T> = numbers.replace(' ', "").split(',')
        .map(|s| s.parse().unwrap_or_else(|_| quit_with_error(
            &format!("failed to parse '{s}' as a node number")))).collect();
    numbers.sort();
    numbers
}


pub fn sign_at_end(num: i32) -> String {
    format!("{}{}", num.abs(), if num >= 0 { "+" } else { "-" })
}


pub fn sign_at_end_vec(nums: &[i32]) -> String {
    nums.iter().map(|&n| sign_at_end(n)).collect::<Vec<_>>().join(",")
}


pub fn up_to_first_space(string: &str) -> String {
    string.split_whitespace().next().unwrap_or("").to_string()
}


pub fn after_first_space(string: &str) -> String {
    string.split_once(char::is_whitespace).map_or("", |(_, rest)| rest).to_string()
}


pub fn first_char_in_file(filename: &Path) -> io::Result<char> {
    first_non_empty_char(text_file_reader(filename)?)
}


fn first_non_empty_char<R: BufRead>(reader: R) -> io::Result<char> {
    for line in reader.lines() {
        let text = line?;
        if let Some(first_char) = text.chars().next() {
            return Ok(first_char);
        }
    }
    Err(io::Error::new(io::ErrorKind::UnexpectedEof, "No non-empty lines found"))
}


pub fn find_replace_i32_tuple(tuple: (i32, i32), find: i32, replace: i32) -> (i32, i32) {
    let (x, y) = tuple;
    (
        if x.abs() == find { replace * x.signum() } else { x },
        if y.abs() == find { replace * y.signum() } else { y },
    )
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;
    use tempfile::tempdir;

    use crate::tests::{make_test_file, make_gzipped_test_file, make_test_bam};

    #[test]
    fn test_is_assembly_file() {
        let dir = tempdir().unwrap();
        for (extension, expected) in [
            ("fasta", true), ("fna", true), ("fa", true),
            ("fasta.gz", true), ("fna.gz", true), ("fa.gz", true),
            ("fasta.bak", false), ("fna.bak", false), ("fa.bak", false),
            ("fasta.txt", false), ("fna.txt", false), ("fa.txt", false),
            ("fa.gz.bak", false), ("gz", false), ("FASTA", false), ("FA.gz", false),
        ] {
            let path = dir.path().join(format!("sample.{extension}"));
            make_test_file(&path, ">a\nACGT\n");
            assert_eq!(is_assembly_file(&path), expected, "{}", path.display());
        }
        let directory = dir.path().join("directory.fasta.gz");
        fs::create_dir(&directory).unwrap();
        assert!(!is_assembly_file(&directory));
        assert!(!is_assembly_file(&dir.path().join("missing.fasta")));
    }

    #[test]
    fn test_decompress_if_gzipped() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("reads.fastq.gz");
        make_gzipped_test_file(&filename, "@read\nACGT\n+\nIIII\n");

        let temp_file = decompress_if_gzipped(&filename).unwrap();
        let temp_path = temp_file.0.clone();
        assert_eq!(temp_path.file_name().unwrap(), "reads.fastq");
        assert_eq!(std::fs::read_to_string(&temp_path).unwrap(),
                   "@read\nACGT\n+\nIIII\n");
        drop(temp_file);
        assert!(!temp_path.exists());

        let filename = dir.path().join("reads.fastq");
        make_test_file(&filename, "@read\nACGT\n+\nIIII\n");
        assert!(decompress_if_gzipped(&filename).is_none());
    }

    #[test]
    fn test_format_duration() {
        let d1 = std::time::Duration::from_micros(123456789);
        let d2 = std::time::Duration::from_micros(3661000001);
        let d3 = std::time::Duration::from_micros(360959000001);
        assert_eq!(format_duration(d1), "0:02:03.456789");
        assert_eq!(format_duration(d2), "1:01:01.000001");
        assert_eq!(format_duration(d3), "100:15:59.000001");
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"GGTATCACTCAGGAAGC"), b"GCTTCCTGAGTGATACC");
        assert_eq!(reverse_complement(b"XYZ"), b"NNN");
        assert_eq!(reverse_complement(b".acgtnACGT."), b".ACGTNNNNN.");
        assert!(reverse_complement(b"").is_empty());
    }

    #[test]
    fn test_usize_division_rounded() {
        assert_eq!(usize_division_rounded(0, 3), 0);
        assert_eq!(usize_division_rounded(1, 3), 0);
        assert_eq!(usize_division_rounded(2, 3), 1);
        assert_eq!(usize_division_rounded(3, 3), 1);
        assert_eq!(usize_division_rounded(4, 3), 1);
        assert_eq!(usize_division_rounded(5, 3), 2);
        assert_eq!(usize_division_rounded(6, 3), 2);
        assert_eq!(usize_division_rounded(7, 3), 2);
        assert_eq!(usize_division_rounded(8, 3), 3);
        assert_eq!(usize_division_rounded(9, 3), 3);
        assert_eq!(usize_division_rounded(10, 3), 3);
        assert_eq!(usize_division_rounded(0, 1), 0);
        assert_eq!(usize_division_rounded(1, 1), 1);
        assert_eq!(usize_division_rounded(10, 1), 10);
    }

    #[test]
    fn test_format_float() {
        assert_eq!(format_float(0.0), "0");
        assert_eq!(format_float(0.1), "0.1");
        assert_eq!(format_float(0.11), "0.11");
        assert_eq!(format_float(0.111), "0.111");
        assert_eq!(format_float(0.1111), "0.1111");
        assert_eq!(format_float(0.11111), "0.11111");
        assert_eq!(format_float(0.111111), "0.111111");
        assert_eq!(format_float(0.1111111), "0.111111");
        assert_eq!(format_float(0.11111111), "0.111111");
        assert_eq!(format_float(10.0), "10");
    }

    #[test]
    fn test_format_float_sigfigs() {
        assert_eq!(format_float_sigfigs(0.0, 1), "0");
        assert_eq!(format_float_sigfigs(0.0, 2), "0.0");
        assert_eq!(format_float_sigfigs(0.0, 5), "0.0000");

        assert_eq!(format_float_sigfigs(0.1, 1), "0.1");
        assert_eq!(format_float_sigfigs(0.1, 3), "0.100");

        assert_eq!(format_float_sigfigs(0.12345678, 1), "0.1");
        assert_eq!(format_float_sigfigs(0.12345678, 2), "0.12");
        assert_eq!(format_float_sigfigs(0.12345678, 3), "0.123");
        assert_eq!(format_float_sigfigs(0.12345678, 4), "0.1235");
        assert_eq!(format_float_sigfigs(0.12345678, 5), "0.12346");
        assert_eq!(format_float_sigfigs(0.12345678, 6), "0.123457");
        assert_eq!(format_float_sigfigs(0.12345678, 7), "0.1234568");
        assert_eq!(format_float_sigfigs(0.12345678, 8), "0.12345678");
        assert_eq!(format_float_sigfigs(0.12345678, 9), "0.123456780");

        assert_eq!(format_float_sigfigs(87.654321, 1), "90");
        assert_eq!(format_float_sigfigs(87.654321, 2), "88");
        assert_eq!(format_float_sigfigs(87.654321, 3), "87.7");
        assert_eq!(format_float_sigfigs(87.654321, 4), "87.65");
        assert_eq!(format_float_sigfigs(87.654321, 5), "87.654");
        assert_eq!(format_float_sigfigs(87.654321, 6), "87.6543");
        assert_eq!(format_float_sigfigs(87.654321, 7), "87.65432");
        assert_eq!(format_float_sigfigs(87.654321, 8), "87.654321");
        assert_eq!(format_float_sigfigs(87.654321, 9), "87.6543210");

        assert_eq!(format_float_sigfigs(-0.12345678, 1), "-0.1");
        assert_eq!(format_float_sigfigs(-0.12345678, 2), "-0.12");
        assert_eq!(format_float_sigfigs(-0.12345678, 3), "-0.123");
        assert_eq!(format_float_sigfigs(-0.12345678, 4), "-0.1235");
        assert_eq!(format_float_sigfigs(-0.12345678, 5), "-0.12346");
        assert_eq!(format_float_sigfigs(-0.12345678, 6), "-0.123457");
        assert_eq!(format_float_sigfigs(-0.12345678, 7), "-0.1234568");
        assert_eq!(format_float_sigfigs(-0.12345678, 8), "-0.12345678");
        assert_eq!(format_float_sigfigs(-0.12345678, 9), "-0.123456780");

        assert_eq!(format_float_sigfigs(-87.654321, 1), "-90");
        assert_eq!(format_float_sigfigs(-87.654321, 2), "-88");
        assert_eq!(format_float_sigfigs(-87.654321, 3), "-87.7");
        assert_eq!(format_float_sigfigs(-87.654321, 4), "-87.65");
        assert_eq!(format_float_sigfigs(-87.654321, 5), "-87.654");
        assert_eq!(format_float_sigfigs(-87.654321, 6), "-87.6543");
        assert_eq!(format_float_sigfigs(-87.654321, 7), "-87.65432");
        assert_eq!(format_float_sigfigs(-87.654321, 8), "-87.654321");
        assert_eq!(format_float_sigfigs(-87.654321, 9), "-87.6543210");

        assert_eq!(format_float_sigfigs(0.0005182, 1), "0.0005");
        assert_eq!(format_float_sigfigs(0.0005182, 2), "0.00052");
        assert_eq!(format_float_sigfigs(0.0005182, 3), "0.000518");
        assert_eq!(format_float_sigfigs(0.0005182, 4), "0.0005182");
        assert_eq!(format_float_sigfigs(0.0005182, 5), "0.00051820");
        assert_eq!(format_float_sigfigs(0.0005182, 6), "0.000518200");

        assert_eq!(format_float_sigfigs(907.001, 1), "900");
        assert_eq!(format_float_sigfigs(907.001, 2), "910");
        assert_eq!(format_float_sigfigs(907.001, 3), "907");
        assert_eq!(format_float_sigfigs(907.001, 4), "907.0");
        assert_eq!(format_float_sigfigs(907.001, 5), "907.00");
        assert_eq!(format_float_sigfigs(907.001, 6), "907.001");
        assert_eq!(format_float_sigfigs(907.001, 7), "907.0010");
        assert_eq!(format_float_sigfigs(907.001, 8), "907.00100");
    }

    #[test]
    fn test_median() {
        assert_eq!(median_usize(&[]), 0);
        assert_eq!(median_usize(&[0, 1, 2, 3, 4]), 2);
        assert_eq!(median_usize(&[4, 3, 2, 1, 0]), 2);
        assert_eq!(median_usize(&[0, 1, 2, 3, 4, 5]), 2);
        assert_eq!(median_usize(&[5, 4, 3, 2, 1, 0]), 2);
        assert_eq!(median_usize(&[0, 2, 4, 6, 8, 10]), 5);
        assert_eq!(median_usize(&[10, 8, 6, 4, 2, 0]), 5);

        assert_eq!(median_isize(&[]), 0);
        assert_eq!(median_isize(&[-4, -1, -2]), -2);
        assert_eq!(median_isize(&[-4, -1, -2, -3]), -2);
        assert_eq!(median_isize(&[-2, 1]), 0);
        assert_eq!(median_isize(&[0, 1, 2, 3, 4]), 2);
        assert_eq!(median_isize(&[4, 3, 2, 1, 0]), 2);
        assert_eq!(median_isize(&[0, 1, 2, 3, 4, 5]), 2);
        assert_eq!(median_isize(&[5, 4, 3, 2, 1, 0]), 2);
        assert_eq!(median_isize(&[0, 2, 4, 6, 8, 10]), 5);
        assert_eq!(median_isize(&[10, 8, 6, 4, 2, 0]), 5);
    }

    #[test]
    fn test_median_absolute_deviation() {
        assert_eq!(mad_usize(&[]), 0);
        assert_eq!(mad_usize(&[1, 1, 2, 2, 4, 6, 9]), 1);
        assert_eq!(mad_usize(&[4, 1, 9, 6, 1, 2, 2]), 1);

        assert_eq!(mad_isize(&[]), 0);
        assert_eq!(mad_isize(&[1, 1, 2, 2, 4, 6, 9]), 1);
        assert_eq!(mad_isize(&[4, 1, 9, 6, 1, 2, 2]), 1);
    }

    #[test]
    fn test_parse_node_numbers() {
        assert_eq!(parse_node_numbers::<u16>(None), Vec::<u16>::new());
        assert_eq!(parse_node_numbers::<u32>(None), Vec::<u32>::new());
        for (input, expected) in [("1,2,3", vec![1, 2, 3]), ("4, 5, 6", vec![4, 5, 6]),
                                   ("  5 , 10 ,15 ", vec![5, 10, 15]), ("3,1,3,0", vec![0, 1, 3, 3])] {
            assert_eq!(parse_node_numbers::<u32>(Some(input.to_string())), expected);
            assert_eq!(parse_node_numbers::<u16>(Some(input.to_string())),
                       expected.iter().map(|&n| n as u16).collect::<Vec<_>>());
        }
        for input in ["", "ABC", "1,X,3", "x,y,z", "^&%^*", "1,,2", "1,\t2", "-1", "4294967296"] {
            assert!(panic::catch_unwind(|| parse_node_numbers::<u16>(Some(input.to_string()))).is_err());
            assert!(panic::catch_unwind(|| parse_node_numbers::<u32>(Some(input.to_string()))).is_err());
        }
        assert_panics_with(|| { parse_node_numbers::<u16>(Some("65536".to_string())); },
                           "failed to parse '65536' as a node number");
        assert_eq!(parse_node_numbers::<u32>(Some("65536,4294967295".to_string())), vec![65536, u32::MAX]);
    }

    #[test]
    fn test_reverse_path() {
        assert_eq!(reverse_path(&[1, -2]), vec![2, -1]);
        assert_eq!(reverse_path(&[4, 8, -3]), vec![3, -8, -4]);
    }

    #[test]
    fn test_sign_at_end() {
        assert_eq!(sign_at_end(123), "123+".to_string());
        assert_eq!(sign_at_end(-321), "321-".to_string());
    }

    #[test]
    fn test_sign_at_end_vec() {
        assert_eq!(sign_at_end_vec(&[8]), "8+".to_string());
        assert_eq!(sign_at_end_vec(&[123, -321]), "123+,321-".to_string());
        assert_eq!(sign_at_end_vec(&[-4, -5, 67, 34345, 1]), "4-,5-,67+,34345+,1+".to_string());
    }

    #[test]
    fn test_up_to_first_space() {
        assert_eq!(up_to_first_space("xyz"), "xyz".to_string());
        assert_eq!(up_to_first_space("1 2 3 4"), "1".to_string());
        assert_eq!(up_to_first_space("abc def"), "abc".to_string());
    }

    #[test]
    fn test_after_first_space() {
        assert_eq!(after_first_space("xyz"), "".to_string());
        assert_eq!(after_first_space("1 2 3 4"), "2 3 4".to_string());
        assert_eq!(after_first_space("abc def"), "def".to_string());
    }

    #[test]
    fn test_first_char_in_file() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");
        for write_file in [make_test_file, make_gzipped_test_file] {
            for (text, expected) in [(">a\nACGT\n", '>'), ("XYZ", 'X'),
                                     ("\n\r\n>seq", '>'), ("\n é", ' '), ("é", 'é')] {
                write_file(&filename, text);
                assert_eq!(first_char_in_file(&filename).unwrap(), expected);
            }
            for text in ["", "\n\r\n"] {
                write_file(&filename, text);
                let error = first_char_in_file(&filename).unwrap_err();
                assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
                assert_eq!(error.to_string(), "No non-empty lines found");
            }
        }
    }

    #[test]
    fn test_load_fasta() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");
        let expected = vec![("a".to_string(), "a".to_string(), "ACGT".to_string()),
                            ("b".to_string(), "b xyz".to_string(), "ACGTACGT".to_string())];
        for write_file in [make_test_file, make_gzipped_test_file] {
            write_file(&filename, "\n>a\r\naCgT\r\n\n>b xyz\nacgt\nACGT");
            assert_eq!(load_fasta(&filename), expected);
            assert_eq!(load_fasta_allow_empty(&filename), expected);

            write_file(&filename, ">  a\tinfo\nac gté\n");
            assert_eq!(load_fasta(&filename),
                       vec![("a".into(), "  a\tinfo".into(), "AC GTé".into())]);
        }
    }

    fn assert_panics_with(action: impl FnOnce() + panic::UnwindSafe, expected: &str) {
        let error = panic::catch_unwind(action).unwrap_err();
        assert_eq!(error.downcast_ref::<String>().unwrap(), expected);
    }

    #[test]
    fn test_load_fasta_validation() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");
        for write_file in [make_test_file, make_gzipped_test_file] {
            for (text, error) in [("\n\n", "contains no sequences"),
                                  (">a\n", "has an empty sequence"),
                                  (">a\nA\n>a\nT\n", "has a duplicate name: a"),
                                  (">a\nA\n>a\nT\n>b\n", "has an empty sequence")] {
                write_file(&filename, text);
                load_fasta_allow_empty(&filename);
                assert_panics_with(|| { load_fasta(&filename); },
                                   &format!("{} {error}", filename.display()));
            }
            for text in ["ACGT\n", ">\n", "> \t\n", ">a\nA\n>\n"] {
                write_file(&filename, text);
                let error = format!("{} is not correctly formatted", filename.display());
                assert_panics_with(|| { load_fasta(&filename); }, &error);
                assert_panics_with(|| { load_fasta_allow_empty(&filename); }, &error);
            }
            write_file(&filename, "");
            assert!(load_fasta_allow_empty(&filename).is_empty());
            let error = if is_file_empty(&filename) { "is an empty file" }
                                              else { "contains no sequences" };
            assert_panics_with(|| { load_fasta(&filename); },
                               &format!("{} {error}", filename.display()));
        }
    }

    #[test]
    fn test_load_fasta_concatenated_gzip() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");
        make_gzipped_test_file(&filename, ">a\nac");
        let mut contents = fs::read(&filename).unwrap();
        make_gzipped_test_file(&filename, "gt\n>b\ntt\n");
        contents.extend(fs::read(&filename).unwrap());
        fs::write(&filename, contents).unwrap();
        assert_eq!(load_fasta(&filename),
                   vec![("a".into(), "a".into(), "ACGT".into()),
                        ("b".into(), "b".into(), "TT".into())]);
    }

    #[test]
    fn test_is_file_empty() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");

        make_test_file(&filename, "");
        assert!(is_file_empty(&filename));

        make_gzipped_test_file(&filename, "");
        assert!(!is_file_empty(&filename));

        make_test_file(&filename, "x");
        assert!(!is_file_empty(&filename));
    }

    #[test]
    fn test_find_replace_i32_tuple() {
        assert_eq!(find_replace_i32_tuple((1, 2), 1, 8), (8, 2));
        assert_eq!(find_replace_i32_tuple((1, 2), 2, 8), (1, 8));
        assert_eq!(find_replace_i32_tuple((-1, 2), 1, 8), (-8, 2));
        assert_eq!(find_replace_i32_tuple((-1, 2), 2, 8), (-1, 8));
        assert_eq!(find_replace_i32_tuple((1, -2), 1, 8), (8, -2));
        assert_eq!(find_replace_i32_tuple((1, -2), 2, 8), (1, -8));
    }

    #[test]
    fn test_total_fasta_length() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");

        make_test_file(&filename, ">a\nACGT\n");
        assert_eq!(total_fasta_length(&filename), 4);

        make_test_file(&filename, ">a\nACGT\n>b xyz\nACGT\nACGT\n");
        assert_eq!(total_fasta_length(&filename), 12);

        make_test_file(&filename, "");
        assert_eq!(total_fasta_length(&filename), 0);
    }

    #[test]
    fn test_is_fasta_empty() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");

        make_test_file(&filename, "");
        assert!(is_fasta_empty(&filename));

        make_gzipped_test_file(&filename, "");
        assert!(is_fasta_empty(&filename));

        make_test_file(&filename, ">a\nACGT\n");
        assert!(!is_fasta_empty(&filename));

        make_test_file(&filename, ">a\n\n");
        assert!(is_fasta_empty(&filename));
    }

    fn read_tuples(files: &[PathBuf]) -> Vec<(String, String, String)> {
        let s = |b: &[u8]| String::from_utf8(b.to_vec()).unwrap();
        read_iter(files).map(|r| (s(&r.head), s(&r.seq), s(&r.qual))).collect()
    }

    #[test]
    fn test_read_iter() {
        let dir = tempdir().unwrap();
        let (fq, gz, bam) = (dir.path().join("a.fastq"), dir.path().join("b.fastq.gz"),
                             dir.path().join("c.bam"));
        make_test_file(&fq, "@a1 info\nACGT\n+\nIIII\n");
        make_gzipped_test_file(&gz, "@b1\nGGA\n+\n!!5\n");
        make_test_bam(&bam, &[("c1", 4, "TTAC", &[0, 10, 20, 40]), ("c2", 4, "A", &[93])]);
        assert!(!is_file_bam(&fq) && !is_file_bam(&gz) && is_file_bam(&bam));
        assert_eq!(read_tuples(&[fq, gz, bam]),
                   [("a1 info", "ACGT", "IIII"), ("b1", "GGA", "!!5"), ("c1", "TTAC", "!+5I"),
                    ("c2", "A", "~")].map(|(a, b, c)| (a.into(), b.into(), c.into())));
    }

    #[test]
    fn test_read_iter_bam_errors() {
        let dir = tempdir().unwrap();
        let bam = dir.path().join("reads.bam");
        for (flags, qual) in [(0, &[30u8][..]), (1 | 4 | 64, &[30]), (4, &[])] {
            make_test_bam(&bam, &[("r1", flags, "A", qual)]);
            let files = vec![bam.clone()];
            assert!(panic::catch_unwind(|| read_tuples(&files)).is_err());
        }
    }

    #[test]
    fn test_unique_read_iter() {
        let dir = tempdir().unwrap();
        let (fq, bam) = (dir.path().join("reads.fastq"), dir.path().join("reads.bam"));
        make_test_file(&fq, "@r1 x\nA\n+\nI\n@r2\nA\n+\nI\n");
        make_test_bam(&bam, &[("r3", 4, "A", &[30])]);
        assert_eq!(unique_read_iter(&[fq.clone(), bam.clone()]).count(), 3);
        assert!(panic::catch_unwind(|| unique_read_iter(&[fq.clone(), fq.clone()]).count()).is_err());
        make_test_file(&fq, "@r1 x\nA\n+\nI\n@r1 y\nA\n+\nI\n");
        let files = vec![fq];
        assert_eq!(read_iter(&files).count(), 2);
        assert!(panic::catch_unwind(|| unique_read_iter(&files).count()).is_err());
    }

    #[test]
    fn test_read_batches() {
        let dir = tempdir().unwrap();
        let files = vec![dir.path().join("reads.fastq")];
        make_test_file(&files[0], "@r1\nAAA\n+\nIII\n@r2\nAA\n+\nII\n@r3\nA\n+\nI\n");
        let sizes = |n| read_batches(&files, n).map(|b| b.len()).collect::<Vec<_>>();
        assert_eq!(sizes(1), [1, 1, 1]);
        assert_eq!(sizes(4), [2, 1]);
        assert_eq!(sizes(100), [3]);
    }

    #[test]
    fn test_single_fastq() {
        let dir = tempdir().unwrap();
        let (fq, bam) = (dir.path().join("reads.fastq"), dir.path().join("reads.bam"));
        make_test_file(&fq, "@r1 x\nACGT\n+\nIIII\n");
        make_test_bam(&bam, &[("r2", 4, "GGA", &[0, 10, 20])]);
        assert!(matches!(single_fastq(vec![fq.clone()]), (p, None) if p == fq));
        for (files, expected) in [(vec![bam.clone()], "@r2\nGGA\n+\n!+5\n"),
                                  (vec![fq.clone(), bam.clone()],
                                   "@r1 x\nACGT\n+\nIIII\n@r2\nGGA\n+\n!+5\n")] {
            let (temp_fq, temp_dir) = single_fastq(files);
            assert_eq!(std::fs::read_to_string(&temp_fq).unwrap(), expected);
            assert!(TEMP_DIRS.lock().unwrap().contains(&temp_dir.unwrap().path().to_path_buf()));
        }
        assert!(panic::catch_unwind(|| single_fastq(vec![fq.clone(), fq.clone()])).is_err());
    }
}
