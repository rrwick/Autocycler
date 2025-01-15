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
use seq_io::fastq::Reader;
use std::collections::HashSet;
use std::fs;
use std::fs::{File, read_dir, create_dir_all, remove_dir_all};
use std::io;
use std::io::{prelude::*, BufReader, Read};
use std::path::{Path, PathBuf};
use std::time::Duration;


pub mod strand {
    // This module lets me use strand::FORWARD for true and strand::REVERSE for false.
    pub const FORWARD: bool = true;
    pub const REVERSE: bool = false;
}


pub fn create_dir(dir_path: &Path) {
    match create_dir_all(dir_path) {
        Ok(_) => {},
        Err(e) => quit_with_error(&format!("failed to create directory {}\n{}", dir_path.display(), e)),
    }
}


pub fn delete_dir_if_exists(dir_path: &Path) {
    if dir_path.exists() && dir_path.is_dir() {
        match remove_dir_all(dir_path) {
            Ok(_) => {},
            Err(e) => quit_with_error(&format!("failed to delete directory {}\n{}", dir_path.display(), e)),
        }
    }
}


pub fn load_file_lines(filename: &Path) -> Vec<String> {
    let file = File::open(filename).unwrap_or_else(|e| {
        quit_with_error(&format!("failed to open file {}\n{}", filename.display(), e));
    });
    let reader = BufReader::new(file);
    reader.lines().map(|line_result| {
        line_result.unwrap_or_else(|e| {
            quit_with_error(&format!("failed to read line\n{}", e));
        })
    }).collect()
}


pub fn find_all_assemblies(in_dir: &Path) -> Vec<PathBuf> {
    let paths = match read_dir(in_dir) {
        Ok(paths) => paths,
        Err(e) => {
            quit_with_error(&format!("unable to read directory {}\n{}", in_dir.display(), e));
        },
    };
    let mut all_assemblies: Vec<PathBuf> = Vec::new();
    for path in paths {
        let path = path.unwrap().path();
        if is_assembly_file(&path) {
            all_assemblies.push(path);
        }
    }
    all_assemblies.sort_unstable();
    if all_assemblies.is_empty() {
        quit_with_error(&format!("no assemblies found in {}", in_dir.display()));
    }
    all_assemblies
}


fn is_assembly_file(path: &Path) -> bool {
    path.is_file() && 
        (path.extension().unwrap_or_default() == "fasta" ||
         path.extension().unwrap_or_default() == "fna" ||
         path.extension().unwrap_or_default() == "fa" ||
         (path.extension().unwrap_or_default() == "gz" &&
          path.file_stem().unwrap_or_default().to_str().unwrap_or_default().ends_with(".fasta") ||
          path.file_stem().unwrap_or_default().to_str().unwrap_or_default().ends_with(".fna") ||
          path.file_stem().unwrap_or_default().to_str().unwrap_or_default().ends_with(".fa")))
}


pub fn check_if_file_exists(filename: &Path) {
    // Quits with an error if the given path is not an existing file.
    let path = Path::new(filename);
    if !path.exists() {
        quit_with_error(&format!("file does not exist: {}", path.display()));
    }
    if !path.is_file() {
        quit_with_error(&format!("{} is not a file", path.display()));
    }
}


pub fn check_if_dir_exists(dir: &Path) {
    // Quits with an error if the given path is not an existing directory.
    let path = Path::new(dir);
    if !path.exists() {
        quit_with_error(&format!("directory does not exist: {}", path.display()));
    }
    if !path.is_dir() {
        quit_with_error(&format!("{} is not a directory", path.display()));
    }
}


pub fn check_if_dir_is_not_dir(dir: &Path) {
    // Quits with an error if the given path exists but is not a directory (not existing is okay).
    if dir.exists() && !dir.is_dir() {
        quit_with_error(&format!("{} exists but is not a directory", dir.display()));
    }
}


#[cfg(not(test))]
pub fn quit_with_error(text: &str) -> ! {
    // For friendly error messages, this function normally just prints the error and quits.
    eprintln!();
    eprintln!("Error: {}", text);
    std::process::exit(1);
}
#[cfg(test)]
pub fn quit_with_error(text: &str) -> ! {
    // But when running unit tests, this function instead panics so I can catch it for the test.
    panic!("{}", text);
}


pub fn load_fasta(filename: &Path) -> Vec<(String, String, String)> {
    // This function loads a FASTA file and runs a few checks on the result. If everything looks
    // good, it returns a vector of name+sequence tuples.
    if is_file_empty(filename) {
        quit_with_error(&format!("{} is an empty file", filename.display()));
    }
    let load_result = if is_file_gzipped(filename) {
        load_fasta_gzipped(filename)
    } else {
        load_fasta_not_gzipped(filename)
    };
    match load_result {
        Ok(_)  => (),
        Err(e) => quit_with_error(&format!("unable to load {}\n{}", filename.display(), e)),
    }
    let fasta_seqs = load_result.unwrap();
    check_load_fasta(&fasta_seqs, filename);
    fasta_seqs
}


fn check_load_fasta(fasta_seqs: &Vec<(String, String, String)>, filename: &Path) {
    // This function looks at the result of the load_fasta function and does some checks to make
    // sure everything looks okay. If any problems are found, it will quit with an error message.
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
    let mut set = HashSet::new();
    for (name, _, _) in fasta_seqs {
        if !set.insert(name) {
            quit_with_error(&format!("{} has a duplicate name: {}", filename.display(), name));
        }
    }
}


pub fn fastq_reader(fastq_file: &Path)
        -> seq_io::fastq::Reader<BufReader<Box<dyn std::io::Read>>> {
    // Returns a reader for a FASTQ file that works on both unzipped and gzipped files.
    let file = File::open(fastq_file).expect("Error opening file");
    let reader: Box<dyn Read> = if is_file_gzipped(fastq_file) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Reader::new(BufReader::new(reader))
}


fn is_file_empty(filename: &Path) -> bool {
    match fs::metadata(filename) {
        Ok(metadata) => metadata.len() == 0,
        Err(_) => false,
    }
}


fn is_file_gzipped(filename: &Path) -> bool {
    // This function returns true if the file appears to be gzipped (based on the first two bytes)
    // and false if not. If it can't open the file or read the first two bytes, it will quit with
    // an error message.
    let open_result = File::open(filename);
    match open_result {
        Ok(_)  => (),
        Err(e) => quit_with_error(&format!("unable to open {}\n{}", filename.display(), e)),
    }
    let file = open_result.unwrap();
    let mut reader = BufReader::new(file);
    let mut buf = vec![0u8; 2];
    let read_result = reader.read_exact(&mut buf);
    match read_result {
        Ok(_)  => (),
        Err(e) => quit_with_error(&format!("{} is too small\n{}", filename.display(), e)),
    }
    buf[0] == 31 && buf[1] == 139
}


fn load_fasta_not_gzipped(filename: &Path) -> io::Result<Vec<(String, String, String)>> {
    let mut fasta_seqs = Vec::new();
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut name = String::new();
    let mut header = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let text = line?;
        if text.is_empty() {continue;}
        if let Some(text) = text.strip_prefix('>') {
            if !name.is_empty() {
                sequence.make_ascii_uppercase();
                fasta_seqs.push((name, header, sequence));
                sequence = String::new();
            }
            header = text.to_string();
            let first_piece = text.split_whitespace().next();
            match first_piece {
                Some(_) => (),
                None    => quit_with_error(&format!("{} is not correctly formatted", filename.display())),
            }
            name = first_piece.unwrap().to_string();
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


fn load_fasta_gzipped(filename: &Path) -> io::Result<Vec<(String, String, String)>> {
    let mut fasta_seqs = Vec::new();
    let file = File::open(filename)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let mut name = String::new();
    let mut header = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let text = line?;
        if text.is_empty() {continue;}
        if let Some(text) = text.strip_prefix('>') {
            if !name.is_empty() {
                sequence.make_ascii_uppercase();
                fasta_seqs.push((name, header, sequence));
                sequence = String::new();
            }
            header = text.to_string();
            let first_piece = text.split_whitespace().next();
            match first_piece {
                Some(_) => (),
                None    => quit_with_error(&format!("{} is not correctly formatted", filename.display())),
            }
            name = first_piece.unwrap().to_string();
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
    let mut rev_seq: Vec<u8> = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        rev_seq.push(complement_base(b));
    }
    rev_seq
}


pub fn format_duration(duration: std::time::Duration) -> String {
    let microseconds = duration.as_micros() % 1000000;
    let seconds =      duration.as_micros() / 1000000 % 60;
    let minutes =      duration.as_micros() / 1000000 / 60 % 60;
    let hours =        duration.as_micros() / 1000000 / 60 / 60;
    format!("{}:{:02}:{:02}.{:06}", hours, minutes, seconds, microseconds)
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
    let mut formatted = format!("{:.6}", num);
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
        format!("{}", rounded_value)
    }
}


pub fn median_usize(values: &[usize]) -> usize {
    if values.is_empty() { return 0; }
    let mut sorted_values = values.to_vec();
    sorted_values.sort();
    let len = sorted_values.len();
    if len % 2 == 0 { (sorted_values[len / 2 - 1] + sorted_values[len / 2]) / 2 }
               else { sorted_values[len / 2] }
}


pub fn median_isize(values: &[isize]) -> isize {
    if values.is_empty() { return 0; }
    let mut sorted_values = values.to_vec();
    sorted_values.sort();
    let len = sorted_values.len();
    if len % 2 == 0 { (sorted_values[len / 2 - 1] + sorted_values[len / 2]) / 2 }
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
        pb.set_message(message.to_string().clone());
        pb
    }
}


pub fn reverse_path(path: &[i32]) -> Vec<i32> {
    path.iter().rev().map(|&num| -num).collect()
}


pub fn sign_at_end(num: i32) -> String {
    if num >= 0 {
        format!("{}+", num.abs())
    } else {
        format!("{}-", num.abs())
    }
}


pub fn sign_at_end_vec(nums: &[i32]) -> String {
    nums.iter().map(|&n| sign_at_end(n)).collect::<Vec<_>>().join(",")
}


pub fn up_to_first_space(string: &str) -> String {
    string.split_whitespace().next().unwrap_or("").to_string()
}


pub fn after_first_space(string: &str) -> String {
    let mut parts = string.splitn(2, char::is_whitespace);
    parts.next();
    parts.next().unwrap_or("").to_string()
}


pub fn first_char_in_file(filename: &Path) -> io::Result<char> {
    if is_file_gzipped(filename) {
        first_char_in_file_gzipped(filename)
    } else {
        first_char_in_file_not_gzipped(filename)
    }
}


fn first_char_in_file_not_gzipped(filename: &Path) -> io::Result<char> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    first_non_empty_char(reader)
}


fn first_char_in_file_gzipped(filename: &Path) -> io::Result<char> {
    let file = File::open(filename)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    first_non_empty_char(reader)
}


fn first_non_empty_char<R: BufRead>(reader: R) -> io::Result<char> {
    for line in reader.lines() {
        let text = line?;
        if !text.is_empty() {
            return text.chars().next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Empty line"));
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
    use tempfile::tempdir;

    use crate::tests::{make_test_file, make_gzipped_test_file};

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

        make_test_file(&filename, ">a\nACGT\n");
        assert_eq!(first_char_in_file(&filename).unwrap(), '>');

        make_test_file(&filename, "XYZ");
        assert_eq!(first_char_in_file(&filename).unwrap(), 'X');
    }

    #[test]
    fn test_first_char_in_file_gzipped() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");

        make_gzipped_test_file(&filename, ">a\nACGT\n");
        assert_eq!(first_char_in_file(&filename).unwrap(), '>');

        make_gzipped_test_file(&filename, "XYZ");
        assert_eq!(first_char_in_file(&filename).unwrap(), 'X');
    }

    #[test]
    fn test_load_fasta() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");

        make_test_file(&filename, ">a\nACGT\n>b xyz\nACGT\nACGT\n");
        let fasta = load_fasta_not_gzipped(&filename).unwrap();
        assert_eq!(fasta, vec![("a".to_string(), "a".to_string(), "ACGT".to_string()),
                               ("b".to_string(), "b xyz".to_string(), "ACGTACGT".to_string())]);

        make_gzipped_test_file(&filename, ">a\nACGT\n>b xyz\nACGT\nACGT\n");
        let fasta = load_fasta_gzipped(&filename).unwrap();
        assert_eq!(fasta, vec![("a".to_string(), "a".to_string(), "ACGT".to_string()),
                               ("b".to_string(), "b xyz".to_string(), "ACGTACGT".to_string())]);
    }

    #[test]
    fn test_is_file_empty() {
        let dir = tempdir().unwrap();
        let filename = dir.path().join("temp.fasta");

        make_test_file(&filename, "");
        assert!(is_file_empty(&filename));

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
}
