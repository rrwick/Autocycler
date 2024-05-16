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
use flate2::read::GzDecoder;
use std::collections::HashSet;
use std::fs::{File, read_dir, create_dir_all};
use std::io;
use std::io::{prelude::*, BufReader};
use std::path::{Path, PathBuf};
use std::time::Duration;


// This module lets me use strand::FORWARD for true and strand::REVERSE for false.
pub mod strand {
    pub const FORWARD: bool = true;
    pub const REVERSE: bool = false;
}


pub fn create_dir(out_dir: &PathBuf) {
    match create_dir_all(&out_dir) {
        Ok(_) => {},
        Err(e) => quit_with_error(&format!("failed to create output directory\n{:?}", e)),
    }
}


pub fn find_all_assemblies(in_dir: &PathBuf) -> Vec<PathBuf> {
    let paths = match read_dir(in_dir) {
        Ok(paths) => paths,
        Err(_) => {
            quit_with_error(&format!("unable to read directory {}", in_dir.display()));
            unreachable!()
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
        quit_with_error(&format!("Error: no assemblies found in {}", in_dir.display()));
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


pub fn check_if_file_exists(filename: &PathBuf) {
    // Quits with an error if the given path is not an existing file.
    let path = Path::new(filename);
    if !path.exists() {
        quit_with_error(&format!("{:?} file does not exist", path));
    }
    if !path.is_file() {
        quit_with_error(&format!("{:?} is not a file", path));
    }
}


pub fn check_if_dir_exists(dir: &PathBuf) {
    // Quits with an error if the given path is not an existing directory.
    let path = Path::new(dir);
    if !path.exists() {
        quit_with_error(&format!("{:?} directory does not exist", path));
    }
    if !path.is_dir() {
        quit_with_error(&format!("{:?} is not a directory", path));
    }
}


pub fn check_if_dir_is_not_dir(dir: &PathBuf) {
    // Quits with an error if the given path exists but is not a directory (not existing is okay).
    if dir.exists() && !dir.is_dir() {
        quit_with_error(&format!("{:?} exists but is not a directory", dir));
    }
}


pub fn quit_with_error(text: &str) {
    eprintln!();
    eprintln!("Error: {}", text);
    std::process::exit(1);
}


/// This function loads a FASTA file and runs a few checks on the result. If everything looks good,
/// it returns a vector of name+sequence tuples.
pub fn load_fasta(filename: &PathBuf) -> Vec<(String, String, String)> {
    let load_result = if is_file_gzipped(&filename) {
        load_fasta_gzipped(&filename)
    } else {
        load_fasta_not_gzipped(&filename)
    };
    match load_result {
        Ok(_)  => (),
        Err(_) => quit_with_error(&format!("unable to load {:?}", filename)),
    }
    let fasta_seqs = load_result.unwrap();
    check_load_fasta(&fasta_seqs, &filename);
    fasta_seqs
}


/// This function looks at the result of the load_fasta function and does some checks to make sure
/// everything looks okay. If any problems are found, it will quit with an error message.
fn check_load_fasta(fasta_seqs: &Vec<(String, String, String)>, filename: &PathBuf) {
    if fasta_seqs.len() == 0 {
        quit_with_error(&format!("{:?} contains no sequences", filename));
    }
    for (name, _, sequence) in fasta_seqs {
        if name.len() == 0 {
            quit_with_error(&format!("{:?} has an unnamed sequence", filename));
        }
        if sequence.len() == 0 {
            quit_with_error(&format!("{:?} has an empty sequence", filename));
        }
    }
    let mut set = HashSet::new();
    for (name, _, _) in fasta_seqs {
        set.insert(name);
    }
    if set.len() < fasta_seqs.len() {
        quit_with_error(&format!("{:?} has a duplicated name", filename));
    }
}


/// This function returns true if the file appears to be gzipped (based on the first two bytes) and
/// false if not. If it can't open the file or read the first two bytes, it will quit with an error
/// message.
fn is_file_gzipped(filename: &PathBuf) -> bool {
    let open_result = File::open(&filename);
    match open_result {
        Ok(_)  => (),
        Err(_) => quit_with_error(&format!("unable to open {:?}", filename)),
    }
    let file = open_result.unwrap();

    let mut reader = BufReader::new(file);
    let mut buf = vec![0u8; 2];

    let read_result = reader.read_exact(&mut buf);
    match read_result {
        Ok(_)  => (),
        Err(_) => quit_with_error(&format!("{:?} is too small", filename)),
    }

    buf[0] == 31 && buf[1] == 139
}


fn load_fasta_not_gzipped(filename: &PathBuf) -> io::Result<Vec<(String, String, String)>> {
    let mut fasta_seqs = Vec::new();
    let file = File::open(&filename)?;
    let reader = BufReader::new(file);
    let mut name = String::new();
    let mut header = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let text = line?;
        if text.len() == 0 {continue;}
        if text.starts_with('>') {
            if name.len() > 0 {
                sequence.make_ascii_uppercase();
                fasta_seqs.push((name, header, sequence));
                sequence = String::new();
            }
            header = (&text[1..]).to_string();
            let first_piece = text[1..].split_whitespace().next();
            match first_piece {
                Some(_) => (),
                None    => quit_with_error(&format!("{:?} is not correctly formatted", filename)),
            }
            name = first_piece.unwrap().to_string();
        } else {
            if name.len() == 0 {
                quit_with_error(&format!("{:?} is not correctly formatted", filename));
            }
            sequence.push_str(&text);
        }
    }
    if name.len() > 0 {
        sequence.make_ascii_uppercase();
        fasta_seqs.push((name, header, sequence));
    }
    Ok(fasta_seqs)
}


fn load_fasta_gzipped(filename: &PathBuf) -> io::Result<Vec<(String, String, String)>> {
    let mut fasta_seqs = Vec::new();
    let file = File::open(&filename)?;
    let reader = BufReader::new(GzDecoder::new(file));
    let mut name = String::new();
    let mut header = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let text = line?;
        if text.len() == 0 {continue;}
        if text.starts_with('>') {
            if name.len() > 0 {
                sequence.make_ascii_uppercase();
                fasta_seqs.push((name, header, sequence));
                sequence = String::new();
            }
            header = (&text[1..]).to_string();
            let first_piece = header.split_whitespace().next();
            match first_piece {
                Some(_) => (),
                None    => quit_with_error(&format!("{:?} is not correctly formatted", filename)),
            }
            name = first_piece.unwrap().to_string();
        } else {
            if name.len() == 0 {
                quit_with_error(&format!("{:?} is not correctly formatted", filename));
            }
            sequence.push_str(&text);
        }
    }
    if name.len() > 0 {
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
    if !formatted.contains('.') {
        return formatted
    }
    while formatted.chars().last().unwrap() == '0' {
        formatted.pop();
    }
    if formatted.chars().last().unwrap() == '.' {
        formatted.pop();
    }
    formatted
}


pub fn round_float(num: f64, digits: u32) -> f64 {
    // Rounds a float to the given number of digits.
    let multiplier = 10f64.powi(digits as i32);
    (num * multiplier).round() / multiplier
}


pub fn median_usize(values: &[usize]) -> usize {
    if values.is_empty() {
        return 0;
    }
    let mut sorted_values = values.to_vec();
    sorted_values.sort();
    let len = sorted_values.len();
    if len % 2 == 0 {
        (sorted_values[len / 2 - 1] + sorted_values[len / 2]) / 2
    } else {
        sorted_values[len / 2]
    }
}


pub fn median_f64(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted_values = values.to_vec();
    sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let len = sorted_values.len();
    if len % 2 == 0 {
        (sorted_values[len / 2 - 1] + sorted_values[len / 2]) / 2.0
    } else {
        sorted_values[len / 2]
    }
}


pub fn spinner(message: &str) -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    pb.enable_steady_tick(Duration::from_millis(100));
    pb.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&vec!["⠋", "⠙", "⠚", "⠞", "⠖", "⠦", "⠴", "⠲", "⠳", "⠓"])  // dots3 from github.com/sindresorhus/cli-spinners 
            .template("{spinner} {msg}").unwrap(),
    );
    pb.set_message(message.to_string().clone());
    pb
}


#[cfg(test)]
mod tests {
    use super::*;

    fn assert_almost_eq(a: f64, b: f64, epsilon: f64) {
        assert!((a - b).abs() < epsilon, "Numbers are not within {:?} of each other: {} vs {}", epsilon, a, b);
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
    fn test_round_float() {
        assert_almost_eq(round_float(0.1234, 1), 0.1, 1e-8);
        assert_almost_eq(round_float(0.1234, 2), 0.12, 1e-8);
        assert_almost_eq(round_float(0.1234, 3), 0.123, 1e-8);
        assert_almost_eq(round_float(0.1234, 4), 0.1234, 1e-8);
        assert_almost_eq(round_float(0.1234, 5), 0.1234, 1e-8);
        assert_almost_eq(round_float(0.1234, 6), 0.1234, 1e-8);
        assert_almost_eq(round_float(0.9876, 1), 1.0, 1e-8);
        assert_almost_eq(round_float(0.9876, 2), 0.99, 1e-8);
        assert_almost_eq(round_float(0.9876, 3), 0.988, 1e-8);
        assert_almost_eq(round_float(0.9876, 4), 0.9876, 1e-8);
        assert_almost_eq(round_float(0.9876, 5), 0.9876, 1e-8);
        assert_almost_eq(round_float(0.9876, 6), 0.9876, 1e-8);
    }

    #[test]
    fn test_median_usize() {
        assert_eq!(median_usize(&mut vec![]), 0);
        assert_eq!(median_usize(&mut vec![0, 1, 2, 3, 4]), 2);
        assert_eq!(median_usize(&mut vec![4, 3, 2, 1, 0]), 2);
        assert_eq!(median_usize(&mut vec![0, 1, 2, 3, 4, 5]), 2);
        assert_eq!(median_usize(&mut vec![5, 4, 3, 2, 1, 0]), 2);
        assert_eq!(median_usize(&mut vec![0, 2, 4, 6, 8, 10]), 5);
        assert_eq!(median_usize(&mut vec![10, 8, 6, 4, 2, 0]), 5);
    }

    #[test]
    fn test_median_f64() {
        assert_almost_eq(median_f64(&mut vec![]), 0.0, 1e-8);
        assert_almost_eq(median_f64(&mut vec![0.0, 1.0, 2.0, 3.0, 4.0]), 2.0, 1e-8);
        assert_almost_eq(median_f64(&mut vec![4.0, 3.0, 2.0, 1.0, 0.0]), 2.0, 1e-8);
        assert_almost_eq(median_f64(&mut vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]), 2.5, 1e-8);
        assert_almost_eq(median_f64(&mut vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0]), 2.5, 1e-8);
        assert_almost_eq(median_f64(&mut vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]), 5.0, 1e-8);
        assert_almost_eq(median_f64(&mut vec![10.0, 8.0, 6.0, 4.0, 2.0, 0.0]), 5.0, 1e-8);
    }
}
