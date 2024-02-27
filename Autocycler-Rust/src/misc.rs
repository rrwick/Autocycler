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

use flate2::read::GzDecoder;

use std::collections::HashSet;
use std::fs::{File, read_dir};
use std::io;
use std::io::{prelude::*, BufReader};
use std::path::{Path, PathBuf};


// This module lets me use strand::FORWARD for true and strand::REVERSE for false.
pub mod strand {
    pub const FORWARD: bool = true;
    pub const REVERSE: bool = false;
}


pub fn find_all_assemblies(in_dir: &PathBuf) -> Vec<PathBuf> {
    let paths = match read_dir(in_dir) {
        Ok(paths) => paths,
        Err(_) => {
            quit_with_error(&format!("unable to read directory {}", in_dir.display()));
            return Vec::new();  // unreachable
        },
    };

    let mut all_assemblies: Vec<PathBuf> = Vec::new();
    for path in paths {
        let path = path.unwrap().path();
        if path.is_file() && 
           (path.extension().unwrap_or_default() == "fasta" ||
            path.extension().unwrap_or_default() == "fna" ||
            path.extension().unwrap_or_default() == "fa" ||
            path.extension().unwrap_or_default() == "gz" &&
            path.file_stem().unwrap_or_default().to_str().unwrap_or_default().ends_with(".fasta") ||
            path.file_stem().unwrap_or_default().to_str().unwrap_or_default().ends_with(".fna") ||
            path.file_stem().unwrap_or_default().to_str().unwrap_or_default().ends_with(".fa")) {
                all_assemblies.push(path);
        }
    }

    all_assemblies.sort_unstable();

    if all_assemblies.is_empty() {
        quit_with_error(&format!("Error: no assemblies found in {}", in_dir.display()));
    }

    all_assemblies
}


pub fn check_if_file_exists(filename: &PathBuf) {
    if !Path::new(filename).exists() {
        quit_with_error(&format!("{:?} file does not exist", filename));
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


fn complement_base_u8(base: u8) -> u8 {
    match base {
        b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G',
        _ => b'N'
    }
}


pub fn reverse_complement_u8(seq: &[u8]) -> Vec<u8> {
    let mut rev_seq: Vec<u8> = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        rev_seq.push(complement_base_u8(b));
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


#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_reverse_complement_u8() {
        assert_eq!(reverse_complement_u8(b"GGTATCACTCAGGAAGC"), b"GCTTCCTGAGTGATACC");
    }
}
