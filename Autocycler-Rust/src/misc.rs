// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

// This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use flate2::read::GzDecoder;

use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::io::{prelude::*, BufReader};
use std::path::{Path, PathBuf};


pub fn check_if_file_exists(filename: &PathBuf) {
    if !Path::new(filename).exists() {
        let error_message = format!("{:?} file does not exist", filename);
        quit_with_error(&error_message);
    }
}


pub fn quit_with_error(text: &str) {
    eprintln!();
    eprintln!("Error: {}", text);
    std::process::exit(1);
}


/// This function loads a FASTA file and runs a few checks on the result. If everything looks good,
/// it returns a vector of name+sequence tuples.
pub fn load_fasta(filename: &PathBuf) -> Vec<(String, String)> {
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
fn check_load_fasta(fasta_seqs: &Vec<(String, String)>, filename: &PathBuf) {
    if fasta_seqs.len() == 0 {
        quit_with_error(&format!("{:?} contains no sequences", filename));
    }
    for (name, sequence) in fasta_seqs {
        if name.len() == 0 {
            quit_with_error(&format!("{:?} has an unnamed sequence", filename));
        }
        if sequence.len() == 0 {
            quit_with_error(&format!("{:?} has an empty sequence", filename));
        }
    }
    let mut set = HashSet::new();
    for (name, _) in fasta_seqs {
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


fn load_fasta_not_gzipped(filename: &PathBuf) -> io::Result<Vec<(String, String)>> {
    let mut fasta_seqs = Vec::new();
    let file = File::open(&filename)?;
    let reader = BufReader::new(file);
    let mut name = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let text = line?;
        if text.len() == 0 {continue;}
        if text.starts_with('>') {
            if name.len() > 0 {
                sequence.make_ascii_uppercase();
                fasta_seqs.push((name, sequence));
                sequence = String::new();
            }
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
        fasta_seqs.push((name, sequence));
    }
    Ok(fasta_seqs)
}


fn load_fasta_gzipped(filename: &PathBuf) -> io::Result<Vec<(String, String)>> {
    let mut fasta_seqs = Vec::new();
    let file = File::open(&filename)?;
    let reader = BufReader::new(GzDecoder::new(file));
    let mut name = String::new();
    let mut sequence = String::new();
    for line in reader.lines() {
        let text = line?;
        if text.len() == 0 {continue;}
        if text.starts_with('>') {
            if name.len() > 0 {
                sequence.make_ascii_uppercase();
                fasta_seqs.push((name, sequence));
                sequence = String::new();
            }
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
        fasta_seqs.push((name, sequence));
    }
    Ok(fasta_seqs)
}


fn complement_base(base: char) -> char {
    match base {
        'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G',
        'a' => 't', 't' => 'a', 'g' => 'c', 'c' => 'g',
        'N' => 'N', 'n' => 'n',
        'R' => 'Y', 'Y' => 'R', 'S' => 'S', 'W' => 'W', 'K' => 'M', 'M' => 'K',
        'B' => 'V', 'V' => 'B', 'D' => 'H', 'H' => 'D',
        'r' => 'y', 'y' => 'r', 's' => 's', 'w' => 'w', 'k' => 'm', 'm' => 'k',
        'b' => 'v', 'v' => 'b', 'd' => 'h', 'h' => 'd',
        '.' => '.', '-' => '-', '?' => '?',
        _ => 'N'
    }
}


pub fn reverse_complement(seq: &str) -> String {
    let mut rev_seq: String = String::with_capacity(seq.len());
    for c in seq.chars().rev() {
        rev_seq.push(complement_base(c))
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


/// This function implements banker's rounding (i.e. round-half-to-even) for positive numbers. I
/// wrote it so I could replicate Python's rounding behaviour, because Rust's round function has
/// round-half-up behaviour. I had tried using math::round::half_to_even, but that didn't seem to
/// work correctly (rounded 42.55 to 42).
pub fn bankers_rounding(float: f64) -> u32 {
    let fractional_part = float - float.floor();
    let rounded_down = float as u32;
    if fractional_part < 0.5 {
        return rounded_down;
    } else if fractional_part > 0.5 {
        return rounded_down + 1;
    } else {  // fractional_part == 0.5
        if rounded_down % 2 == 0 {  // is even
            return rounded_down;
        } else {
            return rounded_down + 1;
        }
    }
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
    fn test_bankers_rounding() {
        assert_eq!(bankers_rounding(0.0), 0);
        assert_eq!(bankers_rounding(123.0), 123);
        assert_eq!(bankers_rounding(98765.0), 98765);

        assert_eq!(bankers_rounding(0.4999), 0);
        assert_eq!(bankers_rounding(0.5),    0);
        assert_eq!(bankers_rounding(0.5001), 1);

        assert_eq!(bankers_rounding(42.45), 42);
        assert_eq!(bankers_rounding(42.5),  42);
        assert_eq!(bankers_rounding(42.55), 43);

        assert_eq!(bankers_rounding(12345.4998), 12345);
        assert_eq!(bankers_rounding(12345.5),    12346);
        assert_eq!(bankers_rounding(12345.5002), 12346);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("GGTATCACTCAGGAAGC"), "GCTTCCTGAGTGATACC");
        assert_eq!(reverse_complement("GGGGaaaaaaaatttatatat"), "atatataaattttttttCCCC");
        assert_eq!(reverse_complement("atatataaattttttttCCCC"), "GGGGaaaaaaaatttatatat");
        assert_eq!(reverse_complement("ACGT123"), "NNNACGT");
    }
}
