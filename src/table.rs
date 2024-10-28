// This file contains the code for the autocycler table subcommand.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use std::fs;
use std::path::PathBuf;

use crate::misc::quit_with_error;


pub fn table(autocycler_dir: Option<PathBuf>, name: String, fields: String) {
    let fields = parse_fields(fields);
    if autocycler_dir.is_none() {
        print_header(fields);
    } else {
        print_values(autocycler_dir.unwrap(), name, fields);
    }
}


fn parse_fields(comma_delimited_fields: String) -> Vec<String> {
    let fields = comma_delimited_fields.split(',').map(|s| s.to_string()).collect();

    // TODO: check to make sure each field is valid

    fields
}


fn print_header(fields: Vec<String>) {
    println!("name\t{}", fields.join("\t"));
}


fn print_values(autocycler_dir: PathBuf, name: String, fields: Vec<String>) {
    if name.contains('\t') {
        quit_with_error("--name cannot contain tab characters")
    }
    print!("{}", name);

    let yaml_files = find_all_yaml_files(&autocycler_dir);
    eprintln!("{:?}", yaml_files);  // TEMP

    // TODO
    // TODO
    // TODO
    // TODO

    println!();
}


fn find_all_yaml_files(autocycler_dir: &PathBuf) -> Vec<PathBuf> {
    let mut yaml_files = Vec::new();
    visit_dirs(autocycler_dir, &mut yaml_files);
    yaml_files
}


fn visit_dirs(dir: &PathBuf, yaml_files: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                visit_dirs(&path, yaml_files);
            } else if path.extension().map_or(false, |ext| ext == "yaml") {
                yaml_files.push(path);
            }
        }
    }
}
