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

use serde_yaml::Value;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use crate::metrics::{ClusteringMetrics, CombineMetrics, InputAssemblyMetrics, SubsampleMetrics,
                     TrimmedClusterMetrics, UntrimmedClusterMetrics};
use crate::misc::{check_if_dir_exists, quit_with_error, format_float_sigfigs};


pub fn table(autocycler_dir: Option<PathBuf>, name: String, fields: String, sigfigs: usize) {
    check_settings(&autocycler_dir, sigfigs);
    let fields = parse_fields(fields);
    if let Some(autocycler_dir) = autocycler_dir {
        print_values(&autocycler_dir, name, fields, sigfigs);
    } else {
        print_header(fields);
    }
}


fn check_settings(autocycler_dir: &Option<PathBuf>, sigfigs: usize) {
    if let Some(dir) = autocycler_dir.as_ref() {
        check_if_dir_exists(dir);
    }
    if sigfigs == 0 {
        quit_with_error("--sigfigs must be 1 or greater");
    }
}


fn parse_fields(comma_delimited_fields: String) -> Vec<String> {
    let fields = comma_delimited_fields.replace(" ", "").split(',')
                                       .map(|s| s.to_string()).collect();
    let mut valid_fields = HashSet::new();
    valid_fields.extend(SubsampleMetrics::get_field_names());
    valid_fields.extend(InputAssemblyMetrics::get_field_names());
    valid_fields.extend(ClusteringMetrics::get_field_names());
    valid_fields.extend(CombineMetrics::get_field_names());
    valid_fields.extend(UntrimmedClusterMetrics::get_field_names());
    valid_fields.extend(TrimmedClusterMetrics::get_field_names());
    for field in &fields {
        if !valid_fields.contains(field) {
            quit_with_error(&format!("{} is not a valid field name", field));
        }
    }
    fields
}


fn print_header(fields: Vec<String>) {
    println!("name\t{}", fields.join("\t"));
}


fn print_values(autocycler_dir: &Path, name: String, fields: Vec<String>, sigfigs: usize) {
    if name.contains('\t') {
        quit_with_error("--name cannot contain tab characters")
    }
    print!("{}", name);

    let yaml_files = find_all_yaml_files(autocycler_dir);
    let subsample_yaml = get_one_copy_yaml(&yaml_files, "subsample.yaml");
    let input_assemblies_yaml = get_one_copy_yaml(&yaml_files, "input_assemblies.yaml");
    let clustering_yaml = get_one_copy_yaml(&yaml_files, "clustering.yaml");
    let consensus_assembly_yaml = get_one_copy_yaml(&yaml_files, "consensus_assembly.yaml");
    let untrimmed_yamls = get_multi_copy_yaml(&yaml_files, "1_untrimmed.yaml");
    let trimmed_yamls = get_multi_copy_yaml(&yaml_files, "2_trimmed.yaml");

    let mut map: HashMap<String, Value> = HashMap::new();
    if let Some(path) = subsample_yaml          { map.extend(load_single_yaml_to_map(&path)); }
    if let Some(path) = input_assemblies_yaml   { map.extend(load_single_yaml_to_map(&path)); }
    if let Some(path) = clustering_yaml         { map.extend(load_single_yaml_to_map(&path)); }
    if let Some(path) = consensus_assembly_yaml { map.extend(load_single_yaml_to_map(&path)); }
    if !untrimmed_yamls.is_empty() { map.extend(load_multi_yaml_to_map(&untrimmed_yamls)); }
    if !trimmed_yamls.is_empty()   { map.extend(load_multi_yaml_to_map(&trimmed_yamls)); }

    for field in fields {
        print!("\t");
        if let Some(value) = map.get(&field) {
            print!("{}", format_value(value, sigfigs));
        }
    }
    println!();
}


fn load_single_yaml_to_map(yaml_path: &Path) -> HashMap<String, Value> {
    let content = fs::read_to_string(yaml_path)
        .unwrap_or_else(|_| quit_with_error("Could not read YAML file"));
    serde_yaml::from_str(&content)
        .unwrap_or_else(|_| quit_with_error("Failed to parse YAML file"))
}


fn load_multi_yaml_to_map(yaml_paths: &Vec<PathBuf>) -> HashMap<String, Value> {
    let mut combined_map: HashMap<String, Vec<Value>> = HashMap::new();
    for yaml_path in yaml_paths {
        let file_map = load_single_yaml_to_map(yaml_path);
        for (key, value) in file_map {
            combined_map.entry(key).or_default().push(value);
        }
    }
    combined_map.into_iter().map(|(key, values)| (key, Value::Sequence(values))).collect()
}


fn find_all_yaml_files(autocycler_dir: &Path) -> Vec<PathBuf> {
    let mut yaml_files = Vec::new();
    visit_dirs_for_yaml_files(autocycler_dir, &mut yaml_files);
    yaml_files.sort();
    yaml_files
}


fn visit_dirs_for_yaml_files(dir: &Path, yaml_files: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                visit_dirs_for_yaml_files(&path, yaml_files);
            } else if path.extension().is_some_and(|ext| ext == "yaml") {
                yaml_files.push(path);
            }
        }
    }
}


fn get_one_copy_yaml(yaml_files: &[PathBuf], filename: &str) -> Option<PathBuf> {
    // Returns the YAML file in the given path with a matching filename. No match is okay and one
    // match is okay, but multiple matches will result in an error.
    let found_files = yaml_files.iter()
        .filter(|path| path.file_name().is_some_and(|name| name == filename)).collect::<Vec<_>>();
    match found_files.len() {
        0 => None,
        1 => Some(found_files[0].clone()),
        _ => quit_with_error(&format!("Multiple {} files found", filename)),
    }
}


fn get_multi_copy_yaml(yaml_files: &[PathBuf], filename: &str) -> Vec<PathBuf> {
    // Returns all YAML files in the given path with a matching filename, excluding those that are
    // in a qc_fail directory.
    yaml_files.iter().filter(|path| {
                         path.file_name().is_some_and(|name| name == filename) &&
                         !path.to_string_lossy().contains("/qc_fail/")
                     }).cloned().collect()
}


fn format_value(value: &Value, sigfigs: usize) -> String {
    // This function formats serde_yaml::Value types. Sequences are formatted with square brackets
    // and commas (no spaces). Mappings are formatted with curly brackets, colons and commas (no
    // spaces).
    match value {
        Value::Number(n) => format_number(n, sigfigs),
        Value::String(s) => s.clone(),
        Value::Bool(b) => b.to_string(),
        Value::Sequence(s) => format_sequence(s, sigfigs),
        Value::Mapping(m) => format_mapping(m, sigfigs),
        _ => String::new(),
    }
}


fn format_number(n: &serde_yaml::Number, sigfigs: usize) -> String {
    if n.is_i64() || n.is_u64()      { n.to_string() }
    else if let Some(f) = n.as_f64() { format_float_sigfigs(f, sigfigs) }
    else                             { n.to_string() }
}


fn format_sequence(s: &[Value], sigfigs: usize) -> String {
    format!("[{}]", s.iter().map(|v| format_value(v, sigfigs)).collect::<Vec<_>>().join(","))
}


fn format_mapping(m: &serde_yaml::Mapping, sigfigs: usize) -> String {
    format!("{{{}}}",
            m.iter().map(|(k, v)| format!("{}:{}",
                                          format_value(k, sigfigs),
                                          format_value(v, sigfigs))).collect::<Vec<_>>().join(","))
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;

    #[test]
    fn test_get_one_copy_yaml() {
        let yaml_files = vec![PathBuf::from("dir/a.yaml"), PathBuf::from("dir/b.yaml"),
                              PathBuf::from("dir/c.yaml"), PathBuf::from("dir2/b.yaml")];
        assert_eq!(get_one_copy_yaml(&yaml_files, "a.yaml"), Some(PathBuf::from("dir/a.yaml")));
        assert_eq!(get_one_copy_yaml(&yaml_files, "c.yaml"), Some(PathBuf::from("dir/c.yaml")));
        assert_eq!(get_one_copy_yaml(&yaml_files, "d.yaml"), None);
        assert!(panic::catch_unwind(|| {
            get_one_copy_yaml(&yaml_files, "b.yaml");
        }).is_err());
    }

    #[test]
    fn test_get_multi_copy_yaml() {
        let yaml_files = vec![PathBuf::from("dir/a.yaml"), PathBuf::from("dir/b.yaml"),
                              PathBuf::from("dir/qc_fail/c.yaml"), PathBuf::from("dir2/b.yaml")];
        assert_eq!(get_multi_copy_yaml(&yaml_files, "a.yaml"), vec![PathBuf::from("dir/a.yaml")]);
        assert_eq!(get_multi_copy_yaml(&yaml_files, "b.yaml"),
                   vec![PathBuf::from("dir/b.yaml"), PathBuf::from("dir2/b.yaml")]);
        let empty_vec: Vec<PathBuf> = Vec::new();
        assert_eq!(get_multi_copy_yaml(&yaml_files, "c.yaml"), empty_vec);

    }

    #[test]
    fn test_parse_fields() {
        assert_eq!(parse_fields("input_assemblies_count,fail_contig_count".to_string()),
                   vec!["input_assemblies_count", "fail_contig_count"]);
        assert_eq!(parse_fields("pass_contig_count,consensus_assembly_fully_resolved".to_string()),
                   vec!["pass_contig_count", "consensus_assembly_fully_resolved"]);
        assert!(panic::catch_unwind(|| {
            parse_fields("input_assemblies_count,abc".to_string());
        }).is_err());
    }

    #[test]
    fn test_format_value_simple() {
        assert_eq!(format_value(&Value::Number(serde_yaml::Number::from(12)), 2), "12");
        assert_eq!(format_value(&Value::Number(serde_yaml::Number::from(1.2)), 1), "1");
        assert_eq!(format_value(&Value::Number(serde_yaml::Number::from(1.2)), 2), "1.2");
        assert_eq!(format_value(&Value::Number(serde_yaml::Number::from(1.2)), 4), "1.200");
        assert_eq!(format_value(&Value::String("abc".to_string()), 2), "abc");
        assert_eq!(format_value(&Value::Bool(true), 2), "true");
        assert_eq!(format_value(&Value::Bool(false), 2), "false");
    }

    #[test]
    fn test_format_value_sequence() {
        let v1 = Value::Number(serde_yaml::Number::from(12));
        let v2 = Value::Number(serde_yaml::Number::from(1.2));
        let v3 = Value::String("abc".to_string());
        let v4 = Value::Bool(true);
        let seq = Value::Sequence(vec![v1, v2, v3, v4]);
        assert_eq!(format_value(&seq, 2), "[12,1.2,abc,true]");
    }

    #[test]
    fn test_format_value_mapping() {
        let v1 = Value::Number(serde_yaml::Number::from(12));
        let v2 = Value::Number(serde_yaml::Number::from(1.2));
        let v3 = Value::String("abc".to_string());
        let v4 = Value::Bool(true);
        let mut map = serde_yaml::Mapping::new();
        map.insert(v1, v2);
        map.insert(v3, v4);
        assert_eq!(format_value(&Value::Mapping(map), 2), "{12:1.2,abc:true}");
    }
}
