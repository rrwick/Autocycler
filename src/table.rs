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
use std::path::PathBuf;

use crate::metrics::{InputAssemblyMetrics, ClusteringMetrics, CombineMetrics};
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
    let fields = comma_delimited_fields.replace(" ", "").split(',')
                                       .map(|s| s.to_string()).collect();
    let mut valid_fields = HashSet::new();
    valid_fields.extend(InputAssemblyMetrics::get_field_names());
    valid_fields.extend(ClusteringMetrics::get_field_names());
    valid_fields.extend(CombineMetrics::get_field_names());
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


fn print_values(autocycler_dir: PathBuf, name: String, fields: Vec<String>) {
    if name.contains('\t') {
        quit_with_error("--name cannot contain tab characters")
    }
    print!("{}", name);

    let yaml_files = find_all_yaml_files(&autocycler_dir);
    let input_assemblies_yaml = get_one_copy_yaml(&yaml_files, "input_assemblies.yaml");
    let clustering_yaml = get_one_copy_yaml(&yaml_files, "clustering.yaml");
    let consensus_assembly_yaml = get_one_copy_yaml(&yaml_files, "consensus_assembly.yaml");

    let mut metrics_map: HashMap<String, Value> = HashMap::new();
    if let Some(path) = input_assemblies_yaml   { metrics_map.extend(load_yaml_to_map(&path)); }
    if let Some(path) = clustering_yaml         { metrics_map.extend(load_yaml_to_map(&path)); }
    if let Some(path) = consensus_assembly_yaml { metrics_map.extend(load_yaml_to_map(&path)); }

    for field in fields {
        print!("\t");
        if let Some(value) = metrics_map.get(&field) {
            print!("{}", format_value(value));
        }
    }
    println!();
}


fn load_yaml_to_map(yaml_path: &PathBuf) -> HashMap<String, Value> {
    let content = fs::read_to_string(yaml_path)
        .unwrap_or_else(|_| quit_with_error("Could not read YAML file"));
    serde_yaml::from_str(&content)
        .unwrap_or_else(|_| quit_with_error("Failed to parse YAML file"))
}


fn find_all_yaml_files(autocycler_dir: &PathBuf) -> Vec<PathBuf> {
    let mut yaml_files = Vec::new();
    visit_dirs_for_yaml_files(autocycler_dir, &mut yaml_files);
    yaml_files
}


fn visit_dirs_for_yaml_files(dir: &PathBuf, yaml_files: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                visit_dirs_for_yaml_files(&path, yaml_files);
            } else if path.extension().map_or(false, |ext| ext == "yaml") {
                yaml_files.push(path);
            }
        }
    }
}


fn get_one_copy_yaml(yaml_files: &Vec<PathBuf>, filename: &str) -> Option<PathBuf> {
    let found_files = yaml_files.iter()
        .filter(|path| path.file_name().map_or(false, |name| name == filename))
        .collect::<Vec<_>>();
    match found_files.len() {
        0 => None,
        1 => Some(found_files[0].clone()),
        _ => quit_with_error(&format!("Multiple {} files found", filename)),
    }
}


fn format_value(value: &Value) -> String {
    // This function formats serde_yaml::Value types. Sequences are formatted with square brackets
    // and commas (no spaces). Mappings are formatted with curly brackets, colons and commas (no
    // spaces).
    match value {
        Value::Number(n) => n.to_string(),
        Value::String(s) => s.clone(),
        Value::Bool(b) => b.to_string(),
        Value::Sequence(s) =>
            format!("[{}]", s.iter().map(format_value).collect::<Vec<_>>().join(",")),
        Value::Mapping(m) =>
            format!("{{{}}}",
                    m.iter().map(|(k, v)| format!("{}:{}", format_value(k),
                                                  format_value(v))).collect::<Vec<_>>().join(",")),
        _ => String::new(),
    }
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
        assert_eq!(format_value(&Value::Number(serde_yaml::Number::from(12))), "12");
        assert_eq!(format_value(&Value::Number(serde_yaml::Number::from(1.2))), "1.2");
        assert_eq!(format_value(&Value::String("abc".to_string())), "abc");
        assert_eq!(format_value(&Value::Bool(true)), "true");
        assert_eq!(format_value(&Value::Bool(false)), "false");
    }

    #[test]
    fn test_format_value_sequence() {
        let v1 = Value::Number(serde_yaml::Number::from(12));
        let v2 = Value::Number(serde_yaml::Number::from(1.2));
        let v3 = Value::String("abc".to_string());
        let v4 = Value::Bool(true);
        let seq = Value::Sequence(vec![v1, v2, v3, v4]);
        assert_eq!(format_value(&seq), "[12,1.2,abc,true]");
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
        assert_eq!(format_value(&Value::Mapping(map)), "{12:1.2,abc:true}");
    }
}
