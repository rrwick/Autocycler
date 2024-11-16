// This file contains functions for outputting formatted text to stderr.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use chrono::prelude::*;
use colored::Colorize;


pub fn section_header(text: &str) {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
    let date = format!("({})", now);
    eprintln!();
    eprintln!("{} {}", text.bold().bright_yellow().underline(), date.dimmed());
}


pub fn explanation(text: &str) {
    let mut term_width = 80;
    if let Some((w, _)) = term_size::dimensions_stderr() {
        term_width = w;
    }
    let indented_text = format!("    {}", text);
    eprintln!("{}", textwrap::fill(&indented_text, term_width).dimmed());
    eprintln!();
}
