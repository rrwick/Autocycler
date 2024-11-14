// This file contains the code for making a cluster's all-vs-all dotplots.

// Copyright 2024 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Autocycler

// This file is part of Autocycler. Autocycler is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Autocycler
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Autocycler. If not, see <http://www.gnu.org/licenses/>.

use ab_glyph::{Font, FontArc, PxScale, ScaleFont};
use image::{RgbImage, Rgb, ImageBuffer};
use imageproc::drawing::{draw_filled_rect_mut, draw_text_mut, draw_hollow_rect_mut};
use imageproc::rect::Rect;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::log::{section_header, explanation};
use crate::misc::{first_char_in_file, quit_with_error, reverse_complement, spinner,
                  find_all_assemblies, load_fasta};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


// Some hard-coded settings. Values between 0 and 1 are relative to full image resolution.
static INITIAL_TOP_LEFT_GAP: f64 = 0.1;
static BORDER_GAP: f64 = 0.015;
static BETWEEN_SEQ_GAP: f64 = 0.01;
static TOTAL_BETWEEN_SEQ_GAP: f64 = 0.1;
static TEXT_GAP: f64 = 0.0025;
static MAX_FONT_SIZE: f64 = 0.025;
static BACKGROUND_COLOUR: image::Rgb<u8> = Rgb([255, 255, 255]);     // white
static SELF_VS_SELF_COLOUR: image::Rgb<u8> = Rgb([211, 211, 211]);   // lightgrey
static SELF_VS_OTHER_COLOUR: image::Rgb<u8> = Rgb([245, 245, 245]);  // whitesmoke
static TEXT_COLOUR: image::Rgb<u8> = Rgb([0, 0, 0]);                 // black
static OUTLINE_COLOUR: image::Rgb<u8> = Rgb([0, 0, 0]);              // black
static FORWARD_DOT_COLOUR: image::Rgb<u8> = Rgb([0, 0, 205]);        // mediumblue
static REVERSE_DOT_COLOUR: image::Rgb<u8> = Rgb([178, 34, 34]);      // firebrick


pub fn dotplot(input: PathBuf, out_png: PathBuf, res: u32, kmer: u32) {
    check_settings(res, kmer);
    let input_type = determine_input_type(&input);
    starting_message();
    print_settings(&input, res, kmer);
    let seqs = load_sequences(&input, input_type);
    create_dotplot(&seqs, &out_png, res, kmer);
    finished_message(&out_png);
}


fn check_settings(res: u32, kmer: u32) {
    if res < 500     { quit_with_error("--res cannot be less than 500"); }
    if res > 10000   { quit_with_error("--res cannot be greater than 10000"); }
    if kmer < 10     { quit_with_error("--kmer cannot be less than 10"); }
    if kmer > 100    { quit_with_error("--kmer cannot be greater than 100"); }
}


enum InputType { Gfa, Fasta, Directory }

fn determine_input_type(input: &Path) -> InputType {
    if input.is_dir() {
        return InputType::Directory;
    }
    if !input.is_file() {
        quit_with_error("--input is neither a file nor a directory");
    }
    let first_char = first_char_in_file(input).unwrap();
    match first_char {
        '>' => InputType::Fasta,
        'H' | 'S' => InputType::Gfa,
        _ => {
            quit_with_error("--input is neither GFA or FASTA");
        }
    }
}


fn starting_message() {
    section_header("Starting autocycler dotplot");
    explanation("This command will take a unitig graph (either before or after trimming) and \
                 generate a dotplot image containing all pairwise comparisons of the sequences.");
}


fn print_settings(input: &Path, res: u32, kmer: u32) {
    eprintln!("Settings:");
    eprintln!("  --input {}", input.display());
    eprintln!("  --res {}", res);
    eprintln!("  --kmer {}", kmer);
    eprintln!();
}


fn finished_message(out_png: &Path) {
    section_header("Finished!");
    eprintln!("Pairwise dotplots: {}", out_png.display());
    eprintln!();
}


#[derive(Eq, Hash, PartialEq, Clone)]
struct FileSeqName {
    filename: String,
    seqname: String,
}


fn load_sequences(input: &Path, input_type: InputType) -> Vec<(FileSeqName, Vec<u8>)> {
    let seqs = match input_type {
        InputType::Gfa => {
            let (graph, sequences) = load_from_graph(input);
            let reconstructed = graph.reconstruct_original_sequences_u8(&sequences);
            reconstructed.into_iter().map(|((filename, seqname), seq)|
                                          (FileSeqName { filename, seqname }, seq)).collect()
        },
        InputType::Fasta => {
            load_from_fasta(input)
        }
        InputType::Directory => {
            load_from_directory(input)
        }
    };
    if seqs.is_empty() {
        quit_with_error("no sequences were loaded")
    }
    seqs
}


fn load_from_graph(gfa: &Path) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading sequences");
    explanation("Sequences are now loaded from the provided unitig graph.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(gfa);
    for s in &sequences {
        eprintln!("{}", s);
    }
    eprintln!();
    (unitig_graph, sequences)
}


fn load_from_fasta(filename: &Path) -> Vec<(FileSeqName, Vec<u8>)> {
    section_header("Loading sequences");
    explanation("Sequences are now loaded from the provided FASTA file.");
    let mut seqs = Vec::new();
    for (name, _, seq) in load_fasta(filename) {
        eprintln!("{} ({} bp)", name, seq.len());
        seqs.push((FileSeqName { filename: String::new(), seqname: name },
                   seq.as_bytes().to_owned()));
    }
    eprintln!();
    seqs
}


fn load_from_directory(dir: &Path) -> Vec<(FileSeqName, Vec<u8>)> {
    section_header("Loading sequences");
    explanation("Sequences are now loaded from FASTA files in the provided directory.");
    let mut seqs = Vec::new();
    let assemblies = find_all_assemblies(dir);
    for assembly in &assemblies {
        let filename = assembly.file_name().and_then(|name| name.to_str()).map(|s| s.to_string()).unwrap();
        for (name, _, seq) in load_fasta(assembly) {
            eprintln!("{} {} ({} bp)", filename, name, seq.len());
            seqs.push((FileSeqName { filename: filename.clone(), seqname: name },
                       seq.as_bytes().to_owned()));
        }
    }
    eprintln!();
    seqs
}


fn create_dotplot(seqs: &Vec<(FileSeqName, Vec<u8>)>, png_filename: &Path, res: u32, kmer: u32) {
    section_header("Creating dotplot");
    explanation("K-mers common between sequences are now used to build the dotplot image.");
    let pb = spinner("creating dotplot...");

    // We create an initial image to test the label sizes.
    let (top_left_gap, border_gap, between_seq_gap, text_gap, max_font_size) =
        get_sizes(res, seqs.len());
    let (start_positions, end_positions, _) =
        get_positions(seqs, res, kmer, top_left_gap, border_gap, between_seq_gap);
    let mut img = ImageBuffer::from_pixel(res, res, BACKGROUND_COLOUR);
    let font = FontArc::try_from_slice(include_bytes!("assets/DejaVuSans.ttf")).unwrap();
    let (text_height, _, _) = reduce_scale(seqs, &start_positions, &end_positions, &font,
                                           max_font_size);

    // Now that we know the values for text_height, we start over, this time readjusting the
    // top-left gap (so it isn't bigger than necessary).
    let top_left_gap = (2.0 * text_height) as u32 + border_gap;
    let (start_positions, end_positions, bp_per_pixel) =
        get_positions(seqs, res, kmer, top_left_gap, border_gap, between_seq_gap);
    draw_sequence_boxes(&mut img, seqs, &start_positions, &end_positions, true);
    draw_labels(&mut img, seqs, &start_positions, &end_positions, text_gap, &font, max_font_size);

    let mut count = 0;
    for (name_a, seq_a) in seqs {
        let rev_comp_seq_a = reverse_complement(seq_a);
        let kmers = get_all_kmer_positions(kmer as usize, seq_a, &rev_comp_seq_a);
        for (name_b, seq_b) in seqs {
            draw_dots(&mut img, start_positions[name_a], start_positions[name_b], seq_a, seq_b,
                      &kmers, bp_per_pixel);
            count += 1;
        }
    }

    // The boxes are drawn once more, this time with no fill, to overwrite any dots which leaked
    // into the outline.
    draw_sequence_boxes(&mut img, seqs, &start_positions, &end_positions, false);

    img.save(png_filename).unwrap();
    pb.finish_and_clear();
    eprintln!("{} pairwise dotplot{} drawn to image", count, match count { 1 => "", _ => "s" });
    eprintln!();
}


fn get_sizes(res: u32, seq_count: usize) -> (u32, u32, u32, u32, u32) {
    let res = res as f64;
    let top_left_gap = (INITIAL_TOP_LEFT_GAP * res).round() as u32;
    let border_gap = ((BORDER_GAP * res).round() as u32).max(2);
    let between_seq_gap = ((between_seq_gap(BETWEEN_SEQ_GAP, TOTAL_BETWEEN_SEQ_GAP,
                                            seq_count) * res).round() as u32).max(2);
    let text_gap = ((TEXT_GAP * res).round() as u32).max(1);
    let max_font_size = ((MAX_FONT_SIZE * res).round() as u32).max(1);
    (top_left_gap, border_gap, between_seq_gap, text_gap, max_font_size)
}


fn between_seq_gap(gap: f64, max_total_gap: f64, seq_count: usize) -> f64 {
    // The amount of space to use for gaps between sequences is set by BETWEEN_SEQ_GAP, unless
    // there are so many sequences that these gaps would exceed TOTAL_BETWEEN_SEQ_GAP, in which
    // case the gap size is scaled down.
    if seq_count <= 1 { return gap; }
    if (seq_count - 1) as f64 * gap > max_total_gap{
        max_total_gap / (seq_count - 1) as f64
    } else {
        gap
    }
}


fn get_positions(seqs: &Vec<(FileSeqName, Vec<u8>)>, res: u32, kmer: u32, top_left_gap: u32,
                 bottom_right_gap: u32, mut between_seq_gap: u32) ->
        (HashMap<FileSeqName, u32>, HashMap<FileSeqName, u32>, f64) {
    // This function returns the image coordinates that start/end each sequence. Since the dotplot
    // is symmetrical, there is only one start/end per sequence (used for both x and y coordinates).
    let mut seq_lengths: HashMap<FileSeqName, u32> = HashMap::new();
    for (key, seq) in seqs {
        seq_lengths.insert(key.clone(), (seq.len() as u32).saturating_sub(kmer).saturating_add(1));
    }
    let mut all_gaps = top_left_gap + bottom_right_gap + between_seq_gap * (seqs.len() as u32 - 1);
    let mut pixels_for_sequence = res.saturating_sub(all_gaps);

    // If there isn't enough room for the dotplots, reduce between_seq_gap.
    if all_gaps > pixels_for_sequence && seqs.len() > 1 {
        between_seq_gap = ((res / 2) - top_left_gap - bottom_right_gap) / (seqs.len() as u32 - 1);
        all_gaps = top_left_gap + bottom_right_gap + between_seq_gap * (seqs.len() as u32 - 1);
        pixels_for_sequence = res.saturating_sub(all_gaps);
    }

    let total_seq_length: u32 = seq_lengths.values().sum();
    let bp_per_pixel = total_seq_length as f64 / pixels_for_sequence as f64;

    let mut start_positions = HashMap::new();
    let mut end_positions = HashMap::new();
    let mut current_pos = top_left_gap;

    for (name, _) in seqs {
        start_positions.insert(name.clone(), current_pos);
        let rect_size = (seq_lengths[name] as f64 / bp_per_pixel).round() as u32;
        current_pos += rect_size;
        end_positions.insert(name.clone(), current_pos);
        current_pos += between_seq_gap;
    }
    (start_positions, end_positions, bp_per_pixel)
}


fn draw_sequence_boxes(img: &mut RgbImage, seqs: &Vec<(FileSeqName, Vec<u8>)>,
                       start_positions: &HashMap<FileSeqName, u32>,
                       end_positions: &HashMap<FileSeqName, u32>, fill: bool) {
    for (name_a, _) in seqs {
        let start_a = start_positions[name_a] - 1;
        let end_a = end_positions[name_a] + 2;
        for (name_b, _) in seqs {
            let start_b = start_positions[name_b] - 1;
            let end_b = end_positions[name_b] + 2;
            let rect = Rect::at(start_a as i32, start_b as i32)
                .of_size(end_a - start_a, end_b - start_b);
            if fill {
                let fill_colour = if name_a == name_b { SELF_VS_SELF_COLOUR }
                                                 else { SELF_VS_OTHER_COLOUR };
                draw_filled_rect_mut(img, rect, fill_colour);
            }
            draw_hollow_rect_mut(img, rect, OUTLINE_COLOUR);
        }
    }
}


fn reduce_scale(seqs: &Vec<(FileSeqName, Vec<u8>)>, start_positions: &HashMap<FileSeqName, u32>,
                end_positions: &HashMap<FileSeqName, u32>, font: &FontArc, max_font_size: u32)
        -> (f32, f32, PxScale) {
    // Reduces the scale as needed to ensure that all labels will fit in their available space.
    let mut text_height = max_font_size as f32;
    let mut scale = PxScale::from(text_height);
    let mut available_width: f32 = 1.0;
    for (name, _) in seqs {
        let start = start_positions[name];
        let end = end_positions[name];
        available_width = (end - start) as f32;
        let text_width = calculate_text_width(&name.filename, scale, font)
                             .max(calculate_text_width(&name.seqname, scale, font));
        if text_width > available_width {
            text_height *= available_width / text_width;
            scale = PxScale::from(text_height);
        }
    }
    (text_height, available_width, scale)
}


struct TextDimensions {
    width: u32,
    height: u32,
}


fn draw_labels(img: &mut RgbImage, seqs: &Vec<(FileSeqName, Vec<u8>)>,
               start_positions: &HashMap<FileSeqName, u32>,
               end_positions: &HashMap<FileSeqName, u32>,
               text_gap: u32, font: &FontArc, max_font_size: u32) {
    let min_pos = start_positions.values().min().cloned().unwrap_or_default();
    let (text_height, available_width, scale) =
        reduce_scale(seqs, start_positions, end_positions, font, max_font_size);
    let dim = TextDimensions { width: (available_width.ceil()) as u32, height: text_height as u32 };
    for (name, _) in seqs {
        let start = start_positions[name];
        let end = end_positions[name];
        let pos_1 = min_pos - text_gap - dim.height;
        let pos_2 = pos_1 - dim.height;

        // Horizontal labels on the top side.
        draw_text_mut(img, TEXT_COLOUR, start as i32, pos_1 as i32, scale, &font, &name.seqname);
        draw_text_mut(img, TEXT_COLOUR, start as i32, pos_2 as i32, scale, &font, &name.filename);

        // Vertical labels on the left side.
        draw_vertical_text(img, &dim, &name.seqname, pos_1, end, &scale, font);
        draw_vertical_text(img, &dim, &name.filename, pos_2, end, &scale, font);
    }
}


fn calculate_text_width(text: &str, scale: PxScale, font: &FontArc) -> f32 {
    let scaled_font = font.as_scaled(scale);
    text.chars().map(|c| {
        let glyph_id = scaled_font.glyph_id(c);
        scaled_font.h_advance(glyph_id)
    }).sum()
}


fn draw_vertical_text(img: &mut RgbImage, dim: &TextDimensions, text: &str, x: u32, y: u32,
                      scale: &PxScale, font: &FontArc) {
    // Draws text onto the image rotated 90 degrees counterclockwise. Does this by creating a temp
    // image with the text and then copying it over, pixel-by-pixel, with the appropriate
    // transformation.
    let (full_width, full_height) = (img.width(), img.height());
    let mut temp_img = ImageBuffer::from_pixel(dim.width, dim.height, BACKGROUND_COLOUR);
    draw_text_mut(&mut temp_img, TEXT_COLOUR, 0, 0, *scale, font, text);
    for i in 0..dim.width {
        let new_y = y - i;
        if new_y >= full_height { continue; }
        for j in 0..dim.height {
            let new_x = x + j;
            if new_x < full_width {
                let pixel = temp_img.get_pixel(i, j);
                if pixel != &BACKGROUND_COLOUR {
                    img.put_pixel(new_x, new_y, *pixel);
                }
            }
        }
    }
}


fn draw_dots(img: &mut RgbImage, a_start_pos: u32, b_start_pos: u32,
             seq_a: &[u8], seq_b: &[u8], a_kmers: &Kmers, bp_per_pixel: f64) {
    let (width, height) = (img.width(), img.height());
    if seq_a.len() < a_kmers.size || seq_b.len() < a_kmers.size {
        return;
    }
    for j in 0..(seq_b.len() - a_kmers.size + 1) {
        let j_pixel = (j as f64 / bp_per_pixel).round() as u32 + b_start_pos;
        let k = &seq_b[j..j + a_kmers.size];
        if let Some(a_reverse_positions) = a_kmers.reverse.get(k) {
            for &i in a_reverse_positions {
                let i_pixel = (i as f64 / bp_per_pixel).round() as u32 + a_start_pos;
                draw_dot(img, i_pixel, j_pixel, width, height, REVERSE_DOT_COLOUR);
            }
        }
        if let Some(a_forward_positions) = a_kmers.forward.get(k) {
            for &i in a_forward_positions {
                let i_pixel = (i as f64 / bp_per_pixel).round() as u32 + a_start_pos;
                draw_dot(img, i_pixel, j_pixel, width, height, FORWARD_DOT_COLOUR);
            }
        }
    }
}


fn draw_dot(img: &mut RgbImage, i: u32, j: u32, width: u32, height: u32, colour: Rgb<u8>) {
    if i < width && j < height {
        img.put_pixel(i, j, colour);
    }
}


struct Kmers<'a> {
    forward: HashMap<&'a [u8], Vec<u32>>,
    reverse: HashMap<&'a [u8], Vec<u32>>,
    size: usize
}


fn get_all_kmer_positions<'a>(kmer_size: usize, seq: &'a [u8], rev_comp_seq: &'a [u8])
        -> Kmers<'a> {
    let mut forward_kmers: HashMap<&[u8], Vec<u32>> = HashMap::new();
    let mut reverse_kmers: HashMap<&[u8], Vec<u32>> = HashMap::new();
    if seq.len() < kmer_size {
        return Kmers { forward: forward_kmers, reverse: reverse_kmers, size: kmer_size };
    }
    let seq_len = seq.len() - kmer_size + 1;
    for i in 0..seq_len {
        let forward_kmer = &seq[i..i+kmer_size];
        forward_kmers.entry(forward_kmer).or_default().push(i as u32);
        let reverse_kmer = &rev_comp_seq[i..i+kmer_size];
        reverse_kmers.entry(reverse_kmer).or_default().push((seq_len - i - 1) as u32);
    }
    assert!(forward_kmers.len() < seq.len());
    assert!(reverse_kmers.len() < seq.len());
    Kmers { forward: forward_kmers, reverse: reverse_kmers, size: kmer_size }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::assert_almost_eq;

    #[test]
    fn test_get_all_kmer_positions() {
        let seq = b"ACGACTGACATCAGCACTGA".to_vec();
        let rev_seq = reverse_complement(&seq);

        let expected_forward: HashMap<&[u8], Vec<u32>> = {
            let mut map = HashMap::new();
            map.insert(&b"ACGA"[..], vec![0]);     map.insert(&b"CGAC"[..], vec![1]);
            map.insert(&b"GACT"[..], vec![2]);     map.insert(&b"ACTG"[..], vec![3, 15]);
            map.insert(&b"CTGA"[..], vec![4, 16]); map.insert(&b"TGAC"[..], vec![5]);
            map.insert(&b"GACA"[..], vec![6]);     map.insert(&b"ACAT"[..], vec![7]);
            map.insert(&b"CATC"[..], vec![8]);     map.insert(&b"ATCA"[..], vec![9]);
            map.insert(&b"TCAG"[..], vec![10]);    map.insert(&b"CAGC"[..], vec![11]);
            map.insert(&b"AGCA"[..], vec![12]);    map.insert(&b"GCAC"[..], vec![13]);
            map.insert(&b"CACT"[..], vec![14]);
            map
        };

        let expected_reverse: HashMap<&[u8], Vec<u32>> = {
            let mut map = HashMap::new();
            map.insert(&b"TCAG"[..], vec![16, 4]); map.insert(&b"CAGT"[..], vec![15, 3]);
            map.insert(&b"AGTG"[..], vec![14]);    map.insert(&b"GTGC"[..], vec![13]);
            map.insert(&b"TGCT"[..], vec![12]);    map.insert(&b"GCTG"[..], vec![11]);
            map.insert(&b"CTGA"[..], vec![10]);    map.insert(&b"TGAT"[..], vec![9]);
            map.insert(&b"GATG"[..], vec![8]);     map.insert(&b"ATGT"[..], vec![7]);
            map.insert(&b"TGTC"[..], vec![6]);     map.insert(&b"GTCA"[..], vec![5]);
            map.insert(&b"AGTC"[..], vec![2]);     map.insert(&b"GTCG"[..], vec![1]);
            map.insert(&b"TCGT"[..], vec![0]);
            map
        };

        let kmers = get_all_kmer_positions(4, &seq, &rev_seq);
        assert_eq!(kmers.forward, expected_forward);
        assert_eq!(kmers.reverse, expected_reverse);
        assert_eq!(kmers.size, 4);
    }

    #[test]
    fn test_between_seq_gap() {
        assert_almost_eq(between_seq_gap(0.01, 0.1, 0), 0.01, 1e-8);
        assert_almost_eq(between_seq_gap(0.01, 0.1, 1), 0.01, 1e-8);
        assert_almost_eq(between_seq_gap(0.01, 0.1, 2), 0.01, 1e-8);
        assert_almost_eq(between_seq_gap(0.01, 0.1, 5), 0.01, 1e-8);
        assert_almost_eq(between_seq_gap(0.01, 0.1, 10), 0.01, 1e-8);
        assert_almost_eq(between_seq_gap(0.01, 0.1, 21), 0.005, 1e-8);
        assert_almost_eq(between_seq_gap(0.01, 0.1, 51), 0.002, 1e-8);
    }
}
