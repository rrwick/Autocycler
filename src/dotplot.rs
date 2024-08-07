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
use std::path::PathBuf;

use crate::log::{section_header, explanation};
use crate::misc::{check_if_file_exists, quit_with_error, reverse_complement, spinner};
use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


// Some hard-coded settings. Values between 0 and 1 are relative to full image resolution.
static INITIAL_TOP_LEFT_GAP: f64 = 0.1;
static BORDER_GAP: f64 = 0.015;
static BETWEEN_SEQ_GAP: f64 = 0.01;
static TEXT_GAP: f64 = 0.0025;
static MAX_FONT_SIZE: f64 = 0.025;
static BACKGROUND_COLOUR: image::Rgb<u8> = Rgb([255, 255, 255]);     // white
static SELF_VS_SELF_COLOUR: image::Rgb<u8> = Rgb([211, 211, 211]);   // lightgrey
static SELF_VS_OTHER_COLOUR: image::Rgb<u8> = Rgb([245, 245, 245]);  // whitesmoke
static TEXT_COLOUR: image::Rgb<u8> = Rgb([0, 0, 0]);                 // black
static OUTLINE_COLOUR: image::Rgb<u8> = Rgb([0, 0, 0]);              // black
static FORWARD_DOT_COLOUR: image::Rgb<u8> = Rgb([0, 0, 205]);        // mediumblue
static REVERSE_DOT_COLOUR: image::Rgb<u8> = Rgb([178, 34, 34]);      // firebrick


pub fn dotplot(in_gfa: PathBuf, out_png: PathBuf, res: u32, kmer: u32) {
    check_settings(&in_gfa, res, kmer);
    starting_message();
    print_settings(&in_gfa, res, kmer);
    let (graph, sequences) = load_graph(&in_gfa);
    create_dotplot(&graph, &sequences, &out_png, res, kmer);
    finished_message(&out_png);
}


fn check_settings(in_gfa: &PathBuf, res: u32, kmer: u32) {
    check_if_file_exists(&in_gfa);
    if res < 100     { quit_with_error("--res cannot be less than 100"); }
    if res > 10000   { quit_with_error("--res cannot be greater than 10000"); }
    if kmer < 10     { quit_with_error("--kmer cannot be less than 10"); }
    if kmer > 100    { quit_with_error("--kmer cannot be greater than 100"); }
}


fn starting_message() {
    section_header("Starting autocycler dotplot");
    explanation("This command will take a unitig graph (either before or after trimming) and \
                 generate a dot plot image containing all pairwise comparisons of the sequences.");
}


fn print_settings(in_gfa: &PathBuf, res: u32, kmer: u32) {
    eprintln!("Settings:");
    eprintln!("  --in_gfa {}", in_gfa.display());
    eprintln!("  --res {}", res);
    eprintln!("  --kmer {}", kmer);
    eprintln!();
}


fn finished_message(out_png: &PathBuf) {
    section_header("Finished!");
    eprintln!("Pairwise dotplots: {}", out_png.display());
    eprintln!();
}


fn load_graph(gfa: &PathBuf) -> (UnitigGraph, Vec<Sequence>) {
    section_header("Loading graph");
    explanation("The unitig graph is now loaded into memory.");
    let (unitig_graph, sequences) = UnitigGraph::from_gfa_file(&gfa);
    unitig_graph.print_basic_graph_info();
    (unitig_graph, sequences)
}


fn create_dotplot(graph: &UnitigGraph, sequences: &Vec<Sequence>, png_filename: &PathBuf, res: u32,
               kmer: u32) {
    let seqs = graph.reconstruct_original_sequences_u8(sequences);
    eprintln!();
    let pb = spinner("creating dot plot...");

    // We create an initial image to test the label sizes.
    let (top_left_gap, border_gap, between_seq_gap, text_gap, max_font_size) = get_sizes(res);
    let (start_positions, end_positions, _) =
        get_positions(&seqs, res, kmer, top_left_gap, border_gap, between_seq_gap);
    let mut img = ImageBuffer::from_pixel(res, res, BACKGROUND_COLOUR);
    let font = FontArc::try_from_slice(include_bytes!("assets/DejaVuSans.ttf")).expect("Error loading font");
    let text_height = draw_labels(&mut img, &seqs, &start_positions, &end_positions, text_gap,
                                  &font, max_font_size, true);

    // Now that we know the values for text_height, we start over, this time readjusting the
    // top-left gap (so it isn't bigger than necessary).
    let top_left_gap = (2.0 * text_height) as u32 + border_gap;
    let (start_positions, end_positions, bp_per_pixel) =
        get_positions(&seqs, res, kmer, top_left_gap, border_gap, between_seq_gap);
    draw_sequence_boxes(&mut img, &seqs, &start_positions, &end_positions, true);
    for (name_a, seq_a) in &seqs {
        let rev_comp_seq_a = reverse_complement(seq_a);
        let (forward_kmers, reverse_kmers) = get_all_kmer_positions(kmer as usize, seq_a,
                                                                    &rev_comp_seq_a);
        for (name_b, seq_b) in &seqs {
            draw_dots(&mut img, name_a, name_b, seq_a, seq_b, &start_positions,
                      &forward_kmers, &reverse_kmers, bp_per_pixel, kmer as usize);
        }
    }
    draw_labels(&mut img, &seqs, &start_positions, &end_positions, text_gap, &font, max_font_size,
                false);

    // The boxes are drawn once more, this time with no fill, to overwrite any dots which leaked
    // into the outline, which would look messy.
    draw_sequence_boxes(&mut img, &seqs, &start_positions, &end_positions, false);

    img.save(png_filename).unwrap();
    pb.finish_and_clear();
}


fn get_sizes(res: u32) -> (u32, u32, u32, u32, u32) {
    let res = res as f64;
    let top_left_gap = (INITIAL_TOP_LEFT_GAP * res).round() as u32;
    let border_gap = ((BORDER_GAP * res).round() as u32).max(2);
    let between_seq_gap = ((BETWEEN_SEQ_GAP * res).round() as u32).max(2);
    let text_gap = ((TEXT_GAP * res).round() as u32).max(1);
    let max_font_size = ((MAX_FONT_SIZE * res).round() as u32).max(1);
    (top_left_gap, border_gap, between_seq_gap, text_gap, max_font_size)
}


fn get_positions(seqs: &Vec<((String, String), Vec<u8>)>, res: u32, kmer: u32, top_left_gap: u32,
                 bottom_right_gap: u32, between_seq_gap: u32) ->
        (HashMap<(String, String), u32>, HashMap<(String, String), u32>, f64) {
    // This function returns the image coordinates that start/end each sequence. Since the dot plot
    // is symmetrical, there is only one start/end per sequence (used for both x and y coordinates).
    let mut seq_lengths: HashMap<(String, String), u32> = HashMap::new();
    for (key, seq) in seqs {
        seq_lengths.insert(key.clone(), (seq.len() as u32).saturating_sub(kmer).saturating_add(1));
    }
    let all_gaps = top_left_gap + bottom_right_gap + between_seq_gap * (seqs.len() as u32 - 1);
    let pixels_for_sequence = res.saturating_sub(all_gaps);
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


fn draw_sequence_boxes(img: &mut RgbImage, seqs: &Vec<((String, String), Vec<u8>)>,
                       start_positions: &HashMap<(String, String), u32>,
                       end_positions: &HashMap<(String, String), u32>, fill: bool) {
    for (name_a, _) in seqs {
        let start_a = start_positions[name_a] - 1;
        let end_a = end_positions[name_a] + 2;
        for (name_b, _) in seqs {
            let start_b = start_positions[name_b] - 1;
            let end_b = end_positions[name_b] + 2;
            let rect = Rect::at(start_a as i32, start_b as i32)
                .of_size((end_a - start_a) as u32, (end_b - start_b) as u32);
            if fill {
                let fill_colour = if name_a == name_b { SELF_VS_SELF_COLOUR }
                                                 else { SELF_VS_OTHER_COLOUR };
                draw_filled_rect_mut(img, rect, fill_colour);
            }
            draw_hollow_rect_mut(img, rect, OUTLINE_COLOUR);
        }
    }
}


fn draw_labels(img: &mut RgbImage, seqs: &Vec<((String, String), Vec<u8>)>,
               start_positions: &HashMap<(String, String), u32>,
               end_positions: &HashMap<(String, String), u32>,
               text_gap: u32, font: &FontArc, max_font_size: u32, skip_draw: bool) -> f32 {
    let min_pos = start_positions.values().min().cloned().unwrap_or_default();
    let mut text_height_f32 = max_font_size as f32;
    let mut scale = PxScale::from(text_height_f32);
    let mut available_width = 1.0;

    // Reduce the scale as needed to ensure that all labels will fit in their available space.
    for (name, _) in seqs {
        let (filename, contig_name) = name;
        let start = start_positions[name];
        let end = end_positions[name];
        available_width = (end - start) as f32;
        let text_width = calculate_text_width(&filename, scale, &font).max(
                         calculate_text_width(&contig_name, scale, &font));
        if text_width > available_width {
            text_height_f32 *= available_width / text_width;
            scale = PxScale::from(text_height_f32);
        }
    }
    if skip_draw { return text_height_f32; }

    let text_height = text_height_f32 as u32;
    let text_width = (available_width.ceil()) as u32;
    for (name, _) in seqs {
        let (filename, contig_name) = name;
        let start = start_positions[name];
        let end = end_positions[name];

        // Horizontal labels on the top side.
        let pos_1 = min_pos - text_gap - text_height;
        let pos_2 = pos_1 - text_height;
        draw_text_mut(img, TEXT_COLOUR, start as i32, pos_1 as i32, scale, &font, contig_name);
        draw_text_mut(img, TEXT_COLOUR, start as i32, pos_2 as i32, scale, &font, filename);

        // Vertical labels on the left side.
        draw_vertical_text(img, text_width, text_height, contig_name, pos_1, end, &scale, &font);
        draw_vertical_text(img, text_width, text_height, filename, pos_2, end, &scale, &font);
    }
    text_height_f32
}


fn calculate_text_width(text: &str, scale: PxScale, font: &FontArc) -> f32 {
    let scaled_font = font.as_scaled(scale);
    text.chars()
        .filter_map(|c| {
            let glyph_id = scaled_font.glyph_id(c);
            Some(scaled_font.h_advance(glyph_id))
        }).sum()
}


fn draw_vertical_text(img: &mut RgbImage, width: u32, height: u32, text: &str, x: u32, y: u32,
                      scale: &PxScale, font: &FontArc) {
    // Draws text onto the image rotated 90 degrees counterclockwise. Does this by creating a temp
    // image with the text and then copying it over, pixel-by-pixel, with the appropriate
    // transformation.
    let (full_width, full_height) = (img.width(), img.height());
    let mut temp_img = ImageBuffer::from_pixel(width, height, BACKGROUND_COLOUR);
    draw_text_mut(&mut temp_img, TEXT_COLOUR, 0, 0, *scale, font, text);
    for i in 0..width {
        let new_y = y - i;
        if new_y >= full_height { continue; }
        for j in 0..height {
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


fn draw_dots(img: &mut RgbImage, name_a: &(String, String), name_b: &(String, String),
             seq_a: &[u8], seq_b: &[u8], start_positions: &HashMap<(String, String), u32>,
             a_forward_kmers: &HashMap<&[u8], Vec<u32>>,
             a_reverse_kmers: &HashMap<&[u8], Vec<u32>>, bp_per_pixel: f64, kmer_size: usize) {
    let a_start_pos = start_positions[name_a];
    let b_start_pos = start_positions[name_b];
    let (width, height) = (img.width(), img.height());
    if seq_a.len() < kmer_size || seq_b.len() < kmer_size {
        return;
    }
    for j in 0..(seq_b.len() - kmer_size + 1) {
        let j_pixel = (j as f64 / bp_per_pixel).round() as u32 + b_start_pos;
        let k = &seq_b[j..j + kmer_size];
        if let Some(a_reverse_positions) = a_reverse_kmers.get(k) {
            for &i in a_reverse_positions {
                let i_pixel = (i as f64 / bp_per_pixel).round() as u32 + a_start_pos;
                draw_dot(img, i_pixel, j_pixel, width, height, REVERSE_DOT_COLOUR);
            }
        }
        if let Some(a_forward_positions) = a_forward_kmers.get(k) {
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


fn get_all_kmer_positions<'a>(kmer_size: usize, seq: &'a [u8], rev_comp_seq: &'a [u8]) ->
        (HashMap<&'a [u8], Vec<u32>>, HashMap<&'a [u8], Vec<u32>>) {
    let mut forward_kmers: HashMap<&[u8], Vec<u32>> = HashMap::new();
    let mut reverse_kmers: HashMap<&[u8], Vec<u32>> = HashMap::new();
    if seq.len() < kmer_size {
        return (forward_kmers, reverse_kmers);
    }
    let seq_len = seq.len() - kmer_size + 1;
    for i in 0..seq_len {
        let forward_kmer = &seq[i..i+kmer_size];
        forward_kmers.entry(forward_kmer).or_insert_with(Vec::new).push(i as u32);
        let reverse_kmer = &rev_comp_seq[i..i+kmer_size];
        reverse_kmers.entry(reverse_kmer).or_insert_with(Vec::new).push((seq_len - i - 1) as u32);
    }
    assert!(forward_kmers.len() < seq.len());
    assert!(reverse_kmers.len() < seq.len());
    (forward_kmers, reverse_kmers)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_all_kmer_positions() {
        let seq = b"ACGACTGACATCAGCACTGA".to_vec();
        let rev_seq = reverse_complement(&seq);

        let expected_forward: HashMap<&[u8], Vec<u32>> = {
            let mut map = HashMap::new();
            map.insert(&b"ACGA"[..], vec![0]);
            map.insert(&b"CGAC"[..], vec![1]);
            map.insert(&b"GACT"[..], vec![2]);
            map.insert(&b"ACTG"[..], vec![3, 15]);
            map.insert(&b"CTGA"[..], vec![4, 16]);
            map.insert(&b"TGAC"[..], vec![5]);
            map.insert(&b"GACA"[..], vec![6]);
            map.insert(&b"ACAT"[..], vec![7]);
            map.insert(&b"CATC"[..], vec![8]);
            map.insert(&b"ATCA"[..], vec![9]);
            map.insert(&b"TCAG"[..], vec![10]);
            map.insert(&b"CAGC"[..], vec![11]);
            map.insert(&b"AGCA"[..], vec![12]);
            map.insert(&b"GCAC"[..], vec![13]);
            map.insert(&b"CACT"[..], vec![14]);
            map
        };

        let expected_reverse: HashMap<&[u8], Vec<u32>> = {
            let mut map = HashMap::new();
            map.insert(&b"TCAG"[..], vec![16, 4]);
            map.insert(&b"CAGT"[..], vec![15, 3]);
            map.insert(&b"AGTG"[..], vec![14]);
            map.insert(&b"GTGC"[..], vec![13]);
            map.insert(&b"TGCT"[..], vec![12]);
            map.insert(&b"GCTG"[..], vec![11]);
            map.insert(&b"CTGA"[..], vec![10]);
            map.insert(&b"TGAT"[..], vec![9]);
            map.insert(&b"GATG"[..], vec![8]);
            map.insert(&b"ATGT"[..], vec![7]);
            map.insert(&b"TGTC"[..], vec![6]);
            map.insert(&b"GTCA"[..], vec![5]);
            map.insert(&b"AGTC"[..], vec![2]);
            map.insert(&b"GTCG"[..], vec![1]);
            map.insert(&b"TCGT"[..], vec![0]);
            map
        };

        let (forward_kmers, reverse_kmers) = get_all_kmer_positions(4, &seq, &rev_seq);
        assert_eq!(forward_kmers, expected_forward);
        assert_eq!(reverse_kmers, expected_reverse);
    }
}
