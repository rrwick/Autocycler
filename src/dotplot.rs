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

use ab_glyph::{FontArc, PxScale, ScaleFont};
use image::{RgbImage, Rgb};
use imageproc::drawing::{draw_filled_rect_mut, draw_text_mut, draw_hollow_rect_mut,
                         draw_line_segment_mut};
use imageproc::rect::Rect;
use std::collections::HashMap;
use std::path::PathBuf;

use crate::sequence::Sequence;
use crate::unitig_graph::UnitigGraph;


// Some hard-coded settings for the dotplot generation. Values between 0 and 1 are relative to full
// image resolution.
static INITIAL_TOP_LEFT_GAP: f64 = 0.1;
static BORDER_GAP: f64 = 0.015;
static BETWEEN_SEQ_GAP: f64 = 0.01;
static OUTLINE_WIDTH: f64 = 0.0015;
static TEXT_GAP: f64 = 0.005;
static MAX_FONT_SIZE: f64 = 0.025;
static BACKGROUND_COLOUR: image::Rgb<u8> = Rgb([255, 255, 255]);  // white
static SELF_VS_SELF_COLOUR: image::Rgb<u8> = Rgb([211, 211, 211]);  // lightgrey
static SELF_VS_OTHER_COLOUR: image::Rgb<u8> = Rgb([245, 245, 245]);  // whitesmoke
static TEXT_COLOUR: image::Rgb<u8> = Rgb([0, 0, 0]);  // black
static OUTLINE_COLOUR: image::Rgb<u8> = Rgb([0, 0, 0]);  // black
static FORWARD_STRAND_DOT_COLOUR: image::Rgb<u8> = Rgb([0, 0, 205]);  // mediumblue
static REVERSE_STRAND_DOT_COLOUR: image::Rgb<u8> = Rgb([178, 34, 34]);  // firebrick


pub fn dotplot(graph: &UnitigGraph, sequences: &Vec<Sequence>, png_filename: &PathBuf, res: u32,
               kmer: u32) {
    let seqs = graph.reconstruct_original_sequences_vec(sequences);
    let font = FontArc::try_from_slice(include_bytes!("assets/DejaVuSans.ttf")).expect("Error loading font");
    let mut img = RgbImage::new(res, res);
    draw_filled_rect_mut(&mut img, Rect::at(0, 0).of_size(res, res), BACKGROUND_COLOUR);
    let (top_left_gap, border_gap, between_seq_gap,
         text_gap, outline_width, max_font_size) = get_sizes(res);

    let (start_positions, end_positions, bp_per_pixel) =
        get_positions(&seqs, res, kmer, top_left_gap, border_gap, between_seq_gap);
    draw_sequence_boxes(&mut img, &seqs, &start_positions, &end_positions, outline_width, true);

    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO

    img.save(png_filename).unwrap();
}


fn get_sizes(res: u32) -> (u32, u32, u32, u32, u32, u32) {
    let res = res as f64;
    let top_left_gap = (INITIAL_TOP_LEFT_GAP * res).round() as u32;
    let border_gap = ((BORDER_GAP * res).round() as u32).max(2);
    let between_seq_gap = ((BETWEEN_SEQ_GAP * res).round() as u32).max(2);
    let text_gap = ((TEXT_GAP * res).round() as u32).max(1);
    let outline_width = ((OUTLINE_WIDTH * res).round() as u32).max(1);
    let max_font_size = ((MAX_FONT_SIZE * res).round() as u32).max(1);
    (top_left_gap, border_gap, between_seq_gap, text_gap, outline_width, max_font_size)
}


fn get_positions(seqs: &Vec<((String, String), String)>, res: u32, kmer: u32, top_left_gap: u32,
                 bottom_right_gap: u32, between_seq_gap: u32) -> (HashMap<(String, String), u32>, HashMap<(String, String), u32>, f64) {
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


fn draw_sequence_boxes(img: &mut RgbImage, seqs: &Vec<((String, String), String)>,
                       start_positions: &HashMap<(String, String), u32>,
                       end_positions: &HashMap<(String, String), u32>,
                       outline_width: u32, fill: bool) {
    for (name_a, _) in seqs {
        let start_a = start_positions[name_a] - outline_width;
        let end_a = end_positions[name_a] + outline_width;
        for (name_b, _) in seqs {
            let start_b = start_positions[name_b] - outline_width;
            let end_b = end_positions[name_b] + outline_width;
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
