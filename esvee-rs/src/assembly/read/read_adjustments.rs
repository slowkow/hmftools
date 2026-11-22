//! ReadAdjustments: poly-G and low base quality trimming, LINE tail marking.

use crate::assembly::common_utils::{below_min_qual, is_uncertain_base_qual};
use crate::assembly::constants::*;
use crate::assembly::line_utils::has_line_tail;
use crate::assembly::read::read::Read;

pub const LOW_QUAL_SCORE: f64 = 1.0 / LOW_BASE_TRIM_PERC - 1.0;

pub fn trim_poly_g_sequences(read: &mut Read) -> bool {
    let mut trailing_g_count = 0;
    if read.positive_strand() {
        for i in (0..read.bases_length() as usize).rev() {
            if read.get_bases()[i] == b'G' {
                trailing_g_count += 1;
            } else {
                break;
            }
        }
    } else {
        for i in 0..read.bases_length() as usize {
            if read.get_bases()[i] == b'C' {
                trailing_g_count += 1;
            } else {
                break;
            }
        }
    }
    if trailing_g_count < POLY_G_TRIM_LENGTH {
        return false;
    }
    let from_start = read.negative_strand();
    read.trim_bases(trailing_g_count, from_start);
    true
}

pub fn mark_line_soft_clips(read: &mut Read) {
    for from_start in [true, false] {
        let sc_base_count = if from_start { read.left_clip_length() } else { read.right_clip_length() };
        if sc_base_count == 0 {
            continue;
        }
        let line_base = if from_start { LINE_BASE_A } else { LINE_BASE_T };
        let soft_clip_index = if from_start { sc_base_count - 1 } else { read.bases_length() - sc_base_count };
        if has_line_tail(read.get_bases(), soft_clip_index, from_start, line_base, Some(read.get_base_quality())) {
            read.mark_line_tail(from_start);
            return;
        }
    }
}

pub fn trim_low_qual_soft_clip_bases(read: &mut Read) -> bool {
    let from_start = read.negative_strand();
    let sc_base_count = if from_start { read.left_clip_length() } else { read.right_clip_length() };
    if sc_base_count == 0 {
        return false;
    }
    let mut line_exclusion_length = 0;
    if sc_base_count >= LINE_POLY_AT_REQ {
        let sc_index_start;
        let sc_index_end;
        let line_test_length = sc_base_count.min(LINE_POLY_AT_TEST_LEN);
        if from_start {
            sc_index_end = sc_base_count - 1;
            sc_index_start = sc_index_end - line_test_length + 1;
        } else {
            sc_index_start = read.bases_length() - sc_base_count;
            sc_index_end = sc_index_start + line_test_length - 1;
        }
        let line_base = if from_start { LINE_BASE_A } else { LINE_BASE_T };
        let has_tail = if from_start { read.has_line_tail_start() } else { read.has_line_tail_end() };
        if has_tail {
            line_exclusion_length = line_test_length;
            let mut base_index = if from_start { sc_index_start } else { sc_index_end };
            let mut base_check = line_test_length - LINE_POLY_AT_REQ;
            while base_check > 0 {
                if read.get_bases()[base_index as usize] != line_base {
                    line_exclusion_length -= 1;
                }
                if from_start {
                    base_index += 1;
                } else {
                    base_index -= 1;
                }
                base_check -= 1;
            }
        }
    }
    let read_index_start;
    let read_index_end;
    let sc_check_length = sc_base_count - line_exclusion_length;
    if from_start {
        read_index_start = 0;
        read_index_end = read_index_start + sc_check_length - 1;
    } else {
        read_index_end = read.bases_length() - 1;
        read_index_start = read_index_end - sc_check_length + 1;
    }
    let trim_count = find_low_base_qual_trim_count(read, read_index_start, read_index_end);
    if trim_count <= 0 {
        return false;
    }
    read.trim_bases(trim_count, from_start);
    true
}

pub fn find_low_base_qual_trim_count(read: &Read, index_start: i32, index_end: i32) -> i32 {
    find_low_base_qual_trim_count_ext(read, index_start, index_end, read.negative_strand(), false)
}

pub fn find_low_base_qual_trim_count_ext(read: &Read, index_start: i32, index_end: i32, from_start: bool, uncertain_only: bool) -> i32 {
    let mut lowest_score = 0.0f64;
    let mut current_score = 0.0f64;
    let mut last_lowest_score_index = -1;
    let mut base_index = if from_start { index_start } else { index_end };
    let quals = read.get_base_quality();
    while base_index >= index_start && base_index <= index_end {
        let qual = quals[base_index as usize];
        let is_low_qual = if uncertain_only { is_uncertain_base_qual(qual) } else { below_min_qual(qual) };
        if is_low_qual {
            current_score -= LOW_QUAL_SCORE;
            if current_score <= lowest_score {
                last_lowest_score_index = base_index;
                lowest_score = current_score;
            }
        } else {
            current_score += 1.0;
        }
        if from_start {
            base_index += 1;
        } else {
            base_index -= 1;
        }
    }
    if last_lowest_score_index < 0 {
        return 0;
    }
    if from_start { last_lowest_score_index - index_start + 1 } else { index_end - last_lowest_score_index + 1 }
}
