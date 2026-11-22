//! LineUtils: poly-A/T tails of LINE insertions.

use std::collections::BTreeMap;

use crate::common::compat;
use crate::assembly::assembly_utils::reverse_complement_str;
use crate::assembly::common_utils::{above_min_qual, is_uncertain_base_qual};
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_utils::{get_read_index_at_reference_position, INVALID_INDEX};
use crate::assembly::sam::{CigarOp, SamRecord};
use crate::assembly::types::assembly_link::{AssemblyLink, LinkP};
use crate::assembly::types::junction::Junction;
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::{LinkType, SupportType};

pub const MAX_NON_LINE_BASES: i32 = LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ;

pub fn find_line_sequence_count(bases: &[u8], index_start: i32, index_end: i32, line_base: u8) -> i32 {
    if index_start < 0 || index_end >= bases.len() as i32 {
        return 0;
    }
    if index_end - index_start + 1 < LINE_POLY_AT_REQ {
        return 0;
    }
    let mut line_base_count = 0;
    let mut other_count = 0;
    for i in index_start..=index_end {
        if bases[i as usize] == line_base {
            line_base_count += 1;
        } else {
            other_count += 1;
            if other_count > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ {
                break;
            }
        }
    }
    if line_base_count >= LINE_POLY_AT_REQ { line_base_count } else { 0 }
}

pub fn has_line_tail_record(record: &SamRecord) -> bool {
    let cigar = &record.cigar;
    if cigar.len() <= 1 {
        return false;
    }
    for from_start in [true, false] {
        let element = if from_start { cigar[0] } else { cigar[cigar.len() - 1] };
        if element.op != CigarOp::S {
            continue;
        }
        let sc_base_count = element.len;
        if sc_base_count == 0 {
            continue;
        }
        let line_base = if from_start { LINE_BASE_A } else { LINE_BASE_T };
        let soft_clip_index = if from_start { sc_base_count - 1 } else { record.read_bases.len() as i32 - sc_base_count };
        if has_line_tail(&record.read_bases, soft_clip_index, from_start, line_base, Some(&record.base_qualities)) {
            return true;
        }
    }
    false
}

pub fn has_line_tail(bases: &[u8], soft_clip_index: i32, left_clipped: bool, line_base: u8, base_quals: Option<&[u8]>) -> bool {
    let len = bases.len() as i32;
    let sc_length = if left_clipped { soft_clip_index + 1 } else { len - soft_clip_index };
    if sc_length < LINE_POLY_AT_REQ {
        return false;
    }
    let mut index = soft_clip_index;
    let mut total_base_count = 0;
    let mut non_line_base_count = 0;
    let mut line_base_count = 0;
    while index >= 0 && index < len {
        total_base_count += 1;
        if bases[index as usize] == line_base {
            line_base_count += 1;
        } else if base_quals.is_none() || !is_uncertain_base_qual(base_quals.unwrap()[index as usize]) {
            if total_base_count > LINE_POLY_AT_TEST_LEN {
                break;
            }
            non_line_base_count += 1;
            if non_line_base_count > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ {
                return false;
            }
        }
        index += if left_clipped { -1 } else { 1 };
    }
    if line_base_count < LINE_POLY_AT_REQ || non_line_base_count > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ {
        return false;
    }
    let ref_index_start = if left_clipped { soft_clip_index + 1 } else { soft_clip_index - 1 };
    let ref_line_base_count = find_ref_base_repeat_count(bases, ref_index_start, left_clipped, line_base);
    if ref_line_base_count > LINE_REF_BASE_MAX_LENGTH {
        return false;
    }
    line_base_count as f64 >= ref_line_base_count as f64 * LINE_REF_BASE_REPEAT_FACTOR
}

pub fn find_ref_base_repeat_count(bases: &[u8], mut index: i32, move_forward: bool, test_base: u8) -> i32 {
    let mut match_count = 0;
    let mut initial_match_count = 0;
    let mut mismatch_count = 0;
    let mut base_count = 0;
    while index >= 0 && index < bases.len() as i32 {
        base_count += 1;
        if bases[index as usize] == test_base {
            match_count += 1;
            if mismatch_count == 0 {
                initial_match_count += 1;
            }
        } else {
            mismatch_count += 1;
        }
        if base_count >= 10 {
            if mismatch_count as f64 / base_count as f64 > 0.1 {
                return initial_match_count;
            }
        } else if mismatch_count >= 2 {
            return initial_match_count;
        }
        index += if move_forward { 1 } else { -1 };
    }
    if mismatch_count == 0 { initial_match_count } else { match_count }
}

pub fn find_consensus_line_extension(reads: &[ReadP], junction: &Junction) -> i32 {
    let mut length_frequency = find_line_extension_frequency(reads, junction, true, true);
    if !length_frequency.is_empty() {
        return calc_median(&length_frequency);
    }
    length_frequency = find_line_extension_frequency(reads, junction, false, true);
    if !length_frequency.is_empty() {
        return calc_median(&length_frequency);
    }
    length_frequency = find_line_extension_frequency(reads, junction, false, false);
    length_frequency.keys().copied().max().unwrap_or(0)
}

pub fn find_line_extension_frequency(reads: &[ReadP], junction: &Junction, five_prime_only: bool, terminating_only: bool) -> BTreeMap<i32, i32> {
    let line_base = if junction.is_forward() { LINE_BASE_T } else { LINE_BASE_A };
    let is_forward = junction.is_forward();
    let mut length_frequency: BTreeMap<i32, i32> = BTreeMap::new();
    for read in reads {
        if !read.has_line_tail(is_forward) {
            continue;
        }
        if five_prime_only && read.positive_strand() == is_forward {
            continue;
        }
        let read_junction_index = get_read_index_at_reference_position(read, junction.position, true);
        if read_junction_index == INVALID_INDEX {
            continue;
        }
        let first_extension_index = read_junction_index + if is_forward { 1 } else { -1 };
        let line_extension_index = find_line_extension_end_index(read, line_base, first_extension_index, is_forward);
        let line_length = (first_extension_index - line_extension_index).abs() + 1;
        if line_length < LINE_POLY_AT_REQ {
            continue;
        }
        if terminating_only {
            let soft_clip_length = if is_forward { read.right_clip_length() } else { read.left_clip_length() };
            if soft_clip_length == line_length {
                continue;
            }
        }
        *length_frequency.entry(line_length).or_insert(0) += 1;
    }
    length_frequency
}

pub fn find_line_extension_end_index(read: &Read, line_base: u8, extension_index: i32, is_forward: bool) -> i32 {
    let mut index = extension_index;
    let mut line_length = 0;
    let mut non_line_count = 0;
    while index >= 0 && index < read.bases_length() {
        if read.get_bases()[index as usize] == line_base {
            line_length += 1;
        } else {
            if line_length >= LINE_POLY_AT_REQ {
                break;
            }
            if above_min_qual(read.get_base_quality()[index as usize]) {
                non_line_count += 1;
                if non_line_count > MAX_NON_LINE_BASES {
                    break;
                }
            }
        }
        index += if is_forward { 1 } else { -1 };
    }
    index += if is_forward { -1 } else { 1 };
    index
}

fn calc_median(length_frequency: &BTreeMap<i32, i32>) -> i32 {
    if length_frequency.is_empty() {
        return 0;
    }
    let mut lengths: Vec<i32> = Vec::new();
    for (k, v) in length_frequency {
        for _ in 0..*v {
            lengths.push(*k);
        }
    }
    lengths.sort();
    lengths[lengths.len() / 2]
}

pub fn is_line_with_local_aligned_insert(assembly: &JunctionAssembly) -> bool {
    if assembly.junction().indel_based() || !assembly.has_line_sequence() {
        return false;
    }
    let poly_at_length = calc_line_sequence_length(&assembly.form_junction_sequence(), assembly.is_forward_junction());
    let mut max_post_poly_at_length = 0;
    let java_a14 = compat::java("A14");
    for read in assembly.support() {
        if read.type_() != SupportType::Indel {
            continue;
        }
        let Some(indel_coords) = read.cached_read().expect("support read without cached read").indel_coords() else { continue };
        if indel_coords.is_delete() {
            continue;
        }
        if (indel_coords.length - poly_at_length).abs() <= 2 {
            let post_insert_aligned_length =
                if assembly.is_forward_junction() { read.alignment_end() - indel_coords.pos_end } else { indel_coords.pos_start - read.alignment_start() };
            if java_a14 {
                // Java: assigned for each read (then floored at 0), so the last read decides
                max_post_poly_at_length = post_insert_aligned_length.max(0);
            } else {
                max_post_poly_at_length = max_post_poly_at_length.max(post_insert_aligned_length);
            }
        }
    }
    max_post_poly_at_length >= MIN_VARIANT_LENGTH
}

pub fn has_line_source_sequence(assembly: &JunctionAssembly) -> bool {
    let ext_index_start;
    let ext_index_end;
    if assembly.is_forward_junction() {
        ext_index_start = assembly.junction_index() + 1;
        ext_index_end = ext_index_start + LINE_POLY_AT_TEST_LEN - 1;
    } else {
        ext_index_end = assembly.junction_index() - 1;
        ext_index_start = ext_index_end - LINE_POLY_AT_TEST_LEN + 1;
    }
    let line_base = if assembly.is_forward_junction() { LINE_BASE_A } else { LINE_BASE_T };
    find_line_sequence_count(assembly.bases(), ext_index_start, ext_index_end, line_base) >= LINE_POLY_AT_REQ
}

pub fn try_line_sequence_link(first: AsmP, second: AsmP, first_reversed: bool, second_reversed: bool) -> Option<LinkP> {
    if !first.has_line_sequence() && !second.has_line_sequence() {
        return None;
    }
    if first.has_line_sequence() && second.has_line_sequence() {
        return None;
    }
    let inserted_bases;
    if first.has_line_sequence() {
        let line_extension_bases = first.form_junction_sequence();
        let poly_at_length = calc_line_sequence_length(&line_extension_bases, true) as usize;
        let line_bases = &line_extension_bases[..poly_at_length];
        let mut other_bases = find_line_insertion_sequence(&second)?;
        if second_reversed {
            other_bases = reverse_complement_str(&other_bases);
        }
        inserted_bases = format!("{line_bases}{other_bases}");
    } else {
        let line_extension_bases = second.form_junction_sequence();
        let poly_at_length = calc_line_sequence_length(&line_extension_bases, false) as usize;
        let line_bases = &line_extension_bases[line_extension_bases.len() - poly_at_length..];
        let mut other_bases = find_line_insertion_sequence(&first)?;
        if first_reversed {
            other_bases = reverse_complement_str(&other_bases);
        }
        inserted_bases = format!("{other_bases}{line_bases}");
    }
    let link = AssemblyLink::new(first, second, LinkType::Split, &inserted_bases, "");
    link.m().mark_insert_site();
    Some(link)
}

fn find_line_insertion_sequence(assembly: &JunctionAssembly) -> Option<String> {
    let extension_bases = assembly.form_junction_sequence();
    let eb = extension_bases.as_bytes();
    let ext_base_length = eb.len() as i32;
    let is_forward = assembly.is_forward_junction();
    let mut i: i32 = if is_forward { ext_base_length - 1 } else { 0 };
    while i >= 0 && i < ext_base_length {
        let base = eb[i as usize];
        if base == LINE_BASE_A || base == LINE_BASE_T {
            let mut line_bases = 1;
            let mut non_line_bases = 0;
            let mut last_line_base_index = 0i32;
            let mut j = if is_forward { i - 1 } else { i + 1 };
            while j >= 0 && j < ext_base_length {
                if eb[j as usize] == base {
                    line_bases += 1;
                    last_line_base_index = j;
                } else {
                    non_line_bases += 1;
                    if non_line_bases > MAX_NON_LINE_BASES || line_bases >= LINE_POLY_AT_REQ {
                        break;
                    }
                }
                j += if is_forward { -1 } else { 1 };
            }
            if line_bases >= LINE_POLY_AT_REQ {
                return Some(if is_forward { extension_bases[..last_line_base_index as usize].to_string() } else { extension_bases[(last_line_base_index + 1) as usize..].to_string() });
            }
        }
        i += if is_forward { -1 } else { 1 };
    }
    None
}

pub fn calc_line_sequence_length(sequence: &str, from_start: bool) -> i32 {
    let s = sequence.as_bytes();
    let sequence_length = s.len() as i32;
    // an empty sequence is a StringIndexOutOfBoundsException in Java
    let line_base = if from_start { s[0] } else { s[(sequence_length - 1) as usize] };
    if line_base != LINE_BASE_A && line_base != LINE_BASE_T {
        return 0;
    }
    let mut count = 1;
    let mut index = if from_start { 1 } else { sequence_length - 2 };
    while index >= 0 && index < sequence_length {
        if s[index as usize] != line_base {
            break;
        }
        count += 1;
        index += if from_start { 1 } else { -1 };
    }
    count
}

/// hmf-common LineElements.isMobileLineElement
pub fn is_mobile_line_element(orientation: Orientation, insert_sequence: &str) -> bool {
    let bases = insert_sequence.as_bytes();
    let ins_seq_length = bases.len() as i32;
    if ins_seq_length < LINE_POLY_AT_REQ {
        return false;
    }
    let poly_at_char = if orientation.is_forward() { b'T' } else { b'A' };
    let test_length = LINE_POLY_AT_TEST_LEN.min(ins_seq_length);
    let mut allowed_non_required_chars = test_length - LINE_POLY_AT_REQ;
    for i in 0..test_length {
        let base = if orientation.is_forward() { bases[i as usize] } else { bases[(ins_seq_length - i - 1) as usize] };
        if base != poly_at_char {
            allowed_non_required_chars -= 1;
        }
        if allowed_non_required_chars < 0 {
            return false;
        }
    }
    true
}
