//! RepeatInfo: short tandem repeats found in a base sequence.

use crate::assembly::common_utils::{above_min_qual, positions_overlap};
use crate::assembly::constants::REPEAT_MIN_COUNT;
use crate::assembly::orientation::Orientation;
use crate::assembly::types::junction_assembly::JunctionAssembly;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatInfo {
    pub index: i32,
    pub bases: String,
    pub count: i32,
    direction: Option<Orientation>,
}

const DUAL_LENGTH: i32 = 2;
const THREE_LENGTH: i32 = 3;
const FOUR_LENGTH: i32 = 4;
const FIVE_LENGTH: i32 = 5;
const DUAL_DUAL_LENGTH: i32 = 4;

impl RepeatInfo {
    pub fn new(index: i32, bases: &str, count: i32) -> Self {
        RepeatInfo { index, bases: bases.to_string(), count, direction: None }
    }
    pub fn with_direction(index: i32, bases: &str, count: i32, orientation: Option<Orientation>) -> Self {
        RepeatInfo { index, bases: bases.to_string(), count, direction: orientation }
    }
    pub fn last_index(&self) -> i32 {
        if self.direction.is_none_or(|d| d.is_forward()) { self.index_end() } else { self.index_start() }
    }
    pub fn post_repeat_index(&self) -> i32 {
        if self.direction.is_none_or(|d| d.is_forward()) { self.index + self.total_length() } else { self.index - self.total_length() }
    }
    pub fn index_start(&self) -> i32 {
        if self.direction().is_forward() { self.index } else { self.index - self.total_length() + 1 }
    }
    pub fn index_end(&self) -> i32 {
        if self.direction().is_forward() { self.index + self.total_length() - 1 } else { self.index }
    }
    pub fn total_length(&self) -> i32 {
        self.count * self.bases.len() as i32
    }
    pub fn repeat_length(&self) -> i32 {
        self.bases.len() as i32
    }
    pub fn direction(&self) -> Orientation {
        self.direction.unwrap_or(Orientation::Forward)
    }
    pub fn matches_type(&self, other: &RepeatInfo) -> bool {
        self.bases == other.bases
    }
    pub fn overlaps(&self, other: &RepeatInfo) -> bool {
        positions_overlap(self.index_start(), self.index_end(), other.index_start(), other.index_end())
    }
}

impl std::fmt::Display for RepeatInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}-{}: {}-{}", self.index_start(), self.index_end(), self.bases, self.count)
    }
}

pub fn find_repeats(bases: &[u8]) -> Vec<RepeatInfo> {
    find_repeats_min(bases, REPEAT_MIN_COUNT)
}

pub fn find_repeats_min(bases: &[u8], min_repeat_count: i32) -> Vec<RepeatInfo> {
    let mut repeats = Vec::new();
    let mut index = 0i32;
    let len = bases.len() as i32;
    while index <= len - REPEAT_MIN_COUNT {
        let repeat = find_single_base_repeat(bases, index, min_repeat_count)
            .or_else(|| find_dual_base_repeat(bases, index, min_repeat_count))
            .or_else(|| find_multi_base_repeat(bases, index, THREE_LENGTH, min_repeat_count))
            .or_else(|| find_multi_base_repeat(bases, index, FOUR_LENGTH, min_repeat_count))
            .or_else(|| find_multi_base_repeat(bases, index, FIVE_LENGTH, min_repeat_count))
            .or_else(|| find_dual_dual_repeat(bases, index, min_repeat_count));
        match repeat {
            Some(r) => {
                index += r.total_length();
                repeats.push(r);
            }
            None => index += 1,
        }
    }
    repeats
}

pub fn find_repeats_offset(bases: &[u8], start_offset: i32, end_offset: i32) -> Vec<RepeatInfo> {
    if start_offset == 0 && end_offset == 0 {
        return find_repeats(bases);
    }
    let extension_bases = &bases[start_offset as usize..bases.len() - end_offset as usize];
    let repeats = find_repeats(extension_bases);
    if start_offset == 0 {
        return repeats;
    }
    repeats.into_iter().map(|x| RepeatInfo::new(x.index + start_offset, &x.bases, x.count)).collect()
}

fn base_str(bases: &[u8], index: i32, len: i32) -> String {
    bases[index as usize..(index + len) as usize].iter().map(|b| *b as char).collect()
}

pub fn find_single_base_repeat(bases: &[u8], index: i32, min_repeat_count: i32) -> Option<RepeatInfo> {
    let len = bases.len() as i32;
    if index + min_repeat_count - 1 >= len {
        return None;
    }
    let mut i = index + 1;
    while i < len {
        if bases[i as usize] != bases[index as usize] {
            break;
        }
        i += 1;
    }
    let repeat_count = i - index;
    if repeat_count >= min_repeat_count { Some(RepeatInfo::new(index, &base_str(bases, index, 1), repeat_count)) } else { None }
}

pub fn find_dual_base_repeat(bases: &[u8], index: i32, min_repeat_count: i32) -> Option<RepeatInfo> {
    let len = bases.len() as i32;
    if index + min_repeat_count * 2 - 1 >= len {
        return None;
    }
    let mut i = index + DUAL_LENGTH;
    while i < len - 1 {
        if bases[i as usize] != bases[index as usize] || bases[i as usize + 1] != bases[index as usize + 1] {
            break;
        }
        i += DUAL_LENGTH;
    }
    let repeat_count = (i - index) / DUAL_LENGTH;
    if repeat_count < min_repeat_count {
        return None;
    }
    Some(RepeatInfo::new(index, &base_str(bases, index, 2), repeat_count))
}

pub fn find_multi_base_repeat(bases: &[u8], index: i32, repeat_length: i32, min_repeat_count: i32) -> Option<RepeatInfo> {
    let len = bases.len() as i32;
    if index + min_repeat_count * repeat_length - 1 >= len {
        return None;
    }
    let mut i = index + repeat_length;
    while i < len - (repeat_length - 1) {
        let mut matched_bases = 0;
        for j in 0..repeat_length {
            if bases[(i + j) as usize] == bases[(index + j) as usize] {
                matched_bases += 1;
            } else {
                break;
            }
        }
        if matched_bases != repeat_length {
            break;
        }
        i += repeat_length;
    }
    let repeat_count = (i - index) / repeat_length;
    if repeat_count < min_repeat_count {
        return None;
    }
    Some(RepeatInfo::new(index, &base_str(bases, index, repeat_length), repeat_count))
}

pub fn find_dual_dual_repeat(bases: &[u8], index: i32, min_repeat_count: i32) -> Option<RepeatInfo> {
    let len = bases.len() as i32;
    if index + min_repeat_count * DUAL_DUAL_LENGTH - 1 >= len {
        return None;
    }
    let ix = index as usize;
    if bases[ix] != bases[ix + 1] || bases[ix + 2] != bases[ix + 3] {
        return None;
    }
    let mut i = index + DUAL_DUAL_LENGTH;
    while i < len - 3 {
        let iu = i as usize;
        if bases[iu] != bases[ix] || bases[iu + 1] != bases[ix + 1] || bases[iu + 2] != bases[ix + 2] || bases[iu + 3] != bases[ix + 3] {
            break;
        }
        i += 4;
    }
    let repeat_length = (i - index) / DUAL_DUAL_LENGTH;
    if repeat_length < min_repeat_count {
        return None;
    }
    Some(RepeatInfo::new(index, &base_str(bases, index, 4), repeat_length))
}

pub fn get_repeat_count(bases: &[u8], repeat: &str, read_index_start: i32, search_forward: bool) -> i32 {
    get_repeat_count_ext(bases, repeat, read_index_start, search_forward, 0, None)
}

pub fn get_repeat_count_ext(bases: &[u8], repeat: &str, read_index_start: i32, search_forward: bool, previous_repeat_count: i32, base_quals: Option<&[u8]>) -> i32 {
    let mut repeat_count = 0;
    let repeat_length = repeat.len() as i32;
    let len = bases.len() as i32;
    let mut read_index = read_index_start;
    if !search_forward {
        read_index -= repeat_length - 1;
    }
    if read_index < 0 || read_index >= len - repeat_length + 1 {
        return -1;
    }
    let repeat_bases = repeat.as_bytes();
    let mut low_qual_mismatch_count = 0;
    let mut total_repeat_span = previous_repeat_count * repeat_length;
    let mut mismatch_repeat_count = 0;
    while read_index >= 0 && read_index < len - repeat_length + 1 {
        let mut has_mismatch = false;
        for j in 0..repeat_length {
            if bases[(read_index + j) as usize] != repeat_bases[j as usize] {
                if base_quals.is_none() || above_min_qual(base_quals.unwrap()[(read_index + j) as usize]) {
                    has_mismatch = true;
                    break;
                }
                low_qual_mismatch_count += 1;
            }
        }
        if has_mismatch {
            break;
        }
        if low_qual_mismatch_count > 0 {
            mismatch_repeat_count += 1;
        } else {
            repeat_count += 1;
        }
        read_index += if search_forward { repeat_length } else { -repeat_length };
        total_repeat_span += repeat_length;
    }
    if low_qual_mismatch_count > 0 && low_qual_mismatch_count as f64 <= 0.1 * total_repeat_span as f64 {
        repeat_count += mismatch_repeat_count;
    }
    repeat_count
}

pub fn build_trimmed_ref_base_sequence(assembly: &JunctionAssembly, max_sequence_length: i32) -> String {
    let repeats = assembly.repeat_info();
    if repeats.is_empty() {
        return String::new();
    }
    let ref_base_length = assembly.ref_base_length();
    let seq_start = if assembly.is_forward_junction() { assembly.junction_index() - ref_base_length + 1 } else { assembly.junction_index() };
    let seq_end = if assembly.is_forward_junction() { assembly.junction_index() } else { assembly.junction_index() + ref_base_length - 1 };
    let mut out = String::new();
    let mut current_index = seq_start;
    let mut current_repeat_index = 0usize;
    let mut trimmed_bases_length = 0;
    let mut current_repeat: Option<&RepeatInfo> = Some(&repeats[0]);
    while current_index <= seq_end {
        while let Some(r) = current_repeat {
            if r.index >= current_index {
                break;
            }
            current_repeat_index += 1;
            current_repeat = repeats.get(current_repeat_index);
        }
        match current_repeat {
            Some(r) if r.index == current_index => {
                if trimmed_bases_length != 0 {
                    out.push('_');
                }
                out.push_str(&format!("{}{}_", r.bases, r.count));
                current_index += r.total_length();
                trimmed_bases_length += r.bases.len() as i32 * 2;
            }
            _ => {
                out.push(assembly.bases()[current_index as usize] as char);
                current_index += 1;
                trimmed_bases_length += 1;
            }
        }
        if trimmed_bases_length >= max_sequence_length {
            break;
        }
    }
    out
}

pub fn calc_trimmed_base_length(seq_start: i32, seq_end: i32, repeats: &[RepeatInfo]) -> i32 {
    if repeats.is_empty() {
        return seq_end - seq_start + 1;
    }
    let mut current_index = seq_start;
    let mut current_repeat_index = 0usize;
    let mut trimmed_bases_length = 0;
    let mut current_repeat: Option<&RepeatInfo> = Some(&repeats[0]);
    while current_index <= seq_end {
        while let Some(r) = current_repeat {
            if r.index >= current_index {
                break;
            }
            current_repeat_index += 1;
            current_repeat = repeats.get(current_repeat_index);
        }
        match current_repeat {
            Some(r) if r.index == current_index => {
                current_index += r.total_length();
                trimmed_bases_length += r.bases.len() as i32 * 2;
            }
            _ => {
                current_index += 1;
                trimmed_bases_length += 1;
            }
        }
    }
    trimmed_bases_length
}

pub fn repeats_as_string(repeats: &[RepeatInfo]) -> String {
    if repeats.is_empty() {
        return String::new();
    }
    let items: Vec<String> = repeats.iter().map(|x| format!("{}x{}", x.count, x.bases)).collect();
    format!("{} {}", repeats.len(), items.join(" "))
}
