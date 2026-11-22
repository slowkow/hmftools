//! HomologyData: the homology at a breakend, from overlapping alignments or from the reference around an indel.

use crate::assembly::alignment::align_data::AlignData;
use crate::assembly::alignment::md_tag::{MdTag, MATCH_BYTE};
use crate::assembly::assembly_utils::reverse_complement_str;
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HomologyData {
    pub homology: String,
    pub exact_start: i32,
    pub exact_end: i32,
    pub inexact_start: i32,
    pub inexact_end: i32,
    /// true for the shared NO_HOMOLOGY instance of the Java class (which is compared by identity)
    no_homology: bool,
}

impl HomologyData {
    pub fn no_homology() -> HomologyData {
        HomologyData { homology: String::new(), exact_start: 0, exact_end: 0, inexact_start: 0, inexact_end: 0, no_homology: true }
    }
    pub fn is_no_homology(&self) -> bool {
        self.no_homology
    }
    pub fn new(homology: &str, exact_start: i32, exact_end: i32, inexact_start: i32, inexact_end: i32) -> HomologyData {
        HomologyData { homology: homology.to_string(), exact_start, exact_end, inexact_start, inexact_end, no_homology: false }
    }
    pub fn length(&self) -> i32 {
        self.inexact_end.abs() + self.inexact_start.abs()
    }
    pub fn exists(&self) -> bool {
        self.inexact_end != 0 || self.inexact_start != 0
    }
    pub fn is_symmetrical(&self) -> bool {
        self.inexact_start.abs() == self.inexact_end
    }
    pub fn position_adjustment(&self, orientation: Orientation) -> i32 {
        if orientation.is_forward() { -self.inexact_end } else { self.inexact_start.abs() }
    }
    pub fn invert(&self, reverse_positions: bool, reverse_bases: bool) -> HomologyData {
        HomologyData::new(
            &if reverse_bases { reverse_complement_str(&self.homology) } else { self.homology.clone() },
            if reverse_positions { -self.exact_end } else { self.exact_start },
            if reverse_positions { self.exact_start.abs() } else { self.exact_end },
            if reverse_positions { -self.inexact_end } else { self.inexact_start },
            if reverse_positions { self.inexact_start.abs() } else { self.inexact_end },
        )
    }

    pub fn determine_homology(full_sequence: &str, left_alignment: &AlignData, right_alignment: &AlignData) -> Option<HomologyData> {
        let overlap_length = left_alignment.sequence_end() - right_alignment.sequence_start() + 1;
        if overlap_length <= 0 {
            return None;
        }
        let left_index_end = left_alignment.segment_length() - 1;
        let left_index_start = left_index_end - overlap_length + 1;
        let right_index_start = 0;
        let right_index_end = overlap_length - 1;
        let mut overlap_bases = full_sequence[right_alignment.sequence_start() as usize..(left_alignment.sequence_end() + 1) as usize].to_string();
        if left_alignment.is_reverse() {
            overlap_bases = reverse_complement_str(&overlap_bases);
        }
        let left_md_tag = MdTag::new(left_alignment.md_tag());
        let right_md_tag = MdTag::new(right_alignment.md_tag());
        if !left_md_tag.has_mismatches() && !right_md_tag.has_mismatches() {
            let half_overlap = overlap_length / 2;
            let exact_start = if overlap_length % 2 == 0 { half_overlap } else { half_overlap + 1 };
            let exact_end = overlap_length - exact_start;
            return Some(HomologyData::new(&overlap_bases, -exact_start, exact_end, -exact_start, exact_end));
        }
        let left_md_seq = left_md_tag.extract_sub_sequence(left_index_start, left_index_end, left_alignment.is_reverse());
        let right_md_seq = right_md_tag.extract_sub_sequence(right_index_start, right_index_end, right_alignment.is_reverse());
        if left_md_seq.len() != right_md_seq.len() {
            return None;
        }
        let n = left_md_seq.len();
        let mut index_mismatches = vec![0i32; n + 1];
        let mut min_mismatch_count = -1;
        let mut min_mismatch_end_index: i32 = -1;
        for end_index in 0..=n {
            let mut mismatch_count = 0;
            for i in 0..n {
                if i < end_index {
                    if left_md_seq[i] != MATCH_BYTE {
                        mismatch_count += 1;
                    }
                } else if right_md_seq[i] != MATCH_BYTE {
                    mismatch_count += 1;
                }
            }
            index_mismatches[end_index] = mismatch_count;
            if min_mismatch_count == -1 || mismatch_count <= min_mismatch_count {
                min_mismatch_count = mismatch_count;
                min_mismatch_end_index = end_index as i32;
            }
        }
        if min_mismatch_count == 0 && min_mismatch_end_index == overlap_length - 1 {
            let half_overlap = overlap_length / 2;
            let exact_start = if overlap_length % 2 == 0 { half_overlap } else { half_overlap + 1 };
            let exact_end = overlap_length - exact_start;
            return Some(HomologyData::new(&overlap_bases, -exact_start, exact_end, -exact_start, exact_end));
        }
        let mut range_start: i32 = -1;
        let mut range_end: i32 = -1;
        let mut longest_range_start: i32 = -1;
        let mut longest_range_end: i32 = -1;
        for (i, m) in index_mismatches.iter().enumerate() {
            let i = i as i32;
            if *m == min_mismatch_count {
                if range_start < 0 {
                    range_start = i;
                }
            } else if range_start >= 0 {
                range_end = i - 1;
                if longest_range_start < 0 || (range_end - range_start + 1 > longest_range_end - longest_range_start + 1) {
                    longest_range_start = range_start;
                    longest_range_end = range_end;
                }
                range_start = -1;
                range_end = -1;
            }
        }
        if range_start >= 0 && range_end < 0 {
            longest_range_start = range_start;
            longest_range_end = index_mismatches.len() as i32 - 1;
        }
        let range_length = longest_range_end - longest_range_start + 1;
        let half_range = range_length / 2;
        let mut exact_start = 0;
        let mut exact_end = 0;
        let mut exact_homology = String::new();
        if longest_range_start >= 0 && longest_range_end > longest_range_start {
            exact_start = half_range.max(1);
            exact_end = range_length - exact_start - 1;
            exact_homology = overlap_bases[longest_range_start as usize..longest_range_end as usize].to_string();
        }
        let inexact_start = longest_range_start + half_range;
        let inexact_end = overlap_length - inexact_start;
        Some(HomologyData::new(&exact_homology, -exact_start, exact_end, -inexact_start, inexact_end))
    }

    pub fn determine_indel_homology(ref_bases_start: &str, ref_bases_end: &str, overlap: i32) -> Option<HomologyData> {
        if ref_bases_start.len() != ref_bases_end.len() {
            return None;
        }
        if ref_bases_start == ref_bases_end {
            let exact_start = (overlap + 1) / 2;
            let exact_end = overlap - exact_start;
            return Some(HomologyData::new(ref_bases_start, -exact_start, exact_end, -exact_start, exact_end));
        }
        let (s, e) = (ref_bases_start.as_bytes(), ref_bases_end.as_bytes());
        let mut exact_match = 0;
        while exact_match < overlap && s[exact_match as usize] == e[exact_match as usize] {
            exact_match += 1;
        }
        if exact_match == 0 {
            return Some(HomologyData::no_homology());
        }
        let exact_start = (exact_match + 1) / 2;
        let exact_end = exact_match - exact_start;
        let mut score_acc = exact_match * INDEL_HOMOLOGY_MATCH_SCORE;
        let mut max_score = score_acc + if exact_match == overlap { 0 } else { INDEL_HOMOLOGY_CLIP_SCORE };
        let mut max_inexact_match = exact_match;
        for i in exact_match..overlap {
            if s[i as usize] == e[i as usize] {
                score_acc += INDEL_HOMOLOGY_MATCH_SCORE;
            } else {
                score_acc += INDEL_HOMOLOGY_MISMATCH_SCORE;
            }
            let remaining = overlap - i - 1;
            let score = score_acc + if remaining > 0 { INDEL_HOMOLOGY_CLIP_SCORE } else { 0 };
            if score >= max_score {
                max_score = score;
                max_inexact_match = i + 1;
            } else if score_acc + remaining * INDEL_HOMOLOGY_MATCH_SCORE < max_score {
                break;
            }
        }
        let inexact_start = (max_inexact_match + 1) / 2;
        let inexact_end = max_inexact_match - inexact_start;
        Some(HomologyData::new(&ref_bases_start[..max_inexact_match as usize], -exact_start, exact_end, -inexact_start, inexact_end))
    }
}

impl std::fmt::Display for HomologyData {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} exact({},{}) inexact({},{})",
            if self.homology.is_empty() { "none" } else { &self.homology },
            self.exact_start,
            self.exact_end,
            self.inexact_start,
            self.inexact_end
        )
    }
}
