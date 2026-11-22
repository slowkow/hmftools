//! SequenceCompare: compares two base sequences allowing for low quality bases, short indels and repeat
//! count differences.

use crate::assembly::assembly_utils::bases_match_qual;
use crate::assembly::common_utils::{below_min_qual, is_medium_base_qual};
use crate::assembly::constants::*;
use crate::assembly::globals::{read_mismatch_long_repeat_count, read_mismatch_medium_repeat_count};
use crate::assembly::sequence_builder::permitted_read_mismatches;
use crate::assembly::sequence_diff::BaseQualType;
use crate::assembly::types::junction_assembly::JunctionAssembly;
use crate::assembly::types::repeat_info::RepeatInfo;

const MISMATCH_RECOVERY_BASE_COUNT: i32 = 3;

pub fn matched_assembly_sequences(first: &JunctionAssembly, second: &JunctionAssembly) -> bool {
    let mut first_index_start;
    let mut first_index_end;
    let mut second_index_start;
    let mut second_index_end;
    let junction_diff = first.junction().position - second.junction().position;
    if first.junction().is_forward() {
        if junction_diff > 0 {
            first_index_start = first.junction_index() - junction_diff;
            second_index_start = second.junction_index();
        } else {
            first_index_start = first.junction_index();
            second_index_start = second.junction_index() + junction_diff;
        }
        let min_distance_from_junction = first.upper_distance_from_junction().min(second.upper_distance_from_junction());
        first_index_end = first_index_start + min_distance_from_junction;
        second_index_end = second_index_start + min_distance_from_junction;
    } else {
        let min_distance_from_junction = first.lower_distance_from_junction().min(second.lower_distance_from_junction());
        if junction_diff > 0 {
            first_index_end = first.junction_index();
            second_index_end = second.junction_index() + junction_diff;
        } else {
            first_index_end = first.junction_index() - junction_diff;
            second_index_end = second.junction_index();
        }
        first_index_start = first_index_end - min_distance_from_junction;
        second_index_start = second_index_end - min_distance_from_junction;
    }
    if first_index_start < 0 || second_index_start < 0 {
        return false;
    }
    first_index_end = first_index_end.min(first.bases().len() as i32 - 1);
    second_index_end = second_index_end.min(second.bases().len() as i32 - 1);
    if first_index_end <= first_index_start || second_index_end <= second_index_start {
        return false;
    }
    let first_full = first.bases();
    let first_overlap = &first_full[first_index_start as usize..=first_index_end as usize];
    let second_full = second.bases();
    if let Some(second_ext_base_index) = index_of(second_full, first_overlap) {
        let second_ext_base_end_index = second_ext_base_index + first_overlap.len() as i32 - 1;
        if second_ext_base_index < second.junction_index() && second_ext_base_end_index > second.junction_index() {
            return true;
        }
    }
    let second_overlap = &second_full[second_index_start as usize..=second_index_end as usize];
    if let Some(first_ext_base_index) = index_of(first_full, second_overlap) {
        let first_ext_base_end_index = first_ext_base_index + second_overlap.len() as i32 - 1;
        if first_ext_base_index < first.junction_index() && first_ext_base_end_index > first.junction_index() {
            return true;
        }
    }
    let overlap_length = (first_index_end - first_index_start + 1).min(second_index_end - second_index_start + 1);
    let max_mismatch_penalty = permitted_read_mismatches(overlap_length).max(PRIMARY_ASSEMBLY_MERGE_MISMATCH as f64);
    let mismatch_penalty = compare_sequences(
        first.bases(), first.base_quals(), first_index_start, first_index_end, first.repeat_info(),
        second.bases(), second.base_quals(), second_index_start, second_index_end, second.repeat_info(), max_mismatch_penalty,
    );
    mismatch_penalty <= max_mismatch_penalty
}

/// String.indexOf for byte sequences
pub fn index_of(haystack: &[u8], needle: &[u8]) -> Option<i32> {
    if needle.is_empty() {
        return Some(0);
    }
    memchr::memmem::find(haystack, needle).map(|i| i as i32)
}

#[allow(clippy::too_many_arguments)]
pub fn compare_sequences(
    first_bases: &[u8], first_base_quals: &[u8], first_index_start: i32, first_index_end: i32, first_repeats: &[RepeatInfo],
    second_bases: &[u8], second_base_quals: &[u8], second_index_start: i32, second_index_end: i32, second_repeats: &[RepeatInfo], max_mismatch_penalty: f64,
) -> f64 {
    compare_sequences_mismatch_spec(
        first_bases, first_base_quals, first_index_start, first_index_end, first_repeats, second_bases, second_base_quals, second_index_start,
        second_index_end, second_repeats, max_mismatch_penalty, true, true,
    )
}

#[allow(clippy::too_many_arguments)]
fn compare_sequences_mismatch_spec(
    first_bases: &[u8], first_base_quals: &[u8], first_index_start: i32, first_index_end: i32, first_repeats: &[RepeatInfo],
    second_bases: &[u8], second_base_quals: &[u8], second_index_start: i32, second_index_end: i32, second_repeats: &[RepeatInfo],
    max_mismatch_penalty: f64, check_forwards: bool, repeats_are_mismatches: bool,
) -> f64 {
    if first_index_start < 0
        || first_index_end >= first_bases.len() as i32
        || first_index_end >= first_base_quals.len() as i32
        || second_index_start < 0
        || second_index_end >= second_bases.len() as i32
        || second_index_end >= second_base_quals.len() as i32
    {
        return if max_mismatch_penalty < 0.0 { max_mismatch_penalty } else { max_mismatch_penalty + 1.0 };
    }
    let mut mismatch_penalty = 0.0;
    let mut first_index = if check_forwards { first_index_start } else { first_index_end };
    let mut second_index = if check_forwards { second_index_start } else { second_index_end };
    let mut last_repeat_skip_bases = 0;
    let step = if check_forwards { 1 } else { -1 };
    while first_index >= 0 && first_index <= first_index_end && second_index >= 0 && second_index <= second_index_end {
        let (fi, si) = (first_index as usize, second_index as usize);
        if bases_match_qual(first_bases[fi], second_bases[si], first_base_quals[fi], second_base_quals[si]) {
            first_index += step;
            second_index += step;
            continue;
        }
        let recovery_skip_bases = check_skip_short_indel(first_bases, first_index, second_bases, second_index, check_forwards);
        if recovery_skip_bases != 0 {
            let skip_count = recovery_skip_bases.abs();
            let first_range_start = if check_forwards { first_index } else { first_index - skip_count };
            let first_range_end = if check_forwards { first_index + skip_count } else { first_index };
            let first_qual_type = range_qual_type(first_base_quals, first_range_start, first_range_end);
            let second_range_start = if check_forwards { second_index } else { second_index - skip_count };
            let second_range_end = if check_forwards { second_index + skip_count } else { second_index };
            let second_qual_type = range_qual_type(second_base_quals, second_range_start, second_range_end);
            if recovery_skip_bases > 0 {
                first_index += if check_forwards { recovery_skip_bases } else { -recovery_skip_bases };
            } else {
                second_index += if check_forwards { -recovery_skip_bases } else { recovery_skip_bases };
            }
            if last_repeat_skip_bases == 0 {
                if (recovery_skip_bases == 1 && matches_repeat(second_index, second_repeats, check_forwards))
                    || (recovery_skip_bases == -1 && matches_repeat(first_index, first_repeats, check_forwards))
                {
                    last_repeat_skip_bases = recovery_skip_bases;
                }
            } else if last_repeat_skip_bases + recovery_skip_bases == 0 {
                last_repeat_skip_bases = 0;
                continue;
            }
            mismatch_penalty += calc_non_snv_mismatch_penalty(first_qual_type, second_qual_type);
            continue;
        }
        let expected_repeat_base_diff = check_repeat_difference(first_index, first_repeats, second_index, second_repeats, check_forwards);
        if expected_repeat_base_diff != 0 {
            last_repeat_skip_bases = 0;
            if expected_repeat_base_diff > 0 {
                // as in Java, the same step in both directions
                first_index += expected_repeat_base_diff;
            } else {
                second_index += if check_forwards { -expected_repeat_base_diff } else { expected_repeat_base_diff };
            }
            if repeats_are_mismatches {
                let skip_count = expected_repeat_base_diff.abs();
                let first_range_start = if check_forwards { first_index } else { first_index - skip_count };
                let first_range_end = if check_forwards { first_index + skip_count } else { first_index };
                let first_qual_type = range_qual_type(first_base_quals, first_range_start, first_range_end);
                let second_range_start = if check_forwards { second_index } else { second_index - skip_count };
                let second_range_end = if check_forwards { second_index + skip_count } else { second_index };
                let second_qual_type = range_qual_type(second_base_quals, second_range_start, second_range_end);
                mismatch_penalty += calc_non_snv_mismatch_penalty(first_qual_type, second_qual_type);
            }
            continue;
        }
        mismatch_penalty += calc_snv_mismatch_penalty(first_base_quals[fi], second_base_quals[si]);
        if max_mismatch_penalty >= 0.0 && mismatch_penalty > max_mismatch_penalty {
            return mismatch_penalty;
        }
        first_index += step;
        second_index += step;
    }
    mismatch_penalty
}

fn calc_non_snv_mismatch_penalty(first: BaseQualType, second: BaseQualType) -> f64 {
    let lower = BaseQualType::select_lower(first, second);
    if lower == BaseQualType::Low {
        return 0.0;
    }
    if lower == BaseQualType::Medium { READ_MISMATCH_MED_QUAL_NON_SNV_PENALTY } else { READ_MISMATCH_PENALTY }
}

fn calc_snv_mismatch_penalty(first_qual: u8, second_qual: u8) -> f64 {
    if below_min_qual(first_qual) || below_min_qual(second_qual) {
        return 0.0;
    }
    if is_medium_base_qual(first_qual) || is_medium_base_qual(second_qual) {
        return READ_MISMATCH_MED_QUAL_SNV_PENALTY;
    }
    READ_MISMATCH_PENALTY
}

fn range_qual_type(base_quals: &[u8], range_start: i32, range_end: i32) -> BaseQualType {
    let mut has_medium = false;
    let range_start = range_start.max(0);
    let range_end = range_end.min(base_quals.len() as i32 - 1);
    let mut i = range_start;
    while i <= range_end {
        let q = base_quals[i as usize];
        if below_min_qual(q) {
            return BaseQualType::Low;
        }
        has_medium |= is_medium_base_qual(q);
        i += 1;
    }
    if has_medium { BaseQualType::Medium } else { BaseQualType::High }
}

pub fn permitted_repeat_count(repeat_count: i32) -> i32 {
    if repeat_count >= read_mismatch_long_repeat_count() {
        return 3;
    }
    if repeat_count >= read_mismatch_medium_repeat_count() {
        return 2;
    }
    1
}

fn check_repeat_difference(first_index: i32, first_repeats: &[RepeatInfo], second_index: i32, second_repeats: &[RepeatInfo], check_forwards: bool) -> i32 {
    let first_repeat = find_repeat(first_index, first_repeats, check_forwards);
    let second_repeat = find_repeat(second_index, second_repeats, check_forwards);
    let (Some(first_repeat), Some(second_repeat)) = (first_repeat, second_repeat) else { return 0 };
    if !first_repeat.matches_type(second_repeat) {
        if first_repeat.bases.contains(&second_repeat.bases) && first_repeat.count == 2 {
            let ratio = first_repeat.repeat_length() / second_repeat.repeat_length();
            if second_repeat.count != first_repeat.count * ratio {
                return 0;
            }
            return first_repeat.bases.len() as i32;
        } else if second_repeat.bases.contains(&first_repeat.bases) && second_repeat.count == 2 {
            let ratio = second_repeat.repeat_length() / first_repeat.repeat_length();
            if first_repeat.count != second_repeat.count * ratio {
                return 0;
            }
            return second_repeat.bases.len() as i32;
        }
        return 0;
    }
    if first_repeat.count == second_repeat.count {
        return 0;
    }
    let permitted_repeat_diff = permitted_repeat_count(first_repeat.count.max(second_repeat.count));
    if (first_repeat.count - second_repeat.count).abs() > permitted_repeat_diff {
        return 0;
    }
    (first_repeat.count - second_repeat.count) * first_repeat.bases.len() as i32
}

fn find_repeat(current_index: i32, repeats: &[RepeatInfo], check_forwards: bool) -> Option<&RepeatInfo> {
    for repeat in repeats {
        if repeat.index < current_index && repeat.post_repeat_index() > current_index {
            return Some(repeat);
        }
        if check_forwards {
            if repeat.post_repeat_index() == current_index {
                return Some(repeat);
            }
        } else if repeat.index == current_index {
            return Some(repeat);
        }
    }
    None
}

fn matches_repeat(current_index: i32, repeats: &[RepeatInfo], check_forwards: bool) -> bool {
    if check_forwards {
        repeats.iter().any(|x| x.index == current_index || x.index == current_index + 1)
    } else {
        repeats.iter().any(|x| x.post_repeat_index() == current_index || x.post_repeat_index() == current_index - 1)
    }
}

fn check_skip_short_indel(first_bases: &[u8], first_index: i32, second_bases: &[u8], second_index: i32, check_forwards: bool) -> i32 {
    for diff in 1..=2 {
        let adjusted_first_index = first_index + if check_forwards { diff } else { -diff };
        if can_recover_match(first_bases, adjusted_first_index, second_bases, second_index, check_forwards) {
            return diff;
        }
    }
    for diff in 1..=2 {
        let adjusted_second_index = second_index + if check_forwards { diff } else { -diff };
        if can_recover_match(first_bases, first_index, second_bases, adjusted_second_index, check_forwards) {
            return -diff;
        }
    }
    0
}

fn can_recover_match(first_bases: &[u8], first_current_index: i32, second_bases: &[u8], second_current_start: i32, check_forwards: bool) -> bool {
    let mut first_index = first_current_index;
    let mut second_index = second_current_start;
    let first_index_end;
    let second_index_end;
    if check_forwards {
        first_index_end = first_current_index + MISMATCH_RECOVERY_BASE_COUNT;
        second_index_end = second_current_start + MISMATCH_RECOVERY_BASE_COUNT;
        if first_index_end >= first_bases.len() as i32 || second_index_end >= second_bases.len() as i32 {
            return false;
        }
    } else {
        first_index_end = first_current_index - MISMATCH_RECOVERY_BASE_COUNT;
        second_index_end = second_current_start - MISMATCH_RECOVERY_BASE_COUNT;
        if first_index_end < 0 || second_index_end < 0 {
            return false;
        }
    }
    loop {
        if first_bases[first_index as usize] != second_bases[second_index as usize] {
            return false;
        }
        if first_index == first_index_end && second_index == second_index_end {
            break;
        }
        first_index += if check_forwards { 1 } else { -1 };
        second_index += if check_forwards { 1 } else { -1 };
    }
    true
}

pub fn compare_sequences_mismatch_no_quals(
    first_bases: &[u8], first_index_start: i32, first_index_end: i32, second_bases: &[u8], second_index_start: i32, second_index_end: i32, max_mismatch_penalty: i32,
) -> i32 {
    if first_index_start < 0 || first_index_end >= first_bases.len() as i32 || second_index_start < 0 || second_index_end >= second_bases.len() as i32 {
        return if max_mismatch_penalty < 0 { max_mismatch_penalty } else { max_mismatch_penalty + 1 };
    }
    let mut mismatch_total = 0;
    let mut first_index = first_index_start;
    let mut second_index = second_index_start;
    let last_repeat_skip_bases = 0;
    while first_index >= 0 && first_index <= first_index_end && second_index >= 0 && second_index <= second_index_end {
        if first_bases[first_index as usize] == second_bases[second_index as usize] {
            first_index += 1;
            second_index += 1;
            continue;
        }
        let recovery_skip_bases = check_skip_short_indel(first_bases, first_index, second_bases, second_index, true);
        if recovery_skip_bases != 0 {
            if recovery_skip_bases > 0 {
                first_index += recovery_skip_bases;
            } else {
                second_index += -recovery_skip_bases;
            }
            if last_repeat_skip_bases + recovery_skip_bases == 0 {
                continue;
            }
            mismatch_total += 1;
            continue;
        }
        mismatch_total += 1;
        if max_mismatch_penalty >= 0 && mismatch_total > max_mismatch_penalty {
            return mismatch_total;
        }
        first_index += 1;
        second_index += 1;
    }
    mismatch_total
}
