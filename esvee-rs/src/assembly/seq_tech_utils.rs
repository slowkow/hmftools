//! SeqTechUtils: behaviour specific to a sequencing technology.

use crate::assembly::constants::*;
use crate::assembly::globals::*;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_adjustments::find_low_base_qual_trim_count_ext;
use crate::assembly::read::read_utils::read_junction_extension_length;
use crate::assembly::types::junction::Junction;
use crate::assembly::types::support_read::SupP;
use crate::common::java_sort::java_sort;

pub const ULTIMA_READ_MISMATCH_MEDIUM_REPEAT_COUNT: i32 = 7;
pub const ULTIMA_READ_MISMATCH_LONG_REPEAT_COUNT: i32 = 11;
pub const SBX_PRIME_POSITION_RANGE_THRESHOLD: i32 = 6;
pub const SBX_MEDIUM_QUAL_DESYNC_COUNT: i32 = 150;
pub const UNPAIRED_READ_PROXIMATE_JUNCTION_DISTANCE: i32 = 100;

const ILLUMINA_HIGH_BASE_QUAL: u8 = 30;
const ULTIMA_MAX_QUAL: u8 = 35;
const SBX_MEDIUM_QUAL_LOWER: u8 = 27;
const SBX_MEDIUM_QUAL_UPPER: u8 = 29;

pub fn illumina_is_high_base_qual(qual: u8) -> bool {
    qual >= ILLUMINA_HIGH_BASE_QUAL
}
pub fn ultima_is_high_base_qual(qual: u8) -> bool {
    qual >= ULTIMA_MAX_QUAL
}
pub fn sbx_is_high_base_qual(qual: u8) -> bool {
    qual > SBX_MEDIUM_QUAL_UPPER
}
pub fn sbx_is_medium_base_qual(qual: u8) -> bool {
    (SBX_MEDIUM_QUAL_LOWER..=SBX_MEDIUM_QUAL_UPPER).contains(&qual)
}

pub fn set_seq_tech_specifics() {
    if is_ultima() {
        set_read_mismatch_medium_repeat_count(ULTIMA_READ_MISMATCH_MEDIUM_REPEAT_COUNT);
        set_read_mismatch_long_repeat_count(ULTIMA_READ_MISMATCH_LONG_REPEAT_COUNT);
    }
    set_proximate_junction_distance(UNPAIRED_READ_PROXIMATE_JUNCTION_DISTANCE);
}

pub fn trim_illumina_adapter_bases(read: &ReadP) {
    if !read.is_paired_read() || read.mate_read().is_none() {
        return;
    }
    let mate_read = read.mate_read().unwrap().clone();
    if read.orientation() == mate_read.orientation() {
        return;
    }
    if read.chromosome() != mate_read.chromosome() {
        return;
    }
    if read.orientation().is_forward() {
        let three_prime_end = read.unclipped_end();
        let mate_five_prime_end = mate_read.unclipped_end();
        if three_prime_end > mate_five_prime_end {
            read.m().trim_bases(three_prime_end - mate_five_prime_end, false);
        }
    } else {
        let three_prime_start = read.unclipped_start();
        let mate_five_prime_start = mate_read.unclipped_start();
        if three_prime_start < mate_five_prime_start {
            read.m().trim_bases(mate_five_prime_start - three_prime_start, true);
        }
    }
}

pub fn trim_sbx_uncertain_bases(read: &mut Read) -> bool {
    let mut last_index = read.bases_length() - 1;
    let mut trim_count_start = find_low_base_qual_trim_count_ext(read, 0, last_index, true, true);
    let mut trim_count_end = find_low_base_qual_trim_count_ext(read, 0, last_index, false, true);
    if trim_count_start > 0 && trim_count_end > 0 && trim_count_start + trim_count_end >= read.bases_length() {
        if trim_count_start < trim_count_end {
            read.trim_bases(trim_count_start, true);
            last_index = read.bases_length() - 1;
            trim_count_end = find_low_base_qual_trim_count_ext(read, 0, last_index, false, true);
            if trim_count_end > 0 {
                read.trim_bases(trim_count_end, false);
            }
        } else {
            read.trim_bases(trim_count_end, false);
            last_index = read.bases_length() - 1;
            // as in Java: the start count is recalculated from the end
            trim_count_start = find_low_base_qual_trim_count_ext(read, 0, last_index, false, true);
            if trim_count_start > 0 {
                read.trim_bases(trim_count_start, true);
            }
        }
    } else {
        if trim_count_start > 0 {
            read.trim_bases(trim_count_start, true);
        }
        if trim_count_end > 0 {
            read.trim_bases(trim_count_end, false);
        }
    }
    trim_count_start > 0 || trim_count_end > 0
}

pub fn pass_sbx_distinct_prime_positions_filter(support: &[SupP]) -> bool {
    let mut min_five = -1;
    let mut max_five = -1;
    let mut min_three = -1;
    let mut max_three = -1;
    for read in support {
        let (five, three) = if read.orientation().is_forward() { (read.untrimmed_start(), read.untrimmed_end()) } else { (read.untrimmed_end(), read.untrimmed_start()) };
        if min_five < 0 || five < min_five {
            min_five = five;
        }
        if min_three < 0 || three < min_three {
            min_three = three;
        }
        max_five = max_five.max(five);
        max_three = max_three.max(three);
    }
    (max_five - min_five) + (max_three - min_three) > SBX_PRIME_POSITION_RANGE_THRESHOLD
}

fn within_likely_duplicate_range(first: &Read, second: &Read) -> bool {
    if first.orientation() != second.orientation() {
        return false;
    }
    let first_start = first.unclipped_start() - first.trim_count_start();
    let first_end = first.unclipped_end() + first.trim_count_end();
    let second_start = second.unclipped_start() - second.trim_count_start();
    let second_end = second.unclipped_end() + second.trim_count_end();
    (first_start - second_start).abs() + (first_end - second_end).abs() <= SBX_PRIME_POSITION_RANGE_THRESHOLD
}

pub fn find_sbx_possible_duplicates(junction: &Junction, reads: &[ReadP]) -> Vec<ReadP> {
    let extension_lengths: Vec<i32> = reads.iter().map(|r| read_junction_extension_length(r, junction)).collect();
    let mut sorted_lengths = extension_lengths.clone();
    java_sort(&mut sorted_lengths, |a, b| b.cmp(a));
    let mut min_long_length = -1;
    let mut reads_above_long_length = 0;
    for i in 0..sorted_lengths.len().saturating_sub(1) {
        let ext_length = sorted_lengths[i];
        let next_length = sorted_lengths[i + 1];
        if ext_length as f64 > next_length as f64 * WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR {
            min_long_length = ext_length;
            reads_above_long_length = i + 1;
            break;
        }
    }
    if reads_above_long_length < 2 {
        return Vec::new();
    }
    let mut candidate_duplicates: Vec<ReadP> = Vec::new();
    for r in 0..reads.len() {
        if extension_lengths[r] >= min_long_length {
            candidate_duplicates.push(reads[r].clone());
        }
    }
    let mut likely_duplicates: Vec<ReadP> = Vec::new();
    for i in 0..candidate_duplicates.len().saturating_sub(1) {
        let read = &candidate_duplicates[i];
        if likely_duplicates.contains(read) {
            continue;
        }
        for next_read in candidate_duplicates.iter().skip(i + 1) {
            if within_likely_duplicate_range(read, next_read) {
                likely_duplicates.push(next_read.clone());
            }
        }
    }
    likely_duplicates
}

pub fn set_sbx_medium_qual_count(read: &Read) -> i32 {
    read.get_base_quality().iter().filter(|q| crate::assembly::common_utils::is_medium_base_qual(**q)).count() as i32
}
