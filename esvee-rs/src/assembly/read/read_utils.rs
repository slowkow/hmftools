//! ReadUtils.

use crate::assembly::common_utils;
use crate::assembly::globals::{max_concordant_fragment_length, max_observed_concordant_frag_length};
use crate::assembly::orientation::Orientation;
use crate::assembly::read::read::Read;
use crate::assembly::types::junction::Junction;
use crate::assembly::types::repeat_info::{calc_trimmed_base_length, find_repeats};

pub const INVALID_INDEX: i32 = -1;

pub fn is_discordant_fragment(read: &Read) -> bool {
    common_utils::is_discordant_fragment(read.bam_record(), max_concordant_fragment_length(max_observed_concordant_frag_length()), read.supplementary_data())
}

pub fn is_valid_support_coords_vs_junction(read: &Read, is_forward_junction: bool, junction_position: i32) -> bool {
    if is_forward_junction {
        if read.negative_strand() || read.is_right_clipped() {
            return false;
        }
        if read.alignment_end() > junction_position {
            return false;
        }
    } else {
        if read.positive_strand() || read.is_left_clipped() {
            return false;
        }
        if read.alignment_start() < junction_position {
            return false;
        }
    }
    true
}

pub fn get_read_index_at_reference_position(read: &Read, ref_position: i32, allow_extrapolation: bool) -> i32 {
    get_read_index_at_reference_position_orient(read, ref_position, allow_extrapolation, None)
}

pub fn get_read_index_at_junction(read: &Read, junction: &Junction, allow_extrapolation: bool) -> i32 {
    if junction.indel_based() && read.indel_coords().is_some() {
        return read_index_from_position(read, junction.position);
    }
    get_read_index_at_reference_position_orient(read, junction.position, allow_extrapolation, Some(junction.orient))
}

pub fn get_read_index_at_reference_position_orient(read: &Read, ref_position: i32, allow_extrapolation: bool, required_orientation: Option<Orientation>) -> i32 {
    let alignment_start = read.alignment_start();
    let alignment_end = read.alignment_end();
    if required_orientation.is_none_or(|o| o.is_reverse()) && read.has_indel_implied_unclipped_start() {
        if read.indel_implied_unclipped_start() > ref_position {
            return INVALID_INDEX;
        }
        return ref_position - read.indel_implied_unclipped_start();
    }
    if required_orientation.is_none_or(|o| o.is_forward()) && read.has_indel_implied_unclipped_end() {
        if read.indel_implied_unclipped_end() < ref_position {
            return INVALID_INDEX;
        }
        return read.bases_length() - 1 - (read.indel_implied_unclipped_end() - ref_position);
    }
    if ref_position <= alignment_start {
        if !allow_extrapolation && ref_position < alignment_start {
            return INVALID_INDEX;
        }
        let base_diff = alignment_start - ref_position;
        let soft_clip_bases = alignment_start - read.unclipped_start();
        return if base_diff <= soft_clip_bases { soft_clip_bases - base_diff } else { INVALID_INDEX };
    } else if ref_position >= alignment_end {
        if !allow_extrapolation && ref_position > alignment_end {
            return INVALID_INDEX;
        }
        let base_diff = ref_position - alignment_end;
        let soft_clip_bases = read.unclipped_end() - alignment_end;
        return if base_diff <= soft_clip_bases { read.bases_length() - (soft_clip_bases - base_diff) - 1 } else { INVALID_INDEX };
    }
    read_index_from_position(read, ref_position)
}

pub fn read_index_from_position(read: &Read, ref_position: i32) -> i32 {
    let mut read_index = 0;
    let mut current_pos = read.alignment_start();
    for element in read.cigar_elements() {
        if !element.op.consumes_reference_bases() {
            if element.op.consumes_read_bases() {
                read_index += element.len;
            }
            continue;
        }
        if current_pos == ref_position {
            break;
        }
        if !element.op.consumes_read_bases() {
            if ref_position >= current_pos && ref_position < current_pos + element.len {
                return read_index - 1;
            }
            current_pos += element.len;
        } else {
            if ref_position >= current_pos && ref_position < current_pos + element.len {
                return read_index + ref_position - current_pos;
            }
            current_pos += element.len;
            read_index += element.len;
        }
    }
    read_index
}

pub fn record_soft_clips_at_junction(read: &Read, junction: &Junction) -> bool {
    if junction.is_forward() {
        read.is_right_clipped() && read.alignment_end() == junction.position
    } else {
        read.is_left_clipped() && read.alignment_start() == junction.position
    }
}

pub fn read_junction_extension_length(read: &Read, junction: &Junction) -> i32 {
    if junction.is_forward() {
        if read.is_right_clipped() { (read.max_unclipped_end() - junction.position).max(0) } else { 0 }
    } else if read.is_left_clipped() {
        (junction.position - read.min_unclipped_start()).max(0)
    } else {
        0
    }
}

pub fn avg_base_quality_read(read: &Read) -> i32 {
    avg_base_quality(read.get_base_quality(), 0, read.bases_length() - 1)
}

pub fn avg_base_quality(base_qualities: &[u8], index_start: i32, index_end: i32) -> i32 {
    if index_start > index_end || index_start < 0 || index_end >= base_qualities.len() as i32 {
        return -1;
    }
    let mut quality_sum = 0i32;
    for i in index_start..=index_end {
        quality_sum += base_qualities[i as usize] as i8 as i32;
    }
    quality_sum / (index_end - index_start + 1)
}

pub fn calc_trimmed_read_base_length(read: &Read, index_start: i32, index_end: i32) -> i32 {
    let from = index_start.max(0) as usize;
    let to = (index_end + 1).min(read.bases_length()) as usize;
    let read_bases = &read.get_bases()[from..to];
    let repeats = find_repeats(read_bases);
    calc_trimmed_base_length(0, read_bases.len() as i32 - 1, &repeats)
}
