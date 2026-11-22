//! PhaseSetMerger: merges the assembly alignments of phase sets whose full sequences overlap.

use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::assembly_utils::{create_base_quals_above_min_threshold, mismatches_per_comparison_length, reverse_complement_str};
use crate::assembly::common_utils::{is_line_insert_pair, position_within};
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::sequence_compare::{compare_sequences, index_of};
use crate::assembly::types::junction_sequence::PHASED_ASSEMBLY_MATCH_SEQ_LENGTH;
use crate::assembly::types::phase_set::PhaseSetP;
use crate::assembly::types::repeat_info::{find_repeats, RepeatInfo};

pub const MATCH_SUBSEQUENCE_LENGTH: i32 = 50;

pub fn merge_phase_sets(phase_sets: &mut Vec<PhaseSetP>) {
    if phase_sets.len() == 1 {
        return;
    }
    let n = phase_sets.len();
    for i in 0..n.saturating_sub(1) {
        let first = phase_sets[i];
        if !first.has_valid_assembly_alignment() || first.is_secondary_line_link() {
            continue;
        }
        for second in phase_sets.iter().skip(i + 1).copied() {
            if !second.has_valid_assembly_alignment() || second.is_secondary_line_link() {
                continue;
            }
            try_merge_phase_sets(first, second);
        }
    }
}

pub fn try_merge_phase_sets(first_phase_set: PhaseSetP, second_phase_set: PhaseSetP) {
    let first_assembly = first_phase_set.assembly_alignment().expect("assembly alignment");
    let second_assembly = second_phase_set.assembly_alignment().expect("assembly alignment");
    let first_is_primary = first_assembly.full_sequence_length() >= second_assembly.full_sequence_length();
    let (primary_assembly, other_assembly, primary_phase_set, other_phase_set) =
        if first_is_primary { (first_assembly, second_assembly, first_phase_set, second_phase_set) } else { (second_assembly, first_assembly, second_phase_set, first_phase_set) };
    let primary_sequence = primary_assembly.full_sequence().expect("full sequence").to_string();
    let other_full_sequence = other_assembly.full_sequence().expect("full sequence").to_string();
    for other_orientation in [Orientation::Forward, Orientation::Reverse] {
        let other_sequence = if other_orientation.is_forward() { other_full_sequence.clone() } else { reverse_complement_str(&other_full_sequence) };
        if let Some(match_indices) = find_merge_indices(primary_assembly, &primary_sequence, other_assembly, &other_sequence, other_orientation.is_reverse()) {
            merge_assembly_alignments(primary_assembly, other_assembly, &primary_sequence, &other_sequence, other_orientation, match_indices);
            primary_phase_set.m().merge_phase_set(other_phase_set);
            other_phase_set.m().set_merged_phase_set_id(primary_phase_set.id());
            return;
        }
    }
}

fn find_merge_indices(first_assembly: AsmAlignP, first_sequence: &str, second_assembly: AsmAlignP, second_sequence: &str, second_reversed: bool) -> Option<[i32; 2]> {
    let (first_bytes, second_bytes) = (first_sequence.as_bytes(), second_sequence.as_bytes());
    if let Some(match_index) = index_of(first_bytes, second_bytes) {
        return Some([match_index, 0]);
    }
    let first_length = first_bytes.len() as i32;
    let second_length = second_bytes.len() as i32;
    let mut prepared: Option<(Vec<u8>, Vec<RepeatInfo>, Vec<u8>, Vec<RepeatInfo>)> = None;
    let sub_seq_iterations = first_length / MATCH_SUBSEQUENCE_LENGTH;
    for i in 0..sub_seq_iterations {
        let match_seq_start_index = i * MATCH_SUBSEQUENCE_LENGTH;
        let match_seq_end_index = match_seq_start_index + MATCH_SUBSEQUENCE_LENGTH;
        if match_seq_end_index >= first_length + 1 {
            break;
        }
        let first_sub_sequence = &first_bytes[match_seq_start_index as usize..match_seq_end_index as usize];
        let Some(second_sub_seq_index) = index_of(second_bytes, first_sub_sequence) else { continue };
        let min_start_distance = match_seq_start_index.min(second_sub_seq_index);
        let first_match_index_start = match_seq_start_index - min_start_distance;
        let second_match_index_start = second_sub_seq_index - min_start_distance;
        let overlap_length = (first_length - first_match_index_start - 1).min(second_length - second_match_index_start - 1);
        if overlap_length < PHASED_ASSEMBLY_MATCH_SEQ_LENGTH {
            continue;
        }
        let first_match_index_end = first_match_index_start + overlap_length - 1;
        let second_match_index_end = second_match_index_start + overlap_length - 1;
        let (first_quals, first_repeats, second_quals, second_repeats) = prepared.get_or_insert_with(|| {
            (create_base_quals_above_min_threshold(first_bytes.len()), find_repeats(first_bytes), create_base_quals_above_min_threshold(second_bytes.len()), find_repeats(second_bytes))
        });
        let permitted_mismatches = mismatches_per_comparison_length(overlap_length);
        let mismatch_penalty = compare_sequences(
            first_bytes, first_quals, first_match_index_start, first_match_index_end, first_repeats,
            second_bytes, second_quals, second_match_index_start, second_match_index_end, second_repeats, permitted_mismatches as f64,
        );
        if mismatch_penalty > permitted_mismatches as f64 {
            continue;
        }
        if !can_merge_phase_sets(first_assembly, first_match_index_start, first_match_index_end, second_assembly, second_match_index_start, second_match_index_end, second_reversed) {
            continue;
        }
        return Some([first_match_index_start, second_match_index_start]);
    }
    None
}

fn can_merge_phase_sets(first: AsmAlignP, first_overlap_start: i32, first_overlap_end: i32, second: AsmAlignP, second_overlap_start: i32, second_overlap_end: i32, second_reversed: bool) -> bool {
    if match_sequence_overlaps_junctions(first, first_overlap_start, first_overlap_end, false) {
        return true;
    }
    if match_sequence_overlaps_junctions(second, second_overlap_start, second_overlap_end, second_reversed) {
        return true;
    }
    for first_assembly in first.assemblies() {
        let first_is_line_site = first_assembly.has_line_sequence();
        for second_assembly in second.assemblies() {
            if (first_is_line_site || second_assembly.has_line_sequence()) && is_line_insert_pair(&first_assembly, &second_assembly) {
                return true;
            }
        }
    }
    false
}

fn match_sequence_overlaps_junctions(assembly_alignment: AsmAlignP, overlap_start: i32, overlap_end: i32, reversed: bool) -> bool {
    let junc_boundary_start = overlap_start + ASSEMBLY_LINK_OVERLAP_BASES;
    let junc_boundary_end = overlap_end - ASSEMBLY_LINK_OVERLAP_BASES;
    for link_index in assembly_alignment.link_indices() {
        let adjusted_link_index = if !reversed { *link_index } else { assembly_alignment.full_sequence_length() - link_index };
        if position_within(adjusted_link_index, junc_boundary_start, junc_boundary_end) {
            return true;
        }
    }
    false
}

pub fn merge_assembly_alignments(primary_assembly: AsmAlignP, other_assembly: AsmAlignP, primary_sequence: &str, other_sequence: &str, other_orientation: Orientation, match_indices: [i32; 2]) {
    let primary_match_start = match_indices[0];
    let other_match_start = match_indices[1];
    let primary_length = primary_sequence.len() as i32;
    let other_length = other_assembly.full_sequence_length();
    let mut new_sequence: Option<String> = None;
    let mut primary_offset_adjust = 0;
    let mut other_offset_adjust = 0;
    if other_match_start > primary_match_start {
        primary_offset_adjust = other_match_start - primary_match_start;
        new_sequence = Some(format!("{}{}", &other_sequence[..primary_offset_adjust as usize], primary_sequence));
    } else {
        other_offset_adjust = primary_match_start - other_match_start;
    }
    if primary_match_start + other_length > primary_length {
        let mut s = new_sequence.unwrap_or_else(|| primary_sequence.to_string());
        let extension_length = primary_match_start + other_length - primary_length;
        let other_seq_index = other_length - extension_length;
        s.push_str(&other_sequence[other_seq_index as usize..]);
        new_sequence = Some(s);
    }
    if let Some(s) = new_sequence {
        primary_assembly.m().update_sequence_info(s);
    }
    if primary_offset_adjust > 0 {
        for assembly in primary_assembly.assemblies() {
            for read in assembly.support() {
                let new_assembly_index = read.full_assembly_index_start() + primary_offset_adjust;
                let orientation = read.full_assembly_orientation();
                read.m().set_full_assembly_info(new_assembly_index, orientation);
            }
        }
    }
    if other_orientation.is_reverse() || other_offset_adjust > 0 {
        for assembly in other_assembly.assemblies() {
            for read in assembly.support() {
                let read_orientation;
                let mut new_assembly_index;
                if other_orientation.is_reverse() {
                    let orig_assembly_index = read.full_assembly_index_start();
                    new_assembly_index = other_length - orig_assembly_index - read.base_length();
                    read_orientation = read.full_assembly_orientation().opposite();
                } else {
                    new_assembly_index = read.full_assembly_index_start();
                    read_orientation = read.full_assembly_orientation();
                }
                new_assembly_index += other_offset_adjust;
                read.m().set_full_assembly_info(new_assembly_index, read_orientation);
            }
        }
    }
}
