//! AssemblyLinker: tests two assemblies for a split link (overlapping extension sequences), an indel link or
//! a facing link.

use std::collections::HashSet;

use crate::assembly::assembly_utils::{extract_insert_sequence, is_local_assembly_candidate};
use crate::assembly::constants::*;
use crate::assembly::line_utils::try_line_sequence_link;
use crate::assembly::sequence_compare::{compare_sequences, index_of};
use crate::assembly::types::assembly_link::{AssemblyLink, LinkP};
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::junction_sequence::{JunctionSequence, PHASED_ASSEMBLY_MATCH_SEQ_LENGTH};
use crate::assembly::types::{AssemblyOutcome, LinkType, SupportType};

/// String.indexOf(str, fromIndex) on byte sequences
fn index_of_from(haystack: &[u8], needle: &[u8], from: i32) -> i32 {
    let from = from.max(0) as usize;
    if from > haystack.len() {
        return -1;
    }
    index_of(&haystack[from..], needle).map(|i| i + from as i32).unwrap_or(-1)
}

pub fn try_assembly_facing(first: AsmP, second: AsmP, split_links: &[LinkP]) -> Option<LinkP> {
    if split_links.iter().any(|x| x.matches(first, second)) {
        return None;
    }
    let (j1, j2) = (first.junction(), second.junction());
    if j1.chromosome != j2.chromosome {
        return None;
    }
    if j1.orient == j2.orient {
        return None;
    }
    let lower = if j1.position < j2.position { first } else { second };
    let upper = if first == lower { second } else { first };
    if !lower.junction().is_reverse() {
        return None;
    }
    let link_distance = upper.junction().position - lower.junction().position + 1;
    if !(PHASED_ASSEMBLY_MIN_TI..=PHASED_ASSEMBLY_MAX_TI).contains(&link_distance) {
        return None;
    }
    let mut require_shared_read = true;
    if !first.ref_side_soft_clips().is_empty()
        && !second.ref_side_soft_clips().is_empty()
        && ref_side_soft_clip_matches_junction(&lower, upper.junction().position)
        && ref_side_soft_clip_matches_junction(&upper, lower.junction().position)
    {
        require_shared_read = false;
    }
    if require_shared_read {
        if first.indel() || second.indel() {
            let first_indel_coords = first.indel_coords();
            let second_indel_coords = second.indel_coords();
            let mut matched = false;
            for support in first.support() {
                if !support.type_().is_split_support() {
                    continue;
                }
                for second_support in second.support() {
                    if !second_support.type_().is_split_support() {
                        continue;
                    }
                    if !second_support.matches_fragment(support, true) {
                        continue;
                    }
                    // a null pointer exception in Java if neither assembly has indel coords
                    let target = match first_indel_coords {
                        Some(c) => c,
                        None => second_indel_coords.expect("indel coords"),
                    };
                    if support.indel_coords().is_some_and(|c| c.matches(target)) || second_support.indel_coords().is_some_and(|c| c.matches(target)) {
                        matched = true;
                    }
                    if matched {
                        break;
                    }
                }
            }
            if !matched {
                return None;
            }
        } else {
            if lower.ref_base_length() > link_distance || upper.ref_base_length() > link_distance {
                return None;
            }
            let mut matched = false;
            for is_start in [true, false] {
                let split_support = if is_start { first.support() } else { second.support() };
                let other_support = if is_start { second.support() } else { first.support() };
                for support in split_support {
                    if !support.type_().is_split_support() {
                        continue;
                    }
                    for other in other_support {
                        if other.type_() != SupportType::Junction && other.type_() != SupportType::Discordant {
                            continue;
                        }
                        if other.id() != support.id() {
                            continue;
                        }
                        if other.flags() == support.flags() {
                            matched = true;
                        } else if other.orientation() != support.orientation() {
                            if other.orientation().is_forward() && other.alignment_start() < support.alignment_start() {
                                matched = true;
                            } else if support.orientation().is_forward() && support.alignment_start() < other.alignment_start() {
                                matched = true;
                            }
                        }
                        if matched {
                            break;
                        }
                    }
                    if matched {
                        break;
                    }
                }
            }
            if !matched {
                return None;
            }
        }
    }
    let upper_pos = upper.junction().position;
    let lower_pos = lower.junction().position;
    lower.m().trim_ref_base_position(upper_pos);
    upper.m().trim_ref_base_position(lower_pos);
    Some(AssemblyLink::new(lower, upper, LinkType::Facing, "", ""))
}

/// `in_facing` tells whether the assembly is already in the set of facing assemblies
pub fn is_facing_assembly_candidate(assembly: AsmP, in_facing: bool, split_links: &[LinkP]) -> bool {
    if in_facing {
        return false;
    }
    if assembly.outcome() == AssemblyOutcome::LocalIndel {
        return false;
    }
    if assembly.discordant_only() && !split_links.iter().any(|x| x.has_assembly(assembly)) {
        return false;
    }
    true
}

fn ref_side_soft_clip_matches_junction(assembly: &JunctionAssembly, other_junction_position: i32) -> bool {
    if !assembly.ref_side_soft_clips().iter().any(|x| x.has_proximate_match(other_junction_position)) {
        return false;
    }
    (assembly.ref_base_position() - other_junction_position).abs() <= PROXIMATE_REF_SIDE_SOFT_CLIPS
}

pub fn is_assembly_indel_link(assembly1: &JunctionAssembly, assembly2: &JunctionAssembly) -> bool {
    assembly1.indel()
        && assembly2.indel()
        && assembly1.is_forward_junction() != assembly2.is_forward_junction()
        // a null pointer exception in Java if either has no indel coords
        && assembly1.indel_coords().expect("indel coords").matches(assembly2.indel_coords().expect("indel coords"))
}

pub fn try_assembly_indel(assembly1: AsmP, assembly2: AsmP) -> Option<LinkP> {
    if !is_assembly_indel_link(&assembly1, &assembly2) {
        return None;
    }
    let inserted_bases = assembly1.indel_coords().unwrap().inserted_bases().to_string();
    if assembly1.junction().position < assembly2.junction().position {
        Some(AssemblyLink::new(assembly1, assembly2, LinkType::Indel, &inserted_bases, ""))
    } else {
        Some(AssemblyLink::new(assembly2, assembly1, LinkType::Indel, &inserted_bases, ""))
    }
}

/// the test overload tryAssemblyOverlap(assembly1, assembly2)
pub fn try_assembly_overlap_basic(assembly1: AsmP, assembly2: AsmP) -> Option<LinkP> {
    let is_local_link = is_local_assembly_candidate(&assembly1, &assembly2, false, false);
    try_assembly_overlap(assembly1, assembly2, true, is_local_link)
}

pub fn try_assembly_overlap(assembly1: AsmP, assembly2: AsmP, allow_mismatches: bool, is_local_indel: bool) -> Option<LinkP> {
    let first;
    let second;
    let mut first_reversed = false;
    let mut second_reversed = false;
    if assembly1.junction().orient != assembly2.junction().orient {
        first = if assembly1.junction().is_forward() { assembly1 } else { assembly2 };
        second = if assembly1.junction().is_reverse() { assembly1 } else { assembly2 };
    } else {
        first = assembly1;
        second = assembly2;
        if assembly1.junction().is_forward() {
            second_reversed = true;
        } else {
            first_reversed = true;
        }
    }
    if let Some(line_link) = try_line_sequence_link(first, second, first_reversed, second_reversed) {
        return Some(line_link);
    }
    let require_ref_base_overlap = is_local_indel;
    let mut first_seq;
    let second_seq;
    if require_ref_base_overlap {
        first_seq = JunctionSequence::form_straddling_match_sequence(&first);
        second_seq = JunctionSequence::form_straddling_match_sequence(&second);
    } else {
        first_seq = JunctionSequence::form_outer_extension_match_sequence(&first, first_reversed);
        second_seq = JunctionSequence::form_outer_extension_match_sequence(&second, second_reversed);
    }
    let second_full = second_seq.full_sequence.as_bytes();
    if let Some(first_seq_index_in_second) = index_of(second_full, first_seq.match_sequence().as_bytes()) {
        return form_link(first, second, &first_seq, &second_seq, first_seq.match_seq_start_index(), first_seq_index_in_second, 0);
    }
    if !allow_mismatches {
        return None;
    }
    if first.extension_length() > 2 * PHASED_ASSEMBLY_MATCH_SEQ_LENGTH {
        first_seq = JunctionSequence::form_full_extension_match_sequence(&first, first_reversed);
    }
    let first_match_sequence = first_seq.match_sequence().as_bytes();
    let first_match_seq_length = first_match_sequence.len() as i32;

    let mut alternative_index_starts: Vec<[i32; 2]> = Vec::new();
    let mut tested_sequences: HashSet<&[u8]> = HashSet::new();
    let sub_seq_iterations = first_match_seq_length / MATCH_SUBSEQUENCE_LENGTH;
    for i in 0..sub_seq_iterations {
        let match_seq_start_index = i * MATCH_SUBSEQUENCE_LENGTH;
        let match_seq_end_index = match_seq_start_index + MATCH_SUBSEQUENCE_LENGTH;
        if match_seq_end_index >= first_match_seq_length + 1 {
            break;
        }
        let first_sub_sequence = &first_match_sequence[match_seq_start_index as usize..match_seq_end_index as usize];
        if !tested_sequences.insert(first_sub_sequence) {
            continue;
        }
        let mut second_sub_seq_index = index_of_from(second_full, first_sub_sequence, 0);
        if second_sub_seq_index < 0 {
            continue;
        }
        let mut implied_second_seq_match_seq_start = second_sub_seq_index - match_seq_start_index;
        let implied_second_seq_match_seq_end = implied_second_seq_match_seq_start + first_match_seq_length - 1;
        if implied_second_seq_match_seq_end >= second_seq.base_length {
            let second_match_length = second_seq.base_length - implied_second_seq_match_seq_start - 1;
            if second_match_length < ASSEMBLY_LINK_OVERLAP_BASES {
                continue;
            }
        }
        alternative_index_starts.push([match_seq_start_index, second_sub_seq_index]);
        second_sub_seq_index = index_of_from(second_full, first_sub_sequence, second_sub_seq_index + MATCH_SUBSEQUENCE_LENGTH);
        while second_sub_seq_index >= 0 {
            implied_second_seq_match_seq_start = second_sub_seq_index - match_seq_start_index;
            if implied_second_seq_match_seq_start + first_match_seq_length > second_seq.base_length {
                break;
            }
            alternative_index_starts.push([match_seq_start_index, second_sub_seq_index]);
            second_sub_seq_index = index_of_from(second_full, first_sub_sequence, second_sub_seq_index + MATCH_SUBSEQUENCE_LENGTH);
        }
    }
    if alternative_index_starts.is_empty() {
        return None;
    }
    let mut min_overlap_length = first.extension_length().min(second.extension_length()).min(ASSEMBLY_LINK_OVERLAP_BASES);
    if first.discordant_only() && second.discordant_only() {
        min_overlap_length = min_overlap_length.min(ASSEMBLY_LINK_DISC_ONLY_OVERLAP_BASES);
    }
    let top_match_indices = if require_ref_base_overlap {
        find_local_sequence_match(&first_seq, &second_seq, min_overlap_length, &alternative_index_starts)
    } else {
        find_best_sequence_match(&first_seq, &second_seq, min_overlap_length, &alternative_index_starts)
    };
    match top_match_indices {
        Some(t) => form_link(first, second, &first_seq, &second_seq, t[0], t[1], t[2]),
        None => None,
    }
}

pub fn find_local_sequence_match(first_seq: &JunctionSequence, second_seq: &JunctionSequence, min_overlap_length: i32, alternative_index_starts: &[[i32; 2]]) -> Option<[i32; 3]> {
    if alternative_index_starts.is_empty() {
        return None;
    }
    let mut top_match_length = 0;
    let mut top_match_mismatch_penalty = 0.0;
    let mut top_match_indices: Option<[i32; 3]> = None;
    let mut tested_offsets: HashSet<i32> = HashSet::new();
    for index_starts in alternative_index_starts {
        let first_match_seq_match_index = index_starts[0];
        let second_match_index = index_starts[1];
        let match_offset = second_match_index - first_match_seq_match_index;
        if !tested_offsets.insert(match_offset) {
            continue;
        }
        let first_match_index = first_match_seq_match_index + first_seq.match_seq_start_index();
        let min_lower_extension = second_match_index.min(first_match_index);
        let max_lower_extension = (first_seq.base_length - first_match_index - 1).min(second_seq.base_length - second_match_index - 1);
        let first_index_start = first_match_index - min_lower_extension;
        let first_index_end = first_match_index + max_lower_extension;
        let second_index_start = second_match_index - min_lower_extension;
        let second_index_end = second_match_index + max_lower_extension;
        let overlap_length = first_index_end - first_index_start + 1;
        if overlap_length < min_overlap_length {
            continue;
        }
        let mismatch_penalty = compare_sequences(
            first_seq.bases(), first_seq.base_quals(), first_index_start, first_index_end, first_seq.repeat_info(),
            second_seq.bases(), second_seq.base_quals(), second_index_start, second_index_end, second_seq.repeat_info(), PRIMARY_ASSEMBLY_MERGE_MISMATCH as f64,
        );
        if mismatch_penalty > PRIMARY_ASSEMBLY_MERGE_MISMATCH as f64 {
            continue;
        }
        if overlap_length > top_match_length || (overlap_length == top_match_length && mismatch_penalty < top_match_mismatch_penalty) {
            let min_ref_overlap_length = MIN_VARIANT_LENGTH / 2;
            let first_ref_base_overlap = first_seq.junction_index() - first_index_start + 1;
            let first_ext_base_overlap = first_index_end - first_seq.junction_index();
            let second_ref_base_overlap = second_index_end - second_seq.junction_index() + 1;
            let second_ext_base_overlap = second_seq.junction_index() - second_index_start;
            if first_ref_base_overlap < min_ref_overlap_length && second_ref_base_overlap < min_ref_overlap_length {
                continue;
            }
            if first_ext_base_overlap < min_ref_overlap_length && second_ext_base_overlap < min_ref_overlap_length {
                continue;
            }
        }
        // as in Java, taken as the top match whether or not it beat the previous one
        top_match_length = overlap_length;
        top_match_indices = Some([first_index_start, second_index_start, 0]);
        top_match_mismatch_penalty = mismatch_penalty;
    }
    top_match_indices
}

pub fn find_best_sequence_match(first_seq: &JunctionSequence, second_seq: &JunctionSequence, min_overlap_length: i32, alternative_index_starts: &[[i32; 2]]) -> Option<[i32; 3]> {
    if alternative_index_starts.is_empty() {
        return None;
    }
    let mut top_match_length = 0;
    let mut top_match_mismatch_penalty = 0.0;
    let mut top_match_indices: Option<[i32; 3]> = None;
    let mut tested_offsets: HashSet<i32> = HashSet::new();
    let first_junction_seq_length = (first_seq.match_sequence().len() as i32).min(PHASED_ASSEMBLY_MATCH_SEQ_LENGTH);
    for index_starts in alternative_index_starts {
        let first_match_seq_match_index = index_starts[0];
        let second_match_index = index_starts[1];
        let match_offset = second_match_index - first_match_seq_match_index;
        if !tested_offsets.insert(match_offset) {
            continue;
        }
        let mut second_index_start = second_match_index - first_match_seq_match_index;
        let mut second_index_end = second_index_start + first_junction_seq_length - 1;
        let mut first_match_index_start = 0;
        let mut first_match_index_end = first_junction_seq_length - 1;
        if second_index_start < 0 {
            first_match_index_start += -second_index_start;
            first_match_index_end += -second_index_start;
            second_index_end += -second_index_start;
            second_index_start = 0;
        }
        if second_index_end >= second_seq.base_length {
            let end_reduction = second_index_end - second_seq.base_length + 1;
            second_index_end -= end_reduction;
            first_match_index_end -= end_reduction;
        }
        let mut first_index_start = first_match_index_start + first_seq.match_seq_start_index();
        let first_index_end = (first_match_index_end + first_seq.match_seq_start_index()).min(first_seq.base_length - 1);
        let mut overlap_length = (first_index_end - first_index_start + 1).min(second_index_end - second_index_start + 1);
        if overlap_length < min_overlap_length {
            continue;
        }
        let mismatch_penalty = compare_sequences(
            first_seq.bases(), first_seq.base_quals(), first_index_start, first_index_end, first_seq.repeat_info(),
            second_seq.bases(), second_seq.base_quals(), second_index_start, second_index_end, second_seq.repeat_info(), PRIMARY_ASSEMBLY_MERGE_MISMATCH as f64,
        );
        if mismatch_penalty > PRIMARY_ASSEMBLY_MERGE_MISMATCH as f64 {
            continue;
        }
        let match_extensions = try_extend_matched_sequence(first_seq.bases(), first_index_start, first_index_end, second_seq.bases(), second_index_start, second_index_end);
        first_index_start -= match_extensions[0];
        second_index_start -= match_extensions[0];
        overlap_length += match_extensions[0] + match_extensions[1];
        if overlap_length > top_match_length || (overlap_length == top_match_length && mismatch_penalty < top_match_mismatch_penalty) {
            top_match_length = overlap_length;
            top_match_indices = Some([first_index_start, second_index_start, 0]);
            top_match_mismatch_penalty = mismatch_penalty;
        }
    }
    top_match_indices
}

fn try_extend_matched_sequence(first_bases: &[u8], first_index_start: i32, first_index_end: i32, second_bases: &[u8], second_index_start: i32, second_index_end: i32) -> [i32; 2] {
    let mut earlier_matched_bases = 0;
    let mut latter_matched_bases = 0;
    let mut first_index = first_index_start - 1;
    let mut second_index = second_index_start - 1;
    while first_index >= 0 && second_index >= 0 {
        if first_bases[first_index as usize] != second_bases[second_index as usize] {
            break;
        }
        earlier_matched_bases += 1;
        first_index -= 1;
        second_index -= 1;
    }
    first_index = first_index_end + 1;
    second_index = second_index_end + 1;
    while first_index < first_bases.len() as i32 && second_index < second_bases.len() as i32 {
        if first_bases[first_index as usize] != second_bases[second_index as usize] {
            break;
        }
        latter_matched_bases += 1;
        first_index += 1;
        second_index += 1;
    }
    [earlier_matched_bases, latter_matched_bases]
}

pub fn form_link(first: AsmP, second: AsmP, first_seq: &JunctionSequence, second_seq: &JunctionSequence, first_index_start: i32, second_index_start: i32, first_mismatch_diff: i32) -> Option<LinkP> {
    let first_junction_offset = first_seq.junction_index() - first_index_start - first_mismatch_diff;
    let second_junction_offset = second_seq.junction_index() - second_index_start;
    let mut first_junction_index_in_second = second_index_start + first_junction_offset;
    let junction_offset_diff;
    if !first_seq.reversed && !second_seq.reversed {
        junction_offset_diff = second_junction_offset - first_junction_offset - 1;
    } else if second_seq.reversed {
        first_junction_index_in_second = second_seq.index_reverted(first_junction_index_in_second);
        junction_offset_diff = first_junction_index_in_second - second.junction_index() - 1;
    } else {
        junction_offset_diff = second.junction_index() - first_junction_index_in_second - 1;
    }
    let mut inserted_bases = String::new();
    let mut overlap_bases = String::new();
    if junction_offset_diff > 0 {
        let second_reversed = first.is_forward_junction() == second.is_forward_junction();
        inserted_bases = extract_insert_sequence(&first, false, &second, second_reversed, junction_offset_diff);
    } else if junction_offset_diff < 0 {
        let overlap_length = junction_offset_diff.abs();
        if overlap_length >= first.ref_base_length() || overlap_length >= second.ref_base_length() {
            return None;
        }
        overlap_bases = match extract_overlap_bases(&first, overlap_length) {
            Some(b) => b,
            None => match extract_overlap_bases(&second, junction_offset_diff) {
                Some(b) => b,
                None => {
                    crate::debug!(
                        "asm({} & {}) invalid insert/overlap junctOffsetDiff({} firstIndex={} secIndex={}) on firstSeq({}) secSeq({})",
                        first.junction().coords(), second.junction().coords(), junction_offset_diff, first_index_start, second_index_start, first_seq, second_seq
                    );
                    return None;
                }
            },
        };
    }
    Some(AssemblyLink::new(first, second, LinkType::Split, &inserted_bases, &overlap_bases))
}

fn extract_overlap_bases(assembly: &JunctionAssembly, overlap_length: i32) -> Option<String> {
    let extra_bases_start_index;
    let extra_bases_end_index;
    if assembly.is_forward_junction() {
        extra_bases_start_index = assembly.junction_index() - overlap_length + 1;
        extra_bases_end_index = assembly.junction_index();
    } else {
        extra_bases_start_index = assembly.junction_index();
        extra_bases_end_index = assembly.junction_index() + overlap_length - 1;
    }
    if extra_bases_start_index >= 0 && extra_bases_end_index >= extra_bases_start_index && extra_bases_end_index < assembly.base_length() {
        return Some(assembly.form_sequence(extra_bases_start_index, extra_bases_end_index));
    }
    None
}
