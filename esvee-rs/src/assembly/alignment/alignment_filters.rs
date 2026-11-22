//! AlignmentFilters: decides which alignments of an assembly sequence are used to form breakends.

use std::cmp::Ordering;

use crate::common::compat;
use crate::assembly::alignment::align_data::AlignP;
use crate::assembly::alignment::alternative_alignment::AlternativeAlignment;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::alignment::breakend_builder::{segment_orientation, segment_orientation_of};
use crate::assembly::alignment::md_tag::MdTag;
use crate::assembly::common_utils::positions_within;
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::sam::{calc_cigar_aligned_length, CigarOp};
use crate::common::java_sort::java_sort;

const SSX2_REGIONS_V37: [(&str, i32, i32); 2] = [("X", 52729628, 52731680), ("X", 52784877, 52786929)];
const SSX2_REGIONS_V38: [(&str, i32, i32); 2] = [("chrX", 52700578, 52702630), ("chrX", 52755800, 52757852)];
const SSX2_GENE_ORIENT: Orientation = Orientation::Reverse;
const SSX2_MAX_MAP_QUAL: i32 = 20;
// caller FilterConstants
const PMS2_V37: (&str, i32, i32) = ("7", 6002870, 6058756);
const PMS2_V38: (&str, i32, i32) = ("chr7", 5960925, 6019106);
const PMS2CL_V37: (&str, i32, i32) = ("7", 6759759, 6803493);
const PMS2CL_V38: (&str, i32, i32) = ("chr7", 6725305, 6761392);

pub fn filter_alignments(assembly_alignment: AsmAlignP, alignments: &[AlignP], valid_alignments: &mut Vec<AlignP>, low_qual_alignments: &mut Vec<AlignP>) {
    let full_sequence = assembly_alignment.full_sequence().expect("full sequence").to_string();
    for x in alignments {
        x.m().set_full_sequence_data(&full_sequence, assembly_alignment.full_sequence_length());
    }
    let mut candidate_alignments: Vec<AlignP> = Vec::new();
    for alignment in alignments {
        if exceeds_initial_map_qual_threshold(*alignment) || alignment.has_alt_alignments() {
            candidate_alignments.push(*alignment);
        } else {
            low_qual_alignments.push(*alignment);
        }
    }
    // AlignmentOrderComparator. Java: never reports equality, so the order of equal starts follows the sort's
    // steps. Corrected: equal starts compare equal and keep their order (the sort is stable)
    let equal_starts = if compat::java("A25") { Ordering::Less } else { Ordering::Equal };
    java_sort(&mut candidate_alignments, |first, second| {
        if first.sequence_start() == second.sequence_start() {
            return equal_starts;
        }
        if first.sequence_start() < second.sequence_start() { Ordering::Less } else { Ordering::Greater }
    });
    for i in 0..candidate_alignments.len() {
        let alignment = candidate_alignments[i];
        let mut overlap_start = 0;
        let mut overlap_end = 0;
        if i > 0 {
            let prev_alignment = candidate_alignments[i - 1];
            if prev_alignment.sequence_end() >= alignment.sequence_start() {
                overlap_start = prev_alignment.sequence_end() - alignment.sequence_start() + 1;
            }
        }
        if i + 1 < candidate_alignments.len() {
            let next_alignment = candidate_alignments[i + 1];
            if alignment.sequence_end() >= next_alignment.sequence_start() {
                overlap_end = alignment.sequence_end() - next_alignment.sequence_start() + 1;
            }
        }
        alignment.m().set_adjusted_alignment(&full_sequence, overlap_start, overlap_end);
    }
    if assembly_alignment.is_simple_indel_assembly() && candidate_alignments.len() == 2 {
        valid_alignments.extend(candidate_alignments);
        return;
    }
    remove_short_adjusted_alignments(&mut candidate_alignments, low_qual_alignments);
    remove_low_modified_map_quals(&mut candidate_alignments, low_qual_alignments);
    check_specific_low_map_pairings(&candidate_alignments);
    let valid_count = candidate_alignments.iter().filter(|x| exceeds_mod_map_qual_threshold(**x)).count();
    if candidate_alignments.len() == valid_count {
        valid_alignments.extend(candidate_alignments);
        return;
    }
    check_local_variants(&candidate_alignments, valid_alignments, low_qual_alignments);
}

fn exceeds_initial_map_qual_threshold(alignment: AlignP) -> bool {
    exceeds_map_qual_threshold(alignment, false, ALIGNMENT_MIN_MOD_MAP_QUAL_NO_XA)
}
fn exceeds_mod_map_qual_threshold(alignment: AlignP) -> bool {
    exceeds_map_qual_threshold(alignment, true, ALIGNMENT_MIN_MOD_MAP_QUAL)
}
fn exceeds_no_alts_mod_map_qual_threshold(alignment: AlignP) -> bool {
    exceeds_map_qual_threshold(alignment, true, ALIGNMENT_MIN_MOD_MAP_QUAL_NO_XA)
}
fn exceeds_map_qual_threshold(alignment: AlignP, use_modified: bool, threshold: i32) -> bool {
    let map_qual = if use_modified { alignment.modified_map_qual() } else { alignment.map_qual() as f64 };
    map_qual >= threshold as f64
}

fn remove_short_adjusted_alignments(candidate_alignments: &mut Vec<AlignP>, low_qual_alignments: &mut Vec<AlignP>) {
    let mut index = 0;
    while index < candidate_alignments.len() {
        let alignment = candidate_alignments[index];
        let mut keep_alignment = alignment.adjusted_alignment() >= ALIGNMENT_MIN_ADJUST_ALIGN_LENGTH;
        if !keep_alignment {
            if index > 0 {
                let prev_alignment = candidate_alignments[index - 1];
                if prev_alignment.adjusted_alignment() >= ALIGNMENT_MIN_ADJUST_ALIGN_LENGTH && are_close_alignments(prev_alignment, alignment) {
                    keep_alignment = true;
                }
            }
            if !keep_alignment && index + 1 < candidate_alignments.len() {
                let next_alignment = candidate_alignments[index + 1];
                if next_alignment.adjusted_alignment() >= ALIGNMENT_MIN_ADJUST_ALIGN_LENGTH && are_close_alignments(next_alignment, alignment) {
                    keep_alignment = true;
                }
            }
            if keep_alignment {
                let md_tag = MdTag::new(alignment.md_tag());
                let insert_count = alignment.cigar_elements().iter().filter(|x| x.op == CigarOp::I).count() as i32;
                if md_tag.mismatch_count() + insert_count > ALIGNMENT_RECOVERY_MAX_MD_ERRORS {
                    keep_alignment = false;
                }
            }
        }
        if !keep_alignment {
            candidate_alignments.remove(index);
            low_qual_alignments.push(alignment);
        } else {
            index += 1;
        }
    }
}

fn are_close_alignments(first: AlignP, second: AlignP) -> bool {
    if first.chromosome() != second.chromosome() {
        return false;
    }
    let min_distance = (first.position_start() - second.position_end()).abs().min((first.position_end() - second.position_start()).abs());
    min_distance <= ALIGNMENT_PROXIMATE_DISTANCE
}

fn remove_low_modified_map_quals(candidate_alignments: &mut Vec<AlignP>, low_qual_alignments: &mut Vec<AlignP>) {
    let mut index = 0;
    while index < candidate_alignments.len() {
        let alignment = candidate_alignments[index];
        if !exceeds_no_alts_mod_map_qual_threshold(alignment) && !alignment.has_alt_alignments() {
            candidate_alignments.remove(index);
            low_qual_alignments.push(alignment);
        } else {
            index += 1;
        }
    }
}

fn check_local_variants(candidate_alignments: &[AlignP], valid_alignments: &mut Vec<AlignP>, low_qual_alignments: &mut Vec<AlignP>) {
    if candidate_alignments.len() < 2 {
        for alignment in candidate_alignments {
            if exceeds_mod_map_qual_threshold(*alignment) {
                valid_alignments.push(*alignment);
            } else {
                low_qual_alignments.push(*alignment);
            }
        }
        return;
    }
    for i in 0..candidate_alignments.len() - 1 {
        let alignment = candidate_alignments[i];
        let next_alignment = candidate_alignments[i + 1];
        if exceeds_mod_map_qual_threshold(alignment) == exceeds_mod_map_qual_threshold(next_alignment) {
            continue;
        }
        find_shortest_local_variant(alignment, next_alignment);
    }
    for i in 0..candidate_alignments.len() - 1 {
        let alignment = candidate_alignments[i];
        let next_alignment = candidate_alignments[i + 1];
        if exceeds_mod_map_qual_threshold(alignment) || exceeds_mod_map_qual_threshold(next_alignment) {
            continue;
        }
        if alignment.has_low_map_qual_alignment() && next_alignment.has_low_map_qual_alignment() {
            continue;
        }
        find_shortest_local_variant(alignment, next_alignment);
    }
    for alignment in candidate_alignments {
        if exceeds_mod_map_qual_threshold(*alignment) || alignment.has_low_map_qual_alignment() {
            valid_alignments.push(*alignment);
        } else {
            low_qual_alignments.push(*alignment);
        }
    }
}

fn collect_alignment_candidates(alignment: AlignP, links_at_end: bool) -> Vec<AlternativeAlignment> {
    if let Some(initial) = alignment.selected_alt_alignment() {
        return vec![initial.clone()];
    }
    let orientation = segment_orientation(&alignment, links_at_end);
    let use_end_position = alignment.is_forward() == links_at_end;
    let position = if use_end_position { alignment.position_end() } else { alignment.position_start() };
    let initial_alignment = AlternativeAlignment::new(alignment.chromosome(), position, orientation, "", alignment.map_qual());
    if exceeds_mod_map_qual_threshold(alignment) {
        return vec![initial_alignment];
    }
    let mut alignments = vec![initial_alignment];
    for alt_alignment in alignment.m().raw_alt_alignments() {
        let alt_orientation = segment_orientation_of(alt_alignment.orient, links_at_end);
        let use_end_position = alt_alignment.orient.is_forward() == links_at_end;
        let mut position = alt_alignment.position;
        if use_end_position {
            position += calc_cigar_aligned_length(&alt_alignment.cigar) - 1;
        }
        alignments.push(AlternativeAlignment::new(&alt_alignment.chromosome, position, alt_orientation, &alt_alignment.cigar, alt_alignment.num_mutations));
    }
    alignments
}

fn find_shortest_local_variant(first: AlignP, second: AlignP) -> bool {
    let first_alignments = collect_alignment_candidates(first, true);
    let second_alignments = collect_alignment_candidates(second, false);
    if first_alignments.is_empty() || second_alignments.is_empty() {
        return false;
    }
    let mut shortest_length = -1;
    // indices of the selected alternatives in the two lists; index zero is the initial alignment
    let mut first_selected_alt: Option<usize> = None;
    let mut second_selected_alt: Option<usize> = None;
    for (fi, first_alt) in first_alignments.iter().enumerate() {
        for (si, second_alt) in second_alignments.iter().enumerate() {
            if first_alt.chromosome != second_alt.chromosome {
                continue;
            }
            let sv_length = (first_alt.position - second_alt.position).abs();
            if sv_length > ALIGNMENT_LOW_MOD_MQ_VARIANT_LENGTH {
                continue;
            }
            if shortest_length < 0 || sv_length < shortest_length {
                first_selected_alt = Some(fi);
                second_selected_alt = Some(si);
                shortest_length = sv_length;
            }
        }
    }
    if first_selected_alt.is_none() && (exceeds_no_alts_mod_map_qual_threshold(first) || (first.has_alt_alignments() && exceeds_map_qual_threshold(first, true, 1))) {
        first_selected_alt = Some(0);
    }
    if second_selected_alt.is_none() && (exceeds_no_alts_mod_map_qual_threshold(second) || (second.has_alt_alignments() && exceeds_map_qual_threshold(second, true, 1))) {
        second_selected_alt = Some(0);
    }
    let has_short_sv_link = shortest_length >= 0;
    mark_alt_alignment(first, first_selected_alt, &first_alignments, has_short_sv_link);
    mark_alt_alignment(second, second_selected_alt, &second_alignments, has_short_sv_link);
    first_selected_alt.is_some() && second_selected_alt.is_some()
}

fn mark_alt_alignment(alignment: AlignP, selected: Option<usize>, all_alignments: &[AlternativeAlignment], has_short_sv_link: bool) {
    let Some(selected) = selected else { return };
    if exceeds_mod_map_qual_threshold(alignment) {
        return;
    }
    let unselected: Vec<AlternativeAlignment> = all_alignments.iter().enumerate().filter(|(i, _)| *i != selected).map(|(_, x)| x.clone()).collect();
    alignment.m().set_selected_alt_alignments(if selected != 0 { Some(all_alignments[selected].clone()) } else { None }, Some(unselected), has_short_sv_link);
}

fn check_specific_low_map_pairings(alignments: &[AlignP]) {
    check_ssx2_low_map_pairings(alignments);
    check_pms2_low_map_pairings(alignments);
}

fn matches_ssx2_region(ssx2_regions: &[(&str, i32, i32); 2], chromosome: &str, position: i32, orientation: Orientation, map_qual: i32) -> bool {
    if map_qual > SSX2_MAX_MAP_QUAL {
        return false;
    }
    if orientation != SSX2_GENE_ORIENT {
        return false;
    }
    ssx2_regions.iter().any(|(c, s, e)| *c == chromosome && position >= *s && position <= *e)
}

fn check_ssx2_low_map_pairings(alignments: &[AlignP]) {
    if alignments.len() < 2 {
        return;
    }
    let mut low_mapped_alignment: Option<AlignP> = None;
    let ssx2_regions = if alignments[0].chromosome().starts_with("chr") { &SSX2_REGIONS_V38 } else { &SSX2_REGIONS_V37 };
    let ssx2_region = ssx2_regions[0];
    for alignment in alignments {
        if matches_ssx2_region(ssx2_regions, alignment.chromosome(), alignment.position_start(), alignment.orientation(), alignment.map_qual()) {
            low_mapped_alignment = Some(*alignment);
        } else {
            for alt_alignment in alignment.m().raw_alt_alignments() {
                if matches_ssx2_region(ssx2_regions, &alt_alignment.chromosome, alt_alignment.position, alt_alignment.orient, alignment.map_qual()) {
                    low_mapped_alignment = Some(*alignment);
                    break;
                }
            }
        }
        if let Some(ssx2_alignment) = low_mapped_alignment {
            if alignments.iter().filter(|x| **x != ssx2_alignment).any(|x| x.chromosome() != ssx2_region.0) {
                let alt_alignment = AlternativeAlignment::new(ssx2_region.0, ssx2_region.1, SSX2_GENE_ORIENT.opposite(), ssx2_alignment.cigar(), 0);
                ssx2_alignment.m().set_selected_low_map_qual_alt_alignment(Some(alt_alignment));
                return;
            }
        }
    }
}

fn matches_pms2_region(region: (&str, i32, i32), chromosome: &str, position_start: i32, position_end: i32) -> bool {
    region.0 == chromosome && positions_within(position_start, position_end, region.1, region.2)
}

fn check_pms2_low_map_pairings(alignments: &[AlignP]) {
    if alignments.is_empty() {
        return;
    }
    let v38 = alignments[0].chromosome().starts_with("chr");
    let pms2 = if v38 { PMS2_V38 } else { PMS2_V37 };
    let pms2cl = if v38 { PMS2CL_V38 } else { PMS2CL_V37 };
    for alignment in alignments {
        if alignment.map_qual() >= ALIGNMENT_MIN_MOD_MAP_QUAL {
            continue;
        }
        if matches_pms2_region(pms2, alignment.chromosome(), alignment.position_start(), alignment.position_end()) {
            alignment.m().set_modified_map_qual(ALIGNMENT_MIN_MOD_MAP_QUAL);
        } else if matches_pms2_region(pms2cl, alignment.chromosome(), alignment.position_start(), alignment.position_end()) {
            for alt_alignment in alignment.m().raw_alt_alignments() {
                if matches_pms2_region(pms2, &alt_alignment.chromosome, alt_alignment.position, alt_alignment.position_end()) {
                    if alignment.map_qual() >= ALIGNMENT_MIN_MOD_MAP_QUAL {
                        alignment.m().set_selected_low_map_qual_alt_alignment(Some(alt_alignment.clone()));
                    } else {
                        let min_qual_alt_alignment = AlternativeAlignment::new(&alt_alignment.chromosome, alt_alignment.position, alt_alignment.orient, &alt_alignment.cigar, ALIGNMENT_MIN_MOD_MAP_QUAL);
                        alignment.m().set_selected_low_map_qual_alt_alignment(Some(min_qual_alt_alignment));
                    }
                    break;
                }
            }
        }
    }
}
