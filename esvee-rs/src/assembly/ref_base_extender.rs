//! RefBaseExtender: extends an assembly's reference bases with junction mates and discordant reads, and
//! branches an assembly where a ref-side soft clip has enough support.

use std::collections::HashSet;

use crate::assembly::assembly_utils::{above_low_base_qual_threshold, bases_match_qual};
use crate::assembly::common_utils::{above_min_qual, is_duplication_fragment};
use crate::assembly::constants::*;
use crate::assembly::globals::has_paired_reads;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_utils::{is_discordant_fragment, is_valid_support_coords_vs_junction};
use crate::assembly::ref_genome::RefGenomeInterface;
use crate::assembly::remote_region_finder::find_remote_regions;
use crate::assembly::sam::inferred_insert_size_abs;
use crate::assembly::sequence_compare::index_of;
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::phase_group::PhaseGroup;
use crate::assembly::types::read_assembly_indices::ReadAssemblyIndices;
use crate::assembly::types::ref_side_soft_clip::RefSideSoftClip;
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::{AssemblyOutcome, SupportType};
use crate::common::java_sort::java_sort;

const REF_READ_SEARCH_LENGTH: i32 = 20;

pub fn find_assembly_candidate_extensions(assembly: AsmP, unfiltered_non_junction_reads: &[ReadP]) {
    let mut new_ref_base_position = assembly.ref_base_position();
    let is_forward_junction = assembly.junction().is_forward();
    let junction_position = assembly.junction().position;
    let is_indel_junction = assembly.indel();
    // (read, type) as the Java NonJunctionRead
    let mut candidate_reads: Vec<(ReadP, SupportType)> = Vec::new();
    let discordant_reads: Vec<ReadP>;
    if is_indel_junction || !has_paired_reads() {
        discordant_reads = Vec::new();
    } else {
        let concordant_read_ids: HashSet<String> = unfiltered_non_junction_reads.iter().filter(|x| is_concordant_read(x)).map(|x| x.id().to_string()).collect();
        discordant_reads = unfiltered_non_junction_reads
            .iter()
            .filter(|x| is_discordant_candidate(x, is_forward_junction, junction_position, &assembly, &concordant_read_ids))
            .cloned()
            .collect();
        for x in &discordant_reads {
            candidate_reads.push((x.clone(), SupportType::Discordant));
        }
        for x in &discordant_reads {
            if x.is_mate_unmapped() {
                if let Some(mate) = x.mate_read().cloned() {
                    if !filter_unmapped(&mate, true) {
                        assembly.m().add_unmapped_read(mate);
                    }
                }
            }
        }
    }
    let mut remote_junction_mates: Vec<ReadP> = Vec::new();
    let mut supp_junction_reads: Vec<ReadP> = Vec::new();
    for read in assembly.support().clone() {
        let cached = read.cached_read().expect("support read without cached read").clone();
        if cached.has_supplementary() {
            supp_junction_reads.push(cached.clone());
        }
        if !is_indel_junction && read.is_discordant() {
            remote_junction_mates.push(cached);
            continue;
        }
        let Some(mate_read) = cached.mate_read().cloned() else { continue };
        mate_read.m().mark_junction_mate();
        if discordant_reads.contains(&mate_read) {
            continue;
        }
        if !mate_read.is_unmapped() {
            let is_past_junction = (is_forward_junction && mate_read.alignment_end() >= junction_position) || (!is_forward_junction && mate_read.alignment_start() <= junction_position);
            if is_past_junction {
                if !mate_read.is_unmapped() && !mate_read.is_left_clipped() && !mate_read.is_right_clipped() {
                    assembly.m().add_concordant_candidate(mate_read);
                }
            } else {
                candidate_reads.push((mate_read, SupportType::JunctionMate));
            }
        } else {
            if filter_unmapped(&mate_read, false) {
                continue;
            }
            assembly.m().add_unmapped_read(mate_read);
        }
    }
    if !candidate_reads.is_empty() {
        for (read, _) in &candidate_reads {
            assembly.m().add_candidate_support(read.clone());
        }
        let mut sorted_candidate_reads = candidate_reads.clone();
        java_sort(&mut sorted_candidate_reads, |a, b| {
            let ka = if is_forward_junction { -a.0.alignment_end() } else { a.0.alignment_start() };
            let kb = if is_forward_junction { -b.0.alignment_end() } else { b.0.alignment_start() };
            ka.cmp(&kb)
        });
        for (read, _) in &sorted_candidate_reads {
            if is_forward_junction {
                if read.has_junction_mate() || read.alignment_end() >= new_ref_base_position - ASSEMBLY_REF_BASE_MAX_GAP {
                    new_ref_base_position = new_ref_base_position.min(read.alignment_start());
                }
            } else if read.has_junction_mate() || read.alignment_start() <= new_ref_base_position + ASSEMBLY_REF_BASE_MAX_GAP {
                new_ref_base_position = new_ref_base_position.max(read.alignment_end());
            }
        }
    }
    if !is_indel_junction {
        find_remote_regions(assembly, &discordant_reads, &remote_junction_mates, &supp_junction_reads);
    }
    RefSideSoftClip::purge_ref_side_soft_clips(assembly.m().ref_side_soft_clips_mut(), new_ref_base_position);
}

fn is_concordant_read(read: &Read) -> bool {
    if !read.is_paired_read() || !read.is_mate_mapped() {
        return false;
    }
    if Some(read.chromosome()) != read.mate_chromosome() || read.orientation() == read.mate_orientation() {
        return false;
    }
    let fragment_size = inferred_insert_size_abs(read.bam_record());
    if fragment_size > DEFAULT_MAX_CONCORDANT_FRAG_LENGTH {
        return false;
    }
    !is_duplication_fragment(read.bam_record(), fragment_size)
}

fn is_discordant_candidate(read: &ReadP, is_forward_junction: bool, junction_position: i32, assembly: &JunctionAssembly, concordant_read_ids: &HashSet<String>) -> bool {
    if !is_valid_support_coords_vs_junction(read, is_forward_junction, junction_position) {
        return false;
    }
    if !is_discordant_fragment(read) {
        if !read.is_mate_unmapped() {
            return false;
        }
        if let Some(mate) = read.mate_read().cloned() {
            if filter_unmapped(&mate, true) {
                return false;
            }
        }
    } else if concordant_read_ids.contains(read.id()) {
        return false;
    }
    !assembly.has_read_support(read.mate_read())
}

fn filter_unmapped(read: &ReadP, is_discordant: bool) -> bool {
    if !read.is_unmapped() {
        return false;
    }
    Read::trim_low_qual_bases(read);
    if is_discordant && read.mate_read().is_some_and(|m| m.mapping_quality() == 0) {
        return true;
    }
    read.bases_length() - read.trim_count() < UNMAPPED_TRIM_THRESHOLD
}

pub fn extend_ref_bases(assembly: AsmP, candidate_support: &mut Vec<ReadP>, ref_genome: Option<&dyn RefGenomeInterface>, allow_branching: bool) {
    if candidate_support.is_empty() {
        return;
    }
    let is_forward_junction = assembly.junction().is_forward();
    let initial_ref_position = assembly.ref_base_position();
    let mut new_ref_base_position = initial_ref_position;
    java_sort(candidate_support, |a, b| {
        let ka = if is_forward_junction { -a.alignment_end() } else { a.alignment_start() };
        let kb = if is_forward_junction { -b.alignment_end() } else { b.alignment_start() };
        ka.cmp(&kb)
    });
    let junction = assembly.junction();
    let mut non_junction_support: Vec<ReadP> = Vec::with_capacity(candidate_support.len());
    for read in candidate_support.iter() {
        if is_forward_junction {
            new_ref_base_position = new_ref_base_position.min(read.alignment_start());
        } else {
            new_ref_base_position = new_ref_base_position.max(read.alignment_end());
        }
        non_junction_support.push(read.clone());
        RefSideSoftClip::check_add_ref_side_soft_clip(assembly.m().ref_side_soft_clips_mut(), &junction, read);
    }
    if non_junction_support.is_empty() {
        return;
    }
    let non_soft_clip_ref_position = new_ref_base_position;
    RefSideSoftClip::purge_ref_side_soft_clips(assembly.m().ref_side_soft_clips_mut(), non_soft_clip_ref_position);
    if assembly.ref_side_soft_clips().is_empty() {
        extend_assembly_ref_bases(assembly, non_soft_clip_ref_position, &non_junction_support, ref_genome, false);
        return;
    }
    java_sort(assembly.m().ref_side_soft_clips_mut(), |a, b| (-a.read_count()).cmp(&(-b.read_count())));
    let candidate_ref_side_soft_clip = assembly.ref_side_soft_clips()[0];
    let sc_pos = candidate_ref_side_soft_clip.position;
    let mut junction_soft_clipped = 0;
    let mut junction_non_soft_clipped = 0;
    for read in assembly.support() {
        if is_forward_junction {
            if read.is_left_clipped() && read.alignment_start() == sc_pos {
                junction_soft_clipped += 1;
            } else if !read.is_left_clipped() && read.alignment_start() < sc_pos {
                junction_non_soft_clipped += 1;
            }
        } else if read.is_right_clipped() && read.alignment_end() == sc_pos {
            junction_soft_clipped += 1;
        } else if !read.is_right_clipped() && read.alignment_end() > sc_pos {
            junction_non_soft_clipped += 1;
        }
    }
    let non_soft_clip_support = if is_forward_junction {
        non_junction_support.iter().filter(|x| !x.is_left_clipped()).filter(|x| x.alignment_start() < sc_pos).count() as i32
    } else {
        non_junction_support.iter().filter(|x| !x.is_right_clipped()).filter(|x| x.alignment_end() > sc_pos).count() as i32
    };
    let primary_ref_position;
    let primary_ref_position_support;
    let second_ref_position_support;
    let mut uses_soft_clipped_position = false;
    let majority_junc_reads_soft_clipped = junction_soft_clipped > junction_non_soft_clipped;
    if majority_junc_reads_soft_clipped || non_soft_clip_support <= candidate_ref_side_soft_clip.read_count() {
        primary_ref_position = sc_pos;
        primary_ref_position_support = candidate_ref_side_soft_clip.read_count();
        second_ref_position_support = non_soft_clip_support;
        uses_soft_clipped_position = true;
    } else {
        primary_ref_position = non_soft_clip_ref_position;
        primary_ref_position_support = non_soft_clip_support;
        second_ref_position_support = candidate_ref_side_soft_clip.read_count();
    }
    let second_ref_position_support_perc = second_ref_position_support as f64 / primary_ref_position_support as f64;
    let has_sufficient_second_ref_support =
        second_ref_position_support >= ASSEMBLY_SPLIT_MIN_READ_SUPPORT && second_ref_position_support_perc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
    if !allow_branching || !has_sufficient_second_ref_support {
        extend_assembly_ref_bases(assembly, primary_ref_position, &non_junction_support, ref_genome, uses_soft_clipped_position);
        return;
    }
    branch_assemblies_from_ref_bases(assembly, non_soft_clip_ref_position, &non_junction_support, primary_ref_position_support, ref_genome);
}

fn extend_assembly_ref_bases(assembly: AsmP, new_ref_position: i32, non_junction_support: &[ReadP], ref_genome: Option<&dyn RefGenomeInterface>, is_soft_clipped: bool) {
    if new_ref_position != assembly.ref_base_position() {
        if assembly.is_forward_junction() == (new_ref_position < assembly.ref_base_position()) {
            assembly.m().extend_ref_bases(new_ref_position, ref_genome);
        } else if is_soft_clipped {
            assembly.m().trim_ref_bases(new_ref_position);
        }
    }
    check_add_ref_base_support(assembly, non_junction_support, &HashSet::new());
}

fn branch_assemblies_from_ref_bases(original_assembly: AsmP, non_soft_clip_ref_position: i32, non_junction_support: &[ReadP], max_ref_side_support: i32, ref_genome: Option<&dyn RefGenomeInterface>) {
    let excluded_read_ids_list = allocate_excluded_reads(&original_assembly, non_junction_support);
    let initial_support: Vec<SupP> = original_assembly.support().clone();
    let mut branched_assemblies: Vec<AsmP> = vec![original_assembly];
    for (i, excluded_reads) in excluded_read_ids_list.iter().enumerate() {
        let non_excluded_non_junction_support: Vec<ReadP> = non_junction_support.iter().filter(|x| !excluded_reads.iter().any(|y| x.id().contains(y.as_str()))).cloned().collect();
        let junction_assembly;
        if i == 0 {
            junction_assembly = original_assembly;
            original_assembly.m().remove_support_reads(excluded_reads);
            extend_assembly_ref_bases(original_assembly, non_soft_clip_ref_position, &non_excluded_non_junction_support, ref_genome, false);
        } else {
            let ref_side_soft_clip = original_assembly.ref_side_soft_clips()[i - 1];
            if ref_side_soft_clip.matches_original() || ref_side_soft_clip.has_proximate_match(original_assembly.ref_base_position()) {
                continue;
            }
            let new_ref_base_length = (ref_side_soft_clip.position - original_assembly.junction().position).abs() + 1;
            let origin_ref_base_length = original_assembly.ref_base_length();
            let ref_base_difference = origin_ref_base_length - new_ref_base_length;
            if new_ref_base_length <= 0 || ref_base_difference <= 0 {
                break;
            }
            let new_support: Vec<SupP> = initial_support.iter().filter(|x| !excluded_reads.contains(x.id())).cloned().collect();
            junction_assembly = JunctionAssembly::new_branched(original_assembly, ref_side_soft_clip, new_ref_base_length, &new_support);
            extend_assembly_ref_bases(junction_assembly, ref_side_soft_clip.position, &non_excluded_non_junction_support, ref_genome, true);
        }
        if junction_assembly != original_assembly {
            let total_support = junction_assembly.support_count() as i32;
            let support_perc = total_support as f64 / max_ref_side_support as f64;
            let has_sufficient_second_ref_support = total_support >= ASSEMBLY_SPLIT_MIN_READ_SUPPORT && support_perc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
            if !has_sufficient_second_ref_support {
                continue;
            }
            branched_assemblies.push(junction_assembly);
            junction_assembly.m().set_outcome(AssemblyOutcome::DupBranched);
            break;
        }
    }
    for junction_assembly in branched_assemblies {
        purge_unrelated_ref_side_soft_clips(junction_assembly);
        if junction_assembly != original_assembly {
            PhaseGroup::add_derived_assembly(original_assembly.phase_group().expect("phase group of a branched assembly"), junction_assembly);
        }
    }
}

fn purge_unrelated_ref_side_soft_clips(assembly: AsmP) {
    let ref_base_position = assembly.ref_base_position();
    assembly.m().ref_side_soft_clips_mut().retain(|sc| sc.has_proximate_match(ref_base_position));
}

fn check_add_ref_base_support(assembly: AsmP, non_junction_reads: &[ReadP], excluded_read_ids: &HashSet<String>) {
    let mut secondary_support_reads: Vec<ReadP> = Vec::new();
    let permitted_mismatches = ASSEMBLY_EXTENSION_BASE_MISMATCH;
    for read in non_junction_reads {
        if excluded_read_ids.contains(read.id()) {
            continue;
        }
        if read.has_junction_mate() {
            let junc_mate = assembly.support().iter().find(|x| x.type_() == SupportType::Junction && x.matches_fragment_read(read, false)).cloned();
            if let Some(jm) = junc_mate {
                if jm.has_reference_mismatches() && jm.reference_mismatches() <= permitted_mismatches {
                    check_add_ref_base_read_overlap(assembly, read, SupportType::JunctionMate, 0);
                    continue;
                }
            }
        }
        secondary_support_reads.push(read.clone());
    }
    let sort_on_read_start = assembly.is_reverse_junction();
    java_sort(&mut secondary_support_reads, |a, b| {
        let ka = if sort_on_read_start { a.alignment_start() } else { -a.alignment_end() };
        let kb = if sort_on_read_start { b.alignment_start() } else { -b.alignment_end() };
        ka.cmp(&kb)
    });
    for read in &secondary_support_reads {
        let type_ = if read.has_junction_mate() { SupportType::JunctionMate } else { SupportType::Discordant };
        check_add_ref_base_read(assembly, read, type_);
    }
}

pub fn check_add_ref_base_read(assembly: AsmP, read: &ReadP, support_type: SupportType) -> bool {
    check_add_ref_base_read_overlap(assembly, read, support_type, ASSEMBLY_READ_OVERLAP_BASES)
}

fn check_add_ref_base_read_overlap(assembly: AsmP, read: &ReadP, support_type: SupportType, required_overlap: i32) -> bool {
    let Some(mut read_assembly_indices) = ReadAssemblyIndices::get_ref_read_indices(&assembly, assembly.ref_base_position(), read) else { return false };
    let read_start_index = read_assembly_indices.read_index_start;
    let permitted_mismatches = ASSEMBLY_EXTENSION_BASE_MISMATCH;
    let mut can_add_read = can_add_ref_base_read(assembly.bases(), assembly.base_quals(), read, &read_assembly_indices, required_overlap, permitted_mismatches, false);
    let read_len = read.get_bases().len() as i32;
    if !can_add_read && read_start_index < read_len - REF_READ_SEARCH_LENGTH {
        let read_test_end_index = read_start_index + REF_READ_SEARCH_LENGTH;
        let length = read_test_end_index - read_start_index + 1;
        if read_start_index < 0 || read_start_index >= read_len || read_start_index + length > read_len {
            crate::error!("refAssembly({}) invalid indices({} - {}) vs readBases({}) for ref extension read search", *assembly, read_start_index, read_test_end_index, read_len);
            return false;
        }
        let read_bases = &read.get_bases()[read_start_index as usize..(read_start_index + length) as usize];
        let assembly_bases = assembly.bases();
        if let Some(assembly_start_index) = index_of(assembly_bases, read_bases) {
            let from = (assembly_start_index + REF_READ_SEARCH_LENGTH) as usize;
            let next_start_index = if from <= assembly_bases.len() { index_of(&assembly_bases[from..], read_bases).map(|i| i + from as i32).unwrap_or(-1) } else { -1 };
            read_assembly_indices = ReadAssemblyIndices::new(read_start_index, read_assembly_indices.read_index_end, assembly_start_index);
            can_add_read = can_add_ref_base_read(assembly_bases, assembly.base_quals(), read, &read_assembly_indices, required_overlap, permitted_mismatches, next_start_index > 0);
        }
    }
    if !can_add_read {
        if support_type == SupportType::JunctionMate {
            let junction_read_start_distance = read_assembly_indices.junction_read_start_distance(assembly.junction_index());
            let support_read = SupportRead::new(read, support_type, junction_read_start_distance, 0, ASSEMBLY_EXTENSION_BASE_MISMATCH + 1);
            assembly.m().support_mut().push(support_read);
        }
        return false;
    }
    assembly.m().add_ref_base_support_read(read, &read_assembly_indices, support_type);
    true
}

fn can_add_ref_base_read(
    assembly_bases: &[u8], assembly_base_quals: &[u8], read: &Read, read_index_info: &ReadAssemblyIndices, required_overlap: i32, permitted_mismatches: i32, require_qual_matches: bool,
) -> bool {
    let mut mismatch_count = 0;
    let mut overlapped_base_count = 0;
    let mut assembly_index = read_index_info.assembly_index_start;
    let read_bases = read.get_bases();
    let read_quals = read.get_base_quality();
    let mut i = read_index_info.read_index_start;
    while i <= read_index_info.read_index_end {
        if i < 0 || assembly_index >= assembly_bases.len() as i32 || i >= read_bases.len() as i32 {
            break;
        }
        // a negative assembly index is an ArrayIndexOutOfBoundsException in Java
        let ai = assembly_index as usize;
        if assembly_bases[ai] == 0 {
            i += 1;
            assembly_index += 1;
            continue;
        }
        if above_min_qual(assembly_base_quals[ai]) {
            overlapped_base_count += 1;
        }
        let mut ref_base_qual = assembly_base_quals[ai];
        if ref_base_qual == 0 {
            ref_base_qual = if require_qual_matches {
                above_low_base_qual_threshold()
            } else {
                crate::common::base_qual::low_base_qual_threshold()
            };
        }
        if !bases_match_qual(read_bases[i as usize], assembly_bases[ai], read_quals[i as usize], ref_base_qual) {
            mismatch_count += 1;
            if mismatch_count > permitted_mismatches {
                return false;
            }
        }
        i += 1;
        assembly_index += 1;
    }
    overlapped_base_count >= required_overlap
}

fn allocate_excluded_reads(assembly: &JunctionAssembly, non_junction_reads: &[ReadP]) -> Vec<HashSet<String>> {
    let ref_side_soft_clips = assembly.ref_side_soft_clips();
    let is_forward_junction = assembly.is_forward_junction();
    let mut nsc_reads: HashSet<String> = HashSet::new();
    let positions = ref_side_soft_clips.iter().filter(|x| !x.matches_original()).map(|x| x.position);
    let max_ref_side_soft_clip_position = if is_forward_junction { positions.min().unwrap_or(0) } else { positions.max().unwrap_or(0) };
    let mut soft_clipped_read_ids: HashSet<String> = HashSet::new();
    for x in ref_side_soft_clips.iter().filter(|x| !x.matches_original()) {
        for id in x.read_ids().iter_unordered() {
            soft_clipped_read_ids.insert(id.clone());
        }
    }
    for read in assembly.support() {
        if soft_clipped_read_ids.contains(read.id()) {
            continue;
        }
        let add_read = if is_forward_junction {
            read.alignment_start() < max_ref_side_soft_clip_position || (!read.is_discordant() && read.is_mate_mapped() && read.mate_alignment_start() < max_ref_side_soft_clip_position)
        } else {
            read.alignment_end() > max_ref_side_soft_clip_position || (!read.is_discordant() && read.is_mate_mapped() && read.mate_alignment_end() > max_ref_side_soft_clip_position)
        };
        if add_read {
            nsc_reads.insert(read.id().to_string());
        }
    }
    for read in non_junction_reads {
        if soft_clipped_read_ids.contains(read.id()) {
            continue;
        }
        let add_read = if is_forward_junction {
            read.alignment_start() < max_ref_side_soft_clip_position || (read.has_junction_mate() && read.is_mate_mapped() && read.mate_alignment_start() < max_ref_side_soft_clip_position)
        } else {
            read.alignment_end() > max_ref_side_soft_clip_position || (read.has_junction_mate() && read.is_mate_mapped() && read.mate_alignment_end() > max_ref_side_soft_clip_position)
        };
        if add_read {
            nsc_reads.insert(read.id().to_string());
        }
    }
    let mut excluded_reads_list: Vec<HashSet<String>> = Vec::with_capacity(1 + ref_side_soft_clips.len());
    excluded_reads_list.push(soft_clipped_read_ids);
    for ref_side_soft_clip in ref_side_soft_clips.iter() {
        let mut excluded_reads: HashSet<String> = HashSet::new();
        if !ref_side_soft_clip.matches_original() {
            excluded_reads.extend(nsc_reads.iter().cloned());
        }
        for other in ref_side_soft_clips.iter() {
            if ref_side_soft_clip == other {
                continue;
            }
            for id in other.read_ids().iter_unordered() {
                excluded_reads.insert(id.clone());
            }
        }
        excluded_reads_list.push(excluded_reads);
    }
    excluded_reads_list
}

pub fn trim_assembly_ref_bases(assembly: AsmP, max_gap: i32) {
    let assembly_base_quals = assembly.base_quals();
    let ref_base_position = assembly.ref_base_position();
    let mut new_ref_base_position = ref_base_position;
    if assembly.is_forward_junction() {
        let mut trim_index = 0;
        let mut last_set_base = -1;
        for i in 0..assembly.junction_index() {
            if assembly_base_quals[i as usize] > 0 {
                if last_set_base >= 0 && i - last_set_base > max_gap {
                    trim_index = i;
                }
                last_set_base = i;
            } else if last_set_base < 0 {
                trim_index = i + 1;
            }
        }
        new_ref_base_position += trim_index;
    } else {
        let last_index = assembly.bases().len() as i32 - 1;
        let mut trim_index = last_index;
        let mut last_set_base = -1;
        let mut i = last_index;
        while i > assembly.junction_index() {
            if assembly_base_quals[i as usize] > 0 {
                if last_set_base >= 0 && last_set_base - i > max_gap {
                    trim_index = i;
                }
                last_set_base = i;
            } else if last_set_base < 0 {
                trim_index = i - 1;
            }
            i -= 1;
        }
        new_ref_base_position -= last_index - trim_index;
    }
    if new_ref_base_position != ref_base_position {
        assembly.m().trim_ref_bases(new_ref_base_position);
    }
}
