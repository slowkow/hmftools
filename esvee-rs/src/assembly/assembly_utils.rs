//! AssemblyUtils, with the hmf-common Nucleotides helpers it uses.

use crate::assembly::common_utils::below_min_qual;
use crate::assembly::constants::*;
use crate::assembly::read::read::Read;
use crate::assembly::types::junction::Junction;
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::phase_group::PhaseGroup;
use crate::assembly::types::repeat_info::calc_trimmed_base_length;
use crate::assembly::types::support_read::SupportRead;
use crate::assembly::types::AssemblyOutcome;

pub const DNA_BASE_BYTES: [u8; 5] = [65, 67, 71, 84, 78];
pub const DNA_BASE_COUNT: usize = 5;
pub const NO_BASE: u8 = 0;
pub const DNA_N_BYTE: u8 = b'N';
/// The quality given to bases that must count as not low quality (reference bases, for one): one above
/// the low band (Java: LOW_BASE_QUAL_THRESHOLD + 1).
#[inline]
pub fn above_low_base_qual_threshold() -> u8 {
    crate::common::base_qual::low_base_qual_threshold() + 1
}

/// Nucleotides.swapDnaBase
#[inline]
pub fn swap_dna_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base,
    }
}

pub fn reverse_complement_bases(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|b| swap_dna_base(*b)).collect()
}

pub fn reverse_complement_str(bases: &str) -> String {
    bases.bytes().rev().map(|b| swap_dna_base(b) as char).collect()
}

/// AssemblyUtils.baseIndex: A, C, G, T, then N for anything else
#[inline]
pub fn base_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

pub fn mismatches_per_comparison_length(sequence_length: i32) -> i32 {
    if sequence_length < 15 {
        return 0;
    }
    if sequence_length <= 100 {
        return 1;
    }
    (sequence_length as f64 / 200.0).ceil() as i32 + 1
}

pub fn read_qual_from_junction(read: &Read, junction: &Junction) -> i32 {
    let read_junction_index = read.get_read_index_at_reference_position_ext(junction.position, true);
    let (mut read_index_start, mut read_index_end) = if junction.is_forward() { (read_junction_index, read.bases_length() - 1) } else { (0, read_junction_index) };
    read_index_start = read_index_start.max(0);
    read_index_end = read_index_end.min(read.get_base_quality().len() as i32 - 1);
    let mut base_qual_total = 0;
    let quals = read.get_base_quality();
    let mut i = read_index_start;
    while i <= read_index_end {
        base_qual_total += quals[i as usize] as i8 as i32;
        i += 1;
    }
    base_qual_total
}

#[inline]
pub fn bases_match_qual(first: u8, second: u8, first_qual: u8, second_qual: u8) -> bool {
    first == second || first == DNA_N_BYTE || second == DNA_N_BYTE || below_min_qual(first_qual) || below_min_qual(second_qual)
}

pub fn bases_match(bases1: &[u8], bases2: &[u8]) -> bool {
    bases1 == bases2
}

pub fn bases_match_range(bases1: &[u8], first_index_start: i32, first_index_end: i32, bases2: &[u8], second_index_start: i32, second_index_end: i32) -> bool {
    if first_index_end - first_index_start != second_index_end - second_index_start {
        return false;
    }
    let mut first_index = first_index_start as usize;
    let mut second_index = second_index_start as usize;
    while first_index < bases1.len() && second_index < bases2.len() {
        if bases1[first_index] != bases2[second_index] {
            return false;
        }
        first_index += 1;
        second_index += 1;
    }
    true
}

pub fn is_local_assembly_candidate(first: &JunctionAssembly, second: &JunctionAssembly, check_concordant_reads: bool, check_line_insertion: bool) -> bool {
    let (j1, j2) = (first.junction(), second.junction());
    if j1.chromosome != j2.chromosome {
        return false;
    }
    if first.is_forward_junction() == second.is_forward_junction() {
        return false;
    }
    let first_is_lower = j1.position <= j2.position;
    let is_del_type = first_is_lower == first.is_forward_junction();
    let junction_distance = (j1.position - j2.position).abs();
    if (is_del_type && junction_distance > PROXIMATE_DEL_LENGTH) || (!is_del_type && junction_distance > PROXIMATE_DUP_LENGTH) {
        return false;
    }
    if !check_concordant_reads && !check_line_insertion {
        return true;
    }
    if check_line_insertion && (first.has_line_sequence() || second.has_line_sequence()) {
        return true;
    }
    if check_concordant_reads {
        let lower_assembly = if first_is_lower { first } else { second };
        let upper_assembly = if !first_is_lower { first } else { second };
        let lower_junction = if first_is_lower { j1 } else { j2 };
        let upper_junction = if !first_is_lower { j1 } else { j2 };
        if lower_assembly.support().iter().any(|x| is_crossing_concordant_read(x, &upper_junction, false))
            || upper_assembly.support().iter().any(|x| is_crossing_concordant_read(x, &lower_junction, true))
        {
            return true;
        }
    }
    false
}

fn is_crossing_concordant_read(read: &SupportRead, junction: &Junction, require_lower: bool) -> bool {
    if read.is_discordant() || !read.is_paired_read() || read.is_mate_unmapped() {
        return false;
    }
    if require_lower { read.mate_alignment_end() <= junction.position } else { read.mate_alignment_start() >= junction.position }
}

pub fn is_supplementary_only(assembly: &JunctionAssembly) -> bool {
    assembly.support().iter().all(|x| x.is_supplementary())
}

pub fn create_base_quals_above_min_threshold(length: usize) -> Vec<u8> {
    vec![above_low_base_qual_threshold(); length]
}

pub fn extract_insert_sequence(first: &JunctionAssembly, first_reversed: bool, second: &JunctionAssembly, second_reversed: bool, insert_length: i32) -> String {
    let mut ext_base_index_start;
    let mut ext_base_index_end;
    if first.is_forward_junction() {
        ext_base_index_start = first.junction_index() + 1;
        ext_base_index_end = (ext_base_index_start + insert_length - 1).min(first.base_length() - 1);
    } else {
        ext_base_index_end = first.junction_index() - 1;
        ext_base_index_start = (ext_base_index_end - insert_length + 1).max(0);
    }
    let mut insert_sequence = first.form_sequence(ext_base_index_start, ext_base_index_end);
    if first_reversed {
        insert_sequence = reverse_complement_str(&insert_sequence);
    }
    if insert_length <= first.extension_length() {
        return insert_sequence;
    }
    let remaining_insert_length = insert_length - first.extension_length();
    if second.is_forward_junction() {
        ext_base_index_start = second.junction_index() + 1;
        ext_base_index_end = (ext_base_index_start + remaining_insert_length - 1).min(second.base_length() - 1);
    } else {
        ext_base_index_end = second.junction_index() - 1;
        ext_base_index_start = (ext_base_index_end - remaining_insert_length + 1).max(0);
    }
    let mut extra_insert_sequence = second.form_sequence(ext_base_index_start, ext_base_index_end);
    if second_reversed {
        extra_insert_sequence = reverse_complement_str(&extra_insert_sequence);
    }
    insert_sequence + &extra_insert_sequence
}

pub fn calc_trimmed_ref_base_length(assembly: &JunctionAssembly) -> i32 {
    let ref_base_length = assembly.ref_base_length();
    if assembly.repeat_info().is_empty() {
        return ref_base_length;
    }
    let seq_start = if assembly.is_forward_junction() { assembly.junction_index() - ref_base_length + 1 } else { assembly.junction_index() };
    let seq_end = if assembly.is_forward_junction() { assembly.junction_index() } else { assembly.junction_index() + ref_base_length - 1 };
    calc_trimmed_base_length(seq_start, seq_end, assembly.repeat_info())
}

pub fn calc_trimmed_extension_base_length(assembly: &JunctionAssembly) -> i32 {
    let ext_base_length = assembly.extension_length();
    if assembly.repeat_info().is_empty() {
        return ext_base_length;
    }
    let seq_start = if assembly.is_forward_junction() { assembly.junction_index() + 1 } else { 0 };
    let seq_end = if assembly.is_forward_junction() { seq_start + ext_base_length - 1 } else { ext_base_length - 1 };
    calc_trimmed_base_length(seq_start, seq_end, assembly.repeat_info())
}

pub fn find_matching_assembly(assemblies: &[AsmP], assembly: AsmP, require_extension_match: bool) -> Option<AsmP> {
    assemblies
        .iter()
        .filter(|x| **x != assembly)
        .filter(|x| x.junction().compare_to(&assembly.junction()) == std::cmp::Ordering::Equal)
        .find(|x| !require_extension_match || x.extension_length() == assembly.extension_length())
        .copied()
}

pub fn set_assembly_outcome(assembly: AsmP) {
    if assembly.outcome() != AssemblyOutcome::Unset {
        return;
    }
    let Some(phase_group) = assembly.phase_group() else {
        assembly.m().set_outcome(AssemblyOutcome::NoLink);
        return;
    };
    let secondary_split_links = PhaseGroup::find_secondary_split_links(&phase_group, assembly);
    if !secondary_split_links.is_empty() {
        assembly.m().set_outcome(AssemblyOutcome::Secondary);
        return;
    }
    if is_supplementary_only(&assembly) {
        assembly.m().set_outcome(AssemblyOutcome::SuppOnly);
        return;
    }
    assembly.m().set_outcome(AssemblyOutcome::NoLink);
}

pub fn has_unset_bases(assembly: &JunctionAssembly) -> bool {
    !find_unset_bases(assembly.bases()).is_empty()
}

pub fn find_unset_bases(bases: &[u8]) -> Vec<[i32; 2]> {
    let mut empty_ranges: Vec<[i32; 2]> = Vec::new();
    let mut open = false;
    for (i, b) in bases.iter().enumerate() {
        if *b == 0 {
            if !open {
                empty_ranges.push([i as i32, -1]);
                open = true;
            }
        } else if open {
            empty_ranges.last_mut().unwrap()[1] = i as i32 - 1;
            open = false;
        }
    }
    if open {
        empty_ranges.last_mut().unwrap()[1] = bases.len() as i32 - 1;
    }
    empty_ranges
}

pub fn non_null_base_str(bases: &[u8]) -> String {
    bases.iter().filter(|b| **b != 0).map(|b| *b as char).collect()
}
