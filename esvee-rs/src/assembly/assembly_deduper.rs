//! AssemblyDeduper: removes one of two assemblies at close junctions with matching sequence.

use crate::assembly::constants::*;
use crate::assembly::globals::proximate_junction_distance;
use crate::assembly::indel_builder::indel_assemblies_match;
use crate::assembly::read::read_utils::record_soft_clips_at_junction;
use crate::assembly::sequence_compare::matched_assembly_sequences;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::SupportType;

pub fn dedup_proximate_assemblies(existing_assemblies: &mut Vec<AsmP>, new_assemblies: &mut Vec<AsmP>, deduped_indels: &mut Vec<AsmP>) {
    if new_assemblies.is_empty() || existing_assemblies.is_empty() {
        return;
    }
    let mut index = existing_assemblies.len() as isize - 1;
    let min_position = new_assemblies[0].junction().position - proximate_junction_distance();
    let mut remove_existing: Vec<AsmP> = Vec::new();
    while index >= 0 && !new_assemblies.is_empty() {
        let assembly = existing_assemblies[index as usize];
        if assembly.junction().position < min_position {
            break;
        }
        let mut new_index = 0usize;
        while new_index < new_assemblies.len() {
            let new_assembly = new_assemblies[new_index];
            if new_assembly.junction() == assembly.junction() || new_assembly.junction().orient != assembly.junction().orient {
                new_index += 1;
                continue;
            }
            if new_assembly.indel() && assembly.indel() {
                if check_dedup_indel_pair(new_assembly, assembly, new_index, new_assemblies, existing_assemblies, &mut remove_existing, deduped_indels) {
                    new_index += 1;
                }
            } else if check_dedup_non_indel(new_assembly, assembly, new_index, new_assemblies, &mut remove_existing) {
                new_index += 1;
            }
        }
        index -= 1;
    }
    for x in remove_existing {
        if let Some(i) = existing_assemblies.iter().position(|a| *a == x) {
            existing_assemblies.remove(i);
        }
    }
}

fn check_dedup_non_indel(new_assembly: AsmP, assembly: AsmP, new_index: usize, new_assemblies: &mut Vec<AsmP>, remove_existing: &mut Vec<AsmP>) -> bool {
    if !can_dedup_assemblies(assembly, new_assembly) {
        return true;
    }
    if select_first_assembly(assembly, new_assembly) {
        assembly.m().add_merged_assembly();
        new_assemblies.remove(new_index);
        false
    } else {
        new_assembly.m().add_merged_assembly();
        remove_existing.push(assembly);
        true
    }
}

fn check_dedup_indel_pair(
    new_assembly: AsmP, assembly: AsmP, new_index: usize, new_assemblies: &mut Vec<AsmP>, existing_assemblies: &[AsmP], remove_existing: &mut Vec<AsmP>,
    deduped_indels: &mut Vec<AsmP>,
) -> bool {
    let matched_indel: Option<AsmP> = if new_assembly.indel() { deduped_indels.iter().find(|x| indel_assemblies_match(x, &new_assembly)).copied() } else { None };
    let dedup_new_assembly;
    match matched_indel {
        None => {
            if !can_dedup_assemblies(assembly, new_assembly) {
                return true;
            }
            dedup_new_assembly = select_first_assembly(assembly, new_assembly);
        }
        Some(m) => {
            dedup_new_assembly = true;
            if let Some(i) = deduped_indels.iter().position(|a| *a == m) {
                deduped_indels.remove(i);
            }
        }
    }
    if dedup_new_assembly {
        assembly.m().add_merged_assembly();
        new_assemblies.remove(new_index);
        if matched_indel.is_none() {
            let matched_existing_indel = existing_assemblies.iter().find(|x| indel_assemblies_match(x, &new_assembly)).copied();
            match matched_existing_indel {
                Some(m) => {
                    if existing_assemblies.iter().any(|x| indel_assemblies_match(x, &assembly)) {
                        remove_existing.push(m);
                    }
                }
                None => deduped_indels.push(new_assembly),
            }
        }
        false
    } else {
        new_assembly.m().add_merged_assembly();
        remove_existing.push(assembly);
        let matched_existing_indel = existing_assemblies.iter().find(|x| indel_assemblies_match(x, &assembly)).copied();
        if let Some(m) = matched_existing_indel {
            if existing_assemblies.iter().any(|x| indel_assemblies_match(x, &new_assembly)) {
                remove_existing.push(m);
            }
        }
        true
    }
}

fn can_dedup_assemblies(first: AsmP, second: AsmP) -> bool {
    if dedup_on_support_ratio(first, second) {
        return true;
    }
    if !matched_assembly_sequences(&first, &second) {
        return false;
    }
    if first.indel() == second.indel() {
        return true;
    }
    let junction_distance = (first.junction().position - second.junction().position).abs();
    if junction_distance <= ASSEMBLY_DEDUP_JITTER_MAX_DIST {
        return true;
    }
    if let Some(second_indel_coords) = second.indel_coords() {
        for first_read in first.support() {
            if first_read.type_() != SupportType::Junction {
                continue;
            }
            if first_read.indel_coords().is_some_and(|c| c.matches(second_indel_coords)) {
                return false;
            }
        }
    }
    if let Some(first_indel_coords) = first.indel_coords() {
        for second_read in second.support() {
            if second_read.type_() != SupportType::Junction {
                continue;
            }
            if second_read.indel_coords().is_some_and(|c| c.matches(first_indel_coords)) {
                return false;
            }
        }
    }
    true
}

fn dedup_on_support_ratio(first: AsmP, second: AsmP) -> bool {
    let junction_distance = (first.junction().position - second.junction().position).abs();
    if junction_distance > ASSEMBLY_DEDUP_JITTER_MAX_DIST {
        return false;
    }
    assembly_exceeds_support_ratio(first, second) || assembly_exceeds_support_ratio(second, first)
}

fn assembly_exceeds_support_ratio(first: AsmP, second: AsmP) -> bool {
    first.support_count() as i32 > second.support_count() as i32 * ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO
}

fn select_first_assembly(first: AsmP, second: AsmP) -> bool {
    if assembly_exceeds_support_ratio(first, second) {
        return true;
    }
    if assembly_exceeds_support_ratio(second, first) {
        return false;
    }
    if first.indel() != second.indel() {
        return first.indel();
    }
    if first.is_saga_matched() != second.is_saga_matched() {
        return first.is_saga_matched();
    }
    let first_has_germline_support = first.support().iter().any(|x| x.is_reference());
    let second_has_germline_support = second.support().iter().any(|x| x.is_reference());
    if first_has_germline_support != second_has_germline_support {
        return first_has_germline_support;
    }
    if first.discordant_only() != second.discordant_only() {
        return !first.discordant_only();
    }
    if first.support_count() != second.support_count() {
        return first.support_count() > second.support_count();
    }
    let precise = |a: AsmP| a.support().iter().filter(|x| record_soft_clips_at_junction(x.cached_read().expect("support read without cached read"), &a.junction())).count();
    precise(first) >= precise(second)
}
