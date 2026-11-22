//! PhaseGroupBuilder with LocalGroupBuilder and RemoteGroupBuilder: puts assemblies that are close together
//! or share reads through remote regions into phase groups.
//!
//! The Java classes run these two stages on several threads, and the grouping then depends on the order in
//! which threads link assemblies. Here both stages run in the order of a single Java thread (they take a
//! negligible part of the run time), so the result does not depend on the thread count.

use std::collections::BTreeMap;
use std::io::Write;

use crate::assembly::assembly_utils::is_local_assembly_candidate;
use crate::assembly::common_utils::positions_overlap;
use crate::assembly::config::{AssemblyConfig, WriteType};
use crate::assembly::constants::*;
use crate::assembly::globals::{dev_debug, has_paired_reads};
use crate::assembly::phase::assembly_linker::is_assembly_indel_link;
use crate::assembly::read::read::ReadP;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::junction_group::{JuncGroupP, JunctionGroup};
use crate::assembly::perf::PerfCounter;
use crate::assembly::types::phase_group::{PhaseGroup, PhaseGroupP};
use crate::assembly::types::remote_region::RemoteRegion;
use crate::assembly::types::support_read::SupP;
use crate::common::java_hash::{JavaHash, JavaHashSet};
use crate::common::java_id_hash::IdHashSet;
use crate::common::java_sort::java_sort;

/// A junction group as a hash set key. JunctionGroup inherits equals and hashCode from BaseRegion, which use
/// only the start and end positions, so two groups with the same bounds count as one (also across chromosomes).
#[derive(Clone)]
struct JunctionGroupKey(JuncGroupP);

impl PartialEq for JunctionGroupKey {
    fn eq(&self, other: &Self) -> bool {
        self.0.start() == other.0.start() && self.0.end() == other.0.end()
    }
}
impl Eq for JunctionGroupKey {}
impl std::hash::Hash for JunctionGroupKey {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        (self.0.start(), self.0.end()).hash(state)
    }
}
impl JavaHash for JunctionGroupKey {
    fn java_hash(&self) -> i32 {
        let result = 31i32.wrapping_add(self.0.start());
        result.wrapping_mul(31).wrapping_add(self.0.end())
    }
}

/// PhaseGroupBuildWriter
pub struct PhaseGroupBuildWriter {
    writer: Option<std::io::BufWriter<std::fs::File>>,
}

impl PhaseGroupBuildWriter {
    pub fn new(config: &AssemblyConfig) -> PhaseGroupBuildWriter {
        let mut writer = None;
        if config.write_types.contains(&WriteType::PhaseGroupBuilding) {
            match std::fs::File::create(config.output_filename(WriteType::PhaseGroupBuilding)) {
                Ok(f) => {
                    let mut w = std::io::BufWriter::new(f);
                    let _ = writeln!(w, "Assembly\tLinkingAssembly\tLinkType\tPhaseGroupSize");
                    writer = Some(w);
                }
                Err(e) => crate::error!("failed to initialise phase group writer: {}", e),
            }
        }
        PhaseGroupBuildWriter { writer }
    }
    pub fn disabled() -> PhaseGroupBuildWriter {
        PhaseGroupBuildWriter { writer: None }
    }
    pub fn enabled(&self) -> bool {
        self.writer.is_some()
    }
    pub fn write_phase_group_build(&mut self, assembly: AsmP, linking_assembly: Option<AsmP>, link_type: &str, phase_group_size: usize) {
        let Some(w) = self.writer.as_mut() else { return };
        let assembly_info = assembly.junction().to_string();
        let linking_assembly_info = linking_assembly.map(|a| a.junction().to_string()).unwrap_or_default();
        let _ = writeln!(w, "{assembly_info}\t{linking_assembly_info}\t{link_type}\t{phase_group_size}");
    }
    pub fn close(&mut self) {
        if let Some(w) = self.writer.as_mut() {
            let _ = w.flush();
        }
    }
}

pub struct PhaseGroupBuilder<'a> {
    config: &'a AssemblyConfig,
    junction_group_map: &'a BTreeMap<String, Vec<JuncGroupP>>,
    phase_groups: Vec<PhaseGroupP>,
    perf_counters: Vec<PerfCounter>,
}

impl<'a> PhaseGroupBuilder<'a> {
    pub fn new(config: &'a AssemblyConfig, junction_group_map: &'a BTreeMap<String, Vec<JuncGroupP>>) -> PhaseGroupBuilder<'a> {
        for junction_groups in junction_group_map.values() {
            for (index, g) in junction_groups.iter().enumerate() {
                g.m().set_index(index as i32);
            }
        }
        PhaseGroupBuilder { config, junction_group_map, phase_groups: Vec::new(), perf_counters: Vec::new() }
    }

    pub fn phase_groups(&self) -> &Vec<PhaseGroupP> {
        &self.phase_groups
    }
    /// the LocalPhaseGroups and RemotePhaseGroups counters of build_groups
    pub fn perf_counters(&self) -> &[PerfCounter] {
        &self.perf_counters
    }

    pub fn build_groups(&mut self, writer: &mut PhaseGroupBuildWriter) {
        let mut all_junction_groups: Vec<JuncGroupP> = Vec::new();
        for junction_groups in self.junction_group_map.values() {
            all_junction_groups.extend(junction_groups.iter().filter(|x| !x.junction_assemblies().is_empty()).copied());
        }
        // stage one: local phase groups
        let mut local_phase_groups: IdHashSet<PhaseGroupP> = IdHashSet::new();
        let mut local_perf = PerfCounter::new("LocalPhaseGroups");
        for junction_group in &all_junction_groups {
            let start = std::time::Instant::now();
            form_local_phase_groups(*junction_group, &mut local_phase_groups, writer);
            local_perf.stop_check_log(start, || format!("juncGroup({})", **junction_group));
        }
        self.perf_counters.push(local_perf);
        self.phase_groups.extend(local_phase_groups.iter().copied());
        crate::info!("building remote phase groups, current group count({})", self.phase_groups.len());
        // stage two: groups linked through remote regions
        let mut remote_perf = PerfCounter::new("RemotePhaseGroups");
        let (remote_phase_groups, removed_phase_groups) = self.build_remote_groups_timed(&all_junction_groups, writer, Some(&mut remote_perf));
        self.perf_counters.push(remote_perf);
        self.phase_groups.extend(remote_phase_groups.iter().copied());
        for y in &removed_phase_groups {
            if let Some(i) = self.phase_groups.iter().position(|x| x == y) {
                self.phase_groups.remove(i);
            }
        }
        crate::info!("phase group building complete, final group count({})", self.phase_groups.len());
        for junction_group in &all_junction_groups {
            for assembly in junction_group.junction_assemblies() {
                if assembly.phase_group().is_none() {
                    self.phase_groups.push(PhaseGroup::new(*assembly, None));
                }
            }
        }
        java_sort(&mut self.phase_groups, |a, b| a.assemblies()[0].junction().compare_to(&b.assemblies()[0].junction()));
        for (i, pg) in self.phase_groups.iter().enumerate() {
            pg.m().set_id(i as i32);
        }
        if dev_debug() {
            for junction_groups in self.junction_group_map.values() {
                for g in junction_groups {
                    for assembly in g.junction_assemblies() {
                        let count = self.phase_groups.iter().filter(|x| x.assemblies().contains(assembly)).count();
                        if count > 1 {
                            crate::error!("asm({}) in {} phase groups", **assembly, count);
                        }
                    }
                }
            }
        }
    }

    /// RemoteGroupBuilder.run over the junction groups in order: the phase groups formed and those merged away
    pub fn build_remote_groups(&self, junction_groups: &[JuncGroupP], writer: &mut PhaseGroupBuildWriter) -> (IdHashSet<PhaseGroupP>, Vec<PhaseGroupP>) {
        self.build_remote_groups_timed(junction_groups, writer, None)
    }

    fn build_remote_groups_timed(&self, junction_groups: &[JuncGroupP], writer: &mut PhaseGroupBuildWriter, mut perf: Option<&mut PerfCounter>) -> (IdHashSet<PhaseGroupP>, Vec<PhaseGroupP>) {
        let mut remote_phase_groups: IdHashSet<PhaseGroupP> = IdHashSet::new();
        let mut removed_phase_groups: Vec<PhaseGroupP> = Vec::new();
        for junction_group in junction_groups {
            let start = std::time::Instant::now();
            for assembly in junction_group.junction_assemblies().clone() {
                self.find_remote_phased_assemblies(*junction_group, assembly, &mut remote_phase_groups, &mut removed_phase_groups, writer);
            }
            if let Some(perf) = perf.as_deref_mut() {
                perf.stop_check_log(start, || format!("juncGroup({})", **junction_group));
            }
        }
        (remote_phase_groups, removed_phase_groups)
    }

    fn find_remote_phased_assemblies(&self, junction_group: JuncGroupP, assembly: AsmP, phase_groups: &mut IdHashSet<PhaseGroupP>, removed_phase_groups: &mut Vec<PhaseGroupP>, writer: &mut PhaseGroupBuildWriter) {
        if assembly.remote_regions().is_empty() {
            return;
        }
        let mut linked_junction_groups: JavaHashSet<JunctionGroupKey> = JavaHashSet::new();
        for region in assembly.remote_regions() {
            if self.is_filtered(region) {
                continue;
            }
            for g in self.find_overlapping_junction_groups(junction_group, region) {
                linked_junction_groups.insert(JunctionGroupKey(g));
            }
        }
        for other_junction_group in linked_junction_groups.to_vec() {
            for other_assembly in other_junction_group.0.junction_assemblies().clone() {
                if assembly == other_assembly {
                    continue;
                }
                if assemblies_already_phased(assembly, other_assembly) {
                    return;
                }
                if !self.can_phase_assemblies(assembly, other_assembly) {
                    continue;
                }
                if !assemblies_already_phased(assembly, other_assembly) {
                    link_to_phase_groups(assembly.phase_group(), assembly, other_assembly, phase_groups, Some(removed_phase_groups), writer, "Remote");
                }
            }
        }
    }

    fn can_phase_assemblies(&self, assembly: AsmP, other_assembly: AsmP) -> bool {
        let Some(overlapping_region) = assembly.remote_regions().iter().find(|x| x.overlaps_assembly(&other_assembly)) else { return false };
        if has_paired_reads() && !assemblies_share_reads(overlapping_region, other_assembly, REMOTE_PHASING_MIN_READS, self.config.apply_remote_phasing_read_check_threshold) {
            return false;
        }
        true
    }

    fn find_overlapping_junction_groups(&self, assembly_junction_group: JuncGroupP, region: &RemoteRegion) -> Vec<JuncGroupP> {
        let Some(junction_groups) = self.junction_group_map.get(region.chromosome) else { return Vec::new() };
        let start_index: i32 = if assembly_junction_group.chromosome() == region.chromosome && assembly_junction_group.overlaps_remote_region(region) {
            assembly_junction_group.index()
        } else {
            JunctionGroup::binary_search(region.start(), junction_groups) as i32
        };
        let mut overlap_groups: Vec<JuncGroupP> = Vec::new();
        for search_up in [true, false] {
            let mut index = if search_up { start_index + 1 } else { start_index };
            while index >= 0 && (index as usize) < junction_groups.len() {
                let g = junction_groups[index as usize];
                if !g.overlaps_remote_region(region) {
                    break;
                }
                overlap_groups.push(g);
                if search_up {
                    index += 1;
                } else {
                    index -= 1;
                }
            }
        }
        overlap_groups
    }

    fn is_filtered(&self, region: &RemoteRegion) -> bool {
        let specific = &self.config.specific_chr_regions;
        if !specific.has_filters() {
            return false;
        }
        if !specific.include_chromosome(region.chromosome) {
            return true;
        }
        if specific.regions.is_empty() {
            return false;
        }
        !specific.regions.iter().any(|x| x.chromosome == region.chromosome && positions_overlap(x.start, x.end, region.start(), region.end()))
    }
}

fn assemblies_already_phased(first: AsmP, second: AsmP) -> bool {
    first.phase_group().is_some() && first.phase_group() == second.phase_group()
}

const NO_FRAG_CHECK_THRESHOLD: i32 = -1;

fn calculate_max_fragment_check_threshold(first_read_count: i32, second_read_count: i32) -> i32 {
    let max_checks: i64 = first_read_count.wrapping_mul(second_read_count) as i64;
    let round = |v: f64| crate::assembly::types::assembly_stats::java_round(v) as i32;
    if max_checks < 500 {
        return NO_FRAG_CHECK_THRESHOLD;
    }
    if max_checks <= 1000 {
        return round(0.6 * max_checks as f64);
    }
    if max_checks <= 2000 {
        return round(0.5 * max_checks as f64);
    }
    if max_checks <= 5000 {
        return round(0.4 * max_checks as f64);
    }
    if max_checks <= 10000 {
        return round(0.3 * max_checks as f64);
    }
    if max_checks <= 250000 {
        return round(0.1 * max_checks as f64);
    }
    50000
}

fn has_matching_support_read(support: &[SupP], read_id: &str) -> bool {
    support.iter().any(|x| x.id() == read_id)
}

fn has_matching_read(support: &[ReadP], read_id: &str) -> bool {
    support.iter().any(|x| x.id() == read_id)
}

fn assemblies_share_reads(first_region: &RemoteRegion, second: AsmP, min_shared_reads: i32, apply_match_threshold: bool) -> bool {
    let first_read_count = first_region.read_ids().len() as i32;
    let second_read_count = (second.support_count() + second.candidate_support().len()) as i32;
    if first_read_count < min_shared_reads || second_read_count < min_shared_reads {
        return false;
    }
    let max_match_checks = if apply_match_threshold { calculate_max_fragment_check_threshold(first_read_count, second_read_count) } else { NO_FRAG_CHECK_THRESHOLD };
    let mut current_checks: i32 = 0;
    let mut matched_count = 0;
    for read_id in first_region.read_ids().iter() {
        if has_matching_support_read(second.support(), read_id) {
            matched_count += 1;
            if matched_count >= min_shared_reads {
                return true;
            }
            continue;
        }
        current_checks += second.support_count() as i32;
        if max_match_checks != NO_FRAG_CHECK_THRESHOLD && current_checks > max_match_checks {
            return false;
        }
        if has_matching_read(second.candidate_support(), read_id) {
            matched_count += 1;
            if matched_count >= min_shared_reads {
                return true;
            }
            continue;
        }
        current_checks += second.candidate_support().len() as i32;
        if max_match_checks != NO_FRAG_CHECK_THRESHOLD && current_checks > max_match_checks {
            return false;
        }
    }
    false
}

/// PhaseGroupBuilder.linkToPhaseGroups
pub fn link_to_phase_groups(
    phase_group: Option<PhaseGroupP>, assembly: AsmP, other_assembly: AsmP, phase_groups: &mut IdHashSet<PhaseGroupP>, removed_phase_groups: Option<&mut Vec<PhaseGroupP>>,
    writer: &mut PhaseGroupBuildWriter, link_type: &str,
) {
    match phase_group {
        None => match other_assembly.phase_group() {
            Some(other_pg) => {
                if writer.enabled() {
                    writer.write_phase_group_build(assembly, Some(other_assembly), link_type, other_pg.assembly_count());
                }
                PhaseGroup::add_assembly(other_pg, assembly);
            }
            None => {
                phase_groups.insert(PhaseGroup::new(assembly, Some(other_assembly)));
                if writer.enabled() {
                    writer.write_phase_group_build(assembly, Some(other_assembly), link_type, 0);
                    writer.write_phase_group_build(other_assembly, Some(assembly), link_type, 1);
                }
            }
        },
        Some(phase_group) => match other_assembly.phase_group() {
            Some(other_pg) => {
                let dest_phase_group = if phase_group.assembly_count() >= other_pg.assembly_count() { phase_group } else { other_pg };
                let src_phase_group = if dest_phase_group == phase_group { other_pg } else { phase_group };
                if writer.enabled() {
                    writer.write_phase_group_build(assembly, Some(other_assembly), link_type, dest_phase_group.assembly_count());
                }
                PhaseGroup::transfer_assemblies(dest_phase_group, src_phase_group);
                match removed_phase_groups {
                    Some(removed) => removed.push(src_phase_group),
                    None => {
                        phase_groups.remove(&src_phase_group);
                    }
                }
            }
            None => {
                if writer.enabled() {
                    writer.write_phase_group_build(other_assembly, Some(assembly), link_type, phase_group.assembly_count());
                }
                PhaseGroup::add_assembly(phase_group, other_assembly);
            }
        },
    }
}

/// LocalGroupBuilder.formLocalPhaseGroups
pub fn form_local_phase_groups(junction_group: JuncGroupP, phase_groups: &mut IdHashSet<PhaseGroupP>, writer: &mut PhaseGroupBuildWriter) {
    if junction_group.junction_assemblies().len() <= 1 {
        return;
    }
    let mut pos_junction_assemblies: Vec<AsmP> = junction_group.junction_assemblies().iter().filter(|x| x.is_forward_junction() && !x.discordant_only()).copied().collect();
    let mut neg_junction_assemblies: Vec<AsmP> = junction_group.junction_assemblies().iter().filter(|x| x.is_reverse_junction() && !x.discordant_only()).copied().collect();
    if pos_junction_assemblies.is_empty() || neg_junction_assemblies.is_empty() {
        return;
    }
    java_sort(&mut pos_junction_assemblies, |a, b| a.junction().position.cmp(&b.junction().position));
    java_sort(&mut neg_junction_assemblies, |a, b| a.junction().position.cmp(&b.junction().position));
    let mut neg_index = 0usize;
    let mut start_neg_assembly = neg_junction_assemblies[neg_index];
    for pos_assembly in pos_junction_assemblies {
        while start_neg_assembly.junction().position < pos_assembly.junction().position - PROXIMATE_DEL_LENGTH {
            neg_index += 1;
            if neg_index >= neg_junction_assemblies.len() {
                return;
            }
            start_neg_assembly = neg_junction_assemblies[neg_index];
        }
        let mut next_neg_index = neg_index;
        let mut neg_assembly = start_neg_assembly;
        while neg_assembly.junction().position < pos_assembly.junction().position + PROXIMATE_DEL_LENGTH {
            if (pos_assembly.phase_group().is_none() || neg_assembly.phase_group() != pos_assembly.phase_group())
                && (is_local_assembly_candidate(&pos_assembly, &neg_assembly, true, true) || is_assembly_indel_link_safe(pos_assembly, neg_assembly) || is_local_facing_link_candidate(pos_assembly, neg_assembly))
            {
                link_to_phase_groups(pos_assembly.phase_group(), pos_assembly, neg_assembly, phase_groups, None, writer, "Local");
            }
            next_neg_index += 1;
            if next_neg_index >= neg_junction_assemblies.len() {
                break;
            }
            neg_assembly = neg_junction_assemblies[next_neg_index];
        }
    }
}

fn is_assembly_indel_link_safe(a: AsmP, b: AsmP) -> bool {
    is_assembly_indel_link(&a, &b)
}

fn is_local_facing_link_candidate(first: AsmP, second: AsmP) -> bool {
    if first.junction().chromosome != second.junction().chromosome || first.junction().orient == second.junction().orient {
        return false;
    }
    let lower = if first.junction().position < second.junction().position { first } else { second };
    let upper = if first == lower { second } else { first };
    if !lower.junction().is_reverse() {
        return false;
    }
    let link_distance = upper.junction().position - lower.junction().position + 1;
    if !(PHASED_ASSEMBLY_MIN_TI..=PHASED_ASSEMBLY_MAX_TI).contains(&link_distance) {
        return false;
    }
    for is_start in [true, false] {
        let assembly = if is_start { lower } else { upper };
        let other_assembly = if is_start { upper } else { lower };
        let other_is_upper = is_start;
        let other_junction_position = other_assembly.junction().position;
        for read in assembly.support() {
            if other_is_upper {
                if read.is_right_clipped() && read.alignment_end() == other_junction_position {
                    return true;
                }
                if read.alignment_end() > other_junction_position + PROXIMATE_REF_SIDE_SOFT_CLIPS {
                    return false;
                }
            } else {
                if read.is_left_clipped() && read.alignment_start() == other_junction_position {
                    return true;
                }
                if read.alignment_start() < other_junction_position - PROXIMATE_REF_SIDE_SOFT_CLIPS {
                    return false;
                }
            }
        }
        for read in assembly.candidate_support() {
            if !read.has_junction_mate() {
                continue;
            }
            if other_is_upper {
                if read.is_right_clipped() && read.alignment_end() == other_junction_position {
                    return true;
                }
                if read.alignment_end() > other_junction_position + PROXIMATE_REF_SIDE_SOFT_CLIPS {
                    return false;
                }
            } else {
                if read.is_left_clipped() && read.alignment_start() == other_junction_position {
                    return true;
                }
                if read.alignment_start() < other_junction_position - PROXIMATE_REF_SIDE_SOFT_CLIPS {
                    return false;
                }
            }
        }
    }
    true
}
