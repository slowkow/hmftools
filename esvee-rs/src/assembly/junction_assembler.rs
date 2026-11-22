//! JunctionAssembler: builds the assemblies for one junction from the reads around it.

use std::collections::HashSet;

use crate::assembly::constants::*;
use crate::assembly::discordant_junction_builder::DiscordantJunctionBuilder;
use crate::assembly::extension_seq_builder::{ExtensionSeqBuilder, JunctionThresholdState};
use crate::assembly::globals::is_sbx;
use crate::assembly::indel_builder::{find_indel_extension_reads, set_indel_status};
use crate::assembly::junction_read_types::{calc_proximate_junction_read_ratio, JunctionReadTypes};
use crate::assembly::line_utils::{has_line_source_sequence, is_line_with_local_aligned_insert};
use crate::assembly::read::read::ReadP;
use crate::assembly::ref_base_seq_builder::RefBaseSeqBuilder;
use crate::assembly::ref_genome::RefGenomeInterface;
use crate::assembly::saga::SagaSequenceMatcher;
use crate::assembly::seq_tech_utils::{find_sbx_possible_duplicates, pass_sbx_distinct_prime_positions_filter};
use crate::assembly::sequence_compare::matched_assembly_sequences;
use crate::assembly::types::junction::{JuncP, Junction};
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::ref_side_soft_clip::RefSideSoftClip;
use crate::assembly::types::support_read::SupP;
use crate::assembly::types::SupportType;

/// PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT
pub const MIN_HOTSPOT_JUNCTION_SUPPORT: i32 = 1;

pub struct JunctionAssembler<'a> {
    junction: JuncP,
    saga_matcher: Option<&'a mut SagaSequenceMatcher>,
    ref_genome: Option<&'a dyn RefGenomeInterface>,
    non_junction_reads: Vec<ReadP>,
}

pub fn min_read_threshold(junction: &Junction) -> i32 {
    if junction.hotspot { MIN_HOTSPOT_JUNCTION_SUPPORT } else { ASSEMBLY_MIN_READ_SUPPORT }
}

fn remove_first(list: &mut Vec<ReadP>, read: &ReadP) {
    if let Some(i) = list.iter().position(|x| x == read) {
        list.remove(i);
    }
}

impl<'a> JunctionAssembler<'a> {
    pub fn new(junction: JuncP, ref_genome: Option<&'a dyn RefGenomeInterface>, saga_matcher: Option<&'a mut SagaSequenceMatcher>) -> JunctionAssembler<'a> {
        JunctionAssembler { junction, saga_matcher, ref_genome, non_junction_reads: Vec::new() }
    }
    pub fn new_basic(junction: JuncP) -> JunctionAssembler<'a> {
        Self::new(junction, None, None)
    }
    pub fn non_junction_reads(&self) -> &Vec<ReadP> {
        &self.non_junction_reads
    }

    pub fn process_junction(&mut self, raw_reads: &[ReadP]) -> Vec<AsmP> {
        if !self.junction.discordant_only {
            set_indel_status(self.junction, raw_reads);
        }
        let mut junction_reads: Vec<ReadP> = Vec::new();
        let mut extension_reads: Vec<ReadP> = Vec::new();
        let mut junc_threshold_state = JunctionThresholdState::new();
        junc_threshold_state.min_required_reads = min_read_threshold(&self.junction);
        junc_threshold_state.min_extension_length = ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
        junc_threshold_state.min_second_extension_length = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

        if self.junction.indel_based() {
            find_indel_extension_reads(self.junction, raw_reads, &mut extension_reads, &mut junction_reads, &mut self.non_junction_reads);
            if extension_reads.is_empty() {
                return Vec::new();
            }
            junc_threshold_state.extension_length_valid = true;
            junc_threshold_state.second_extension_length_valid = true;
        } else if self.junction.discordant_only {
            let mut discordant_junction_builder = DiscordantJunctionBuilder::new(self.junction);
            discordant_junction_builder.assess_discordant_junction(raw_reads, &mut extension_reads, &mut junction_reads, &mut self.non_junction_reads);
            if extension_reads.is_empty() {
                return Vec::new();
            }
            if let Some(nj) = discordant_junction_builder.new_junction() {
                self.junction = nj;
            }
            junc_threshold_state.extension_length_valid = true;
            junc_threshold_state.second_extension_length_valid = true;
        } else {
            let junc_read_types = JunctionReadTypes::new(self.junction, raw_reads, self.ref_genome);
            junction_reads.extend(junc_read_types.junction_reads());
            junc_read_types.establish_extension_reads(&mut extension_reads, &mut junc_threshold_state);
            if !junc_threshold_state.extension_length_valid {
                return Vec::new();
            }
            if !junc_threshold_state.second_extension_length_valid {
                if junc_threshold_state.is_line {
                    return Vec::new();
                }
                junc_threshold_state.require_non_extension_support = true;
            }
            self.non_junction_reads.extend(junc_read_types.non_junction_reads().iter().cloned());
        }

        junc_threshold_state.min_reads_valid = extension_reads.len() as i32 >= junc_threshold_state.min_required_reads;
        if !junc_threshold_state.min_reads_valid && !junc_threshold_state.require_non_extension_support {
            return Vec::new();
        }
        let duplicate_long_extension_reads: Vec<ReadP> = if is_sbx() { find_sbx_possible_duplicates(&self.junction, &extension_reads) } else { Vec::new() };
        for x in &duplicate_long_extension_reads {
            remove_first(&mut extension_reads, x);
        }
        for x in &duplicate_long_extension_reads {
            remove_first(&mut junction_reads, x);
        }

        let mut extension_seq_builder = ExtensionSeqBuilder::new(self.junction, &extension_reads);
        extension_seq_builder.check_validity(&mut junc_threshold_state);
        if !junc_threshold_state.extension_length_valid {
            return Vec::new();
        }
        if (!junc_threshold_state.second_extension_length_valid || !junc_threshold_state.min_reads_valid) && !junc_threshold_state.require_non_extension_support {
            return Vec::new();
        }
        let mut assembly_support: Vec<SupP> = extension_seq_builder.form_assembly_support();
        let first_assembly = JunctionAssembly::new(
            self.junction,
            extension_seq_builder.extension_bases().to_vec(),
            extension_seq_builder.base_qualities().to_vec(),
            &assembly_support,
            extension_seq_builder.repeats().clone(),
        );
        if !first_assembly.indel() && has_line_source_sequence(&first_assembly) {
            return Vec::new();
        }
        first_assembly.m().set_ext_base_build_info(extension_seq_builder.build_information());
        if extension_seq_builder.has_line_sequence() {
            first_assembly.m().mark_line_sequence();
        }
        let initial_assembly_support = assembly_support.len();
        self.add_junction_reads(first_assembly, &mut extension_seq_builder, &junction_reads);

        if junc_threshold_state.require_non_extension_support {
            for read in first_assembly.support().clone() {
                if assembly_support.contains(&read) {
                    continue;
                }
                if read.extension_base_mismatches() > 0 {
                    continue;
                }
                let base_matches = read.extension_base_matches();
                let extension_length = read.extension_base_overlap();
                let match_perc = if extension_length > 0 { base_matches as f64 / extension_length as f64 } else { 0.0 };
                if match_perc < ASSEMBLY_MAX_JUNCTION_READ_LOW_QUAL_MISMATCH_PERC {
                    continue;
                }
                if extension_length >= junc_threshold_state.min_second_extension_length {
                    junc_threshold_state.second_extension_length_valid = true;
                    assembly_support.push(read);
                    first_assembly.m().stats_mut().non_extension_support_reads += 1;
                }
            }
            junc_threshold_state.min_reads_valid = assembly_support.len() as i32 >= junc_threshold_state.min_required_reads;
        }
        if !self.meets_min_support_threshold(&assembly_support, junc_threshold_state.min_required_reads) {
            return Vec::new();
        }
        if first_assembly.has_line_sequence() && is_line_with_local_aligned_insert(&first_assembly) {
            first_assembly.m().unmark_line_sequence();
            if first_assembly.extension_length() < ASSEMBLY_MIN_SOFT_CLIP_LENGTH {
                return Vec::new();
            }
        }

        let second_assembly = self.check_second_assembly(&extension_seq_builder.mismatch_reads(), first_assembly, &junction_reads);
        let mut assemblies: Vec<AsmP> = vec![first_assembly];
        if let Some(second_assembly) = second_assembly {
            assemblies.push(second_assembly);
            if !self.keep_second_assembly(second_assembly.support_count() as i32, initial_assembly_support as i32) {
                assemblies.retain(|a| *a != first_assembly);
            }
        }

        for assembly in assemblies.iter().copied() {
            if assemblies.len() > 1 {
                RefSideSoftClip::check_support_vs_ref_side_soft_clip(assembly);
            }
            let ref_base_seq_builder = RefBaseSeqBuilder::new(assembly);
            assembly.m().set_ref_bases(&ref_base_seq_builder);
            assembly.m().build_repeat_info();
            let saga_matched = match self.saga_matcher.as_deref_mut() {
                Some(m) => assembly.m().match_to_saga(m),
                None => false,
            };
            let is_valid_prep_junction = if self.junction.soft_clip_based() {
                self.junction.max_soft_clip_length >= MIN_VARIANT_LENGTH || assembly.has_line_sequence() || self.junction.hotspot
            } else {
                true
            };
            if (junc_threshold_state.uses_lower_saga_limits || !is_valid_prep_junction) && !saga_matched {
                return Vec::new();
            }
            assembly.m().stats_mut().proximate_junc_read_ratio = calc_proximate_junction_read_ratio(&self.junction, raw_reads);
        }
        assemblies
    }

    fn check_second_assembly(&mut self, extension_reads: &[ReadP], first_assembly: AsmP, junction_reads: &[ReadP]) -> Option<AsmP> {
        if extension_reads.is_empty() || self.junction.discordant_only {
            return None;
        }
        if first_assembly.has_line_sequence() {
            return None;
        }
        let second_support = extension_reads.len() as i32;
        let second_support_perc = second_support as f64 / first_assembly.support_count() as f64;
        let mut junc_threshold_state = JunctionThresholdState::new();
        junc_threshold_state.min_required_reads = ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
        junc_threshold_state.min_extension_length = ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
        junc_threshold_state.min_second_extension_length = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
        if second_support < junc_threshold_state.min_required_reads || second_support_perc < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC {
            return None;
        }
        let mut extension_seq_builder = ExtensionSeqBuilder::new(self.junction, extension_reads);
        extension_seq_builder.check_validity(&mut junc_threshold_state);
        if !junc_threshold_state.is_valid() {
            return None;
        }
        let assembly_support = extension_seq_builder.form_assembly_support();
        if !self.keep_second_assembly(first_assembly.support_count() as i32, assembly_support.len() as i32) {
            return None;
        }
        if !self.pass_distinct_fragments_filter(&assembly_support) {
            return None;
        }
        let new_assembly = JunctionAssembly::new(
            self.junction,
            extension_seq_builder.extension_bases().to_vec(),
            extension_seq_builder.base_qualities().to_vec(),
            &assembly_support,
            extension_seq_builder.repeats().clone(),
        );
        if new_assembly.extension_length() < ASSEMBLY_MIN_SOFT_CLIP_LENGTH {
            return None;
        }
        if matched_assembly_sequences(&first_assembly, &new_assembly) {
            return None;
        }
        if has_line_source_sequence(&new_assembly) {
            return None;
        }
        new_assembly.m().set_ext_base_build_info(extension_seq_builder.build_information());
        self.add_junction_reads(new_assembly, &mut extension_seq_builder, junction_reads);
        Some(new_assembly)
    }

    pub fn keep_second_assembly(&self, first_support_count: i32, second_support_count: i32) -> bool {
        let second_support_perc = second_support_count as f64 / first_support_count as f64;
        second_support_count >= ASSEMBLY_SPLIT_MIN_READ_SUPPORT && second_support_perc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC
    }

    fn add_junction_reads(&self, assembly: AsmP, extension_seq_builder: &mut ExtensionSeqBuilder, junction_reads: &[ReadP]) {
        let mut mismatch_read_count = 0;
        for read in junction_reads {
            if assembly.support().iter().any(|x| x.cached_read() == Some(read)) {
                continue;
            }
            match extension_seq_builder.check_add_junction_read(read) {
                None => mismatch_read_count += 1,
                Some(read_info) => {
                    if read_info.mismatched() {
                        mismatch_read_count += 1;
                    } else if extension_seq_builder.sufficient_qual_matches(&read_info) {
                        assembly.m().add_junction_support(read, SupportType::Junction, read_info.start_index(), &read_info);
                    }
                }
            }
        }
        assembly.m().add_mismatch_read_count(mismatch_read_count);
    }

    fn meets_min_support_threshold(&self, support: &[SupP], min_required_read_count: i32) -> bool {
        if !self.above_min_support_threshold(support, min_required_read_count) {
            return false;
        }
        self.pass_distinct_fragments_filter(support)
    }

    fn above_min_support_threshold(&self, support: &[SupP], min_required_read_count: i32) -> bool {
        if support.len() as i32 >= min_required_read_count * 2 {
            return true;
        }
        let unique_read_ids: HashSet<&str> = support.iter().map(|x| x.id()).collect();
        unique_read_ids.len() as i32 >= min_required_read_count
    }

    fn pass_distinct_fragments_filter(&self, support: &[SupP]) -> bool {
        if self.junction.hotspot {
            return true;
        }
        if is_sbx() && !pass_sbx_distinct_prime_positions_filter(support) {
            return false;
        }
        let mut read_positions: HashSet<i32> = HashSet::new();
        let mut read_end_positions: HashSet<i32> = HashSet::new();
        for read in support {
            if read.is_paired_read() {
                read_positions.insert(read.cached_read().expect("support read without cached read").five_prime_fragment_position());
                if read_positions.len() as i32 >= ASSEMBLY_MIN_DISTINCT_FRAGS {
                    return true;
                }
            } else {
                read_positions.insert(read.unclipped_start());
                read_end_positions.insert(read.unclipped_end());
                if read_positions.len() as i32 >= ASSEMBLY_MIN_DISTINCT_FRAGS && read_end_positions.len() as i32 >= ASSEMBLY_MIN_DISTINCT_FRAGS {
                    return true;
                }
            }
        }
        false
    }
}
