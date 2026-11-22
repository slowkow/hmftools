//! JunctionGroupAssembler: loads the reads of a junction group and builds its junction assemblies.

use std::collections::{HashMap, HashSet};

use crate::assembly::alignment::alignment_checker::AlignmentChecker;
use crate::assembly::assembly_deduper::dedup_proximate_assemblies;
use crate::assembly::common_utils::{position_within, positions_overlap};
use crate::assembly::config::AssemblyConfig;
use crate::assembly::constants::*;
use crate::assembly::globals::{dev_debug, is_illumina, is_sbx};
use crate::assembly::indel_builder::calc_indel_inferred_unclipped_positions;
use crate::assembly::junction_assembler::JunctionAssembler;
use crate::assembly::read::bam_reader::BamReader;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_adjustments::{mark_line_soft_clips, trim_low_qual_soft_clip_bases, trim_poly_g_sequences};
use crate::assembly::read::read_stats::ReadStats;
use crate::assembly::ref_base_extender::find_assembly_candidate_extensions;
use crate::assembly::remote_region_finder::filter_low_qual_read;
use crate::assembly::saga::SagaSequenceMatcher;
use crate::assembly::sam::*;
use crate::assembly::seq_tech_utils::{find_illumina_adapter_trims, trim_illumina_adapter_bases, trim_sbx_uncertain_bases};
use crate::common::compat;
use crate::assembly::mem_stats;
use crate::assembly::types::junction_assembly::{start_tracking_created, take_created, AsmP};
use crate::assembly::types::junction_group::JuncGroupP;

pub struct JunctionGroupAssembler<'a> {
    config: &'a AssemblyConfig,
    bam_reader: BamReader,
    alignment_checker: AlignmentChecker,
    read_stats: ReadStats,
    decoy_assemblies: Vec<AsmP>,
    saga_matcher: Option<SagaSequenceMatcher>,
}

/// reads of one fragment within a junction group
struct ReadGroup {
    reads: Vec<ReadP>,
}

impl ReadGroup {
    fn form_read_links(&self) {
        if self.reads.len() == 1 {
            return;
        }
        for i in 0..self.reads.len() - 1 {
            for j in i + 1..self.reads.len() {
                Read::make_read_links(&self.reads[i], &self.reads[j]);
            }
        }
    }
}

impl<'a> JunctionGroupAssembler<'a> {
    pub fn new(config: &'a AssemblyConfig, bam_reader: BamReader, alignment_checker: AlignmentChecker, saga_matcher: Option<SagaSequenceMatcher>) -> JunctionGroupAssembler<'a> {
        JunctionGroupAssembler { config, bam_reader, alignment_checker, read_stats: ReadStats::default(), decoy_assemblies: Vec::new(), saga_matcher }
    }
    pub fn read_stats(&self) -> &ReadStats {
        &self.read_stats
    }
    pub fn decoy_assemblies(&self) -> &Vec<AsmP> {
        &self.decoy_assemblies
    }
    pub fn take_decoy_assemblies(&mut self) -> Vec<AsmP> {
        std::mem::take(&mut self.decoy_assemblies)
    }
    /// The decoy match rows formed since the last call (for the junction group just processed).
    pub fn take_decoy_match_rows(&mut self) -> Vec<String> {
        self.alignment_checker.take_pending_rows()
    }

    pub fn process_junction_group(&mut self, junction_group: JuncGroupP) {
        let slice_start = junction_group.read_range_start();
        let slice_end = junction_group.read_range_end();
        let mut read_group_map: HashMap<String, ReadGroup> = HashMap::new();
        let mut supplementary_repeats: HashMap<String, SamRecord> = HashMap::new();
        let config = self.config;
        let read_stats = &mut self.read_stats;
        {
            self.bam_reader.slice_bam(config, junction_group.chromosome(), slice_start, slice_end, |record, is_reference| {
                process_record(record, is_reference, junction_group, read_stats, &mut read_group_map, &mut supplementary_repeats);
            });
            for group in read_group_map.values() {
                group.form_read_links();
            }
            drop(read_group_map);
            drop(supplementary_repeats);
            if is_illumina() && !compat::java("A12") {
                finish_reads_after_linking(junction_group, read_stats);
            }
        }

        start_tracking_created();
        let mut junction_group_assemblies: Vec<AsmP> = Vec::new();
        let mut deduped_indels: Vec<AsmP> = Vec::new();
        for junction in junction_group.junctions().iter().copied() {
            let (junction_boundary_start, junction_boundary_end) =
                if junction.is_forward() { (junction.position - BAM_READ_JUNCTION_BUFFER, junction.position) } else { (junction.position, junction.position + BAM_READ_JUNCTION_BUFFER) };
            let candidate_reads: Vec<ReadP> = junction_group
                .candidate_reads()
                .iter()
                .filter(|x| positions_overlap(junction_boundary_start, junction_boundary_end, x.unclipped_start(), x.unclipped_end()))
                .cloned()
                .collect();
            if candidate_reads.is_empty() {
                continue;
            }
            if junction.discordant_only && config.discordant_only_disabled {
                continue;
            }
            let mut junction_assembler = JunctionAssembler::new(junction, config.ref_genome(), self.saga_matcher.as_mut());
            let mut candidate_assemblies = junction_assembler.process_junction(&candidate_reads);
            dedup_proximate_assemblies(&mut junction_group_assemblies, &mut candidate_assemblies, &mut deduped_indels);
            for assembly in candidate_assemblies {
                if !assembly.is_saga_matched() && self.alignment_checker.matches_decoy(&assembly) {
                    self.read_stats.decoy_sequences += 1;
                    self.decoy_assemblies.push(assembly);
                    continue;
                }
                if self.alignment_checker.fails_mappability(&assembly) {
                    self.read_stats.ref_base_alignment_fails += 1;
                    continue;
                }
                find_assembly_candidate_extensions(assembly, junction_assembler.non_junction_reads());
                junction_group_assemblies.push(assembly);
            }
        }
        junction_group.m().add_junction_assemblies(&junction_group_assemblies);
        for x in &junction_group_assemblies {
            x.m().clear_support_cached_reads();
            x.m().set_read_source_group(junction_group);
        }
        junction_group.m().set_assemblies_holding_reads(junction_group_assemblies.len());
        // assemblies that were built and then rejected or deduplicated are not referred to again: empty them
        // so the reads they hold can be freed (the decoy assemblies are kept for the assembly TSV)
        let kept: HashSet<AsmP> = junction_group_assemblies.iter().chain(self.decoy_assemblies.iter()).copied().collect();
        for assembly in take_created() {
            if !kept.contains(&assembly) {
                assembly.m().release_discarded();
            }
        }
        // Java clears the list before the cached reads of the assemblies; the order makes no difference to
        // the results and clearing it last lets the reads only used as junction support be freed here
        junction_group.m().clear_candidate_reads();
        mem_stats::READS_RETAINED.fetch_add(junction_group.retained_read_count() as i64, std::sync::atomic::Ordering::Relaxed);
    }
}

fn process_record(
    mut record: SamRecord, is_reference: bool, junction_group: JuncGroupP, read_stats: &mut ReadStats, read_group_map: &mut HashMap<String, ReadGroup>,
    supplementary_repeats: &mut HashMap<String, SamRecord>,
) {
    let java_a12 = compat::java("A12");
    if ignore_identical_supplementary(&record, read_stats, supplementary_repeats) {
        return;
    }
    if !record.read_unmapped() && !record.cigar.iter().any(|x| x.op == CigarOp::M) {
        return;
    }
    if filter_low_qual_read(&mut record) {
        read_stats.low_base_qual_filtered += 1;
        return;
    }
    let read = Read::new(record);
    read_stats.total_reads += 1;
    if is_reference {
        read.m().mark_reference();
    }
    if is_illumina() && trim_poly_g_sequences(read.m()) {
        read_stats.poly_g_trimmed += 1;
    }
    mark_line_soft_clips(read.m());
    if is_illumina() {
        if trim_low_qual_soft_clip_bases(read.m()) {
            read_stats.low_base_qual_trimmed += 1;
        }
        if java_a12 {
            // Java: runs before mates are linked, so it has no effect. Corrected: done for the whole junction
            // group once the links exist (finish_reads_after_linking), with the steps that follow it here
            trim_illumina_adapter_bases(&read);
        }
    } else if is_sbx() && trim_sbx_uncertain_bases(read.m()) {
        read_stats.low_base_qual_trimmed += 1;
    }
    if read.trim_count() > 0 && !has_valid_cigar(read.cigar_elements()) {
        read_stats.low_base_qual_filtered += 1;
        return;
    }
    if (java_a12 || !is_illumina()) && calc_indel_inferred_unclipped_positions(read.m()) {
        read_stats.indel_soft_clip_converted += 1;
    }
    junction_group.m().add_candidate_read(read.clone());
    if let Some(read_group) = read_group_map.get_mut(read.id()) {
        read_group.reads.push(read);
        return;
    }
    let has_local_mate = read.is_mate_unmapped()
        || (read.is_mate_mapped()
            && read.mate_chromosome() == Some(read.chromosome())
            && positions_overlap(read.mate_alignment_start(), read.mate_alignment_end(), junction_group.read_range_start(), junction_group.read_range_end()));
    let has_local_supplementary = read
        .supplementary_data()
        .is_some_and(|s| s.chromosome == read.chromosome() && position_within(s.position, junction_group.read_range_start(), junction_group.read_range_end()));
    if !has_local_mate && !has_local_supplementary {
        return;
    }
    let mut read_group = ReadGroup { reads: vec![read.clone()] };
    if read.is_unmapped() {
        for cached_read in junction_group.candidate_reads() {
            if *cached_read != read && cached_read.is_mate_mapped() && cached_read.id() == read.id() && !read_group.reads.contains(cached_read) {
                read_group.reads.push(cached_read.clone());
                break;
            }
        }
    }
    read_group_map.insert(read.id().to_string(), read_group);
}

/// The corrected form of finding A12, for Illumina reads: the adapter trim needs the mate, so it is done here,
/// after the reads of the junction group are linked, followed by the two steps that Java has after it in
/// processRecord: a read left without a valid cigar is dropped, and the indel-inferred unclipped positions
/// are worked out on the final cigar.
fn finish_reads_after_linking(junction_group: JuncGroupP, read_stats: &mut ReadStats) {
    let trims = find_illumina_adapter_trims(junction_group.candidate_reads());
    let mut dropped: Vec<ReadP> = Vec::new();
    for (index, count, from_start) in trims {
        let read = &junction_group.candidate_reads()[index];
        read.m().trim_bases(count, from_start);
        read_stats.adapter_trimmed += 1;
        if !has_valid_cigar(read.cigar_elements()) {
            read_stats.low_base_qual_filtered += 1;
            dropped.push(read.clone());
        }
    }
    if !dropped.is_empty() {
        for read in &dropped {
            // as if the read had not been loaded: its mate has no link to it
            if let Some(mate_read) = read.mate_read().cloned() {
                if mate_read.mate_read() == Some(read) {
                    mate_read.m().clear_mate_read();
                }
            }
            read.m().clear_mate_read();
        }
        junction_group.m().remove_candidate_reads(&dropped);
    }
    for read in junction_group.candidate_reads() {
        if calc_indel_inferred_unclipped_positions(read.m()) {
            read_stats.indel_soft_clip_converted += 1;
        }
    }
}

fn ignore_identical_supplementary(read: &SamRecord, read_stats: &mut ReadStats, supplementary_repeats: &mut HashMap<String, SamRecord>) -> bool {
    if !dev_debug() {
        return false;
    }
    if !read.supplementary() {
        return false;
    }
    let Some(previous_read) = supplementary_repeats.get(&read.read_name) else {
        supplementary_repeats.insert(read.read_name.clone(), read.clone());
        return false;
    };
    if previous_read.alignment_start != read.alignment_start || previous_read.flags != read.flags {
        return false;
    }
    if previous_read.cigar_string() != read.cigar_string() {
        return false;
    }
    read_stats.identical_supplementaries += 1;
    supplementary_repeats.remove(&read.read_name);
    true
}
