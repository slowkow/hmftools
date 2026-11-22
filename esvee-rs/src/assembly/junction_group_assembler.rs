//! JunctionGroupAssembler: loads the reads of a junction group and builds its junction assemblies.

use std::collections::HashMap;

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
use crate::assembly::seq_tech_utils::{trim_illumina_adapter_bases, trim_sbx_uncertain_bases};
use crate::assembly::types::junction_assembly::AsmP;
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
        self.bam_reader.slice_bam(config, junction_group.chromosome(), slice_start, slice_end, |record, is_reference| {
            process_record(record, is_reference, junction_group, read_stats, &mut read_group_map, &mut supplementary_repeats);
        });
        for group in read_group_map.values() {
            group.form_read_links();
        }
        drop(read_group_map);

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
        junction_group.m().clear_candidate_reads();
        for x in &junction_group_assemblies {
            x.m().clear_support_cached_reads();
        }
    }
}

fn process_record(
    mut record: SamRecord, is_reference: bool, junction_group: JuncGroupP, read_stats: &mut ReadStats, read_group_map: &mut HashMap<String, ReadGroup>,
    supplementary_repeats: &mut HashMap<String, SamRecord>,
) {
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
        // as in Java this runs before mates are linked, so it has no effect
        trim_illumina_adapter_bases(&read);
    } else if is_sbx() && trim_sbx_uncertain_bases(read.m()) {
        read_stats.low_base_qual_trimmed += 1;
    }
    if read.trim_count() > 0 && !has_valid_cigar(read.cigar_elements()) {
        read_stats.low_base_qual_filtered += 1;
        return;
    }
    if calc_indel_inferred_unclipped_positions(read.m()) {
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
