//! AssemblyAligner and Alignment: aligns each assembly alignment's full sequence to the reference, forms
//! breakends and allocates read support to them.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::common::compat;
use crate::assembly::alignment::align_data::{AlignData, AlignP};
use crate::assembly::alignment::alignment_fragments::AlignmentFragments;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::alignment::breakend_builder::BreakendBuilder;
use crate::assembly::alignment::bwa::IBwaMemAligner;
use crate::assembly::config::{AssemblyConfig, WriteType};
use crate::assembly::constants::*;
use crate::assembly::globals::{is_illumina, sample_discordant_rate};
use crate::assembly::orientation::Orientation;
use crate::assembly::output::alignment_writer::AlignmentWriter;
use crate::assembly::perf::PerfCounter;
use crate::assembly::ptr::P;
use crate::assembly::saga::{SagaJunctionInfo, SagaMatcherFactory, SagaSequenceMatcher};
use crate::assembly::seq_tech_utils::{SBX_MEDIUM_QUAL_DESYNC_COUNT, SBX_PRIME_POSITION_RANGE_THRESHOLD};
use crate::assembly::types::repeat_info::find_repeats_min;
use crate::assembly::types::support_read::SupP;
use crate::assembly::types::{StructuralVariantType, SupportType};

/// The alignments kept for an assembly alignment, and the supplementaries replaced by a requery.
#[derive(Default)]
pub struct AlignmentInfo {
    pub alignments: Vec<AlignP>,
    pub requeried_alignments: Vec<AlignP>,
}

pub struct AssemblyAligner<'a> {
    config: &'a AssemblyConfig,
    aligner: &'a dyn IBwaMemAligner,
    requeried_supp_count: i32,
    requeried_soft_clip_count: i32,
    saga_matcher: Option<SagaSequenceMatcher>,
}

impl<'a> AssemblyAligner<'a> {
    pub fn new(config: &'a AssemblyConfig, aligner: &'a dyn IBwaMemAligner) -> Self {
        Self::new_saga(config, aligner, None)
    }
    pub fn new_saga(config: &'a AssemblyConfig, aligner: &'a dyn IBwaMemAligner, saga_matcher_factory: Option<&SagaMatcherFactory>) -> Self {
        AssemblyAligner { config, aligner, requeried_supp_count: 0, requeried_soft_clip_count: 0, saga_matcher: saga_matcher_factory.map(|f| f.create_sequence_matcher()) }
    }
    pub fn requeried_supp_count(&self) -> i32 {
        self.requeried_supp_count
    }
    pub fn requeried_soft_clip_count(&self) -> i32 {
        self.requeried_soft_clip_count
    }

    /// Returns None where Java writes nothing (invalid), otherwise the alignments to write.
    pub fn process_assembly(&mut self, assembly_alignment: AsmAlignP) -> Option<AlignmentInfo> {
        if !assembly_alignment.is_valid() {
            crate::warn!("assembly alignment({}) invalid, skipping", *assembly_alignment);
            return None;
        }
        if assembly_alignment.is_merged() {
            return Some(AlignmentInfo::default());
        }
        // try to form breakends from a matched SAGA variant first, otherwise from alignment to the ref genome
        let mut alignment_info: Option<AlignmentInfo> = None;
        if self.try_match_assembly_alignment_to_saga(assembly_alignment) {
            alignment_info = self.process_saga_matched_assembly(assembly_alignment);
        }
        Some(match alignment_info {
            Some(info) => info,
            None => self.align_assembly(assembly_alignment),
        })
    }

    fn try_match_assembly_alignment_to_saga(&mut self, assembly_alignment: AsmAlignP) -> bool {
        let Some(saga_matcher) = self.saga_matcher.as_mut() else { return false };
        // SAGA is only used for the simple case of 0 or 1 links
        let is_chained = assembly_alignment.phase_set().is_some_and(|p| p.assembly_links().len() > 1);
        if is_chained {
            return false;
        }
        let junction_infos: Vec<SagaJunctionInfo> = assembly_alignment.link_indices().iter().map(|i| SagaJunctionInfo { offset: *i }).collect();
        // a LINE site has a lower extension sequence requirement, so the junction overlap requirement is lowered too
        let lower_junction_overlap = assembly_alignment.assemblies().iter().any(|x| x.has_line_sequence());
        let saga_match = saga_matcher.match_by_sequence(assembly_alignment.full_sequence().unwrap_or("").as_bytes(), &junction_infos, lower_junction_overlap, false);
        let matched = saga_match.is_some();
        assembly_alignment.m().set_saga_match(saga_match);
        matched
    }

    fn process_saga_matched_assembly(&mut self, assembly_alignment: AsmAlignP) -> Option<AlignmentInfo> {
        let saga_alignment = assembly_alignment.saga_match()?.alignment.clone();
        let align_data = P::new(AlignData::from_saga(&saga_alignment));
        let full_sequence = assembly_alignment.full_sequence().unwrap_or("").to_string();
        align_data.m().set_full_sequence_data(&full_sequence, assembly_alignment.full_sequence_length());
        // a single alignment, so there is no neighbouring-alignment overlap to account for
        align_data.m().set_adjusted_alignment(&full_sequence, 0, 0);
        let breakend_builder = BreakendBuilder::new(self.config.ref_genome(), assembly_alignment);
        if !breakend_builder.form_breakends_from_saga(&saga_alignment, align_data) {
            return None;
        }
        Some(AlignmentInfo { alignments: vec![align_data], requeried_alignments: Vec::new() })
    }

    /// The second half of Java's processAssembly. Support reads are shared by assembly alignments whose phase
    /// sets have secondary links to the same assemblies, and each allocation both reads and overwrites their
    /// fragment lengths, so this is run in assembly alignment order (the single-threaded Java order).
    pub fn allocate_support(config: &AssemblyConfig, assembly_alignment: AsmAlignP, combined_sample_ids: &[String]) {
        let _ = config;
        let alignment_fragments = AlignmentFragments::new(assembly_alignment, combined_sample_ids);
        alignment_fragments.allocate_breakend_support();
    }

    fn align(&self, bases: &str) -> Vec<AlignP> {
        self.aligner.align_sequence(bases.as_bytes()).iter().filter_map(|x| AlignData::from_ref(x, self.config.ref_gen_version_v38)).map(P::new).collect()
    }

    fn align_assembly(&mut self, assembly_alignment: AsmAlignP) -> AlignmentInfo {
        let full_sequence = assembly_alignment.full_sequence().unwrap_or("").to_string();
        let mut alignments = self.align(&full_sequence);
        for x in &alignments {
            x.m().set_full_sequence_data(&full_sequence, assembly_alignment.full_sequence_length());
        }
        let mut requeried_alignments: Vec<AlignP> = Vec::new();
        alignments = self.requery_soft_clip_alignments(assembly_alignment, alignments);
        alignments = self.requery_supplementary_alignments(assembly_alignment, alignments, &mut requeried_alignments);
        if !alignments.is_empty() {
            let breakend_builder = BreakendBuilder::new(self.config.ref_genome(), assembly_alignment);
            breakend_builder.form_breakends(&alignments);
            if is_weak_single_read_extension_assembly(assembly_alignment) {
                assembly_alignment.m().breakends_mut().clear();
            }
        }
        AlignmentInfo { alignments, requeried_alignments }
    }

    fn requery_soft_clip_alignments(&mut self, assembly_alignment: AsmAlignP, alignments: Vec<AlignP>) -> Vec<AlignP> {
        let full_length = assembly_alignment.full_sequence_length();
        let min_seq_start = alignments.iter().map(|x| x.sequence_start()).min().unwrap_or(0);
        let max_seq_end = alignments.iter().map(|x| x.sequence_end()).max().unwrap_or(full_length);
        let missing_bases_on_end = full_length - max_seq_end - 1;
        let relevant_alignment: AlignP;
        let is_left_clip: bool;
        let first_bases_missing: bool;
        // Java dereferences a null alignment here if the list is empty (not reachable: both tests fail then)
        if min_seq_start >= ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH && min_seq_start >= missing_bases_on_end {
            relevant_alignment = *alignments.iter().find(|x| x.sequence_start() == min_seq_start).expect("alignment");
            is_left_clip = relevant_alignment.orientation().is_forward();
            first_bases_missing = true;
        } else if missing_bases_on_end >= ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH {
            relevant_alignment = *alignments.iter().find(|x| x.sequence_end() == max_seq_end).expect("alignment");
            is_left_clip = relevant_alignment.orientation().is_reverse();
            first_bases_missing = false;
        } else {
            return alignments;
        }
        let soft_clip_length = if is_left_clip { relevant_alignment.left_soft_clip_length() } else { relevant_alignment.right_soft_clip_length() };
        let full_sequence = assembly_alignment.full_sequence().unwrap_or("");
        let soft_clip_bases = if first_bases_missing { &full_sequence[..soft_clip_length as usize] } else { &full_sequence[(full_length - soft_clip_length) as usize..] };
        self.requeried_soft_clip_count += 1;
        let mut new_alignments = self.align(soft_clip_bases);
        if new_alignments.is_empty() {
            return alignments;
        }
        let soft_clip_seq_index_start = if is_left_clip == relevant_alignment.orientation().is_forward() {
            (relevant_alignment.sequence_start() - soft_clip_length).max(0)
        } else {
            relevant_alignment.sequence_end() + 1
        };
        for rq_alignment in &new_alignments {
            let (adj_sequence_start, adj_sequence_end);
            if rq_alignment.orientation().is_forward() {
                adj_sequence_start = soft_clip_seq_index_start + rq_alignment.sequence_start();
                adj_sequence_end = soft_clip_seq_index_start + rq_alignment.sequence_end();
            } else {
                let new_sequence_start = (soft_clip_length - 1) - (rq_alignment.raw_sequence_end() - 1);
                let new_sequence_end = new_sequence_start + rq_alignment.segment_length() - 1;
                adj_sequence_start = soft_clip_seq_index_start + new_sequence_start;
                adj_sequence_end = soft_clip_seq_index_start + new_sequence_end;
            }
            rq_alignment.m().set_requeried_sequence_coords(adj_sequence_start, adj_sequence_end);
        }
        if is_left_clip {
            new_alignments.extend(alignments);
            new_alignments
        } else {
            let mut combined = alignments;
            combined.extend(new_alignments);
            combined
        }
    }

    fn requery_supplementary_alignments(&mut self, assembly_alignment: AsmAlignP, alignments: Vec<AlignP>, requeried_alignments: &mut Vec<AlignP>) -> Vec<AlignP> {
        if !alignments.iter().any(|x| x.is_supplementary()) || alignments.iter().all(|x| x.is_supplementary()) {
            return alignments;
        }
        let mut new_alignments: Vec<AlignP> = Vec::new();
        for align_data in alignments {
            if !align_data.is_supplementary() {
                new_alignments.push(align_data);
                continue;
            }
            requeried_alignments.push(align_data);
            new_alignments.extend(self.requery_supplementary_alignment(assembly_alignment, align_data));
        }
        new_alignments
    }

    fn requery_supplementary_alignment(&mut self, assembly_alignment: AsmAlignP, align_data: AlignP) -> Vec<AlignP> {
        self.requeried_supp_count += 1;
        let full_sequence = assembly_alignment.full_sequence().unwrap_or("");
        align_data.m().set_full_sequence_data(full_sequence, assembly_alignment.full_sequence_length());
        let alignment_sequence = &full_sequence[align_data.sequence_start() as usize..(align_data.sequence_end() + 1) as usize];
        let requery_alignments = self.align(alignment_sequence);
        let mut converted_alignments: Vec<AlignP> = Vec::new();
        for rq_alignment in requery_alignments {
            rq_alignment.m().set_full_sequence_data(alignment_sequence, alignment_sequence.len() as i32);
            let converted_alignment = P::new(AlignData::new(
                rq_alignment.chromosome(),
                rq_alignment.position_start(),
                rq_alignment.position_end(),
                rq_alignment.raw_sequence_start(),
                rq_alignment.raw_sequence_end(),
                rq_alignment.map_qual(),
                rq_alignment.score(),
                rq_alignment.flags(),
                rq_alignment.cigar(),
                rq_alignment.n_matches(),
                Some(rq_alignment.xa_tag()),
                Some(rq_alignment.md_tag()),
            ));
            let rq_seq_offset_start = rq_alignment.sequence_start();
            let adj_sequence_start = align_data.sequence_start() + rq_seq_offset_start;
            let rq_seq_offset_end = alignment_sequence.len() as i32 - 1 - rq_alignment.sequence_end();
            let adj_sequence_end = align_data.sequence_end() - rq_seq_offset_end;
            converted_alignment.m().set_requeried_sequence_coords(adj_sequence_start, adj_sequence_end);
            let left = converted_alignment.left_soft_clip_length() + align_data.left_soft_clip_length();
            let right = converted_alignment.right_soft_clip_length() + align_data.right_soft_clip_length();
            converted_alignment.m().set_soft_clip_lengths(left, right);
            converted_alignments.push(converted_alignment);
            align_data.m().mark_dropped_on_requery();
        }
        converted_alignments
    }
}

fn is_weak_single_read_extension_assembly(assembly_alignment: AsmAlignP) -> bool {
    if assembly_alignment.breakends().is_empty() {
        return false;
    }
    let assemblies = assembly_alignment.assemblies();
    if assemblies.len() != 1 && !assembly_alignment.is_local_indel_assembly() {
        return false;
    }
    if assembly_alignment.breakends().len() == 2 {
        let breakend = assembly_alignment.breakends()[0];
        if breakend.is_short_local_del_dup_ins() {
            return false;
        }
    }
    let assembly = assemblies[0];
    if assembly.stats().max_ext_base_match_count < ASSEMBLY_UNLINKED_WEAK_ASSEMBLY_EXTENSION_LENGTH {
        return false;
    }
    if assembly.has_line_sequence() || assembly.junction().discordant_only || assembly.junction().indel_based() {
        return false;
    }
    let has_suspect_extension = has_suspect_extension(assembly_alignment);
    if is_illumina() {
        if assembly.stats().soft_clip_second_max_length < ASSEMBLY_MIN_SOFT_CLIP_LENGTH {
            return true;
        }
        if !has_suspect_extension {
            return false;
        }
    }
    has_longer_minority_extensions(assembly.support(), assembly.junction().orient, has_suspect_extension)
}

pub fn has_suspect_extension(assembly_alignment: AsmAlignP) -> bool {
    let breakend = assembly_alignment.breakends()[0];
    if breakend.sv_type() != StructuralVariantType::Bnd && breakend.sv_length() < WEAK_ASSEMBLY_LONG_LENGTH {
        return false;
    }
    if breakend.homology.length() >= INV_SHORT_MAX_HOMOLOGY_HIGHER && sample_discordant_rate() > WEAK_ASSEMBLY_DISC_RATE {
        return true;
    }
    if breakend.inserted_bases.len() as i32 >= WEAK_ASSEMBLY_MIN_INS_LENGTH {
        return true;
    }
    // Java: the loop over the assemblies returns on its first one, so only that one is examined
    let examined = if compat::java("A26") { 1 } else { usize::MAX };
    for assembly in assembly_alignment.assemblies().iter().take(examined) {
        let ref_bases = assembly.form_ref_base_sequence();
        let min_ref_length = (WEAK_ASSEMBLY_MIN_REF_REPEAT * WEAK_ASSEMBLY_MAX_REF_REPEAT_LENGTH) as usize;
        let adjacent_ref_bases: &str = if ref_bases.len() <= min_ref_length {
            &ref_bases
        } else if assembly.is_forward_junction() {
            &ref_bases[ref_bases.len() - min_ref_length..]
        } else {
            &ref_bases[..min_ref_length]
        };
        let ref_repeats = find_repeats_min(adjacent_ref_bases.as_bytes(), WEAK_ASSEMBLY_MIN_REF_REPEAT);
        let adjacent_repeat = if assembly.is_forward_junction() {
            ref_repeats.iter().any(|x| x.last_index() >= adjacent_ref_bases.len() as i32 - x.repeat_length())
        } else {
            ref_repeats.iter().any(|x| x.index == 0)
        };
        if adjacent_repeat {
            return true;
        }
    }
    false
}

pub fn has_longer_minority_extensions(reads: &[SupP], assembly_orientation: Orientation, use_increased_limits: bool) -> bool {
    let mut extension_lengths: Vec<i32> = Vec::with_capacity(reads.len());
    for read in reads {
        if read.type_() != SupportType::Junction {
            continue;
        }
        extension_lengths.push(read.junction_extension_length(assembly_orientation));
    }
    extension_lengths.sort_by(|a, b| b.cmp(a));
    let mut min_long_read_count = 1;
    if extension_lengths.len() > 2 {
        let second_longest_length = extension_lengths[1];
        let top_reads: Vec<&SupP> =
            reads.iter().filter(|x| x.type_() == SupportType::Junction && x.junction_extension_length(assembly_orientation) >= second_longest_length).collect();
        let read1 = top_reads[0];
        let read2 = top_reads[1];
        if reads_close_match(read1, read2) && (read1.medium_qual_count() >= SBX_MEDIUM_QUAL_DESYNC_COUNT || read2.medium_qual_count() >= SBX_MEDIUM_QUAL_DESYNC_COUNT) {
            min_long_read_count = 2;
        }
    }
    let mut max_long_read_count: f64 = (reads.len() as f64 / WEAK_ASSEMBLY_UNPAIRED_READ_FACTOR).max(min_long_read_count as f64);
    if use_increased_limits {
        max_long_read_count += 1.0;
    }
    max_long_read_count = max_long_read_count.min(WEAK_ASSEMBLY_UNPAIRED_MAX_READS as f64);
    let upper_factor = WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR;
    let lower_factor = if use_increased_limits { WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR_LOWER } else { upper_factor };
    for i in 0..extension_lengths.len().saturating_sub(1) {
        let ext_length = extension_lengths[i];
        let next_length = extension_lengths[i + 1];
        let read_count = (i + 1) as f64;
        if read_count > max_long_read_count {
            return false;
        }
        let long_ext_factor = if read_count == max_long_read_count { upper_factor } else { lower_factor };
        if ext_length as f64 > next_length as f64 * long_ext_factor {
            return true;
        }
    }
    false
}

fn reads_close_match(read1: &SupP, read2: &SupP) -> bool {
    if read1.orientation() != read2.orientation() {
        return false;
    }
    let three_prime1 = if read1.orientation().is_forward() { read1.untrimmed_end() } else { read1.untrimmed_start() };
    let three_prime2 = if read2.orientation().is_forward() { read2.untrimmed_end() } else { read2.untrimmed_start() };
    (three_prime1 - three_prime2).abs() <= SBX_PRIME_POSITION_RANGE_THRESHOLD
}

/// Alignment.run: aligns all assembly alignments on the configured number of threads. The phased assembly and
/// alignment TSV rows are written in assembly alignment order, which is the order of a single-threaded Java run.
pub fn run_alignment(config: &AssemblyConfig, aligner: &dyn IBwaMemAligner, saga_matcher_factory: Option<&SagaMatcherFactory>, assembly_alignments: &[AsmAlignP]) -> PerfCounter {
    let single_assemblies = assembly_alignments.iter().filter(|x| x.assemblies().len() == 1).count();
    let linked_assemblies = assembly_alignments.len() - single_assemblies;
    crate::info!("running alignment for {} assemblies, linked({}) single({})", assembly_alignments.len(), linked_assemblies, single_assemblies);
    let task_count = (config.threads.max(0) as usize).min(assembly_alignments.len());
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<(usize, AlignmentInfo)>> = Mutex::new(Vec::with_capacity(assembly_alignments.len()));
    let counts: Mutex<(i32, i32)> = Mutex::new((0, 0));
    let perf: Mutex<PerfCounter> = Mutex::new(PerfCounter::new("AssemblerAlignment"));
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..task_count {
            let (next, results, counts, perf) = (&next, &results, &counts, &perf);
            handles.push(scope.spawn(move || {
                let mut task = AssemblyAligner::new_saga(config, aligner, saga_matcher_factory);
                let mut local: Vec<(usize, AlignmentInfo)> = Vec::new();
                let mut perf_counter = PerfCounter::new("AssemblerAlignment");
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= assembly_alignments.len() {
                        break;
                    }
                    let start = std::time::Instant::now();
                    if let Some(info) = task.process_assembly(assembly_alignments[i]) {
                        local.push((i, info));
                    }
                    perf_counter.stop_check_log(start, || format!("alignment count({}) assemblies({})", assembly_alignments[i].assemblies().len(), assembly_alignments[i].info()));
                }
                perf.lock().unwrap().merge(&perf_counter);
                results.lock().unwrap().extend(local);
                let mut c = counts.lock().unwrap();
                c.0 += task.requeried_supp_count();
                c.1 += task.requeried_soft_clip_count();
            }));
        }
        for h in handles {
            if h.join().is_err() {
                std::process::exit(1);
            }
        }
    });
    let (requeried_supp_count, requeried_soft_clip_count) = counts.into_inner().unwrap();
    if requeried_supp_count > 0 || requeried_soft_clip_count > 0 {
        crate::debug!("requeried supp alignments({}) soft-clips({})", requeried_supp_count, requeried_soft_clip_count);
    }
    let mut results = results.into_inner().unwrap();
    results.sort_by_key(|x| x.0);
    let combined_sample_ids = config.combined_sample_ids();
    for (i, _) in &results {
        if !assembly_alignments[*i].is_merged() {
            AssemblyAligner::allocate_support(config, assembly_alignments[*i], &combined_sample_ids);
        }
    }
    crate::info!("alignment complete");
    let mut writer = AlignmentWriter::new(config);
    for (i, info) in &results {
        write_assembly_data(&mut writer, config, assembly_alignments[*i], &info.alignments, &info.requeried_alignments);
    }
    writer.close();
    perf.into_inner().unwrap()
}

/// Alignment.writeAssemblyData
pub fn write_assembly_data(writer: &mut AlignmentWriter, config: &AssemblyConfig, assembly_alignment: AsmAlignP, alignments: &[AlignP], requeried_alignments: &[AlignP]) {
    if !assembly_alignment.is_valid() {
        return;
    }
    if config.write_types.contains(&WriteType::PhasedAssembly) {
        writer.write_phased_assembly(assembly_alignment);
    }
    if config.write_types.contains(&WriteType::Alignment) {
        let mut alignments_to_write: Vec<AlignP> = alignments.to_vec();
        alignments_to_write.extend(requeried_alignments.iter().copied());
        writer.write_alignment_details(assembly_alignment, &alignments_to_write);
    }
}
