//! JunctionAssembly: the consensus sequence built from the reads at one junction, with its support.

use std::collections::HashSet;

use crate::assembly::assembly_utils::{calc_trimmed_ref_base_length, read_qual_from_junction};
use crate::assembly::common_utils::{above_min_qual, below_min_qual, IndelCoords};
use crate::assembly::constants::*;
use crate::assembly::globals::clear_cached_read;
use crate::assembly::indel_builder::{converted_indel_crosses_junction, find_inserted_bases};
use crate::assembly::ptr::{A, P};
use crate::assembly::read::read::{subset_array, ReadP};
use crate::assembly::read_parse_state::ReadParseState;
use crate::assembly::ref_base_seq_builder::RefBaseSeqBuilder;
use crate::assembly::saga::{SagaJunctionInfo, SagaSequenceMatch, SagaSequenceMatcher};
use crate::assembly::sam::{cigar_elements_to_str, CigarElement};
use crate::assembly::types::assembly_stats::AssemblyStats;
use crate::assembly::types::junction_group::JuncGroupP;
use crate::assembly::types::junction::JuncP;
use crate::assembly::types::phase_group::PhaseGroupP;
use crate::assembly::types::phase_set::PhaseSetP;
use crate::assembly::types::read_assembly_indices::ReadAssemblyIndices;
use crate::assembly::types::ref_side_soft_clip::{RefSideSoftClip, RefSideSoftClipP};
use crate::assembly::types::remote_region::RemoteRegionP;
use crate::assembly::types::repeat_info::{build_trimmed_ref_base_sequence, find_repeats, RepeatInfo};
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::{AssemblyOutcome, SupportType};
use crate::assembly::ref_genome::RefGenomeInterface;

pub type AsmP = P<JunctionAssembly>;

thread_local! {
    /// assemblies created on this thread since `start_tracking_created`, see `take_created`
    static CREATED: std::cell::RefCell<Option<Vec<AsmP>>> = const { std::cell::RefCell::new(None) };
}

fn new_tracked(assembly: JunctionAssembly) -> AsmP {
    let assembly = P::new(assembly);
    crate::assembly::mem_stats::ASSEMBLIES_CREATED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    CREATED.with(|c| {
        if let Some(created) = c.borrow_mut().as_mut() {
            created.push(assembly);
        }
    });
    assembly
}

/// Starts a record of the assemblies created on this thread. Assemblies live for the whole run, so the
/// ones a task creates and then discards are found through this record and emptied (`release_discarded`).
pub fn start_tracking_created() {
    CREATED.with(|c| *c.borrow_mut() = Some(Vec::new()));
}

/// Ends the record and returns the assemblies created since `start_tracking_created`.
pub fn take_created() -> Vec<AsmP> {
    CREATED.with(|c| c.borrow_mut().take()).unwrap_or_default()
}

pub struct JunctionAssembly {
    assembly_id: i32,
    junction: JuncP,
    /// the junction group whose reads this assembly's candidate lists refer to, until they are dropped
    /// (for the release of the group's reads; not in the Java class)
    read_source_group: Option<JuncGroupP>,
    junction_index: i32,
    ref_base_position: i32,
    ref_base_cigar_elements: Vec<CigarElement>,
    indel_coords: Option<A<IndelCoords>>,
    has_line_sequence: bool,
    bases: Vec<u8>,
    base_quals: Vec<u8>,
    support: Vec<SupP>,
    candidate_support: Vec<ReadP>,
    concordant_candidates: Vec<ReadP>,
    unmapped_candidates: Vec<ReadP>,
    ref_side_soft_clips: Vec<RefSideSoftClipP>,
    repeat_info: Vec<RepeatInfo>,
    ref_bases_repeated_trimmed: Option<String>,
    ref_base_trim_length: i32,
    remote_regions: Vec<RemoteRegionP>,
    phase_group: Option<PhaseGroupP>,
    non_phased_id: i32,
    outcome: AssemblyOutcome,
    assembly_alignment_info: Option<String>,
    saga_match: Option<SagaSequenceMatch>,
    initial_read_id: Option<String>,
    merged_assemblies: i32,
    mismatch_read_count: i32,
    ext_base_build_info: Option<String>,
    stats: AssemblyStats,
}

impl JunctionAssembly {
    pub fn new(junction: JuncP, bases: Vec<u8>, base_qualities: Vec<u8>, assembly_support: &[SupP], repeat_info: Vec<RepeatInfo>) -> AsmP {
        let mut max_junction_base_qual_read: Option<ReadP> = None;
        let mut max_junction_base_qual_total = 0;
        let mut indel_coords: Option<A<IndelCoords>> = None;
        for support in assembly_support {
            let read = support.cached_read().expect("support read without cached read").clone();
            let junction_base_qual_total = read_qual_from_junction(&read, &junction);
            if junction_base_qual_total > max_junction_base_qual_total {
                max_junction_base_qual_total = junction_base_qual_total;
                max_junction_base_qual_read = Some(read.clone());
            }
            if indel_coords.is_none() {
                if let Some(ic) = read.indel_coords() {
                    if ic.length >= MIN_INDEL_LENGTH {
                        if ic.is_insert() {
                            let inserted = find_inserted_bases(&read);
                            ic.m().set_inserted_bases(inserted);
                        }
                        indel_coords = Some(ic);
                    }
                }
            }
        }
        let initial_read_id = match &max_junction_base_qual_read {
            Some(r) => r.id().to_string(),
            None => match assembly_support.first() {
                Some(s) => s.id().to_string(),
                None => "null".to_string(),
            },
        };
        let junction_index = if junction.is_forward() { 0 } else { bases.len() as i32 - 1 };
        new_tracked(JunctionAssembly {
            assembly_id: -1,
            junction,
            read_source_group: None,
            junction_index,
            ref_base_position: junction.position,
            ref_base_cigar_elements: Vec::new(),
            indel_coords,
            has_line_sequence: false,
            bases,
            base_quals: base_qualities,
            support: assembly_support.to_vec(),
            candidate_support: Vec::new(),
            concordant_candidates: Vec::new(),
            unmapped_candidates: Vec::new(),
            ref_side_soft_clips: Vec::new(),
            repeat_info,
            ref_bases_repeated_trimmed: None,
            ref_base_trim_length: 0,
            remote_regions: Vec::new(),
            phase_group: None,
            non_phased_id: -1,
            outcome: AssemblyOutcome::Unset,
            assembly_alignment_info: None,
            saga_match: None,
            initial_read_id: Some(initial_read_id),
            merged_assemblies: 0,
            mismatch_read_count: 0,
            ext_base_build_info: None,
            stats: AssemblyStats::new(),
        })
    }

    /// An assembly branched from an initial one, with its ref bases cut at a ref-side soft clip.
    pub fn new_branched(initial_assembly: AsmP, ref_side_soft_clip: RefSideSoftClipP, ref_base_length: i32, initial_support: &[SupP]) -> AsmP {
        let junction = initial_assembly.junction();
        let extension_length = initial_assembly.extension_length();
        let new_base_length = extension_length + ref_base_length;
        let initial_base_length = initial_assembly.base_length();
        let base_length_diff = initial_base_length - new_base_length;
        let mut bases = vec![0u8; new_base_length as usize];
        let mut base_quals = vec![0u8; new_base_length as usize];
        let mut assembly_index_offset = 0;
        let junction_index;
        if junction.is_forward() {
            assembly_index_offset = base_length_diff;
            junction_index = new_base_length - extension_length - 1;
        } else {
            junction_index = initial_assembly.junction_index();
        }
        let mut assembly_index = assembly_index_offset;
        for i in 0..bases.len() {
            if assembly_index < 0 {
                bases[i] = 0;
                base_quals[i] = 0;
            } else {
                if assembly_index as usize >= initial_assembly.bases().len() {
                    break;
                }
                bases[i] = initial_assembly.bases()[assembly_index as usize];
                base_quals[i] = initial_assembly.base_quals()[assembly_index as usize];
            }
            assembly_index += 1;
        }
        let mut stats = AssemblyStats::new();
        let mut support: Vec<SupP> = Vec::new();
        let mut initial_read: Option<SupP> = None;
        for s in initial_support {
            support.push(s.clone());
            stats.add_read(s, &junction, s.cached_read().map(|r| &**r));
            if Some(s.id()) == initial_assembly.initial_read_id() {
                initial_read = Some(s.clone());
            }
        }
        let initial_read_id = match initial_read {
            Some(r) => r.id().to_string(),
            None => support.first().map(|s| s.id().to_string()).unwrap_or_default(),
        };
        new_tracked(JunctionAssembly {
            assembly_id: -1,
            junction,
            read_source_group: None,
            junction_index,
            ref_base_position: ref_side_soft_clip.position,
            ref_base_cigar_elements: initial_assembly.ref_base_cigar_elements.clone(),
            indel_coords: initial_assembly.indel_coords.clone(),
            has_line_sequence: false,
            bases,
            base_quals,
            support,
            candidate_support: Vec::new(),
            concordant_candidates: Vec::new(),
            unmapped_candidates: Vec::new(),
            ref_side_soft_clips: vec![ref_side_soft_clip],
            repeat_info: Vec::new(),
            ref_bases_repeated_trimmed: Some(String::new()),
            ref_base_trim_length: 0,
            remote_regions: Vec::new(),
            phase_group: None,
            non_phased_id: 0,
            outcome: AssemblyOutcome::Unset,
            assembly_alignment_info: None,
            saga_match: None,
            initial_read_id: Some(initial_read_id),
            merged_assemblies: 0,
            mismatch_read_count: 0,
            ext_base_build_info: None,
            stats,
        })
    }

    /// the test constructor JunctionAssembly(junction, bases, quals, junctionIndex)
    pub fn new_test(junction: JuncP, bases: &[u8], quals: &[u8], junction_index: i32) -> AsmP {
        let ref_base_position = if junction.is_forward() { junction.position - junction_index } else { junction.position + (bases.len() as i32 - junction_index) - 1 };
        new_tracked(JunctionAssembly {
            assembly_id: -1,
            junction,
            read_source_group: None,
            junction_index,
            ref_base_position,
            ref_base_cigar_elements: Vec::new(),
            indel_coords: None,
            has_line_sequence: false,
            bases: bases.to_vec(),
            base_quals: quals.to_vec(),
            support: Vec::new(),
            candidate_support: Vec::new(),
            concordant_candidates: Vec::new(),
            unmapped_candidates: Vec::new(),
            ref_side_soft_clips: Vec::new(),
            repeat_info: Vec::new(),
            ref_bases_repeated_trimmed: Some(String::new()),
            ref_base_trim_length: 0,
            remote_regions: Vec::new(),
            phase_group: None,
            non_phased_id: 0,
            outcome: AssemblyOutcome::Unset,
            assembly_alignment_info: None,
            saga_match: None,
            initial_read_id: None,
            merged_assemblies: 0,
            mismatch_read_count: 0,
            ext_base_build_info: None,
            stats: AssemblyStats::new(),
        })
    }

    pub fn set_id(&mut self, id: i32) {
        self.assembly_id = id;
    }
    pub fn id(&self) -> i32 {
        self.assembly_id
    }
    #[inline]
    pub fn junction(&self) -> JuncP {
        self.junction
    }
    pub fn is_forward_junction(&self) -> bool {
        self.junction.is_forward()
    }
    pub fn is_reverse_junction(&self) -> bool {
        self.junction.is_reverse()
    }
    pub fn indel(&self) -> bool {
        self.junction.indel_based()
    }
    pub fn discordant_only(&self) -> bool {
        self.junction.discordant_only
    }
    pub fn indel_coords(&self) -> Option<&A<IndelCoords>> {
        self.indel_coords.as_ref()
    }
    pub fn set_indel_coords(&mut self, coords: Option<A<IndelCoords>>) {
        self.indel_coords = coords;
    }
    pub fn merged_assembly_count(&self) -> i32 {
        self.merged_assemblies
    }
    pub fn add_merged_assembly(&mut self) {
        self.merged_assemblies += 1;
    }
    #[inline]
    pub fn junction_index(&self) -> i32 {
        self.junction_index
    }
    pub fn set_junction_index(&mut self, index: i32) {
        self.junction_index = index;
    }
    pub fn lower_distance_from_junction(&self) -> i32 {
        self.junction_index
    }
    pub fn upper_distance_from_junction(&self) -> i32 {
        self.bases.len() as i32 - self.junction_index - 1
    }
    pub fn ref_base_length(&self) -> i32 {
        (if self.junction.is_forward() { self.lower_distance_from_junction() } else { self.upper_distance_from_junction() }) + 1
    }
    pub fn extension_length(&self) -> i32 {
        if self.junction.is_forward() { self.upper_distance_from_junction() } else { self.lower_distance_from_junction() }
    }
    pub fn ref_base_position(&self) -> i32 {
        self.ref_base_position
    }
    pub fn ref_base_cigar(&self) -> String {
        cigar_elements_to_str(&self.ref_base_cigar_elements)
    }
    #[inline]
    pub fn base_length(&self) -> i32 {
        self.bases.len() as i32
    }
    #[inline]
    /// the Java unit tests change assembly bases through the array that bases() returns
    #[cfg(test)]
    pub fn bases_mut(&mut self) -> &mut Vec<u8> {
        &mut self.bases
    }
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
    #[inline]
    pub fn base_quals(&self) -> &[u8] {
        &self.base_quals
    }
    pub fn initial_read_id(&self) -> Option<&str> {
        self.initial_read_id.as_deref()
    }
    pub fn has_line_sequence(&self) -> bool {
        self.has_line_sequence
    }
    pub fn mark_line_sequence(&mut self) {
        self.has_line_sequence = true;
    }
    pub fn unmark_line_sequence(&mut self) {
        self.has_line_sequence = false;
    }
    pub fn support(&self) -> &Vec<SupP> {
        &self.support
    }
    pub fn support_mut(&mut self) -> &mut Vec<SupP> {
        &mut self.support
    }
    pub fn support_count(&self) -> usize {
        self.support.len()
    }
    pub fn add_candidate_support(&mut self, read: ReadP) {
        self.candidate_support.push(read);
        self.stats.candidate_support_count += 1;
    }
    pub fn candidate_support(&self) -> &Vec<ReadP> {
        &self.candidate_support
    }
    pub fn candidate_support_mut(&mut self) -> &mut Vec<ReadP> {
        &mut self.candidate_support
    }
    pub fn add_concordant_candidate(&mut self, read: ReadP) {
        self.concordant_candidates.push(read);
    }
    pub fn concordant_candidates(&self) -> &Vec<ReadP> {
        &self.concordant_candidates
    }
    pub fn add_unmapped_read(&mut self, read: ReadP) {
        self.unmapped_candidates.push(read);
        self.stats.unmapped_read_count += 1;
    }
    pub fn unmapped_reads(&self) -> &Vec<ReadP> {
        &self.unmapped_candidates
    }
    pub fn unmapped_reads_mut(&mut self) -> &mut Vec<ReadP> {
        &mut self.unmapped_candidates
    }
    /// Empties an assembly that was created and then discarded, so that the reads and buffers it holds
    /// are freed (Java leaves this to the garbage collector). The assembly must not be used afterwards.
    pub fn release_discarded(&mut self) {
        self.support = Vec::new();
        self.candidate_support = Vec::new();
        self.concordant_candidates = Vec::new();
        self.unmapped_candidates = Vec::new();
        self.bases = Vec::new();
        self.base_quals = Vec::new();
        self.repeat_info = Vec::new();
        self.ref_side_soft_clips = Vec::new();
        self.remote_regions = Vec::new();
        crate::assembly::mem_stats::ASSEMBLIES_RELEASED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
    pub fn set_read_source_group(&mut self, junction_group: JuncGroupP) {
        self.read_source_group = Some(junction_group);
    }
    /// To be called once the candidate reads and the cached reads of the support have been dropped at the
    /// end of phase set building: tells the junction group, which frees its reads when all of its
    /// assemblies have reported.
    pub fn candidate_reads_dropped(&mut self) {
        if let Some(junction_group) = self.read_source_group.take() {
            junction_group.m().assembly_reads_dropped();
        }
    }
    /// Drops the reads an assembly still refers to once phase set building is over (not in Java, where
    /// assemblies that phase set building skipped keep them to the end): nothing looks at the candidate
    /// reads or the cached reads of the support after that stage, and no statistics are taken here.
    pub fn drop_reads_after_phasing(&mut self) {
        self.clear_candidate_support();
        for read in &self.support {
            read.m().clear_cached_read();
        }
        self.read_source_group = None;
    }
    pub fn clear_candidate_support(&mut self) {
        self.candidate_support = Vec::new();
        self.unmapped_candidates = Vec::new();
        self.concordant_candidates = Vec::new();
    }
    pub fn stats(&self) -> &AssemblyStats {
        &self.stats
    }
    pub fn stats_mut(&mut self) -> &mut AssemblyStats {
        &mut self.stats
    }
    pub fn mismatch_read_count(&self) -> i32 {
        self.mismatch_read_count
    }
    pub fn add_mismatch_read_count(&mut self, count: i32) {
        self.mismatch_read_count += count;
    }

    pub fn add_ref_base_support_read(&mut self, read: &ReadP, read_assembly_indices: &ReadAssemblyIndices, type_: SupportType) {
        if read_assembly_indices.is_invalid_indices() {
            return;
        }
        let mut mismatch_count = 0;
        let mut high_qual_match_count = 0;
        let mut assembly_index = read_assembly_indices.assembly_index_start;
        let read_bases = read.get_bases();
        let read_quals = read.get_base_quality();
        let mut i = read_assembly_indices.read_index_start;
        while i <= read_assembly_indices.read_index_end {
            if assembly_index < 0 {
                i += 1;
                assembly_index += 1;
                continue;
            }
            if assembly_index as usize >= self.bases.len() || i as usize >= read_bases.len() {
                break;
            }
            let base = read_bases[i as usize];
            let qual = read_quals[i as usize];
            let ai = assembly_index as usize;
            if self.bases[ai] == 0 {
                self.bases[ai] = base;
                self.base_quals[ai] = qual;
                if above_min_qual(qual) {
                    high_qual_match_count += 1;
                }
            } else if self.bases[ai] == base || below_min_qual(qual) {
                if qual > self.base_quals[ai] {
                    self.base_quals[ai] = qual;
                }
                if above_min_qual(qual) {
                    high_qual_match_count += 1;
                }
            } else if below_min_qual(self.base_quals[ai]) || qual > self.base_quals[ai] {
                self.bases[ai] = base;
                self.base_quals[ai] = qual;
            } else {
                mismatch_count += 1;
            }
            i += 1;
            assembly_index += 1;
        }
        let junction_read_start_distance = read_assembly_indices.junction_read_start_distance(self.junction_index);
        let support = SupportRead::new(read, type_, junction_read_start_distance, high_qual_match_count, mismatch_count);
        self.support.push(support);
    }

    pub fn add_junction_support(&mut self, read: &ReadP, type_: SupportType, junction_read_start_distance: i32, read_info: &ReadParseState) {
        let is_indel_crossing_junction = converted_indel_crosses_junction(&self.junction, read);
        let adjusted_type = if type_ == SupportType::Junction && is_indel_crossing_junction { SupportType::Indel } else { type_ };
        let matched_bases = read_info.matched_bases();
        let mismatch_count = read_info.mismatch_count(true);
        let support = SupportRead::new(read, adjusted_type, junction_read_start_distance, matched_bases, mismatch_count);
        read_info.read().m().set_extension_mismatches(read_info.mismatches_shared());
        support.m().set_extension_match_info(read_info);
        self.support.push(support);
    }

    pub fn set_ref_bases(&mut self, ref_base_seq_builder: &RefBaseSeqBuilder) {
        self.ref_base_position = ref_base_seq_builder.ref_base_position();
        self.ref_base_cigar_elements.extend_from_slice(ref_base_seq_builder.cigar_elements());
        let existing_bases = std::mem::take(&mut self.bases);
        let existing_quals = std::mem::take(&mut self.base_quals);
        let ref_base_extension = ref_base_seq_builder.ref_base_length() - 1;
        let new_base_length = existing_bases.len() + ref_base_extension as usize;
        let is_forward_junction = self.junction.is_forward();
        let mut ref_base_index = if is_forward_junction { 0usize } else { 1 };
        let base_offset = if is_forward_junction { ref_base_extension as usize } else { 0 };
        if is_forward_junction {
            self.junction_index += ref_base_extension;
        }
        self.bases = vec![0u8; new_base_length];
        self.base_quals = vec![0u8; new_base_length];
        for i in 0..new_base_length {
            if is_forward_junction {
                if i < base_offset {
                    self.bases[i] = ref_base_seq_builder.bases()[ref_base_index];
                    self.base_quals[i] = ref_base_seq_builder.base_qualities()[ref_base_index];
                    ref_base_index += 1;
                } else {
                    self.bases[i] = existing_bases[i - base_offset];
                    self.base_quals[i] = existing_quals[i - base_offset];
                }
            } else if i < existing_bases.len() {
                self.bases[i] = existing_bases[i];
                self.base_quals[i] = existing_quals[i];
            } else {
                self.bases[i] = ref_base_seq_builder.bases()[ref_base_index];
                self.base_quals[i] = ref_base_seq_builder.base_qualities()[ref_base_index];
                ref_base_index += 1;
            }
        }
        for i in 0..self.support.len() {
            let read = self.support[i].clone();
            let read_state = &ref_base_seq_builder.reads()[i];
            if read_state.is_valid() && !read_state.mismatched() {
                read.m().set_reference_mismatches(read_state.mismatch_count(true));
                let cached = read.cached_read().expect("support read without cached read").clone();
                self.check_add_ref_side_soft_clip(&cached);
            }
        }
    }

    pub fn trim_ref_base_position(&mut self, new_ref_base_position: i32) {
        let len = self.bases.len() as i32;
        if self.is_forward_junction() {
            let trim_length = new_ref_base_position - self.ref_base_position;
            if trim_length <= 0 || trim_length >= len - 1 {
                return;
            }
            self.bases = subset_array(&self.bases, trim_length as usize, (len - 1) as usize);
            self.base_quals = subset_array(&self.base_quals, trim_length as usize, (len - 1) as usize);
            self.junction_index -= trim_length;
        } else {
            let trim_length = self.ref_base_position - new_ref_base_position;
            if trim_length <= 0 || trim_length >= len - 1 {
                return;
            }
            self.bases = subset_array(&self.bases, 0, (len - 1 - trim_length) as usize);
            self.base_quals = subset_array(&self.base_quals, 0, (len - 1 - trim_length) as usize);
        }
        self.ref_base_position = new_ref_base_position;
    }

    pub fn extend_ref_bases(&mut self, new_ref_base_position: i32, ref_genome: Option<&dyn RefGenomeInterface>) {
        let existing_bases = std::mem::take(&mut self.bases);
        let existing_quals = std::mem::take(&mut self.base_quals);
        let ref_base_extension = (new_ref_base_position - self.ref_base_position).abs();
        let new_base_length = existing_bases.len() + ref_base_extension as usize;
        let existing_ref_base_position = self.ref_base_position;
        self.ref_base_position = new_ref_base_position;
        let is_forward_junction = self.junction.is_forward();
        if is_forward_junction {
            self.junction_index += ref_base_extension;
        }
        let base_offset = if is_forward_junction { ref_base_extension as usize } else { 0 };
        self.bases = vec![0u8; new_base_length];
        self.base_quals = vec![0u8; new_base_length];
        let ref_genome_bases: Option<Vec<u8>> = ref_genome.and_then(|g| {
            let (ref_base_start, ref_base_end) =
                if is_forward_junction { (new_ref_base_position, existing_ref_base_position - 1) } else { (existing_ref_base_position + 1, new_ref_base_position) };
            g.get_bases(self.junction.chromosome, ref_base_start, ref_base_end)
        });
        let mut ref_base_index = 0usize;
        for i in 0..new_base_length {
            let from_ref = if is_forward_junction { i < base_offset } else { i >= existing_bases.len() };
            if from_ref {
                self.bases[i] = match &ref_genome_bases {
                    Some(b) if ref_base_index < b.len() => b[ref_base_index],
                    _ => 0,
                };
                ref_base_index += 1;
                self.base_quals[i] = 0;
            } else if is_forward_junction {
                self.bases[i] = existing_bases[i - base_offset];
                self.base_quals[i] = existing_quals[i - base_offset];
            } else {
                self.bases[i] = existing_bases[i];
                self.base_quals[i] = existing_quals[i];
            }
        }
    }

    pub fn expand_extension_bases(&mut self, extension_bases: &[u8], extension_base_quals: &[u8], support_reads: &[SupP]) {
        let existing_bases = std::mem::take(&mut self.bases);
        let existing_quals = std::mem::take(&mut self.base_quals);
        let new_extension_length = extension_bases.len() as i32;
        let existing_base_length = existing_bases.len() as i32;
        let existing_extension_length = if self.junction.is_forward() { existing_base_length - self.junction_index - 1 } else { self.junction_index };
        let base_length_change = new_extension_length - existing_extension_length;
        let new_base_length = existing_base_length + base_length_change;
        if self.junction.is_reverse() {
            self.junction_index += base_length_change;
        }
        self.bases = vec![0u8; new_base_length as usize];
        self.base_quals = vec![0u8; new_base_length as usize];
        let mut new_ext_base_index = 0usize;
        for i in 0..new_base_length {
            let iu = i as usize;
            if self.junction.is_forward() {
                if i < existing_base_length {
                    self.bases[iu] = existing_bases[iu];
                    self.base_quals[iu] = existing_quals[iu];
                } else {
                    self.bases[iu] = extension_bases[new_ext_base_index];
                    self.base_quals[iu] = extension_base_quals[new_ext_base_index];
                }
                if i > self.junction_index {
                    new_ext_base_index += 1;
                }
            } else if i < base_length_change {
                self.bases[iu] = extension_bases[new_ext_base_index];
                self.base_quals[iu] = extension_base_quals[new_ext_base_index];
                new_ext_base_index += 1;
            } else {
                self.bases[iu] = existing_bases[(i - base_length_change) as usize];
                self.base_quals[iu] = existing_quals[(i - base_length_change) as usize];
            }
        }
        for support in support_reads {
            self.support.push(support.clone());
        }
    }

    pub fn remove_support_reads(&mut self, read_ids: &HashSet<String>) {
        self.support.retain(|s| !read_ids.contains(s.id()));
    }

    pub fn trim_ref_bases(&mut self, new_ref_position: i32) {
        let existing_bases = std::mem::take(&mut self.bases);
        let existing_quals = std::mem::take(&mut self.base_quals);
        let existing_ref_position = self.ref_base_position;
        let new_base_length = existing_bases.len() as i32 - (existing_ref_position - new_ref_position).abs();
        let base_length_change = existing_bases.len() as i32 - new_base_length;
        if self.junction.is_forward() {
            self.junction_index -= base_length_change;
        }
        self.ref_base_position = new_ref_position;
        let base_offset = if self.junction.is_forward() { base_length_change } else { 0 };
        self.bases = vec![0u8; new_base_length as usize];
        self.base_quals = vec![0u8; new_base_length as usize];
        for i in 0..new_base_length {
            self.bases[i as usize] = existing_bases[(i + base_offset) as usize];
            self.base_quals[i as usize] = existing_quals[(i + base_offset) as usize];
        }
        self.build_repeat_info();
    }

    pub fn has_read_support(&self, read: Option<&ReadP>) -> bool {
        match read {
            None => false,
            Some(r) => self.support.iter().any(|x| x.cached_read() == Some(r)),
        }
    }

    pub fn clear_support_cached_reads(&mut self) {
        if !clear_cached_read() {
            return;
        }
        for read in &self.support {
            if let Some(cached) = read.cached_read().cloned() {
                self.stats.add_read(read, &self.junction, Some(&cached));
                read.m().clear_cached_read();
            }
        }
    }

    pub fn repeat_info(&self) -> &Vec<RepeatInfo> {
        &self.repeat_info
    }

    pub fn build_repeat_info(&mut self) {
        self.repeat_info.clear();
        let repeats = find_repeats(&self.bases);
        if !repeats.is_empty() {
            self.repeat_info = repeats;
            self.ref_bases_repeated_trimmed = Some(build_trimmed_ref_base_sequence(self, MIN_VARIANT_LENGTH));
        } else {
            self.ref_bases_repeated_trimmed = Some(self.form_ref_base_sequence_max(MIN_VARIANT_LENGTH));
        }
        self.ref_base_trim_length = calc_trimmed_ref_base_length(self);
    }

    /// null in Java until repeat info has been built
    pub fn ref_bases_repeated_trimmed(&self) -> Option<&str> {
        self.ref_bases_repeated_trimmed.as_deref()
    }
    pub fn ref_base_trim_length(&self) -> i32 {
        self.ref_base_trim_length
    }
    pub fn ref_side_soft_clips(&self) -> &Vec<RefSideSoftClipP> {
        &self.ref_side_soft_clips
    }
    pub fn ref_side_soft_clips_mut(&mut self) -> &mut Vec<RefSideSoftClipP> {
        &mut self.ref_side_soft_clips
    }
    pub fn check_add_ref_side_soft_clip(&mut self, read: &ReadP) -> bool {
        RefSideSoftClip::check_add_ref_side_soft_clip(&mut self.ref_side_soft_clips, &self.junction, read)
    }
    pub fn remote_regions(&self) -> &Vec<RemoteRegionP> {
        &self.remote_regions
    }
    pub fn remote_regions_mut(&mut self) -> &mut Vec<RemoteRegionP> {
        &mut self.remote_regions
    }
    pub fn add_remote_regions(&mut self, regions: &[RemoteRegionP]) {
        self.remote_regions.extend_from_slice(regions);
    }
    pub fn phase_group(&self) -> Option<PhaseGroupP> {
        self.phase_group
    }
    pub fn set_phase_group(&mut self, phase_group: Option<PhaseGroupP>) {
        self.phase_group = phase_group;
    }
    /// JunctionAssembly.phaseSet()
    pub fn phase_set(this: AsmP) -> Option<PhaseSetP> {
        this.phase_group.and_then(|pg| pg.find_phase_set(this))
    }
    pub fn non_phase_set_id(&self) -> i32 {
        self.non_phased_id
    }
    pub fn set_non_phase_set_id(&mut self, id: i32) {
        self.non_phased_id = id;
    }
    pub fn outcome(&self) -> AssemblyOutcome {
        self.outcome
    }
    pub fn set_outcome(&mut self, outcome: AssemblyOutcome) {
        self.set_outcome_override(outcome, false)
    }
    pub fn set_outcome_override(&mut self, outcome: AssemblyOutcome, override_: bool) {
        if !override_ && self.outcome.ordinal() <= outcome.ordinal() {
            return;
        }
        self.outcome = outcome;
    }
    pub fn set_assembly_alignment_info(&mut self, info: Option<String>) {
        self.assembly_alignment_info = info;
    }
    pub fn assembly_alignment_info(&self) -> String {
        match &self.assembly_alignment_info {
            Some(i) => i.clone(),
            None => self.junction.coords(),
        }
    }
    pub fn set_ext_base_build_info(&mut self, info: String) {
        self.ext_base_build_info = Some(info);
    }
    pub fn ext_base_build_info(&self) -> &str {
        self.ext_base_build_info.as_deref().unwrap_or("")
    }

    pub fn form_full_sequence(&self) -> String {
        bytes_to_string(&self.bases)
    }
    pub fn form_junction_sequence(&self) -> String {
        self.form_junction_sequence_ref(0)
    }
    pub fn form_junction_sequence_ref(&self, include_ref_base_count: i32) -> String {
        let (seq_index_start, seq_index_end) = if self.junction.is_forward() {
            ((self.junction_index + 1 - include_ref_base_count).max(0), self.bases.len() as i32 - 1)
        } else {
            (0, (self.junction_index - 1 + include_ref_base_count).min(self.bases.len() as i32 - 1))
        };
        self.form_sequence(seq_index_start, seq_index_end)
    }
    pub fn form_ref_base_sequence(&self) -> String {
        self.form_ref_base_sequence_max(-1)
    }
    pub fn form_ref_base_sequence_max(&self, max_ref_base_count: i32) -> String {
        let seq_index_start;
        let seq_index_end;
        if self.junction.is_forward() {
            seq_index_end = self.junction_index;
            seq_index_start = if max_ref_base_count > 0 { seq_index_end - max_ref_base_count + 1 } else { 0 };
        } else {
            seq_index_start = self.junction_index;
            seq_index_end = if max_ref_base_count > 0 { seq_index_start + max_ref_base_count - 1 } else { self.bases.len() as i32 - 1 };
        }
        self.form_sequence(seq_index_start, seq_index_end)
    }
    pub fn form_sequence(&self, seq_index_start: i32, seq_index_end: i32) -> String {
        let from = seq_index_start.max(0);
        let to = seq_index_end.min(self.bases.len() as i32 - 1);
        if from > to {
            return String::new();
        }
        bytes_to_string(&self.bases[from as usize..=to as usize])
    }

    pub fn min_aligned_position(&self) -> i32 {
        if self.junction.is_forward() { self.ref_base_position } else { self.junction.position }
    }
    pub fn max_aligned_position(&self) -> i32 {
        if self.junction.is_reverse() { self.ref_base_position } else { self.junction.position }
    }

    pub fn is_saga_matched(&self) -> bool {
        self.saga_match.is_some()
    }
    pub fn saga_match(&self) -> Option<&SagaSequenceMatch> {
        self.saga_match.as_ref()
    }
    pub fn match_to_saga(&mut self, saga_matcher: &mut SagaSequenceMatcher) -> bool {
        let junction_offset = if self.junction.is_forward() { self.ref_base_length() } else { self.base_length() - self.ref_base_length() };
        let junction_infos = vec![SagaJunctionInfo { offset: junction_offset }];
        let lower_junction_overlap = self.has_line_sequence();
        self.saga_match = saga_matcher.match_by_sequence(&self.bases, &junction_infos, lower_junction_overlap, true);
        self.saga_match.is_some()
    }

    /// test helper addJunctionRead(read)
    pub fn add_junction_read(&mut self, read: &ReadP) {
        let junction_read_start_distance = self.junction.position - read.unclipped_start();
        let extension_bases = if self.junction.is_forward() { read.unclipped_end() - self.junction.position } else { self.junction.position - read.unclipped_start() };
        let support = SupportRead::new(read, SupportType::Junction, junction_read_start_distance, extension_bases, 0);
        support.m().set_reference_mismatches(0);
        self.support.push(support);
    }
    /// test helper addJunctionRead(read, type, junctionReadStartDistance)
    pub fn add_junction_read_typed(&mut self, read: &ReadP, type_: SupportType, junction_read_start_distance: i32) {
        let is_indel_crossing_junction = converted_indel_crosses_junction(&self.junction, read);
        let adjusted_type = if type_ == SupportType::Junction && is_indel_crossing_junction { SupportType::Indel } else { type_ };
        let support = SupportRead::new(read, adjusted_type, junction_read_start_distance, junction_read_start_distance, 0);
        self.support.push(support);
    }
}

/// new String(byte[]) for base sequences (ASCII, possibly with unset zero bytes)
pub fn bytes_to_string(bytes: &[u8]) -> String {
    bytes.iter().map(|b| *b as char).collect()
}

impl std::fmt::Display for JunctionAssembly {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "junc({}) coords(extLen={} refLen={} refBasePos={} len={} juncIndex={}) support({}) mismatches({})",
            self.junction.coords_typed(),
            self.extension_length(),
            self.ref_base_length(),
            self.ref_base_position(),
            self.base_length(),
            self.junction_index,
            self.support.len(),
            self.mismatch_read_count
        )
    }
}
