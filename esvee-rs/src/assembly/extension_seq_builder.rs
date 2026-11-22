//! ExtensionSeqBuilder: the consensus of the bases that extend past a junction, and the reads that support it.

use crate::assembly::common_utils::above_min_qual;
use crate::assembly::constants::*;
use crate::assembly::globals::assembly_build_debug;
use crate::assembly::line_utils::{find_consensus_line_extension, has_line_tail};
use crate::assembly::ptr::P;
use crate::assembly::read::read::ReadP;
use crate::assembly::read::read_utils::{get_read_index_at_junction, INVALID_INDEX};
use crate::assembly::read_parse_state::{ReadParseState, RpsP};
use crate::assembly::sequence_builder::{get_read_repeat_count, SequenceBuilder, NEXT_BASE_CHECK_COUNT};
use crate::assembly::sequence_compare::permitted_repeat_count;
use crate::assembly::sequence_diff::{SequenceDiffInfo, SequenceDiffType};
use crate::assembly::types::junction::JuncP;
use crate::assembly::types::repeat_info::{get_repeat_count, RepeatInfo};
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::SupportType;
use crate::assembly::assembly_utils::bases_match;
use crate::common::java_sort::java_sort;

/// JunctionThresholdState
#[derive(Default, Clone, Debug)]
pub struct JunctionThresholdState {
    pub extension_length_valid: bool,
    pub second_extension_length_valid: bool,
    pub min_reads_valid: bool,
    pub uses_lower_saga_limits: bool,
    pub is_line: bool,
    pub require_non_extension_support: bool,
    pub min_required_reads: i32,
    pub min_extension_length: i32,
    pub min_second_extension_length: i32,
}

impl JunctionThresholdState {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn is_valid(&self) -> bool {
        self.extension_length_valid && self.second_extension_length_valid && self.min_reads_valid
    }
}

impl std::fmt::Display for JunctionThresholdState {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "minReads({} req={}) mainLength({} req={}) secondLength({} req={}) line({}) allowJS({}) sagaLimits({})",
            self.min_reads_valid,
            self.min_required_reads,
            self.extension_length_valid,
            self.min_extension_length,
            self.second_extension_length_valid,
            self.min_second_extension_length,
            self.is_line,
            self.require_non_extension_support,
            self.uses_lower_saga_limits
        )
    }
}

/// ReadJunctionInfo
pub struct ReadJunctionInfo {
    pub read: ReadP,
    pub extension_length: i32,
    pub matches_junction: bool,
}

pub struct ExtensionSeqBuilder {
    junction: JuncP,
    reads: Vec<RpsP>,
    build_forwards: bool,
    sequence_builder: SequenceBuilder,
    has_line_sequence: bool,
    ref_repeat_count: i32,
}

pub type ExtBuilderP = P<ExtensionSeqBuilder>;

impl ExtensionSeqBuilder {
    pub fn new(junction: JuncP, reads: &[ReadP]) -> ExtensionSeqBuilder {
        let build_forwards = junction.is_forward();
        let mut rps: Vec<RpsP> = Vec::with_capacity(reads.len());
        let mut max_extension = 0;
        let mut has_line_reads = false;
        for read in reads {
            let read_junction_index = get_read_index_at_junction(read, &junction, true);
            if read_junction_index == INVALID_INDEX {
                continue;
            }
            has_line_reads |= read.has_line_tail(build_forwards);
            let extension_length = if junction.is_forward() { read.bases_length() - read_junction_index - 1 } else { read_junction_index };
            if extension_length < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH {
                continue;
            }
            max_extension = max_extension.max(extension_length);
            rps.push(ReadParseState::new(build_forwards, read, read_junction_index));
        }
        let base_length = max_extension + 1;
        let sequence_builder = SequenceBuilder::new(&rps, build_forwards, base_length);
        let mut b = ExtensionSeqBuilder { junction, reads: rps, build_forwards, sequence_builder, has_line_sequence: false, ref_repeat_count: 0 };
        b.has_line_sequence = if has_line_reads { b.check_line_sequence() } else { false };
        b.check_repeat_in_ref_bases();
        if assembly_build_debug() {
            crate::debug!("junc({}) initial extension bases sequence {}", b.junction.coords(), b.sequence_builder.base_string());
        }
        b
    }

    pub fn extension_bases(&self) -> &[u8] {
        self.sequence_builder.bases()
    }
    pub fn extension_length(&self) -> i32 {
        self.sequence_builder.bases().len() as i32 - 1
    }
    pub fn base_qualities(&self) -> &[u8] {
        self.sequence_builder.base_quals()
    }
    pub fn repeats(&self) -> &Vec<RepeatInfo> {
        self.sequence_builder.repeats()
    }
    pub fn ref_base_repeat_count(&self) -> i32 {
        self.ref_repeat_count
    }
    pub fn reads(&self) -> &Vec<RpsP> {
        &self.reads
    }

    pub fn form_assembly_support(&self) -> Vec<SupP> {
        let mut support_reads = Vec::new();
        for read in &self.reads {
            if read.mismatched() {
                continue;
            }
            if !self.sufficient_qual_matches(read) {
                continue;
            }
            let support_read = SupportRead::new(read.read(), SupportType::Junction, read.start_index(), read.matched_bases(), read.mismatch_count(true));
            read.read().m().set_extension_mismatches(read.mismatches_shared());
            support_read.m().set_extension_match_info(read);
            support_reads.push(support_read);
        }
        support_reads
    }

    pub fn sufficient_qual_matches(&self, read: &ReadParseState) -> bool {
        let overlap_bases = (read.overlap_base_count() - 1) as f64;
        if overlap_bases > 0.0 {
            let mismatch_perc = read.mismatches().len() as f64 / overlap_bases;
            if mismatch_perc > ASSEMBLY_MAX_EXTENSION_READ_LOW_QUAL_MISMATCH_PERC {
                return false;
            }
        }
        let ref_base_repeat_buffer = if self.ref_base_repeat_count() > 0 { 2 } else { 0 };
        let high_qual_extension_matches = read.high_qual_matches() - 1;
        high_qual_extension_matches >= ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH + ref_base_repeat_buffer
    }

    pub fn mismatch_reads(&self) -> Vec<ReadP> {
        self.reads.iter().filter(|x| x.mismatched()).map(|x| x.read().clone()).collect()
    }
    pub fn mismatches(&self) -> i32 {
        self.reads.iter().filter(|x| x.mismatched()).count() as i32
    }

    pub fn check_add_junction_read(&mut self, read: &ReadP) -> Option<RpsP> {
        let read_junction_index = get_read_index_at_junction(read, &self.junction, true);
        if read_junction_index == INVALID_INDEX {
            return None;
        }
        let read_parse_state = ReadParseState::new(self.build_forwards, read, read_junction_index);
        let ext_len = self.sequence_builder.bases().len() as i32;
        let mut extension_index: i32 = if self.build_forwards { 0 } else { self.sequence_builder.base_quals().len() as i32 - 1 };
        while extension_index >= 0 && extension_index < ext_len {
            let consensus_base = self.sequence_builder.bases()[extension_index as usize];
            let base = read_parse_state.current_base();
            let qual = read_parse_state.current_qual();
            let hq = above_min_qual(qual);
            let mut ext_base_move = 1;
            if base == consensus_base {
                read_parse_state.m().add_base_match(hq);
                read_parse_state.m().move_next();
            } else {
                ext_base_move = self.assess_read_mismatch(&read_parse_state, extension_index);
                if read_parse_state.mismatched() {
                    return Some(read_parse_state);
                }
            }
            if read_parse_state.exhausted() {
                break;
            }
            extension_index += if self.build_forwards { ext_base_move } else { -ext_base_move };
        }
        if self.sufficient_qual_matches(&read_parse_state) {
            let mut replaced = false;
            for i in 0..self.reads.len() {
                if self.reads[i].read() == read {
                    self.reads[i] = read_parse_state.clone();
                    replaced = true;
                    break;
                }
            }
            if !replaced {
                self.reads.push(read_parse_state.clone());
            }
        }
        Some(read_parse_state)
    }

    fn assess_read_mismatch(&self, read: &RpsP, extension_index: i32) -> i32 {
        let mut seq_diff_info = SequenceDiffInfo::unset();
        let read_next_bases = read.get_next_bases(NEXT_BASE_CHECK_COUNT);
        let mut consensus_next_bases = vec![0u8; NEXT_BASE_CHECK_COUNT as usize];
        let mut consensus_current_and_next = vec![0u8; NEXT_BASE_CHECK_COUNT as usize];
        self.sequence_builder.populate_bases(&mut consensus_current_and_next, extension_index);
        self.sequence_builder.populate_bases(&mut consensus_next_bases, extension_index + if self.build_forwards { 1 } else { -1 });

        let mut consensus_repeat: Option<RepeatInfo> = None;
        let mut previous_repeat_length = 0;
        if read_next_bases.is_none() || bases_match(read_next_bases.as_ref().unwrap(), &consensus_next_bases) {
            seq_diff_info = SequenceDiffInfo::from_snv(read, extension_index);
        } else {
            let prior_ext_index = extension_index + if self.build_forwards { -1 } else { 1 };
            consensus_repeat = self.sequence_builder.find_repeat(prior_ext_index).cloned();
            if let Some(cr) = &consensus_repeat {
                let prev_ext_base_length = (cr.index - extension_index).abs();
                let repeat_base_length = cr.repeat_length();
                let prev_repeat_count = prev_ext_base_length / repeat_base_length;
                let repeat_ext_remainder = prev_ext_base_length % repeat_base_length;
                let mut read_repeat_count = 0;
                if prev_repeat_count > 0 && repeat_ext_remainder == 0 {
                    read_repeat_count = get_read_repeat_count(read, &cr.bases, prev_repeat_count, self.build_forwards, false);
                }
                let permitted_repeat_diff = permitted_repeat_count(cr.count);
                let repeat_count_diff = (read_repeat_count - cr.count).abs();
                if read_repeat_count > 0 && repeat_count_diff <= permitted_repeat_diff {
                    previous_repeat_length = prev_ext_base_length;
                    let read_repeat_start = SequenceDiffInfo::repeat_index_static(read.read_index(), read_repeat_count, cr, previous_repeat_length, self.build_forwards, true);
                    let read_repeat_end = SequenceDiffInfo::repeat_index_static(read.read_index(), read_repeat_count, cr, previous_repeat_length, self.build_forwards, false);
                    let read_repeat_index_begin = if self.build_forwards { read_repeat_start } else { read_repeat_end };
                    let base_qual_type = read.range_min_qual_type(read_repeat_start, read_repeat_end);
                    seq_diff_info = SequenceDiffInfo::new_full(
                        read.read_index(),
                        extension_index,
                        Some(cr.bases.clone()),
                        SequenceDiffType::Repeat,
                        base_qual_type,
                        read_repeat_count,
                        read_repeat_index_begin,
                        0,
                    );
                }
            }
            if seq_diff_info.type_ == SequenceDiffType::Unset {
                seq_diff_info = self.sequence_builder.assess_novel_indel_differences(read, seq_diff_info, read_next_bases.as_deref(), &consensus_next_bases, &consensus_current_and_next);
            }
        }
        if seq_diff_info.is_unset() {
            seq_diff_info = SequenceDiffInfo::from_snv(read, extension_index);
        }
        self.sequence_builder.apply_read_mismatches(read, &mut seq_diff_info, consensus_repeat.as_ref(), previous_repeat_length);
        let mut ext_base_move = 1;
        if seq_diff_info.type_ == SequenceDiffType::Repeat {
            let cr = consensus_repeat.as_ref().unwrap();
            if cr.count > seq_diff_info.repeat_count {
                ext_base_move = (cr.count - seq_diff_info.repeat_count) * cr.repeat_length();
            } else if seq_diff_info.repeat_count > cr.count {
                ext_base_move = 0;
            }
        }
        ext_base_move
    }

    pub fn check_validity(&self, junc_threshold_state: &mut JunctionThresholdState) {
        let mut max_valid_extension_length = 0;
        let mut valid_extension_read_count = 0;
        let mut has_extension_length_read = false;
        for read in &self.reads {
            if read.mismatched() {
                continue;
            }
            let r = read.read();
            let sc_length = if self.junction.discordant_only {
                if self.junction.is_forward() { r.unclipped_end() - self.junction.position } else { self.junction.position - r.unclipped_start() }
            } else if self.junction.is_forward() {
                r.right_clip_length()
            } else {
                r.left_clip_length()
            };
            if sc_length >= junc_threshold_state.min_second_extension_length {
                has_extension_length_read |= sc_length >= junc_threshold_state.min_extension_length;
                valid_extension_read_count += 1;
            } else if r.indel_coords().is_some_and(|c| c.length >= MIN_INDEL_LENGTH) {
                valid_extension_read_count += 1;
                has_extension_length_read = true;
            } else {
                continue;
            }
            max_valid_extension_length = read.overlap_base_count().max(max_valid_extension_length);
        }
        junc_threshold_state.extension_length_valid = has_extension_length_read && max_valid_extension_length > 0;
        junc_threshold_state.second_extension_length_valid |= valid_extension_read_count >= 2;
        junc_threshold_state.min_reads_valid |= valid_extension_read_count >= junc_threshold_state.min_required_reads;
    }

    fn check_repeat_in_ref_bases(&mut self) {
        let Some(first_repeat) = self.sequence_builder.repeats().first().cloned() else { return };
        let repeat_junction_distance = if self.build_forwards { first_repeat.index } else { self.sequence_builder.base_length() - first_repeat.index - 1 };
        if repeat_junction_distance > first_repeat.repeat_length() {
            return;
        }
        let Some(sample_read) = self.reads.iter().find(|x| !x.mismatched()) else { return };
        let junction_index = sample_read.start_index();
        let repeat_index_start = junction_index + if self.build_forwards { repeat_junction_distance - 1 } else { -repeat_junction_distance + 1 };
        let previous_repeat_count = get_repeat_count(sample_read.read().get_bases(), &first_repeat.bases, repeat_index_start, !self.build_forwards);
        if previous_repeat_count > 0 {
            self.ref_repeat_count = previous_repeat_count;
        }
    }

    pub fn has_line_sequence(&self) -> bool {
        self.has_line_sequence
    }

    fn check_line_sequence(&self) -> bool {
        let valid_read_count = self.sequence_builder.reads().iter().filter(|x| !x.mismatched()).count();
        let valid_line_reads: Vec<ReadP> =
            self.sequence_builder.reads().iter().filter(|x| !x.mismatched() && x.read().has_line_tail(self.build_forwards)).map(|x| x.read().clone()).collect();
        if (valid_line_reads.len() as f64) < valid_read_count as f64 * 0.5 {
            return false;
        }
        let line_extension_length = find_consensus_line_extension(&valid_line_reads, &self.junction);
        if line_extension_length >= LINE_POLY_AT_REQ {
            return true;
        }
        let first_ext_index = if self.build_forwards { 1 } else { self.sequence_builder.base_length() - 2 };
        let line_base = if self.junction.is_forward() { LINE_BASE_A } else { LINE_BASE_T };
        has_line_tail(self.sequence_builder.bases(), first_ext_index, !self.build_forwards, line_base, None)
    }

    pub fn build_information(&self) -> String {
        let mut sj: Vec<String> = Vec::new();
        let mut exact_match = 0;
        let mut mismatched_reads = 0;
        let mut low_qual_mismatches = 0;
        let mut snvs = 0;
        let mut indels = 0;
        let mut homopolymers = 0;
        let mut other_repeats = 0;
        for read in &self.reads {
            if !self.sufficient_qual_matches(read) {
                continue;
            }
            if read.mismatch_count(false) == 0 {
                exact_match += 1;
            } else if read.mismatched() {
                mismatched_reads += 1;
            } else {
                for mismatch in read.mismatches() {
                    if mismatch.mismatch_penalty > 0.0 {
                        match mismatch.type_ {
                            SequenceDiffType::Base => snvs += 1,
                            SequenceDiffType::Repeat => {
                                if mismatch.bases.as_ref().is_some_and(|b| b.len() == 1) {
                                    homopolymers += 1;
                                } else {
                                    other_repeats += 1;
                                }
                            }
                            SequenceDiffType::Insert | SequenceDiffType::Delete => indels += 1,
                            _ => {}
                        }
                    } else {
                        low_qual_mismatches += 1;
                    }
                }
            }
        }
        sj.push(format!("RC={}", self.reads.len()));
        sj.push(format!("EM={exact_match}"));
        if low_qual_mismatches > 0 {
            sj.push(format!("LQ={low_qual_mismatches}"));
        }
        if mismatched_reads > 0 {
            sj.push(format!("MM={mismatched_reads}"));
        }
        if snvs > 0 {
            sj.push(format!("SNV={snvs}"));
        }
        if homopolymers > 0 {
            sj.push(format!("HP={homopolymers}"));
        }
        if other_repeats > 0 {
            sj.push(format!("OR={other_repeats}"));
        }
        if indels > 0 {
            sj.push(format!("ID={indels}"));
        }
        let mut repeats: Vec<RepeatInfo> = self.sequence_builder.repeats().clone();
        if repeats.len() > 5 {
            java_sort(&mut repeats, |a, b| (-a.count).cmp(&(-b.count)));
        }
        for repeat in repeats.iter().take(5) {
            sj.push(format!("REP:{}x{}", repeat.count, repeat.bases));
        }
        sj.join(";")
    }

    pub fn junction_sequence(&self) -> String {
        self.sequence_builder.base_string()
    }
}

impl std::fmt::Display for ExtensionSeqBuilder {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "junc({}) reads({}) baseLength({})", self.junction.coords_typed(), self.reads.len(), self.sequence_builder.bases().len())
    }
}
