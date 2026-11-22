//! SequenceBuilder: builds a consensus sequence base by base from a set of reads, classifying each read's
//! differences (SNV, one or two base indel, repeat expansion or contraction) and dropping reads whose
//! mismatch penalty grows too large.

use std::cmp::Ordering;

use crate::assembly::assembly_utils::{self, bases_match, bases_match_range, DNA_BASE_BYTES, DNA_BASE_COUNT, NO_BASE};
use crate::assembly::common_utils::{above_min_qual, below_min_qual, is_high_base_qual};
use crate::assembly::constants::*;
use crate::assembly::globals::{read_mismatch_long_repeat_count, read_mismatch_max_repeat_diff, read_mismatch_medium_repeat_count};
use crate::assembly::orientation::Orientation;
use crate::assembly::read_parse_state::{ReadParseState, RpsP};
use crate::assembly::sequence_diff::{BaseQualType, SequenceDiffInfo, SequenceDiffType};
use crate::assembly::types::repeat_info::{find_repeats, get_repeat_count, get_repeat_count_ext, RepeatInfo};
use crate::common::java_sort::java_sort;

pub const NEXT_BASE_CHECK_COUNT: i32 = 3;
const MIN_PREVIOUS_REPEAT_COUNT: i32 = READ_MISMATCH_LOW_REPEAT_COUNT - 1;

pub struct SequenceBuilder {
    build_forwards: bool,
    reads: Vec<RpsP>,
    mismatch_penalty_max_read_length: i32,
    disable_mismatch_penalty: bool,
    current_index: i32,
    bases: Vec<u8>,
    base_quals: Vec<u8>,
    repeats: Vec<RepeatInfo>,
}

fn read_comparator(first: &RpsP, second: &RpsP) -> Ordering {
    let first_overlap = first.overlap_base_count();
    let second_overlap = second.overlap_base_count();
    if first_overlap != second_overlap {
        return if first_overlap > second_overlap { Ordering::Less } else { Ordering::Greater };
    }
    first.read().id().cmp(second.read().id())
}

#[inline]
fn max_qual(a: u8, b: u8) -> u8 {
    a.max(b)
}

impl SequenceBuilder {
    pub fn new(reads: &[RpsP], build_forwards: bool, initial_base_length: i32) -> SequenceBuilder {
        Self::new_ext(reads, build_forwards, initial_base_length, false)
    }

    pub fn new_ext(reads: &[RpsP], build_forwards: bool, initial_base_length: i32, disable_mismatch_penalty: bool) -> SequenceBuilder {
        let mut sorted: Vec<RpsP> = reads.to_vec();
        java_sort(&mut sorted, read_comparator);
        let mismatch_penalty_max_read_length = if sorted.len() >= 2 { sorted[1].overlap_base_count() } else { -1 };
        let mut sb = SequenceBuilder {
            build_forwards,
            reads: sorted,
            mismatch_penalty_max_read_length,
            disable_mismatch_penalty,
            current_index: -1,
            bases: vec![0u8; initial_base_length as usize],
            base_quals: vec![0u8; initial_base_length as usize],
            repeats: Vec::new(),
        };
        sb.build_sequence();
        sb.trim_final_sequence();
        sb.find_non_mismatch_repeats();
        sb
    }

    /// a builder without reads or bases, for tests of the protected methods
    pub fn new_for_test(build_forwards: bool, bases: Vec<u8>, current_index: i32) -> SequenceBuilder {
        let n = bases.len();
        SequenceBuilder { build_forwards, reads: Vec::new(), mismatch_penalty_max_read_length: -1, disable_mismatch_penalty: false, current_index, bases, base_quals: vec![0u8; n], repeats: Vec::new() }
    }

    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
    pub fn base_length(&self) -> i32 {
        self.bases.len() as i32
    }
    pub fn base_string(&self) -> String {
        self.bases.iter().map(|b| *b as char).collect()
    }
    pub fn base_quals(&self) -> &[u8] {
        &self.base_quals
    }
    pub fn reads(&self) -> &Vec<RpsP> {
        &self.reads
    }
    pub fn repeats(&self) -> &Vec<RepeatInfo> {
        &self.repeats
    }

    pub fn find_repeat(&self, base_index: i32) -> Option<&RepeatInfo> {
        self.repeats.iter().find(|repeat| {
            if self.build_forwards {
                base_index >= repeat.index && base_index <= repeat.last_index()
            } else {
                base_index >= repeat.last_index() && base_index <= repeat.index
            }
        })
    }

    fn build_sequence(&mut self) {
        let len = self.bases.len() as i32;
        self.current_index = if self.build_forwards { 0 } else { len - 1 };
        let mut active_reads: Vec<RpsP> = self.reads.clone();
        let mut hp_count = 0;
        let mut hp_base = NO_BASE;

        while self.current_index >= 0 && self.current_index < len {
            let mut consensus_base: u8 = 0;
            let mut consensus_max_qual: u8 = 0;
            let mut consensus_high_qual_count = 0;
            let mut consensus_med_qual_count = 0;
            let mut has_arrays = false;
            let mut max_quals = [0u8; DNA_BASE_COUNT];
            let mut high_quals = [0i32; DNA_BASE_COUNT];
            let mut medium_quals = [0i32; DNA_BASE_COUNT];

            for read in active_reads.iter() {
                let base = read.current_base();
                let qual = read.current_qual();
                let is_high_qual = is_high_base_qual(qual);
                let mut has_mismatch = has_arrays;
                if !has_arrays {
                    if consensus_base == 0 || (base != consensus_base && below_min_qual(consensus_max_qual) && above_min_qual(qual)) {
                        consensus_base = base;
                        consensus_max_qual = qual;
                        if is_high_qual {
                            consensus_high_qual_count += 1;
                        } else {
                            consensus_med_qual_count += 1;
                        }
                    } else if base == consensus_base {
                        consensus_max_qual = max_qual(qual, consensus_max_qual);
                        if is_high_qual {
                            consensus_high_qual_count += 1;
                        } else {
                            consensus_med_qual_count += 1;
                        }
                    } else {
                        has_mismatch = true;
                        has_arrays = true;
                        let base_index = assembly_utils::base_index(consensus_base);
                        max_quals[base_index] = consensus_max_qual;
                        high_quals[base_index] = consensus_high_qual_count;
                        medium_quals[base_index] = consensus_med_qual_count;
                    }
                }
                if has_mismatch {
                    let base_index = assembly_utils::base_index(base);
                    max_quals[base_index] = max_qual(max_quals[base_index], qual);
                    if is_high_qual {
                        high_quals[base_index] += 1;
                    } else {
                        medium_quals[base_index] += 1;
                    }
                }
            }

            let mut has_consensus_repeat = false;
            if has_arrays {
                let mut max_high_qual = 0;
                let mut max_medium_qual = 0;
                let mut max_base_index = 0usize;
                let mut tie_break = false;
                for b in 0..DNA_BASE_COUNT {
                    if high_quals[b] > max_high_qual {
                        max_high_qual = high_quals[b];
                        max_medium_qual = medium_quals[b];
                        max_base_index = b;
                        tie_break = false;
                    } else if high_quals[b] == max_high_qual && medium_quals[b] > max_medium_qual {
                        max_high_qual = high_quals[b];
                        max_medium_qual = medium_quals[b];
                        max_base_index = b;
                        tie_break = false;
                    } else if high_quals[b] == max_high_qual && medium_quals[b] == max_medium_qual {
                        tie_break = true;
                    }
                }
                if tie_break {
                    let min_mismatches = active_reads.iter().map(|x| x.mismatch_penalty()).fold(f64::INFINITY, f64::min);
                    if let Some(read) = active_reads.iter().find(|x| x.mismatch_penalty() == min_mismatches) {
                        max_base_index = assembly_utils::base_index(read.current_base());
                    }
                }
                consensus_base = DNA_BASE_BYTES[max_base_index];
                consensus_max_qual = max_quals[max_base_index];
                let mut seq_diff_infos = self.assess_differences(consensus_base, &active_reads, hp_base, hp_count);
                has_consensus_repeat = self.apply_read_mismatches_and_consensus(&mut active_reads, consensus_base, consensus_max_qual, &mut seq_diff_infos);
            } else {
                self.bases[self.current_index as usize] = consensus_base;
                self.base_quals[self.current_index as usize] = consensus_max_qual;
                self.current_index += if self.build_forwards { 1 } else { -1 };
                for read in active_reads.iter() {
                    if read.current_base() == consensus_base {
                        let hq = above_min_qual(read.current_qual());
                        read.m().add_base_match(hq);
                    } else {
                        let info = SequenceDiffInfo::from_snv(read, self.current_index);
                        read.m().add_mismatch_info(info);
                    }
                    read.m().move_next();
                }
            }

            if has_consensus_repeat {
                hp_base = NO_BASE;
                hp_count = 0;
            } else if hp_base == consensus_base {
                hp_count += 1;
            } else {
                hp_base = consensus_base;
                hp_count = 1;
            }

            active_reads.retain(|read| !(read.exhausted() || read.mismatched()));
            if active_reads.is_empty() {
                break;
            }
        }
    }

    fn assess_differences(&mut self, consensus_base: u8, active_reads: &[RpsP], hp_base: u8, hp_count: i32) -> Vec<SequenceDiffInfo> {
        let mut seq_diff_infos: Vec<SequenceDiffInfo> = vec![SequenceDiffInfo::unset(); active_reads.len()];
        let read_count = active_reads.len();
        let mut all_read_next_bases: Vec<Option<Vec<u8>>> = Vec::with_capacity(read_count);
        let mut consensus_read_next_bases: Vec<Option<Vec<u8>>> = Vec::with_capacity(read_count);
        let match_info = SequenceDiffInfo::from_match(-1, self.current_index);

        for (read_index, read) in active_reads.iter().enumerate() {
            let read_next_bases = read.get_next_bases(NEXT_BASE_CHECK_COUNT);
            if read.current_base() == consensus_base {
                consensus_read_next_bases.push(read_next_bases.clone());
                seq_diff_infos[read_index] = match_info.clone();
            } else if read_next_bases.is_none() {
                seq_diff_infos[read_index] = SequenceDiffInfo::from_snv(read, self.current_index);
            }
            all_read_next_bases.push(read_next_bases);
        }

        let consensus_next_bases = find_consensus_next_bases(&consensus_read_next_bases);
        if let Some(cnb) = &consensus_next_bases {
            let mut post_consensus_mismatch_count = 0;
            for (read_index, read) in active_reads.iter().enumerate() {
                if seq_diff_infos[read_index].type_ == SequenceDiffType::Match {
                    continue;
                }
                let Some(read_next_bases) = &all_read_next_bases[read_index] else { continue };
                if bases_match(read_next_bases, cnb) {
                    seq_diff_infos[read_index] = SequenceDiffInfo::from_snv(read, self.current_index);
                } else {
                    post_consensus_mismatch_count += 1;
                }
            }
            if post_consensus_mismatch_count == 0 {
                return seq_diff_infos;
            }
        }

        if hp_count >= MIN_PREVIOUS_REPEAT_COUNT {
            self.assess_repeat_differences(active_reads, &mut seq_diff_infos, hp_base, hp_count);
            if !seq_diff_infos.iter().any(|x| x.is_unset()) {
                return seq_diff_infos;
            }
        }

        if let Some(cnb) = &consensus_next_bases {
            self.assess_novel_indel_differences_all(consensus_base, cnb, active_reads, &all_read_next_bases, &mut seq_diff_infos);
            if !seq_diff_infos.iter().any(|x| x.is_unset()) {
                return seq_diff_infos;
            }
        }

        let remaining_unset = seq_diff_infos.iter().filter(|x| x.is_unset()).count();
        if remaining_unset as f64 > TOTAL_READ_MISMATCH_RATE * read_count as f64 {
            self.assess_repeat_differences(active_reads, &mut seq_diff_infos, NO_BASE, 0);
        }
        for (read_index, read) in active_reads.iter().enumerate() {
            if seq_diff_infos[read_index].is_unset() {
                seq_diff_infos[read_index] = SequenceDiffInfo::from_snv(read, self.current_index);
            }
        }
        seq_diff_infos
    }

    fn assess_novel_indel_differences_all(
        &self, consensus_base: u8, consensus_next_bases: &[u8], active_reads: &[RpsP], all_read_next_bases: &[Option<Vec<u8>>], seq_diff_infos: &mut [SequenceDiffInfo],
    ) {
        let consensus_current_and_next = self.current_and_next_base(consensus_base, consensus_next_bases);
        for read_index in 0..active_reads.len() {
            let t = seq_diff_infos[read_index].type_;
            if t == SequenceDiffType::Match || t == SequenceDiffType::Base {
                continue;
            }
            let info = self.assess_novel_indel_differences(
                &active_reads[read_index],
                seq_diff_infos[read_index].clone(),
                all_read_next_bases[read_index].as_deref(),
                consensus_next_bases,
                &consensus_current_and_next,
            );
            seq_diff_infos[read_index] = info;
        }
    }

    pub fn assess_novel_indel_differences(
        &self, read: &ReadParseState, seq_diff_info: SequenceDiffInfo, read_next_bases: Option<&[u8]>, consensus_next_bases: &[u8], consensus_current_and_next: &[u8],
    ) -> SequenceDiffInfo {
        if seq_diff_info.type_ == SequenceDiffType::Match || seq_diff_info.type_ == SequenceDiffType::Base {
            return seq_diff_info;
        }
        let Some(read_next_bases) = read_next_bases else { return seq_diff_info };
        if bases_match(consensus_current_and_next, read_next_bases) {
            let qual_type = read.qual_type(read.read_index());
            let base = (read.current_base() as char).to_string();
            return SequenceDiffInfo::from_indel(read.read_index(), self.current_index, &base, SequenceDiffType::Insert, qual_type, 1);
        }
        let n = read_next_bases.len() as i32;
        let read_next_index_start = if self.build_forwards { 1 } else { 0 };
        let read_next_index_end = if self.build_forwards { n - 1 } else { n - 2 };
        let consensus_next_index_start = if self.build_forwards { 0 } else { 1 };
        let consensus_next_index_end = if self.build_forwards { n - 2 } else { n - 1 };
        if bases_match_range(consensus_current_and_next, consensus_next_index_start, consensus_next_index_end, read_next_bases, read_next_index_start, read_next_index_end) {
            let qual_type = read.range_min_qual_type(read.read_index(), read.read_index() + 1);
            let bases = format!("{}{}", read.current_base() as char, read_next_bases[0] as char);
            return SequenceDiffInfo::from_indel(read.read_index(), self.current_index, &bases, SequenceDiffType::Insert, qual_type, 2);
        }
        let read_current_and_next = self.current_and_next_base(read.current_base(), read_next_bases);
        if bases_match(consensus_next_bases, &read_current_and_next) {
            let previous_index = read.read_index() + if self.build_forwards { -1 } else { 1 };
            let qual_type = read.straddling_qual_type(read.read_index(), previous_index);
            return SequenceDiffInfo::from_indel(read.read_index(), self.current_index, "", SequenceDiffType::Delete, qual_type, -1);
        }
        seq_diff_info
    }

    fn assess_repeat_differences(&mut self, active_reads: &[RpsP], seq_diff_infos: &mut [SequenceDiffInfo], hp_base: u8, hp_count: i32) {
        let mut selected_repeat_bases: Option<String> = None;
        let mut count_of_reads_with_selected_repeat = 0;
        let mut selected_repeat_read_repeat_counts: Vec<i32> = Vec::new();
        let mut selected_repeat_previous_repeat_count = 0;
        let len = self.bases.len() as i32;

        if hp_count > 0 {
            let bases = (hp_base as char).to_string();
            selected_repeat_previous_repeat_count = hp_count;
            selected_repeat_read_repeat_counts = vec![0; active_reads.len()];
            for (read_index, read) in active_reads.iter().enumerate() {
                let repeat_count = self.get_read_repeat_count(read, &bases, hp_count);
                selected_repeat_read_repeat_counts[read_index] = repeat_count;
                if repeat_count > 0 {
                    count_of_reads_with_selected_repeat += 1;
                }
            }
            selected_repeat_bases = Some(bases);
        } else {
            let last_index = if self.build_forwards { self.current_index - 1 } else { self.current_index + 1 };
            if let Some(last_repeat) = self.repeats.last() {
                if (self.build_forwards && last_repeat.index_end() >= last_index) || (!self.build_forwards && last_repeat.index_start() <= last_index) {
                    return;
                }
            }
            for repeat_length in 2..=REPEAT_MAX_BASE_LENGTH {
                let repeat_bases: String;
                if self.build_forwards {
                    if self.current_index - repeat_length < 0 {
                        break;
                    }
                    repeat_bases = self.bases[(self.current_index - repeat_length) as usize..self.current_index as usize].iter().map(|b| *b as char).collect();
                } else {
                    if self.current_index + repeat_length >= len {
                        break;
                    }
                    repeat_bases = self.bases[(self.current_index + 1) as usize..(self.current_index + repeat_length + 1) as usize].iter().map(|b| *b as char).collect();
                }
                let repeat_index_start = last_index;
                let previous_repeat_count = get_repeat_count(&self.bases, &repeat_bases, repeat_index_start, !self.build_forwards);
                if previous_repeat_count < MIN_PREVIOUS_REPEAT_COUNT {
                    continue;
                }
                let mut repeat_read_repeat_counts = vec![0; active_reads.len()];
                let mut non_zero_repeat_read_count = 0;
                for (read_index, read) in active_reads.iter().enumerate() {
                    let repeat_count = self.get_read_repeat_count(read, &repeat_bases, previous_repeat_count);
                    repeat_read_repeat_counts[read_index] = repeat_count;
                    if repeat_count > 0 {
                        non_zero_repeat_read_count += 1;
                    }
                }
                if non_zero_repeat_read_count > count_of_reads_with_selected_repeat {
                    count_of_reads_with_selected_repeat = non_zero_repeat_read_count;
                    selected_repeat_read_repeat_counts = repeat_read_repeat_counts;
                    selected_repeat_bases = Some(repeat_bases);
                    selected_repeat_previous_repeat_count = previous_repeat_count;
                    if non_zero_repeat_read_count == active_reads.len() {
                        break;
                    }
                }
            }
            if selected_repeat_bases.is_none() {
                return;
            }
        }
        let selected_repeat_bases = selected_repeat_bases.unwrap();

        let max_observed_repeat_count = selected_repeat_read_repeat_counts.iter().copied().max().unwrap_or(0);
        let reads_with_min_repeat_count = selected_repeat_read_repeat_counts.iter().filter(|x| **x >= READ_MISMATCH_LOW_REPEAT_COUNT).count();
        if max_observed_repeat_count < READ_MISMATCH_LOW_REPEAT_COUNT {
            return;
        }
        let min_required_reads_with_repeat = active_reads.len() as f64 * 0.5;
        if (count_of_reads_with_selected_repeat as f64) < min_required_reads_with_repeat || (reads_with_min_repeat_count as f64) < min_required_reads_with_repeat {
            return;
        }
        let consensus_repeat_count = find_consensus_repeat_count(&selected_repeat_read_repeat_counts);
        let repeat_base_length = selected_repeat_bases.len() as i32;
        let consensus_repeat_offset = selected_repeat_previous_repeat_count * repeat_base_length;
        let consensus_repeat_index_start = if self.build_forwards { self.current_index - consensus_repeat_offset } else { self.current_index + consensus_repeat_offset };
        let consensus_repeat = RepeatInfo::with_direction(
            consensus_repeat_index_start,
            &selected_repeat_bases,
            consensus_repeat_count,
            Some(if self.build_forwards { Orientation::Forward } else { Orientation::Reverse }),
        );
        self.repeats.push(consensus_repeat.clone());

        for (read_index, read) in active_reads.iter().enumerate() {
            let read_repeat_count = selected_repeat_read_repeat_counts[read_index];
            let diff_type = if read_repeat_count == consensus_repeat_count { SequenceDiffType::Match } else { SequenceDiffType::Repeat };
            let read_repeat_start = SequenceDiffInfo::repeat_index_static(read.read_index(), read_repeat_count, &consensus_repeat, consensus_repeat_offset, self.build_forwards, true);
            let read_repeat_end = SequenceDiffInfo::repeat_index_static(read.read_index(), read_repeat_count, &consensus_repeat, consensus_repeat_offset, self.build_forwards, false);
            let read_repeat_index_begin = if self.build_forwards { read_repeat_start } else { read_repeat_end };
            let base_qual_type = read.range_min_qual_type(read_repeat_start, read_repeat_end);
            let read_repeat_base_count_diff = (read_repeat_count - consensus_repeat_count) * repeat_base_length;
            seq_diff_infos[read_index] = SequenceDiffInfo::new_full(
                read.read_index(),
                self.current_index,
                Some(selected_repeat_bases.clone()),
                diff_type,
                base_qual_type,
                read_repeat_count,
                read_repeat_index_begin,
                read_repeat_base_count_diff,
            );
        }
    }

    fn get_read_repeat_count(&self, read: &ReadParseState, repeat: &str, previous_repeat_count: i32) -> i32 {
        get_read_repeat_count(read, repeat, previous_repeat_count, self.build_forwards, true)
    }

    fn current_and_next_base(&self, current_base: u8, next_bases: &[u8]) -> Vec<u8> {
        let n = next_bases.len();
        let mut current_and_next = vec![0u8; n];
        if self.build_forwards {
            current_and_next[0] = current_base;
            current_and_next[1..n].copy_from_slice(&next_bases[..n - 1]);
        } else {
            current_and_next[n - 1] = current_base;
            // as in Java: the remaining positions take next bases n-1 down to 1
            let mut index = n as isize - 2;
            for i in (1..n).rev() {
                current_and_next[index as usize] = next_bases[i];
                index -= 1;
            }
        }
        current_and_next
    }

    /// Sets the consensus, moves each read on and registers its mismatch; returns true if a repeat was applied.
    fn apply_read_mismatches_and_consensus(&mut self, active_reads: &mut Vec<RpsP>, consensus_base: u8, consensus_qual: u8, seq_diff_infos: &mut [SequenceDiffInfo]) -> bool {
        let consensus_repeat: Option<RepeatInfo> = if seq_diff_infos.iter().any(|x| x.repeat_count > 0) { Some(self.repeats[self.repeats.len() - 1].clone()) } else { None };
        let previous_read_bases = consensus_repeat.as_ref().map(|r| (self.current_index - r.index).abs()).unwrap_or(0);
        self.set_consensus_from_mismatches(consensus_base, consensus_qual, consensus_repeat.as_ref());
        let mut has_mismatched = false;
        for (read_index, read) in active_reads.iter().enumerate() {
            self.apply_read_mismatches(read, &mut seq_diff_infos[read_index], consensus_repeat.as_ref(), previous_read_bases);
            has_mismatched |= read.mismatched();
        }
        if has_mismatched {
            active_reads.retain(|r| !r.mismatched());
        }
        consensus_repeat.is_some()
    }

    fn set_consensus_from_mismatches(&mut self, consensus_base: u8, consensus_qual: u8, consensus_repeat: Option<&RepeatInfo>) {
        let Some(consensus_repeat) = consensus_repeat else {
            self.bases[self.current_index as usize] = consensus_base;
            self.base_quals[self.current_index as usize] = consensus_qual;
            self.current_index += if self.build_forwards { 1 } else { -1 };
            return;
        };
        let repeat_length = consensus_repeat.repeat_length();
        let previous_read_bases = (self.current_index - consensus_repeat.index).abs();
        let remaining_repeat_bases = consensus_repeat.total_length() - previous_read_bases;
        let mut repeat_base_index = if self.build_forwards { 0 } else { repeat_length - 1 };
        let repeat_bases = consensus_repeat.bases.as_bytes();
        for _ in 0..remaining_repeat_bases {
            let repeat_base = repeat_bases[repeat_base_index as usize];
            self.bases[self.current_index as usize] = repeat_base;
            self.base_quals[self.current_index as usize] = consensus_qual;
            if self.build_forwards {
                self.current_index += 1;
                repeat_base_index += 1;
                if repeat_base_index == repeat_length {
                    repeat_base_index = 0;
                }
            } else {
                self.current_index -= 1;
                repeat_base_index -= 1;
                if repeat_base_index < 0 {
                    repeat_base_index = repeat_length - 1;
                }
            }
            if self.current_index < 0 || self.current_index >= self.bases.len() as i32 {
                break;
            }
        }
    }

    pub fn apply_read_mismatches(&self, read: &RpsP, seq_diff_info: &mut SequenceDiffInfo, consensus_repeat: Option<&RepeatInfo>, previous_repeat_length: i32) {
        if seq_diff_info.type_ != SequenceDiffType::Match {
            seq_diff_info.mismatch_penalty = calc_mismatch_penalty(seq_diff_info, consensus_repeat);
            read.m().add_mismatch_info(seq_diff_info.clone());
            if self.check_read_mismatches(read) {
                read.m().mark_mismatched();
                return;
            }
        } else {
            let hq = above_min_qual(read.current_qual());
            read.m().add_base_match(hq);
        }
        read.m().move_on_match_type(seq_diff_info);
        if let Some(consensus_repeat) = consensus_repeat {
            self.process_read_repeat(read, seq_diff_info, consensus_repeat, previous_repeat_length);
        }
    }

    pub fn process_read_repeat(&self, read: &RpsP, seq_diff_info: &SequenceDiffInfo, consensus_repeat: &RepeatInfo, previous_repeat_length: i32) {
        let remaining_read_repeat_bases = (seq_diff_info.repeat_count * consensus_repeat.repeat_length() - previous_repeat_length).max(0);
        read.m().move_next_bases(remaining_read_repeat_bases);
        let read_repeat_start = seq_diff_info.repeat_index(consensus_repeat, self.build_forwards, true);
        let read_repeat_end = seq_diff_info.repeat_index(consensus_repeat, self.build_forwards, false);
        let high_qual_repeat_bases = read.non_low_qual_base_count(read_repeat_start, read_repeat_end);
        read.m().add_base_matches(remaining_read_repeat_bases, high_qual_repeat_bases);
    }

    fn check_read_mismatches(&self, read: &ReadParseState) -> bool {
        if self.disable_mismatch_penalty {
            return false;
        }
        exceeds_read_mismatches_max(read, self.mismatch_penalty_max_read_length)
    }

    pub fn populate_bases(&self, bases: &mut [u8], index_start: i32) {
        let mut base_index = index_start;
        let n = bases.len();
        for i in 0..n {
            if base_index < 0 || base_index >= self.bases.len() as i32 {
                break;
            }
            let dest_index = if self.build_forwards { i } else { n - i - 1 };
            bases[dest_index] = self.bases[base_index as usize];
            base_index += if self.build_forwards { 1 } else { -1 };
        }
    }

    fn find_non_mismatch_repeats(&mut self) {
        let mut all_repeats = find_repeats(&self.bases);
        if all_repeats.is_empty() {
            return;
        }
        java_sort(&mut all_repeats, |a, b| (-a.count).cmp(&(-b.count)));
        // as in Java: a repeat that is kept is seen again on the next pass, now overlapping itself, and removed
        let index = 0;
        while index < all_repeats.len() {
            let repeat = all_repeats[index].clone();
            if self.repeats.iter().any(|x| x.overlaps(&repeat)) {
                all_repeats.remove(index);
            } else if self.build_forwards {
                self.repeats.push(repeat);
            } else {
                self.repeats.push(RepeatInfo::with_direction(repeat.last_index(), &repeat.bases, repeat.count, Some(Orientation::Reverse)));
            }
        }
    }

    fn trim_final_sequence(&mut self) {
        let len = self.bases.len() as i32;
        let mut current_index = if self.build_forwards { 0 } else { len - 1 };
        let mut valid_length = 0;
        while current_index >= 0 && current_index < len {
            if self.bases[current_index as usize] == 0 {
                break;
            }
            valid_length += 1;
            current_index += if self.build_forwards { 1 } else { -1 };
        }
        if valid_length == len {
            return;
        }
        let reduction = len - valid_length;
        let start_index = if self.build_forwards { 0 } else { reduction };
        let end_index = start_index + valid_length;
        self.bases = self.bases[start_index as usize..end_index as usize].to_vec();
        self.base_quals = self.base_quals[start_index as usize..end_index as usize].to_vec();
        if !self.build_forwards {
            for repeat in self.repeats.iter_mut() {
                *repeat = RepeatInfo::with_direction(repeat.index - reduction, &repeat.bases, repeat.count, Some(repeat.direction()));
            }
        }
    }
}

pub fn get_read_repeat_count(read: &ReadParseState, repeat: &str, previous_repeat_count: i32, build_forwards: bool, use_base_quals: bool) -> i32 {
    let r = read.read();
    let count_forward =
        get_repeat_count_ext(r.get_bases(), repeat, read.read_index(), build_forwards, previous_repeat_count, if use_base_quals { Some(r.get_base_quality()) } else { None }).max(0);
    previous_repeat_count + count_forward
}

fn find_consensus_repeat_count(read_repeat_counts: &[i32]) -> i32 {
    let mut repeat_frequencies: Vec<(i32, i32)> = Vec::new();
    let mut all_counts: Vec<i32> = Vec::new();
    for &read_repeat_count in read_repeat_counts {
        if read_repeat_count < READ_MISMATCH_LOW_REPEAT_COUNT {
            continue;
        }
        match repeat_frequencies.iter_mut().find(|(k, _)| *k == read_repeat_count) {
            Some((_, v)) => *v += 1,
            None => repeat_frequencies.push((read_repeat_count, 1)),
        }
        all_counts.push(read_repeat_count);
    }
    // highest frequency, then the lower repeat count: the same whatever the map order
    let mut max_repeat_freq = 0;
    let mut max_repeat_count = 0;
    for (key, value) in &repeat_frequencies {
        if *value > max_repeat_freq || (*value == max_repeat_freq && *key < max_repeat_count) {
            max_repeat_count = *key;
            max_repeat_freq = *value;
        }
    }
    if (max_repeat_freq as f64) < 0.33 * all_counts.len() as f64 {
        return median_integer(&mut all_counts);
    }
    max_repeat_count
}

/// Doubles.medianInteger
pub fn median_integer(values: &mut Vec<i32>) -> i32 {
    values.sort();
    let count = values.len();
    if count == 0 {
        return 0;
    }
    if count % 2 == 0 { (values[count / 2] + values[count / 2 - 1]) / 2 } else { values[count / 2] }
}

/// The next bases shared by most of the reads that match the consensus base, if any.
/// (The Java map is keyed by array identity; its iteration order cannot change the result: at most one key
/// has the same content as a read's bases, and a tie for the highest count is at most half of the reads.)
fn find_consensus_next_bases(consensus_read_next_bases: &[Option<Vec<u8>>]) -> Option<Vec<u8>> {
    if consensus_read_next_bases.is_empty() {
        return None;
    }
    if consensus_read_next_bases.len() == 1 {
        return consensus_read_next_bases[0].clone();
    }
    let mut consensus_next_bases: Option<&Vec<u8>> = None;
    let mut freq: Option<Vec<(&Vec<u8>, i32)>> = None;
    let mut consensus_read_next_base_match_count = 0;
    let consensus_read_count = consensus_read_next_bases.len();
    let required_match_count = 0.5 * consensus_read_count as f64;
    for read_next_bases in consensus_read_next_bases.iter().flatten() {
        if consensus_next_bases.is_none() {
            consensus_next_bases = Some(read_next_bases);
            consensus_read_next_base_match_count = 1;
        } else if freq.is_none() {
            if bases_match(read_next_bases, consensus_next_bases.unwrap()) {
                consensus_read_next_base_match_count += 1;
                if consensus_read_next_base_match_count as f64 > required_match_count {
                    return consensus_next_bases.cloned();
                }
            } else {
                freq = Some(vec![(consensus_next_bases.unwrap(), consensus_read_next_base_match_count), (read_next_bases, 1)]);
            }
        } else {
            let f = freq.as_mut().unwrap();
            let mut found = false;
            for entry in f.iter_mut() {
                if bases_match(read_next_bases, entry.0) {
                    let new_count = entry.1 + 1;
                    if new_count as f64 > required_match_count {
                        // as in Java, the first bases seen are returned, not this entry's
                        return consensus_next_bases.cloned();
                    }
                    entry.1 = new_count;
                    found = true;
                    break;
                }
            }
            if !found {
                f.push((read_next_bases, 1));
            }
        }
    }
    if let Some(f) = &freq {
        for entry in f {
            if entry.1 > consensus_read_next_base_match_count {
                consensus_read_next_base_match_count = entry.1;
                consensus_next_bases = Some(entry.0);
            }
        }
    }
    let matched_fraction = consensus_read_next_base_match_count as f64 / consensus_read_count as f64;
    if matched_fraction > 0.5 { consensus_next_bases.cloned() } else { None }
}

pub fn exceeds_read_mismatches(read: &ReadParseState) -> bool {
    exceeds_read_mismatches_max(read, -1)
}

fn exceeds_read_mismatches_max(read: &ReadParseState, max_read_length: i32) -> bool {
    let mismatch_penalty = read.mismatch_penalty();
    let read_overlap_length = read.overlap_base_count();
    let capped_read_overlap = if max_read_length > 0 { max_read_length.min(read_overlap_length) } else { read_overlap_length };
    let permitted_penalty = permitted_read_mismatches(capped_read_overlap);
    if mismatch_penalty > permitted_penalty {
        return true;
    }
    let evaluated_bases = read.evaluated_base_count();
    if evaluated_bases >= READ_MISMATCH_RATE_MIN_BASE {
        let mismatch_count = read.mismatch_count(false);
        if mismatch_count as f64 >= READ_MAX_MISMATCH_RATE * evaluated_bases as f64 {
            return true;
        }
    }
    false
}

pub fn permitted_read_mismatches(read_overlap: i32) -> f64 {
    if read_overlap <= READ_MISMATCH_PENALTY_LENGTH_1 {
        return READ_MISMATCH_PENALTY_PENALTY_1;
    }
    if read_overlap <= READ_MISMATCH_PENALTY_LENGTH_2 {
        return READ_MISMATCH_PENALTY_PENALTY_2;
    }
    if read_overlap <= READ_MISMATCH_PENALTY_LENGTH_3 {
        return READ_MISMATCH_PENALTY_PENALTY_3;
    }
    if read_overlap <= READ_MISMATCH_PENALTY_LENGTH_4 {
        return READ_MISMATCH_PENALTY_PENALTY_4;
    }
    READ_MISMATCH_PENALTY_PENALTY_LONG
}

pub fn calc_mismatch_penalty(seq_diff_info: &SequenceDiffInfo, consensus_repeat: Option<&RepeatInfo>) -> f64 {
    use SequenceDiffType::*;
    if seq_diff_info.type_ == Match {
        return 0.0;
    }
    if seq_diff_info.qual_type == BaseQualType::Low && seq_diff_info.type_ != Repeat {
        return 0.0;
    }
    if seq_diff_info.type_ == Base {
        return if seq_diff_info.qual_type == BaseQualType::Medium { READ_MISMATCH_MED_QUAL_SNV_PENALTY } else { READ_MISMATCH_PENALTY };
    }
    if seq_diff_info.type_ == Insert || seq_diff_info.type_ == Delete {
        return if seq_diff_info.qual_type == BaseQualType::Medium { READ_MISMATCH_MED_QUAL_NON_SNV_PENALTY } else { READ_MISMATCH_PENALTY };
    }
    let Some(consensus_repeat) = consensus_repeat else { return READ_MISMATCH_PENALTY };
    let repeat_count = consensus_repeat.count;
    let repeat_count_diff = (seq_diff_info.repeat_count - consensus_repeat.count).abs();
    if seq_diff_info.qual_type == BaseQualType::Low {
        if repeat_count >= read_mismatch_long_repeat_count() && seq_diff_info.repeat_count < consensus_repeat.count - read_mismatch_max_repeat_diff() {
            return (consensus_repeat.count - seq_diff_info.repeat_count) as f64 * (1.0 / read_mismatch_max_repeat_diff() as f64);
        }
        return 0.0;
    }
    let penalty_count = if repeat_count < read_mismatch_medium_repeat_count() {
        (repeat_count_diff - 1).max(0)
    } else if repeat_count < read_mismatch_long_repeat_count() {
        (repeat_count_diff - 2).max(0)
    } else {
        0
    };
    if seq_diff_info.qual_type == BaseQualType::Medium {
        return penalty_count as f64 * READ_MISMATCH_MED_QUAL_NON_SNV_PENALTY;
    }
    penalty_count as f64
}
