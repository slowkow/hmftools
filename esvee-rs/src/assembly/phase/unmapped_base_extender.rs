//! UnmappedBaseExtender: extends an assembly's extension bases with unmapped mates and remote reads that
//! overlap the existing extension sequence.

use std::cmp::Ordering;
use std::collections::HashSet;

use crate::common::compat;
use crate::assembly::assembly_utils::{self, mismatches_per_comparison_length, reverse_complement_bases, swap_dna_base, DNA_BASE_BYTES, DNA_BASE_COUNT};
use crate::assembly::common_utils::below_min_qual;
use crate::assembly::constants::*;
use crate::assembly::globals::assembly_build_debug;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_utils::calc_trimmed_read_base_length;
use crate::assembly::sequence_compare::{compare_sequences, index_of};
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::repeat_info::{find_repeats, RepeatInfo};
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::SupportType;
use crate::common::java_sort::java_sort;

pub struct UnmappedBaseExtender {
    junction_assembly: AsmP,
    is_forward_junction: bool,
    bases: Vec<u8>,
    base_quals: Vec<u8>,
    repeats: Vec<RepeatInfo>,
    support_reads: Vec<SupP>,
    /// the bases at the start of the current pass, which reads are searched against
    extension_bases: Vec<u8>,
    ref_base_count: i32,
}

#[derive(Clone)]
struct ReadSequenceMatch {
    read: ReadP,
    read_seq_start: i32,
    extension_base_seq_start: i32,
    overlap: i32,
    mismatch_penalty: f64,
    consensus_mismatches: i32,
}

fn reverse_read_bases(read: &Read) -> bool {
    read.orientation() == read.mate_orientation()
}

fn permitted_read_mismatches(read_base_overlap: i32) -> i32 {
    mismatches_per_comparison_length(read_base_overlap) + 1
}

fn index_of_from(haystack: &[u8], needle: &[u8], from: i32) -> i32 {
    let from = from.max(0) as usize;
    if from > haystack.len() {
        return -1;
    }
    index_of(&haystack[from..], needle).map(|i| i + from as i32).unwrap_or(-1)
}

impl UnmappedBaseExtender {
    pub fn new(assembly: AsmP) -> UnmappedBaseExtender {
        let is_forward_junction = assembly.junction().is_forward();
        let ref_base_count = ASSEMBLY_LINK_OVERLAP_BASES.min(assembly.ref_base_length());
        let (extension_index_start, extension_index_end) =
            if is_forward_junction { (assembly.junction_index() - ref_base_count + 1, assembly.base_length() - 1) } else { (0, assembly.junction_index() + ref_base_count - 1) };
        UnmappedBaseExtender {
            junction_assembly: assembly,
            is_forward_junction,
            bases: assembly.bases()[extension_index_start as usize..=extension_index_end as usize].to_vec(),
            base_quals: assembly.base_quals()[extension_index_start as usize..=extension_index_end as usize].to_vec(),
            repeats: Vec::new(),
            support_reads: Vec::new(),
            extension_bases: Vec::new(),
            ref_base_count,
        }
    }

    pub fn extension_bases(&self) -> Vec<u8> {
        let n = self.bases.len();
        let r = self.ref_base_count as usize;
        if self.is_forward_junction { self.bases[r..n].to_vec() } else { self.bases[0..n - r].to_vec() }
    }
    pub fn extension_base_length(&self) -> i32 {
        self.bases.len() as i32 - self.ref_base_count
    }
    pub fn base_qualities(&self) -> Vec<u8> {
        let n = self.base_quals.len();
        let r = self.ref_base_count as usize;
        if self.is_forward_junction { self.base_quals[r..n].to_vec() } else { self.base_quals[0..n - r].to_vec() }
    }
    pub fn support_reads(&self) -> &Vec<SupP> {
        &self.support_reads
    }

    /// reads that are matched are removed from the list
    pub fn process_reads(&mut self, reads: &mut Vec<ReadP>) {
        let mut iteration = 0;
        loop {
            if !self.try_add_reads(reads, iteration) {
                break;
            }
            iteration += 1;
        }
    }

    fn max_read_extension(&self, m: &ReadSequenceMatch) -> i32 {
        if self.is_forward_junction { m.extension_base_seq_start + m.read.bases_length() - m.read_seq_start } else { self.bases.len() as i32 - m.extension_base_seq_start + m.read_seq_start }
    }

    pub fn try_add_reads(&mut self, reads: &mut Vec<ReadP>, iteration: i32) -> bool {
        self.extension_bases = self.bases.clone();
        self.repeats = find_repeats(&self.bases);
        let mut read_sequence_matches: Vec<ReadSequenceMatch> = Vec::new();
        let mut read_index = 0;
        while read_index < reads.len() {
            let read = reads[read_index].clone();
            match self.find_read_sequence_overlap(&read) {
                Some(m) => {
                    if self.max_read_extension(&m) > 0 {
                        read_sequence_matches.push(m);
                    }
                    reads.remove(read_index);
                }
                None => read_index += 1,
            }
        }
        if read_sequence_matches.is_empty() {
            return false;
        }
        if iteration == 0 && assembly_build_debug() {
            crate::debug!("junc({}) initial sequence {}", self.junction_assembly.junction().coords(), self.bases.iter().map(|b| *b as char).collect::<String>());
        }
        self.check_extension_consensus(&mut read_sequence_matches);
        let new_extension_length = read_sequence_matches.iter().map(|x| self.max_read_extension(x)).max().unwrap_or(0);
        let mut base_offset = 0;
        if new_extension_length > self.bases.len() as i32 {
            if !self.is_forward_junction {
                base_offset = new_extension_length - self.bases.len() as i32;
            }
            self.extend_bases(new_extension_length);
        }
        let mut added_reads = 0;
        java_sort(&mut read_sequence_matches, |a, b| {
            if a.consensus_mismatches != b.consensus_mismatches {
                return if a.consensus_mismatches < b.consensus_mismatches { Ordering::Less } else { Ordering::Greater };
            }
            if a.overlap != b.overlap {
                return if a.overlap > b.overlap { Ordering::Less } else { Ordering::Greater };
            }
            if a.mismatch_penalty != b.mismatch_penalty {
                return if a.mismatch_penalty < b.mismatch_penalty { Ordering::Less } else { Ordering::Greater };
            }
            Ordering::Equal
        });
        for m in &read_sequence_matches {
            if self.add_read(m, base_offset) {
                added_reads += 1;
            }
        }
        self.trim_unset_bases();
        if added_reads > 0 && assembly_build_debug() {
            crate::debug!(
                "junc({}) iteration({}) added {} unmapped reads, new sequence {}",
                self.junction_assembly.junction().coords(), iteration, added_reads, self.bases.iter().map(|b| *b as char).collect::<String>()
            );
        }
        added_reads > 0
    }

    fn extend_bases(&mut self, new_extension_length: i32) {
        let existing_bases = std::mem::take(&mut self.bases);
        let existing_quals = std::mem::take(&mut self.base_quals);
        let n = new_extension_length as usize;
        let base_offset = n - existing_bases.len();
        self.bases = vec![0u8; n];
        self.base_quals = vec![0u8; n];
        if self.is_forward_junction {
            self.bases[..existing_bases.len()].copy_from_slice(&existing_bases);
            self.base_quals[..existing_quals.len()].copy_from_slice(&existing_quals);
        } else {
            self.bases[base_offset..].copy_from_slice(&existing_bases);
            self.base_quals[base_offset..].copy_from_slice(&existing_quals);
        }
    }

    fn trim_unset_bases(&mut self) {
        let is_forward = self.junction_assembly.is_forward_junction();
        let len = self.bases.len() as i32;
        let mut current_index = if is_forward { 0 } else { len - 1 };
        let mut valid_length = 0;
        while current_index >= 0 && current_index < len {
            if self.bases[current_index as usize] == 0 {
                break;
            }
            valid_length += 1;
            current_index += if is_forward { 1 } else { -1 };
        }
        if valid_length == len {
            return;
        }
        let reduction = len - valid_length;
        let start_index = if is_forward { 0 } else { reduction } as usize;
        let end_index = start_index + valid_length as usize;
        self.bases = self.bases[start_index..end_index].to_vec();
        self.base_quals = self.base_quals[start_index..end_index].to_vec();
    }

    fn check_extension_consensus(&self, read_sequence_matches: &mut [ReadSequenceMatch]) {
        if read_sequence_matches.len() <= 2 {
            return;
        }
        let read_count = read_sequence_matches.len();
        let mut read_index_starts = vec![0i32; read_count];
        let mut read_index_ends = vec![0i32; read_count];
        let mut max_extension_length = -1;
        for i in 0..read_count {
            let m = &read_sequence_matches[i];
            if self.is_forward_junction {
                read_index_starts[i] = m.read_seq_start + m.overlap;
                read_index_ends[i] = m.read.bases_length() - 1;
            } else {
                read_index_starts[i] = 0;
                read_index_ends[i] = m.read_seq_start - 1;
            }
            let read_extension = read_index_ends[i] - read_index_starts[i] + 1;
            max_extension_length = max_extension_length.max(read_extension);
        }
        let read_base_at = |m: &ReadSequenceMatch, i: usize, base_index: i32| -> Option<u8> {
            let read = &m.read;
            let mut read_base_index = if self.is_forward_junction { read_index_starts[i] + base_index } else { read_index_ends[i] - base_index };
            if read_base_index < 0 || read_base_index >= read.bases_length() {
                return None;
            }
            if reverse_read_bases(read) {
                read_base_index = read.bases_length() - read_base_index - 1;
                Some(swap_dna_base(read.get_bases()[read_base_index as usize]))
            } else {
                Some(read.get_bases()[read_base_index as usize])
            }
        };
        let java_a23 = compat::java("A23");
        for base_index in 0..max_extension_length {
            let mut consensus_base: u8 = 0;
            let mut base_counts: Option<[i32; DNA_BASE_COUNT]> = None;
            // the reads so far that have a base here, all of them the first base seen
            let mut agreeing_reads = 0;
            for i in 0..read_count {
                let Some(read_base) = read_base_at(&read_sequence_matches[i], i, base_index) else { continue };
                if base_counts.is_none() {
                    if consensus_base == 0 {
                        consensus_base = read_base;
                        agreeing_reads += 1;
                        continue;
                    } else if consensus_base == read_base {
                        agreeing_reads += 1;
                        continue;
                    }
                    let mut counts = [0i32; DNA_BASE_COUNT];
                    // Java: the count for the first base is the index of this read, which also counts the reads
                    // that do not reach this base
                    counts[assembly_utils::base_index(consensus_base)] = if java_a23 { i as i32 } else { agreeing_reads };
                    base_counts = Some(counts);
                }
                base_counts.as_mut().unwrap()[assembly_utils::base_index(read_base)] += 1;
            }
            let Some(base_counts) = base_counts else { continue };
            let mut max_count = 0;
            let mut max_base: u8 = 0;
            for b in 0..DNA_BASE_COUNT {
                if base_counts[b] > max_count {
                    max_count = base_counts[b];
                    max_base = DNA_BASE_BYTES[b];
                }
            }
            for i in 0..read_count {
                let Some(read_base) = read_base_at(&read_sequence_matches[i], i, base_index) else { continue };
                if read_base != max_base {
                    read_sequence_matches[i].consensus_mismatches += 1;
                }
            }
        }
    }

    fn add_read(&mut self, m: &ReadSequenceMatch, base_offset: i32) -> bool {
        let mut mismatch_count = 0;
        let mut read_index_start;
        let read_index_end;
        let mut junction_read_start_distance;
        let mut ext_base_index_start;
        let read = &m.read;
        if self.is_forward_junction {
            read_index_start = m.read_seq_start;
            read_index_end = read.bases_length() - 1;
            ext_base_index_start = m.extension_base_seq_start;
            junction_read_start_distance = -m.extension_base_seq_start - m.read_seq_start;
            junction_read_start_distance += self.ref_base_count;
        } else {
            read_index_start = 0;
            read_index_end = m.read_seq_start - 1;
            ext_base_index_start = base_offset + m.extension_base_seq_start - m.read_seq_start;
            if ext_base_index_start < 0 {
                read_index_start += ext_base_index_start.abs();
                ext_base_index_start = 0;
            }
            junction_read_start_distance = self.bases.len() as i32 - m.extension_base_seq_start - m.read_seq_start;
            junction_read_start_distance += base_offset;
            junction_read_start_distance -= self.ref_base_count;
        }
        let reverse_bases = reverse_read_bases(read);
        let read_base_length = read.bases_length();
        let read_bases = read.get_bases();
        let read_quals = read.get_base_quality();
        let len = self.bases.len() as i32;
        let mut ext_base_index = ext_base_index_start;
        let mut i = read_index_start;
        while i <= read_index_end {
            let read_index = (if reverse_bases { read_base_length - 1 - i } else { i }) as usize;
            let mut base = read_bases[read_index];
            if reverse_bases {
                base = swap_dna_base(base);
            }
            let qual = read_quals[read_index];
            if ext_base_index >= len {
                break;
            }
            let e = ext_base_index as usize;
            if !(self.bases[e] == 0 || self.bases[e] == base || below_min_qual(qual) || below_min_qual(self.base_quals[e])) {
                mismatch_count += 1;
            }
            i += 1;
            ext_base_index += 1;
        }
        let read_base_overlap = read_index_end - read_index_start + 1;
        let permitted_mismatches = permitted_read_mismatches(read_base_overlap);
        if mismatch_count > permitted_mismatches {
            return false;
        }
        mismatch_count = 0;
        ext_base_index = ext_base_index_start;
        i = read_index_start;
        while i <= read_index_end {
            let read_index = (if reverse_bases { read_base_length - 1 - i } else { i }) as usize;
            let mut base = read_bases[read_index];
            if reverse_bases {
                base = swap_dna_base(base);
            }
            let qual = read_quals[read_index];
            if ext_base_index >= len {
                break;
            }
            let e = ext_base_index as usize;
            if self.bases[e] == 0 {
                self.bases[e] = base;
                self.base_quals[e] = qual;
            } else if self.bases[e] == base || below_min_qual(qual) {
                if qual > self.base_quals[e] {
                    self.base_quals[e] = qual;
                }
            } else if below_min_qual(self.base_quals[e]) {
                self.bases[e] = base;
                self.base_quals[e] = qual;
            } else {
                mismatch_count += 1;
            }
            i += 1;
            ext_base_index += 1;
        }
        let matched_count = read_base_overlap - mismatch_count;
        let support_read = SupportRead::new(read, SupportType::Extension, junction_read_start_distance, matched_count, mismatch_count);
        self.support_reads.push(support_read);
        true
    }

    fn find_read_sequence_overlap(&self, read: &ReadP) -> Option<ReadSequenceMatch> {
        let subsequence_length = MATCH_SUBSEQUENCE_LENGTH;
        Read::trim_low_qual_bases(read);
        let read_base_count = read.get_bases().len() as i32;
        let reverse_bases = reverse_read_bases(read);
        let read_bases: Vec<u8> = if reverse_bases { reverse_complement_bases(read.get_bases()) } else { read.get_bases().to_vec() };
        let mut read_base_quals: Option<Vec<u8>> = None;
        let mut read_sequence_match: Option<ReadSequenceMatch> = None;
        let mut tested_matches: HashSet<(i32, i32, i32)> = HashSet::new();
        let mut read_index = 0;
        while read_index + subsequence_length < read.bases_length() {
            let read_seq_bases = &read_bases[read_index as usize..(read_index + subsequence_length) as usize];
            let mut ext_base_match_index = index_of_from(&self.extension_bases, read_seq_bases, 0);
            while ext_base_match_index >= 0 {
                if read_base_quals.is_none() {
                    read_base_quals = Some(if reverse_bases { read.get_base_quality().iter().rev().copied().collect() } else { read.get_base_quality().to_vec() });
                }
                let new_match = self.find_read_sequence_match(read, &read_bases, read_base_quals.as_ref().unwrap(), read_index, read_base_count, ext_base_match_index, &mut tested_matches);
                if let Some(nm) = new_match {
                    if read_sequence_match.is_some() {
                        return None;
                    }
                    read_sequence_match = Some(nm);
                }
                ext_base_match_index = index_of_from(&self.extension_bases, read_seq_bases, ext_base_match_index + subsequence_length);
            }
            read_index += subsequence_length;
        }
        read_sequence_match
    }

    #[allow(clippy::too_many_arguments)]
    fn find_read_sequence_match(
        &self, read: &ReadP, read_bases: &[u8], read_base_quals: &[u8], read_index: i32, read_base_count: i32, ext_base_match_index: i32, tested_matches: &mut HashSet<(i32, i32, i32)>,
    ) -> Option<ReadSequenceMatch> {
        let read_lower_base_count = read_index;
        let read_upper_base_count = read_base_count - read_index;
        let ext_base_lower_base_count = ext_base_match_index;
        let ext_base_upper_base_count = self.bases.len() as i32 - ext_base_match_index;
        let lower_overlap = read_lower_base_count.min(ext_base_lower_base_count);
        let upper_overlap = read_upper_base_count.min(ext_base_upper_base_count);
        let total_overlap = lower_overlap + upper_overlap;
        let read_index_start = read_index - lower_overlap;
        let read_index_end = read_index + upper_overlap - 1;
        if total_overlap < ASSEMBLY_LINK_OVERLAP_BASES {
            let trimmed_read_length = calc_trimmed_read_base_length(read, read_index_start, read_index_end);
            if trimmed_read_length < ASSEMBLY_READ_TRIMMED_OVERLAP_BASES {
                return None;
            }
        }
        let ext_base_index_start = ext_base_match_index - lower_overlap;
        let ext_base_index_end = ext_base_match_index + upper_overlap - 1;
        if !tested_matches.insert((read_index_start, ext_base_index_start, total_overlap)) {
            return None;
        }
        let permitted_mismatches = permitted_read_mismatches(total_overlap);
        let mismatch_penalty = compare_sequences(
            &self.bases, &self.base_quals, ext_base_index_start, ext_base_index_end, &self.repeats,
            read_bases, read_base_quals, read_index_start, read_index_end, &[], permitted_mismatches as f64,
        );
        if mismatch_penalty <= permitted_mismatches as f64 {
            Some(ReadSequenceMatch { read: read.clone(), read_seq_start: read_index_start, extension_base_seq_start: ext_base_index_start, overlap: total_overlap, mismatch_penalty, consensus_mismatches: 0 })
        } else {
            None
        }
    }
}
