//! ReadParseState: a cursor over a read's bases (and optionally its CIGAR) with its matches and mismatches
//! against a consensus sequence being built.

use std::cell::OnceCell;

use crate::assembly::assembly_utils::NO_BASE;
use crate::assembly::common_utils::{above_min_qual, below_min_qual, is_medium_base_qual};
use crate::assembly::globals::is_ultima;
use crate::assembly::ptr::A;
use crate::assembly::read::read::ReadP;
use crate::assembly::sam::{CigarElement, CigarOp, SamRecord, BASE_QUAL_TAG_ATTRIBUTE};
use crate::assembly::sequence_builder::calc_mismatch_penalty;
use crate::assembly::sequence_diff::{BaseQualType, SequenceDiffInfo, SequenceDiffType};
use crate::common::fmt::java_fixed;

pub type RpsP = A<ReadParseState>;

pub struct ReadParseState {
    read: ReadP,
    move_forward: bool,
    base_length: i32,
    start_index: i32,
    read_index: i32,
    exhausted: bool,
    valid: bool,
    ref_position: i32,
    element_count: i32,
    /// the read's CIGAR elements when the CIGAR is tracked
    elements: Option<Vec<CigarElement>>,
    cigar_index: i32,
    element: Option<CigarElement>,
    element_index: i32,
    mismatched: bool,
    pub base_matches: i32,
    pub high_qual_matches: i32,
    mismatches: Option<A<Vec<SequenceDiffInfo>>>,
    low_qual_indices: OnceCell<Vec<i32>>,
}

/// UltimaBamUtils.extractLowQualIndices
pub fn ultima_extract_low_qual_indices(record: &SamRecord) -> Vec<i32> {
    let Some(qual_tag) = record.get_string_attribute(BASE_QUAL_TAG_ATTRIBUTE) else { return Vec::new() };
    let Some((_, indices)) = qual_tag.split_once('=') else { panic!("invalid UQ tag: {qual_tag}") };
    let mut values = Vec::new();
    // String.split drops trailing empty items
    let mut items: Vec<&str> = indices.split(',').collect();
    while items.last().is_some_and(|x| x.is_empty()) {
        items.pop();
    }
    for item in items {
        if let Some((s, e)) = item.split_once('-') {
            let (s, e): (i32, i32) = (s.parse().expect("UQ index"), e.parse().expect("UQ index"));
            for i in s..=e {
                values.push(i);
            }
        } else {
            values.push(item.parse().expect("UQ index"));
        }
    }
    values
}

impl ReadParseState {
    pub fn new(move_forward: bool, read: &ReadP, start_index: i32) -> RpsP {
        Self::new_cigar(move_forward, read, start_index, false)
    }

    pub fn new_cigar(move_forward: bool, read: &ReadP, start_index: i32, track_cigar: bool) -> RpsP {
        let elements = if track_cigar { Some(read.cigar_elements().to_vec()) } else { None };
        let element_count = elements.as_ref().map(|e| e.len() as i32).unwrap_or(0);
        let mut s = ReadParseState {
            read: read.clone(),
            move_forward,
            base_length: read.bases_length(),
            start_index,
            read_index: 0,
            exhausted: false,
            valid: true,
            ref_position: 0,
            element_count,
            elements,
            cigar_index: 0,
            element: None,
            element_index: 0,
            mismatched: false,
            base_matches: 0,
            high_qual_matches: 0,
            mismatches: None,
            low_qual_indices: OnceCell::new(),
        };
        s.reset_index();
        A::new(s)
    }

    pub fn read(&self) -> &ReadP {
        &self.read
    }
    pub fn mismatched(&self) -> bool {
        self.mismatched
    }
    pub fn mark_mismatched(&mut self) {
        self.mismatched = true;
    }
    pub fn is_valid(&self) -> bool {
        self.valid
    }
    pub fn mark_invalid(&mut self) {
        self.valid = false;
    }
    pub fn exhausted(&self) -> bool {
        self.exhausted
    }
    pub fn start_index(&self) -> i32 {
        self.start_index
    }
    #[inline]
    pub fn read_index(&self) -> i32 {
        self.read_index
    }
    #[inline]
    pub fn current_base(&self) -> u8 {
        self.read.get_bases()[self.read_index as usize]
    }
    #[inline]
    pub fn current_qual(&self) -> u8 {
        self.read.get_base_quality()[self.read_index as usize]
    }
    pub fn ref_position(&self) -> i32 {
        self.ref_position
    }
    pub fn operator(&self) -> Option<CigarOp> {
        self.element.map(|e| e.op)
    }
    pub fn element_length(&self) -> i32 {
        self.element.map(|e| e.len).unwrap_or(0)
    }
    pub fn base_length(&self) -> i32 {
        self.base_length
    }
    pub fn overlap_base_count(&self) -> i32 {
        if self.move_forward { self.base_length - self.start_index } else { self.start_index + 1 }
    }
    pub fn evaluated_base_count(&self) -> i32 {
        (if self.move_forward { self.read_index.min(self.base_length - 1) - self.start_index } else { self.start_index - self.read_index.max(0) }) + 1
    }

    pub fn move_on_match_type(&mut self, seq_diff_info: &SequenceDiffInfo) {
        if seq_diff_info.repeat_count != 0 {
            return;
        }
        match seq_diff_info.type_ {
            SequenceDiffType::Match | SequenceDiffType::Base => self.move_next(),
            SequenceDiffType::Delete => {}
            SequenceDiffType::Insert => self.move_next_bases(seq_diff_info.indel_length + 1),
            _ => {}
        }
    }

    pub fn move_next_bases(&mut self, base_count: i32) {
        if self.exhausted {
            return;
        }
        for _ in 0..base_count {
            self.move_next();
        }
    }

    pub fn move_next(&mut self) {
        if self.exhausted {
            return;
        }
        if self.elements.is_some() {
            self.move_next_cigar_aware();
        } else {
            self.move_next_read_index();
        }
    }

    fn move_next_read_index(&mut self) {
        if self.move_forward {
            self.read_index += 1;
        } else {
            self.read_index -= 1;
        }
        if self.read_index < 0 || self.read_index >= self.base_length {
            self.exhausted = true;
        }
    }

    fn move_next_cigar_aware(&mut self) {
        let elements = self.elements.as_ref().unwrap();
        let element = self.element.expect("cigar element");
        if self.move_forward {
            self.element_index += 1;
            if self.element_index >= element.len {
                let was_clipping = element.op.is_clipping();
                self.cigar_index += 1;
                if self.cigar_index >= self.element_count {
                    self.exhausted = true;
                    return;
                }
                self.element_index = 0;
                let element = elements[self.cigar_index as usize];
                self.element = Some(element);
                if self.cigar_index > 0 && element.op.is_clipping() {
                    self.exhausted = true;
                } else {
                    if !was_clipping && element.op != CigarOp::I {
                        self.ref_position += 1;
                    }
                    if element.op.consumes_read_bases() {
                        self.read_index += 1;
                    }
                }
            } else {
                if element.op.consumes_reference_bases() {
                    self.ref_position += 1;
                }
                if element.op.consumes_read_bases() {
                    self.read_index += 1;
                }
            }
        } else {
            self.element_index -= 1;
            if self.element_index < 0 || element.op == CigarOp::H {
                let was_clipping = element.op.is_clipping();
                self.cigar_index -= 1;
                if self.cigar_index < 0 {
                    self.exhausted = true;
                    return;
                }
                let element = elements[self.cigar_index as usize];
                self.element = Some(element);
                self.element_index = element.len - 1;
                if self.cigar_index < self.element_count - 1 && element.op.is_clipping() {
                    self.exhausted = true;
                } else {
                    if !was_clipping && element.op != CigarOp::I {
                        self.ref_position -= 1;
                    }
                    if element.op.consumes_read_bases() {
                        self.read_index -= 1;
                    }
                }
            } else {
                if element.op.consumes_reference_bases() {
                    self.ref_position -= 1;
                }
                if element.op.consumes_read_bases() {
                    self.read_index -= 1;
                }
            }
        }
        if self.read_index < 0 || self.read_index >= self.read.get_bases().len() as i32 {
            self.exhausted = true;
        }
    }

    pub fn skip_insert(&mut self) {
        while !self.exhausted() && self.operator() == Some(CigarOp::I) {
            self.move_next();
        }
    }

    pub fn reset_index(&mut self) {
        self.exhausted = false;
        self.read_index = if self.move_forward { 0 } else { self.base_length - 1 };
        if let Some(elements) = &self.elements {
            if self.move_forward {
                self.ref_position = self.read.alignment_start();
                self.cigar_index = 0;
                self.element = Some(elements[0]);
                self.element_index = 0;
            } else {
                self.ref_position = self.read.alignment_end();
                self.cigar_index = self.element_count - 1;
                let e = elements[self.cigar_index as usize];
                self.element = Some(e);
                self.element_index = e.len - 1;
            }
            self.check_skip_hard_clips();
        }
        if self.start_index >= 0 {
            while self.read_index != self.start_index && !self.exhausted {
                self.move_next();
            }
        }
    }

    fn get_adjacent_base(&self, get_upper: bool) -> u8 {
        let bases = self.read.get_bases();
        if self.move_forward == get_upper {
            if self.read_index >= self.base_length - 1 { NO_BASE } else { bases[(self.read_index + 1) as usize] }
        } else if self.read_index < 1 {
            NO_BASE
        } else {
            bases[(self.read_index - 1) as usize]
        }
    }
    pub fn get_next_bases(&self, base_count: i32) -> Option<Vec<u8>> {
        self.get_adjacent_bases(base_count, self.move_forward)
    }
    pub fn get_previous_bases(&self, base_count: i32) -> Option<Vec<u8>> {
        self.get_adjacent_bases(base_count, !self.move_forward)
    }
    pub fn get_previous_base(&self) -> u8 {
        self.get_adjacent_base(false)
    }
    pub fn get_next_base(&self) -> u8 {
        self.get_adjacent_base(true)
    }
    fn get_adjacent_bases(&self, base_count: i32, get_upper: bool) -> Option<Vec<u8>> {
        let mut bases = vec![0u8; base_count as usize];
        let read_bases = self.read.get_bases();
        let mut read_index = self.read_index + if get_upper { 1 } else { -1 };
        if get_upper {
            for b in bases.iter_mut() {
                if read_index >= self.base_length {
                    return None;
                }
                *b = read_bases[read_index as usize];
                read_index += 1;
            }
        } else {
            for i in (0..base_count as usize).rev() {
                if read_index < 0 {
                    return None;
                }
                bases[i] = read_bases[read_index as usize];
                read_index -= 1;
            }
        }
        Some(bases)
    }

    pub fn non_low_qual_base_count(&self, read_index_start: i32, read_index_end: i32) -> i32 {
        let quals = self.read.get_base_quality();
        let mut count = 0;
        let mut i = read_index_start;
        while i <= read_index_end {
            if i >= 0 && (i as usize) < quals.len() && above_min_qual(quals[i as usize]) {
                count += 1;
            }
            i += 1;
        }
        count
    }

    pub fn qual_type(&self, index: i32) -> BaseQualType {
        self.range_min_qual_type(index, index)
    }

    pub fn straddling_qual_type(&self, read_index_start: i32, read_index_end: i32) -> BaseQualType {
        let start = self.qual_type(read_index_start);
        let end = self.qual_type(read_index_end);
        if start == end {
            return start;
        }
        if start > end { start } else { end }
    }

    pub fn range_min_qual_type(&self, read_index_start: i32, read_index_end: i32) -> BaseQualType {
        if read_index_start < 0 {
            return BaseQualType::Low;
        }
        if is_ultima() {
            let low = self.low_qual_indices.get_or_init(|| ultima_extract_low_qual_indices(self.read.bam_record()));
            let mut i = read_index_start;
            while i <= read_index_end {
                if low.contains(&i) {
                    return BaseQualType::Low;
                }
                i += 1;
            }
            BaseQualType::High
        } else {
            let base_quals = self.read.get_base_quality();
            if read_index_end >= base_quals.len() as i32 {
                return BaseQualType::Low;
            }
            let mut has_medium = false;
            let mut i = read_index_start;
            while i <= read_index_end {
                let q = base_quals[i as usize];
                if below_min_qual(q) {
                    return BaseQualType::Low;
                }
                has_medium |= is_medium_base_qual(q);
                i += 1;
            }
            if has_medium { BaseQualType::Medium } else { BaseQualType::High }
        }
    }

    pub fn matched_bases(&self) -> i32 {
        self.base_matches
    }
    pub fn add_base_match(&mut self, above_min_qual: bool) {
        self.base_matches += 1;
        if above_min_qual {
            self.high_qual_matches += 1;
        }
    }
    pub fn add_base_matches(&mut self, count: i32, above_min_qual_count: i32) {
        self.base_matches += count;
        self.high_qual_matches += above_min_qual_count;
    }
    pub fn high_qual_matches(&self) -> i32 {
        self.high_qual_matches
    }
    pub fn mismatches(&self) -> &[SequenceDiffInfo] {
        match &self.mismatches {
            Some(m) => m,
            None => &[],
        }
    }
    /// The list object itself, which Java hands to the read (Read.setExtensionMismatches) so later changes
    /// are seen through both; an immutable empty list when there is none yet.
    pub fn mismatches_shared(&self) -> A<Vec<SequenceDiffInfo>> {
        match &self.mismatches {
            Some(m) => m.clone(),
            None => A::new(Vec::new()),
        }
    }
    pub fn mismatch_count(&self, exclude_low_qual: bool) -> i32 {
        match &self.mismatches {
            None => 0,
            Some(m) => {
                if exclude_low_qual {
                    m.iter().filter(|x| x.mismatch_penalty != 0.0).count() as i32
                } else {
                    m.len() as i32
                }
            }
        }
    }
    pub fn add_mismatch_info(&mut self, seq_diff_info: SequenceDiffInfo) {
        if self.mismatches.is_none() {
            self.mismatches = Some(A::new(Vec::new()));
        }
        self.mismatches.as_ref().unwrap().m().push(seq_diff_info);
    }
    pub fn add_mismatch(&mut self, consensus_index: i32, type_: SequenceDiffType) {
        let mut seq_diff_info = if type_ == SequenceDiffType::Base {
            SequenceDiffInfo::from_snv(self, consensus_index)
        } else {
            SequenceDiffInfo::new(self.read_index(), consensus_index, Some(String::new()), type_, self.qual_type(self.read_index()))
        };
        seq_diff_info.mismatch_penalty = calc_mismatch_penalty(&seq_diff_info, None);
        self.add_mismatch_info(seq_diff_info);
    }
    pub fn reset_matches(&mut self) {
        self.high_qual_matches = 0;
        self.base_matches = 0;
        if let Some(m) = &self.mismatches {
            m.m().clear();
        }
    }
    pub fn mismatch_penalty(&self) -> f64 {
        match &self.mismatches {
            None => 0.0,
            // DoubleStream.sum() in Java, which is a compensated sum
            Some(m) => crate::common::java_math::double_stream_sum(m.iter().map(|x| x.mismatch_penalty)),
        }
    }

    pub fn mismatch_info(&self) -> String {
        let mismatches = self.mismatches();
        if mismatches.is_empty() {
            return "EXACT".to_string();
        }
        if self.mismatched {
            return "MISMATCHED".to_string();
        }
        let mut sj: Vec<String> = vec![format!("HQ={}", self.high_qual_matches)];
        let mut low_qual_mismatches = 0;
        let mut snvs = 0;
        let mut indels = 0;
        let mut homopolymers = 0;
        let mut other_repeats = 0;
        for mismatch in mismatches {
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
        if low_qual_mismatches > 0 {
            sj.push(format!("LQ={low_qual_mismatches}"));
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
        sj.join(";")
    }

    pub fn reset_all(&mut self) {
        self.reset_index();
        self.reset_matches();
    }

    pub fn check_skip_hard_clips(&mut self) {
        let elements = self.elements.as_ref().unwrap();
        if self.element.map(|e| e.op) != Some(CigarOp::H) {
            return;
        }
        if self.move_forward {
            self.cigar_index += 1;
            if self.cigar_index >= self.element_count {
                self.exhausted = true;
                return;
            }
            self.element_index = 0;
            self.element = Some(elements[self.cigar_index as usize]);
        } else {
            self.cigar_index -= 1;
            if self.cigar_index < 0 {
                self.exhausted = true;
                return;
            }
            let e = elements[self.cigar_index as usize];
            self.element = Some(e);
            self.element_index = e.len - 1;
        }
    }
}

impl std::fmt::Display for ReadParseState {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let cigar_info = match (&self.elements, &self.element) {
            (Some(_), Some(e)) => format!(" refPos({}) element({}/{}{} {}/{}) ", self.ref_position, self.element_index + 1, e.len, e.op.as_char(), self.cigar_index + 1, self.element_count),
            _ => " ".to_string(),
        };
        write!(
            f,
            "{}: index({}/{}){}state({}) match({} hq={} pen={} mm={})",
            self.read.id(),
            self.read_index,
            self.base_length - 1,
            cigar_info,
            if self.mismatched { "mismatched" } else if self.exhausted { "exhausted" } else { "active" },
            self.base_matches,
            self.high_qual_matches,
            java_fixed(self.mismatch_penalty(), 1),
            self.mismatches().len()
        )
    }
}
