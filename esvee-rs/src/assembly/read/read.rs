//! Read: a BAM record with the cached and adjusted state the assembler keeps for it (Java assembly.read.Read).

use std::cell::Cell;
use std::sync::atomic::Ordering::Relaxed;
use std::sync::{Mutex, OnceLock};

use crate::assembly::common_utils::IndelCoords;
use crate::assembly::constants::*;
use crate::assembly::globals::READ_ID_TRIMMER;
use crate::assembly::mem_stats;
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::{A, H};
use crate::assembly::read::read_heap;
use crate::assembly::read::read_adjustments::find_low_base_qual_trim_count;
use crate::assembly::read::read_utils;
use crate::assembly::sam::*;
use crate::assembly::sequence_diff::SequenceDiffInfo;

pub type ReadP = H<Read>;

pub struct Read {
    record: SamRecord,
    orig_cigar_string: String,
    cigar_elements: Vec<CigarElement>,
    id: OnceLock<String>,
    alignment_start: i32,
    alignment_end: i32,
    unclipped_start: i32,
    unclipped_end: i32,
    mate_alignment_end: Cell<Option<i32>>,
    /// set when bases have been trimmed
    trimmed: Option<Box<TrimmedRead>>,

    mate_read: Option<ReadP>,
    /// what this read reports about its mate once the mate itself has been released
    mate_snapshot: Option<Box<MateSnapshot>>,
    has_junction_mate: bool,
    supplementary_data: OnceLock<Option<Box<SupplementaryReadData>>>,

    checked_indel_coords: Cell<bool>,
    indel_coords: Cell<Option<A<IndelCoords>>>,
    invalid_indel: bool,

    /// NO_INFERRED_POSITION when not set
    indel_inferred_unclipped_start: i32,
    indel_inferred_unclipped_end: i32,

    is_reference: bool,
    has_line_tail_start: bool,
    has_line_tail_end: bool,
    trim_count_start: i32,
    trim_count_end: i32,
    low_qual_trimmed: bool,
    mismatches: Option<A<Vec<SequenceDiffInfo>>>,
    /// estimate of the bytes held, for the memory diagnostics
    approx_bytes: u32,
}

impl Drop for Read {
    fn drop(&mut self) {
        mem_stats::READS_LIVE.fetch_sub(1, Relaxed);
        mem_stats::READ_BYTES_LIVE.fetch_sub(self.approx_bytes as i64, Relaxed);
    }
}

/// The properties of a linked mate that a read reports (mate chromosome, alignment, fragment end, mapped
/// state). When a junction group is finished and only this read still refers to its mate, the link is
/// replaced by this copy so the mate can be freed; nothing can change the mate after that point, so the
/// read reports the same values as it would through the link.
#[derive(Clone, Copy)]
struct MateSnapshot {
    chromosome: &'static str,
    alignment_start: i32,
    alignment_end: i32,
    fragment_end: i32,
    unmapped: bool,
}

/// the bases, qualities and cigar of a read after trimming (the record keeps the original ones)
struct TrimmedRead {
    bases: Vec<u8>,
    base_quals: Vec<u8>,
    cigar_string: String,
}

/// value of the indel-inferred unclipped positions when there are none (Java null)
const NO_INFERRED_POSITION: i32 = i32::MIN;

static TRIM_LOCK: Mutex<()> = Mutex::new(());

impl Read {
    pub fn new(record: SamRecord) -> ReadP {
        let orig_cigar_string = read_heap::string_from(&record.cigar_string());
        let cigar_elements = read_heap::vec_from(&cigar_elements_from_str(&orig_cigar_string));
        let start = record.alignment_start;
        let mut read = Read {
            record,
            orig_cigar_string,
            cigar_elements,
            id: OnceLock::new(),
            alignment_start: 0,
            alignment_end: 0,
            unclipped_start: 0,
            unclipped_end: 0,
            mate_alignment_end: Cell::new(None),
            trimmed: None,
            mate_read: None,
            mate_snapshot: None,
            has_junction_mate: false,
            supplementary_data: OnceLock::new(),
            checked_indel_coords: Cell::new(false),
            indel_coords: Cell::new(None),
            invalid_indel: false,
            indel_inferred_unclipped_start: NO_INFERRED_POSITION,
            indel_inferred_unclipped_end: NO_INFERRED_POSITION,
            is_reference: false,
            has_line_tail_start: false,
            has_line_tail_end: false,
            trim_count_start: 0,
            trim_count_end: 0,
            low_qual_trimmed: false,
            mismatches: None,
            approx_bytes: 0,
        };
        let r = &read.record;
        read.approx_bytes = (std::mem::size_of::<Read>()
            + 16
            + r.read_name.capacity()
            + r.cigar.capacity() * std::mem::size_of::<CigarElement>()
            + r.read_bases.capacity()
            + r.base_qualities.capacity()
            + r.aux.capacity()
            + read.orig_cigar_string.capacity()
            + read.cigar_elements.capacity() * std::mem::size_of::<CigarElement>()) as u32;
        mem_stats::READS_CREATED.fetch_add(1, Relaxed);
        mem_stats::READS_LIVE.fetch_add(1, Relaxed);
        mem_stats::READ_BYTES_LIVE.fetch_add(read.approx_bytes as i64, Relaxed);
        read.set_boundaries(start);
        H::new(read)
    }

    fn set_boundaries(&mut self, new_read_start: i32) {
        self.alignment_start = new_read_start;
        self.unclipped_start = self.alignment_start;
        if self.cigar_elements.is_empty() {
            self.alignment_end = self.alignment_start;
            self.unclipped_end = self.alignment_start;
            return;
        }
        let mut current_position = self.alignment_start;
        let n = self.cigar_elements.len();
        for i in 0..n {
            let element = self.cigar_elements[i];
            if i == 0 && element.op == CigarOp::S {
                self.unclipped_start -= element.len;
            }
            if element.op.consumes_reference_bases() {
                current_position += element.len;
            }
            if i == n - 1 {
                self.alignment_end = current_position - 1;
                self.unclipped_end = if element.op == CigarOp::S { self.alignment_end + element.len } else { self.alignment_end };
            }
        }
    }

    pub fn bam_record(&self) -> &SamRecord {
        &self.record
    }
    pub fn bam_record_mut(&mut self) -> &mut SamRecord {
        &mut self.record
    }

    pub fn set_mate_read(&mut self, mate: &ReadP) {
        self.mate_alignment_end.set(Some(mate.alignment_end()));
        self.mate_read = Some(mate.clone());
    }
    pub fn mate_read(&self) -> Option<&ReadP> {
        self.mate_read.as_ref()
    }
    /// drops the reference to the mate (used to release reads; not part of the Java class)
    pub fn clear_mate_read(&mut self) {
        self.mate_read = None;
    }
    /// Frees what a read no longer needs once its junction group has been assembled (not part of the Java
    /// class): the full read name, which is only used to form the id.
    pub fn release_unused_record_data(&mut self) {
        let _ = self.id();
        let name_bytes = self.record.read_name.capacity();
        self.record.read_name = String::new();
        let freed = (name_bytes as u32).min(self.approx_bytes);
        self.approx_bytes -= freed;
        mem_stats::READ_BYTES_LIVE.fetch_sub(freed as i64, Relaxed);
    }
    /// Replaces the link to the mate by a copy of what this read reports about it (not part of the Java
    /// class). Only for a mate that nothing else refers to, and after the junction group's assemblies are
    /// built: from then on the mate is only reached through the mate accessors of this read.
    pub fn replace_mate_read_with_snapshot(&mut self) {
        if let Some(m) = &self.mate_read {
            self.mate_snapshot = Some(Box::new(MateSnapshot {
                chromosome: m.chromosome(),
                alignment_start: m.alignment_start(),
                alignment_end: m.alignment_end(),
                fragment_end: if m.orientation().is_forward() { m.unclipped_start() } else { m.unclipped_end() },
                unmapped: m.is_unmapped(),
            }));
            self.mate_read = None;
        }
    }
    pub fn has_junction_mate(&self) -> bool {
        self.has_junction_mate
    }
    pub fn mark_junction_mate(&mut self) {
        self.has_junction_mate = true;
    }

    pub fn id(&self) -> &str {
        self.id.get_or_init(|| read_heap::string_from(&READ_ID_TRIMMER.trim(&self.record.read_name)))
    }
    pub fn chromosome(&self) -> &'static str {
        self.record.reference_name
    }
    pub fn cigar_elements(&self) -> &[CigarElement] {
        &self.cigar_elements
    }
    pub fn cigar_string(&self) -> &str {
        match &self.trimmed {
            Some(t) => &t.cigar_string,
            None => &self.orig_cigar_string,
        }
    }
    pub fn original_cigar_string(&self) -> &str {
        &self.orig_cigar_string
    }

    #[inline]
    pub fn alignment_start(&self) -> i32 {
        self.alignment_start
    }
    #[inline]
    pub fn alignment_end(&self) -> i32 {
        self.alignment_end
    }
    #[inline]
    pub fn unclipped_start(&self) -> i32 {
        self.unclipped_start
    }
    #[inline]
    pub fn unclipped_end(&self) -> i32 {
        self.unclipped_end
    }
    pub fn five_prime_fragment_position(&self) -> i32 {
        if self.positive_strand() { self.unclipped_start } else { self.unclipped_end }
    }

    pub fn is_left_clipped(&self) -> bool {
        self.unclipped_start != self.alignment_start || self.indel_inferred_unclipped_start != NO_INFERRED_POSITION
    }
    pub fn is_right_clipped(&self) -> bool {
        self.unclipped_end != self.alignment_end || self.indel_inferred_unclipped_end != NO_INFERRED_POSITION
    }
    pub fn left_clip_length(&self) -> i32 {
        (self.alignment_start - self.unclipped_start).max(0)
    }
    pub fn right_clip_length(&self) -> i32 {
        (self.unclipped_end - self.alignment_end).max(0)
    }

    #[inline]
    pub fn get_bases(&self) -> &[u8] {
        match &self.trimmed {
            Some(t) => &t.bases,
            None => &self.record.read_bases,
        }
    }
    #[inline]
    pub fn get_base_quality(&self) -> &[u8] {
        match &self.trimmed {
            Some(t) => &t.base_quals,
            None => &self.record.base_qualities,
        }
    }
    #[inline]
    pub fn bases_length(&self) -> i32 {
        self.get_bases().len() as i32
    }

    pub fn get_flags(&self) -> i32 {
        self.record.flags
    }
    pub fn is_unmapped(&self) -> bool {
        self.record.read_unmapped()
    }
    pub fn is_paired_read(&self) -> bool {
        self.record.read_paired()
    }
    pub fn positive_strand(&self) -> bool {
        !self.record.read_negative_strand()
    }
    pub fn negative_strand(&self) -> bool {
        self.record.read_negative_strand()
    }
    pub fn orientation(&self) -> Orientation {
        if self.record.read_negative_strand() { Orientation::Reverse } else { Orientation::Forward }
    }
    pub fn first_in_pair(&self) -> bool {
        self.record.read_paired() && self.record.first_of_pair()
    }
    pub fn mapping_quality(&self) -> i32 {
        self.record.mapping_quality
    }

    pub fn mate_chromosome(&self) -> Option<&'static str> {
        if let Some(m) = &self.mate_read {
            return Some(m.chromosome());
        }
        if let Some(m) = &self.mate_snapshot {
            return Some(m.chromosome);
        }
        if self.is_mate_mapped() { Some(self.record.mate_reference_name) } else { None }
    }
    pub fn mate_alignment_start(&self) -> i32 {
        if let Some(m) = &self.mate_read {
            return m.alignment_start();
        }
        if let Some(m) = &self.mate_snapshot {
            return m.alignment_start;
        }
        self.record.mate_alignment_start
    }
    pub fn mate_alignment_end(&self) -> i32 {
        if let Some(m) = &self.mate_read {
            return m.alignment_end();
        }
        if let Some(m) = &self.mate_snapshot {
            return m.alignment_end;
        }
        if let Some(e) = self.mate_alignment_end.get() {
            return e;
        }
        let mut end = get_mate_alignment_end(&self.record);
        if end == NO_POSITION {
            end = self.mate_alignment_start() + self.bases_length() - 1;
        }
        self.mate_alignment_end.set(Some(end));
        end
    }
    pub fn mate_fragment_end(&self) -> i32 {
        if let Some(m) = &self.mate_read {
            return if m.orientation().is_forward() { m.unclipped_start() } else { m.unclipped_end() };
        }
        if let Some(m) = &self.mate_snapshot {
            return m.fragment_end;
        }
        match self.record.get_string_attribute(MATE_CIGAR_ATTRIBUTE) {
            None => NO_POSITION,
            Some(mc) => get_five_prime_unclipped_position(self.mate_alignment_start(), &mc, !self.record.mate_negative_strand()),
        }
    }
    pub fn is_mate_mapped(&self) -> bool {
        if let Some(m) = &self.mate_read {
            return !m.is_unmapped();
        }
        if let Some(m) = &self.mate_snapshot {
            return !m.unmapped;
        }
        self.record.read_paired() && !self.record.mate_unmapped()
    }
    pub fn is_mate_unmapped(&self) -> bool {
        !self.is_mate_mapped()
    }
    pub fn mate_orientation(&self) -> Orientation {
        if !self.record.read_paired() {
            return Orientation::Forward;
        }
        if self.record.mate_negative_strand() { Orientation::Reverse } else { Orientation::Forward }
    }

    pub fn has_supplementary(&self) -> bool {
        self.supplementary_data().is_some()
    }
    pub fn is_supplementary(&self) -> bool {
        self.record.supplementary()
    }

    /// Read.makeReadLinks
    pub fn make_read_links(this: &ReadP, other: &ReadP) {
        if this.record.supplementary() == other.record.supplementary() && this.first_in_pair() != other.first_in_pair() {
            this.m().mate_read = Some(other.clone());
            other.m().set_mate_read(this);
            if !this.record.read_paired() {
                this.m().record.set_flag(FLAG_READ_PAIRED, true);
            }
            if !other.record.read_paired() {
                other.m().record.set_flag(FLAG_READ_PAIRED, true);
            }
            if this.is_unmapped() && !other.record.mate_unmapped() {
                other.m().record.set_flag(FLAG_MATE_UNMAPPED, true);
            } else if other.is_unmapped() && !this.record.mate_unmapped() {
                this.m().record.set_flag(FLAG_MATE_UNMAPPED, true);
            }
        }
    }

    pub fn supplementary_data(&self) -> Option<&SupplementaryReadData> {
        self.supplementary_data
            .get_or_init(|| SupplementaryReadData::extract_alignment(self.record.get_string_attribute(SUPPLEMENTARY_ATTRIBUTE).as_deref()).map(Box::new))
            .as_deref()
    }

    pub fn get_read_index_at_reference_position(&self, ref_position: i32) -> i32 {
        self.get_read_index_at_reference_position_ext(ref_position, false)
    }
    pub fn get_read_index_at_reference_position_ext(&self, ref_position: i32, allow_extrapolation: bool) -> i32 {
        read_utils::get_read_index_at_reference_position(self, ref_position, allow_extrapolation)
    }

    pub fn indel_coords(&self) -> Option<A<IndelCoords>> {
        if !self.checked_indel_coords.get() && !self.invalid_indel {
            self.checked_indel_coords.set(true);
            self.indel_coords.set(IndelCoords::find_indel_coords(self.alignment_start, &self.cigar_elements, MIN_INDEL_SUPPORT_LENGTH).map(A::new));
        }
        let v = self.indel_coords.take();
        self.indel_coords.set(v.clone());
        v
    }
    pub fn mark_invalid_indel(&mut self) {
        self.invalid_indel = true;
    }
    pub fn invalid_indel(&self) -> bool {
        self.invalid_indel
    }

    pub fn matches_fragment(&self, other: &Read, allow_read_match: bool) -> bool {
        if self.id() != other.id() {
            return false;
        }
        allow_read_match || self.get_flags() != other.get_flags()
    }

    pub fn find_matching_fragment_support(support: &[ReadP], read: &Read) -> Vec<ReadP> {
        support.iter().filter(|x| x.matches_fragment(read, false)).cloned().collect()
    }

    pub fn sample_index(&self) -> i32 {
        self.record.sample_index
    }
    pub fn is_reference(&self) -> bool {
        self.is_reference
    }
    pub fn mark_reference(&mut self) {
        self.is_reference = true;
    }
    pub fn has_line_tail_start(&self) -> bool {
        self.has_line_tail_start
    }
    pub fn has_line_tail_end(&self) -> bool {
        self.has_line_tail_end
    }
    pub fn has_line_tail(&self, is_forward_junction: bool) -> bool {
        if is_forward_junction { self.has_line_tail_end } else { self.has_line_tail_start }
    }
    pub fn mark_line_tail(&mut self, is_start: bool) {
        if is_start {
            self.has_line_tail_start = true;
        } else {
            self.has_line_tail_end = true;
        }
    }

    pub fn get_bases_string(&self) -> String {
        String::from_utf8_lossy(self.get_bases()).into_owned()
    }

    pub fn trim_bases(&mut self, count: i32, from_start: bool) {
        let mut remaining_bases = count;
        let new_base_length = (self.bases_length() - count).max(1) as usize;
        let mut new_bases = vec![0u8; new_base_length];
        let mut new_base_quals = vec![0u8; new_base_length];
        let mut new_read_start = self.alignment_start;
        if from_start {
            while !self.cigar_elements.is_empty() {
                let element = self.cigar_elements[0];
                if remaining_bases == 0 {
                    if element.op.consumes_reference_bases() && !element.op.consumes_read_bases() {
                        self.cigar_elements.remove(0);
                        new_read_start += element.len;
                        continue;
                    }
                    break;
                }
                if !element.op.consumes_read_bases() {
                    self.cigar_elements.remove(0);
                    if element.op.consumes_reference_bases() {
                        new_read_start += element.len;
                    }
                    continue;
                }
                if element.len <= remaining_bases {
                    self.cigar_elements.remove(0);
                    remaining_bases -= element.len;
                    if element.op.consumes_reference_bases() {
                        new_read_start += element.len;
                    }
                } else {
                    self.cigar_elements[0] = CigarElement::new(element.len - remaining_bases, element.op);
                    if element.op.consumes_reference_bases() {
                        new_read_start += remaining_bases;
                    }
                    remaining_bases = 0;
                }
            }
            copy_array(self.get_bases(), &mut new_bases, count as usize, 0);
            copy_array(self.get_base_quality(), &mut new_base_quals, count as usize, 0);
        } else {
            while !self.cigar_elements.is_empty() {
                let last_index = self.cigar_elements.len() - 1;
                let element = self.cigar_elements[last_index];
                if remaining_bases == 0 {
                    if element.op.consumes_reference_bases() && !element.op.consumes_read_bases() {
                        self.cigar_elements.remove(last_index);
                        continue;
                    }
                    break;
                }
                if !element.op.consumes_read_bases() {
                    self.cigar_elements.remove(last_index);
                    continue;
                }
                if element.len <= remaining_bases {
                    self.cigar_elements.remove(last_index);
                    remaining_bases -= element.len;
                } else {
                    self.cigar_elements[last_index] = CigarElement::new(element.len - remaining_bases, element.op);
                    remaining_bases = 0;
                }
            }
            copy_array(self.get_bases(), &mut new_bases, 0, 0);
            copy_array(self.get_base_quality(), &mut new_base_quals, 0, 0);
        }
        self.trimmed = Some(Box::new(TrimmedRead { bases: new_bases, base_quals: new_base_quals, cigar_string: cigar_elements_to_str(&self.cigar_elements) }));
        if from_start {
            self.trim_count_start += count;
        } else {
            self.trim_count_end += count;
        }
        self.set_boundaries(new_read_start);
    }

    pub fn trim_count(&self) -> i32 {
        self.trim_count_start + self.trim_count_end
    }
    pub fn trim_count_start(&self) -> i32 {
        self.trim_count_start
    }
    pub fn trim_count_end(&self) -> i32 {
        self.trim_count_end
    }
    pub fn mark_low_qual_trimmed(&mut self) {
        self.low_qual_trimmed = true;
    }
    pub fn low_qual_trimmed(&self) -> bool {
        self.low_qual_trimmed
    }

    pub fn set_indel_inferred_unclipped_positions(&mut self, inferred_start: Option<i32>, inferred_end: Option<i32>) {
        self.indel_inferred_unclipped_start = inferred_start.unwrap_or(NO_INFERRED_POSITION);
        self.indel_inferred_unclipped_end = inferred_end.unwrap_or(NO_INFERRED_POSITION);
    }
    pub fn has_indel_implied_unclipped_start(&self) -> bool {
        self.indel_inferred_unclipped_start != NO_INFERRED_POSITION
    }
    pub fn has_indel_implied_unclipped_end(&self) -> bool {
        self.indel_inferred_unclipped_end != NO_INFERRED_POSITION
    }
    pub fn indel_implied_unclipped_start(&self) -> i32 {
        if self.indel_inferred_unclipped_start == NO_INFERRED_POSITION { 0 } else { self.indel_inferred_unclipped_start }
    }
    pub fn indel_implied_unclipped_end(&self) -> i32 {
        if self.indel_inferred_unclipped_end == NO_INFERRED_POSITION { 0 } else { self.indel_inferred_unclipped_end }
    }
    pub fn min_unclipped_start(&self) -> i32 {
        if self.indel_inferred_unclipped_start == NO_INFERRED_POSITION { self.unclipped_start } else { self.unclipped_start.min(self.indel_inferred_unclipped_start) }
    }
    pub fn max_unclipped_end(&self) -> i32 {
        if self.indel_inferred_unclipped_end == NO_INFERRED_POSITION { self.unclipped_end } else { self.unclipped_end.max(self.indel_inferred_unclipped_end) }
    }

    /// the list is shared with the ReadParseState it came from, as in Java
    pub fn set_extension_mismatches(&mut self, seq_mismatches: A<Vec<SequenceDiffInfo>>) {
        self.mismatches = Some(seq_mismatches);
    }
    pub fn mismatches(&self) -> &[SequenceDiffInfo] {
        match &self.mismatches {
            Some(m) => m,
            None => &[],
        }
    }

    /// Read.trimLowQualBases (synchronized in Java: reads can be reached from several phase group tasks)
    pub fn trim_low_qual_bases(this: &ReadP) {
        let _guard = TRIM_LOCK.lock().unwrap();
        if this.low_qual_trimmed() {
            return;
        }
        let read_index_start = 0;
        let read_index_end = this.bases_length() - 1;
        let trim_count = find_low_base_qual_trim_count(this, read_index_start, read_index_end);
        if trim_count > 0 {
            let from_start = this.negative_strand();
            this.m().trim_bases(trim_count, from_start);
            this.m().mark_low_qual_trimmed();
        }
    }
}

impl std::fmt::Display for Read {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "id({}) coords({}:{}-{}) cigar({}) mate({}:{}) flags({})",
            self.id(),
            self.chromosome(),
            self.alignment_start,
            self.alignment_end,
            self.cigar_string(),
            self.mate_chromosome().unwrap_or("null"),
            self.mate_alignment_start(),
            self.record.flags
        )
    }
}

/// hmf-common Arrays.copyArray(source, dest, sourceIndexStart, destIndexStart)
pub fn copy_array(source: &[u8], dest: &mut [u8], source_index_start: usize, dest_index_start: usize) {
    let mut d = dest_index_start;
    let mut s = source_index_start;
    while s < source.len() && d < dest.len() {
        dest[d] = source[s];
        s += 1;
        d += 1;
    }
}

/// hmf-common Arrays.subsetArray (inclusive indices)
pub fn subset_array(source: &[u8], start_index: usize, end_index: usize) -> Vec<u8> {
    source[start_index..=end_index].to_vec()
}
