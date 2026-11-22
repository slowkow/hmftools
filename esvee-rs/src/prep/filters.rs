//! Read filters (ReadFilters, ReadFilterConfig), with the helpers they use from other packages: LINE tails
//! (LineUtils), repeats (RepeatInfo) and discordant fragments (CommonUtils).

use super::rec::*;
use super::types::*;
use rust_htslib::bam;

// PrepConstants and SvConstants
pub const MIN_ALIGNMENT_BASES: i32 = 50;
pub const MIN_CALC_ALIGNMENT_SCORE: i32 = 40;
pub const MIN_SUPP_CALC_ALIGNMENT_SCORE: i32 = 10;
pub const MIN_CALC_ALIGNMENT_LOWER_SCORE: i32 = 35;
pub const MIN_ALIGNMENT_SCORE_DIFF: i32 = 15;
pub const MIN_INSERT_ALIGNMENT_OVERLAP: i32 = 5;
pub const MIN_VARIANT_LENGTH: i32 = 32;
pub const MIN_VARIANT_LENGTH_LOWER: i32 = 24;
pub const MIN_SOFT_CLIP_LENGTH: i32 = MIN_VARIANT_LENGTH_LOWER;
pub const LINE_POLY_AT_TEST_LEN: i32 = 18;
pub const LINE_POLY_AT_REQ: i32 = 16;
pub const MIN_LINE_SOFT_CLIP_LENGTH: i32 = LINE_POLY_AT_TEST_LEN;
pub const MIN_SOFT_CLIP_HIGH_QUAL_PERC: f64 = 0.75;
pub const MAX_SOFT_CLIP_LOW_QUAL_COUNT: i32 = 5;
pub const MIN_INDEL_ANCHOR_LENGTH: i32 = 10;
pub const SUPP_SHORT_ALIGNMENT_MIN_MAP_QUAL: i32 = 10;
pub const REPEAT_BREAK_CHECK_LENGTH: i32 = 9;
pub const REPEAT_BREAK_MIN_MAP_QUAL: i32 = 40;
pub const REPEAT_BREAK_MIN_SC_LENGTH: i32 = 50;
pub const MAX_JUNC_POSITION_DIFF: i32 = 15;
pub const MAX_JUNC_SC_VS_INDEL_POSITION_DIFF: i32 = 50;
pub const UNPAIRED_READ_JUNCTION_DISTANCE: i32 = 5;
pub const MAX_SUPPORT_FRAGMENT_DISTANCE: i32 = 1000;
pub const MAX_HIGH_QUAL_BASE_MISMATCHES: i32 = 1;
pub const MIN_EXACT_BASE_PERC: f64 = 0.75;
pub const MIN_MAP_QUALITY: i32 = 20;
pub const MIN_INDEL_LENGTH: i32 = MIN_VARIANT_LENGTH;
pub const MIN_JUNCTION_SUPPORT: i32 = 2;
pub const MIN_HOTSPOT_JUNCTION_SUPPORT: i32 = 1;
pub const DEFAULT_MAX_FRAGMENT_LENGTH: i32 = 1100;
pub const MIN_UPPER_FRAGMENT_LENGTH: i32 = 800;
pub const DEFAULT_READ_LENGTH: i32 = 151;
pub const LINE_REF_BASE_REPEAT_FACTOR: f64 = 1.5;
pub const LINE_REF_BASE_MAX_LENGTH: i32 = 20;
pub const POLY_G_LENGTH: i32 = 16;
pub const LOW_BASE_QUAL_THRESHOLD: u8 = 26;
pub const REPEAT_MIN_COUNT: usize = 4;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SequencingType {
    Illumina,
    Ultima,
    Sbx,
}

impl SequencingType {
    pub fn parse(s: &str) -> Option<SequencingType> {
        match s {
            "ILLUMINA" => Some(SequencingType::Illumina),
            "ULTIMA" => Some(SequencingType::Ultima),
            "SBX" => Some(SequencingType::Sbx),
            _ => None,
        }
    }
    pub fn is_illumina(self) -> bool {
        self == SequencingType::Illumina
    }
    pub fn is_sbx(self) -> bool {
        self == SequencingType::Sbx
    }
    /// SvConstants.hasPairedReads
    pub fn has_paired_reads(self) -> bool {
        self.is_illumina()
    }
}

/// SvConstants.maxConcordantFragmentLength
pub fn max_concordant_fragment_length(observed_max: i32) -> i32 {
    observed_max.max(MIN_UPPER_FRAGMENT_LENGTH)
}

#[inline]
pub fn below_min_qual(q: u8) -> bool {
    // Java bytes are signed: a quality above 127 is negative there and so low. One comparison covers both
    // ranges: 0..=26 wraps to 229..=255 and 128..=255 becomes 101..=228.
    q.wrapping_sub(LOW_BASE_QUAL_THRESHOLD + 1) > 100
}

#[inline]
pub fn above_min_qual(q: u8) -> bool {
    !below_min_qual(q)
}

/// ReadFilterConfig
#[derive(Clone, Debug)]
pub struct ReadFilterConfig {
    pub min_alignment_bases: i32,
    pub min_map_quality: i32,
    pub min_insert_alignment_overlap: i32,
    pub min_indel_length: i32,
    pub min_soft_clip_length: i32,
    pub min_soft_clip_high_qual_perc: f64,
    pub min_junction_support: i32,
    pub observed_frag_length_min: i32,
    pub observed_frag_length_max: i32,
    pub max_supporting_fragment_distance: i32,
    pub sequencing: SequencingType,
}

impl ReadFilterConfig {
    pub fn new(min_alignment_bases: i32, min_junction_support: i32) -> Self {
        ReadFilterConfig {
            min_alignment_bases,
            min_map_quality: MIN_MAP_QUALITY,
            min_insert_alignment_overlap: MIN_INSERT_ALIGNMENT_OVERLAP,
            min_indel_length: MIN_INDEL_LENGTH,
            min_soft_clip_length: MIN_SOFT_CLIP_LENGTH,
            min_soft_clip_high_qual_perc: MIN_SOFT_CLIP_HIGH_QUAL_PERC,
            min_junction_support,
            observed_frag_length_min: DEFAULT_READ_LENGTH,
            observed_frag_length_max: DEFAULT_MAX_FRAGMENT_LENGTH,
            max_supporting_fragment_distance: DEFAULT_MAX_FRAGMENT_LENGTH,
            sequencing: SequencingType::Illumina,
        }
    }
    pub fn set_fragment_lengths(&mut self, min_length: i32, max_length: i32) {
        self.observed_frag_length_min = min_length;
        self.observed_frag_length_max = max_length;
        self.max_supporting_fragment_distance = max_length.min(MAX_SUPPORT_FRAGMENT_DISTANCE);
    }
}

impl Default for ReadFilterConfig {
    fn default() -> Self {
        Self::new(MIN_ALIGNMENT_BASES, MIN_JUNCTION_SUPPORT)
    }
}

/// CommonUtils.isDiscordantFragment without supplementary data
pub fn is_discordant_fragment(rec: &bam::Record, max_concordant_fragment_length: i32) -> bool {
    let flags = rec.flags();
    if flags & FLAG_UNMAPPED != 0 || flags & FLAG_PAIRED == 0 || flags & FLAG_MATE_UNMAPPED != 0 {
        return false;
    }
    if rec.tid() != rec.mtid() {
        return true;
    }
    if (flags & FLAG_REVERSE != 0) == (flags & FLAG_MATE_REVERSE != 0) {
        return true;
    }
    let fragment_size = (rec.insert_size() as i32).wrapping_abs();
    if fragment_size == 0 || (max_concordant_fragment_length > 0 && fragment_size >= max_concordant_fragment_length) {
        return true;
    }
    // a duplication orientation of fragments that do not overlap
    if fragment_size > 2 * rec.seq_len() as i32 {
        if flags & FLAG_REVERSE == 0 {
            if rec.mpos() < rec.pos() {
                return true;
            }
        } else if rec.pos() < rec.mpos() {
            return true;
        }
    }
    false
}

/// ReadFilters.isChimericRead: discordant, or with an unmapped mate
pub fn is_chimeric_read(rec: &bam::Record, config: &ReadFilterConfig) -> bool {
    if mate_unmapped(rec) {
        return true;
    }
    is_discordant_fragment(rec, max_concordant_fragment_length(config.observed_frag_length_max))
}

/// ReadFilters.filterLowQualRead: more than half of the bases are low quality. A base that is not A, C, G, T
/// or N becomes N in the record; in Java only for as far as the check reads (finding A17).
pub fn filter_low_qual_read(rec: &mut bam::Record, config: &ReadFilterConfig) -> bool {
    filter_low_qual_read_perc(rec, config, 0.5)
}

pub fn filter_low_qual_read_perc(rec: &mut bam::Record, config: &ReadFilterConfig, max_permitted_low_percent: f64) -> bool {
    if !config.sequencing.is_illumina() {
        return false;
    }
    let base_length = rec.seq_len();
    let threshold = (base_length as f64 * max_permitted_low_percent).floor() as usize + 1;
    let mut non_standard = has_non_standard_bases(rec);
    if non_standard && !crate::common::compat::java("A17") {
        // Java: rewrites inside the counting loop, so the bases after its early exit keep their letters
        replace_non_standard_bases(rec, base_length);
        non_standard = false;
    }
    if !non_standard && max_permitted_low_percent == 0.5 {
        // With the threshold above half the read only one of the two counts can reach it, so the answer is
        // whether the low-quality bases do, and without bases to rewrite the order of the bases is of no
        // account. This is the path of nearly every read.
        let low_qual_count: usize = rec.qual().chunks(128).map(|c| c.iter().map(|&q| below_min_qual(q) as u8).sum::<u8>() as usize).sum();
        return low_qual_count >= threshold;
    }
    let mut low_qual_count = 0usize;
    let mut result = false;
    let mut checked = base_length;
    {
        let quals = rec.qual();
        for (i, &q) in quals.iter().enumerate().take(base_length) {
            if below_min_qual(q) {
                low_qual_count += 1;
                if low_qual_count >= threshold {
                    result = true;
                    checked = i + 1;
                    break;
                }
            } else if i + 1 - low_qual_count >= threshold {
                checked = i + 1;
                break;
            }
        }
    }
    if non_standard {
        replace_non_standard_bases(rec, checked);
    }
    result
}

/// Whether both bases of a packed sequence byte are one of A, C, G, T, N.
static STANDARD_PAIR: [bool; 256] = {
    let mut t = [false; 256];
    let mut i = 0;
    while i < 256 {
        let (hi, lo) = (i >> 4, i & 0xf);
        t[i] = matches!(hi, 1 | 2 | 4 | 8 | 15) && matches!(lo, 1 | 2 | 4 | 8 | 15);
        i += 1;
    }
    t
};

fn packed_seq(rec: &bam::Record) -> &[u8] {
    let b = &rec.inner;
    let off = b.core.l_qname as usize + b.core.n_cigar as usize * 4;
    unsafe { std::slice::from_raw_parts(b.data.add(off), (b.core.l_qseq as usize).div_ceil(2)) }
}

fn has_non_standard_bases(rec: &bam::Record) -> bool {
    let seq = packed_seq(rec);
    let n = rec.seq_len();
    if n == 0 {
        return false;
    }
    let full = n / 2;
    // without an early exit: nearly every read is all standard bases, and the loop then has no branch to miss
    let mut standard = true;
    for &b in &seq[..full] {
        standard &= STANDARD_PAIR[b as usize];
    }
    // an odd length leaves a last base in the high half of the last byte
    if n % 2 == 1 {
        standard &= STANDARD_PAIR[(seq[full] | 0x0f) as usize];
    }
    !standard
}

/// Sets the first `count` bases that are none of A, C, G, T, N to N.
fn replace_non_standard_bases(rec: &mut bam::Record, count: usize) {
    let b = &rec.inner;
    let off = b.core.l_qname as usize + b.core.n_cigar as usize * 4;
    let seq = unsafe { std::slice::from_raw_parts_mut(b.data.add(off), (b.core.l_qseq as usize).div_ceil(2)) };
    let standard = |code: u8| matches!(code, 1 | 2 | 4 | 8 | 15);
    let full = count / 2;
    for byte in seq.iter_mut().take(full) {
        if !standard(*byte >> 4) {
            *byte |= 0xf0;
        }
        if !standard(*byte & 0xf) {
            *byte |= 0x0f;
        }
    }
    if count % 2 == 1 && !standard(seq[full] >> 4) {
        seq[full] |= 0xf0;
    }
}

/// ReadFilters.filterPoorMappedSupplementary
pub fn filter_poor_mapped_supplementary(rec: &bam::Record) -> bool {
    if !has(rec, FLAG_SUPPLEMENTARY) {
        return false;
    }
    let aligned_length: i32 = rec.raw_cigar().iter().filter(|c| op(**c) == OP_M).map(|c| oplen(*c)).sum();
    let has_low_map_qual = (rec.mapq() as i32) < SUPP_SHORT_ALIGNMENT_MIN_MAP_QUAL;
    if aligned_length < MIN_CALC_ALIGNMENT_SCORE && has_low_map_qual {
        return true;
    }
    // Java: asks for the score with a minimum length of 40; under that the answer is -1 ("not calculated"),
    // which the tests below take for a score, so every supplementary with under 40 aligned bases is dropped
    // whatever its mapping quality. Corrected: the score is calculated for any length.
    let min_length = if crate::common::compat::java("P10") { MIN_CALC_ALIGNMENT_SCORE } else { 1 };
    if aligned_length < min_length {
        return true;
    }
    let bases = seq_bytes(rec);
    let adjusted = calc_repeat_trimmed_alignment_score(rec, &bases, aligned_length, min_length, true);
    if adjusted < MIN_CALC_ALIGNMENT_SCORE as f64 && has_low_map_qual {
        return true;
    }
    adjusted < MIN_SUPP_CALC_ALIGNMENT_SCORE as f64
}

/// LineUtils.hasLineTail(record): a poly-A soft clip on the left or a poly-T soft clip on the right
pub fn has_line_tail(rec: &bam::Record) -> bool {
    let cigar = rec.raw_cigar();
    if cigar.len() <= 1 {
        return false;
    }
    let mut bases: Option<Vec<u8>> = None;
    for from_start in [true, false] {
        let element = if from_start { cigar[0] } else { cigar[cigar.len() - 1] };
        if op(element) != OP_S {
            continue;
        }
        let sc_base_count = oplen(element);
        if sc_base_count == 0 {
            continue;
        }
        if sc_base_count < LINE_POLY_AT_REQ {
            continue; // too short to test, as the inner routine finds
        }
        let line_base = if from_start { b'A' } else { b'T' };
        let bases = bases.get_or_insert_with(|| seq_bytes(rec));
        let soft_clip_index = if from_start { sc_base_count - 1 } else { bases.len() as i32 - sc_base_count };
        if has_line_tail_bases(bases, soft_clip_index, from_start, line_base, Some(rec.qual())) {
            return true;
        }
    }
    false
}

/// LineUtils.hasLineTail over bases
pub fn has_line_tail_bases(bases: &[u8], soft_clip_index: i32, left_clipped: bool, line_base: u8, base_quals: Option<&[u8]>) -> bool {
    let n = bases.len() as i32;
    let sc_length = if left_clipped { soft_clip_index + 1 } else { n - soft_clip_index };
    if sc_length < LINE_POLY_AT_REQ {
        return false;
    }
    let mut index = soft_clip_index;
    let (mut total_base_count, mut non_line_base_count, mut line_base_count) = (0, 0, 0);
    while index >= 0 && index < n {
        total_base_count += 1;
        if bases[index as usize] == line_base {
            line_base_count += 1;
        } else if base_quals.is_none_or(|q| !is_uncertain_base_qual(q[index as usize])) {
            if total_base_count > LINE_POLY_AT_TEST_LEN {
                break;
            }
            non_line_base_count += 1;
            if non_line_base_count > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ {
                return false;
            }
        }
        index += if left_clipped { -1 } else { 1 };
    }
    if line_base_count < LINE_POLY_AT_REQ || non_line_base_count > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ {
        return false;
    }
    // the reference bases next to the soft clip must not hold a similar run of the same base
    let ref_index_start = if left_clipped { soft_clip_index + 1 } else { soft_clip_index - 1 };
    let ref_line_base_count = find_ref_base_repeat_count(bases, ref_index_start, left_clipped, line_base);
    if ref_line_base_count > LINE_REF_BASE_MAX_LENGTH {
        return false;
    }
    line_base_count as f64 >= ref_line_base_count as f64 * LINE_REF_BASE_REPEAT_FACTOR
}

/// BaseQualAdjustment.isUncertainBaseQual (a signed byte comparison in Java)
#[inline]
fn is_uncertain_base_qual(q: u8) -> bool {
    q <= 1 || q > 127
}

/// LineUtils.findRefBaseRepeatCount
pub fn find_ref_base_repeat_count(bases: &[u8], mut index: i32, move_forward: bool, test_base: u8) -> i32 {
    let (mut match_count, mut initial_match_count, mut mismatch_count, mut base_count) = (0, 0, 0, 0);
    while index >= 0 && index < bases.len() as i32 {
        base_count += 1;
        if bases[index as usize] == test_base {
            match_count += 1;
            if mismatch_count == 0 {
                initial_match_count += 1;
            }
        } else {
            mismatch_count += 1;
        }
        if base_count >= 10 {
            if mismatch_count as f64 / base_count as f64 > 0.1 {
                return initial_match_count;
            }
        } else if mismatch_count >= 2 {
            return initial_match_count;
        }
        index += if move_forward { 1 } else { -1 };
    }
    if mismatch_count == 0 { initial_match_count } else { match_count }
}

/// ReadFilters.ignoreRead: low-quality reads, poorly mapped supplementaries and ordinary reads of a proper
/// local pair take no further part.
pub fn ignore_read(rec: &mut bam::Record, has_line_tail: bool, config: &ReadFilterConfig) -> bool {
    if !has_line_tail && filter_low_qual_read(rec, config) {
        return true;
    }
    if filter_poor_mapped_supplementary(rec) {
        return true;
    }
    let flags = rec.flags();
    if flags & FLAG_PAIRED == 0 {
        return false;
    }
    let cigar = rec.raw_cigar();
    if flags & FLAG_SUPPLEMENTARY != 0 || cigar.len() > 1 {
        return false;
    }
    if !cigar.iter().any(|c| op(*c) == OP_M) {
        return false;
    }
    if flags & (FLAG_UNMAPPED | FLAG_MATE_UNMAPPED) != 0 {
        return false;
    }
    if is_chimeric_read(rec, config) {
        return false;
    }
    match aux_str(rec, b"MC") {
        Some(mc) => mate_cigar_is_aligned(mc, rec.seq_len()),
        None => false,
    }
}

/// Whether a cigar string is "<length>M" as String.format("%dM") writes it.
pub fn mate_cigar_is_aligned(cigar: &[u8], length: usize) -> bool {
    let n = cigar.len();
    if n < 2 || cigar[n - 1] != b'M' || cigar[0] == b'0' && n > 2 {
        return false;
    }
    let mut v: usize = 0;
    for &c in &cigar[..n - 1] {
        if !c.is_ascii_digit() {
            return false;
        }
        v = v.wrapping_mul(10).wrapping_add((c - b'0') as usize);
    }
    v == length
}

/// ReadFilters.isCandidateSupportingRead
pub fn is_candidate_supporting_read(read: &PrepRead, config: &ReadFilterConfig) -> bool {
    if is_chimeric_read(&read.rec, config) {
        return true;
    }
    read.is_left_clipped() || read.is_right_clipped() || read.max_indel_length >= MIN_INDEL_SUPPORT_LENGTH
}

fn contains_run(bases: &[u8], base: u8, length: usize) -> bool {
    let mut run = 0;
    for &b in bases {
        if b == base {
            run += 1;
            if run >= length {
                return true;
            }
        } else {
            run = 0;
        }
    }
    false
}

/// ReadFilters.checkFilters: every filter that stops a read from establishing a junction
pub fn check_filters(read: &mut PrepRead, config: &ReadFilterConfig) {
    let aligned_bases = read.aligned_base_length;
    if aligned_bases < config.min_alignment_bases {
        read.add_filter(F_MIN_ALIGN_MATCH);
    }
    if read.map_quality() < config.min_map_quality {
        read.add_filter(F_MIN_MAP_QUAL);
    }
    let flags = read.rec.flags();
    if flags & FLAG_PAIRED != 0 && !mate_unmapped(&read.rec) && flags & FLAG_UNMAPPED == 0 {
        let insert_alignment_overlap = (inferred_insert_size_abs(&read.rec) - aligned_bases).wrapping_abs();
        if insert_alignment_overlap < config.min_insert_alignment_overlap {
            read.add_filter(F_INSERT_MAP_OVERLAP);
        }
    }
    let (sc_left, sc_right) = (read.sc_left, read.sc_right);
    if !read.has_line_tail && sc_left < config.min_soft_clip_length && sc_right < config.min_soft_clip_length && read.max_indel_length < config.min_indel_length {
        read.add_filter(F_SOFT_CLIP_LENGTH);
    }
    if sc_left > 0 || sc_right > 0 {
        let quals = read.rec.qual();
        let n = quals.len() as i32;
        let use_left_clip = sc_left > sc_right;
        let sc_range_start = if use_left_clip { 0 } else { n - sc_right };
        let sc_range_end = if use_left_clip { sc_left } else { n };
        let sc_length = if use_left_clip { sc_left } else { sc_right };
        let above_qual = quals[sc_range_start as usize..sc_range_end as usize].iter().filter(|q| above_min_qual(**q)).count() as i32;
        let mut add = 0;
        if !read.has_line_tail && (above_qual as f64 / sc_length as f64) < config.min_soft_clip_high_qual_perc {
            add |= F_SOFT_CLIP_BASE_QUAL;
            if sc_length - above_qual >= MAX_SOFT_CLIP_LOW_QUAL_COUNT {
                add |= F_SOFT_CLIP_LOW_BASE_QUAL;
            }
        }
        // poly-G or poly-C inserts, then a break in a repetitive section
        let sc_bases = &read.bases[sc_range_start as usize..sc_range_end as usize];
        let poly_g = sc_length >= POLY_G_LENGTH && (contains_run(sc_bases, b'G', POLY_G_LENGTH as usize) || contains_run(sc_bases, b'C', POLY_G_LENGTH as usize));
        if poly_g {
            add |= F_POLY_G_SC;
        } else if !read.has_filter(F_SOFT_CLIP_LENGTH)
            && (read.map_quality() < REPEAT_BREAK_MIN_MAP_QUAL || sc_length < REPEAT_BREAK_MIN_SC_LENGTH)
            && is_repetitive_section_break(&read.bases, use_left_clip, sc_length)
        {
            add |= F_BREAK_IN_REPEAT;
        }
        read.add_filter(add);
    }
}

/// ReadFilters.isRepetitiveSectionBreak: the 18 bases around the soft clip start are one repeat of 1 to 3 bases
pub fn is_repetitive_section_break(read_bases: &[u8], left_clipped: bool, sc_length: i32) -> bool {
    let read_length = read_bases.len() as i32;
    let start_index = if left_clipped { (sc_length - REPEAT_BREAK_CHECK_LENGTH).max(0) } else { (read_length - sc_length - REPEAT_BREAK_CHECK_LENGTH).max(0) };
    let end_index = (start_index + 2 * REPEAT_BREAK_CHECK_LENGTH).min(read_length);
    if end_index - start_index < 2 * REPEAT_BREAK_CHECK_LENGTH {
        return false;
    }
    let (s, e) = (start_index as usize, end_index as usize);
    let b = read_bases;
    let first = b[s];
    if (s + 1..e).all(|i| b[i] == first) {
        return true;
    }
    let second = b[s + 1];
    let mut all_match = true;
    let mut i = s + 2;
    while i < e - 1 {
        if b[i] != first || b[i + 1] != second {
            all_match = false;
            break;
        }
        i += 2;
    }
    if all_match {
        return true;
    }
    let third = b[s + 2];
    let mut i = s + 3;
    while i < e - 2 {
        if b[i] != first || b[i + 1] != second || b[i + 2] != third {
            return false;
        }
        i += 3;
    }
    true
}

/// A repeat found by RepeatInfo.findRepeats: unit length and count
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatInfo {
    pub index: usize,
    pub bases: Vec<u8>,
    pub count: usize,
}

impl RepeatInfo {
    pub fn total_length(&self) -> usize {
        self.count * self.bases.len()
    }
}

fn find_single_base_repeat(bases: &[u8], index: usize, min_repeat_count: usize) -> Option<RepeatInfo> {
    if index + min_repeat_count - 1 >= bases.len() {
        return None;
    }
    let mut i = index + 1;
    while i < bases.len() && bases[i] == bases[index] {
        i += 1;
    }
    let count = i - index;
    (count >= min_repeat_count).then(|| RepeatInfo { index, bases: vec![bases[index]], count })
}

fn find_dual_base_repeat(bases: &[u8], index: usize, min_repeat_count: usize) -> Option<RepeatInfo> {
    if index + min_repeat_count * 2 - 1 >= bases.len() {
        return None;
    }
    let mut i = index + 2;
    while i + 1 < bases.len() {
        if bases[i] != bases[index] || bases[i + 1] != bases[index + 1] {
            break;
        }
        i += 2;
    }
    let count = (i - index) / 2;
    (count >= min_repeat_count).then(|| RepeatInfo { index, bases: bases[index..index + 2].to_vec(), count })
}

fn find_multi_base_repeat(bases: &[u8], index: usize, repeat_length: usize, min_repeat_count: usize) -> Option<RepeatInfo> {
    if index + min_repeat_count * repeat_length - 1 >= bases.len() {
        return None;
    }
    let mut i = index + repeat_length;
    while i + (repeat_length - 1) < bases.len() {
        if bases[i..i + repeat_length] != bases[index..index + repeat_length] {
            break;
        }
        i += repeat_length;
    }
    let count = (i - index) / repeat_length;
    (count >= min_repeat_count).then(|| RepeatInfo { index, bases: bases[index..index + repeat_length].to_vec(), count })
}

fn find_dual_dual_repeat(bases: &[u8], index: usize, min_repeat_count: usize) -> Option<RepeatInfo> {
    if index + min_repeat_count * 4 - 1 >= bases.len() {
        return None;
    }
    if bases[index] != bases[index + 1] || bases[index + 2] != bases[index + 3] {
        return None;
    }
    let mut i = index + 4;
    while i + 3 < bases.len() {
        if bases[i..i + 4] != bases[index..index + 4] {
            break;
        }
        i += 4;
    }
    let count = (i - index) / 4;
    (count >= min_repeat_count).then(|| RepeatInfo { index, bases: bases[index..index + 4].to_vec(), count })
}

/// RepeatInfo.findRepeats: repeats of one to five bases with at least four units, the shorter unit first
pub fn find_repeats(bases: &[u8]) -> Vec<RepeatInfo> {
    let min = REPEAT_MIN_COUNT;
    let mut repeats = Vec::new();
    let mut index = 0usize;
    while index + REPEAT_MIN_COUNT <= bases.len() {
        let repeat = find_single_base_repeat(bases, index, min)
            .or_else(|| find_dual_base_repeat(bases, index, min))
            .or_else(|| find_multi_base_repeat(bases, index, 3, min))
            .or_else(|| find_multi_base_repeat(bases, index, 4, min))
            .or_else(|| find_multi_base_repeat(bases, index, 5, min))
            .or_else(|| find_dual_dual_repeat(bases, index, min));
        match repeat {
            Some(r) => {
                index += r.total_length();
                repeats.push(r);
            }
            None => index += 1,
        }
    }
    repeats
}

/// The bases removed from the alignment score for repeats: each repeat counts twice at most.
fn repeat_reduction(bases: &[u8]) -> i32 {
    find_repeats(bases).iter().map(|r| r.count.saturating_sub(2) * r.bases.len()).sum::<usize>() as i32
}

/// ReadFilters.calcRepeatTrimmedAlignmentScore: the alignment score less the repeated bases, scaled by the
/// square of score over aligned length. -1 when the aligned length is below the minimum.
pub fn calc_repeat_trimmed_alignment_score(rec: &bam::Record, read_bases: &[u8], aligned_length: i32, min_length: i32, apply_power_adjust: bool) -> f64 {
    if aligned_length < min_length {
        return -1.0;
    }
    let mut aligned: Vec<u8> = Vec::with_capacity(aligned_length as usize);
    let mut read_index = 0usize;
    for &c in rec.raw_cigar() {
        let (o, len) = (op(c), oplen(c) as usize);
        if o == OP_M {
            aligned.extend_from_slice(&read_bases[read_index..read_index + len]);
        }
        if consumes_read(o) {
            read_index += len;
        }
    }
    let reduction = repeat_reduction(&aligned);
    let alignment_score = match aux_int(rec, b"AS") {
        Some(v) => v as i32,
        None => aligned_length - aux_int(rec, b"NM").unwrap_or(0) as i32,
    };
    let mut adjusted = (alignment_score - reduction) as f64;
    if apply_power_adjust {
        let ratio = alignment_score as f64 / aligned_length as f64;
        adjusted *= ratio * ratio;
    }
    adjusted
}

pub fn calc_repeat_trimmed_alignment_score_read(read: &PrepRead, min_length: i32, apply_power_adjust: bool) -> f64 {
    calc_repeat_trimmed_alignment_score(&read.rec, &read.bases, read.aligned_base_length, min_length, apply_power_adjust)
}

pub fn above_repeat_trimmed_alignment_threshold(read: &PrepRead, min_length: i32, apply_power_adjust: bool) -> bool {
    calc_repeat_trimmed_alignment_score_read(read, min_length, apply_power_adjust) >= min_length as f64
}
