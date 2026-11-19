//! Evidence stage: ReadContextCounter and its per-read helpers (RawContext, QualityCalculator,
//! JitterMatch, Realignment, strand bias, fragment coordinates, edge distance).
//!
//! Not ported: the visualiser hooks.

use crate::candidate::{calc_event_penalty, calc_soft_clip_adjustment, calc_with_mnv_raw};
use crate::config::Config;
use crate::constants::*;
use crate::matcher::*;
use crate::quality::*;
use crate::read_context::*;
use crate::refseq::RefSequence;
use crate::repeat::RepeatInfo;
use crate::types::*;
use std::sync::Arc as Rc;

pub const INVALID_BASE_QUAL: f64 = -1.0;
pub const INVALID_INDEX: i32 = -1;

/// A read formed from two overlapping mates, with the originals retained.
pub struct FragmentData {
    pub first: Read,
    pub second: Read,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PosType {
    None,
    Aligned,
    SoftClip,
    LowQual,
    Deleted,
    Skipped,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ReadMatchType {
    AltSupportExact,
    AltSupportLowQualMismatches,
    RefSupport,
    NoAltRefMatch,
    Unrelated,
    NonCore,
    SoftClip,
    InSplit,
    MapQual,
    Chimeric,
    MaxCoverage,
}

/// RawContext.exceedsSoftClipLowBaseQual
pub(crate) fn exceeds_soft_clip_low_base_qual(quals: &[u8], start: i32, sc_len: i32) -> bool {
    if sc_len < MAX_SOFT_CLIP_LOW_QUAL_COUNT {
        return false;
    }
    let required = MIN_SOFT_CLIP_HIGH_QUAL_PERC * sc_len as f64;
    let mut above = 0;
    for i in start..start + sc_len {
        if quals[i as usize] as i32 >= MIN_SOFT_CLIP_MIN_BASE_QUAL {
            above += 1;
            if above as f64 >= required {
                return false;
            }
        }
    }
    sc_len - above >= MAX_SOFT_CLIP_LOW_QUAL_COUNT
}

/// RawContext.createFromRead: (readVariantIndex, positionType)
pub fn raw_context(variant: &SimpleVariant, read: &Read) -> (i32, PosType) {
    let n = read.cigar.len();
    let mut read_index = 0;
    let mut ref_base = read.start;
    for (i, el) in read.cigar.iter().enumerate() {
        let mut result: Option<(i32, PosType)> = None;
        match el.op {
            Op::H | Op::P => {}
            Op::S => {
                if i == 0 {
                    if variant.pos < read.start {
                        let idx = el.len - (read.start - variant.pos);
                        let t = if exceeds_soft_clip_low_base_qual(&read.quals, 0, el.len) {
                            PosType::LowQual
                        } else {
                            PosType::SoftClip
                        };
                        result = Some((idx, t));
                    }
                } else if i == n - 1 {
                    let unclipped_end = ref_base + el.len - 1;
                    if variant.pos <= unclipped_end {
                        let sc_start = read.quals.len() as i32 - el.len;
                        let idx = sc_start + (variant.pos - read.end) - 1;
                        let t = if exceeds_soft_clip_low_base_qual(&read.quals, sc_start, el.len) {
                            PosType::LowQual
                        } else {
                            PosType::SoftClip
                        };
                        result = Some((idx, t));
                    }
                }
                read_index += el.len;
            }
            Op::N => {
                let ref_position = ref_base - 1;
                if ref_position + el.len >= variant.pos {
                    result = Some((read_index - 1, PosType::Skipped));
                }
                ref_base += el.len;
            }
            Op::D => {
                let ref_position = ref_base - 1;
                if ref_position == variant.pos {
                    result = Some((read_index - 1, PosType::Aligned));
                } else if position_within(variant.pos, ref_position, ref_position + el.len) {
                    result = Some((read_index - 1, PosType::Deleted));
                }
                ref_base += el.len;
            }
            Op::I => {
                let ri = if read_index == 0 { 0 } else { read_index - 1 };
                if ref_base - 1 == variant.pos {
                    result = Some((ri, PosType::Aligned));
                }
                read_index += el.len;
            }
            Op::M | Op::Eq | Op::X => {
                let mut ref_end = ref_base + el.len - 1;
                if variant.is_indel() {
                    ref_end -= 1;
                }
                if ref_base <= variant.pos && variant.pos <= ref_end {
                    result = Some((read_index + variant.pos - ref_base, PosType::Aligned));
                }
                read_index += el.len;
                ref_base += el.len;
            }
        }
        if let Some(r) = result {
            return r;
        }
    }
    (-1, PosType::None)
}

/// CigarUtils.getReadIndexFromPosition(start, cigar, position, lastIfGapped=false, inferFromSoftClip=true)
pub fn read_index_from_position(alignment_start: i32, cigar: &[CigarEl], position: i32) -> i32 {
    let mut ref_pos = alignment_start;
    let mut index = 0;
    for el in cigar {
        if el.op == Op::S && ref_pos == alignment_start {
            ref_pos -= el.len;
        }
        let consumes_ref = el.op.consumes_ref() || el.op == Op::S;
        if consumes_ref && ref_pos + el.len - 1 >= position {
            if el.op.consumes_read() {
                index += position - ref_pos;
            } else {
                return INVALID_INDEX;
            }
            return index;
        }
        if consumes_ref {
            ref_pos += el.len;
        }
        if el.op.consumes_read() {
            index += el.len;
        }
    }
    INVALID_INDEX
}

/// SplitReadSegment: the part of a split (RNA) read, between `N` elements, which covers the variant.
pub struct SplitReadSegment {
    pub read_var_index: i32,
    /// indices in the original read
    pub segment_index_start: i32,
    #[cfg_attr(not(test), allow(dead_code))]
    pub segment_index_end: i32,
    pub read_bases: Vec<u8>,
    pub read_quals: Vec<u8>,
}

impl SplitReadSegment {
    pub fn length(&self) -> i32 {
        self.read_bases.len() as i32
    }

    /// SplitReadSegment.formSegment
    pub fn form_segment(read: &Read, position: i32, read_var_index: i32) -> Option<SplitReadSegment> {
        let mut read_index = 0;
        let mut ref_position = read.start - read.left_soft_clip();
        let mut segment_index_start = -1;
        let mut segment_index_end = -1;
        for el in &read.cigar {
            if el.op.consumes_ref() || el.op == Op::S {
                let ref_end_position = ref_position + el.len - 1;
                let position_within_section = position_within(position, ref_position, ref_end_position);
                if position_within_section && (el.op == Op::N || el.op == Op::D) {
                    return None;
                }
                if el.op == Op::N {
                    if ref_position > position {
                        break;
                    }
                    // reset base capture after any split
                    segment_index_start = -1;
                    segment_index_end = -1;
                } else if (el.op == Op::M || el.op == Op::S) && segment_index_start < 0 {
                    segment_index_start = read_index;
                }
                ref_position += el.len;
            }
            if el.op.consumes_read() {
                if segment_index_start >= 0 {
                    if segment_index_end == -1 {
                        segment_index_end = segment_index_start + el.len - 1;
                    } else {
                        segment_index_end += el.len;
                    }
                }
                read_index += el.len;
            }
        }
        if segment_index_start < 0 {
            return None;
        }
        // Arrays.copyOfRange pads with zeros past the end of the array
        let copy = |src: &[u8]| -> Vec<u8> {
            (segment_index_start..=segment_index_end).map(|i| src.get(i as usize).copied().unwrap_or(0)).collect()
        };
        Some(SplitReadSegment {
            read_var_index: read_var_index - segment_index_start,
            segment_index_start,
            segment_index_end,
            read_bases: copy(&read.bases),
            read_quals: copy(&read.quals),
        })
    }
}

/// Realignment.realignedReadIndexPosition
pub(crate) fn realigned_read_index_position(ctx: &VariantReadContext, read: &Read) -> i32 {
    let core_end_index = read_index_from_position(read.start, &read.cigar, ctx.core_position_end);
    if core_end_index == INVALID_INDEX {
        return INVALID_INDEX;
    }
    core_end_index - ctx.right_core_length()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Jitter {
    Shortened,
    Lengthened,
    None,
}

/// JitterMatch.hasValidBaseQuals (Illumina: no medium-qual band)
pub fn has_valid_base_quals(ctx: &VariantReadContext, alt_index_upper: i32, read: &Read, read_var_index: i32) -> bool {
    let flank_start = read_var_index - ctx.left_length();
    let flank_end = read_var_index + ctx.right_length() - 1;
    let core_start = flank_start + ctx.left_flank_length();
    let core_end = flank_end - ctx.right_flank_length();
    if core_start < 0 || core_end >= read.quals.len() as i32 {
        return false;
    }
    let lower = ctx.var_index;
    let upper = alt_index_upper + 1;
    let mut ctx_index = 0;
    for read_index in flank_start..=flank_end {
        let within_core = ctx_index >= ctx.core_index_start && ctx_index <= ctx.core_index_end;
        if within_core && ctx_index >= lower && ctx_index <= upper && is_uncertain_base_qual(read.quals[read_index as usize]) {
            return false;
        }
        if within_core && is_medium_base_qual(read.quals[read_index as usize] as i8) {
            return false;
        }
        ctx_index += 1;
    }
    true
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn has_jitter_match_type(
    repeat: &RepeatInfo, ctx: &VariantReadContext, alt_index_upper: i32, read_var_index: i32, read_bases: &[u8],
    read_quals: &[u8], jitter: Jitter, jitter_at_start: bool, jitter_count: i32,
) -> bool {
    let repeat_length = repeat.repeat_length() * jitter_count;
    let repeat_end_index = repeat.end_index();
    let offset = read_var_index - ctx.var_index;
    let n = read_bases.len() as i32;
    if repeat.index + offset < 0 || repeat_end_index + offset >= n {
        return false;
    }
    let mut flank_start = read_var_index - ctx.left_length();
    let flank_end = read_var_index + ctx.right_length() - 1;
    if jitter_at_start {
        if jitter == Jitter::Shortened {
            flank_start += repeat_length;
        } else {
            flank_start -= repeat_length;
        }
    }
    let lengthened_extra = if jitter == Jitter::Lengthened { repeat_length } else { 0 };
    let core_start = flank_start + ctx.left_flank_length();
    let core_end = flank_end - ctx.right_flank_length() + lengthened_extra;
    if core_start < 0 || core_end >= n {
        return false;
    }
    let mut ctx_index = 0;
    let mut read_index = flank_start;
    let mut index_adjusted = false;
    let mut permitted_low_qual_mismatches = 1;
    let lower = ctx.var_index;
    let upper = alt_index_upper + lengthened_extra;
    // checkJitter passes null base quals for Ultima
    let use_quals = !is_ultima();
    // Java: for(; readIndex <= flankReadIndexEnd; ++readIndex, ++readContextIndex)
    while read_index <= flank_end {
        if read_index < 0 {
            read_index += 1;
            ctx_index += 1;
            continue;
        }
        if !index_adjusted {
            if jitter == Jitter::Shortened && ctx_index == repeat_end_index - repeat_length {
                index_adjusted = true;
                read_index -= repeat_length;
            } else if jitter == Jitter::Lengthened && ctx_index == repeat_end_index + 1 {
                index_adjusted = true;
                ctx_index -= repeat_length;
            }
            if read_index < 0 || ctx_index < 0 {
                return false;
            }
        }
        if read_index >= n || ctx_index >= ctx.read_bases.len() as i32 {
            break;
        }
        let ctx_base = ctx.read_bases[ctx_index as usize];
        let within_core = ctx_index >= ctx.core_index_start && ctx_index <= ctx.core_index_end;
        let within_permitted = within_core && ctx_index >= lower && ctx_index <= upper;
        let q = read_quals[read_index as usize];
        if use_quals && within_permitted && is_uncertain_base_qual(q) {
            return false;
        }
        if use_quals && within_core && is_medium_base_qual(q as i8) {
            return false;
        }
        if read_bases[read_index as usize] != ctx_base {
            let mut permitted = false;
            if use_quals && q < matching_base_quality() {
                if within_core {
                    if permitted_low_qual_mismatches > 0 && !within_permitted {
                        permitted_low_qual_mismatches -= 1;
                        permitted = true;
                    }
                } else {
                    permitted = true;
                }
            }
            if !permitted {
                return false;
            }
        }
        read_index += 1;
        ctx_index += 1;
    }
    true
}

/// JitterMatch.checkJitter
pub fn check_jitter(ctx: &VariantReadContext, matcher: &Matcher, read: &Read, read_var_index: i32, is_msi: bool) -> Jitter {
    if ctx.all_repeats.is_empty() {
        return Jitter::None;
    }
    let double_count = if is_msi { MSI_DOUBLE_JITTER_REPEAT_COUNT } else { DEFAULT_DOUBLE_JITTER_REPEAT_COUNT };
    let check_double = ctx.max_repeat.as_ref().map_or(false, |r| r.count >= double_count);
    for repeat in &ctx.all_repeats {
        for jitter in [Jitter::Shortened, Jitter::Lengthened] {
            for at_start in [true, false] {
                let counts: &[i32] = if check_double { &[1, 2] } else { &[1] };
                for &count in counts {
                    if has_jitter_match_type(
                        repeat, ctx, matcher.alt_index_upper, read_var_index, &read.bases, &read.quals, jitter, at_start,
                        count,
                    ) {
                        return jitter;
                    }
                }
            }
        }
    }
    Jitter::None
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Realigned {
    None,
    Exact,
    LowQualMismatches,
    Shortened,
    Lengthened,
}

/// Realignment.checkRealignment
pub(crate) fn check_realignment(
    ctx: &VariantReadContext, matcher: &mut Matcher, read: &Read, read_index: i32, realigned_read_index: i32,
    split: Option<&SplitReadSegment>, is_msi: bool,
) -> Realigned {
    let mut realigned_read_index = realigned_read_index;
    if read_index == realigned_read_index {
        return Realigned::None;
    }
    if read_index < 0 && !ctx.variant.is_delete() {
        return Realigned::None;
    }
    let offset = if ctx.variant.is_delete() { 0 } else { realigned_read_index - read_index };
    let (bases, quals): (&[u8], &[u8]) = match split {
        Some(seg) => {
            // the index stays relative to the segment for the jitter check below, which Java makes
            // against the whole read
            realigned_read_index -= seg.segment_index_start;
            (&seg.read_bases, &seg.read_quals)
        }
        None => (&read.bases, &read.quals),
    };
    if realigned_read_index < 0 || realigned_read_index >= bases.len() as i32 {
        return Realigned::None;
    }
    matcher.set_realignment_index_offset(offset);
    let info = matcher.determine(ctx, bases, Some(quals), realigned_read_index, true);
    matcher.clear_realignment_index_offset();
    if matches!(info.kind, MatchType::Full | MatchType::PartialCore) {
        return if info.exact { Realigned::Exact } else { Realigned::LowQualMismatches };
    }
    match check_jitter(ctx, matcher, read, realigned_read_index, is_msi) {
        Jitter::Shortened => Realigned::Shortened,
        Jitter::Lengthened => Realigned::Lengthened,
        Jitter::None => Realigned::None,
    }
}

#[derive(Default, Clone, Debug)]
pub struct SupportCounts {
    pub full: i32,
    pub partial_core: i32,
    pub core: i32,
    pub realigned: i32,
    pub ref_: i32,
    pub total: i32,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Support {
    Full,
    PartialCore,
    Core,
    Realigned,
    Ref,
}

impl SupportCounts {
    fn add(&mut self, support: Option<Support>, count: i32) {
        self.total += count;
        match support {
            Some(Support::Full) => self.full += count,
            Some(Support::PartialCore) => self.partial_core += count,
            Some(Support::Core) => self.core += count,
            Some(Support::Realigned) => self.realigned += count,
            Some(Support::Ref) => self.ref_ += count,
            None => {}
        }
    }
    pub fn alt_support(&self) -> i32 {
        self.full + self.partial_core + self.core + self.realigned
    }
    pub fn strong_support(&self) -> i32 {
        self.full + self.partial_core + self.realigned
    }
    pub fn to_array(&self) -> [i32; 6] {
        [self.full, self.partial_core, self.core, self.realigned, self.ref_, self.total]
    }
}

#[derive(Default, Clone, Debug)]
pub struct QualCounters {
    pub alt_seq_tech_base_quality_total: i32,
    pub recalibrated_base_quality_total: i32,
    pub alt_recalibrated_base_quality_total: i32,
    pub strong_alt_recalibrated_base_quality_total: f64,
    pub strong_alt_medium_recalibrated_base_quality_total: f64,
    pub alt_final_base_quality_total: f64,
    pub alt_medium_final_base_quality_total: f64,
    pub map_quality_total: i32,
    pub alt_map_quality_total: i32,
    pub alt_final_map_quality_total: f64,
    pub low_qual_alt_support_count: i32,
}

impl QualCounters {
    fn update(&mut self, q: &QualityScores, map_quality: i32, kind: MatchType) {
        self.recalibrated_base_quality_total += jround(q.recalibrated_base_quality) as i32;
        self.map_quality_total += map_quality;
        if kind.full_alt_support() {
            self.alt_final_base_quality_total += q.final_base_quality;
            self.strong_alt_recalibrated_base_quality_total += q.recalibrated_base_quality;
            if q.is_medium_qual {
                self.strong_alt_medium_recalibrated_base_quality_total += q.recalibrated_base_quality;
                self.alt_medium_final_base_quality_total += q.final_base_quality;
            }
            self.alt_final_map_quality_total += q.final_map_quality as f64;
        }
        if kind.supports_alt() {
            self.alt_recalibrated_base_quality_total += jround(q.recalibrated_base_quality) as i32;
            self.alt_seq_tech_base_quality_total += jround(q.seq_tech_base_quality) as i32;
            self.alt_map_quality_total += map_quality;
        }
    }
}

#[derive(Default, Clone, Debug)]
pub struct StrandBias {
    pub forward: i32,
    pub reverse: i32,
}

impl StrandBias {
    fn add(&mut self, fwd: bool) {
        if fwd {
            self.forward += 1;
        } else {
            self.reverse += 1;
        }
    }
    pub fn depth(&self) -> i32 {
        self.forward + self.reverse
    }
    pub fn bias(&self) -> f64 {
        let depth = (self.forward + self.reverse) as f64;
        if depth > 0.0 {
            self.forward as f64 / depth
        } else {
            0.5
        }
    }
    pub fn min_bias(&self) -> f64 {
        let b = self.bias();
        b.min(1.0 - b)
    }
    /// registerFragment
    fn register_fragment(&mut self, read: &Read) {
        let fwd = !read.negative_strand();
        if !read.paired() {
            self.add(fwd);
            return;
        }
        if read.consensus == ConsensusType::Dual {
            self.forward += 1;
            self.reverse += 1;
            return;
        }
        let mate_fwd = !read.flag(FLAG_MATE_REVERSE);
        let first = read.flag(FLAG_FIRST);
        let first_fwd = if first { fwd } else { mate_fwd };
        let second_fwd = if !first { fwd } else { mate_fwd };
        if first_fwd != second_fwd {
            self.add(first_fwd);
        }
    }
}

/// FragmentCoords: unique 5' fragment coordinates of alt reads, capped.
#[derive(Clone, Debug)]
pub struct FragmentCoords {
    pub(crate) max: usize,
    pub(crate) lower: Vec<(i32, i32, bool)>,
    pub(crate) upper: Vec<(i32, i32, bool)>,
}

fn add_unique(v: &mut Vec<(i32, i32, bool)>, x: (i32, i32, bool)) {
    if !v.contains(&x) {
        v.push(x);
    }
}

/// SamRecordUtils.getFivePrimeUnclippedPosition(readStart, cigarStr, forwardStrand)
fn five_prime_unclipped_from_cigar(start: i32, cigar: &str, forward: bool) -> i32 {
    let mut els: Vec<(i32, u8)> = Vec::new();
    let mut num = 0;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            num = num * 10 + (b - b'0') as i32;
        } else {
            els.push((num, b));
            num = 0;
        }
    }
    if forward {
        match els.first() {
            Some(&(l, b'S')) => start - l,
            _ => start,
        }
    } else {
        let ref_len: i32 = els.iter().filter(|(_, o)| matches!(o, b'M' | b'D' | b'N' | b'=' | b'X')).map(|(l, _)| l).sum();
        let end = start + ref_len - 1;
        match els.last() {
            Some(&(l, b'S')) => end + l,
            _ => end,
        }
    }
}

impl FragmentCoords {
    fn new(max: usize) -> Self {
        FragmentCoords { max, lower: Vec::new(), upper: Vec::new() }
    }
    pub fn min_count(&self) -> i32 {
        self.lower.len().min(self.upper.len()) as i32
    }
    fn add_read(&mut self, read: &Read, mate: Option<&Read>) {
        if self.lower.len() >= self.max && self.upper.len() >= self.max {
            return;
        }
        let fwd = !read.negative_strand();
        let five = if fwd { read.start - read.left_soft_clip() } else { read.end + read.right_soft_clip() };
        let five_coord = (-1, five, fwd);
        if !read.paired() || read.flag(FLAG_MATE_UNMAPPED) {
            let three = if fwd { read.end + read.right_soft_clip() } else { read.start - read.left_soft_clip() };
            let three_coord = (-1, three, fwd);
            add_unique(&mut self.lower, if five < three { five_coord } else { three_coord });
            add_unique(&mut self.upper, if five < three { three_coord } else { five_coord });
            return;
        }
        let (mate_pos, mate_fwd) = match mate {
            Some(m) => {
                let mf = !m.negative_strand();
                (if mf { m.start - m.left_soft_clip() } else { m.end + m.right_soft_clip() }, mf)
            }
            None => {
                let mf = !read.flag(FLAG_MATE_REVERSE);
                let Some(mc) = &read.mate_cigar else {
                    add_unique(&mut self.lower, five_coord);
                    return;
                };
                (five_prime_unclipped_from_cigar(read.mate_start, mc, mf), mf)
            }
        };
        let same_chr = read.tid == read.mate_tid;
        let mate_coord = if same_chr { (-1, mate_pos, mate_fwd) } else { (read.mate_tid, mate_pos, mate_fwd) };
        if (same_chr && five <= mate_pos) || (!same_chr && read.tid < read.mate_tid) {
            add_unique(&mut self.lower, five_coord);
            add_unique(&mut self.upper, mate_coord);
        } else {
            add_unique(&mut self.lower, mate_coord);
            add_unique(&mut self.upper, five_coord);
        }
    }
}

#[derive(Clone, Debug)]
pub struct ReadEdgeDistance {
    variant_position: i32,
    pub max_distance_from_edge: f64,
    pub max_alt_distance_from_edge: f64,
    total: f64,
    total_alt: f64,
    updates: i32,
    updates_alt: i32,
}

fn calc_distance_from_read_edge(pos: i32, read: &Read) -> i32 {
    if !position_within(pos, read.start, read.end) {
        return -1;
    }
    (pos - read.start).min(read.end - pos)
}

impl ReadEdgeDistance {
    fn new(variant_position: i32) -> Self {
        ReadEdgeDistance {
            variant_position,
            max_distance_from_edge: 0.0,
            max_alt_distance_from_edge: 0.0,
            total: 0.0,
            total_alt: 0.0,
            updates: 0,
            updates_alt: 0,
        }
    }
    pub fn avg_distance_from_edge(&self) -> f64 {
        if self.updates > 0 {
            self.total / self.updates as f64
        } else {
            0.0
        }
    }
    pub fn avg_alt_distance_from_edge(&self) -> f64 {
        if self.updates_alt > 0 {
            self.total_alt / self.updates_alt as f64
        } else {
            0.0
        }
    }
    fn update(&mut self, read: &Read, fragment: Option<&FragmentData>, alt_support: bool) {
        if read.has_op(Op::N) {
            return;
        }
        let pos = self.variant_position;
        let mut min_distance = -1;
        let read_length;
        if let Some(f) = fragment {
            let d1 = calc_distance_from_read_edge(pos, &f.first);
            let d2 = calc_distance_from_read_edge(pos, &f.second);
            if d1 > -1 && d2 > -1 {
                min_distance = calc_distance_from_read_edge(pos, read);
            } else if d1 > -1 {
                min_distance = d1;
            } else if d2 > -1 {
                min_distance = d2;
            }
            read_length = f.first.len().max(f.second.len());
        } else {
            min_distance = calc_distance_from_read_edge(pos, read);
            read_length = read.len();
        }
        let perc = min_distance as f64 / read_length as f64;
        self.updates += 1;
        self.max_distance_from_edge = jmax(perc, self.max_distance_from_edge);
        if alt_support {
            self.max_alt_distance_from_edge = jmax(perc, self.max_alt_distance_from_edge);
            self.total_alt += perc;
            self.updates_alt += 1;
        }
        self.total += perc;
    }
}

#[derive(Default, Clone, Debug)]
pub struct JitterData {
    pub lengthened: i32,
    pub shortened: i32,
    pub valid_qual_full_support: i32,
    pub qual_boost: f64,
    pub filter_on_noise: bool,
    pub hard_filter_on_noise: bool,
}

impl JitterData {
    fn update(&mut self, j: Jitter) {
        match j {
            Jitter::Lengthened => self.lengthened += 1,
            Jitter::Shortened => self.shortened += 1,
            Jitter::None => {}
        }
    }
}

pub struct ReadContextCounter {
    pub tier: Tier,
    pub ctx: Rc<VariantReadContext>,
    pub matcher: Matcher,
    pub max_coverage: i32,
    pub min_number_of_events: i32,
    pub is_mnv: bool,
    pub is_indel: bool,
    pub qual_cache: QualCache,
    pub qualities: SupportCounts,
    pub counts: SupportCounts,
    pub simple_alt_matches: i32,
    pub high_qual_strong_support: i32,
    pub medium_qual_strong_support: i32,
    pub allow_uncertain_core_bases: bool,
    pub uncertain_core_base_count: i32,
    pub alt_fragment_strand_bias: StrandBias,
    pub non_alt_fragment_strand_bias: StrandBias,
    pub alt_read_strand_bias: StrandBias,
    pub non_alt_read_strand_bias: StrandBias,
    pub jitter: JitterData,
    pub improper_pair_count: i32,
    pub qual_counters: QualCounters,
    pub max_candidate_delete_length: i32,
    pub read_edge_distance: ReadEdgeDistance,
    pub max_position_vs_read_start: i32,
    pub non_alt_nm_count_total: i64,
    pub local_phase_sets: Vec<i32>,
    pub lps_counts: Vec<i32>,
    pub consensus_type_counts: Option<[i32; 6]>,
    pub fragment_coords: FragmentCoords,
    /// mFragmentLengthData (-write_frag_lengths): (length, [ref count, alt count]) in insertion order
    pub fragment_length_data: Option<Vec<(i32, [i32; 2])>>,
    pub tumor_qual_probability: f64,
    pub map_qual_factor: f64,
    /// mUltimaData: only for Ultima
    pub ultima: Option<crate::ultima::UltimaVariantData>,
    /// the counter ReadContextCounterFactory builds when the Ultima models cannot be built for a variant
    /// near the start of a chromosome or on a non-standard one: it has no matcher and counts nothing
    pub no_matcher: bool,
}

/// FragmentLengthCounts.addLength: lengths of 1 to 500 are counted, as (length, [ref, alt]) in the
/// order the lengths were first seen
pub fn add_fragment_length(lengths: &mut Vec<(i32, [i32; 2])>, length: i32, is_alt: bool) {
    if length <= 0 || length > 500 {
        return;
    }
    let index = match lengths.iter().position(|(l, _)| *l == length) {
        Some(i) => i,
        None => {
            lengths.push((length, [0, 0]));
            lengths.len() - 1
        }
    };
    lengths[index].1[is_alt as usize] += 1;
}

impl ReadContextCounter {
    pub fn new(
        ctx: Rc<VariantReadContext>, tier: Tier, max_coverage: i32, min_number_of_events: i32, qc: &QualityCalculator,
        is_reference_sample: bool,
    ) -> Self {
        let variant = &ctx.variant;
        // Java leaves the matcher null for a context marked invalid (append mode only)
        let matcher =
            if ctx.marked_invalid { Matcher::placeholder(is_reference_sample) } else { Matcher::new(&ctx, true, is_reference_sample) };
        let qual_cache = QualCache::new(&ctx, qc);
        let is_indel = variant.is_indel();
        let mut allow_uncertain = false;
        if is_indel {
            let alt_base_diff = match (&qual_cache.msi_indel_repeat, qual_cache.uses_msi_indel_error_qual()) {
                (Some(repeat), true) => {
                    let alt_bases: &[u8] = if variant.is_delete() { &variant.ref_[1..] } else { &variant.alt[1..] };
                    implied_alt_change(variant.is_delete(), alt_bases, repeat).abs()
                }
                _ => variant.indel_length_abs(),
            };
            allow_uncertain = alt_base_diff >= INDEL_UNCERTAIN_BASE_REPEAT_MIN;
        }
        let max_position_vs_read_start =
            if variant.is_delete() { variant.pos + variant.indel_length().abs() } else { variant.pos };
        // calcAdjustedVariantPosition(position, indelLengthAbs + 1)
        let indel_len = variant.indel_length_abs() + 1;
        let adjusted_pos = variant.pos + indel_len / 2;
        ReadContextCounter {
            tier,
            is_mnv: variant.is_mnv(),
            is_indel,
            matcher,
            max_coverage,
            min_number_of_events,
            qual_cache,
            qualities: SupportCounts::default(),
            counts: SupportCounts::default(),
            simple_alt_matches: 0,
            high_qual_strong_support: 0,
            medium_qual_strong_support: 0,
            allow_uncertain_core_bases: allow_uncertain,
            uncertain_core_base_count: 0,
            alt_fragment_strand_bias: StrandBias::default(),
            non_alt_fragment_strand_bias: StrandBias::default(),
            alt_read_strand_bias: StrandBias::default(),
            non_alt_read_strand_bias: StrandBias::default(),
            jitter: JitterData { qual_boost: 1.0, ..Default::default() },
            improper_pair_count: 0,
            qual_counters: QualCounters::default(),
            max_candidate_delete_length: 0,
            read_edge_distance: ReadEdgeDistance::new(adjusted_pos),
            max_position_vs_read_start,
            non_alt_nm_count_total: 0,
            local_phase_sets: Vec::new(),
            lps_counts: Vec::new(),
            consensus_type_counts: None,
            fragment_coords: FragmentCoords::new(REQUIRED_UNIQUE_FRAG_COORDS_2 as usize),
            fragment_length_data: None,
            tumor_qual_probability: 0.0,
            map_qual_factor: 0.0,
            ultima: if is_ultima() && !ctx.marked_invalid { Some(crate::ultima::UltimaVariantData::new(&ctx)) } else { None },
            no_matcher: ctx.marked_invalid,
            ctx,
        }
    }

    pub fn variant(&self) -> &SimpleVariant {
        &self.ctx.variant
    }
    pub fn position(&self) -> i32 {
        self.ctx.variant.pos
    }
    pub fn alt_support(&self) -> i32 {
        self.counts.alt_support()
    }
    pub fn strong_alt_support(&self) -> i32 {
        self.counts.strong_support()
    }
    pub fn depth(&self) -> i32 {
        self.counts.total
    }
    pub fn vaf(&self) -> f64 {
        if self.counts.total == 0 {
            0.0
        } else {
            self.counts.alt_support() as f64 / self.counts.total as f64
        }
    }
    /// tumorQuality: `int *= double` truncates toward zero
    pub fn tumor_quality(&self) -> f64 {
        let mut total = self.qualities.full + self.qualities.partial_core + self.qualities.realigned;
        if !self.jitter.filter_on_noise {
            total = (total as f64 * self.jitter.qual_boost) as i32;
        }
        total as f64
    }
    pub fn use_msi_error_rate(&self) -> bool {
        self.qual_cache.uses_msi_indel_error_qual()
    }
    pub fn exceeds_max_coverage(&self) -> bool {
        self.counts.total >= self.max_coverage
    }
    pub fn average_alt_seq_tech_base_quality(&self) -> f64 {
        let n = self.counts.alt_support() + self.qual_counters.low_qual_alt_support_count;
        if n > 0 {
            self.qual_counters.alt_seq_tech_base_quality_total as f64 / n as f64
        } else {
            0.0
        }
    }
    pub fn average_alt_recalibrated_base_quality(&self) -> f64 {
        let n = self.counts.alt_support();
        if n > 0 {
            self.qual_counters.alt_recalibrated_base_quality_total as f64 / n as f64
        } else {
            0.0
        }
    }
    pub fn is_long_indel(&self) -> bool {
        self.ctx.variant.indel_length_abs() >= LONG_INSERT_LENGTH
    }

    fn has_uncertain_core_bases(&mut self, read: &Read, read_index: i32, split: Option<&SplitReadSegment>) -> bool {
        if self.allow_uncertain_core_bases {
            return false;
        }
        let uncertain = match split {
            Some(seg) => self.matcher.has_uncertain_core_bases(&self.ctx, &seg.read_quals, seg.read_var_index),
            None => self.matcher.has_uncertain_core_bases(&self.ctx, &read.quals, read_index),
        };
        if uncertain {
            self.uncertain_core_base_count += 1;
            return true;
        }
        false
    }

    fn count_consensus_type(&mut self, read: &Read, read_var_index: i32, supports_variant: bool) {
        if self.consensus_type_counts.is_none() {
            if !read.has_consensus_tag {
                return;
            }
            let mut c = [0; 6];
            c[0] = self.depth();
            c[3] = self.strong_alt_support();
            self.consensus_type_counts = Some(c);
        }
        let consensus =
            if is_sbx() { sbx_consensus_type(&self.ctx, read_var_index, read) } else { read.consensus };
        let c = self.consensus_type_counts.as_mut().unwrap();
        let t = consensus as usize;
        c[t] += 1;
        if supports_variant {
            c[t + 3] += 1;
        }
    }

    fn register_read_support(
        &mut self, read: &Read, support: Option<Support>, read_var_index: i32, combined_quality: f64, seq_tech_qual: f64,
        is_medium: bool,
    ) {
        let strong = matches!(support, Some(Support::Full | Support::PartialCore | Support::Realigned));
        self.count_consensus_type(read, read_var_index, strong);
        self.counts.add(support, 1);
        self.qualities.add(support, combined_quality as i32);
        if strong {
            if is_high_base_qual(seq_tech_qual) {
                self.high_qual_strong_support += 1;
            } else if is_medium {
                self.medium_qual_strong_support += 1;
            }
        }
        if let Some(lengths) = self.fragment_length_data.as_mut() {
            let supports_alt = strong || support == Some(Support::Core);
            if support == Some(Support::Ref) || supports_alt {
                // SamRecordUtils.inferredInsertSizeAbs
                let length = if read.paired() { read.insert_size.abs() } else { (read.end - read.start + 1).abs() };
                add_fragment_length(lengths, length, supports_alt);
            }
        }
        if matches!(support, Some(Support::Full | Support::PartialCore)) && (read.improper_pair() || read.supplementary()) {
            self.improper_pair_count += 1;
        }
        if support == Some(Support::Full) {
            if let Some(ultima) = self.ultima.as_mut() {
                ultima.add_read_support_info(&self.ctx, read, read_var_index);
            }
        }
    }

    /// StrandBiasData.registerRead(record, fragment, counter)
    fn register_read_strand(&self, is_alt: bool, read: &Read, fragment: Option<&FragmentData>) -> (i32, i32) {
        let mut fwd_rev = (0, 0);
        let mut add = |r: &Read| {
            if r.negative_strand() {
                fwd_rev.1 += 1;
            } else {
                fwd_rev.0 += 1;
            }
        };
        let Some(f) = fragment else {
            add(read);
            return fwd_rev;
        };
        if f.first.negative_strand() == f.second.negative_strand() {
            return fwd_rev;
        }
        let pos = self.position();
        let within_first = position_within(pos, f.first.start, f.first.end);
        let within_second = position_within(pos, f.second.start, f.second.end);
        let supports = |r: &Read| -> bool {
            let (idx, _) = raw_context(self.variant(), r);
            if idx < 0 || !self.matcher.covers_variant(&self.ctx, r.len(), idx) {
                return false;
            }
            let kind = self.matcher.determine_read_match(&self.ctx, r, idx).kind;
            if is_alt {
                kind.supports_alt()
            } else {
                kind == MatchType::Ref
            }
        };
        if within_first && !within_second {
            if supports(&f.first) {
                add(&f.first);
            }
        } else if !within_first && within_second {
            if supports(&f.second) {
                add(&f.second);
            }
        } else {
            if supports(&f.first) {
                add(&f.first);
            }
            if supports(&f.second) {
                add(&f.second);
            }
        }
        fwd_rev
    }

    fn count_alt_support_metrics(&mut self, read: &Read, fragment: Option<&FragmentData>) {
        let (f, r) = self.register_read_strand(true, read, fragment);
        self.alt_read_strand_bias.forward += f;
        self.alt_read_strand_bias.reverse += r;
        self.alt_fragment_strand_bias.register_fragment(read);
        match fragment {
            Some(fd) => self.fragment_coords.add_read(&fd.first, Some(&fd.second)),
            None => self.fragment_coords.add_read(read, None),
        }
    }

    /// ReadContextCounter.processRead
    pub fn process_read(
        &mut self, read: &Read, number_of_events: i32, fragment: Option<&FragmentData>, qc: &QualityCalculator,
    ) -> ReadMatchType {
        use ReadMatchType::*;
        if self.no_matcher {
            return Unrelated;
        }
        if self.exceeds_max_coverage() {
            return MaxCoverage;
        }
        if self.tier != Tier::Hotspot && read.mapq < EVIDENCE_MIN_MAP_QUAL {
            return MapQual;
        }
        let high_depth = qc.cfg.high_depth_mode;
        if high_depth && fragment.is_none() && is_chimeric_read(read) {
            return Chimeric;
        }
        // belowQualThreshold
        let below_qual_threshold =
            |qual_cache: &QualCache, q: f64| high_depth && !qual_cache.uses_msi_indel_error_qual() && !is_high_base_qual(q);
        let ctx = self.ctx.clone();
        let variant = &ctx.variant;
        let (read_var_index, pos_type) = raw_context(variant, read);
        let split_segment =
            if read.has_op(Op::N) { SplitReadSegment::form_segment(read, variant.pos, read_var_index) } else { None };
        let split = split_segment.as_ref();
        // coversVariant: with a split segment Java ignores the index passed in
        let covers = |matcher: &Matcher, index: i32| match split {
            Some(seg) => matcher.covers_variant(&ctx, seg.length(), seg.read_var_index),
            None => matcher.covers_variant(&ctx, read.len(), index),
        };
        // determineReadContextMatch
        let determine = |matcher: &Matcher| match split {
            Some(seg) => matcher.determine(&ctx, &seg.read_bases, Some(&seg.read_quals), seg.read_var_index, false),
            None => matcher.determine_read_match(&ctx, read, read_var_index),
        };
        if pos_type == PosType::LowQual {
            return SoftClip;
        }
        let mut check_realigned = false;
        let mut realigned_read_index: Option<i32> = None;
        if read_var_index < 0 {
            if variant.is_delete() && self.max_position_vs_read_start >= read.start - read.left_soft_clip() {
                let idx = realigned_read_index_position(&ctx, read);
                realigned_read_index = Some(idx);
                check_realigned = idx != INVALID_INDEX;
            }
            if !check_realigned {
                return Unrelated;
            }
        }
        if pos_type == PosType::Skipped {
            return InSplit;
        }
        let variant_covered = covers(&self.matcher, read_var_index);
        if !variant_covered && !check_realigned {
            return NonCore;
        }

        let mut adjusted_events = number_of_events as f64;
        if self.is_mnv {
            adjusted_events = calc_with_mnv_raw(number_of_events, &variant.ref_, &variant.alt) as f64;
        }
        if variant.ref_.len().max(variant.alt.len()) as f64 <= SC_READ_EVENTS_FACTOR {
            adjusted_events += calc_soft_clip_adjustment(read.left_soft_clip() + read.right_soft_clip());
        }
        adjusted_events = jmax(self.min_number_of_events as f64, adjusted_events);

        let mut match_info = NO_MATCH;
        let mut match_type = MatchType::None;
        let mut seq_tech_qual = 0.0;
        let mut combined_quality = 0.0;
        let mut scores: Option<QualityScores> = None;

        if variant_covered {
            if self.has_uncertain_core_bases(read, read_var_index, split) {
                return NonCore;
            }
            let q = qc.calculate_quality_scores(self, read_var_index, read, adjusted_events);
            if !q.valid() {
                return Unrelated;
            }
            seq_tech_qual = q.seq_tech_base_quality;
            if below_qual_threshold(&self.qual_cache, seq_tech_qual) {
                if pos_type != PosType::Deleted {
                    let info = determine(&self.matcher);
                    if info.kind.supports_alt() {
                        self.count_alt_support_metrics(read, fragment);
                        self.qual_counters.alt_seq_tech_base_quality_total += jround(q.seq_tech_base_quality) as i32;
                        self.qual_counters.low_qual_alt_support_count += 1;
                    }
                }
                return Unrelated;
            }
            combined_quality = q.combined_quality;
            if pos_type != PosType::Deleted {
                match_info = determine(&self.matcher);
                match_type = match_info.kind;
            }
            if match_type.supports_alt() {
                let support = match match_type {
                    MatchType::Full => Support::Full,
                    MatchType::PartialCore => Support::PartialCore,
                    MatchType::Realigned => Support::Realigned,
                    _ => Support::Core,
                };
                self.register_read_support(read, Some(support), read_var_index, combined_quality, seq_tech_qual, q.is_medium_qual);
                self.qual_counters.update(&q, read.mapq, match_type);
                self.read_edge_distance.update(read, fragment, match_type.full_alt_support());
                self.count_alt_support_metrics(read, fragment);
                if support == Support::Full
                    && has_valid_base_quals(&ctx, self.matcher.alt_index_upper, read, read_var_index)
                {
                    self.jitter.valid_qual_full_support += 1;
                }
                return if match_info.exact && match_type.full_alt_support() {
                    AltSupportExact
                } else {
                    AltSupportLowQualMismatches
                };
            }
            scores = Some(q);
        }

        let mut realigned_type = Realigned::None;
        if match_type != MatchType::Ref || (is_sbx() && !self.is_indel && ctx.has_indel_in_core) {
            let idx = *realigned_read_index.get_or_insert_with(|| realigned_read_index_position(&ctx, read));
            let can_realign = idx != INVALID_INDEX && covers(&self.matcher, idx);
            if can_realign {
                let is_msi = self.qual_cache.is_msi_sample_and_variant;
                realigned_type = check_realignment(&ctx, &mut self.matcher, read, read_var_index, idx, split, is_msi);
            }
            if realigned_type != Realigned::None {
                if self.has_uncertain_core_bases(read, idx, split) {
                    return Unrelated;
                }
                let q = qc.calculate_quality_scores(self, idx, read, adjusted_events);
                if !q.valid() {
                    return Unrelated;
                }
                seq_tech_qual = q.seq_tech_base_quality;
                combined_quality = q.combined_quality;
                if below_qual_threshold(&self.qual_cache, seq_tech_qual) {
                    return Unrelated;
                }
                if matches!(realigned_type, Realigned::Exact | Realigned::LowQualMismatches) {
                    self.register_read_support(read, Some(Support::Realigned), read_var_index, combined_quality, seq_tech_qual, q.is_medium_qual);
                    self.qual_counters.update(&q, read.mapq, MatchType::Realigned);
                    self.read_edge_distance.update(read, fragment, true);
                    self.count_alt_support_metrics(read, fragment);
                    return if realigned_type == Realigned::Exact { AltSupportExact } else { AltSupportLowQualMismatches };
                }
                if realigned_type == Realigned::Shortened {
                    self.jitter.update(Jitter::Shortened);
                } else if realigned_type == Realigned::Lengthened {
                    self.jitter.update(Jitter::Lengthened);
                }
                scores = Some(q);
            } else if !variant_covered || read_var_index < 0 {
                return if read_var_index < 0 { Unrelated } else { NonCore };
            }
        }

        if realigned_type == Realigned::None {
            let j = check_jitter(&ctx, &self.matcher, read, read_var_index, self.qual_cache.is_msi_sample_and_variant);
            self.jitter.update(j);
        }

        let mut support = None;
        if match_type == MatchType::Ref {
            support = Some(Support::Ref);
        } else if match_type == MatchType::SimpleAlt && pos_type != PosType::SoftClip {
            self.simple_alt_matches += 1;
        }

        if match_type == MatchType::None && variant.is_insert() && read_var_index - variant.indel_length() < 0 {
            return Unrelated;
        }

        // Java dereferences qualityScores here; it is always set on the paths that reach this point
        let q = scores.expect("quality scores");
        self.qual_counters.update(&q, read.mapq, match_type);
        self.non_alt_fragment_strand_bias.register_fragment(read);
        let (f, r) = self.register_read_strand(false, read, fragment);
        self.non_alt_read_strand_bias.forward += f;
        self.non_alt_read_strand_bias.reverse += r;
        self.register_read_support(read, support, read_var_index, combined_quality, seq_tech_qual, q.is_medium_qual);
        self.read_edge_distance.update(read, fragment, false);
        self.non_alt_nm_count_total += number_of_events as i64;

        if match_type == MatchType::Ref {
            return RefSupport;
        }
        if !is_simple_alt_match(variant, read, read_var_index) {
            return NoAltRefMatch;
        }
        Unrelated
    }

    /// addLocalPhaseSet
    pub fn add_local_phase_set(&mut self, lps: i32, read_count: i32, alloc_count: f64) {
        let mut index = 0;
        while index < self.lps_counts.len() {
            if read_count as f64 + alloc_count > self.lps_counts[index] as f64 {
                break;
            }
            index += 1;
        }
        self.local_phase_sets.insert(index, lps);
        self.lps_counts.insert(index, read_count + alloc_count as i32);
    }
}

/// QualityCalculator.calculateQualityScores and helpers that need the counter.
impl<'a> QualityCalculator<'a> {
    pub fn calculate_quality_scores(
        &self, rc: &ReadContextCounter, read_var_index: i32, read: &Read, number_of_events: f64,
    ) -> QualityScores {
        let ctx = &rc.ctx;
        let uses_msi = rc.qual_cache.uses_msi_indel_error_qual();
        let has_valid_msi_jitter_quals =
            uses_msi && has_valid_base_quals(ctx, rc.matcher.alt_index_upper, read, read_var_index);

        // calcSequencingTechBaseQuality
        let artefact_qual: f64 =
            ctx.artefact.as_ref().map_or(INVALID_BASE_QUAL, |a| a.find_applicable_base_qual(read, read_var_index) as f64);
        let seq_tech = if is_ultima() {
            let models = rc.ultima.as_ref().and_then(|u| u.qual_models.as_ref()).expect("Ultima qual models");
            models.calculate_qual(uses_msi, rc.qual_cache.msi_indel_error_qual, rc.is_indel, read_var_index, read)
        } else if rc.is_indel {
            let mut msi_qual = rc.qual_cache.msi_indel_error_qual as f64;
            let mut early: Option<f64> = None;
            if uses_msi {
                if !has_valid_msi_jitter_quals {
                    msi_qual /= 2.0;
                }
                if rc.allow_uncertain_core_bases {
                    early = Some(msi_qual);
                }
            }
            if let Some(v) = early {
                v
            } else if uses_msi && artefact_qual != INVALID_BASE_QUAL {
                jmin(artefact_qual, msi_qual)
            } else {
                let avg = if is_sbx() {
                    sbx_indel_core_quality(ctx, read, read_var_index)
                } else {
                    average_core_quality(ctx, read, read_var_index)
                };
                if uses_msi {
                    jmin(avg, msi_qual)
                } else if artefact_qual != INVALID_BASE_QUAL {
                    jmin(avg, artefact_qual)
                } else {
                    avg
                }
            }
        } else if artefact_qual != INVALID_BASE_QUAL {
            artefact_qual
        } else if is_sbx() && ctx.has_indel_in_core {
            sbx_indel_core_quality(ctx, read, read_var_index)
        } else {
            let var_len = ctx.variant.ref_len();
            min_base_qual_across_range(read_var_index, read_var_index + var_len - 1, read)
        };
        if seq_tech == INVALID_BASE_QUAL {
            return QualityScores::invalid();
        }

        let mut recalibrated = seq_tech;
        if is_ultima() {
            if ctx.variant.is_snv() {
                let bqr_qual = rc.qual_cache.get_qual(
                    self, crate::ultima::max_raw_qual() as u8, read.consensus, 0, !read.negative_strand(),
                );
                recalibrated = jmin(seq_tech, bqr_qual);
            }
        } else if !rc.is_indel && ctx.artefact.is_none() {
            recalibrated = self.recalibrated_base_quality(rc, read_var_index, read);
            // don't boost the recalibrated qual of SBX artefacts
            if is_sbx()
                && ctx.has_indel_in_core
                && seq_tech < min_base_qual_across_range(read_var_index, read_var_index + ctx.variant.ref_len() - 1, read)
            {
                recalibrated = jmin(recalibrated, seq_tech);
            }
        }

        // calcMapQuality
        let penalised_map_quality = if self.is_highly_polymorphic(ctx.variant.pos) {
            HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY.min(read.mapq - self.cfg.fixed_map_qual_penalty)
        } else {
            let improper = if read.improper_pair() { self.cfg.improper_pair_penalty } else { 0 };
            let event_penalty = calc_event_penalty(number_of_events, read.len(), self.cfg.read_events_penalty);
            read.mapq - self.cfg.fixed_map_qual_penalty - improper - event_penalty
        };

        let mut penalised_base_quality = recalibrated;
        if !rc.use_msi_error_rate() {
            penalised_base_quality -= DEFAULT_BASE_QUAL_FIXED_PENALTY as f64;
        }

        // readEdgeDistancePenalty
        let n = read.len();
        let min_distance = if rc.is_indel {
            let lower = read_var_index - (ctx.var_index - rc.matcher.alt_index_lower);
            let upper = read_var_index + (rc.matcher.alt_index_upper - ctx.var_index);
            0.max(lower.min(n - 1 - upper))
        } else {
            let upper = read_var_index + ctx.variant.alt_len() - 1;
            0.max(read_var_index.min(n - 1 - upper))
        };
        let edge_penalty = if min_distance <= 0 {
            READ_EDGE_PENALTY_0
        } else if min_distance == 1 {
            READ_EDGE_PENALTY_1
        } else {
            0
        };
        penalised_base_quality -= edge_penalty as f64;

        let combined = jmax(0.0, jmin(penalised_map_quality as f64, penalised_base_quality));
        let is_medium = if uses_msi { !has_valid_msi_jitter_quals } else { is_medium_base_qual(seq_tech as i8) };
        QualityScores {
            seq_tech_base_quality: seq_tech,
            recalibrated_base_quality: recalibrated,
            final_map_quality: penalised_map_quality.max(0),
            final_base_quality: jmax(0.0, penalised_base_quality),
            combined_quality: combined,
            is_medium_qual: is_medium,
        }
    }

    fn recalibrated_base_quality(&self, rc: &ReadContextCounter, start_read_index: i32, read: &Read) -> f64 {
        let variant = &rc.ctx.variant;
        let pos_orient = !read.negative_strand();
        let mut consensus = read.consensus;
        if is_sbx() && consensus == ConsensusType::Dual && !in_duplex_region(pos_orient, read.duplex_index, start_read_index) {
            consensus = ConsensusType::Single;
        }
        if variant.is_snv() {
            return rc.qual_cache.get_qual(self, read.quals[start_read_index as usize], consensus, 0, pos_orient);
        }
        let max_index = (start_read_index + variant.ref_len()).min(read.quals.len() as i32) - 1;
        let max_length = max_index - start_read_index + 1;
        let mut quality = i32::MAX as f64;
        for i in 0..max_length {
            let raw = read.quals[(start_read_index + i) as usize];
            quality = jmin(quality, rc.qual_cache.get_qual(self, raw, consensus, i as usize, pos_orient));
        }
        quality
    }
}

fn min_base_qual_across_range(start: i32, end: i32, read: &Read) -> f64 {
    let mut min = read.quals[start as usize];
    for i in start + 1..=end {
        min = min.min(read.quals[i as usize]);
    }
    min as f64
}

pub(crate) fn average_core_quality(ctx: &VariantReadContext, read: &Read, read_var_index: i32) -> f64 {
    let s = (read_var_index - ctx.left_core_length()).max(0);
    let e = (read_var_index + ctx.right_core_length()).min(read.len() - 1);
    let len = e - s + 1;
    if len <= 0 {
        return 0.0;
    }
    let total: f64 = (s..=e).map(|i| read.quals[i as usize] as f64).sum();
    jround(total / len as f64) as i32 as f64
}

/// FilterUtils.isChimericRead
pub fn is_chimeric_read(read: &Read) -> bool {
    if read.paired() {
        if read.flag(FLAG_MATE_UNMAPPED) {
            return true;
        }
        if read.tid != read.mate_tid {
            return true;
        }
        if read.negative_strand() == read.flag(FLAG_MATE_REVERSE) {
            return true;
        }
    }
    read.insert_size.abs() > 1000
}

pub fn is_high_base_qual(q: f64) -> bool {
    is_high_base_qual_byte(q as i8)
}

/// SbxBamUtils.inDuplexRegion
fn in_duplex_region(pos_orientation_read: bool, duplex_base_index: i32, base_index: i32) -> bool {
    if duplex_base_index < 0 {
        return false;
    }
    if pos_orientation_read {
        base_index >= duplex_base_index
    } else {
        base_index <= duplex_base_index
    }
}

/// SbxUtils.determineConsensusType: a duplex read only counts as DUAL if the whole core is in its duplex region
fn sbx_consensus_type(ctx: &VariantReadContext, read_var_index: i32, read: &Read) -> ConsensusType {
    if read.consensus != ConsensusType::Dual {
        return read.consensus;
    }
    let start = read_var_index - ctx.left_core_length();
    let end = read_var_index + ctx.right_core_length();
    for i in start..=end {
        if !in_duplex_region(!read.negative_strand(), read.duplex_index, i) {
            return ConsensusType::Single;
        }
    }
    ConsensusType::Dual
}

/// SbxUtils.indelCoreQuality: the minimum base quality across the core
fn sbx_indel_core_quality(ctx: &VariantReadContext, read: &Read, read_var_index: i32) -> f64 {
    let start = (read_var_index - ctx.left_core_length()).max(0);
    let end = (read_var_index + ctx.right_core_length()).min(read.len() - 1);
    if end - start + 1 <= 0 {
        return 0.0;
    }
    min_base_qual_across_range(start, end, read)
}

/// Unused config placeholder so the counter's constructor signature stays stable.
pub fn max_coverage(cfg: &Config, tier: Tier) -> i32 {
    if matches!(tier, Tier::Hotspot | Tier::Panel) {
        cfg.max_read_depth_panel
    } else {
        cfg.max_read_depth
    }
}

/// Reference bases are only needed by the quality calculator; re-exported for callers.
pub type RefSeq = RefSequence;
