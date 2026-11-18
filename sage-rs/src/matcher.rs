//! ReadContextMatcher: classifies a read against a variant's read context.

use crate::constants::*;
use crate::read_context::VariantReadContext;
use crate::types::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchType {
    None,
    Core,
    PartialCore,
    Full,
    Realigned,
    Ref,
    SimpleAlt,
}

impl MatchType {
    pub fn supports_alt(self) -> bool {
        matches!(self, MatchType::Core | MatchType::PartialCore | MatchType::Full | MatchType::Realigned)
    }
    pub fn full_alt_support(self) -> bool {
        matches!(self, MatchType::PartialCore | MatchType::Full | MatchType::Realigned)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MatchInfo {
    pub kind: MatchType,
    pub exact: bool,
}

pub const NO_MATCH: MatchInfo = MatchInfo { kind: MatchType::None, exact: false };

#[derive(Clone, Copy, PartialEq, Eq)]
enum BaseMatch {
    Match,
    LowQualMismatches,
    Incomplete,
    Mismatch,
}

fn valid(m: BaseMatch) -> bool {
    matches!(m, BaseMatch::Match | BaseMatch::LowQualMismatches)
}

#[derive(Clone, Debug)]
struct LowQualExclusion {
    indices: Vec<i32>,
    is_range: bool,
    offset: i32,
}

impl LowQualExclusion {
    fn covers_with(&self, index: i32, offset: i32) -> bool {
        let index = index + offset;
        if self.is_range {
            index >= self.indices[0] && index <= self.indices[1]
        } else {
            self.indices.contains(&index)
        }
    }
    fn covers(&self, index: i32) -> bool {
        let mut c = self.covers_with(index, self.offset);
        if self.offset != 0 {
            c |= self.covers_with(index, 0);
        }
        c
    }
}

#[derive(Clone, Debug)]
pub struct Matcher {
    max_core_low_qual: i32,
    pub is_reference: bool,
    pub alt_index_lower: i32,
    pub alt_index_upper: i32,
    excl_read: Option<LowQualExclusion>,
    excl_ref: Option<LowQualExclusion>,
}

fn calc_max_low_qual_mismatches(len: i32) -> i32 {
    1 + len / CORE_LOW_QUAL_MISMATCH_FACTOR
}

fn index_of(hay: &[u8], needle: &[u8], from: i32) -> i32 {
    // String.indexOf(str, fromIndex): negative from is treated as 0
    let from = from.max(0) as usize;
    if needle.is_empty() {
        return if from <= hay.len() { from as i32 } else { hay.len() as i32 };
    }
    if hay.len() < needle.len() {
        return -1;
    }
    let mut i = from;
    while i + needle.len() <= hay.len() {
        if &hay[i..i + needle.len()] == needle {
            return i as i32;
        }
        i += 1;
    }
    -1
}

pub fn java_index_of(hay: &[u8], needle: &[u8], from: i32) -> i32 {
    index_of(hay, needle, from)
}

impl Matcher {
    pub fn new(ctx: &VariantReadContext, allow_mismatches: bool, is_reference: bool) -> Matcher {
        let max_core_low_qual = if allow_mismatches {
            let mut core_length = ctx.core_length();
            if let Some(mr) = &ctx.max_repeat {
                let core = ctx.core();
                let mut count = 0;
                let mut idx = index_of(core, &mr.bases, mr.index - ctx.left_flank_length());
                while idx >= 0 {
                    count += 1;
                    idx = index_of(core, &mr.bases, idx + 1);
                }
                if count > 2 {
                    core_length -= (count - 2) * mr.repeat_length();
                }
            }
            calc_max_low_qual_mismatches(core_length)
        } else {
            0
        };

        let v = &ctx.variant;
        let hom_len = ctx.homology.as_ref().map_or(0, |h| h.length);
        let mut lower = ctx.var_index;
        let mut upper = if v.is_insert() {
            ctx.var_index + hom_len + 1 + v.indel_length().abs()
        } else if v.is_delete() {
            ctx.var_index + hom_len + 1
        } else {
            ctx.var_index + v.alt_len() - 1
        };
        if is_sbx() && !v.is_indel() && ctx.has_indel_in_core {
            lower = ctx.core_index_start;
            upper = ctx.core_index_end;
        } else if !is_reference && !ctx.all_repeats.is_empty() {
            let mut min_r = ctx.all_repeats[0].index - 1;
            let mut max_r = ctx.all_repeats[0].end_index() + 1;
            for r in &ctx.all_repeats[1..] {
                min_r = min_r.min(r.index - 1);
                max_r = max_r.max(r.end_index() + 1);
            }
            lower = lower.min(min_r);
            upper = upper.max(max_r);
        }

        let (excl_read, excl_ref) = if allow_mismatches {
            if v.is_indel() {
                let low = if v.is_insert() { ctx.var_index + v.indel_length() } else { ctx.var_index };
                // determineIndelLowQualRefReadDiffIndex
                let mut ref_index = ctx.variant_ref_index();
                let mut read_index = ctx.var_index;
                while read_index <= ctx.core_index_end && ref_index < ctx.ref_bases.len() as i32 {
                    if ctx.ref_bases[ref_index as usize] != ctx.read_bases[read_index as usize] {
                        break;
                    }
                    read_index += 1;
                    ref_index += 1;
                }
                let high = read_index.min(ctx.core_index_end);
                // Java collects a HashSet<Integer> into a list; only membership is used
                let mut set = vec![low, high];
                if ctx.homology.is_some() {
                    set.extend([ctx.var_index, ctx.var_index + 1, upper - 1, upper]);
                }
                let diff = ctx.left_flank_length();
                (
                    Some(LowQualExclusion { indices: set, is_range: false, offset: 0 }),
                    Some(LowQualExclusion { indices: vec![low - diff, high - diff], is_range: false, offset: 0 }),
                )
            } else {
                let alt_range = v.alt_len() - 1;
                let ref_index = ctx.left_core_length();
                (
                    Some(LowQualExclusion { indices: vec![ctx.var_index, ctx.var_index + alt_range], is_range: true, offset: 0 }),
                    Some(LowQualExclusion { indices: vec![ref_index, ref_index + alt_range], is_range: true, offset: 0 }),
                )
            }
        } else {
            (None, None)
        };

        Matcher { max_core_low_qual, is_reference, alt_index_lower: lower, alt_index_upper: upper, excl_read, excl_ref }
    }

    /// Stands in for the null matcher of a counter whose read context is marked invalid; never used.
    pub fn placeholder(is_reference: bool) -> Matcher {
        Matcher { max_core_low_qual: 0, is_reference, alt_index_lower: 0, alt_index_upper: 0, excl_read: None, excl_ref: None }
    }

    pub fn set_realignment_index_offset(&mut self, offset: i32) {
        if let Some(e) = self.excl_read.as_mut() {
            e.offset = offset;
        }
    }
    pub fn clear_realignment_index_offset(&mut self) {
        self.set_realignment_index_offset(0);
    }
    fn read_excl_offset(&self) -> i32 {
        self.excl_read.as_ref().map_or(0, |e| e.offset)
    }

    pub fn covers_variant(&self, ctx: &VariantReadContext, read_len: i32, read_var_index: i32) -> bool {
        if read_var_index < 0 {
            return false;
        }
        let lower = read_var_index - ctx.var_index + self.alt_index_lower;
        let upper = read_var_index - ctx.var_index + self.alt_index_upper;
        !(lower < 0 || upper >= read_len)
    }

    pub fn has_uncertain_core_bases(&self, ctx: &VariantReadContext, quals: &[u8], read_var_index: i32) -> bool {
        let lower = read_var_index - ctx.var_index + self.alt_index_lower;
        let upper = read_var_index - ctx.var_index + self.alt_index_upper;
        let mut i = lower.max(0);
        while i <= upper.min(quals.len() as i32 - 1) {
            if is_uncertain_base_qual(quals[i as usize]) {
                return true;
            }
            i += 1;
        }
        false
    }

    /// determineReadMatchInfo(SAMRecord, readVarIndex)
    pub fn determine_read_match(&self, ctx: &VariantReadContext, read: &Read, read_var_index: i32) -> MatchInfo {
        let info = self.determine(ctx, &read.bases, Some(&read.quals), read_var_index, false);
        if info.kind == MatchType::None && self.is_reference && is_simple_alt_match(&ctx.variant, read, read_var_index) {
            return MatchInfo { kind: MatchType::SimpleAlt, exact: true };
        }
        info
    }

    /// determineReadMatchInfo(readBases, readQuals, readVarIndex, skipRefMatch)
    pub fn determine(
        &self, ctx: &VariantReadContext, read_bases: &[u8], read_quals: Option<&[u8]>, read_var_index: i32,
        skip_ref_match: bool,
    ) -> MatchInfo {
        if !skip_ref_match {
            let m = self.core_matches_ref(ctx, read_bases, read_quals, read_var_index);
            if m.kind != MatchType::None {
                return m;
            }
        }
        let core = self.determine_core_match(ctx, read_bases, read_quals, read_var_index);
        if core.kind == MatchType::None {
            if self.is_reference && ctx.variant.is_long_insert() {
                if let Some(ext) = &ctx.extended_ref_bases {
                    let ext_start = ctx.variant_ref_index() + ctx.left_flank_length() + self.read_excl_offset();
                    let (ref_m, read_m) =
                        count_ref_and_read_differences(read_bases, &ctx.read_bases, ext, read_var_index, ctx.var_index, ext_start);
                    if read_m - ref_m >= LONG_GERMLINE_INSERT_READ_VS_REF_DIFF {
                        return MatchInfo { kind: MatchType::PartialCore, exact: true };
                    }
                }
            }
            return NO_MATCH;
        }
        let left = self.determine_flank_match(ctx, read_bases, read_quals, read_var_index, true);
        let right = self.determine_flank_match(ctx, read_bases, read_quals, read_var_index, false);
        if right == BaseMatch::Mismatch || left == BaseMatch::Mismatch {
            return MatchInfo { kind: MatchType::Core, exact: core.exact };
        }
        if !valid(left) && !valid(right) {
            return MatchInfo { kind: MatchType::Core, exact: core.exact };
        }
        let exact = core.exact && left != BaseMatch::LowQualMismatches && right != BaseMatch::LowQualMismatches;
        if core.kind == MatchType::PartialCore {
            MatchInfo { kind: MatchType::PartialCore, exact }
        } else {
            MatchInfo { kind: MatchType::Full, exact }
        }
    }

    fn core_matches_ref(&self, ctx: &VariantReadContext, read_bases: &[u8], quals: Option<&[u8]>, read_var_index: i32) -> MatchInfo {
        let mut ref_index_start = 0;
        let ref_index = ctx.variant_ref_index();
        let mut read_index_start = read_var_index - ref_index;
        let mut read_index_end = read_index_start + ctx.ref_bases.len() as i32 - 1;
        if read_index_start < 0 {
            ref_index_start += read_index_start.abs();
            read_index_start = 0;
        }
        read_index_end = read_index_end.min(read_bases.len() as i32 - 1);
        let compare_length = read_index_end - read_index_start + 1;
        let m = match_type(
            &ctx.ref_bases, read_bases, quals, ref_index_start, read_index_start, compare_length, self.max_core_low_qual,
            self.excl_ref.as_ref(),
        );
        if valid(m) {
            MatchInfo { kind: MatchType::Ref, exact: m == BaseMatch::Match }
        } else {
            NO_MATCH
        }
    }

    fn determine_core_match(&self, ctx: &VariantReadContext, read_bases: &[u8], quals: Option<&[u8]>, read_var_index: i32) -> MatchInfo {
        let mut read_index_start = read_var_index - ctx.left_core_length();
        let read_index_end = read_var_index + ctx.right_core_length();
        let n = read_bases.len() as i32;
        if read_index_start >= 0 && read_index_end < n {
            let m = match_type(
                &ctx.read_bases, read_bases, quals, ctx.core_index_start, read_index_start, ctx.core_length(),
                self.max_core_low_qual, self.excl_read.as_ref(),
            );
            if valid(m) {
                MatchInfo { kind: MatchType::Core, exact: m == BaseMatch::Match }
            } else {
                NO_MATCH
            }
        } else {
            let left_trim = if read_index_start < 0 { read_index_start.abs() } else { 0 };
            let right_trim = if read_index_end >= n { (n - read_index_end - 1).abs() } else { 0 };
            read_index_start += left_trim;
            let core_index_start = ctx.core_index_start + left_trim;
            let compare_length = ctx.core_length() - left_trim - right_trim;
            let max_lq = self.max_core_low_qual.min(calc_max_low_qual_mismatches(compare_length));
            let m = match_type(
                &ctx.read_bases, read_bases, quals, core_index_start, read_index_start, compare_length, max_lq,
                self.excl_read.as_ref(),
            );
            if valid(m) {
                MatchInfo { kind: MatchType::PartialCore, exact: m == BaseMatch::Match }
            } else {
                NO_MATCH
            }
        }
    }

    fn determine_flank_match(
        &self, ctx: &VariantReadContext, read_bases: &[u8], quals: Option<&[u8]>, read_var_index: i32, is_left: bool,
    ) -> BaseMatch {
        let (mut read_index_start, read_index_end, mut flank_start_index, mut flank_length);
        if is_left {
            read_index_start = read_var_index - ctx.left_length();
            read_index_end = read_var_index - ctx.left_core_length() - 1;
            flank_start_index = 0;
            flank_length = ctx.left_flank_length();
        } else {
            read_index_start = read_var_index + ctx.right_core_length() + 1;
            read_index_end = read_var_index + ctx.right_length() - 1;
            flank_start_index = ctx.core_index_end + 1;
            flank_length = ctx.right_flank_length();
        }
        let n = read_bases.len() as i32;
        let mut is_partial = false;
        if read_index_start < 0 || read_index_end >= n {
            is_partial = true;
            let left_trim = if read_index_start < 0 { read_index_start.abs() } else { 0 };
            let right_trim = if read_index_end >= n { (n - read_index_end - 1).abs() } else { 0 };
            flank_length -= left_trim + right_trim;
            read_index_start += left_trim;
            flank_start_index += left_trim;
        }
        if flank_length <= 0 {
            return BaseMatch::Incomplete;
        }
        let m = match_type(
            &ctx.read_bases, read_bases, quals, flank_start_index, read_index_start, flank_length,
            FLANK_LOW_QUAL_MISMATCHES, self.excl_read.as_ref(),
        );
        if !valid(m) {
            return BaseMatch::Mismatch;
        }
        if is_partial {
            BaseMatch::Incomplete
        } else {
            m
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn match_type(
    bases: &[u8], read_bases: &[u8], quals: Option<&[u8]>, base_index_start: i32, read_index_start: i32,
    compare_length: i32, max_low_qual_mismatches: i32, excl: Option<&LowQualExclusion>,
) -> BaseMatch {
    if compare_length <= 0 {
        return BaseMatch::Mismatch;
    }
    let base_index_end = base_index_start + compare_length - 1;
    let read_index_end = read_index_start + compare_length - 1;
    if base_index_end >= bases.len() as i32 || read_index_end >= read_bases.len() as i32 {
        return BaseMatch::Incomplete;
    }
    // standard base quals cannot be used for Ultima
    let quals = if is_ultima() { None } else { quals };
    let mut mismatches = 0;
    let (mut i, mut j) = (base_index_start, read_index_start);
    while i <= base_index_end && j <= read_index_end {
        if bases[i as usize] != read_bases[j as usize] {
            let Some(q) = quals else { return BaseMatch::Mismatch };
            if excl.map_or(false, |e| e.covers(i)) {
                return BaseMatch::Mismatch;
            }
            if q[j as usize] < matching_base_quality() {
                mismatches += 1;
                if mismatches > max_low_qual_mismatches {
                    return BaseMatch::Mismatch;
                }
            } else {
                return BaseMatch::Mismatch;
            }
        }
        i += 1;
        j += 1;
    }
    if mismatches == 0 {
        BaseMatch::Match
    } else {
        BaseMatch::LowQualMismatches
    }
}

/// ReadContextMatcher.isSimpleAltMatch
pub fn is_simple_alt_match(variant: &SimpleVariant, read: &Read, read_var_index: i32) -> bool {
    if variant.is_indel() {
        let mut read_index = 0;
        for el in &read.cigar {
            if read_index == read_var_index + 1 {
                if el.len != variant.indel_length_abs() {
                    return false;
                }
                if variant.is_delete() {
                    return el.op == Op::D;
                }
                if el.op != Op::I {
                    return false;
                }
                for i in 0..el.len {
                    if read.bases[(read_index + i) as usize] != variant.alt[(i + 1) as usize] {
                        return false;
                    }
                }
                return true;
            }
            if el.op.consumes_read() {
                read_index += el.len;
            }
        }
        false
    } else {
        for i in 0..variant.alt_len() {
            if read.bases[(read_var_index + i) as usize] != variant.alt[i as usize] {
                return false;
            }
        }
        true
    }
}

fn count_ref_and_read_differences(
    read_bases: &[u8], ctx_bases: &[u8], ref_bases: &[u8], read_index_start: i32, ctx_index_start: i32, ref_index_start: i32,
) -> (i32, i32) {
    let (mut ctx_matches, mut ref_matches) = (0, 0);
    let mut i = (-ref_index_start).max(0);
    loop {
        let rb = read_index_start + i;
        if rb >= read_bases.len() as i32 {
            break;
        }
        let ci = ctx_index_start + i;
        if ci >= ctx_bases.len() as i32 {
            break;
        }
        let ri = ref_index_start + i;
        if ri >= ref_bases.len() as i32 {
            break;
        }
        if read_bases[rb as usize] == ctx_bases[ci as usize] {
            ctx_matches += 1;
        }
        if read_bases[rb as usize] == ref_bases[ri as usize] {
            ref_matches += 1;
        }
        i += 1;
    }
    (ref_matches, ctx_matches)
}
