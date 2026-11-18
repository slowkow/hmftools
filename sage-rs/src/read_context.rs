//! VariantReadContext, its builder, ReadCigarInfo, Microhomology and IlluminaArtefactContext.

use crate::constants::*;
use crate::refseq::RefSequence;
use crate::repeat::*;
use crate::types::*;

#[derive(Clone, Debug)]
pub struct Microhomology {
    pub bases: Vec<u8>,
    pub length: i32,
}

/// Microhomology.findHomology
pub fn find_homology(variant: &SimpleVariant, read_bases: &[u8], var_read_index: i32, apply_indel_length: bool) -> Option<Microhomology> {
    if !variant.is_indel() {
        return None;
    }
    let indel_alt_length = variant.indel_length_abs();
    let indel_bases: &[u8] = if variant.is_insert() { &variant.alt[1..] } else { &variant.ref_[1..] };
    let mut hom_read_index = var_read_index + 1;
    if apply_indel_length {
        hom_read_index += indel_alt_length;
    }
    let n = read_bases.len() as i32;
    let mut homology: Vec<u8> = Vec::new();
    let mut i = 0;
    while i < indel_bases.len() && hom_read_index < n {
        // a negative index throws in Java; the builder treats that as "no context"
        let indel_base = indel_bases[i];
        if indel_base != read_bases[hom_read_index as usize] {
            break;
        }
        homology.push(indel_base);
        i += 1;
        hom_read_index += 1;
    }
    if homology.is_empty() {
        return None;
    }
    let mut homology_length = homology.len() as i32;
    if homology_length == indel_alt_length && hom_read_index < n {
        let mut matched = true;
        while matched {
            let mut match_count = 0;
            let mut i = 0;
            while i < indel_alt_length && hom_read_index < n {
                if homology[i as usize] != read_bases[hom_read_index as usize] {
                    matched = false;
                    break;
                }
                match_count += 1;
                i += 1;
                hom_read_index += 1;
            }
            if match_count > 0 {
                homology_length += match_count;
            } else {
                break;
            }
        }
    }
    Some(Microhomology { bases: homology, length: homology_length })
}

/// Microhomology.findLeftHomologyShift
pub fn find_left_homology_shift(variant: &SimpleVariant, ref_seq: &RefSequence, read_bases: &[u8], var_read_index: i32) -> i32 {
    if !variant.is_indel() {
        return 0;
    }
    let mut offset = 0;
    let mut pos = variant.pos;
    let mut last_ref = *variant.ref_.last().unwrap();
    let mut alt = variant.alt.clone();
    let alt_len = alt.len();
    let mut last_alt = alt[alt_len - 1];
    while last_ref == last_alt {
        pos -= 1;
        offset += 1;
        if pos < ref_seq.start {
            break;
        }
        last_ref = ref_seq.base(pos);
        if var_read_index - offset < 0 {
            break;
        }
        let prev = read_bases[(var_read_index - offset) as usize];
        alt.pop();
        alt.insert(0, prev);
        last_alt = alt[alt_len - 1];
    }
    offset
}

#[derive(Clone, Debug)]
pub struct ArtefactContext {
    pub offset: [i32; 2],
    pub base: [i8; 2],
}

pub const ARTEFACT_NO_INDEX: i32 = -10;
const HOMOPOLYMER_BASE_SEARCH: i32 = 3;
const HOMOPOLYMER_REPEAT_LENGTH: i32 = 8;

fn has_homopolymer_repeat(bases: &[u8], start: i32, down: bool) -> Option<bool> {
    let n = bases.len() as i32;
    if start < 0 || start >= n {
        return None; // StringIndexOutOfBounds in Java
    }
    let mut count = 1;
    let last = bases[start as usize];
    let mut i = start + if down { -1 } else { 1 };
    while i >= 0 && i < n {
        if bases[i as usize] == last {
            count += 1;
            if count >= HOMOPOLYMER_REPEAT_LENGTH {
                return Some(true);
            }
        } else {
            return Some(false);
        }
        i += if down { -1 } else { 1 };
    }
    Some(false)
}

impl ArtefactContext {
    /// IlluminaArtefactContext.buildContext. Err(()) stands for the exception Java throws on an
    /// out-of-range charAt, which the caller turns into "no read context".
    pub fn build(rc: &VariantReadContext) -> Result<Option<ArtefactContext>, ()> {
        let bases = &rc.read_bases;
        let n = bases.len() as i32;
        let at = |i: i32| -> Result<u8, ()> {
            if i < 0 || i >= n {
                Err(())
            } else {
                Ok(bases[i as usize])
            }
        };
        let index_in_bases = rc.var_index;
        let variant = &rc.variant;
        let mut hp_start_offset = ARTEFACT_NO_INDEX;
        let mut hp_end_offset = ARTEFACT_NO_INDEX;
        let mut hp_start_base: i8 = ARTEFACT_NO_INDEX as i8;
        let mut hp_end_base: i8 = ARTEFACT_NO_INDEX as i8;
        let alt_length = variant.alt_len();
        for var_index in index_in_bases..index_in_bases + alt_length {
            if has_homopolymer_repeat(bases, var_index, true).ok_or(())? {
                hp_start_base = at(var_index)? as i8;
                hp_start_offset = index_in_bases - var_index;
                if at(var_index + 1)? == at(var_index)? {
                    hp_start_offset -= 1;
                }
            }
            if hp_end_offset != ARTEFACT_NO_INDEX {
                break;
            }
            if has_homopolymer_repeat(bases, var_index, false).ok_or(())? {
                hp_end_base = at(var_index)? as i8;
                hp_end_offset = var_index - index_in_bases;
                if at(var_index - 1)? == at(var_index)? {
                    hp_end_offset -= 1;
                }
            } else if variant.is_delete() {
                let deleted = &variant.ref_[1..];
                let del_len = deleted.len() as i32;
                let mut del_start = index_in_bases + 1;
                while del_start + del_len < n && deleted == &bases[del_start as usize..(del_start + del_len) as usize] {
                    del_start += del_len;
                }
                if has_homopolymer_repeat(bases, del_start, false).ok_or(())? {
                    hp_end_base = at(del_start)? as i8;
                    hp_end_offset = del_start - index_in_bases;
                }
            }
        }
        if hp_start_offset == ARTEFACT_NO_INDEX && hp_end_offset == ARTEFACT_NO_INDEX {
            return Ok(None);
        }
        Ok(Some(ArtefactContext { offset: [hp_start_offset, hp_end_offset], base: [hp_start_base, hp_end_base] }))
    }

    /// findApplicableBaseQual; -1 (INVALID_BASE_QUAL) when not applicable.
    pub fn find_applicable_base_qual(&self, read: &Read, var_read_index: i32) -> i32 {
        let side = if read.negative_strand() { SE_END } else { SE_START };
        if self.offset[side] == ARTEFACT_NO_INDEX {
            return -1;
        }
        let search_down = side == SE_START;
        let hp_offset = self.offset[side];
        let hp_base = self.base[side];
        let hp_start = if search_down { var_read_index - hp_offset } else { var_read_index + hp_offset };
        let mut min_qual: i32 = -1;
        for i in 0..HOMOPOLYMER_BASE_SEARCH {
            let hp_index = hp_start + if search_down { -i } else { i };
            if hp_index < 0 || hp_index >= read.quals.len() as i32 {
                break;
            }
            if read.bases[hp_index as usize] as i8 != hp_base {
                continue;
            }
            let q = read.quals[hp_index as usize] as i32;
            min_qual = if min_qual < 0 { q } else { min_qual.min(q) };
        }
        min_qual as i8 as i32
    }
}

#[derive(Clone, Debug)]
pub struct VariantReadContext {
    pub variant: SimpleVariant,
    pub alignment_start: i32,
    pub alignment_end: i32,
    pub ref_bases: Vec<u8>,
    pub read_bases: Vec<u8>,
    pub core_index_start: i32,
    pub var_index: i32,
    pub core_index_end: i32,
    pub homology: Option<Microhomology>,
    pub max_repeat: Option<RepeatInfo>,
    pub all_repeats: Vec<RepeatInfo>,
    pub core_position_start: i32,
    pub core_position_end: i32,
    pub read_cigar: String,
    pub has_indel_in_core: bool,
    pub artefact: Option<ArtefactContext>,
    pub max_ref_repeat: Option<RepeatInfo>,
    pub extended_ref_bases: Option<Vec<u8>>,
    pub marked_invalid: bool,
}

impl VariantReadContext {
    pub fn core_length(&self) -> i32 {
        self.core_index_end - self.core_index_start + 1
    }
    pub fn left_core_length(&self) -> i32 {
        self.var_index - self.core_index_start
    }
    pub fn right_core_length(&self) -> i32 {
        self.core_index_end - self.var_index
    }
    pub fn left_flank_length(&self) -> i32 {
        self.core_index_start
    }
    pub fn right_flank_length(&self) -> i32 {
        self.read_bases.len() as i32 - self.core_index_end - 1
    }
    pub fn left_length(&self) -> i32 {
        self.var_index
    }
    pub fn right_length(&self) -> i32 {
        self.read_bases.len() as i32 - self.var_index
    }
    pub fn total_length(&self) -> i32 {
        self.read_bases.len() as i32
    }
    pub fn variant_ref_index(&self) -> i32 {
        self.variant.pos - self.core_position_start
    }
    pub fn core(&self) -> &[u8] {
        &self.read_bases[self.core_index_start as usize..=self.core_index_end as usize]
    }
    pub fn core_str(&self) -> &str {
        std::str::from_utf8(self.core()).unwrap()
    }
    pub fn left_flank_str(&self) -> &str {
        std::str::from_utf8(&self.read_bases[..self.core_index_start as usize]).unwrap()
    }
    pub fn right_flank_str(&self) -> &str {
        std::str::from_utf8(&self.read_bases[(self.core_index_end + 1) as usize..]).unwrap()
    }
    pub fn max_repeat_count(&self) -> i32 {
        self.max_repeat.as_ref().map_or(0, |r| r.count)
    }
    pub fn trinucleotide(&self) -> &[u8] {
        let i = self.variant_ref_index() as usize;
        &self.ref_bases[i - 1..=i + 1]
    }
    pub fn is_valid(&self) -> bool {
        if self.marked_invalid {
            return false;
        }
        if self.core_index_start <= 0 || self.core_index_end >= self.read_bases.len() as i32 - 1 {
            return false;
        }
        let min_core = if self.variant.is_indel() { MIN_CORE_DISTANCE * 2 } else { MIN_CORE_DISTANCE * 2 + 1 };
        if self.core_length() < min_core {
            return false;
        }
        if self.var_index <= self.core_index_start || self.core_index_end <= self.var_index {
            return false;
        }
        if self.variant.pos <= self.core_position_start || self.variant.pos >= self.core_position_end {
            return false;
        }
        true
    }
}

pub struct ReadCigarInfo {
    pub cigar: Vec<CigarEl>,
    pub read_alignment_start: i32,
    pub flank_position_start: i32,
    pub flank_position_end: i32,
    pub core_position_start: i32,
    pub core_position_end: i32,
    pub flank_index_start: i32,
    pub flank_index_end: i32,
}

impl ReadCigarInfo {
    pub fn is_valid(&self) -> bool {
        let read_bases: i32 = self.cigar.iter().filter(|x| x.op.consumes_read()).map(|x| x.len).sum();
        self.flank_index_start >= 0
            && self.flank_index_end > self.flank_index_start
            && self.core_position_start < self.core_position_end
            && self.core_position_start >= self.flank_position_start
            && self.core_position_end <= self.flank_position_end
            && read_bases == self.flank_index_end - self.flank_index_start + 1
    }

    /// ReadCigarInfo.buildReadCigar
    pub fn build(
        read_start: i32, elements: &[CigarEl], left_flank_index: i32, left_core_index: i32, right_core_index: i32,
        right_flank_index: i32,
    ) -> ReadCigarInfo {
        let mut cigar: Vec<CigarEl> = Vec::new();
        let mut read_index = 0;
        let mut ref_position = read_start;
        let (mut flank_pos_start, mut flank_pos_end, mut core_pos_start, mut core_pos_end) = (0, 0, 0, 0);
        let mut final_index_start = left_flank_index;
        let mut final_index_end = right_flank_index;

        for &element in elements {
            let op = element.op;
            if read_index == 0 && op == Op::S {
                ref_position -= element.len;
            }
            let element_end_index = read_index + element.len - 1;
            let is_read_bases = op.consumes_read();

            if is_read_bases && read_index <= left_flank_index && element_end_index >= right_flank_index {
                let el_start = read_index.max(left_flank_index);
                let el_end = element_end_index.min(right_flank_index);
                cigar.push(CigarEl::new(el_end - el_start + 1, op));
                flank_pos_start = ref_position + (left_flank_index - read_index);
                flank_pos_end = ref_position + (right_flank_index - read_index);
                core_pos_start = ref_position + (left_core_index - read_index);
                core_pos_end = ref_position + (right_core_index - read_index);
                break;
            }

            let mut added = false;
            if element_end_index >= left_flank_index {
                let read_end_index = if is_read_bases { element_end_index } else { read_index };
                if core_pos_start == 0 && position_within(left_core_index, read_index, read_end_index) {
                    core_pos_start = match op {
                        Op::I => ref_position - 1,
                        Op::D | Op::N => ref_position + element.len,
                        _ => ref_position + (left_core_index - read_index),
                    };
                }
                if core_pos_end == 0 && position_within(right_core_index, read_index, read_end_index) {
                    core_pos_end = match op {
                        Op::I => ref_position,
                        Op::D | Op::N => ref_position + element.len,
                        _ => ref_position + (right_core_index - read_index),
                    };
                }
                if is_read_bases {
                    if flank_pos_start == 0 {
                        let el_start = read_index.max(left_flank_index);
                        let el_len =
                            if op.is_indel() || op == Op::N { element.len } else { element_end_index - el_start + 1 };
                        cigar.push(CigarEl::new(el_len, op));
                        if op.is_indel() {
                            if cigar.len() > 1 {
                                cigar.remove(0);
                            }
                            let extra_index_start = left_flank_index - read_index + 1;
                            cigar.insert(0, CigarEl::new(1, Op::M));
                            final_index_start -= extra_index_start;
                            flank_pos_start = (ref_position - 1).max(read_start);
                        } else {
                            flank_pos_start = ref_position + (left_flank_index - read_index);
                        }
                        added = true;
                    } else if flank_pos_end == 0 && element_end_index >= right_flank_index {
                        if op.is_indel() || op == Op::N {
                            cigar.push(element);
                            cigar.push(CigarEl::new(1, Op::M));
                            if op == Op::I {
                                flank_pos_end = ref_position;
                                final_index_end += element_end_index + 1 - right_flank_index;
                            } else {
                                flank_pos_end = ref_position + element.len;
                            }
                        } else {
                            let el_end = element_end_index.min(right_flank_index);
                            cigar.push(CigarEl::new(el_end - read_index + 1, op));
                            flank_pos_end = ref_position + (right_flank_index - read_index);
                        }
                        break;
                    }
                }
            }
            if !added {
                let past_flank_start = read_index - if op == Op::D { 1 } else { 0 } >= left_flank_index;
                if past_flank_start {
                    cigar.push(element);
                }
            }
            if op.consumes_read() {
                read_index += element.len;
            }
            if op.consumes_ref() || op == Op::S {
                ref_position += element.len;
            }
            if read_index > right_flank_index {
                break;
            }
        }
        ReadCigarInfo {
            cigar,
            read_alignment_start: read_start,
            flank_position_start: flank_pos_start,
            flank_position_end: flank_pos_end,
            core_position_start: core_pos_start,
            core_position_end: core_pos_end,
            flank_index_start: final_index_start,
            flank_index_end: final_index_end,
        }
    }
}

/// ReadCigarInfo.hasIndelInCore (the read index advances on every element, as in Java)
pub fn has_indel_in_core(cigar: &[CigarEl], core_index_start: i32, core_index_end: i32) -> bool {
    let mut read_index = 0;
    for el in cigar {
        if read_index > core_index_end {
            break;
        }
        if el.op.consumes_read() {
            if el.op == Op::I && read_index + el.len - 1 >= core_index_start {
                return true;
            }
        } else if el.op == Op::D && read_index >= core_index_start {
            return true;
        }
        read_index += el.len;
    }
    false
}

struct SoftClipAdjustment {
    alignment_start: i32,
    cigar: Vec<CigarEl>,
}

fn check_indel_soft_clip_adjustment(read: &Read, variant: &SimpleVariant, var_read_index: i32) -> Option<SoftClipAdjustment> {
    if !variant.is_long_insert() {
        return None;
    }
    let orig = &read.cigar;
    if orig[0].op == Op::S && var_read_index < orig[0].len {
        let left_sc = orig[0].len;
        let mut converted = orig.clone();
        converted.remove(0);
        converted.insert(0, CigarEl::new(var_read_index + 1, Op::M));
        converted.insert(1, CigarEl::new(variant.indel_length(), Op::I));
        let remaining = left_sc - (var_read_index + 1) - variant.indel_length();
        if remaining > 0 {
            let post = converted[2].len;
            converted[2] = CigarEl::new(remaining + post, Op::M);
        }
        return Some(SoftClipAdjustment { alignment_start: variant.pos - var_read_index, cigar: converted });
    }
    let last = orig.len() - 1;
    if orig[last].op == Op::S {
        let right_sc = orig[last].len;
        let right_sc_start = read.len() - right_sc;
        if var_read_index >= right_sc_start {
            let mut converted = orig.clone();
            converted.remove(last);
            converted.push(CigarEl::new(variant.indel_length(), Op::I));
            converted.push(CigarEl::new(right_sc - variant.indel_length(), Op::M));
            return Some(SoftClipAdjustment { alignment_start: read.start, cigar: converted });
        }
    }
    None
}

/// VariantReadContextBuilder.buildContext
pub fn build_context(
    variant: &SimpleVariant, read: &Read, var_index_in_read: i32, ref_seq: &RefSequence, flank_size: i32,
) -> Option<VariantReadContext> {
    let read_bases = &read.bases;
    let n = read_bases.len() as i32;
    let mut read_core_start;
    let mut read_core_end;
    let mut homology: Option<Microhomology> = None;
    let mut sc_adjust: Option<SoftClipAdjustment> = None;

    if variant.is_indel() {
        read_core_start = var_index_in_read - MIN_CORE_DISTANCE + 1;
        read_core_end =
            var_index_in_read + if variant.is_insert() { variant.indel_length() + 1 } else { 1 } + MIN_CORE_DISTANCE - 1;
        let ref_base_length = n - var_index_in_read + if variant.is_delete() { variant.indel_length_abs() } else { 0 };
        let homology_ref_bases = ref_seq.base_range(variant.pos, variant.pos + ref_base_length)?;
        sc_adjust = check_indel_soft_clip_adjustment(read, variant, var_index_in_read);
        let mut ref_homology = None;
        if sc_adjust.is_some() {
            homology = find_homology(variant, homology_ref_bases, 0, false);
        } else {
            if var_index_in_read + 1 + if variant.is_insert() { variant.indel_length_abs() } else { 0 } < 0 {
                return None;
            }
            homology = find_homology(variant, read_bases, var_index_in_read, variant.is_insert());
            ref_homology = find_homology(variant, homology_ref_bases, 0, variant.is_delete());
        }
        if let Some(rh) = ref_homology {
            if homology.as_ref().map_or(true, |h| rh.length > h.length) {
                homology = Some(rh);
            }
        }
        if let Some(h) = &homology {
            read_core_end += h.length;
        }
    } else {
        read_core_start = var_index_in_read - MIN_CORE_DISTANCE;
        read_core_end = var_index_in_read + variant.alt_len() - 1 + MIN_CORE_DISTANCE;
    }

    if read_core_start < 0 || read_core_end >= n {
        return None;
    }

    let boundaries = find_repeat_boundaries(read_bases, read_core_start, read_core_end, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);
    read_core_start = read_core_start.min(boundaries.lower_index);
    read_core_end = read_core_end.max(boundaries.upper_index);

    let mut read_flank_start = read_core_start - flank_size;
    let mut read_flank_end = read_core_end + flank_size;
    if read_flank_start < 0 || read_flank_end >= n {
        return None;
    }

    let (cigar_start, cigar_elements): (i32, &[CigarEl]) = match &sc_adjust {
        Some(a) => (a.alignment_start, &a.cigar),
        None => (read.start, &read.cigar),
    };
    let info =
        ReadCigarInfo::build(cigar_start, cigar_elements, read_flank_start, read_core_start, read_core_end, read_flank_end);
    if !info.is_valid() {
        return None;
    }

    let read_position_start = info.read_alignment_start;
    read_flank_start = info.flank_index_start;
    read_flank_end = info.flank_index_end;
    if read_flank_start < 0 || read_flank_end >= n {
        return None;
    }
    let core_position_start = info.core_position_start;
    let core_position_end = info.core_position_end;
    let alignment_start = read_position_start.max(info.flank_position_start);
    let alignment_end = read.end.min(info.flank_position_end);

    let context_read_bases = read_bases[read_flank_start as usize..=read_flank_end as usize].to_vec();
    let offset = read_flank_start;
    let core_index_start = read_core_start - offset;
    let read_var_index = var_index_in_read - offset;
    let core_index_end = read_core_end - offset;

    let ref_bases = ref_seq.base_range(core_position_start, core_position_end)?.to_vec();

    let max_repeat = boundaries.max_repeat.as_ref().map(|r| RepeatInfo::new(r.index - offset, r.bases.clone(), r.count));
    let all_repeats: Vec<RepeatInfo> = boundaries
        .all_repeats
        .as_ref()
        .map(|l| l.iter().map(|r| RepeatInfo::new(r.index - offset, r.bases.clone(), r.count)).collect())
        .unwrap_or_default();

    let has_indel = has_indel_in_core(&info.cigar, core_index_start, core_index_end);
    Some(VariantReadContext {
        variant: variant.clone(),
        alignment_start,
        alignment_end,
        ref_bases,
        read_bases: context_read_bases,
        core_index_start,
        var_index: read_var_index,
        core_index_end,
        homology,
        max_repeat,
        all_repeats,
        core_position_start,
        core_position_end,
        read_cigar: cigar_str(&info.cigar),
        has_indel_in_core: has_indel,
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    })
}

fn above_min_base_qual(rc: &VariantReadContext, read: &Read, var_read_index: i32) -> bool {
    let start = var_read_index - rc.left_core_length() - rc.left_flank_length();
    let end = var_read_index + rc.right_core_length() + rc.right_flank_length();
    if start < 0 || end >= read.quals.len() as i32 {
        return false;
    }
    (start..=end).all(|i| !is_uncertain_base_qual(read.quals[i as usize]))
}

/// VariantReadContextBuilder.createContext (Illumina)
pub fn create_context(
    variant: &SimpleVariant, read: &Read, var_read_index: i32, ref_seq: &RefSequence, flank_size: i32,
) -> Option<VariantReadContext> {
    let mut rc = build_context(variant, read, var_read_index, ref_seq, flank_size)?;
    if rc.left_flank_length() < flank_size || rc.right_flank_length() < flank_size {
        return None;
    }
    if !above_min_base_qual(&rc, read, var_read_index) {
        return None;
    }
    if rc.variant_ref_index() <= 0 || rc.variant_ref_index() >= rc.ref_bases.len() as i32 - 1 {
        return None;
    }
    if variant.is_long_insert() {
        if let Some(b) = ref_seq.base_range(
            rc.core_position_start - rc.left_flank_length(),
            rc.core_position_end + rc.right_flank_length() + 1,
        ) {
            rc.extended_ref_bases = Some(b.to_vec());
        }
    }
    let ref_repeat_index = if variant.is_indel() { rc.variant_ref_index() + 1 } else { rc.variant_ref_index() };
    rc.max_ref_repeat = find_max_repeat(
        &rc.ref_bases, ref_repeat_index, ref_repeat_index, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1, false,
        ref_repeat_index,
    );
    rc.artefact = ArtefactContext::build(&rc).ok()?;
    Some(rc)
}
