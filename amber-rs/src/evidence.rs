//! The counting of reads at sites: RegionTask (with checkFragmentBoundaries), ReadCigarState and
//! PositionEvidenceChecker, over BAM records as htslib holds them.

use rust_htslib::bam;

use crate::types::AmberBase;

/// The counts at one site of one sample (the mutable part of PositionEvidence).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub read_depth: i32,
    pub indel_count: i32,
    pub ref_support: i32,
    pub alt_support: i32,
    pub base_qual_filtered: i32,
    pub map_qual_filtered: i32,
    pub seq_tech_filtered: i32,
}

/// A site and its counts: PositionEvidence without the LOH data.
#[derive(Clone, Copy, Debug)]
pub struct Evidence {
    pub chr: &'static str,
    pub pos: i32,
    pub ref_: AmberBase,
    pub alt: AmberBase,
    pub c: Counts,
}

impl Evidence {
    pub fn new(chr: &'static str, pos: i32, ref_: AmberBase, alt: AmberBase) -> Evidence {
        Evidence { chr, pos, ref_, alt, c: Counts::default() }
    }
    pub fn copy_site(&self) -> Evidence {
        Evidence::new(self.chr, self.pos, self.ref_, self.alt)
    }
    /// PositionEvidence.isValid
    pub fn is_valid(&self) -> bool {
        self.c.indel_count == 0
    }
    pub fn to_position_evidence(&self) -> crate::types::PositionEvidence {
        let mut pe = crate::types::PositionEvidence::new(self.chr, self.pos, self.ref_, self.alt);
        pe.read_depth = self.c.read_depth;
        pe.indel_count = self.c.indel_count;
        pe.ref_support = self.c.ref_support;
        pe.alt_support = self.c.alt_support;
        pe.base_qual_filtered = self.c.base_qual_filtered;
        pe.map_qual_filtered = self.c.map_qual_filtered;
        pe.seq_tech_filtered = self.c.seq_tech_filtered;
        pe
    }
}

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// CIGAR operations as BAM codes them
pub const OP_M: u8 = 0;
pub const OP_I: u8 = 1;
pub const OP_D: u8 = 2;
pub const OP_N: u8 = 3;
pub const OP_S: u8 = 4;
pub const OP_H: u8 = 5;
pub const OP_P: u8 = 6;
pub const OP_EQ: u8 = 7;
pub const OP_X: u8 = 8;

#[inline]
fn consumes_ref(op: u8) -> bool {
    matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X)
}
#[inline]
fn consumes_read(op: u8) -> bool {
    matches!(op, OP_M | OP_I | OP_S | OP_EQ | OP_X)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CigarElem {
    pub len: i32,
    pub op: u8,
}

/// CigarUtils.cigarElementsFromStr: digits then an operator character, "*" for none. None for an
/// operator htsjdk does not know (CigarOperator.characterToEnum throws).
pub fn cigar_from_str(s: &[u8]) -> Option<Vec<CigarElem>> {
    if s == b"*" {
        return Some(Vec::new());
    }
    let mut out = Vec::with_capacity(4);
    let mut len: i32 = 0;
    for &c in s {
        if c.is_ascii_digit() {
            len = len.wrapping_mul(10).wrapping_add((c - b'0') as i32);
        } else {
            let op = match c {
                b'M' => OP_M,
                b'I' => OP_I,
                b'D' => OP_D,
                b'N' => OP_N,
                b'S' => OP_S,
                b'H' => OP_H,
                b'P' => OP_P,
                b'=' => OP_EQ,
                b'X' => OP_X,
                _ => return None,
            };
            out.push(CigarElem { len, op });
            len = 0;
        }
    }
    Some(out)
}

/// A BAM record's fields, read from htslib's buffer without copying.
pub struct RecView<'a> {
    /// 1-based alignment start
    pub start: i32,
    pub flags: u16,
    pub mapq: u8,
    pub tid: i32,
    pub mtid: i32,
    /// 1-based mate alignment start
    pub mate_start: i32,
    cigar: &'a [u8],
    n_cigar: usize,
    seq: &'a [u8],
    pub seq_len: usize,
    quals: &'a [u8],
    aux: &'a [u8],
}

const SEQ_NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

impl<'a> RecView<'a> {
    pub fn of(rec: &'a bam::Record) -> RecView<'a> {
        let b = rec.inner();
        let c = &b.core;
        // SAFETY: htslib guarantees `data` holds `l_data` bytes, and the slice borrows the record
        let data: &[u8] = if b.l_data > 0 { unsafe { std::slice::from_raw_parts(b.data, b.l_data as usize) } } else { &[] };
        let name_len = c.l_qname as usize;
        let n_cigar = c.n_cigar as usize;
        let seq_len = c.l_qseq.max(0) as usize;
        let cigar_end = name_len + n_cigar * 4;
        let seq_end = cigar_end + seq_len.div_ceil(2);
        let qual_end = seq_end + seq_len;
        RecView {
            start: c.pos as i32 + 1,
            flags: c.flag,
            mapq: c.qual,
            tid: c.tid,
            mtid: c.mtid,
            mate_start: c.mpos as i32 + 1,
            cigar: &data[name_len..cigar_end],
            n_cigar,
            seq: &data[cigar_end..seq_end],
            seq_len,
            quals: &data[seq_end..qual_end],
            aux: data.get(qual_end..).unwrap_or(&[]),
        }
    }
    #[inline]
    pub fn cigar(&self, i: usize) -> CigarElem {
        let c = &self.cigar[i * 4..i * 4 + 4];
        let v = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
        CigarElem { len: (v >> 4) as i32, op: (v & 0xf) as u8 }
    }
    pub fn n_cigar(&self) -> usize {
        self.n_cigar
    }
    /// SAMRecord.getAlignmentEnd: the last reference base, start - 1 for no reference bases
    pub fn alignment_end(&self) -> i32 {
        let mut len = 0;
        for i in 0..self.n_cigar {
            let e = self.cigar(i);
            if consumes_ref(e.op) {
                len += e.len;
            }
        }
        self.start + len - 1
    }
    /// the read base at an index, as htsjdk decodes it
    #[inline]
    pub fn base(&self, i: usize) -> u8 {
        let b = self.seq[i / 2];
        SEQ_NT16[(if i % 2 == 0 { b >> 4 } else { b & 0xf }) as usize]
    }
    /// the base quality at an index, as a Java byte
    #[inline]
    pub fn qual(&self, i: usize) -> i32 {
        self.quals[i] as i8 as i32
    }
    pub fn has_flag(&self, f: u16) -> bool {
        self.flags & f != 0
    }
    /// A string (Z) tag's value
    pub fn string_tag(&self, tag: &[u8; 2]) -> Option<&'a [u8]> {
        let mut a = self.aux;
        while a.len() >= 3 {
            let (key, ty) = (&a[0..2], a[2]);
            let rest = &a[3..];
            let size = match ty {
                b'A' | b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                b'Z' | b'H' => {
                    let end = rest.iter().position(|&x| x == 0)?;
                    if key == tag {
                        return if ty == b'Z' { Some(&rest[..end]) } else { None };
                    }
                    a = &rest[end + 1..];
                    continue;
                }
                b'B' => {
                    if rest.len() < 5 {
                        return None;
                    }
                    let sub = match rest[0] {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        _ => 4,
                    };
                    let n = u32::from_le_bytes([rest[1], rest[2], rest[3], rest[4]]) as usize;
                    5 + sub * n
                }
                _ => return None,
            };
            if rest.len() < size {
                return None;
            }
            a = &rest[size..];
        }
        None
    }
}

/// ReadCigarState: walks the CIGAR of a read to the reference positions of the sites in order.
///
/// AM2: Java starts at the first element (or the second after a leading soft clip) and adds the length
/// of every element that does not consume the reference to the read index, so a hard clip (or padding)
/// moves the read index although it holds no read bases, and a leading insertion, or a soft clip after a
/// hard clip, leaves the walk one base out: the base read for a site is then another base of the read.
/// Corrected: elements that hold no read bases do not move the read index, and the walk starts at the
/// first element that consumes the reference.
pub struct ReadCigarState {
    /// None for a single aligned element
    single: bool,
    pub base_length: i32,
    pub read_index: i32,
    pub ref_position: i32,
    pub cigar_index: i32,
    pub element_index: i32,
    java: bool,
}

impl ReadCigarState {
    pub fn new(elements: &[CigarElem], read_start: i32, base_length: i32) -> ReadCigarState {
        let java = crate::compat::java("AM2");
        let single = elements.len() == 1 && elements[0].op == OP_M;
        let mut s =
            ReadCigarState { single, base_length, read_index: 0, ref_position: read_start, cigar_index: 0, element_index: 0, java };
        if java {
            if elements.first().is_some_and(|e| e.op == OP_S) {
                s.cigar_index = 1;
                s.read_index = elements[0].len;
            }
        } else {
            while (s.cigar_index as usize) < elements.len() && !consumes_ref(elements[s.cigar_index as usize].op) {
                let e = elements[s.cigar_index as usize];
                if consumes_read(e.op) {
                    s.read_index += e.len;
                }
                s.cigar_index += 1;
            }
        }
        s
    }

    pub fn move_to_position(&mut self, elements: &[CigarElem], position: i32) {
        if position <= self.ref_position {
            return;
        }
        if self.single {
            let pos_difference = position - self.ref_position;
            let remaining = self.base_length - self.element_index - 1;
            if pos_difference > remaining {
                self.read_index = self.base_length;
                return;
            }
            self.ref_position += pos_difference;
            self.read_index += pos_difference;
            self.element_index += pos_difference;
            return;
        }
        while (self.cigar_index as usize) < elements.len() {
            let e = elements[self.cigar_index as usize];
            if !consumes_ref(e.op) {
                if self.java || consumes_read(e.op) {
                    self.read_index += e.len;
                }
                self.element_index = -1;
                self.cigar_index += 1;
                continue;
            }
            let pos_difference = position - self.ref_position;
            let remaining = e.len - self.element_index - 1;
            let pos_move = pos_difference.min(remaining);
            self.ref_position += pos_move;
            if consumes_read(e.op) {
                self.read_index += pos_move;
            }
            self.element_index += pos_move;
            if self.read_index >= self.base_length {
                return;
            }
            if self.ref_position == position {
                return;
            }
            self.element_index = -1;
            self.cigar_index += 1;
        }
    }

    pub fn is_valid(&self, elements: &[CigarElem]) -> bool {
        if self.single {
            return self.read_index < self.base_length;
        }
        self.cigar_index >= 0 && (self.cigar_index as usize) < elements.len() && self.read_index < self.base_length
    }
    pub fn is_delete(&self, elements: &[CigarElem]) -> bool {
        if !self.is_valid(elements) || self.single {
            return false;
        }
        elements[self.cigar_index as usize].op == OP_D
    }
    /// at the end of an aligned element with an insert following
    pub fn is_insert(&self, elements: &[CigarElem]) -> bool {
        if !self.is_valid(elements) || self.single {
            return false;
        }
        let ci = self.cigar_index as usize;
        let e = elements[ci];
        if e.op == OP_M && self.element_index == e.len - 1 && ci < elements.len() - 1 {
            return elements[ci + 1].op == OP_I;
        }
        false
    }
    pub fn is_indel(&self, elements: &[CigarElem]) -> bool {
        self.is_delete(elements) || self.is_insert(elements)
    }
}

/// RegionTask.checkFragmentBoundaries: for concordant overlapping pairs, the part of this read whose
/// bases are counted, so that bases both reads cover are counted once and adapter bases not at all.
pub fn check_fragment_boundaries(rec: &RecView, read_end: i32) -> Option<(i32, i32)> {
    if !rec.has_flag(FLAG_PAIRED) || rec.has_flag(FLAG_MATE_UNMAPPED) {
        return None;
    }
    if rec.tid != rec.mtid {
        return None;
    }
    let reverse = rec.has_flag(FLAG_REVERSE);
    if reverse == rec.has_flag(FLAG_MATE_REVERSE) {
        return None;
    }
    if read_end < rec.mate_start {
        return None;
    }
    let mate_cigar = rec.string_tag(b"MC")?;
    let mate_elements = match cigar_from_str(mate_cigar) {
        Some(e) => e,
        None => {
            crate::error!("thread exception: invalid mate cigar {}", String::from_utf8_lossy(mate_cigar));
            std::process::exit(1);
        }
    };
    fragment_bounds(rec.start, read_end, rec.mate_start, &mate_elements, reverse)
}

fn fragment_bounds(read_start: i32, read_end: i32, mate_start: i32, mate: &[CigarElem], reverse: bool) -> Option<(i32, i32)> {
    let (mut read_start, mut read_end) = (read_start, read_end);
    let mate_end = mate_start + mate.iter().filter(|e| consumes_ref(e.op)).map(|e| e.len).sum::<i32>() - 1;
    if mate_end < read_start {
        return None;
    }
    if mate.is_empty() {
        // Java takes the first element of the empty list and fails
        crate::error!("thread exception: java.lang.IndexOutOfBoundsException: empty mate cigar");
        std::process::exit(1);
    }
    let is_first_read = if read_start < mate_start {
        true
    } else if read_start == mate_start {
        if read_end == mate_end { !reverse } else { read_end < mate_end }
    } else {
        false
    };
    if !is_first_read {
        read_start = read_start.max(mate_end + 1);
    }
    if reverse {
        let mut unclipped_mate_start = mate_start;
        if mate[0].op == OP_S {
            unclipped_mate_start -= mate[0].len;
        }
        read_start = unclipped_mate_start.max(read_start);
    } else {
        let mut unclipped_mate_end = mate_end;
        let last = mate[mate.len() - 1];
        if last.op == OP_S {
            unclipped_mate_end += last.len;
        }
        read_end = unclipped_mate_end.min(read_end);
    }
    Some((read_start, read_end))
}

/// The Ultima low-quality read indices of the UQ tag ("count=a,b-c,..."), parsed once per read.
pub enum LowQual {
    Unread,
    None,
    Ranges(Vec<(i32, i32)>),
}

fn parse_low_qual(tag: &[u8]) -> Option<Vec<(i32, i32)>> {
    let s = std::str::from_utf8(tag).ok()?;
    let (count, items) = s.split_once('=')?;
    crate::config::java_parse_int(count)?;
    let mut out = Vec::new();
    for item in items.split(',') {
        if let Some((a, b)) = item.split_once('-') {
            out.push((crate::config::java_parse_int(a)?, crate::config::java_parse_int(b)?));
        } else {
            let v = crate::config::java_parse_int(item)?;
            out.push((v, v));
        }
    }
    Some(out)
}

impl LowQual {
    fn contains(&mut self, rec: &RecView, index: i32) -> bool {
        if let LowQual::Unread = self {
            *self = match rec.string_tag(b"UQ") {
                None => LowQual::None,
                Some(tag) => match parse_low_qual(tag) {
                    Some(r) => LowQual::Ranges(r),
                    None => {
                        crate::error!("thread exception: invalid UQ tag {}", String::from_utf8_lossy(tag));
                        std::process::exit(1);
                    }
                },
            };
        }
        match self {
            LowQual::Ranges(r) => r.iter().any(|&(a, b)| index >= a && index <= b),
            _ => false,
        }
    }
}

/// PositionEvidenceChecker
#[derive(Clone, Copy)]
pub struct EvidenceChecker {
    pub min_map_quality: i32,
    pub min_base_quality: i32,
    pub ultima: bool,
}

impl EvidenceChecker {
    /// addEvidence: one read at one site
    pub fn add_evidence(&self, ev: &mut Evidence, rec: &RecView, elements: &[CigarElem], low_qual: &mut LowQual, read_index: i32, is_indel: bool) {
        let base_quality = if is_indel { indel_base_quality(ev.pos, rec, elements) } else { rec.qual(read_index as usize) };
        let mut filtered = false;
        if (rec.mapq as i32) < self.min_map_quality {
            ev.c.map_qual_filtered += 1;
            filtered = true;
        }
        if base_quality < self.min_base_quality {
            ev.c.base_qual_filtered += 1;
            filtered = true;
        }
        if self.ultima && low_qual.contains(rec, read_index) {
            ev.c.seq_tech_filtered += 1;
            filtered = true;
        }
        ev.c.read_depth += 1;
        if filtered {
            return;
        }
        if read_index >= 0 {
            if !is_indel {
                let base = rec.base(read_index as usize);
                match AmberBase::from_char(base) {
                    Some(b) if b == ev.ref_ => ev.c.ref_support += 1,
                    Some(b) if b == ev.alt => ev.c.alt_support += 1,
                    Some(_) => {}
                    None => {
                        // AM3: AmberBase.valueOf throws for a read base other than G, A, T, C or N (an IUPAC
                        // code or '='), and the uncaught exception ends the run. Corrected: the base
                        // supports neither allele.
                        if crate::compat::java("AM3") {
                            crate::error!(
                                "thread exception: java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.amber.AmberBase.{}",
                                base as char
                            );
                            std::process::exit(1);
                        }
                    }
                }
            } else {
                ev.c.indel_count += 1;
            }
        }
    }
}

/// getBaseQuality: the quality of the first aligned base at or after the position (the base after a
/// deletion), 0 when there is none; SAMRecord.getReadPositionAtReferencePosition over the alignment blocks
fn indel_base_quality(position: i32, rec: &RecView, elements: &[CigarElem]) -> i32 {
    let end = rec.alignment_end();
    let mut pos = position;
    while pos <= end {
        let rp = read_position_at_reference_position(rec.start, elements, pos);
        if rp > 0 {
            return rec.qual((rp - 1) as usize);
        }
        pos += 1;
    }
    0
}

/// htsjdk SAMRecord.getReadPositionAtReferencePosition(pos): 1-based read position of an aligned base,
/// 0 for a deleted or uncovered position
pub fn read_position_at_reference_position(start: i32, elements: &[CigarElem], pos: i32) -> i32 {
    if pos <= 0 {
        return 0;
    }
    let (mut read_base, mut ref_base) = (1, start);
    for e in elements {
        match e.op {
            OP_M | OP_EQ | OP_X => {
                // an alignment block from ref_base, read_base, of e.len bases
                if ref_base + e.len - 1 >= pos {
                    return if pos < ref_base { 0 } else { pos - ref_base + read_base };
                }
                read_base += e.len;
                ref_base += e.len;
            }
            OP_S | OP_I => read_base += e.len,
            OP_D | OP_N => ref_base += e.len,
            _ => {}
        }
    }
    0
}

/// RegionTask: the sites of one region, the index below which every site lies before the reads still
/// to come, and whether all are passed.
pub struct RegionTask {
    pub current_index: usize,
    pub complete: bool,
}

impl RegionTask {
    pub fn new() -> RegionTask {
        RegionTask { current_index: 0, complete: false }
    }

    /// processRecord
    pub fn process_record(&mut self, positions: &mut [Evidence], rec: &RecView, checker: &EvidenceChecker, elements: &mut Vec<CigarElem>) {
        let read_alignment_start = rec.start;
        let read_end = rec.alignment_end();
        let (alignment_start, alignment_end) = match check_fragment_boundaries(rec, read_end) {
            Some((s, e)) => {
                if s > e {
                    return;
                }
                (s, e)
            }
            None => (read_alignment_start, read_end),
        };
        let mut state: Option<ReadCigarState> = None;
        let mut low_qual = LowQual::Unread;
        let mut index = self.current_index;
        while index < positions.len() {
            let pos = positions[index].pos;
            if read_alignment_start > pos {
                self.current_index += 1;
                index += 1;
                continue;
            }
            if alignment_start > pos {
                index += 1;
                continue;
            }
            if alignment_end < pos {
                break;
            }
            if state.is_none() {
                elements.clear();
                elements.extend((0..rec.n_cigar()).map(|i| rec.cigar(i)));
                state = Some(ReadCigarState::new(elements, read_alignment_start, rec.seq_len as i32));
            }
            let s = state.as_mut().unwrap();
            s.move_to_position(elements, pos);
            if !s.is_valid(elements) {
                break;
            }
            let (read_index, is_indel) = (s.read_index, s.is_indel(elements));
            checker.add_evidence(&mut positions[index], rec, elements, &mut low_qual, read_index, is_indel);
            index += 1;
        }
        if self.current_index >= positions.len() {
            self.complete = true;
        }
    }
}

impl Default for RegionTask {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, CigarString};

    fn elems(s: &str) -> Vec<CigarElem> {
        cigar_from_str(s.as_bytes()).unwrap()
    }

    /// A synthetic paired read (SamRecordTestUtils.createSamRecord): chromosome index 0 or 1
    #[allow(clippy::too_many_arguments)]
    pub fn read(
        tid: i32, start: i32, bases: &str, cigar: &str, mate_tid: i32, mate_start: i32, reverse: bool, mate_reverse: bool,
        mate_cigar: Option<&str>,
    ) -> bam::Record {
        let mut r = bam::Record::new();
        let cigar = CigarString::try_from(cigar.as_bytes()).unwrap();
        let quals = vec![37u8; bases.len()];
        r.set(b"READ_001", Some(&cigar), bases.as_bytes(), &quals);
        r.set_tid(tid);
        r.set_pos((start - 1) as i64);
        r.set_mtid(mate_tid);
        r.set_mpos((mate_start - 1) as i64);
        r.set_mapq(60);
        let mut flags = FLAG_PAIRED;
        if reverse {
            flags |= FLAG_REVERSE;
        }
        if mate_reverse {
            flags |= FLAG_MATE_REVERSE;
        }
        r.set_flags(flags);
        if let Some(mc) = mate_cigar {
            r.push_aux(b"MC", Aux::String(mc)).unwrap();
        }
        r
    }

    fn bounds(r: &bam::Record) -> Option<(i32, i32)> {
        let v = RecView::of(r);
        check_fragment_boundaries(&v, v.alignment_end())
    }

    #[test]
    fn read_cigar_state() {
        crate::compat::with_java(true, || {
            let e = elems("10M");
            let mut s = ReadCigarState::new(&e, 100, 10);
            assert_eq!((s.cigar_index, s.element_index, s.read_index, s.ref_position), (0, 0, 0, 100));
            assert!(s.is_valid(&e) && !s.is_delete(&e) && !s.is_insert(&e));
            s.move_to_position(&e, 109);
            assert_eq!((s.cigar_index, s.element_index, s.read_index, s.ref_position), (0, 9, 9, 109));
            s.move_to_position(&e, 110);
            assert!(!s.is_valid(&e));

            let e = elems("5M2D5M2I5M");
            let mut s = ReadCigarState::new(&e, 100, 20);
            assert_eq!((s.cigar_index, s.element_index, s.read_index, s.ref_position), (0, 0, 0, 100));
            s.move_to_position(&e, 104);
            assert_eq!((s.ref_position, s.cigar_index, s.element_index, s.read_index), (104, 0, 4, 4));
            assert!(!s.is_delete(&e));
            s.move_to_position(&e, 105);
            assert_eq!((s.ref_position, s.cigar_index, s.element_index, s.read_index), (105, 1, 0, 4));
            assert!(s.is_delete(&e));
            s.move_to_position(&e, 106);
            assert_eq!((s.ref_position, s.cigar_index, s.element_index, s.read_index), (106, 1, 1, 4));
            assert!(s.is_delete(&e));
            s.move_to_position(&e, 107);
            assert_eq!((s.ref_position, s.cigar_index, s.element_index, s.read_index), (107, 2, 0, 5));
            assert!(!s.is_delete(&e));
            s.move_to_position(&e, 111);
            assert_eq!((s.ref_position, s.cigar_index, s.element_index, s.read_index), (111, 2, 4, 9));
            assert!(!s.is_delete(&e) && s.is_insert(&e));
            s.move_to_position(&e, 112);
            assert_eq!((s.ref_position, s.cigar_index, s.element_index, s.read_index), (112, 4, 0, 12));
            assert!(!s.is_insert(&e));
            s.move_to_position(&e, 116);
            assert_eq!((s.ref_position, s.cigar_index, s.element_index, s.read_index), (116, 4, 4, 16));

            let e = elems("5S10M10I10M10D10M5S");
            let mut s = ReadCigarState::new(&e, 100, 50);
            assert_eq!((s.cigar_index, s.element_index, s.read_index, s.ref_position), (1, 0, 5, 100));
            s.move_to_position(&e, 105);
            assert_eq!((s.read_index, s.ref_position), (10, 105));
            s.move_to_position(&e, 109);
            assert_eq!((s.read_index, s.ref_position), (14, 109));
            assert!(!s.is_delete(&e) && s.is_insert(&e));
            s.move_to_position(&e, 110);
            assert_eq!((s.read_index, s.ref_position), (25, 110));
            assert!(!s.is_delete(&e) && !s.is_insert(&e));
            s.move_to_position(&e, 119);
            assert_eq!((s.read_index, s.ref_position), (34, 119));
            s.move_to_position(&e, 120);
            assert_eq!((s.read_index, s.ref_position), (34, 120));
            assert!(s.is_delete(&e) && !s.is_insert(&e));
            s.move_to_position(&e, 129);
            assert_eq!((s.read_index, s.ref_position), (34, 129));
            assert!(s.is_delete(&e));
            s.move_to_position(&e, 130);
            assert_eq!((s.read_index, s.ref_position), (35, 130));
            assert!(!s.is_delete(&e) && !s.is_insert(&e));
            s.move_to_position(&e, 139);
            assert_eq!((s.read_index, s.ref_position), (44, 139));
            s.move_to_position(&e, 140);
            assert!(!s.is_valid(&e));
        });
    }

    /// AM2: a leading hard clip or insertion moves Java's read index off the base at the site
    #[test]
    fn read_cigar_state_clips() {
        // 5H10M: reference 100-109 are read bases 0-9
        let e = elems("5H10M");
        let walk = |java: bool, pos: i32| {
            crate::compat::with_java(java, || {
                let mut s = ReadCigarState::new(&e, 100, 10);
                s.move_to_position(&e, pos);
                (s.read_index, s.is_valid(&e))
            })
        };
        assert_eq!(walk(false, 103), (3, true));
        assert_eq!(walk(true, 103), (8, true));
        assert_eq!(walk(false, 100), (0, true));
        // 2I10M: reference 100-109 are read bases 2-11
        let e = elems("2I10M");
        let walk = |java: bool, pos: i32| {
            crate::compat::with_java(java, || {
                let mut s = ReadCigarState::new(&e, 100, 12);
                s.move_to_position(&e, pos);
                s.read_index
            })
        };
        assert_eq!(walk(false, 100), 2);
        assert_eq!(walk(false, 105), 7);
        // Java reads the inserted base at the start, and its element index is one short
        assert_eq!(walk(true, 100), 0);
        assert_eq!(walk(true, 105), 7);
        // 2I5M3D5M: reference 105 is deleted; Java lands on the base after the deletion
        let e = elems("2I5M3D5M");
        let at = |java: bool, pos: i32| {
            crate::compat::with_java(java, || {
                let mut s = ReadCigarState::new(&e, 100, 12);
                s.move_to_position(&e, pos);
                (s.read_index, s.is_delete(&e))
            })
        };
        assert_eq!(at(false, 105), (6, true));
        assert_eq!(at(true, 105), (7, false));
        // 3H2S10M: reference 100 is read base 2
        let e = elems("3H2S10M");
        let s = crate::compat::with_java(false, || ReadCigarState::new(&e, 100, 12));
        assert_eq!(s.read_index, 2);
        // ordinary CIGARs walk the same either way
        for cig in ["10M", "3S7M", "4M2D6M", "4M2I4M2S"] {
            let e = elems(cig);
            for pos in 100..112 {
                let a = crate::compat::with_java(true, || {
                    let mut s = ReadCigarState::new(&e, 100, 10);
                    s.move_to_position(&e, pos);
                    (s.read_index, s.is_valid(&e), s.is_indel(&e))
                });
                let b = crate::compat::with_java(false, || {
                    let mut s = ReadCigarState::new(&e, 100, 10);
                    s.move_to_position(&e, pos);
                    (s.read_index, s.is_valid(&e), s.is_indel(&e))
                });
                assert_eq!(a, b, "{cig} {pos}");
            }
        }
    }

    #[test]
    fn fragment_read_sync() {
        let bases = "A".repeat(20);
        let r = read(0, 100, &bases, "20M", 0, 110, false, true, Some("20M"));
        assert_eq!(bounds(&r), Some((100, 119)));
        // the mate only processes new / non-overlapping bases
        let r = read(0, 110, &bases, "20M", 0, 100, true, false, Some("20M"));
        assert_eq!(bounds(&r), Some((120, 129)));
        // two reads overlapping at the start
        let r = read(0, 100, &bases, "20M", 0, 100, false, true, Some("25M"));
        assert_eq!(bounds(&r), Some((100, 119)));
        let r = read(0, 100, &"A".repeat(25), "25M", 0, 100, true, false, Some("20M"));
        assert_eq!(bounds(&r), Some((120, 124)));
        // 3' adapter trimming
        let r = read(0, 100, &bases, "20M", 0, 102, true, false, Some("20M"));
        assert_eq!(bounds(&r), Some((102, 119)));
        let r = read(0, 100, &bases, "20M", 0, 110, true, false, Some("5S20M"));
        assert_eq!(bounds(&r), Some((105, 119)));
        let r = read(0, 100, &bases, "20M", 0, 90, false, true, Some("20M5S"));
        assert_eq!(bounds(&r), Some((110, 114)));
        // no sync: BND, INV, unmapped mate, beyond the other read
        assert_eq!(bounds(&read(0, 100, &bases, "20M", 1, 110, false, true, Some("20M"))), None);
        assert_eq!(bounds(&read(0, 100, &bases, "20M", 0, 110, false, false, Some("20M"))), None);
        let mut r = read(0, 100, &bases, "20M", 0, 110, false, true, Some("20M"));
        r.set_flags(r.flags() | FLAG_MATE_UNMAPPED);
        assert_eq!(bounds(&r), None);
        assert_eq!(bounds(&read(0, 150, &bases, "20M", 0, 110, false, true, Some("20M"))), None);
        assert_eq!(bounds(&read(0, 100, &bases, "20M", 0, 110, false, true, None)), None);
    }

    fn checker() -> EvidenceChecker {
        EvidenceChecker { min_map_quality: 50, min_base_quality: 18, ultima: false }
    }

    /// TumorEvidenceTest.useQualityOfBaseAfterDel
    #[test]
    fn quality_of_base_after_del() {
        let min_quality = (b'J' - 33) as i32;
        let mut low_base_qual_del = read(0, 1000, "CT", "1M1D1M", 0, 1200, false, true, None);
        let cig = CigarString::try_from("1M1D1M".as_bytes()).unwrap();
        low_base_qual_del.set(b"READ_01", Some(&cig), b"CT", &[(b'F' - 33), (b'I' - 33)]);
        let mut low_map_qual = read(0, 1000, "CT", "1M1D1M", 0, 1200, false, true, None);
        low_map_qual.set_mapq(40);
        let mut ev = Evidence::new("5", 1001, AmberBase::A, AmberBase::T);
        let check = EvidenceChecker { min_map_quality: 50, min_base_quality: min_quality, ultima: false };
        for (r, depth, bq, mq) in [(&low_base_qual_del, 1, Some(1), 0), (&low_map_qual, 2, None, 1)] {
            let v = RecView::of(r);
            let e: Vec<CigarElem> = (0..v.n_cigar()).map(|i| v.cigar(i)).collect();
            let mut s = ReadCigarState::new(&e, v.start, v.seq_len as i32);
            s.move_to_position(&e, 1001);
            check.add_evidence(&mut ev, &v, &e, &mut LowQual::Unread, s.read_index, s.is_indel(&e));
            assert_eq!((ev.c.read_depth, ev.c.map_qual_filtered), (depth, mq));
            if let Some(bq) = bq {
                assert_eq!(ev.c.base_qual_filtered, bq);
            }
        }
    }

    /// TumorEvidenceTest.testFragmentSync
    #[test]
    fn fragment_sync() {
        let mut positions = vec![Evidence::new("chr1", 1020, AmberBase::A, AmberBase::C), Evidence::new("chr1", 1025, AmberBase::A, AmberBase::C)];
        let mut task = RegionTask::new();
        let mut scratch = Vec::new();
        let bases = "A".repeat(30);
        // both reads overlap the same 2 positions but are only counted once
        for r in [read(0, 1000, &bases, "30M", 0, 1010, false, true, Some("30M")), read(0, 1010, &bases, "30M", 0, 1000, true, false, Some("30M"))] {
            task.process_record(&mut positions, &RecView::of(&r), &checker(), &mut scratch);
        }
        assert_eq!((positions[0].c.ref_support, positions[0].c.read_depth), (1, 1));
        assert_eq!((positions[1].c.ref_support, positions[1].c.read_depth), (1, 1));
        for r in [read(0, 1010, &bases, "30M", 0, 1021, false, true, Some("30M")), read(0, 1021, &bases, "30M", 0, 1010, true, false, Some("30M"))] {
            task.process_record(&mut positions, &RecView::of(&r), &checker(), &mut scratch);
        }
        assert_eq!((positions[0].c.ref_support, positions[0].c.read_depth), (2, 2));
        assert_eq!((positions[1].c.ref_support, positions[1].c.read_depth), (2, 2));
        assert!(!task.complete);
        // a read past both sites completes the task
        task.process_record(&mut positions, &RecView::of(&read(0, 1030, &bases, "30M", 0, 1000, true, false, None)), &checker(), &mut scratch);
        assert!(task.complete && task.current_index == 2);
    }

    #[test]
    fn evidence_counts() {
        // ref, alt, other base, low base quality, insert after the site, deletion at the site
        let mut ev = Evidence::new("chr1", 105, AmberBase::C, AmberBase::T);
        let mut task = RegionTask::new();
        let mut scratch = Vec::new();
        let reads = [("AAAAACAAAA", "10M"), ("AAAAATAAAA", "10M"), ("AAAAAGAAAA", "10M"), ("AAAAACAAAAA", "6M1I4M"), ("AAAAAAAAA", "5M1D4M")];
        for (bases, cigar) in reads {
            let r = read(0, 100, bases, cigar, 0, 500, false, true, None);
            task.process_record(std::slice::from_mut(&mut ev), &RecView::of(&r), &checker(), &mut scratch);
        }
        let mut lowq = read(0, 100, "AAAAACAAAA", "10M", 0, 500, false, true, None);
        let cig = CigarString::try_from("10M".as_bytes()).unwrap();
        lowq.set(b"r", Some(&cig), b"AAAAACAAAA", &[37, 37, 37, 37, 37, 10, 37, 37, 37, 37]);
        task.process_record(std::slice::from_mut(&mut ev), &RecView::of(&lowq), &checker(), &mut scratch);
        assert_eq!(ev.c, Counts { read_depth: 6, indel_count: 2, ref_support: 1, alt_support: 1, base_qual_filtered: 1, map_qual_filtered: 0, seq_tech_filtered: 0 });
    }

    #[test]
    fn ultima_low_qual_and_iupac() {
        let mut ev = Evidence::new("chr1", 102, AmberBase::C, AmberBase::T);
        let mut r = read(0, 100, "AACAA", "5M", 0, 500, false, true, None);
        r.push_aux(b"UQ", Aux::String("2=0-2,4")).unwrap();
        let ultima = EvidenceChecker { ultima: true, ..checker() };
        let mut task = RegionTask::new();
        task.process_record(std::slice::from_mut(&mut ev), &RecView::of(&r), &ultima, &mut Vec::new());
        assert_eq!((ev.c.read_depth, ev.c.seq_tech_filtered, ev.c.ref_support), (1, 1, 0));
        let mut task = RegionTask::new();
        task.process_record(std::slice::from_mut(&mut ev), &RecView::of(&r), &checker(), &mut Vec::new());
        assert_eq!((ev.c.read_depth, ev.c.ref_support), (2, 1));
        // an IUPAC base supports neither allele (AM3, corrected)
        let r = read(0, 100, "AAMAA", "5M", 0, 500, false, true, None);
        let mut task = RegionTask::new();
        task.process_record(std::slice::from_mut(&mut ev), &RecView::of(&r), &checker(), &mut Vec::new());
        assert_eq!((ev.c.read_depth, ev.c.ref_support, ev.c.alt_support), (3, 1, 0));
        assert_eq!(parse_low_qual(b"3=1,5-7"), Some(vec![(1, 1), (5, 7)]));
        assert_eq!(parse_low_qual(b"3"), None);
    }

    #[test]
    fn read_positions() {
        let e = elems("2S3M2D3M");
        assert_eq!(read_position_at_reference_position(100, &e, 100), 3);
        assert_eq!(read_position_at_reference_position(100, &e, 103), 0);
        assert_eq!(read_position_at_reference_position(100, &e, 105), 6);
        assert_eq!(read_position_at_reference_position(100, &e, 108), 0);
        assert_eq!(cigar_from_str(b"*"), Some(vec![]));
        assert_eq!(cigar_from_str(b"5S10M"), Some(vec![CigarElem { len: 5, op: OP_S }, CigarElem { len: 10, op: OP_M }]));
        assert_eq!(cigar_from_str(b"5Q"), None);
    }
}
