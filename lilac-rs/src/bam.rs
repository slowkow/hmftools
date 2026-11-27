//! BAM records as LILAC sees them (the htsjdk SAMRecord methods it calls), BamSlicer, and read.Read.

use std::sync::Arc;

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as _};

use crate::config::Params;
use crate::constants::*;
use crate::indel::Indel;
use crate::refdata::RefBases;

/// The fields of a BAM record LILAC uses
#[derive(Clone, Debug)]
pub struct Rec {
    pub name: String,
    pub flags: u16,
    pub contig: String,
    /// 1-based alignment start
    pub start: i32,
    /// 1-based alignment end (last reference base)
    pub end: i32,
    pub cigar: Vec<(u8, u32)>,
    pub seq: Vec<u8>,
    pub quals: Vec<u8>,
    pub mapq: u8,
    /// the mate's contig name ("*" for none)
    pub mate_contig: String,
    pub mate_tid: i32,
    /// 1-based, 0 for none
    pub mate_start: i32,
    pub isize: i64,
    /// the Ultima UQ tag
    pub uq: Option<String>,
}

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_FIRST: u16 = 0x40;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

fn consumes_ref(op: u8) -> bool {
    matches!(op, b'M' | b'D' | b'N' | b'=' | b'X')
}
fn consumes_read(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'S' | b'=' | b'X')
}

impl Rec {
    pub fn from_bam(r: &bam::Record, header: &bam::HeaderView) -> Rec {
        let contig = if r.tid() >= 0 { String::from_utf8_lossy(header.tid2name(r.tid() as u32)).into_owned() } else { "*".to_string() };
        let mate_contig = if r.mtid() >= 0 { String::from_utf8_lossy(header.tid2name(r.mtid() as u32)).into_owned() } else { "*".to_string() };
        let cigar: Vec<(u8, u32)> = r.cigar().iter().map(|c| (c.char() as u8, c.len())).collect();
        let start = r.pos() as i32 + 1;
        let ref_len: i32 = cigar.iter().filter(|c| consumes_ref(c.0)).map(|c| c.1 as i32).sum();
        let end = start + ref_len - 1;
        let uq = match r.aux(b"UQ") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => Some(s.to_string()),
            _ => None,
        };
        Rec {
            name: String::from_utf8_lossy(r.qname()).into_owned(),
            flags: r.flags(),
            contig,
            start,
            end,
            cigar,
            seq: r.seq().as_bytes(),
            quals: r.qual().to_vec(),
            mapq: r.mapq(),
            mate_contig,
            mate_tid: r.mtid(),
            mate_start: r.mpos() as i32 + 1,
            isize: r.insert_size(),
            uq,
        }
    }

    #[inline]
    pub fn flag(&self, f: u16) -> bool {
        self.flags & f != 0
    }
    pub fn paired(&self) -> bool {
        self.flag(FLAG_PAIRED)
    }
    pub fn reverse(&self) -> bool {
        self.flag(FLAG_REVERSE)
    }
    pub fn first_in_pair(&self) -> bool {
        self.flag(FLAG_FIRST)
    }
    pub fn supplementary(&self) -> bool {
        self.flag(FLAG_SUPPLEMENTARY)
    }
    pub fn read_len(&self) -> usize {
        self.seq.len()
    }

    pub fn cigar_string(&self) -> String {
        if self.cigar.is_empty() {
            return "*".to_string();
        }
        let mut s = String::new();
        for (op, len) in &self.cigar {
            s.push_str(&len.to_string());
            s.push(*op as char);
        }
        s
    }

    pub fn left_soft_clip(&self) -> i32 {
        match self.cigar.first() {
            Some((b'S', l)) => *l as i32,
            _ => 0,
        }
    }
    pub fn right_soft_clip(&self) -> i32 {
        match self.cigar.last() {
            Some((b'S', l)) => *l as i32,
            _ => 0,
        }
    }
    pub fn has_op(&self, op: u8) -> bool {
        self.cigar.iter().any(|c| c.0 == op)
    }

    /// SAMUtils.getAlignmentBlocks: (1-based read start, reference start, length) of each M, = and X
    pub fn alignment_blocks(&self) -> Vec<(i32, i32, i32)> {
        let mut out = Vec::new();
        let mut read_base = 1i32;
        let mut ref_base = self.start;
        for &(op, len) in &self.cigar {
            let len = len as i32;
            match op {
                b'H' | b'P' => {}
                b'S' | b'I' => read_base += len,
                b'N' | b'D' => ref_base += len,
                b'M' | b'=' | b'X' => {
                    out.push((read_base, ref_base, len));
                    read_base += len;
                    ref_base += len;
                }
                _ => {}
            }
        }
        out
    }

    /// SAMRecord.getReferencePositionAtReadPosition (1-based read position; 0 when not aligned)
    pub fn ref_pos_at_read_pos(&self, offset: i32) -> i32 {
        if offset == 0 {
            return 0;
        }
        for (read_start, ref_start, len) in self.alignment_blocks() {
            if read_start + len - 1 < offset {
                continue;
            }
            if offset < read_start {
                return 0;
            }
            return ref_start + offset - read_start;
        }
        0
    }

    /// SAMRecord.getReadPositionAtReferencePosition (1-based; 0 when not aligned, or with `last_if_deleted`
    /// the last aligned read position before it)
    pub fn read_pos_at_ref_pos(&self, pos: i32, last_if_deleted: bool) -> i32 {
        if pos <= 0 {
            return 0;
        }
        let mut last = 0;
        for (read_start, ref_start, len) in self.alignment_blocks() {
            if ref_start + len - 1 >= pos {
                if pos < ref_start {
                    return if last_if_deleted { last } else { 0 };
                }
                return pos - ref_start + read_start;
            }
            last = read_start + len - 1;
        }
        0
    }

    /// CigarUtils.getReadIndexFromPosition(alignmentStart, elements, position, gapMode 0, inferFromSoftClip)
    pub fn read_index_from_position(&self, position: i32, infer_from_soft_clip: bool) -> i32 {
        if position < self.start && !infer_from_soft_clip {
            return -1;
        }
        let mut ref_pos = self.start;
        let mut index = 0;
        for &(op, len) in &self.cigar {
            let len = len as i32;
            if op == b'S' && ref_pos == self.start && infer_from_soft_clip {
                ref_pos -= len;
            }
            let consumes = consumes_ref(op) || (op == b'S' && infer_from_soft_clip);
            if consumes && ref_pos + len - 1 >= position {
                if consumes_read(op) {
                    index += position - ref_pos;
                } else {
                    return -1;
                }
                return index;
            }
            if consumes {
                ref_pos += len;
            }
            if consumes_read(op) {
                index += len;
            }
        }
        -1
    }

    /// SamRecordUtils.generateMappedCoords: reference ranges of the M runs, joined over D and I
    pub fn mapped_coords(&self) -> Vec<(i32, i32)> {
        let mut out: Vec<(i32, i32)> = Vec::new();
        let mut offset = 0;
        let mut cont = false;
        for &(op, len) in &self.cigar {
            let len = len as i32;
            match op {
                b'D' => {
                    offset += len;
                    cont = true;
                }
                b'I' => cont = true,
                b'N' => {
                    offset += len;
                    cont = false;
                }
                b'M' => {
                    let s = self.start + offset;
                    let e = s + len - 1;
                    if cont && !out.is_empty() {
                        out.last_mut().unwrap().1 = e;
                    } else {
                        out.push((s, e));
                    }
                    offset += len;
                    cont = false;
                }
                _ => {}
            }
        }
        out
    }

    /// UltimaBamUtils.extractLowQualIndices
    pub fn ultima_low_qual_indices(&self) -> Vec<i32> {
        let Some(tag) = &self.uq else { return Vec::new() };
        let mut parts = tag.splitn(2, '=');
        let _total = parts.next();
        let Some(rest) = parts.next() else { return Vec::new() };
        let mut out = Vec::new();
        for item in hmf_common::java_text::split(rest, ',') {
            if item.contains('-') {
                let mut se = item.splitn(2, '-');
                let s = hmf_common::java_text::parse_int(se.next().unwrap_or("")).unwrap_or(0);
                let e = hmf_common::java_text::parse_int(se.next().unwrap_or("")).unwrap_or(-1);
                for i in s..=e {
                    out.push(i);
                }
            } else if let Some(v) = hmf_common::java_text::parse_int(item) {
                out.push(v);
            }
        }
        out
    }

    /// UltimaBamUtils.extractLowQualCount
    pub fn ultima_low_qual_count(&self) -> i32 {
        let Some(tag) = &self.uq else { return 0 };
        hmf_common::java_text::parse_int(tag.splitn(2, '=').next().unwrap_or("")).unwrap_or(0)
    }

    /// readInfo: contig:start-end cigar
    pub fn read_info(&self) -> String {
        format!("{}:{}-{} {}", self.contig, self.start, self.end, self.cigar_string())
    }
}

/// BaseQualAdjustment.isLowBaseQual on the Java (signed) byte: quals of 128 and above are negative, so low
pub fn is_low_base_qual(q: u8) -> bool {
    (q as i8 as i32) <= LOW_BASE_QUAL_THRESHOLD as i32
}

pub fn is_medium_base_qual(q: u8, params: &Params) -> bool {
    params.is_sbx() && (SBX_MEDIUM_QUAL_LOWER..=SBX_MEDIUM_QUAL_UPPER).contains(&q)
}

pub fn is_uncertain_base_qual(q: i8) -> bool {
    q <= BASE_QUAL_MINIMUM as i8
}

/// BamSlicerFilter(0, false, false, false).passesFilters
pub fn passes_filters(r: &Rec) -> bool {
    !(r.flag(FLAG_UNMAPPED) || r.flag(FLAG_SECONDARY) || r.flag(FLAG_SUPPLEMENTARY) || r.flag(FLAG_DUPLICATE))
}

/// An indexed BAM or CRAM
pub struct Bam {
    reader: bam::IndexedReader,
    pub path: String,
}

impl Bam {
    pub fn open(path: &str, ref_genome: &str) -> Result<Bam> {
        let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("failed to open BAM {path}"))?;
        if !ref_genome.is_empty() {
            let _ = reader.set_reference(ref_genome);
        }
        Ok(Bam { reader, path: path.to_string() })
    }

    /// BamSlicer.slice: the records overlapping the region that pass the filters, in file order
    pub fn slice(&mut self, chr: &str, start: i32, end: i32) -> Result<Vec<Arc<Rec>>> {
        let header = self.reader.header().clone();
        let Some(tid) = header.tid(chr.as_bytes()) else {
            crate::warn!("cannot find sequence index for chromosome {chr} in bam header");
            return Ok(Vec::new());
        };
        self.reader.fetch((tid as i32, (start - 1).max(0) as i64, end as i64))?;
        let mut out = Vec::new();
        let mut rec = bam::Record::new();
        while let Some(r) = self.reader.read(&mut rec) {
            r?;
            let x = Rec::from_bam(&rec, &header);
            if passes_filters(&x) {
                out.push(Arc::new(x));
            }
        }
        Ok(out)
    }

    /// BamSlicer.queryMate: the first record at the mate's start with the other first-in-pair flag, the same
    /// supplementary flag and the same name
    pub fn query_mate(&mut self, r: &Rec) -> Result<Option<Arc<Rec>>> {
        if r.mate_tid < 0 {
            // queryUnmapped: not supported here (the mates looked up are those of mapped reads near a variant)
            return Ok(None);
        }
        let header = self.reader.header().clone();
        self.reader.fetch((r.mate_tid, (r.mate_start - 1).max(0) as i64, r.mate_start as i64))?;
        let mut rec = bam::Record::new();
        while let Some(res) = self.reader.read(&mut rec) {
            res?;
            if rec.pos() as i32 + 1 != r.mate_start {
                continue;
            }
            let x = Rec::from_bam(&rec, &header);
            if !x.paired() {
                if x.name == r.name {
                    crate::error!("read({}) loc({}:{}) isFirst({}) mate not paired", r.name, r.contig, r.start, r.first_in_pair());
                    return Ok(None);
                }
                continue;
            }
            if x.first_in_pair() == r.first_in_pair() || x.supplementary() != r.supplementary() {
                continue;
            }
            if x.name == r.name {
                return Ok(if passes_filters(&x) { Some(Arc::new(x)) } else { None });
            }
        }
        Ok(None)
    }
}

/// read.Read
#[derive(Clone, Debug)]
pub struct Read {
    pub id: Arc<str>,
    pub soft_clipped_start: i32,
    pub soft_clipped_end: i32,
    pub position_start: i32,
    pub position_end: i32,
    pub read_index_start: i32,
    pub read_index_end: i32,
    pub valid_indels: Vec<Indel>,
    pub ignored_indels: Vec<Indel>,
    pub record: Arc<Rec>,
    pub trimmed_bases: i32,
}

/// The context of read creation: sequencing type and the HLA region's reference bases
pub struct ReadCtx<'a> {
    pub params: &'a Params,
    pub ref_bases: Option<&'a RefBases>,
}

impl Read {
    pub fn contains_soft_clip(&self) -> bool {
        self.soft_clipped_start > 0 || self.soft_clipped_end > 0
    }
    pub fn contains_valid_indel(&self) -> bool {
        !self.valid_indels.is_empty()
    }
    pub fn max_valid_indel_size(&self) -> i32 {
        self.valid_indels.iter().map(|x| x.length.abs()).max().unwrap_or(0)
    }
    pub fn read_info(&self) -> String {
        self.record.read_info()
    }

    /// populateCodingRegion: bases (reverse-complemented for a reverse-strand gene), low-qual flags and quals
    pub fn coding_region(&self, params: &Params, reverse: bool) -> (Vec<u8>, Vec<bool>, Vec<i8>) {
        let n = (self.read_index_end - self.read_index_start + 1).max(0) as usize;
        let mut bases = vec![0u8; n];
        let mut low = vec![false; n];
        let mut quals = vec![0i8; n];
        let ultima: Option<Vec<i32>> = if params.is_ultima() { Some(self.record.ultima_low_qual_indices()) } else { None };
        let r = &self.record;
        for (k, i) in (self.read_index_start..=self.read_index_end).enumerate() {
            let idx = if reverse { n - 1 - k } else { k };
            let b = r.seq[i as usize];
            bases[idx] = if reverse { swap_dna_base(b) } else { b };
            let (q, l) = match &ultima {
                None => {
                    let q = r.quals[i as usize];
                    (q as i8, is_low_base_qual(q))
                }
                Some(u) => (MISSING_BASE_QUAL, u.contains(&i)),
            };
            low[idx] = l;
            quals[idx] = q;
        }
        (bases, low, quals)
    }

    /// Read.createRead
    pub fn create(ctx: &ReadCtx, region: (i32, i32), record: &Arc<Rec>, include_soft_clips: bool, include_indels: bool) -> Result<Read> {
        let r = record.as_ref();
        let sc_left = r.left_soft_clip();
        let sc_right = r.right_soft_clip();
        let alignment_start = r.start;
        let alignment_end = r.end;
        let mut position_start = alignment_start;
        let mut position_end = alignment_end;
        let mut read_index_start;
        let mut read_index_end = r.read_len() as i32 - sc_right - 1;
        read_index_start = sc_left;
        let mut soft_clip_start = sc_left;
        let mut soft_clip_end = sc_right;
        if r.paired() && r.isize.abs() <= r.read_len() as i64 {
            if r.reverse() {
                soft_clip_start = 0;
            } else {
                soft_clip_end = 0;
            }
        }
        let trimmed = if ctx.params.is_illumina() { find_low_qual_trimmed_bases(r) } else { 0 };
        if trimmed > 0 {
            if r.reverse() {
                soft_clip_start = (soft_clip_start - trimmed).max(0);
                read_index_start += trimmed - sc_left;
                position_start = r.ref_pos_at_read_pos(read_index_start + 1);
            } else {
                soft_clip_end = (soft_clip_end - trimmed).max(0);
                read_index_end -= trimmed - sc_right;
                position_end = r.ref_pos_at_read_pos(read_index_end + 1);
            }
        }
        let unclipped_start = if soft_clip_start > 0 { (alignment_start - soft_clip_start).max(region.0) } else { 0 };
        let unclipped_end = if soft_clip_end > 0 { (alignment_end + soft_clip_end).min(region.1) } else { 0 };
        if include_soft_clips {
            if soft_clip_start > 0 {
                position_start = unclipped_start;
            }
            if soft_clip_end > 0 {
                position_end = unclipped_end;
            }
        }
        position_start = position_start.max(region.0);
        position_end = position_end.min(region.1);
        read_index_start = r.read_index_from_position(position_start, true);
        read_index_end = r.read_index_from_position(position_end, true);
        if r.has_op(b'N') {
            for (s, e) in r.mapped_coords() {
                if e < region.0 {
                    continue;
                }
                if s > region.1 {
                    break;
                }
                position_start = position_start.max(s).max(region.0);
                position_end = position_end.min(e).min(region.1);
                read_index_start = r.read_pos_at_ref_pos(position_start, true) - 1;
                read_index_end = r.read_pos_at_ref_pos(position_end, true) - 1;
                break;
            }
        }
        let soft_clipped_start = (alignment_start - position_start).max(0);
        let soft_clipped_end = (position_end - alignment_end).max(0);
        let indels = if include_indels { read_indels(ctx, position_start, position_end, r)? } else { Vec::new() };
        let (ignored, valid): (Vec<Indel>, Vec<Indel>) = indels.into_iter().partition(|x| x.ignore);
        Ok(Read {
            id: Arc::from(r.name.as_str()),
            soft_clipped_start,
            soft_clipped_end,
            position_start,
            position_end,
            read_index_start,
            read_index_end,
            valid_indels: valid,
            ignored_indels: ignored,
            record: record.clone(),
            trimmed_bases: trimmed,
        })
    }

    /// alignmentsOnly: one Read per alignment block overlapping the read's coding range. The block ends are
    /// taken one past the block (refStart + length), as in Java; the callers only count them.
    pub fn alignments_only(&self, ctx: &ReadCtx) -> Result<Vec<Read>> {
        let mut out = Vec::new();
        for (_, ref_start, len) in self.record.alignment_blocks() {
            let (s, e) = (ref_start, ref_start + len);
            if s <= self.position_end && e >= self.position_start {
                let region = (self.position_start.max(s), self.position_end.min(e));
                out.push(Read::create(ctx, region, &self.record, false, false)?);
            }
        }
        Ok(out)
    }
}

pub fn swap_dna_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        x => x,
    }
}

/// Read.findLowQualTrimmedBases: the bases at the read's 3' end to drop for low base quality
fn find_low_qual_trimmed_bases(r: &Rec) -> i32 {
    let from_start = r.reverse();
    let n = r.read_len() as i32;
    let mut idx: i32 = if from_start { 0 } else { n - 1 };
    let low_qual_score = 1.0 / LOW_BASE_TRIM_PERC - 1.0;
    let mut lowest = 0.0f64;
    let mut current = 0.0f64;
    let mut last_lowest: i32 = -1;
    while idx >= 0 && idx < n {
        if is_low_base_qual(r.quals[idx as usize]) {
            current -= low_qual_score;
            if current <= lowest {
                lowest = current;
                last_lowest = idx;
            }
        } else {
            current += 1.0;
        }
        if from_start {
            idx += 1;
        } else {
            idx -= 1;
        }
    }
    if last_lowest < 0 {
        return 0;
    }
    if from_start { last_lowest + 1 } else { n - last_lowest }
}

/// Read.indels: the inserts and deletes whose reference position lies in [start, end], via CigarHandler
fn read_indels(ctx: &ReadCtx, start: i32, end: i32, r: &Rec) -> Result<Vec<Indel>> {
    let mut out = Vec::new();
    let params = ctx.params;
    let ultima: Option<Vec<i32>> = if params.is_ultima() { Some(r.ultima_low_qual_indices()) } else { None };
    let low_or_medium = |i: i32| -> bool {
        match &ultima {
            Some(u) => u.contains(&i),
            None => {
                let q = r.quals[i as usize];
                is_low_base_qual(q) || is_medium_base_qual(q, params)
            }
        }
    };
    let n_ops = r.cigar.len();
    let mut read_index: i32 = 0;
    let mut ref_base: i32 = r.start;
    for (k, &(op, len)) in r.cigar.iter().enumerate() {
        let len = len as i32;
        match op {
            b'H' | b'P' => {}
            b'S' => {
                let _ = (k, n_ops);
                read_index += len;
            }
            b'N' => ref_base += len,
            b'D' => {
                let ri = read_index - 1;
                let rp = ref_base - 1;
                if start <= rp && rp <= end {
                    let base = *r.seq.get(ri as usize).context("java.lang.StringIndexOutOfBoundsException")? as char;
                    let del_start = rp + 1;
                    let del_end = del_start + len - 1;
                    let refb = ctx.ref_bases.context("java.lang.NullPointerException (no reference bases)")?.get(del_start, del_end)?;
                    let del_bases = format!("{base}{refb}");
                    let mut ignore = false;
                    if params.is_ultima() || params.is_sbx() {
                        ignore = low_or_medium(ri);
                        if !ignore {
                            ignore = low_or_medium(ri + 1);
                        }
                    }
                    out.push(Indel::new(&r.contig, rp, ri, &del_bases, &base.to_string(), ignore));
                }
                ref_base += len;
            }
            b'I' => {
                let (ri, rp) = if read_index == 0 { (0, ref_base - 1) } else { (read_index - 1, ref_base - 1) };
                if start <= rp && rp <= end {
                    let base = r.seq[ri as usize] as char;
                    let mut ignore = false;
                    if params.is_ultima() || params.is_sbx() {
                        for i in ri + 1..=ri + len {
                            ignore = low_or_medium(i);
                            if ignore {
                                break;
                            }
                        }
                    }
                    let alt = String::from_utf8_lossy(&r.seq[ri as usize..(ri + len + 1) as usize]).into_owned();
                    out.push(Indel::new(&r.contig, rp, ri, &base.to_string(), &alt, ignore));
                }
                read_index += len;
            }
            b'M' | b'=' | b'X' => {
                read_index += len;
                ref_base += len;
            }
            _ => {}
        }
    }
    Ok(out)
}

/// BamRecordReader.filterLowQualRead: half or more of the bases low-qual
pub fn filter_low_qual_read(r: &Rec, params: &Params) -> bool {
    let n = r.read_len();
    let threshold = (n as f64 * MAX_LOW_BASE_PERC).floor() as i32 + 1;
    if params.is_ultima() {
        return r.ultima_low_qual_count() >= threshold;
    }
    let mut low = 0;
    for i in 0..n {
        if is_low_base_qual(r.quals[i]) {
            low += 1;
            if low >= threshold {
                return true;
            }
        } else {
            let high = i as i32 + 1 - low;
            if high >= threshold {
                return false;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    pub fn rec(start: i32, cigar: &[(u8, u32)], seq: &[u8]) -> Rec {
        let ref_len: i32 = cigar.iter().filter(|c| consumes_ref(c.0)).map(|c| c.1 as i32).sum();
        Rec {
            name: "R1".into(),
            flags: 0,
            contig: "chr6".into(),
            start,
            end: start + ref_len - 1,
            cigar: cigar.to_vec(),
            seq: seq.to_vec(),
            quals: vec![37; seq.len()],
            mapq: 60,
            mate_contig: "*".into(),
            mate_tid: -1,
            mate_start: 0,
            isize: 0,
            uq: None,
        }
    }

    #[test]
    fn positions() {
        // 2S 5M 2I 3M 1D 4M
        let r = rec(100, &[(b'S', 2), (b'M', 5), (b'I', 2), (b'M', 3), (b'D', 1), (b'M', 4)], b"AAACCCGGGTTTACGT");
        assert_eq!(r.end, 100 + 5 + 3 + 1 + 4 - 1);
        assert_eq!(r.ref_pos_at_read_pos(1), 0);
        assert_eq!(r.ref_pos_at_read_pos(3), 100);
        assert_eq!(r.ref_pos_at_read_pos(8), 0);
        assert_eq!(r.ref_pos_at_read_pos(10), 105);
        assert_eq!(r.read_pos_at_ref_pos(108, false), 0);
        assert_eq!(r.read_pos_at_ref_pos(108, true), 12);
        assert_eq!(r.read_pos_at_ref_pos(109, false), 13);
        assert_eq!(r.read_index_from_position(98, true), 0);
        assert_eq!(r.read_index_from_position(100, true), 2);
        assert_eq!(r.read_index_from_position(108, true), -1);
        assert_eq!(r.mapped_coords(), vec![(100, 112)]);
    }
}
