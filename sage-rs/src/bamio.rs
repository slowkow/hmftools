//! BAM slicing: the equivalent of SamSlicer / BamSlicer over an indexed BAM.

use crate::types::*;
use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::htslib;
use std::cell::RefCell;

pub struct BamSlicer {
    reader: bam::IndexedReader,
    record: bam::Record,
    /// reused by `slice_ref`
    scratch: Read,
    /// the records of the last `slice_ref_caching` call
    cache: ReadCache,
}

/// A BAM record as stored in the file: the fixed fields and the variable-length data (name, cigar,
/// packed bases, quals, tags).
#[derive(Clone, Copy)]
struct RawRecord<'a> {
    core: RawCore,
    data: &'a [u8],
}

#[derive(Clone, Copy)]
struct RawCore {
    tid: i32,
    /// 0-based
    pos: i32,
    flags: u16,
    mapq: u8,
    /// without the terminating and padding NULs
    name_len: u16,
    /// offset of the cigar in the data
    cigar_offset: u16,
    n_cigar: u32,
    seq_len: u32,
    mate_tid: i32,
    mate_pos: i32,
    insert_size: i32,
}

impl<'a> RawRecord<'a> {
    fn of(rec: &'a bam::Record) -> Self {
        let b = rec.inner();
        let c = &b.core;
        // SAFETY: htslib guarantees `data` holds `l_data` bytes, and the slice borrows the record
        let data: &[u8] = if b.l_data > 0 { unsafe { std::slice::from_raw_parts(b.data, b.l_data as usize) } } else { &[] };
        RawRecord {
            core: RawCore {
                tid: c.tid,
                pos: c.pos as i32,
                flags: c.flag,
                mapq: c.qual,
                name_len: (c.l_qname as u16).saturating_sub(1 + c.l_extranul as u16),
                cigar_offset: c.l_qname as u16,
                n_cigar: c.n_cigar,
                seq_len: c.l_qseq as u32,
                mate_tid: c.mtid,
                mate_pos: c.mpos as i32,
                insert_size: c.isize_ as i32,
            },
            data,
        }
    }
    fn name(&self) -> &'a [u8] {
        &self.data[..self.core.name_len as usize]
    }
    fn cigar(&self) -> impl Iterator<Item = u32> + 'a {
        let start = self.core.cigar_offset as usize;
        self.data[start..start + self.core.n_cigar as usize * 4]
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
    }
    fn seq_offset(&self) -> usize {
        self.core.cigar_offset as usize + self.core.n_cigar as usize * 4
    }
    fn packed_bases(&self) -> &'a [u8] {
        let start = self.seq_offset();
        &self.data[start..start + (self.core.seq_len as usize).div_ceil(2)]
    }
    fn quals(&self) -> &'a [u8] {
        let start = self.seq_offset() + (self.core.seq_len as usize).div_ceil(2);
        &self.data[start..start + self.core.seq_len as usize]
    }
    fn aux(&self) -> &'a [u8] {
        let start = self.seq_offset() + (self.core.seq_len as usize).div_ceil(2) + self.core.seq_len as usize;
        self.data.get(start..).unwrap_or(&[])
    }
    /// htslib's bam_endpos: the 0-based exclusive end on the reference, one past the start for a
    /// record that covers no reference bases
    fn end_pos(&self) -> i32 {
        let mut ref_len = 0;
        for c in self.cigar() {
            // M, D, N, = and X consume the reference
            if matches!(c & 0xf, 0 | 2 | 3 | 7 | 8) {
                ref_len += (c >> 4) as i32;
            }
        }
        self.core.pos + ref_len.max(1)
    }
}

/// Whether an indexed query of the 1-based closed `region` returns a record at `pos` (0-based) with
/// exclusive end `end_pos`.
fn query_returns(pos: i32, end_pos: i32, region: &BaseRegion) -> bool {
    pos < region.end && end_pos > region.start - 1
}

/// Whether a record is returned by the queries of `regions`, which are sorted and do not overlap.
/// `next` is the first region that may still overlap; records arrive in coordinate order.
fn in_regions(pos: i32, end_pos: i32, regions: &[BaseRegion], next: &mut usize) -> bool {
    while *next < regions.len() && regions[*next].end <= pos {
        *next += 1;
    }
    *next < regions.len() && query_returns(pos, end_pos, &regions[*next])
}

fn sorted_and_disjoint(regions: &[BaseRegion]) -> bool {
    regions.iter().all(|r| r.start <= r.end) && regions.windows(2).all(|w| w[0].end < w[1].start)
}

/// The most record data (MB) kept per sample for one region partition, by default; a partition
/// with more is read from the file a second time instead.
pub const DEFAULT_READ_CACHE_MB: i32 = 32;

/// The size (KB) of the buffer each alignment file handle reads through (-file_buffer_kb; 0 leaves
/// htslib's default).
pub const DEFAULT_FILE_BUFFER_KB: i32 = 1024;
pub static FILE_BUFFER_KB: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(DEFAULT_FILE_BUFFER_KB);
const READ_CACHE_GROWTH_BYTES: usize = 4 << 20;

/// The records one sample has within a set of regions, kept as stored in the BAM, so that the
/// evidence stage can go over the reads the candidate stage has just read without reading and
/// decompressing them from the file again.
#[derive(Default)]
struct ReadCache {
    tid: i32,
    /// what was queried: every record these queries return is held, in file order
    regions: Vec<BaseRegion>,
    entries: Vec<(RawCore, usize, i32)>,
    data: Vec<u8>,
    usable: bool,
    max_bytes: usize,
}

impl ReadCache {
    fn reset(&mut self, tid: i32, regions: &[BaseRegion]) {
        self.tid = tid;
        self.regions.clear();
        self.regions.extend_from_slice(regions);
        // after an unusually deep partition, go back to about what the last one needed
        let recent = self.data.len() + READ_CACHE_GROWTH_BYTES;
        if self.data.capacity() > recent * 2 {
            self.data = Vec::with_capacity(recent);
            self.entries = Vec::new();
        }
        self.entries.clear();
        self.data.clear();
        self.usable = self.max_bytes > 0;
    }
    fn invalidate(&mut self) {
        self.usable = false;
        // a partition this deep is rare: give the memory back
        self.entries = Vec::new();
        self.data = Vec::new();
    }
    fn push(&mut self, raw: &RawRecord, end_pos: i32) {
        if !self.usable {
            return;
        }
        let needed = self.data.len() + raw.data.len();
        if needed > self.max_bytes {
            self.invalidate();
            return;
        }
        if needed > self.data.capacity() {
            // grow in steps rather than by doubling, which would hold up to twice the limit
            self.data.reserve_exact((needed - self.data.len()).max(READ_CACHE_GROWTH_BYTES));
        }
        self.entries.push((raw.core, self.data.len(), end_pos));
        // everything up to the tags, then only the tags that are read
        let aux = raw.aux();
        self.data.extend_from_slice(&raw.data[..raw.data.len() - aux.len()]);
        let data = &mut self.data;
        let vis = crate::vis::is_enabled();
        for_each_aux_span(aux, |tag, _, span| {
            if matches!(tag, b"NM" | b"MC" | b"YX" | b"UT" | b"tp" | b"t0" | b"UQ") || (vis && tag == b"CR") {
                data.extend_from_slice(span);
            }
        });
    }
    /// Whether every record the queries of `regions` return is held.
    fn covers(&self, tid: i32, regions: &[BaseRegion]) -> bool {
        self.usable
            && self.tid == tid
            && sorted_and_disjoint(regions)
            && regions.iter().all(|r| self.regions.iter().any(|c| c.start <= r.start && r.end <= c.end))
    }
    fn record(&self, i: usize) -> (RawRecord<'_>, i32) {
        let (core, offset, end_pos) = self.entries[i];
        let end = self.entries.get(i + 1).map_or(self.data.len(), |e| e.1);
        (RawRecord { core, data: &self.data[offset..end] }, end_pos)
    }
}

/// The 16 four-bit base codes of the BAM format, as htslib's seq_nt16_str.
const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Both bases of every packed sequence byte.
static BASE_PAIRS: [[u8; 2]; 256] = {
    let mut t = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = [NT16[i >> 4], NT16[i & 15]];
        i += 1;
    }
    t
};

fn decode_bases(encoded: &[u8], len: usize, out: &mut Vec<u8>) {
    out.clear();
    out.resize(len.div_ceil(2) * 2, 0);
    for (pair, b) in out.chunks_exact_mut(2).zip(encoded) {
        pair.copy_from_slice(&BASE_PAIRS[*b as usize]);
    }
    out.truncate(len);
}

/// The values of the tags SAGE reads; other types are skipped.
pub(crate) enum AuxValue<'a> {
    Int(i32),
    Str(&'a [u8]),
    /// B:c or B:C
    ByteArray(&'a [u8]),
    Other,
}

/// Walks the tags of a BAM record without validating or allocating; stops at a malformed tag.
pub(crate) fn for_each_aux<'a, F: FnMut(&'a [u8], AuxValue<'a>)>(aux: &'a [u8], mut f: F) {
    for_each_aux_span(aux, |tag, value, _| f(tag, value))
}

/// As `for_each_aux`, also giving the tag as stored (name, type and value).
fn for_each_aux_span<'a, F: FnMut(&'a [u8], AuxValue<'a>, &'a [u8])>(mut aux: &'a [u8], mut f: F) {
    while aux.len() >= 3 {
        let tag = &aux[..2];
        let ty = aux[2];
        let body = &aux[3..];
        let fixed = |n: usize| if body.len() >= n { Some(n) } else { None };
        let (value, used) = match ty {
            b'c' => match fixed(1) {
                Some(n) => (AuxValue::Int(body[0] as i8 as i32), n),
                None => return,
            },
            b'C' => match fixed(1) {
                Some(n) => (AuxValue::Int(body[0] as i32), n),
                None => return,
            },
            b's' => match fixed(2) {
                Some(n) => (AuxValue::Int(i16::from_le_bytes([body[0], body[1]]) as i32), n),
                None => return,
            },
            b'S' => match fixed(2) {
                Some(n) => (AuxValue::Int(u16::from_le_bytes([body[0], body[1]]) as i32), n),
                None => return,
            },
            b'i' => match fixed(4) {
                Some(n) => (AuxValue::Int(i32::from_le_bytes([body[0], body[1], body[2], body[3]])), n),
                None => return,
            },
            b'I' => match fixed(4) {
                Some(n) => (AuxValue::Int(u32::from_le_bytes([body[0], body[1], body[2], body[3]]) as i32), n),
                None => return,
            },
            b'A' => match fixed(1) {
                Some(n) => (AuxValue::Other, n),
                None => return,
            },
            b'f' => match fixed(4) {
                Some(n) => (AuxValue::Other, n),
                None => return,
            },
            b'd' => match fixed(8) {
                Some(n) => (AuxValue::Other, n),
                None => return,
            },
            b'Z' | b'H' => match memchr::memchr(0, body) {
                Some(end) => (if ty == b'Z' { AuxValue::Str(&body[..end]) } else { AuxValue::Other }, end + 1),
                None => return,
            },
            b'B' => {
                if body.len() < 5 {
                    return;
                }
                let width = match body[0] {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return,
                };
                let count = u32::from_le_bytes([body[1], body[2], body[3], body[4]]) as usize;
                let end = 5 + count * width;
                if body.len() < end {
                    return;
                }
                (if width == 1 { AuxValue::ByteArray(&body[5..end]) } else { AuxValue::Other }, end)
            }
            _ => return,
        };
        f(tag, value, &aux[..3 + used]);
        aux = &body[used..];
    }
}

pub(crate) fn blank_read() -> Read {
    Read {
        name: String::new(),
        flags: 0,
        tid: -1,
        start: 0,
        end: 0,
        mapq: 0,
        cigar: Vec::new(),
        bases: Vec::new(),
        quals: Vec::new(),
        mate_tid: -1,
        mate_start: 0,
        insert_size: 0,
        nm: None,
        mate_cigar: None,
        consensus: ConsensusType::None,
        has_consensus_tag: false,
        duplex_index: -1,
        tp: Vec::new(),
        t0: Vec::new(),
        low_qual_tag: None,
        consensus_read_tag: None,
    }
}

/// Fill `read` from a BAM record, reusing its buffers.
fn fill(rec: &RawRecord, read: &mut Read) {
    read.cigar.clear();
    let mut ref_len = 0;
    for c in rec.cigar() {
        let len = (c >> 4) as i32;
        let op = match c & 0xf {
            0 => Op::M,
            1 => Op::I,
            2 => Op::D,
            3 => Op::N,
            4 => Op::S,
            5 => Op::H,
            6 => Op::P,
            7 => Op::Eq,
            _ => Op::X,
        };
        if op.consumes_ref() {
            ref_len += len;
        }
        read.cigar.push(CigarEl::new(len, op));
    }
    let start = rec.core.pos + 1;
    // one pass over the tags: NM, MC and UT (SamRecordUtils.extractConsensusType; UT is set by REDUX)
    read.nm = None;
    read.consensus = ConsensusType::None;
    read.has_consensus_tag = false;
    read.duplex_index = -1;
    read.tp.clear();
    read.t0.clear();
    read.low_qual_tag = None;
    read.consensus_read_tag = None;
    let ultima = crate::constants::is_ultima();
    let vis = crate::vis::is_enabled();
    let mut has_mate_cigar = false;
    for_each_aux(rec.aux(), |tag, value| match (tag, value) {
        (b"NM", AuxValue::Int(v)) => read.nm = Some(v),
        (b"MC", AuxValue::Str(s)) => {
            if let Ok(s) = std::str::from_utf8(s) {
                has_mate_cigar = true;
                match read.mate_cigar.as_mut() {
                    Some(mc) => {
                        mc.clear();
                        mc.push_str(s);
                    }
                    None => read.mate_cigar = Some(s.to_string()),
                }
            }
        }
        (b"YX", AuxValue::Int(v)) => read.duplex_index = v,
        // UltimaBamUtils: tp (B:c), t0 (a quality string) and UQ (low-qual indices, set by REDUX)
        (b"tp", AuxValue::ByteArray(a)) if ultima => read.tp.extend(a.iter().map(|v| *v as i8)),
        (b"t0", AuxValue::Str(s)) if ultima => read.t0.extend(s.iter().map(|b| (*b as i32 - 33) as i8)),
        (b"UQ", AuxValue::Str(s)) if ultima => read.low_qual_tag = Some(String::from_utf8_lossy(s).into_owned()),
        (b"CR", AuxValue::Str(s)) if vis => read.consensus_read_tag = Some(String::from_utf8_lossy(s).into_owned()),
        (b"UT", AuxValue::Str(s)) => {
            read.has_consensus_tag = true;
            read.consensus = match s {
                b"DUAL" => ConsensusType::Dual,
                b"SINGLE" => ConsensusType::Single,
                _ => ConsensusType::None,
            };
        }
        _ => {}
    });
    if !has_mate_cigar {
        read.mate_cigar = None;
    }
    read.name.clear();
    let name = rec.name();
    match std::str::from_utf8(name) {
        Ok(name) => read.name.push_str(name),
        Err(_) => read.name.push_str(&String::from_utf8_lossy(name)),
    }
    read.flags = rec.core.flags;
    read.tid = rec.core.tid;
    read.start = start;
    read.end = start + ref_len - 1;
    read.mapq = rec.core.mapq as i32;
    decode_bases(rec.packed_bases(), rec.core.seq_len as usize, &mut read.bases);
    read.quals.clear();
    read.quals.extend_from_slice(rec.quals());
    read.mate_tid = rec.core.mate_tid;
    read.mate_start = rec.core.mate_pos + 1;
    read.insert_size = rec.core.insert_size;
}

/// Reads whose buffers can be reused.
#[derive(Default)]
pub struct ReadPool(RefCell<Vec<Read>>);

impl ReadPool {
    pub fn take(&self) -> Read {
        self.0.borrow_mut().pop().unwrap_or_else(blank_read)
    }
    pub fn put(&self, read: Read) {
        self.0.borrow_mut().push(read);
    }
}

fn is_cram(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".cram")
}

impl BamSlicer {
    /// Opens an indexed BAM or CRAM; `ref_genome` is the FASTA a CRAM is decoded against.
    pub fn open(path: &str, ref_genome: &str, threads: usize) -> Result<Self> {
        Self::open_with_cache(path, ref_genome, threads, DEFAULT_READ_CACHE_MB as usize)
    }

    /// As `open`, with the size limit of the read cache in MB (0 for none).
    pub fn open_with_cache(path: &str, ref_genome: &str, threads: usize, read_cache_mb: usize) -> Result<Self> {
        let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("opening BAM {path}"))?;
        if is_cram(path) {
            reader.set_reference(ref_genome).with_context(|| format!("setting CRAM reference {ref_genome}"))?;
            // htslib by default regenerates the MD and NM tags a CRAM leaves out; htsjdk does not, and
            // Java SAGE then counts the mismatches itself (NumberEvents.getOrCalcNm), as is done here
            // SAFETY: the option takes one int argument
            let ret = unsafe {
                htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_CRAM_OPT_DECODE_MD, 0 as std::os::raw::c_int)
            };
            if ret != 0 {
                anyhow::bail!("setting CRAM options for {path}");
            }
        }
        let buffer_kb = FILE_BUFFER_KB.load(std::sync::atomic::Ordering::Relaxed);
        if buffer_kb > 0 {
            // htslib reads the file in blocks of 32 KB or the file system's preferred size; with many
            // threads each reading its own part of a large file, bigger blocks mean far fewer requests
            // SAFETY: the option takes one int argument
            let ret = unsafe {
                htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_HTS_OPT_BLOCK_SIZE, (buffer_kb << 10) as std::os::raw::c_int)
            };
            if ret != 0 {
                anyhow::bail!("setting the file buffer size for {path}");
            }
        }
        if threads > 1 {
            reader.set_threads(threads)?;
        }
        Ok(BamSlicer { reader, record: bam::Record::new(), scratch: blank_read(), cache: ReadCache { max_bytes: read_cache_mb << 20, ..ReadCache::default() } })
    }

    pub fn sequence_names(&self) -> Vec<String> {
        let h = self.reader.header();
        (0..h.target_count()).map(|i| String::from_utf8_lossy(h.tid2name(i)).into_owned()).collect()
    }

    /// Reads overlapping any of the sorted regions, each returned once, in coordinate order.
    /// Filters as BamSlicerFilter(minMapQual=0, keepDuplicates=false, keepSupplementaries, keepSecondaries=false).
    /// The reads are taken from the cache where the last `slice_ref_caching` call holds them all;
    /// each read's buffers come from `pool` where it has any, and the consumer returns the reads it
    /// is finished with to the pool, so a steady state allocates nothing per read.
    pub fn slice_pooled<F: FnMut(Read)>(
        &mut self, chr: &str, regions: &[BaseRegion], keep_supplementaries: bool, pool: &ReadPool, mut consumer: F,
    ) -> Result<()> {
        let tid = self.tid(chr)?;
        if self.cache.covers(tid as i32, regions) {
            let mut next = 0;
            for i in 0..self.cache.entries.len() {
                let (raw, end_pos) = self.cache.record(i);
                if raw.core.flags & FLAG_SUPPLEMENTARY != 0 && !keep_supplementaries {
                    continue;
                }
                if in_regions(raw.core.pos, end_pos, regions, &mut next) {
                    let mut read = pool.take();
                    fill(&raw, &mut read);
                    consumer(read);
                }
            }
            return Ok(());
        }
        self.slice_impl(tid, regions, keep_supplementaries, |rec, _| {
            let mut read = pool.take();
            fill(&RawRecord::of(rec), &mut read);
            consumer(read);
        })
    }

    /// As `slice_pooled`, but lends each read to the consumer from one reused buffer.
    pub fn slice_ref<F: FnMut(&Read)>(
        &mut self, chr: &str, regions: &[BaseRegion], keep_supplementaries: bool, mut consumer: F,
    ) -> Result<()> {
        let tid = self.tid(chr)?;
        self.slice_impl(tid, regions, keep_supplementaries, |rec, scratch| {
            fill(&RawRecord::of(rec), scratch);
            consumer(scratch);
        })
    }

    /// As `slice_ref` keeping supplementary reads, and also keeps every record within `margin` bases
    /// of the regions, so that a following `slice_pooled` of regions inside that span does not read
    /// the file. The consumer is given exactly the reads `slice_ref` would give it.
    pub fn slice_ref_caching<F: FnMut(&Read)>(
        &mut self, chr: &str, regions: &[BaseRegion], margin: i32, mut consumer: F,
    ) -> Result<()> {
        let tid = self.tid(chr)?;
        if !sorted_and_disjoint(regions) {
            self.cache.invalidate();
            return self.slice_ref(chr, regions, true, consumer);
        }
        // widen each region, merging those that then meet
        let mut widened: Vec<BaseRegion> = Vec::with_capacity(regions.len());
        for r in regions {
            let (start, end) = ((r.start - margin).max(1), r.end.saturating_add(margin));
            match widened.last_mut() {
                Some(last) if start <= last.end + 1 => last.end = end,
                _ => widened.push(BaseRegion::new(start, end)),
            }
        }
        let BamSlicer { reader, record, scratch, cache } = self;
        cache.reset(tid as i32, &widened);
        let mut next = 0;
        slice_records(reader, record, tid, &widened, true, |rec| {
            let raw = RawRecord::of(rec);
            let end_pos = raw.end_pos();
            cache.push(&raw, end_pos);
            if in_regions(raw.core.pos, end_pos, regions, &mut next) {
                fill(&raw, scratch);
                consumer(scratch);
            }
        })
    }

    fn tid(&self, chr: &str) -> Result<u32> {
        self.reader.header().tid(chr.as_bytes()).with_context(|| format!("unknown contig {chr}"))
    }

    fn slice_impl<F: FnMut(&bam::Record, &mut Read)>(
        &mut self, tid: u32, regions: &[BaseRegion], keep_supplementaries: bool, mut consumer: F,
    ) -> Result<()> {
        let BamSlicer { reader, record, scratch, .. } = self;
        slice_records(reader, record, tid, regions, keep_supplementaries, |rec| consumer(rec, scratch))
    }
}

/// The records an indexed query of each region returns, each once, in file order, without the
/// unmapped, secondary and duplicate ones.
fn slice_records<F: FnMut(&bam::Record)>(
    reader: &mut bam::IndexedReader, record: &mut bam::Record, tid: u32, regions: &[BaseRegion],
    keep_supplementaries: bool, mut consumer: F,
) -> Result<()> {
    let mut prev_end = 0;
    for region in regions {
        reader.fetch((tid, (region.start - 1) as i64, region.end as i64))?;
        while let Some(r) = reader.read(record) {
            r?;
            let flags = record.flags();
            if flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_DUPLICATE) != 0 {
                continue;
            }
            if flags & FLAG_SUPPLEMENTARY != 0 && !keep_supplementaries {
                continue;
            }
            // already returned by an earlier interval
            if prev_end > 0 && record.pos() as i32 + 1 <= prev_end {
                continue;
            }
            consumer(record);
        }
        prev_end = region.end;
    }
    Ok(())
}
