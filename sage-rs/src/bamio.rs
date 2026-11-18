//! BAM slicing: the equivalent of SamSlicer / BamSlicer over an indexed BAM.

use crate::types::*;
use anyhow::{Context, Result};
use rust_htslib::bam::{self, record::Aux, Read as BamRead};

pub struct BamSlicer {
    reader: bam::IndexedReader,
    record: bam::Record,
    /// reused by `slice_ref`
    scratch: Read,
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
    out.reserve(len + 1);
    for b in &encoded[..len.div_ceil(2)] {
        out.extend_from_slice(&BASE_PAIRS[*b as usize]);
    }
    out.truncate(len);
}

fn blank_read() -> Read {
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
    }
}

/// Fill `read` from a BAM record, reusing its buffers.
fn fill(rec: &bam::Record, read: &mut Read) {
    read.cigar.clear();
    let mut ref_len = 0;
    for c in rec.raw_cigar() {
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
    let start = rec.pos() as i32 + 1;
    // one pass over the tags: NM, MC and UT (SamRecordUtils.extractConsensusType; UT is set by REDUX)
    read.nm = None;
    read.consensus = ConsensusType::None;
    read.has_consensus_tag = false;
    read.duplex_index = -1;
    read.tp.clear();
    read.t0.clear();
    read.low_qual_tag = None;
    let ultima = crate::constants::is_ultima();
    let mut has_mate_cigar = false;
    for (tag, value) in rec.aux_iter().flatten() {
        match tag {
            b"NM" => {
                read.nm = match value {
                    Aux::U8(v) => Some(v as i32),
                    Aux::I8(v) => Some(v as i32),
                    Aux::U16(v) => Some(v as i32),
                    Aux::I16(v) => Some(v as i32),
                    Aux::U32(v) => Some(v as i32),
                    Aux::I32(v) => Some(v),
                    _ => None,
                }
            }
            b"MC" => {
                if let Aux::String(s) = value {
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
            b"YX" => {
                read.duplex_index = match value {
                    Aux::U8(v) => v as i32,
                    Aux::I8(v) => v as i32,
                    Aux::U16(v) => v as i32,
                    Aux::I16(v) => v as i32,
                    Aux::U32(v) => v as i32,
                    Aux::I32(v) => v,
                    _ => -1,
                }
            }
            // UltimaBamUtils: tp (B:c), t0 (a quality string) and UQ (low-qual indices, set by REDUX)
            b"tp" if ultima => match value {
                Aux::ArrayI8(a) => read.tp.extend(a.iter()),
                Aux::ArrayU8(a) => read.tp.extend(a.iter().map(|v| v as i8)),
                _ => {}
            },
            b"t0" if ultima => {
                if let Aux::String(s) = value {
                    read.t0.extend(s.bytes().map(|b| (b as i32 - 33) as i8));
                }
            }
            b"UQ" if ultima => {
                if let Aux::String(s) = value {
                    read.low_qual_tag = Some(s.to_string());
                }
            }
            b"UT" => {
                if let Aux::String(s) = value {
                    read.has_consensus_tag = true;
                    read.consensus = match s {
                        "DUAL" => ConsensusType::Dual,
                        "SINGLE" => ConsensusType::Single,
                        _ => ConsensusType::None,
                    };
                }
            }
            _ => {}
        }
    }
    if !has_mate_cigar {
        read.mate_cigar = None;
    }
    read.name.clear();
    match std::str::from_utf8(rec.qname()) {
        Ok(name) => read.name.push_str(name),
        Err(_) => read.name.push_str(&String::from_utf8_lossy(rec.qname())),
    }
    read.flags = rec.flags();
    read.tid = rec.tid();
    read.start = start;
    read.end = start + ref_len - 1;
    read.mapq = rec.mapq() as i32;
    let seq = rec.seq();
    decode_bases(seq.encoded, seq.len(), &mut read.bases);
    read.quals.clear();
    read.quals.extend_from_slice(rec.qual());
    read.mate_tid = rec.mtid();
    read.mate_start = rec.mpos() as i32 + 1;
    read.insert_size = rec.insert_size() as i32;
}

impl BamSlicer {
    pub fn open(path: &str, threads: usize) -> Result<Self> {
        let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("opening BAM {path}"))?;
        if threads > 1 {
            reader.set_threads(threads)?;
        }
        Ok(BamSlicer { reader, record: bam::Record::new(), scratch: blank_read() })
    }

    pub fn sequence_names(&self) -> Vec<String> {
        let h = self.reader.header();
        (0..h.target_count()).map(|i| String::from_utf8_lossy(h.tid2name(i)).into_owned()).collect()
    }

    /// Reads overlapping any of the sorted regions, each returned once, in coordinate order.
    /// Filters as BamSlicerFilter(minMapQual=0, keepDuplicates=false, keepSupplementaries, keepSecondaries=false).
    pub fn slice<F: FnMut(Read)>(
        &mut self, chr: &str, regions: &[BaseRegion], keep_supplementaries: bool, mut consumer: F,
    ) -> Result<()> {
        self.slice_impl(chr, regions, keep_supplementaries, |rec, _| {
            let mut read = blank_read();
            fill(rec, &mut read);
            consumer(read);
        })
    }

    /// As `slice`, but lends each read to the consumer from one reused buffer, so nothing is
    /// allocated per read.
    pub fn slice_ref<F: FnMut(&Read)>(
        &mut self, chr: &str, regions: &[BaseRegion], keep_supplementaries: bool, mut consumer: F,
    ) -> Result<()> {
        self.slice_impl(chr, regions, keep_supplementaries, |rec, scratch| {
            fill(rec, scratch);
            consumer(scratch);
        })
    }

    fn slice_impl<F: FnMut(&bam::Record, &mut Read)>(
        &mut self, chr: &str, regions: &[BaseRegion], keep_supplementaries: bool, mut consumer: F,
    ) -> Result<()> {
        let tid = self.reader.header().tid(chr.as_bytes()).with_context(|| format!("unknown contig {chr}"))?;
        let mut prev_end = 0;
        for region in regions {
            self.reader.fetch((tid, (region.start - 1) as i64, region.end as i64))?;
            while let Some(r) = self.reader.read(&mut self.record) {
                r?;
                let rec = &self.record;
                let flags = rec.flags();
                if flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_DUPLICATE) != 0 {
                    continue;
                }
                if flags & FLAG_SUPPLEMENTARY != 0 && !keep_supplementaries {
                    continue;
                }
                // already returned by an earlier interval
                if prev_end > 0 && rec.pos() as i32 + 1 <= prev_end {
                    continue;
                }
                consumer(rec, &mut self.scratch);
            }
            prev_end = region.end;
        }
        Ok(())
    }
}
