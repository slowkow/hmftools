//! BAM access: BamSlicer's region queries and filters, and the conversion of records to `Read`s.

use crate::read::{CigarElem, Read};
use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::{self, record::Aux, Read as BamRead};
use std::sync::OnceLock;

/// contig names of the BAM header, kept for the life of the program
static CONTIGS: OnceLock<Vec<&'static str>> = OnceLock::new();

pub fn contig_name(tid: i32) -> &'static str {
    if tid < 0 {
        return "*";
    }
    CONTIGS.get().and_then(|c| c.get(tid as usize).copied()).unwrap_or("*")
}

/// BamSlicerFilter
#[derive(Clone, Copy, Debug)]
pub struct SliceFilter {
    pub min_map_quality: i32,
    pub keep_duplicates: bool,
    pub keep_supplementaries: bool,
    pub keep_secondaries: bool,
}

impl SliceFilter {
    pub fn passes(&self, rec: &bam::Record) -> bool {
        if (rec.mapq() as i32) < self.min_map_quality {
            return false;
        }
        if rec.is_unmapped() {
            return false;
        }
        if rec.is_secondary() && !self.keep_secondaries {
            return false;
        }
        if rec.is_supplementary() && !self.keep_supplementaries {
            return false;
        }
        !rec.is_duplicate() || self.keep_duplicates
    }
}

pub struct BamReader {
    reader: bam::IndexedReader,
    record: bam::Record,
    tid_by_name: std::collections::HashMap<String, i32>,
}

impl BamReader {
    pub fn open(path: &str, ref_genome: Option<&str>) -> Result<BamReader> {
        let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("failed to open BAM {path}"))?;
        if let Some(r) = ref_genome {
            // CRAM input needs the reference; harmless for BAM
            let _ = reader.set_reference(r);
        }
        let header = reader.header().clone();
        let names: Vec<String> = (0..header.target_count()).map(|t| String::from_utf8_lossy(header.tid2name(t)).to_string()).collect();
        let _ = CONTIGS.get_or_init(|| names.iter().map(|n| &*Box::leak(n.clone().into_boxed_str())).collect());
        let tid_by_name = names.iter().enumerate().map(|(i, n)| (n.clone(), i as i32)).collect();
        Ok(BamReader { reader, record: bam::Record::new(), tid_by_name })
    }

    pub fn tid(&self, chromosome: &str) -> Option<i32> {
        self.tid_by_name.get(chromosome).copied()
    }

    /// BamSlicer.slice over one region: the records overlapping it (1-based inclusive) that pass the
    /// filter; the callback returns false to halt
    pub fn slice<F: FnMut(&bam::Record) -> bool>(&mut self, chromosome: &str, start: i32, end: i32, filter: &SliceFilter, mut f: F) -> Result<()> {
        let Some(tid) = self.tid(chromosome) else {
            hmf_common::warn!("cannot find sequence index for chromosome {} in bam header", chromosome);
            return Ok(());
        };
        // htsjdk's QueryInterval with an end of 0 or less (or past the contig) means to the end
        let fetch_end = if end <= 0 { i64::MAX } else { end as i64 };
        self.reader.fetch((tid, (start.max(1) - 1) as i64, fetch_end)).map_err(|e| anyhow!("BAM query failed: {e}"))?;
        loop {
            match self.reader.read(&mut self.record) {
                None => break,
                Some(Err(e)) => return Err(anyhow!("BAM read failed: {e}")),
                Some(Ok(())) => {}
            }
            if !filter.passes(&self.record) {
                continue;
            }
            if !f(&self.record) {
                break;
            }
        }
        Ok(())
    }
}

pub fn aux_string(rec: &bam::Record, tag: &[u8]) -> Option<String> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

pub fn record_cigar(rec: &bam::Record) -> Vec<CigarElem> {
    rec.raw_cigar().iter().map(|c| CigarElem { op: (c & 0xf) as u8, len: c >> 4 }).collect()
}

pub fn record_name(rec: &bam::Record) -> &str {
    std::str::from_utf8(rec.qname()).unwrap_or("")
}

/// SAMRecord.getAlignmentStart / getMateAlignmentStart: 1-based, 0 where unset
pub fn rec_start(rec: &bam::Record) -> i32 {
    (rec.pos() + 1) as i32
}
pub fn rec_mate_start(rec: &bam::Record) -> i32 {
    (rec.mpos() + 1) as i32
}
/// SAMRecord.getAlignmentEnd
pub fn rec_end(rec: &bam::Record) -> i32 {
    let len: i64 = rec.raw_cigar().iter().filter(|c| matches!(*c & 0xf, 0 | 2 | 3 | 7 | 8)).map(|c| (c >> 4) as i64).sum();
    (rec.pos() + len) as i32
}

pub fn to_read(rec: &bam::Record) -> Read {
    let sa = aux_string(rec, b"SA");
    let xa = aux_string(rec, b"XA");
    let cr = aux_string(rec, b"CR");
    Read::new(
        record_name(rec),
        contig_name(rec.tid()),
        rec_start(rec),
        rec.flags(),
        record_cigar(rec),
        rec.seq().as_bytes(),
        contig_name(rec.mtid()),
        rec_mate_start(rec),
        rec.insert_size() as i32,
        rec.mapq(),
        sa.as_deref(),
        xa.as_deref(),
        cr.as_deref(),
    )
}
