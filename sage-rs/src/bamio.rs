//! BAM slicing: the equivalent of SamSlicer / BamSlicer over an indexed BAM.

use crate::types::*;
use anyhow::{Context, Result};
use rust_htslib::bam::{self, record::Aux, record::Cigar, Read as BamRead};

pub struct BamSlicer {
    reader: bam::IndexedReader,
    record: bam::Record,
}

fn convert(rec: &bam::Record) -> Read {
    let cigar: Vec<CigarEl> = rec
        .cigar()
        .iter()
        .map(|c| match *c {
            Cigar::Match(l) => CigarEl::new(l as i32, Op::M),
            Cigar::Ins(l) => CigarEl::new(l as i32, Op::I),
            Cigar::Del(l) => CigarEl::new(l as i32, Op::D),
            Cigar::RefSkip(l) => CigarEl::new(l as i32, Op::N),
            Cigar::SoftClip(l) => CigarEl::new(l as i32, Op::S),
            Cigar::HardClip(l) => CigarEl::new(l as i32, Op::H),
            Cigar::Pad(l) => CigarEl::new(l as i32, Op::P),
            Cigar::Equal(l) => CigarEl::new(l as i32, Op::Eq),
            Cigar::Diff(l) => CigarEl::new(l as i32, Op::X),
        })
        .collect();
    let ref_len: i32 = cigar.iter().filter(|c| c.op.consumes_ref()).map(|c| c.len).sum();
    let start = rec.pos() as i32 + 1;
    let nm = match rec.aux(b"NM") {
        Ok(Aux::U8(v)) => Some(v as i32),
        Ok(Aux::I8(v)) => Some(v as i32),
        Ok(Aux::U16(v)) => Some(v as i32),
        Ok(Aux::I16(v)) => Some(v as i32),
        Ok(Aux::U32(v)) => Some(v as i32),
        Ok(Aux::I32(v)) => Some(v),
        _ => None,
    };
    let mate_cigar = match rec.aux(b"MC") {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    };
    // SamRecordUtils.extractConsensusType: UT tag, set by REDUX
    let (consensus, has_consensus_tag) = match rec.aux(b"UT") {
        Ok(Aux::String(s)) if s == "DUAL" => (ConsensusType::Dual, true),
        Ok(Aux::String(s)) if s == "SINGLE" => (ConsensusType::Single, true),
        Ok(Aux::String(_)) => (ConsensusType::None, true),
        _ => (ConsensusType::None, false),
    };
    Read {
        name: String::from_utf8_lossy(rec.qname()).into_owned(),
        flags: rec.flags(),
        tid: rec.tid(),
        start,
        end: start + ref_len - 1,
        mapq: rec.mapq() as i32,
        cigar,
        bases: rec.seq().as_bytes(),
        quals: rec.qual().to_vec(),
        mate_tid: rec.mtid(),
        mate_start: rec.mpos() as i32 + 1,
        insert_size: rec.insert_size() as i32,
        nm,
        mate_cigar,
        consensus,
        has_consensus_tag,
    }
}

impl BamSlicer {
    pub fn open(path: &str, threads: usize) -> Result<Self> {
        let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("opening BAM {path}"))?;
        if threads > 1 {
            reader.set_threads(threads)?;
        }
        Ok(BamSlicer { reader, record: bam::Record::new() })
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
                let read = convert(rec);
                // already returned by an earlier interval
                if prev_end > 0 && read.start <= prev_end {
                    continue;
                }
                consumer(read);
            }
            prev_end = region.end;
        }
        Ok(())
    }
}
