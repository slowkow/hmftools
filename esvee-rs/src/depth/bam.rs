//! Reads for a region from an indexed BAM or CRAM: Java BamSlicer(0, false, true, false) over
//! htsjdk's queryOverlapping.

use super::task::{ReadInfo, ReadSource, SamRead, FLAG_DUPLICATE, FLAG_SECONDARY, FLAG_UNMAPPED};
use anyhow::{Context, Result};
use rust_htslib::bam::{self, record::Aux, Read as BamRead};

pub struct BamSource {
    reader: bam::IndexedReader,
    record: bam::Record,
}

struct Record<'a> {
    record: &'a bam::Record,
    info: ReadInfo,
}

impl SamRead for Record<'_> {
    fn name(&self) -> &[u8] {
        self.record.qname()
    }
    fn info(&self) -> ReadInfo {
        self.info
    }
    fn supplementary_data(&self) -> Option<&[u8]> {
        match self.record.aux(b"SA") {
            Ok(Aux::String(s)) => Some(s.as_bytes()),
            _ => None,
        }
    }
}

const CIGAR_SOFT_CLIP: u32 = 4;

impl BamSource {
    pub fn open(path: &str, ref_genome: &str) -> Result<BamSource> {
        let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("opening BAM {path}"))?;
        if path.ends_with(".cram") {
            reader.set_reference(ref_genome).with_context(|| format!("setting CRAM reference {ref_genome}"))?;
        }
        Ok(BamSource { reader, record: bam::Record::new() })
    }
}

impl ReadSource for BamSource {
    fn slice(&mut self, chromosome: &str, start: i32, end: i32, consumer: &mut dyn FnMut(&dyn SamRead) -> Result<bool>) -> Result<()> {
        let Some(tid) = self.reader.header().tid(chromosome.as_bytes()) else {
            warn!("cannot find sequence index for chromosome {} in bam header", chromosome);
            return Ok(());
        };
        self.reader.fetch((tid, (start as i64 - 1).max(0), end as i64))?;

        while let Some(result) = self.reader.read(&mut self.record) {
            result?;
            let record = &self.record;
            let flags = record.flags();
            // BamSlicerFilter: unmapped reads, secondaries and duplicates are dropped
            if flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_DUPLICATE) != 0 {
                continue;
            }
            let cigar = record.raw_cigar();
            let mut ref_length = 0i32;
            for c in cigar {
                // M, D, N, = and X consume the reference
                if matches!(c & 0xf, 0 | 2 | 3 | 7 | 8) {
                    ref_length += (c >> 4) as i32;
                }
            }
            let read_start = record.pos() as i32 + 1;
            let read_end = read_start + ref_length - 1;
            // htsjdk takes a read as overlapping when its alignment end, which is one before the start for
            // a read covering no reference bases, is not before the region; htslib counts such a read as
            // one base long
            if read_end < start || read_start > end {
                continue;
            }
            let info = ReadInfo {
                start: read_start,
                end: read_end,
                flags,
                ref_index: record.tid(),
                mate_ref_index: record.mtid(),
                mate_start: record.mpos() as i32 + 1,
                insert_size: record.insert_size() as i32,
                read_length: record.seq_len() as i32,
                left_soft_clipped: cigar.first().is_some_and(|c| c & 0xf == CIGAR_SOFT_CLIP),
                right_soft_clipped: cigar.last().is_some_and(|c| c & 0xf == CIGAR_SOFT_CLIP),
            };
            if !consumer(&Record { record, info })? {
                break;
            }
        }
        Ok(())
    }
}
