//! BamWriter.captureReadInfo: the reads that feed the jitter analysis and the base quality recalibration.

use crate::bqr::BqrRegionReader;
use crate::jitter::MsJitterProcessor;
use rust_htslib::bam;

pub const BQR_MIN_MAP_QUAL: u8 = 50;

pub struct ReadCapture<'a> {
    pub jitter: Option<MsJitterProcessor<'a>>,
    pub bqr: BqrRegionReader<'a>,
}

impl ReadCapture<'_> {
    #[inline]
    pub fn capture(&mut self, rec: &bam::Record, chromosome: &str) {
        // unmapped, duplicate, supplementary, secondary
        if rec.flags() & (0x4 | 0x400 | 0x800 | 0x100) != 0 {
            return;
        }
        if rec.mapq() < BQR_MIN_MAP_QUAL {
            return;
        }
        if let Some(j) = &mut self.jitter {
            j.process_read(rec, chromosome);
        }
        if self.bqr.is_active() {
            self.bqr.process_record(rec);
        }
    }
}
