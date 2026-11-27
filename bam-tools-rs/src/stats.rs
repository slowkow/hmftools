//! FlagStats, FragmentLengths, ReadCounts and PartitionStats.

use hmf_common::metrics::{FlagQcStats, FlagStatType as T, ValueFrequency, FLAG_STAT_TYPES};

use crate::compat;
use crate::record::{Rec, OP_N};

/// FlagStats: samtools-flagstat-like counts, QC-passed and failed
#[derive(Clone, Default)]
pub struct FlagStats {
    pub stats: [FlagQcStats; FLAG_STAT_TYPES],
}

const MAP_QUAL_THRESHOLD: u8 = 5;

impl FlagStats {
    #[inline]
    pub fn increment(&mut self, t: T, passes: bool) {
        let s = &mut self.stats[t as usize];
        if passes {
            s.passed += 1;
        } else {
            s.failed += 1;
        }
    }

    #[inline]
    fn decrement(&mut self, t: T, passes: bool) {
        let s = &mut self.stats[t as usize];
        if passes {
            s.passed -= 1;
        } else {
            s.failed -= 1;
        }
    }

    #[cfg(test)]
    pub fn pass_count(&self, t: T) -> i64 {
        self.stats[t as usize].passed
    }

    pub fn merge(&mut self, other: &FlagStats) {
        for i in 0..FLAG_STAT_TYPES {
            self.stats[i].passed += other.stats[i].passed;
            self.stats[i].failed += other.stats[i].failed;
        }
    }

    /// FlagStats.processRead
    pub fn process_read(&mut self, rec: &Rec, is_consensus: bool, expect_duplicates: bool) {
        let passes = !rec.qc_fail();
        let is_supplementary = rec.supplementary();
        // BT4: Java takes the supplementary flag for isSecondary. Corrected: the secondary flag.
        let is_secondary = if compat::java("BT4") { rec.supplementary() } else { rec.secondary() };
        if is_consensus && expect_duplicates {
            // consensus reads take back the extra duplicate; they add to no other count
            self.decrement(T::Duplicate, passes);
            if !is_supplementary && !is_secondary {
                self.decrement(T::PrimaryDuplicate, passes);
            }
            return;
        }
        self.increment(T::Total, passes);
        let is_duplicate = rec.duplicate();
        let is_mapped = !rec.unmapped();
        if is_duplicate {
            self.increment(T::Duplicate, passes);
        }
        if is_mapped {
            self.increment(T::Mapped, passes);
        }
        if rec.secondary() {
            self.increment(T::Secondary, passes);
            return;
        }
        if is_supplementary {
            self.increment(T::Supplementary, passes);
            return;
        }
        self.increment(T::Primary, passes);
        if is_duplicate {
            self.increment(T::PrimaryDuplicate, passes);
        }
        if is_mapped {
            self.increment(T::PrimaryMapped, passes);
        }
        if !rec.paired() {
            return;
        }
        self.increment(T::Paired, passes);
        if rec.first_of_pair() {
            self.increment(T::Read1, passes);
        }
        if rec.second_of_pair() {
            self.increment(T::Read2, passes);
        }
        if !is_mapped {
            return;
        }
        if rec.flag(crate::record::FLAG_PROPER_PAIR) {
            self.increment(T::ProperlyPaired, passes);
        }
        if rec.mate_unmapped() {
            self.increment(T::Singleton, passes);
            return;
        }
        self.increment(T::PairMapped, passes);
        if rec.tid != rec.mtid {
            self.increment(T::InterChrPairMapped, passes);
            if rec.mapq >= MAP_QUAL_THRESHOLD {
                self.increment(T::InterChrPairMapQualGe5, passes);
            }
        }
    }
}

pub const FRAG_LENGTH_DIST_MAX_LENGTH: usize = 1500;

/// FragmentLengths: fragments by length bucket
#[derive(Clone)]
pub struct FragmentLengths {
    pub counts: Vec<i64>,
}

impl Default for FragmentLengths {
    fn default() -> Self {
        FragmentLengths { counts: vec![0; FRAG_LENGTH_DIST_MAX_LENGTH + 1] }
    }
}

/// FragmentLengths.getLengthBucket: units up to 1000, tens up to 3000, then hundreds
#[inline]
pub fn length_bucket(length: i32) -> i32 {
    if length < 1000 {
        return length;
    }
    // Math.round: floor(x + 0.5)
    if length < 3000 {
        return 10 * (length as f64 / 10.0 + 0.5).floor() as i32;
    }
    100 * (length as f64 / 100.0 + 0.5).floor() as i32
}

impl FragmentLengths {
    /// FragmentLengths.processRead: primary non-duplicates, first of pair with the mate on the same
    /// chromosome and the other strand, no N in the CIGAR
    #[inline]
    pub fn process_read(&mut self, rec: &Rec) {
        if rec.duplicate() || rec.secondary() || rec.supplementary() {
            return;
        }
        if rec.paired() {
            if rec.second_of_pair() {
                return;
            }
            if rec.tid != rec.mtid {
                return;
            }
            if rec.mate_reverse() == rec.reverse() {
                return;
            }
        }
        if rec.has_op(OP_N) {
            return;
        }
        // SamRecordUtils.inferredInsertSizeAbs
        let size = if rec.paired() { rec.tlen } else { rec.end.wrapping_sub(rec.start).wrapping_add(1) };
        let size = size.wrapping_abs();
        let length = length_bucket(size);
        if length <= 0 || length as usize > FRAG_LENGTH_DIST_MAX_LENGTH {
            return;
        }
        self.counts[length as usize] += 1;
    }

    pub fn merge(&mut self, other: &FragmentLengths) {
        for (a, b) in self.counts.iter_mut().zip(&other.counts) {
            *a += *b;
        }
    }

    /// the lengths seen, ascending
    pub fn rows(&self) -> Vec<ValueFrequency> {
        self.counts.iter().enumerate().filter(|(_, c)| **c != 0).map(|(l, c)| ValueFrequency::new(l as i32, *c)).collect()
    }
}

/// ReadCounts
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct ReadCounts {
    pub total: i64,
    pub duplicates: i64,
    pub dual_strand: i64,
    pub off_target: i64,
}

impl ReadCounts {
    pub fn merge(&mut self, o: &ReadCounts) {
        self.total += o.total;
        self.duplicates += o.duplicates;
        self.dual_strand += o.dual_strand;
        self.off_target += o.off_target;
    }
}

/// PartitionStats
#[derive(Clone, Copy, Default, Debug)]
pub struct PartitionStats {
    pub total_reads: i64,
    pub duplicate_reads: i64,
    pub chimeric_reads: i64,
    pub inter_partition: i64,
    pub unmapped_reads: i64,
}

impl PartitionStats {
    /// PartitionStats.processRead; `partition` is (tid, start, end)
    #[inline]
    pub fn process_read(&mut self, rec: &Rec, partition: (i32, i32, i32), is_consensus: bool, expect_duplicates: bool) {
        if is_consensus {
            if expect_duplicates {
                self.duplicate_reads -= 1;
            } else {
                self.total_reads += 1;
            }
        } else {
            self.total_reads += 1;
            if rec.duplicate() {
                self.duplicate_reads += 1;
            }
        }
        if rec.unmapped() {
            self.unmapped_reads += 1;
        } else if rec.is_discordant() {
            self.chimeric_reads += 1;
            let (tid, s, e) = partition;
            if !(rec.mtid == tid && rec.mate_start >= s && rec.mate_start <= e) {
                self.inter_partition += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn buckets() {
        assert_eq!(length_bucket(999), 999);
        assert_eq!(length_bucket(1004), 1000);
        assert_eq!(length_bucket(1005), 1010);
        assert_eq!(length_bucket(1495), 1500);
        assert_eq!(length_bucket(1505), 1510);
        assert_eq!(length_bucket(3049), 3000);
        assert_eq!(length_bucket(3050), 3100);
    }
}
