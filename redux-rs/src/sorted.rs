//! SortedBamWriter: holds records until no earlier one can arrive, then writes them in coordinate order.
//! A record that can no longer be placed in order is dropped and counted, as in Java.

use crate::out::RawBamOut;
use crate::rec::*;
use anyhow::Result;
use rust_htslib::bam;
use std::cmp::Ordering;
use std::collections::BTreeSet;

const READ_POS_CACHE_BUFFER: i32 = 30;
const MIN_WRITE_COUNT: usize = 10;
const POS_BUFFER_SIZE: i32 = 500;

struct Entry {
    start: i32,
    counter: u64,
    rec: bam::Record,
}

/// htsjdk SAMRecordCoordinateComparator within one reference, then the order of arrival
impl Ord for Entry {
    fn cmp(&self, o: &Self) -> Ordering {
        let c = self.start.cmp(&o.start);
        if c != Ordering::Equal {
            return c;
        }
        let (a, b) = (&self.rec, &o.rec);
        let (ar, br) = (has(a, FLAG_REVERSE), has(b, FLAG_REVERSE));
        if ar != br {
            return if ar { Ordering::Greater } else { Ordering::Less };
        }
        a.qname()
            .cmp(b.qname())
            .then(a.flags().cmp(&b.flags()))
            .then(a.mapq().cmp(&b.mapq()))
            .then(a.mtid().cmp(&b.mtid()))
            .then(a.mpos().cmp(&b.mpos()))
            .then(a.insert_size().cmp(&b.insert_size()))
            .then(self.counter.cmp(&o.counter))
    }
}
impl PartialOrd for Entry {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl PartialEq for Entry {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for Entry {}

pub struct SortedBamWriter {
    records: BTreeSet<Entry>,
    current_tid: i32,
    last_written_position: i32,
    upper_bound_position: i32,
    upper_writable_position: i32,
    min_cached_position: i32,
    read_count: u64,
    pub reads_written: u64,
    /// written records, kept for reuse so that each read does not allocate a new one
    pub pool: Vec<bam::Record>,
    position_buffer: i32,
    read_pos_cache_buffer: i32,
    min_write_count: usize,
}

const POOL_MAX: usize = 4096;

impl SortedBamWriter {
    pub fn new() -> Self {
        Self::with_config(POS_BUFFER_SIZE, READ_POS_CACHE_BUFFER, MIN_WRITE_COUNT)
    }

    /// SortedBamConfig(positionBuffer, readPosCacheBuffer, minWriteCount)
    pub fn with_config(position_buffer: i32, read_pos_cache_buffer: i32, min_write_count: usize) -> Self {
        SortedBamWriter {
            records: BTreeSet::new(),
            current_tid: -2,
            last_written_position: 0,
            upper_bound_position: 0,
            upper_writable_position: 0,
            min_cached_position: 0,
            read_count: 0,
            reads_written: 0,
            pool: Vec::new(),
            position_buffer,
            read_pos_cache_buffer,
            min_write_count,
        }
    }

    pub fn initialise_start_position(&mut self, tid: i32, start_position: i32) {
        if self.current_tid != tid {
            self.current_tid = tid;
            self.last_written_position = 0;
        }
        self.min_cached_position = start_position;
        self.upper_bound_position = start_position + self.position_buffer;
        self.upper_writable_position = start_position;
    }
    pub fn set_upper_bound_position(&mut self, position: i32) {
        self.upper_bound_position = self.upper_bound_position.max(position);
    }
    pub fn set_upper_writable_position(&mut self, position: i32) {
        self.upper_writable_position = position - self.read_pos_cache_buffer;
    }

    #[cfg(test)]
    pub fn cached(&self) -> usize {
        self.records.len()
    }

    pub fn can_write_record(&self, rec: &bam::Record) -> bool {
        if rec.tid() != self.current_tid {
            return false;
        }
        let s = start(rec);
        s >= self.last_written_position && s <= self.upper_bound_position
    }

    pub fn add_record(&mut self, rec: bam::Record, out: &mut RawBamOut) -> Result<()> {
        let s = start(&rec);
        self.records.insert(Entry { start: s, counter: self.read_count, rec });
        self.read_count += 1;
        self.min_cached_position = if self.min_cached_position == 0 { s } else { self.min_cached_position.min(s) };
        if self.upper_writable_position < self.min_cached_position {
            return Ok(());
        }
        let upper = self.upper_writable_position;
        let writable = self.records.iter().take_while(|e| e.start <= upper).take(self.min_write_count).count();
        if writable < self.min_write_count {
            return Ok(());
        }
        self.write_records(Some(upper), out)?;
        self.min_cached_position = match self.records.first() {
            Some(e) => e.start,
            None => self.upper_writable_position + 1,
        };
        Ok(())
    }

    /// Writes the records starting at or before `upper`, or all of them.
    fn write_records(&mut self, upper: Option<i32>, out: &mut RawBamOut) -> Result<()> {
        let Some(first) = self.records.first() else { return Ok(()) };
        if let Some(u) = upper {
            if first.start > u {
                return Ok(());
            }
        }
        if first.start < self.last_written_position {
            anyhow::bail!(
                "sorted BAM cache writing earlier(readStart={} vs last={}) from {} records",
                first.start, self.last_written_position, self.records.len()
            );
        }
        let mut last = self.last_written_position;
        while let Some(e) = self.records.first() {
            if upper.is_some_and(|u| e.start > u) {
                break;
            }
            let e = self.records.pop_first().unwrap();
            out.write(&e.rec)?;
            last = e.start;
            self.reads_written += 1;
            if self.pool.len() < POOL_MAX {
                self.pool.push(e.rec);
            }
        }
        self.last_written_position = last;
        Ok(())
    }

    pub fn flush(&mut self, out: &mut RawBamOut) -> Result<()> {
        self.write_records(None, out)
    }
}
