//! SortedBamWriter: holds records until no earlier one can arrive, then writes them in coordinate order.
//! A record that can no longer be placed in order is dropped and counted, as in Java.

use crate::out::RawBamOut;
use crate::rec::*;
use anyhow::Result;
use rust_htslib::bam;
use std::cmp::Ordering;
use std::collections::BTreeMap;

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
    /// the held records by start position; each list is in order of arrival and is sorted when written
    records: BTreeMap<i32, Vec<Entry>>,
    cached: usize,
    /// emptied lists, kept for reuse
    spare: Vec<Vec<Entry>>,
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
            records: BTreeMap::new(),
            cached: 0,
            spare: Vec::new(),
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
        self.cached
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
        let entry = Entry { start: s, counter: self.read_count, rec };
        match self.records.entry(s) {
            std::collections::btree_map::Entry::Occupied(mut o) => o.get_mut().push(entry),
            std::collections::btree_map::Entry::Vacant(v) => {
                let mut list = self.spare.pop().unwrap_or_default();
                list.push(entry);
                v.insert(list);
            }
        }
        self.cached += 1;
        self.read_count += 1;
        self.min_cached_position = if self.min_cached_position == 0 { s } else { self.min_cached_position.min(s) };
        if self.upper_writable_position < self.min_cached_position {
            return Ok(());
        }
        let upper = self.upper_writable_position;
        let mut writable = 0;
        for (_, list) in self.records.range(..=upper) {
            writable += list.len();
            if writable >= self.min_write_count {
                break;
            }
        }
        if writable < self.min_write_count {
            return Ok(());
        }
        self.write_records(Some(upper), out)?;
        self.min_cached_position = match self.records.first_key_value() {
            Some((&start, _)) => start,
            None => self.upper_writable_position + 1,
        };
        Ok(())
    }

    /// Writes the records starting at or before `upper`, or all of them.
    fn write_records(&mut self, upper: Option<i32>, out: &mut RawBamOut) -> Result<()> {
        let Some((&first_start, _)) = self.records.first_key_value() else { return Ok(()) };
        if let Some(u) = upper {
            if first_start > u {
                return Ok(());
            }
        }
        if first_start < self.last_written_position {
            anyhow::bail!(
                "sorted BAM cache writing earlier(readStart={} vs last={}) from {} records",
                first_start, self.last_written_position, self.cached
            );
        }
        let mut last = self.last_written_position;
        while let Some((&start, _)) = self.records.first_key_value() {
            if upper.is_some_and(|u| start > u) {
                break;
            }
            let (_, mut list) = self.records.pop_first().unwrap();
            if list.len() > 1 {
                list.sort();
            }
            for e in list.drain(..) {
                out.write(&e.rec)?;
                self.reads_written += 1;
                self.cached -= 1;
                if self.pool.len() < POOL_MAX {
                    self.pool.push(e.rec);
                }
            }
            last = start;
            if self.spare.len() < POOL_MAX {
                self.spare.push(list);
            }
        }
        self.last_written_position = last;
        Ok(())
    }

    pub fn flush(&mut self, out: &mut RawBamOut) -> Result<()> {
        self.write_records(None, out)
    }
}
