//! Region slices over one or more input BAMs, merged by position as the Java BamReader does.

use anyhow::Result;
use rust_htslib::bam::{self, Read};

struct FileReader {
    reader: bam::IndexedReader,
    current: Option<bam::Record>,
    active: bool,
}

impl FileReader {
    fn slice(&mut self, chr: &str, start: i32, end: i32) {
        self.current = None;
        self.active = false;
        let Some(tid) = self.reader.header().tid(chr.as_bytes()) else { return };
        // htsjdk queryOverlapping with 1-based inclusive start and end
        if self.reader.fetch((tid, (start - 1) as i64, end as i64)).is_err() {
            return;
        }
        self.active = true;
        self.move_next();
    }
    fn move_next(&mut self) {
        let mut rec = self.current.take().unwrap_or_else(bam::Record::new);
        match self.reader.read(&mut rec) {
            Some(Ok(())) => self.current = Some(rec),
            _ => {
                self.current = None;
                self.active = false;
            }
        }
    }
    /// 1-based alignment start
    fn position(&self) -> i32 {
        self.current.as_ref().map(|r| r.pos() as i32 + 1).unwrap_or(-1)
    }
}

/// A coordinate-sorted BAM of placed reads without alignments, read front to back and merged in as one more
/// reader after the input files.
struct ExtraStream {
    path: String,
    reader: bam::Reader,
    next: Option<bam::Record>,
    /// (reference index, 1-based position) of the last slice start, to notice a slice that goes back
    last_start: (i32, i32),
}

impl ExtraStream {
    fn open(path: &str) -> Result<Self> {
        let mut e = ExtraStream { path: path.to_string(), reader: bam::Reader::from_path(path)?, next: None, last_start: (-1, 0) };
        e.advance();
        Ok(e)
    }
    fn advance(&mut self) {
        let mut rec = self.next.take().unwrap_or_else(bam::Record::new);
        self.next = match self.reader.read(&mut rec) {
            Some(Ok(())) => Some(rec),
            _ => None,
        };
    }
    fn key(&self) -> Option<(i32, i32)> {
        self.next.as_ref().map(|r| (r.tid(), r.pos() as i32 + 1))
    }
    /// Moves to the first record at or after (tid, start).
    fn seek(&mut self, tid: i32, start: i32) -> Result<()> {
        if (tid, start) < self.last_start {
            *self = ExtraStream::open(&self.path)?;
        }
        self.last_start = (tid, start);
        while self.key().is_some_and(|k| k < (tid, start)) {
            self.advance();
        }
        Ok(())
    }
}

pub struct BamReader {
    readers: Vec<FileReader>,
    extra: Option<ExtraStream>,
}

impl BamReader {
    pub fn open(files: &[String], ref_genome: &str, decompress_threads: usize) -> Result<Self> {
        let mut readers = Vec::new();
        for f in files {
            let mut reader = bam::IndexedReader::from_path(f)?;
            reader.set_reference(ref_genome)?;
            if decompress_threads > 1 {
                reader.set_threads(decompress_threads)?;
            }
            readers.push(FileReader { reader, current: None, active: false });
        }
        Ok(BamReader { readers, extra: None })
    }

    pub fn header(&self) -> &bam::HeaderView {
        self.readers[0].reader.header()
    }

    /// Sets the extra sorted BAM, or none.
    pub fn set_extra(&mut self, path: Option<&str>) -> Result<()> {
        self.extra = match path {
            Some(p) => Some(ExtraStream::open(p)?),
            None => None,
        };
        Ok(())
    }

    /// Calls `consumer` with every record overlapping the region, in the order of the Java merge:
    /// by alignment start, a reader keeps its turn while its next start is not above the next reader's.
    pub fn slice_region<F: FnMut(&mut bam::Record)>(&mut self, chr: &str, start: i32, end: i32, mut consumer: F) {
        if self.extra.is_some() {
            return self.slice_region_with_extra(chr, start, end, consumer);
        }
        if self.readers.len() == 1 {
            let r = &mut self.readers[0];
            r.slice(chr, start, end);
            while r.current.is_some() {
                consumer(r.current.as_mut().unwrap());
                r.move_next();
            }
            return;
        }
        // indices of active readers ordered by position
        let mut active: Vec<usize> = Vec::new();
        for i in 0..self.readers.len() {
            self.readers[i].slice(chr, start, end);
            if self.readers[i].current.is_some() {
                Self::add_in_position(&self.readers, &mut active, i);
            }
        }
        while !active.is_empty() {
            let top = active[0];
            consumer(self.readers[top].current.as_mut().unwrap());
            self.readers[top].move_next();
            if self.readers[top].current.is_some() && active.len() > 1 && self.readers[top].position() <= self.readers[active[1]].position() {
                continue;
            }
            active.remove(0);
            if self.readers[top].current.is_some() {
                Self::add_in_position(&self.readers, &mut active, top);
            }
        }
    }

    /// The same merge with the extra BAM as the last reader. Its records are placed reads without an alignment,
    /// so each is in the slice when its position is.
    fn slice_region_with_extra<F: FnMut(&mut bam::Record)>(&mut self, chr: &str, start: i32, end: i32, mut consumer: F) {
        let n = self.readers.len();
        let tid = self.header().tid(chr.as_bytes()).map(|t| t as i32).unwrap_or(-2);
        let extra = self.extra.as_mut().unwrap();
        if extra.seek(tid, start).is_err() {
            return;
        }
        let readers = &mut self.readers;
        // the extra reader's position, or -1 past the slice
        let extra_pos = |extra: &ExtraStream| -> i32 {
            match extra.key() {
                Some((t, p)) if t == tid && p <= end => p,
                _ => -1,
            }
        };
        let position = |readers: &[FileReader], extra: &ExtraStream, i: usize| -> i32 { if i < n { readers[i].position() } else { extra_pos(extra) } };
        let mut active: Vec<usize> = Vec::new();
        let add = |active: &mut Vec<usize>, readers: &[FileReader], extra: &ExtraStream, i: usize| {
            let p = position(readers, extra, i);
            let mut index = 0;
            while index < active.len() && p > position(readers, extra, active[index]) {
                index += 1;
            }
            active.insert(index, i);
        };
        for i in 0..n {
            readers[i].slice(chr, start, end);
            if readers[i].current.is_some() {
                add(&mut active, readers, extra, i);
            }
        }
        if extra_pos(extra) >= 0 {
            add(&mut active, readers, extra, n);
        }
        while !active.is_empty() {
            let top = active[0];
            let finished;
            if top < n {
                consumer(readers[top].current.as_mut().unwrap());
                readers[top].move_next();
                finished = readers[top].current.is_none();
            } else {
                consumer(extra.next.as_mut().unwrap());
                extra.advance();
                finished = extra_pos(extra) < 0;
            }
            if !finished && active.len() > 1 && position(readers, extra, top) <= position(readers, extra, active[1]) {
                continue;
            }
            active.remove(0);
            if !finished {
                add(&mut active, readers, extra, top);
            }
        }
    }

    fn add_in_position(readers: &[FileReader], active: &mut Vec<usize>, i: usize) {
        let mut index = 0;
        while index < active.len() && readers[i].position() > readers[active[index]].position() {
            index += 1;
        }
        active.insert(index, i);
    }
}
