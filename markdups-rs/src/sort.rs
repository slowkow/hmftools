//! The coordinate sort of samtools 1.21 (`samtools sort`), for the one-step mode: STAR's unsorted BAM in,
//! the records in the order samtools sort would write them out, without writing the sorted BAM.
//!
//! samtools sorts by (reference index as an unsigned number, so unplaced reads come last; position + 1;
//! reverse strand), and keeps the input order among equal keys (its in-memory sort is a merge sort and
//! its merge of the temporary files breaks ties by file and record number). Its header is the input's
//! with @HD SO:coordinate (GO and SS removed) and a @PG line for itself.
//!
//! Records are collected up to a memory budget, sorted by that key and the input number, and written
//! to a temporary run (BGZF, deflate level 1) in TMP_DIR; the last run stays in memory. Each pass of
//! MarkDuplicates merges the runs again.

use anyhow::{Context, Result, bail};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::bam::{IdRecord, Rec, SamHeader, read_header};
use crate::bgzf::{BgzfReader, BgzfWriter};
use crate::input::{InputSet, RecordSource};

/// samtools' sort key: (tid as u32, pos + 1, reverse)
#[inline]
pub fn sort_key(rec: &[u8]) -> (u32, u64, bool) {
    let r = Rec(rec);
    (r.ref_id() as u32, (r.pos() as i64 + 1) as u64, r.is_reverse())
}

/// How samtools sort changes the header: SO:coordinate on @HD (added as `@HD VN:1.6 SO:coordinate` if
/// there is none), GO and SS removed, and its @PG line(s) when its command line is known.
pub fn samtools_sorted_header(mut h: SamHeader, had_hd: bool, pg_command_line: Option<&str>, version: &str) -> SamHeader {
    if !had_hd {
        // sam_hdr_add_line(HD, VN, SAM_FORMAT_VERSION, SO, coordinate)
        h.hd = vec![(b"VN".to_vec(), b"1.6".to_vec())];
    }
    h.set_hd(b"SO", Some(b"coordinate"));
    h.set_hd(b"GO", None);
    h.set_hd(b"SS", None);
    if let Some(cl) = pg_command_line {
        // sam_hdr_add_pg: one line per end of a PP chain, each with a free id (samtools, samtools.1, ...)
        let referenced: Vec<Vec<u8>> = h.pg.iter().filter_map(|p| p.get(b"PP").map(|v| v.to_vec())).collect();
        let ends: Vec<Vec<u8>> = h.pg.iter().filter(|p| !referenced.contains(&p.id)).map(|p| p.id.clone()).collect();
        let mut counter = 1;
        let mut new_id = |h: &SamHeader| -> Vec<u8> {
            if !h.pg.iter().any(|p| p.id == b"samtools") {
                return b"samtools".to_vec();
            }
            loop {
                let id = format!("samtools.{counter}").into_bytes();
                counter += 1;
                if !h.pg.iter().any(|p| p.id == id) {
                    return id;
                }
            }
        };
        let make = |id: Vec<u8>, pp: Option<&[u8]>| {
            let mut r = IdRecord { id, attrs: vec![(b"PN".to_vec(), b"samtools".to_vec())] };
            if let Some(pp) = pp {
                r.set(b"PP", Some(pp));
            }
            r.set(b"VN", Some(version.as_bytes()));
            r.set(b"CL", Some(cl.as_bytes()));
            r
        };
        if ends.is_empty() {
            let id = new_id(&h);
            h.pg.push(make(id, None));
        } else {
            for e in ends {
                let id = new_id(&h);
                h.pg.push(make(id, Some(&e)));
            }
        }
    }
    h
}

struct MemRun {
    data: Vec<u8>,
    /// (offset of the record's length prefix) in sorted order
    order: Vec<usize>,
}

/// The sorted input: temporary runs on disk and the last one in memory.
pub struct SortedInput {
    header: SamHeader,
    runs: Vec<String>,
    mem: MemRun,
    pub records: u64,
}

impl Drop for SortedInput {
    fn drop(&mut self) {
        for r in &self.runs {
            let _ = std::fs::remove_file(r);
        }
    }
}

/// Reads `path` (any order) and sorts it as samtools sort would.
pub fn sort_input(
    path: &str,
    tmp_dir: &str,
    memory_bytes: usize,
    reader_threads: usize,
    compress_threads: usize,
    pg_command_line: Option<&str>,
    samtools_version: &str,
) -> Result<SortedInput> {
    let mut r = BgzfReader::open(path, reader_threads)?;
    let h = read_header(&mut r)?;
    let had_hd = h.sam.has_hd;
    let header = samtools_sorted_header(h.sam, had_hd, pg_command_line, samtools_version);
    let mut runs = Vec::new();
    // reserved, not touched: the pages are only used as records come
    let mut data: Vec<u8> = Vec::with_capacity(memory_bytes + (1 << 20));
    let mut keys: Vec<(u32, u64, bool, u64, usize)> = Vec::new();
    let mut rec = Vec::new();
    let mut n: u64 = 0;
    let pid = std::process::id();
    while r.read_record(&mut rec)? {
        Rec::check(&rec)?;
        let (tid, pos, rev) = sort_key(&rec);
        keys.push((tid, pos, rev, n, data.len()));
        data.extend_from_slice(&(rec.len() as u32).to_le_bytes());
        data.extend_from_slice(&rec);
        n += 1;
        if data.len() + keys.len() * 32 >= memory_bytes {
            keys.sort_unstable();
            let run = format!("{tmp_dir}/markdups-rs.{pid}.sort.{:04}.tmp", runs.len());
            let mut w = BgzfWriter::create(&run, 1, compress_threads).with_context(|| format!("temporary file in {tmp_dir}"))?;
            for &(.., off) in &keys {
                let len = u32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as usize;
                w.write(&data[off..off + 4 + len])?;
            }
            w.finish()?;
            runs.push(run);
            data.clear();
            keys.clear();
        }
    }
    keys.sort_unstable();
    let order = keys.iter().map(|k| k.4).collect();
    crate::log_info!("Sorted {n} records ({} temporary runs)", runs.len());
    Ok(SortedInput { header, runs, mem: MemRun { data, order }, records: n })
}

struct RunReader {
    r: BgzfReader,
    rec: Vec<u8>,
}

struct Merger<'a> {
    files: Vec<RunReader>,
    mem: &'a MemRun,
    mem_at: usize,
    /// (key, run number): the smallest first; the memory run is the last
    heap: BinaryHeap<Reverse<((u32, u64, bool), usize)>>,
}

impl<'a> Merger<'a> {
    fn mem_rec(&self, i: usize) -> &'a [u8] {
        let off = self.mem.order[i];
        let len = u32::from_le_bytes(self.mem.data[off..off + 4].try_into().unwrap()) as usize;
        &self.mem.data[off + 4..off + 4 + len]
    }
}

impl RecordSource for Merger<'_> {
    fn next(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        let Some(Reverse((_, run))) = self.heap.pop() else { return Ok(false) };
        if run == self.files.len() {
            buf.clear();
            buf.extend_from_slice(self.mem_rec(self.mem_at));
            self.mem_at += 1;
            if self.mem_at < self.mem.order.len() {
                let k = sort_key(self.mem_rec(self.mem_at));
                self.heap.push(Reverse((k, run)));
            }
        } else {
            let f = &mut self.files[run];
            std::mem::swap(buf, &mut f.rec);
            if f.r.read_record(&mut f.rec)? {
                let k = sort_key(&f.rec);
                self.heap.push(Reverse((k, run)));
            }
        }
        Ok(true)
    }
}

impl InputSet for SortedInput {
    fn open(&self, threads: usize) -> Result<(SamHeader, Box<dyn RecordSource + '_>)> {
        let per = if self.runs.is_empty() { 0 } else { threads / self.runs.len() };
        let mut files = Vec::new();
        let mut heap = BinaryHeap::new();
        for (i, p) in self.runs.iter().enumerate() {
            let mut r = BgzfReader::open(p, per)?;
            let mut rec = Vec::new();
            if !r.read_record(&mut rec)? {
                bail!("{p}: empty temporary run");
            }
            heap.push(Reverse((sort_key(&rec), i)));
            files.push(RunReader { r, rec });
        }
        let mut m = Merger { files, mem: &self.mem, mem_at: 0, heap };
        if !self.mem.order.is_empty() {
            let k = sort_key(m.mem_rec(0));
            m.heap.push(Reverse((k, m.files.len())));
        }
        Ok((self.header.clone(), Box::new(m)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::SamHeader;

    #[test]
    fn header_as_samtools_sort() {
        let h = SamHeader::parse(b"@HD\tVN:1.4\tGO:query\n@SQ\tSN:c\tLN:10\n@PG\tID:STAR\tPN:STAR\n").unwrap();
        let s = samtools_sorted_header(h, true, Some("samtools sort -o x.bam in.bam"), "1.21");
        assert_eq!(
            String::from_utf8(s.encode()).unwrap(),
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:10\n@PG\tID:STAR\tPN:STAR\n@PG\tID:samtools\tPN:samtools\tPP:STAR\tVN:1.21\tCL:samtools sort -o x.bam in.bam\n"
        );
        let h = SamHeader::parse(b"@SQ\tSN:c\tLN:10\n@PG\tID:samtools\tPN:samtools\n").unwrap();
        let s = samtools_sorted_header(h, false, Some("samtools sort"), "1.21");
        assert_eq!(s.pg[1].id, b"samtools.1");
        assert_eq!(s.pg[1].get(b"PP"), Some(&b"samtools"[..]));
    }

    #[test]
    fn unplaced_reads_sort_last() {
        let mut a = vec![0u8; 36];
        a[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        a[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        let mut b = vec![0u8; 36];
        b[0..4].copy_from_slice(&5i32.to_le_bytes());
        b[4..8].copy_from_slice(&100i32.to_le_bytes());
        assert!(sort_key(&b) < sort_key(&a));
        assert_eq!(sort_key(&a), (u32::MAX, 0, false));
    }
}
