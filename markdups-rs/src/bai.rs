//! The .bai index as htsjdk's BAMIndexer builds it while BAMFileWriter writes (CREATE_INDEX): bins with
//! chunks coalesced when they start in the same or the next block address, the metadata pseudo-bin,
//! the linear index with empty windows filled in, and the count of records without a position.
//! A record's chunk is the virtual offset before and after it, an offset at the end of a full block
//! being the start of the next block (htsjdk's writer deflates a block as soon as it is full).

use std::collections::BTreeMap;

const MAX_BINS: u32 = 37450;
const UNINITIALIZED: i64 = -1;

/// One record as the indexer sees it.
#[derive(Clone, Copy, Debug)]
pub struct IndexEntry {
    pub ref_id: i32,
    /// 1-based alignment start; 0 for none
    pub start: i32,
    /// SAMRecord.getAlignmentEnd: 0 for unmapped reads
    pub end: i32,
    pub bin: u32,
    pub unmapped: bool,
    pub chunk_start: u64,
    pub chunk_end: u64,
}

struct RefBuilder {
    bins: BTreeMap<u32, Vec<(u64, u64)>>,
    index: Vec<i64>,
    largest: i32,
    first_offset: i64,
    last_offset: u64,
    aligned: u64,
    unaligned: u64,
}

impl RefBuilder {
    fn new() -> Self {
        RefBuilder { bins: BTreeMap::new(), index: Vec::new(), largest: -1, first_offset: -1, last_offset: 0, aligned: 0, unaligned: 0 }
    }
}

pub struct BaiBuilder {
    n_ref: usize,
    current: i32,
    builder: RefBuilder,
    no_coordinate: u64,
    out: Vec<u8>,
}

/// LinearIndex.convertToLinearIndexOffset
fn window(pos: i32) -> i32 {
    let p = if pos <= 0 { 0 } else { pos - 1 };
    p >> 14
}

impl BaiBuilder {
    pub fn new(n_ref: usize) -> Self {
        let mut out = Vec::new();
        out.extend_from_slice(b"BAI\x01");
        out.extend_from_slice(&(n_ref as i32).to_le_bytes());
        BaiBuilder { n_ref, current: 0, builder: RefBuilder::new(), no_coordinate: 0, out }
    }

    /// BAMIndexer.processAlignment
    pub fn add(&mut self, e: &IndexEntry) -> anyhow::Result<()> {
        if e.ref_id != -1 && e.ref_id != self.current {
            if e.ref_id < self.current {
                anyhow::bail!("Exception creating BAM index: records out of order (reference {} after {})", e.ref_id, self.current);
            }
            self.advance_to(e.ref_id as usize);
        }
        // BAMIndexMetaData.recordMetaData
        if e.start == 0 {
            self.no_coordinate += 1;
            return Ok(());
        }
        let b = &mut self.builder;
        if e.unmapped {
            b.unaligned += 1;
        } else {
            b.aligned += 1;
        }
        if (e.chunk_start as i64) <= b.first_offset || b.first_offset == -1 {
            b.first_offset = e.chunk_start as i64;
        }
        if b.last_offset <= e.chunk_end {
            b.last_offset = e.chunk_end;
        }
        if e.ref_id != self.current {
            anyhow::bail!("Unexpected reference {} when constructing index for {}", e.ref_id, self.current);
        }
        // BinningIndexBuilder.processFeature
        let chunks = b.bins.entry(e.bin).or_default();
        match chunks.last_mut() {
            Some(last) if (last.1 >> 16) == (e.chunk_start >> 16) || (last.1 >> 16) + 1 == (e.chunk_start >> 16) => last.1 = e.chunk_end,
            _ => chunks.push((e.chunk_start, e.chunk_end)),
        }
        let (start_window, end_window) = if e.end == 0 {
            let w = window(e.start - 1);
            (w, w)
        } else {
            (window(e.start), window(e.end))
        };
        if end_window > b.largest {
            b.largest = end_window;
        }
        if end_window >= 0 && b.index.len() <= end_window as usize {
            b.index.resize(end_window as usize + 1, UNINITIALIZED);
        }
        for w in start_window.max(0)..=end_window {
            let slot = &mut b.index[w as usize];
            if *slot == UNINITIALIZED || (e.chunk_start as i64) < *slot {
                *slot = e.chunk_start as i64;
            }
        }
        Ok(())
    }

    fn advance_to(&mut self, next: usize) {
        while (self.current as usize) < next {
            self.write_reference();
            self.current += 1;
            self.builder = RefBuilder::new();
        }
    }

    fn write_reference(&mut self) {
        let b = &mut self.builder;
        if b.bins.is_empty() {
            self.out.extend_from_slice(&0u64.to_le_bytes());
            return;
        }
        let o = &mut self.out;
        o.extend_from_slice(&((b.bins.len() + 1) as i32).to_le_bytes());
        for (bin, chunks) in &b.bins {
            if *bin == MAX_BINS {
                continue;
            }
            o.extend_from_slice(&bin.to_le_bytes());
            o.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
            for (s, e) in chunks {
                o.extend_from_slice(&s.to_le_bytes());
                o.extend_from_slice(&e.to_le_bytes());
            }
        }
        o.extend_from_slice(&MAX_BINS.to_le_bytes());
        o.extend_from_slice(&2i32.to_le_bytes());
        o.extend_from_slice(&b.first_offset.to_le_bytes());
        o.extend_from_slice(&b.last_offset.to_le_bytes());
        o.extend_from_slice(&b.aligned.to_le_bytes());
        o.extend_from_slice(&b.unaligned.to_le_bytes());
        // linear index up to the largest window seen, empty windows filled with the last value before them
        let n = (b.largest + 1).max(0) as usize;
        o.extend_from_slice(&(n as i32).to_le_bytes());
        let mut last = 0i64;
        for w in 0..n {
            let v = b.index.get(w).copied().unwrap_or(UNINITIALIZED);
            let v = if v == UNINITIALIZED { last } else { v };
            last = v;
            o.extend_from_slice(&v.to_le_bytes());
        }
    }

    /// BAMIndexer.finish: the remaining references and the count of records without a position.
    pub fn finish(mut self) -> Vec<u8> {
        if self.n_ref > 0 {
            self.advance_to(self.n_ref);
        }
        self.out.extend_from_slice(&self.no_coordinate.to_le_bytes());
        self.out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(ref_id: i32, start: i32, end: i32, cs: u64, ce: u64) -> IndexEntry {
        let bin = crate::bam::region_to_bin(start - 1, if end <= 0 { start } else { end }) as u32;
        IndexEntry { ref_id, start, end, bin, unmapped: end == 0, chunk_start: cs, chunk_end: ce }
    }

    /// The index GATK 4.6.1.0 wrote for its own output of the tiny synthetic case (records in block 0 from
    /// offset 0x65f; the last chunk ending at 0x0be5 in the same block).
    #[test]
    fn as_htsjdk_writer() {
        let mut b = BaiBuilder::new(2);
        let offs = [0x65fu64, 0x6f0, 0x781, 0x812, 0x89b, 0x924, 0xa07, 0xaf6, 0xbe5];
        let starts = [100, 100, 100, 200, 200, 300, 300, 300];
        for i in 0..8 {
            b.add(&entry(0, starts[i], starts[i] + 49, offs[i], offs[i + 1])).unwrap();
        }
        let out = b.finish();
        let want: Vec<u8> = vec![
            0x42, 0x41, 0x49, 0x01, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x49, 0x12, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x5f, 0x06, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0xe5, 0x0b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4a, 0x92, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x5f, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0xe5, 0x0b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00,
            0x00, 0x00, 0x5f, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        assert_eq!(out, want);
    }

    #[test]
    fn unplaced_records_count_and_empty_references() {
        let mut b = BaiBuilder::new(3);
        b.add(&entry(1, 20000, 20100, 1 << 16, 2 << 16)).unwrap();
        b.add(&IndexEntry { ref_id: -1, start: 0, end: 0, bin: 4680, unmapped: true, chunk_start: 2 << 16, chunk_end: 3 << 16 }).unwrap();
        let out = b.finish();
        // reference 0 and 2 empty (8 zero bytes each); no-coordinate count 1 at the end
        assert_eq!(&out[8..16], &[0u8; 8]);
        assert_eq!(&out[out.len() - 8..], &1u64.to_le_bytes());
        assert_eq!(&out[out.len() - 16..out.len() - 8], &[0u8; 8]);
    }
}
