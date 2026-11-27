//! VCFWriter.writeBaseDepths: the SNP-check sites of the reference as a bgzipped VCF, with the tabix index
//! htsjdk builds on the fly (TabixIndexCreator, BinningIndexBuilder, TabixIndex.write).
//!
//! The BGZF blocks hold 65498 uncompressed bytes each and are deflated by zlib at level 5, as htsjdk's
//! BlockCompressedOutputStream does, so the compressed blocks, their addresses and so the index are those
//! Java writes (as long as the zlib in use deflates as Java's; checked on the small pair).

use anyhow::{Context, Result};
use std::io::Write;

use crate::evidence::Evidence;

/// htsjdk BlockCompressedStreamConstants.DEFAULT_UNCOMPRESSED_BLOCK_SIZE (seen in its files: 65498 bytes)
const BLOCK_SIZE: usize = 65498;

/// The empty block that ends a BGZF file
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
];

/// A BGZF writer that reports htsjdk's virtual file pointers.
pub struct BgzfWriter<W: Write> {
    out: W,
    buffer: Vec<u8>,
    /// compressed bytes written so far: the address of the block being filled
    block_address: u64,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(out: W) -> Self {
        BgzfWriter { out, buffer: Vec::with_capacity(BLOCK_SIZE), block_address: 0 }
    }
    /// BlockCompressedOutputStream.getFilePointer
    pub fn file_pointer(&self) -> u64 {
        (self.block_address << 16) | self.buffer.len() as u64
    }
    pub fn write_all(&mut self, mut data: &[u8]) -> Result<()> {
        while !data.is_empty() {
            let n = (BLOCK_SIZE - self.buffer.len()).min(data.len());
            self.buffer.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.buffer.len() == BLOCK_SIZE {
                self.deflate_block()?;
            }
        }
        Ok(())
    }
    fn deflate_block(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let mut enc = flate2::write::DeflateEncoder::new(Vec::with_capacity(self.buffer.len()), flate2::Compression::new(5));
        enc.write_all(&self.buffer)?;
        let mut data = enc.finish()?;
        if data.len() + 26 > 65536 {
            // incompressible: store
            let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::none());
            enc.write_all(&self.buffer)?;
            data = enc.finish()?;
        }
        let mut crc = flate2::Crc::new();
        crc.update(&self.buffer);
        let bsize = (data.len() + 25) as u16;
        let mut block = Vec::with_capacity(data.len() + 26);
        block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00]);
        block.extend_from_slice(&bsize.to_le_bytes());
        block.extend_from_slice(&data);
        block.extend_from_slice(&crc.sum().to_le_bytes());
        block.extend_from_slice(&(self.buffer.len() as u32).to_le_bytes());
        self.out.write_all(&block)?;
        self.block_address += block.len() as u64;
        self.buffer.clear();
        Ok(())
    }
    pub fn finish(mut self) -> Result<W> {
        self.deflate_block()?;
        self.out.write_all(&BGZF_EOF)?;
        self.out.flush()?;
        Ok(self.out)
    }
}

/// GenomicIndexUtil.regionToBin (reg2bin), zero-based half-open
fn region_to_bin(beg: i32, end: i32) -> usize {
    let end = end - 1;
    let b = if beg >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (beg >> 26)
    } else {
        0
    };
    b as usize
}

/// LinearIndex.convertToLinearIndexOffset
fn linear_index_offset(pos: i32) -> usize {
    ((if pos <= 0 { 0 } else { pos - 1 }) >> 14) as usize
}

/// BinningIndexBuilder for one reference sequence
struct SequenceIndex {
    /// bin number to chunks (start, end); kept in bin order when written
    bins: std::collections::BTreeMap<usize, Vec<(u64, u64)>>,
    linear: Vec<u64>,
    largest_window: i64,
}

impl SequenceIndex {
    fn new() -> Self {
        SequenceIndex { bins: Default::default(), linear: Vec::new(), largest_window: -1 }
    }
    /// processFeature: a record at [start, end] (1-based) from file pointer chunk_start to chunk_end
    fn process(&mut self, start: i32, end: i32, chunk_start: u64, chunk_end: u64) {
        let bin = region_to_bin(start - 1, if end <= 0 { start } else { end });
        let chunks = self.bins.entry(bin).or_default();
        match chunks.last_mut() {
            // coalesce chunks in the same or adjacent blocks (BlockCompressedFilePointerUtil.areInSameOrAdjacentBlocks)
            Some(last) if (last.1 >> 16) == (chunk_start >> 16) || (last.1 >> 16) + 1 == (chunk_start >> 16) => last.1 = chunk_end,
            _ => chunks.push((chunk_start, chunk_end)),
        }
        let start_window = linear_index_offset(start);
        let end_window = if end == 0 { start_window } else { linear_index_offset(end) };
        if end_window as i64 > self.largest_window {
            self.largest_window = end_window as i64;
        }
        if self.linear.len() <= end_window {
            self.linear.resize(end_window + 1, 0);
        }
        for w in start_window..=end_window {
            if self.linear[w] == 0 || chunk_start < self.linear[w] {
                self.linear[w] = chunk_start;
            }
        }
    }
    /// generateIndexContent: the linear index up to the largest window, zeros after the first entry
    /// filled with the one before
    fn write(&mut self, out: &mut Vec<u8>) {
        out.extend_from_slice(&(self.bins.len() as i32).to_le_bytes());
        for (bin, chunks) in &self.bins {
            out.extend_from_slice(&(*bin as u32).to_le_bytes());
            out.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
            for (s, e) in chunks {
                out.extend_from_slice(&s.to_le_bytes());
                out.extend_from_slice(&e.to_le_bytes());
            }
        }
        let n = (self.largest_window + 1) as usize;
        self.linear.resize(n, 0);
        let mut last = 0u64;
        for v in self.linear.iter_mut() {
            if *v == 0 {
                *v = last;
            } else {
                last = *v;
            }
        }
        out.extend_from_slice(&(n as i32).to_le_bytes());
        for v in &self.linear {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
}

/// The header VCFWriter.header gives, as htsjdk writes it (meta lines sorted)
pub fn header_text(samples: &[&str]) -> String {
    let mut s = String::from("##fileformat=VCFv4.2\n");
    s.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    s.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n");
    s.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n");
    s.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    s.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for sample in samples {
        s.push('\t');
        s.push_str(sample);
    }
    s.push('\n');
    s
}

/// VCFWriter.writeBaseDepths: the sites sorted by position, each with GT 0/1, AD and DP; the index beside
pub fn write_base_depths(path: &str, sites: &[&Evidence], sample: &str) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let mut w = BgzfWriter::new(std::io::BufWriter::new(file));
    w.write_all(header_text(&[sample]).as_bytes())?;
    let mut names: Vec<&str> = Vec::new();
    let mut indexes: Vec<SequenceIndex> = Vec::new();
    // the record waiting for the next one's start as its chunk end
    let mut pending: Option<(i32, i32, u64)> = None;
    for e in sites {
        let pos = w.file_pointer();
        if names.last() != Some(&e.chr) {
            if names.contains(&e.chr) {
                anyhow::bail!("Sequence {} added out sequence of order", e.chr);
            }
            if let Some((s, en, start)) = pending.take() {
                indexes.last_mut().unwrap().process(s, en, start, pos);
            }
            names.push(e.chr);
            indexes.push(SequenceIndex::new());
        } else if let Some((s, en, start)) = pending.take() {
            indexes.last_mut().unwrap().process(s, en, start, pos);
        }
        let line = format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT:AD:DP\t0/1:{},{}:{}\n",
            e.chr,
            e.pos,
            e.ref_.as_str(),
            e.alt.as_str(),
            e.c.ref_support,
            e.c.alt_support,
            e.c.read_depth
        );
        w.write_all(line.as_bytes())?;
        pending = Some((e.pos, e.pos, pos));
    }
    // TabixIndexCreator.finalizeIndex takes the file pointer once the last block is written: the start of
    // the next block
    w.deflate_block()?;
    let end = w.file_pointer();
    if let Some((s, en, start)) = pending.take() {
        indexes.last_mut().unwrap().process(s, en, start, end);
    }
    w.finish()?;

    // TabixIndex.write, bgzipped
    let mut idx: Vec<u8> = Vec::new();
    idx.extend_from_slice(b"TBI\x01");
    for v in [names.len() as i32, 2, 1, 2, 0, b'#' as i32, 0] {
        idx.extend_from_slice(&v.to_le_bytes());
    }
    let name_block: usize = names.iter().map(|n| n.len() + 1).sum();
    idx.extend_from_slice(&(name_block as i32).to_le_bytes());
    for n in &names {
        idx.extend_from_slice(n.as_bytes());
        idx.push(0);
    }
    for s in indexes.iter_mut() {
        s.write(&mut idx);
    }
    let tbi = format!("{path}.tbi");
    let file = std::fs::File::create(&tbi).with_context(|| format!("creating {tbi}"))?;
    let mut w = BgzfWriter::new(std::io::BufWriter::new(file));
    w.write_all(&idx)?;
    w.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::AmberBase;

    #[test]
    fn bins() {
        assert_eq!(region_to_bin(17592104, 17592105), 5754);
        assert_eq!(region_to_bin(0, 1), 4681);
        assert_eq!(linear_index_offset(17592105), 1073);
    }

    /// The index answers queries through htslib
    #[test]
    fn vcf_and_index_read_back() {
        let dir = std::env::temp_dir().join(format!("amber_rs_vcf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.vcf.gz").to_str().unwrap().to_string();
        let mut sites = Vec::new();
        for (chr, pos) in [("chr1", 100), ("chr1", 20000), ("chr1", 70000), ("chr2", 5)] {
            let mut e = Evidence::new(chr, pos, AmberBase::A, AmberBase::G);
            e.c.read_depth = 10;
            e.c.ref_support = 6;
            e.c.alt_support = 4;
            sites.push(e);
        }
        let refs: Vec<&Evidence> = sites.iter().collect();
        write_base_depths(&path, &refs, "S").unwrap();
        use rust_htslib::bcf::{IndexedReader, Read};
        let mut r = IndexedReader::from_path(&path).unwrap();
        let rid = r.header().name2rid(b"chr1").unwrap();
        r.fetch(rid, 19000, Some(80000)).unwrap();
        let positions: Vec<i64> = r.records().map(|x| x.unwrap().pos() + 1).collect();
        assert_eq!(positions, vec![20000, 70000]);
        let text = {
            use std::io::Read as _;
            let mut s = String::new();
            flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap()).read_to_string(&mut s).unwrap();
            s
        };
        assert!(text.ends_with("chr2\t5\t.\tA\tG\t.\t.\t.\tGT:AD:DP\t0/1:6,4:10\n"));
        std::fs::remove_dir_all(&dir).ok();
    }
}
