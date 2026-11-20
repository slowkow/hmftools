//! BAM output: block-compressed record streams written per partition and joined at the end without recompression.

use anyhow::{bail, Result};
use rust_htslib::bam;
use rust_htslib::htslib;
use std::ffi::CString;
use std::io::{Read, Seek, SeekFrom, Write};

/// A BGZF stream of BAM records, with or without the BAM header.
pub struct RawBamOut {
    fp: *mut htslib::BGZF,
    pub path: String,
    pub written: u64,
    /// index built while writing (see `enable_index`); offsets are relative to this file
    idx: *mut htslib::hts_idx_t,
    idx_failed: bool,
}

unsafe impl Send for RawBamOut {}

/// The 28-byte empty block that ends a BGZF file.
pub const BGZF_EOF_LEN: u64 = 28;

impl RawBamOut {
    pub fn create(path: &str, level: u32) -> Result<Self> {
        let cpath = CString::new(path)?;
        let mode = CString::new(format!("w{level}"))?;
        let fp = unsafe { htslib::bgzf_open(cpath.as_ptr(), mode.as_ptr()) };
        if fp.is_null() {
            bail!("failed to open {path} for writing");
        }
        Ok(RawBamOut { fp, path: path.to_string(), written: 0, idx: std::ptr::null_mut(), idx_failed: false })
    }
    pub fn write_header(&mut self, header: &bam::HeaderView) -> Result<()> {
        if unsafe { htslib::bam_hdr_write(self.fp, header.inner_ptr()) } < 0 {
            bail!("failed to write the BAM header to {}", self.path);
        }
        Ok(())
    }
    #[inline]
    pub fn write(&mut self, rec: &bam::Record) -> Result<()> {
        if unsafe { htslib::bam_write1(self.fp, &rec.inner) } < 0 {
            bail!("failed to write a record to {}", self.path);
        }
        self.written += 1;
        if !self.idx.is_null() && !self.idx_failed {
            let b = &rec.inner;
            let r = unsafe {
                htslib::hts_idx_push(self.idx, b.core.tid, b.core.pos, htslib::bam_endpos(b), self.tell(), (b.core.flag & 4 == 0) as i32)
            };
            // records out of order: no index from this file, the caller indexes the final BAM instead
            self.idx_failed = r < 0;
        }
        Ok(())
    }
    /// The virtual offset of the next byte (bgzf_tell)
    fn tell(&self) -> u64 {
        unsafe { (((*self.fp).block_address as u64) << 16) | ((*self.fp).block_offset as u64 & 0xFFFF) }
    }
    /// Builds a .bai for this file while it is written, saved as `<path>.bai` at `close`. The file must
    /// hold only records, in coordinate order; `merge_part_indexes` turns the indexes of the parts into
    /// the index of the final BAM, which saves reading the final BAM again.
    pub fn enable_index(&mut self) {
        self.idx = unsafe { htslib::hts_idx_init(0, htslib::HTS_FMT_BAI as i32, 0, 14, 5) };
    }
    /// Ends the current block, so that the next record starts a new one.
    pub fn flush(&mut self) -> Result<()> {
        if unsafe { htslib::bgzf_flush(self.fp) } < 0 {
            bail!("failed to flush {}", self.path);
        }
        Ok(())
    }
    /// Compressed file offset; valid after `flush`.
    pub fn compressed_offset(&self) -> u64 {
        (unsafe { (*self.fp).block_address }) as u64
    }
    /// Appends already compressed blocks.
    pub fn write_raw(&mut self, data: &[u8]) -> Result<()> {
        let n = unsafe { htslib::bgzf_raw_write(self.fp, data.as_ptr() as *const std::ffi::c_void, data.len()) };
        if n != data.len() as isize {
            bail!("failed to write to {}", self.path);
        }
        Ok(())
    }
    pub fn close(mut self) -> Result<()> {
        if !self.idx.is_null() {
            let idx = std::mem::replace(&mut self.idx, std::ptr::null_mut());
            let ok = !self.idx_failed && unsafe { htslib::bgzf_flush(self.fp) } >= 0 && unsafe { htslib::hts_idx_finish(idx, self.tell()) } >= 0;
            if ok {
                let cpath = CString::new(self.path.as_str())?;
                let cidx = CString::new(format!("{}.bai", self.path))?;
                unsafe { htslib::hts_idx_save_as(idx, cpath.as_ptr(), cidx.as_ptr(), htslib::HTS_FMT_BAI as i32) };
            }
            unsafe { htslib::hts_idx_destroy(idx) };
        }
        let fp = std::mem::replace(&mut self.fp, std::ptr::null_mut());
        if unsafe { htslib::bgzf_close(fp) } < 0 {
            bail!("failed to close {}", self.path);
        }
        Ok(())
    }
}

impl Drop for RawBamOut {
    fn drop(&mut self) {
        if !self.idx.is_null() {
            unsafe { htslib::hts_idx_destroy(self.idx) };
        }
        if !self.fp.is_null() {
            unsafe { htslib::bgzf_close(self.fp) };
        }
    }
}

/// Appends bytes `[start, end)` of a block-compressed file to `out`.
pub fn append_file_range(out: &mut RawBamOut, path: &str, start: u64, end: u64) -> Result<()> {
    if end <= start {
        return Ok(());
    }
    let mut f = std::fs::File::open(path)?;
    f.seek(SeekFrom::Start(start))?;
    let mut remaining = end - start;
    let mut buf = vec![0u8; 4 << 20];
    while remaining > 0 {
        let n = (remaining as usize).min(buf.len());
        f.read_exact(&mut buf[..n])?;
        out.write_raw(&buf[..n])?;
        remaining -= n as u64;
    }
    Ok(())
}

/// Appends a whole record stream written by `RawBamOut`, without its end-of-file block.
pub fn append_file(out: &mut RawBamOut, path: &str) -> Result<()> {
    let len = std::fs::metadata(path)?.len();
    if len < BGZF_EOF_LEN {
        bail!("{path} is not a complete block-compressed file");
    }
    append_file_range(out, path, 0, len - BGZF_EOF_LEN)
}

/// One reference of a .bai: bins with their chunks, the mapped/unmapped pseudo-bin, the linear index
#[derive(Default)]
struct BaiRef {
    bins: std::collections::BTreeMap<u32, Vec<(u64, u64)>>,
    /// file range of the reference, then mapped and unmapped read counts
    meta: Option<[u64; 4]>,
    linear: Vec<u64>,
}

const BAI_META_BIN: u32 = 37450;

fn read_bai(path: &str) -> Result<Vec<BaiRef>> {
    let data = std::fs::read(path)?;
    let mut at = 0usize;
    let mut take = |n: usize| -> Result<&[u8]> {
        if at + n > data.len() {
            bail!("{path} is truncated");
        }
        at += n;
        Ok(&data[at - n..at])
    };
    if take(4)? != b"BAI\x01" {
        bail!("{path} is not a BAM index");
    }
    let u32_of = |b: &[u8]| u32::from_le_bytes(b.try_into().unwrap());
    let u64_of = |b: &[u8]| u64::from_le_bytes(b.try_into().unwrap());
    let n_ref = u32_of(take(4)?) as usize;
    let mut refs = Vec::with_capacity(n_ref);
    for _ in 0..n_ref {
        let mut r = BaiRef::default();
        for _ in 0..u32_of(take(4)?) {
            let bin = u32_of(take(4)?);
            let n_chunk = u32_of(take(4)?) as usize;
            let mut chunks = Vec::with_capacity(n_chunk);
            for _ in 0..n_chunk {
                chunks.push((u64_of(take(8)?), u64_of(take(8)?)));
            }
            if bin == BAI_META_BIN {
                if chunks.len() != 2 {
                    bail!("{path}: unexpected pseudo-bin");
                }
                r.meta = Some([chunks[0].0, chunks[0].1, chunks[1].0, chunks[1].1]);
            } else {
                r.bins.insert(bin, chunks);
            }
        }
        for _ in 0..u32_of(take(4)?) {
            r.linear.push(u64_of(take(8)?));
        }
        refs.push(r);
    }
    Ok(refs)
}

/// Writes the index of a BAM made by appending `parts` (path of a part's own index, and the file offset
/// at which that part starts in the BAM), in order. The parts together are coordinate sorted.
pub fn merge_part_indexes(parts: &[(String, u64)], n_ref: usize, no_coor: u64, out_path: &str) -> Result<()> {
    let mut merged: Vec<BaiRef> = (0..n_ref).map(|_| BaiRef::default()).collect();
    for (path, start) in parts {
        let shift = |v: u64| v + (start << 16);
        for (tid, r) in read_bai(path)?.into_iter().enumerate() {
            if tid >= n_ref {
                bail!("{path} has more references than the BAM");
            }
            let m = &mut merged[tid];
            for (bin, chunks) in r.bins {
                m.bins.entry(bin).or_default().extend(chunks.into_iter().map(|(a, b)| (shift(a), shift(b))));
            }
            if let Some(x) = r.meta {
                let (a, b) = (shift(x[0]), shift(x[1]));
                m.meta = Some(match m.meta {
                    Some(y) => [y[0].min(a), y[1].max(b), y[2] + x[2], y[3] + x[3]],
                    None => [a, b, x[2], x[3]],
                });
            }
            // a window already covered by an earlier part keeps that part's offset, which is the lower one
            for (w, v) in r.linear.into_iter().enumerate() {
                if w >= m.linear.len() {
                    m.linear.push(shift(v));
                }
            }
        }
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(out_path)?);
    w.write_all(b"BAI\x01")?;
    w.write_all(&(n_ref as u32).to_le_bytes())?;
    for m in &merged {
        w.write_all(&((m.bins.len() + m.meta.is_some() as usize) as u32).to_le_bytes())?;
        for (bin, chunks) in &m.bins {
            w.write_all(&bin.to_le_bytes())?;
            w.write_all(&(chunks.len() as u32).to_le_bytes())?;
            for (a, b) in chunks {
                w.write_all(&a.to_le_bytes())?;
                w.write_all(&b.to_le_bytes())?;
            }
        }
        if let Some(x) = m.meta {
            w.write_all(&BAI_META_BIN.to_le_bytes())?;
            w.write_all(&2u32.to_le_bytes())?;
            for v in x {
                w.write_all(&v.to_le_bytes())?;
            }
        }
        w.write_all(&(m.linear.len() as u32).to_le_bytes())?;
        for v in &m.linear {
            w.write_all(&v.to_le_bytes())?;
        }
    }
    w.write_all(&no_coor.to_le_bytes())?;
    w.flush()?;
    Ok(())
}

/// Builds the .bai index of a coordinate-sorted BAM.
pub fn index_bam(path: &str, threads: usize) -> Result<()> {
    let cpath = CString::new(path)?;
    let r = unsafe { htslib::sam_index_build3(cpath.as_ptr(), std::ptr::null(), 0, threads as i32) };
    if r < 0 {
        bail!("failed to index {path}");
    }
    Ok(())
}

/// The output header: the first input's header as coordinate sorted, with the read groups and the first program
/// record of the other inputs added (BamUtils.buildCombinedHeader).
pub fn build_combined_header(bam_files: &[String], ref_genome: &str) -> Result<bam::Header> {
    use rust_htslib::bam::Read as _;
    let mut texts = Vec::new();
    for f in bam_files {
        let mut r = bam::Reader::from_path(f)?;
        r.set_reference(ref_genome)?;
        texts.push(String::from_utf8_lossy(r.header().as_bytes()).to_string());
    }
    let mut hd: Option<String> = None;
    let (mut sq, mut rg, mut pg, mut co) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let id_of = |line: &str| line.split('\t').find_map(|f| f.strip_prefix("ID:")).map(str::to_string);
    for line in texts[0].lines() {
        match line.get(..3) {
            Some("@HD") => hd = Some(line.to_string()),
            Some("@SQ") => sq.push(line.to_string()),
            Some("@RG") => rg.push(line.to_string()),
            Some("@PG") => pg.push(line.to_string()),
            Some("@CO") => co.push(line.to_string()),
            _ => {}
        }
    }
    for (i, text) in texts.iter().enumerate().skip(1) {
        let mut first_pg = true;
        for line in text.lines() {
            if line.starts_with("@RG") {
                let id = id_of(line);
                if !rg.iter().any(|l| id_of(l) == id) {
                    rg.push(line.to_string());
                }
            } else if line.starts_with("@PG") && first_pg {
                first_pg = false;
                let Some(id) = id_of(line) else { continue };
                let new_id = format!("{id}.{i}");
                if !pg.iter().any(|l| id_of(l).as_deref() == Some(&new_id)) {
                    let fields: Vec<String> =
                        line.split('\t').map(|f| if f.starts_with("ID:") { format!("ID:{new_id}") } else { f.to_string() }).collect();
                    pg.push(fields.join("\t"));
                }
            }
        }
    }
    let hd_line = match hd {
        Some(l) => {
            let mut fields: Vec<String> = l.split('\t').filter(|f| !f.starts_with("SO:")).map(str::to_string).collect();
            fields.insert(fields.len().min(2), "SO:coordinate".to_string());
            fields.join("\t")
        }
        None => "@HD\tVN:1.6\tSO:coordinate".to_string(),
    };
    let mut text = hd_line;
    text.push('\n');
    for l in sq.iter().chain(rg.iter()).chain(pg.iter()).chain(co.iter()) {
        text.push_str(l);
        text.push('\n');
    }
    let mut header = bam::Header::new();
    // Header::new() is empty; the text is kept verbatim
    let _ = &mut header;
    Ok(bam::Header::from_template(&bam::HeaderView::from_bytes(text.as_bytes())))
}

/// Writes text to a file through a buffer.
pub fn write_text(path: &str, text: &str) -> Result<()> {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    f.write_all(text.as_bytes())?;
    Ok(())
}
