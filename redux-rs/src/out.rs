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
        Ok(RawBamOut { fp, path: path.to_string(), written: 0 })
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
        Ok(())
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
        let fp = std::mem::replace(&mut self.fp, std::ptr::null_mut());
        if unsafe { htslib::bgzf_close(fp) } < 0 {
            bail!("failed to close {}", self.path);
        }
        Ok(())
    }
}

impl Drop for RawBamOut {
    fn drop(&mut self) {
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
