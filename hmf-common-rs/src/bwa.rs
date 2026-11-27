//! BWA-MEM alignment through the vendored bwa library (vendor/bwa, built with the `bwa` feature), as the
//! Java tools align through gatk-bwamem-jni: an index image (`<ref genome>.img`, or one held in memory),
//! the `mem_opt_t` option block of `mem_opt_init()` changed field by field (`MemOpt`, the struct the Java
//! BwaMemAligner writes through a direct ByteBuffer), and the alignment records of BwaMemAlignment.
//!
//! The C code is esvee-rs's (`esvee_bwa.c`, here `hmf_bwa.c`). esvee-rs still builds its own copy; it can
//! take this module instead (its BwaMemAlignParams map onto `MemOpt` plus `set_scoring_matrix`).

use std::ffi::{c_char, c_int, c_void, CStr, CString};

use anyhow::{bail, Result};

unsafe extern "C" {
    fn hmf_bwa_open_index(path: *const c_char) -> *mut c_void;
    fn hmf_bwa_open_index_mem(mem: *mut u8, size: i64) -> *mut c_void;
    fn hmf_bwa_on_network_fs(path: *const c_char) -> c_int;
    fn hmf_bwa_destroy_index(idx: *mut c_void, mapped: c_int) -> c_int;
    fn hmf_bwa_contig_count(idx: *mut c_void) -> c_int;
    fn hmf_bwa_contig_name(idx: *mut c_void, i: c_int) -> *const c_char;
    fn hmf_bwa_contig_length(idx: *mut c_void, i: c_int) -> c_int;
    fn hmf_bwa_default_opts() -> *mut MemOpt;
    fn hmf_bwa_destroy_opts(opts: *mut MemOpt);
    fn hmf_bwa_fill_scmat(a: c_int, b: c_int, mat: *mut i8);
    fn hmf_bwa_opts_layout(which: c_int) -> usize;
    fn hmf_bwa_create_alignments(idx: *mut c_void, opts: *mut MemOpt, seqs: *mut c_char, buf_size: *mut usize) -> *mut c_void;
    fn hmf_bwa_free(p: *mut c_void);
}

/// bwa 0.7.16a's mem_opt_t, field for field (the layout is checked against the C struct when options are
/// created). The Java BwaMemAligner setters write these fields at their byte offsets.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct MemOpt {
    /// match score (setMatchScoreOption)
    pub a: i32,
    /// mismatch penalty (setMismatchPenaltyOption)
    pub b: i32,
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    pub pen_unpaired: i32,
    pub pen_clip5: i32,
    pub pen_clip3: i32,
    /// band width
    pub w: i32,
    pub zdrop: i32,
    pub max_mem_intv: u64,
    /// output score threshold
    pub t: i32,
    /// MEM_F_* flags
    pub flag: i32,
    pub min_seed_len: i32,
    pub min_chain_weight: i32,
    pub max_chain_extend: i32,
    pub split_factor: f32,
    pub split_width: i32,
    pub max_occ: i32,
    pub max_chain_gap: i32,
    pub n_threads: i32,
    pub chunk_size: i32,
    pub mask_level: f32,
    pub drop_ratio: f32,
    pub xa_drop_ratio: f32,
    pub mask_level_redun: f32,
    pub map_q_coef_len: f32,
    pub map_q_coef_fac: i32,
    pub max_ins: i32,
    pub max_matesw: i32,
    pub max_xa_hits: i32,
    pub max_xa_hits_alt: i32,
    /// the 5x5 scoring matrix; mem_opt_init fills it from a and b (bwa_fill_scmat), and the Java setters
    /// of a and b leave it as it is
    pub mat: [i8; 25],
}

/// BwaMemAligner.MEM_F_ALL: output all alignments
pub const MEM_F_ALL: i32 = 0x8;

/// An option block from mem_opt_init(), to be changed through `opts_mut`.
pub struct BwaOptions {
    ptr: *mut MemOpt,
}

unsafe impl Send for BwaOptions {}
unsafe impl Sync for BwaOptions {}

impl BwaOptions {
    /// mem_opt_init(), as BwaMemIndex.createDefaultOptions
    pub fn new() -> Result<BwaOptions> {
        unsafe {
            let layout = [
                std::mem::size_of::<MemOpt>(),
                std::mem::offset_of!(MemOpt, max_mem_intv),
                std::mem::offset_of!(MemOpt, n_threads),
                std::mem::offset_of!(MemOpt, drop_ratio),
                std::mem::offset_of!(MemOpt, mat),
            ];
            for (i, &v) in layout.iter().enumerate() {
                let c = hmf_bwa_opts_layout(i as c_int);
                if c != v {
                    bail!("mem_opt_t layout differs from the Rust mirror (item {i}: C {c}, Rust {v})");
                }
            }
            let ptr = hmf_bwa_default_opts();
            if ptr.is_null() {
                bail!("mem_opt_init failed");
            }
            Ok(BwaOptions { ptr })
        }
    }
    pub fn opts(&self) -> &MemOpt {
        unsafe { &*self.ptr }
    }
    pub fn opts_mut(&mut self) -> &mut MemOpt {
        unsafe { &mut *self.ptr }
    }
    /// The matrix as hmf-common's BwaMemAligner sets it: match and mismatch scores, and `uncertain` as the
    /// penalty for any N.
    pub fn set_scoring_matrix(&mut self, match_score: i32, mismatch: i32, uncertain: i32) {
        let o = self.opts_mut();
        let mut k = 0;
        for i in 0..4 {
            for j in 0..4 {
                o.mat[k] = if i == j { match_score as i8 } else { -mismatch as i8 };
                k += 1;
            }
            o.mat[k] = -uncertain as i8;
            k += 1;
        }
        for _ in 0..5 {
            o.mat[k] = -uncertain as i8;
            k += 1;
        }
    }
    /// bwa_fill_scmat(a, b): the matrix of `bwa mem` (-1 for any N)
    pub fn fill_scoring_matrix(&mut self, match_score: i32, mismatch: i32) {
        unsafe { hmf_bwa_fill_scmat(match_score, mismatch, self.opts_mut().mat.as_mut_ptr()) };
    }
}

impl Drop for BwaOptions {
    fn drop(&mut self) {
        unsafe { hmf_bwa_destroy_opts(self.ptr) };
    }
}

/// The image is memory-mapped, as in gatk-bwamem-jni, and the aligner then touches pages all over it.
/// On a network file system (a bucket mounted with gcsfuse on Google Batch) each of those pages is a
/// separate small read, so there the file is first read through once in large sequential reads, which
/// leaves it in the page cache that the mapping then uses. Not a Java behaviour; the output is the same.
/// `<env_var>`=1 or 0 forces it on or off (esvee-rs: ESVEE_IMAGE_PRELOAD).
fn preload_if_remote(path: &str, c_path: &CString, env_var: &str, always: bool) {
    let wanted = match std::env::var(env_var).ok().as_deref() {
        Some("1") => true,
        Some("0") => false,
        _ => always || (unsafe { hmf_bwa_on_network_fs(c_path.as_ptr()) }) == 1,
    };
    if !wanted {
        return;
    }
    let start = std::time::Instant::now();
    let read = (|| -> std::io::Result<u64> {
        use std::io::Read;
        let mut file = std::fs::File::open(path)?;
        let mut buf = vec![0u8; 16 << 20];
        let mut total = 0u64;
        loop {
            let n = file.read(&mut buf)?;
            if n == 0 {
                return Ok(total);
            }
            total += n as u64;
        }
    })();
    match read {
        Ok(n) => crate::debug!("read the genome image {path} through once ({:.1} GB, {:.1}s)", n as f64 / 1e9, start.elapsed().as_secs_f64()),
        // not fatal: the mapping still works, only more slowly
        Err(e) => crate::debug!("could not read the genome image {path} through: {e}"),
    }
}

/// BwaMemIndex: an index image, mapped from a file or copied into memory.
pub struct BwaIndex {
    idx: *mut c_void,
    mapped: bool,
    // the image of `from_image_bytes`, 8-byte aligned; the index points into it
    _mem: Option<Vec<u64>>,
}

unsafe impl Send for BwaIndex {}
unsafe impl Sync for BwaIndex {}

impl BwaIndex {
    /// Maps the image file (BwaMemIndex(path)); `preload_env` names the variable that forces reading it
    /// through first (see preload_if_remote), which happens by default on network file systems only.
    pub fn open(path: &str, preload_env: &str) -> Result<BwaIndex> {
        Self::open_with_preload(path, preload_env, false)
    }

    /// As `open`; with `always` the image is read through first wherever it is (unless `preload_env` is 0).
    /// On a Mac with little free memory, page faults into the mapped GRCh38 image made a CIDER alignment
    /// of a few hundred sequences take 80 to 180 s, against 6 s after a read-through (1 s when cached).
    pub fn open_with_preload(path: &str, preload_env: &str, always: bool) -> Result<BwaIndex> {
        let c_path = CString::new(path)?;
        if !std::path::Path::new(path).is_file() {
            bail!("BWA index image not found: {path}");
        }
        preload_if_remote(path, &c_path, preload_env, always);
        let idx = unsafe { hmf_bwa_open_index(c_path.as_ptr()) };
        if idx.is_null() {
            bail!("failed to open BWA index image: {path}");
        }
        Ok(BwaIndex { idx, mapped: true, _mem: None })
    }

    /// An image held in memory (for example a resource compiled into the program), copied so that it is
    /// aligned as the index needs.
    pub fn from_image_bytes(image: &[u8]) -> Result<BwaIndex> {
        if image.len() < 64 {
            bail!("BWA index image too short");
        }
        let mut mem = vec![0u64; image.len().div_ceil(8)];
        unsafe { std::ptr::copy_nonoverlapping(image.as_ptr(), mem.as_mut_ptr() as *mut u8, image.len()) };
        let idx = unsafe { hmf_bwa_open_index_mem(mem.as_mut_ptr() as *mut u8, image.len() as i64) };
        if idx.is_null() {
            bail!("failed to open BWA index image from memory");
        }
        Ok(BwaIndex { idx, mapped: false, _mem: Some(mem) })
    }

    /// reference contigs of the index in order: (name, length)
    pub fn contigs(&self) -> Vec<(String, i32)> {
        let n = unsafe { hmf_bwa_contig_count(self.idx) };
        (0..n)
            .map(|i| unsafe { (CStr::from_ptr(hmf_bwa_contig_name(self.idx, i)).to_string_lossy().into_owned(), hmf_bwa_contig_length(self.idx, i)) })
            .collect()
    }

    /// BwaMemAligner.alignSeqs: the alignments of each sequence, in order. The option block's n_threads
    /// sets the threads BWA uses inside the call; the results do not depend on it.
    pub fn align(&self, opts: &BwaOptions, sequences: &[&[u8]]) -> Vec<Vec<BwaMemAlignment>> {
        if sequences.is_empty() {
            return Vec::new();
        }
        // a 4-byte count, then each sequence with a trailing NUL
        let mut buf: Vec<u8> = Vec::with_capacity(4 + sequences.iter().map(|s| s.len() + 1).sum::<usize>());
        buf.extend_from_slice(&(sequences.len() as u32).to_ne_bytes());
        for s in sequences {
            buf.extend_from_slice(s);
            buf.push(0);
        }
        let mut size: usize = 0;
        let results = unsafe { hmf_bwa_create_alignments(self.idx, opts.ptr, buf.as_mut_ptr() as *mut c_char, &mut size) };
        let ints: &[i32] = if size == 0 || results.is_null() { &[] } else { unsafe { std::slice::from_raw_parts(results as *const i32, size / 4) } };
        let parsed = parse_alignments(ints, sequences.len());
        unsafe { hmf_bwa_free(results) };
        parsed
    }
}

impl Drop for BwaIndex {
    fn drop(&mut self) {
        unsafe {
            hmf_bwa_destroy_index(self.idx, if self.mapped { 1 } else { 0 });
        }
    }
}

/// gatk-bwamem-jni BwaMemAlignment
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BwaMemAlignment {
    pub sam_flag: i32,
    pub ref_id: i32,
    /// 0-based inclusive
    pub ref_start: i32,
    /// 0-based exclusive
    pub ref_end: i32,
    pub seq_start: i32,
    pub seq_end: i32,
    pub map_qual: i32,
    /// the edit distance of the aligned part (NM)
    pub n_mismatches: i32,
    pub aligner_score: i32,
    pub suboptimal_score: i32,
    pub cigar: String,
    pub md_tag: Option<String>,
    pub xa_tag: Option<String>,
    pub mate_ref_id: i32,
    pub mate_ref_start: i32,
    pub template_len: i32,
}

/// Decodes the packed result buffer as gatk-bwamem-jni's BwaMemAligner.alignSeqs does.
pub fn parse_alignments(buf: &[i32], n_sequences: usize) -> Vec<Vec<BwaMemAlignment>> {
    const CIGAR_OPS: &[u8; 16] = b"MID?S???????????";
    struct Cursor<'a> {
        buf: &'a [i32],
        pos: usize,
    }
    impl Cursor<'_> {
        fn next(&mut self) -> i32 {
            let v = self.buf[self.pos];
            self.pos += 1;
            v
        }
        fn tag(&mut self) -> Option<String> {
            let tag_len = self.next() as usize;
            if tag_len == 0 {
                return None;
            }
            let words = tag_len.div_ceil(4);
            let bytes: Vec<u8> = self.buf[self.pos..self.pos + words].iter().flat_map(|w| w.to_ne_bytes()).take(tag_len).collect();
            self.pos += words;
            Some(String::from_utf8_lossy(&bytes).into_owned())
        }
    }
    let mut c = Cursor { buf, pos: 0 };
    let mut all = Vec::with_capacity(n_sequences);
    for _ in 0..n_sequences {
        let n_aligns = c.next();
        let mut alignments = Vec::with_capacity(n_aligns.max(0) as usize);
        for _ in 0..n_aligns {
            let flag_map_q = c.next();
            let flags = ((flag_map_q as u32) >> 16) as i32;
            let map_qual = flag_map_q & 0xff;
            let (ref_id, ref_start, ref_end, seq_start, seq_end, n_mismatches, aligner_score, suboptimal_score);
            let mut cigar = String::new();
            let (md_tag, xa_tag);
            if flags & 0x4 != 0 {
                ref_id = -1;
                ref_start = -1;
                ref_end = -1;
                seq_start = -1;
                seq_end = -1;
                n_mismatches = 0;
                aligner_score = 0;
                suboptimal_score = 0;
                md_tag = None;
                xa_tag = None;
            } else {
                ref_id = c.next();
                ref_start = c.next();
                n_mismatches = c.next();
                aligner_score = c.next();
                suboptimal_score = c.next();
                let n_cigar_ops = c.next();
                if n_cigar_ops <= 0 {
                    seq_start = 0;
                    seq_end = 0;
                    ref_end = ref_start;
                } else {
                    let mut ref_len = 0;
                    let mut seq_len = 0;
                    let mut start = 0;
                    for i in 0..n_cigar_ops {
                        let len_op = c.next() as u32;
                        let len = (len_op >> 4) as i32;
                        let op = CIGAR_OPS[(len_op & 0x0f) as usize];
                        cigar.push_str(&len.to_string());
                        cigar.push(op as char);
                        if i == 0 {
                            start = if op == b'S' { len } else { 0 };
                        }
                        if op == b'M' || op == b'D' {
                            ref_len += len;
                        }
                        if op == b'M' || op == b'I' {
                            seq_len += len;
                        }
                    }
                    seq_start = start;
                    ref_end = ref_start + ref_len;
                    seq_end = seq_start + seq_len;
                }
                md_tag = c.tag();
                xa_tag = c.tag();
            }
            let (mate_ref_id, mate_ref_start, template_len);
            if flags & 0x1 == 0 || flags & 0x8 != 0 {
                mate_ref_id = -1;
                mate_ref_start = -1;
                template_len = 0;
            } else {
                mate_ref_id = c.next();
                mate_ref_start = c.next();
                template_len = c.next();
            }
            alignments.push(BwaMemAlignment {
                sam_flag: flags,
                ref_id,
                ref_start,
                ref_end,
                seq_start,
                seq_end,
                map_qual,
                n_mismatches,
                aligner_score,
                suboptimal_score,
                cigar,
                md_tag,
                xa_tag,
                mate_ref_id,
                mate_ref_start,
                template_len,
            });
        }
        all.push(alignments);
    }
    all
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_options_are_mem_opt_init() {
        let o = BwaOptions::new().unwrap();
        let m = o.opts();
        assert_eq!((m.a, m.b, m.o_del, m.e_del, m.o_ins, m.e_ins), (1, 4, 6, 1, 6, 1));
        assert_eq!((m.w, m.t, m.zdrop, m.pen_clip5, m.pen_clip3, m.min_seed_len), (100, 30, 100, 5, 5, 19));
        assert_eq!((m.max_mem_intv, m.split_width, m.max_occ, m.max_chain_gap, m.n_threads), (20, 10, 500, 10000, 1));
        assert_eq!((m.split_factor, m.drop_ratio, m.mask_level), (1.5, 0.5, 0.5));
        // bwa_fill_scmat(1, 4): 1 on the diagonal, -4 off it, -1 for N
        assert_eq!(&m.mat[..6], &[1, -4, -4, -4, -1, -4]);
        assert_eq!(&m.mat[20..], &[-1, -1, -1, -1, -1]);
    }

    #[test]
    fn matrix_with_uncertain_penalty() {
        let mut o = BwaOptions::new().unwrap();
        o.set_scoring_matrix(2, 3, 1);
        assert_eq!(&o.opts().mat[..10], &[2, -3, -3, -3, -1, -3, 2, -3, -3, -1]);
        o.fill_scoring_matrix(1, 4);
        assert_eq!(&o.opts().mat[..5], &[1, -4, -4, -4, -1]);
    }

    #[test]
    fn parse_unmapped_and_mapped() {
        // one sequence with an unmapped record, one with a mapped one: rid 2, pos 10, NM 1, AS 30, XS 0,
        // cigar 5S20M, no MD, no XA
        let buf = vec![1, 4 << 16, 1, 0, 2, 10, 1, 30, 0, 2, (5 << 4) | 4, 20 << 4, 0, 0];
        let a = parse_alignments(&buf, 2);
        assert_eq!(a[0][0].sam_flag, 4);
        assert_eq!(a[0][0].ref_id, -1);
        let m = &a[1][0];
        assert_eq!((m.ref_id, m.ref_start, m.ref_end, m.seq_start, m.seq_end), (2, 10, 30, 5, 25));
        assert_eq!((m.n_mismatches, m.aligner_score, m.cigar.as_str()), (1, 30, "5S20M"));
    }
}
