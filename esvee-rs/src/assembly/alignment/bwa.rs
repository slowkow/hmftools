//! BWA-MEM alignment through the vendored bwa library (vendor/bwa), with the interface of the hmf-common
//! BwaMemAligner and the alignment record of gatk-bwamem-jni's BwaMemAlignment.

use std::ffi::{c_char, c_int, c_void, CStr, CString};

use anyhow::{bail, Result};

unsafe extern "C" {
    fn esvee_bwa_open_index(path: *const c_char) -> *mut c_void;
    fn esvee_bwa_on_network_fs(path: *const c_char) -> c_int;
    fn esvee_bwa_destroy_index(idx: *mut c_void) -> c_int;
    fn esvee_bwa_contig_count(idx: *mut c_void) -> c_int;
    fn esvee_bwa_contig_name(idx: *mut c_void, i: c_int) -> *const c_char;
    fn esvee_bwa_contig_length(idx: *mut c_void, i: c_int) -> c_int;
    #[allow(clippy::too_many_arguments)]
    fn esvee_bwa_create_opts(
        match_: c_int, mismatch: c_int, uncertain: c_int, gap_open: c_int, gap_extend: c_int, clip: c_int, min_seed_len: c_int, split_width: c_int,
        max_mem_intv: c_int, max_occ: c_int, split_factor: f32, drop_ratio: f32, band_width: c_int, zdrop: c_int, min_score: c_int, all_alignments: c_int,
    ) -> *mut c_void;
    fn esvee_bwa_destroy_opts(opts: *mut c_void);
    fn esvee_bwa_create_alignments(idx: *mut c_void, opts: *mut c_void, seqs: *mut c_char, buf_size: *mut usize) -> *mut c_void;
    fn esvee_bwa_free(p: *mut c_void);
}

/// The image is memory-mapped, as in gatk-bwamem-jni, and the aligner then touches pages all over it.
/// On a network file system (a bucket mounted with gcsfuse on Google Batch) each of those pages is a
/// separate small read: on a 2M-read-pair test the two alignment passes took 28 minutes at 1% CPU. So
/// there the file is first read through once in large sequential reads, which leaves it in the page
/// cache that the mapping then uses (5.6 GB for GRCh38, about a minute from a bucket). Not a Java
/// behaviour; the output is the same. ESVEE_IMAGE_PRELOAD=1 or 0 forces it on or off.
fn preload_if_remote(path: &str, c_path: &CString) {
    let wanted = match std::env::var("ESVEE_IMAGE_PRELOAD").ok().as_deref() {
        Some("1") => true,
        Some("0") => false,
        _ => (unsafe { esvee_bwa_on_network_fs(c_path.as_ptr()) }) == 1,
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

/// hmf-common BwaMemAlignParams
#[derive(Clone, Copy, Debug)]
pub struct BwaMemAlignParams {
    pub match_reward: i32,
    pub mismatch_penalty: i32,
    pub uncertain_base_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
    pub clip_penalty: i32,
    pub seed_length_min: i32,
    pub seed2_max_occurrence: i32,
    pub seed3_max_occurrence: i32,
    pub mem_max_occurrence: i32,
    pub mem_reseed_factor: f32,
    pub chain_overlap_factor: f32,
    pub band_width: i32,
    pub z_dropoff: i32,
    pub min_align_score: i32,
}

impl BwaMemAlignParams {
    pub const DEFAULT: BwaMemAlignParams = BwaMemAlignParams {
        match_reward: 1,
        mismatch_penalty: 4,
        uncertain_base_penalty: 1,
        gap_open_penalty: 6,
        gap_extend_penalty: 1,
        clip_penalty: 5,
        seed_length_min: 19,
        seed2_max_occurrence: 10,
        seed3_max_occurrence: 20,
        mem_max_occurrence: 500,
        mem_reseed_factor: 1.5,
        chain_overlap_factor: 0.5,
        band_width: 100,
        z_dropoff: 100,
        min_align_score: 30,
    };
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

/// hmf-common IBwaMemAligner
pub trait IBwaMemAligner: Send + Sync {
    fn align_sequences(&self, sequences: &[&[u8]]) -> Vec<Vec<BwaMemAlignment>>;
    fn align_sequence(&self, sequence: &[u8]) -> Vec<BwaMemAlignment> {
        self.align_sequences(&[sequence]).pop().unwrap_or_default()
    }
}

pub struct BwaMemAligner {
    idx: *mut c_void,
    opts: *mut c_void,
}

unsafe impl Send for BwaMemAligner {}
unsafe impl Sync for BwaMemAligner {}

impl BwaMemAligner {
    /// Opens the BWA index image (`<ref genome>.img`) with the given alignment parameters.
    pub fn new(index_path: &str, params: &BwaMemAlignParams) -> Result<BwaMemAligner> {
        Self::new_all_alignments(index_path, params, false)
    }

    /// As `new`, with the BWA "all alignments" option (BwaMemAlignerConfig.allAlignments) on or off.
    pub fn new_all_alignments(index_path: &str, params: &BwaMemAlignParams, all_alignments: bool) -> Result<BwaMemAligner> {
        let c_path = CString::new(index_path)?;
        preload_if_remote(index_path, &c_path);
        let idx = unsafe { esvee_bwa_open_index(c_path.as_ptr()) };
        if idx.is_null() {
            bail!("failed to open BWA index image: {index_path}");
        }
        let opts = unsafe {
            esvee_bwa_create_opts(
                params.match_reward,
                params.mismatch_penalty,
                params.uncertain_base_penalty,
                params.gap_open_penalty,
                params.gap_extend_penalty,
                params.clip_penalty,
                params.seed_length_min,
                params.seed2_max_occurrence,
                params.seed3_max_occurrence,
                params.mem_max_occurrence,
                params.mem_reseed_factor,
                params.chain_overlap_factor,
                params.band_width,
                params.z_dropoff,
                params.min_align_score,
                if all_alignments { 1 } else { 0 },
            )
        };
        Ok(BwaMemAligner { idx, opts })
    }

    /// reference contigs of the index in order: (name, length)
    pub fn contigs(&self) -> Vec<(String, i32)> {
        let n = unsafe { esvee_bwa_contig_count(self.idx) };
        (0..n)
            .map(|i| unsafe { (CStr::from_ptr(esvee_bwa_contig_name(self.idx, i)).to_string_lossy().into_owned(), esvee_bwa_contig_length(self.idx, i)) })
            .collect()
    }
}

impl Drop for BwaMemAligner {
    fn drop(&mut self) {
        unsafe {
            esvee_bwa_destroy_opts(self.opts);
            esvee_bwa_destroy_index(self.idx);
        }
    }
}

impl IBwaMemAligner for BwaMemAligner {
    fn align_sequences(&self, sequences: &[&[u8]]) -> Vec<Vec<BwaMemAlignment>> {
        // a 4-byte count, then each sequence with a trailing NUL
        let mut buf: Vec<u8> = Vec::with_capacity(4 + sequences.iter().map(|s| s.len() + 1).sum::<usize>());
        buf.extend_from_slice(&(sequences.len() as u32).to_ne_bytes());
        for s in sequences {
            buf.extend_from_slice(s);
            buf.push(0);
        }
        let mut size: usize = 0;
        let results = unsafe { esvee_bwa_create_alignments(self.idx, self.opts, buf.as_mut_ptr() as *mut c_char, &mut size) };
        let ints: &[i32] = if size == 0 || results.is_null() { &[] } else { unsafe { std::slice::from_raw_parts(results as *const i32, size / 4) } };
        let parsed = parse_alignments(ints, sequences.len());
        unsafe { esvee_bwa_free(results) };
        parsed
    }
}

/// Decodes the packed result buffer as gatk-bwamem-jni's BwaMemAligner.alignSeqs does.
fn parse_alignments(buf: &[i32], n_sequences: usize) -> Vec<Vec<BwaMemAlignment>> {
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
