//! Indexed FASTA through htslib's faidx, in place of rust-htslib's `faidx::Reader` (same names and
//! results).
//!
//! rust-htslib 1.0.2 `Reader::fetch_seq` hands the buffer that htslib allocated with C `malloc` to a Rust
//! `Vec`, which the Rust global allocator then frees. With mimalloc as the global allocator that crashes
//! on Linux (seen in the unit tests); with the system allocator it works by accident. It also does not
//! check for a failed fetch (a NULL buffer), nor `fai_load` for a failed load. Here the bases are copied
//! out and htslib's buffer is returned with `hts_free`, and both failures are errors.

use anyhow::{bail, Result};
use rust_htslib::htslib;
use std::ffi::{CStr, CString};
use std::path::Path;

pub struct Reader {
    inner: *mut htslib::faidx_t,
}

// as rust-htslib's Reader: one thread at a time (htslib's faidx keeps file state between fetches)
unsafe impl Send for Reader {}

impl Reader {
    /// fai_load: reads `<path>.fai` (and `.gzi` for a bgzipped file), building the index if absent
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        if !path.exists() {
            bail!("file not found: {}", path.display());
        }
        let Some(text) = path.to_str() else { bail!("path is not valid UTF-8: {}", path.display()) };
        let cpath = CString::new(text)?;
        let inner = unsafe { htslib::fai_load(cpath.as_ptr()) };
        if inner.is_null() {
            bail!("failed to load the FASTA index of {text}");
        }
        Ok(Reader { inner })
    }

    /// Bases from `begin` to `end`, zero-based and inclusive, cut at the contig end (faidx_fetch_seq64)
    pub fn fetch_seq<N: AsRef<str>>(&self, name: N, begin: usize, end: usize) -> Result<Vec<u8>> {
        if begin > i64::MAX as usize || end > i64::MAX as usize {
            bail!("FASTA position too large: {begin}-{end}");
        }
        let cname = CString::new(name.as_ref().as_bytes())?;
        let mut len: htslib::hts_pos_t = 0;
        let ptr = unsafe { htslib::faidx_fetch_seq64(self.inner, cname.as_ptr(), begin as htslib::hts_pos_t, end as htslib::hts_pos_t, &mut len) };
        if ptr.is_null() {
            bail!("failed to fetch {}:{}-{} from the FASTA (htslib code {len})", name.as_ref(), begin, end);
        }
        let bases = if len > 0 { unsafe { std::slice::from_raw_parts(ptr as *const u8, len as usize) }.to_vec() } else { Vec::new() };
        unsafe { htslib::hts_free(ptr as *mut std::ffi::c_void) };
        Ok(bases)
    }

    /// faidx_seq_len, as a u64 as in rust-htslib (an unknown name gives u64::MAX)
    pub fn fetch_seq_len<N: AsRef<str>>(&self, name: N) -> u64 {
        let cname = CString::new(name.as_ref().as_bytes()).unwrap();
        (unsafe { htslib::faidx_seq_len(self.inner, cname.as_ptr()) }) as u64
    }

    /// The contig names in index order
    pub fn seq_names(&self) -> Result<Vec<String>> {
        let n = unsafe { htslib::faidx_nseq(self.inner) };
        let mut names = Vec::with_capacity(n.max(0) as usize);
        for i in 0..n {
            let name = unsafe { CStr::from_ptr(htslib::faidx_iseq(self.inner, i)) };
            names.push(name.to_str()?.to_string());
        }
        Ok(names)
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        unsafe { htslib::fai_destroy(self.inner) };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fasta(name: &str, text: &str) -> (std::path::PathBuf, String) {
        let dir = std::env::temp_dir().join(format!("faidx_{name}_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.fa");
        std::fs::write(&path, text).unwrap();
        let p = path.to_str().unwrap().to_string();
        (dir, p)
    }

    #[test]
    fn same_as_rust_htslib() {
        let (dir, path) = fasta("same", ">a x\nACGTACGTAC\nGTA\n>b\nNNNN\n");
        let ours = Reader::from_path(&path).unwrap();
        let theirs = rust_htslib::faidx::Reader::from_path(&path).unwrap();
        assert_eq!(ours.seq_names().unwrap(), theirs.seq_names().unwrap());
        for name in ["a", "b", "nope"] {
            assert_eq!(ours.fetch_seq_len(name), theirs.fetch_seq_len(name));
        }
        assert_eq!(ours.fetch_seq("a", 0, 12).unwrap(), b"ACGTACGTACGTA");
        assert_eq!(ours.fetch_seq("a", 8, 100).unwrap(), b"ACGTA");
        assert_eq!(ours.fetch_seq("a", 2, 9).unwrap(), b"GTACGTAC");
        assert_eq!(ours.fetch_seq("b", 1, 2).unwrap(), b"NN");
        // rust-htslib's fetch_seq is not called: its Vec would free htslib's buffer with the Rust allocator
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn failures_are_errors() {
        let (dir, path) = fasta("fail", ">a\nACGT\n");
        let ours = Reader::from_path(&path).unwrap();
        assert!(ours.fetch_seq("nope", 0, 2).is_err());
        assert!(Reader::from_path(dir.join("missing.fa")).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
