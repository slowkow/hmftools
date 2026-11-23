//! Indexed FASTA access. Positions are 1-based and inclusive, as in the Java RefGenomeInterface.

use anyhow::Result;
use crate::common::faidx;

/// An indexed FASTA. An uncompressed file is read with positioned reads at offsets worked out from the
/// index; htslib reads it a byte at a time, which showed in profiles. Anything else goes through htslib.
pub struct FastaReader {
    reader: faidx::Reader,
    file: Option<std::fs::File>,
    /// per contig: length, offset of the first base, bases per line, bytes per line
    index: std::collections::HashMap<String, (u64, u64, u64, u64)>,
}

impl FastaReader {
    pub fn open(path: &str) -> Result<Self> {
        let reader = faidx::Reader::from_path(path)?;
        let mut index = std::collections::HashMap::new();
        let mut file = None;
        let compressed = path.ends_with(".gz") || path.ends_with(".bgz") || std::path::Path::new(&format!("{path}.gzi")).exists();
        if !compressed {
            if let Ok(text) = std::fs::read_to_string(format!("{path}.fai")) {
                for line in text.lines() {
                    let f: Vec<&str> = line.split('\t').collect();
                    if f.len() < 5 {
                        continue;
                    }
                    if let (Ok(len), Ok(off), Ok(lb), Ok(lw)) = (f[1].parse::<u64>(), f[2].parse::<u64>(), f[3].parse::<u64>(), f[4].parse::<u64>()) {
                        if lb > 0 && lw >= lb {
                            index.insert(f[0].to_string(), (len, off, lb, lw));
                        }
                    }
                }
                file = std::fs::File::open(path).ok();
            }
        }
        Ok(FastaReader { reader, file, index })
    }
    pub fn seq_names(&self) -> Result<Vec<String>> {
        Ok(self.reader.seq_names()?)
    }
    pub fn seq_len(&self, name: &str) -> u64 {
        self.reader.fetch_seq_len(name)
    }
    /// Bases from `begin` to `end`, zero-based and inclusive, cut at the contig end (as faidx_fetch_seq)
    pub fn fetch(&self, chr: &str, begin: usize, end: usize) -> Result<Vec<u8>> {
        if let (Some(file), Some(&(len, off, lb, lw))) = (&self.file, self.index.get(chr)) {
            use std::os::unix::fs::FileExt;
            let begin = begin as u64;
            let end = (end as u64).min(len.saturating_sub(1));
            if len == 0 || begin > end {
                return Ok(Vec::new());
            }
            let from = off + begin / lb * lw + begin % lb;
            let to = off + end / lb * lw + end % lb + 1;
            let mut buf = vec![0u8; (to - from) as usize];
            file.read_exact_at(&mut buf, from)?;
            buf.retain(|b| b.is_ascii_graphic());
            return Ok(buf);
        }
        Ok(self.reader.fetch_seq(chr, begin, end)?)
    }
}

pub struct RefGenome {
    reader: FastaReader,
    names: Vec<String>,
    lengths: Vec<i32>,
}

impl RefGenome {
    pub fn open(path: &str) -> Result<Self> {
        let reader = FastaReader::open(path)?;
        let names = reader.seq_names()?;
        let lengths = names.iter().map(|n| reader.seq_len(n) as i32).collect();
        Ok(RefGenome { reader, names, lengths })
    }
    /// contigs in dictionary order
    pub fn sequence_names(&self) -> &[String] {
        &self.names
    }
    pub fn chromosome_length(&self, chr: &str) -> i32 {
        self.names.iter().position(|n| n == chr).map(|i| self.lengths[i]).unwrap_or(0)
    }
    pub fn bases(&self, chr: &str, start: i32, end: i32) -> Result<Vec<u8>> {
        self.reader.fetch(chr, (start - 1) as usize, (end - 1) as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The positioned reads give the same bases as htslib over line ends, contig ends and a short last line
    #[test]
    fn positioned_fetch_matches_htslib() {
        let dir = std::env::temp_dir().join(format!("redux_rs_fasta_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.fa");
        let mut text = String::new();
        let mut lens = Vec::new();
        for (c, len) in [("a", 1usize), ("b", 60), ("c", 61), ("d", 187), ("e", 7)] {
            text.push_str(&format!(">{c} some description\n"));
            let seq: Vec<u8> = (0..len).map(|i| b"ACGTNacgtn"[(i * 7 + len) % 10]).collect();
            for line in seq.chunks(60) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
            lens.push((c, len));
        }
        std::fs::write(&path, text).unwrap();
        let path = path.to_str().unwrap().to_string();
        let plain = faidx::Reader::from_path(&path).unwrap(); // builds the index
        let fast = FastaReader::open(&path).unwrap();
        assert!(fast.file.is_some() && fast.index.len() == 5);
        for (c, len) in lens {
            for begin in 0..len + 2 {
                for end in begin..len + 3 {
                    assert_eq!(fast.fetch(c, begin, end).unwrap(), plain.fetch_seq(c, begin, end).unwrap(), "{c} {begin} {end}");
                }
            }
        }
        std::fs::remove_dir_all(&dir).ok();
    }
}
