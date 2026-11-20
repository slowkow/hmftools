//! Indexed FASTA access. Positions are 1-based and inclusive, as in the Java RefGenomeInterface.

use anyhow::Result;
use rust_htslib::faidx;

pub struct RefGenome {
    reader: faidx::Reader,
    names: Vec<String>,
    lengths: Vec<i32>,
}

impl RefGenome {
    pub fn open(path: &str) -> Result<Self> {
        let reader = faidx::Reader::from_path(path)?;
        let names = reader.seq_names()?;
        let lengths = names.iter().map(|n| reader.fetch_seq_len(n) as i32).collect();
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
        Ok(self.reader.fetch_seq(chr, (start - 1) as usize, (end - 1) as usize)?)
    }
}
