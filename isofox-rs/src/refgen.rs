//! Reference bases, as RefGenomeSource.getBaseString gives them (1-based, inclusive), with a window cached
//! for the GC ratios of the reads of a gene collection.

use anyhow::Result;
use hmf_common::refgenome::RefGenome;

pub struct RefSource {
    fasta: Option<RefGenome>,
    /// tests: sequences by chromosome, 0-based as MockRefGenome
    mock: std::collections::HashMap<String, Vec<u8>>,
    window_chr: String,
    window_start: i32,
    window: Vec<u8>,
    /// a window to load on first use: chromosome, start, end
    pending: Option<(String, i32, i32)>,
}

impl RefSource {
    pub fn open(path: &str) -> Result<RefSource> {
        Ok(RefSource { fasta: Some(RefGenome::open(path)?), mock: Default::default(), window_chr: String::new(), window_start: 0, window: Vec::new(), pending: None })
    }

    /// an in-memory reference (MockRefGenome: positions index the sequence from 0)
    #[cfg(test)]
    pub fn mock(seqs: Vec<(&str, String)>) -> RefSource {
        RefSource { fasta: None, mock: seqs.into_iter().map(|(c, s)| (c.to_string(), s.into_bytes())).collect(), window_chr: String::new(), window_start: 0, window: Vec::new(), pending: None }
    }

    fn fetch(&self, chr: &str, start: i32, end: i32) -> Vec<u8> {
        match &self.fasta {
            Some(f) => f.bases(chr, start, end).unwrap_or_default(),
            None => match self.mock.get(chr) {
                Some(s) if start >= 0 && (end as usize) < s.len() => s[start as usize..=end as usize].to_vec(),
                _ => Vec::new(),
            },
        }
    }

    pub fn chromosome_length(&self, chr: &str) -> i32 {
        match &self.fasta {
            Some(f) => f.chromosome_length(chr),
            None => self.mock.get(chr).map(|s| s.len() as i32).unwrap_or(0),
        }
    }

    /// loads the pending window if the region falls in it
    fn check_pending(&mut self, chr: &str, start: i32, end: i32) {
        if let Some((c, s, e)) = &self.pending {
            if c == chr && start >= *s && end <= *e {
                let (c, s, e) = self.pending.take().unwrap();
                self.set_window(&c, s, e);
            }
        }
    }

    /// the bases, or empty where the region is not in the reference (Java throws there)
    pub fn bases(&mut self, chr: &str, start: i32, end: i32) -> Vec<u8> {
        if end < start {
            return Vec::new();
        }
        self.check_pending(chr, start, end);
        if chr == self.window_chr && start >= self.window_start && end < self.window_start + self.window.len() as i32 {
            let o = (start - self.window_start) as usize;
            return self.window[o..o + (end - start + 1) as usize].to_vec();
        }
        self.fetch(chr, start, end)
    }

    pub fn base_string(&mut self, chr: &str, start: i32, end: i32) -> String {
        String::from_utf8(self.bases(chr, start, end)).unwrap_or_default()
    }

    /// a window to cache when first used (a gene collection's exons and reads)
    pub fn set_window_lazy(&mut self, chr: &str, start: i32, end: i32) {
        self.pending = Some((chr.to_string(), start.max(1), end));
    }

    /// caches a window of the reference
    pub fn set_window(&mut self, chr: &str, start: i32, end: i32) {
        let start = start.max(1);
        let len = self.chromosome_length(chr);
        let end = end.min(len);
        if end < start {
            self.window.clear();
            return;
        }
        self.window = self.fetch(chr, start, end);
        self.window_chr = chr.to_string();
        self.window_start = start;
    }

    /// GcRatioCounts.calcGcRatioFromReadRegions
    pub fn gc_ratio_of_regions(&mut self, chr: &str, regions: &[hmf_common::regions::BaseRegion]) -> f64 {
        let mut gc_total = 0.0;
        let mut bases_total = 0;
        for r in regions {
            let (count, len) = self.gc_count(chr, r.start, r.end);
            bases_total += len;
            gc_total += (count as f64 / len as f64) * len as f64;
        }
        gc_total / bases_total as f64
    }

    fn gc_count(&mut self, chr: &str, start: i32, end: i32) -> (i32, i32) {
        self.check_pending(chr, start, end);
        if chr == self.window_chr && start >= self.window_start && end < self.window_start + self.window.len() as i32 && end >= start {
            let o = (start - self.window_start) as usize;
            let s = &self.window[o..o + (end - start + 1) as usize];
            return (crate::adjusts::calc_gc_count(s), s.len() as i32);
        }
        let b = self.bases(chr, start, end);
        (crate::adjusts::calc_gc_count(&b), b.len() as i32)
    }
}
