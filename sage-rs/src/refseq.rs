//! RefSequence: reference bases for a region plus a 1000 base buffer either side.

use anyhow::Result;
use crate::faidx;

pub struct RefGenome {
    reader: faidx::Reader,
}

impl RefGenome {
    pub fn open(path: &str) -> Result<Self> {
        Ok(RefGenome { reader: faidx::Reader::from_path(path)? })
    }
    pub fn sequence_names(&self) -> Result<Vec<String>> {
        Ok(self.reader.seq_names()?)
    }
    pub fn chromosome_length(&self, chr: &str) -> i32 {
        self.reader.fetch_seq_len(chr) as i32
    }
    /// 1-based inclusive
    pub fn bases(&self, chr: &str, start: i32, end: i32) -> Result<Vec<u8>> {
        Ok(self.reader.fetch_seq(chr, (start - 1) as usize, (end - 1) as usize)?)
    }
}

pub struct RefSequence {
    pub start: i32,
    pub end: i32,
    pub bases: Vec<u8>,
}

const BUFFER: i32 = 1000;

impl RefSequence {
    pub fn new(genome: &RefGenome, chr: &str, region_start: i32, region_end: i32) -> Result<Self> {
        let len = genome.chromosome_length(chr);
        let start = 1.max(region_start - BUFFER);
        let end = len.min(region_end + BUFFER);
        Ok(RefSequence { start, end, bases: genome.bases(chr, start, end)? })
    }
    pub fn from_bases(start: i32, bases: Vec<u8>) -> Self {
        RefSequence { start, end: start + bases.len() as i32 - 1, bases }
    }
    pub fn index(&self, pos: i32) -> i32 {
        pos - self.start
    }
    pub fn base(&self, pos: i32) -> u8 {
        self.bases[(pos - self.start) as usize]
    }
    /// baseRange: inclusive positions, None when out of range
    pub fn base_range(&self, pos_start: i32, pos_end: i32) -> Option<&[u8]> {
        let s = self.index(pos_start);
        let e = self.index(pos_end);
        if s < 0 || e < s || e >= self.bases.len() as i32 {
            return None;
        }
        Some(&self.bases[s as usize..=e as usize])
    }
    /// positionBases: note the Java bound is `end >= length` -> null
    pub fn position_bases(&self, pos_start: i32, pos_end: i32) -> Option<&[u8]> {
        self.base_range(pos_start, pos_end)
    }
}
