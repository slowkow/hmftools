//! The `common` package: enums, transcript-exon references, the merging of regions, fragment type counts
//! and the per-gene-collection base depth.

use hmf_common::ensembl::TranscriptData;
use hmf_common::regions::BaseRegion;

pub const SE_START: usize = 0;
pub const SE_END: usize = 1;
pub const FS_UP: usize = 0;
pub const FS_DOWN: usize = 1;
pub const ORIENT_FWD: i8 = 1;
pub const ORIENT_REV: i8 = -1;
pub const NO_GENE_ID: i32 = -1;
pub const NO_EXON: i32 = -1;

pub const TRANS_COUNT: usize = 0;
pub const UNIQUE_TRANS_COUNT: usize = 1;

#[inline]
pub fn switch_index(i: usize) -> usize {
    1 - i
}

#[inline]
pub fn position_within(pos: i32, start: i32, end: i32) -> bool {
    pos >= start && pos <= end
}
#[inline]
pub fn positions_within(inner_start: i32, inner_end: i32, outer_start: i32, outer_end: i32) -> bool {
    inner_start <= inner_end && inner_start >= outer_start && inner_end <= outer_end
}
#[inline]
pub fn positions_overlap(start1: i32, end1: i32, start2: i32, end2: i32) -> bool {
    !(start1 > end2 || end1 < start2)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RegionMatchType {
    None,
    ExonBoundary,
    WithinExon,
    ExonIntron,
    Intron,
}

impl RegionMatchType {
    #[allow(dead_code)]
    pub fn name(self) -> &'static str {
        match self {
            RegionMatchType::None => "NONE",
            RegionMatchType::ExonBoundary => "EXON_BOUNDARY",
            RegionMatchType::WithinExon => "WITHIN_EXON",
            RegionMatchType::ExonIntron => "EXON_INTRON",
            RegionMatchType::Intron => "INTRON",
        }
    }
    pub fn valid_exon_match(self) -> bool {
        self == RegionMatchType::ExonBoundary || self == RegionMatchType::WithinExon
    }
    pub fn exon_boundary(self) -> bool {
        self == RegionMatchType::ExonBoundary
    }
    pub fn match_rank(self) -> i32 {
        match self {
            RegionMatchType::ExonBoundary => 4,
            RegionMatchType::WithinExon => 3,
            RegionMatchType::ExonIntron => 2,
            RegionMatchType::Intron => 1,
            RegionMatchType::None => 0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransMatchType {
    Unknown,
    SpliceJunction,
    Exonic,
    Unspliced,
    OtherTrans,
    Alt,
}

impl TransMatchType {
    #[allow(dead_code)]
    pub fn name(self) -> &'static str {
        match self {
            TransMatchType::Unknown => "UNKNOWN",
            TransMatchType::SpliceJunction => "SPLICE_JUNCTION",
            TransMatchType::Exonic => "EXONIC",
            TransMatchType::Unspliced => "UNSPLICED",
            TransMatchType::OtherTrans => "OTHER_TRANS",
            TransMatchType::Alt => "ALT",
        }
    }
}

/// ReadTranscriptUtils.validTranscriptType
pub fn valid_transcript_type(t: TransMatchType) -> bool {
    t == TransMatchType::Exonic || t == TransMatchType::SpliceJunction
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FragmentType {
    Total = 0,
    TransSupporting,
    Alt,
    Unspliced,
    Chimeric,
    Duplicate,
    MultiMapped,
    ForwardStrand,
    ReverseStrand,
}

impl FragmentType {
    #[allow(dead_code)]
    pub fn name(self) -> &'static str {
        match self {
            FragmentType::Total => "TOTAL",
            FragmentType::TransSupporting => "TRANS_SUPPORTING",
            FragmentType::Alt => "ALT",
            FragmentType::Unspliced => "UNSPLICED",
            FragmentType::Chimeric => "CHIMERIC",
            FragmentType::Duplicate => "DUPLICATE",
            FragmentType::MultiMapped => "MULTI_MAPPED",
            FragmentType::ForwardStrand => "FORWARD_STRAND",
            FragmentType::ReverseStrand => "REVERSE_STRAND",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub enum FragmentMatchType {
    Spliced = 0,
    Short,
    Long,
    Unspliced,
    Discordant,
}

/// FragmentTypeCounts
#[derive(Clone, Debug, Default)]
pub struct FragmentTypeCounts {
    pub counts: [i64; 9],
}

impl FragmentTypeCounts {
    pub fn type_count(&self, t: FragmentType) -> i64 {
        self.counts[t as usize]
    }
    pub fn add_count(&mut self, t: FragmentType, count: i64) {
        self.counts[t as usize] += count;
    }
    pub fn combine(&mut self, other: &FragmentTypeCounts) {
        for i in 0..9 {
            self.counts[i] += other.counts[i];
        }
    }
}

/// TransExonRef: a transcript's exon
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransExonRef {
    pub gene_id: &'static str,
    pub trans_id: i32,
    pub trans_name: &'static str,
    pub exon_rank: i32,
    pub is_canonical: Option<bool>,
}

impl TransExonRef {
    pub fn new(trans: &'static TranscriptData, exon_rank: i32) -> Self {
        TransExonRef { gene_id: &trans.gene_id, trans_id: trans.trans_id, trans_name: &trans.trans_name, exon_rank, is_canonical: Some(trans.is_canonical) }
    }
    pub fn matches(&self, other: &TransExonRef) -> bool {
        other.trans_id == self.trans_id && other.exon_rank == self.exon_rank
    }
    fn matches_within(&self, other: &TransExonRef, max_diff_vs_other: i32) -> bool {
        if self.trans_id != other.trans_id {
            return false;
        }
        if max_diff_vs_other >= 0 {
            other.exon_rank >= self.exon_rank && other.exon_rank <= self.exon_rank + max_diff_vs_other
        } else {
            other.exon_rank >= self.exon_rank + max_diff_vs_other && other.exon_rank <= self.exon_rank
        }
    }
    pub fn is_canonical(&self) -> bool {
        self.is_canonical == Some(true)
    }
    pub fn has_match(list1: &[TransExonRef], list2: &[TransExonRef]) -> bool {
        list1.iter().any(|x| list2.iter().any(|y| x.matches(y)))
    }
    pub fn has_match_within_range(list1: &[TransExonRef], list2: &[TransExonRef], max_list2_diff: i32) -> bool {
        list1.iter().any(|x| list2.iter().any(|y| x.matches_within(y, max_list2_diff)))
    }
}

/// CommonUtils.deriveCommonRegions: merges overlapping (or adjacent) regions of two sorted lists
pub fn derive_common_regions(regions1: &[BaseRegion], regions2: &[BaseRegion]) -> Vec<BaseRegion> {
    if regions1.is_empty() || regions2.is_empty() {
        return if regions1.is_empty() { regions2.to_vec() } else { regions1.to_vec() };
    }
    let mut new_regions: Vec<BaseRegion> = Vec::with_capacity(regions1.len() + regions2.len());
    if regions1[regions1.len() - 1].end < regions2[0].start {
        new_regions.extend_from_slice(regions1);
        new_regions.extend_from_slice(regions2);
        return new_regions;
    } else if regions2[regions2.len() - 1].end < regions1[0].start {
        new_regions.extend_from_slice(regions2);
        new_regions.extend_from_slice(regions1);
        return new_regions;
    }
    let mut i1 = 0;
    let mut i2 = 0;
    while i1 < regions1.len() || i2 < regions2.len() {
        let r1 = regions1.get(i1).copied();
        let r2 = regions2.get(i2).copied();
        match (r1, r2) {
            (Some(a), Some(b)) => {
                if a.end < b.start - 1 {
                    new_regions.push(a);
                    i1 += 1;
                    continue;
                } else if b.end < a.start - 1 {
                    new_regions.push(b);
                    i2 += 1;
                    continue;
                }
                let mut nr = BaseRegion::new(a.start.min(b.start), a.end.max(b.end));
                i1 += 1;
                i2 += 1;
                loop {
                    if i1 >= regions1.len() && i2 >= regions2.len() {
                        break;
                    }
                    let mut merged = false;
                    if let Some(a) = regions1.get(i1) {
                        if a.start <= nr.end + 1 {
                            nr.end = a.end.max(nr.end);
                            i1 += 1;
                            merged = true;
                        }
                    }
                    if let Some(b) = regions2.get(i2) {
                        if b.start <= nr.end + 1 {
                            nr.end = b.end.max(nr.end);
                            i2 += 1;
                            merged = true;
                        }
                    }
                    if !merged {
                        break;
                    }
                }
                new_regions.push(nr);
            }
            (Some(a), None) => {
                new_regions.push(a);
                i1 += 1;
            }
            (None, Some(b)) => {
                new_regions.push(b);
                i2 += 1;
            }
            (None, None) => break,
        }
    }
    new_regions
}

/// impliedSvType
pub fn implied_sv_type(chromosomes: &[&str; 2], orientations: [i8; 2]) -> &'static str {
    if chromosomes[0] == chromosomes[1] {
        if orientations[0] == orientations[1] {
            "INV"
        } else if orientations[0] == ORIENT_FWD {
            "DEL"
        } else {
            "DUP"
        }
    } else {
        "BND"
    }
}

pub fn canonical_acceptor(context: &str, strand: i8) -> bool {
    if strand == 1 { context == "AG" } else { context == "CT" }
}

pub fn canonical_donor(context: &str, strand: i8) -> bool {
    if strand == 1 { context == "GT" || context == "GC" } else { context == "AC" || context == "GC" }
}

/// HumanChromosome.lowerChromosome
pub fn lower_chromosome(chr: &str, other: &str) -> bool {
    hmf_common::regions::chromosome_rank_java(chr) < hmf_common::regions::chromosome_rank_java(other)
}

#[allow(dead_code)]
const BD_MIN_DEPTH_COUNT: i32 = 2;
const BD_DEFAULT_SIZE: usize = 3_000_000;

/// BaseDepth: read depth over a gene collection's range
pub struct BaseDepth {
    pub range: [i32; 2],
    depth: Vec<i32>,
    /// the entries the last range used, the only ones to clear
    used: usize,
}

impl Default for BaseDepth {
    fn default() -> Self {
        Self::new()
    }
}

impl BaseDepth {
    pub fn new() -> Self {
        BaseDepth { range: [0, 0], depth: vec![0; BD_DEFAULT_SIZE], used: BD_DEFAULT_SIZE }
    }
    pub fn length(&self) -> i32 {
        self.range[1] - self.range[0] + 1
    }
    pub fn initialise(&mut self, range: [i32; 2]) {
        self.range = range;
        let len = self.length().max(0) as usize;
        if len > self.depth.len() {
            hmf_common::info!("growing base depth: old({}) new({}) range({}-{})", self.depth.len(), len, range[0], range[1]);
            self.depth = vec![0; len];
        } else {
            // Java clears the whole array; only the entries of the last range can be set
            let n = self.used.min(self.depth.len());
            self.depth[..n].iter_mut().for_each(|d| *d = 0);
        }
        self.used = len.min(self.depth.len());
    }
    pub fn process_read(&mut self, alignments: &[BaseRegion]) {
        for a in alignments {
            let (s, e) = (a.start, a.end);
            if s > self.range[1] || e < self.range[0] {
                continue;
            }
            let idx = if s > self.range[0] { s - self.range[0] } else { 0 };
            let overlap = e.min(self.range[1]) - s.max(self.range[0]) + 1;
            if idx + overlap > self.length() {
                hmf_common::error!(
                    "baseDepth(region({}-{})) read coords({} -> {}) regionBaseIndex({}) overlap({}) regionLength({})",
                    self.range[0],
                    self.range[1],
                    s,
                    e,
                    idx,
                    overlap,
                    self.length()
                );
                return;
            }
            for d in &mut self.depth[idx as usize..(idx + overlap) as usize] {
                *d += 1;
            }
        }
    }
    pub fn has_position(&self, pos: i32) -> bool {
        position_within(pos, self.range[0], self.range[1])
    }
    pub fn depth_at_base(&self, pos: i32) -> i32 {
        if !self.has_position(pos) {
            return 0;
        }
        self.depth[(pos - self.range[0]) as usize]
    }
    #[allow(dead_code)]
    pub fn bases_with_depth(&self) -> i32 {
        self.depth[..self.used].iter().filter(|d| **d >= BD_MIN_DEPTH_COUNT).count() as i32
    }
}

/// GeneRegionFilters' gene id list loader is in config; this is TranscriptData.exons as a slice
#[inline]
#[allow(dead_code)]
pub fn exons(t: &TranscriptData) -> &[hmf_common::ensembl::ExonData] {
    &t.exons
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(s: i32, e: i32) -> BaseRegion {
        BaseRegion::new(s, e)
    }

    // CommonUtilsTest.testCommonRegions
    #[test]
    fn common_regions() {
        let l1 = vec![r(10, 20), r(30, 40), r(50, 60)];
        let l2 = vec![r(5, 15), r(25, 35), r(55, 65)];
        let c = derive_common_regions(&l1, &l2);
        assert_eq!(c, vec![r(5, 20), r(25, 40), r(50, 65)]);
        // no overlaps
        let c = derive_common_regions(&[r(10, 20)], &[r(30, 40)]);
        assert_eq!(c, vec![r(10, 20), r(30, 40)]);
        let c = derive_common_regions(&[r(30, 40)], &[r(10, 20)]);
        assert_eq!(c, vec![r(10, 20), r(30, 40)]);
        // one spanning all
        let c = derive_common_regions(&[r(10, 100)], &[r(20, 30), r(40, 50), r(90, 110)]);
        assert_eq!(c, vec![r(10, 110)]);
        // regions wholly before the other list's are kept apart, even if adjacent (the early exit)
        let c = derive_common_regions(&[r(10, 20)], &[r(21, 30)]);
        assert_eq!(c, vec![r(10, 20), r(21, 30)]);
        // adjacent regions merge once the lists interleave
        let c = derive_common_regions(&[r(10, 20), r(40, 50)], &[r(21, 30)]);
        assert_eq!(c, vec![r(10, 30), r(40, 50)]);
        let c = derive_common_regions(&[], &[r(1, 2)]);
        assert_eq!(c, vec![r(1, 2)]);
    }

    #[test]
    fn base_depth() {
        let mut bd = BaseDepth::new();
        bd.initialise([100, 200]);
        bd.process_read(&[r(90, 110), r(150, 160)]);
        bd.process_read(&[r(105, 120)]);
        assert_eq!(bd.depth_at_base(100), 1);
        assert_eq!(bd.depth_at_base(106), 2);
        assert_eq!(bd.depth_at_base(99), 0);
        assert_eq!(bd.depth_at_base(155), 1);
        assert_eq!(bd.bases_with_depth(), 6);
    }

    #[test]
    fn sv_types() {
        assert_eq!(implied_sv_type(&["1", "1"], [1, -1]), "DEL");
        assert_eq!(implied_sv_type(&["1", "1"], [-1, 1]), "DUP");
        assert_eq!(implied_sv_type(&["1", "1"], [1, 1]), "INV");
        assert_eq!(implied_sv_type(&["1", "2"], [1, -1]), "BND");
    }
}
