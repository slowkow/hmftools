//! LohIdentifier (com.hartwig.hmftools.amber.LohIdentifier) and the logic of BafLohData: tumor-only LOH
//! segments, whose sites are kept ("recovered") although they fail the tumor-only support filter.
//!
//! Per autosome: the valid sites (depth at least the tumor minimum depth, less than 7.5% of reads filtered
//! for base quality, mapping quality or sequencing technology, outside HLA-A/B/C) each get a BafLohData,
//! which says whether the site's major allele fraction is below 0.96, 0.90 and 0.80 (with at least two reads
//! of each allele). A window of 101 valid sites centred on each site gives the rates of such sites; a site
//! is an LOH candidate where the rates are below 0.15, 0.04 and 0. Runs of candidates spanning at least
//! 20Mb, less any centromere gaps they overlap, are LOH segments: every fifth site of a segment, widened by
//! 50 sites on each side, is marked recovered.
//!
//! Findings (REPORT.md): AM51 (the first full window is never used, so the first valid site of a chromosome
//! can never be recovered), AM52 (an overlapping centromere gap is subtracted whole, not only its overlap).

#![allow(dead_code)]

use crate::purity::genome;
use crate::types::{BafLohData, PositionEvidence, TumorBaf};
use hmf_common::regions::ChrBaseRegion;

pub const LOH_QUAL_FILTERED_THRESHOLD: f64 = 0.075;
pub const LOH_REF_ALT_MIN_COUNT: i32 = 2;
pub const LOH_MOVING_AVG_COUNT: usize = 50;
pub const LOH_MOVING_ITEM_COUNT: usize = LOH_MOVING_AVG_COUNT * 2 + 1;
pub const LOH_REF_ALT_THRESHOLD_HIGH: f64 = 0.96;
pub const LOH_REF_ALT_THRESHOLD_MED: f64 = 0.90;
pub const LOH_REF_ALT_THRESHOLD_LOW: f64 = 0.80;
pub const LOH_MIN_LENGTH: i32 = 20_000_000;
/// every fifth point of a segment is recovered
pub const LOH_POINT_DENSITY: usize = 5;
pub const VALID_BAF_RATE_HIGH: f64 = 0.15;
pub const VALID_BAF_RATE_MED: f64 = 0.04;
pub const VALID_BAF_RATE_LOW: f64 = 0.0;

impl BafLohData {
    /// BafLohData.setBafStatus
    pub fn set_baf_status(&mut self, ref_support: i32, alt_support: i32) {
        if ref_support.min(alt_support) >= LOH_REF_ALT_MIN_COUNT {
            // Java: double totalCount = refSupport + altSupport (an int sum)
            let total_count = ref_support.wrapping_add(alt_support) as f64;
            let baf = ref_support.max(alt_support) as f64 / total_count;
            if baf < LOH_REF_ALT_THRESHOLD_HIGH {
                self.valid_baf_high = true;
                if baf < LOH_REF_ALT_THRESHOLD_MED {
                    self.valid_baf_medium = true;
                    if baf < LOH_REF_ALT_THRESHOLD_LOW {
                        self.valid_baf_low = true;
                    }
                }
            }
        }
    }

    /// BafLohData.setValidBafData: the rates of valid sites in the window around this one
    pub fn set_valid_baf_data(&mut self, valid_count_high: i32, valid_count_med: i32, valid_count_low: i32) {
        self.valid_rate_high = valid_count_high as f64 / LOH_MOVING_ITEM_COUNT as f64;
        self.valid_rate_medium = valid_count_med as f64 / LOH_MOVING_ITEM_COUNT as f64;
        self.valid_rate_low = valid_count_low as f64 / LOH_MOVING_ITEM_COUNT as f64;
        if self.valid_rate_high < VALID_BAF_RATE_HIGH && self.valid_rate_medium < VALID_BAF_RATE_MED && self.valid_rate_low <= VALID_BAF_RATE_LOW {
            self.loh_candidate = true;
        }
    }
}

impl PositionEvidence {
    /// PositionEvidence.setLohData: a new BafLohData with the BAF status of this site's counts
    pub fn set_loh_data(&mut self) {
        let mut data = BafLohData::default();
        data.set_baf_status(self.ref_support, self.alt_support);
        self.loh_data = Some(data);
    }
}

/// AmberUtils.aboveQualFilter: the share of reads filtered for base quality, mapping quality or sequencing
/// technology is below `threshold` (false with no depth)
pub fn above_qual_filter(evidence: &PositionEvidence, threshold: f64) -> bool {
    if evidence.read_depth == 0 {
        return false;
    }
    let filtered = evidence.base_qual_filtered.wrapping_add(evidence.map_qual_filtered).wrapping_add(evidence.seq_tech_filtered);
    let filtered_percent = filtered as f64 / evidence.read_depth as f64;
    filtered_percent < threshold
}

pub struct LohIdentifier {
    v38: bool,
    tumor_min_depth: i32,
    hla_regions: Vec<ChrBaseRegion>,
}

impl LohIdentifier {
    pub fn new(v38: bool, tumor_min_depth: i32) -> LohIdentifier {
        LohIdentifier { v38, tumor_min_depth, hla_regions: genome::hla_regions(v38) }
    }

    /// LohIdentifier.isValidSite
    fn is_valid_site(&self, evidence: &PositionEvidence) -> bool {
        if evidence.read_depth < self.tumor_min_depth {
            return false;
        }
        if !above_qual_filter(evidence, LOH_QUAL_FILTERED_THRESHOLD) {
            return false;
        }
        !self.hla_regions.iter().any(|r| r.chromosome == evidence.chromosome && r.contains_position(evidence.position))
    }

    /// findLohSegments for one chromosome: sets `tumor.loh_data` of the sites in `bafs`, which are in the
    /// chromosome's site order. Java skips X and Y; so does this (the caller may pass them).
    pub fn find_chromosome_loh_segments(&self, chromosome: &str, bafs: &mut [TumorBaf]) {
        // Both X and Y are excluded (Java: chromosome.isAllosome())
        if genome::is_allosome(chromosome) {
            return;
        }
        let valid: Vec<usize> = (0..bafs.len()).filter(|&i| self.is_valid_site(&bafs[i].tumor)).collect();
        for &i in &valid {
            bafs[i].tumor.set_loh_data();
        }
        let flags = |b: &TumorBaf| {
            let d = b.tumor.loh_data.as_ref().unwrap();
            (d.valid_baf_high as i32, d.valid_baf_medium as i32, d.valid_baf_low as i32)
        };

        // a moving count over windows of LOH_MOVING_ITEM_COUNT valid sites, set on the window's centre
        let java_window = crate::compat::java("AM51");
        let (mut high, mut med, mut low) = (0i32, 0i32, 0i32);
        for i in 0..valid.len() {
            let (h, m, l) = flags(&bafs[valid[i]]);
            high += h;
            med += m;
            low += l;
            if i >= LOH_MOVING_ITEM_COUNT {
                let (h, m, l) = flags(&bafs[valid[i - LOH_MOVING_ITEM_COUNT]]);
                high -= h;
                med -= m;
                low -= l;
                if high < 0 || med < 0 || low < 0 {
                    crate::error!("chromosome({}) invalid valid BAF counts", java_chromosome_name(chromosome));
                }
            }
            // Java sets the centre only once a site has been dropped from the window (i > 100), so the
            // first full window, sites 0 to 100, never reaches site 50
            let full = if java_window { i >= LOH_MOVING_ITEM_COUNT } else { i + 1 >= LOH_MOVING_ITEM_COUNT };
            if full {
                let target = valid[i - LOH_MOVING_AVG_COUNT];
                bafs[target].tumor.loh_data.as_mut().unwrap().set_valid_baf_data(high, med, low);
            }
        }

        // find continuous regions of LOH
        let gaps = genome::centromere_gaps(self.v38, chromosome);
        let java_gaps = crate::compat::java("AM52");
        let mut loh_start: Option<usize> = None;
        let mut loh_end = 0usize;
        for i in 0..valid.len() {
            let candidate = bafs[valid[i]].tumor.loh_data.as_ref().unwrap().loh_candidate;
            if candidate {
                if loh_start.is_some() {
                    loh_end = i;
                } else {
                    loh_start = Some(i);
                    loh_end = i;
                }
            } else if let Some(start) = loh_start {
                // Java leaves a run still open after the last site unprocessed; the last 50 sites are never
                // candidates, so no run can be open there
                self.mark_segment(chromosome, bafs, &valid, start, loh_end, &gaps, java_gaps);
                loh_start = None;
            }
        }
    }

    fn mark_segment(&self, chromosome: &str, bafs: &mut [TumorBaf], valid: &[usize], start: usize, end: usize, gaps: &[(i32, i32)], java_gaps: bool) {
        let start_pos = bafs[valid[start]].tumor.position;
        let end_pos = bafs[valid[end]].tumor.position;
        let mut loh_length = end_pos.wrapping_sub(start_pos);
        // BaseRegion.positionsOverlap; Java subtracts the whole base length of each overlapping gap
        let mut gap_length = 0i32;
        for &(gap_start, gap_end) in gaps {
            if !(gap_start > end_pos || gap_end < start_pos) {
                gap_length = gap_length.wrapping_add(if java_gaps { gap_end - gap_start + 1 } else { gap_end.min(end_pos) - gap_start.max(start_pos) + 1 });
            }
        }
        loh_length = loh_length.wrapping_sub(gap_length);
        if loh_length < LOH_MIN_LENGTH {
            return;
        }
        crate::debug!("chromosome({}) LOH segment({} - {}) points({})", java_chromosome_name(chromosome), start_pos, end_pos, end - start + 1);

        // mark these points as LOH, plus the ones either side used in the moving average
        let widened_start = start.saturating_sub(LOH_MOVING_AVG_COUNT);
        let widened_end = (end + LOH_MOVING_AVG_COUNT).min(valid.len() - 1);
        let mut last_recovered: Option<usize> = None;
        for j in widened_start..=widened_end {
            let evidence = &mut bafs[valid[j]].tumor;
            if evidence.read_depth >= self.tumor_min_depth && last_recovered.is_none_or(|last| j - last >= LOH_POINT_DENSITY) {
                evidence.loh_data.as_mut().unwrap().loh_recovered = true;
                last_recovered = Some(j);
            }
        }
    }
}

/// HumanChromosome.toString, as Java logs the chromosome
fn java_chromosome_name(chromosome: &str) -> &str {
    if chromosome.len() >= 3 && chromosome[..3].eq_ignore_ascii_case("chr") { &chromosome[3..] } else { chromosome }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compat::with_java;
    use crate::types::AmberBase;

    fn status(ref_support: i32, alt_support: i32) -> (bool, bool, bool) {
        let mut d = BafLohData::default();
        d.set_baf_status(ref_support, alt_support);
        (d.valid_baf_high, d.valid_baf_medium, d.valid_baf_low)
    }

    #[test]
    fn baf_status() {
        assert_eq!(status(1, 50), (false, false, false)); // fewer than two reads of one allele
        assert_eq!(status(50, 1), (false, false, false));
        assert_eq!(status(2, 48), (false, false, false)); // 0.96 is not below 0.96
        assert_eq!(status(3, 47), (true, false, false));
        assert_eq!(status(5, 45), (true, false, false)); // 0.90
        assert_eq!(status(6, 44), (true, true, false));
        assert_eq!(status(40, 10), (true, true, false)); // 0.80
        assert_eq!(status(39, 11), (true, true, true));
        assert_eq!(status(25, 25), (true, true, true));
        assert_eq!(status(0, 0), (false, false, false));
    }

    #[test]
    fn valid_baf_data() {
        let rates = |h, m, l| {
            let mut d = BafLohData::default();
            d.set_valid_baf_data(h, m, l);
            d
        };
        assert!(rates(15, 4, 0).loh_candidate); // 15/101 < 0.15, 4/101 < 0.04
        assert!(!rates(16, 4, 0).loh_candidate);
        assert!(!rates(15, 5, 0).loh_candidate);
        assert!(!rates(15, 4, 1).loh_candidate);
        let d = rates(10, 2, 0);
        assert_eq!(d.valid_rate_high, 10.0 / 101.0);
        assert_eq!(d.valid_rate_medium, 2.0 / 101.0);
        assert_eq!(d.valid_rate_low, 0.0);
    }

    #[test]
    fn qual_filter() {
        let mut e = PositionEvidence::new("chr1", 100, AmberBase::A, AmberBase::C);
        assert!(!above_qual_filter(&e, 0.075));
        e.read_depth = 100;
        e.base_qual_filtered = 4;
        e.map_qual_filtered = 3;
        assert!(above_qual_filter(&e, 0.075));
        e.seq_tech_filtered = 1;
        assert!(!above_qual_filter(&e, 0.075));
        assert!(above_qual_filter(&e, 0.15));
    }

    /// A chromosome of sites `spacing` apart from `first`; sites in `loh` (index ranges) have one allele
    /// only (LOH), the others are heterozygous.
    fn chromosome(chr: &str, n: usize, first: i32, spacing: i32, loh: &[std::ops::Range<usize>]) -> Vec<TumorBaf> {
        (0..n)
            .map(|i| {
                let mut normal = PositionEvidence::new(chr, first + i as i32 * spacing, AmberBase::A, AmberBase::G);
                normal.read_depth = 40;
                let mut baf = TumorBaf::from_normal(&normal);
                let t = &mut baf.tumor;
                t.read_depth = 40;
                if loh.iter().any(|r| r.contains(&i)) {
                    t.ref_support = 40;
                } else {
                    t.ref_support = 20;
                    t.alt_support = 20;
                }
                baf
            })
            .collect()
    }

    fn recovered(bafs: &[TumorBaf]) -> Vec<usize> {
        (0..bafs.len()).filter(|&i| bafs[i].tumor.loh_recovered()).collect()
    }

    /// The expected result worked out directly: candidates are the sites whose full window has at most 15
    /// high, 4 medium and 0 low valid sites, among the centres Java sets (`first_centre` to n - 51).
    fn expected(bafs: &[TumorBaf], first_centre: usize, min_length: i32) -> Vec<usize> {
        let n = bafs.len();
        let valid_low: Vec<bool> = bafs.iter().map(|b| status(b.tumor.ref_support, b.tumor.alt_support).2).collect();
        let valid_high: Vec<bool> = bafs.iter().map(|b| status(b.tumor.ref_support, b.tumor.alt_support).0).collect();
        let mut candidate = vec![false; n];
        for c in first_centre..n.saturating_sub(50) {
            let h = (c - 50..=c + 50).filter(|&j| valid_high[j]).count();
            let l = (c - 50..=c + 50).filter(|&j| valid_low[j]).count();
            candidate[c] = h <= 15 && l == 0;
        }
        let mut out = Vec::new();
        let mut i = 0;
        while i < n {
            if !candidate[i] {
                i += 1;
                continue;
            }
            let s = i;
            while i < n && candidate[i] {
                i += 1;
            }
            let e = i - 1;
            if bafs[e].tumor.position - bafs[s].tumor.position >= min_length {
                let (ws, we) = (s.saturating_sub(50), (e + 50).min(n - 1));
                out.extend((ws..=we).step_by(5));
            }
        }
        out
    }

    #[test]
    fn segment_in_the_middle() {
        // 1000 sites 100kb apart: an LOH block of 400 sites (40Mb) from site 300; the candidates are the
        // 300 sites whose window lies inside it
        for java in [true, false] {
            with_java(java, || {
                let mut bafs = chromosome("chr3", 1000, 1_000_000, 100_000, &[300..700]);
                LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr3", &mut bafs);
                let got = recovered(&bafs);
                assert_eq!(got, expected(&bafs, if java { 51 } else { 50 }, LOH_MIN_LENGTH));
                assert!(!got.is_empty());
                assert!(bafs.iter().all(|b| b.tumor.loh_data.is_some()));
            });
        }
    }

    #[test]
    fn short_segment_is_not_recovered() {
        // 150 LOH sites 100kb apart: 14.9Mb, below 20Mb
        let mut bafs = chromosome("chr3", 1000, 1_000_000, 100_000, &[400..550]);
        LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr3", &mut bafs);
        assert!(recovered(&bafs).is_empty());
        assert!(bafs.iter().any(|b| b.tumor.loh_data.as_ref().unwrap().loh_candidate));
    }

    #[test]
    fn segment_at_chromosome_start() {
        // LOH from the first site: Java's first centre is site 51, so site 0 is never recovered (AM51)
        let run = |java| {
            with_java(java, || {
                let mut bafs = chromosome("chr5", 1000, 1_000_000, 100_000, &[0..400]);
                LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr5", &mut bafs);
                assert!(!bafs[50].tumor.loh_data.as_ref().unwrap().loh_candidate || !java);
                (recovered(&bafs), expected(&bafs, if java { 51 } else { 50 }, LOH_MIN_LENGTH))
            })
        };
        let (got, want) = run(true);
        assert_eq!(got, want);
        assert_eq!(got[0], 1);
        let (got, want) = run(false);
        assert_eq!(got, want);
        assert_eq!(got[0], 0);
    }

    #[test]
    fn segment_at_chromosome_end() {
        let mut bafs = chromosome("chr5", 1000, 1_000_000, 100_000, &[600..1000]);
        LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr5", &mut bafs);
        let got = recovered(&bafs);
        assert_eq!(got, expected(&bafs, 51, LOH_MIN_LENGTH));
        // the segment ends at site 948 (the last centre); widened by 50 to the last site
        assert_eq!((got[0], *got.last().unwrap()), (600, 995));
    }

    #[test]
    fn allosomes_are_skipped() {
        for chr in ["chrX", "chrY", "X"] {
            let mut bafs = chromosome(chr, 1000, 1_000_000, 100_000, &[400..700]);
            LohIdentifier::new(true, 25).find_chromosome_loh_segments(chr, &mut bafs);
            assert!(bafs.iter().all(|b| b.tumor.loh_data.is_none()));
        }
    }

    #[test]
    fn few_sites() {
        // 101 valid sites, all LOH: Java sets no centre (none is ever dropped from the window); the
        // corrected window sets site 50 only, a run of one site, too short
        for java in [true, false] {
            with_java(java, || {
                let mut bafs = chromosome("chr2", 101, 1_000_000, 100_000, &[0..101]);
                LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr2", &mut bafs);
                let candidates = bafs.iter().filter(|b| b.tumor.loh_data.as_ref().unwrap().loh_candidate).count();
                assert_eq!(candidates, if java { 0 } else { 1 });
                assert!(recovered(&bafs).is_empty());
            });
        }
        let mut bafs = chromosome("chr2", 10, 1_000_000, 100_000, &[]);
        LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr2", &mut bafs);
        assert!(recovered(&bafs).is_empty());
    }

    #[test]
    fn invalid_sites_are_left_out() {
        // shallow, quality-filtered and HLA sites get no LOH data and do not count in the windows
        let mut bafs = chromosome("chr6", 1000, 1_000_000, 50_000, &[200..900]);
        bafs[10].tumor.read_depth = 24;
        bafs[11].tumor.base_qual_filtered = 3; // 3/40 = 0.075, not below
        let hla = bafs.iter().position(|b| b.tumor.position >= 29942532 - 5000).unwrap();
        LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr6", &mut bafs);
        assert!(bafs[10].tumor.loh_data.is_none() && bafs[11].tumor.loh_data.is_none());
        assert!(bafs[hla].tumor.loh_data.is_none());
        assert!(bafs[12].tumor.loh_data.is_some());
        // the same without the invalid sites
        let mut kept: Vec<TumorBaf> = chromosome("chr6", 1000, 1_000_000, 50_000, &[200..900]);
        kept.remove(hla);
        kept.remove(11);
        kept.remove(10);
        LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr6", &mut kept);
        let a: Vec<i32> = bafs.iter().filter(|b| b.tumor.loh_recovered()).map(|b| b.tumor.position).collect();
        let b: Vec<i32> = kept.iter().filter(|b| b.tumor.loh_recovered()).map(|b| b.tumor.position).collect();
        assert_eq!(a, b);
        assert!(!a.is_empty());
    }

    #[test]
    fn centromere_gap_is_subtracted() {
        // chr10 (gap 40340102-41140102, 800,001 bases): sites every 100kb from 21Mb, with none in the gap.
        // A run spanning the gap is 20.6Mb of positions, 19.8Mb without the gap: not a segment.
        let mut bafs: Vec<TumorBaf> = chromosome("chr10", 700, 20_000_000, 100_000, &[]);
        bafs.retain(|b| b.tumor.position < 40_340_102 || b.tumor.position > 41_140_102);
        let n = bafs.len();
        let first = bafs.iter().position(|b| b.tumor.position >= 30_000_000).unwrap();
        // the candidate run covers sites first..=last, widened by the window: LOH sites from first - 50
        let last = bafs.iter().position(|b| b.tumor.position >= 30_000_000 + 20_600_000).unwrap();
        for b in bafs[first - 50..=last + 50].iter_mut() {
            b.tumor.ref_support = 40;
            b.tumor.alt_support = 0;
        }
        assert_eq!(bafs[last].tumor.position - bafs[first].tumor.position, 20_600_000);
        LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr10", &mut bafs);
        let candidates: Vec<usize> = (0..n).filter(|&i| bafs[i].tumor.loh_data.as_ref().unwrap().loh_candidate).collect();
        assert_eq!(candidates.first(), Some(&first));
        assert_eq!(candidates.last(), Some(&last));
        assert!(recovered(&bafs).is_empty());
        // without the gap the same run is long enough
        let mut no_gap: Vec<TumorBaf> = bafs.iter().cloned().collect();
        for b in no_gap.iter_mut() {
            b.tumor.loh_data = None;
            b.tumor.chromosome = "chr10_alt".to_string();
        }
        LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr10_alt", &mut no_gap);
        assert!(!recovered(&no_gap).is_empty());
    }

    #[test]
    fn partly_overlapping_gap() {
        // AM52: with a site inside the gap (not in the standard loci file) Java subtracts the whole gap
        // chr21 (v38 gap 11590184-12390184): run from 11,990,184 (in the gap) to 32,290,184: 20.3Mb, the
        // part of the gap inside it 400,001 bases, the whole gap 800,001
        let mut bafs: Vec<TumorBaf> = chromosome("chr21", 600, 1_990_184, 100_000, &[]);
        let first = bafs.iter().position(|b| b.tumor.position == 11_990_184).unwrap();
        let last = bafs.iter().position(|b| b.tumor.position == 32_290_184).unwrap();
        for b in bafs[first - 50..=last + 50].iter_mut() {
            b.tumor.ref_support = 40;
            b.tumor.alt_support = 0;
        }
        for java in [true, false] {
            with_java(java, || {
                let mut run = bafs.clone();
                LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr21", &mut run);
                // Java: 20.3Mb - 800,001 < 20Mb; corrected: 20.3Mb - 400,001 < 20Mb too
                assert!(recovered(&run).is_empty());
            });
        }
        // a longer run: 20.7Mb - 800,001 < 20Mb in Java, 20.7Mb - 400,001 >= 20Mb corrected
        let last = bafs.iter().position(|b| b.tumor.position == 32_690_184).unwrap();
        for b in bafs[first - 50..=last + 50].iter_mut() {
            b.tumor.ref_support = 40;
            b.tumor.alt_support = 0;
        }
        for java in [true, false] {
            with_java(java, || {
                let mut run = bafs.clone();
                LohIdentifier::new(true, 25).find_chromosome_loh_segments("chr21", &mut run);
                assert_eq!(recovered(&run).is_empty(), java);
            });
        }
    }
}
