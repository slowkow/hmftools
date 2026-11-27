//! ObservedRegion, ObservedRegionFactory (BAF and COBALT ratio accumulation per segment, germline status,
//! germline amp/del candidate regions, min-start extension) and GermlineStatusCalcs.

use std::collections::BTreeMap;

use crate::constants::*;
use crate::inputs::{AmberBaf, CobaltChromosomes, CobaltRatio, WINDOW_SIZE};
use crate::segment::SupportSegment;
use crate::types::*;

#[derive(Clone, Debug, PartialEq)]
pub struct ObservedRegion {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub ratio_support: bool,
    pub support: SegmentSupport,
    pub baf_count: i32,
    pub observed_baf: f64,
    pub depth_window_count: i32,
    pub observed_tumor_ratio: f64,
    pub observed_normal_ratio: f64,
    pub unnormalised_observed_normal_ratio: f64,
    pub germline_status: GermlineStatus,
    pub sv_cluster: bool,
    pub gc_content: f64,
    pub min_start: i32,
    pub max_start: i32,
    pub minor_allele_copy_number_deviation: f64,
    pub major_allele_copy_number_deviation: f64,
    pub deviation_penalty: f64,
    pub event_penalty: f64,
    pub ref_normalised_copy_number: f64,
    pub tumor_copy_number: f64,
    pub tumor_baf: f64,
    pub fitted_tumor_copy_number: f64,
    pub fitted_baf: f64,
}

impl ObservedRegion {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        chromosome: &str,
        start: i32,
        end: i32,
        ratio_support: bool,
        support: SegmentSupport,
        baf_count: i32,
        observed_baf: f64,
        depth_window_count: i32,
        observed_tumor_ratio: f64,
        observed_normal_ratio: f64,
        unnormalised_observed_normal_ratio: f64,
        germline_status: GermlineStatus,
        sv_cluster: bool,
        gc_content: f64,
        min_start: i32,
        max_start: i32,
    ) -> ObservedRegion {
        ObservedRegion {
            chromosome: chromosome.to_string(),
            start,
            end,
            ratio_support,
            support,
            baf_count,
            observed_baf,
            depth_window_count,
            observed_tumor_ratio,
            observed_normal_ratio,
            unnormalised_observed_normal_ratio,
            germline_status,
            sv_cluster,
            gc_content,
            min_start,
            max_start,
            minor_allele_copy_number_deviation: 0.0,
            major_allele_copy_number_deviation: 0.0,
            deviation_penalty: 0.0,
            event_penalty: 0.0,
            ref_normalised_copy_number: 0.0,
            tumor_copy_number: 0.0,
            tumor_baf: 0.0,
            fitted_tumor_copy_number: 0.0,
            fitted_baf: 0.0,
        }
    }

    pub fn major_allele_copy_number(&self) -> f64 {
        self.tumor_baf * self.tumor_copy_number
    }
    pub fn minor_allele_copy_number(&self) -> f64 {
        self.tumor_copy_number - self.major_allele_copy_number()
    }
    pub fn is_diploid_cn(cn: f64) -> bool {
        greater_or_equal(cn, MIN_DIPLOID_COPY_NUMBER) && less_or_equal(cn, MAX_DIPLOID_COPY_NUMBER)
    }
    pub fn is_diploid(&self) -> bool {
        Self::is_diploid_cn(self.major_allele_copy_number()) && Self::is_diploid_cn(self.minor_allele_copy_number())
    }
    pub fn bases(&self) -> i32 {
        1 + self.end - self.start
    }
    pub fn chr(&self) -> Chromosome {
        chr(&self.chromosome)
    }

    /// PurpleSegment.toTsv of ObservedRegion.toSegment (svCluster is written as false)
    pub fn to_tsv(&self) -> String {
        let f = |v: f64| hmf_common::fmt::decimal_format_min(v, 4, 4);
        [
            self.chromosome.clone(),
            self.start.to_string(),
            self.end.to_string(),
            self.germline_status.name().to_string(),
            "false".to_string(),
            self.baf_count.to_string(),
            f(self.observed_baf),
            f(self.minor_allele_copy_number_deviation),
            f(self.observed_tumor_ratio),
            f(self.observed_normal_ratio),
            f(self.unnormalised_observed_normal_ratio),
            f(self.major_allele_copy_number_deviation),
            f(self.deviation_penalty),
            f(self.tumor_copy_number),
            f(self.fitted_tumor_copy_number),
            f(self.fitted_baf),
            f(self.ref_normalised_copy_number),
            self.ratio_support.to_string(),
            self.support.name().to_string(),
            self.depth_window_count.to_string(),
            f(self.tumor_baf),
            f(self.gc_content),
            f(self.event_penalty),
            self.min_start.to_string(),
            self.max_start.to_string(),
        ]
        .join("\t")
    }
}

pub const SEGMENT_HEADER: &str = "chromosome\tstart\tend\tgermlineStatus\tsvCluster\tbafCount\tobservedBAF\tminorAlleleCopyNumberDeviation\tobservedTumorRatio\tobservedNormalRatio\tunnormalisedObservedNormalRatio\tmajorAlleleCopyNumberDeviation\tdeviationPenalty\ttumorCopyNumber\tfittedTumorCopyNumber\tfittedBAF\trefNormalisedCopyNumber\tratioSupport\tsupport\tdepthWindowCount\ttumorBAF\tgcContent\teventPenalty\tminStart\tmaxStart";

/// PurpleSegment.write
pub fn write_segments(path: &str, regions: &[ObservedRegion]) -> anyhow::Result<()> {
    let mut out = String::with_capacity(regions.len() * 200 + 400);
    out.push_str(SEGMENT_HEADER);
    out.push('\n');
    for r in regions {
        out.push_str(&r.to_tsv());
        out.push('\n');
    }
    std::fs::write(path, out).map_err(|e| anyhow::anyhow!("{path}: {e}"))
}

/// GermlineStatusCalcs
pub fn calc_germline_status(cobalt: &CobaltChromosomes, contig: &str, normal_ratio: f64, tumor_ratio: f64, depth_window_count: i32) -> GermlineStatus {
    let Some(c) = cobalt.get(contig) else { return GermlineStatus::Unknown };
    if depth_window_count == 0 {
        return GermlineStatus::Unknown;
    }
    let adjustment = c.actual_ratio;
    let hom = 0.1 * adjustment;
    if less_than(normal_ratio, hom) && less_than(tumor_ratio, hom) {
        return GermlineStatus::HomDeletion;
    }
    if less_than(normal_ratio, (1.0 - 0.3) * adjustment) {
        return GermlineStatus::HetDeletion;
    }
    if less_than(normal_ratio, (1.0 - 0.15) * adjustment) {
        return GermlineStatus::LikelyDiploid;
    }
    if less_than(normal_ratio, (1.0 + 0.15) * adjustment) {
        return GermlineStatus::Diploid;
    }
    if less_than(normal_ratio, (1.0 + 0.3) * adjustment) {
        return GermlineStatus::LikelyDiploid;
    }
    if less_than(normal_ratio, 2.2 * adjustment) {
        return GermlineStatus::Amplification;
    }
    GermlineStatus::Noise
}

/// The static region lists of ObservedRegionFactory.setSpecificRegions
#[derive(Clone, Debug, Default)]
pub struct SpecialRegions {
    pub excluded_immune: Vec<ChrRegion>,
    pub centromeric: Vec<ChrRegion>,
    pub germline_amp_del_exclusions: Vec<ChrRegion>,
}

/// ImmuneRegions: IG, TR and HLA regions
pub fn immune_regions(v38: bool) -> Vec<ChrRegion> {
    let flank = 5000;
    let r = ChrRegion::new;
    if v38 {
        vec![
            r("chr2", 88857161, 90315836),
            r("chr14", 105586437, 106879844),
            r("chr22", 22026076, 22922913),
            r("chr7", 38239580, 38367882),
            r("chr7", 142299177, 142812869),
            r("chr14", 21622293, 22552156),
            r("chr6", 29942532 - flank, 29945870 + flank),
            r("chr6", 31353872 - flank, 31367067 + flank),
            r("chr6", 31268749 - flank, 31272092 + flank),
        ]
    } else {
        vec![
            r("2", 89156674, 90538397),
            r("14", 106032614, 107288051),
            r("22", 22380474, 23265085),
            r("7", 38279181, 38407483),
            r("7", 141999017, 142510554),
            r("14", 22090446, 23021099),
            r("6", 29909037 - flank, 29913661 + flank),
            r("6", 31321649 - flank, 31324964 + flank),
            r("6", 31236526 - flank, 31239863 + flank),
        ]
    }
}

impl SpecialRegions {
    pub fn new(v38: bool) -> SpecialRegions {
        let coords = crate::segment::ReferenceCoords { v38 };
        let half = CENTROMERIC_WIDTH / 2;
        let mut centromeric = Vec::new();
        for c in Chromosome::ALL {
            let centromere = coords.centromere(c);
            centromeric.push(ChrRegion::new(&c.versioned(v38), centromere - half, centromere + half));
        }
        let vc = |s: &str| if v38 { format!("chr{s}") } else { s.to_string() };
        let germline_amp_del_exclusions = vec![
            ChrRegion::new(&vc("1"), 1, GERMLINE_AMP_DEL_EXCLUSION_CHR_1),
            ChrRegion::new(&vc("9"), GERMLINE_AMP_DEL_EXCLUSION_CHR_9, coords.length(chr("9"))),
            ChrRegion::new(&vc("17"), GERMLINE_AMP_DEL_EXCLUSION_CHR_17, coords.length(chr("17"))),
            ChrRegion::new(&vc("19"), 1, GERMLINE_AMP_DEL_EXCLUSION_CHR_19),
        ];
        SpecialRegions { excluded_immune: immune_regions(v38), centromeric, germline_amp_del_exclusions }
    }
}

pub struct ObservedRegionFactory<'a> {
    pub cobalt: &'a CobaltChromosomes,
    pub regions: &'a SpecialRegions,
}

/// RatioAccumulator
#[derive(Default)]
struct RatioAccumulator {
    sum: f64,
    count: i32,
    values: Vec<f64>,
}

impl RatioAccumulator {
    fn mean(&self) -> f64 {
        if self.count > 0 { self.sum / self.count as f64 } else { 0.0 }
    }
    /// the values were inserted in order (after equal ones); sorted at the end here
    fn median(&mut self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mid = self.values.len() / 2;
        if self.values.len() % 2 == 0 { (self.values[mid - 1] + self.values[mid]) * 0.5 } else { self.values[mid] }
    }
    fn add(&mut self, ratio: f64, keep: bool) -> bool {
        if !greater_than(ratio, -1.0) {
            return false;
        }
        self.count += 1;
        self.sum += ratio;
        if keep {
            self.values.push(ratio);
        }
        true
    }
}

impl<'a> ObservedRegionFactory<'a> {
    fn germline_status(&self, region: &SupportSegment, normal_ratio: f64, tumor_ratio: f64, depth_window_count: i32) -> GermlineStatus {
        let within = |list: &[ChrRegion]| list.iter().any(|x| x.chromosome == region.chromosome && positions_within(region.start, region.end, x.start, x.end));
        if within(&self.regions.excluded_immune) {
            return GermlineStatus::Excluded;
        }
        if within(&self.regions.centromeric) {
            return GermlineStatus::Centrometic;
        }
        calc_germline_status(self.cobalt, &region.chromosome, normal_ratio, tumor_ratio, depth_window_count)
    }

    fn is_amp_del_ratio(raw_normal_ratio: f64, normal_ratio: f64) -> bool {
        if normal_ratio == 0.0 {
            return false;
        }
        let ratio = raw_normal_ratio / normal_ratio;
        ratio < GERMLINE_DEL_RATIO || ratio > GERMLINE_AMP_RATIO
    }

    fn is_amp_del_candidate(&self, region: &SupportSegment, status: GermlineStatus, raw_normal_ratio: f64, normal_ratio: f64) -> bool {
        if status != GermlineStatus::Diploid || normal_ratio <= 0.0 || !is_autosome(chr(&region.chromosome)) {
            return false;
        }
        if self
            .regions
            .germline_amp_del_exclusions
            .iter()
            .any(|x| x.chromosome == region.chromosome && positions_within(region.start, region.end, x.start, x.end))
        {
            return false;
        }
        Self::is_amp_del_ratio(raw_normal_ratio, normal_ratio)
    }

    /// ObservedRegionFactory.formObservedRegions
    pub fn form_observed_regions(
        &self,
        segments: Vec<SupportSegment>,
        bafs: &BTreeMap<Chromosome, Vec<AmberBaf>>,
        ratios: &BTreeMap<Chromosome, Vec<CobaltRatio>>,
    ) -> Vec<ObservedRegion> {
        let refined = crate::segment::refine_excluded(segments, &self.regions.excluded_immune);
        let mut observed: Vec<ObservedRegion> = Vec::with_capacity(refined.len());
        let mut candidates: Vec<usize> = Vec::new();
        let window = Window { size: WINDOW_SIZE };
        // forward cursors per chromosome (the selectors only move forward through sorted positions)
        let mut baf_cursor: BTreeMap<Chromosome, usize> = BTreeMap::new();
        let mut ratio_cursor: BTreeMap<Chromosome, usize> = BTreeMap::new();
        for region in &refined {
            let c = chr(&region.chromosome);
            // BAFs
            let mut baf_values: Vec<f64> = Vec::new();
            if let Some(list) = bafs.get(&c) {
                let cur = baf_cursor.entry(c).or_insert(0);
                while *cur < list.len() && list[*cur].position < region.start {
                    *cur += 1;
                }
                while *cur < list.len() && list[*cur].position <= region.end {
                    let baf = &list[*cur];
                    if let Some(cc) = self.cobalt.get(&baf.chromosome) {
                        if cc.is_normal() && cc.is_diploid() && !baf.tumor_modified_baf().is_nan() {
                            baf_values.push(baf.tumor_modified_baf());
                        }
                    }
                    *cur += 1;
                }
            }
            let baf_count = baf_values.len() as i32;
            let median_baf = if baf_count > 0 {
                baf_values.sort_by(|a, b| dcmp(*a, *b));
                let n = baf_values.len();
                if n % 2 == 0 { (baf_values[n / 2] + baf_values[n / 2 - 1]) / 2.0 } else { baf_values[n / 2] }
            } else {
                0.0
            };
            // COBALT ratios
            let mut reference = RatioAccumulator::default();
            let mut unnormalised = RatioAccumulator::default();
            let mut tumor_gc = RatioAccumulator::default();
            let mut tumor_content = RatioAccumulator::default();
            if let Some(list) = ratios.get(&c) {
                let cur = ratio_cursor.entry(c).or_insert(0);
                // GenomePositionSelectorListImpl: back to the first position at or after the start
                while *cur > 0 && list[*cur].position >= region.start {
                    *cur -= 1;
                }
                while *cur + 1 < list.len() && list[*cur].position < region.start {
                    *cur += 1;
                }
                while *cur < list.len() && list[*cur].position >= region.start && list[*cur].position <= region.end {
                    let r = &list[*cur];
                    if window.end(r.position) <= region.end && r.reference_gc_diploid_ratio >= 0.0 {
                        reference.add(r.reference_gc_diploid_ratio, false);
                        unnormalised.add(r.reference_gc_ratio, false);
                        if tumor_gc.add(r.tumor_gc_ratio, true) {
                            tumor_content.add(r.tumor_gc_content, true);
                        }
                    }
                    *cur += 1;
                }
                *cur = (*cur).min(list.len().saturating_sub(1));
            }
            let tumor_ratio = tumor_gc.median();
            let tumor_mean_content = tumor_content.mean();
            let normal_ratio = reference.mean();
            let depth_window_count = tumor_gc.count;
            let status = self.germline_status(region, normal_ratio, tumor_ratio, depth_window_count);
            let o = ObservedRegion::new(
                &region.chromosome,
                region.start,
                region.end,
                region.ratio_support,
                region.support,
                baf_count,
                median_baf,
                depth_window_count,
                tumor_ratio,
                normal_ratio,
                unnormalised.mean(),
                status,
                region.sv_cluster,
                tumor_mean_content,
                region.min_start,
                region.max_start,
            );
            if o.start > o.end || !positions_within(region.min_start, region.max_start, o.start, o.end) {
                crate::error!(
                    "invalid observed region: loc({}:{}-{}) ratio({}) seg({}) status({}) start(min={} max={})",
                    o.chromosome,
                    o.start,
                    o.end,
                    o.ratio_support,
                    o.support,
                    o.germline_status,
                    o.min_start,
                    o.max_start
                );
            }
            let raw_normal_ratio = o.unnormalised_observed_normal_ratio;
            observed.push(o);
            if self.is_amp_del_candidate(region, status, raw_normal_ratio, normal_ratio) {
                candidates.push(observed.len() - 1);
            }
        }
        self.find_germline_amp_del_regions(&mut observed, &candidates);
        extend_min_support(&mut observed);
        observed
    }

    fn find_germline_amp_del_regions(&self, regions: &mut [ObservedRegion], candidates: &[usize]) {
        if candidates.is_empty() {
            return;
        }
        let is_candidate = |r: &ObservedRegion| Self::is_amp_del_ratio(r.unnormalised_observed_normal_ratio, r.observed_normal_ratio);
        let mut processed: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &ci in candidates {
            if processed.contains(&ci) {
                continue;
            }
            let chromosome = regions[ci].chromosome.clone();
            let (mut range_min, mut range_max, mut index_min, mut index_max) = (0, 0, 0usize, 0usize);
            let mut j = ci as i64 - 1;
            while j >= 0 {
                let r = &regions[j as usize];
                if r.chromosome != chromosome {
                    break;
                }
                if r.support == SegmentSupport::Centromere || (r.germline_status == GermlineStatus::Diploid && !is_candidate(r)) {
                    range_min = r.end;
                    break;
                } else if r.support == SegmentSupport::Telomere {
                    range_min = r.start;
                    index_min = j as usize;
                    break;
                }
                index_min = j as usize;
                j -= 1;
            }
            for j in ci + 1..regions.len() {
                let r = &regions[j];
                if r.chromosome != chromosome {
                    index_max = j - 1;
                    range_max = regions[index_max].end;
                    break;
                }
                if r.support == SegmentSupport::Centromere || (r.germline_status == GermlineStatus::Diploid && !is_candidate(r)) {
                    range_max = r.start;
                    break;
                } else if r.support == SegmentSupport::Telomere {
                    range_max = r.end;
                    index_max = j;
                    break;
                }
                index_max = j;
            }
            if range_max - range_min >= GERMLINE_DEL_MIN_LENGTH {
                for j in index_min..=index_max {
                    let r = &mut regions[j];
                    r.observed_normal_ratio = r.unnormalised_observed_normal_ratio;
                    if r.germline_status != GermlineStatus::Excluded && r.germline_status != GermlineStatus::Centrometic {
                        r.germline_status = calc_germline_status(self.cobalt, &r.chromosome, r.observed_normal_ratio, r.observed_tumor_ratio, r.depth_window_count);
                    }
                    processed.insert(j);
                }
            }
        }
    }
}

/// ObservedRegionFactory.extendMinSupport
pub fn extend_min_support(regions: &mut [ObservedRegion]) {
    for i in 0..regions.len() {
        if regions[i].support == SegmentSupport::None && regions[i].germline_status == GermlineStatus::Diploid {
            let mut j = i as i64 - 1;
            while j >= 0 {
                let (prior_status, prior_start, prior_support) = {
                    let p = &regions[j as usize];
                    (p.germline_status, p.start, p.support)
                };
                if prior_status == GermlineStatus::Diploid {
                    break;
                }
                regions[i].min_start = regions[i].min_start.min(prior_start);
                if prior_support != SegmentSupport::None {
                    break;
                }
                j -= 1;
            }
        }
    }
}

#[cfg(test)]
mod germline_status_tests {
    use super::*;
    use crate::test_support::cobalt_chromosomes;

    fn status(female: bool, chromosome: &str, ratio: f64) -> GermlineStatus {
        calc_germline_status(&cobalt_chromosomes(female), chromosome, ratio, 0.01, 1)
    }

    fn check(chromosome: &str, ratio: f64, male: GermlineStatus, female: GermlineStatus) {
        assert_eq!(status(false, chromosome, ratio), male, "male {chromosome} {ratio}");
        assert_eq!(status(true, chromosome, ratio), female, "female {chromosome} {ratio}");
    }

    // GermlineStatusCalcsTest
    #[test]
    fn autosome() {
        use GermlineStatus::*;
        for (r, s) in [
            (0.0, HomDeletion),
            (0.09, HomDeletion),
            (0.10, HetDeletion),
            (0.69, HetDeletion),
            (0.70, LikelyDiploid),
            (0.84, LikelyDiploid),
            (0.8499, LikelyDiploid),
            (0.85, Diploid),
            (0.85001, Diploid),
            (1.1499, Diploid),
            (1.15, LikelyDiploid),
            (1.15001, LikelyDiploid),
            (1.16, LikelyDiploid),
            (1.299, LikelyDiploid),
            (1.3, Amplification),
            (1.3001, Amplification),
            (2.199, Amplification),
            (2.20, Noise),
            (2.21, Noise),
        ] {
            check("1", r, s, s);
        }
    }

    #[test]
    fn x() {
        use GermlineStatus::*;
        check("X", 0.0, HomDeletion, HomDeletion);
        for (r, m, f) in [
            (0.049, HomDeletion, HomDeletion),
            (0.05, HetDeletion, HomDeletion),
            (0.09, HetDeletion, HomDeletion),
            (0.10, HetDeletion, HetDeletion),
            (0.424, LikelyDiploid, HetDeletion),
            (0.425, Diploid, HetDeletion),
            (0.574, Diploid, HetDeletion),
            (0.575, LikelyDiploid, HetDeletion),
            (0.576, LikelyDiploid, HetDeletion),
            (0.84, Amplification, LikelyDiploid),
            (0.85, Amplification, Diploid),
            (1.10, Noise, Diploid),
            (1.15, Noise, LikelyDiploid),
            (1.16, Noise, LikelyDiploid),
            (1.8, Noise, Amplification),
            (2.20, Noise, Noise),
            (2.21, Noise, Noise),
        ] {
            check("X", r, m, f);
        }
    }

    #[test]
    fn y() {
        use GermlineStatus::*;
        for (r, m) in [(0.0, HomDeletion), (0.049, HomDeletion), (0.05, HetDeletion), (0.424, LikelyDiploid), (0.425, Diploid), (0.575, LikelyDiploid), (0.576, LikelyDiploid), (0.9, Amplification), (1.10, Noise), (1.11, Noise)] {
            check("Y", r, m, Unknown);
        }
    }
}
