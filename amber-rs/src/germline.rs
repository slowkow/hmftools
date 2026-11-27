//! GermlineAnalysis: the reference sample(s) at every site, and from that the heterozygous and homozygous
//! sites, the SNP-check sites, the regions of homozygosity, consanguinity and uniparental disomy.

use anyhow::Result;
use hmf_common::java_hash::{JavaHash, JavaHashSet};

use crate::config::AmberConfig;
use crate::constants::{doubles, median_positive_value};
use crate::evidence::Evidence;
use crate::roh::{consanguinity_proportion, find_uniparental_disomy, RegionOfHomozygosity, RegionOfHomozygosityFinder};
use crate::sites::{human_chromosome, ByChr, Locus};
use crate::types::AmberBase;

/// Java's Math.round(double) as an int: half up
pub fn java_round(x: f64) -> i32 {
    let f = x.floor();
    (if x - f >= 0.5 { f + 1.0 } else { f }) as i64 as i32
}

/// BaseDepthFilter: depth within a proportion of the median positive depth (off when both are 0)
pub struct BaseDepthFilter {
    min_depth: i32,
    max_depth: i32,
}

impl BaseDepthFilter {
    pub fn new<'a>(min_depth_percentage: f64, max_depth_percentage: f64, evidence: impl Iterator<Item = &'a Evidence>) -> Self {
        if min_depth_percentage == 0.0 && max_depth_percentage == 0.0 {
            return BaseDepthFilter { min_depth: 0, max_depth: 0 };
        }
        let median = median_positive_value(evidence.map(|e| e.c.read_depth));
        let min_depth = java_round(median as f64 * min_depth_percentage);
        let max_depth = java_round(median as f64 * max_depth_percentage);
        crate::info!("median normal depth({median}) reads, filtered(min={min_depth} max={max_depth})");
        BaseDepthFilter { min_depth, max_depth }
    }
    pub fn test(&self, e: &Evidence) -> bool {
        if self.min_depth == 0 && self.max_depth == 0 {
            return true;
        }
        e.c.read_depth > 0 && e.c.read_depth >= self.min_depth && e.c.read_depth <= self.max_depth
    }
}

/// NormalHeterozygousCheck through NormalHeterozygousFilter
#[derive(Clone, Copy)]
pub struct HeterozygousCheck {
    pub min_het_af_percentage: f64,
    pub max_het_af_percentage: f64,
}

impl HeterozygousCheck {
    pub fn test(&self, e: &Evidence) -> bool {
        let c = &e.c;
        c.indel_count == 0 && c.alt_support > 0 && c.ref_support > 0 && {
            let between = |count: i32, min: f64, max: f64| doubles::greater_or_equal(count as f64, min) && doubles::less_or_equal(count as f64, max);
            let depth = c.read_depth as f64;
            between(c.ref_support, (1.0 - self.max_het_af_percentage) * depth, (1.0 - self.min_het_af_percentage) * depth)
                && between(c.alt_support, self.min_het_af_percentage * depth, self.max_het_af_percentage * depth)
        }
    }
}

/// NormalHomozygousFilter (no alt support)
pub fn is_homozygous(e: &Evidence) -> bool {
    e.is_valid() && e.c.alt_support == 0
}

/// An AmberSite as BaseDepthIntersectFilter keeps it (AmberUtils.depthAsSite: no SNP check flag)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SiteKey {
    pub chr: &'static str,
    pub pos: i32,
    pub ref_: AmberBase,
    pub alt: AmberBase,
}

impl SiteKey {
    pub fn of(e: &Evidence) -> SiteKey {
        SiteKey { chr: e.chr, pos: e.pos, ref_: e.ref_, alt: e.alt }
    }
}

impl JavaHash for SiteKey {
    /// AmberSite.hashCode with the SNP check flag false
    fn java_hash(&self) -> i32 {
        crate::types::AmberSite::new(self.chr, self.pos, self.ref_.as_str(), self.alt.as_str(), false).java_hash()
    }
}

/// BaseDepthIntersectFilter: the sites heterozygous in every reference
pub struct IntersectFilter {
    additional: bool,
    intersection: JavaHashSet<SiteKey>,
}

impl IntersectFilter {
    pub fn new() -> Self {
        IntersectFilter { additional: false, intersection: JavaHashSet::new() }
    }
    pub fn additional<'a>(&mut self, evidence: impl Iterator<Item = &'a Evidence>) {
        // Collectors.toSet(): a new HashSet, added to in stream order
        let mut set: JavaHashSet<SiteKey> = JavaHashSet::new();
        for e in evidence {
            set.insert(SiteKey::of(e));
        }
        if self.additional {
            self.intersection.retain(|k| set.contains(k));
        } else {
            self.additional = true;
            for k in set.iter() {
                self.intersection.insert(k.clone());
            }
        }
    }
    pub fn test(&self, e: &Evidence) -> bool {
        !self.additional || self.intersection.contains(&SiteKey::of(e))
    }
    /// sites(): the intersection by chromosome, sorted by position (a stable sort of the set's order)
    pub fn sites(&self) -> ByChr<Locus> {
        let mut result: ByChr<Locus> = ByChr::new();
        for k in self.intersection.iter() {
            result.put(
                human_chromosome(k.chr),
                Locus { chr: k.chr, pos: k.pos, ref_: k.ref_, alt: k.alt, snp_check: false, gnomad_frequency: 0.5 },
            );
        }
        for list in result.lists.iter_mut() {
            list.sort_by_key(|l| l.pos);
        }
        result
    }
}

impl Default for IntersectFilter {
    fn default() -> Self {
        Self::new()
    }
}

pub struct GermlineAnalysis {
    pub snp_checked: ByChr<Evidence>,
    pub homozygous: ByChr<Evidence>,
    pub heterozygous: ByChr<Evidence>,
    pub regions_of_homozygosity: Vec<RegionOfHomozygosity>,
    pub consanguinity_proportion: f64,
    pub uniparental_disomy: Option<usize>,
}

/// germlineDepth: the reads of one reference BAM at the sites
pub fn germline_depth(config: &AmberConfig, bam: &str, sites: &ByChr<Locus>) -> Result<ByChr<Evidence>> {
    crate::info!("processing {} Amber sites in reference bam({bam})", sites.len());
    let mut evidence = sites.map(|s| Evidence::new(s.chr, s.pos, s.ref_, s.alt));
    crate::bam::process_bam(config, bam, &mut evidence)?;
    Ok(evidence)
}

impl GermlineAnalysis {
    pub fn run(config: &AmberConfig, sites: &ByChr<Locus>) -> Result<GermlineAnalysis> {
        let het_check = HeterozygousCheck { min_het_af_percentage: config.min_het_af_percent, max_het_af_percentage: config.max_het_af_percent };
        // SnpCheckFilter: the sites the loci file marks for the SNP check
        let snp_sites: std::collections::HashSet<SiteKey> =
            sites.values().filter(|s| s.snp_check).map(|s| SiteKey { chr: s.chr, pos: s.pos, ref_: s.ref_, alt: s.alt }).collect();

        let unfiltered = germline_depth(config, &config.reference_bams[0], sites)?;
        let depth_filter = BaseDepthFilter::new(config.min_depth_percent, config.max_depth_percent, unfiltered.values());
        let snp_checked = unfiltered.filter(|e| snp_sites.contains(&SiteKey::of(e)));
        let homozygous = unfiltered.filter(|e| depth_filter.test(e) && is_homozygous(e) && e.is_valid());
        let primary_het = unfiltered.filter(|e| depth_filter.test(e) && het_check.test(e) && e.is_valid());
        let mut intersect = IntersectFilter::new();
        intersect.additional(primary_het.values());

        for bam in config.reference_bams.iter().skip(1) {
            let additional = germline_depth(config, bam, &intersect.sites())?;
            let filter = BaseDepthFilter::new(config.min_depth_percent, config.max_depth_percent, additional.values());
            intersect.additional(additional.values().filter(|e| filter.test(e) && het_check.test(e) && e.is_valid()));
        }

        let heterozygous = if config.write_unfiltered_germline { unfiltered.clone() } else { primary_het.filter(|e| intersect.test(e)) };
        crate::info!("{} heterozygous, {} homozygous in reference bams", heterozygous.len(), homozygous.len());

        let roh_finder = RegionOfHomozygosityFinder::new(config.v38, config.min_depth_percent, config.max_depth_percent);
        let regions_of_homozygosity = roh_finder.find_regions(&unfiltered);
        let consanguinity_proportion = consanguinity_proportion(&regions_of_homozygosity);
        let uniparental_disomy = find_uniparental_disomy(&regions_of_homozygosity);
        Ok(GermlineAnalysis { snp_checked, homozygous, heterozygous, regions_of_homozygosity, consanguinity_proportion, uniparental_disomy })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::Counts;

    fn ev(chr: &'static str, pos: i32, depth: i32, r: i32, a: i32) -> Evidence {
        let mut e = Evidence::new(chr, pos, AmberBase::A, AmberBase::T);
        e.c = Counts { read_depth: depth, ref_support: r, alt_support: a, ..Counts::default() };
        e
    }

    /// BaseDepthIntersectFilterTest
    #[test]
    fn intersect_filter() {
        let (first, second, third) = (ev("1", 11, 5, 3, 2), ev("2", 22, 5, 3, 2), ev("3", 33, 5, 3, 2));
        let mut f = IntersectFilter::new();
        f.additional([&first, &second, &third].into_iter());
        assert!(f.test(&first) && f.test(&second) && f.test(&third));

        let mut f = IntersectFilter::new();
        f.additional([&first, &second].into_iter());
        f.additional([&third, &second].into_iter());
        assert!(!f.test(&first) && f.test(&second) && !f.test(&third));

        let mut f = IntersectFilter::new();
        f.additional([&first, &second].into_iter());
        f.additional([&third, &second].into_iter());
        f.additional([&first, &third].into_iter());
        assert!(!f.test(&first) && !f.test(&second) && !f.test(&third));
        assert!(IntersectFilter::new().test(&first));
    }

    #[test]
    fn depth_filter() {
        let evs = [ev("1", 1, 10, 0, 0), ev("1", 2, 20, 0, 0), ev("1", 3, 30, 0, 0), ev("1", 4, 0, 0, 0)];
        let f = BaseDepthFilter::new(0.5, 1.5, evs.iter());
        // median 20: 10 to 30
        assert!(evs[..3].iter().all(|e| f.test(e)) && !f.test(&evs[3]) && !f.test(&ev("1", 5, 31, 0, 0)));
        let off = BaseDepthFilter::new(0.0, 0.0, evs.iter());
        assert!(off.test(&evs[3]));
        assert_eq!(java_round(2.5), 3);
        assert_eq!(java_round(0.49999999999999994), 0);
        assert_eq!(java_round(-2.5), -2);
    }

    #[test]
    fn heterozygous_check() {
        let c = HeterozygousCheck { min_het_af_percentage: 0.4, max_het_af_percentage: 0.65 };
        assert!(c.test(&ev("1", 1, 100, 50, 50)));
        assert!(c.test(&ev("1", 1, 100, 60, 40)));
        assert!(!c.test(&ev("1", 1, 100, 61, 39)));
        assert!(c.test(&ev("1", 1, 100, 35, 65)));
        assert!(!c.test(&ev("1", 1, 100, 34, 66)));
        assert!(!c.test(&ev("1", 1, 100, 100, 0)));
        let mut e = ev("1", 1, 100, 50, 50);
        e.c.indel_count = 1;
        assert!(!c.test(&e));
    }
}
