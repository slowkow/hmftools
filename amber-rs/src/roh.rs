//! RegionOfHomozygosityFinder, RegionOfHomozygosity and ConsanguinityAnalyser: stretches of germline sites
//! without heterozygous ones, and the consanguinity and uniparental disomy measures taken from them.

use std::collections::HashMap;

use crate::constants::*;
use crate::evidence::Evidence;
use crate::germline::BaseDepthFilter;
use crate::poisson::Poisson;
use crate::sites::{bed_resource, human_chromosome_name, strip_chr, ByChr, CENTROMERES_37, CENTROMERES_38, ROH_EXCLUDED_37, ROH_EXCLUDED_38};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Zygosity {
    Heterozygous,
    Homozygous,
}

#[derive(Clone, Copy, Debug)]
pub struct LocusZygosity {
    pub position: i32,
    pub zygosity: Zygosity,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionOfHomozygosity {
    /// HumanChromosome index (1-24); written without the chr prefix, as Java writes it
    pub chromosome: usize,
    pub start: i32,
    pub end: i32,
    pub num_homozygous: i32,
    pub num_heterozygous: i32,
    pub num_unclear: i32,
}

impl RegionOfHomozygosity {
    pub fn snp_count(&self) -> i32 {
        self.num_homozygous + self.num_heterozygous + self.num_unclear
    }
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
    pub fn chromosome_name(&self) -> &'static str {
        human_chromosome_name(self.chromosome)
    }
}

const EXCLUDED_REGION_EXPAND_SIZE: i32 = 1_000_000;
const HALF_CENTROMERE_SIZE: i32 = 1_500_000;

pub struct RegionOfHomozygosityFinder {
    /// chromosome name, start, end
    excluded: Vec<(String, i32, i32)>,
    min_depth_percent: f64,
    max_depth_percent: f64,
    min_homozygous_region_size: i32,
    min_snp_loci_count: i32,
    site_window_size: i32,
    max_het_in_window: i32,
}

/// isAlleleHomozygous, the definition PURPLE uses
pub fn is_allele_homozygous(total_count: i32, allele_count: i32) -> bool {
    if total_count == allele_count {
        return true;
    }
    if allele_count as f64 > 0.75 * total_count as f64 {
        let p = Poisson::new(total_count as f64 / 2.0).cumulative_probability(total_count - allele_count);
        return p < 0.005;
    }
    false
}

impl RegionOfHomozygosityFinder {
    pub fn new(v38: bool, min_depth_percent: f64, max_depth_percent: f64) -> Self {
        Self::with_params(
            v38,
            min_depth_percent,
            max_depth_percent,
            HOMOZYGOUS_REGION_MIN_SIZE,
            HOMOZYGOUS_REGION_MIN_SNP_LOCI_COUNT,
            HOMOZYGOUS_REGION_WINDOW_SIZE,
            HOMOZYGOUS_REGION_MAX_HET_IN_WINDOW,
        )
    }

    pub fn with_params(
        v38: bool, min_depth_percent: f64, max_depth_percent: f64, min_homozygous_region_size: i32, min_snp_loci_count: i32,
        site_window_size: i32, max_het_in_window: i32,
    ) -> Self {
        // loadExcludedRegions: the resource BED and 1.5 Mb either side of each centromere
        let mut excluded = bed_resource(if v38 { ROH_EXCLUDED_38 } else { ROH_EXCLUDED_37 });
        for line in (if v38 { CENTROMERES_38 } else { CENTROMERES_37 }).lines().filter(|l| !l.is_empty()) {
            let (chr, pos) = line.split_once('\t').unwrap();
            let pos: i32 = pos.trim().parse().unwrap();
            // HumanChromosome.toString: no prefix
            excluded.push((strip_chr(chr).to_string(), pos - HALF_CENTROMERE_SIZE, pos + HALF_CENTROMERE_SIZE));
        }
        RegionOfHomozygosityFinder {
            excluded,
            min_depth_percent,
            max_depth_percent,
            min_homozygous_region_size,
            min_snp_loci_count,
            site_window_size,
            max_het_in_window,
        }
    }

    /// findRegions: the depth-filtered sites of each autosome
    pub fn find_regions(&self, base_depths: &ByChr<Evidence>) -> Vec<RegionOfHomozygosity> {
        let depth_filter = BaseDepthFilter::new(self.min_depth_percent, self.max_depth_percent, base_depths.values());
        let mut cache: HashMap<(i32, i32), bool> = HashMap::new();
        let mut is_hom = |total: i32, allele: i32| *cache.entry((total, allele)).or_insert_with(|| is_allele_homozygous(total, allele));
        let mut regions = Vec::new();
        for (ci, list) in base_depths.chromosomes() {
            if ci >= 23 {
                continue; // allosomes
            }
            let sites: Vec<LocusZygosity> = list
                .iter()
                .filter(|e| depth_filter.test(e))
                .map(|e| {
                    let hom = is_hom(e.c.read_depth, e.c.ref_support) || is_hom(e.c.read_depth, e.c.alt_support);
                    LocusZygosity { position: e.pos, zygosity: if hom { Zygosity::Homozygous } else { Zygosity::Heterozygous } }
                })
                .collect();
            regions.extend(self.find_regions_for_chromosome(ci, &sites));
        }
        regions.sort_by(|a, b| a.chromosome.cmp(&b.chromosome).then(a.start.cmp(&b.start)).then(a.end.cmp(&b.end)));
        regions
    }

    pub fn find_regions_for_chromosome(&self, chromosome: usize, sites: &[LocusZygosity]) -> Vec<RegionOfHomozygosity> {
        let chr_name = human_chromosome_name(chromosome);
        let mut regions = Vec::new();
        let mut sliding_window_start: i64 = -1;
        let mut roh_start: i64 = -1;
        let mut roh_end: i64 = -1;
        let mut num_het_in_window = 0;
        for (i, site) in sites.iter().enumerate() {
            let i = i as i64;
            let position = site.position;
            if roh_start != -1 {
                let from = if roh_end != -1 { sites[roh_end as usize].position } else { position };
                let crossed_excluded = self.overlaps_excluded(chr_name, from, position);
                if site.zygosity == Zygosity::Heterozygous {
                    num_het_in_window += 1;
                }
                if sliding_window_start <= i - self.site_window_size as i64 {
                    if sites[sliding_window_start as usize].zygosity == Zygosity::Heterozygous {
                        num_het_in_window -= 1;
                    }
                    sliding_window_start += 1;
                }
                if num_het_in_window > self.max_het_in_window || crossed_excluded {
                    if let Some(r) = self.create_region_if_pass(sites, chromosome, roh_start, roh_end) {
                        regions.push(r);
                    }
                    roh_start = -1;
                    roh_end = -1;
                }
            }
            if site.zygosity == Zygosity::Homozygous && !self.overlaps_excluded(chr_name, position, position) {
                if roh_start != -1 {
                    roh_end = i;
                } else {
                    roh_start = i;
                    roh_end = i;
                    sliding_window_start = i;
                    num_het_in_window = 0;
                }
            }
        }
        if let Some(r) = self.create_region_if_pass(sites, chromosome, roh_start, roh_end) {
            regions.push(r);
        }
        regions
    }

    /// overlapWithExcludedRegions, with each region widened by 1 Mb.
    /// AM1: Java compares the chromosome name of HumanChromosome (no prefix) with the BED's names, which
    /// for 38 carry the chr prefix, so in 38 none of the BED regions (the chr1 and chr9 heterochromatin,
    /// the chr16 gap and the short arms of 13, 14, 15, 21, 22) is ever excluded; only the centromere
    /// regions, named without a prefix, are. Corrected: names are compared without the prefix.
    fn overlaps_excluded(&self, chromosome: &str, p1: i32, p2: i32) -> bool {
        let java = crate::compat::java("AM1");
        self.excluded.iter().any(|(chr, start, end)| {
            let same = if java { chr == chromosome } else { strip_chr(chr) == strip_chr(chromosome) };
            same && p1 <= end + EXCLUDED_REGION_EXPAND_SIZE && p2 >= start - EXCLUDED_REGION_EXPAND_SIZE
        })
    }

    /// createRegionIfPassFilter.
    /// AM4: Java returns no region when the start index is 0 (`regionStartIndex <= 0`), so a stretch that
    /// begins at the first site of a chromosome is never reported, however long; -1 marks no region.
    /// Corrected: only -1 means no region.
    fn create_region_if_pass(&self, sites: &[LocusZygosity], chromosome: usize, start_index: i64, end_index: i64) -> Option<RegionOfHomozygosity> {
        let none = if crate::compat::java("AM4") { start_index <= 0 || end_index <= 0 } else { start_index < 0 || end_index < 0 };
        if none {
            return None;
        }
        let (mut hom, mut het) = (0, 0);
        for s in &sites[start_index as usize..=end_index as usize] {
            match s.zygosity {
                Zygosity::Homozygous => hom += 1,
                Zygosity::Heterozygous => het += 1,
            }
        }
        if hom + het < self.min_snp_loci_count {
            return None;
        }
        let (start, end) = (sites[start_index as usize].position, sites[end_index as usize].position);
        if end - start >= self.min_homozygous_region_size {
            return Some(RegionOfHomozygosity { chromosome, start, end, num_homozygous: hom, num_heterozygous: het, num_unclear: 0 });
        }
        None
    }
}

/// ConsanguinityAnalyser.calcConsanguinityProportion: the length of the regions over 3 Mb over 2.88 Gb
pub fn consanguinity_proportion(rohs: &[RegionOfHomozygosity]) -> f64 {
    let mut sum = 0.0;
    for r in rohs {
        if r.length() > HOMOZYGOUS_REGION_LONG_SIZE {
            sum += r.length() as f64;
        }
    }
    sum / 2_880_000_000i64 as f64
}

/// ConsanguinityAnalyser.findUniparentalDisomy: the one chromosome holding every region of 3 Mb or more,
/// when they add up to 10 Mb
pub fn find_uniparental_disomy(rohs: &[RegionOfHomozygosity]) -> Option<usize> {
    let mut chromosome: Option<usize> = None;
    let mut sum: i64 = 0;
    for r in rohs {
        if r.length() >= HOMOZYGOUS_REGION_LONG_SIZE {
            match chromosome {
                None => chromosome = Some(r.chromosome),
                Some(c) if c != r.chromosome => return None,
                _ => {}
            }
            sum += r.length() as i64;
        }
    }
    if sum >= UNIPARENTAL_DISOMY_MIN_LENGTH { chromosome } else { None }
}

#[cfg(test)]
mod tests {
    use super::*;
    use Zygosity::{Heterozygous as HET, Homozygous as HOM};

    fn site(position: i32, zygosity: Zygosity) -> LocusZygosity {
        LocusZygosity { position, zygosity }
    }

    fn finder(window: i32) -> RegionOfHomozygosityFinder {
        RegionOfHomozygosityFinder::with_params(false, 0.5, 1.5, 200, 5, window, 3)
    }

    #[test]
    fn find_region_sample() {
        let f = finder(50);
        let mut sites = vec![
            site(10, HET),
            site(13, HOM),
            site(100, HOM),
            site(140, HOM),
            site(141, HOM),
            site(200, HOM),
            site(219, HOM),
            site(222, HET),
            site(223, HET),
            site(225, HET),
        ];
        let r = f.find_regions_for_chromosome(1, &sites);
        assert_eq!(r.len(), 1);
        assert_eq!((r[0].start, r[0].end), (13, 219));
        sites.extend([
            site(226, HET),
            site(250, HOM),
            site(260, HOM),
            site(270, HOM),
            site(280, HOM),
            site(470, HOM),
            site(550, HET),
            site(570, HET),
            site(600, HET),
        ]);
        let r = f.find_regions_for_chromosome(1, &sites);
        assert_eq!(r.len(), 2);
        assert_eq!((r[0].start, r[0].end, r[1].start, r[1].end), (13, 219, 250, 470));
    }

    #[test]
    fn find_region_hets_in_window() {
        let f = finder(100);
        let mut sites = vec![
            site(10, HET),
            site(13, HOM),
            site(100, HOM),
            site(110, HET),
            site(140, HOM),
            site(150, HET),
            site(200, HOM),
            site(209, HET),
            site(219, HOM),
            site(222, HET),
            site(223, HET),
            site(225, HET),
        ];
        let r = f.find_regions_for_chromosome(1, &sites);
        assert_eq!(r.len(), 1);
        assert_eq!((r[0].start, r[0].end), (13, 219));
        sites.extend([
            site(226, HET),
            site(250, HOM),
            site(255, HET),
            site(260, HOM),
            site(265, HET),
            site(270, HOM),
            site(275, HET),
            site(280, HOM),
            site(470, HOM),
            site(550, HET),
            site(570, HET),
            site(600, HET),
        ]);
        let r = f.find_regions_for_chromosome(1, &sites);
        assert_eq!(r.len(), 2);
        assert_eq!((r[0].start, r[0].end, r[1].start, r[1].end), (13, 219, 250, 470));
        let sites = vec![
            site(10, HET),
            site(13, HOM),
            site(100, HOM),
            site(110, HET),
            site(120, HET),
            site(140, HOM),
            site(150, HET),
            site(200, HOM),
            site(209, HET),
            site(219, HOM),
            site(222, HET),
            site(223, HET),
            site(225, HET),
        ];
        assert!(f.find_regions_for_chromosome(1, &sites).is_empty());
    }

    /// AM4: a stretch from the first site is a region once corrected
    #[test]
    fn region_from_first_site() {
        let f = finder(50);
        let sites: Vec<_> = [10, 100, 200, 300, 400, 500].iter().map(|&p| site(p, HOM)).collect();
        assert_eq!(f.find_regions_for_chromosome(1, &sites).len(), 1);
        assert!(crate::compat::with_java(true, || f.find_regions_for_chromosome(1, &sites)).is_empty());
    }

    /// AM1: in 38 the BED regions (here the short arm of chr21) are excluded once corrected
    #[test]
    fn excluded_regions_in_38() {
        let f = RegionOfHomozygosityFinder::new(true, 0.5, 1.5);
        assert!(f.overlaps_excluded("21", 5_000_000, 5_000_000));
        assert!(!crate::compat::with_java(true, || f.overlaps_excluded("21", 5_000_000, 5_000_000)));
        // the centromere regions apply either way
        assert!(crate::compat::with_java(true, || f.overlaps_excluded("1", 123_605_523, 123_605_523)));
        let f37 = RegionOfHomozygosityFinder::new(false, 0.5, 1.5);
        assert!(crate::compat::with_java(true, || f37.overlaps_excluded("21", 5_000_000, 5_000_000)));
    }

    #[test]
    fn uniparental_disomy() {
        let roh = |c, s, e| RegionOfHomozygosity { chromosome: c, start: s, end: e, num_homozygous: 0, num_heterozygous: 0, num_unclear: 0 };
        let rohs = vec![roh(1, 1, 5_000_000), roh(1, 1_000_000, 6_000_000)];
        assert_eq!(find_uniparental_disomy(&rohs), Some(1));
        let rohs = vec![roh(1, 1, 5_000_000), roh(1, 1_000_000, 6_000_000), roh(2, 1, 5_000_000)];
        assert_eq!(find_uniparental_disomy(&rohs), None);
        assert!((consanguinity_proportion(&[roh(1, 1, 5_760_000)]) - 0.002).abs() < 1e-15);
    }

    #[test]
    fn allele_homozygous() {
        assert!(is_allele_homozygous(0, 0) && is_allele_homozygous(30, 30));
        assert!(!is_allele_homozygous(30, 15));
        assert!(is_allele_homozygous(40, 39));
        assert!(!is_allele_homozygous(8, 7));
    }
}
