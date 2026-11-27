//! TumorContamination and TumorContaminationModel: from the tumor reads at germline homozygous sites,
//! the level of contamination whose Poisson model best fits the alt read counts.

use hmf_common::java_hash::JavaHashMap;

use crate::constants::*;
use crate::evidence::Counts;
use crate::poisson::Poisson;
use crate::types::AmberBase;

/// A germline homozygous site with tumor alt reads (TumorContamination; the counts as BaseDepthData).
#[derive(Clone, Debug)]
pub struct TumorContamination {
    pub chr: &'static str,
    pub pos: i32,
    pub ref_: AmberBase,
    pub alt: AmberBase,
    pub normal: Counts,
    pub tumor: Counts,
}

impl TumorContamination {
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn tumor_vaf(&self) -> f64 {
        if self.tumor.read_depth == 0 {
            return -1.0;
        }
        self.tumor.alt_support as f64 / self.tumor.read_depth as f64
    }
}

/// TumorContaminationModel.calcContamination: 0 unless enough sites have three or more alt reads.
/// The sites come in the order the alt-count map is built from, which decides only the order of a float
/// sum (see REPORT.md).
pub fn calc_contamination(sites: &[TumorContamination], amber_site_count: i64) -> f64 {
    let (mut three_plus, mut three_plus_low_vaf, mut two_plus) = (0i64, 0i64, 0i64);
    for site in sites {
        if site.tumor.alt_support >= 2 {
            two_plus += 1;
            if site.tumor.alt_support >= 3 {
                three_plus += 1;
                let vaf = site.tumor.alt_support as f64 / site.tumor.read_depth as f64;
                if vaf < CONTAMINATON_THREE_PLUS_READS_VAF_MIN {
                    three_plus_low_vaf += 1;
                }
            }
        }
    }
    let mut calc = false;
    if three_plus >= CONTAMINATON_THREE_PLUS_READS_MIN {
        if three_plus as f64 > CONTAMINATON_THREE_PLUS_READS_SITE_PERC * amber_site_count as f64
            || three_plus_low_vaf as f64 > CONTAMINATON_THREE_PLUS_READS_SITE_LOW_VAF_PERC * amber_site_count as f64
        {
            calc = true;
        }
    }
    if !calc {
        crate::debug!("siteCount({amber_site_count}) altSites(3={three_plus} 2={two_plus} 3-lowVaf={three_plus_low_vaf})");
        return 0.0;
    }
    // Collectors.groupingBy into a HashMap, in the order of the sites
    let mut alt_support_frequencies: JavaHashMap<i32, i64> = JavaHashMap::new();
    for site in sites {
        let k = site.tumor.alt_support;
        let n = alt_support_frequencies.get(&k).copied().unwrap_or(0);
        alt_support_frequencies.insert(k, n + 1);
    }
    let frequencies: Vec<(i32, i64)> = alt_support_frequencies.iter().map(|(k, v)| (*k, *v)).collect();
    let median_tumor_read_depth = median_positive_value(sites.iter().map(|s| s.tumor.read_depth));
    let level = calc_contamination_level(median_tumor_read_depth, two_plus, &frequencies);
    crate::info!("contamination level identified({}) ", hmf_common::fmt::java_fixed(level, 3));
    crate::debug!(
        "contamination({}) siteCount({amber_site_count}) medianTumorDepth({median_tumor_read_depth}) altSites(3={three_plus} 2={two_plus} 3-lowVaf={three_plus_low_vaf})",
        hmf_common::fmt::java_fixed(level, 3)
    );
    level
}

/// calcContamination(medianTumorReadDepth, twoPlusReadCount, altSupportFrequencies): the level, in steps of
/// 0.001 up to 1, with the lowest score. `alt_support_frequencies` is in the Java map's iteration order.
pub fn calc_contamination_level(median_tumor_read_depth: i32, two_plus_read_count: i64, alt_support_frequencies: &[(i32, i64)]) -> f64 {
    const INCREMENT: f64 = 0.001;
    let mut contamination = 0.0;
    let mut lowest_score = f64::MAX;
    let mut i = INCREMENT;
    while doubles::less_or_equal(i, 1.0) {
        let score = contamination_score(i, two_plus_read_count, median_tumor_read_depth as i64, alt_support_frequencies);
        if doubles::less_than(score, lowest_score) {
            lowest_score = score;
            contamination = i;
        }
        i += INCREMENT;
    }
    contamination
}

fn unadjusted_likelihood(het: &Poisson, hom_alt: &Poisson, alt_support: i32) -> f64 {
    0.5 * het.probability(alt_support) + 0.25 * hom_alt.probability(alt_support)
}

fn contamination_score(contamination: f64, two_plus_read_count: i64, median_tumor_read_depth: i64, alt_support: &[(i32, i64)]) -> f64 {
    let het = Poisson::new(0.5 * contamination * median_tumor_read_depth as f64);
    let hom_alt = Poisson::new(contamination * median_tumor_read_depth as f64);
    let alt_support_0 = unadjusted_likelihood(&het, &hom_alt, 0) + 0.25;
    let alt_support_1 = unadjusted_likelihood(&het, &hom_alt, 1);
    let alt_support_2_plus = 1.0 - alt_support_0 - alt_support_1;
    let mut total_model_percentage = 0.0;
    let mut total_difference = 0.0;
    for &(support, count) in alt_support {
        if support > 1 {
            let model_percentage = unadjusted_likelihood(&het, &hom_alt, support) / alt_support_2_plus;
            let actual_percentage = count as f64 * 1.0 / two_plus_read_count as f64;
            total_difference += (actual_percentage - model_percentage).abs();
            total_model_percentage += model_percentage;
        }
    }
    total_difference + (1.0 - total_model_percentage)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// TumorContaminationModelTest.fromFile: lines of "count altSupport", into a HashMap
    fn from_file(text: &str) -> Vec<(i32, i64)> {
        let mut map: JavaHashMap<i32, i64> = JavaHashMap::new();
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let f: Vec<&str> = line.split_whitespace().collect();
            map.insert(f[1].parse().unwrap(), f[0].parse().unwrap());
        }
        map.iter().map(|(k, v)| (*k, *v)).collect()
    }

    fn above(min: i32, map: &[(i32, i64)]) -> i64 {
        map.iter().filter(|(k, _)| *k >= min).map(|(_, v)| v).sum()
    }

    #[test]
    fn high_contamination() {
        let map = from_file(include_str!("../testdata/contamination.high"));
        assert_eq!(above(2, &map), 249329);
        assert_eq!(above(3, &map), 248936);
        let c = calc_contamination_level(107, above(2, &map), &map);
        assert!((c - 1.0).abs() < 1e-10, "{c}");
    }

    #[test]
    fn low_contamination() {
        let map = from_file(include_str!("../testdata/contamination.low"));
        assert_eq!(above(2, &map), 73058);
        assert_eq!(above(3, &map), 31170);
        let c = calc_contamination_level(99, above(2, &map), &map);
        assert!((c - 0.02).abs() < 1e-10, "{c}");
    }

    /// Against TumorContaminationModel itself: the file harness/ContaminationRef.java writes (path in
    /// AMBER_CONTAMINATION_REF). Run with `cargo test --release -- --ignored contamination_matches_java`.
    #[test]
    #[ignore]
    fn contamination_matches_java() {
        let path = std::env::var("AMBER_CONTAMINATION_REF").expect("AMBER_CONTAMINATION_REF");
        let (mut n, mut differ, mut differ_bits) = (0, 0, 0);
        for line in std::fs::read_to_string(path).unwrap().lines() {
            let f: Vec<&str> = line.split(' ').collect();
            let median: i32 = f[0].parse().unwrap();
            let two_plus: i64 = f[1].parse().unwrap();
            let map: Vec<(i32, i64)> = f[2]
                .split(',')
                .map(|kv| {
                    let (k, v) = kv.split_once(':').unwrap();
                    (k.parse().unwrap(), v.parse().unwrap())
                })
                .collect();
            let java = f64::from_bits(u64::from_str_radix(f[3], 16).unwrap());
            let rust = calc_contamination_level(median, two_plus, &map);
            n += 1;
            if rust.to_bits() != java.to_bits() {
                differ_bits += 1;
                if hmf_common::fmt::java_fixed(rust, 4) != hmf_common::fmt::java_fixed(java, 4) {
                    differ += 1;
                    eprintln!("differs: median {median} two-plus {two_plus}: java {java} rust {rust}");
                }
            }
        }
        eprintln!("{n} cases: {differ_bits} differ in any bit, {differ} in the level");
        assert_eq!(differ_bits, 0);
    }

    #[test]
    fn tumor_vaf() {
        let t = Counts { read_depth: 100, ref_support: 95, alt_support: 5, ..Counts::default() };
        let tc = TumorContamination { chr: "chr1", pos: 1000, ref_: AmberBase::T, alt: AmberBase::G, normal: Counts::default(), tumor: t };
        assert!((tc.tumor_vaf() - 0.05).abs() < 0.001);
    }

    #[test]
    fn no_contamination_below_minimum_sites() {
        let site = |alt| TumorContamination {
            chr: "chr1",
            pos: 1,
            ref_: AmberBase::A,
            alt: AmberBase::C,
            normal: Counts::default(),
            tumor: Counts { read_depth: 100, alt_support: alt, ..Counts::default() },
        };
        let sites: Vec<_> = (0..9).map(|_| site(5)).collect();
        assert_eq!(calc_contamination(&sites, 100), 0.0);
    }
}
