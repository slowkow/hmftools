//! GnomadFrequencySupplier, DefaultGnomadFrequencySupplier, GnomadFrequencyAnalysis and
//! PeakGnomadFrequenciesChecker: a peak whose captured points with VAF below 0.5, or those at 0.5 and
//! above, have a mean gnomAD frequency more than 0.2 from that of the heterozygous baseline is no
//! contamination or copy number peak.
//!
//! Findings (REPORT.md): AM54 (Java sums the frequencies in the iteration order of a HashSet of
//! PositionEvidence, whose hash codes use enum identity hashes), AM55 (the frequency is looked up by position
//! only, so at a position with two sites either allele's frequency may be taken).

use crate::types::{AmberSite, PositionEvidence};
use std::collections::HashMap;

pub const GNOMAD_FREQUENCY_TOLERANCE: f64 = 0.20;

/// GnomadFrequencySupplier
pub trait GnomadFrequencySupplier {
    fn frequency(&self, evidence: &PositionEvidence) -> f64;
}

/// a lambda (chromosome, position) -> frequency, as the Java tests use
impl<F: Fn(&str, i32) -> f64> GnomadFrequencySupplier for F {
    fn frequency(&self, evidence: &PositionEvidence) -> f64 {
        self(&evidence.chromosome, evidence.position)
    }
}

/// the sites of one chromosome, sorted by position (stably, as Arrays.sort on objects)
enum ChromosomeSites {
    /// sites[start..end] are this chromosome's, already in position order
    Range(usize, usize),
    /// indices into sites
    Indices(Vec<u32>),
}

/// DefaultGnomadFrequencySupplier: the loci file's frequencies by versioned chromosome and position
pub struct DefaultGnomadFrequencySupplier<'a> {
    sites: &'a [AmberSite],
    chromosomes: HashMap<String, ChromosomeSites>,
}

impl<'a> DefaultGnomadFrequencySupplier<'a> {
    /// `sites` in loci-file order (Java: the ListMultimap<Chromosome, AmberSite>, whose lists keep file order)
    pub fn new(sites: &'a [AmberSite], v38: bool) -> Self {
        // runs of one chromosome name; one run per chromosome, sorted, is the usual case and needs no copy
        let mut runs: Vec<(String, usize, usize)> = Vec::new();
        let mut i = 0;
        while i < sites.len() {
            let chr = &sites[i].chromosome;
            let mut j = i + 1;
            while j < sites.len() && sites[j].chromosome == *chr {
                j += 1;
            }
            runs.push((hmf_common::regions::versioned_chromosome(v38, chr), i, j));
            i = j;
        }
        let mut grouped: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
        for (name, start, end) in runs {
            grouped.entry(name).or_default().push((start, end));
        }
        let mut chromosomes = HashMap::new();
        for (name, ranges) in grouped {
            let sorted = |s: usize, e: usize| sites[s..e].windows(2).all(|w| w[0].position <= w[1].position);
            let entry = if ranges.len() == 1 && sorted(ranges[0].0, ranges[0].1) {
                ChromosomeSites::Range(ranges[0].0, ranges[0].1)
            } else {
                let mut indices: Vec<u32> = ranges.iter().flat_map(|&(s, e)| (s as u32)..(e as u32)).collect();
                indices.sort_by_key(|&k| sites[k as usize].position);
                ChromosomeSites::Indices(indices)
            };
            chromosomes.insert(name, entry);
        }
        DefaultGnomadFrequencySupplier { sites, chromosomes }
    }

    fn site_at(&self, entry: &ChromosomeSites, k: usize) -> usize {
        match entry {
            ChromosomeSites::Range(start, _) => start + k,
            ChromosomeSites::Indices(indices) => indices[k] as usize,
        }
    }

    /// the frequency at a position, as Java's binary search finds it; panics where Java throws
    pub fn frequency_at(&self, chromosome: &str, position: i32) -> f64 {
        self.lookup(chromosome, position, None)
    }

    fn lookup(&self, chromosome: &str, position: i32, alleles: Option<(&str, &str)>) -> f64 {
        if let Some(entry) = self.chromosomes.get(chromosome) {
            let len = match entry {
                ChromosomeSites::Range(s, e) => e - s,
                ChromosomeSites::Indices(v) => v.len(),
            };
            let (mut low, mut high) = (0i64, len as i64 - 1);
            while low <= high {
                let mid = ((low + high) as u64 >> 1) as i64;
                let mid_site = &self.sites[self.site_at(entry, mid as usize)];
                if mid_site.position < position {
                    low = mid + 1;
                } else if mid_site.position > position {
                    high = mid - 1;
                } else {
                    if let Some((r, a)) = alleles {
                        if !(mid_site.ref_ == r && mid_site.alt == a) {
                            // AM55: another site at this position with the evidence's alleles
                            let at = |k: i64| &self.sites[self.site_at(entry, k as usize)];
                            let mut k = mid - 1;
                            while k >= 0 && at(k).position == position {
                                if at(k).ref_ == r && at(k).alt == a {
                                    return at(k).gnomad_frequency;
                                }
                                k -= 1;
                            }
                            let mut k = mid + 1;
                            while k < len as i64 && at(k).position == position {
                                if at(k).ref_ == r && at(k).alt == a {
                                    return at(k).gnomad_frequency;
                                }
                                k += 1;
                            }
                        }
                    }
                    return mid_site.gnomad_frequency;
                }
            }
        }
        panic!("java.lang.IllegalArgumentException: No Gnomad data for: {chromosome}:{position}");
    }
}

impl GnomadFrequencySupplier for DefaultGnomadFrequencySupplier<'_> {
    fn frequency(&self, evidence: &PositionEvidence) -> f64 {
        if crate::compat::java("AM55") {
            self.lookup(&evidence.chromosome, evidence.position, None)
        } else {
            self.lookup(&evidence.chromosome, evidence.position, Some((evidence.ref_.as_str(), evidence.alt.as_str())))
        }
    }
}

/// GnomadFrequencyAnalysis.getMeanFrequency: DoubleStream.average(), 0 for no positions
pub fn mean_frequency<'e>(supplier: &dyn GnomadFrequencySupplier, positions: impl IntoIterator<Item = &'e PositionEvidence>) -> f64 {
    let mut count = 0usize;
    let sum = hmf_common::java_math::double_stream_sum(positions.into_iter().map(|p| {
        count += 1;
        supplier.frequency(p)
    }));
    if count == 0 { 0.0 } else { sum / count as f64 }
}

/// DescriptiveStatistics.getMean (commons-math3 Mean.evaluate): the sum divided by n, plus the mean
/// residual; NaN for no values
pub fn descriptive_mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let n = values.len() as f64;
    let mut sum = 0.0;
    for v in values {
        sum += v;
    }
    let xbar = sum / n;
    let mut correction = 0.0;
    for v in values {
        correction += v - xbar;
    }
    xbar + correction / n
}

/// The outcome of PeakGnomadFrequenciesChecker.checkGnomadFrequencies, with the two means it compares.
#[derive(Clone, Copy, Debug, Default)]
pub struct GnomadCheck {
    pub low_mean: f64,
    pub high_mean: f64,
    pub low_count: usize,
    pub high_count: usize,
    pub pass: bool,
}

/// PeakGnomadFrequenciesChecker.checkGnomadFrequencies over a peak's captured points, in the order given
/// (Java: the HashSet's)
pub fn check_gnomad_frequencies<'e>(supplier: &dyn GnomadFrequencySupplier, points: impl IntoIterator<Item = &'e PositionEvidence>, expected_mean: f64) -> GnomadCheck {
    let mut low = Vec::new();
    let mut high = Vec::new();
    for evidence in points {
        let frequency = supplier.frequency(evidence);
        if evidence.vaf() < 0.5 {
            low.push(frequency);
        } else {
            high.push(frequency);
        }
    }
    let low_mean = descriptive_mean(&low);
    let high_mean = descriptive_mean(&high);
    // Java: return false if either is empty, or if abs(mean - expected) > TOLERANCE (which a NaN passes)
    let pass = !(low.is_empty() || high.is_empty())
        && !((low_mean - expected_mean).abs() > GNOMAD_FREQUENCY_TOLERANCE)
        && !((high_mean - expected_mean).abs() > GNOMAD_FREQUENCY_TOLERANCE);
    GnomadCheck { low_mean, high_mean, low_count: low.len(), high_count: high.len(), pass }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compat::with_java;
    use crate::types::AmberBase;

    fn site(chromosome: &str, position: i32, frequency: f64) -> AmberSite {
        AmberSite::with_frequency(chromosome, position, "A", "G", false, frequency)
    }

    fn at(chromosome: &str, position: i32) -> PositionEvidence {
        PositionEvidence::new(chromosome, position, AmberBase::A, AmberBase::G)
    }

    // DefaultGnomadFrequencySupplierTest
    #[test]
    fn get_frequency() {
        let sites = vec![
            site("chr1", 1_000, 0.01),
            site("chr1", 2_000, 0.02),
            site("chr1", 3_000, 0.03),
            site("chr1", 4_000, 0.04),
            site("chr1", 40_000, 0.041),
            site("chr1", 50_000, 0.042),
            site("chr1", 60_000, 0.043),
            site("chr1", 70_000, 0.044),
            site("chr1", 80_000, 0.043),
            site("chr2", 2_000, 0.02),
            site("chr3", 3_000, 0.24),
        ];
        let supplier = DefaultGnomadFrequencySupplier::new(&sites, true);
        assert_eq!(supplier.frequency(&at("chr1", 1000)), 0.01);
        assert_eq!(supplier.frequency(&at("chr1", 2000)), 0.02);
        assert_eq!(supplier.frequency(&at("chr1", 3000)), 0.03);
        assert_eq!(supplier.frequency(&at("chr1", 50_000)), 0.042);
        assert_eq!(supplier.frequency(&at("chr1", 60_000)), 0.043);
        assert_eq!(supplier.frequency(&at("chr1", 80_000)), 0.043);
        assert_eq!(supplier.frequency(&at("chr3", 3000)), 0.24);
    }

    #[test]
    fn unsorted_and_split_chromosomes() {
        // the sites of a chromosome need not be together or sorted; v37 names are looked up without "chr"
        let sites = vec![site("1", 3_000, 0.3), site("2", 1_000, 0.9), site("1", 1_000, 0.1), site("1", 2_000, 0.2)];
        let supplier = DefaultGnomadFrequencySupplier::new(&sites, false);
        assert_eq!(supplier.frequency_at("1", 1000), 0.1);
        assert_eq!(supplier.frequency_at("1", 2000), 0.2);
        assert_eq!(supplier.frequency_at("1", 3000), 0.3);
        assert_eq!(supplier.frequency_at("2", 1000), 0.9);
    }

    #[test]
    #[should_panic(expected = "No Gnomad data for: chr1:1000")]
    fn unknown_position_panics() {
        let sites: Vec<AmberSite> = Vec::new();
        DefaultGnomadFrequencySupplier::new(&sites, true).frequency_at("chr1", 1000);
    }

    #[test]
    #[should_panic(expected = "No Gnomad data for: chr1:2000")]
    fn position_not_in_sites_panics() {
        let sites = vec![site("chr1", 1_000, 0.01), site("chr1", 3_000, 0.03)];
        DefaultGnomadFrequencySupplier::new(&sites, true).frequency_at("chr1", 2000);
    }

    #[test]
    #[should_panic(expected = "No Gnomad data for: chrY:1000")]
    fn unknown_chromosome_panics() {
        let sites = vec![site("chr1", 1_000, 0.01)];
        DefaultGnomadFrequencySupplier::new(&sites, true).frequency_at("chrY", 1000);
    }

    #[test]
    fn two_sites_at_one_position() {
        // AM55: Java's binary search lands on the middle site, whatever its alleles
        let sites = vec![
            site("chr1", 1_000, 0.1),
            AmberSite::with_frequency("chr1", 2_000, "A", "C", false, 0.7),
            AmberSite::with_frequency("chr1", 2_000, "A", "T", false, 0.2),
            site("chr1", 3_000, 0.3),
        ];
        let supplier = DefaultGnomadFrequencySupplier::new(&sites, true);
        let e = PositionEvidence::new("chr1", 2_000, AmberBase::A, AmberBase::T);
        assert_eq!(with_java(true, || supplier.frequency(&e)), 0.7);
        assert_eq!(with_java(false, || supplier.frequency(&e)), 0.2);
    }

    // GnomadFrequencyAnalysisTest
    #[test]
    fn mean_test() {
        let frequencies = |_: &str, position: i32| position as f64 / 1000.0;
        let data = [at("chr1", 100), at("chr2", 100), at("chr2", 400)];
        assert!((mean_frequency(&frequencies, data.iter()) - 0.2).abs() < 0.001);
        assert_eq!(mean_frequency(&frequencies, [].iter()), 0.0);
    }

    #[test]
    fn descriptive_mean_two_passes() {
        assert!(descriptive_mean(&[]).is_nan());
        assert_eq!(descriptive_mean(&[0.5]), 0.5);
        let values = [0.1, 0.2, 0.3, 0.4];
        let naive = values.iter().sum::<f64>() / 4.0;
        let mean = descriptive_mean(&values);
        assert!((mean - naive).abs() < 1e-15);
    }
}
