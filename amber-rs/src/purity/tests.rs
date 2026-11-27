//! Ports of TumorOnlyPurityAnalysisTest and PeakGnomadFrequenciesCheckerTest, and tests of the pieces
//! that join them (the Java HashSet order of a peak's points).

use super::*;
use crate::types::AmberBase;

/// a small deterministic generator in place of the Java test's Math.random and Collections.shuffle
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn bases(&mut self) -> (AmberBase, AmberBase) {
        let mut bases = [AmberBase::A, AmberBase::C, AmberBase::G, AmberBase::T];
        for i in (1..4).rev() {
            let j = (self.next() * (i + 1) as f64) as usize;
            bases.swap(i, j);
        }
        (bases[0], bases[1])
    }
}

const BASELINE_READINGS_PER_ARM: i32 = 100;
const BASELINE_DEPTH: i32 = 1000;

struct Scenario {
    rng: Rng,
    /// per chromosome ordinal: (evidence, gnomad frequency), kept sorted by position as Java's TreeSet
    evidence: Vec<Vec<(PositionEvidence, f64)>>,
    sites: Vec<AmberSite>,
}

fn chr(ordinal: usize) -> String {
    format!("chr{}", hmf_common::regions::HUMAN_CHROMOSOMES[ordinal])
}

impl Scenario {
    fn new(seed: u64) -> Scenario {
        Scenario { rng: Rng(seed), evidence: vec![Vec::new(); 24], sites: Vec::new() }
    }

    fn alt_depth_for_germline_point(&mut self, gnomad_frequency: f64, depth: i32) -> i32 {
        let alt_count = (0..2).filter(|_| self.rng.next() < gnomad_frequency).count();
        match alt_count {
            0 => 0,
            2 => depth,
            _ => depth / 2,
        }
    }

    fn add(&mut self, ordinal: usize, e: PositionEvidence, frequency: f64) {
        self.sites.push(AmberSite::with_frequency(&e.chromosome, e.position, e.ref_.as_str(), e.alt.as_str(), false, frequency));
        let list = &mut self.evidence[ordinal];
        // TreeSet: one entry per position
        if let Err(at) = list.binary_search_by_key(&e.position, |(x, _)| x.position) {
            list.insert(at, (e, frequency));
        }
    }

    fn start_position(ordinal: usize, q: bool) -> i32 {
        if !q {
            return 1_000_000;
        }
        // as the Java test (a division where an addition was presumably meant)
        (genome::centromere(true, &chr(ordinal)) + 3_000_000) / 1_000_000
    }

    fn create_baseline_readings(&mut self) {
        for ordinal in 0..24 {
            let short_arm = matches!(ordinal, 12 | 13 | 14 | 20 | 21);
            let arms: &[bool] = if short_arm { &[true] } else { &[false, true] };
            for &q in arms {
                let mut position = Scenario::start_position(ordinal, q);
                for _ in 0..BASELINE_READINGS_PER_ARM {
                    let (r, a) = self.rng.bases();
                    let mut e = PositionEvidence::new(&chr(ordinal), position, r, a);
                    let alt = self.alt_depth_for_germline_point(0.5, BASELINE_DEPTH);
                    e.read_depth = BASELINE_DEPTH;
                    e.ref_support = BASELINE_DEPTH - alt;
                    e.alt_support = alt;
                    self.add(ordinal, e, 0.5);
                    position += 1000;
                }
            }
        }
    }

    fn create_copy_number_event_points(&mut self, vaf: f64, num_points: i32) {
        let mut position = Scenario::start_position(0, false) + 50;
        let mut low = false;
        for _ in 0..num_points {
            let (r, a) = self.rng.bases();
            let mut e = PositionEvidence::new(&chr(0), position, r, a);
            let alt = if low { (BASELINE_DEPTH as f64 * vaf) as i32 } else { (BASELINE_DEPTH as f64 * (1.0 - vaf)) as i32 };
            e.read_depth = BASELINE_DEPTH;
            e.ref_support = BASELINE_DEPTH - alt;
            e.alt_support = alt;
            self.add(0, e, 0.5);
            position += BASELINE_DEPTH;
            low = !low;
        }
    }

    fn add_in_contamination_data(&mut self, contamination: f64) {
        let depth = (BASELINE_DEPTH as f64 * contamination) as i32;
        for ordinal in 0..24 {
            for k in 0..self.evidence[ordinal].len() {
                let frequency = self.evidence[ordinal][k].1;
                let alt = self.alt_depth_for_germline_point(frequency, depth);
                let e = &mut self.evidence[ordinal][k].0;
                e.read_depth += depth;
                e.ref_support += depth - alt;
                e.alt_support += alt;
            }
        }
    }

    fn analysis(&self) -> TumorOnlyPurityAnalysis {
        let evidence: Vec<PositionEvidence> = self.evidence.iter().flat_map(|l| l.iter().map(|(e, _)| e.clone())).collect();
        let config = PurityAnalysisConfig { tumor_id: "tumor123".into(), v38: true, output_dir: String::new(), threads: 2 };
        TumorOnlyPurityAnalysis::new(&evidence, &self.sites, &config)
    }
}

#[test]
fn no_peaks() {
    let mut s = Scenario::new(1);
    s.create_baseline_readings();
    let a = s.analysis();
    assert!((a.cutoff() - 0.04).abs() < 0.001);
    assert!(a.contamination_peaks().is_empty());
}

#[test]
fn one_peak_above_minimum_cutoff() {
    let mut s = Scenario::new(2);
    s.create_baseline_readings();
    let vaf = 0.21;
    s.create_copy_number_event_points(vaf - 0.01, 20);
    s.create_copy_number_event_points(vaf, 40);
    s.create_copy_number_event_points(vaf + 0.01, 20);
    let a = s.analysis();
    assert!((a.cutoff() - MIN_CUTOFF).abs() < 0.001);
    assert_eq!(a.contamination_peaks().len(), 0);
    assert_eq!(a.copy_number_peaks().len(), 1);
    assert!((a.copy_number_peaks()[0].vaf() - vaf).abs() < 0.01);
}

#[test]
fn one_cn_event_below_minimum_cutoff() {
    let mut s = Scenario::new(3);
    s.create_baseline_readings();
    let vaf = 0.1;
    s.create_copy_number_event_points(vaf - 0.01, 20);
    s.create_copy_number_event_points(vaf, 40);
    s.create_copy_number_event_points(vaf + 0.01, 20);
    let a = s.analysis();
    assert!((a.cutoff() - vaf / 3.0).abs() < 0.01);
    assert_eq!(a.contamination_peaks().len(), 0);
    assert_eq!(a.copy_number_peaks().len(), 2);
    let peak0 = a.copy_number_peaks()[0].stats().unwrap();
    assert!((peak0.homozygous_proportion - 1.0).abs() < 0.01);
    assert!((a.copy_number_peaks()[0].vaf() - vaf).abs() < 0.01);
    assert_eq!((peak0.homozygous_count, peak0.heterozygous_count), (40, 0));
    let peak1 = a.copy_number_peaks()[1].stats().unwrap();
    assert!((a.copy_number_peaks()[1].vaf() - vaf * 2.0).abs() < 0.1);
    assert!(peak1.homozygous_proportion.abs() < 0.01);
    assert_eq!((peak1.homozygous_count, peak1.heterozygous_count), (0, 40));
}

#[test]
fn contamination() {
    let mut s = Scenario::new(4);
    s.create_baseline_readings();
    s.add_in_contamination_data(0.25);
    let a = s.analysis();
    assert_eq!(a.contamination_peaks().len(), 1);
    assert!(a.contamination() > 0.0);
}

#[test]
fn compat_and_corrected_agree_on_synthetic_data() {
    // AM53 to AM55 cannot change these small cases
    let mut s = Scenario::new(4);
    s.create_baseline_readings();
    s.add_in_contamination_data(0.25);
    let summary = |a: &TumorOnlyPurityAnalysis| (a.cutoff(), a.contamination(), a.copy_number_peaks().len(), a.contamination_peaks().len());
    let java = crate::compat::with_java(true, || summary(&s.analysis()));
    let fixed = crate::compat::with_java(false, || summary(&s.analysis()));
    assert_eq!(java, fixed);
}

// PeakGnomadFrequenciesCheckerTest
fn checker_points(points: &[(i32, i32, i32)]) -> (Vec<PositionEvidence>, CandidatePeak) {
    let evidence: Vec<PositionEvidence> = points
        .iter()
        .map(|&(position, depth, alt)| {
            let mut e = PositionEvidence::new("chr1", position, AmberBase::A, AmberBase::T);
            e.read_depth = depth;
            e.ref_support = depth - alt;
            e.alt_support = alt;
            e
        })
        .collect();
    let mut peak = CandidatePeak::new(0.1);
    for (i, e) in evidence.iter().enumerate() {
        peak.test(i as u32, e);
    }
    (evidence, peak)
}

fn check(points: &[(i32, i32, i32)], supplier: impl Fn(&str, i32) -> f64, expected: f64) -> bool {
    let (evidence, peak) = checker_points(points);
    let refs: Vec<&PositionEvidence> = evidence.iter().collect();
    let keys = EvidenceKeys::new(&refs);
    let captured = peak.all_captured_points(&keys);
    gnomad::check_gnomad_frequencies(&supplier, captured.iter().map(|&i| refs[i as usize]), expected).pass
}

const LEVEL1: [(i32, i32, i32); 4] = [(1000, 100, 90), (2000, 100, 91), (3000, 100, 10), (4000, 100, 9)];

#[test]
fn high_allele_depth_variants_are_skewed() {
    assert!(!check(&LEVEL1, |_, p| if p < 2500 { 0.9 } else { 0.5 }, 0.45));
}

#[test]
fn low_allele_depth_variants_are_skewed() {
    assert!(!check(&LEVEL1, |_, p| if p > 2500 { 0.9 } else { 0.5 }, 0.45));
}

#[test]
fn low_and_high_variants_are_skewed() {
    assert!(!check(&LEVEL1, |_, _| 0.9, 0.45));
}

#[test]
fn variant_means_match_expected() {
    assert!(check(&LEVEL1, |_, _| 0.5, 0.45));
}

#[test]
fn check_tolerance() {
    assert!(check(&LEVEL1, |_, _| 0.5, 0.36));
    assert!(check(&LEVEL1, |_, _| 0.5, 0.64));
}

#[test]
fn include_het_and_hom_captures() {
    let points = [(1000, 1000, 900), (2000, 1000, 50)];
    let supplier = |_: &str, p: i32| if p < 1500 { 0.95 } else { 0.35 };
    assert!(!check(&points, supplier, 0.95));
    assert!(!check(&points, supplier, 0.35));
}

#[test]
fn handle_no_high_allele_depth_variants() {
    assert!(!check(&[(2000, 1000, 50)], |_, _| 0.5, 0.5));
}

#[test]
fn handle_no_low_allele_depth_variants() {
    assert!(!check(&[(2000, 1000, 950)], |_, _| 0.5, 0.5));
}

#[test]
fn position_evidence_hash_under_constant_identity_hashes() {
    // Objects.hash("chr1", 1000, <1>, <1>) = 31 * (31 * (31 * (31 + "chr1".hashCode()) + 1000) + 1) + 1
    let e = PositionEvidence::new("chr1", 1000, AmberBase::A, AmberBase::T);
    let s = hmf_common::java_hash::java_string_hash("chr1");
    let expected = 31i32.wrapping_mul(31i32.wrapping_mul(31i32.wrapping_mul(31i32.wrapping_add(s)).wrapping_add(1000)).wrapping_add(1)).wrapping_add(1);
    assert_eq!(position_evidence_hash(&e), expected);
}
