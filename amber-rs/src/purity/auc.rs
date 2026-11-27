//! AucCalculator, CategoryEvidence, CategoryEvidenceIntegral, CanonicalSnvType and the chromosome arm
//! (hmf-common ChrArm / ChrArmLocator.defaultLocator) of the tumor-only purity analysis.
//!
//! An AUC compares how a peak's captured points spread over categories (chromosome arms, or the six
//! canonical SNV types) with how the baseline (heterozygous) points do: the categories are ordered by their
//! ratio of peak points to baseline points, and the area under the cumulative curve is divided by the
//! area of the triangle; 1 means the peak's points are spread as the baseline's, as contamination would be.
//!
//! Finding AM53 (REPORT.md): Java multiplies the two counts of a category as ints, which overflows on
//! whole-genome data.

use crate::types::{AmberBase, PositionEvidence};
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// CanonicalSnvType, in its ordinal order
#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CanonicalSnvType {
    C_A,
    C_G,
    C_T,
    T_A,
    T_C,
    T_G,
}

impl CanonicalSnvType {
    /// CanonicalSnvType.type: the SNV as a change from C or T. Java throws IllegalArgumentException for
    /// equal bases or N; this panics with the same exception text (the standard loci file has neither).
    pub fn of(ref_: AmberBase, alt: AmberBase) -> CanonicalSnvType {
        use AmberBase::*;
        if ref_ == alt || ref_ == N || alt == N {
            panic!("java.lang.IllegalArgumentException");
        }
        let (r, a) = if ref_ == A || ref_ == G { (ref_.complement(), alt.complement()) } else { (ref_, alt) };
        match (r, a) {
            (C, A) => CanonicalSnvType::C_A,
            (C, G) => CanonicalSnvType::C_G,
            (C, T) => CanonicalSnvType::C_T,
            (T, A) => CanonicalSnvType::T_A,
            (T, C) => CanonicalSnvType::T_C,
            (T, G) => CanonicalSnvType::T_G,
            // unreachable: r is C or T and differs from a, neither N
            _ => panic!("java.lang.NullPointerException"),
        }
    }
}

/// TumorOnlyPurityAnalysis.snvTypeClassifier: the type of the change to the minor allele, taking the alt
/// as minor up to a VAF of 0.65 (HET_VAF_UPPER_BOUND)
pub fn snv_type(evidence: &PositionEvidence) -> CanonicalSnvType {
    if evidence.vaf() <= super::HET_VAF_UPPER_BOUND {
        CanonicalSnvType::of(evidence.ref_, evidence.alt)
    } else {
        CanonicalSnvType::of(evidence.alt, evidence.ref_)
    }
}

/// ChrArm(HumanChromosome, Arm): ordered by chromosome ordinal, then P before Q
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ChrArm {
    pub chromosome: usize,
    /// false: P, true: Q
    pub q: bool,
}

/// ChrArmLocator.defaultLocator(version).map(chromosome, position): P below the centromere position
pub fn chr_arm(v38: bool, chromosome: &str, position: i32) -> ChrArm {
    let centre = super::genome::centromere(v38, chromosome);
    let Some(ordinal) = super::genome::chromosome_ordinal(chromosome) else {
        panic!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._{chromosome}");
    };
    ChrArm { chromosome: ordinal, q: position >= centre }
}

/// CategoryEvidence: of the points of one category, how many are evidence (here: in the peak)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CategoryEvidence<T> {
    pub category: T,
    total_points: i32,
    evidence_points: i32,
}

impl<T: Ord + Copy> CategoryEvidence<T> {
    pub fn new(category: T) -> Self {
        CategoryEvidence { category, total_points: 0, evidence_points: 0 }
    }
    pub fn register(&mut self, is_evidence: bool) {
        self.total_points += 1;
        if is_evidence {
            self.evidence_points += 1;
        }
    }
    pub fn set(&mut self, total_points: i32, evidence_points: i32) {
        self.total_points = total_points;
        self.evidence_points = evidence_points;
    }
    pub fn total_points(&self) -> i32 {
        self.total_points
    }
    pub fn evidence_points(&self) -> i32 {
        self.evidence_points
    }
    /// Double.MAX_VALUE when there are no points
    pub fn ratio(&self) -> f64 {
        if self.total_points == 0 {
            return f64::MAX;
        }
        self.evidence_points as f64 / self.total_points as f64
    }
    /// compareTo: by ratio (Double.compare), then by category
    pub fn compare(&self, other: &Self) -> Ordering {
        self.ratio().total_cmp(&other.ratio()).then(self.category.cmp(&other.category))
    }
}

/// CategoryEvidenceIntegral, over the categories sorted by CategoryEvidence.compareTo
pub struct CategoryEvidenceIntegral<T> {
    evidence: Vec<CategoryEvidence<T>>,
}

impl<T: Ord + Copy> CategoryEvidenceIntegral<T> {
    /// the categories are distinct (Java takes a Set), so the order is total and any sort gives it
    pub fn new(mut evidence: Vec<CategoryEvidence<T>>) -> Self {
        evidence.sort_by(|a, b| a.compare(b));
        CategoryEvidenceIntegral { evidence }
    }
    pub fn total_hits(&self) -> i32 {
        self.evidence.iter().fold(0i32, |s, e| s.wrapping_add(e.evidence_points))
    }
    pub fn total_points(&self) -> i32 {
        self.evidence.iter().fold(0i32, |s, e| s.wrapping_add(e.total_points))
    }
    /// the area under the cumulative curve: per category a rectangle under the height so far and a triangle
    pub fn value(&self) -> f64 {
        let java = crate::compat::java("AM53");
        let mut height_so_far = 0.0f64;
        let mut area = 0.0f64;
        for e in &self.evidence {
            area += e.total_points as f64 * height_so_far;
            // Java: evidencePoints() * totalPoints() / 2.0, an int product that overflows past 2^31
            let triangle = if java { e.evidence_points.wrapping_mul(e.total_points) as f64 / 2.0 } else { e.evidence_points as f64 * e.total_points as f64 / 2.0 };
            area += triangle;
            height_so_far += e.evidence_points as f64;
        }
        area
    }
    /// whether value() overflows Java's int product (AM53)
    pub fn int_product_overflows(&self) -> bool {
        self.evidence.iter().any(|e| e.evidence_points.checked_mul(e.total_points).is_none())
    }
}

/// AucCalculator, with the classifier applied by the caller: the baseline's category counts
pub struct AucCalculator<T> {
    baseline: BTreeMap<T, i32>,
}

impl<T: Ord + Copy> AucCalculator<T> {
    pub fn new(baseline: impl IntoIterator<Item = T>) -> Self {
        AucCalculator { baseline: counts(baseline) }
    }

    /// calculateAuc: NaN when the samples or the baseline are empty (0 / 0)
    pub fn calculate_auc(&self, samples: impl IntoIterator<Item = T>) -> f64 {
        self.integral(samples).0
    }

    /// the AUC and whether Java's int product overflowed in computing it (AM53)
    pub fn integral(&self, samples: impl IntoIterator<Item = T>) -> (f64, bool) {
        let counts = counts(samples);
        let mut categories: Vec<T> = counts.keys().copied().collect();
        categories.extend(self.baseline.keys().copied().filter(|c| !counts.contains_key(c)));
        let evidence: Vec<CategoryEvidence<T>> = categories
            .into_iter()
            .map(|c| {
                let mut e = CategoryEvidence::new(c);
                e.set(self.baseline.get(&c).copied().unwrap_or(0), counts.get(&c).copied().unwrap_or(0));
                e
            })
            .collect();
        let integral = CategoryEvidenceIntegral::new(evidence);
        let total_hits = integral.total_hits();
        let total_points = integral.total_points();
        let possible_max = 0.5 * total_hits as f64 * total_points as f64;
        (integral.value() / possible_max, integral.int_product_overflows())
    }
}

fn counts<T: Ord + Copy>(samples: impl IntoIterator<Item = T>) -> BTreeMap<T, i32> {
    let mut counts = BTreeMap::new();
    for s in samples {
        *counts.entry(s).or_insert(0) += 1;
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compat::with_java;
    use AmberBase::*;
    use CanonicalSnvType::*;

    fn arm(chromosome: usize, q: bool) -> ChrArm {
        ChrArm { chromosome, q }
    }

    // CanonicalSnvTypeTest
    #[test]
    fn type_test() {
        assert_eq!(C_A, CanonicalSnvType::of(C, A));
        assert_eq!(C_A, CanonicalSnvType::of(G, T));
        assert_eq!(C_G, CanonicalSnvType::of(C, G));
        assert_eq!(C_G, CanonicalSnvType::of(G, C));
        assert_eq!(C_T, CanonicalSnvType::of(C, T));
        assert_eq!(C_T, CanonicalSnvType::of(G, A));
        assert_eq!(T_A, CanonicalSnvType::of(T, A));
        assert_eq!(T_A, CanonicalSnvType::of(A, T));
        assert_eq!(T_C, CanonicalSnvType::of(T, C));
        assert_eq!(T_C, CanonicalSnvType::of(A, G));
        assert_eq!(T_G, CanonicalSnvType::of(T, G));
        assert_eq!(T_G, CanonicalSnvType::of(A, C));
    }

    #[test]
    #[should_panic]
    fn type_of_n_panics() {
        CanonicalSnvType::of(N, A);
    }

    // TumorOnlyPurityAnalysisTest.mutationTypeClassifierTest
    #[test]
    fn mutation_type_classifier() {
        let pe = |r: AmberBase, a: AmberBase, vaf: f64| {
            let mut e = PositionEvidence::new("1", 1000, r, a);
            e.read_depth = 1000;
            e.alt_support = (1000.0 * vaf) as i32;
            e.ref_support = 1000 - e.alt_support;
            e
        };
        assert_eq!(C_A, snv_type(&pe(C, A, 0.5)));
        assert_eq!(C_A, snv_type(&pe(C, A, 0.65)));
        assert_eq!(T_G, snv_type(&pe(C, A, 0.66)));
    }

    fn with_values(category: ChrArm, hits: i32, total: i32) -> CategoryEvidence<ChrArm> {
        let mut e = CategoryEvidence::new(category);
        for _ in 0..hits {
            e.register(true);
        }
        for _ in hits..total {
            e.register(false);
        }
        e
    }

    // CategoryEvidenceTest
    #[test]
    fn category_evidence() {
        let mut e = CategoryEvidence::new(arm(0, false));
        assert_eq!(e.ratio(), f64::MAX);
        e.register(false);
        assert_eq!(e.ratio(), 0.0);
        e.register(true);
        assert_eq!(e.ratio(), 0.5);
        let mut e = CategoryEvidence::new(arm(0, false));
        e.set(12, 6);
        assert_eq!((e.total_points(), e.evidence_points()), (12, 6));

        let ae1p = with_values(arm(0, false), 1, 10);
        let ae2q = with_values(arm(1, true), 10, 100);
        let ae1q = with_values(arm(0, true), 10, 20);
        let ae2p = with_values(arm(1, false), 20, 20);
        let ae3p = with_values(arm(2, false), 5, 200);
        let ae3q = with_values(arm(2, true), 5, 200);
        assert_eq!(ae1p.compare(&ae2p), Ordering::Less);
        assert_eq!(ae1p.compare(&ae2q), Ordering::Less);
        assert_eq!(ae1p.compare(&ae1q), Ordering::Less);
        assert_eq!(ae1q.compare(&ae2p), Ordering::Less);
        assert_eq!(ae3p.compare(&ae3q), Ordering::Less);
    }

    // CategoryEvidenceIntegralTest
    #[test]
    fn integral() {
        let three = || vec![with_values(arm(0, false), 10, 20), with_values(arm(1, false), 5, 10), with_values(arm(2, false), 2, 20)];
        assert_eq!(CategoryEvidenceIntegral::new(three()).total_hits(), 17);
        assert_eq!(CategoryEvidenceIntegral::new(three()).total_points(), 50);
        assert_eq!(CategoryEvidenceIntegral::new(vec![with_values(arm(0, false), 10, 20)]).value(), 100.0);
        assert_eq!(CategoryEvidenceIntegral::new(three()).value(), 305.0);
        let five = vec![
            with_values(arm(0, false), 2, 10),
            with_values(arm(0, true), 3, 10),
            with_values(arm(1, false), 1, 8),
            with_values(arm(1, true), 2, 12),
            with_values(arm(2, false), 1, 10),
        ];
        assert_eq!(CategoryEvidenceIntegral::new(five).value(), 178.0);
    }

    #[test]
    fn integral_int_overflow() {
        // AM53: 50,000 x 50,000 overflows an int; Java's product wraps to -1,794,967,296
        let big = || vec![with_values_set(arm(0, false), 50_000, 50_000), with_values_set(arm(1, false), 0, 100)];
        let java = with_java(true, || CategoryEvidenceIntegral::new(big()).value());
        let fixed = with_java(false, || CategoryEvidenceIntegral::new(big()).value());
        assert_eq!(fixed, 0.0 * 100.0 + 1.25e9);
        assert_eq!(java, (50_000i32.wrapping_mul(50_000)) as f64 / 2.0);
        assert!(CategoryEvidenceIntegral::new(big()).int_product_overflows());
    }

    fn with_values_set(category: ChrArm, hits: i32, total: i32) -> CategoryEvidence<ChrArm> {
        let mut e = CategoryEvidence::new(category);
        e.set(total, hits);
        e
    }

    // AucCalculatorTest
    const CENTRE: i32 = 10_000_000;

    fn arm_of(e: &PositionEvidence) -> ChrArm {
        arm(super::super::genome::chromosome_ordinal(&e.chromosome).unwrap(), e.position >= CENTRE)
    }

    fn pe(chromosome: usize, position: i32, ref_support: i32, alt_support: i32, depth: i32) -> PositionEvidence {
        let chr = format!("chr{}", hmf_common::regions::HUMAN_CHROMOSOMES[chromosome]);
        let mut e = PositionEvidence::new(&chr, position, A, C);
        e.read_depth = depth;
        e.ref_support = ref_support;
        e.alt_support = alt_support;
        e
    }

    #[test]
    fn uneven_distribution_cost() {
        let mut baseline = Vec::new();
        let mut hits = Vec::new();
        for i in 1..=8 {
            baseline.push(pe(0, i * 1000, 5, 5, 10));
        }
        hits.push(pe(0, 9000, 3, 7, 10));
        hits.push(pe(0, 10000, 3, 7, 10));
        for i in 1..=9 {
            baseline.push(pe(1, i * 1000, 4, 6, 10));
        }
        hits.push(pe(1, 10000, 3, 7, 10));
        for i in 1..=4 {
            baseline.push(pe(2, i * 1000, 9, 1, 10));
        }
        hits.push(pe(2, 5000, 3, 7, 10));
        for i in 1..=10 {
            baseline.push(pe(3, i * 1000, 2, 8, 10));
        }
        for i in 11..=15 {
            hits.push(pe(3, i * 1000, 3, 7, 10));
        }
        let p_arm = |e: &PositionEvidence| arm(super::super::genome::chromosome_ordinal(&e.chromosome).unwrap(), false);
        let calculator = AucCalculator::new(baseline.iter().map(p_arm));
        assert!((calculator.calculate_auc(hits.iter().map(p_arm)) - 99.5 / 139.5).abs() < 1e-4);
    }

    /// createEvenlySpacedDataPointsWithContamination
    fn spaced(chromosome: usize, start: i32, blocks: i32, contaminated_per_block: i32, per_block: i32, alt: i32, total: i32) -> (Vec<PositionEvidence>, Vec<PositionEvidence>) {
        let (mut baseline, mut contamination) = (Vec::new(), Vec::new());
        let mut position = start;
        for _ in 0..blocks {
            for point in 0..per_block {
                if point < contaminated_per_block {
                    contamination.push(pe(chromosome, position, total - alt, alt, total));
                } else {
                    baseline.push(pe(chromosome, position, total, 0, total));
                }
                position += 1000;
            }
        }
        (baseline, contamination)
    }

    fn auc_of(baseline: &[PositionEvidence], contamination: &[PositionEvidence]) -> f64 {
        AucCalculator::new(baseline.iter().map(arm_of)).calculate_auc(contamination.iter().map(arm_of))
    }

    #[test]
    fn same_contamination_on_all_chromosome_arms() {
        let (mut b, mut c) = (Vec::new(), Vec::new());
        for chromosome in 0..24 {
            let (pb, pc) = spaced(chromosome, 1000, 50, 10, 100, 20, 100);
            b.extend(pb);
            c.extend(pc);
            let (qb, qc) = spaced(chromosome, CENTRE + 1000, 30, 10, 100, 20, 100);
            b.extend(qb);
            c.extend(qc);
        }
        assert!((auc_of(&b, &c) - 1.0).abs() < 1e-4);
    }

    #[test]
    fn all_contamination_in_1p() {
        let (mut b, mut c) = (Vec::new(), Vec::new());
        for chromosome in 0..24 {
            if chromosome == 0 {
                let (pb, pc) = spaced(chromosome, 1000, 50, 10, 100, 20, 100);
                b.extend(pb);
                c.extend(pc);
                b.extend(spaced(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100).0);
            } else {
                b.extend(spaced(chromosome, 1000, 50, 0, 100, 20, 100).0);
                b.extend(spaced(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100).0);
            }
        }
        assert!((auc_of(&b, &c) - 4500.0 / 191500.0).abs() < 1e-4);
    }

    #[test]
    fn high_contamination_in_p_arms_only() {
        let (mut b, mut c) = (Vec::new(), Vec::new());
        for chromosome in 0..24 {
            let (pb, pc) = spaced(chromosome, 1000, 10, 20, 100, 10, 50);
            b.extend(pb);
            c.extend(pc);
            let (qb, qc) = spaced(chromosome, CENTRE + 1000, 10, 0, 100, 0, 40);
            b.extend(qb);
            c.extend(qc);
        }
        assert!((auc_of(&b, &c) - 4.0 / 9.0).abs() < 1e-4);
    }

    #[test]
    fn no_contamination() {
        let (mut b, mut c) = (Vec::new(), Vec::new());
        for chromosome in 0..24 {
            let (pb, pc) = spaced(chromosome, 1000, 50, 0, 100, 20, 100);
            b.extend(pb);
            c.extend(pc);
            let (qb, qc) = spaced(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100);
            b.extend(qb);
            c.extend(qc);
        }
        assert!(auc_of(&b, &c).is_nan());
    }

    #[test]
    fn default_locator() {
        assert_eq!(chr_arm(true, "chr1", 123605522), arm(0, false));
        assert_eq!(chr_arm(true, "chr1", 123605523), arm(0, true));
        assert_eq!(chr_arm(false, "X", 60132012), arm(22, true));
    }
}
