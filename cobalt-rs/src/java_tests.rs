//! The Java unit tests of the COBALT classes ported here (cobalt/src/test/java), adapted to this crate's
//! types: one Rust test per Java test, named after it. Tests of Java object plumbing (equals, hashCode,
//! toString) and of class wiring with Mockito mocks (BamCalculationTest, TumorCalculationTest,
//! ReferenceCalculationTest, DoNothingNormaliserTest, NoOpConsolidatorTest, ReadDepthAccumulatorTest's
//! accessors) are not ported: the port has no such classes. ReadCountMedianTest tests hmf-common's
//! Doubles.interpolatedMedian, which COBALT does not use.

use crate::calc::*;
use crate::chromosome::{Chr, PerChr};
use crate::count::{ChrCounts, SampleCounts};
use crate::gc::{GcProfile, WindowStatuses};
use crate::targeted::{Scope, TargetRegions};

fn c(n: &str) -> Chr {
    Chr::from_name(n).unwrap()
}

fn close(expected: f64, actual: f64, eps: f64) {
    assert!((expected - actual).abs() <= eps, "expected {expected}, got {actual}");
}

/// CalculationsTestBase.br: a BamRatio of a DepthReading (GC 0 for depth 0)
fn br(pos: i32, depth: f64, gc: f64, on_target: bool) -> BamRatio {
    BamRatio::new(pos, depth, if depth == 0.0 { 0.0 } else { gc }, on_target)
}

// BamRatioTest

#[test]
fn bam_ratio_depth0() {
    let mut r = BamRatio::new(1001, 0.0, 0.0, true);
    assert_eq!((r.read_depth, r.ratio(), r.gc_content, r.diploid_adjusted_ratio()), (0.0, 0.0, 0.0, -1.0));
    for f in [3.4, 45.89, 34.56] {
        r.normalise(f);
        assert_eq!((r.read_depth, r.ratio(), r.gc_content, r.diploid_adjusted_ratio()), (0.0, 0.0, 0.0, -1.0));
    }
}

#[test]
fn bam_ratio_constructors() {
    let r = BamRatio::simple(19_001, 123.4, 0.55);
    assert_eq!((r.position, r.ratio(), r.read_depth, r.gc_content, r.diploid_adjusted_ratio()), (19_001, 123.4, 123.4, 0.55, -1.0));
}

#[test]
fn bam_ratio_override() {
    let mut r = BamRatio::new(1001, 82.0, 0.49, true);
    assert_eq!((r.read_depth, r.ratio(), r.gc_content), (82.0, 82.0, 0.49));
    r.override_ratio(1.8);
    assert_eq!(r.ratio(), 1.8);
    // overrideRatioSetsIncludedStatus
    let mut r = BamRatio::new(1001, 82.0, 0.49, false);
    assert_eq!(r.ratio(), -1.0);
    r.override_ratio(1.8);
    assert_eq!(r.ratio(), 1.8);
    // overrideWithNegativeRatioSetsIncludedStatus
    let mut r = BamRatio::new(1001, 82.0, 0.49, false);
    r.override_ratio(-1.0);
    assert_eq!(r.ratio(), -1.0);
}

#[test]
fn bam_ratio_normalise_by_mean() {
    let mut r = BamRatio::new(1001, 82.0, 0.49, true);
    r.normalise(0.5);
    assert_eq!((r.read_depth, r.ratio(), r.gc_content), (82.0, 164.0, 0.49));
    // normalise0ByMean
    let mut r = BamRatio::new(1001, 0.0, 0.0, true);
    r.normalise(0.5);
    assert_eq!((r.read_depth, r.ratio(), r.gc_content), (0.0, 0.0, 0.0));
    // normaliseNaNByMean, handleNaNDepthInConstructor
    let mut r = BamRatio::new(1001, f64::NAN, 0.49, true);
    assert_eq!((r.ratio(), r.read_depth, r.gc_content), (-1.0, -1.0, 0.49));
    r.normalise(0.5);
    assert_eq!((r.position, r.ratio(), r.read_depth, r.gc_content), (1001, -1.0, -1.0, 0.49));
    // handleNaNGcInConstructor
    let r = BamRatio::new(1001, 34.0, f64::NAN, true);
    assert_eq!((r.ratio(), r.read_depth, r.gc_content, r.diploid_adjusted_ratio()), (34.0, 34.0, -1.0, -1.0));
}

#[test]
fn bam_ratio_diploid_adjusted() {
    // normaliseDiploidRatioWhenNotSet
    let mut r = BamRatio::new(1001, 82.0, 0.49, true);
    r.normalise_diploid_adjusted_ratio(0.5);
    assert_eq!((r.ratio(), r.diploid_adjusted_ratio()), (82.0, -1.0));
    // normaliseDiploidRatio
    let mut r = BamRatio::new(1001, 82.0, 0.49, true);
    r.set_diploid_adjusted_ratio(1.8);
    assert_eq!(r.diploid_adjusted_ratio(), 1.8);
    r.normalise_diploid_adjusted_ratio(2.0);
    assert_eq!((r.read_depth, r.diploid_adjusted_ratio()), (82.0, 0.9));
    // normaliseDiploidRatioByNaN, normaliseDiploidRatioByZero
    for f in [f64::NAN, 0.0] {
        let mut r = BamRatio::new(1001, 82.0, 0.49, true);
        r.set_diploid_adjusted_ratio(1.8);
        r.normalise_diploid_adjusted_ratio(f);
        assert_eq!(r.diploid_adjusted_ratio(), -1.0);
    }
    // diploidRatioWhenReadRatioIsZero
    let mut r = BamRatio::simple(1001, 0.0, 0.0);
    assert_eq!(r.diploid_adjusted_ratio(), -1.0);
    r.set_diploid_adjusted_ratio(1.8);
    assert_eq!(r.diploid_adjusted_ratio(), 0.0);
    r.normalise_diploid_adjusted_ratio(10.0);
    assert_eq!(r.diploid_adjusted_ratio(), 0.0);
}

#[test]
fn bam_ratio_target_region_and_gc() {
    // inTargetRegion, outsideTargetRegion
    assert_eq!(BamRatio::new(1001, 82.0, 0.49, true).ratio(), 82.0);
    assert_eq!(BamRatio::new(1001, 82.0, 0.49, false).ratio(), -1.0);
    // gcNormalise
    let mut r = BamRatio::new(1001, 82.0, 0.49, true);
    r.normalise(100.0);
    close(0.82, r.ratio(), 1e-12);
    // gcNormaliseByZero, gcNormaliseNegativeValue
    for f in [0.0, -1.0] {
        let mut r = BamRatio::new(1001, 82.0, 0.49, true);
        r.normalise(f);
        assert_eq!((r.ratio(), r.read_depth, r.gc_content), (-1.0, 82.0, 0.49));
    }
    // applyEnrichment
    let mut r = BamRatio::new(1001, 82.0, 0.49, true);
    r.normalise(0.5);
    assert_eq!(r.ratio(), 164.0);
    // handleNaNEnrichment, handleNaNEnrichmentForDepth0Window
    let mut r = BamRatio::new(1001, 82.0, 0.49, true);
    r.normalise(f64::NAN);
    assert_eq!(r.ratio(), -1.0);
    let mut r = BamRatio::simple(1001, 0.0, 0.0);
    r.normalise(f64::NAN);
    assert_eq!(r.ratio(), -1.0);
}

/// CO4: a ratio of 0 under a factor that masks (Java keeps the 0)
#[test]
fn zero_ratio_under_a_masking_factor() {
    for (java, expected) in [(true, 0.0), (false, -1.0)] {
        crate::compat::with_java(java, || {
            let mut r = BamRatio::new(1001, 0.0, 0.0, true);
            r.normalise(-1.0);
            assert_eq!(r.ratio(), expected);
            let mut r = BamRatio::new(1001, 0.0, 0.0, true);
            r.normalise(0.0);
            assert_eq!(r.ratio(), expected);
            // a factor that applies keeps the 0 in both modes
            let mut r = BamRatio::new(1001, 0.0, 0.0, true);
            r.normalise(2.0);
            assert_eq!(r.ratio(), 0.0);
        });
    }
}

// CobaltWindowTest

#[test]
fn cobalt_window_construct_and_include() {
    let w = CobaltWindow::new(1001, 82.0, 0.49, true, true);
    assert!(w.excluded && w.in_target && w.bucket.is_none());
    assert!(!CobaltWindow::new(5001, 1892.0, 0.52, true, true).include());
    assert!(!CobaltWindow::new(5001, 1892.0, 0.52, true, false).include());
    assert!(CobaltWindow::new(5001, 1892.0, 0.52, false, true).include());
    assert!(!CobaltWindow::new(5001, 1892.0, 0.52, false, false).include());
}

#[test]
fn cobalt_window_to_bam_ratio() {
    let r = CobaltWindow::new(5001, 1892.0, 0.52, true, true).to_bam_ratio();
    assert_eq!((r.read_depth, r.ratio(), r.gc_content), (1892.0, -1.0, 0.52));
    let r = CobaltWindow::new(5001, 1892.0, 0.52, true, false).to_bam_ratio();
    assert_eq!(r.ratio(), -1.0);
    let r = CobaltWindow::new(5001, 1892.0, 0.52, false, true).to_bam_ratio();
    assert_eq!(r.ratio(), 1892.0);
    let r = CobaltWindow::new(5001, 1892.0, 0.52, false, false).to_bam_ratio();
    assert_eq!(r.ratio(), -1.0);
}

#[test]
fn cobalt_window_corrected_by_reference_value() {
    let mut asked = 0;
    let w = CobaltWindow::new(1001, 82.0, 0.49, false, true);
    let corrected = w.clone().corrected_by_reference_value(|| {
        asked += 1;
        Some(0.89)
    });
    assert_eq!((corrected, asked), (w, 0));
    let w = CobaltWindow::new(1001, -1.0, 0.49, false, true);
    assert_eq!(w.clone().corrected_by_reference_value(|| Some(0.89)), w);
    let w = CobaltWindow::new(1001, 0.99, 0.49, false, true);
    assert_eq!(w.corrected_by_reference_value(|| Some(0.89)).gc_content, 0.89);
}

#[test]
fn cobalt_window_bucketed() {
    let mut pails: [Vec<f64>; 101] = std::array::from_fn(|_| Vec::new());
    // bucketedTest
    let b = CobaltWindow::new(5001, 1892.0, 0.52, false, true).bucketed(c("2"), &mut pails).unwrap();
    assert_eq!((b.bucket, b.excluded, b.in_target), (Some(52), false, true));
    assert_eq!(pails[52], vec![1892.0]);
    // bucketedAllosomeTest
    let b = CobaltWindow::new(5001, 1892.0, 0.52, false, true).bucketed(c("Y"), &mut pails).unwrap();
    assert_eq!(b.bucket, Some(52));
    assert_eq!(pails[52].len(), 1);
    // bucketedOffTargetRegionTest
    let b = CobaltWindow::new(5001, 1892.0, 0.52, false, false).bucketed(c("2"), &mut pails).unwrap();
    assert_eq!((b.bucket, b.in_target), (Some(52), false));
    assert_eq!(pails[52].len(), 1);
    // bucketedDepth0Test
    pails[0].push(99.0);
    let b = CobaltWindow::new(5001, 0.0, 0.0, false, true).bucketed(c("2"), &mut pails).unwrap();
    assert_eq!(b.bucket, Some(0));
    assert_eq!(pails[0], vec![99.0]);
    // an excluded window has no bucket
    assert_eq!(CobaltWindow::new(5001, 10.0, 0.52, true, true).bucketed(c("2"), &mut pails).unwrap().bucket, None);
}

// GCPailTest, GCPailsListTest

#[test]
fn gc_pail_median() {
    use crate::stats::median;
    let mut v: Vec<f64> = vec![];
    assert_eq!(if v.is_empty() { 0.0 } else { median(&mut v) }, 0.0);
    v.push(34.5);
    assert_eq!(median(&mut v.clone()), 34.5);
    v.push(34.5);
    v.push(38.5);
    assert_eq!(median(&mut v.clone()), 34.5);
    v.push(38.5);
    assert_eq!(median(&mut v.clone()), (34.5 + 38.5) / 2.0);
}

#[test]
fn gc_pails_list() {
    for (gc, bucket) in [(0.0001, 0), (0.0002, 0), (0.003, 0), (0.2001, 20), (0.21, 21), (0.56999, 57), (0.99499, 99), (0.999, 100), (1.0, 100)] {
        assert_eq!(gc_bucket_index(gc), bucket, "{gc}");
    }
}

// GcBucketStatisticsTest

fn bucket_stats(readings: &[(f64, f64)], min: i32, max: i32) -> GcBucketStatistics {
    let mut pails: [Vec<f64>; 101] = std::array::from_fn(|_| Vec::new());
    for &(gc, depth) in readings {
        pails[gc_bucket_index(gc) as usize].push(depth);
    }
    GcBucketStatistics::new(&mut pails, min, max).unwrap()
}

fn check_buckets(s: &GcBucketStatistics, expected: &[(usize, f64)]) {
    for &(bucket, e) in expected {
        close(e, s.median_read_depth(Some(bucket)), 0.001);
    }
}

#[test]
fn gc_bucket_statistics_no_data() {
    let s = bucket_stats(&[], 20, 70);
    for i in 0..100 {
        close(if !(21..=69).contains(&i) { -1.0 } else { 0.0 }, s.median_read_depth(Some(i)), 0.001);
    }
}

#[test]
fn gc_bucket_statistics_one_bucket() {
    let s = bucket_stats(&[(0.5, 10.0)], 23, 67);
    check_buckets(&s, &[(48, 0.0), (49, 3.3333), (50, 3.3333), (51, 3.3333), (52, 0.0)]);
    let s = bucket_stats(&[(0.5, 10.0), (0.5, 10.0), (0.5, 12.0)], 23, 67);
    check_buckets(&s, &[(48, 0.0), (49, 3.3333), (50, 3.3333), (51, 3.3333), (52, 0.0)]);
}

#[test]
fn gc_bucket_statistics_two_buckets() {
    let s = bucket_stats(&[(0.50, 8.0), (0.51, 9.0)], 23, 67);
    check_buckets(&s, &[(48, 0.0), (49, 8.0 / 3.0), (50, 17.0 / 3.0), (51, 17.0 / 3.0), (52, 3.0), (53, 0.0)]);
    let s = bucket_stats(&[(0.50, 8.0), (0.52, 10.0)], 23, 67);
    check_buckets(&s, &[(48, 0.0), (49, 8.0 / 3.0), (50, 8.0 / 3.0), (51, 6.0), (52, 10.0 / 3.0), (53, 10.0 / 3.0), (54, 0.0)]);
}

#[test]
fn gc_bucket_statistics_many_buckets() {
    let readings = [
        (0.50, 8.0),
        (0.50, 8.0),
        (0.50, 10.0),
        (0.51, 12.0),
        (0.51, 8.0),
        (0.52, 12.0),
        (0.52, 12.0),
        (0.53, 15.0),
        (0.53, 14.0),
        (0.53, 13.0),
        (0.54, 16.0),
        (0.55, 18.0),
    ];
    let s = bucket_stats(&readings, 23, 67);
    check_buckets(
        &s,
        &[(48, 0.0), (49, 8.0 / 3.0), (50, 6.0), (51, 10.0), (52, 12.0), (53, 14.0), (54, 16.0), (55, 34.0 / 3.0), (56, 6.0), (57, 0.0)],
    );
}

#[test]
fn gc_bucket_statistics_all_buckets() {
    let readings: Vec<(f64, f64)> = (0..100).map(|i| (0.01 * i as f64, i as f64)).collect();
    let s = bucket_stats(&readings, 19, 80);
    check_buckets(&s, &[(0, -1.0), (10, -1.0), (19, -1.0), (20, 20.0), (50, 50.0), (79, 79.0), (80, -1.0), (81, -1.0), (99, -1.0), (100, -1.0)]);
    assert_eq!(s.median_read_depth(None), -1.0);
}

#[test]
fn gc_bucket_statistics_allowed() {
    // isAllowedTest: 25-68 allowed with limits 24 and 68 (68 has no mean: bucket 69 is not allowed)
    let s = bucket_stats(&[], 24, 68);
    for i in 0..100 {
        let allowed = s.median_read_depth(Some(i)) != -1.0 || i == 68;
        assert_eq!(allowed, (25..=68).contains(&i), "{i}");
    }
}

// ResultsConsolidatorTest, LowCoverageConsolidatorTest

#[test]
fn consolidation_count() {
    for (depth, count) in [(16.0, 1), (8.1, 1), (8.0, 10), (3.0, 30), (0.8, 100), (0.3, 300), (0.15, 500), (0.08, 1000), (0.008, 1000), (0.00015, 1000)] {
        assert_eq!(calc_consolidation_count(depth), count, "{depth}");
    }
}

#[test]
fn consolidate_boundaries() {
    let mut positions = vec![1_001, 2_001, 5_001, 6_001, 10_001, 3_020_001, 3_021_001, 3_024_001, 3_029_001];
    let b = consolidate_into_buckets(&positions, 4);
    let se: Vec<(i32, i32)> = b.iter().map(|b| (b.start, b.end)).collect();
    assert_eq!(se, vec![(1_001, 8_001), (9_001, 11_001), (3_020_001, 3_030_001)]);
    positions.push(3_034_001);
    let b = consolidate_into_buckets(&positions, 4);
    let se: Vec<(i32, i32)> = b.iter().map(|b| (b.start, b.end)).collect();
    assert_eq!(se, vec![(1_001, 8_001), (9_001, 11_001), (3_020_001, 3_031_001), (3_032_001, 3_035_001)]);
}

fn chr_ratios(chr: &str, list: Vec<BamRatio>) -> ChrRatios {
    ChrRatios { chr: c(chr), ratios: list }
}

fn forty_windows() -> Vec<ChrRatios> {
    ["1", "2"]
        .iter()
        .map(|chr| chr_ratios(chr, (0..40).map(|i| br(i * 1000 + 1, 10.0 + i as f64 * 0.01, 0.40 + i as f64 * 0.01, true)).collect()))
        .collect()
}

#[test]
fn ratio_and_gc_are_averaged_in_consolidation() {
    let ratios = forty_windows();
    let mut consolidator = Consolidator::LowCoverage { count: calc_consolidation_count(8.0), boundaries: Vec::new() };
    let consolidated = consolidator.consolidated_for_test(&ratios);
    assert_eq!(consolidated.len(), 2);
    for (_, list) in &consolidated {
        let got: Vec<(i32, f64, f64)> = list.iter().map(|r| (r.position, r.ratio(), r.gc_content)).collect();
        let expected = [(4001, 10.045, 0.445), (14001, 10.145, 0.545), (24001, 10.245, 0.645), (35001, 10.345, 0.745)];
        assert_eq!(got.len(), 4);
        for (g, e) in got.iter().zip(expected) {
            assert_eq!(g.0, e.0);
            close(e.1, g.1, 1e-4);
            close(e.2, g.2, 1e-4);
        }
    }
}

#[test]
fn use_consolidation_count() {
    let ratios = forty_windows();
    let mut consolidator = Consolidator::LowCoverage { count: calc_consolidation_count(4.0), boundaries: Vec::new() };
    let consolidated = consolidator.consolidated_for_test(&ratios);
    for (_, list) in &consolidated {
        let got: Vec<(i32, f64, f64)> = list.iter().map(|r| (r.position, r.ratio(), r.gc_content)).collect();
        assert_eq!(got.len(), 2);
        assert_eq!((got[0].0, got[1].0), (9001, 30001));
        close(10.095, got[0].1, 1e-4);
        close(0.495, got[0].2, 1e-4);
        close(10.295, got[1].1, 1e-4);
        close(0.695, got[1].2, 1e-4);
    }
}

#[test]
fn masked_ratios_are_skipped_in_consolidation() {
    let list: Vec<BamRatio> = (0..200).map(|i| if i % 5 == 0 { br(i * 1000 + 1, 10.0, 0.5, true) } else { br(i * 1000 + 1, -1.0, 0.4, true) }).collect();
    let ratios = vec![chr_ratios("1", list)];
    let mut consolidator = Consolidator::LowCoverage { count: calc_consolidation_count(8.0), boundaries: Vec::new() };
    let consolidated = consolidator.consolidated_for_test(&ratios);
    assert_eq!(consolidated.len(), 1);
    let got: Vec<i32> = consolidated[0].1.iter().map(|r| r.position).collect();
    assert_eq!(got, vec![23001, 72001, 122001, 172001]);
    close(10.0, consolidated[0].1[0].ratio(), 1e-4);
    close(0.5, consolidated[0].1[0].gc_content, 1e-4);
}

#[test]
fn consolidated_regions_are_limited_in_extent() {
    let list: Vec<BamRatio> =
        (0..200).map(|i| if i % 10 == 0 { br(i * 100_000 + 1, 10.0, 0.5, true) } else { br(i * 100_000 + 1, -1.0, 0.4, true) }).collect();
    let ratios = vec![chr_ratios("1", list)];
    let mut consolidator = Consolidator::LowCoverage { count: calc_consolidation_count(8.0), boundaries: Vec::new() };
    let consolidated = consolidator.consolidated_for_test(&ratios);
    let got: Vec<i32> = consolidated[0].1.iter().map(|r| r.position).collect();
    assert_eq!(got, vec![1_000_001, 4_000_001, 7_000_001, 10_000_001, 13_000_001, 16_000_001, 18_500_001]);
}

#[test]
fn boundaries_are_reused() {
    let tumor: Vec<ChrRatios> = [("1", 10.0, 0.5), ("2", 12.0, 0.52)]
        .iter()
        .map(|(chr, d, gc)| chr_ratios(chr, (0..1000).map(|i| if i % 10 == 0 { br(i * 1000 + 1, *d, *gc, true) } else { br(i * 1000 + 1, -1.0, 0.4, true) }).collect()))
        .collect();
    let reference: Vec<ChrRatios> =
        [("1", 11.0, 0.51), ("2", 13.0, 0.53)].iter().map(|(chr, d, gc)| chr_ratios(chr, (0..1000).map(|i| br(i * 1000 + 1, *d, *gc, true)).collect())).collect();
    let mut consolidator = Consolidator::LowCoverage { count: calc_consolidation_count(8.0), boundaries: Vec::new() };
    let t = consolidator.consolidated_for_test(&tumor);
    let r = consolidator.consolidated_for_test(&reference);
    let positions = |x: &Vec<(Chr, Vec<BamRatio>)>| -> Vec<Vec<i32>> { x.iter().map(|(_, l)| l.iter().map(|b| b.position).collect()).collect() };
    assert_eq!(positions(&t), positions(&r));
}

// BamRatiosTest (consolidation applied to the windows)

fn hundred_windows() -> Vec<ChrRatios> {
    let mut seed = 11u64;
    let mut next = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    ["1", "2", "3"]
        .iter()
        .map(|chr| {
            chr_ratios(
                chr,
                (0..100)
                    .map(|i| {
                        let mut r = BamRatio::simple(i * 1000 + 1, next() * 10.0, (next() * 40.0 + 25.0) / 100.1);
                        r.override_ratio(next());
                        r
                    })
                    .collect(),
            )
        })
        .collect()
}

#[test]
fn bam_ratios_consolidation_identical() {
    let mut ratios = hundred_windows();
    let before = ratios.clone();
    let mut consolidator = Consolidator::NoOp;
    crate::calc::apply_consolidation(&mut ratios, &mut consolidator).unwrap();
    for (a, b) in ratios.iter().zip(&before) {
        assert_eq!(a.ratios, b.ratios);
    }
}

#[test]
fn bam_ratios_consolidate_by_10_and_30() {
    for step in [10, 30] {
        let mut ratios = hundred_windows();
        let before = ratios.clone();
        let consolidated: Vec<(Chr, Vec<BamRatio>)> = ratios
            .iter()
            .map(|c| (c.chr, (1..=(100 / step)).filter(|i| i * step < 100).map(|i| BamRatio::simple(i * step * 1000 + 1, 0.1 * i as f64, 0.5)).collect()))
            .collect();
        apply_consolidated(&mut ratios, &consolidated).unwrap();
        for (ci, c) in ratios.iter().enumerate() {
            for i in 1..100usize {
                let r = &c.ratios[i];
                assert_eq!((r.position, r.read_depth, r.gc_content), (before[ci].ratios[i].position, before[ci].ratios[i].read_depth, before[ci].ratios[i].gc_content));
                if i % step as usize == 0 {
                    close(consolidated[ci].1[i / step as usize - 1].ratio(), r.ratio(), 1e-12);
                } else {
                    assert_eq!(r.ratio(), -1.0);
                }
            }
        }
    }
}

#[test]
fn bam_ratios_consolidate_empty() {
    let mut ratios: Vec<ChrRatios> = Vec::new();
    let consolidated = vec![(c("1"), vec![BamRatio::simple(10_001, 1.0, 0.5)])];
    apply_consolidated(&mut ratios, &consolidated).unwrap();
    assert!(ratios.is_empty());
}

// DiploidRatioNormaliserTest

/// Java's expectations hold in Java mode only: from the tenth valid ratio on its rolling window keeps a 0.51
/// that was removed (CO2). The windows of the last two ratios are {0.49, 0.49, 0.58, 0.59} (median 0.535)
/// and {0.49, 0.49, 0.59} (median 0.49).
#[test]
fn diploid_ratio_normaliser() {
    let input = [-1.0, 0.50, 0.51, 0.51, 0.51, 0.49, 0.58, 0.51, 0.58, 0.49, 0.49, -1.0, 0.59];
    let java = crate::compat::with_java(true, || diploid_ratio_normalise(&input, 1.0, 2, 2));
    let expected = [-1.0, 0.5 / 0.51, 1.0, 1.0, 1.0, 0.49 / 0.51, 0.58 / 0.51, 1.0, 0.58 / 0.51, 0.49 / 0.51, 0.49 / 0.51, -1.0, 0.59 / 0.50];
    for (e, o) in expected.iter().zip(&java) {
        close(*e, *o, 1e-10);
    }
    let fixed = crate::compat::with_java(false, || diploid_ratio_normalise(&input, 1.0, 2, 2));
    assert_eq!(&fixed[..10], &java[..10]);
    close(0.49 / 0.535, fixed[10], 1e-10);
    close(0.59 / 0.49, fixed[12], 1e-10);
    let mut values: Vec<f64> = input.iter().copied().filter(|v| *v > 0.0).collect();
    assert_eq!(values.len(), 11);
    close(0.51, crate::stats::median(&mut values), 1e-4);
}

#[test]
fn diploid_ratio_normaliser_close_to_zero() {
    let input = [0.0, 0.0, 0.002, 0.0, 0.0];
    assert_eq!(diploid_ratio_normalise(&input, 1.0, 5, 5), input.to_vec());
}

#[test]
fn diploid_ratio_normaliser_max_window_distance() {
    let out = diploid_ratio_normalise(&[1.0, 1.5, -1.0, 1.1, 1.2], 1.15, 2, 1);
    let expected = [1.15 / 1.1, 1.5 * 1.15 / 1.15, -1.0, 1.1 * 1.15 / 1.15, 1.2 * 1.15 / 1.2];
    for (e, o) in expected.iter().zip(&out) {
        close(*e, *o, 1e-10);
    }
}

#[test]
fn diploid_ratio_normaliser_min_coverage() {
    assert_eq!(diploid_ratio_normalise(&[1.0, 1.5, 2.0, -1.0, -1.0], 1.5, 1, 3), vec![1.0, 1.5, 2.0, -1.0, -1.0]);
}

// DiploidNormaliserTest

/// DiploidNormaliserTest.addRatios: ratios of depth/100, and the median and count of the positive ones
fn add_ratios(chr: &str, depths: &[f64]) -> (ChrRatios, f64, i32) {
    let mut list = Vec::new();
    for (i, &d) in depths.iter().enumerate() {
        let mut r = br(i as i32 * 1000 + 1, d, 0.45, true);
        if d > 0.0 {
            r.normalise(100.0);
        }
        list.push(r);
    }
    let mut positive: Vec<f64> = depths.iter().copied().filter(|d| *d > 0.0).collect();
    let count = positive.len() as i32;
    (chr_ratios(chr, list), crate::constants::doubles::median(&mut positive) / 100.0, count)
}

#[test]
fn diploid_normaliser_median_ratios() {
    let (c1, m1, n1) = add_ratios("1", &[-1.0, 22.0, 23.0, 24.0, -1.0, 34.0, 27.0, 25.0, -1.0, 30.0, 27.0, 27.0, 25.0, 23.0, 21.0, 29.0, 23.0, 25.0, 31.0, 80.0, 32.0]);
    let (c4, m4, n4) = add_ratios("4", &[19.0, 19.0, 20.0, 22.0, 23.0, 19.0, 17.0, 20.0, 27.0, 27.0, 23.0, 33.0, 21.0, 29.0, 23.0, 25.0, 31.0, 23.0, 23.0]);
    let (c2, m2, n2) = add_ratios("2", &[-1.0, 32.0, 33.0, 34.0, 36.0, 27.0, 25.0, 30.0, 27.0, 37.0, 35.0, 33.0, 31.0, 29.0, 23.0, 25.0, 31.0, 80.0, 23.0]);
    let mut ratios = vec![c1, c4, c2];
    let medians = diploid_normaliser_with(&mut ratios, true, 5, 7);
    let got: Vec<(&str, f64, i32)> = medians.iter().map(|m| (m.chromosome.as_str(), m.median_ratio, m.count)).collect();
    assert_eq!(got.len(), 3);
    for (g, e) in got.iter().zip([("chr1", m1, n1), ("chr2", m2, n2), ("chr4", m4, n4)]) {
        assert_eq!((g.0, g.2), (e.0, e.2));
        close(e.1, g.1, 0.001);
    }
}

#[test]
fn diploid_normaliser_normalise() {
    let (c1, ..) = add_ratios("1", &[20.0, 20.0, 20.0, 25.0, 25.0, 25.0]);
    let (c4, ..) = add_ratios("4", &[20.0, 20.0, 20.0, 20.0, 30.0, 30.0, 30.0]);
    let (cx, ..) = add_ratios("X", &[-1.0, 30.0, 30.0, 50.0, 30.0, 30.0, 40.0, 40.0]);
    let (cy, ..) = add_ratios("Y", &[10.0, 10.0, 15.0, -1.0, 10.0, 15.0]);
    let mut ratios = vec![c1, c4, cx, cy];
    let medians = diploid_normaliser_with(&mut ratios, true, 3, 3);
    let d = |ci: usize| -> Vec<f64> { ratios[ci].ratios.iter().map(|r| r.diploid_adjusted_ratio()).collect() };
    let check = |got: Vec<f64>, expected: &[f64]| {
        assert_eq!(got.len(), expected.len());
        for (g, e) in got.iter().zip(expected) {
            close(*e, *g, 0.001);
        }
    };
    check(d(0), &[1.0, 1.0, 20.0 / 22.5, 25.0 / 22.5, 1.0, 1.0]);
    check(d(1), &[1.0, 1.0, 1.0, 1.0, 30.0 / 25.0, 1.0, 1.0]);
    check(d(2), &[-1.0, 0.5, 0.5, 0.5 * 50.0 / 30.0, 0.5, 0.5 * 30.0 / 35.0, 0.5, 0.5 * 40.0 / 35.0]);
    check(d(3), &[0.1, 0.1, 0.15, -1.0, 0.1, 0.15]);
    let names: Vec<&str> = medians.iter().map(|m| m.chromosome.as_str()).collect();
    assert_eq!(names, vec!["chr1", "chr4", "chrX", "chrY"]);
}

// ReadDepthStatisticsNormaliserTest

fn depth_factor(list: Vec<(&str, f64, f64)>) -> f64 {
    let mut by_chr: Vec<ChrRatios> = Vec::new();
    for (chr, depth, gc) in list {
        let r = br(1, depth, gc, true);
        match by_chr.iter_mut().find(|x| x.chr == c(chr)) {
            Some(x) => x.ratios.push(r),
            None => by_chr.push(chr_ratios(chr, vec![r])),
        }
    }
    let (mean, median) = read_depth_statistics_pub(&by_chr, 0.24, 0.68);
    mean / median
}

#[test]
fn read_depth_statistics_normaliser() {
    let mut r = br(1, 1.0, 0.5, true);
    r.normalise(depth_factor(vec![("1", 10.0, 0.5), ("1", 10.0, 0.5), ("1", 100.0, 0.5)]));
    close(0.25, r.ratio(), 1e-4);
    // valuesWithExtremeGcAreIgnored
    let f = depth_factor(vec![
        ("1", 100.0, 0.23),
        ("1", 100.0, 0.23999),
        ("1", 10.0, 0.24),
        ("1", 10.0, 0.25),
        ("1", 10.0, 0.26),
        ("1", 10.0, 0.34),
        ("1", 10.0, 0.54),
        ("1", 10.0, 0.67),
        ("1", 10.0, 0.68),
        ("1", 1000.0, 0.680001),
        ("1", 1000.0, 0.69),
    ]);
    close(1.0, f, 1e-4);
    // sexChromosomeValuesAreIgnored
    let f = depth_factor(vec![("1", 10.0, 0.44), ("1", 10.0, 0.45), ("1", 10.0, 0.46), ("1", 10.0, 0.34), ("1", 10.0, 0.54), ("X", 1000.0, 0.67), ("Y", 10000.0, 0.48)]);
    close(1.0, f, 1e-4);
    // nonPositiveValuesIgnored
    let f = depth_factor(vec![("1", 10.0, 0.44), ("1", 10.0, 0.45), ("1", 10.0, 0.46), ("1", 0.0, 0.34), ("1", -1.0, 0.54), ("1", -0.0001, 0.54)]);
    close(1.0, f, 1e-4);
}

// UnityNormaliserTest

fn unity(list: Vec<(&str, BamRatio)>) -> Vec<BamRatio> {
    let mut by_chr: Vec<ChrRatios> = Vec::new();
    let mut order = Vec::new();
    for (chr, r) in list {
        let i = match by_chr.iter().position(|x| x.chr == c(chr)) {
            Some(i) => i,
            None => {
                by_chr.push(chr_ratios(chr, vec![]));
                by_chr.len() - 1
            }
        };
        order.push((i, by_chr[i].ratios.len()));
        by_chr[i].ratios.push(r);
    }
    unity_normaliser_pub(&mut by_chr);
    order.iter().map(|&(i, j)| by_chr[i].ratios[j].clone()).collect()
}

#[test]
fn unity_normaliser_record() {
    let out = unity(vec![("1", br(1001, 10.0, 0.45, true)), ("1", br(2001, 20.0, 0.46, false)), ("3", br(2001, 30.0, 0.46, true))]);
    close(0.5, out[0].ratio(), 0.001);
    assert_eq!(out[1].ratio(), -1.0);
    close(1.5, out[2].ratio(), 0.001);
}

#[test]
fn unity_normaliser_megabase_scale_values() {
    let mut a = br(1001, 10.0, 0.45, true);
    a.set_diploid_adjusted_ratio(10.0);
    let mut b = br(2001, 20.0, 0.46, false);
    b.set_diploid_adjusted_ratio(-1.0);
    let mut d = br(2001, 30.0, 0.46, true);
    d.set_diploid_adjusted_ratio(12.0);
    let out = unity(vec![("1", a), ("1", b), ("3", d)]);
    close(10.0 / 11.0, out[0].diploid_adjusted_ratio(), 0.001);
    assert_eq!(out[1].diploid_adjusted_ratio(), -1.0);
    close(12.0 / 11.0, out[2].diploid_adjusted_ratio(), 0.001);
}

#[test]
fn unity_normaliser_allosomes_and_zeros() {
    let out = unity(vec![("1", br(1001, 1.0, 0.45, true)), ("X", br(2001, 2.0, 0.46, true)), ("Y", br(2001, 3.0, 0.46, true))]);
    let got: Vec<f64> = out.iter().map(|r| r.ratio()).collect();
    assert_eq!(got, vec![1.0, 2.0, 3.0]);
    // handle0s
    let mut masked = br(2001, 10.0, 0.46, true);
    masked.normalise(-1.0);
    let out = unity(vec![("1", br(1001, 10.0, 0.45, true)), ("1", br(2001, 0.0, 0.46, true)), ("3", br(2001, 50.0, 0.46, true)), ("3", masked)]);
    close(10.0 / 30.0, out[0].ratio(), 0.001);
    assert_eq!(out[1].ratio(), 0.0);
    close(50.0 / 30.0, out[2].ratio(), 0.001);
    assert_eq!(out[3].ratio(), -1.0);
}

// ResultsCollatorTest

fn collator_inputs() -> (Vec<ChrRatios>, Vec<ChrRatios>) {
    let t = |pos: i32, depth: f64, gc: f64| {
        let mut r = br(pos, depth, gc, true);
        if depth > 0.0 {
            r.normalise(100.0);
        }
        r
    };
    let r = |pos: i32, depth: f64, gc: f64, diploid: f64| {
        let mut x = t(pos, depth, gc);
        x.set_diploid_adjusted_ratio(diploid);
        x
    };
    let tumor = vec![
        chr_ratios("1", vec![t(1, -1.0, -1.0), t(1001, -1.0, -1.0), t(2001, 20.0, 0.40), t(3001, 25.0, 0.40), t(4001, 26.0, 0.40), t(5001, 24.0, 0.39), t(6001, -1.0, 0.39)]),
        chr_ratios("2", vec![t(1, -1.0, -1.0), t(1001, 35.0, 0.53), t(2001, 32.0, 0.54), t(3001, 25.0, 0.40), t(4001, -1.0, 0.40)]),
    ];
    let reference = vec![
        chr_ratios(
            "1",
            vec![
                r(1, -1.0, -1.0, -1.0),
                r(1001, -1.0, -1.0, -1.0),
                r(2001, -1.0, -1.0, -1.0),
                r(3001, 10.0, 0.40, 0.12),
                r(4001, 9.0, 0.40, 0.13),
                r(5001, 11.0, 0.40, 0.10),
                r(6001, 15.0, 0.48, 0.11),
            ],
        ),
        chr_ratios("2", vec![r(1, -1.0, -1.0, -1.0), r(1001, -1.0, -1.0, -1.0), r(2001, 21.0, 0.32, -1.0), r(3001, -1.0, -1.0, 0.12), r(4001, 28.0, 0.36, 0.28)]),
    ];
    (tumor, reference)
}

fn row_values(r: &CobaltRatio) -> [f64; 7] {
    [r.reference_read_depth, r.tumor_read_depth, r.reference_gc_ratio, r.tumor_gc_ratio, r.reference_gc_diploid_ratio, r.reference_gc_content, r.tumor_gc_content]
}

fn check_rows(collated: &[CollatedChr], expected: &[(&str, i32, [f64; 7])]) {
    let rows: Vec<CobaltRatio> = collated.iter().flat_map(|c| (0..c.len()).map(|i| c.row(i)).collect::<Vec<_>>()).collect();
    assert_eq!(rows.len(), expected.len());
    for (r, (chr, pos, values)) in rows.iter().zip(expected) {
        assert_eq!((r.chr, r.position), (c(chr), *pos));
        for (g, e) in row_values(r).iter().zip(values) {
            close(*e, *g, 1e-4);
        }
    }
}

#[test]
fn results_collator_collate() {
    let (tumor, reference) = collator_inputs();
    let collated = collate(&tumor, &reference).unwrap();
    let m = -1.0;
    check_rows(
        &collated,
        &[
            ("1", 1, [m, m, m, m, m, m, m]),
            ("1", 1001, [m, m, m, m, m, m, m]),
            ("1", 2001, [m, 20.0, m, 0.20, m, m, 0.40]),
            ("1", 3001, [10.0, 25.0, 0.10, 0.25, 0.12, 0.40, 0.40]),
            ("1", 4001, [9.0, 26.0, 0.09, 0.26, 0.13, 0.40, 0.40]),
            ("1", 5001, [11.0, 24.0, 0.11, 0.24, 0.10, 0.40, 0.39]),
            ("1", 6001, [15.0, m, 0.15, m, 0.11, 0.48, 0.39]),
            ("2", 1, [m, m, m, m, m, m, m]),
            ("2", 1001, [m, 35.0, m, 0.35, m, m, 0.53]),
            ("2", 2001, [21.0, 32.0, 0.21, 0.32, m, 0.32, 0.54]),
            ("2", 3001, [m, 25.0, m, 0.25, 0.12, m, 0.40]),
            ("2", 4001, [28.0, m, 0.28, m, 0.28, 0.36, 0.40]),
        ],
    );
}

#[test]
fn results_collator_reference_and_tumor_only() {
    let (tumor, reference) = collator_inputs();
    let m = -1.0;
    let collated = collate(&[], &reference).unwrap();
    check_rows(
        &collated,
        &[
            ("1", 1, [m, m, m, m, m, m, m]),
            ("1", 1001, [m, m, m, m, m, m, m]),
            ("1", 2001, [m, m, m, m, m, m, m]),
            ("1", 3001, [10.0, m, 0.10, m, 0.12, 0.40, m]),
            ("1", 4001, [9.0, m, 0.09, m, 0.13, 0.40, m]),
            ("1", 5001, [11.0, m, 0.11, m, 0.10, 0.40, m]),
            ("1", 6001, [15.0, m, 0.15, m, 0.11, 0.48, m]),
            ("2", 1, [m, m, m, m, m, m, m]),
            ("2", 1001, [m, m, m, m, m, m, m]),
            ("2", 2001, [21.0, m, 0.21, m, m, 0.32, m]),
            ("2", 3001, [m, m, m, m, 0.12, m, m]),
            ("2", 4001, [28.0, m, 0.28, m, 0.28, 0.36, m]),
        ],
    );
    let collated = collate(&tumor, &[]).unwrap();
    check_rows(
        &collated,
        &[
            ("1", 1, [m, m, m, m, m, m, m]),
            ("1", 1001, [m, m, m, m, m, m, m]),
            ("1", 2001, [m, 20.0, m, 0.20, m, m, 0.40]),
            ("1", 3001, [m, 25.0, m, 0.25, m, m, 0.40]),
            ("1", 4001, [m, 26.0, m, 0.26, m, m, 0.40]),
            ("1", 5001, [m, 24.0, m, 0.24, m, m, 0.39]),
            ("1", 6001, [m, m, m, m, m, m, 0.39]),
            ("2", 1, [m, m, m, m, m, m, m]),
            ("2", 1001, [m, 35.0, m, 0.35, m, m, 0.53]),
            ("2", 2001, [m, 32.0, m, 0.32, m, m, 0.54]),
            ("2", 3001, [m, 25.0, m, 0.25, m, m, 0.40]),
            ("2", 4001, [m, m, m, m, m, m, 0.40]),
        ],
    );
}

// WindowStatusTest, WindowStatusesTest

fn profile(start: i32, gc: f64, mappable: f64) -> GcProfile {
    GcProfile { start, gc_content: gc, mappable_percentage: mappable }
}

#[test]
fn window_status_masked_by_mappability() {
    for (m, masked) in [(0.99, false), (0.87, false), (0.86, false), (0.85, false), (0.84, true), (0.83, true), (0.01, true)] {
        assert_eq!(!profile(1001, 0.47, m).is_mappable(), masked, "{m}");
    }
}

#[test]
fn window_statuses() {
    let mut profiles: PerChr<Vec<GcProfile>> = PerChr::default();
    *profiles.get_mut(c("1")) =
        vec![profile(1, 0.47, 0.70), profile(1001, 0.48, 0.50), profile(2001, 0.49, 1.0), profile(3001, 0.50, 1.0), profile(4001, 0.51, 1.0), profile(5001, 0.52, 1.0)];
    *profiles.get_mut(c("2")) = vec![profile(1, 0.57, 1.0), profile(1001, 0.58, 1.0), profile(2001, 0.59, 0.50), profile(3001, 0.57, 1.0)];
    *profiles.get_mut(c("3")) = [0.47, 0.48, 0.49, 0.50, 0.47, 0.48, 0.49, 0.50, 0.47, 0.48, 0.49, 0.50]
        .iter()
        .enumerate()
        .map(|(i, gc)| profile(i as i32 * 1000 + 1, *gc, if i == 3 { 0.8 } else { 1.0 }))
        .collect();
    let regions: Vec<crate::gc::Region> = [("chr1", 1100, 1200), ("chr1", 1300, 1400), ("chr1", 4100, 4200), ("chr3", 1100, 1200), ("chr3", 7100, 7200)]
        .iter()
        .map(|(chr, s, e)| crate::gc::Region { chromosome: chr.to_string(), start: *s, end: *e })
        .collect();
    let excluded = crate::gc::excluded_windows(&regions, &profiles).unwrap();
    let mut diploid: PerChr<Vec<bool>> = PerChr::default();
    *diploid.get_mut(c("1")) = vec![true; 7];
    *diploid.get_mut(c("2")) = vec![true; 4];
    *diploid.get_mut(c("3")) = vec![true, true, true, false, true, true, false, false, true, true, true];
    let s = WindowStatuses::new(profiles, &excluded, Some(&diploid));
    // gcReferenceValueTest
    assert_eq!(s.reference_gc(c("1"), 0), Some(0.47));
    assert_eq!(s.reference_gc(c("1"), 1000), Some(0.48));
    assert_eq!(s.reference_gc(c("2"), 1), Some(0.57));
    assert_eq!(s.reference_gc(c("3"), 1000), Some(0.48));
    // excludeTest
    let excluded_at = [
        ("1", 1, true),
        ("1", 1001, true),
        ("1", 2001, false),
        ("1", 3001, false),
        ("1", 4001, true),
        ("1", 5001, false),
        ("2", 1, false),
        ("2", 1001, false),
        ("2", 2001, true),
        ("2", 3001, false),
        ("3", 1, false),
        ("3", 1001, true),
        ("3", 2001, false),
        ("3", 3001, true),
        ("3", 4001, false),
        ("3", 5001, false),
        ("3", 6001, true),
        ("3", 7001, true),
        ("3", 8001, false),
        ("3", 9001, false),
        ("3", 10001, false),
        ("3", 11001, true),
    ];
    for (chr, pos, e) in excluded_at {
        assert_eq!(s.exclude(c(chr), pos).unwrap(), e, "{chr}:{pos}");
    }
}

// ExcludedRegionsFileTest, SuppliedExcludedRegionsTest

#[test]
fn supplied_excluded_regions_find_intersections() {
    let mut profiles: PerChr<Vec<GcProfile>> = PerChr::default();
    for i in 0..10 {
        profiles.get_mut(c("1")).push(profile(i * 1000 + 1, 0.5, 0.9));
        profiles.get_mut(c("3")).push(profile(i * 1000 + 1, 0.5, 0.9));
    }
    for i in 10..20 {
        profiles.get_mut(c("1")).push(profile(i * 1000 + 1, 0.5, 0.9));
        profiles.get_mut(c("2")).push(profile(i * 1000 + 1, 0.5, 0.9));
    }
    let regions: Vec<crate::gc::Region> = [("1", 100, 120), ("1", 200, 220), ("1", 1200, 1220), ("1", 13000, 15000), ("2", 1100, 1220), ("2", 11000, 12020)]
        .iter()
        .map(|(chr, s, e)| crate::gc::Region { chromosome: chr.to_string(), start: *s, end: *e })
        .collect();
    let ex = crate::gc::excluded_windows(&regions, &profiles).unwrap();
    let starts = |chr: &str| -> Vec<i32> {
        ex.get(c(chr)).iter().zip(profiles.get(c(chr))).filter(|(e, _)| **e).map(|(_, p)| p.start).collect()
    };
    assert_eq!(starts("1"), vec![1, 1_001, 12_001, 13_001, 14_001]);
    // chromosome 2's profiles start at 10,001: the region 1100-1220 overlaps none
    assert_eq!(starts("2"), vec![10_001, 11_001, 12_001]);
    assert!(starts("3").is_empty());
}

// DiploidRegionLoaderTest

#[test]
fn diploid_region_loader() {
    let dir = std::env::temp_dir().join(format!("cobalt_rs_diploid_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("diploid.bed");
    // DiploidFileSection writes "chr<n>\tstart\tstop"
    std::fs::write(&path, "chr1\t1000\t3000\nchr1\t5000\t6000\nchr2\t1000\t3000\nchr2\t7000\t13000").unwrap();
    let r = crate::gc::load_diploid_regions(path.to_str().unwrap()).unwrap();
    assert_eq!(r.get(c("1")), &vec![false, true, true, false, false, true]);
    let r2 = r.get(c("2"));
    assert_eq!(r2.len(), 13);
    assert_eq!(&r2[..9], &[false, true, true, false, false, false, false, true, true]);
    assert!(r2[12]);
}

// TargetRegionsTest, WholeGenomeTest

#[test]
fn target_regions_enrichment() {
    let mut slots: PerChr<Option<Vec<Option<f64>>>> = PerChr::default();
    let mut one = vec![None; 15];
    for (pos, e) in [(5_001, 0.5), (6_001, 0.6), (7_001, 0.5), (10_001, 0.5), (11_001, 0.4)] {
        one[(pos - 1) / 1000] = Some(e);
    }
    *slots.get_mut(c("1")) = Some(one);
    let mut two = vec![None; 10];
    two[3] = Some(0.52);
    two[4] = Some(0.62);
    *slots.get_mut(c("2")) = Some(two);
    let t = TargetRegions::from_slots(slots);
    let s = Scope::Targeted(&t);
    let check = |e: f64, chr: &str, pos: i32| assert_eq!(s.region_enrichment(c(chr), pos).unwrap(), e, "{chr}:{pos}");
    for pos in [1, 1001, 2001, 3001, 4001, 8001, 9001, 12001] {
        check(-1.0, "1", pos);
    }
    check(0.5, "1", 5001);
    check(0.6, "1", 6001);
    check(0.4, "1", 11001);
    check(0.52, "2", 3001);
    check(-1.0, "2", 5001);
    check(-1.0, "3", 1);
    assert!(s.on_target(c("1"), 5001).unwrap() && !s.on_target(c("1"), 1).unwrap() && !s.on_target(c("3"), 1).unwrap());
    // WholeGenomeTest
    let w = Scope::WholeGenome;
    assert!(w.on_target(c("1"), 1).unwrap());
    assert_eq!(w.region_enrichment(c("5"), 1001).unwrap(), 1.0);
    assert!(matches!(Consolidator::for_whole_genome(f64::NAN), Consolidator::NoOp));
    assert!(matches!(Consolidator::for_whole_genome(10.0), Consolidator::NoOp));
    assert!(matches!(Consolidator::for_whole_genome(3.0), Consolidator::LowCoverage { count: 30, .. }));
}

#[test]
fn merge_tumor_ratios() {
    let mut all = vec![chr_ratios("1", vec![br(1, -1.0, 0.5, true), br(1001, 2.0, 0.5, true), br(2001, -1.0, 0.5, true)])];
    let panel = vec![chr_ratios("1", vec![br(1, 3.0, 0.5, true), br(1001, 4.0, 0.5, true), br(2001, -1.0, 0.5, true)])];
    merge_panel_ratios(&mut all, panel).unwrap();
    let got: Vec<f64> = all[0].ratios.iter().map(|r| r.ratio()).collect();
    assert_eq!(got, vec![3.0, 2.0, -1.0]);
}

// CobaltCalculatorTest

const TD1: [f64; 15] = [-1.0, -1.0, -1.0, 25.0, 25.0, 24.0, 22.0, 23.0, 24.0, 22.0, 26.0, 26.0, 27.0, 23.0, 25.0];
const RD1: [f64; 15] = [-1.0, -1.0, -1.0, 14.0, 12.0, 12.0, 16.0, 16.0, 17.0, 18.0, 18.0, 19.0, 12.0, 12.0, 15.0];
const TGC1: [f64; 15] = [0.60, 0.40, 0.41, 0.40, 0.40, 0.40, 0.43, 0.42, 0.42, 0.43, 0.41, 0.41, 0.41, 0.42, 0.43];
const RGC1: [f64; 15] = [0.60, 0.50, 0.51, 0.50, 0.50, 0.50, 0.51, 0.51, 0.51, 0.52, 0.52, 0.52, 0.53, 0.53, 0.53];
const TD2: [f64; 6] = [26.0, 30.0, 30.0, 34.0, 28.0, 32.0];
const RD2: [f64; 6] = [16.0, 20.0, 18.0, 22.0, 20.0, 24.0];
const TGC2: [f64; 6] = [0.44, 0.44, 0.46, 0.46, 0.45, 0.45];
const RGC2: [f64; 6] = [0.54, 0.54, 0.55, 0.55, 0.56, 0.56];

/// counts giving each window the depth and GC content (bases = depth * 1000, G/C bases = GC * bases)
fn counts(chrs: &[(&str, &[f64], &[f64])]) -> SampleCounts {
    SampleCounts {
        chromosomes: chrs
            .iter()
            .map(|(chr, depth, gc)| {
                let bases: Vec<i32> = depth.iter().map(|d| (d * 1000.0).round() as i32).collect();
                let gcs: Vec<i32> = bases.iter().zip(gc.iter()).map(|(b, g)| (*b as f64 * g).round() as i32).collect();
                ChrCounts { chr: c(chr), contig: format!("chr{chr}"), bases, gc: gcs }
            })
            .collect(),
    }
}

fn calculator_statuses() -> WindowStatuses {
    let mut profiles: PerChr<Vec<GcProfile>> = PerChr::default();
    // the test's profiles start at 0, 1000, ...: the window index is start / 1000 either way
    *profiles.get_mut(c("1")) = (0..15).map(|i| profile(i * 1000, 0.5, if i == 1 { 0.80 } else { 1.0 })).collect();
    *profiles.get_mut(c("2")) = (0..6).map(|i| profile(i * 1000, 0.5, 1.0)).collect();
    let regions = vec![
        crate::gc::Region { chromosome: "chr1".into(), start: 2010, end: 2160 },
        crate::gc::Region { chromosome: "chr1".into(), start: 2210, end: 2360 },
    ];
    let excluded = crate::gc::excluded_windows(&regions, &profiles).unwrap();
    WindowStatuses::new(profiles, &excluded, None)
}

fn cfg() -> CalcConfig {
    CalcConfig { gc_bucket_min: 24, gc_bucket_max: 68, gc_ratio_min: 0.24, gc_ratio_max: 0.68 }
}

fn gc_ratio(map: &[(f64, f64)], gc: f64) -> f64 {
    map.iter().find(|(g, _)| *g == gc).unwrap().1
}

const REF_GC: [(f64, f64); 9] =
    [(0.49, 12.0 / 3.0), (0.50, 28.0 / 3.0), (0.51, 46.0 / 3.0), (0.52, 46.0 / 3.0), (0.53, 48.0 / 3.0), (0.54, 50.0 / 3.0), (0.55, 60.0 / 3.0), (0.56, 42.0 / 3.0), (0.57, 22.0 / 3.0)];
const TUM_GC: [(f64, f64); 9] =
    [(0.39, 25.0 / 3.0), (0.40, 51.0 / 3.0), (0.41, 74.0 / 3.0), (0.42, 71.0 / 3.0), (0.43, 73.0 / 3.0), (0.44, 80.0 / 3.0), (0.45, 90.0 / 3.0), (0.46, 62.0 / 3.0), (0.47, 32.0 / 3.0)];

fn mean_median_ratio(a: &[f64], b: &[f64]) -> f64 {
    let v: Vec<f64> = a[3..].iter().chain(b.iter()).copied().collect();
    crate::stats::mean(&v) / crate::stats::median(&mut v.clone())
}

fn run_calculator(tumor: &SampleCounts, reference: Option<&SampleCounts>, scope: &Scope) -> (Vec<ChrRatios>, Vec<ChrRatios>, Vec<MedianRatio>) {
    let statuses = calculator_statuses();
    let t = calculate(tumor, &statuses, scope, &cfg(), Sample::Tumor).unwrap();
    let mut consolidator = t.consolidator.unwrap();
    let targeted = matches!(scope, Scope::Targeted(_));
    let (mut reference_ratios, mut medians) = (Vec::new(), Vec::new());
    if let Some(r) = reference {
        let rc = calculate(r, &statuses, scope, &cfg(), Sample::Reference { consolidator: &mut consolidator, diploid_normalise: !targeted, v38: true }).unwrap();
        reference_ratios = rc.ratios;
        medians = if targeted { panel_median_ratios(&mut reference_ratios, true) } else { rc.median_ratios.unwrap() };
    }
    (t.ratios, reference_ratios, medians)
}

fn check_ratio(r: &CobaltRatio, chr: &str, pos: i32, values: [f64; 7]) {
    assert_eq!((r.chr, r.position), (c(chr), pos));
    for (g, e) in row_values(r).iter().zip(values) {
        close(e, *g, 1e-4);
    }
}

#[test]
fn cobalt_calculator_tumor_and_reference_whole_genome() {
    let tumor = counts(&[("1", &TD1, &TGC1), ("2", &TD2, &TGC2)]);
    let reference = counts(&[("1", &RD1, &RGC1), ("2", &RD2, &RGC2)]);
    let (t, r, medians) = run_calculator(&tumor, Some(&reference), &Scope::WholeGenome);
    let collated = collate(&t, &r).unwrap();
    assert_eq!(collated.len(), 2);
    let (tm, rm) = (mean_median_ratio(&TD1, &TD2), mean_median_ratio(&RD1, &RD2));
    let rows1: Vec<CobaltRatio> = (0..collated[0].len()).map(|i| collated[0].row(i)).collect();
    assert_eq!(rows1.len(), 15);
    for i in 0..3 {
        check_ratio(&rows1[i], "1", i as i32 * 1000 + 1, [RD1[i], TD1[i], -1.0, -1.0, -1.0, RGC1[i], TGC1[i]]);
    }
    let mut ref1 = Vec::new();
    for i in 3..15 {
        let ref_ratio = RD1[i] / gc_ratio(&REF_GC, RGC1[i]) / rm;
        ref1.push(ref_ratio);
        let tum_ratio = TD1[i] / gc_ratio(&TUM_GC, TGC1[i]) / tm;
        check_ratio(&rows1[i], "1", i as i32 * 1000 + 1, [RD1[i], TD1[i], ref_ratio, tum_ratio, ref_ratio, RGC1[i], TGC1[i]]);
    }
    let mut ref2 = Vec::new();
    for i in 0..6 {
        let ref_ratio = RD2[i] / gc_ratio(&REF_GC, RGC2[i]) / rm;
        ref2.push(ref_ratio);
        let tum_ratio = TD2[i] / gc_ratio(&TUM_GC, TGC2[i]) / tm;
        check_ratio(&collated[1].row(i), "2", i as i32 * 1000 + 1, [RD2[i], TD2[i], ref_ratio, tum_ratio, ref_ratio, RGC2[i], TGC2[i]]);
    }
    assert_eq!(medians.len(), 2);
    assert_eq!((medians[0].chromosome.as_str(), medians[0].count), ("chr1", 12));
    close(crate::constants::doubles::median(&mut ref1), medians[0].median_ratio, 0.001);
    assert_eq!((medians[1].chromosome.as_str(), medians[1].count), ("chr2", 6));
    close(crate::constants::doubles::median(&mut ref2), medians[1].median_ratio, 0.001);
}

#[test]
fn cobalt_calculator_tumor_only_and_reference_only_whole_genome() {
    let tumor = counts(&[("1", &TD1, &TGC1), ("2", &TD2, &TGC2)]);
    let (t, _, medians) = run_calculator(&tumor, None, &Scope::WholeGenome);
    let collated = collate(&t, &[]).unwrap();
    let tm = mean_median_ratio(&TD1, &TD2);
    for i in 0..3 {
        check_ratio(&collated[0].row(i), "1", i as i32 * 1000 + 1, [-1.0, TD1[i], -1.0, -1.0, -1.0, -1.0, TGC1[i]]);
    }
    for i in 3..15 {
        let tum_ratio = TD1[i] / gc_ratio(&TUM_GC, TGC1[i]) / tm;
        check_ratio(&collated[0].row(i), "1", i as i32 * 1000 + 1, [-1.0, TD1[i], -1.0, tum_ratio, -1.0, -1.0, TGC1[i]]);
    }
    assert!(medians.is_empty());
    // referenceOnlyWholeGenomeTest
    let reference = counts(&[("1", &RD1, &RGC1), ("2", &RD2, &RGC2)]);
    let empty = SampleCounts { chromosomes: Vec::new() };
    let (t, r, medians) = run_calculator(&empty, Some(&reference), &Scope::WholeGenome);
    let collated = collate(&t, &r).unwrap();
    let rm = mean_median_ratio(&RD1, &RD2);
    for i in 0..3 {
        check_ratio(&collated[0].row(i), "1", i as i32 * 1000 + 1, [RD1[i], -1.0, -1.0, -1.0, -1.0, RGC1[i], -1.0]);
    }
    for i in 3..15 {
        let ref_ratio = RD1[i] / gc_ratio(&REF_GC, RGC1[i]) / rm;
        check_ratio(&collated[0].row(i), "1", i as i32 * 1000 + 1, [RD1[i], -1.0, ref_ratio, -1.0, ref_ratio, RGC1[i], -1.0]);
    }
    assert_eq!(medians.iter().map(|m| m.count).collect::<Vec<_>>(), vec![12, 6]);
}

fn calculator_panel() -> TargetRegions {
    // CobaltCalculatorTest's scope: chromosome 1 on target from 4001 (enrichment 1.1), chromosome 2 everywhere (1.2)
    let mut slots: PerChr<Option<Vec<Option<f64>>>> = PerChr::default();
    *slots.get_mut(c("1")) = Some((0..15).map(|i| if i >= 4 { Some(1.1) } else { None }).collect());
    *slots.get_mut(c("2")) = Some(vec![Some(1.2); 6]);
    TargetRegions::from_slots(slots)
}

#[test]
fn cobalt_calculator_targeted() {
    let tumor = counts(&[("1", &TD1, &TGC1), ("2", &TD2, &TGC2)]);
    let reference = counts(&[("1", &RD1, &RGC1), ("2", &RD2, &RGC2)]);
    let panel = calculator_panel();
    let scope = Scope::Targeted(&panel);
    let mut tum_gc_targeted = TUM_GC.to_vec();
    for (gc, v) in [(0.39, 24.5 / 3.0), (0.40, (24.5 + 26.0) / 3.0), (0.41, (24.5 + 26.0 + 23.0) / 3.0)] {
        tum_gc_targeted.iter_mut().find(|x| x.0 == gc).unwrap().1 = v;
    }
    let (ef1, ef2) = (1.1, 1.2);
    let raw_tumor: Vec<f64> = (4..15)
        .map(|i| TD1[i] / (gc_ratio(&tum_gc_targeted, TGC1[i]) * ef1))
        .chain((0..6).map(|i| TD2[i] / (gc_ratio(&tum_gc_targeted, TGC2[i]) * ef2)))
        .collect();
    let raw_reference: Vec<f64> =
        (4..15).map(|i| RD1[i] / (gc_ratio(&REF_GC, RGC1[i]) * ef1)).chain((0..6).map(|i| RD2[i] / (gc_ratio(&REF_GC, RGC2[i]) * ef2))).collect();
    let tf = raw_tumor.iter().sum::<f64>() / raw_tumor.len() as f64;
    let rf = raw_reference.iter().sum::<f64>() / raw_reference.len() as f64;
    let (t, r, _) = run_calculator(&tumor, Some(&reference), &scope);
    let collated = collate(&t, &r).unwrap();
    for i in 0..4 {
        check_ratio(&collated[0].row(i), "1", i as i32 * 1000 + 1, [RD1[i], TD1[i], -1.0, -1.0, -1.0, RGC1[i], TGC1[i]]);
    }
    for i in 4..15 {
        let ref_ratio = RD1[i] / (gc_ratio(&REF_GC, RGC1[i]) * ef1) / rf;
        let tum_ratio = TD1[i] / (gc_ratio(&tum_gc_targeted, TGC1[i]) * ef1 * tf);
        check_ratio(&collated[0].row(i), "1", i as i32 * 1000 + 1, [RD1[i], TD1[i], ref_ratio, tum_ratio, ref_ratio, RGC1[i], TGC1[i]]);
    }
    for i in 0..6 {
        let ref_ratio = RD2[i] / (gc_ratio(&REF_GC, RGC2[i]) * ef2) / rf;
        let tum_ratio = TD2[i] / (gc_ratio(&tum_gc_targeted, TGC2[i]) * ef2 * tf);
        check_ratio(&collated[1].row(i), "2", i as i32 * 1000 + 1, [RD2[i], TD2[i], ref_ratio, tum_ratio, ref_ratio, RGC2[i], TGC2[i]]);
    }
    // tumorOnlyTargetedTest
    let (t, _, medians) = run_calculator(&tumor, None, &scope);
    let collated = collate(&t, &[]).unwrap();
    for i in 4..15 {
        let tum_ratio = TD1[i] / (gc_ratio(&tum_gc_targeted, TGC1[i]) * ef1 * tf);
        check_ratio(&collated[0].row(i), "1", i as i32 * 1000 + 1, [-1.0, TD1[i], -1.0, tum_ratio, -1.0, -1.0, TGC1[i]]);
    }
    assert!(medians.is_empty());
}

// CobaltRatioSegmenterTest (the shared segmentation with COBALT's values and a locator putting positions
// above 20,000 on the q arm)

fn segment_test(points: &[(&str, i32, f64)]) -> std::collections::BTreeMap<hmf_common::segmentation::ChrArm, Vec<hmf_common::segmentation::PcfSegment>> {
    use hmf_common::segmentation::*;
    let pts: Vec<SegmentPoint> = points.iter().map(|(chr, pos, v)| SegmentPoint { chromosome: chr.to_string(), position: *pos, value: *v }).collect();
    let locator = |chr: &str, pos: i32| -> anyhow::Result<ChrArm> {
        Ok(ChrArm { chromosome: HumanChromosome::parse(chr)?, arm: if pos > 20_000 { Arm::Q } else { Arm::P } })
    };
    segment_with_locator(&pts, &locator, 100.0, WINDOWED, 2).unwrap()
}

/// SegmentationTestBase.ratioForValue: 2^value, -1 for a negative value
fn rv(v: f64) -> f64 {
    if v >= 0.0 { 2f64.powf(v) } else { -1.0 }
}

#[test]
fn cobalt_ratio_segmenter_write_segmentation() {
    let values = [
        (1, 0.5),
        (1001, 0.51),
        (2001, 0.49),
        (3001, -1.0),
        (4001, -1.0),
        (5001, -1.0),
        (6001, -1.0),
        (13001, 10.09),
        (14001, 9.92),
        (15001, 10.02),
        (16001, -1.0),
        (20_001, 0.5),
        (21_001, 0.51),
        (22_001, 0.49),
        (23_001, 0.52),
    ];
    let mut points = Vec::new();
    for chr in ["1", "2", "3"] {
        for (p, v) in values {
            points.push((chr, p, rv(v)));
        }
    }
    let segments = segment_test(&points);
    let chr1: Vec<&hmf_common::segmentation::PcfSegment> = segments.iter().filter(|(a, _)| a.chromosome.short_name() == "1").flat_map(|(_, l)| l).collect();
    assert_eq!(chr1.len(), 3);
    assert_eq!((chr1[0].start, chr1[0].end), (1, 3000));
    close((rv(0.5) + rv(0.51) + rv(0.49)) / 3.0, chr1[0].mean_ratio, 1e-4);
    assert_eq!((chr1[1].start, chr1[1].end), (13001, 16000));
    close((rv(10.09) + rv(9.92) + rv(10.02)) / 3.0, chr1[1].mean_ratio, 1e-4);
    assert_eq!((chr1[2].start, chr1[2].end), (20001, 24000));
}

#[test]
fn cobalt_ratio_segmenter_masked_windows_and_zeros() {
    let values = [
        (1, 0.5),
        (1001, -1.0),
        (2001, -1.0),
        (3001, 0.51),
        (4001, 0.49),
        (5001, -1.0),
        (6001, -1.0),
        (10001, 10.09),
        (11001, -1.0),
        (12001, -1.0),
        (13001, 9.92),
        (14001, -1.0),
        (15001, 10.02),
        (16001, -1.0),
        (20_001, -1.0),
        (21_001, 0.5),
        (22_001, 0.51),
        (23_001, 0.49),
        (24_001, -1.0),
        (25_001, -1.0),
        (26_001, 0.52),
    ];
    let points: Vec<(&str, i32, f64)> = values.iter().map(|(p, v)| ("1", *p, rv(*v))).collect();
    let segments = segment_test(&points);
    let s: Vec<(i32, i32)> = segments.values().flatten().map(|x| (x.start, x.end)).collect();
    assert_eq!(s, vec![(1, 5000), (10001, 16000), (21001, 27000)]);
    // handleZeroes
    let values = [
        (1, 0.5),
        (1001, -1.0),
        (2001, -1.0),
        (3001, 0.51),
        (4001, 0.49),
        (5001, -1.0),
        (6001, -1.0),
        (10001, 100.09),
        (11001, 120.0),
        (12001, -1.0),
        (13001, 119.92),
        (14001, -1.0),
        (15001, 0.02),
        (16001, 0.0),
        (17001, 0.0),
    ];
    let points: Vec<(&str, i32, f64)> = values.iter().map(|(p, v)| ("1", *p, *v)).collect();
    let segments = segment_test(&points);
    let s: Vec<f64> = segments.values().flatten().map(|x| x.mean_ratio).collect();
    assert_eq!(s.len(), 4);
    close(0.0, s[3], 1e-4);
}

#[test]
fn cobalt_ratio_segmenter_values() {
    use hmf_common::segmentation::cobalt_segmentation_value;
    // base2LogsOfValuesAreUsed, lowerBoundOnLoggedValues
    close(0.5, cobalt_segmentation_value(rv(0.5)), 0.001);
    close(0.7, cobalt_segmentation_value(rv(0.7)), 0.001);
    close(-9.966, cobalt_segmentation_value(0.0009999), 0.01);
    close(0.001f64.ln() / 2f64.ln(), cobalt_segmentation_value(0.001), 0.001);
    close(0.01f64.ln() / 2f64.ln(), cobalt_segmentation_value(0.01), 0.001);
}

// ReadDepthAccumulatorTest

#[test]
fn read_depth_accumulator_window_read_depth() {
    // two windows; a 1,000-base block from 501 with G/C at read indices below 400 and from 800, then a
    // block from 1 with G at every second base
    let mut bases = vec![b'A'; 1000];
    for (i, b) in bases.iter_mut().enumerate() {
        *b = if !(400..800).contains(&i) {
            if i % 2 == 0 { b'G' } else { b'C' }
        } else if i % 2 == 0 {
            b'A'
        } else {
            b'T'
        };
    }
    let (mut counts, mut gc) = (vec![0i32; 2], vec![0i32; 2]);
    crate::count::add_block_for_test(&mut counts, &mut gc, 501, 1000, &bases);
    assert_eq!((counts[0] as f64 / 1000.0, gc[0] as f64 / counts[0] as f64), (0.5, 0.8));
    assert_eq!((counts[1] as f64 / 1000.0, gc[1] as f64 / counts[1] as f64), (0.5, 0.4));
    let bases2: Vec<u8> = (0..1000).map(|i| if i % 2 == 0 { b'G' } else { b'A' }).collect();
    crate::count::add_block_for_test(&mut counts, &mut gc, 1, 1000, &bases2);
    close(1.5, counts[0] as f64 / 1000.0, 1e-10);
    close(0.6, gc[0] as f64 / counts[0] as f64, 1e-10);
    close(0.4, gc[1] as f64 / counts[1] as f64, 1e-10);
}
