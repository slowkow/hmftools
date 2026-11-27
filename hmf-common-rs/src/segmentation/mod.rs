//! hmf-common `segmentation`: PerArmSegmenter, ChrArmLocator, ChrArm, the absolute and windowed segments and
//! SegmentsFile, over the piecewise constant fit of `copynumber`. AMBER segments its tumor BAFs with it
//! (`write_absolute_segments`, as BAFSegmenter), COBALT its GC-corrected ratios (`write_windowed_segments`,
//! as CobaltRatioSegmenter); PURPLE reads the files back with `pcf::read_positions`.
//!
//! The segmentation is per chromosome arm. Values below 0 are dropped. One segment penalty serves every
//! arm: `gamma` times the squared median absolute deviation of all values (in arm order) from their running
//! median. Each arm is then fitted on its own, in parallel; the result does not depend on the thread count.
//! A segment runs from the position of its first point to that of its last (plus the 1,000-base window less
//! one for windowed values) and its value is the plain mean of the raw values in it.
//!
//! Java faults corrected here (see NOTES.md; `-java_compat` restores Java's behaviour):
//! - SG1: with no value left to segment, Java's Gamma throws and the tool fails; corrected, the file has
//!   the header only, as AMBER writes for an empty BAF list.
//! - SG2: infinite values pass Java's `value >= 0` filter and make the costs NaN (a wrong segmentation, or
//!   no end to the backtracking when the first value of an arm is infinite); corrected, they are dropped as
//!   invalid like negative values.
//! - SG3 (copynumber::windowed_median): Java's running-median heaps lose track of 0.0 and -0.0 in one
//!   window; corrected, the heaps take both as 0.0. Not reachable from AMBER or COBALT.

pub mod copynumber;
pub mod pcf;

use crate::compat::java;
use anyhow::{anyhow, bail, Context, Result};
use std::collections::BTreeMap;
use std::io::Write;

/// One value to segment, at a genome position.
#[derive(Clone, Debug)]
pub struct SegmentPoint {
    pub chromosome: String,
    pub position: i32,
    pub value: f64,
}

/// GCProfileFactory.WINDOW_SIZE: the window of a COBALT ratio, from its position on.
pub const WINDOW_SIZE: i32 = 1000;

/// The value COBALT segments in place of a ratio below 0.001 (log2 of 0.001, as the R script had it).
pub const COBALT_FLOOR_LOG2: f64 = -9.965784;

/// HumanChromosome: 1-22, X and Y, in that order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct HumanChromosome(u8);

const SHORT_NAMES: [&str; 24] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y",
];

impl HumanChromosome {
    pub const ALL: [HumanChromosome; 24] = {
        let mut a = [HumanChromosome(0); 24];
        let mut i = 0;
        while i < 24 {
            a[i] = HumanChromosome(i as u8);
            i += 1;
        }
        a
    };

    /// HumanChromosome.fromString: an optional "chr" prefix in any case, then exactly "1".."22", "X" or "Y"
    /// (Java throws for anything else).
    pub fn from_name(name: &str) -> Option<HumanChromosome> {
        let b = name.as_bytes();
        let bare = if b.len() >= 3 && b[..3].eq_ignore_ascii_case(b"chr") { &b[3..] } else { b };
        let digit = |c: u8| c.is_ascii_digit();
        let i = match bare {
            [b'X'] => 22,
            [b'Y'] => 23,
            [d] if digit(*d) && *d != b'0' => d - b'1',
            [t, u] if digit(*t) && digit(*u) && *t != b'0' => {
                let n = (t - b'0') * 10 + (u - b'0');
                if n > 22 {
                    return None;
                }
                n - 1
            }
            _ => return None,
        };
        Some(HumanChromosome(i))
    }

    /// HumanChromosome.fromString, with Java's exception as an error.
    pub fn parse(name: &str) -> Result<HumanChromosome> {
        Self::from_name(name).ok_or_else(|| anyhow!("No enum constant HumanChromosome for chromosome '{name}'"))
    }

    /// HumanChromosome.contains: a lower-case "chr" prefix stripped, then a number 1-22, X or Y. (Java throws
    /// on an empty name; that is false here.)
    pub fn contains(name: &str) -> bool {
        let bare = name.strip_prefix("chr").unwrap_or(name);
        if !bare.is_empty() && bare.bytes().all(|b| b.is_ascii_digit()) {
            return bare.parse::<i64>().is_ok_and(|n| (1..=22).contains(&n));
        }
        bare == "X" || bare == "Y"
    }

    /// The ordinal, 0 for chromosome 1 up to 23 for Y.
    pub fn ordinal(self) -> usize {
        self.0 as usize
    }

    pub fn short_name(self) -> &'static str {
        SHORT_NAMES[self.0 as usize]
    }

    /// RefGenomeVersion.versionedChromosome: "chr1" for 38, "1" for 37.
    pub fn versioned(self, v38: bool) -> String {
        if v38 { format!("chr{}", self.short_name()) } else { self.short_name().to_string() }
    }
}

/// RefGenomeCoordinates centromeres (hmf-common resources centromeres.37.tsv and centromeres.38.tsv).
const CENTROMERES_37: [i32; 24] = [
    123035434, 93826171, 92004854, 51160117, 47905641, 60330166, 59554331, 45338887, 48867679, 40754935, 53144205, 36356694, 17500000,
    17500000, 18500000, 36835801, 23763006, 16960898, 26181782, 27869569, 12788129, 14500000, 60132012, 11604553,
];
const CENTROMERES_38: [i32; 24] = [
    123605523, 93139351, 92214016, 50726026, 48272854, 59191911, 59498944, 44955505, 44377363, 40640102, 52751711, 35977330, 17025624,
    17086762, 18362627, 37295920, 24849830, 18161053, 25844927, 28237290, 11890184, 14004553, 60509061, 10430492,
];

/// RefGenomeCoordinates.centromere(String): 0 for a chromosome that HumanChromosome.contains rejects.
pub fn centromere(chromosome: &str, v38: bool) -> Result<i32> {
    if !HumanChromosome::contains(chromosome) {
        return Ok(0);
    }
    let c = HumanChromosome::parse(chromosome)?;
    Ok(if v38 { CENTROMERES_38[c.ordinal()] } else { CENTROMERES_37[c.ordinal()] })
}

/// Arm: P before Q.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Arm {
    P,
    Q,
}

/// ChrArm, ordered by chromosome and then arm.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ChrArm {
    pub chromosome: HumanChromosome,
    pub arm: Arm,
}

/// ChrArmLocator.defaultLocator: the P arm before the centromere, the Q arm from it on.
pub fn default_locator(v38: bool) -> impl Fn(&str, i32) -> Result<ChrArm> + Sync {
    move |chromosome: &str, position: i32| {
        let centre = centromere(chromosome, v38)?;
        let c = HumanChromosome::parse(chromosome)?;
        Ok(ChrArm { chromosome: c, arm: if position < centre { Arm::P } else { Arm::Q } })
    }
}

/// DataForSegmentation: the values that are segmented and those that are averaged, of equal length.
#[derive(Clone, Debug, Default)]
pub struct DataForSegmentation {
    pub values_for_segmentation: Vec<f64>,
    pub raw_values: Vec<f64>,
}

impl DataForSegmentation {
    pub fn new(values_for_segmentation: Vec<f64>, raw_values: Vec<f64>) -> Result<Self> {
        if values_for_segmentation.len() != raw_values.len() {
            bail!("DataForSegmentation: arrays of different lengths");
        }
        Ok(DataForSegmentation { values_for_segmentation, raw_values })
    }
    pub fn count(&self) -> usize {
        self.values_for_segmentation.len()
    }
    pub fn is_empty(&self) -> bool {
        self.count() == 0
    }
}

/// PcfSegment: a segment and the mean of its raw values. `chromosome` is the short name, as Java's.
pub use pcf::PcfSegment;

/// How the values of an arm are segmented: AMBER's absolute positions or COBALT's 1,000-base windows, and
/// the value of each raw value that is segmented (identity for AMBER, the floored log2 for COBALT).
#[derive(Clone, Copy)]
pub struct SegmentationKind {
    pub windowed: bool,
    pub segmentation_value: fn(f64) -> f64,
}

/// AMBER's BAFSegmenter: absolute positions, the values themselves.
pub const ABSOLUTE: SegmentationKind = SegmentationKind { windowed: false, segmentation_value: identity };
/// COBALT's CobaltRatioSegmenter: windows, the ratio floored at 0.001 and taken as a float log2.
pub const WINDOWED: SegmentationKind = SegmentationKind { windowed: true, segmentation_value: cobalt_segmentation_value };

fn identity(v: f64) -> f64 {
    v
}

/// CobaltRatioSegmenter.buildSegmentationData: -9.965784 below 0.001, else `(float) FastMath.log(2, v)`.
pub fn cobalt_segmentation_value(v: f64) -> f64 {
    if v < 0.001 { COBALT_FLOOR_LOG2 } else { crate::java_fastmath::log_base(2.0, v) as f32 as f64 }
}

/// Whether PerArmSegmenter keeps a value: `value >= 0` (SG2: and finite).
fn keep_value(v: f64, keep_infinite: bool) -> bool {
    v >= 0.0 && (v.is_finite() || keep_infinite)
}

/// PerArmSegmenter with the default locator, then getSegmentation: the segments of each arm.
pub fn segment(points: &[SegmentPoint], v38: bool, gamma: f64, kind: SegmentationKind, threads: usize) -> Result<BTreeMap<ChrArm, Vec<PcfSegment>>> {
    segment_with_locator(points, &default_locator(v38), gamma, kind, threads)
}

/// PerArmSegmenter with a locator of one's own (the unit tests put whole chromosomes on one arm).
pub fn segment_with_locator(
    points: &[SegmentPoint],
    locator: &(dyn Fn(&str, i32) -> Result<ChrArm> + Sync),
    gamma: f64,
    kind: SegmentationKind,
    threads: usize,
) -> Result<BTreeMap<ChrArm, Vec<PcfSegment>>> {
    // ArmToRatios: the points of each arm in input order (an arm lies on one chromosome, and Java keeps the
    // order of each chromosome's list)
    let mut arms: BTreeMap<ChrArm, (Vec<i32>, DataForSegmentation)> = BTreeMap::new();
    let keep_infinite = java("SG2");
    for p in points {
        if !keep_value(p.value, keep_infinite) {
            continue;
        }
        let arm = locator(&p.chromosome, p.position)?;
        let (positions, data) = arms.entry(arm).or_default();
        positions.push(p.position);
        data.values_for_segmentation.push((kind.segmentation_value)(p.value));
        data.raw_values.push(p.value);
    }
    // the penalty from all values, arms in order
    let all: Vec<f64> = arms.values().flat_map(|(_, d)| d.values_for_segmentation.iter().copied()).collect();
    crate::debug!("using uniform segmentation penalty, number of ratios({})", all.len());
    if all.is_empty() {
        if java("SG1") {
            bail!("Input array must not be empty");
        }
        return Ok(BTreeMap::new());
    }
    let penalty = copynumber::segment_penalty(&all, gamma, true)?;
    crate::debug!("uniform segmentation penalty({})", crate::fmt::java_fixed(penalty, 3));
    drop(all);

    // the arms, largest first, over the threads
    let work: Vec<(ChrArm, &(Vec<i32>, DataForSegmentation))> = {
        let mut w: Vec<_> = arms.iter().map(|(a, d)| (*a, d)).collect();
        w.sort_by(|a, b| b.1 .0.len().cmp(&a.1 .0.len()).then(a.0.cmp(&b.0)));
        w
    };
    let threads = threads.max(1).min(work.len());
    let next = std::sync::atomic::AtomicUsize::new(0);
    let results: std::sync::Mutex<Vec<(ChrArm, Result<Vec<PcfSegment>>)>> = std::sync::Mutex::new(Vec::new());
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let Some((arm, (positions, data))) = work.get(i) else { break };
                let r = arm_segments(*arm, positions, data, penalty, kind.windowed);
                results.lock().unwrap().push((*arm, r));
            });
        }
    });
    // every arm is done, so the error reported (the first in arm order) does not depend on the threads
    let mut out = BTreeMap::new();
    let mut results = results.into_inner().unwrap();
    results.sort_by_key(|r| r.0);
    for (arm, r) in results {
        out.insert(arm, r?);
    }
    Ok(out)
}

/// calculateSegmentsForChromosome and ChromosomeArmSegments: fit the arm, then each segment's first and last
/// position and the mean of its raw values.
fn arm_segments(arm: ChrArm, positions: &[i32], data: &DataForSegmentation, penalty: f64, windowed: bool) -> Result<Vec<PcfSegment>> {
    let lengths = copynumber::least_cost_segment_lengths(&data.values_for_segmentation, penalty)
        .with_context(|| format!("segmenting chromosome {} arm {:?}", arm.chromosome.short_name(), arm.arm))?;
    let mut segments = Vec::with_capacity(lengths.len());
    let mut start = 0;
    for len in lengths {
        let end_position = positions[start + len - 1];
        segments.push(PcfSegment {
            chromosome: arm.chromosome.short_name().to_string(),
            start: positions[start],
            end: if windowed { end_position + WINDOW_SIZE - 1 } else { end_position },
            mean_ratio: copynumber::mean(&data.raw_values[start..start + len]),
        });
        start += len;
    }
    Ok(segments)
}

/// SegmentsFile.write: a header, then the segments of each arm in arm order, the chromosome named for the
/// genome version and the mean written as DecimalFormat("#.####"). (Java gzips a name ending in .gz; that is
/// not done here, and no AMBER or COBALT file name does.)
pub fn write_segments_file(path: &str, v38: bool, segments: &BTreeMap<ChrArm, Vec<PcfSegment>>) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("failed to create {path}"))?;
    let mut w = std::io::BufWriter::new(file);
    w.write_all(b"Chromosome\tStart\tEnd\tMeanRatio\n")?;
    for (arm, list) in segments {
        let chromosome = arm.chromosome.versioned(v38);
        for s in list {
            writeln!(w, "{chromosome}\t{}\t{}\t{}", s.start, s.end, crate::fmt::decimal_format(s.mean_ratio, 4))?;
        }
    }
    w.flush()?;
    Ok(())
}

/// PerArmSegmenter for values that are not windowed (AbsoluteSegments), as AMBER's BAFSegmenter, then
/// SegmentsFile.write to `path`. `points` are in genome order; Java drops values below 0. With no points
/// the file holds the header only.
pub fn write_absolute_segments(points: &[SegmentPoint], v38: bool, gamma: f64, threads: usize, path: &str) -> anyhow::Result<()> {
    if points.is_empty() {
        return write_segments_file(path, v38, &BTreeMap::new());
    }
    // BAFSegmenter keys every BAF by HumanChromosome first, whatever its value
    for p in points {
        HumanChromosome::parse(&p.chromosome)?;
    }
    let segments = segment(points, v38, gamma, ABSOLUTE, threads)?;
    write_segments_file(path, v38, &segments)
}

/// CobaltRatioSegmenter.writeTumorSegments / writeReferenceSegments: the ratios (tumorGCRatio or
/// referenceGCDiploidRatio) of the windows starting at each position, segmented as floored float log2
/// values, averaged as ratios, written as windows (WindowSegments). COBALT's gamma is `-pcf_gamma`
/// (default 100). Java has no case for an empty list (SG1).
pub fn write_windowed_segments(points: &[SegmentPoint], v38: bool, gamma: f64, threads: usize, path: &str) -> anyhow::Result<()> {
    let segments = segment(points, v38, gamma, WINDOWED, threads)?;
    write_segments_file(path, v38, &segments)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// hmf-common's test resources (the Java tests' segmentation and pcf files), next to this crate in the
    /// hmftools repository, or where HMF_COMMON_TEST_RESOURCES points.
    pub(crate) fn resources_dir() -> Option<std::path::PathBuf> {
        let candidates = [
            std::env::var("HMF_COMMON_TEST_RESOURCES").ok().map(std::path::PathBuf::from),
            Some(std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../hmf-common/src/test/resources")),
        ];
        let found = candidates.into_iter().flatten().find(|d| d.join("segmentation").is_dir());
        if found.is_none() {
            eprintln!("hmf-common test resources not found: test skipped");
        }
        found
    }

    fn arm(c: &str, a: Arm) -> ChrArm {
        ChrArm { chromosome: HumanChromosome::parse(c).unwrap(), arm: a }
    }

    // ---- ChrArmLocatorTest
    #[test]
    fn default_locator_test() {
        let locate = default_locator(true);
        assert_eq!(arm("14", Arm::P), locate("chr14", 10_000_000).unwrap());
        assert_eq!(arm("14", Arm::P), locate("chr14", 17086761).unwrap());
        assert_eq!(arm("14", Arm::Q), locate("chr14", 17086762).unwrap());
        let locate37 = default_locator(false);
        assert_eq!(arm("1", Arm::P), locate37("1", 123035433).unwrap());
        assert_eq!(arm("1", Arm::Q), locate37("1", 123035434).unwrap());
        assert_eq!(arm("Y", Arm::Q), locate37("Y", 11604553).unwrap());
        // a v38 run given a name without the prefix still finds the arm
        assert_eq!(arm("X", Arm::P), locate("X", 60509060).unwrap());
        assert!(locate("chrM", 100).is_err() && locate("chr1_KI270706v1_random", 100).is_err());
        // "CHR1" passes fromString but not contains: no centromere, so always Q
        assert_eq!(arm("1", Arm::Q), locate("CHR1", 100).unwrap());
    }

    // ---- ChrArmTest
    #[test]
    fn chr_arm_compare() {
        use std::cmp::Ordering::*;
        assert_eq!(Equal, arm("1", Arm::P).cmp(&arm("1", Arm::P)));
        assert_eq!(Less, arm("1", Arm::P).cmp(&arm("1", Arm::Q)));
        assert_eq!(Greater, arm("2", Arm::P).cmp(&arm("1", Arm::Q)));
        assert_eq!(Greater, arm("2", Arm::P).cmp(&arm("1", Arm::P)));
        assert_eq!(Less, arm("22", Arm::Q).cmp(&arm("X", Arm::P)));
        assert_eq!(Less, arm("9", Arm::Q).cmp(&arm("10", Arm::P)));
    }

    // ---- DataForSegmentationTest
    #[test]
    fn data_for_segmentation() {
        let d = DataForSegmentation::new(vec![1.0, 2.0, 3.0], vec![1.1, 2.1, 3.1]).unwrap();
        assert_eq!(3, d.count());
        assert!(!d.is_empty());
        assert!(DataForSegmentation::new(vec![], vec![]).unwrap().is_empty());
        assert!(DataForSegmentation::new(vec![1.0], vec![]).is_err());
    }

    #[test]
    fn chromosome_names() {
        assert_eq!(HumanChromosome::from_name("chr1"), HumanChromosome::from_name("1"));
        assert_eq!(HumanChromosome::from_name("Chr22").map(|c| c.ordinal()), Some(21));
        assert!(HumanChromosome::from_name("x").is_none() && HumanChromosome::from_name("MT").is_none());
        assert!(HumanChromosome::from_name("01").is_none() && HumanChromosome::contains("01"));
        assert!(HumanChromosome::contains("chrY") && !HumanChromosome::contains("chr23") && !HumanChromosome::contains("CHR1"));
        assert_eq!(HumanChromosome::parse("X").unwrap().versioned(true), "chrX");
        assert_eq!(HumanChromosome::parse("chrX").unwrap().versioned(false), "X");
    }

    // ---- PerArmSegmenterTest: the input order of the chromosomes does not change the result
    #[test]
    fn segmentation_is_deterministic() {
        let locator = |c: &str, _p: i32| -> Result<ChrArm> { Ok(ChrArm { chromosome: HumanChromosome::parse(c)?, arm: Arm::P }) };
        let block_size = 20;
        let shift = 300.0;
        let mut blocks: Vec<Vec<SegmentPoint>> = Vec::new();
        for hc in HumanChromosome::ALL {
            let multiplier = if hc.ordinal() % 2 == 0 { 1.0 } else { 100.0 };
            let value = (hc.ordinal() + 1) as f64 * multiplier;
            let mut position = 1000;
            let mut block = Vec::new();
            for i in 0..block_size {
                let v = value + if hc.ordinal() == 0 && i >= block_size / 2 { shift } else { 0.0 };
                block.push(SegmentPoint { chromosome: hc.short_name().to_string(), position, value: v });
                position += 1000;
            }
            blocks.push(block);
        }
        let identity_windowed = SegmentationKind { windowed: true, segmentation_value: identity };
        let first = segment_with_locator(&blocks.concat(), &locator, 100.0, identity_windowed, 8).unwrap();
        assert_eq!(first.len(), 24);
        let mut state = 7u64;
        for _ in 0..100 {
            // shuffle the chromosome blocks
            for i in (1..blocks.len()).rev() {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                blocks.swap(i, (state >> 33) as usize % (i + 1));
            }
            let threads = 1 + (state >> 60) as usize % 4;
            let again = segment_with_locator(&blocks.concat(), &locator, 100.0, identity_windowed, threads).unwrap();
            assert_eq!(first, again);
        }
    }

    /// The Java setup of PerArmSegmenterTest, three points on each of 1, 2 and X, on one arm each
    #[test]
    fn small_arms_and_windowed_ends() {
        let mut points = Vec::new();
        for c in ["1", "2", "X"] {
            for (pos, v) in [(1, 0.5), (1001, 0.51), (2001, 0.49)] {
                points.push(SegmentPoint { chromosome: c.into(), position: pos, value: v });
            }
        }
        let locator = |c: &str, _p: i32| -> Result<ChrArm> { Ok(ChrArm { chromosome: HumanChromosome::parse(c)?, arm: Arm::P }) };
        let kind = SegmentationKind { windowed: true, segmentation_value: identity };
        let result = segment_with_locator(&points, &locator, 100.0, kind, 2).unwrap();
        assert_eq!(result.len(), 3);
        for segments in result.values() {
            assert_eq!(segments.len(), 1);
            assert_eq!((segments[0].start, segments[0].end), (1, 2001 + WINDOW_SIZE - 1));
            assert!((segments[0].mean_ratio - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn negative_and_infinite_values() {
        let pts = |vals: &[f64]| -> Vec<SegmentPoint> {
            vals.iter().enumerate().map(|(i, v)| SegmentPoint { chromosome: "chr3".into(), position: 1000 + i as i32, value: *v }).collect()
        };
        // nothing left: header only (SG1 corrected)
        let r = segment(&pts(&[-1.0, f64::NAN, -0.5]), true, 100.0, ABSOLUTE, 1).unwrap();
        assert!(r.is_empty());
        // an infinite value is dropped (SG2 corrected)
        let r = segment(&pts(&[f64::INFINITY, 0.5, 0.5, 0.5]), true, 100.0, ABSOLUTE, 1).unwrap();
        let s = &r[&arm("3", Arm::P)];
        assert_eq!((s.len(), s[0].start, s[0].end), (1, 1001, 1003));
    }

    #[test]
    fn cobalt_values() {
        assert_eq!(cobalt_segmentation_value(0.0009), COBALT_FLOOR_LOG2);
        assert_eq!(cobalt_segmentation_value(1.0), 0.0);
        assert_eq!(cobalt_segmentation_value(2.0), 1.0);
        assert_eq!(cobalt_segmentation_value(3.0), 1.5849625007211563f64 as f32 as f64);
    }

    #[test]
    fn file_layout() {
        let dir = std::env::temp_dir().join(format!("hmf_common_seg_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.pcf");
        let path = path.to_str().unwrap();
        write_absolute_segments(&[], true, 100.0, 2, path).unwrap();
        assert_eq!(std::fs::read_to_string(path).unwrap(), "Chromosome\tStart\tEnd\tMeanRatio\n");
        let mut points = Vec::new();
        for i in 0..40 {
            let v = if i < 20 { 0.52 + (i % 3) as f64 * 0.01 } else { 0.95 + (i % 2) as f64 * 0.01 };
            points.push(SegmentPoint { chromosome: "chr2".into(), position: 100_000_000 + i * 1000, value: v });
        }
        write_absolute_segments(&points, true, 100.0, 2, path).unwrap();
        let text = std::fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines, ["Chromosome\tStart\tEnd\tMeanRatio", "chr2\t100000000\t100019000\t0.5295", "chr2\t100020000\t100039000\t0.955"]);
        write_absolute_segments(&points, false, 100.0, 1, path).unwrap();
        assert!(std::fs::read_to_string(path).unwrap().contains("\n2\t100000000\t"));
        let bad = [SegmentPoint { chromosome: "chrM".into(), position: 1, value: -1.0 }];
        assert!(write_absolute_segments(&bad, true, 100.0, 1, path).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
