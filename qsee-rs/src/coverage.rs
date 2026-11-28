//! qsee `coverage` (upstream bdceb2302a; `lod` before): LimitOfDetectionCalcs (the missed variant likelihood and
//! its inverse, the LoD VAF), ExonCoverageQc (per exon of the reportable genes: median depth, LoD VAF and a
//! detectability status against the EXON_VARIANT_LOD threshold), ExonCoverageQcFile
//! (`<tumor>.qsee.exon_coverage_qc.tsv.gz`), ExonCoverageQcFileBuilder (the file, and its own entry point),
//! ExonExclusions (`-exon_exclusions_file`) and ExonExclusionsFileBuilder (the exclusions from a cohort, an entry
//! point of its own).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use anyhow::{anyhow, Result};

use hmf_common::commons_math::binomial_cdf;
use hmf_common::delim_file::{write_delim_file, DelimFileReader};
use hmf_common::metrics::{depth_range_field_name, DepthRange, ExonDepth, GeneDepthFile};

use crate::common::{decimal_8, generate_filename, generate_filename_opt, SampleType, COL_SAMPLE_ID, COL_SAMPLE_TYPE, MULTISAMPLE_SAMPLE_ID};
use crate::config::QseeConfig;
use crate::feature::FeatureType;
use crate::prep::{depth_counts_map, read_input, reportable_genes_of, PrepError};
use crate::status::{QcStatusType, ThresholdKey, ThresholdRegistry};

// ---- LimitOfDetectionCalcs ----

pub const DEFAULT_LOD_MIN_ALT_READS: i32 = 3;
/// the LoD VAF's floor, an approximate noise floor below which variants are more likely artefacts
pub const LOD_VAF_LOWER_BOUND: f64 = 0.01;
const LOD_VAF_BISECTION_ITERATIONS: usize = 20;
const LOD_TARGET_MISSED_VARIANT_LIKELIHOOD: f64 = 0.05;
pub const REF_SAMPLE_PURITY: f64 = 1.0;

/// LimitOfDetectionCalcs.calcMissedVariantLikelihood: the probability of fewer than `min_alt_reads` alt reads at
/// each depth range's midpoint, weighted by the range's base count (ranges in file order)
pub fn missed_variant_likelihood(depth_counts: &[(DepthRange, i32)], vaf: f64, min_alt_reads: i32) -> f64 {
    let total = depth_counts.iter().fold(0i32, |s, (_, c)| s.wrapping_add(*c));
    if total == 0 {
        return 0.0;
    }
    let mut mean = 0.0;
    for (range, count) in depth_counts {
        let weight = *count as f64 / total as f64;
        if weight == 0.0 {
            continue;
        }
        let mut at_depth = 1.0;
        if range.midpoint() > 0 && vaf > 0.0 {
            at_depth = binomial_cdf(range.midpoint(), vaf, min_alt_reads.wrapping_sub(1));
        }
        mean += at_depth * weight;
    }
    mean
}

/// LimitOfDetectionCalcs.calcMinSampleVafToDetectVariant: the VAF at which the missed variant likelihood is 5%,
/// by 20 bisections, floored at 0.01; 1 when even a VAF of 1 misses more
pub fn min_sample_vaf_to_detect_variant(depth_counts: &[(DepthRange, i32)], min_alt_reads: i32) -> f64 {
    if missed_variant_likelihood(depth_counts, 1.0, min_alt_reads) > LOD_TARGET_MISSED_VARIANT_LIKELIHOOD {
        return 1.0;
    }
    let (mut lower, mut upper, mut mid) = (0.0f64, 1.0f64, 0.5f64);
    for _ in 0..LOD_VAF_BISECTION_ITERATIONS {
        if missed_variant_likelihood(depth_counts, mid, min_alt_reads) <= LOD_TARGET_MISSED_VARIANT_LIKELIHOOD {
            upper = mid;
        } else {
            lower = mid;
        }
        mid = (lower + upper) / 2.0;
    }
    // Math.max
    if mid.is_nan() { mid } else { mid.max(LOD_VAF_LOWER_BOUND) }
}

// ---- ExonExclusions ----

/// ExonExclusions: exons by gene name
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ExonExclusions {
    exclusions: HashMap<String, BTreeSet<i32>>,
}

impl ExonExclusions {
    pub fn empty() -> ExonExclusions {
        ExonExclusions::default()
    }
    pub fn add(&mut self, gene: &str, exon_rank: i32) {
        self.exclusions.entry(gene.to_string()).or_default().insert(exon_rank);
    }
    pub fn is_excluded(&self, gene: &str, exon_rank: i32) -> bool {
        self.exclusions.get(gene).is_some_and(|s| s.contains(&exon_rank))
    }

    /// ExonExclusions.read: columns GeneName and ExonRank
    pub fn read(path: &str) -> Result<ExonExclusions> {
        hmf_common::info!("Reading exon exclusions from: {path}");
        let reader = DelimFileReader::open(path)?;
        let mut out = ExonExclusions::empty();
        for row in reader.rows() {
            let gene = row.get("GeneName")?;
            let rank = row.get_int("ExonRank")?;
            out.add(gene, rank);
        }
        Ok(out)
    }

    /// ExonExclusions.write: by gene name, then exon rank
    pub fn write(&self, path: &str) -> Result<()> {
        hmf_common::info!("Writing exon exclusions to: {path}");
        let sorted: BTreeMap<&String, &BTreeSet<i32>> = self.exclusions.iter().collect();
        let mut rows: Vec<(String, i32)> = Vec::new();
        for (gene, ranks) in sorted {
            for r in ranks.iter() {
                rows.push((gene.clone(), *r));
            }
        }
        let columns = vec!["GeneName".to_string(), "ExonRank".to_string()];
        write_delim_file(path, &columns, &rows, |(g, r), row| {
            row.set("GeneName", g.as_str())?;
            row.set("ExonRank", r.to_string())
        })
    }
}

/// ExonExclusions.generateFilename: `qsee.exon_exclusions[.<outputId>].tsv`
pub fn exon_exclusions_filename(output_dir: &str, output_id: Option<&str>) -> String {
    generate_filename_opt(output_dir, None, "exon_exclusions", output_id, "tsv")
}

// ---- ExonCoverageQc ----

/// VariantDetectabilityStatus
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantDetectabilityStatus {
    Pass,
    Warn,
    Excluded,
}

impl VariantDetectabilityStatus {
    pub fn name(self) -> &'static str {
        match self {
            VariantDetectabilityStatus::Pass => "PASS",
            VariantDetectabilityStatus::Warn => "WARN",
            VariantDetectabilityStatus::Excluded => "EXCLUDED",
        }
    }
}

/// ExonCoverageQc.LOD_VAF: the feature name of the EXON_VARIANT_LOD threshold
pub const LOD_VAF: &str = "LOD_VAF";

/// ExonCoverageQc
#[derive(Clone, Debug)]
pub struct ExonCoverageQc {
    pub sample_id: String,
    pub sample_type: SampleType,
    pub gene_name: String,
    pub exon_rank: i32,
    pub median_depth: f64,
    pub lod_vaf: f64,
    pub status: VariantDetectabilityStatus,
    pub depth_counts: Vec<(DepthRange, i32)>,
}

/// ExonCoverageQc.determineVariantDetectabilityStatus: EXCLUDED, else WARN when the LoD VAF meets the
/// EXON_VARIANT_LOD LOD_VAF threshold of the sample type, else PASS
fn detectability_status(sample_type: SampleType, gene: &str, exon_rank: i32, lod_vaf: f64, exclusions: &ExonExclusions, thresholds: &ThresholdRegistry) -> Result<VariantDetectabilityStatus> {
    if exclusions.is_excluded(gene, exon_rank) {
        return Ok(VariantDetectabilityStatus::Excluded);
    }
    let threshold = thresholds.get(&ThresholdKey::new(sample_type, FeatureType::ExonVariantLod, LOD_VAF, QcStatusType::Warn))?;
    Ok(if threshold.sample_qc_status(lod_vaf).status_type == QcStatusType::Warn { VariantDetectabilityStatus::Warn } else { VariantDetectabilityStatus::Pass })
}

/// ExonCoverageQc.fromExonDepths: the exons of reportable genes, in file order
pub fn coverage_qc_per_exon(
    sample_id: &str, sample_type: SampleType, exons: &[ExonDepth], exclusions: &ExonExclusions, reportable: &[String], thresholds: &ThresholdRegistry, min_alt_reads: i32,
) -> Result<Vec<ExonCoverageQc>> {
    let mut out = Vec::new();
    for e in exons.iter().filter(|e| reportable.contains(&e.depth.gene)) {
        let counts = depth_counts_map(&e.depth.depth_counts);
        let lod_vaf = min_sample_vaf_to_detect_variant(&counts, min_alt_reads);
        let status = detectability_status(sample_type, &e.depth.gene, e.exon_rank, lod_vaf, exclusions, thresholds)?;
        out.push(ExonCoverageQc {
            sample_id: sample_id.to_string(),
            sample_type,
            gene_name: e.depth.gene.clone(),
            exon_rank: e.exon_rank,
            median_depth: e.median_depth,
            lod_vaf,
            status,
            depth_counts: counts,
        });
    }
    Ok(out)
}

/// ExonCoverageQc.calcLowCoverageExonProportion: WARN exons over WARN and PASS exons (0 without any)
pub fn low_coverage_exon_proportion(exons: &[ExonCoverageQc]) -> f64 {
    let warn = exons.iter().filter(|x| x.status == VariantDetectabilityStatus::Warn).count() as i64;
    let pass = exons.iter().filter(|x| x.status == VariantDetectabilityStatus::Pass).count() as i64;
    let assessed = warn + pass;
    if assessed > 0 { warn as f64 / assessed as f64 } else { 0.0 }
}

// ---- ExonCoverageQcFile, ExonCoverageQcFileBuilder ----

/// ExonCoverageQcFile.generateFilename
pub fn exon_coverage_qc_filename(config: &QseeConfig) -> String {
    let sample = if config.is_single_patient() { config.sample_ids(SampleType::Tumor)[0].as_str() } else { MULTISAMPLE_SAMPLE_ID };
    generate_filename(&config.output_dir, sample, "exon_coverage_qc", config.output_id.as_deref(), "tsv.gz")
}

/// ExonCoverageQcFileBuilder.loadExonCoverage then ExonCoverageQc.fromExonDepths: an IOException (a missing
/// input) is `PrepError::Missing`; anything else ends the run, as Java's uncaught exceptions do
fn sample_coverage_qc(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> std::result::Result<Vec<ExonCoverageQc>, PrepError> {
    let path = GeneDepthFile::generate_exon_coverage_filename(&config.bam_metrics_dir(sample_id, sample_type)?, sample_id);
    read_input(&path)?;
    let exons = GeneDepthFile::read_exon_coverage(&path)?;
    let reportable = reportable_genes_of(&config.driver_genes, sample_type);
    Ok(coverage_qc_per_exon(sample_id, sample_type, &exons, &config.exon_exclusions, &reportable, &config.thresholds, config.lod_min_alt_reads)?)
}

/// ExonCoverageQcFileBuilder.run: every tumor then every normal sample; the file is written even when every
/// input is missing (with -allow_missing_input). Returns the exons, for ExonExclusionsFileBuilder.
pub fn run(config: &QseeConfig) -> Result<Vec<ExonCoverageQc>> {
    let mut exons = Vec::new();
    for sample_type in SampleType::ALL {
        for sample_id in config.sample_ids(sample_type) {
            match sample_coverage_qc(config, sample_id, sample_type) {
                Ok(l) => exons.extend(l),
                Err(PrepError::Missing(m)) => {
                    if config.allow_missing_input {
                        hmf_common::warn!("sampleType({}) exon variant LoD - missing input file(s): {m}", sample_type.name());
                    } else {
                        hmf_common::error!("sampleType({}) exon variant LoD - missing input file(s): {m}", sample_type.name());
                        return Err(anyhow::Error::new(crate::app::Exit));
                    }
                }
                Err(PrepError::Fatal(e)) => return Err(e),
                Err(PrepError::Exit) => return Err(anyhow::Error::new(crate::app::Exit)),
            }
        }
    }
    if !config.is_single_patient() {
        exons.sort_by(|a, b| hmf_common::java_text::compare(&a.sample_id, &b.sample_id));
    }
    write_exon_coverage_qc(&exon_coverage_qc_filename(config), &exons, config.write_lod_details)?;
    Ok(exons)
}

/// ExonCoverageQcFile.write: with -write_lod_details, a column per depth range of the first exon
pub fn write_exon_coverage_qc(path: &str, exons: &[ExonCoverageQc], details: bool) -> Result<()> {
    hmf_common::info!("Writing exon coverage QC to: {path}");
    let ranges: Vec<DepthRange> = if details && !exons.is_empty() { exons[0].depth_counts.iter().map(|(r, _)| *r).collect() } else { Vec::new() };
    let mut columns: Vec<String> =
        [COL_SAMPLE_ID, COL_SAMPLE_TYPE, "GeneName", "ExonRank", "MedianDepth", "LodVaf", "VariantDetectabilityStatus"].iter().map(|s| s.to_string()).collect();
    columns.extend(ranges.iter().map(depth_range_field_name));
    write_delim_file(path, &columns, exons, |g, row| {
        row.set(COL_SAMPLE_ID, g.sample_id.as_str())?;
        row.set(COL_SAMPLE_TYPE, g.sample_type.name())?;
        row.set("GeneName", g.gene_name.as_str())?;
        row.set("ExonRank", g.exon_rank.to_string())?;
        row.set("MedianDepth", decimal_8(g.median_depth))?;
        row.set("LodVaf", decimal_8(g.lod_vaf))?;
        row.set("VariantDetectabilityStatus", g.status.name())?;
        for r in &ranges {
            // a range the exon does not have: Java unboxes a null Integer (NullPointerException)
            let c = g.depth_counts.iter().find(|(x, _)| x == r).map(|(_, c)| *c).ok_or_else(|| anyhow!("java.lang.NullPointerException: no count for {}", depth_range_field_name(r)))?;
            row.set(&depth_range_field_name(r), c.to_string())?;
        }
        Ok(())
    })
}

// ---- ExonExclusionsFileBuilder ----

pub const DEFAULT_MIN_LOW_COVERAGE_LOD_VAF: f64 = 0.07;
pub const DEFAULT_MIN_LOW_COVERAGE_COHORT_FRACTION: f64 = 0.3;

/// ExonExclusionsFileBuilder.calcExonExclusions: an exon is excluded when at least ceil(n x fraction) of the n
/// tumor samples with it have a LoD VAF at or above `min_lod_vaf`
pub fn calc_exon_exclusions(exons: &[ExonCoverageQc], min_lod_vaf: f64, min_cohort_fraction: f64) -> ExonExclusions {
    let mut by_exon: HashMap<(&str, i32), Vec<&ExonCoverageQc>> = HashMap::new();
    for e in exons.iter().filter(|e| e.sample_type == SampleType::Tumor) {
        by_exon.entry((e.gene_name.as_str(), e.exon_rank)).or_default().push(e);
    }
    let mut out = ExonExclusions::empty();
    for ((gene, rank), samples) in by_exon {
        let low = samples.iter().filter(|x| x.lod_vaf >= min_lod_vaf).count() as i64;
        // (int) Math.ceil(size * fraction): Rust's cast saturates and takes NaN to 0, as Java's
        let min_count = (samples.len() as f64 * min_cohort_fraction).ceil() as i32 as i64;
        if low >= min_count {
            out.add(gene, rank);
        }
    }
    out
}

/// ExonExclusionsFileBuilder.run: needs more than one tumor sample; writes the cohort's exon coverage QC file,
/// then `qsee.exon_exclusions.tsv`
pub fn run_exclusions_builder(config: &QseeConfig, min_lod_vaf: f64, min_cohort_fraction: f64) -> Result<()> {
    if config.sample_ids(SampleType::Tumor).len() <= 1 {
        hmf_common::error!("Building exon exclusions requires a cohort of tumor samples");
        return Err(anyhow::Error::new(crate::app::Exit));
    }
    let exons = run(config)?;
    let exclusions = calc_exon_exclusions(&exons, min_lod_vaf, min_cohort_fraction);
    exclusions.write(&exon_exclusions_filename(&config.output_dir, config.output_id.as_deref()))
}
