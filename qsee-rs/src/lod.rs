//! qsee `lod`: LimitOfDetectionCalcs (missed variant likelihood and its inverse), GeneLodCalculator and
//! GeneLodFile (`<tumor>.qsee.gene_lod.tsv.gz`).

use anyhow::{anyhow, Result};

use hmf_common::commons_math::binomial_cdf;
use hmf_common::delim_file::write_delim_file;
use hmf_common::metrics::{depth_range_field_name, DepthRange, ExonDepth, GeneDepthFile};
use hmf_common::purple_files::{purple_purity_filename, purple_purity_from_lines};

use crate::common::{decimal_8, generate_filename, SampleType, COL_SAMPLE_ID, COL_SAMPLE_TYPE, MULTISAMPLE_SAMPLE_ID};
use crate::config::QseeConfig;
use crate::prep::{depth_counts_map, read_input, reportable_genes, PrepError};

const MIN_ALT_READS: i32 = 3;
const LOD_VAF_BISECTION_ITERATIONS: usize = 20;
const LOD_TARGET_MISSED_VARIANT_LIKELIHOOD: f64 = 0.05;
pub const REF_SAMPLE_PURITY: f64 = 1.0;

const TUMOR_VAF_CLONAL_HOM_VARIANT_DETECTABLE: f64 = 1.0;
const TUMOR_VAF_CLONAL_HET_VARIANT_DETECTABLE: f64 = 0.5;

/// LimitOfDetectionCalcs.calcMissedVariantLikelihood: the probability of at most two alt reads at each
/// depth range's midpoint, weighted by the range's base count (ranges in file order)
pub fn missed_variant_likelihood(depth_counts: &[(DepthRange, i32)], vaf: f64) -> f64 {
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
            at_depth = binomial_cdf(range.midpoint(), vaf, MIN_ALT_READS - 1);
        }
        mean += at_depth * weight;
    }
    mean
}

/// LimitOfDetectionCalcs.calcMinSampleVafToDetectVariant: the VAF at which the missed variant
/// likelihood is 5%, by 20 bisections; 1 when even a VAF of 1 misses more
pub fn min_sample_vaf_to_detect_variant(depth_counts: &[(DepthRange, i32)]) -> f64 {
    if missed_variant_likelihood(depth_counts, 1.0) > LOD_TARGET_MISSED_VARIANT_LIKELIHOOD {
        return 1.0;
    }
    let (mut lower, mut upper, mut mid) = (0.0f64, 1.0f64, 0.5f64);
    for _ in 0..LOD_VAF_BISECTION_ITERATIONS {
        if missed_variant_likelihood(depth_counts, mid) <= LOD_TARGET_MISSED_VARIANT_LIKELIHOOD {
            upper = mid;
        } else {
            lower = mid;
        }
        mid = (lower + upper) / 2.0;
    }
    mid
}

/// VariantDetectabilityStatus
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantDetectabilityStatus {
    Pass,
    WarnHetUnreliable,
    Fail,
}

impl VariantDetectabilityStatus {
    pub fn name(self) -> &'static str {
        match self {
            VariantDetectabilityStatus::Pass => "PASS",
            VariantDetectabilityStatus::WarnHetUnreliable => "WARN_HET_UNRELIABLE",
            VariantDetectabilityStatus::Fail => "FAIL",
        }
    }
    fn from_tumor_vaf_lod(lod: f64) -> VariantDetectabilityStatus {
        if lod >= TUMOR_VAF_CLONAL_HOM_VARIANT_DETECTABLE {
            VariantDetectabilityStatus::Fail
        } else if lod > TUMOR_VAF_CLONAL_HET_VARIANT_DETECTABLE {
            VariantDetectabilityStatus::WarnHetUnreliable
        } else {
            VariantDetectabilityStatus::Pass
        }
    }
}

/// GeneLod
#[derive(Clone, Debug)]
pub struct GeneLod {
    pub sample_id: String,
    pub sample_type: SampleType,
    pub gene_name: String,
    pub exon_rank: i32,
    pub median_depth: f64,
    pub sample_vaf_lod: f64,
    pub tumor_vaf_lod: f64,
    pub status: VariantDetectabilityStatus,
    pub depth_counts: Vec<(DepthRange, i32)>,
}

/// GeneLodCalculator.calcGeneLods: the exons of reportable genes, in file order
pub fn calc_gene_lods(sample_id: &str, sample_type: SampleType, exons: &[ExonDepth], reportable: &[String], tumor_purity: f64) -> Vec<GeneLod> {
    exons
        .iter()
        .filter(|e| reportable.contains(&e.depth.gene))
        .map(|e| {
            let counts = depth_counts_map(&e.depth.depth_counts);
            let sample_vaf_lod = min_sample_vaf_to_detect_variant(&counts);
            let tumor_vaf_lod = sample_vaf_lod / tumor_purity;
            GeneLod {
                sample_id: sample_id.to_string(),
                sample_type,
                gene_name: e.depth.gene.clone(),
                exon_rank: e.exon_rank,
                median_depth: e.median_depth,
                sample_vaf_lod,
                tumor_vaf_lod,
                status: VariantDetectabilityStatus::from_tumor_vaf_lod(tumor_vaf_lod),
                depth_counts: counts,
            }
        })
        .collect()
}

/// calcSampleGeneLods: an IOException (a missing input) is `PrepError::Missing`; anything else ends
/// the run, as Java's uncaught exceptions do
fn calc_sample_gene_lods(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> std::result::Result<Vec<GeneLod>, PrepError> {
    let path = GeneDepthFile::generate_exon_coverage_filename(&config.bam_metrics_dir(sample_id, sample_type)?, sample_id);
    read_input(&path)?;
    let exons = GeneDepthFile::read_exon_coverage(&path)?;
    let reportable = reportable_genes(config, sample_type);
    let purity = if sample_type == SampleType::Tumor {
        let p = purple_purity_filename(&config.purple_dir(sample_id)?, sample_id);
        let lines = read_input(&p)?;
        purple_purity_from_lines(&lines)?.purity
    } else {
        REF_SAMPLE_PURITY
    };
    Ok(calc_gene_lods(sample_id, sample_type, &exons, &reportable, purity))
}

/// GeneLodFile.generateFilename
pub fn gene_lod_filename(config: &QseeConfig) -> String {
    let sample = if config.is_single_patient() { config.sample_ids(SampleType::Tumor)[0].as_str() } else { MULTISAMPLE_SAMPLE_ID };
    generate_filename(&config.output_dir, sample, "gene_lod", config.output_id.as_deref(), "tsv.gz")
}

/// GeneLodCalculator.run: every tumor then every normal sample; the file is written even when every
/// input is missing (with -allow_missing_input)
pub fn run(config: &QseeConfig) -> Result<()> {
    let mut lods = Vec::new();
    for sample_type in SampleType::ALL {
        for sample_id in config.sample_ids(sample_type) {
            match calc_sample_gene_lods(config, sample_id, sample_type) {
                Ok(l) => lods.extend(l),
                Err(PrepError::Missing(m)) => {
                    if config.allow_missing_input {
                        hmf_common::warn!("sampleType({}) gene LoD - missing input file(s): {m}", sample_type.name());
                    } else {
                        hmf_common::error!("sampleType({}) gene LoD - missing input file(s): {m}", sample_type.name());
                        return Err(anyhow!("missing input"));
                    }
                }
                Err(PrepError::Fatal(e)) => return Err(e),
                Err(PrepError::Exit) => return Err(anyhow!("exit")),
            }
        }
    }
    if !config.is_single_patient() {
        lods.sort_by(|a, b| hmf_common::java_text::compare(&a.sample_id, &b.sample_id));
    }
    write_gene_lods(&gene_lod_filename(config), &lods, config.write_lod_details)
}

/// GeneLodFile.write: with -write_lod_details, a column per depth range of the first exon
pub fn write_gene_lods(path: &str, lods: &[GeneLod], details: bool) -> Result<()> {
    hmf_common::info!("Writing gene limit of detection to: {path}");
    let ranges: Vec<DepthRange> = if details && !lods.is_empty() { lods[0].depth_counts.iter().map(|(r, _)| *r).collect() } else { Vec::new() };
    let mut columns: Vec<String> =
        [COL_SAMPLE_ID, COL_SAMPLE_TYPE, "GeneName", "ExonRank", "MedianDepth", "SampleVafLod", "TumorVafLod", "VariantDetectabilityStatus"].iter().map(|s| s.to_string()).collect();
    columns.extend(ranges.iter().map(depth_range_field_name));
    write_delim_file(path, &columns, lods, |g, row| {
        row.set(COL_SAMPLE_ID, g.sample_id.as_str())?;
        row.set(COL_SAMPLE_TYPE, g.sample_type.name())?;
        row.set("GeneName", g.gene_name.as_str())?;
        row.set("ExonRank", g.exon_rank.to_string())?;
        row.set("MedianDepth", decimal_8(g.median_depth))?;
        row.set("SampleVafLod", decimal_8(g.sample_vaf_lod))?;
        row.set("TumorVafLod", decimal_8(g.tumor_vaf_lod))?;
        row.set("VariantDetectabilityStatus", g.status.name())?;
        for r in &ranges {
            // a range the exon does not have: Java unboxes a null Integer (NullPointerException)
            let c = g.depth_counts.iter().find(|(x, _)| x == r).map(|(_, c)| *c).ok_or_else(|| anyhow!("java.lang.NullPointerException: no count for {}", depth_range_field_name(r)))?;
            row.set(&depth_range_field_name(r), c.to_string())?;
        }
        Ok(())
    })
}
