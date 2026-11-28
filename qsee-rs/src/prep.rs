//! qsee `feature.category`: the ten categories of features, each from one tool's files, and
//! CategoryPrepTask's handling of missing inputs.
//!
//! Java tells a missing input (an IOException: no such file, unreadable, not UTF-8) from any other
//! failure (an input it cannot parse), which ends the run. `PrepError` keeps the two apart: the inputs
//! are read through `read_input`, whose errors are missing inputs, and parsed after, whose errors are
//! fatal.

use anyhow::{anyhow, Result};

use hmf_common::driver_panel::DriverGene;
use hmf_common::metrics::{BamFlagStats, BamMetricCoverage, BamMetricFragmentLength, BamMetricSummary, DepthRange, ExonDepth, GeneDepthFile, ValueFrequency};
use hmf_common::purple_files::{purple_purity_filename, purple_purity_from_lines, purple_qc_filename, read_purple_qc, PurplePurity, PurpleQc};
use hmf_common::redux_files::{self, BqrRecord, JitterTableRow};
use hmf_common::{cobalt_files, esvee_files, java_strictmath};

use crate::common::{form_multi_field, form_single_field, BinnedFrequencies, SampleType};
use crate::compat;
use crate::config::{PrepCategory, QseeConfig, SequencingType};
use crate::feature::{
    BaseQualBin, DiscordantFragGroup, Feature, FeatureKey, FeatureMetadata, FeatureType, RepeatType, SourceTool, SummaryTableFeature,
};
use crate::coverage::{self, ExonCoverageQc, ExonExclusions};
use crate::status::{purple_to_qsee_status, QcStatus, QcStatusType, PURPLE_QC_STATUSES};

/// How a category's preparation ends when it does not give features
#[derive(Debug)]
pub enum PrepError {
    /// an IOException: the message Java logs (usually the path)
    Missing(String),
    /// any other exception: the run ends
    Fatal(anyhow::Error),
    /// System.exit(1) after the error was logged
    Exit,
}

impl From<anyhow::Error> for PrepError {
    fn from(e: anyhow::Error) -> PrepError {
        PrepError::Fatal(e)
    }
}

pub type PrepResult = std::result::Result<Vec<Feature>, PrepError>;

/// Files.readAllLines as an input check: the lines, or the IOException Java would throw (no such
/// file, a directory, bytes that are not UTF-8) as a missing input
pub fn read_input(path: &str) -> std::result::Result<Vec<String>, PrepError> {
    let bytes = std::fs::read(path).map_err(|e| match e.kind() {
        std::io::ErrorKind::NotFound => PrepError::Missing(path.to_string()),
        _ => PrepError::Missing(format!("{path}: {e}")),
    })?;
    let text = String::from_utf8(bytes).map_err(|_| PrepError::Missing("Input length = 1".to_string()))?;
    Ok(hmf_common::java_text::lines_of(&text))
}

/// CategoryPrepTask.missingInputFilesError: a warning with -allow_missing_input, else an error and exit
pub fn missing_input_files_error(allow: bool, category: PrepCategory, sample_type: SampleType, message: &str) -> std::result::Result<(), PrepError> {
    if allow {
        hmf_common::warn!("sampleType({}) category({}) - missing input file(s): {message}", sample_type.name(), category.name());
        Ok(())
    } else {
        hmf_common::error!("sampleType({}) category({}) - missing input file(s): {message}", sample_type.name(), category.name());
        Err(PrepError::Exit)
    }
}

/// CategoryPrepTask.run: the category's features; a missing input is logged and gives none (or ends the
/// run without -allow_missing_input); `missing` counts the samples with a missing input (training)
pub fn run_category_task(config: &QseeConfig, category: PrepCategory, sample_id: &str, sample_type: SampleType, missing: Option<&std::sync::atomic::AtomicUsize>) -> Result<Vec<Feature>, PrepError> {
    match extract_sample_data(config, category, sample_id, sample_type) {
        Ok(f) => Ok(f),
        Err(PrepError::Missing(msg)) => {
            missing_input_files_error(config.allow_missing_input, category, sample_type, &msg)?;
            if let Some(m) = missing {
                m.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            Ok(Vec::new())
        }
        Err(PrepError::Fatal(e)) => {
            hmf_common::error!("sampleType({}) category({}) - Failed to run prep for sample({}): {e:#}", sample_type.name(), category.name(), sample_id);
            Err(PrepError::Exit)
        }
        Err(PrepError::Exit) => Err(PrepError::Exit),
    }
}

/// CategoryPrep.extractSampleData
pub fn extract_sample_data(config: &QseeConfig, category: PrepCategory, sample_id: &str, sample_type: SampleType) -> PrepResult {
    match category {
        PrepCategory::SummaryTableBamMetrics => summary_table_bam_metrics(config, sample_id, sample_type),
        PrepCategory::SummaryTablePurple => summary_table_purple(config, sample_id, sample_type),
        PrepCategory::CoverageDistribution => coverage_distribution(config, sample_id, sample_type),
        PrepCategory::FragLengthDistribution => frag_length_distribution(config, sample_id, sample_type),
        PrepCategory::MissedGeneVariant => missed_gene_variant(config, sample_id, sample_type),
        PrepCategory::DuplicateFreq => duplicate_freq(config, sample_id, sample_type),
        PrepCategory::GcBias => gc_bias(config, sample_id),
        PrepCategory::DiscordantFragFreq => discordant_frag_freq(config, sample_id, sample_type),
        PrepCategory::BaseQualRecalibration => base_qual_recalibration(config, sample_id, sample_type),
        PrepCategory::MsIndelError => ms_indel_error(config, sample_id, sample_type),
    }
}

// ---- summary table: BamMetrics ----

/// BamMetricsData: the four files (the exon coverage file since upstream bdceb2302a), each that can be read
pub struct BamMetricsData {
    pub summary: Option<BamMetricSummary>,
    pub coverage: Option<Vec<ValueFrequency>>,
    pub flag_stats: Option<BamFlagStats>,
    pub exon_depths: Option<Vec<ExonDepth>>,
    pub missing: Vec<String>,
}

fn load_bam_metrics_data(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> std::result::Result<BamMetricsData, PrepError> {
    let dir = config.bam_metrics_dir(sample_id, sample_type)?;
    let mut missing = Vec::new();
    let summary_path = BamMetricSummary::generate_filename(&dir, sample_id);
    let summary = match read_input(&summary_path) {
        Ok(lines) => Some(BamMetricSummary::from_lines(&lines).map_err(|e| anyhow!("{summary_path}: {e:#}"))?),
        Err(PrepError::Missing(_)) => {
            missing.push(summary_path);
            None
        }
        Err(e) => return Err(e),
    };
    let coverage_path = BamMetricCoverage::generate_filename(&dir, sample_id);
    let coverage = match read_input(&coverage_path) {
        Ok(_) => Some(BamMetricCoverage::read(&coverage_path)?),
        Err(PrepError::Missing(_)) => {
            missing.push(coverage_path);
            None
        }
        Err(e) => return Err(e),
    };
    let flag_path = BamFlagStats::generate_filename(&dir, sample_id);
    let flag_stats = match read_input(&flag_path) {
        Ok(lines) => Some(BamFlagStats::from_lines(&lines).map_err(|e| anyhow!("{flag_path}: {e:#}"))?),
        Err(PrepError::Missing(_)) => {
            missing.push(flag_path);
            None
        }
        Err(e) => return Err(e),
    };
    let exon_path = GeneDepthFile::generate_exon_coverage_filename(&dir, sample_id);
    let exon_depths = match read_input(&exon_path) {
        Ok(_) => Some(GeneDepthFile::read_exon_coverage(&exon_path)?),
        Err(PrepError::Missing(_)) => {
            missing.push(exon_path);
            None
        }
        Err(e) => return Err(e),
    };
    Ok(BamMetricsData { summary, coverage, flag_stats, exon_depths, missing })
}

/// SummaryTableBamMetricsPrep.putFeatures for the coverage file: the proportion of bases at or above
/// each coverage level
pub fn coverage_above_features(coverage: &[ValueFrequency], sample_type: SampleType, thresholds: &crate::status::ThresholdRegistry) -> Result<Vec<(SummaryTableFeature, Feature)>> {
    use SummaryTableFeature::*;
    let mut out = Vec::new();
    let total = coverage.iter().fold(0i64, |s, v| s.wrapping_add(v.count));
    for (f, level) in [(CoverageAbove10, 10), (CoverageAbove20, 20), (CoverageAbove30, 30), (CoverageAbove60, 60), (CoverageAbove100, 100), (CoverageAbove250, 250)] {
        let above = coverage.iter().filter(|v| v.value >= level).fold(0i64, |s, v| s.wrapping_add(v.count));
        let prop = above as f64 / total as f64;
        let status = thresholds.get_summary(sample_type, f, QcStatusType::Warn)?.sample_qc_status(prop);
        out.push((f, f.feature(prop, status)));
    }
    Ok(out)
}

/// SummaryTableBamMetricsPrep.createFeatures: the features of the files there are, in enum order; the low
/// coverage exons from the exons' coverage QC, when the exon coverage file was read
pub fn bam_metrics_features(
    data: &BamMetricsData, coverage_qc: Option<&[ExonCoverageQc]>, sample_type: SampleType, thresholds: &crate::status::ThresholdRegistry,
) -> Result<Vec<Feature>> {
    use SummaryTableFeature::*;
    let mut map: Vec<(SummaryTableFeature, Feature)> = Vec::new();
    let warn = |f: SummaryTableFeature, v: f64| -> Result<(SummaryTableFeature, Feature)> {
        Ok((f, f.feature(v, thresholds.get_summary(sample_type, f, QcStatusType::Warn)?.sample_qc_status(v))))
    };
    if let Some(s) = &data.summary {
        map.push(warn(MeanCoverage, s.mean_coverage)?);
        map.push(warn(LowMapQual, s.low_map_qual_percent)?);
        map.push(warn(LowBaseQual, s.low_base_qual_percent)?);
        map.push(warn(DuplicateReads, s.duplicate_reads as f64 / s.total_reads as f64)?);
        map.push(warn(DualStrandReads, s.dual_strand_reads as f64 / s.total_reads as f64)?);
    }
    if let Some(c) = &data.coverage {
        map.extend(coverage_above_features(c, sample_type, thresholds)?);
    }
    if let Some(fs) = &data.flag_stats {
        let status = thresholds.get_summary(sample_type, MappedProportion, QcStatusType::Fail)?.sample_qc_status(fs.mapped_proportion);
        map.push((MappedProportion, MappedProportion.feature(fs.mapped_proportion, status)));
    }
    if let Some(exons) = coverage_qc {
        map.push(warn(LowCoverageExons, coverage::low_coverage_exon_proportion(exons))?);
    }
    // an EnumMap: values in enum order
    map.sort_by_key(|(f, _)| *f);
    Ok(map.into_iter().map(|(_, f)| f).collect())
}

fn summary_table_bam_metrics(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    let data = load_bam_metrics_data(config, sample_id, sample_type)?;
    if !data.missing.is_empty() {
        let names: Vec<&str> = data.missing.iter().map(|p| p.rsplit('/').next().unwrap_or(p)).collect();
        let dir = config.bam_metrics_dir(sample_id, sample_type)?;
        missing_input_files_error(config.allow_missing_input, PrepCategory::SummaryTableBamMetrics, sample_type, &format!("{} in dir: {dir}", names.join(", ")))?;
    }
    // SummaryTableBamMetricsPrep.calcCoverageQcPerExon
    let coverage_qc = match &data.exon_depths {
        Some(exons) => {
            let reportable = reportable_genes(config, sample_type);
            Some(coverage::coverage_qc_per_exon(sample_id, sample_type, exons, &config.exon_exclusions, &reportable, &config.thresholds, config.lod_min_alt_reads)?)
        }
        None => None,
    };
    Ok(bam_metrics_features(&data, coverage_qc.as_deref(), sample_type, &config.thresholds)?)
}

// ---- summary table: PURPLE ----

/// PurityContext: what the PURPLE features take from the QC and purity files
pub struct PurityContext {
    pub qc: PurpleQc,
    pub purity: PurplePurity,
}

/// PurityContextFile.readWithQC: the QC file, then the purity file
fn load_purity_context(config: &QseeConfig, sample_id: &str) -> std::result::Result<PurityContext, PrepError> {
    let dir = config.purple_dir(sample_id)?;
    let qc_path = purple_qc_filename(&dir, sample_id);
    read_input(&qc_path)?;
    let qc = read_purple_qc(&qc_path)?;
    let purity_path = purple_purity_filename(&dir, sample_id);
    let lines = read_input(&purity_path)?;
    let purity = purple_purity_from_lines(&lines)?;
    Ok(PurityContext { qc, purity })
}

/// SummaryTablePurplePrep.createQcStatusMappings: each PURPLE status's qsee status, PASS (with the
/// threshold) for those PURPLE did not set
pub fn purple_qc_status_mappings(purple_statuses: &[String], thresholds: &crate::status::ThresholdRegistry) -> Result<Vec<(&'static str, QcStatus)>> {
    let mut out = Vec::new();
    for s in PURPLE_QC_STATUSES {
        let qsee = purple_to_qsee_status(thresholds, s)?;
        let status = if purple_statuses.iter().any(|p| p == s) { qsee } else { QcStatus::pass(&qsee) };
        out.push((s, status));
    }
    Ok(out)
}

fn mapping(m: &[(&'static str, QcStatus)], s: &str) -> QcStatus {
    m.iter().find(|(k, _)| *k == s).map(|(_, v)| *v).expect("every PURPLE status is mapped")
}

/// SummaryTablePurplePrep.getTincQcStatus: FAIL_TINC if set, else WARN_TINC if set, else PASS with the
/// fail threshold
pub fn tinc_qc_status(m: &[(&'static str, QcStatus)]) -> QcStatus {
    worst_status(m, "FAIL_TINC", "WARN_TINC")
}

fn worst_status(m: &[(&'static str, QcStatus)], fail: &str, warn: &str) -> QcStatus {
    let maybe_fail = mapping(m, fail);
    if maybe_fail.status_type != QcStatusType::Pass {
        return maybe_fail;
    }
    let maybe_warn = mapping(m, warn);
    if maybe_warn.status_type != QcStatusType::Pass {
        return maybe_warn;
    }
    QcStatus::pass(&maybe_fail)
}

/// SummaryTablePurplePrep.createFeatures
pub fn purple_features(ctx: &PurityContext, thresholds: &crate::status::ThresholdRegistry) -> Result<Vec<Feature>> {
    use SummaryTableFeature::*;
    let m = purple_qc_status_mappings(&ctx.qc.status, thresholds)?;
    // QS5: Java takes only FAIL_CONTAMINATION, so PURPLE's WARN_CONTAMINATION shows as PASS
    let contamination = if compat::java("QS5") { mapping(&m, "FAIL_CONTAMINATION") } else { worst_status(&m, "FAIL_CONTAMINATION", "WARN_CONTAMINATION") };
    Ok(vec![
        Purity.feature(ctx.purity.purity, mapping(&m, "WARN_LOW_PURITY")),
        Ploidy.feature(ctx.purity.ploidy, QcStatus::empty()),
        LohPercent.feature(ctx.qc.loh_percent, QcStatus::empty()),
        DeletedGenes.feature(ctx.qc.deleted_genes as f64, mapping(&m, "WARN_DELETED_GENES")),
        UnsupportedCnSegments.feature(ctx.qc.unsupported_copy_number_segments as f64, mapping(&m, "WARN_HIGH_COPY_NUMBER_NOISE")),
        Tinc.feature(ctx.qc.tinc_level, tinc_qc_status(&m)),
        Contamination.feature(ctx.qc.contamination, contamination),
        TmbSmallVariants.feature(ctx.purity.tmb_per_mb, QcStatus::empty()),
        TmbMsIndels.feature(ctx.purity.ms_indels_per_mb, QcStatus::empty()),
        // RefGenomeConstants.MB_PER_GENOME
        TmbStructuralVariants.feature(ctx.purity.sv_tumor_mutational_burden as f64 / 2859.0, QcStatus::empty()),
    ])
}

fn summary_table_purple(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    if sample_type != SampleType::Tumor {
        return Ok(Vec::new());
    }
    let ctx = load_purity_context(config, sample_id)?;
    Ok(purple_features(&ctx, &config.thresholds)?)
}

// ---- distributions ----

fn distribution_features(values: &[ValueFrequency], title: &str, feature_type: FeatureType, tool: SourceTool) -> Result<Vec<Feature>> {
    let b = BinnedFrequencies::from_value_frequencies(values)?;
    Ok(b.proportional_density_features(title).into_iter().map(|(n, v)| Feature::new(FeatureKey::new(n, feature_type, tool), v)).collect())
}

fn coverage_distribution(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    let path = BamMetricCoverage::generate_filename(&config.bam_metrics_dir(sample_id, sample_type)?, sample_id);
    read_input(&path)?;
    let coverage = BamMetricCoverage::read(&path)?;
    Ok(distribution_features(&coverage, "ReadDepth", FeatureType::CoverageDistribution, SourceTool::BamMetrics)?)
}

fn frag_length_distribution(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    let path = BamMetricFragmentLength::generate_filename(&config.bam_metrics_dir(sample_id, sample_type)?, sample_id);
    read_input(&path)?;
    let lengths = BamMetricFragmentLength::read(&path)?;
    Ok(distribution_features(&lengths, "FragLength", FeatureType::FragLengthDistribution, SourceTool::BamMetrics)?)
}

fn duplicate_freq(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    if sample_type != SampleType::Tumor {
        return Ok(Vec::new());
    }
    let path = redux_files::duplicate_freq_filename(&config.redux_dir_unjoined(sample_id, sample_type)?, sample_id);
    let lines = read_input(&path)?;
    let freqs = redux_files::duplicate_frequencies_from_lines(&lines)?;
    let b = BinnedFrequencies::new(freqs.iter().map(|f| f.read_count as i64).collect(), freqs.iter().map(|f| f.count as f64).collect())?;
    Ok(b.proportional_density_features("ReadCount").into_iter().map(|(n, v)| Feature::new(FeatureKey::new(n, FeatureType::DuplicateFreq, SourceTool::Redux), v)).collect())
}

// ---- missed variant likelihood ----

/// ReportableGenes.forSampleType (MissedGeneVariantPrep.getReportableGenes before upstream bdceb2302a): the
/// panel's genes reported somatically (tumor) or in the germline (normal), in panel order
pub fn reportable_genes(config: &QseeConfig, sample_type: SampleType) -> Vec<String> {
    reportable_genes_of(&config.driver_genes, sample_type)
}

pub fn reportable_genes_of(driver_genes: &[DriverGene], sample_type: SampleType) -> Vec<String> {
    driver_genes
        .iter()
        .filter(|g| if sample_type == SampleType::Tumor { g.report_somatic() } else { g.report_germline() })
        .map(|g| g.gene.clone())
        .collect()
}

/// a GeneDepth's depth counts as Java's LinkedHashMap holds them: a repeated range keeps its first place
/// and its last count
pub fn depth_counts_map(counts: &[(DepthRange, i32)]) -> Vec<(DepthRange, i32)> {
    let mut out: Vec<(DepthRange, i32)> = Vec::with_capacity(counts.len());
    for (r, c) in counts {
        match out.iter_mut().find(|(x, _)| x == r) {
            Some(e) => e.1 = *c,
            None => out.push((*r, *c)),
        }
    }
    out
}

/// HumanChromosome.fromString(..).isDiploid(gender): an error for a name that is not 1-22, X or Y
pub fn is_diploid(chromosome: &str, gender: &str) -> Result<bool> {
    let name = if chromosome.to_lowercase().starts_with("chr") { &chromosome[3..] } else { chromosome };
    let autosome = matches!(name.parse::<u32>(), Ok(1..=22)) && !name.starts_with('+') && !name.starts_with('0');
    if autosome {
        return Ok(true);
    }
    match name {
        "X" => Ok(gender != "MALE"),
        "Y" => Ok(false),
        _ => Err(anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._{name}")),
    }
}

/// MissedGeneVariantPrep.getMissedVariantLikelihoods (upstream bdceb2302a): per reportable gene, the depth
/// counts of its exons that are not excluded, summed (genes in first-seen order, ranges in the first exon's
/// order; the chromosome of the gene's last exon), and their missed variant likelihood. A gene with every exon
/// excluded has none.
pub fn missed_variant_likelihoods(
    exons: &[ExonDepth], reportable: &[String], exclusions: &ExonExclusions, purity: f64, gender: &str, min_alt_reads: i32,
) -> Result<Vec<Feature>> {
    // LinkedHashMap<gene, LinkedHashMap<range, count>> and HashMap<gene, chromosome>
    let mut genes: Vec<(String, Vec<(DepthRange, i32)>, String)> = Vec::new();
    for e in exons {
        if !reportable.contains(&e.depth.gene) || exclusions.is_excluded(&e.depth.gene, e.exon_rank) {
            continue;
        }
        let i = match genes.iter().position(|(g, _, _)| *g == e.depth.gene) {
            Some(i) => i,
            None => {
                genes.push((e.depth.gene.clone(), Vec::new(), String::new()));
                genes.len() - 1
            }
        };
        let counts = &mut genes[i].1;
        // exonDepth.DepthCounts.forEach(merge(range, count, Integer::sum))
        for (r, c) in depth_counts_map(&e.depth.depth_counts) {
            match counts.iter_mut().find(|(x, _)| *x == r) {
                Some(x) => x.1 = x.1.wrapping_add(c),
                None => counts.push((r, c)),
            }
        }
        genes[i].2 = e.depth.chromosome.clone();
    }
    let mut out = Vec::new();
    for (gene, counts, chromosome) in &genes {
        let key = FeatureKey::new(form_single_field("Gene", gene), FeatureType::MissedVariantLikelihood, SourceTool::BamMetrics);
        let alleles = if is_diploid(chromosome, gender)? { 2 } else { 1 };
        let vaf = purity * (1.0 / alleles as f64);
        out.push(Feature::new(key, coverage::missed_variant_likelihood(counts, vaf, min_alt_reads)));
    }
    Ok(out)
}

fn load_purple_purity(config: &QseeConfig, sample_id: &str) -> std::result::Result<PurplePurity, PrepError> {
    let path = purple_purity_filename(&config.purple_dir(sample_id)?, sample_id);
    let lines = read_input(&path)?;
    Ok(purple_purity_from_lines(&lines)?)
}

fn missed_gene_variant(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    let path = GeneDepthFile::generate_exon_coverage_filename(&config.bam_metrics_dir(sample_id, sample_type)?, sample_id);
    read_input(&path)?;
    let exons = GeneDepthFile::read_exon_coverage(&path)?;
    let reportable = reportable_genes(config, sample_type);
    let (purity, gender) = if sample_type == SampleType::Tumor {
        let p = load_purple_purity(config, sample_id)?;
        (p.purity, p.sex)
    } else {
        // the gender of the paired tumor sample's PURPLE output
        let i = config.samples.reference_ids.iter().position(|r| r == sample_id).ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: Index -1"))?;
        let tumor = config.samples.tumor_ids.get(i).ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: Index {i}"))?;
        (coverage::REF_SAMPLE_PURITY, load_purple_purity(config, tumor)?.sex)
    };
    Ok(missed_variant_likelihoods(&exons, &reportable, &config.exon_exclusions, purity, &gender, config.lod_min_alt_reads)?)
}

// ---- GC bias ----

/// GcBiasPrep.normaliseMedianReadDepths: by bucket, the median over the sample median (-1 read as 0)
pub fn gc_bias_features(g: &cobalt_files::GcMedianReadDepth) -> Vec<Feature> {
    let mut buckets = g.median_per_bucket.clone();
    buckets.sort_by_key(|(b, _)| *b);
    buckets
        .into_iter()
        .map(|(b, m)| {
            let m = if m == cobalt_files::NO_READ_DEPTH_VALUE { 0.0 } else { m };
            Feature::new(FeatureKey::new(form_single_field("GCBucket", &b.to_string()), FeatureType::GcBias, SourceTool::Cobalt), m / g.median)
        })
        .collect()
}

fn gc_bias(config: &QseeConfig, sample_id: &str) -> PrepResult {
    let path = cobalt_files::cobalt_gc_median_filename(&config.cobalt_dir(sample_id)?, sample_id);
    let lines = match read_input(&path) {
        Ok(l) => l,
        Err(e) => {
            hmf_common::error!("Failed to load Cobalt GC median read depth file for sample: {sample_id}");
            return Err(e);
        }
    };
    let g = cobalt_files::gc_median_from_lines(&lines)?;
    Ok(gc_bias_features(&g))
}

// ---- discordant fragments ----

/// DiscordantFragFreqPrep.calcDiscordantProportions and formFeatures: each group's count over the total
/// reads, in group order
pub fn discordant_features(stats: &esvee_files::EsveeDiscordantStats) -> Vec<Feature> {
    let mut groups: Vec<(DiscordantFragGroup, i64)> = Vec::new();
    for (i, c) in stats.type_counts.iter().enumerate() {
        let g = DiscordantFragGroup::from_type_index(i);
        match groups.iter_mut().find(|(x, _)| *x == g) {
            Some(e) => e.1 = e.1.wrapping_add(*c),
            None => groups.push((g, *c)),
        }
    }
    groups.sort_by_key(|(g, _)| *g);
    groups
        .into_iter()
        .map(|(g, c)| {
            let key = FeatureKey::new(form_single_field("DiscordantFragType", g.name()), FeatureType::DiscordantFragFreq, SourceTool::Esvee);
            let metadata = FeatureMetadata { display_name: g.display_name().to_string(), ..Default::default() };
            Feature::with_metadata(key, c as f64 / stats.total_reads as f64, metadata)
        })
        .collect()
}

fn discordant_frag_freq(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    if sample_type != SampleType::Tumor {
        return Ok(Vec::new());
    }
    let dir = config.esvee_dir(sample_id)?;
    let new_path = esvee_files::esvee_disc_stats_filename(&dir, sample_id);
    let old_path = format!("{dir}/{sample_id}.esvee.prep.disc_stats.tsv");
    let path = if std::path::Path::new(&new_path).is_file() {
        new_path
    } else if std::path::Path::new(&old_path).is_file() {
        old_path
    } else {
        let name = |p: &str| p.rsplit('/').next().unwrap_or(p).to_string();
        return Err(PrepError::Missing(format!("{} or {}", name(&new_path), name(&old_path))));
    };
    let lines = read_input(&path)?;
    let stats = esvee_files::esvee_disc_stats_from_lines(&lines)?;
    Ok(discordant_features(&stats))
}

// ---- base quality recalibration ----

/// BaseQualBinRange
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaseQualBinRange {
    pub bin: BaseQualBin,
    pub lower: i8,
    pub upper: i8,
}

impl BaseQualBinRange {
    /// toString: "HIGH (30+)" or "LOW (0-29)"
    pub fn display(&self) -> String {
        if self.bin == BaseQualBin::High { format!("{} ({}+)", self.bin.name(), self.lower) } else { format!("{} ({}-{})", self.bin.name(), self.lower, self.upper) }
    }
}

/// BaseQualBinner: the quality bins of a sequencing type (BaseQualAdjustment's low, medium and high
/// qualities), each from its lowest quality to the one before the next bin's
#[derive(Clone, Debug)]
pub struct BaseQualBinner {
    pub ranges: Vec<BaseQualBinRange>,
}

impl BaseQualBinner {
    pub fn new(sequencing_type: SequencingType) -> BaseQualBinner {
        let within = |bin: BaseQualBin, q: i8| -> bool {
            match bin {
                // BaseQualAdjustment.isLowBaseQual: at or below 26
                BaseQualBin::Low => q <= 26,
                BaseQualBin::Medium => sequencing_type == SequencingType::Sbx && (27..=29).contains(&q),
                BaseQualBin::High => match sequencing_type {
                    SequencingType::Illumina => q >= 30,
                    SequencingType::Sbx => q > 29,
                    SequencingType::Ultima => q >= 35,
                },
            }
        };
        let mut lowers: Vec<(BaseQualBin, i8)> = Vec::new();
        for q in 0..i8::MAX {
            for bin in [BaseQualBin::Low, BaseQualBin::Medium, BaseQualBin::High] {
                if within(bin, q) && !lowers.iter().any(|(b, _)| *b == bin) {
                    lowers.push((bin, q));
                }
            }
        }
        lowers.sort_by_key(|(b, _)| *b);
        let mut ranges = Vec::new();
        for (i, (bin, lower)) in lowers.iter().enumerate() {
            let upper = if *bin == BaseQualBin::High {
                i8::MAX
            } else {
                match lowers.get(i + 1) {
                    Some((_, next)) => next.wrapping_sub(1),
                    // no upper bound set (cannot happen: HIGH exists for every type); Java would fail
                    None => i8::MAX,
                }
            };
            ranges.push(BaseQualBinRange { bin: *bin, lower: *lower, upper });
        }
        BaseQualBinner { ranges }
    }

    pub fn range(&self, bin: BaseQualBin) -> Option<&BaseQualBinRange> {
        self.ranges.iter().find(|r| r.bin == bin)
    }

    /// binFrom: an error (RuntimeException) for a quality in no bin
    pub fn bin_from(&self, q: i8) -> Result<&BaseQualBinRange> {
        self.ranges.iter().find(|r| q >= r.lower && q <= r.upper).ok_or_else(|| anyhow!("java.lang.RuntimeException: No bin found for base quality: {q}"))
    }
}

/// Java's (char) of a byte: sign-extended, so a byte above 127 is U+FF80 to U+FFFF
fn java_char(b: u8) -> char {
    char::from_u32((b as i8 as i32 as u32) & 0xFFFF).unwrap_or('\u{FFFD}')
}

/// Nucleotides.swapDnaBase(char): A, T, C and G complemented, any other character kept
fn swap_dna_char(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => c,
    }
}

/// Orientation of a standard context's read (upstream 883ab11bba)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadOrientation {
    Forward,
    Reverse,
}

impl ReadOrientation {
    pub fn name(self) -> &'static str {
        match self {
            ReadOrientation::Forward => "FORWARD",
            ReadOrientation::Reverse => "REVERSE",
        }
    }
    /// Orientation.asChar
    #[allow(dead_code)] // Java's, for mutationContextString
    pub fn as_char(self) -> char {
        match self {
            ReadOrientation::Forward => '+',
            ReadOrientation::Reverse => '-',
        }
    }
}

/// StandardTrinucleotideContext (upstream 883ab11bba): the mutation with a C or T reference. REDUX's BQR keys are in
/// the read's own orientation, so a G or A reference (complemented here, its context reverse complemented) is a read
/// of the strand opposite the standard one: REVERSE.
#[derive(Clone, Debug, PartialEq)]
pub struct StandardTrinucleotideContext {
    pub ref_base: char,
    pub alt: char,
    pub context: String,
    pub orientation: ReadOrientation,
}

impl StandardTrinucleotideContext {
    /// fromNucleotides
    pub fn from_nucleotides(ref_base: char, alt: char, context: &str) -> StandardTrinucleotideContext {
        let standard = ref_base == 'C' || ref_base == 'T';
        let (r, a) = if standard { (ref_base, alt) } else { (swap_dna_char(ref_base), swap_dna_char(alt)) };
        let orientation = if ref_base == r { ReadOrientation::Forward } else { ReadOrientation::Reverse };
        let context = match orientation {
            ReadOrientation::Forward => context.to_string(),
            // Nucleotides.reverseComplementBases(String)
            ReadOrientation::Reverse => context.chars().rev().map(swap_dna_char).collect(),
        };
        StandardTrinucleotideContext { ref_base: r, alt: a, context, orientation }
    }
    /// mutationContextString: e.g. "G[T>G]A +"
    #[allow(dead_code)] // Java's; its tests use it
    pub fn mutation_context_string(&self) -> Result<String> {
        let c: Vec<char> = self.context.chars().collect();
        if c.len() < 3 {
            anyhow::bail!("java.lang.StringIndexOutOfBoundsException: context '{}'", self.context);
        }
        Ok(format!("{}[{}>{}]{} {}", c[0], self.ref_base, self.alt, c[2], self.orientation.as_char()))
    }
    /// mutationString: e.g. "C>A"
    pub fn mutation_string(&self) -> String {
        format!("{}>{}", self.ref_base, self.alt)
    }
    /// mutationStringSplitCpG: C>T as "C>T CpG" (a G after the C) or "C>T other"
    pub fn mutation_string_split_cpg(&self) -> Result<String> {
        let mut m = self.mutation_string();
        if self.ref_base == 'C' && self.alt == 'T' {
            let third = self.context.chars().nth(2).ok_or_else(|| anyhow!("java.lang.StringIndexOutOfBoundsException: context '{}'", self.context))?;
            m.push_str(if third == 'G' { " CpG" } else { " other" });
        }
        Ok(m)
    }
}

/// StandardBqrRecord (upstream 883ab11bba)
#[derive(Clone, Debug)]
pub struct StandardBqrRecord {
    pub context: StandardTrinucleotideContext,
    pub original_quality: i8,
    pub read_type: redux_files::ConsensusType,
    pub count: i32,
    pub recalibrated_quality: f64,
}

impl StandardBqrRecord {
    /// StandardBqrRecord.from
    pub fn from(r: &BqrRecord) -> StandardBqrRecord {
        // new String(TrinucleotideContext) (UTF-8) and (char) of the bases
        let context = String::from_utf8_lossy(&r.trinucleotide_context).into_owned();
        StandardBqrRecord {
            context: StandardTrinucleotideContext::from_nucleotides(java_char(r.ref_base), java_char(r.alt), &context),
            original_quality: r.quality,
            read_type: r.read_type,
            count: r.count,
            recalibrated_quality: r.recalibrated_quality,
        }
    }
}

/// loadSnvBqrRecords after reading: SNVs only, standardised, sorted (stable) by read type, mutation string and
/// standard context
pub fn prepare_bqr_records(records: Vec<BqrRecord>) -> Vec<StandardBqrRecord> {
    let mut out: Vec<StandardBqrRecord> = records.iter().filter(|r| r.ref_base != r.alt).map(StandardBqrRecord::from).collect();
    out.sort_by(|a, b| {
        a.read_type
            .cmp(&b.read_type)
            .then_with(|| hmf_common::java_text::compare(&a.context.mutation_string(), &b.context.mutation_string()))
            .then_with(|| hmf_common::java_text::compare(&a.context.context, &b.context.context))
    });
    out
}

const PSEUDOCOUNT: i32 = 5;
/// groups with fewer bases than this have no value (NaN), upstream 883ab11bba
const MIN_COUNT: i32 = 5;

/// calcWeightedMeanChangeInQualPerGroup: NaN for a group of fewer than 5 counts; else commons-math3
/// Mean.evaluate(values, weights) of the change in quality, weighted by the count plus a pseudocount
fn weighted_mean_change(records: &[&StandardBqrRecord]) -> Result<f64> {
    let total = records.iter().fold(0i32, |s, r| s.wrapping_add(r.count));
    if total < MIN_COUNT {
        return Ok(f64::NAN);
    }
    let values: Vec<f64> = records.iter().map(|r| r.recalibrated_quality - r.original_quality as f64).collect();
    let weights: Vec<f64> = records.iter().map(|r| r.count.wrapping_add(PSEUDOCOUNT) as f64).collect();
    weighted_mean(&values, &weights)
}

/// commons-math3 Mean.evaluate(values, weights): the definitional mean with a second-pass correction.
/// A NaN, infinite or negative weight, or no positive weight, is an error (MathIllegalArgumentException);
/// no values give NaN.
pub fn weighted_mean(values: &[f64], weights: &[f64]) -> Result<f64> {
    if weights.iter().any(|w| w.is_nan() || w.is_infinite() || *w < 0.0) || !weights.iter().any(|w| *w > 0.0) {
        anyhow::bail!("org.apache.commons.math3.exception.MathIllegalArgumentException: invalid weights");
    }
    if values.is_empty() {
        return Ok(f64::NAN);
    }
    let mut sumw = 0.0;
    for w in weights {
        sumw += w;
    }
    let mut sum = 0.0;
    for (v, w) in values.iter().zip(weights) {
        sum += v * w;
    }
    let xbarw = sum / sumw;
    let mut correction = 0.0;
    for (v, w) in values.iter().zip(weights) {
        correction += w * (v - xbarw);
    }
    Ok(xbarw + correction / sumw)
}

fn grouped_features(records: &[&StandardBqrRecord], feature_type: FeatureType, name: impl Fn(&StandardBqrRecord) -> Result<String>) -> Result<Vec<Feature>> {
    let mut groups: Vec<(String, Vec<&StandardBqrRecord>)> = Vec::new();
    for r in records {
        let n = name(r)?;
        match groups.iter_mut().find(|(k, _)| *k == n) {
            Some(g) => g.1.push(r),
            None => groups.push((n, vec![r])),
        }
    }
    groups.into_iter().map(|(n, rs)| Ok(Feature::new(FeatureKey::new(n, feature_type, SourceTool::Redux), weighted_mean_change(&rs)?))).collect()
}

/// the value of a grouping field when records are not grouped by it
const ANY: &str = "ANY";

/// formGroupingKey
fn grouping_key(read_type: &str, mutation: &str, context: &str, orientation: &str, qual_bin: &str) -> String {
    form_multi_field(&[
        ("ReadType", read_type),
        ("StandardMutation", mutation),
        ("StandardTrinucContext", context),
        ("ReadOrientation", orientation),
        ("OriginalQualBin", qual_bin),
    ])
}

/// aggregateChangeInQualPerOriginalQual: by read type, mutation (C>T split into CpG and other), read orientation
/// and quality bin
pub fn bqr_per_original_qual(records: &[StandardBqrRecord], binner: &BaseQualBinner) -> Result<Vec<Feature>> {
    let refs: Vec<&StandardBqrRecord> = records.iter().collect();
    grouped_features(&refs, FeatureType::BqrPerOrigQual, |r| {
        Ok(grouping_key(r.read_type.name(), &r.context.mutation_string_split_cpg()?, ANY, r.context.orientation.name(), &binner.bin_from(r.original_quality)?.display()))
    })
}

/// aggregateChangeInQualPerTrinucContext: high qualities only, by read type, mutation, context, read orientation
/// and quality bin
pub fn bqr_per_trinuc_context(records: &[StandardBqrRecord], binner: &BaseQualBinner) -> Result<Vec<Feature>> {
    let high = binner.range(BaseQualBin::High).ok_or_else(|| anyhow!("java.lang.NullPointerException: no HIGH bin"))?.lower;
    let refs: Vec<&StandardBqrRecord> = records.iter().filter(|r| r.original_quality >= high).collect();
    grouped_features(&refs, FeatureType::BqrPerSnv96Context, |r| {
        Ok(grouping_key(r.read_type.name(), &r.context.mutation_string(), &r.context.context, r.context.orientation.name(), &binner.bin_from(r.original_quality)?.display()))
    })
}

fn base_qual_recalibration(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    let path = redux_files::bqr_filename(&config.redux_dir(sample_id, sample_type)?, sample_id);
    // QS1: BqrFile.read returns null for a file it cannot read, and Java's NullPointerException then
    // ends the run, -allow_missing_input or not
    let lines = match read_input(&path) {
        Ok(l) => l,
        Err(PrepError::Missing(m)) => {
            if compat::java("QS1") {
                hmf_common::error!("failed to read BQR file({path}) record index(0): java.nio.file.NoSuchFileException: {m}");
                return Err(PrepError::Fatal(anyhow!("java.lang.NullPointerException: no BQR records")));
            }
            return Err(PrepError::Missing(m));
        }
        Err(e) => return Err(e),
    };
    let records = redux_files::bqr_records_from_lines(&lines).map_err(|e| {
        hmf_common::error!("failed to read BQR file({path}): {e:#}");
        anyhow!("java.lang.NullPointerException: no BQR records")
    })?;
    let records = prepare_bqr_records(records);
    let binner = BaseQualBinner::new(config.sequencing_type);
    let mut features = bqr_per_original_qual(&records, &binner)?;
    features.extend(bqr_per_trinuc_context(&records, &binner)?);
    Ok(features)
}

// ---- microsatellite indel errors ----

const MAX_REF_UNIT_COUNT: i32 = 15;
pub const INSERTIONS_KEY: i32 = 1;
pub const DELETIONS_KEY: i32 = -1;

/// aggregateByJitter: the insertions (jitter above 0) and deletions (below 0) of each row, summed as ints
pub fn aggregate_by_jitter(table: &[JitterTableRow]) -> Vec<JitterTableRow> {
    table
        .iter()
        .map(|row| {
            let ins = row.jitter_counts.iter().filter(|(j, _)| *j > 0).fold(0i32, |s, (_, c)| s.wrapping_add(*c));
            let del = row.jitter_counts.iter().filter(|(j, _)| *j < 0).fold(0i32, |s, (_, c)| s.wrapping_add(*c));
            let mut r = JitterTableRow::new(row.ref_num_units, &row.repeat_unit, row.consensus_type);
            r.total_read_count = row.total_read_count;
            r.set_jitter_read_count(INSERTIONS_KEY, ins);
            r.set_jitter_read_count(DELETIONS_KEY, del);
            r
        })
        .collect()
}

/// aggregateByRepeatType: rows of the same consensus type, repeat type and number of units summed (in
/// the order the groups first appear), the repeat unit replaced by the repeat type's name
pub fn aggregate_by_repeat_type(table: &[JitterTableRow]) -> Result<Vec<JitterTableRow>> {
    let mut groups: Vec<(String, Vec<&JitterTableRow>)> = Vec::new();
    for row in table {
        let repeat = RepeatType::from_repeat_unit(&row.repeat_unit)?.display_name();
        let name = form_multi_field(&[("ConsensusType", row.consensus_type.name()), ("RepeatUnitType", repeat), ("RefNumUnits", &row.ref_num_units.to_string())]);
        match groups.iter_mut().find(|(k, _)| *k == name) {
            Some(g) => g.1.push(row),
            None => groups.push((name, vec![row])),
        }
    }
    let mut out = Vec::new();
    for (_, rows) in groups {
        let first = rows[0];
        let mut r = JitterTableRow::new(first.ref_num_units, RepeatType::from_repeat_unit(&first.repeat_unit)?.display_name(), first.consensus_type);
        let (mut total, mut ins, mut del) = (0i32, 0i32, 0i32);
        for x in &rows {
            total = total.wrapping_add(x.total_read_count);
            ins = ins.wrapping_add(x.jitter_read_count(INSERTIONS_KEY));
            del = del.wrapping_add(x.jitter_read_count(DELETIONS_KEY));
        }
        r.total_read_count = total;
        r.set_jitter_read_count(INSERTIONS_KEY, ins);
        r.set_jitter_read_count(DELETIONS_KEY, del);
        out.push(r);
    }
    Ok(out)
}

/// calcPhredScore: -10 log10(count / total), NaN for a zero count or total. Java's Math.log10 on
/// aarch64 is fdlibm's, which StrictMath.log10 is.
fn phred(count: i32, total: i32) -> f64 {
    if count == 0 || total == 0 {
        return f64::NAN;
    }
    -10.0 * java_strictmath::log10(count as f64 / total as f64)
}

/// calcPhredScores: the error rates of every row (where defined), then the deletion-insertion biases
pub fn ms_indel_phred_scores(table: &[JitterTableRow]) -> Vec<Feature> {
    let mut rates = Vec::new();
    let mut biases = Vec::new();
    for row in table {
        let name = form_multi_field(&[("ConsensusType", row.consensus_type.name()), ("RepeatUnitType", &row.repeat_unit), ("RefNumUnits", &row.ref_num_units.to_string())]);
        let del = row.jitter_read_count(DELETIONS_KEY);
        let ins = row.jitter_read_count(INSERTIONS_KEY);
        let indel = del.wrapping_add(ins);
        let rate = phred(indel, row.total_read_count);
        if !rate.is_nan() {
            rates.push(Feature::new(FeatureKey::new(name.clone(), FeatureType::MsIndelErrorRates, SourceTool::Redux), rate));
        }
        let bias = phred(del, row.total_read_count) - phred(ins, row.total_read_count);
        if !bias.is_nan() {
            biases.push(Feature::new(FeatureKey::new(name, FeatureType::MsIndelErrorBias, SourceTool::Redux), bias));
        }
    }
    rates.extend(biases);
    rates
}

fn ms_indel_error(config: &QseeConfig, sample_id: &str, sample_type: SampleType) -> PrepResult {
    let dir = config.redux_dir(sample_id, sample_type)?;
    let new_path = redux_files::jitter_counts_filename(&dir, sample_id);
    let old_path = format!("{dir}/{sample_id}.ms_table.tsv.gz");
    let path = if std::path::Path::new(&new_path).is_file() {
        new_path
    } else if std::path::Path::new(&old_path).is_file() {
        old_path
    } else {
        let name = |p: &str| p.rsplit('/').next().unwrap_or(p).to_string();
        return Err(PrepError::Missing(format!("sample({sample_id}) - could not determine the REDUX({}) or SAGE({}) jitter file", name(&new_path), name(&old_path))));
    };
    let tables = redux_files::read_jitter_counts_tables(&path)?;
    let rows: Vec<JitterTableRow> = tables.into_iter().flat_map(|t| t.rows).filter(|r| r.ref_num_units <= MAX_REF_UNIT_COUNT).collect();
    let rows = aggregate_by_repeat_type(&aggregate_by_jitter(&rows))?;
    Ok(ms_indel_phred_scores(&rows))
}

