//! QseeConfig, QseePlotConfig and TrainConfig: the options as Java registers them, and the configuration
//! made from them.

use anyhow::{anyhow, bail, Result};

use hmf_common::config_builder::ConfigBuilder;
use hmf_common::driver_panel::DriverGene;
use hmf_common::metrics::check_add_dir_separator;

use crate::common::{SampleIds, SampleType};
use crate::compat;
use crate::status::{read_threshold_overrides, ThresholdRegistry};

/// SequencingType, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequencingType {
    Illumina,
    Ultima,
    Sbx,
}

impl SequencingType {
    pub fn value_of(s: &str) -> Result<SequencingType> {
        match s {
            "ILLUMINA" => Ok(SequencingType::Illumina),
            "ULTIMA" => Ok(SequencingType::Ultima),
            "SBX" => Ok(SequencingType::Sbx),
            _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.sequencing.SequencingType.{s}"),
        }
    }
}

/// PrepCategory, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PrepCategory {
    BaseQualRecalibration,
    CoverageDistribution,
    DiscordantFragFreq,
    DuplicateFreq,
    FragLengthDistribution,
    GcBias,
    MissedGeneVariant,
    MsIndelError,
    SummaryTableBamMetrics,
    SummaryTablePurple,
}

impl PrepCategory {
    pub const ALL: [PrepCategory; 10] = [
        PrepCategory::BaseQualRecalibration,
        PrepCategory::CoverageDistribution,
        PrepCategory::DiscordantFragFreq,
        PrepCategory::DuplicateFreq,
        PrepCategory::FragLengthDistribution,
        PrepCategory::GcBias,
        PrepCategory::MissedGeneVariant,
        PrepCategory::MsIndelError,
        PrepCategory::SummaryTableBamMetrics,
        PrepCategory::SummaryTablePurple,
    ];

    /// FeaturePrep.createCategoryPreps: the order the categories are prepared in
    pub const PREP_ORDER: [PrepCategory; 10] = [
        PrepCategory::SummaryTableBamMetrics,
        PrepCategory::SummaryTablePurple,
        PrepCategory::CoverageDistribution,
        PrepCategory::FragLengthDistribution,
        PrepCategory::MissedGeneVariant,
        PrepCategory::DuplicateFreq,
        PrepCategory::GcBias,
        PrepCategory::DiscordantFragFreq,
        PrepCategory::BaseQualRecalibration,
        PrepCategory::MsIndelError,
    ];

    pub fn name(self) -> &'static str {
        match self {
            PrepCategory::BaseQualRecalibration => "BASE_QUAL_RECALIBRATION",
            PrepCategory::CoverageDistribution => "COVERAGE_DISTRIBUTION",
            PrepCategory::DiscordantFragFreq => "DISCORDANT_FRAG_FREQ",
            PrepCategory::DuplicateFreq => "DUPLICATE_FREQ",
            PrepCategory::FragLengthDistribution => "FRAG_LENGTH_DISTRIBUTION",
            PrepCategory::GcBias => "GC_BIAS",
            PrepCategory::MissedGeneVariant => "MISSED_GENE_VARIANT",
            PrepCategory::MsIndelError => "MS_INDEL_ERROR",
            PrepCategory::SummaryTableBamMetrics => "SUMMARY_TABLE_BAM_METRICS",
            PrepCategory::SummaryTablePurple => "SUMMARY_TABLE_PURPLE",
        }
    }
}

/// QseeConfig.parseCategories: every category but those in the comma-separated list (in any case); the
/// two summary-table categories cannot be skipped
pub fn parse_categories(skip: Option<&str>) -> Result<Vec<PrepCategory>> {
    let Some(skip) = skip else { return Ok(PrepCategory::ALL.to_vec()) };
    let mut skipped = Vec::new();
    for s in hmf_common::java_text::split(skip, ',') {
        let upper = s.to_uppercase();
        let c = PrepCategory::ALL
            .into_iter()
            .find(|c| c.name() == upper)
            .ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.qsee.feature.category.PrepCategory.{upper}"))?;
        if c == PrepCategory::SummaryTablePurple || c == PrepCategory::SummaryTableBamMetrics {
            bail!("Skipping category({}) is not supported", c.name());
        }
        skipped.push(c);
    }
    Ok(PrepCategory::ALL.into_iter().filter(|c| !skipped.contains(c)).collect())
}

/// QseeConfig.registerConfig (with QseePlotConfig's and TrainConfig's items when asked)
pub fn register(cb: &mut ConfigBuilder, plot: bool, train: bool) {
    let plurals = "(s), separated by ','";
    cb.add_str("tumor", false, &format!("Tumor ID{plurals}"));
    cb.add_str("reference", false, &format!("Reference ID{plurals}"));
    cb.add_path("sample_id_file", false, "Sample ID TSV file containing columns TumorId and optionally ReferenceId");
    cb.add_path("sample_data_dir", false, "Path to sample pipeline files");
    cb.add_path("bam_metrics_tumor_dir", false, "Path to BamMetrics tumor pipeline files");
    cb.add_path("bam_metrics_ref_dir", false, "Path to BamMetrics reference pipeline files");
    cb.add_path("redux_tumor_dir", false, "Path to Redux tumor pipeline files");
    cb.add_path("redux_ref_dir", false, "Path to Redux reference pipeline files");
    cb.add_path("cobalt_dir", false, "Path to Cobalt pipeline files");
    cb.add_path("esvee_dir", false, "Path to Esvee pipeline files");
    cb.add_path("purple_dir", false, "Path to Purple pipeline files");
    cb.add_str_default("sequencing_type", false, "Sequencing types: ILLUMINA, ULTIMA, SBX", "ILLUMINA");
    cb.add_path("driver_gene_panel", false, "Path to driver gene panel");
    cb.add_path("cohort_percentiles_file", false, "Path to the cohort percentiles file");
    cb.add_path("threshold_overrides_file", false, "Path to threshold overrides file");
    cb.add_str("skip_categories", false, "Comma-separated list of categories to skip");
    cb.add_flag("targeted_mode", "Run in targeted mode (use different threshold defaults for targeted mode)");
    cb.add_flag("allow_missing_input", "Continue sample data extraction even if some input files are missing");
    cb.add_flag("skip_plot", "Write output files only, without plotting");
    cb.add_flag("write_lod_details", "Write exon depth distribution columns to gene LoD file");
    cb.add_path("output_dir", true, "Output directory");
    cb.add_str("output_id", false, "Output directory");
    cb.add_str_default("threads", false, "Number of threads, 0 or 1 not multi-threaded", "1");
    cb.add_logging_options();
    if plot {
        register_plot(cb);
    }
    if train {
        cb.add_integer("num_percentiles", "Number of percentiles (e.g. 5 would give percentiles: 0, 25, 50, 75, 100)", 101);
        cb.add_decimal("percentile_increment", "Percentile increment (e.g. 25 would give percentiles: 0, 25, 50, 75, 100)", f64::NAN);
        cb.add_flag("write_cohort_features", "Write extracted features to file");
    }
}

/// QseePlotConfig.registerConfig
pub fn register_plot(cb: &mut ConfigBuilder) {
    cb.add_path("vis_data_file", false, "Path to the vis data file");
    cb.add_str("tumor", false, "Tumor ID(s), separated by ','");
    cb.add_str("reference", false, "Reference ID(s), separated by ','");
    cb.add_path("sample_id_file", false, "Sample ID CSV file with 'SampleId' column");
    cb.add_flag("show_plot_warnings", "Show ggplot2 warnings");
    cb.add_flag("merge_plots", "In multisample mode, merge plots into one file and delete the individual plot files");
    cb.add_path("cohort_percentiles_file", false, "Path to the cohort percentiles file");
    cb.add_path("output_dir", true, "Output directory");
    cb.add_str("output_id", false, "Output file suffix");
    cb.add_str_default("threads", false, "Number of threads, 0 or 1 not multi-threaded", "1");
    cb.add_logging_options();
}

/// TaskExecutor.parseThreads: Integer.parseInt of the option
pub fn parse_threads(cb: &ConfigBuilder) -> Result<i32> {
    let v = cb.value("threads").unwrap_or_else(|| "1".into());
    hmf_common::java_text::parse_int(&v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
}

/// ConfigUtils.convertWildcardSamplePath: every '*' replaced by the sample id. Java's replaceAll takes the
/// id as a regex replacement, so a '$' or '\' in an id is read as a group reference or an escape (QS11).
pub fn convert_wildcard_sample_path(path: &str, sample_id: &str) -> Result<String> {
    if !path.contains('*') {
        return Ok(path.to_string());
    }
    if compat::java("QS11") && (sample_id.contains('$') || sample_id.contains('\\')) {
        // Matcher.appendReplacement: "\x" is the character x; "$n" a group (none exists) or "${name}"
        let mut rep = String::new();
        let chars: Vec<char> = sample_id.chars().collect();
        let mut i = 0;
        while i < chars.len() {
            match chars[i] {
                '\\' => {
                    i += 1;
                    let c = chars.get(i).ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: character to be escaped is missing"))?;
                    rep.push(*c);
                }
                '$' => match chars.get(i + 1) {
                    // group 0 is the match itself
                    Some('0') => {
                        rep.push('*');
                        i += 1;
                    }
                    Some(c) if c.is_ascii_digit() => bail!("java.lang.IndexOutOfBoundsException: No group {c}"),
                    _ => bail!("java.lang.IllegalArgumentException: Illegal group reference"),
                },
                c => rep.push(c),
            }
            i += 1;
        }
        return Ok(path.replace('*', &rep));
    }
    Ok(path.replace('*', sample_id))
}

/// QseeConfig
#[derive(Clone, Debug)]
pub struct QseeConfig {
    pub samples: SampleIds,
    #[allow(dead_code)] // Java's field: the default of each tool directory
    pub sample_data_dir: String,
    pub bam_metrics_dir_tumor: String,
    pub bam_metrics_dir_ref: String,
    pub redux_dir_tumor: String,
    pub redux_dir_ref: String,
    pub cobalt_dir: String,
    pub esvee_dir: String,
    pub purple_dir: String,
    pub driver_genes: Vec<DriverGene>,
    pub cohort_percentiles_file: Option<String>,
    #[allow(dead_code)] // Java's field: read into `thresholds`
    pub thresholds_file: Option<String>,
    pub thresholds: ThresholdRegistry,
    pub sequencing_type: SequencingType,
    pub categories: Vec<PrepCategory>,
    #[allow(dead_code)] // Java's field: applied in `thresholds`
    pub targeted_mode: bool,
    pub allow_missing_input: bool,
    pub skip_plot: bool,
    pub write_lod_details: bool,
    pub output_dir: String,
    pub output_id: Option<String>,
    pub threads: i32,
}

impl QseeConfig {
    /// new QseeConfig(configBuilder): in Java's order, so that the first error is Java's
    pub fn from(cb: &ConfigBuilder) -> Result<QseeConfig> {
        let samples = SampleIds::from_config(cb)?;
        let sample_data_dir = cb.value_or("sample_data_dir", "");
        let dir = |name: &str| cb.value(name).unwrap_or_else(|| sample_data_dir.clone());
        let bam_metrics_dir_tumor = dir("bam_metrics_tumor_dir");
        let bam_metrics_dir_ref = dir("bam_metrics_ref_dir");
        let redux_dir_tumor = dir("redux_tumor_dir");
        let redux_dir_ref = dir("redux_ref_dir");
        let cobalt_dir = dir("cobalt_dir");
        let esvee_dir = dir("esvee_dir");
        let purple_dir = dir("purple_dir");
        let driver_genes = hmf_common::driver_panel::load_driver_genes(cb.value("driver_gene_panel").as_deref())?;
        let cohort_percentiles_file = cb.value("cohort_percentiles_file");
        let thresholds_file = cb.value("threshold_overrides_file");
        let targeted_mode = cb.has_flag("targeted_mode");
        let thresholds = match &thresholds_file {
            None => ThresholdRegistry::create_default(targeted_mode),
            Some(f) => read_threshold_overrides(f, targeted_mode)?,
        };
        let sequencing_type = SequencingType::value_of(&cb.value("sequencing_type").unwrap_or_default())?;
        let categories = parse_categories(cb.value("skip_categories").as_deref())?;
        let output_dir = check_add_dir_separator(&cb.value("output_dir").unwrap_or_default());
        let output_id = cb.value("output_id");
        let threads = parse_threads(cb)?;
        Ok(QseeConfig {
            samples,
            sample_data_dir,
            bam_metrics_dir_tumor,
            bam_metrics_dir_ref,
            redux_dir_tumor,
            redux_dir_ref,
            cobalt_dir,
            esvee_dir,
            purple_dir,
            driver_genes,
            cohort_percentiles_file,
            thresholds_file,
            thresholds,
            sequencing_type,
            categories,
            targeted_mode,
            allow_missing_input: cb.has_flag("allow_missing_input"),
            skip_plot: cb.has_flag("skip_plot"),
            write_lod_details: cb.has_flag("write_lod_details"),
            output_dir,
            output_id,
            threads,
        })
    }

    pub fn sample_ids(&self, sample_type: SampleType) -> &[String] {
        self.samples.get(sample_type)
    }
    pub fn is_single_patient(&self) -> bool {
        self.samples.is_single_patient()
    }

    /// a tool directory's base, as Java's tool file names join it: an empty directory (no -<tool>_dir
    /// and no -sample_data_dir) is the working directory here, where Java's basePath + "/" + file names
    /// point to the file system root for PURPLE, COBALT, ESVEE and REDUX's BQR and jitter files (QS7)
    fn base(&self, dir: &str, sample_id: &str) -> Result<String> {
        let d = convert_wildcard_sample_path(dir, sample_id)?;
        Ok(if d.is_empty() && !compat::java("QS7") { ".".to_string() } else { d })
    }
    pub fn cobalt_dir(&self, sample_id: &str) -> Result<String> {
        self.base(&self.cobalt_dir, sample_id)
    }
    pub fn esvee_dir(&self, sample_id: &str) -> Result<String> {
        self.base(&self.esvee_dir, sample_id)
    }
    pub fn purple_dir(&self, sample_id: &str) -> Result<String> {
        self.base(&self.purple_dir, sample_id)
    }
    /// getBamMetricsDir; the BamMetrics file names add a separator only where one is missing, so an
    /// empty directory is the working directory in Java too
    pub fn bam_metrics_dir(&self, sample_id: &str, sample_type: SampleType) -> Result<String> {
        let d = if sample_type == SampleType::Tumor { &self.bam_metrics_dir_tumor } else { &self.bam_metrics_dir_ref };
        convert_wildcard_sample_path(d, sample_id)
    }
    pub fn redux_dir(&self, sample_id: &str, sample_type: SampleType) -> Result<String> {
        let d = if sample_type == SampleType::Tumor { &self.redux_dir_tumor } else { &self.redux_dir_ref };
        self.base(d, sample_id)
    }
    /// the REDUX directory as DuplicateFrequency's name joins it (checkAddDirSeparator: "" stays "")
    pub fn redux_dir_unjoined(&self, sample_id: &str, sample_type: SampleType) -> Result<String> {
        let d = if sample_type == SampleType::Tumor { &self.redux_dir_tumor } else { &self.redux_dir_ref };
        convert_wildcard_sample_path(d, sample_id)
    }
}

/// QseePlotConfig
#[derive(Clone, Debug)]
pub struct QseePlotConfig {
    pub vis_data_file: Option<String>,
    pub samples: SampleIds,
    pub cohort_percentiles_file: Option<String>,
    pub show_plot_warnings: bool,
    pub merge_plots: bool,
    pub output_dir: String,
    pub output_id: Option<String>,
    /// Java plots several patients on this many threads; the port runs R for them one after the other
    #[allow(dead_code)]
    pub threads: i32,
}

impl QseePlotConfig {
    pub fn from(cb: &ConfigBuilder) -> Result<QseePlotConfig> {
        let vis_data_file = cb.value("vis_data_file");
        let samples = SampleIds::from_config(cb)?;
        Ok(QseePlotConfig {
            vis_data_file,
            samples,
            cohort_percentiles_file: cb.value("cohort_percentiles_file"),
            show_plot_warnings: cb.has_flag("show_plot_warnings"),
            merge_plots: cb.has_flag("merge_plots"),
            output_dir: check_add_dir_separator(&cb.value("output_dir").unwrap_or_default()),
            output_id: cb.value("output_id"),
            threads: parse_threads(cb)?,
        })
    }
}

/// the log level the R script is given: QC_LOGGER.getLevel() (INFO unless -log_debug or -log_level)
pub fn log_level_name(cb: &ConfigBuilder) -> String {
    if cb.is_registered("log_debug") && cb.has_flag("log_debug") {
        return "DEBUG".into();
    }
    if cb.is_registered("log_level") && cb.has_value("log_level") {
        return cb.value("log_level").unwrap_or_default().to_uppercase();
    }
    "INFO".into()
}

/// Level.valueOf for -log_level: one of log4j's levels, in any case
pub fn check_log_level(cb: &ConfigBuilder) -> Result<()> {
    if cb.is_registered("log_debug") && cb.has_flag("log_debug") {
        return Ok(());
    }
    if cb.is_registered("log_level") && cb.has_value("log_level") {
        let v = cb.value("log_level").unwrap_or_default();
        let levels = ["OFF", "FATAL", "ERROR", "WARN", "INFO", "DEBUG", "TRACE", "ALL"];
        if !levels.contains(&v.to_uppercase().as_str()) {
            bail!("java.lang.IllegalArgumentException: Unknown level constant [{}].", v.to_uppercase());
        }
    }
    Ok(())
}
