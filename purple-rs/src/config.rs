//! PurpleConfig, FittingConfig, SomaticFitConfig, ChartConfig and SampleDataFiles, with the options
//! registered as PurpleApplication registers them.

use anyhow::Result;
use hmf_common::config_builder::{check_add_dir_separator, ConfigBuilder, ItemType};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::constants::*;
use crate::types::RunMode;

pub const VERSION: &str = "4.5";
pub const APP_NAME: &str = "Purple";

/// PurpleConstants.AMBIGUOUS_BAF_THRESHOLD: a static that FittingConfig overwrites from -ambiguous_baf_limit
static AMBIGUOUS_BAF_THRESHOLD_BITS: AtomicU64 = AtomicU64::new(0x3FE0CCCCCCCCCCCD); // 0.525

pub fn ambiguous_baf_threshold() -> f64 {
    f64::from_bits(AMBIGUOUS_BAF_THRESHOLD_BITS.load(Ordering::Relaxed))
}
pub fn set_ambiguous_baf_threshold(v: f64) {
    AMBIGUOUS_BAF_THRESHOLD_BITS.store(v.to_bits(), Ordering::Relaxed);
}

#[derive(Clone, Debug)]
pub struct FittingConfig {
    pub min_purity: f64,
    pub max_purity: f64,
    pub purity_increment: f64,
    pub min_ploidy: f64,
    pub max_ploidy: f64,
    pub min_diploid_tumor_ratio_count: i32,
    pub min_diploid_tumor_ratio_count_at_centromere: i32,
    pub ploidy_penalty_factor: f64,
    pub ploidy_penalty_standard_deviation: f64,
    pub ploidy_penalty_min_standard_deviation_per_ploidy: f64,
    pub ploidy_penalty_major_allele_sub_one_multiplier: f64,
    pub ploidy_penalty_major_allele_sub_one_additional: f64,
    pub ploidy_penalty_min_deviation: f64,
    pub deviation_penalty_gc_min_adjust: f64,
    pub gc_ratio_exponent: f64,
}

impl FittingConfig {
    pub fn from(cb: &ConfigBuilder, targeted: bool) -> FittingConfig {
        let c = FittingConfig {
            min_purity: cb.decimal("min_purity"),
            max_purity: cb.decimal("max_purity"),
            purity_increment: cb.decimal("purity_increment"),
            min_ploidy: cb.decimal("min_ploidy"),
            max_ploidy: cb.decimal("max_ploidy"),
            min_diploid_tumor_ratio_count: cb.integer_or(
                "min_diploid_tumor_ratio_count",
                if targeted { TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT } else { MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT },
            ),
            min_diploid_tumor_ratio_count_at_centromere: cb.integer_or(
                "min_diploid_tumor_ratio_count_centromere",
                if targeted { TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT } else { MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT },
            ),
            ploidy_penalty_factor: cb.decimal("ploidy_penalty_factor"),
            ploidy_penalty_standard_deviation: cb.decimal_or(
                "ploidy_penalty_standard_deviation",
                if targeted { TARGETED_PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT } else { PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT },
            ),
            ploidy_penalty_min_standard_deviation_per_ploidy: cb.decimal("ploidy_penalty_min_standard_deviation_per_ploidy"),
            ploidy_penalty_major_allele_sub_one_multiplier: cb.decimal_or(
                "ploidy_penalty_sub_one_major_allele_multiplier",
                PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT,
            ),
            ploidy_penalty_major_allele_sub_one_additional: cb.decimal("ploidy_penalty_sub_min_additional"),
            ploidy_penalty_min_deviation: cb.decimal_or("ploidy_penalty_min", PLOIDY_PENALTY_MIN_DEFAULT),
            deviation_penalty_gc_min_adjust: cb
                .decimal_or("deviation_penalty_gc_min_adjust", if targeted { TARGETED_DEVIATION_PENALTY_GC_MIN_ADJUST_DEFAULT } else { 0.0 }),
            gc_ratio_exponent: cb.decimal_or("gc_ratio_exponent", if targeted { TARGETED_GC_RATIO_EXPONENT_DEFAULT } else { 0.0 }),
        };
        if cb.has_value("ambiguous_baf_limit") {
            set_ambiguous_baf_threshold(cb.decimal("ambiguous_baf_limit"));
        }
        c
    }

    pub fn has_valid_values(&self) -> bool {
        self.min_purity > 0.0
            && self.max_purity >= self.min_purity
            && self.max_purity <= 1.0
            && self.min_ploidy > 0.0
            && self.max_ploidy >= self.min_ploidy
            && self.purity_increment > 0.0
    }

    /// default values, as the Java tests construct FittingConfig
    #[allow(dead_code)]
    pub fn defaults() -> FittingConfig {
        FittingConfig {
            min_purity: MIN_PURITY_DEFAULT,
            max_purity: MAX_PURITY_DEFAULT,
            purity_increment: PURITY_INCREMENT_DEFAULT,
            min_ploidy: MIN_PLOIDY_DEFAULT,
            max_ploidy: MAX_PLOIDY_DEFAULT,
            min_diploid_tumor_ratio_count: MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT,
            min_diploid_tumor_ratio_count_at_centromere: MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT,
            ploidy_penalty_factor: PLOIDY_PENALTY_FACTOR_DEFAULT,
            ploidy_penalty_standard_deviation: PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT,
            ploidy_penalty_min_standard_deviation_per_ploidy: PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT,
            ploidy_penalty_major_allele_sub_one_multiplier: PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT,
            ploidy_penalty_major_allele_sub_one_additional: PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT,
            ploidy_penalty_min_deviation: PLOIDY_PENALTY_MIN_DEFAULT,
            deviation_penalty_gc_min_adjust: 0.0,
            gc_ratio_exponent: 0.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SomaticFitConfig {
    pub min_total_variants: i32,
    pub min_peak_variants: i32,
    pub min_somatic_purity: f64,
    pub min_somatic_purity_spread: f64,
    pub somatic_penalty_weight: f64,
    pub highly_diploid_percentage: f64,
}

#[derive(Clone, Debug)]
pub struct ChartConfig {
    pub disabled: bool,
    pub circos_binary: Option<String>,
    pub plot_directory: String,
    pub circos_directory: String,
    /// run the Circos binary given by -circos, as Java does, instead of drawing the Circos plots natively
    pub external_circos: bool,
    /// keep plot/<tumor>.somatic_data.tsv (the rainfall plot's data, which Java deletes)
    pub write_plot_data: bool,
}

#[derive(Clone, Debug, Default)]
pub struct SampleDataFiles {
    pub sample_data_dir: Option<String>,
    pub somatic_sv_vcf_file: Option<String>,
    pub germline_sv_vcf_file: Option<String>,
    pub somatic_vcf_file: Option<String>,
    pub germline_vcf_file: Option<String>,
    pub amber_directory: Option<String>,
    pub cobalt_directory: Option<String>,
    pub redux_tumor_directory: Option<String>,
}

#[derive(Clone, Debug)]
pub struct PurpleConfig {
    pub reference_id: Option<String>,
    pub tumor_id: String,
    pub version: String,
    pub output_dir: String,
    pub run_drivers: bool,
    pub sample_files: SampleDataFiles,
    pub fitting: FittingConfig,
    pub somatic_fitting: SomaticFitConfig,
    pub charting: ChartConfig,
    pub target_regions_mode: bool,
    pub ignore_plot_errors: bool,
    pub sv_qual_filter: i32,
    pub threads: i32,
    pub write_all_somatics: bool,
    pub skip_circos_genes: bool,
    pub specific_chromosomes: Vec<String>,
    pub specific_regions: Vec<(String, i32, i32)>,
}

impl PurpleConfig {
    pub fn tumor_only_mode(&self) -> bool {
        self.reference_id.is_none()
    }
    pub fn germline_mode(&self) -> bool {
        self.reference_id.as_deref() == Some(self.tumor_id.as_str())
    }
    pub fn run_tumor(&self) -> bool {
        !self.germline_mode()
    }
    pub fn run_germline(&self) -> bool {
        !self.tumor_only_mode()
    }
    pub fn run_mode(&self) -> RunMode {
        if self.tumor_only_mode() {
            RunMode::Tumor
        } else if self.germline_mode() {
            RunMode::Germline
        } else {
            RunMode::TumorGermline
        }
    }

    /// PurpleConfig.excludeOnSpecificRegion
    pub fn exclude_on_specific_region(&self, chromosome: &str, position: i32) -> bool {
        if self.specific_chromosomes.is_empty() && self.specific_regions.is_empty() {
            return false;
        }
        if !self.specific_regions.is_empty() {
            return !self.specific_regions.iter().any(|(c, s, e)| c == chromosome && position >= *s && position <= *e);
        }
        !self.specific_chromosomes.iter().any(|c| c == chromosome)
    }
}

/// The options PurpleConfig.registerConfig and addLoggingOptions register, in that order.
pub fn register(cb: &mut ConfigBuilder) {
    cb.add_str("reference", false, "Reference ID");
    cb.add_str("tumor", false, "Tumor ID");
    cb.add_str("output_dir", false, "Path to the output directory. If <sample_dir> is set, then is sample_dir/output_dir/.");
    cb.add_flag("write_all_somatics", "Write all variants regardless of filters");
    cb.add_flag("ignore_plot_errors", "Run to completion if plotting fails");
    cb.add_flag("skip_circos_genes", "Skip plotting driver genes on Circos");
    cb.add_integer("sv_qual_filter", "Min SV qual, 0 is disabled", 0);

    // FittingConfig
    cb.add_decimal("min_purity", "Minimum purity", MIN_PURITY_DEFAULT);
    cb.add_decimal("max_purity", "Maximum purity", MAX_PURITY_DEFAULT);
    cb.add_decimal("purity_increment", "Purity increment", PURITY_INCREMENT_DEFAULT);
    cb.add_decimal("min_ploidy", "Minimum ploidy", MIN_PLOIDY_DEFAULT);
    cb.add_decimal("max_ploidy", "Maximum ploidy", MAX_PLOIDY_DEFAULT);
    let targeted_int = |cb: &mut ConfigBuilder, name: &str, desc: &str, d: i32, t: i32| {
        cb.add(ItemType::Integer, name, false, &format!("{desc}, default={d} targeted default={t}"), Some(d.to_string()));
    };
    targeted_int(
        cb,
        "min_diploid_tumor_ratio_count",
        "Minimum ratio count while smoothing before diploid regions become suspect",
        MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT,
        TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT,
    );
    targeted_int(
        cb,
        "min_diploid_tumor_ratio_count_centromere",
        "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere",
        MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT,
        TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT,
    );
    cb.add_decimal("ploidy_penalty_factor", "Penalty factor to apply to the number of copy number events", PLOIDY_PENALTY_FACTOR_DEFAULT);
    cb.add_decimal(
        "ploidy_penalty_standard_deviation",
        "Standard deviation of normal distribution modelling ploidy deviation from whole number",
        PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT,
    );
    cb.add_decimal(
        "ploidy_penalty_min_standard_deviation_per_ploidy",
        "Minimum ploidy penalty standard deviation to be applied",
        PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT,
    );
    cb.add_decimal(
        "ploidy_penalty_sub_min_additional",
        "Additional penalty to apply to major allele < 1 or minor allele < 0",
        PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT,
    );
    cb.add_decimal(
        "ploidy_penalty_sub_one_major_allele_multiplier",
        "Penalty multiplier applied to major allele < 1",
        PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT,
    );
    cb.add_decimal("ploidy_penalty_min", "Minimum ploidy penalty", PLOIDY_PENALTY_MIN_DEFAULT);
    cb.add_decimal("ambiguous_baf_limit", "Amber ambiguous BAF limit", DEFAULT_AMBIGUOUS_BAF_THRESHOLD);
    let targeted_dec = |cb: &mut ConfigBuilder, name: &str, desc: &str, d: f64, t: f64| {
        let desc = format!("{desc}, default={} targeted default={}", java_g3(d), java_g3(t));
        cb.add(ItemType::Decimal, name, false, &desc, Some(hmf_common::fmt::java_double_to_string(d)));
    };
    targeted_dec(cb, "deviation_penalty_gc_min_adjust", "Adjust deviation penalty by tumor GC Ratio", 0.0, TARGETED_DEVIATION_PENALTY_GC_MIN_ADJUST_DEFAULT);
    targeted_dec(cb, "gc_ratio_exponent", "Adjust GC Ratio by exponent in penalty calc", 0.0, TARGETED_GC_RATIO_EXPONENT_DEFAULT);

    // SomaticFitConfig
    cb.add_integer("somatic_min_peak", "Minimum number of somatic variants to consider a peak", SOMATIC_MIN_PEAK_DEFAULT);
    cb.add_integer("somatic_min_variants", "Minimum number of somatic variants required to assist highly diploid fits", SOMATIC_MIN_VARIANTS_DEFAULT);
    cb.add_decimal("somatic_min_purity", "Minimum spread within candidate purities before somatics can be used", SOMATIC_MIN_PURITY_DEFAULT);
    cb.add_decimal(
        "somatic_min_purity_spread",
        "Minimum spread within candidate purities before somatics can be used",
        SOMATIC_MIN_PURITY_SPREAD_DEFAULT,
    );
    cb.add_decimal("somatic_penalty_weight", "Proportion of somatic deviation to include in fitted purity score", SOMATIC_PENALTY_WEIGHT_DEFAULT);
    cb.add_decimal("highly_diploid_percentage", "Proportion of genome that must be diploid before using somatic fit", HIGHLY_DIPLOID_PERCENTAGE_DEFAULT);

    // ReferenceData
    cb.add_str_default("ref_genome_version", false, "Ref genome version, 37 or 38", "V37");
    cb.add_path("ref_genome", false, "Path to reference genome fasta files");
    cb.add_str_default("somatic_hotspots", false, "Path to somatic hotspot VCF", "");
    cb.add_str_default("germline_hotspots", false, "Path to germline hotspot VCF", "");
    cb.add_path("germline_amp_del_freq_file", false, "Path to cohort germline deletions frequency file");
    cb.add_path("target_regions_bed", false, "Target regions BED file");
    cb.add_integer("target_coding_factor", "Targeted panel coding base factor", DEFAULT_CODING_BASE_FACTOR);
    cb.add_decimal("target_tmb_ratio", "Targeted panel TMB adjustment factor", DEFAULT_TARGETED_TMB_RATIO);
    cb.add_decimal("target_tml_ratio", "Targeted panel TML adjustment factor", DEFAULT_TARGETED_TML_RATIO);
    cb.add_path("ensembl_data_dir", true, "Ensembl data file directory");
    cb.add_path("driver_gene_panel", false, "Path to driver gene panel");

    // ChartConfig
    cb.add_path("circos", false, "Location of circos binary");
    cb.add_flag("no_charts", "Disable charts");
    cb.add_flag("external_circos", "Run the Circos binary given by -circos instead of drawing the Circos plots natively");
    cb.add_flag("circos_native", "Draw the Circos plots natively, as -circos does, without naming a Circos binary");
    cb.add_flag("write_plot_data", "Keep the rainfall plot's data, plot/<tumor>.somatic_data.tsv");

    // SampleDataFiles
    cb.add_path("sample_dir", false, "Sample data directory containing Cobalt, Amber, Pave and Esvee files");
    cb.add_path("cobalt_dir", false, "Path to Cobalt pipeline files");
    cb.add_path("cobalt", false, "Path to Cobalt pipeline files deprecated, use: cobalt_dir");
    cb.add_path("amber", false, "Path to Amber pipeline files deprecated, use: amber_dir");
    cb.add_path("amber_dir", false, "Path to Amber pipeline files");
    cb.add_path("esvee_dir", false, "Path to Esvee pipeline files");
    cb.add_path("pave_somatic_dir", false, "Path to Pave somatic pipeline files");
    cb.add_path("pave_germline_dir", false, "Path to Pave germline pipeline files");
    cb.add_path("somatic_sv_vcf", false, "Somatic SV VCF");
    cb.add_path("germline_sv_vcf", false, "Germline SV VCF to annotate");
    cb.add_path("germline_vcf", false, "Germline variant VCF");
    cb.add_path("somatic_vcf", false, "Somatic variant VCF");
    cb.add_path("redux_tumor_dir", false, "Path to Redux tumor pipeline files");
    cb.add_str(
        "pipeline_format",
        false,
        "Pipeline output structure: OA_V3_1, OA_V3_0, OA_V2_3, OA_V2_2, OA_V2_0, PIP5_V6_0, DB_V6_0. Default: OA_V3_1",
    );
    cb.add_path("pipeline_format_file", false, "File describing expected tool directory structure.");

    cb.add_thread_options(1);
    cb.add_str("specific_chr", false, "Restrict to chromosome(s) separated by ';'");
    cb.add_str(
        "specific_regions",
        false,
        "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd",
    );
    cb.add_logging_options();
}

/// String.format("%.3g", v) for the two defaults FittingConfig prints (0 and positive values below 1000)
fn java_g3(v: f64) -> String {
    if v == 0.0 {
        return "0.00".into();
    }
    let exp = v.abs().log10().floor() as i32;
    let decimals = (2 - exp).max(0) as usize;
    hmf_common::fmt::java_fixed(v, decimals)
}

fn create_directory(dir: &str) -> bool {
    let p = std::path::Path::new(dir);
    if !p.exists() && std::fs::create_dir_all(p).is_err() {
        crate::error!("unable to create chart directory {dir}");
        return false;
    }
    true
}

/// SampleDataFiles.getFilename
fn sample_file(cb: &ConfigBuilder, file_cfg: &str, dir_cfg: &str, pipeline_dir: &str, sample_dir: Option<&str>, sample: &str, suffix: &str) -> Option<String> {
    if cb.has_value(file_cfg) {
        let f = cb.value(file_cfg).unwrap();
        if std::path::Path::new(&f).exists() {
            return Some(f);
        }
        crate::error!("missing file: {f}");
        return None;
    }
    let base = if cb.has_value(dir_cfg) {
        check_add_dir_separator(&cb.value(dir_cfg).unwrap())
    } else if let Some(sd) = sample_dir {
        check_add_dir_separator(sd) + &check_add_dir_separator(pipeline_dir)
    } else {
        return Some(String::new());
    };
    let filename = format!("{base}{sample}{suffix}");
    let exists = |f: &str| std::path::Path::new(f).exists();
    if exists(&filename) {
        return Some(filename);
    }
    for sub in ["/germline", "/somatic"] {
        if pipeline_dir.contains(sub) {
            let alt = filename.replacen(sub, "", 1);
            if exists(&alt) {
                return Some(alt);
            }
        }
    }
    let filename = format!("{}{sample}{suffix}", sample_dir.unwrap_or("null"));
    Some(if exists(&filename) { filename } else { String::new() })
}

impl SampleDataFiles {
    pub fn from(cb: &ConfigBuilder, sample: &str) -> Result<SampleDataFiles> {
        let sample_dir = if cb.has_value("sample_dir") { Some(check_add_dir_separator(&cb.value("sample_dir").unwrap())) } else { None };
        if cb.has_value("pipeline_format_file") || cb.value("pipeline_format").is_some_and(|f| f != "OA_V3_1" && f != "OA_V3_0") {
            anyhow::bail!("only the default pipeline format (OA_V3_1) is supported by this port");
        }
        let dir = |primary: &str, deprecated: Option<&str>, pipeline: &str| -> Option<String> {
            if cb.has_value(primary) {
                Some(check_add_dir_separator(&cb.value(primary).unwrap()))
            } else if deprecated.is_some_and(|d| cb.has_value(d)) {
                Some(check_add_dir_separator(&cb.value(deprecated.unwrap()).unwrap()))
            } else {
                sample_dir.as_ref().map(|sd| format!("{sd}{pipeline}/"))
            }
        };
        let sd = sample_dir.as_deref();
        Ok(SampleDataFiles {
            amber_directory: dir("amber_dir", Some("amber"), "amber"),
            cobalt_directory: dir("cobalt_dir", Some("cobalt"), "cobalt"),
            somatic_sv_vcf_file: sample_file(cb, "somatic_sv_vcf", "esvee_dir", "esvee", sd, sample, ".esvee.somatic.vcf.gz"),
            germline_sv_vcf_file: sample_file(cb, "germline_sv_vcf", "esvee_dir", "esvee", sd, sample, ".esvee.germline.vcf.gz"),
            somatic_vcf_file: sample_file(cb, "somatic_vcf", "pave_somatic_dir", "pave/somatic", sd, sample, ".pave.somatic.vcf.gz"),
            germline_vcf_file: sample_file(cb, "germline_vcf", "pave_germline_dir", "pave/germline", sd, sample, ".pave.germline.vcf.gz"),
            redux_tumor_directory: cb.value("redux_tumor_dir"),
            sample_data_dir: sample_dir,
        })
    }
}

impl PurpleConfig {
    /// PurpleConfig(version, configBuilder); the bool is mIsValid
    pub fn from(cb: &ConfigBuilder) -> Result<(PurpleConfig, bool)> {
        let mut valid = true;
        let reference_id = cb.value("reference");
        let tumor_id = cb.value("tumor").or_else(|| reference_id.clone());
        let Some(tumor_id) = tumor_id else {
            return Ok((Self::empty(), false));
        };
        let output_dir_raw = cb.value("output_dir").unwrap_or_else(|| "null".to_string());
        let output_dir = if cb.has_value("sample_dir") {
            check_add_dir_separator(&(check_add_dir_separator(&cb.value("sample_dir").unwrap()) + &output_dir_raw))
        } else {
            check_add_dir_separator(&output_dir_raw)
        };
        let sample_files = SampleDataFiles::from(cb, &tumor_id)?;
        let charting = ChartConfig {
            disabled: cb.has_flag("no_charts"),
            plot_directory: format!("{output_dir}plot/"),
            circos_directory: format!("{output_dir}circos/"),
            // -circos_native: the Circos plots without a Circos binary (the binary is never run unless -external_circos)
            circos_binary: cb.value("circos").or_else(|| cb.has_flag("circos_native").then(|| "native".to_string())),
            external_circos: cb.has_flag("external_circos"),
            write_plot_data: cb.has_flag("write_plot_data"),
        };
        valid &= create_directory(&output_dir);
        if !charting.disabled {
            if charting.circos_binary.is_some() {
                valid &= create_directory(&charting.circos_directory);
            }
            valid &= create_directory(&charting.plot_directory);
        }
        let target_regions_mode = cb.has_value("target_regions_bed");
        let fitting = FittingConfig::from(cb, target_regions_mode);
        let somatic_fitting = SomaticFitConfig {
            min_total_variants: cb.integer("somatic_min_variants"),
            min_peak_variants: cb.integer("somatic_min_peak"),
            min_somatic_purity: cb.decimal("somatic_min_purity"),
            min_somatic_purity_spread: cb.decimal("somatic_min_purity_spread"),
            somatic_penalty_weight: cb.decimal("somatic_penalty_weight"),
            highly_diploid_percentage: cb.decimal("highly_diploid_percentage"),
        };
        let specific = hmf_common::regions::SpecificRegions::from_java(cb.value("specific_regions").as_deref(), cb.value("specific_chr").as_deref());
        let (specific_chromosomes, specific_regions) = match specific {
            Ok(s) => (s.chromosomes.clone(), s.regions.iter().map(|r| (r.chromosome.clone(), r.start, r.end)).collect()),
            Err(e) => {
                crate::error!("{e}");
                valid = false;
                (Vec::new(), Vec::new())
            }
        };
        Ok((
            PurpleConfig {
                reference_id,
                tumor_id,
                version: VERSION.to_string(),
                output_dir,
                run_drivers: cb.has_value("driver_gene_panel"),
                sample_files,
                fitting,
                somatic_fitting,
                charting,
                target_regions_mode,
                ignore_plot_errors: cb.has_flag("ignore_plot_errors"),
                sv_qual_filter: cb.integer("sv_qual_filter"),
                threads: cb.integer("threads"),
                write_all_somatics: cb.has_flag("write_all_somatics"),
                skip_circos_genes: cb.has_flag("skip_circos_genes"),
                specific_chromosomes,
                specific_regions,
            },
            valid,
        ))
    }

    fn empty() -> PurpleConfig {
        PurpleConfig {
            reference_id: None,
            tumor_id: String::new(),
            version: VERSION.into(),
            output_dir: String::new(),
            run_drivers: false,
            sample_files: SampleDataFiles::default(),
            fitting: FittingConfig::defaults(),
            somatic_fitting: SomaticFitConfig {
                min_total_variants: 0,
                min_peak_variants: 0,
                min_somatic_purity: 0.0,
                min_somatic_purity_spread: 0.0,
                somatic_penalty_weight: 0.0,
                highly_diploid_percentage: 0.0,
            },
            charting: ChartConfig { disabled: true, circos_binary: None, plot_directory: String::new(), circos_directory: String::new(), external_circos: false, write_plot_data: false },
            target_regions_mode: false,
            ignore_plot_errors: false,
            sv_qual_filter: 0,
            threads: 1,
            write_all_somatics: false,
            skip_circos_genes: false,
            specific_chromosomes: Vec::new(),
            specific_regions: Vec::new(),
        }
    }
}
