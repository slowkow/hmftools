//! SigsConfig and VariantFilters: the options, registered in Java's order.

use anyhow::{anyhow, bail, Result};
use hmf_common::config_builder::ConfigBuilder;

pub const SAMPLE: &str = "sample";
pub const SAMPLE_ID_FILE: &str = "sample_id_file";
pub const SOMATIC_VCF_FILE: &str = "somatic_vcf_file";
pub const SAMPLE_COUNTS_FILE: &str = "sample_counts_file";
pub const SIGNATURES_FILE: &str = "signatures_file";
pub const REF_GENOME_VERSION: &str = "ref_genome_version";
pub const OUTPUT_DIR: &str = "output_dir";
pub const OUTPUT_ID: &str = "output_id";
pub const POSITION_BUCKET_SIZE: &str = "position_bucket_size";
pub const MIN_ALLOC: &str = "min_alloc";
pub const MIN_ALLOC_PERC: &str = "min_alloc_perc";
pub const WRITE_POS_COORDS: &str = "write_pos_freq_coords";
pub const SUBCLONAL_MIN: &str = "subclonal_min";
pub const SUBCLONAL_MAX: &str = "subclonal_max";
pub const PLOIDY_MIN: &str = "ploidy_min";
pub const PLOIDY_MAX: &str = "ploidy_max";
pub const TUMOR_ONLY: &str = "tumor_only";

pub const DEFAULT_MIN_ALLOCATION: f64 = 0.01;
pub const DEFAULT_MIN_ALLOCATION_PERC: f64 = 0.0005;
pub const TUMOR_ONLY_AF_THRESHOLD: f64 = 0.05;

/// SigsConfig.registerConfig, then VariantFilters.registerConfig
pub fn register(cb: &mut ConfigBuilder) {
    cb.add_str(SAMPLE, false, "Sample ID");
    cb.add_path(SAMPLE_ID_FILE, false, "Sample ID CSV file with 'SampleId' column");
    cb.add_path(SOMATIC_VCF_FILE, false, "Somatic variant VCF file");
    cb.add_path(SAMPLE_COUNTS_FILE, false, "Path to the main input file");
    cb.add_path(SIGNATURES_FILE, false, "Signature definitions");
    cb.add_str_default(REF_GENOME_VERSION, false, "Ref genome version, 37 or 38", "V37");
    cb.add_str(OUTPUT_DIR, false, "Output directory");
    cb.add_str(OUTPUT_ID, false, "Output file suffix");
    cb.add_logging_options();
    cb.add_integer(POSITION_BUCKET_SIZE, "Position bucket size", 0);
    cb.add_decimal(MIN_ALLOC_PERC, "Min signature allocation as percentage", DEFAULT_MIN_ALLOCATION_PERC);
    cb.add_decimal(MIN_ALLOC, "Min signature allocation", DEFAULT_MIN_ALLOCATION);
    cb.add_flag(WRITE_POS_COORDS, "Include coordinates with pos-frequency counts");
    cb.add_decimal(SUBCLONAL_MIN, "Subclonal minimum", 0.0);
    cb.add_decimal(SUBCLONAL_MAX, "Subclonal minimum", 1.0);
    cb.add_decimal(PLOIDY_MIN, "Ploidy minimum", 0.0);
    cb.add_decimal(PLOIDY_MAX, "Ploidy minimum", 1.0);
    cb.add_flag(TUMOR_ONLY, "Apply tumor-only filters");
}

pub struct SigsConfig {
    /// RefGenomeVersion V38 (else V37, the default)
    pub v38: bool,
    pub sample_ids: Vec<String>,
    pub snv_counts_file: Option<String>,
    pub signatures_file: Option<String>,
    pub vcf_file: Option<String>,
    /// with a trailing separator; None without -output_dir
    pub output_dir: Option<String>,
    pub output_id: Option<String>,
    pub write_pos_freq_coords: bool,
    pub position_bucket_size: i32,
    pub min_fit: f64,
    pub min_fit_perc: f64,
}

/// RefGenomeVersion.from(String)
pub fn parse_ref_genome_version(version: &str) -> Result<bool> {
    match version {
        "V37" | "37" | "HG37" => Ok(false),
        "V38" | "38" | "HG38" => Ok(true),
        _ => bail!("java.lang.IllegalArgumentException: Cannot resolve ref genome version: {version}"),
    }
}

/// ConfigUtils.loadDelimitedIdFile(filename, "SampleId", ","): the SampleId column if the first line holds
/// that text, else the single column of a file of one column; lines that are empty or start with '#' are
/// skipped. A file that cannot be read gives no ids (an error is logged).
pub fn load_sample_ids_file(filename: &str) -> Result<Vec<String>> {
    let lines = match hmf_common::java_text::read_lines(filename) {
        Ok(l) => l,
        Err(e) => {
            hmf_common::error!("failed to read SampleId file({filename}): {e}");
            return Ok(Vec::new());
        }
    };
    let Some(first) = lines.first() else {
        bail!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0 (empty sample ID file {filename})");
    };
    let column_count = hmf_common::java_text::split_all(first, ',').len();
    let (id_index, body) = if first.contains("SampleId") {
        let fields = hmf_common::java_text::fields_index_map(first, ',');
        let i = *fields
            .get("SampleId")
            .ok_or_else(|| anyhow!("java.lang.NullPointerException: no SampleId column in {filename}"))?;
        (i, &lines[1..])
    } else {
        if column_count > 1 {
            return Ok(Vec::new());
        }
        (0, &lines[..])
    };
    let mut ids = Vec::new();
    for line in body {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let items = hmf_common::java_text::split_all(line, ',');
        let id = items
            .get(id_index)
            .ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: line '{line}' of {filename} has no column {id_index}"))?;
        ids.push(id.to_string());
    }
    Ok(ids)
}

impl SigsConfig {
    pub fn from(cb: &ConfigBuilder) -> Result<SigsConfig> {
        let mut sample_ids = Vec::new();
        if cb.has_value(SAMPLE) {
            sample_ids.push(cb.value(SAMPLE).unwrap());
        } else if cb.has_value(SAMPLE_ID_FILE) {
            sample_ids.extend(load_sample_ids_file(&cb.value(SAMPLE_ID_FILE).unwrap())?);
        }
        let v38 = if cb.has_value(REF_GENOME_VERSION) { parse_ref_genome_version(&cb.value(REF_GENOME_VERSION).unwrap())? } else { false };
        Ok(SigsConfig {
            v38,
            sample_ids,
            snv_counts_file: cb.value(SAMPLE_COUNTS_FILE),
            signatures_file: cb.value(SIGNATURES_FILE),
            vcf_file: cb.value(SOMATIC_VCF_FILE),
            output_dir: cb.value(OUTPUT_DIR).map(|d| check_add_dir_separator(&d)),
            output_id: cb.value(OUTPUT_ID),
            write_pos_freq_coords: cb.has_flag(WRITE_POS_COORDS),
            position_bucket_size: cb.integer(POSITION_BUCKET_SIZE),
            min_fit: cb.decimal(MIN_ALLOC),
            min_fit_perc: cb.decimal(MIN_ALLOC_PERC),
        })
    }

    /// The output directory as Java's file names start with it: "null" without -output_dir (SI5), as Java
    /// concatenates the null; corrected, the working directory.
    pub fn output_prefix(&self) -> String {
        match &self.output_dir {
            Some(d) => d.clone(),
            None if crate::compat::java("SI5") => "null".to_string(),
            None => String::new(),
        }
    }
}

/// FileWriterUtils.checkAddDirSeparator
pub fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

/// VariantFilters: the subclonal-likelihood and variant copy number limits are applied only when given on
/// the command line; their registered defaults are never used.
#[derive(Clone, Debug, Default)]
pub struct VariantFilters {
    pub ploidy_min: Option<f64>,
    pub ploidy_max: Option<f64>,
    pub subclonal_min: Option<f64>,
    pub subclonal_max: Option<f64>,
    pub tumor_only: bool,
}

impl VariantFilters {
    pub fn from(cb: &ConfigBuilder) -> VariantFilters {
        let opt = |name: &str| if cb.has_value(name) { Some(cb.decimal(name)) } else { None };
        VariantFilters {
            subclonal_min: opt(SUBCLONAL_MIN),
            subclonal_max: opt(SUBCLONAL_MAX),
            ploidy_min: opt(PLOIDY_MIN),
            ploidy_max: opt(PLOIDY_MAX),
            tumor_only: cb.has_flag(TUMOR_ONLY),
        }
    }
}
