//! ChordConfig and PredictConfig: the options of ChordApplication (and ChordDataPrep) and of the
//! predict-only entry point ChordModel, registered and parsed as the Java ConfigBuilder does
//! (hmf_common::config_builder).

use anyhow::{anyhow, bail, Result};

use hmf_common::config_builder::{check_add_dir_separator, ConfigBuilder, ItemType};
use hmf_common::java_text;
use hmf_common::error;

pub const APP_NAME: &str = "CHORD";
pub const VERSION: &str = "2.1.2";

pub const SAMPLE: &str = "sample";
pub const SAMPLE_ID_FILE: &str = "sample_id_file";
pub const PURPLE_DIR: &str = "purple_dir";
pub const SNV_INDEL_VCF_FILE: &str = "snv_indel_vcf_file";
pub const SV_VCF_FILE: &str = "sv_vcf_file";
pub const CHORD_MODEL_FILE: &str = "chord_model_file";
pub const REF_GENOME: &str = "ref_genome";
pub const OUTPUT_DIR: &str = "output_dir";
pub const OUTPUT_ID: &str = "output_id";
pub const THREADS: &str = "threads";
pub const INCLUDE_NON_PASS: &str = "include_non_pass";
pub const WRITE_DETAILED_FILES: &str = "write_detailed_files";

pub const MUT_CONTEXTS_FILE: &str = "mut_contexts_file";
pub const OUTPUT_FILE: &str = "output_file";

/// PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX, PURPLE_SV_VCF_SUFFIX
pub const PURPLE_SOMATIC_VCF_SUFFIX: &str = ".purple.somatic.vcf.gz";
pub const PURPLE_SV_VCF_SUFFIX: &str = ".purple.sv.vcf.gz";

#[derive(Clone, Debug, Default)]
pub struct ChordConfig {
    pub sample_ids: Vec<String>,
    pub purple_dir: Option<String>,
    pub snv_indel_vcf_file: Option<String>,
    pub sv_vcf_file: Option<String>,
    pub chord_model_file: Option<String>,
    pub ref_genome_file: Option<String>,
    /// FileWriterUtils.parseOutputDir: with a trailing separator, or None
    pub output_dir: Option<String>,
    pub output_id: Option<String>,
    pub threads: i32,
    pub include_non_pass: bool,
    pub write_detailed_files: bool,
    /// the log level as log4j holds it (passed to the R script by Java; the port only logs it)
    pub log_level: String,
}

/// ChordConfig.registerConfig
pub fn register_config(cb: &mut ConfigBuilder) {
    cb.add_str(SAMPLE, false, "Sample ID");
    cb.add_path(SAMPLE_ID_FILE, false, "Sample ID CSV file with 'SampleId' column");
    cb.add_path(PURPLE_DIR, false, "Path to Purple pipeline files");
    cb.add_path(SNV_INDEL_VCF_FILE, false, "Path to the VCF containing SNVs and indels");
    cb.add_path(SV_VCF_FILE, false, "Path to the VCF containing structural variants");
    cb.add_path(CHORD_MODEL_FILE, false, "Path to the CHORD.rds file");
    cb.add_str(REF_GENOME, false, "Path to reference genome fasta files");
    cb.add_str(OUTPUT_DIR, false, "Output directory");
    cb.add_str(OUTPUT_ID, false, "Output file suffix");
    // registered as a plain item with a default (not ConfigItemType.INTEGER), as ChordConfig does
    cb.add(ItemType::Str, THREADS, false, "Number of threads, 0 or 1 not multi-threaded", Some("1".to_string()));
    cb.add_flag(INCLUDE_NON_PASS, "Include non pass variants when counting mutation types");
    cb.add_flag(WRITE_DETAILED_FILES, "Write mutation context attributes per variant");
    cb.add_logging_options();
}

/// The log level as Java sets it (ConfigUtils.setLogLevel; log4j2.xml's default is INFO)
pub fn log_level(cb: &ConfigBuilder) -> String {
    if cb.is_registered("log_debug") && cb.has_flag("log_debug") {
        return "DEBUG".to_string();
    }
    if cb.is_registered("log_level") && cb.has_value("log_level") {
        return cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    }
    "INFO".to_string()
}

impl ChordConfig {
    /// the ChordConfig(ConfigBuilder) constructor; Err where Java logs an error and exits 1
    pub fn from(cb: &ConfigBuilder) -> Result<ChordConfig> {
        let sample_ids = sample_ids(cb)?;
        let purple_dir = cb.value(PURPLE_DIR);
        let snv_indel_vcf_file = cb.value(SNV_INDEL_VCF_FILE);
        let sv_vcf_file = cb.value(SV_VCF_FILE);
        check_required_input_paths(purple_dir.as_deref(), snv_indel_vcf_file.as_deref(), sv_vcf_file.as_deref())?;
        let threads_str = cb.value(THREADS).unwrap_or_default();
        let threads = java_text::parse_int(&threads_str)
            .ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{threads_str}\""))?;
        Ok(ChordConfig {
            sample_ids,
            purple_dir,
            snv_indel_vcf_file,
            sv_vcf_file,
            chord_model_file: cb.value(CHORD_MODEL_FILE),
            ref_genome_file: cb.value(REF_GENOME),
            output_dir: cb.value(OUTPUT_DIR).map(|d| if d.is_empty() { d } else { check_add_dir_separator(&d) }),
            output_id: cb.value(OUTPUT_ID),
            threads,
            include_non_pass: cb.has_flag(INCLUDE_NON_PASS),
            write_detailed_files: cb.has_flag(WRITE_DETAILED_FILES),
            log_level: log_level(cb),
        })
    }

    pub fn is_multi_sample(&self) -> bool {
        self.sample_ids.len() > 1
    }

    pub fn is_single_sample(&self) -> bool {
        !self.is_multi_sample()
    }

    /// the output directory as Java concatenates it ("null" when not given)
    pub fn output_dir_str(&self) -> String {
        self.output_dir.clone().unwrap_or_else(|| "null".to_string())
    }

    pub fn snv_indel_vcf_file(&self, sample_id: &str) -> Result<String> {
        match &self.snv_indel_vcf_file {
            Some(f) => convert_wildcard_sample_path(f, sample_id),
            None => convert_wildcard_sample_path(
                &format!("{}/*{}", self.purple_dir.as_deref().unwrap_or("null"), PURPLE_SOMATIC_VCF_SUFFIX),
                sample_id,
            ),
        }
    }

    pub fn sv_vcf_file(&self, sample_id: &str) -> Result<String> {
        match &self.sv_vcf_file {
            Some(f) => convert_wildcard_sample_path(f, sample_id),
            None => convert_wildcard_sample_path(&format!("{}/*{}", self.purple_dir.as_deref().unwrap_or("null"), PURPLE_SV_VCF_SUFFIX), sample_id),
        }
    }
}

/// ChordConfig.getSampleIds
fn sample_ids(cb: &ConfigBuilder) -> Result<Vec<String>> {
    if cb.has_value(SAMPLE) {
        let s = cb.value(SAMPLE).unwrap_or_default();
        return Ok(java_text::split(&s, ',').into_iter().map(str::to_string).collect());
    }
    if cb.has_value(SAMPLE_ID_FILE) {
        return load_sample_ids_file(&cb.value(SAMPLE_ID_FILE).unwrap_or_default());
    }
    error!("Either -{SAMPLE} or -{SAMPLE_ID_FILE} must be provided to config");
    bail!("no sample ids");
}

/// ConfigUtils.loadDelimitedIdFile(file, "SampleId", ","): a 'SampleId' column, or a single column
pub fn load_sample_ids_file(path: &str) -> Result<Vec<String>> {
    let mut lines = match java_text::read_lines(path) {
        Ok(l) => l,
        Err(e) => {
            error!("failed to read SampleId file({path}): {e}");
            return Ok(Vec::new());
        }
    };
    let Some(first) = lines.first().cloned() else {
        bail!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0 (empty sample id file)");
    };
    let column_count = java_text::split_all(&first, ',').len();
    let id_index;
    if first.contains("SampleId") {
        let fields = java_text::fields_index_map(&first, ',');
        id_index = *fields.get("SampleId").ok_or_else(|| anyhow!("java.lang.NullPointerException: no SampleId column in {path}"))?;
        lines.remove(0);
    } else {
        if column_count > 1 {
            return Ok(Vec::new());
        }
        id_index = 0;
    }
    let mut ids = Vec::new();
    for line in &lines {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let items = java_text::split_all(line, ',');
        let id = items.get(id_index).ok_or_else(|| {
            anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {id_index} out of bounds for length {}", items.len())
        })?;
        ids.push(id.to_string());
    }
    Ok(ids)
}

/// ChordConfig.checkRequiredInputPaths
fn check_required_input_paths(purple_dir: Option<&str>, snv_indel: Option<&str>, sv: Option<&str>) -> Result<()> {
    if purple_dir.is_some() || (snv_indel.is_some() && sv.is_some()) {
        return Ok(());
    }
    error!("Please provide either 1) {PURPLE_DIR}, or 2) {SNV_INDEL_VCF_FILE} and {SV_VCF_FILE}");
    bail!("missing input paths");
}

/// ConfigUtils.convertWildcardSamplePath: `samplePath.replaceAll("\\*", sampleId)`, with the sample id
/// taken as a Matcher replacement string ('\' escapes, '$n' group references), as Java does
pub fn convert_wildcard_sample_path(sample_path: &str, sample_id: &str) -> Result<String> {
    if !sample_path.contains('*') {
        return Ok(sample_path.to_string());
    }
    let replacement = java_replacement(sample_id, "*")?;
    Ok(sample_path.replace('*', &replacement))
}

/// Matcher.appendReplacement's reading of a replacement string for a pattern without groups whose match
/// is `group0`
fn java_replacement(replacement: &str, group0: &str) -> Result<String> {
    let chars: Vec<char> = replacement.chars().collect();
    let mut out = String::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c == '\\' {
            i += 1;
            if i == chars.len() {
                bail!("java.lang.IllegalArgumentException: character to be escaped is missing");
            }
            out.push(chars[i]);
            i += 1;
        } else if c == '$' {
            i += 1;
            if i == chars.len() {
                bail!("java.lang.IllegalArgumentException: Illegal group reference: group index is missing");
            }
            if chars[i] == '{' {
                let name: String = chars[i + 1..].iter().take_while(|c| **c != '}').collect();
                bail!("java.lang.IllegalArgumentException: No group with name {{{name}}}");
            }
            let Some(d) = chars[i].to_digit(10) else {
                bail!("java.lang.IllegalArgumentException: Illegal group reference");
            };
            i += 1;
            // the group count is 0, so no further digit extends the reference
            if d != 0 {
                bail!("java.lang.IndexOutOfBoundsException: No group {d}");
            }
            out.push_str(group0);
        } else {
            out.push(c);
            i += 1;
        }
    }
    Ok(out)
}

/// PredictConfig: the options of ChordModel.main (predict only)
pub fn register_predict_config(cb: &mut ConfigBuilder) {
    cb.add_path(MUT_CONTEXTS_FILE, true, "Path to *.chord.mutation_contexts.tsv file");
    cb.add_str(OUTPUT_FILE, true, "Path to output predictions");
    cb.add_path(CHORD_MODEL_FILE, false, "Path to the CHORD.rds file");
}

pub struct PredictConfig {
    pub mut_contexts_file: String,
    pub output_file: String,
    pub chord_model_file: Option<String>,
}

impl PredictConfig {
    pub fn from(cb: &ConfigBuilder) -> PredictConfig {
        PredictConfig {
            mut_contexts_file: cb.value(MUT_CONTEXTS_FILE).unwrap_or_default(),
            output_file: cb.value(OUTPUT_FILE).unwrap_or_default(),
            chord_model_file: cb.value(CHORD_MODEL_FILE),
        }
    }
}
