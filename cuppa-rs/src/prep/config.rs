//! PrepConfig, SampleIdsLoader and the categories (CategoryType), with the jar's options and the
//! pipeline's `-virus_dir` (Virus Interpreter's annotated file, as CUPPA 2.5.1 read it; the jar reads
//! Viridian's report from `-viridian_dir`).

use anyhow::{anyhow, bail, Result};
use hmf_common::config_builder::{check_add_dir_separator, ConfigBuilder};
use hmf_common::{debug, error};

pub const APP_NAME: &str = "Cuppa";
pub const VERSION: &str = "2.5";
pub const NO_RNA_SAMPLE_ID: &str = "NO_RNA";
pub const GEN_POS_BUCKET_SIZE: i32 = 500_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum CategoryType {
    Snv,
    Sv,
    SampleTrait,
    Driver,
    GeneExp,
    AltSj,
}

impl CategoryType {
    pub const ALL: [CategoryType; 6] =
        [CategoryType::Snv, CategoryType::Sv, CategoryType::SampleTrait, CategoryType::Driver, CategoryType::GeneExp, CategoryType::AltSj];

    pub fn name(self) -> &'static str {
        match self {
            CategoryType::Snv => "SNV",
            CategoryType::Sv => "SV",
            CategoryType::SampleTrait => "SAMPLE_TRAIT",
            CategoryType::Driver => "DRIVER",
            CategoryType::GeneExp => "GENE_EXP",
            CategoryType::AltSj => "ALT_SJ",
        }
    }
    pub fn is_rna(self) -> bool {
        matches!(self, CategoryType::GeneExp | CategoryType::AltSj)
    }
    fn value_of(s: &str) -> Result<CategoryType> {
        CategoryType::ALL
            .into_iter()
            .find(|c| c.name() == s)
            .ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.cup.prep.CategoryType.{s}"))
    }
}

#[derive(Clone, Debug)]
#[allow(dead_code)] // sample_data_dir: the default of the tool directories, kept as Java keeps it
pub struct PrepConfig {
    pub sample_ids: Vec<String>,
    pub rna_sample_ids: Vec<String>,
    pub sample_data_dir: String,
    pub linx_dir: String,
    pub purple_dir: String,
    /// `-viridian_dir` (or the sample data directory): Viridian's `<sample>.viridian.report.tsv`
    pub viridian_dir: String,
    /// `-virus_dir`: Virus Interpreter's `<sample>.virus.annotated.tsv`, as the pipeline passes it
    pub virus_dir: Option<String>,
    pub isofox_dir: String,
    pub somatic_variants_dir: String,
    pub categories: Vec<CategoryType>,
    pub v38: bool,
    pub alt_sj_sites: Option<String>,
    /// with the trailing separator, as FileWriterUtils.parseOutputDir gives it
    pub output_dir: String,
    pub output_id: Option<String>,
    pub write_by_category: bool,
    pub threads: usize,
}

/// PrepConfig.registerConfig (with PredictionConfig.registerConfig for the PredictionRunner)
pub fn register(cb: &mut ConfigBuilder, runner: bool) {
    cb.add_str("sample", false, "Sample ID");
    cb.add_str("rna_sample", false, "RNA sample ID. Note: Output file prefix uses value from `-sample`");
    cb.add_path(
        "sample_id_file",
        false,
        "TSV file with columns SampleId and/or RnaSampleId. RnaSampleId values can be blank (SampleId value used instead) or NO_RNA (skips extracting RNA features)",
    );
    cb.add_path("sample_data_dir", false, "Path to sample pipeline files");
    cb.add_path("linx_dir", false, "Path to Linx pipeline files");
    cb.add_path("purple_dir", false, "Path to Purple pipeline files");
    cb.add_path("viridian_dir", false, "Path to Viridian pipeline files");
    cb.add_path("virus_dir", false, "Path to Virus Interpreter pipeline files (Cuppa 2.5.1)");
    cb.add_path("isofox_dir", false, "Path to Isofox pipeline files");
    cb.add_path("somatic_variants_dir", false, "Directory containing somatic variant generic tabular files");
    cb.add_str("categories", false, "Categories to build ref data for");
    cb.add_str_default("ref_genome_version", false, "Ref genome version, 37 or 38", "V37");
    cb.add_path("ref_alt_sj_sites", false, "RNA required alternative splice junction sites");
    cb.add_str("output_dir", true, "Output directory");
    cb.add_str("output_id", false, "Output file suffix");
    cb.add_flag("write_by_category", "Cohort mode - write files by category");
    cb.add_str_default("threads", false, "Number of threads to use in multi sample mode", "1");
    if runner {
        cb.add_path("features_path", false, "Path to or directory containing the input features file(s), i.e. *.cuppa_data.tsv file(s)");
        cb.add_path("classifier_path", true, "Path to the CUPPA classifier file (.pickle or .pickle.gz file)");
        // the port needs no interpreter: required and checked only in Java mode (CU5)
        cb.add_path("python_path", crate::compat::java("CU5"), "Path to Python interpreter binary");
        cb.add_str("clf_group", false, "Classifier groups to filter probabilities for. Can be 'all' (default) or 'dna'");
        cb.add_str("cv_predictions_path", false, "Path to a CuppaPrediction tsv file containing the cross-validation predictions");
    }
    cb.add_logging_options();
}

/// FileWriterUtils.readLines then SampleIdsLoader.loadFromLines
fn load_sample_ids_file(path: &str) -> Result<(Vec<String>, Vec<String>)> {
    let mut lines = hmf_common::java_text::read_lines(path).map_err(|e| anyhow!("Failed to load sample IDs file: {e}"))?;
    if lines.is_empty() {
        bail!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0");
    }
    let header = lines.remove(0);
    let fields: Vec<&str> = hmf_common::java_text::split_all(&header, '\t');
    let sample_index = fields.iter().position(|f| *f == "SampleId");
    let rna_index = fields.iter().position(|f| *f == "RnaSampleId");
    if lines.is_empty() {
        error!("No sample IDs found in sample IDs file");
        bail!("no sample ids");
    }
    let mut samples = Vec::new();
    let mut rna = Vec::new();
    for line in &lines {
        let values: Vec<&str> = hmf_common::java_text::split_all(line, '\t');
        let Some(si) = sample_index else { bail!("java.lang.NullPointerException: no SampleId column in {path}") };
        let sample = values.get(si).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: {line}"))?.to_string();
        let mut r = match rna_index {
            Some(ri) => values.get(ri).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: {line}"))?.to_string(),
            None => String::new(),
        };
        if r.is_empty() {
            r = sample.clone();
        }
        samples.push(sample);
        rna.push(r);
    }
    debug!("Loaded {} sample IDs from file", samples.len());
    Ok((samples, rna))
}

/// ConfigUtils.convertWildcardSamplePath: '*' replaced by the sample id as a Java replacement string (the
/// sample ids are checked in PrepConfig::from)
pub fn wildcard(path: &str, sample: &str) -> String {
    hmf_common::config_builder::convert_wildcard_sample_path(path, sample).unwrap_or_else(|_| path.replace('*', sample))
}

impl PrepConfig {
    pub fn from(cb: &ConfigBuilder) -> Result<PrepConfig> {
        let (sample_ids, rna_sample_ids) = if !cb.has_value("sample_id_file") && !cb.has_value("sample") && !cb.has_value("rna_sample") {
            error!("Sample ID(s) must be provided to 1) -sample_id_file or to 2) -sample and/or -rna_sample");
            bail!("no sample ids");
        } else if let Some(f) = cb.value("sample_id_file") {
            load_sample_ids_file(&f)?
        } else {
            // loadFromStrings: a null -sample is added as null (an -rna_sample alone)
            let s = cb.value("sample");
            let r = cb.value("rna_sample").or_else(|| s.clone());
            match s {
                Some(s) => (vec![s], vec![r.unwrap()]),
                None => bail!("java.lang.NullPointerException: no -sample given with -rna_sample"),
            }
        };
        let categories = match cb.value("categories") {
            None => CategoryType::ALL[..4].to_vec(),
            Some(c) => {
                let c = c.to_uppercase();
                match c.as_str() {
                    "ALL" => CategoryType::ALL.to_vec(),
                    "DNA" => CategoryType::ALL[..4].to_vec(),
                    "RNA" => CategoryType::ALL[4..].to_vec(),
                    _ => hmf_common::java_text::split_all(&c, ';').iter().map(|s| CategoryType::value_of(s)).collect::<Result<Vec<_>>>()?,
                }
            }
        };
        let sample_data_dir = cb.value("sample_data_dir").unwrap_or_default();
        let or_data = |n: &str| cb.value(n).unwrap_or_else(|| sample_data_dir.clone());
        let v38 = match cb.value("ref_genome_version") {
            None => false,
            Some(v) => match v.as_str() {
                "V37" | "37" | "HG37" => false,
                "V38" | "38" | "HG38" => true,
                _ => bail!("java.lang.IllegalArgumentException: Cannot resolve ref genome version: {v}"),
            },
        };
        // a sample id that is not a valid Java replacement string fails where a wildcard path is converted
        let dirs = ["sample_data_dir", "linx_dir", "purple_dir", "viridian_dir", "virus_dir", "isofox_dir", "somatic_variants_dir"];
        if dirs.iter().any(|d| cb.value(d).is_some_and(|v| v.contains('*'))) {
            for s in sample_ids.iter().chain(&rna_sample_ids) {
                hmf_common::config_builder::java_replacement(s, "*")?;
            }
        }
        let threads_value = cb.value("threads").unwrap_or_else(|| "1".to_string());
        let threads = hmf_common::java_text::parse_int(&threads_value)
            .ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{threads_value}\""))?;
        Ok(PrepConfig {
            sample_ids,
            rna_sample_ids,
            linx_dir: or_data("linx_dir"),
            purple_dir: or_data("purple_dir"),
            viridian_dir: or_data("viridian_dir"),
            virus_dir: cb.value("virus_dir"),
            isofox_dir: or_data("isofox_dir"),
            somatic_variants_dir: or_data("somatic_variants_dir"),
            sample_data_dir,
            categories,
            v38,
            alt_sj_sites: cb.value("ref_alt_sj_sites"),
            output_dir: check_add_dir_separator(&cb.value("output_dir").unwrap()),
            output_id: cb.value("output_id"),
            write_by_category: cb.has_flag("write_by_category"),
            threads: threads.max(1) as usize,
        })
    }

    pub fn is_single_sample(&self) -> bool {
        self.sample_ids.len() == 1
    }

    pub fn purple_dir(&self, s: &str) -> String {
        wildcard(&self.purple_dir, s)
    }
    pub fn linx_dir(&self, s: &str) -> String {
        wildcard(&self.linx_dir, s)
    }

    /// PurpleCommon.purpleSomaticVcfFile
    pub fn purple_somatic_vcf(&self, s: &str) -> String {
        format!("{}{s}.purple.somatic.vcf.gz", check_add_dir_separator(&self.purple_dir(s)))
    }
    /// SomaticVariant.generateFilename
    pub fn somatic_variants_generic(&self, s: &str) -> String {
        format!("{}{s}.somatic_variants.csv", check_add_dir_separator(&wildcard(&self.somatic_variants_dir, s)))
    }
    pub fn purple_sv_vcf(&self, s: &str) -> String {
        format!("{}{s}.purple.sv.vcf.gz", check_add_dir_separator(&self.purple_dir(s)))
    }
    pub fn purple_purity(&self, s: &str) -> String {
        format!("{}{s}.purple.purity.tsv", check_add_dir_separator(&self.purple_dir(s)))
    }
    /// DriverCatalogFile.generateSomaticFilename (basePath + File.separator + ...)
    pub fn purple_driver_catalog(&self, s: &str) -> String {
        format!("{}/{s}.purple.driver.catalog.somatic.tsv", self.purple_dir(s))
    }
    /// LinxDriver.generateCatalogFilename(dir, sample, true)
    pub fn linx_driver_catalog(&self, s: &str) -> String {
        format!("{}/{s}.linx.driver.catalog.tsv", self.linx_dir(s))
    }
    pub fn linx_clusters(&self, s: &str) -> String {
        format!("{}/{s}.linx.clusters.tsv", self.linx_dir(s))
    }
    pub fn linx_fusions(&self, s: &str) -> String {
        format!("{}/{s}.linx.fusion.tsv", self.linx_dir(s))
    }
    pub fn viridian_report(&self, s: &str) -> String {
        format!("{}{s}.viridian.report.tsv", check_add_dir_separator(&wildcard(&self.viridian_dir, s)))
    }
    pub fn virus_annotated(&self, dir: &str, s: &str) -> String {
        format!("{}{s}.virus.annotated.tsv", check_add_dir_separator(&wildcard(dir, s)))
    }
    pub fn gene_expression(&self, s: &str) -> String {
        format!("{}{s}.isf.gene_data.tsv", check_add_dir_separator(&wildcard(&self.isofox_dir, s)))
    }
    /// PrepConfig.altSpliceJunctionFile: the unfiltered file if it exists
    pub fn alt_splice_junctions(&self, s: &str) -> String {
        let dir = check_add_dir_separator(&wildcard(&self.isofox_dir, s));
        let unfiltered = format!("{dir}{s}.isf.alt_splice_junc_unfiltered.tsv");
        if std::path::Path::new(&unfiltered).exists() { unfiltered } else { format!("{dir}{s}.isf.alt_splice_junc.tsv") }
    }
}
