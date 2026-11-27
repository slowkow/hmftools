//! NeoScorer (`com.hartwig.hmftools.neo.score`): scores NEO_FINDER's neo-epitopes.
//!
//! For each sample: the neo-epitopes (`<sample>.neo.neo_data.tsv`), LILAC's alleles and PURPLE's ploidy;
//! TPMs from ISOFOX (with RNA) or the cohort TPM medians; RNA fragment support from ISOFOX's neo-epitope file
//! (fusions) and an RNA-annotated somatic VCF (point mutations). TpmCalculator derives the peptides of each
//! neo-epitope (8 to 12 amino acids, with 3 flanking ones) and their expected and effective TPMs;
//! peptides found in an upstream transcript's wild-type protein are dropped, and the rest are scored for each
//! of the sample's alleles (BindScorer). Writes `<sample>.neo.peptide_scores.tsv` and
//! `<sample>.neo.neoepitope.tsv` (`neo_cohort.*` for several samples).

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;
use std::io::Write;

use hmf_common::config_builder::ConfigBuilder;
use hmf_common::ensembl::{check_old_filename, ensembl_amino_acid_data_file, load_transcript_amino_acid_data};
use hmf_common::fmt::{java_fixed, java_sci};
use hmf_common::hla::{lilac_allele_filename, read_lilac_alleles};
use hmf_common::java_hash::JavaHashSet;
use hmf_common::java_text::{fields_index_map, parse_bool, parse_double, parse_int, read_lines, split, split_all};
use hmf_common::neo::{
    check_add_dir_separator, extract_transcript_names, neo_epitope_filename, parse_neo_epitope_type, point_mutation_info,
    read_rna_neo_epitopes, rna_neo_epitope_filename, NeoEpitopeType, RnaNeoEpitope,
};
use hmf_common::purple_files::read_purity_context;
use hmf_common::regions::RefGenomeVersion;
use hmf_common::rna::{read_transcript_expression, transcript_expression_filename, RnaExpressionMatrix, EXPRESSION_SCOPE_TRANS, INVALID_EXP};
use hmf_common::variant_context::{Value, VcfReader};
use hmf_common::{debug, error, info, warn};

use crate::bind::*;
use crate::compat;
use crate::epitope::generate_peptides;

pub const FS_UP: usize = 0;
pub const FS_DOWN: usize = 1;

/// NeoRnaData.NO_TPM_VALUE
pub const NO_TPM_VALUE: f64 = INVALID_EXP;

/// TpmCalculator.DEFAULT_PEPTIDE_LENGTH_RANGE
pub const DEFAULT_PEPTIDE_LENGTH_RANGE: [i32; 2] = [MIN_PEPTIDE_LENGTH, REF_PEPTIDE_LENGTH];

/// NeoScorerConfig.addConfig
pub fn register(cb: &mut ConfigBuilder) {
    cb.add_path("sample_id_file", false, "Sample ID CSV file with 'SampleId' column");
    cb.add_str("sample", false, "Sample ID");
    cb.add_str("cancer_type", false, "Cancer type for sample, matching those in cohort TPM medians file");
    cb.add_path("sample_data_dir", false, "Path to sample pipeline files");
    cb.add_path("neo_dir", false, "Path to Neo pipeline files");
    cb.add_path("purple_dir", false, "Path to Purple pipeline files");
    cb.add_path("lilac_dir", false, "Path to Lilac pipeline files");
    cb.add_path("rna_somatic_vcf", false, "Directory for Purple somatic variant RNA-appended files");
    cb.add_str("rna_sample", false, "RNA sample ID in Sage-appended VCF");
    cb.add_path("isofox_dir", false, "Directory for Isofox files (Transcript expresion, Neoepitope coverage)");
    cb.add_path("cancer_tpm_medians_file", false, "TPM medians per cancer type and pan-cancer");
    cb.add_path("cohort_trans_exp_file", false, "Cohort transcript expression matrix, instead of loading per-sample Isofox transcript data");
    cb.add_str("peptide_lengths", false, "Peptide length min-max, separated by '-', eg 8-12");
    register_score_config(cb);
    cb.add_logging_options();
    cb.add_path("ensembl_data_dir", false, "Ensembl data file directory");
    cb.add_str_default("ref_genome_version", false, "Ref genome version, 37 or 38", "V37");
    cb.add_decimal("rank_threshold", "Rank threshold to write full peptide data, default 0 (not applied)", 0.0);
    cb.add_decimal("sim_threshold", "Immunogenic similarity threshold to write full peptide data", 10.0);
    cb.add_str("write_types", false, "Valid types: ALLELE_PEPTIDE, PEPTIDE, NEOEPITOPE, SAMPLE_SUMMARY sep by ';'");
    cb.add_str("output_dir", false, "Output directory");
    cb.add_str("output_id", false, "Output file suffix");
    cb.add_thread_options(1);
}

// ---------------------------------------------------------------------------------------------
// config and samples

/// OutputType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputType {
    SampleSummary,
    Neoepitope,
    Peptide,
    AllelePeptide,
}

impl OutputType {
    fn value_of(s: &str) -> Result<OutputType> {
        Ok(match s {
            "SAMPLE_SUMMARY" => OutputType::SampleSummary,
            "NEOEPITOPE" => OutputType::Neoepitope,
            "PEPTIDE" => OutputType::Peptide,
            "ALLELE_PEPTIDE" => OutputType::AllelePeptide,
            _ => bail!("IllegalArgumentException: No enum constant com.hartwig.hmftools.neo.score.OutputType.{s}"),
        })
    }
}

/// TpmSource
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TpmSource {
    Sample,
    CancerType,
    Cohort,
}

impl TpmSource {
    fn name(self) -> &'static str {
        match self {
            TpmSource::Sample => "SAMPLE",
            TpmSource::CancerType => "CANCER_TYPE",
            TpmSource::Cohort => "COHORT",
        }
    }
}

/// SampleData
#[derive(Clone, Debug)]
pub struct SampleData {
    pub tumor_id: String,
    /// null in Java for samples from a sample id file; file names then carry "null"
    pub rna_sample_id: Option<String>,
    pub cancer_type: String,
    pub has_rna: bool,
}

impl SampleData {
    fn rna_id_str(&self) -> &str {
        self.rna_sample_id.as_deref().unwrap_or("null")
    }
}

/// SampleData.loadFromConfig: the sample id file (SampleId, optional CancerType and HasRna columns), else
/// -sample
pub fn load_samples(cb: &ConfigBuilder) -> Result<Vec<SampleData>> {
    let mut samples = Vec::new();
    if let Some(filename) = cb.value("sample_id_file") {
        let lines = match read_lines(&filename) {
            Ok(l) => l,
            Err(e) => {
                error!("failed to read sample ID file({}): {}", filename, e);
                return Ok(samples);
            }
        };
        let Some(header) = lines.first() else { return Ok(samples) };
        let f = fields_index_map(header, ',');
        let sample_i = *f.get("SampleId").ok_or_else(|| anyhow!("NullPointerException: no SampleId column in {filename}"))?;
        let cancer_i = f.get("CancerType").copied();
        let rna_i = f.get("HasRna").copied();
        for line in &lines[1..] {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let v = split_all(line, ',');
            let g = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("ArrayIndexOutOfBoundsException: Index {i} in {filename}"));
            let tumor_id = g(sample_i)?.to_string();
            let cancer_type = match cancer_i {
                Some(i) => g(i)?.to_string(),
                None => String::new(),
            };
            let has_rna = match rna_i {
                Some(i) => parse_bool(g(i)?),
                None => false,
            };
            samples.push(SampleData { tumor_id, rna_sample_id: None, cancer_type, has_rna });
        }
    } else if let Some(sample) = cb.value("sample") {
        samples.push(SampleData {
            tumor_id: sample,
            rna_sample_id: cb.value("rna_sample"),
            cancer_type: cb.value_or("cancer_type", ""),
            has_rna: cb.has_value("isofox_dir") && cb.has_value("rna_somatic_vcf"),
        });
    }
    Ok(samples)
}

/// NeoScorerConfig
pub struct NeoScorerConfig {
    pub samples: Vec<SampleData>,
    pub neo_dir: String,
    pub lilac_dir: String,
    pub purple_dir: String,
    pub isofox_dir: String,
    pub rna_somatic_vcf: String,
    pub output_dir: String,
    pub output_id: Option<String>,
    pub write_types: Vec<OutputType>,
    pub likelihood_threshold: f64,
    pub similarity_threshold: f64,
    pub threads: i32,
    pub peptide_length_range: [i32; 2],
}

/// `-peptide_lengths min-max`
fn parse_peptide_lengths(s: &str) -> Result<[i32; 2]> {
    let parts: Vec<&str> = s.split('-').collect();
    let bad = || anyhow!("invalid peptide_lengths({s}): expected min-max, eg 8-12");
    if parts.len() != 2 {
        return Err(bad());
    }
    let (a, b) = (parse_int(parts[0].trim()).ok_or_else(bad)?, parse_int(parts[1].trim()).ok_or_else(bad)?);
    if a < 1 || b < a {
        return Err(bad());
    }
    Ok([a, b])
}

impl NeoScorerConfig {
    pub fn new(cb: &ConfigBuilder) -> Result<NeoScorerConfig> {
        let samples = load_samples(cb)?;
        let sample_data_dir = cb.value("sample_data_dir").map(|d| check_add_dir_separator(&d)).unwrap_or_default();
        let dir = |name: &str| check_add_dir_separator(&cb.value(name).unwrap_or_else(|| sample_data_dir.clone()));
        let output_dir = match cb.value("output_dir") {
            Some(d) => check_add_dir_separator(&d),
            None => sample_data_dir.clone(),
        };
        let mut write_types = Vec::new();
        match cb.value("write_types") {
            Some(t) => {
                for x in split(&t, ';') {
                    write_types.push(OutputType::value_of(x)?);
                }
            }
            None => {
                write_types.push(OutputType::Neoepitope);
                write_types.push(OutputType::AllelePeptide);
            }
        }
        // NE9: Java registers -peptide_lengths but never reads it
        let peptide_length_range = match cb.value("peptide_lengths") {
            Some(v) if !compat::java("NE9") => parse_peptide_lengths(&v)?,
            _ => DEFAULT_PEPTIDE_LENGTH_RANGE,
        };
        Ok(NeoScorerConfig {
            samples,
            neo_dir: dir("neo_dir"),
            lilac_dir: dir("lilac_dir"),
            purple_dir: dir("purple_dir"),
            isofox_dir: dir("isofox_dir"),
            // NE11: Java defaults the RNA somatic VCF to the sample data directory, which it then fails to read
            // as a VCF (a null variant list, then a NullPointerException for any point mutation); corrected, the
            // VCF is read only when given (as HasRna already requires)
            rna_somatic_vcf: cb.value("rna_somatic_vcf").unwrap_or_else(|| if compat::java("NE11") { sample_data_dir.clone() } else { String::new() }),
            output_dir,
            output_id: cb.value("output_id"),
            write_types,
            likelihood_threshold: cb.decimal("rank_threshold"),
            similarity_threshold: cb.decimal("sim_threshold"),
            threads: cb.integer("threads"),
            peptide_length_range,
        })
    }

    /// formFilename
    pub fn form_filename(&self, file_id: &str) -> String {
        let mut f = self.output_dir.clone();
        if self.samples.len() == 1 {
            f += &format!("{}.neo", self.samples[0].tumor_id);
        } else {
            f += "neo_cohort";
        }
        if let Some(id) = self.output_id.as_deref().filter(|s| !s.is_empty()) {
            f += &format!(".{id}");
        }
        f += &format!(".{file_id}.tsv");
        f
    }
}

// ---------------------------------------------------------------------------------------------
// reference data

/// TpmMediansCache: pan-cancer ("All") and cancer-type TPM medians by transcript. Only the columns the
/// samples ask for are parsed (Java parses every column of every line).
pub struct TpmMediansCache {
    cancer_types: std::collections::HashSet<String>,
    /// transcript -> the values of `columns`
    values: HashMap<String, Vec<f64>>,
    columns: Vec<String>,
}

const ALL_TYPES: &str = "All";

impl TpmMediansCache {
    pub fn new(filename: Option<&str>, wanted_cancer_types: &[String]) -> Result<TpmMediansCache> {
        let mut cache = TpmMediansCache { cancer_types: Default::default(), values: HashMap::new(), columns: Vec::new() };
        let Some(filename) = filename else { return Ok(cache) };
        if !std::path::Path::new(filename).exists() {
            error!("invalid cohort TPM medians file({})", filename);
            return Ok(cache);
        }
        let mut reader = hmf_common::gz::LineReader::open_threaded(filename)?;
        let Some(header) = reader.next_line()? else { bail!("NullPointerException: empty cohort TPM medians file({filename})") };
        let header = String::from_utf8_lossy(header).into_owned();
        let f = fields_index_map(&header, ',');
        for k in f.keys() {
            if k != "GeneId" && k != "GeneName" && k != "TransName" {
                cache.cancer_types.insert(k.clone());
            }
        }
        let trans_i = *f.get("TransName").ok_or_else(|| anyhow!("NullPointerException: no TransName column in {filename}"))?;
        let mut cols: Vec<usize> = Vec::new();
        let mut wanted: Vec<String> = vec![ALL_TYPES.to_string()];
        for c in wanted_cancer_types {
            if !c.is_empty() && !wanted.contains(c) {
                wanted.push(c.clone());
            }
        }
        for c in &wanted {
            if cache.cancer_types.contains(c) {
                cache.columns.push(c.clone());
                cols.push(f[c]);
            }
        }
        while let Some(line) = reader.next_line()? {
            let line = std::str::from_utf8(line).map_err(|e| anyhow!("{filename}: {e}"))?;
            let items = split_all(line, ',');
            let g = |i: usize| items.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: Index {i}"));
            let trans = g(trans_i)?;
            let mut vals = Vec::with_capacity(cols.len());
            for &c in &cols {
                let s = g(c)?;
                vals.push(parse_double(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))?);
            }
            cache.values.insert(trans.to_string(), vals);
        }
        info!("loaded {} cohort TPM median data", cache.values.len());
        Ok(cache)
    }

    pub fn has_cancer_type(&self, cancer_type: &str) -> bool {
        self.cancer_types.contains(cancer_type)
    }

    /// getTranscriptTpm: [cancer type, pan-cancer], NO_TPM_VALUE where missing
    pub fn transcript_tpm(&self, trans_name: &str, cancer_type: &str) -> [f64; 2] {
        let mut r = [NO_TPM_VALUE, NO_TPM_VALUE];
        if let Some(v) = self.values.get(trans_name) {
            for (i, c) in self.columns.iter().enumerate() {
                if c == ALL_TYPES {
                    r[1] = v[i];
                }
                if !cancer_type.is_empty() && c == cancer_type {
                    r[0] = v[i];
                }
            }
        }
        r
    }
}

/// ReferenceData
pub struct ReferenceData {
    pub scorer: BindScorer,
    /// transcript name -> amino acids
    pub trans_amino_acids: HashMap<String, String>,
    pub transcript_expression: Option<RnaExpressionMatrix>,
    pub tpm_medians: TpmMediansCache,
}

impl ReferenceData {
    pub fn new(cb: &ConfigBuilder, config: &NeoScorerConfig) -> Result<ReferenceData> {
        let scorer = BindScorer::new(ScoreConfig::new(cb));
        let transcript_expression = match cb.value("cohort_trans_exp_file") {
            Some(f) => {
                info!("loading cohort transcript expression");
                Some(RnaExpressionMatrix::load(&f, EXPRESSION_SCOPE_TRANS)?)
            }
            None => None,
        };
        info!("loading cohort transcript medians");
        let cancer_types: Vec<String> = config.samples.iter().map(|s| s.cancer_type.clone()).collect();
        let tpm_medians = TpmMediansCache::new(cb.value("cancer_tpm_medians_file").as_deref(), &cancer_types)?;
        // the result of loadTranscriptAminoAcidData is ignored: without the file no peptide is found in a
        // wild-type protein
        let mut trans_amino_acids = HashMap::new();
        let version = RefGenomeVersion::from_java(&cb.value_or("ref_genome_version", "V37"))?;
        let dir = cb.value("ensembl_data_dir").map(|d| check_add_dir_separator(&d)).unwrap_or_else(|| "null".to_string());
        let filename = check_old_filename(&format!("{dir}{}", ensembl_amino_acid_data_file(version)));
        if std::path::Path::new(&filename).exists() {
            for (name, t) in load_transcript_amino_acid_data(&dir, version, &[], false)? {
                trans_amino_acids.insert(name, t.amino_acids);
            }
        }
        Ok(ReferenceData { scorer, trans_amino_acids, transcript_expression, tpm_medians })
    }

    /// peptideMatchesWildtype
    pub fn peptide_matches_wildtype(&self, peptide: &str, trans_name: &str) -> bool {
        self.trans_amino_acids.get(trans_name).is_some_and(|aas| aas.contains(peptide))
    }
}

// ---------------------------------------------------------------------------------------------
// neo-epitope data

/// NeoRnaData: sample and cohort TPMs, RNA fragment support and depth
#[derive(Clone, Debug)]
pub struct NeoRnaData {
    pub tpm_cancer: [f64; 2],
    pub tpm_pan_cancer: [f64; 2],
    pub trans_expression: [f64; 2],
    has_expression: bool,
    fragment_support: i32,
    base_depth: [i32; 2],
    has_coverage: bool,
}

impl Default for NeoRnaData {
    fn default() -> Self {
        NeoRnaData {
            tpm_cancer: [NO_TPM_VALUE; 2],
            tpm_pan_cancer: [NO_TPM_VALUE; 2],
            trans_expression: [NO_TPM_VALUE; 2],
            has_expression: false,
            fragment_support: 0,
            base_depth: [0, 0],
            has_coverage: false,
        }
    }
}

impl NeoRnaData {
    pub fn has_coverage(&self) -> bool {
        self.has_coverage
    }
    pub fn has_expression(&self) -> bool {
        self.has_expression
    }
    pub fn fragment_support(&self) -> i32 {
        self.fragment_support
    }
    pub fn set_cohort_values(&mut self, tpm_cancer: [f64; 2], tpm_pan_cancer: [f64; 2]) {
        self.tpm_cancer = tpm_cancer;
        self.tpm_pan_cancer = tpm_pan_cancer;
    }
    pub fn set_expression(&mut self, expression: Option<[f64; 2]>) {
        if let Some(e) = expression {
            self.has_expression = true;
            self.trans_expression = e;
        }
    }
    pub fn set_coverage(&mut self, fragments: i32, depth_up: i32, depth_down: i32) {
        self.has_coverage = true;
        self.base_depth = [depth_up, depth_down];
        self.fragment_support = fragments;
    }
    /// getTPM: the sample's TPM, else the cancer type's, else the pan-cancer one
    pub fn tpm(&self, stream: usize) -> f64 {
        if self.trans_expression[stream] != NO_TPM_VALUE {
            return self.trans_expression[stream];
        }
        if self.tpm_cancer[stream] != NO_TPM_VALUE { self.tpm_cancer[stream] } else { self.tpm_pan_cancer[stream] }
    }
    pub fn average_base_depth(&self) -> f64 {
        (self.base_depth[0] + self.base_depth[1]) as f64 * 0.5
    }
}

/// NeoEpitopeData: a neo-epitope as the scorer reads it, with its RNA data, peptides and TPMs
#[derive(Clone, Debug)]
pub struct NeoEpitopeData {
    pub id: i32,
    pub variant_type: NeoEpitopeType,
    pub variant_info: String,
    pub gene_name: String,
    pub up_aas: String,
    pub novel_aas: String,
    pub down_aas: String,
    pub transcripts: [Vec<String>; 2],
    pub nmd_bases: [i32; 2],
    pub coding_bases_length: [i32; 2],
    pub fused_intron_length: i32,
    pub skipped_acceptors_donors: [i32; 2],
    pub variant_copy_number: f64,
    pub copy_number: f64,
    pub subclonal_likelihood: f64,
    pub rna: NeoRnaData,
    /// indices into the sample's peptide list
    pub peptides: Vec<usize>,
    pub raw_effective_tpm: f64,
    pub effective_tpm: f64,
    pub expected_tpm: f64,
    pub copy_number_factor: f64,
}

impl NeoEpitopeData {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: i32, variant_type: NeoEpitopeType, variant_info: &str, gene_name: &str, up_aas: &str, novel_aas: &str, down_aas: &str,
        trans_up: Vec<String>, trans_down: Vec<String>, nmd: [i32; 2], variant_copy_number: f64, copy_number: f64,
        subclonal_likelihood: f64, coding_bases_length: [i32; 2], fused_intron_length: i32, skipped: [i32; 2],
    ) -> NeoEpitopeData {
        NeoEpitopeData {
            id,
            variant_type,
            variant_info: variant_info.to_string(),
            gene_name: gene_name.to_string(),
            up_aas: up_aas.to_string(),
            novel_aas: novel_aas.to_string(),
            down_aas: down_aas.to_string(),
            transcripts: [trans_up, trans_down],
            nmd_bases: nmd,
            coding_bases_length,
            fused_intron_length,
            skipped_acceptors_donors: skipped,
            variant_copy_number,
            copy_number,
            subclonal_likelihood,
            rna: NeoRnaData::default(),
            peptides: Vec::new(),
            raw_effective_tpm: NO_TPM_VALUE,
            effective_tpm: 0.0,
            expected_tpm: NO_TPM_VALUE,
            copy_number_factor: 0.0,
        }
    }

    /// setExpressionData: the sample's TPMs (per sample, or from the cohort matrix) and the cohort medians,
    /// summed over each stream's transcripts. Java starts the cohort sums at the -1 "no value" sentinel
    /// (NE6), so each is 1 too low; corrected, a sum starts at 0 when its first value is added.
    pub fn set_expression_data(
        &mut self, sample: &SampleData, sample_tpms: &HashMap<String, f64>, matrix: Option<&RnaExpressionMatrix>, medians: &TpmMediansCache,
    ) {
        let mut sample_expression: Option<[f64; 2]> = None;
        let mut pan = [NO_TPM_VALUE; 2];
        let mut cancer = [NO_TPM_VALUE; 2];
        let mut pan_set = [false; 2];
        let mut cancer_set = [false; 2];
        let java_ne6 = compat::java("NE6");
        for fs in [FS_UP, FS_DOWN] {
            for trans_name in &self.transcripts[fs] {
                let mut trans_expression = INVALID_EXP;
                if !sample_tpms.is_empty() {
                    match sample_tpms.get(trans_name) {
                        Some(&v) => trans_expression = v,
                        None => warn!("sample({}) missing transcript({}) TPM", sample.tumor_id, trans_name),
                    }
                } else if let Some(m) = matrix {
                    trans_expression = m.expression(trans_name, &sample.tumor_id);
                }
                if trans_expression != INVALID_EXP {
                    let e = sample_expression.get_or_insert([0.0, 0.0]);
                    e[fs] += trans_expression;
                }
                let tpm = medians.transcript_tpm(trans_name, &sample.cancer_type);
                if tpm[1] != NO_TPM_VALUE {
                    if !java_ne6 && !pan_set[fs] {
                        pan[fs] = 0.0;
                        pan_set[fs] = true;
                    }
                    pan[fs] += tpm[1];
                }
                if tpm[0] != NO_TPM_VALUE {
                    if !java_ne6 && !cancer_set[fs] {
                        cancer[fs] = 0.0;
                        cancer_set[fs] = true;
                    }
                    cancer[fs] += tpm[0];
                }
            }
        }
        self.rna.set_expression(sample_expression);
        self.rna.set_cohort_values(cancer, pan);
    }

    /// setFusionRnaSupport
    pub fn set_fusion_rna_support(&mut self, rna_neos: &[RnaNeoEpitope]) {
        if !self.variant_type.is_fusion() {
            return;
        }
        if let Some(m) = rna_neos.iter().find(|x| x.id == self.id && x.variant_info == self.variant_info) {
            self.rna.set_coverage(m.fragment_count, m.base_depth[0], m.base_depth[1]);
        }
    }

    /// setMutationRnaSupport
    pub fn set_mutation_rna_support(&mut self, variants: &[SomaticVariant]) {
        if !self.variant_type.is_point_mutation() {
            return;
        }
        if let Some(m) = variants.iter().find(|x| x.variant_info() == self.variant_info) {
            self.rna.set_coverage(m.rna_fragments, m.rna_depth, m.rna_depth);
        }
    }

    /// matchingTranscripts
    pub fn matching_transcripts(&self, other: &NeoEpitopeData, stream: usize) -> bool {
        if self.transcripts[stream].len() != other.transcripts[stream].len() {
            return false;
        }
        self.transcripts[stream].iter().all(|x| other.transcripts[stream].contains(x))
    }

    /// setCalculatedTpmValues
    pub fn set_calculated_tpm_values(&mut self, norm_factor: f64, sample_ploidy: f64, prob_low: f64, prob_high: f64) {
        self.raw_effective_tpm = if norm_factor != NO_TPM_VALUE { self.rna.fragment_support as f64 * norm_factor } else { NO_TPM_VALUE };
        let sc_likelihood = if self.variant_type.is_point_mutation() { self.subclonal_likelihood } else { 1.0 };
        let variant_cn = self.variant_copy_number.max(0.0);
        let segment_cn = if self.rna.has_expression() { self.copy_number } else { sample_ploidy };
        let adjusted_cn = if segment_cn > 0.0 { variant_cn / segment_cn } else { 1.0 };
        self.copy_number_factor = 1.0 - sc_likelihood * (1.0 - f64::min(1.0, adjusted_cn));
        self.expected_tpm = self.rna.tpm(FS_UP) * self.copy_number_factor;
        if norm_factor == NO_TPM_VALUE || !self.rna.has_expression() || self.raw_effective_tpm == NO_TPM_VALUE {
            self.effective_tpm = self.expected_tpm;
        } else {
            self.effective_tpm = calc_effective_tpm(norm_factor, self.raw_effective_tpm, self.expected_tpm, prob_low, prob_high);
        }
    }
}

/// PeptideScoreData: a peptide of one or more neo-epitopes of a variant, its effective TPM and allele scores
#[derive(Clone, Debug)]
pub struct PeptideScoreData {
    pub peptide: String,
    pub up_flank: String,
    pub down_flank: String,
    pub ne_ids: JavaHashSet<i32>,
    pub effective_tpm: f64,
    pub tpm_calculated: bool,
    pub alleles: Vec<BindData>,
    pub written: bool,
}

impl PeptideScoreData {
    pub fn new(peptide: &str, up_flank: &str, down_flank: &str, ne_id: i32) -> PeptideScoreData {
        PeptideScoreData {
            peptide: peptide.to_string(),
            up_flank: up_flank.to_string(),
            down_flank: down_flank.to_string(),
            // Sets.newHashSet(neId)
            ne_ids: JavaHashSet::of_one(ne_id),
            effective_tpm: 0.0,
            tpm_calculated: false,
            alleles: Vec::new(),
            written: false,
        }
    }
    pub fn add_allele(&mut self, allele: &str) {
        let mut b = BindData::new(allele, &self.peptide, &self.up_flank, &self.down_flank);
        b.set_tpm(self.effective_tpm);
        self.alleles.push(b);
    }
    pub fn set_effective_tpm(&mut self, tpm: f64) {
        self.effective_tpm = tpm;
        self.tpm_calculated = true;
    }
    /// neoIdsStr
    pub fn ne_ids_str(&self) -> String {
        self.ne_ids.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(";")
    }
}

/// AlleleCoverage
#[derive(Clone, Debug)]
pub struct AlleleCoverage {
    pub allele: String,
    pub copy_number: f64,
    pub variant_count: f64,
}

impl AlleleCoverage {
    pub fn is_lost(&self) -> bool {
        self.copy_number < 0.5 || self.variant_count > 0.0
    }
}

/// SomaticVariant: a point mutation's tumor and RNA support from an RNA-annotated VCF (the tumor counts are
/// read, as in Java, but not used)
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct SomaticVariant {
    pub chromosome: String,
    pub position: i32,
    pub ref_bases: String,
    pub alt: String,
    pub tumor_fragments: i32,
    pub tumor_depth: i32,
    pub rna_fragments: i32,
    pub rna_depth: i32,
}

impl SomaticVariant {
    pub fn variant_info(&self) -> String {
        point_mutation_info(&self.chromosome, self.position, &self.ref_bases, &self.alt)
    }
}

// ---------------------------------------------------------------------------------------------
// DataLoader

fn int_of(s: &str) -> Result<i32> {
    parse_int(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
}

fn dbl_of(s: &str) -> Result<f64> {
    parse_double(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
}

/// DataLoader.loadNeoEpitopes: Ok(None) where Java logs an IOException and returns null
pub fn load_neo_epitopes(sample_id: &str, neo_dir: &str) -> Result<Option<Vec<NeoEpitopeData>>> {
    let filename = neo_epitope_filename(neo_dir, sample_id);
    let lines = match read_lines(&filename) {
        Ok(l) => l,
        Err(e) => {
            error!("failed to read sample({}) neo-epitope file: java.nio.file.NoSuchFileException: {} ({})", sample_id, filename, e);
            return Ok(None);
        }
    };
    let header = lines.first().ok_or_else(|| anyhow!("IndexOutOfBoundsException: Index 0 out of bounds for length 0 ({filename})"))?;
    let f = fields_index_map(header, '\t');
    let c = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("NullPointerException: no column {n} in {filename}"));
    let idx = [
        c("NeId")?,
        c("VariantType")?,
        c("VariantInfo")?,
        c("GeneIdDown")?,
        c("GeneNameUp")?,
        c("GeneNameDown")?,
        c("UpstreamAA")?,
        c("DownstreamAA")?,
        c("NovelAA")?,
        c("NmdMin")?,
        c("NmdMax")?,
        c("VariantCopyNumber")?,
        c("CopyNumber")?,
        c("SubclonalLikelihood")?,
        c("CodingBasesLengthMin")?,
        c("CodingBasesLengthMax")?,
        c("FusedIntronLength")?,
        c("SkippedDonors")?,
        c("SkippedAcceptors")?,
        c("UpTranscripts")?,
        c("DownTranscripts")?,
    ];
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let v = split_all(line, '\t');
        let g = |k: usize| v.get(idx[k]).copied().ok_or_else(|| anyhow!("ArrayIndexOutOfBoundsException: Index {} ({filename})", idx[k]));
        let id = int_of(g(0)?)?;
        let (gene_up, gene_down) = (g(4)?, g(5)?);
        let gene_name = if gene_up == gene_down { gene_up.to_string() } else { format!("{gene_up}_{gene_down}") };
        let (trans_up, trans_down) = extract_transcript_names(g(19)?, g(20)?);
        out.push(NeoEpitopeData::new(
            id,
            parse_neo_epitope_type(g(1)?)?,
            g(2)?,
            &gene_name,
            g(6)?,
            g(8)?,
            g(7)?,
            trans_up,
            trans_down,
            [int_of(g(9)?)?, int_of(g(10)?)?],
            dbl_of(g(11)?)?,
            dbl_of(g(12)?)?,
            dbl_of(g(13)?)?,
            [int_of(g(14)?)?, int_of(g(15)?)?],
            int_of(g(16)?)?,
            [int_of(g(17)?)?, int_of(g(18)?)?],
        ));
    }
    Ok(Some(out))
}

/// DataLoader.loadAlleleCoverage: Ok(None) where Java logs an IOException and returns null
pub fn load_allele_coverage(sample_id: &str, lilac_dir: &str) -> Result<Option<Vec<AlleleCoverage>>> {
    let filename = lilac_allele_filename(lilac_dir, sample_id);
    let checked = hmf_common::hla::check_file_extension_rename(&filename);
    if !std::path::Path::new(&checked).exists() {
        error!("failed to read sample({}) Lilac allele coverage file: java.nio.file.NoSuchFileException: {}", sample_id, checked);
        return Ok(None);
    }
    let alleles = read_lilac_alleles(&filename)?;
    Ok(Some(
        alleles
            .iter()
            .map(|a| AlleleCoverage {
                allele: a.allele.replace('*', "").replace(':', ""),
                copy_number: a.tumor_copy_number,
                variant_count: a.somatic_missense + a.somatic_splice + a.somatic_nonsense_or_frameshift + a.somatic_inframe_indel,
            })
            .collect(),
    ))
}

/// DataLoader.loadRnaNeoData: Ok(None) where Java returns null (the file missing with RNA, or unreadable)
pub fn load_rna_neo_data(sample: &SampleData, isofox_dir: &str) -> Result<Option<Vec<RnaNeoEpitope>>> {
    let filename = rna_neo_epitope_filename(isofox_dir, sample.rna_id_str());
    if !std::path::Path::new(&filename).exists() {
        return Ok(if sample.has_rna { None } else { Some(Vec::new()) });
    }
    match read_rna_neo_epitopes(&filename) {
        Ok(v) => {
            let v: Vec<RnaNeoEpitope> = v.into_iter().filter(|x| x.variant_type.is_fusion()).collect();
            debug!("sample({}) loaded {} fusion RNA-annotated neoepitopes", sample.rna_id_str(), v.len());
            Ok(Some(v))
        }
        Err(e) if e.to_string().contains("NoSuchFile") => {
            error!("failed to read sample({}) Isofox neoepitopes file: {}", sample.rna_id_str(), e);
            Ok(None)
        }
        Err(e) => Err(e),
    }
}

/// SomaticVariant.extractVariantCounts: the read context counts FULL to REALIGNED (the first four), and DP
fn extract_variant_counts(g: &hmf_common::variant_context::Genotype) -> Result<(i32, i32)> {
    let raw = match g.attrs.get("RC_CNT") {
        Some(Value::Raw(s)) => s.clone(),
        Some(v) => v.java_string(),
        None => "0".to_string(),
    };
    let counts = split_all(&raw, ',');
    let mut support = 0;
    for i in 0..=3 {
        let s = counts.get(i).ok_or_else(|| anyhow!("ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", counts.len()))?;
        support += int_of(s)?;
    }
    Ok((g.dp.unwrap_or(-1), support))
}

/// DataLoader.loadSomaticVariants: the unfiltered variants of the RNA-annotated VCF that match a point
/// mutation neo-epitope; Ok(None) where Java returns null
pub fn load_somatic_variants(sample: &SampleData, rna_somatic_vcf: &str, point_neos: &[&NeoEpitopeData]) -> Result<Option<Vec<SomaticVariant>>> {
    let vcf = rna_somatic_vcf.replace('*', &sample.tumor_id);
    if !std::path::Path::new(&vcf).exists() {
        return Ok(if sample.has_rna { None } else { Some(Vec::new()) });
    }
    let result = (|| -> Result<Vec<SomaticVariant>> {
        let mut reader = VcfReader::open(&vcf)?;
        let header = reader.header.clone();
        let mut matched = Vec::new();
        while let Some(mut rec) = reader.next_record()? {
            if rec.is_filtered() {
                continue;
            }
            rec.decode(&header)?;
            let Some(tumor) = rec.decoded_genotype(&sample.tumor_id) else { continue };
            let Some(rna_id) = sample.rna_sample_id.as_deref() else { continue };
            let Some(rna) = rec.decoded_genotype(rna_id) else { continue };
            let (tumor_depth, tumor_fragments) = extract_variant_counts(tumor)?;
            let (rna_depth, rna_fragments) = extract_variant_counts(rna)?;
            let alt = if rec.alts.is_empty() { rec.ref_allele.clone() } else { rec.alts[0].clone() };
            let v = SomaticVariant {
                chromosome: rec.chrom.clone(),
                position: rec.pos,
                ref_bases: rec.ref_allele.clone(),
                alt,
                tumor_fragments,
                tumor_depth,
                rna_fragments,
                rna_depth,
            };
            let info = v.variant_info();
            if point_neos.iter().any(|x| x.variant_info == info) {
                matched.push(v);
            }
        }
        Ok(matched)
    })();
    match result {
        Ok(v) => {
            debug!("sample({}) loaded {} somatic RNA-annotated variants", sample.tumor_id, v.len());
            Ok(Some(v))
        }
        Err(e) => {
            error!("failed to read somatic VCF file({}): {:#}", vcf, e);
            Ok(None)
        }
    }
}

/// DataLoader.loadTranscriptExpression: TPM by transcript name (Collectors.toMap: a repeated name throws)
pub fn load_transcript_expression(isofox_dir: &str, rna_sample_id: &str) -> Result<HashMap<String, f64>> {
    let filename = transcript_expression_filename(isofox_dir, rna_sample_id);
    let mut map = HashMap::new();
    for t in read_transcript_expression(&filename)? {
        if map.insert(t.transcript_name.clone(), t.tpm).is_some() {
            bail!("IllegalStateException: Duplicate key {}", t.transcript_name);
        }
    }
    Ok(map)
}

// ---------------------------------------------------------------------------------------------
// TpmCalculator

/// TpmCalculator constants
pub const LOW_PROBABILITY: f64 = 0.25;
pub const HIGH_PROBABILITY: f64 = 1.0 - LOW_PROBABILITY;
pub const LOW_PROBABILITY_ZERO_OBS_MEAN: f64 = 1.39;
pub const EFFECTIVE_TPM_ACTUAL_PERC: f64 = 0.8;

/// the bundled `/ref/expected_fragment_ranges.csv`
const EXPECTED_FRAGMENT_RANGES: &str = include_str!("../resources/expected_fragment_ranges.csv");

/// TpmCalculator.calcEffectiveTpm: 0.8 x the raw effective TPM + 0.2 x the expected TPM clamped to the
/// Poisson range of the fragment count
pub fn calc_effective_tpm(norm_factor: f64, raw_effective_tpm: f64, expected_tpm: f64, prob_low: f64, prob_high: f64) -> f64 {
    let range_low = prob_low * norm_factor;
    let range_high = prob_high * norm_factor;
    EFFECTIVE_TPM_ACTUAL_PERC * raw_effective_tpm + (1.0 - EFFECTIVE_TPM_ACTUAL_PERC) * f64::max(range_high, f64::min(range_low, expected_tpm))
}

/// TpmCalculator.calculateTpmNormalisation: the summed up-stream TPM over the summed RNA depth of the missense
/// neo-epitopes with both
pub fn calculate_tpm_normalisation(neos: &[NeoEpitopeData]) -> f64 {
    let mut tpm_up_total = 0.0;
    let mut depth_total = 0.0;
    let mut count = 0;
    for n in neos {
        if n.variant_type != NeoEpitopeType::Missense {
            continue;
        }
        if n.rna.has_coverage() && n.rna.has_expression() {
            tpm_up_total += n.rna.trans_expression[FS_UP];
            depth_total += n.rna.average_base_depth();
            count += 1;
        }
    }
    if count > 0 { tpm_up_total / depth_total } else { NO_TPM_VALUE }
}

/// TpmCalculator.calcPoissonMean
pub fn calc_poisson_mean(fragment_count: i32, required_prob: f64) -> f64 {
    if fragment_count < 0 {
        return 0.0;
    }
    if fragment_count == 0 {
        return if required_prob == HIGH_PROBABILITY { 0.0 } else { LOW_PROBABILITY_ZERO_OBS_MEAN };
    }
    hmf_common::poisson::calc_poisson_noise_value(fragment_count, required_prob)
}

struct SharedUpstreamGroup {
    ne_ids: JavaHashSet<i32>,
    tpm_up: f64,
    tpm_down_total: f64,
}

/// TpmCalculator
pub struct TpmCalculator {
    flank_count: i32,
    peptide_length_range: [i32; 2],
    poisson_range_values: HashMap<i32, (f64, f64)>,
}

impl TpmCalculator {
    pub fn new(flank_count: i32, peptide_length_range: [i32; 2]) -> TpmCalculator {
        let mut poisson_range_values = HashMap::new();
        for line in EXPECTED_FRAGMENT_RANGES.lines().skip(1) {
            let v: Vec<&str> = line.splitn(3, ',').collect();
            if v.len() == 3 {
                if let (Some(f), Some(l), Some(h)) = (parse_int(v[0]), parse_double(v[1]), parse_double(v[2])) {
                    poisson_range_values.insert(f, (l, h));
                }
            }
        }
        TpmCalculator { flank_count, peptide_length_range, poisson_range_values }
    }

    fn range_values(&mut self, fragment_count: i32) -> (f64, f64) {
        if let Some(&v) = self.poisson_range_values.get(&fragment_count) {
            return v;
        }
        let v = (calc_poisson_mean(fragment_count, LOW_PROBABILITY), calc_poisson_mean(fragment_count, HIGH_PROBABILITY));
        self.poisson_range_values.insert(fragment_count, v);
        v
    }

    /// compute: the peptides of each neo-epitope, shared between the neo-epitopes of a variant, and their
    /// expected and effective TPMs; the peptides go into `peptides`
    pub fn compute(&mut self, sample_id: &str, neos: &mut [NeoEpitopeData], sample_ploidy: f64, peptides: &mut Vec<PeptideScoreData>) {
        let norm_factor = calculate_tpm_normalisation(neos);
        debug!("sample({}) TPM norm-factor({})", sample_id, java_fixed(norm_factor, 3));
        // the neo-epitopes by variant; Java's HashMap order does not matter, as each variant's peptides are
        // its own
        let mut groups: Vec<(String, Vec<usize>)> = Vec::new();
        let mut index: HashMap<String, usize> = HashMap::new();
        for (i, n) in neos.iter().enumerate() {
            match index.get(&n.variant_info) {
                Some(&g) => groups[g].1.push(i),
                None => {
                    index.insert(n.variant_info.clone(), groups.len());
                    groups.push((n.variant_info.clone(), vec![i]));
                }
            }
        }
        let mut total = 0;
        for (_, variant_neos) in &groups {
            total += self.process_variant_neos(neos, variant_neos, norm_factor, sample_ploidy, peptides);
        }
        debug!("sample({}) neoepitopes({}) derived {} peptides", sample_id, neos.len(), total);
    }

    fn process_variant_neos(
        &mut self, neos: &mut [NeoEpitopeData], vn: &[usize], norm_factor: f64, sample_ploidy: f64, all_peptides: &mut Vec<PeptideScoreData>,
    ) -> usize {
        // neo-epitopes with the same upstream transcripts
        let mut shared: Option<Vec<SharedUpstreamGroup>> = None;
        if vn.len() > 1 {
            for i in 0..vn.len() - 1 {
                let n1 = &neos[vn[i]];
                // NE7: Java tests whether a group's ids (integers) contain the neo-epitope object, never true
                if !compat::java("NE7") && shared.as_ref().is_some_and(|g| g.iter().any(|x| x.ne_ids.contains(&n1.id))) {
                    continue;
                }
                let mut group: Option<usize> = None;
                for &j in &vn[i + 1..] {
                    let n2 = &neos[j];
                    if n1.matching_transcripts(n2, FS_UP) {
                        let groups = shared.get_or_insert_with(Vec::new);
                        let gi = match group {
                            Some(g) => g,
                            None => {
                                groups.push(SharedUpstreamGroup {
                                    ne_ids: JavaHashSet::of_one(n1.id),
                                    tpm_up: n1.rna.tpm(FS_UP),
                                    tpm_down_total: n1.rna.tpm(FS_DOWN),
                                });
                                group = Some(groups.len() - 1);
                                groups.len() - 1
                            }
                        };
                        groups[gi].ne_ids.insert(n2.id);
                        groups[gi].tpm_down_total += n2.rna.tpm(FS_DOWN);
                    }
                }
            }
        }

        // the peptides of each neo-epitope, those repeated across the variant's neo-epitopes shared
        let start = all_peptides.len();
        let mut by_peptide: HashMap<String, usize> = HashMap::new();
        for &ni in vn {
            let n = &neos[ni];
            for p in generate_peptides(&n.up_aas, &n.novel_aas, &n.down_aas, self.peptide_length_range, self.flank_count) {
                if p.peptide.contains('U') {
                    continue;
                }
                match by_peptide.get(&p.peptide) {
                    Some(&k) => {
                        all_peptides[k].ne_ids.insert(n.id);
                    }
                    None => {
                        by_peptide.insert(p.peptide.clone(), all_peptides.len());
                        all_peptides.push(PeptideScoreData::new(&p.peptide, &p.up_flank, &p.down_flank, n.id));
                    }
                }
            }
        }
        let peptide_count = all_peptides.len() - start;

        // expected and effective TPMs
        for &ni in vn {
            let id = neos[ni].id;
            let neo_peptides: Vec<usize> = (start..all_peptides.len()).filter(|&k| all_peptides[k].ne_ids.contains(&id)).collect();
            neos[ni].peptides.extend(neo_peptides.iter().copied());
            let (low, high) = self.range_values(neos[ni].rna.fragment_support());
            neos[ni].set_calculated_tpm_values(norm_factor, sample_ploidy, low, high);
            let raw_effective_tpm = neos[ni].raw_effective_tpm;
            for &k in &neo_peptides {
                if all_peptides[k].tpm_calculated {
                    continue;
                }
                let mut tpm_up = 0.0;
                if all_peptides[k].ne_ids.len() == 1 && shared.is_none() {
                    tpm_up = neos[ni].rna.tpm(FS_UP);
                } else {
                    let mut processed: Vec<usize> = Vec::new();
                    for ne_id in all_peptides[k].ne_ids.iter() {
                        let group = shared.as_ref().and_then(|g| g.iter().position(|x| x.ne_ids.contains(ne_id)));
                        match group {
                            Some(gi) => {
                                if !processed.contains(&gi) {
                                    processed.push(gi);
                                    let g = &shared.as_ref().unwrap()[gi];
                                    let mut tpm_down = 0.0;
                                    for group_ne_id in g.ne_ids.iter() {
                                        if !all_peptides[k].ne_ids.contains(group_ne_id) {
                                            continue;
                                        }
                                        let pn = vn.iter().find(|&&x| neos[x].id == *group_ne_id).expect("neo-epitope of the group");
                                        tpm_down += neos[*pn].rna.tpm(FS_DOWN);
                                    }
                                    if g.tpm_down_total > 0.0 {
                                        tpm_up += g.tpm_up * tpm_down / g.tpm_down_total;
                                    }
                                }
                            }
                            None => {
                                let pn = vn.iter().find(|&&x| neos[x].id == *ne_id).expect("neo-epitope of the peptide");
                                tpm_up += neos[*pn].rna.tpm(FS_UP);
                            }
                        }
                    }
                }
                let n = &neos[ni];
                let expected_tpm = tpm_up * n.copy_number_factor;
                let effective_tpm = if norm_factor == NO_TPM_VALUE || !n.rna.has_expression() {
                    expected_tpm
                } else {
                    calc_effective_tpm(norm_factor, raw_effective_tpm, expected_tpm, low, high)
                };
                all_peptides[k].set_effective_tpm(effective_tpm);
            }
        }
        peptide_count
    }
}

// ---------------------------------------------------------------------------------------------
// NeoDataWriter

/// NeoDataWriter.passesThreshold
pub fn passes_threshold(value: f64, threshold: f64) -> bool {
    value != INVALID_CALC && (threshold == 0.0 || threshold > 0.0 && value <= threshold)
}

struct OutFile {
    path: String,
    file: std::fs::File,
    buf: Vec<u8>,
}

/// NeoDataWriter: Java opens both files at the start and writes them through BufferedWriters closed at the
/// end; where Java exits before closing them (no samples, scoring data that fails to load), they stay empty.
/// The port writes each file's text at close.
pub struct NeoDataWriter {
    multi: bool,
    peptide: Option<OutFile>,
    neo: Option<OutFile>,
    likelihood_threshold: f64,
    similarity_threshold: f64,
}

fn open_out(path: String, header: String, what: &str) -> Option<OutFile> {
    match std::fs::File::create(&path) {
        Ok(file) => {
            let mut buf = Vec::with_capacity(1 << 20);
            buf.extend_from_slice(header.as_bytes());
            buf.push(b'\n');
            Some(OutFile { path, file, buf })
        }
        Err(e) => {
            error!("failed to create {} writer: java.io.FileNotFoundException: {} ({})", what, path, e);
            None
        }
    }
}

impl NeoDataWriter {
    pub fn new(config: &NeoScorerConfig) -> NeoDataWriter {
        let multi = config.samples.len() > 1;
        let prefix = if multi { "SampleId\t" } else { "" };
        let peptide = if config.write_types.contains(&OutputType::AllelePeptide) {
            let h = "NeIds\tVariantType\tVariantInfo\tGene\tAllele\tPeptide\tFlankUp\tFlankDown\tScore\tRank\tLikelihood\tLikelihoodRank\t\
                     ExpLikelihood\tExpLikelihoodRank\tEffectiveTpm\tRecogSim\tOtherAlleleRecogSim\tAlleleCN\tAlleleDisrupted";
            open_out(config.form_filename("peptide_scores"), format!("{prefix}{h}"), "peptide")
        } else {
            None
        };
        let neo = if config.write_types.contains(&OutputType::Neoepitope) {
            let h = "NeId\tVariantType\tVariantInfo\tGeneName\tUpAminoAcids\tNovelAminoAcids\tDownAminoAcids\tPeptideCount\tTpmSource\t\
                     RnaFrags\tRnaDepth\tTpmUp\tTpmDown\tExpectedTpm\tRawEffectiveTpm\tEffectiveTpm\tTpmCancerUp\tTpmCancerDown\t\
                     TpmPanCancerUp\tTpmPanCancerDown\tNmdMin\tNmdMax\tCodingBasesLengthMin\tCodingBasesLengthMax\tFusedIntronLength\t\
                     SkippedDonors\tSkippedAcceptors\tTranscriptsUp\tTranscriptsDown\tVariantCopyNumber\tCopyNumber\tSubclonalLikelihood";
            open_out(config.form_filename("neoepitope"), format!("{prefix}{h}"), "neoepitope")
        } else {
            None
        };
        NeoDataWriter {
            multi,
            peptide,
            neo,
            likelihood_threshold: config.likelihood_threshold,
            similarity_threshold: config.similarity_threshold,
        }
    }

    pub fn close(mut self) -> Result<()> {
        for f in [self.neo.take(), self.peptide.take()].into_iter().flatten() {
            let mut file = f.file;
            file.write_all(&f.buf).with_context(|| format!("writing {}", f.path))?;
        }
        Ok(())
    }

    /// logPeptide
    fn log_peptide(&self, b: &BindData) -> bool {
        passes_threshold(b.likelihood_rank, self.likelihood_threshold)
            || passes_threshold(b.recognition_similarity, self.similarity_threshold)
            || passes_threshold(b.other_allele_recognition_similarity, self.similarity_threshold)
    }

    /// writePeptideData: each peptide once, with the first of its neo-epitopes to be written
    pub fn write_peptide_data(&mut self, sample_id: &str, neo: &NeoEpitopeData, peptides: &mut [PeptideScoreData], coverages: &[AlleleCoverage]) {
        if self.peptide.is_none() {
            return;
        }
        let mut out = String::new();
        for &k in &neo.peptides {
            let p = &mut peptides[k];
            if p.written {
                continue;
            }
            p.written = true;
            for b in &p.alleles {
                if !self.log_peptide(b) {
                    continue;
                }
                let cov = coverages.iter().find(|x| x.allele == b.allele).expect("allele coverage of a scored allele");
                if self.multi {
                    out.push_str(sample_id);
                    out.push('\t');
                }
                let fields = [
                    p.ne_ids_str(),
                    neo.variant_type.name().to_string(),
                    neo.variant_info.clone(),
                    neo.gene_name.clone(),
                    b.allele.clone(),
                    b.peptide.clone(),
                    b.up_flank.clone(),
                    b.down_flank.clone(),
                    java_fixed(b.score, 4),
                    java_fixed(b.rank_percentile, 6),
                    java_fixed(b.likelihood, 6),
                    java_fixed(b.likelihood_rank, 6),
                    java_fixed(b.exp_likelihood, 6),
                    java_fixed(b.exp_likelihood_rank, 6),
                    java_sci(p.effective_tpm, 3),
                    java_fixed(b.recognition_similarity, 1),
                    java_fixed(b.other_allele_recognition_similarity, 1),
                    java_fixed(cov.copy_number, 1),
                    cov.is_lost().to_string(),
                ];
                out.push_str(&fields.join("\t"));
                out.push('\n');
            }
        }
        self.peptide.as_mut().unwrap().buf.extend_from_slice(out.as_bytes());
    }

    /// writeNeoData
    pub fn write_neo_data(&mut self, sample_id: &str, tpm_source: TpmSource, neo: &NeoEpitopeData) {
        let Some(w) = self.neo.as_mut() else { return };
        let r = &neo.rna;
        let mut fields: Vec<String> = Vec::with_capacity(33);
        if self.multi {
            fields.push(sample_id.to_string());
        }
        fields.extend([
            neo.id.to_string(),
            neo.variant_type.name().to_string(),
            neo.variant_info.clone(),
            neo.gene_name.clone(),
            neo.up_aas.clone(),
            neo.novel_aas.clone(),
            neo.down_aas.clone(),
            neo.peptides.len().to_string(),
            tpm_source.name().to_string(),
            r.fragment_support().to_string(),
            java_fixed(r.average_base_depth(), 0),
            java_sci(r.trans_expression[FS_UP], 3),
            java_sci(r.trans_expression[FS_DOWN], 3),
            java_sci(neo.expected_tpm, 3),
            java_sci(neo.raw_effective_tpm, 3),
            java_sci(neo.effective_tpm, 3),
            java_sci(r.tpm_cancer[FS_UP], 3),
            java_sci(r.tpm_cancer[FS_DOWN], 3),
            java_sci(r.tpm_pan_cancer[FS_UP], 3),
            java_sci(r.tpm_pan_cancer[FS_DOWN], 3),
            neo.nmd_bases[0].to_string(),
            neo.nmd_bases[1].to_string(),
            neo.coding_bases_length[0].to_string(),
            neo.coding_bases_length[1].to_string(),
            neo.fused_intron_length.to_string(),
            neo.skipped_acceptors_donors[0].to_string(),
            neo.skipped_acceptors_donors[1].to_string(),
            neo.transcripts[FS_UP].join(";"),
            neo.transcripts[FS_DOWN].join(";"),
            java_fixed(neo.variant_copy_number, 4),
            java_fixed(neo.copy_number, 4),
            java_fixed(neo.subclonal_likelihood, 4),
        ]);
        w.buf.extend_from_slice(fields.join("\t").as_bytes());
        w.buf.push(b'\n');
    }
}

// ---------------------------------------------------------------------------------------------
// NeoScorer and NeoScorerTask

/// NeoScorer.run
pub fn run(cb: &ConfigBuilder) -> Result<()> {
    let t0 = std::time::Instant::now();
    let config = NeoScorerConfig::new(cb)?;
    let mut reference = ReferenceData::new(cb, &config)?;
    let mut writers = NeoDataWriter::new(&config);
    debug!("reference data loaded in {:.3}s", t0.elapsed().as_secs_f64());
    if config.samples.is_empty() {
        return Ok(());
    }
    if !reference.scorer.load_scoring_data()? {
        error!("failed to load scoring data");
        std::process::exit(1);
    }
    debug!("scoring data loaded in {:.3}s", t0.elapsed().as_secs_f64());
    if config.samples.len() == 1 {
        info!("running neoepitope scoring for {}", config.samples[0].tumor_id);
    } else {
        info!("running neoepitope scoring for {} samples", config.samples.len());
    }
    // Java splits several samples between -threads tasks, whose rows then interleave in the files as the
    // threads run (NE10); the port takes the samples in order (as one Java thread does) and uses the
    // threads to score each sample's peptides
    let reference = &reference;
    if config.samples.len() > 1 {
        info!("0: processing {} samples", config.samples.len());
    }
    for (i, sample) in config.samples.iter().enumerate() {
        process_sample(sample, &config, reference, &mut writers).with_context(|| format!("sample({}) failed processing", sample.tumor_id))?;
        if i > 0 && i % 10 == 0 {
            info!("0: processed {} samples of {}", i, config.samples.len());
        }
    }
    if config.samples.len() > 1 {
        info!("0: processing complete");
    }
    writers.close()?;
    info!("Neo peptide scoring complete");
    Ok(())
}

/// NeoScorerTask.processSample
pub fn process_sample(sample: &SampleData, config: &NeoScorerConfig, reference: &ReferenceData, writers: &mut NeoDataWriter) -> Result<()> {
    let sample_id = sample.tumor_id.as_str();
    let neos = load_neo_epitopes(sample_id, &config.neo_dir)?;
    let coverages = load_allele_coverage(sample_id, &config.lilac_dir)?;
    let (Some(mut neos), Some(coverages)) = (neos, coverages) else { std::process::exit(1) };
    info!("sample({}) processing {} neoepitopes", sample_id, neos.len());

    // Collectors.toSet(): a HashSet in Java order
    let mut unique_alleles: JavaHashSet<String> = JavaHashSet::new();
    for c in &coverages {
        unique_alleles.insert(c.allele.clone());
    }
    let purple_dir = config.purple_dir.replace('*', sample_id);
    let purity = read_purity_context(&purple_dir, sample_id).map_err(|e| {
        error!("failed to read sample({}) purity data: {:#}", sample_id, e);
        anyhow!("NullPointerException: no purity context")
    })?;
    let sample_ploidy = purity.ploidy;

    let mut sample_tpms: HashMap<String, f64> = HashMap::new();
    let tpm_source;
    if sample.has_rna {
        match &reference.transcript_expression {
            None => {
                debug!("sample({}) loading transcript expression", sample.rna_id_str());
                match load_transcript_expression(&config.isofox_dir, sample.rna_id_str()) {
                    Ok(m) => sample_tpms = m,
                    Err(e) => {
                        error!("failed to load sample({}) transcript expression: {:#}", sample.rna_id_str(), e);
                        std::process::exit(1);
                    }
                }
            }
            Some(m) => {
                if !m.has_sample_id(sample_id) {
                    error!("sample({}) missing from transcript expression matrix", sample_id);
                    std::process::exit(1);
                }
            }
        }
        tpm_source = TpmSource::Sample;
    } else if reference.tpm_medians.has_cancer_type(&sample.cancer_type) {
        tpm_source = TpmSource::CancerType;
    } else {
        tpm_source = TpmSource::Cohort;
    }

    let rna_neos = load_rna_neo_data(sample, &config.isofox_dir)?;
    let mut somatic_variants = Some(Vec::new());
    if !config.rna_somatic_vcf.is_empty() {
        let point_neos: Vec<&NeoEpitopeData> = neos.iter().filter(|x| x.variant_type.is_point_mutation()).collect();
        somatic_variants = load_somatic_variants(sample, &config.rna_somatic_vcf, &point_neos)?;
    }
    if sample.has_rna && (rna_neos.is_none() || somatic_variants.is_none()) {
        error!(
            "sample({}) missing required RNA: fusions({}) variants({})",
            sample_id,
            if rna_neos.is_none() { "missing" } else { "present" },
            if somatic_variants.is_none() { "missing" } else { "present" }
        );
        std::process::exit(1);
    }
    for n in neos.iter_mut() {
        n.set_expression_data(sample, &sample_tpms, reference.transcript_expression.as_ref(), &reference.tpm_medians);
        if n.variant_type.is_fusion() {
            let r = rna_neos.as_deref().ok_or_else(|| anyhow!("NullPointerException: no RNA neo-epitope data"))?;
            n.set_fusion_rna_support(r);
        }
        if n.variant_type.is_point_mutation() {
            let v = somatic_variants.as_deref().ok_or_else(|| anyhow!("NullPointerException: no RNA somatic variants"))?;
            n.set_mutation_rna_support(v);
        }
    }

    let mut peptides: Vec<PeptideScoreData> = Vec::new();
    let mut calc = TpmCalculator::new(FLANK_AA_COUNT, config.peptide_length_range);
    calc.compute(sample_id, &mut neos, sample_ploidy, &mut peptides);

    // drop peptides found in an upstream transcript's protein, give the rest the sample's alleles
    let mut to_score: Vec<usize> = Vec::new();
    for n in neos.iter_mut() {
        let mut i = 0;
        while i < n.peptides.len() {
            let k = n.peptides[i];
            if n.transcripts[FS_UP].iter().any(|t| reference.peptide_matches_wildtype(&peptides[k].peptide, t)) {
                n.peptides.remove(i);
                continue;
            }
            if peptides[k].alleles.is_empty() {
                for a in unique_alleles.iter() {
                    peptides[k].add_allele(a);
                }
                if !peptides[k].alleles.is_empty() {
                    to_score.push(k);
                }
            }
            i += 1;
        }
    }
    let score_count = score_peptides(&reference.scorer, &mut peptides, &to_score, config.threads);
    debug!("sample({}) neoepitopes({}) scored {} allele-peptides", sample_id, neos.len(), score_count);

    if config.write_types.contains(&OutputType::AllelePeptide) {
        for n in &neos {
            writers.write_peptide_data(sample_id, n, &mut peptides, &coverages);
        }
    }
    if config.write_types.contains(&OutputType::Neoepitope) {
        for n in &neos {
            writers.write_neo_data(sample_id, tpm_source, n);
        }
    }
    Ok(())
}

/// BindScorer.calcScoreData for every allele of the given peptides, on up to `threads` threads; the
/// recognition similarities of a peptide are found once for all its alleles
fn score_peptides(scorer: &BindScorer, peptides: &mut [PeptideScoreData], to_score: &[usize], threads: i32) -> usize {
    let mut work: Vec<(usize, Vec<BindData>)> = to_score.iter().map(|&k| (k, std::mem::take(&mut peptides[k].alleles))).collect();
    let count: usize = work.iter().map(|(_, v)| v.len()).sum();
    let score_one = |alleles: &mut Vec<BindData>| {
        let needs_similarity = alleles.iter().any(|b| scorer.has_matrix(&b.allele, b.peptide_length()));
        let sims = if needs_similarity { scorer.peptide_similarities(&alleles[0].peptide) } else { None };
        for b in alleles.iter_mut() {
            scorer.calc_score_data(b, sims.as_ref());
        }
    };
    let threads = threads.max(1) as usize;
    if threads <= 1 || work.len() < 2 {
        for (_, v) in work.iter_mut() {
            score_one(v);
        }
    } else {
        let chunk = work.len().div_ceil(threads);
        std::thread::scope(|s| {
            for part in work.chunks_mut(chunk) {
                s.spawn(move || {
                    for (_, v) in part.iter_mut() {
                        score_one(v);
                    }
                });
            }
        });
    }
    for (k, v) in work {
        peptides[k].alleles = v;
    }
    count
}
