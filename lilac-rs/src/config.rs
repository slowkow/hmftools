//! LilacConfig and the command line, parsed as the Java ConfigBuilder does: `-name value` and `-flag`,
//! required items, path checks, integer and decimal checks, `-version` and `-help` as the sole argument.
//! Also `Params`: the LilacConstants values the config overrides, and the sequencing type.

use anyhow::{bail, Result};
use std::collections::HashMap;

use hmf_common::java_text::{parse_double, parse_int};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ItemType {
    Str,
    Path,
    Integer,
    Decimal,
    Flag,
}

impl ItemType {
    fn name(self) -> &'static str {
        match self {
            ItemType::Str => "STRING",
            ItemType::Path => "PATH",
            ItemType::Integer => "INTEGER",
            ItemType::Decimal => "DECIMAL",
            ItemType::Flag => "FLAG",
        }
    }
}

struct Item {
    name: &'static str,
    ty: ItemType,
    required: bool,
    desc: &'static str,
    default: Option<&'static str>,
}

/// The items LilacConfig.addConfig registers, in that order.
fn items() -> Vec<Item> {
    use ItemType::*;
    let item = |name, ty, required, desc, default| Item { name, ty, required, desc, default };
    vec![
        item("sample", Str, true, "Name of sample", None),
        item("sample_data_dir", Path, false, "Path to sample pipeline files", None),
        item("reference_bam", Path, false, "Reference BAM file", None),
        item("tumor_bam", Path, false, "Tumor BAM file", None),
        item("rna_bam", Path, false, "RNA BAM file", None),
        item("resource_dir", Path, true, "Path to resource files", None),
        item("genes", Str, false, "Gene sets to use", Some("ALL")),
        item("sequencing_type", Str, false, "Sequencing types: ILLUMINA, ULTIMA, SBX", Some("ILLUMINA")),
        item("targeted_panel", Flag, false, "Run with targeted panel settings", None),
        item("min_evidence_factor", Decimal, false, "Min fragment evidence required", Some("0.006")),
        item("max_ref_fragments", Integer, false, "Cap ref fragments in solution search, 0 uses all", Some("1000")),
        item("min_high_qual_evidence_factor", Decimal, false, "Min high-qual fragment evidence factor", Some("0.003")),
        item("min_evidence_support", Integer, false, "Min fragment evidence depth", Some("1")),
        item("min_depth_filter", Integer, false, "Raw depth below which we do not filter", Some("10")),
        item("min_fragments_per_allele", Integer, false, "Min fragments per allele", Some("7")),
        item("min_fragments_to_remove_single", Integer, false, "Min fragments to remote single", Some("40")),
        item(
            "fatal_total_low_coverage_positions",
            Integer,
            false,
            "Fatal total low coverage positions across all HLA loci",
            Some("300"),
        ),
        item("hla_y_threshold", Decimal, false, "HLA-Y percent threshold", Some("0.003")),
        item("freq_score_penalty", Decimal, false, "Population frequency penalty", Some("9.0E-4")),
        item("top_score_threshold", Decimal, false, "Max distance from top score", Some("0.005")),
        item("actual_alleles", Str, false, "Comma separated known actual alleles for the sample", None),
        item("restricted_alleles", Str, false, "Comma separated restricted analysis allele list", None),
        item(
            "max_elim_candidates",
            Integer,
            false,
            "Revert to only common alleles if candidate allele count exceeds this after elimination",
            Some("0"),
        ),
        item("gene_copy_number", Path, false, "Path to gene copy number file", None),
        item("purple_dir", Path, false, "Path to Purple pipeline files", None),
        item("somatic_vcf", Path, false, "Path to sample Purple somatic VCF", None),
        item("debug_phasing", Flag, false, "More detailed logging of phasing", None),
        item("run_validation", Flag, false, "Run validation checks", None),
        item("log_perf", Flag, false, "Log performance metrics", None),
        item("stack_sample_rate", Integer, false, "Stack sampling rate", Some("0")),
        item("write_types", Str, false, "Write types: SUMMARY (default), FRAGMENTS, READS, REF_COUNTS or ALL", Some("SUMMARY")),
        item("bam_validation", Str, false, "BAM validation: STRICT (default), LENIENT or SILENT", Some("SILENT")),
        item("ref_genome_version", Str, true, "Ref genome version, 37 or 38", None),
        item("ref_genome", Path, true, "Path to reference genome fasta files", None),
        item("output_dir", Str, false, "Output directory", None),
        item("threads", Integer, false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        item("log_debug", Flag, false, "Log at DEBUG level", None),
        item("log_level", Str, false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", None),
    ]
}

/// The parsed command line: values by name and the flags given.
pub struct ConfigBuilder {
    items: Vec<Item>,
    values: HashMap<&'static str, String>,
}

pub enum Parsed {
    Run(ConfigBuilder),
    /// -version or -help was printed; exit 0
    Exit,
}

impl ConfigBuilder {
    /// ConfigBuilder.checkAndParseCommandLine: an error (logged) for anything Java rejects.
    pub fn parse(argv: &[String]) -> Result<Parsed> {
        let items = items();
        if argv.len() == 1 && (argv[0] == "-version" || argv[0] == "-help") {
            // a closed pipe (lilac -version | sed ...) is not an error
            use std::io::Write;
            let mut out = std::io::stdout().lock();
            let _ = writeln!(out, "Lilac version {}", crate::VERSION);
            if argv[0] == "-help" {
                for line in item_lines(&items) {
                    let _ = writeln!(out, "{line}");
                }
            }
            return Ok(Parsed::Exit);
        }
        let mut values: HashMap<&'static str, String> = HashMap::new();
        let mut matched: Option<usize> = None;
        let invalid = |items: &[Item], msg: String| -> anyhow::Error {
            crate::error!("{msg}");
            for line in item_lines(items) {
                crate::info!("{line}");
            }
            anyhow::anyhow!("invalid command line")
        };
        for arg in argv {
            if let Some(i) = matched {
                let ty = items[i].ty;
                if arg.starts_with('-') && ty != ItemType::Decimal && ty != ItemType::Integer {
                    return Err(invalid(&items, format!("config item({}) has invalid argument: {arg}", items[i].name)));
                }
                values.insert(items[i].name, arg.clone());
                matched = None;
            } else if let Some(name) = arg.strip_prefix('-') {
                let Some(i) = items.iter().position(|it| it.name == name) else {
                    return Err(invalid(&items, format!("unregistered config item: {arg}")));
                };
                if items[i].ty == ItemType::Flag {
                    values.insert(items[i].name, "true".to_string());
                } else {
                    matched = Some(i);
                }
            } else {
                return Err(invalid(&items, format!("expecting config name but encountered invalid argument: {arg}")));
            }
        }
        let mut ok = true;
        for item in &items {
            let value = values.get(item.name);
            if item.required && value.is_none() {
                crate::error!("missing required config: {}: {}", item.name, item.desc);
                ok = false;
                continue;
            }
            let Some(value) = value else { continue };
            match item.ty {
                ItemType::Path => {
                    if !std::path::Path::new(value).exists() {
                        crate::error!("invalid path for config: {} = {}", item.name, value);
                        ok = false;
                    }
                }
                ItemType::Integer if parse_int(value).is_none() => {
                    crate::error!("invalid type for config: {}: {} = {}", item.name, item.ty.name(), value);
                    ok = false;
                }
                ItemType::Decimal if parse_double(value).is_none() => {
                    crate::error!("invalid type for config: {}: {} = {}", item.name, item.ty.name(), value);
                    ok = false;
                }
                _ => {}
            }
        }
        if !ok {
            bail!("invalid command line");
        }
        Ok(Parsed::Run(ConfigBuilder { items, values }))
    }

    pub fn has_flag(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }
    pub fn has_value(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }
    /// getValue: the value given, or the default
    pub fn value(&self, name: &str) -> Option<String> {
        if let Some(v) = self.values.get(name) {
            return Some(v.clone());
        }
        self.items.iter().find(|i| i.name == name).and_then(|i| i.default.map(|d| d.to_string()))
    }
    pub fn value_or(&self, name: &str, default: &str) -> String {
        self.value(name).unwrap_or_else(|| default.to_string())
    }
    pub fn integer(&self, name: &str) -> Result<i32> {
        let v = self.value(name).unwrap_or_default();
        parse_int(&v).ok_or_else(|| anyhow::anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
    pub fn decimal(&self, name: &str) -> Result<f64> {
        let v = self.value(name).unwrap_or_default();
        parse_double(&v).ok_or_else(|| anyhow::anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
}

/// ConfigBuilder.printItems: one line per item, sorted, after a title line
fn item_lines(items: &[Item]) -> Vec<String> {
    let mut lines: Vec<String> = items
        .iter()
        .map(|item| {
            let mut s = format!(" -{}:", item.name);
            if item.ty != ItemType::Str {
                s.push_str(&format!(" type({})", item.ty.name()));
            }
            if let Some(d) = item.default.filter(|_| item.ty != ItemType::Flag) {
                s.push_str(&format!(" default({d})"));
            }
            if item.required {
                s.push_str(" REQUIRED");
            }
            s.push_str(&format!(" desc({})", item.desc));
            s
        })
        .collect();
    lines.sort();
    lines.insert(0, "registered config items:".to_string());
    lines
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum SequencingType {
    Illumina,
    Ultima,
    Sbx,
}

/// The LilacConstants values that the config sets or overrides, and the global sequencing type.
#[derive(Clone, Debug)]
pub struct Params {
    pub min_evidence_factor: f64,
    pub min_high_qual_evidence_factor: f64,
    pub min_evidence_support: i32,
    pub min_depth_filter: i32,
    pub freq_score_penalty: f64,
    pub fail_low_coverage_threshold: i32,
    pub warn_low_coverage_threshold: f64,
    pub warn_unmatched_haplotype_support: f64,
    pub warn_indel_threshold: f64,
    pub sequencing_type: SequencingType,
    /// HLA_CHR: chr6 or 6
    pub hla_chr: String,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            min_evidence_factor: 0.006,
            min_high_qual_evidence_factor: 0.003,
            min_evidence_support: 1,
            min_depth_filter: 10,
            freq_score_penalty: 0.0009,
            fail_low_coverage_threshold: 360,
            warn_low_coverage_threshold: 180.0,
            warn_unmatched_haplotype_support: 0.01,
            warn_indel_threshold: 0.005,
            sequencing_type: SequencingType::Illumina,
            hla_chr: "6".to_string(),
        }
    }
}

impl Params {
    pub fn is_illumina(&self) -> bool {
        self.sequencing_type == SequencingType::Illumina
    }
    pub fn is_ultima(&self) -> bool {
        self.sequencing_type == SequencingType::Ultima
    }
    pub fn is_sbx(&self) -> bool {
        self.sequencing_type == SequencingType::Sbx
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum WriteType {
    Summary,
    Fragments,
    RefCounts,
}

#[derive(Clone, Debug)]
pub struct LilacConfig {
    pub sample: String,
    pub reference_bam: String,
    pub tumor_bam: String,
    pub rna_bam: String,
    pub resource_dir: String,
    pub ref_genome: String,
    pub v38: bool,
    pub sample_data_dir: String,
    /// with the trailing separator
    pub output_dir: String,
    pub genes: Vec<crate::hla::GeneSelector>,
    pub hla_y_percent_threshold: f64,
    pub min_fragments_per_allele: i32,
    pub min_fragments_to_remove_single: i32,
    pub max_ref_fragments: i32,
    pub threads: i32,
    pub top_score_threshold: f64,
    pub copy_number_file: String,
    pub somatic_variants_file: String,
    pub debug_phasing: bool,
    pub run_validation: bool,
    pub fatal_total_low_coverage_positions: i32,
    pub max_elimination_candidates: i32,
    pub log_perf: bool,
    /// the allele names given (four-digit form taken when the reference data is loaded)
    pub actual_alleles: Vec<String>,
    pub restricted_alleles: Vec<String>,
    pub write_types: Vec<WriteType>,
    pub params: Params,
}

fn check_add_dir_separator(d: &str) -> String {
    if d.is_empty() || d.ends_with('/') { d.to_string() } else { format!("{d}/") }
}

fn file_if_exists(f: String) -> String {
    if std::path::Path::new(&f).exists() { f } else { String::new() }
}

/// FileWriterUtils.pathFromFile: the parent directory with a separator, or "./"
fn path_from_file(f: &str) -> String {
    if !f.contains('/') {
        return "./".to_string();
    }
    match std::path::Path::new(f).parent() {
        Some(p) if !p.as_os_str().is_empty() => format!("{}/", p.display()),
        Some(_) => "/".to_string(),
        None => "null/".to_string(),
    }
}

impl LilacConfig {
    /// The LilacConfig constructor; errors are those Java throws from it
    pub fn from(cb: &ConfigBuilder) -> Result<LilacConfig> {
        let sample = cb.value("sample").unwrap_or_default();
        let sample_data_dir;
        let output_dir;
        let reference_bam;
        let tumor_bam;
        let rna_bam;
        if cb.has_value("sample_data_dir") {
            sample_data_dir = check_add_dir_separator(&cb.value("sample_data_dir").unwrap());
            output_dir = if cb.has_value("output_dir") { check_add_dir_separator(&cb.value("output_dir").unwrap()) } else { sample_data_dir.clone() };
            let Some(t) = sample.rfind('T') else {
                bail!("java.lang.StringIndexOutOfBoundsException: begin 0, end -1, length {}", sample.len());
            };
            reference_bam = format!("{sample_data_dir}{}R.hla.bam", &sample[..t]);
            tumor_bam = file_if_exists(format!("{sample_data_dir}{sample}.hla.bam"));
            rna_bam = file_if_exists(format!("{sample_data_dir}{sample}.rna.hla.bam"));
        } else {
            sample_data_dir = String::new();
            reference_bam = cb.value_or("reference_bam", "");
            tumor_bam = cb.value_or("tumor_bam", "");
            rna_bam = cb.value_or("rna_bam", "");
            output_dir = if cb.has_value("output_dir") { check_add_dir_separator(&cb.value("output_dir").unwrap()) } else { path_from_file(&reference_bam) };
        }
        let (somatic_variants_file, copy_number_file) = if cb.has_value("somatic_vcf") && cb.has_value("gene_copy_number") {
            (cb.value("somatic_vcf").unwrap(), cb.value("gene_copy_number").unwrap())
        } else if cb.has_value("purple_dir") {
            let d = check_add_dir_separator(&cb.value("purple_dir").unwrap());
            (format!("{d}{sample}.purple.somatic.vcf.gz"), format!("{d}{sample}.purple.cnv.gene.tsv"))
        } else if !sample_data_dir.is_empty() {
            let cn = file_if_exists(format!("{sample_data_dir}{sample}.purple.cnv.gene.tsv"));
            let sv = file_if_exists(format!("{sample_data_dir}{sample}.purple.somatic.vcf.gz"));
            (sv, cn)
        } else {
            (String::new(), String::new())
        };
        let resource_dir = check_add_dir_separator(&cb.value_or("resource_dir", ""));
        let ref_genome = cb.value_or("ref_genome", "");
        let v38 = match cb.value("ref_genome_version") {
            None => false,
            Some(v) => match v.as_str() {
                "V37" | "37" | "HG37" => false,
                "V38" | "38" | "HG38" => true,
                _ => bail!("java.lang.IllegalArgumentException: Cannot resolve ref genome version: {v}"),
            },
        };
        let mut params = Params::default();
        let seq_type = cb.value("sequencing_type").unwrap();
        params.sequencing_type = match seq_type.as_str() {
            "ILLUMINA" => SequencingType::Illumina,
            "ULTIMA" => SequencingType::Ultima,
            "SBX" => SequencingType::Sbx,
            _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.sequencing.SequencingType.{seq_type}"),
        };
        if cb.has_flag("targeted_panel") || params.sequencing_type != SequencingType::Illumina {
            params.warn_unmatched_haplotype_support *= 2.0;
            params.warn_indel_threshold *= 2.0;
        }
        let genes = crate::hla::GeneSelector::parse_arg(&cb.value("genes").unwrap()).map_err(|e| anyhow::anyhow!(e))?;
        params.min_evidence_factor = cb.decimal("min_evidence_factor")?;
        params.min_high_qual_evidence_factor = cb.decimal("min_high_qual_evidence_factor")?;
        params.min_evidence_support = cb.integer("min_evidence_support")?;
        params.min_depth_filter = cb.integer("min_depth_filter")?;
        params.freq_score_penalty = cb.decimal("freq_score_penalty")?;
        params.hla_chr = if v38 { "chr6".to_string() } else { "6".to_string() };
        // parseAlleleList: split at ';' (ITEM_DELIM), although the options' descriptions say "comma separated"
        let parse_list = |name: &str| -> Vec<String> {
            match cb.value(name) {
                Some(s) if !s.is_empty() => hmf_common::java_text::split_all(&s, ';').into_iter().map(|x| x.to_string()).collect(),
                _ => Vec::new(),
            }
        };
        let write_types_str = cb.value("write_types").unwrap();
        let write_types = if write_types_str == "ALL" {
            vec![WriteType::Summary, WriteType::Fragments, WriteType::RefCounts]
        } else {
            let mut v = Vec::new();
            for t in hmf_common::java_text::split(&write_types_str, ',') {
                let w = match t {
                    "SUMMARY" => WriteType::Summary,
                    "FRAGMENTS" => WriteType::Fragments,
                    "REF_COUNTS" => WriteType::RefCounts,
                    _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.lilac.ResultsWriter.WriteType.{t}"),
                };
                if !v.contains(&w) {
                    v.push(w);
                }
            }
            v
        };
        let bam_validation = cb.value("bam_validation").unwrap();
        if !matches!(bam_validation.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            bail!("java.lang.IllegalArgumentException: No enum constant htsjdk.samtools.ValidationStringency.{bam_validation}");
        }
        let config = LilacConfig {
            sample,
            reference_bam,
            tumor_bam,
            rna_bam,
            resource_dir,
            ref_genome,
            v38,
            sample_data_dir,
            output_dir,
            genes,
            hla_y_percent_threshold: cb.decimal("hla_y_threshold")?,
            min_fragments_per_allele: cb.integer("min_fragments_per_allele")?,
            min_fragments_to_remove_single: cb.integer("min_fragments_to_remove_single")?,
            max_ref_fragments: cb.integer("max_ref_fragments")?,
            threads: cb.integer("threads")?,
            top_score_threshold: crate::math::java_min(cb.decimal("top_score_threshold")?, 0.5),
            copy_number_file,
            somatic_variants_file,
            debug_phasing: cb.has_flag("debug_phasing"),
            run_validation: cb.has_flag("run_validation"),
            fatal_total_low_coverage_positions: cb.integer("fatal_total_low_coverage_positions")?,
            max_elimination_candidates: cb.integer("max_elim_candidates")?,
            log_perf: cb.has_flag("log_perf"),
            actual_alleles: parse_list("actual_alleles"),
            restricted_alleles: parse_list("restricted_alleles"),
            write_types,
            params,
        };
        cb.integer("stack_sample_rate")?;
        Ok(config)
    }

    /// the config of the Java unit tests (LilacConfig(sampleId))
    #[cfg(test)]
    pub fn for_tests(sample: &str) -> LilacConfig {
        LilacConfig {
            sample: sample.to_string(),
            reference_bam: String::new(),
            tumor_bam: String::new(),
            rna_bam: String::new(),
            resource_dir: String::new(),
            ref_genome: String::new(),
            v38: false,
            sample_data_dir: String::new(),
            output_dir: String::new(),
            genes: vec![crate::hla::GeneSelector::MhcClass1],
            hla_y_percent_threshold: 0.003,
            min_fragments_per_allele: 7,
            min_fragments_to_remove_single: 40,
            max_ref_fragments: 1000,
            threads: 0,
            top_score_threshold: 0.005,
            copy_number_file: String::new(),
            somatic_variants_file: String::new(),
            debug_phasing: false,
            run_validation: true,
            fatal_total_low_coverage_positions: 600,
            max_elimination_candidates: 0,
            log_perf: false,
            actual_alleles: Vec::new(),
            restricted_alleles: Vec::new(),
            write_types: vec![WriteType::Summary],
            params: Params::default(),
        }
    }

    pub fn tumor_only(&self) -> bool {
        self.reference_bam.is_empty() && !self.tumor_bam.is_empty()
    }
    pub fn allele_calling_bam(&self) -> &str {
        if !self.reference_bam.is_empty() { &self.reference_bam } else { &self.tumor_bam }
    }
    /// formFileId: <output dir><sample>.lilac.<id>
    pub fn form_file_id(&self, id: &str) -> String {
        format!("{}{}.lilac.{}", self.output_dir, self.sample, id)
    }
    pub fn writes(&self, w: WriteType) -> bool {
        self.write_types.contains(&w)
    }
}
