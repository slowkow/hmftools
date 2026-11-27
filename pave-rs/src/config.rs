//! PaveConfig and the command line, parsed as the Java ConfigBuilder does: `-name value` and `-flag`,
//! required items, path checks, integer checks, `-version` and `-help` as the sole argument.

use anyhow::{bail, Result};
use std::collections::HashMap;

use hmf_common::java_text;
use hmf_common::regions::{RefGenomeVersion, SpecificRegions};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ItemType {
    Str,
    Path,
    Integer,
    Flag,
}

impl ItemType {
    /// ConfigItemType names as ConfigBuilder prints them
    fn name(self) -> &'static str {
        match self {
            ItemType::Str => "STRING",
            ItemType::Path => "PATH",
            ItemType::Integer => "INTEGER",
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

/// The items PaveConfig.addConfig registers, in that order.
fn items() -> Vec<Item> {
    use ItemType::*;
    let item = |name, ty, required, desc, default| Item { name, ty, required, desc, default };
    vec![
        item("sample", Str, true, "Sample ID", None),
        item("input_vcf", Path, true, "VCF input file", None),
        item("output_vcf", Str, false, "VCF output file (optional), otherwise will append 'pave' suffix to input filename", None),
        item("ref_genome", Path, true, "Path to reference genome fasta files", None),
        item("ref_genome_version", Str, true, "Ref genome version, 37 or 38", None),
        item("ensembl_data_dir", Path, true, "Ensembl data file directory", None),
        item("driver_gene_panel", Path, false, "Path to driver gene panel", None),
        item("sequencing_type", Str, false, "Sequencing types: ILLUMINA, ULTIMA, SBX", Some("ILLUMINA")),
        item("pon_file", Path, false, "PON entries", None),
        item("pon_artefact_file", Path, false, "PON artefacts to filter", None),
        item("pon_filters", Str, false, "PON filters per tier, format: TIER:MAX_SAMPLES:MAX_COUNT separated by ';'", None),
        item("write_transcript_data", Flag, false, "Write variant impacts per transcript to TSV", None),
        item("only_canonical", Flag, false, "Only check canonical transcripts", None),
        item("filter_type", Str, false, "Variant filters: ALL;PASS;PANEL separated by ';'", None),
        item("write_pass_only", Flag, false, "Only annotate passing variants", None),
        item("set_reportable", Flag, false, "Set reportable and hotspot flags", None),
        item("write_detailed", Flag, false, "Write detailed transcript impact info", None),
        item("gnomad_freq_file", Path, false, "Gnomad frequency file", None),
        item("gnomad_freq_dir", Path, false, "Gnomad frequency directory", None),
        item("gnomad_no_filter", Flag, false, "No Gnomad filter is applied", None),
        item("mappability_bed", Path, false, "Mappability BED file", None),
        item("clinvar_vcf", Path, false, "Clinvar annotation VCF", None),
        item("blacklist_bed", Path, false, "Blacklist BED file", None),
        item("blacklist_vcf", Path, false, "Blacklist VCF file", None),
        item("threads", Integer, false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        item("specific_chr", Str, false, "Restrict to chromosome(s) separated by ';'", None),
        item(
            "specific_regions",
            Str,
            false,
            "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd",
            None,
        ),
        item("output_dir", Str, false, "Output directory", None),
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
            println!("Pave version {}", crate::VERSION);
            if argv[0] == "-help" {
                for line in item_lines(&items) {
                    println!("{line}");
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
                if arg.starts_with('-') && ty != ItemType::Integer {
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
        // a trailing name without a value is left unset, as in Java
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
                ItemType::Integer if java_text::parse_int(value).is_none() => {
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

    #[allow(dead_code)]
    pub fn has_value(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }
    pub fn has_flag(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }
    /// getValue: the value given, or the default
    pub fn value(&self, name: &str) -> Option<String> {
        if let Some(v) = self.values.get(name) {
            return Some(v.clone());
        }
        self.items.iter().find(|i| i.name == name).and_then(|i| i.default.map(|d| d.to_string()))
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
            if let Some(d) = &item.default {
                if item.ty != ItemType::Flag {
                    s.push_str(&format!(" default({d})"));
                }
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

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SequencingType {
    Illumina,
    Ultima,
    Sbx,
}

/// FilterType
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FilterType {
    All,
    Pass,
    Panel,
}

/// The options of one run (PaveConfig, with the reference data options ReferenceData reads)
pub struct PaveConfig {
    pub sample_id: String,
    pub vcf_file: String,
    pub ref_gen_version: RefGenomeVersion,
    pub output_dir: String,
    pub output_vcf_file: Option<String>,
    pub write_transcript_file: bool,
    pub only_canonical: bool,
    pub filter: FilterType,
    pub write_pass_only: bool,
    pub write_detailed: bool,
    pub set_reportable: bool,
    pub specific_regions: SpecificRegions,
    pub threads: usize,
    pub sequencing_type: SequencingType,

    pub ref_genome: String,
    pub ensembl_data_dir: String,
    pub driver_gene_panel: Option<String>,
    pub pon_file: Option<String>,
    pub pon_artefact_file: Option<String>,
    pub pon_filters: Option<String>,
    pub gnomad_freq_file: Option<String>,
    pub gnomad_freq_dir: Option<String>,
    pub gnomad_no_filter: bool,
    pub mappability_bed: Option<String>,
    pub clinvar_vcf: Option<String>,
    pub blacklist_bed: Option<String>,
    pub blacklist_vcf: Option<String>,
}

/// FileWriterUtils.checkAddDirSeparator
pub fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

/// FileWriterUtils.pathFromFile: the file's parent directory with a separator, "./" for a bare name
fn path_from_file(file: &str) -> String {
    if !file.contains('/') {
        return "./".to_string();
    }
    // File.getParent: the path before the last separator (trailing separators ignored)
    let trimmed = file.trim_end_matches('/');
    match trimmed.rfind('/') {
        Some(0) => "//".to_string(),
        Some(i) => format!("{}/", &trimmed[..i]),
        None => "null/".to_string(),
    }
}

impl PaveConfig {
    pub fn from(cb: &ConfigBuilder) -> Result<PaveConfig> {
        let sample_id = cb.value("sample").unwrap_or_default();
        let vcf_file = cb.value("input_vcf").unwrap_or_default();
        let output_vcf_file = cb.value("output_vcf");
        let version_str = cb.value("ref_genome_version").unwrap_or_default();
        let ref_gen_version = RefGenomeVersion::from_java(&version_str)?;
        let filter = match cb.value("filter_type") {
            None => FilterType::Pass,
            Some(v) => match v.as_str() {
                "ALL" => FilterType::All,
                "PASS" => FilterType::Pass,
                "PANEL" => FilterType::Panel,
                _ => bail!("No enum constant com.hartwig.hmftools.pave.FilterType.{v}"),
            },
        };
        let seq = cb.value("sequencing_type").unwrap_or_default();
        let sequencing_type = match seq.as_str() {
            "ILLUMINA" => SequencingType::Illumina,
            "ULTIMA" => SequencingType::Ultima,
            "SBX" => SequencingType::Sbx,
            _ => bail!("No enum constant com.hartwig.hmftools.common.sequencing.SequencingType.{seq}"),
        };
        // TaskExecutor.parseThreads: at least 1
        let threads = java_text::parse_int(&cb.value("threads").unwrap_or_default()).unwrap_or(1).max(1) as usize;
        let specific_regions = SpecificRegions::from_java(cb.value("specific_regions").as_deref(), cb.value("specific_chr").as_deref())?;
        let output_dir = if let Some(dir) = cb.value("output_dir") {
            check_add_dir_separator(&dir)
        } else {
            path_from_file(output_vcf_file.as_deref().unwrap_or(&vcf_file))
        };
        Ok(PaveConfig {
            write_transcript_file: cb.has_flag("write_transcript_data"),
            only_canonical: cb.has_flag("only_canonical"),
            filter,
            write_pass_only: cb.has_flag("write_pass_only"),
            write_detailed: cb.has_flag("write_detailed"),
            set_reportable: cb.has_flag("set_reportable"),
            specific_regions,
            threads,
            sequencing_type,
            ref_genome: cb.value("ref_genome").unwrap_or_default(),
            ensembl_data_dir: cb.value("ensembl_data_dir").unwrap_or_default(),
            driver_gene_panel: cb.value("driver_gene_panel"),
            pon_file: cb.value("pon_file"),
            pon_artefact_file: cb.value("pon_artefact_file"),
            pon_filters: cb.value("pon_filters"),
            gnomad_freq_file: cb.value("gnomad_freq_file"),
            gnomad_freq_dir: cb.value("gnomad_freq_dir"),
            gnomad_no_filter: cb.has_flag("gnomad_no_filter"),
            mappability_bed: cb.value("mappability_bed"),
            clinvar_vcf: cb.value("clinvar_vcf"),
            blacklist_bed: cb.value("blacklist_bed"),
            blacklist_vcf: cb.value("blacklist_vcf"),
            sample_id,
            vcf_file,
            ref_gen_version,
            output_dir,
            output_vcf_file,
        })
    }

    /// requireIndex: a .vcf.gz input is read by chromosome through its index
    pub fn require_index(&self) -> bool {
        self.vcf_file.ends_with(".vcf.gz")
    }

    pub fn is_sbx(&self) -> bool {
        self.sequencing_type == SequencingType::Sbx
    }
    pub fn is_ultima(&self) -> bool {
        self.sequencing_type == SequencingType::Ultima
    }

    /// PaveApplication.initialiseVcfWriter: the output VCF, by default the input's name with ".pave"
    /// before ".vcf", always zipped
    pub fn output_vcf(&self) -> String {
        if let Some(f) = &self.output_vcf_file {
            return f.clone();
        }
        let filename = self.vcf_file.rsplit('/').next().unwrap_or(&self.vcf_file);
        let i = filename.find(".vcf").unwrap_or(filename.len());
        let mut out = format!("{}{}.pave{}", self.output_dir, &filename[..i], &filename[i..]);
        if !out.ends_with(".gz") {
            out.push_str(".gz");
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_like_config_builder() {
        let dir = std::env::temp_dir();
        let d = dir.to_str().unwrap();
        let argv = strings(&[
            "-sample", "S", "-input_vcf", d, "-ref_genome", d, "-ref_genome_version", "38", "-ensembl_data_dir", d, "-threads", "-2",
            "-gnomad_no_filter",
        ]);
        let Parsed::Run(cb) = ConfigBuilder::parse(&argv).unwrap() else { panic!() };
        assert!(cb.has_flag("gnomad_no_filter"));
        assert_eq!(cb.value("sequencing_type").unwrap(), "ILLUMINA");
        let c = PaveConfig::from(&cb).unwrap();
        assert_eq!(c.threads, 1);
        assert_eq!(c.filter, FilterType::Pass);
        assert!(ConfigBuilder::parse(&strings(&["-sample", "S"])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-unknown", "S"])).is_err());
    }

    #[test]
    fn output_names() {
        let dir = std::env::temp_dir();
        let d = dir.to_str().unwrap();
        let argv = strings(&["-sample", "S", "-input_vcf", d, "-ref_genome", d, "-ref_genome_version", "38", "-ensembl_data_dir", d]);
        let Parsed::Run(cb) = ConfigBuilder::parse(&argv).unwrap() else { panic!() };
        let mut c = PaveConfig::from(&cb).unwrap();
        c.vcf_file = "/a/b/S.sage.somatic.vcf".to_string();
        c.output_dir = "/o/".to_string();
        assert_eq!(c.output_vcf(), "/o/S.sage.somatic.pave.vcf.gz");
        c.vcf_file = "/a/b/S.vcf.gz".to_string();
        assert_eq!(c.output_vcf(), "/o/S.pave.vcf.gz");
        assert_eq!(path_from_file("/a/b/S.vcf.gz"), "/a/b/");
        assert_eq!(path_from_file("S.vcf.gz"), "./");
    }
}
