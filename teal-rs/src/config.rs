//! The command line, parsed as the Java ConfigBuilder does (`-name value`, `-flag`, required items, path,
//! integer and decimal checks, `-version` and `-help` as the sole argument), and the parameters of the
//! three entry points: TealCommonParams (TealPipelineTelbamApp), TealPipelineParams (TealPipelineApp) and
//! TealParams (TealApplication).

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;

use crate::tellength::SequencingType;

pub const APP_NAME: &str = "Teal";
pub const VERSION: &str = "1.4";

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

/// Which entry point
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum App {
    /// TealPipelineTelbamApp
    Prep,
    /// TealPipelineApp
    Pipeline,
    /// TealApplication (the jar's main class)
    Standalone,
}

fn item(name: &'static str, ty: ItemType, required: bool, desc: &'static str, default: Option<&'static str>) -> Item {
    Item { name, ty, required, desc, default }
}

/// TealCommonParams.registerConfig
fn common_items() -> Vec<Item> {
    use ItemType::*;
    vec![
        item("reference", Str, false, "Reference ID", None),
        item("reference_bam", Path, false, "Reference BAM file", None),
        item("tumor", Str, false, "Tumor ID", None),
        item("tumor_bam", Path, false, "Tumor BAM file", None),
        item("ref_genome", Path, false, "Path to reference genome fasta files", None),
        item("output_dir", Str, false, "Output directory", None),
        item("log_debug", Flag, false, "Log at DEBUG level", None),
        item("log_level", Str, false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", None),
        item("threads", Integer, false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        item("ref_genome_version", Str, false, "Ref genome version, 37 or 38", Some("V37")),
        item("sequencing_type", Str, false, "Sequencing types: ILLUMINA, ULTIMA, SBX", Some("ILLUMINA")),
    ]
}

fn items(app: App) -> Vec<Item> {
    use ItemType::*;
    match app {
        App::Prep => common_items(),
        App::Pipeline => {
            // TE10: -purple_dir is registered as required, though only a tumor needs it
            let purple_required = crate::compat::java("TE10");
            let mut v = vec![
                item("purple_dir", Path, purple_required, "Path to Purple pipeline files", None),
                item("cobalt_dir", Path, true, "Path to Cobalt pipeline files", None),
                item("reference_wgs_metrics", Path, false, "Path to reference BAM METRICS file", None),
                item("tumor_wgs_metrics", Path, false, "Path to tumor BAM METRICS file", None),
            ];
            v.extend(common_items());
            v
        }
        App::Standalone => {
            let mut v = common_items();
            v.extend([
                item("reference_duplicate_proportion", Str, false, "Proportion of reads that are marked duplicates in the reference sample BAM", None),
                item("reference_mean_read_depth", Str, false, "Mean read depth of the reference sample", None),
                item("reference_gc50_read_depth", Str, false, "GC 50 read depth of the reference sample", None),
                item("tumor_purity", Decimal, false, "Purity of the tumor sample", Some("1.0")),
                item("tumor_ploidy", Decimal, false, "Ploidy of the tumor", Some("2.0")),
                item("tumor_duplicate_proportion", Decimal, false, "Proportion of reads that are marked duplicates in the tumor sample BAM", Some("0.0")),
                item("tumor_mean_read_depth", Str, false, "Mean read depth of the tumor sample", None),
                item("tumor_gc50_read_depth", Str, false, "GC 50 read depth of the tumor sample", None),
            ]);
            v
        }
    }
}

/// Integer.parseInt
pub fn java_parse_int(s: &str) -> Option<i32> {
    let digits = s.strip_prefix('-').or_else(|| s.strip_prefix('+')).unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse::<i32>().ok()
}

/// Double.parseDouble for the forms a command line or a TSV holds
pub fn java_parse_double(s: &str) -> Option<f64> {
    let t = s.trim_matches(|c: char| c <= ' ');
    let t = t.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(t);
    match t {
        "NaN" | "+NaN" | "-NaN" => return Some(f64::NAN),
        "Infinity" | "+Infinity" => return Some(f64::INFINITY),
        "-Infinity" => return Some(f64::NEG_INFINITY),
        _ => {}
    }
    if t.is_empty() || t.contains(|c: char| c.is_ascii_alphabetic() && c != 'e' && c != 'E') {
        return None;
    }
    t.parse::<f64>().ok()
}

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
    /// checkAndParseCommandLine
    pub fn parse(app: App, argv: &[String]) -> Result<Parsed> {
        let items = items(app);
        if argv.len() == 1 && (argv[0] == "-version" || argv[0] == "-help") {
            println!("{APP_NAME} version {VERSION}");
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
            anyhow!("invalid command line")
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
                ItemType::Integer if java_parse_int(value).is_none() => {
                    crate::error!("invalid type for config: {}: {} = {}", item.name, item.ty.name(), value);
                    ok = false;
                }
                ItemType::Decimal if java_parse_double(value).is_none() => {
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
    pub fn integer(&self, name: &str) -> Result<i32> {
        let v = self.value(name).unwrap_or_default();
        java_parse_int(&v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
    pub fn decimal(&self, name: &str) -> Result<f64> {
        let v = self.value(name).ok_or_else(|| anyhow!("java.lang.NullPointerException ({name})"))?;
        java_parse_double(&v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
    fn optional_decimal(&self, name: &str) -> Result<Option<f64>> {
        if self.has_value(name) { Ok(Some(self.decimal(name)?)) } else { Ok(None) }
    }
}

/// ConfigBuilder.printItems
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

/// TealCommonParams
#[derive(Clone, Debug)]
pub struct CommonParams {
    pub reference_sample_id: Option<String>,
    pub tumor_sample_id: Option<String>,
    pub reference_bam: Option<String>,
    pub tumor_bam: Option<String>,
    pub ref_genome: Option<String>,
    /// with a trailing separator, as parseOutputDir gives it
    pub output_dir: Option<String>,
    pub threads: i32,
    pub v37: bool,
    pub sequencing_type: SequencingType,
}

impl CommonParams {
    pub fn from(cb: &ConfigBuilder) -> Result<CommonParams> {
        let output_dir = cb.value("output_dir").map(|d| if d.is_empty() || d.ends_with('/') { d } else { format!("{d}/") });
        // RefGenomeVersion.from(configBuilder): V37 unless given
        let v37 = if cb.has_value("ref_genome_version") {
            match cb.value("ref_genome_version").unwrap().as_str() {
                "V37" | "37" | "HG37" => true,
                "V38" | "38" | "HG38" => false,
                v => bail!("java.lang.IllegalArgumentException: Cannot resolve ref genome version: {v}"),
            }
        } else {
            true
        };
        let sequencing_type = if cb.has_value("sequencing_type") {
            let v = cb.value("sequencing_type").unwrap();
            SequencingType::parse(&v)
                .ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.sequencing.SequencingType.{v}"))?
        } else {
            SequencingType::Illumina
        };
        Ok(CommonParams {
            reference_sample_id: cb.value("reference"),
            tumor_sample_id: cb.value("tumor"),
            reference_bam: cb.value("reference_bam"),
            tumor_bam: cb.value("tumor_bam"),
            ref_genome: cb.value("ref_genome"),
            output_dir,
            threads: cb.integer("threads")?,
            v37,
            sequencing_type,
        })
    }

    pub fn reference_only(&self) -> bool {
        self.reference_sample_id.is_some() && self.tumor_sample_id.is_none()
    }
    pub fn tumor_only(&self) -> bool {
        self.reference_sample_id.is_none() && self.tumor_sample_id.is_some()
    }

    /// validate (a ParameterException in Java)
    pub fn validate(&self) -> Result<()> {
        if self.reference_sample_id.is_none() && self.reference_bam.is_some() {
            bail!("com.beust.jcommander.ParameterException: reference is required when reference_bam({}) is specified", self.reference_bam.as_ref().unwrap());
        }
        if self.reference_sample_id.is_some() && self.reference_bam.is_none() {
            bail!("com.beust.jcommander.ParameterException: reference_bam is required when reference({}) is specified", self.reference_sample_id.as_ref().unwrap());
        }
        if self.tumor_sample_id.is_none() && self.tumor_bam.is_some() {
            bail!("com.beust.jcommander.ParameterException: tumor is required when tumor_bam({}) is specified", self.tumor_bam.as_ref().unwrap());
        }
        if self.tumor_sample_id.is_some() && self.tumor_bam.is_none() {
            bail!("com.beust.jcommander.ParameterException: tumor_bam is required when tumor({}) is specified", self.tumor_sample_id.as_ref().unwrap());
        }
        Ok(())
    }

    fn out(&self) -> String {
        // Java builds "${outputDir}/..." (outputDir may be null)
        format!("{}/", self.output_dir.as_deref().unwrap_or("null"))
    }
    pub fn germline_telbam_path(&self) -> String {
        format!("{}{}.teal.telbam.bam", self.out(), self.reference_sample_id.as_deref().unwrap_or("null"))
    }
    pub fn tumor_telbam_path(&self) -> String {
        format!("{}{}.teal.telbam.bam", self.out(), self.tumor_sample_id.as_deref().unwrap_or("null"))
    }
    pub fn germline_tellength_path(&self) -> String {
        format!("{}{}.teal.tellength.tsv", self.out(), self.reference_sample_id.as_deref().unwrap_or("null"))
    }
    pub fn tumor_tellength_path(&self) -> String {
        format!("{}{}.teal.tellength.tsv", self.out(), self.tumor_sample_id.as_deref().unwrap_or("null"))
    }
    pub fn breakend_path(&self) -> String {
        format!("{}{}.teal.breakend.tsv.gz", self.out(), self.tumor_sample_id.as_deref().unwrap_or("null"))
    }
}

/// TealParams: the values TealApplication needs beyond the common ones
#[derive(Clone, Debug)]
pub struct TealParams {
    pub common: CommonParams,
    pub germline_duplicate_proportion: f64,
    pub germline_mean_read_depth: Option<f64>,
    pub germline_gc50_read_depth: Option<f64>,
    pub tumor_purity: f64,
    pub tumor_ploidy: f64,
    pub tumor_duplicate_proportion: f64,
    pub tumor_mean_read_depth: Option<f64>,
    pub tumor_gc50_read_depth: Option<f64>,
}

impl TealParams {
    /// TealParams(commonParams = ...): the defaults
    pub fn with_common(common: CommonParams) -> TealParams {
        TealParams {
            common,
            germline_duplicate_proportion: 0.0,
            germline_mean_read_depth: None,
            germline_gc50_read_depth: None,
            tumor_purity: 1.0,
            tumor_ploidy: 2.0,
            tumor_duplicate_proportion: 0.0,
            tumor_mean_read_depth: None,
            tumor_gc50_read_depth: None,
        }
    }

    /// TealParams(configBuilder)
    pub fn from(cb: &ConfigBuilder) -> Result<TealParams> {
        Ok(TealParams {
            common: CommonParams::from(cb)?,
            germline_duplicate_proportion: cb.optional_decimal("reference_duplicate_proportion")?.unwrap_or(0.0),
            germline_mean_read_depth: cb.optional_decimal("reference_mean_read_depth")?,
            germline_gc50_read_depth: cb.optional_decimal("reference_gc50_read_depth")?,
            tumor_purity: cb.decimal("tumor_purity")?,
            tumor_ploidy: cb.decimal("tumor_ploidy")?,
            tumor_duplicate_proportion: cb.optional_decimal("tumor_duplicate_proportion")?.unwrap_or(0.0),
            tumor_mean_read_depth: cb.optional_decimal("tumor_mean_read_depth")?,
            tumor_gc50_read_depth: cb.optional_decimal("tumor_gc50_read_depth")?,
        })
    }

    pub fn validate(&self) -> Result<()> {
        self.common.validate()?;
        if let Some(r) = &self.common.reference_sample_id {
            if self.germline_mean_read_depth.is_none() {
                bail!("com.beust.jcommander.ParameterException: reference_mean_read_depth is required when reference({r}) is specified");
            }
        }
        if let Some(t) = &self.common.tumor_sample_id {
            if self.tumor_mean_read_depth.is_none() {
                bail!("com.beust.jcommander.ParameterException: tumor_mean_read_depth is required when tumor({t}) is specified");
            }
        }
        Ok(())
    }
}

/// TealPipelineParams
#[derive(Clone, Debug)]
pub struct PipelineParams {
    pub purple: Option<String>,
    pub cobalt: Option<String>,
    pub reference_wgs_metrics: Option<String>,
    pub tumor_wgs_metrics: Option<String>,
    pub common: CommonParams,
}

impl PipelineParams {
    pub fn from(cb: &ConfigBuilder) -> Result<PipelineParams> {
        Ok(PipelineParams {
            purple: cb.value("purple_dir"),
            cobalt: cb.value("cobalt_dir"),
            reference_wgs_metrics: cb.value("reference_wgs_metrics"),
            tumor_wgs_metrics: cb.value("tumor_wgs_metrics"),
            common: CommonParams::from(cb)?,
        })
    }

    pub fn validate(&self) -> Result<()> {
        self.common.validate()?;
        if let Some(r) = &self.common.reference_sample_id {
            if self.reference_wgs_metrics.is_none() {
                bail!("com.beust.jcommander.ParameterException: reference_wgs_metrics is required when reference({r}) is specified");
            }
        }
        if let Some(t) = &self.common.tumor_sample_id {
            if self.purple.is_none() {
                bail!("com.beust.jcommander.ParameterException: purple_dir is required when tumor({t}) is specified");
            }
            if self.tumor_wgs_metrics.is_none() {
                bail!("com.beust.jcommander.ParameterException: tumor_wgs_metrics is required when tumor({t}) is specified");
            }
        }
        Ok(())
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
        let Parsed::Run(cb) = ConfigBuilder::parse(App::Prep, &strings(&["-tumor", "T", "-threads", "2", "-output_dir", "out", "-log_debug"])).unwrap() else {
            panic!()
        };
        let c = CommonParams::from(&cb).unwrap();
        assert_eq!(c.tumor_sample_id.as_deref(), Some("T"));
        assert_eq!(c.output_dir.as_deref(), Some("out/"));
        assert_eq!(c.threads, 2);
        assert!(c.v37);
        assert!(c.tumor_only() && !c.reference_only());
        assert!(c.validate().is_err());
        assert_eq!(c.tumor_telbam_path(), "out//T.teal.telbam.bam");
        assert!(ConfigBuilder::parse(App::Prep, &strings(&["-unknown", "x"])).is_err());
        assert!(ConfigBuilder::parse(App::Prep, &strings(&["-threads", "x"])).is_err());
        assert!(ConfigBuilder::parse(App::Pipeline, &strings(&["-tumor", "T"])).is_err());
    }

    #[test]
    fn standalone_values() {
        let Parsed::Run(cb) = ConfigBuilder::parse(App::Standalone, &strings(&["-reference", "R", "-reference_mean_read_depth", "30.5", "-ref_genome_version", "38"])).unwrap() else {
            panic!()
        };
        let p = TealParams::from(&cb).unwrap();
        assert_eq!(p.germline_mean_read_depth, Some(30.5));
        assert_eq!(p.tumor_purity, 1.0);
        assert!(!p.common.v37);
        assert_eq!(p.germline_duplicate_proportion, 0.0);
    }
}
