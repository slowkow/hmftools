//! The step's options: Java DepthConfig, with the command line read the way ConfigBuilder reads it.

use super::vcf::parse_java_int;
use crate::common::compat;
use anyhow::{bail, Result};

#[derive(Clone, Copy, PartialEq, Eq)]
enum ItemType {
    Str,
    Path,
    Paths,
    Integer,
    Decimal,
    Flag,
}

impl ItemType {
    fn name(self) -> &'static str {
        match self {
            ItemType::Str => "STRING",
            ItemType::Path => "PATH",
            ItemType::Paths => "PATHS",
            ItemType::Integer => "INTEGER",
            ItemType::Decimal => "DECIMAL",
            ItemType::Flag => "FLAG",
        }
    }
}

struct Item {
    ty: ItemType,
    name: &'static str,
    required: bool,
    description: &'static str,
    default: Option<&'static str>,
    value: Option<String>,
}

/// The options in the order Java registers them (DepthConfig.registerConfig).
fn registered_items() -> Vec<Item> {
    let item = |ty, name, required, description, default| Item { ty, name, required, description, default, value: None };
    // D10: Java registers -sample and -bam_file as optional and then fails with a NullPointerException
    // when one is absent, in looking up the tumor and reference options of the other steps. Corrected:
    // both are required, and reported like any other missing required option.
    let samples_required = !compat::java("D10");
    vec![
        item(ItemType::Path, "input_vcf", false, "Input VCF", None),
        item(ItemType::Str, "sample", samples_required, "List of samples separated by ','", None),
        item(ItemType::Paths, "bam_file", samples_required, "BAM file paths separated by ','", None),
        item(ItemType::Str, "vcf_tag_prefix", false, "VCF tag prefix for testing & comparison", None),
        item(ItemType::Str, "ref_genome_version", true, "Ref genome version, 37 or 38", None),
        item(ItemType::Path, "ref_genome", true, "Path to reference genome fasta files", None),
        item(ItemType::Path, "unmap_regions", false, "Genome regions to unmap reads", None),
        item(ItemType::Decimal, "vaf_cap", false, "Ref support depth limit as function of variant fragments", Some("0.001")),
        item(ItemType::Integer, "proximity_distance", false, "Proximity distance to group variants", Some("2000")),
        item(ItemType::Str, "bam_validation", false, "BAM validation: STRICT (default), LENIENT or SILENT", Some("SILENT")),
        item(ItemType::Str, "specific_chr", false, "Restrict to chromosome(s) separated by ';'", None),
        item(ItemType::Str, "specific_regions", false, "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd", None),
        item(ItemType::Integer, "threads", false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        item(ItemType::Str, "output_dir", false, "Output directory", None),
        item(ItemType::Str, "output_id", false, "Output file suffix", None),
        item(ItemType::Flag, "log_debug", false, "Log at DEBUG level", Some("false")),
        item(ItemType::Str, "log_level", false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", None),
        item(ItemType::Decimal, "perf_log_time", false, "Performance log time threshold (seconds)", Some("0.0")),
    ]
}

/// Double.parseDouble, near enough: surrounding whitespace and a trailing d/f type suffix are accepted
fn parse_java_double(s: &str) -> Option<f64> {
    let t = s.trim();
    let t = t.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(t);
    if t.is_empty() || t.eq_ignore_ascii_case("inf") || t.eq_ignore_ascii_case("+inf") || t.eq_ignore_ascii_case("-inf") || t.eq_ignore_ascii_case("nan") {
        return None;
    }
    match t {
        "Infinity" | "+Infinity" => Some(f64::INFINITY),
        "-Infinity" => Some(f64::NEG_INFINITY),
        "NaN" => Some(f64::NAN),
        _ if t.to_ascii_lowercase().contains("infinity") => None,
        _ => t.parse().ok(),
    }
}

/// String.split(delim) without a limit: trailing empty strings are dropped
fn java_split(s: &str, delim: char) -> Vec<&str> {
    let mut parts: Vec<&str> = s.split(delim).collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

struct Items(Vec<Item>);

impl Items {
    fn get(&self, name: &str) -> &Item {
        self.0.iter().find(|x| x.name == name).unwrap()
    }
    fn has_value(&self, name: &str) -> bool {
        self.get(name).value.is_some()
    }
    /// the given value, or else the default
    fn value(&self, name: &str) -> Option<String> {
        let item = self.get(name);
        item.value.clone().or(item.default.map(|x| x.to_string()))
    }

    /// ConfigBuilder.printItems
    fn print(&self, as_log: bool) {
        let mut output: Vec<String> = self
            .0
            .iter()
            .map(|item| {
                let mut s = format!(" -{}:", item.name);
                if item.ty != ItemType::Str {
                    s.push_str(&format!(" type({})", item.ty.name()));
                }
                if let (Some(default), true) = (item.default, item.ty != ItemType::Flag) {
                    s.push_str(&format!(" default({default})"));
                }
                if item.required {
                    s.push_str(" REQUIRED");
                }
                s.push_str(&format!(" desc({})", item.description));
                s
            })
            .collect();
        output.sort();
        output.insert(0, "registered config items:".to_string());
        for line in output {
            if as_log {
                info!("{line}");
            } else {
                println!("{line}");
            }
        }
    }

    /// ConfigBuilder.parseCommandLine; false when Java would stop with exit status 1
    fn parse(&mut self, args: &[String]) -> bool {
        if args.len() == 1 && (args[0] == "-help" || args[0] == "-version") {
            println!("Esvee version {}", crate::VERSION);
            if args[0] == "-help" {
                self.print(false);
            }
            std::process::exit(0);
        }

        let mut incorrect_argument = false;
        let mut matched: Option<usize> = None;
        for argument in args {
            if let Some(index) = matched {
                let item = &mut self.0[index];
                // a leading '-' is only accepted for a number, which may be negative
                if argument.starts_with('-') && item.ty != ItemType::Decimal && item.ty != ItemType::Integer {
                    error!("config item({}) has invalid argument: {}", item.name, argument);
                    incorrect_argument = true;
                    break;
                }
                item.value = Some(argument.clone());
                matched = None;
            } else if let Some(name) = argument.strip_prefix('-') {
                let Some(index) = self.0.iter().position(|x| x.name == name) else {
                    error!("unregistered config item: {}", argument);
                    incorrect_argument = true;
                    break;
                };
                if self.0[index].ty == ItemType::Flag {
                    self.0[index].value = Some("true".to_string());
                } else {
                    matched = Some(index);
                }
            } else {
                error!("expecting config name but encountered invalid argument: {}", argument);
                incorrect_argument = true;
                break;
            }
        }
        if incorrect_argument {
            self.print(true);
            return false;
        }

        // ConfigBuilder.isValid
        let mut valid = true;
        for item in &self.0 {
            let value = item.value.clone().or(item.default.map(|x| x.to_string()));
            if item.required && item.value.is_none() {
                error!("missing required config: {}: {}", item.name, item.description);
                valid = false;
            } else if item.ty == ItemType::Path || item.ty == ItemType::Paths {
                let Some(value) = &item.value else { continue };
                let paths: Vec<&str> = if item.ty == ItemType::Paths { value.split(',').collect() } else { vec![value.as_str()] };
                let exists = paths.iter().all(|p| p.contains('*') || p.contains('$') || std::path::Path::new(p).exists());
                if !exists {
                    error!("invalid path for config: {} = {}", item.name, value);
                    valid = false;
                }
            } else {
                let ok = match (item.ty, &value) {
                    (ItemType::Decimal, Some(v)) => parse_java_double(v).is_some(),
                    (ItemType::Integer, Some(v)) => parse_java_int(v).is_some(),
                    (ItemType::Decimal | ItemType::Integer, None) => false,
                    _ => true,
                };
                if !ok {
                    error!("invalid type for config: {}: {} = {}", item.name, item.ty.name(), value.as_deref().unwrap_or("null"));
                    valid = false;
                }
            }
        }
        valid
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ValidationStringency {
    Strict,
    Lenient,
    Silent,
}

pub struct DepthConfig {
    pub input_vcf: String,
    /// with a trailing separator; None when the option is not given
    pub output_dir: Option<String>,
    pub output_id: Option<String>,
    pub sample_ids: Vec<String>,
    pub bam_files: Vec<String>,
    pub ref_genome: String,
    pub unmap_regions_file: Option<String>,
    pub is_38: bool,
    pub vaf_cap: f64,
    pub proximity_distance: i32,
    pub vcf_tag_prefix: Option<String>,
    /// Accepted and checked as Java does, but htsjdk's record validation is not ported: a record that
    /// STRICT would stop on in Java is read here like any other.
    #[allow(dead_code)]
    pub bam_stringency: ValidationStringency,
    pub threads: i32,
    pub perf_log_time: f64,
    pub specific_regions: Vec<(String, i32, i32)>,
    /// -specific_chr, used when no -specific_regions are given (always empty with Java behaviour, D14)
    pub specific_chromosomes: Vec<String>,
}

/// The restriction given by -specific_regions and -specific_chr: (regions, chromosomes), or the text of
/// the region that cannot be read. As SpecificRegions.loadSpecificChromsomesOrRegions of hmf-common,
/// which the other steps use: regions when given, otherwise the chromosomes.
pub fn parse_specific(regions: Option<&str>, chromosomes: Option<&str>) -> std::result::Result<(Vec<(String, i32, i32)>, Vec<String>), String> {
    if let Some(value) = regions.filter(|v| !v.is_empty()) {
        let mut parsed = Vec::new();
        for region in value.split(';') {
            parsed.push(parse_region(region).ok_or_else(|| region.to_string())?);
        }
        return Ok((parsed, Vec::new()));
    }
    let chromosomes = chromosomes.map(|value| java_split(value, ';').into_iter().filter(|x| !x.is_empty()).map(|x| x.to_string()).collect()).unwrap_or_default();
    Ok((Vec::new(), chromosomes))
}

/// SpecificRegions.parseStandardFormat: chr:start-end or chr:start:end
fn parse_region(region: &str) -> Option<(String, i32, i32)> {
    let items = java_split(region, ':');
    if !(items.len() == 3 || (items.len() == 2 && items[1].contains('-'))) {
        return None;
    }
    let (start, end) = if items.len() == 3 {
        (parse_java_int(items[1])?, parse_java_int(items[2])?)
    } else {
        let (start, end) = items[1].split_once('-')?;
        (parse_java_int(start)?, parse_java_int(end)?)
    };
    // ChrBaseRegion.hasValidPositions
    if start > 0 && end >= start { Some((items[0].to_string(), start, end)) } else { None }
}

/// FileCommon.formOutputFile. D13: without -output_dir Java's string concatenation puts the text "null"
/// in front of the file name; corrected, the file is in the current directory under its own name.
pub fn form_output_file(output_dir: Option<&str>, sample_id: &str, file_type: &str, output_id: Option<&str>) -> String {
    let mut filename = output_dir.unwrap_or(if compat::java("D13") { "null" } else { "" }).to_string();
    filename.push_str(sample_id);
    filename.push_str(".esvee");
    if let Some(id) = output_id {
        filename.push('.');
        filename.push_str(id);
    }
    filename.push('.');
    filename.push_str(file_type);
    filename
}

impl DepthConfig {
    /// Reads the command line as ConfigBuilder.checkAndParseCommandLine and the DepthConfig constructor
    /// do, setting the log level on the way. Exits with status 1 where Java does.
    pub fn from_args(args: &[String]) -> Result<DepthConfig> {
        let mut items = Items(registered_items());
        if !items.parse(args) {
            std::process::exit(1);
        }

        // ConfigUtils.setLogLevel
        if items.has_value("log_debug") {
            super::set_log_level(super::LEVEL_DEBUG);
        } else if let Some(level) = items.get("log_level").value.as_deref() {
            let name = level.trim().to_uppercase();
            let Some(rank) = super::LEVEL_NAMES.iter().position(|x| *x == name) else {
                bail!("java.lang.IllegalArgumentException: Unknown level constant [{name}].");
            };
            super::set_log_level(rank as u8);
        }
        info!("Esvee version {}", crate::VERSION);

        // FileCommon.parseSampleIds and parseBamFiles look up the tumor and reference options of the other
        // steps when -sample or -bam_file is absent; this step does not register them, so Java fails
        let list = |name: &str, fallback: &str| -> Result<Vec<String>> {
            match items.get(name).value.as_deref() {
                Some(value) => Ok(java_split(value, ',').into_iter().map(|x| x.to_string()).collect()),
                None => {
                    warn!("invalid config item({}) requested", fallback);
                    bail!("java.lang.NullPointerException: -{name} is required (config item {fallback} is not registered)")
                }
            }
        };
        let sample_ids = list("sample", "tumor")?;
        let bam_files = list("bam_file", "tumor_bam")?;

        // FileWriterUtils.parseOutputDir
        let output_dir = items.get("output_dir").value.clone().map(|dir| if dir.is_empty() || dir.ends_with('/') { dir } else { format!("{dir}/") });
        let output_id = items.get("output_id").value.clone();

        let input_vcf = match items.get("input_vcf").value.clone() {
            Some(path) => path,
            None => {
                let Some(first) = sample_ids.first() else { bail!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0") };
                form_output_file(output_dir.as_deref(), first, "raw.vcf.gz", output_id.as_deref())
            }
        };

        let ref_genome = items.value("ref_genome").unwrap_or_default();
        let version = items.value("ref_genome_version").unwrap_or_default();
        let is_38 = match version.as_str() {
            "V37" | "37" | "HG37" => false,
            "V38" | "38" | "HG38" => true,
            _ => bail!("java.lang.IllegalArgumentException: Cannot resolve ref genome version: {version}"),
        };
        let vaf_cap = parse_java_double(&items.value("vaf_cap").unwrap()).unwrap();
        let proximity_distance = parse_java_int(&items.value("proximity_distance").unwrap()).unwrap();
        let stringency = items.value("bam_validation").unwrap();
        let bam_stringency = match stringency.as_str() {
            "STRICT" => ValidationStringency::Strict,
            "LENIENT" => ValidationStringency::Lenient,
            "SILENT" => ValidationStringency::Silent,
            _ => bail!("java.lang.IllegalArgumentException: No enum constant htsjdk.samtools.ValidationStringency.{stringency}"),
        };
        let perf_log_time = parse_java_double(&items.value("perf_log_time").unwrap()).unwrap();
        let threads = parse_java_int(&items.value("threads").unwrap()).unwrap();

        // SpecificRegions.loadSpecificRegions. D14: Java registers -specific_chr and never reads it, and
        // after an invalid -specific_regions reports the failure and carries on without any restriction.
        // Corrected: -specific_chr restricts the run, and an invalid region stops it.
        let java_specific = compat::java("D14");
        let specific_chr = if java_specific { None } else { items.get("specific_chr").value.as_deref() };
        let (specific_regions, specific_chromosomes) = match parse_specific(items.get("specific_regions").value.as_deref(), specific_chr) {
            Ok(parsed) => parsed,
            Err(_) if java_specific => {
                error!("failed to load specific regions");
                (Vec::new(), Vec::new())
            }
            Err(region) => {
                error!("invalid region in -specific_regions: {}", region);
                std::process::exit(1);
            }
        };

        Ok(DepthConfig {
            input_vcf,
            output_dir,
            output_id,
            sample_ids,
            bam_files,
            ref_genome,
            unmap_regions_file: items.get("unmap_regions").value.clone(),
            is_38,
            vaf_cap,
            proximity_distance,
            vcf_tag_prefix: items.get("vcf_tag_prefix").value.clone(),
            bam_stringency,
            threads,
            perf_log_time,
            specific_regions,
            specific_chromosomes,
        })
    }

    /// the configuration of the Java unit tests: DepthConfig(vcfCap, proximityDistance)
    #[cfg(test)]
    pub fn for_test(vaf_cap: f64, proximity_distance: i32) -> DepthConfig {
        DepthConfig {
            input_vcf: String::new(),
            output_dir: Some(String::new()),
            output_id: Some(String::new()),
            sample_ids: Vec::new(),
            bam_files: Vec::new(),
            ref_genome: String::new(),
            unmap_regions_file: None,
            is_38: false,
            vaf_cap,
            proximity_distance,
            vcf_tag_prefix: Some(String::new()),
            bam_stringency: ValidationStringency::Strict,
            threads: 0,
            perf_log_time: 0.0,
            specific_regions: Vec::new(),
            specific_chromosomes: Vec::new(),
        }
    }

    /// DepthConfig.getVcfTag
    pub fn vcf_tag(&self, tag: &str) -> String {
        match &self.vcf_tag_prefix {
            Some(prefix) => format!("{prefix}_{tag}"),
            None => tag.to_string(),
        }
    }

    /// RefGenomeVersion.versionedChromosome for a HumanChromosome name
    pub fn versioned_chromosome(&self, name: &str) -> String {
        if self.is_38 { format!("chr{name}") } else { name.to_string() }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions() {
        assert_eq!(parse_region("1:100-200"), Some(("1".to_string(), 100, 200)));
        assert_eq!(parse_region("chr2:5:9"), Some(("chr2".to_string(), 5, 9)));
        assert_eq!(parse_region("1:200-100"), None);
        assert_eq!(parse_region("1:0-100"), None);
        assert_eq!(parse_region("1"), None);
        assert_eq!(parse_region("1:100"), None);
    }

    #[test]
    fn output_names() {
        assert_eq!(form_output_file(Some("out/"), "S", "ref_depth.vcf.gz", None), "out/S.esvee.ref_depth.vcf.gz");
        assert_eq!(form_output_file(Some("out/"), "S", "raw.vcf.gz", Some("x")), "out/S.esvee.x.raw.vcf.gz");
        // D13: Java's name without an output directory
        assert_eq!(compat::with_java(true, || form_output_file(None, "S", "raw.vcf.gz", None)), "nullS.esvee.raw.vcf.gz");
        assert_eq!(compat::with_java(false, || form_output_file(None, "S", "raw.vcf.gz", Some("x"))), "S.esvee.x.raw.vcf.gz");
    }

    #[test]
    fn specific_chromosomes_and_regions() {
        // D14
        assert_eq!(parse_specific(None, Some("chr1;chr7")), Ok((Vec::new(), vec!["chr1".to_string(), "chr7".to_string()])));
        assert_eq!(parse_specific(None, None), Ok((Vec::new(), Vec::new())));
        // regions come first, as in the other steps
        assert_eq!(parse_specific(Some("chr2:5-9"), Some("chr1")), Ok((vec![("chr2".to_string(), 5, 9)], Vec::new())));
        assert_eq!(parse_specific(Some("chr2:5-9;chr3"), None), Err("chr3".to_string()));
        assert_eq!(parse_specific(Some("chr2:5-9;"), None), Err(String::new()));
    }

    #[test]
    fn samples_and_bams_are_required() {
        // D10; -ref_genome only has to be a path that exists
        let dir = std::env::temp_dir().to_string_lossy().to_string();
        let args = |extra: &[&str]| -> Vec<String> {
            ["-ref_genome", dir.as_str(), "-ref_genome_version", "38"].iter().chain(extra).map(|x| x.to_string()).collect()
        };
        let parse = |java: bool, extra: &[&str]| compat::with_java(java, || Items(registered_items()).parse(&args(extra)));
        assert!(parse(false, &["-sample", "A", "-bam_file", dir.as_str()]));
        assert!(!parse(false, &["-sample", "A"]));
        assert!(!parse(false, &["-bam_file", dir.as_str()]));
        // Java accepts the command line and fails later
        assert!(parse(true, &["-sample", "A"]));
    }

    #[test]
    fn doubles() {
        assert_eq!(parse_java_double("0.001"), Some(0.001));
        assert_eq!(parse_java_double(" 2d "), Some(2.0));
        assert_eq!(parse_java_double("x"), None);
        assert_eq!(parse_java_double("inf"), None);
    }
}
