//! CobaltConfig and the command line, parsed as the Java ConfigBuilder does: `-name value` and `-flag`,
//! required items, path checks, integer and decimal checks, `-version` and `-help` as the sole argument.

use anyhow::{bail, Result};
use std::collections::HashMap;

use crate::constants::*;
use hmf_common::regions::{ChrBaseRegion, SpecificRegions};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ItemType {
    Str,
    Path,
    Integer,
    Decimal,
    Flag,
}

impl ItemType {
    /// ConfigItemType names as ConfigBuilder prints them
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

/// The items CobaltConfig.registerConfig registers, in that order.
fn items() -> Vec<Item> {
    use ItemType::*;
    let item = |name, ty, required, desc, default| Item { name, ty, required, desc, default };
    vec![
        item("tumor", Str, false, "Tumor ID", None),
        item("tumor_bam", Path, false, "Tumor BAM file", None),
        item("reference", Str, false, "Reference ID", None),
        item("reference_bam", Path, false, "Reference BAM file", None),
        item("gc_profile", Path, true, "Path to GC profile", None),
        item("ref_genome", Path, false, "Path to reference genome fasta files, required when using CRAM files", None),
        item("ref_genome_version", Str, true, "Ref genome version, 37 or 38", None),
        item("gc_ratio_min", Decimal, false, "Restrict GC ratios to above minimum", Some("0.24")),
        item("gc_ratio_max", Decimal, false, "Restrict GC ratios to below maximum", Some("0.68")),
        item("tumor_only_diploid_bed", Path, false, "Diploid regions for tumor-only mode", None),
        item("target_region_norm_file", Path, false, "Targeted regions normalisation file", None),
        item("min_quality", Integer, false, "Min map quality", Some("10")),
        item("pcf_gamma", Integer, false, "Gamma value for copy number PCF", Some("100")),
        item("include_duplicates", Flag, false, "Include duplicate reads in depth counts", None),
        item("skip_pcf_calc", Flag, false, "Skip final PCF output", None),
        item("specific_chr", Str, false, "Restrict to chromosome(s) separated by ';'", None),
        item(
            "specific_regions",
            Str,
            false,
            "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd",
            None,
        ),
        item("output_dir", Str, false, "Output directory", None),
        item("threads", Integer, false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        item("bam_validation", Str, false, "BAM validation: STRICT (default), LENIENT or SILENT", Some("SILENT")),
        item("log_debug", Flag, false, "Log at DEBUG level", None),
        item("log_level", Str, false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", None),
    ]
}

/// Java's Integer.parseInt: an optional sign and decimal digits, within the int range
pub fn java_parse_int(s: &str) -> Option<i32> {
    let digits = s.strip_prefix('-').or_else(|| s.strip_prefix('+')).unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse::<i32>().ok()
}

/// Java's Double.parseDouble for the forms a command line or a TSV holds (leading and trailing white
/// space is ignored, as Java does; a trailing d/D/f/F is accepted)
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
            println!("Cobalt version {}", crate::VERSION);
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
    pub fn integer(&self, name: &str) -> Result<i32> {
        let v = self.value(name).unwrap_or_default();
        java_parse_int(&v).ok_or_else(|| anyhow::anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
    pub fn decimal(&self, name: &str) -> Result<f64> {
        let v = self.value(name).unwrap_or_default();
        java_parse_double(&v).ok_or_else(|| anyhow::anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
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

/// The options of one run (CobaltConfig).
#[derive(Clone, Debug)]
pub struct CobaltConfig {
    pub tumor_id: Option<String>,
    pub tumor_bam: Option<String>,
    pub reference_id: Option<String>,
    pub reference_bam: Option<String>,
    pub ref_genome_path: Option<String>,
    pub v38: bool,
    pub gc_profile_path: String,
    /// -threads as given (Java builds its thread pool with this count)
    pub threads: i32,
    /// with the trailing separator, as parseOutputDir gives it; None when -output_dir is not given
    pub output_dir: Option<String>,
    pub min_mapping_quality: i32,
    pub pcf_gamma: i32,
    pub include_duplicates: bool,
    pub skip_pcf_calc: bool,
    pub tumor_only_diploid_bed: Option<String>,
    pub target_region_norm_file: Option<String>,
    pub gc_ratio_min: f64,
    pub gc_ratio_max: f64,
    /// None when SpecificRegions.from fails (Java goes on with a null and fails when it reads a BAM)
    pub specific_regions: Option<SpecificRegions>,
}

impl CobaltConfig {
    /// The CobaltConfig constructor. Its errors are those Java throws before validate() (an unknown genome
    /// version or validation stringency, a number that does not parse).
    pub fn from(cb: &ConfigBuilder) -> Result<CobaltConfig> {
        // RefGenomeVersion.from(configBuilder): V37 unless given
        let v38 = match cb.value("ref_genome_version") {
            None => false,
            Some(version) => match version.as_str() {
                "V37" | "37" | "HG37" => false,
                "V38" | "38" | "HG38" => true,
                _ => bail!("java.lang.IllegalArgumentException: Cannot resolve ref genome version: {version}"),
            },
        };
        let gc_ratio_min = cb.decimal("gc_ratio_min")?;
        let gc_ratio_max = cb.decimal("gc_ratio_max")?;
        let min_mapping_quality = cb.integer("min_quality")?;
        let pcf_gamma = cb.integer("pcf_gamma")?;
        let bam_validation = cb.value("bam_validation").unwrap();
        if !matches!(bam_validation.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            bail!("java.lang.IllegalArgumentException: No enum constant htsjdk.samtools.ValidationStringency.{bam_validation}");
        }
        let output_dir = cb.value("output_dir").map(|d| if !d.is_empty() && !d.ends_with('/') { format!("{d}/") } else { d });
        let threads = cb.integer("threads")?;
        let specific_regions = parse_specific_regions(cb.value("specific_regions").as_deref(), cb.value("specific_chr").as_deref());
        Ok(CobaltConfig {
            tumor_id: cb.value("tumor"),
            tumor_bam: cb.value("tumor_bam"),
            reference_id: cb.value("reference"),
            reference_bam: cb.value("reference_bam"),
            ref_genome_path: cb.value("ref_genome"),
            v38,
            gc_profile_path: cb.value("gc_profile").unwrap(),
            threads,
            output_dir,
            min_mapping_quality,
            pcf_gamma,
            include_duplicates: cb.has_flag("include_duplicates"),
            skip_pcf_calc: cb.has_flag("skip_pcf_calc"),
            tumor_only_diploid_bed: cb.value("tumor_only_diploid_bed"),
            target_region_norm_file: cb.value("target_region_norm_file"),
            gc_ratio_min,
            gc_ratio_max,
            specific_regions,
        })
    }

    /// CobaltConfig.validate, which also creates the output directory. The message is that of the
    /// exception Java logs as "config loading failed: {}".
    pub fn validate(&self) -> Result<(), String> {
        if self.reference_id.is_none() {
            if self.reference_bam.is_some() {
                return Err("java.lang.Exception: reference_bam option not allowed in tumor only mode".into());
            }
        } else if self.tumor_only_diploid_bed.is_some() {
            return Err("java.lang.Exception: tumor_only_diploid_bed option is only allowed in tumor only mode".into());
        }
        if self.tumor_id.is_none() && self.tumor_bam.is_some() {
            return Err("java.lang.Exception: tumor_bam option not allowed in germline only mode".into());
        }
        let Some(dir) = &self.output_dir else {
            return Err("java.lang.NullPointerException".into());
        };
        if !std::path::Path::new(dir).exists() && std::fs::create_dir_all(dir).is_err() {
            return Err(format!("java.lang.Exception: failed to create output directory({dir})"));
        }
        if self.gc_profile_path.ends_with(".gz") {
            return Err(format!("java.lang.Exception: invalid GC-profile file({}), must be uncompressed", self.gc_profile_path));
        }
        Ok(())
    }

    pub fn out_dir(&self) -> &str {
        self.output_dir.as_deref().unwrap_or("")
    }

    /// hasReferenceId: a reference id that is not empty
    pub fn has_reference_id(&self) -> bool {
        self.reference_id.as_ref().is_some_and(|r| !r.is_empty())
    }

    /// CobaltRatioFile.generateFilename(OutputDir, TumorId != null ? TumorId : ReferenceId)
    pub fn ratios_file_name(&self) -> String {
        let sample = self.tumor_id.as_deref().or(self.reference_id.as_deref()).unwrap_or("null");
        format!("{}{sample}.cobalt.ratio.tsv.gz", self.out_dir())
    }

    /// PCFFile.generateRatioFilename
    pub fn pcf_file_name(&self, sample: &str) -> String {
        format!("{}{sample}.cobalt.ratio.pcf", self.out_dir())
    }

    /// CobaltMedianRatioFile.generateFilename: basePath + File.separator + sample + extension
    pub fn median_ratios_file_name(&self, sample: &str) -> String {
        format!("{}/{sample}.cobalt.ratio.median.tsv", self.out_dir())
    }

    /// CobaltGcMedianFile.generateFilename: basePath + File.separator + sample + extension
    pub fn gc_median_file_name(&self, sample: &str) -> String {
        format!("{}/{sample}.cobalt.gc.median.tsv", self.out_dir())
    }

    /// The GC bucket limits of the GC normalisation. CO1: Java sets GC_RATIO_MIN/MAX from -gc_ratio_min and
    /// -gc_ratio_max but never uses them: its bucket limits come from the defaults when the class loads.
    pub fn gc_bucket_limits(&self) -> (i32, i32) {
        if crate::compat::java("CO1") {
            (calc_gc_bucket(DEFAULT_GC_RATIO_MIN), calc_gc_bucket(DEFAULT_GC_RATIO_MAX))
        } else {
            (calc_gc_bucket(self.gc_ratio_min), calc_gc_bucket(self.gc_ratio_max))
        }
    }

    /// The GC range of the windows whose read depths give the mean and median (ReadDepthStatisticsNormaliser
    /// uses the defaults; CO1 corrected: the options).
    pub fn gc_ratio_limits(&self) -> (f64, f64) {
        if crate::compat::java("CO1") {
            (DEFAULT_GC_RATIO_MIN, DEFAULT_GC_RATIO_MAX)
        } else {
            (self.gc_ratio_min, self.gc_ratio_max)
        }
    }
}

/// SpecificRegions.from: `-specific_regions` (chr:start-end or chr:start:end, separated by ';') wins over
/// `-specific_chr`; the regions are sorted (ChrBaseRegion.compareTo) and give the chromosome list. None
/// (logged as Java does) when a region does not parse.
pub fn parse_specific_regions(regions: Option<&str>, chromosomes: Option<&str>) -> Option<SpecificRegions> {
    let mut s = SpecificRegions::default();
    if let Some(r) = regions.filter(|r| !r.is_empty()) {
        for item in r.split(';') {
            let parts: Vec<&str> = item.split(':').collect();
            let region = if parts.len() == 3 {
                match (java_parse_int(parts[1]), java_parse_int(parts[2])) {
                    (Some(a), Some(b)) => Some(ChrBaseRegion::new(parts[0], a, b)),
                    _ => None,
                }
            } else if parts.len() == 2 && parts[1].contains('-') {
                let (a, b) = parts[1].split_once('-').unwrap();
                match (java_parse_int(a), java_parse_int(b)) {
                    (Some(a), Some(b)) => Some(ChrBaseRegion::new(parts[0], a, b)),
                    _ => None,
                }
            } else {
                None
            };
            // ChrBaseRegion.isValid(false): positive start, end not before start
            match region {
                Some(reg) if reg.start > 0 && reg.end >= reg.start => s.regions.push(reg),
                _ => {
                    crate::error!("invalid specific config: regions({r}) chromosomes({})", chromosomes.unwrap_or("null"));
                    return None;
                }
            }
        }
        s.regions.sort_by(|a, b| {
            hmf_common::regions::compare_chromosomes(&a.chromosome, &b.chromosome).then(a.start.cmp(&b.start)).then(a.end.cmp(&b.end))
        });
        for r in &s.regions {
            if !s.chromosomes.contains(&r.chromosome) {
                s.chromosomes.push(r.chromosome.clone());
            }
        }
    } else if let Some(c) = chromosomes.filter(|c| !c.is_empty()) {
        s.chromosomes = c.split(';').map(|x| x.to_string()).collect();
    }
    Some(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_like_config_builder() {
        let base = ["-gc_profile", "/", "-ref_genome_version", "38"];
        let args = |extra: &[&str]| strings(&[&base[..], extra].concat());
        let Parsed::Run(cb) = ConfigBuilder::parse(&args(&["-threads", "2", "-include_duplicates", "-gc_ratio_min", "-0.5"])).unwrap() else {
            panic!()
        };
        assert_eq!(cb.integer("threads").unwrap(), 2);
        assert!(cb.has_flag("include_duplicates") && !cb.has_flag("skip_pcf_calc"));
        assert_eq!(cb.decimal("gc_ratio_min").unwrap(), -0.5);
        assert_eq!(cb.decimal("gc_ratio_max").unwrap(), 0.68);
        assert_eq!(cb.integer("min_quality").unwrap(), 10);
        assert_eq!(cb.integer("pcf_gamma").unwrap(), 100);
        assert!(ConfigBuilder::parse(&args(&["-bogus"])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-gc_profile", "/"])).is_err()); // ref_genome_version is required
        assert!(ConfigBuilder::parse(&strings(&["-gc_profile", "/no/such/file", "-ref_genome_version", "38"])).is_err());
        assert!(ConfigBuilder::parse(&args(&["-threads", "x"])).is_err());
        assert!(ConfigBuilder::parse(&args(&["-tumor", "-x"])).is_err());
    }

    #[test]
    fn config_values() {
        let Parsed::Run(cb) =
            ConfigBuilder::parse(&strings(&["-gc_profile", "/", "-ref_genome_version", "V38", "-output_dir", "out", "-tumor", "T"])).unwrap()
        else {
            panic!()
        };
        let c = CobaltConfig::from(&cb).unwrap();
        assert!(c.v38 && c.output_dir.as_deref() == Some("out/") && c.threads == 1 && c.min_mapping_quality == 10);
        assert_eq!(c.ratios_file_name(), "out/T.cobalt.ratio.tsv.gz");
        assert_eq!(c.pcf_file_name("T"), "out/T.cobalt.ratio.pcf");
        assert_eq!(c.gc_median_file_name("T"), "out//T.cobalt.gc.median.tsv");
        assert_eq!(c.gc_bucket_limits(), (24, 68));
        let Parsed::Run(cb) = ConfigBuilder::parse(&strings(&["-gc_profile", "/", "-ref_genome_version", "39"])).unwrap() else { panic!() };
        assert!(CobaltConfig::from(&cb).is_err());
    }

    #[test]
    fn specific_regions() {
        let s = parse_specific_regions(Some("chr2:10-20;chr1:5:8"), Some("chr3")).unwrap();
        assert_eq!(s.chromosomes, vec!["chr1", "chr2"]);
        let c = parse_specific_regions(None, Some("chr3;chr4")).unwrap();
        assert!(c.include_chromosome("chr4") && !c.include_chromosome("chr1"));
        assert!(parse_specific_regions(Some("chr1:20-10"), None).is_none());
        assert!(parse_specific_regions(Some("chr1"), None).is_none());
    }

    #[test]
    fn java_number_parsing() {
        assert_eq!(java_parse_int("+5"), Some(5));
        assert_eq!(java_parse_int("5.0"), None);
        assert_eq!(java_parse_double(" 0.5 "), Some(0.5));
        assert_eq!(java_parse_double("2d"), Some(2.0));
        assert!(java_parse_double("abc").is_none());
    }
}
