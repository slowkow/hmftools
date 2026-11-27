//! AmberConfig and the command line, parsed as the Java ConfigBuilder does: `-name value` and `-flag`,
//! required items, path checks, integer and decimal checks, `-version` and `-help` as the sole argument.

use anyhow::{bail, Result};
use std::collections::HashMap;

use crate::constants::*;
use hmf_common::regions::{ChrBaseRegion, SpecificRegions};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ItemType {
    Str,
    Path,
    Paths,
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
            ItemType::Paths => "PATHS",
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
    desc: String,
    default: Option<String>,
}

/// The items AmberConfig.registerConfig registers, in that order.
fn items() -> Vec<Item> {
    use ItemType::*;
    let item = |name, ty, required, desc: &str, default: Option<&str>| Item {
        name,
        ty,
        required,
        desc: desc.to_string(),
        default: default.map(|s| s.to_string()),
    };
    vec![
        item("tumor", Str, false, "Tumor ID", None),
        item("tumor_bam", Path, false, "Tumor BAM file", None),
        item("reference", Str, false, "Reference ID", None),
        item("reference_bam", Paths, false, "Reference BAM file", None),
        item("loci", Path, true, "Path to BAF loci vcf file", None),
        item("target_regions_bed", Path, false, "Target regions BED file", None),
        item("ref_genome_version", Str, false, "Ref genome version, 37 or 38", Some("V37")),
        item("ref_genome", Path, false, "Path to reference genome fasta files", None),
        item("sequencing_type", Str, false, "Sequencing types: ILLUMINA, ULTIMA, SBX", Some("ILLUMINA")),
        item("tumor_only_min_support", Integer, false, "Min support in ref and alt in tumor-only mode", Some("2")),
        item("tumor_only_min_vaf", Decimal, false, "Min VAF in ref and alt in tumor-only mode", Some("0.025")),
        item(
            "tumor_min_depth",
            Str,
            false,
            &format!("Min tumor depth, default tumor/normal({DEFAULT_TUMOR_MIN_DEPTH}) tumor-only({DEFAULT_TUMOR_ONLY_MIN_DEPTH})"),
            None,
        ),
        item("min_base_quality", Integer, false, "Minimum quality for a base to be considered", Some("18")),
        item("min_map_quality", Integer, false, "Minimum mapping quality for an alignment to be used", Some("50")),
        item("position_gap", Integer, false, "Gap between site for BAM reading", Some("0")),
        item("min_depth_percent", Decimal, false, "Min percentage of median depth", Some("0.5")),
        item("max_depth_percent", Decimal, false, "Max percentage of median depth", Some("1.5")),
        item("min_het_af_percent", Decimal, false, "Max heterozygous AF%", Some("0.4")),
        item("max_het_af_percent", Decimal, false, "Max heterozygous AF%", Some("0.65")),
        item("perf_debug", Flag, false, "Detailed performance tracking and logging", None),
        item("skip_baf_segmentation", Flag, false, "Skip BAF segmentation", None),
        item("write_unfiltered_germline", Flag, false, "Write all (unfiltered) germline points", None),
        item("write_tumor_data", Flag, false, "Write tumor points", None),
        item("output_dir", Str, false, "Output directory", None),
        item("threads", Integer, false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        item("bam_validation", Str, false, "BAM validation: STRICT (default), LENIENT or SILENT", Some("SILENT")),
        item("log_debug", Flag, false, "Log at DEBUG level", None),
        item("log_level", Str, false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", None),
        item("specific_chr", Str, false, "Restrict to chromosome(s) separated by ';'", None),
        item(
            "specific_regions",
            Str,
            false,
            "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd",
            None,
        ),
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
            println!("Amber version {}", crate::VERSION);
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
                ItemType::Path | ItemType::Paths => {
                    let paths: Vec<&str> = if item.ty == ItemType::Paths { value.split(',').collect() } else { vec![value.as_str()] };
                    if paths.iter().any(|p| !std::path::Path::new(p).exists()) {
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
    pub fn has_flag(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }
    /// getValue: the value given, or the default
    pub fn value(&self, name: &str) -> Option<String> {
        if let Some(v) = self.values.get(name) {
            return Some(v.clone());
        }
        self.items.iter().find(|i| i.name == name).and_then(|i| i.default.clone())
    }
    pub fn integer(&self, name: &str) -> Result<i32> {
        let v = self.value(name).unwrap_or_default();
        java_parse_int(&v).ok_or_else(|| anyhow::anyhow!("invalid integer for {name}: {v}"))
    }
    pub fn decimal(&self, name: &str) -> Result<f64> {
        let v = self.value(name).unwrap_or_default();
        java_parse_double(&v).ok_or_else(|| anyhow::anyhow!("invalid decimal for {name}: {v}"))
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

/// The options of one run (AmberConfig).
#[derive(Clone, Debug)]
pub struct AmberConfig {
    pub tumor_id: Option<String>,
    pub tumor_bam: Option<String>,
    pub reference_ids: Vec<String>,
    pub reference_bams: Vec<String>,
    pub baf_loci_path: String,
    pub target_regions_bed: Option<String>,
    pub v38: bool,
    pub ref_genome_file: Option<String>,
    pub sequencing_type: SequencingType,
    pub tumor_only_min_support: i32,
    /// registered and parsed, but used nowhere in AMBER 4.4
    #[allow(dead_code)]
    pub tumor_only_min_vaf: f64,
    pub tumor_min_depth: i32,
    pub min_base_quality: i32,
    pub min_mapping_quality: i32,
    pub min_depth_percent: f64,
    pub max_depth_percent: f64,
    pub min_het_af_percent: f64,
    pub max_het_af_percent: f64,
    pub write_unfiltered_germline: bool,
    pub write_tumor_data: bool,
    pub position_gap: i32,
    /// with the trailing separator, as parseOutputDir gives it
    pub output_dir: String,
    pub threads: i32,
    pub skip_baf_segmentation: bool,
    pub perf_debug: bool,
    pub specific_regions: SpecificRegions,
}

impl AmberConfig {
    pub fn from(cb: &ConfigBuilder) -> Result<AmberConfig> {
        let tumor_id = cb.value("tumor");
        let tumor_bam = cb.value("tumor_bam");
        let mut reference_ids = Vec::new();
        let mut reference_bams = Vec::new();
        if cb.has_value("reference") && cb.has_value("reference_bam") {
            // String.split drops trailing empty strings
            let split = |s: String| -> Vec<String> {
                let mut v: Vec<String> = s.split(',').map(|x| x.to_string()).collect();
                while v.len() > 1 && v.last().is_some_and(|x| x.is_empty()) {
                    v.pop();
                }
                v
            };
            reference_ids = split(cb.value("reference").unwrap());
            reference_bams = split(cb.value("reference_bam").unwrap());
        }
        let version = cb.value("ref_genome_version").unwrap();
        let v38 = match version.as_str() {
            "V37" | "37" | "HG37" => false,
            "V38" | "38" | "HG38" => true,
            _ => bail!("Cannot resolve ref genome version: {version}"),
        };
        let seq = cb.value("sequencing_type").unwrap();
        let sequencing_type = match seq.as_str() {
            "ILLUMINA" => SequencingType::Illumina,
            "ULTIMA" => SequencingType::Ultima,
            "SBX" => SequencingType::Sbx,
            _ => bail!("No enum constant com.hartwig.hmftools.common.sequencing.SequencingType.{seq}"),
        };
        let tumor_min_depth = if cb.has_value("tumor_min_depth") {
            cb.integer("tumor_min_depth")?
        } else if reference_ids.is_empty() {
            DEFAULT_TUMOR_ONLY_MIN_DEPTH
        } else {
            DEFAULT_TUMOR_MIN_DEPTH
        };
        let target_regions_bed = cb.value("target_regions_bed");
        let (min_depth_percent, max_depth_percent) =
            if target_regions_bed.is_none() { (cb.decimal("min_depth_percent")?, cb.decimal("max_depth_percent")?) } else { (0.0, 0.0) };
        let output_dir = match cb.value("output_dir") {
            Some(d) if !d.is_empty() && !d.ends_with('/') => format!("{d}/"),
            Some(d) => d,
            None => bail!("missing output_dir"),
        };
        let bam_validation = cb.value("bam_validation").unwrap();
        if !matches!(bam_validation.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            bail!("No enum constant htsjdk.samtools.ValidationStringency.{bam_validation}");
        }
        let specific_regions = parse_specific_regions(cb.value("specific_regions").as_deref(), cb.value("specific_chr").as_deref())?;
        Ok(AmberConfig {
            tumor_id,
            tumor_bam,
            reference_ids,
            reference_bams,
            baf_loci_path: cb.value("loci").unwrap(),
            target_regions_bed,
            v38,
            ref_genome_file: cb.value("ref_genome"),
            sequencing_type,
            tumor_only_min_support: cb.integer("tumor_only_min_support")?,
            tumor_only_min_vaf: cb.decimal("tumor_only_min_vaf")?,
            tumor_min_depth,
            min_base_quality: cb.integer("min_base_quality")?,
            min_mapping_quality: cb.integer("min_map_quality")?,
            min_depth_percent,
            max_depth_percent,
            min_het_af_percent: cb.decimal("min_het_af_percent")?,
            max_het_af_percent: cb.decimal("max_het_af_percent")?,
            write_unfiltered_germline: cb.has_flag("write_unfiltered_germline"),
            write_tumor_data: cb.has_flag("write_tumor_data"),
            position_gap: cb.integer("position_gap")?,
            output_dir,
            threads: cb.integer("threads")?,
            skip_baf_segmentation: cb.has_flag("skip_baf_segmentation"),
            perf_debug: cb.has_flag("perf_debug"),
            specific_regions,
        })
    }

    pub fn is_ultima(&self) -> bool {
        self.sequencing_type == SequencingType::Ultima
    }
    pub fn primary_reference(&self) -> &str {
        &self.reference_ids[0]
    }
    pub fn is_tumor_only(&self) -> bool {
        self.reference_bams.is_empty() && self.tumor_bam.is_some()
    }
    pub fn is_germline_only(&self) -> bool {
        !self.reference_bams.is_empty() && self.tumor_bam.is_none()
    }
    /// the tumor id if there is one, otherwise the primary reference id
    pub fn sample_id(&self) -> &str {
        match &self.tumor_id {
            Some(t) => t,
            None => self.primary_reference(),
        }
    }
    /// The worker count of the BAM reading (Java: max(threads, 1))
    pub fn reader_threads(&self) -> usize {
        self.threads.max(1) as usize
    }

    /// AmberConfig.isValid, which also creates the output directory
    pub fn is_valid(&self) -> bool {
        if self.reference_ids.len() != self.reference_bams.len() {
            crate::error!("Each reference sample must have matching bam");
            return false;
        }
        if self.tumor_id.is_none() != self.tumor_bam.is_none() {
            crate::error!(
                "Unmatched: TumorId: {} and TumorBamPath: {}",
                self.tumor_id.as_deref().unwrap_or("null"),
                self.tumor_bam.as_deref().unwrap_or("null")
            );
            return false;
        }
        if !std::path::Path::new(&self.output_dir).exists() {
            let _ = std::fs::create_dir_all(&self.output_dir);
        }
        if !std::path::Path::new(&self.baf_loci_path).exists() {
            crate::error!("Unable to locate vcf file {}", self.baf_loci_path);
            return false;
        }
        if let Some(t) = &self.tumor_bam {
            if !std::path::Path::new(t).exists() {
                crate::error!("Unable to locate tumor bam file {t}");
                return false;
            }
        }
        if !self.is_tumor_only() {
            for r in &self.reference_bams {
                if !std::path::Path::new(r).exists() {
                    crate::error!("Unable to locate reference bam {r}");
                    return false;
                }
            }
        }
        true
    }
}

/// SpecificRegions.from: `-specific_regions` (chr:start-end or chr:start:end, separated by ';') wins over
/// `-specific_chr`; the regions are sorted (ChrBaseRegion.compareTo) and give the chromosome list.
pub fn parse_specific_regions(regions: Option<&str>, chromosomes: Option<&str>) -> Result<SpecificRegions> {
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
                    bail!("invalid specific region: {item}");
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
    Ok(s)
}

/// SpecificRegions.includePosition
pub fn include_position(s: &SpecificRegions, chromosome: &str, position: i32) -> bool {
    if s.regions.is_empty() {
        return s.include_chromosome(chromosome);
    }
    s.regions.iter().any(|r| r.chromosome == chromosome && r.contains_position(position))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_like_config_builder() {
        let Parsed::Run(cb) = ConfigBuilder::parse(&strings(&["-loci", "/", "-threads", "2", "-write_tumor_data", "-min_depth_percent", "-1"])).unwrap()
        else {
            panic!()
        };
        assert_eq!(cb.integer("threads").unwrap(), 2);
        assert!(cb.has_flag("write_tumor_data") && !cb.has_flag("perf_debug"));
        assert_eq!(cb.decimal("min_depth_percent").unwrap(), -1.0);
        assert_eq!(cb.value("ref_genome_version").as_deref(), Some("V37"));
        assert!(ConfigBuilder::parse(&strings(&["-loci", "/", "-bogus"])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-threads", "2"])).is_err()); // loci is required
        assert!(ConfigBuilder::parse(&strings(&["-loci", "/no/such/file"])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-loci", "/", "-threads", "x"])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-loci", "/", "-tumor", "-x"])).is_err());
    }

    #[test]
    fn specific_regions() {
        let s = parse_specific_regions(Some("chr2:10-20;chr1:5:8"), Some("chr3")).unwrap();
        assert_eq!(s.chromosomes, vec!["chr1", "chr2"]);
        assert!(include_position(&s, "chr1", 5) && include_position(&s, "chr1", 8) && !include_position(&s, "chr1", 9));
        assert!(!include_position(&s, "chr3", 10));
        let c = parse_specific_regions(None, Some("chr3;chr4")).unwrap();
        assert!(include_position(&c, "chr4", 1) && !include_position(&c, "chr1", 1));
        assert!(parse_specific_regions(Some("chr1:20-10"), None).is_err());
        assert!(parse_specific_regions(Some("chr1"), None).is_err());
    }

    #[test]
    fn java_number_parsing() {
        assert_eq!(java_parse_int("+5"), Some(5));
        assert_eq!(java_parse_int("-5"), Some(-5));
        assert_eq!(java_parse_int("5.0"), None);
        assert_eq!(java_parse_double(" 0.5 "), Some(0.5));
        assert_eq!(java_parse_double("1e-3"), Some(0.001));
        assert_eq!(java_parse_double("2d"), Some(2.0));
        assert!(java_parse_double("abc").is_none());
    }
}
