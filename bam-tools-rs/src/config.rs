//! The command line, parsed as the Java ConfigBuilder does (`-name value` and `-flag`, required items,
//! path checks, integer checks, `-version` and `-help` as the sole argument), and MetricsConfig.

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;

use hmf_common::java_hash::JavaHashMap;
use hmf_common::java_text::parse_int;
use hmf_common::regions::{BaseRegion, ChrBaseRegion, RefGenomeVersion, SpecificRegions};

use crate::compat;

pub const APP_NAME: &str = "BamTools";
pub const VERSION: &str = "1.6.2";

// the defaults of the items below (DEFAULT_CHR_PARTITION_SIZE 1,000,000, map quality 20, base quality 10,
// max coverage 250, partition read count check 1,000,000); a panel's max coverage is 10,000
#[cfg(test)]
pub const DEFAULT_CHR_PARTITION_SIZE: i32 = 1_000_000;
#[cfg(test)]
pub const DEFAULT_MAP_QUAL_THRESHOLD: i32 = 20;
#[cfg(test)]
pub const DEFAULT_BASE_QUAL_THRESHOLD: i32 = 10;
pub const DEFAULT_PANEL_MAX_COVERAGE: i32 = 10000;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ItemType {
    Str,
    Path,
    Integer,
    Flag,
}

impl ItemType {
    fn name(self) -> &'static str {
        match self {
            ItemType::Str => "STRING",
            ItemType::Path => "PATH",
            ItemType::Integer => "INTEGER",
            ItemType::Flag => "FLAG",
        }
    }
}

pub struct Item {
    pub name: &'static str,
    pub ty: ItemType,
    pub required: bool,
    pub desc: &'static str,
    pub default: Option<&'static str>,
}

/// The items MetricsConfig.addConfig registers, in that order
fn items() -> Vec<Item> {
    use ItemType::*;
    let it = |name, ty, required, desc, default| Item { name, ty, required, desc, default };
    vec![
        it("ref_genome", Path, true, "Path to reference genome fasta files", None),
        it("ref_genome_version", Str, false, "Ref genome version, 37 or 38", Some("V37")),
        it("specific_chr", Str, false, "Restrict to chromosome(s) separated by ';'", None),
        it("specific_regions", Str, false, "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd", None),
        it("bam_file", Path, true, "Input BAM file", None),
        it("sample", Str, false, "Sample ID", None),
        it(
            "regions_file",
            Path,
            false,
            "TSV or BED file with regions to analyse, expected columns Chromosome,PositionStart,PositionEnd or no headers",
            None,
        ),
        it("partition_size", Integer, false, "Partition size", Some("1000000")),
        it("map_qual_threshold", Integer, false, "Map quality threshold", Some("20")),
        it("base_qual_threshold", Integer, false, "Base quality threshold", Some("10")),
        it("max_coverage", Integer, false, "Max coverage", Some("250")),
        it("partition_read_count_check", Integer, false, "Partition read count log", Some("1000000")),
        it("only_target", Flag, false, "Only capture metrics within the specific regions file", None),
        it("no_duplicates", Flag, false, "Process BAM with duplicates dropped", None),
        it(
            "off_target_frag_overlap_threshold",
            Integer,
            false,
            "Write regions of high off-target fragment overlap if pile-up above threshold (0=disabled)",
            Some("0"),
        ),
        it("ensembl_data_dir", Path, false, "Ensembl data file directory", None),
        it("driver_gene_panel", Path, false, "Path to driver gene panel", None),
        it("exclude_zero_coverage", Flag, false, "Exclude bases with zero coverage", None),
        it("write_off_target", Flag, false, "Write off-target data", None),
        it("log_read_ids", Str, false, "Log specific read IDs, separated by ';'", None),
        it("perf_debug", Flag, false, "Detailed performance tracking and logging", None),
        it("output_dir", Str, false, "Output directory", None),
        it("output_id", Str, false, "Output file suffix", None),
        it("log_debug", Flag, false, "Log at DEBUG level", None),
        it("log_level", Str, false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", None),
        it("threads", Integer, false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        // port options (not in Java)
        it("file_buffer_kb", Integer, false, "Rust port: BAM file read size in KB (0 = htslib default)", Some("1024")),
    ]
}

/// The parsed command line
pub struct ConfigBuilder {
    items: Vec<Item>,
    values: HashMap<&'static str, String>,
}

pub enum Parsed {
    Run(ConfigBuilder),
    /// -version or -help was printed; exit 0
    Exit,
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
            if let Some(d) = item.default {
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

impl ConfigBuilder {
    /// ConfigBuilder.checkAndParseCommandLine: an error (logged) for anything Java rejects
    pub fn parse(argv: &[String]) -> Result<Parsed> {
        let items = items();
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
                if arg.starts_with('-') && items[i].ty != ItemType::Integer {
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
            match item.ty {
                ItemType::Path => {
                    if let Some(v) = value {
                        if !v.contains('*') && !v.contains('$') && !std::path::Path::new(v).exists() {
                            crate::error!("invalid path for config: {} = {}", item.name, v);
                            ok = false;
                        }
                    }
                }
                ItemType::Integer => {
                    let v = value.map(|s| s.as_str()).or(item.default);
                    if v.and_then(parse_int).is_none() {
                        crate::error!("invalid type for config: {}: {} = {}", item.name, item.ty.name(), v.unwrap_or("null"));
                        ok = false;
                    }
                }
                _ => {}
            }
        }
        if !ok {
            bail!("invalid command line");
        }
        Ok(Parsed::Run(ConfigBuilder { items, values }))
    }

    /// hasValue: given on the command line (a default does not count)
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
        self.items.iter().find(|i| i.name == name).and_then(|i| i.default.map(|s| s.to_string()))
    }
    pub fn integer(&self, name: &str) -> i32 {
        self.value(name).as_deref().and_then(parse_int).expect("checked by parse")
    }
}

/// checkAddDirSeparator
pub fn check_add_dir_separator(path: &str) -> String {
    if path.is_empty() || path.ends_with('/') { path.to_string() } else { format!("{path}/") }
}

/// FileWriterUtils.pathFromFile: the file's directory with a separator, or "./"
fn path_from_file(path: &str) -> String {
    if path.contains('/') {
        match std::path::Path::new(path).parent() {
            Some(p) => format!("{}/", p.to_string_lossy()),
            None => "null/".to_string(),
        }
    } else {
        "./".to_string()
    }
}

/// BamUtils.deriveRefGenomeVersion: 38 when the first sequence of the BAM header starts with "chr"
fn derive_ref_genome_version(bam_file: &str) -> Result<RefGenomeVersion> {
    use rust_htslib::bam::Read as _;
    let reader = rust_htslib::bam::Reader::from_path(bam_file).map_err(|e| anyhow!("opening {bam_file}: {e}"))?;
    let h = reader.header();
    if h.target_count() == 0 {
        bail!("no sequences in the header of {bam_file}");
    }
    Ok(if h.tid2name(0).starts_with(b"chr") { RefGenomeVersion::V38 } else { RefGenomeVersion::V37 })
}

/// MetricsConfig
pub struct MetricsConfig {
    pub sample_id: String,
    pub bam_file: String,
    pub ref_genome_file: String,
    pub ref_gen_version: RefGenomeVersion,
    /// the version the Ensembl data are read for (BT10)
    pub gene_ref_gen_version: RefGenomeVersion,
    pub map_quality_threshold: i32,
    pub base_quality_threshold: i32,
    pub max_coverage: i32,
    pub partition_size: i32,
    /// genome_unmappable_regions.<37|38>.bed, chromosomes as HumanChromosome names them ("1" for "chr1")
    pub unmappable_regions: Vec<ChrBaseRegion>,
    /// -regions_file, by chromosome, as Java's HashMap iterates
    pub target_regions: JavaHashMap<String, Vec<BaseRegion>>,
    pub only_target_regions: bool,
    pub ensembl_data_dir: Option<String>,
    pub driver_gene_panel: Option<String>,
    pub exclude_zero_coverage: bool,
    pub write_off_target: bool,
    pub no_duplicates: bool,
    pub high_fragment_overlap_threshold: i32,
    pub partition_read_count_check: i32,
    pub output_dir: String,
    pub output_id: Option<String>,
    pub threads: i32,
    pub specific_chr_regions: SpecificRegions,
    pub log_read_ids: Vec<String>,
    pub perf_debug: bool,
    pub file_buffer_kb: i32,
}

impl MetricsConfig {
    /// new MetricsConfig(configBuilder), then isValid()
    pub fn from(cb: &ConfigBuilder) -> Result<MetricsConfig> {
        let sample_id = cb.value("sample").unwrap_or_else(|| "null".to_string());
        let bam_file = cb.value("bam_file").ok_or_else(|| anyhow!("missing config: bam"))?;
        let ref_genome_file = cb.value("ref_genome").ok_or_else(|| anyhow!("missing config: refGenome"))?;
        let output_dir = if cb.has_value("output_dir") { check_add_dir_separator(&cb.value("output_dir").unwrap()) } else { path_from_file(&bam_file) };
        let output_id = cb.value("output_id");
        let ref_gen_version = if cb.has_value("ref_genome_version") {
            RefGenomeVersion::from_java(&cb.value("ref_genome_version").unwrap())?
        } else {
            derive_ref_genome_version(&bam_file)?
        };
        // BT10: the Ensembl cache takes RefGenomeVersion.from(configBuilder), V37 when the option is not
        // given, though the BAM's version was derived above. Corrected: the BAM's version.
        let gene_ref_gen_version = if compat::java("BT10") && !cb.has_value("ref_genome_version") {
            RefGenomeVersion::V37
        } else {
            ref_gen_version
        };
        crate::info!("refGenome({}), bam({})", ref_gen_version.name(), bam_file);
        crate::info!("output({})", output_dir);
        let partition_size = cb.integer("partition_size");
        let map_quality_threshold = cb.integer("map_qual_threshold");
        let base_quality_threshold = cb.integer("base_qual_threshold");

        // BT12: BaseRegion.checkMergeOverlaps cuts a region short when the next one lies inside it.
        // Corrected: the larger end is kept.
        let target_regions = match cb.value("regions_file") {
            Some(f) => hmf_common::regions::load_chr_base_regions(&f, !compat::java("BT12"))?,
            None => JavaHashMap::new(),
        };
        let only_target_regions = !target_regions.is_empty() && cb.has_flag("only_target");
        let max_coverage = if !target_regions.is_empty() {
            if cb.has_value("max_coverage") { cb.integer("max_coverage") } else { DEFAULT_PANEL_MAX_COVERAGE }
        } else {
            cb.integer("max_coverage")
        };
        if max_coverage < 0 {
            // Java: NegativeArraySizeException creating the coverage frequencies
            bail!("java.lang.NegativeArraySizeException: {}", max_coverage + 1);
        }
        if partition_size <= 0 {
            // Java: ArithmeticException (/ by zero) or a never-ending partitioning loop
            bail!("invalid partition size: {partition_size}");
        }
        let unmappable_regions = load_unmappable_regions(ref_gen_version)?;
        let specific_chr_regions = SpecificRegions::from_java(cb.value("specific_regions").as_deref(), cb.value("specific_chr").as_deref())
            .inspect_err(|_| {
                crate::error!(
                    "invalid specific config: regions({}) chromosomes({})",
                    cb.value("specific_regions").unwrap_or_else(|| "null".into()),
                    cb.value("specific_chr").unwrap_or_else(|| "null".into())
                );
            })?;
        for r in &specific_chr_regions.regions {
            crate::debug!("filtering for specific region: {r}");
        }
        let log_read_ids = match cb.value("log_read_ids") {
            Some(v) if cb.has_value("log_read_ids") => hmf_common::java_text::split_all(&v, ';').into_iter().map(|s| s.to_string()).collect(),
            _ => Vec::new(),
        };
        let config = MetricsConfig {
            sample_id,
            bam_file,
            ref_genome_file,
            ref_gen_version,
            gene_ref_gen_version,
            map_quality_threshold,
            base_quality_threshold,
            max_coverage,
            partition_size,
            unmappable_regions,
            target_regions,
            only_target_regions,
            ensembl_data_dir: if cb.has_value("ensembl_data_dir") { cb.value("ensembl_data_dir") } else { None },
            driver_gene_panel: if cb.has_value("driver_gene_panel") { cb.value("driver_gene_panel") } else { None },
            exclude_zero_coverage: cb.has_flag("exclude_zero_coverage"),
            write_off_target: cb.has_flag("write_off_target"),
            no_duplicates: cb.has_flag("no_duplicates"),
            high_fragment_overlap_threshold: cb.integer("off_target_frag_overlap_threshold"),
            partition_read_count_check: cb.integer("partition_read_count_check"),
            output_dir,
            output_id,
            threads: cb.integer("threads"),
            specific_chr_regions,
            log_read_ids,
            perf_debug: cb.has_flag("perf_debug"),
            file_buffer_kb: cb.integer("file_buffer_kb").clamp(0, 1 << 20),
        };
        // isValid(): the BAM and reference files exist (checked by the parser as PATH items)
        Ok(config)
    }

    pub fn expect_duplicates(&self) -> bool {
        !self.no_duplicates
    }

    /// formFilename(fileType): <output_dir><sample>.bam_metric.<type>[.<output_id>].tsv
    pub fn form_filename(&self, file_type: &str) -> String {
        let mut f = format!("{}{}{}.{file_type}", self.output_dir, self.sample_id, hmf_common::metrics::BAM_METRICS_FILE_ID);
        if let Some(id) = &self.output_id {
            f.push('.');
            f.push_str(id);
        }
        f + ".tsv"
    }

    /// TargetRegions: the regions of a chromosome
    pub fn target_regions_of(&self, chromosome: &str) -> Option<&Vec<BaseRegion>> {
        self.target_regions.get(&chromosome.to_string())
    }

    /// For tests: MetricsConfig(maxCoverage) of the Java tests
    #[cfg(test)]
    pub fn for_tests(max_coverage: i32) -> MetricsConfig {
        MetricsConfig {
            sample_id: "SAMPLE_ID".into(),
            bam_file: String::new(),
            ref_genome_file: String::new(),
            ref_gen_version: RefGenomeVersion::V37,
            gene_ref_gen_version: RefGenomeVersion::V37,
            map_quality_threshold: DEFAULT_MAP_QUAL_THRESHOLD,
            base_quality_threshold: DEFAULT_BASE_QUAL_THRESHOLD,
            max_coverage,
            partition_size: DEFAULT_CHR_PARTITION_SIZE,
            unmappable_regions: Vec::new(),
            target_regions: JavaHashMap::new(),
            only_target_regions: false,
            ensembl_data_dir: None,
            driver_gene_panel: None,
            exclude_zero_coverage: false,
            write_off_target: false,
            no_duplicates: false,
            high_fragment_overlap_threshold: 0,
            partition_read_count_check: 0,
            output_dir: String::new(),
            output_id: None,
            threads: 0,
            specific_chr_regions: SpecificRegions::default(),
            log_read_ids: Vec::new(),
            perf_debug: false,
            file_buffer_kb: 0,
        }
    }
}

const UNMAPPABLE_37: &str = include_str!("../resources/genome_unmappable_regions.37.bed");
const UNMAPPABLE_38: &str = include_str!("../resources/genome_unmappable_regions.38.bed");

/// MetricsConfig.loadUnmappableRegions: the resource BED through BedFileReader.loadBedFile(lines)
fn load_unmappable_regions(version: RefGenomeVersion) -> Result<Vec<ChrBaseRegion>> {
    let text = if version.is_37() { UNMAPPABLE_37 } else { UNMAPPABLE_38 };
    let lines = hmf_common::java_text::lines_of(text);
    hmf_common::regions::load_bed_lines(&lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_as_java() {
        let dir = std::env::temp_dir();
        let d = dir.to_str().unwrap();
        let Parsed::Run(cb) = ConfigBuilder::parse(&strings(&["-bam_file", d, "-ref_genome", d, "-threads", "-2", "-only_target"])).unwrap() else { panic!() };
        assert_eq!(cb.integer("threads"), -2);
        assert_eq!(cb.integer("max_coverage"), 250);
        assert!(cb.has_flag("only_target") && !cb.has_value("max_coverage"));
        assert!(ConfigBuilder::parse(&strings(&["-bam_file", d])).is_err()); // ref_genome required
        assert!(ConfigBuilder::parse(&strings(&["-bam_file", d, "-ref_genome", d, "-sample", "-x"])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-bam_file", d, "-ref_genome", d, "-threads", "two"])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-bam_file", "/no/such/file", "-ref_genome", d])).is_err());
        assert!(ConfigBuilder::parse(&strings(&["-bam_file", d, "-ref_genome", d, "-unknown", "1"])).is_err());
    }

    #[test]
    fn unmappable_regions() {
        let r38 = load_unmappable_regions(RefGenomeVersion::V38).unwrap();
        assert!(!r38.is_empty() && r38[0].chromosome == "1" && r38[0].start == 1);
        let r37 = load_unmappable_regions(RefGenomeVersion::V37).unwrap();
        assert!(!r37.is_empty() && r37.iter().all(|r| !r.chromosome.starts_with("chr")));
    }

    #[test]
    fn help_lines() {
        let lines = item_lines(&items());
        assert_eq!(lines[0], "registered config items:");
        assert!(lines.contains(&" -bam_file: type(PATH) REQUIRED desc(Input BAM file)".to_string()));
        assert!(lines.contains(&" -ref_genome_version: default(V37) desc(Ref genome version, 37 or 38)".to_string()));
        assert!(lines.contains(&" -only_target: type(FLAG) desc(Only capture metrics within the specific regions file)".to_string()));
    }
}
