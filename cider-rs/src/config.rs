//! The command line, parsed as the Java ConfigBuilder does (`-name value` and `-flag`, required items,
//! path checks, integer and decimal checks, `-version` and `-help` as the sole argument), and CiderParams.

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;

use hmf_common::java_text::{parse_double, parse_int};
use hmf_common::regions::RefGenomeVersion;


pub const APP_NAME: &str = "Cider";
pub const VERSION: &str = "1.2";

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

/// CiderParams.registerConfig and ConfigUtils.addLoggingOptions
fn items() -> Vec<Item> {
    use ItemType::*;
    let it = |name, ty, required, desc, default| Item { name, ty, required, desc, default };
    vec![
        it("sample", Str, false, "Sample ID", None),
        it("bam", Path, true, "input BAM / CRAM file", None),
        it("ref_genome", Path, false, "Path to reference genome fasta files", None),
        it("output_dir", Str, false, "Output directory", None),
        it("threads", Integer, false, "Number of threads, 0 or 1 not multi-threaded", Some("1")),
        it("ref_genome_version", Str, false, "Ref genome version, 37 or 38", Some("V37")),
        it("max_fragment_length", Integer, false, "Approximate maximum fragment length", Some("1000")),
        it("min_base_quality", Integer, false, "Minimum quality for a base to be considered", Some("25")),
        it("write_cider_bam", Flag, false, "Write a output BAM file containing all CDR3 reads", None),
        it("report_match_ref_seq", Flag, false, "Report VDJ sequences that match reference genome", None),
        it("num_trim_bases", Integer, false, "Number of bases to trim on each side of reads", Some("0")),
        it(
            "max_low_qual_base_fraction",
            Decimal,
            false,
            "Maximum fraction of bases in a read that can be low quality. Reads that exceed this limit are discarded",
            Some("0.1"),
        ),
        it(
            "max_reads_per_gene",
            Integer,
            false,
            "Maximum number of reads per gene. If number of reads exceed this limit, they are downsampled.",
            Some("100000"),
        ),
        it("primer_csv", Path, false, "Path to csv file containing primers", None),
        it("primer_mismatch_max", Integer, false, "Maximum number of mismatch bases for matching primer sequence", Some("0")),
        it("bwa_lib", Path, false, "Path to BWA library", None),
        it("bwa_index_image", Path, false, "Reference genome BWA-MEM index GATK image file", None),
        it("log_debug", Flag, false, "Log at DEBUG level", None),
        it("log_level", Str, false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", None),
    ]
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
                if arg.starts_with('-') && items[i].ty != ItemType::Integer && items[i].ty != ItemType::Decimal {
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
                ItemType::Decimal => {
                    let v = value.map(|s| s.as_str()).or(item.default);
                    if v.and_then(parse_double).is_none() {
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
        self.items.iter().find(|i| i.name == name).and_then(|i| i.default.map(|s| s.to_string()))
    }
    pub fn integer(&self, name: &str) -> i32 {
        self.value(name).as_deref().and_then(parse_int).expect("checked by parse")
    }
    pub fn decimal(&self, name: &str) -> f64 {
        self.value(name).as_deref().and_then(parse_double).expect("checked by parse")
    }
}

/// CiderParams
#[derive(Clone, Debug)]
pub struct CiderParams {
    pub sample_id: String,
    pub bam_path: String,
    pub ref_genome_path: Option<String>,
    pub output_dir: String,
    pub thread_count: i32,
    pub ref_genome_version: RefGenomeVersion,
    pub approx_max_fragment_length: i32,
    pub min_base_quality: i32,
    pub write_filtered_bam: bool,
    pub report_match_ref_seq: bool,
    pub num_bases_to_trim: i32,
    pub max_low_qual_base_fraction: f64,
    pub max_read_count_per_gene: i32,
    pub primer_csv: Option<String>,
    pub primer_mismatch_max: i32,
    pub bwa_index_image_path: Option<String>,
}

/// FileWriterUtils.checkAddDirSeparator
pub fn check_add_dir_separator(path: &str) -> String {
    if path.is_empty() || path.ends_with('/') { path.to_string() } else { format!("{path}/") }
}

impl CiderParams {
    /// CiderParams.fromConfigBuilder
    pub fn from(cb: &ConfigBuilder) -> Result<CiderParams> {
        let ref_genome_path = cb.value("ref_genome");
        let default_index = ref_genome_path.as_ref().map(|r| format!("{r}.img"));
        let sample_id = cb.value("sample").ok_or_else(|| anyhow!("java.lang.NullPointerException (-sample)"))?;
        let output_dir = cb.value("output_dir").ok_or_else(|| anyhow!("java.lang.NullPointerException (-output_dir)"))?;
        let ref_genome_version = RefGenomeVersion::from_java(&cb.value("ref_genome_version").unwrap_or_default())?;
        Ok(CiderParams {
            sample_id,
            bam_path: cb.value("bam").expect("required"),
            ref_genome_path,
            output_dir: check_add_dir_separator(&output_dir),
            thread_count: cb.integer("threads"),
            ref_genome_version,
            approx_max_fragment_length: cb.integer("max_fragment_length"),
            min_base_quality: cb.integer("min_base_quality"),
            write_filtered_bam: cb.has_flag("write_cider_bam"),
            report_match_ref_seq: cb.has_flag("report_match_ref_seq"),
            num_bases_to_trim: cb.integer("num_trim_bases"),
            max_low_qual_base_fraction: cb.decimal("max_low_qual_base_fraction"),
            max_read_count_per_gene: cb.integer("max_reads_per_gene"),
            primer_csv: cb.value("primer_csv"),
            primer_mismatch_max: cb.integer("primer_mismatch_max"),
            bwa_index_image_path: cb.value("bwa_index_image").or(default_index),
        })
    }

    pub fn file(&self, suffix: &str) -> String {
        format!("{}{}{}", self.output_dir, self.sample_id, suffix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::*;

    fn args(s: &str) -> Vec<String> {
        s.split_whitespace().map(|x| x.to_string()).collect()
    }

    #[test]
    fn parses_pipeline_command() {
        let bam = std::env::current_exe().unwrap().to_string_lossy().into_owned();
        let a = args(&format!("-sample S -bam {bam} -ref_genome_version 38 -write_cider_bam -threads 2 -output_dir out"));
        let Parsed::Run(cb) = ConfigBuilder::parse(&a).unwrap() else { panic!() };
        let p = CiderParams::from(&cb).unwrap();
        assert_eq!(p.sample_id, "S");
        assert_eq!(p.thread_count, 2);
        assert!(p.write_filtered_bam && !p.report_match_ref_seq);
        assert_eq!(p.ref_genome_version, RefGenomeVersion::V38);
        assert_eq!(p.output_dir, "out/");
        assert_eq!(p.max_low_qual_base_fraction, DEFAULT_MAX_LOW_QUAL_BASES_FRACTION);
        assert_eq!(p.max_read_count_per_gene, DEFAULT_MAX_READ_COUNT_PER_GENE);
        assert_eq!(p.approx_max_fragment_length, DEFAULT_MAX_FRAGMENT_LENGTH);
        assert_eq!(p.min_base_quality, DEFAULT_MIN_BASE_QUALITY);
        assert!(p.bwa_index_image_path.is_none());
        assert_eq!(p.file(".cider.bam"), "out/S.cider.bam");
    }

    #[test]
    fn rejects_bad_options() {
        assert!(ConfigBuilder::parse(&args("-sample S -bam /nonexistent/x.bam -output_dir o")).is_err());
        assert!(ConfigBuilder::parse(&args("-unknown 1")).is_err());
        assert!(ConfigBuilder::parse(&args("-sample S")).is_err());
    }
}
