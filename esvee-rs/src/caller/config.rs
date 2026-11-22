//! Command line of the caller step (CallerConfig with hmf-common's ConfigBuilder rules).

use super::filters::*;
use super::types::*;
use anyhow::{bail, Result};
use std::collections::HashMap;

#[derive(Clone, Copy, PartialEq)]
enum ItemType {
    Flag,
    Integer,
    Text,
    Path,
    Paths,
}

/// every option CallerConfig.registerConfig registers: (name, type, default)
const ITEMS: &[(&str, ItemType, Option<&str>)] = &[
    ("tumor", ItemType::Text, None),
    ("reference", ItemType::Text, None),
    ("input_vcf", ItemType::Path, None),
    ("esvee_prep_dir", ItemType::Paths, None),
    ("manual_ref_depth", ItemType::Integer, Some("0")),
    ("sequencing_type", ItemType::Text, Some("ILLUMINA")),
    ("output_dir", ItemType::Text, None),
    ("output_id", ItemType::Text, None),
    ("log_debug", ItemType::Flag, None),
    ("log_level", ItemType::Text, None),
    ("ref_genome_version", ItemType::Text, None),
    ("specific_chr", ItemType::Text, None),
    ("specific_regions", ItemType::Text, None),
    ("pon_sv_file", ItemType::Path, None),
    ("pon_sgl_file", ItemType::Path, None),
    ("pon_margin", ItemType::Integer, Some("3")),
    ("pon_sgl_margin", ItemType::Integer, Some("10")),
    ("artefact_pon_sv_file", ItemType::Path, None),
    ("artefact_pon_sgl_file", ItemType::Path, None),
    ("write_breakend_tsv", ItemType::Flag, None),
    ("known_hotspot_file", ItemType::Path, None),
    ("filter_sgls", ItemType::Flag, None),
    ("min_qual", ItemType::Integer, Some("30")),
    ("hotspot_min_qual", ItemType::Integer, Some("30")),
    ("min_support", ItemType::Integer, Some("4")),
    ("hotspot_min_support", ItemType::Integer, Some("2")),
    ("sgl_min_support", ItemType::Integer, Some("6")),
    ("repeat_mask_file", ItemType::Path, None),
    ("target_regions_bed", ItemType::Path, None),
];

const REQUIRED: &[&str] = &["ref_genome_version"];

/// The values given on the command line
pub struct ConfigValues {
    values: HashMap<&'static str, String>,
}

impl ConfigValues {
    /// ConfigBuilder.parseCommandLine and isValid. An error is the list of messages Java logs.
    pub fn parse(args: &[String]) -> std::result::Result<ConfigValues, Vec<String>> {
        let mut values: HashMap<&'static str, String> = HashMap::new();
        let mut matched: Option<(&'static str, ItemType)> = None;
        for argument in args {
            if let Some((name, item_type)) = matched {
                if argument.starts_with('-') && item_type != ItemType::Integer {
                    return Err(vec![format!("config item({name}) has invalid argument: {argument}")]);
                }
                values.insert(name, argument.clone());
                matched = None;
            } else if let Some(name) = argument.strip_prefix('-') {
                let Some(item) = ITEMS.iter().find(|i| i.0 == name) else {
                    return Err(vec![format!("unregistered config item: {argument}")]);
                };
                if item.1 == ItemType::Flag {
                    values.insert(item.0, "true".to_string());
                } else {
                    matched = Some((item.0, item.1));
                }
            } else {
                return Err(vec![format!("expecting config name but encountered invalid argument: {argument}")]);
            }
        }

        let mut errors = Vec::new();
        for (name, item_type, default) in ITEMS {
            let value = values.get(name);
            if REQUIRED.contains(name) && value.is_none() {
                errors.push(format!("missing required config: {name}"));
                continue;
            }
            match item_type {
                ItemType::Path | ItemType::Paths => {
                    if let Some(v) = value {
                        let paths: Vec<&str> = if *item_type == ItemType::Paths { v.split(',').collect() } else { vec![v.as_str()] };
                        if paths.iter().any(|p| !p.contains('*') && !p.contains('$') && !std::path::Path::new(p).exists()) {
                            errors.push(format!("invalid path for config: {name} = {v}"));
                        }
                    }
                }
                ItemType::Integer => {
                    let v = value.map(String::as_str).or(*default).unwrap_or("");
                    if super::vcf::java_parse_int(v).is_none() {
                        errors.push(format!("invalid type for config: {name}: INTEGER = {v}"));
                    }
                }
                _ => {}
            }
        }
        if errors.is_empty() { Ok(ConfigValues { values }) } else { Err(errors) }
    }

    pub fn has_value(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }
    /// the value, or the default if there is one
    pub fn value(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str).or_else(|| ITEMS.iter().find(|i| i.0 == name).and_then(|i| i.2))
    }
    pub fn integer(&self, name: &str) -> i32 {
        super::vcf::java_parse_int(self.value(name).unwrap_or("0")).unwrap_or(0)
    }
    pub fn has_flag(&self, name: &str) -> bool {
        self.values.get(name).is_some_and(|v| v == "true")
    }
}

pub struct CallerConfig {
    pub tumor_id: Option<String>,
    pub reference_id: Option<String>,
    pub ref_gen_version: RefGenomeVersion,
    pub vcf_file: String,
    pub prep_dir: String,
    /// Java builds file names by string concatenation, so a missing output directory reads `null`
    pub output_dir: String,
    pub output_id: Option<String>,
    pub specific_chromosomes: Vec<String>,
    /// the regions of -specific_regions; when there are any they are used instead of the chromosomes
    pub specific_regions: Vec<ChrBaseRegion>,
    pub manual_ref_depth: i32,
    pub write_breakend_tsv: bool,
    pub sequencing_type: SequencingType,
}

/// FileWriterUtils.checkAddDirSeparator
fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

/// FileCommon.formOutputFile
pub fn form_output_file(output_dir: &str, sample_id: &str, app_stage: &str, file_type: &str, output_id: Option<&str>) -> String {
    match output_id {
        Some(id) => format!("{output_dir}{sample_id}.{app_stage}.{id}.{file_type}"),
        None => format!("{output_dir}{sample_id}.{app_stage}.{file_type}"),
    }
}

/// FileCommon.formPrepInputFilename: the file with the output id if it exists, else the one without
pub fn form_prep_input_filename(dir: &str, sample_id: &str, file_id: &str, output_id: Option<&str>) -> String {
    let with_id = form_output_file(dir, sample_id, "esvee.prep", file_id, output_id);
    if std::path::Path::new(&with_id).exists() {
        return with_id;
    }
    form_output_file(dir, sample_id, "esvee.prep", file_id, None)
}

/// SpecificRegions.loadSpecificRegions: `chr:start-end` or `chr:start:end` items separated by ';'
pub fn parse_specific_regions(text: &str) -> Result<Vec<ChrBaseRegion>> {
    let mut regions = Vec::new();
    for item in text.split(';') {
        let int = |s: &str| super::vcf::java_parse_int(s);
        let parts: Vec<&str> = item.split(':').collect();
        let region = match parts.as_slice() {
            [chromosome, start, end] => int(start).zip(int(end)).map(|(s, e)| ChrBaseRegion::new(chromosome, s, e)),
            [chromosome, positions] => {
                positions.split_once('-').and_then(|(s, e)| int(s).zip(int(e))).map(|(s, e)| ChrBaseRegion::new(chromosome, s, e))
            }
            _ => None,
        };
        match region {
            Some(r) if !r.chromosome.is_empty() && r.start > 0 && r.end >= r.start => regions.push(r),
            _ => bail!("invalid specific region: {item} (expected chr:start-end, with 1 <= start <= end)"),
        }
    }
    Ok(regions)
}

impl CallerConfig {
    pub fn new(values: &ConfigValues) -> Result<CallerConfig> {
        let tumor_id = values.value("tumor").map(str::to_string);
        let reference_id = values.value("reference").map(str::to_string);
        let output_dir = match values.value("output_dir") {
            Some(dir) => check_add_dir_separator(dir),
            // Java (C8): the missing directory is null, which string concatenation turns into file names
            // that start with "null". Corrected: the current directory.
            None if crate::common::compat::java("C8") => "null".to_string(),
            None => "./".to_string(),
        };
        let output_id = values.value("output_id").map(str::to_string);
        let file_sample_id = tumor_id.clone().or(reference_id.clone()).unwrap_or_else(|| "null".to_string());
        let vcf_file = match values.value("input_vcf") {
            Some(v) if values.has_value("input_vcf") => v.to_string(),
            _ => form_output_file(&output_dir, &file_sample_id, "esvee", "ref_depth.vcf.gz", output_id.as_deref()),
        };
        let prep_dir = match values.value("esvee_prep_dir") {
            Some(dir) => check_add_dir_separator(dir),
            None => output_dir.clone(),
        };
        let version_str = values.value("ref_genome_version").unwrap_or("");
        let Some(ref_gen_version) = RefGenomeVersion::from_str(version_str) else {
            bail!("Cannot resolve ref genome version: {version_str}");
        };
        let sequencing_type = if values.has_value("sequencing_type") {
            let v = values.value("sequencing_type").unwrap();
            match SequencingType::value_of(v) {
                Some(t) => t,
                None => bail!("IllegalArgumentException: No enum constant SequencingType.{v}"),
            }
        } else {
            SequencingType::Illumina
        };
        let specific_chromosomes = match values.value("specific_chr") {
            Some(s) if !s.is_empty() => {
                // String.split drops trailing empty items
                let mut items: Vec<String> = s.split(';').map(str::to_string).collect();
                while items.last().is_some_and(|x| x.is_empty()) {
                    items.pop();
                }
                items
            }
            _ => Vec::new(),
        };
        // Java (C7): -specific_regions is registered but never read. Corrected: as in hmf-common's
        // SpecificRegions, regions take the place of -specific_chr.
        let mut specific_chromosomes = specific_chromosomes;
        let mut specific_regions = Vec::new();
        if !crate::common::compat::java("C7") {
            if let Some(s) = values.value("specific_regions").filter(|s| !s.is_empty()) {
                specific_regions = parse_specific_regions(s)?;
                specific_chromosomes.clear();
            }
        }
        Ok(CallerConfig {
            tumor_id,
            reference_id,
            ref_gen_version,
            vcf_file,
            prep_dir,
            output_dir,
            output_id,
            specific_chromosomes,
            specific_regions,
            manual_ref_depth: values.integer("manual_ref_depth"),
            write_breakend_tsv: values.has_flag("write_breakend_tsv"),
            sequencing_type,
        })
    }

    pub fn has_tumor(&self) -> bool {
        self.tumor_id.is_some()
    }
    pub fn has_reference(&self) -> bool {
        self.reference_id.is_some()
    }
    pub fn germline_only(&self) -> bool {
        self.reference_id.is_some() && self.tumor_id.is_none()
    }
    pub fn file_sample_id(&self) -> String {
        self.tumor_id.clone().or(self.reference_id.clone()).unwrap_or_else(|| "null".to_string())
    }
}

/// FilterConstants.from(configBuilder)
pub fn filter_constants(values: &ConfigValues, version: RefGenomeVersion) -> FilterConstants {
    let targeted_mode = values.has_value("target_regions_bed");
    let min_qual = if values.has_value("min_qual") || !targeted_mode { values.integer("min_qual") } else { DEFAULT_MIN_QUAL_PANEL };
    FilterConstants::new(
        min_qual,
        values.integer("hotspot_min_qual"),
        values.integer("min_support"),
        values.integer("sgl_min_support"),
        values.integer("hotspot_min_support"),
        values.has_flag("filter_sgls"),
        version,
    )
}
