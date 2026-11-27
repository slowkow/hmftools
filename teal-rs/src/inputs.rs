//! The pipeline files TealPipelineApp reads: PURPLE's purity (PurityContextFile, with its QC file), COBALT's
//! GC median read depths (CobaltGcMedianFile) and BamTools' metrics summary (BamMetricSummary). Only the
//! values TEAL uses are kept; the other fields are parsed where Java would fail on them.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;

use crate::config::java_parse_double;

fn read_lines(path: &str) -> Result<Vec<String>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    Ok(text.lines().map(|l| l.to_string()).collect())
}

fn parse_double(s: &str) -> Result<f64> {
    java_parse_double(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

fn parse_int(s: &str) -> Result<i64> {
    crate::config::java_parse_int(s).map(|v| v as i64).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

fn parse_long(s: &str) -> Result<i64> {
    let digits = s.strip_prefix('-').or_else(|| s.strip_prefix('+')).unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        bail!("java.lang.NumberFormatException: For input string: \"{s}\"");
    }
    s.parse::<i64>().map_err(|_| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

/// createFieldsIndexMap
fn fields(header: &str) -> HashMap<String, usize> {
    header.split('\t').enumerate().map(|(i, f)| (f.to_string(), i)).collect()
}

fn column<'a>(index: &HashMap<String, usize>, values: &[&'a str], name: &str) -> Result<&'a str> {
    let i = *index.get(name).ok_or_else(|| anyhow!("java.lang.NullPointerException (column {name} missing)"))?;
    values.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException (column {name})"))
}

/// PurityContextFile.read: the best fit's purity and ploidy
pub fn read_purity(purple_dir: &str, sample: &str) -> Result<(f64, f64)> {
    // PurpleQCFile.read comes first
    let qc_path = format!("{purple_dir}/{sample}.purple.qc");
    let qc = read_lines(&qc_path)?;
    let get = |field: &str| -> Option<String> {
        qc.iter().find_map(|l| {
            let v: Vec<&str> = l.split('\t').collect();
            if v[0] == field { v.get(1).map(|s| s.to_string()) } else { None }
        })
    };
    for f in ["CopyNumberSegments", "UnsupportedCopyNumberSegments", "DeletedGenes", "AmberMeanDepth"] {
        if let Some(v) = get(f) {
            parse_int(&v)?;
        }
    }
    for f in ["Purity", "Contamination", "LohPercent", "TincLevel"] {
        if let Some(v) = get(f) {
            parse_double(&v)?;
        }
    }
    for f in ["AmberGender", "CobaltGender"] {
        if get(f).is_none() {
            bail!("java.lang.NullPointerException: Name is null ({f} missing in {qc_path})");
        }
    }
    let path = format!("{purple_dir}/{sample}.purple.purity.tsv");
    let lines = read_lines(&path)?;
    if lines.len() < 2 {
        bail!("java.lang.IndexOutOfBoundsException (no values in {path})");
    }
    let index = fields(&lines[0]);
    let values: Vec<&str> = lines[1].split('\t').collect();
    for f in [
        "normFactor",
        "score",
        "diploidProportion",
        "polyclonalProportion",
        "minPurity",
        "maxPurity",
        "minPloidy",
        "maxPloidy",
        "minDiploidProportion",
        "maxDiploidProportion",
        "somaticPenalty",
        "msIndelsPerMb",
        "tmbPerMb",
    ] {
        parse_double(column(&index, &values, f)?)?;
    }
    for f in ["tml", "svTumorMutationalBurden"] {
        parse_int(column(&index, &values, f)?)?;
    }
    for f in ["gender", "status", "wholeGenomeDuplication", "msStatus", "tmlStatus", "tmbStatus", "runMode", "targeted"] {
        column(&index, &values, f)?;
    }
    let purity = parse_double(column(&index, &values, "purity")?)?;
    let ploidy = parse_double(column(&index, &values, "ploidy")?)?;
    Ok((purity, ploidy))
}

/// CobaltGcMedianFile.read: the sample mean and the median read depth of GC bucket 50 (-1 when absent)
pub fn read_gc_median(cobalt_dir: &str, sample: &str) -> Result<(f64, f64)> {
    let path = format!("{cobalt_dir}/{sample}.cobalt.gc.median.tsv");
    let lines = read_lines(&path)?;
    let mut i = 1;
    let mut mean = 0.0;
    if lines.len() > i {
        let v: Vec<&str> = lines[i].split('\t').collect();
        i += 1;
        mean = parse_double(v[0])?;
        parse_double(v.get(1).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException ({path})"))?)?;
    }
    i += 1;
    let mut gc50 = -1.0;
    for line in lines.iter().skip(i) {
        let v: Vec<&str> = line.split('\t').collect();
        if v.len() == 2 {
            let bucket = parse_int(v[0])?;
            let depth = parse_double(v[1])?;
            if bucket == 50 {
                gc50 = depth;
            }
        }
    }
    Ok((mean, gc50))
}

/// BamMetricSummary.read: DuplicatePercent
pub fn read_duplicate_percent(path: &str) -> Result<f64> {
    let lines = read_lines(path)?;
    if lines.len() < 2 {
        bail!("java.lang.IndexOutOfBoundsException (no values in {path})");
    }
    let index = fields(&lines[0]);
    let values: Vec<&str> = lines[1].split('\t').collect();
    for (i, c) in lines[0].split('\t').enumerate() {
        if let Some(level) = c.strip_prefix("DepthCoverage_") {
            parse_int(level)?;
            parse_double(values.get(i).copied().unwrap_or(""))?;
        }
    }
    if !index.contains_key("UnmappedPercent") && !index.contains_key("UnpairedPercent") {
        bail!("java.lang.NullPointerException (column UnmappedPercent missing in {path})");
    }
    for f in ["TotalRegionBases", "TotalReads", "DuplicateReads", "DualStrandReads"] {
        parse_long(column(&index, &values, f)?)?;
    }
    if index.contains_key("OffTargetReads") {
        parse_long(column(&index, &values, "OffTargetReads")?)?;
    }
    for f in ["MeanCoverage", "StdDevCoverage", "MedianCoverage", "MadCoverage", "LowMapQualPercent", "LowBaseQualPercent", "OverlappingReadPercent", "CappedCoverage"] {
        parse_double(column(&index, &values, f)?)?;
    }
    parse_double(column(&index, &values, "DuplicatePercent")?)
}
