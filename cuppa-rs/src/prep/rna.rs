//! GeneExpressionPrep (GENE_EXP: log(AdjTPM + 1) per gene of Isofox's gene data) and AltSpliceJunctionPrep
//! (ALT_SJ: fragment counts of Isofox's alternative splice junctions at the reference sites).

use std::collections::HashMap;

use anyhow::{anyhow, bail, Result};
use hmf_common::{error, warn};

use super::config::PrepConfig;
use super::items::{DataItem, ItemType, Source};

/// FileDelimiters.inferFileDelimiter
fn infer_delim(path: &str) -> char {
    if path.ends_with(".csv") || path.ends_with(".csv.gz") { ',' } else { '\t' }
}

/// None where Java returns null (any failure; logged)
pub fn extract_gene_expression(config: &PrepConfig, sample: &str) -> Option<Vec<DataItem>> {
    let path = config.gene_expression(sample);
    let run = || -> Result<Vec<DataItem>> {
        let delim = infer_delim(&path);
        let lines = hmf_common::java_text::read_lines(&path).map_err(|e| anyhow!("java.nio.file.NoSuchFileException: {path} ({e})"))?;
        let Some(header) = lines.first() else { bail!("java.lang.IndexOutOfBoundsException: empty file {path}") };
        let f = hmf_common::java_text::fields_index_map(header, delim);
        let gi = *f.get("GeneName").ok_or_else(|| anyhow!("java.lang.NullPointerException: no GeneName column"))?;
        let ti = *f.get("AdjTPM").ok_or_else(|| anyhow!("java.lang.NullPointerException: no AdjTPM column"))?;
        let mut items = Vec::new();
        for line in &lines[1..] {
            let v: Vec<&str> = hmf_common::java_text::split_all(line, delim);
            let g = v.get(gi).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: {line}"))?;
            let t = v.get(ti).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: {line}"))?;
            let tpm = hmf_common::java_text::parse_double(t).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{t}\""))?;
            items.push(DataItem::sci(Source::Rna, ItemType::Expression, g, (tpm + 1.0).ln()));
        }
        Ok(items)
    };
    match run() {
        Ok(items) => Some(items),
        Err(e) => {
            error!("sample({sample}) failed to extract category(GENE_EXP): {e:#}");
            None
        }
    }
}

/// NovelSpliceJunctionFile.formKey
fn form_key(chr: &str, start: i32, end: i32) -> String {
    format!("{chr};{start};{end}")
}

/// AltSpliceJunctionPrep.loadRefAltSjIndices
fn load_ref_sites(path: &str) -> Result<HashMap<String, usize>> {
    let lines = hmf_common::java_text::read_lines(path).map_err(|e| anyhow!("Failed to read ref_alt_sj_sites file({path}): {e}"))?;
    let delim = infer_delim(path);
    let Some(header) = lines.first() else { bail!("Failed to read ref_alt_sj_sites file({path}): empty") };
    let f = hmf_common::java_text::fields_index_map(header, delim);
    let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("Failed to read ref_alt_sj_sites file({path}): no column {n}"));
    let (ci, si, ei) = (col("Chromosome")?, col("PosStart")?, col("PosEnd")?);
    let mut map = HashMap::new();
    for (i, line) in lines[1..].iter().enumerate() {
        let v: Vec<&str> = hmf_common::java_text::split_all(line, delim);
        let get = |k: usize| v.get(k).copied().ok_or_else(|| anyhow!("Failed to read ref_alt_sj_sites file({path}): short line"));
        let int = |s: &str| hmf_common::java_text::parse_int(s).ok_or_else(|| anyhow!("Failed to read ref_alt_sj_sites file({path}): {s}"));
        map.insert(form_key(get(ci)?, int(get(si)?)?, int(get(ei)?)?), i);
    }
    Ok(map)
}

struct Junction {
    chromosome: String,
    start: i32,
    end: i32,
    fragments: i32,
}

/// NovelSpliceJunctionFile.read (plain text, the delimiter from the name; the fields Java parses checked)
fn read_junctions(path: &str) -> Result<Vec<Junction>> {
    let text = std::fs::read_to_string(path).map_err(|e| anyhow!("{path}: {e}"))?;
    let mut lines = text.lines();
    let delim = infer_delim(path);
    let header = lines.next().unwrap_or("");
    let f = hmf_common::java_text::fields_index_map(header, delim);
    let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {n} in {path}"));
    let mut out = Vec::new();
    const TYPES: [&str; 11] = [
        "SKIPPED_EXONS",
        "NOVEL_5_PRIME",
        "NOVEL_3_PRIME",
        "NOVEL_EXON",
        "NOVEL_INTRON",
        "MIXED_TRANS",
        "INTRONIC",
        "EXON_INTRON",
        "CIRCULAR",
        "CANONICAL",
        "UNKNOWN",
    ];
    const CONTEXTS: [&str; 5] = ["SPLICE_JUNC", "EXONIC", "INTRONIC", "MIXED", "UNKNOWN"];
    for line in lines {
        let v: Vec<&str> = hmf_common::java_text::split_all(line, delim);
        let g = |n: &str| -> Result<&str> {
            let i = col(n)?;
            v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len()))
        };
        let int = |s: &str| hmf_common::java_text::parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""));
        let chromosome = g("Chromosome")?.to_string();
        let start = int(g("SjStart")?)?;
        let end = int(g("SjEnd")?)?;
        g("GeneName")?;
        let t = g("Type")?;
        if !TYPES.contains(&t) {
            bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.rna.AltSpliceJunctionType.{t}");
        }
        for opt in ["ExonStart", "ExonEnd"] {
            if f.contains_key(opt) {
                int(g(opt)?)?;
            }
        }
        let fragments = int(g("FragCount")?)?;
        int(g("DepthStart")?)?;
        int(g("DepthEnd")?)?;
        for c in ["RegionStart", "RegionEnd"] {
            let s = g(c)?;
            if !CONTEXTS.contains(&s) {
                bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.{s}");
            }
        }
        if f.contains_key("CohortFrequency") {
            int(g("CohortFrequency")?)?;
        }
        out.push(Junction { chromosome, start, end, fragments });
    }
    Ok(out)
}

/// AltSpliceJunctionPrep: None where Java returns null (no reference sites, or no junction file)
pub fn extract_alt_sj(config: &PrepConfig, sample: &str) -> Result<Option<Vec<DataItem>>> {
    // the reference sites are loaded for each sample (a new AltSpliceJunctionPrep per sample); Java exits on
    // a failure to load them
    let sites = match &config.alt_sj_sites {
        Some(p) => load_ref_sites(p)?,
        None => bail!("Failed to read ref_alt_sj_sites file(null): java.lang.NullPointerException"),
    };
    if sites.is_empty() {
        return Ok(None);
    }
    let path = config.alt_splice_junctions(sample);
    if !std::path::Path::new(&path).exists() {
        return Ok(None);
    }
    let junctions = read_junctions(&path)?;
    let mut items = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for j in junctions {
        let key = form_key(&j.chromosome, j.start, j.end);
        if !sites.contains_key(&key) || seen.contains(&key) {
            continue;
        }
        items.push(DataItem::int(Source::Rna, ItemType::AltSj, &key, j.fragments as i64));
        seen.insert(key);
    }
    if items.is_empty() {
        warn!("sample({sample}) had no matching alt-SJs of the {} provided in configItem(-ref_alt_sj_sites)", sites.len());
    }
    Ok(Some(items))
}
