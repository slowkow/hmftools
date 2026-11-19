//! Append mode: adds the evidence of further samples to the variants of an existing SAGE VCF.
//!
//! Ported from append/SageAppendApplication, SageAppendConfig, RegionAppendTask, phase/AppendVariantPhaser,
//! phase/LpsReadCounts, vcf/CandidateSerialisation and the append constructor of vcf/VariantVCF. Java
//! reads each record with htsjdk, decodes every field with the header's types (`fullyDecode`) and writes
//! it again, so existing values are re-formatted; `VariantContext` here reproduces that round trip.

use crate::bamio::BamSlicer;
use crate::candidate::Candidate;
use crate::config::Config;
use crate::constants::*;
use crate::evidence::ReadContextCounter;
use crate::evidence_stage::{self, PhaseSink};
use crate::fmt::{java_fixed, vcf_double};
use crate::quality::{BqrMap, MsiJitter};
use crate::read_context::{create_context, VariantReadContext};
use crate::refdata::RefData;
use crate::refseq::{RefGenome, RefSequence};
use crate::sync::FragmentSync;
use crate::types::*;
use crate::vcfout;
use anyhow::{anyhow, bail, Context, Result};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{BufRead, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

pub const LPS_APPEND_INFO: &str = "LPSA";
const LPS_APPEND_HEADER: &str = "##FORMAT=<ID=LPSA,Number=1,Type=String,Description=\"Local Phase Set Append Info\">";
const MIN_PRIOR_VERSION: &str = "2.8";
const READ_CONTEXT_EXTENSION_LENGTH: i32 = 50;
const MISSING: &str = ".";

// --- VCF header field definitions ---

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FieldType {
    Integer,
    Float,
    Flag,
    Character,
    Str,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FieldCount {
    Fixed(i32),
    A,
    R,
    G,
    Unbounded,
}

#[derive(Clone, Copy, Debug)]
pub struct FieldDef {
    pub ty: FieldType,
    pub count: FieldCount,
}

/// The INFO and FORMAT definitions of a header. As in htsjdk, the first line for an id is the one used.
#[derive(Default)]
pub struct HeaderDefs {
    info: HashMap<String, FieldDef>,
    format: HashMap<String, FieldDef>,
}

/// The `key=value` pairs inside the `<...>` of a compound header line.
fn compound_attributes(body: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut key = String::new();
    let mut value = String::new();
    let mut in_value = false;
    let mut in_quotes = false;
    let mut escaped = false;
    for c in body.chars() {
        if in_quotes {
            if escaped {
                value.push(c);
                escaped = false;
            } else if c == '\\' {
                escaped = true;
            } else if c == '"' {
                in_quotes = false;
            } else {
                value.push(c);
            }
        } else if c == '"' && in_value {
            in_quotes = true;
        } else if c == '=' && !in_value {
            in_value = true;
        } else if c == ',' {
            out.push((std::mem::take(&mut key), std::mem::take(&mut value)));
            in_value = false;
        } else if in_value {
            value.push(c);
        } else {
            key.push(c);
        }
    }
    if !key.is_empty() {
        out.push((key, value));
    }
    out
}

/// An INFO or FORMAT header line: (is INFO, id, definition)
fn parse_field_line(line: &str) -> Option<(bool, String, FieldDef)> {
    let (is_info, rest) = if let Some(r) = line.strip_prefix("##INFO=<") {
        (true, r)
    } else if let Some(r) = line.strip_prefix("##FORMAT=<") {
        (false, r)
    } else {
        return None;
    };
    let body = rest.strip_suffix('>').unwrap_or(rest);
    let attrs = compound_attributes(body);
    let get = |k: &str| attrs.iter().find(|(key, _)| key == k).map(|(_, v)| v.as_str());
    let id = get("ID")?;
    let ty = match get("Type") {
        Some("Integer") => FieldType::Integer,
        Some("Float") => FieldType::Float,
        Some("Flag") => FieldType::Flag,
        Some("Character") => FieldType::Character,
        _ => FieldType::Str,
    };
    let count = match get("Number") {
        Some("A") => FieldCount::A,
        Some("R") => FieldCount::R,
        Some("G") => FieldCount::G,
        Some(n) => n.parse().map(FieldCount::Fixed).unwrap_or(FieldCount::Unbounded),
        None => FieldCount::Unbounded,
    };
    Some((is_info, id.to_string(), FieldDef { ty, count }))
}

/// htsjdk VCFStandardHeaderLines: (is INFO, id, number, type, description)
const STANDARD_HEADER_LINES: [(bool, &str, &str, &str, &str); 17] = [
    (false, "GT", "1", "String", "Genotype"),
    (false, "GQ", "1", "Integer", "Genotype Quality"),
    (false, "DP", "1", "Integer", "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"),
    (false, "PL", "G", "Integer", "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"),
    (false, "AD", "R", "Integer", "Allelic depths for the ref and alt alleles in the order listed"),
    (false, "FT", ".", "String", "Genotype-level filter"),
    (false, "PQ", "1", "Float", "Read-backed phasing quality"),
    (true, "END", "1", "Integer", "Stop position of the interval"),
    (true, "DB", "0", "Flag", "dbSNP Membership"),
    (true, "DP", "1", "Integer", "Approximate read depth; some reads may have been filtered"),
    (true, "SB", "1", "Float", "Strand Bias"),
    (true, "AF", "A", "Float", "Allele Frequency, for each ALT allele, in the same order as listed"),
    (true, "AC", "A", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed"),
    (true, "AN", "1", "Integer", "Total number of alleles in called genotypes"),
    (true, "MQ0", "1", "Integer", "Total Mapping Quality Zero Reads"),
    (true, "MQ", "1", "Float", "RMS Mapping Quality"),
    (true, "SOMATIC", "0", "Flag", "Somatic event"),
];

/// VCFStandardHeaderLines.repairStandardHeaderLines, which htsjdk applies when it reads a header: a line
/// for a standard field whose number or type is not the standard one is replaced by the standard line
/// (a different description alone is kept).
fn repair_standard_header_line(line: &str) -> String {
    let Some((is_info, id, def)) = parse_field_line(line) else { return line.to_string() };
    let Some(standard) = STANDARD_HEADER_LINES.iter().find(|s| s.0 == is_info && s.1 == id) else { return line.to_string() };
    let standard_line = format!(
        "##{}=<ID={},Number={},Type={},Description=\"{}\">",
        if is_info { "INFO" } else { "FORMAT" },
        standard.1,
        standard.2,
        standard.3,
        standard.4
    );
    let standard_def = parse_field_line(&standard_line).unwrap().2;
    // a fixed count is compared by value, the other count types by kind
    if def.count != standard_def.count || def.ty != standard_def.ty { standard_line } else { line.to_string() }
}

impl HeaderDefs {
    pub fn add_line(&mut self, line: &str) {
        if let Some((is_info, id, def)) = parse_field_line(line) {
            let map = if is_info { &mut self.info } else { &mut self.format };
            map.entry(id).or_insert(def);
        }
    }

    /// VariantContextUtils.getMetaDataForField: the FORMAT line if there is one, otherwise the INFO line
    fn field(&self, key: &str) -> Option<FieldDef> {
        self.format.get(key).or_else(|| self.info.get(key)).copied()
    }
}

// --- records ---

/// A genotype as htsjdk holds it: GT, the fields with their own accessors (AD, DP, GQ, PL), the genotype
/// filter and the remaining ("extended") attributes. Values are stored as they will be written.
#[derive(Clone, Debug, Default)]
pub struct Genotype {
    pub sample: String,
    pub gt: Option<String>,
    pub standard: BTreeMap<String, String>,
    pub filters: Option<String>,
    pub attributes: BTreeMap<String, String>,
}

/// A VCF record after htsjdk's `fullyDecode`: every value is held in the form the encoder writes.
#[derive(Clone, Debug, Default)]
pub struct VariantContext {
    pub chromosome: String,
    pub position: i32,
    pub id: String,
    pub ref_: String,
    pub alts: Vec<String>,
    /// log10 of the error probability; None when QUAL is missing
    pub log10_p_error: Option<f64>,
    /// None: filters not applied ("."); empty: PASS
    pub filters: Option<Vec<String>>,
    /// a flag has an empty value, a missing value is "."
    pub info: BTreeMap<String, String>,
    pub genotypes: Vec<Genotype>,
}

/// String.split(","): trailing empty strings are dropped
fn java_split(s: &str, delim: char) -> Vec<&str> {
    let mut parts: Vec<&str> = s.split(delim).collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// VariantContext.decodeOne, returning the value as VCFEncoder.formatVCFField writes it
fn decode_one(field: &str, s: &str, def: FieldDef) -> Result<String> {
    if s == MISSING {
        return Ok(MISSING.to_string());
    }
    let failed = || anyhow!("could not decode field {field} with value {s} of declared type {:?}", def.ty);
    Ok(match def.ty {
        FieldType::Character | FieldType::Str => s.to_string(),
        FieldType::Flag => {
            if !(s.eq_ignore_ascii_case("true") || s == "1") {
                bail!("FLAG field {field} cannot contain false values");
            }
            String::new()
        }
        FieldType::Integer => s.parse::<i32>().map_err(|_| failed())?.to_string(),
        FieldType::Float => vcf_double(s.parse::<f64>().map_err(|_| failed())?),
    })
}

fn decode_list(field: &str, parts: &[&str], def: FieldDef) -> Result<String> {
    if parts.is_empty() {
        return Ok(MISSING.to_string());
    }
    let decoded: Result<Vec<String>> = parts.iter().map(|p| decode_one(field, p, def)).collect();
    Ok(decoded?.join(","))
}

fn parse_ints(field: &str, value: &str) -> Result<String> {
    let ints: Result<Vec<String>> = value
        .split(',')
        .map(|x| x.parse::<i32>().map(|v| v.to_string()).map_err(|_| anyhow!("invalid {field} value {value}")))
        .collect();
    Ok(ints?.join(","))
}

/// AbstractVCFCodec.parseFilters
fn parse_filters(s: &str) -> Option<Vec<String>> {
    if s == MISSING {
        None
    } else if s == "PASS" {
        Some(Vec::new())
    } else {
        Some(s.split(';').map(|x| x.to_string()).collect())
    }
}

/// VCFEncoder.isMissingValue
fn is_missing_value(s: &str) -> bool {
    s.bytes().all(|b| b == b'.' || b == b',')
}

const STANDARD_GENOTYPE_KEYS: [&str; 4] = ["AD", "DP", "GQ", "PL"];

impl VariantContext {
    /// One data line, read as VCFCodec does and then fully decoded with the header's field types.
    pub fn parse(line: &str, defs: &HeaderDefs, samples: &[String]) -> Result<VariantContext> {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            bail!("invalid VCF record: {line}");
        }
        let position: i32 = cols[1].parse().with_context(|| format!("invalid VCF position: {}", cols[1]))?;
        let is_symbolic = |a: &str| a.starts_with('<') || a.contains('[') || a.contains(']') || a == "*";
        // Allele.create holds bases in upper case
        let allele = |a: &str| if is_symbolic(a) { a.to_string() } else { a.to_ascii_uppercase() };
        let alts: Vec<String> = if cols[4] == MISSING { Vec::new() } else { cols[4].split(',').map(allele).collect() };

        // parseQual
        let log10_p_error = if cols[5] == MISSING {
            None
        } else {
            let val: f64 = cols[5].parse().with_context(|| format!("invalid QUAL: {}", cols[5]))?;
            if val < 0.0 && (val - -1.0).abs() < 0.00005 { None } else { Some(val / -10.0) }
        };

        let mut info: BTreeMap<String, String> = BTreeMap::new();
        if cols[7] != MISSING {
            for field in cols[7].split(';') {
                if field.is_empty() {
                    continue;
                }
                let info_def = defs.info.get(field.split('=').next().unwrap()).copied();
                let (key, value) = match field.split_once('=') {
                    Some((key, value_string)) => {
                        // ParsingUtils.split keeps empty elements
                        let parts: Vec<&str> = value_string.split(',').collect();
                        let def = defs.field(key).ok_or_else(|| anyhow!("no header line for INFO field {key}"))?;
                        if parts.len() == 1 {
                            if info_def.is_some_and(|d| d.ty == FieldType::Flag) && parts[0] == "0" {
                                // a flag written as KEY=0 is not set
                                continue;
                            }
                            // KEY= is read as missing
                            let single = if parts[0].is_empty() { MISSING } else { parts[0] };
                            (key, decode_one(key, single, def)?)
                        } else {
                            (key, decode_list(key, &parts, def)?)
                        }
                    }
                    None => {
                        if defs.field(field).is_none() {
                            bail!("no header line for INFO field {field}");
                        }
                        // a key without a value is a flag unless the header says otherwise
                        let value = if info_def.is_some_and(|d| d.ty != FieldType::Flag) { MISSING.to_string() } else { String::new() };
                        (field, value)
                    }
                };
                info.insert(key.to_string(), value);
            }
        }

        let mut genotypes = Vec::new();
        if cols.len() > 9 {
            let keys: Vec<&str> = cols[8].split(':').collect();
            for (column, sample) in cols[9..].iter().zip(samples) {
                let values: Vec<&str> = column.split(':').collect();
                if keys.len() < values.len() {
                    bail!("there are too many keys for the sample {sample} at {}:{}", cols[0], position);
                }
                let mut g = Genotype { sample: sample.clone(), ..Default::default() };
                for (i, key) in keys.iter().enumerate() {
                    if *key == "GT" {
                        if i > 0 {
                            bail!("saw GT field at position {i}, but it must be at the first position for genotypes when present");
                        }
                        g.gt = Some(values[0].to_string());
                        continue;
                    }
                    // a value which is absent or missing is not added to the genotype
                    let Some(value) = values.get(i).filter(|v| **v != MISSING) else { continue };
                    match *key {
                        "FT" => {
                            if let Some(mut filters) = parse_filters(value).filter(|f| !f.is_empty()) {
                                filters.sort();
                                g.filters = Some(filters.join(";"));
                            }
                        }
                        "GQ" => {
                            if *value != "-1" {
                                let gq: f64 = value.parse().with_context(|| format!("invalid GQ value {value}"))?;
                                g.standard.insert("GQ".to_string(), (jround(gq) as i32).to_string());
                            }
                        }
                        "AD" | "PL" | "DP" => {
                            g.standard.insert(key.to_string(), parse_ints(key, value)?);
                        }
                        _ => {
                            let def = defs.field(key).ok_or_else(|| anyhow!("no header line for FORMAT field {key}"))?;
                            let decoded = if value.contains(',') {
                                decode_list(key, &java_split(value, ','), def)?
                            } else {
                                decode_one(key, value, def)?
                            };
                            g.attributes.insert(key.to_string(), decoded);
                        }
                    }
                }
                genotypes.push(g);
            }
        }

        Ok(VariantContext {
            chromosome: cols[0].to_string(),
            position,
            id: cols[2].to_string(),
            ref_: allele(cols[3]),
            alts,
            log10_p_error,
            filters: parse_filters(cols[6]),
            info,
            genotypes,
        })
    }

    /// VariantContext.isFiltered() && !getFilters().contains(PASS)
    fn is_filtered(&self) -> bool {
        self.filters.as_ref().is_some_and(|f| !f.is_empty() && !f.iter().any(|x| x == "PASS"))
    }

    /// getAttributeAsIntList
    fn info_ints(&self, key: &str) -> Vec<i32> {
        self.info.get(key).map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_default()
    }

    /// The record as htsjdk's VCFEncoder writes it.
    pub fn encode(&self, defs: &HeaderDefs) -> String {
        let mut line = String::new();
        line.push_str(&self.chromosome);
        line.push('\t');
        line.push_str(&self.position.to_string());
        line.push('\t');
        line.push_str(&self.id);
        line.push('\t');
        line.push_str(&self.ref_);
        line.push('\t');
        if self.alts.is_empty() {
            line.push_str(MISSING);
        } else {
            line.push_str(&self.alts.join(","));
        }
        line.push('\t');
        match self.log10_p_error {
            None => line.push_str(MISSING),
            Some(log10) => {
                let qual = java_fixed(log10 * -10.0 + 0.0, 2);
                line.push_str(qual.strip_suffix(".00").unwrap_or(&qual));
            }
        }
        line.push('\t');
        match &self.filters {
            None => line.push_str(MISSING),
            Some(f) if f.is_empty() => line.push_str("PASS"),
            Some(f) => {
                let mut sorted = f.clone();
                sorted.sort();
                line.push_str(&sorted.join(";"));
            }
        }
        line.push('\t');
        if self.info.is_empty() {
            line.push_str(MISSING);
        } else {
            let mut first = true;
            for (key, value) in &self.info {
                if !first {
                    line.push(';');
                }
                first = false;
                line.push_str(key);
                let is_flag = self.info_def_count(defs, key) == Some(FieldCount::Fixed(0));
                if !value.is_empty() && !is_flag {
                    line.push('=');
                    line.push_str(value);
                }
            }
        }
        if self.genotypes.is_empty() {
            return line;
        }

        // calcVCFGenotypeKeys
        let mut keys: BTreeSet<&str> = BTreeSet::new();
        let mut saw_gt = false;
        for g in &self.genotypes {
            keys.extend(g.attributes.keys().map(|k| k.as_str()));
            keys.extend(g.standard.keys().map(|k| k.as_str()));
            if g.filters.is_some() {
                keys.insert("FT");
            }
            saw_gt |= g.gt.is_some();
        }
        let mut format: Vec<&str> = Vec::new();
        if saw_gt || keys.is_empty() {
            format.push("GT");
        }
        format.extend(keys.iter());
        line.push('\t');
        line.push_str(&format.join(":"));

        let allele_count = 1 + self.alts.len() as i32;
        for g in &self.genotypes {
            line.push('\t');
            let has_gt = format[0] == "GT";
            if has_gt {
                line.push_str(g.gt.as_deref().unwrap_or("./."));
            }
            let mut values: Vec<String> = Vec::new();
            for key in format.iter().filter(|k| **k != "GT") {
                let value = if *key == "FT" {
                    g.filters.clone().unwrap_or_else(|| "PASS".to_string())
                } else if STANDARD_GENOTYPE_KEYS.contains(key) {
                    g.standard.get(*key).cloned().unwrap_or_else(|| MISSING.to_string())
                } else {
                    match g.attributes.get(*key) {
                        Some(v) => v.clone(),
                        None => {
                            // a missing field which expects several values is written with all of them
                            let expected = match defs.format.get(*key).map(|d| d.count) {
                                Some(FieldCount::Fixed(n)) => n,
                                Some(FieldCount::A) => allele_count - 1,
                                Some(FieldCount::R) => allele_count,
                                Some(FieldCount::G) => allele_count * (allele_count + 1) / 2,
                                Some(FieldCount::Unbounded) | None => -1,
                            };
                            if expected > 1 { vec![MISSING; expected as usize].join(",") } else { MISSING.to_string() }
                        }
                    }
                };
                values.push(value);
            }
            // trailing missing values are dropped; an empty string counts as missing
            while values.last().is_some_and(|v| is_missing_value(v)) {
                values.pop();
            }
            for (i, value) in values.iter().enumerate() {
                if i > 0 || has_gt {
                    line.push(':');
                }
                line.push_str(value);
            }
        }
        line
    }

    fn info_def_count(&self, defs: &HeaderDefs, key: &str) -> Option<FieldCount> {
        defs.info.get(key).map(|d| d.count)
    }
}

// --- CandidateSerialisation ---

/// ReadContextVcfInfo
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadContextVcfInfo {
    pub alignment_start: i32,
    pub var_index: i32,
    pub left_flank: String,
    pub core: String,
    pub right_flank: String,
    pub cigar: String,
}

impl ReadContextVcfInfo {
    pub fn from_vcf_tag(tag: &str) -> Option<ReadContextVcfInfo> {
        let values: Vec<&str> = tag.splitn(6, '-').collect();
        if values.len() != 6 {
            return None;
        }
        Some(ReadContextVcfInfo {
            alignment_start: values[0].parse().ok()?,
            var_index: values[1].parse().ok()?,
            left_flank: values[2].to_string(),
            core: values[3].to_string(),
            right_flank: values[4].to_string(),
            cigar: values[5].to_string(),
        })
    }

    pub fn read_bases(&self) -> String {
        format!("{}{}{}", self.left_flank, self.core, self.right_flank)
    }
}

/// CandidateSerialisation.toContext: the record fields which describe a candidate. The caller writes its
/// records directly (vcfout::record), so this is only the counterpart of `to_candidate` for the tests.
#[cfg_attr(not(test), allow(dead_code))]
pub fn to_context(chromosome: &str, candidate: &Candidate) -> VariantContext {
    let ctx = &candidate.ctx;
    let mut info: BTreeMap<String, String> = BTreeMap::new();
    info.insert("TIER".to_string(), candidate.tier.name().to_string());
    info.insert("RC_INFO".to_string(), vcfout::read_context_vcf_info(ctx));
    info.insert("RC_NM".to_string(), candidate.min_events.to_string());
    info.insert("TNC".to_string(), String::from_utf8_lossy(ctx.trinucleotide()).into_owned());
    if let Some(h) = &ctx.homology {
        let bases = String::from_utf8_lossy(&h.bases).into_owned();
        info.insert("RC_MH".to_string(), bases.clone());
        info.insert("MH".to_string(), bases);
    }
    if let Some(r) = &ctx.max_repeat {
        info.insert("RC_REPC".to_string(), r.count.to_string());
        info.insert("RC_REPS".to_string(), r.bases_str().to_string());
    }
    if let Some(r) = &ctx.max_ref_repeat {
        info.insert("REP_C".to_string(), r.count.to_string());
        info.insert("REP_S".to_string(), r.bases_str().to_string());
    }
    VariantContext {
        chromosome: chromosome.to_string(),
        position: candidate.position(),
        id: MISSING.to_string(),
        ref_: candidate.variant().ref_str().to_string(),
        alts: vec![candidate.variant().alt_str().to_string()],
        log10_p_error: None,
        filters: None,
        info,
        genotypes: Vec::new(),
    }
}

fn parse_cigar(cigar: &str) -> Result<Vec<CigarEl>> {
    let mut out = Vec::new();
    let mut len: i32 = 0;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            len = len * 10 + d as i32;
            continue;
        }
        let op = match c {
            'M' => Op::M,
            'I' => Op::I,
            'D' => Op::D,
            'N' => Op::N,
            'S' => Op::S,
            'H' => Op::H,
            'P' => Op::P,
            '=' => Op::Eq,
            'X' => Op::X,
            _ => bail!("invalid cigar {cigar}"),
        };
        out.push(CigarEl::new(len, op));
        len = 0;
    }
    Ok(out)
}

fn parse_tier(name: &str) -> Result<Tier> {
    Ok(match name {
        "HOTSPOT" => Tier::Hotspot,
        "PANEL" => Tier::Panel,
        "HIGH_CONFIDENCE" => Tier::HighConfidence,
        "LOW_CONFIDENCE" => Tier::LowConfidence,
        _ => bail!("invalid variant tier {name}"),
    })
}

/// CandidateSerialisation.toCandidate: rebuilds the read context from the RC_INFO tag by treating it as a
/// read. A context which cannot be rebuilt, or is not valid, is kept but marked invalid.
pub fn to_candidate(context: &VariantContext, ref_sequence: &RefSequence) -> Result<Candidate> {
    let alt = context.alts.first().ok_or_else(|| anyhow!("variant {}:{} has no alt allele", context.chromosome, context.position))?;
    let variant = SimpleVariant::new(context.position, context.ref_.as_bytes(), alt.as_bytes());
    let tier = match context.info.get("TIER").filter(|t| t.as_str() != MISSING) {
        Some(t) => parse_tier(t)?,
        None => Tier::LowConfidence,
    };
    let tag = context.info.get("RC_INFO").map(|s| s.as_str()).unwrap_or("");
    let vcf_info = ReadContextVcfInfo::from_vcf_tag(tag).ok_or_else(|| {
        anyhow!("variant {}:{} has no valid read context info: {tag}", context.chromosome, context.position)
    })?;
    let bases = vcf_info.read_bases().into_bytes();
    let cigar = parse_cigar(&vcf_info.cigar)?;
    let ref_length: i32 = cigar.iter().filter(|e| e.op.consumes_ref()).map(|e| e.len).sum();
    let read = Read {
        name: "CANDIDATE".to_string(),
        flags: 0,
        tid: -1,
        start: vcf_info.alignment_start,
        end: vcf_info.alignment_start + ref_length - 1,
        mapq: 0,
        cigar,
        // maximum across technologies
        quals: vec![40; bases.len()],
        bases,
        mate_tid: -1,
        mate_start: 0,
        insert_size: 0,
        nm: None,
        mate_cigar: None,
        consensus: ConsensusType::None,
        has_consensus_tag: false,
        duplex_index: -1,
        tp: Vec::new(),
        t0: Vec::new(),
        low_qual_tag: None,
        consensus_read_tag: None,
    };
    let read_context = create_context(&variant, &read, vcf_info.var_index, ref_sequence, DEFAULT_FLANK_LENGTH);
    let read_context = match read_context {
        Some(ctx) if ctx.is_valid() => ctx,
        other => {
            eprintln!(
                "variant({}:{} {}>{}) failed to recreate read context from VCF info({tag})",
                context.chromosome,
                variant.pos,
                variant.ref_str(),
                variant.alt_str()
            );
            let mut ctx = other.unwrap_or_else(|| VariantReadContext {
                variant: variant.clone(),
                alignment_start: vcf_info.alignment_start,
                alignment_end: 0,
                ref_bases: Vec::new(),
                read_bases: read.bases.clone(),
                core_index_start: 0,
                var_index: vcf_info.var_index,
                core_index_end: 0,
                homology: None,
                max_repeat: None,
                all_repeats: Vec::new(),
                core_position_start: 0,
                core_position_end: 0,
                read_cigar: String::new(),
                has_indel_in_core: false,
                artefact: None,
                max_ref_repeat: None,
                extended_ref_bases: None,
                marked_invalid: false,
            });
            ctx.marked_invalid = true;
            ctx
        }
    };
    let min_events = context.info.get("RC_NM").and_then(|v| v.parse().ok()).unwrap_or(0);
    Ok(Candidate { tier, ctx: Arc::new(read_context), min_events, full_match_support: 0, core_match_support: 0 })
}

// --- AppendVariantPhaser ---

/// phase/LpsReadCounts
#[derive(Clone, Debug)]
pub struct LpsReadCounts {
    pub id: i32,
    pub variants: Vec<SimpleVariant>,
    pub subset_variant_counts: Option<Vec<LpsReadCounts>>,
    pub depth: i32,
    pub alt_support: i32,
}

fn variants_match(own: &[SimpleVariant], variants: &[SimpleVariant]) -> bool {
    own.len() == variants.len() && own.iter().all(|x| variants.contains(x))
}

impl LpsReadCounts {
    fn new(id: i32, variants: Vec<SimpleVariant>) -> Self {
        LpsReadCounts { id, variants, subset_variant_counts: None, depth: 0, alt_support: 0 }
    }

    fn register_count(&mut self, variants: &[SimpleVariant], has_alt_support: bool) {
        let counts = if variants.len() < self.variants.len() {
            let subsets = self.subset_variant_counts.get_or_insert_with(Vec::new);
            let index = match subsets.iter().position(|x| variants_match(&x.variants, variants)) {
                Some(i) => i,
                None => {
                    subsets.push(LpsReadCounts::new(subsets.len() as i32, variants.to_vec()));
                    subsets.len() - 1
                }
            };
            &mut subsets[index]
        } else {
            self
        };
        counts.depth += 1;
        if has_alt_support {
            counts.alt_support += 1;
        }
    }
}

struct LocalPhaseSet {
    id: i32,
    variants: Vec<SimpleVariant>,
}

/// phase/AppendVariantPhaser: counts, per appended sample, the reads which support or cover the variants of
/// each existing local phase set.
#[derive(Default)]
pub struct AppendVariantPhaser {
    /// in the iteration order of Java's HashMap<Integer, LocalPhaseSet>
    local_phase_sets: Vec<LocalPhaseSet>,
    variant_lps_ids: HashMap<SimpleVariant, Vec<i32>>,
    sample_lps_read_counts: Vec<(String, HashMap<i32, LpsReadCounts>)>,
    current_sample: usize,
}

impl AppendVariantPhaser {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn initialise(&mut self, sample: &str) {
        match self.sample_lps_read_counts.iter().position(|(s, _)| s == sample) {
            Some(i) => {
                self.sample_lps_read_counts[i].1 = HashMap::new();
                self.current_sample = i;
            }
            None => {
                self.sample_lps_read_counts.push((sample.to_string(), HashMap::new()));
                self.current_sample = self.sample_lps_read_counts.len() - 1;
            }
        }
    }

    pub fn sample_lps_read_counts(&self, sample: &str) -> Option<&HashMap<i32, LpsReadCounts>> {
        self.sample_lps_read_counts.iter().find(|(s, _)| s == sample).map(|(_, counts)| counts)
    }

    pub fn register_local_phase_sets(&mut self, candidates: &[Candidate], variant_contexts: &[VariantContext]) {
        // in order of first insertion, which with the key values decides the map's iteration order
        let mut sets: Vec<LocalPhaseSet> = Vec::new();
        for (candidate, context) in candidates.iter().zip(variant_contexts) {
            if context.is_filtered() || !context.info.contains_key("LPS") {
                continue;
            }
            for lps_id in context.info_ints("LPS") {
                match sets.iter_mut().find(|s| s.id == lps_id) {
                    Some(set) => set.variants.push(candidate.variant().clone()),
                    None => sets.push(LocalPhaseSet { id: lps_id, variants: vec![candidate.variant().clone()] }),
                }
            }
        }
        let inserted: Vec<i32> = sets.iter().map(|s| s.id).collect();
        let order = crate::phase::java_int_hash_order(&inserted);
        // purge any single-variant groups; the map keeps its capacity, so the order of the rest is unchanged
        let mut by_id: HashMap<i32, LocalPhaseSet> = sets.into_iter().filter(|s| s.variants.len() != 1).map(|s| (s.id, s)).collect();
        self.local_phase_sets = order.into_iter().filter_map(|id| by_id.remove(&id)).collect();

        // make links from variants to their LPS IDs
        for set in &self.local_phase_sets {
            for variant in &set.variants {
                self.variant_lps_ids.entry(variant.clone()).or_default().push(set.id);
            }
        }
    }

    /// registeredPhasedVariants(posCounters, negCounters)
    pub fn registered_phased_variants(&mut self, pos_counters: &[&ReadContextCounter], neg_counters: &[&ReadContextCounter]) {
        if pos_counters.len() >= 2 {
            let variants: Vec<SimpleVariant> = pos_counters.iter().map(|c| c.variant().clone()).collect();
            self.register_variants(&variants, true);
        }
        if neg_counters.len() >= 2 {
            let variants: Vec<SimpleVariant> = neg_counters.iter().map(|c| c.variant().clone()).collect();
            self.register_variants(&variants, false);
        }
    }

    fn variants_match_lps(&self, variants: &[SimpleVariant]) -> bool {
        self.local_phase_sets.iter().any(|s| variants_match(&s.variants, variants))
    }

    fn register_variants(&mut self, variants: &[SimpleVariant], has_alt_support: bool) {
        let mut processed_lps_ids: Vec<i32> = Vec::new();
        for variant in variants {
            let Some(lps_ids) = self.variant_lps_ids.get(variant).cloned() else { continue };
            for lps_id in lps_ids {
                if processed_lps_ids.contains(&lps_id) {
                    continue;
                }
                processed_lps_ids.push(lps_id);
                let Some(set_index) = self.local_phase_sets.iter().position(|s| s.id == lps_id) else {
                    eprintln!("lpsId({lps_id}) missing");
                    continue;
                };
                // skip if the variants from the read aren't an exact match for this LPS and are for another
                if !variants_match(&self.local_phase_sets[set_index].variants, variants) && self.variants_match_lps(variants) {
                    continue;
                }
                let set_variants = &self.local_phase_sets[set_index].variants;
                let mut lps_variants: Vec<SimpleVariant> = variants.to_vec();
                if variants.iter().any(|x| !set_variants.contains(x)) {
                    // remove any variants not in the LPS
                    lps_variants.retain(|x| set_variants.contains(x));
                    if lps_variants.len() < 2 {
                        continue;
                    }
                }
                let set_variants = set_variants.clone();
                let count_variants =
                    self.sample_lps_read_counts[self.current_sample].1.get(&lps_id).map_or(set_variants.len(), |c| c.variants.len());
                // the counts are created before the subset test, as in Java
                self.sample_lps_read_counts[self.current_sample]
                    .1
                    .entry(lps_id)
                    .or_insert_with(|| LpsReadCounts::new(lps_id, set_variants));
                // only look for subsets if this variant set doesn't match any others precisely
                if lps_variants.len() < count_variants && self.variants_match_lps(&lps_variants) {
                    continue;
                }
                self.sample_lps_read_counts[self.current_sample]
                    .1
                    .get_mut(&lps_id)
                    .unwrap()
                    .register_count(&lps_variants, has_alt_support);
            }
        }
    }

    /// populateLocalPhaseSetInfo: sets the LPSA genotype field of the appended samples
    pub fn populate_local_phase_set_info(&self, candidates: &[Candidate], variant_contexts: &mut [VariantContext]) {
        for (candidate, context) in candidates.iter().zip(variant_contexts.iter_mut()) {
            let Some(lps_ids) = self.variant_lps_ids.get(candidate.variant()) else { continue };
            for genotype in context.genotypes.iter_mut() {
                let Some(sample_read_counts) = self.sample_lps_read_counts(&genotype.sample) else { continue };
                let mut lps_info: Vec<String> = Vec::new();
                for lps_id in lps_ids {
                    let Some(read_counts) = sample_read_counts.get(lps_id) else { continue };
                    lps_info.push(format!("{}={}/{}", read_counts.id, read_counts.alt_support, read_counts.depth));
                    if let Some(subsets) = &read_counts.subset_variant_counts {
                        for subset in subsets.iter().filter(|s| s.variants.contains(candidate.variant())) {
                            lps_info.push(format!("{}_{}={}/{}", read_counts.id, subset.id, subset.alt_support, subset.depth));
                        }
                    }
                }
                genotype.attributes.insert(LPS_APPEND_INFO.to_string(), lps_info.join(";"));
            }
        }
    }
}

impl PhaseSink for AppendVariantPhaser {
    fn register(&mut self, counters: &[ReadContextCounter], positive: &[usize], negative: &[usize]) {
        let pos: Vec<&ReadContextCounter> = positive.iter().map(|&i| &counters[i]).collect();
        let neg: Vec<&ReadContextCounter> = negative.iter().map(|&i| &counters[i]).collect();
        self.registered_phased_variants(&pos, &neg);
    }
}

// --- RegionAppendTask ---

/// VariantContextFactory.createGenotype(counter, sampleId, existingTags)
fn create_genotype(counter: &ReadContextCounter, sample: &str, existing_tags: &BTreeSet<String>) -> Genotype {
    let mut g = Genotype { sample: sample.to_string(), ..Default::default() };
    for (key, value) in vcfout::genotype_fields(counter) {
        match key {
            "GT" => g.gt = Some(value),
            "AD" | "DP" => {
                g.standard.insert(key.to_string(), value);
            }
            _ => {
                g.attributes.insert(key.to_string(), value);
            }
        }
    }
    for tag in existing_tags {
        if !g.attributes.contains_key(tag) {
            g.attributes.insert(tag.clone(), String::new());
        }
    }
    g
}

/// RegionAppendTask.checkGenotypeFields: an existing genotype gets a default for each field it lacks
fn check_genotype_fields(genotype: &mut Genotype, expected_tags: &BTreeSet<String>) {
    for tag in expected_tags {
        if genotype.attributes.contains_key(tag) {
            continue;
        }
        let value = match tag.as_str() {
            "MUC" => "0".to_string(),
            "ASBQ" => "0,0".to_string(),
            "AED" => format!("{},{}", vcf_double(0.0), vcf_double(0.0)),
            "UMI_CNT" => "0,0,0,0,0,0".to_string(),
            _ => String::new(),
        };
        genotype.attributes.insert(tag.clone(), value);
    }
}

/// RegionAppendTask.addGenotype
fn add_genotype(context: &VariantContext, counters: &[&ReadContextCounter], sample_ids: &[String]) -> VariantContext {
    let mut existing_tags: BTreeSet<String> = BTreeSet::new();
    for genotype in &context.genotypes {
        existing_tags.extend(genotype.attributes.keys().cloned());
    }
    let new_genotypes: Vec<Genotype> =
        counters.iter().zip(sample_ids).map(|(counter, sample)| create_genotype(counter, sample, &existing_tags)).collect();
    let expected_tags: BTreeSet<String> =
        new_genotypes.first().map(|g| g.attributes.keys().cloned().collect()).unwrap_or_default();
    let mut updated = context.clone();
    for genotype in updated.genotypes.iter_mut() {
        check_genotype_fields(genotype, &expected_tags);
    }
    updated.genotypes.extend(new_genotypes);
    updated
}

struct Shared<'a> {
    cfg: &'a Config,
    ref_data: &'a RefData,
    bqr: &'a [BqrMap],
    msi: &'a [MsiJitter],
    defs: &'a HeaderDefs,
}

struct Worker {
    genome: RefGenome,
    bams: Vec<BamSlicer>,
}

impl Worker {
    fn open(cfg: &Config) -> Result<Worker> {
        let bams: Result<Vec<BamSlicer>> = cfg.reference_bams.iter().map(|p| BamSlicer::open(p, &cfg.ref_genome, 1)).collect();
        Ok(Worker { genome: RefGenome::open(&cfg.ref_genome)?, bams: bams? })
    }
}

/// RegionAppendTask.call: the region's records with the genotypes of the appended samples added
fn process_region(
    sh: &Shared, w: &mut Worker, chr: &str, region: BaseRegion, variants: &[&VariantContext],
) -> Result<(Vec<String>, Vec<String>)> {
    let cfg = sh.cfg;
    let ref_sequence = RefSequence::new(
        &w.genome,
        chr,
        1.max(region.start - READ_CONTEXT_EXTENSION_LENGTH),
        region.end + READ_CONTEXT_EXTENSION_LENGTH,
    )?;
    let original: Vec<VariantContext> = variants.iter().map(|v| (*v).clone()).collect();
    let candidates: Result<Vec<Candidate>> = original.iter().map(|v| to_candidate(v, &ref_sequence)).collect();
    let candidates = candidates?;

    let mut phaser = AppendVariantPhaser::new();
    phaser.register_local_phase_sets(&candidates, &original);

    // one FragmentSync per region task, shared by the samples
    let mut sync = FragmentSync::new();
    let mut sample_counters: Vec<Vec<ReadContextCounter>> = Vec::new();
    for (i, sample_id) in cfg.reference_ids.iter().enumerate() {
        phaser.initialise(sample_id);
        let mut sample = evidence_stage::Sample { id: sample_id, bam: &mut w.bams[i], bqr: &sh.bqr[i], msi: &sh.msi[i], is_reference: true };
        sample_counters.push(evidence_stage::collect_evidence(
            cfg, chr, &candidates, &w.genome, sh.ref_data, &mut sample, &mut sync, Some(&mut phaser),
        )?);
    }

    let mut final_variants: Vec<VariantContext> = original
        .iter()
        .enumerate()
        .map(|(i, context)| {
            let counters: Vec<&ReadContextCounter> = sample_counters.iter().map(|c| &c[i]).collect();
            add_genotype(context, &counters, &cfg.reference_ids)
        })
        .collect();
    phaser.populate_local_phase_set_info(&candidates, &mut final_variants);
    let mut frag_lines = Vec::new();
    if cfg.write_frag_lengths {
        for (i, candidate) in candidates.iter().enumerate() {
            for (sample_id, counters) in cfg.reference_ids.iter().zip(&sample_counters) {
                crate::fragment_length_lines(sample_id, chr, candidate.variant(), &counters[i], &mut frag_lines);
            }
        }
    }
    Ok((final_variants.iter().map(|v| v.encode(sh.defs)).collect(), frag_lines))
}

// --- SageAppendApplication ---

/// The numeric components of a version, compared as java.lang.module.ModuleDescriptor.Version does for
/// plain dotted versions.
fn version_numbers(version: &str) -> Vec<u64> {
    version.split(|c: char| !c.is_ascii_digit()).filter(|x| !x.is_empty()).map(|x| x.parse().unwrap_or(0)).collect()
}

struct InputVcf {
    meta_lines: Vec<String>,
    samples: Vec<String>,
    records: Vec<String>,
}

fn read_vcf(path: &str) -> Result<InputVcf> {
    let file = std::fs::File::open(path).with_context(|| format!("opening input VCF {path}"))?;
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") || path.ends_with(".bgz") {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };
    let mut vcf = InputVcf { meta_lines: Vec::new(), samples: Vec::new(), records: Vec::new() };
    let mut has_column_line = false;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("##") {
            vcf.meta_lines.push(line);
        } else if line.starts_with('#') {
            vcf.samples = line.split('\t').skip(9).map(|s| s.to_string()).collect();
            has_column_line = true;
        } else if !line.is_empty() {
            vcf.records.push(line);
        }
    }
    if !has_column_line {
        bail!("invalid input VCF({path})");
    }
    Ok(vcf)
}

/// SageCommon.setReadLength with no panel: the specific regions, or else a fixed region of the first
/// reference sequence
fn sample_read_length(cfg: &Config, genome: &RefGenome, bam: &mut BamSlicer) -> Result<i32> {
    let mut regions: Vec<(String, BaseRegion)> =
        cfg.specific_regions.iter().map(|(c, s, e)| (c.clone(), BaseRegion::new(*s, *e))).collect();
    if regions.is_empty() {
        if let Some(first) = genome.sequence_names()?.into_iter().next() {
            regions.push((first, BaseRegion::new(100_000, 200_000)));
        }
    }
    for (chr, region) in regions {
        let mut count = 0;
        let mut max_length = 0;
        let sliced = bam.slice_ref(&chr, &[BaseRegion::new(region.start.max(1), region.end)], false, |read| {
            if count < 1000 {
                count += 1;
                max_length = max_length.max(read.len());
            }
        });
        if sliced.is_ok() && count > 0 && max_length > 0 {
            return Ok(max_length);
        }
    }
    Ok(DEFAULT_READ_LENGTH)
}

struct VcfOutput {
    path: String,
    out: rust_htslib::bgzf::Writer,
}

impl VcfOutput {
    /// VariantVCF(reference, additionalSampleIds, existingHeader, outputVcf)
    fn create(path: &str, meta_lines: Vec<String>, samples: &[String]) -> Result<VcfOutput> {
        let mut out = rust_htslib::bgzf::Writer::from_path(path).with_context(|| format!("creating {path}"))?;
        let mut lines = vec!["##fileformat=VCFv4.2".to_string()];
        lines.extend(meta_lines.into_iter().filter(|l| !l.starts_with("##fileformat=")));
        let mut header = vcfout::sorted_header_lines(lines).join("\n");
        header.push_str("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !samples.is_empty() {
            header.push_str("\tFORMAT");
            for sample in samples {
                header.push('\t');
                header.push_str(sample);
            }
        }
        header.push('\n');
        out.write_all(header.as_bytes())?;
        Ok(VcfOutput { path: path.to_string(), out })
    }

    fn write(&mut self, line: &str) -> Result<()> {
        self.out.write_all(line.as_bytes())?;
        self.out.write_all(b"\n")?;
        Ok(())
    }

    fn close(mut self) -> Result<()> {
        self.out.flush()?;
        let path = self.path.clone();
        drop(self.out);
        vcfout::index_vcf(&path)
    }
}

/// `sage-rs append <args>`: SageAppendApplication.main
pub fn run(args: &[String]) -> Result<()> {
    set_append_mode(true);
    let mut cfg = Config::from_append_args(args)?;
    let input_vcf = cfg.input_vcf.clone().expect("input_vcf");

    // SageConfig.isValid
    if cfg.reference_ids.len() != cfg.reference_bams.len() {
        bail!("each reference sample must have matching BAM file");
    }
    for bam in &cfg.reference_bams {
        if !std::path::Path::new(bam).exists() {
            bail!("unable to locate reference bam({bam})");
        }
    }
    if input_vcf == cfg.output_vcf {
        bail!("input and output VCFs must be different");
    }
    if cfg.reference_ids.is_empty() {
        bail!("missing reference Id must be supplied");
    }
    eprintln!("appending variants with {} sample BAM(s)", cfg.reference_bams.len());
    eprintln!("reading and validating file: {input_vcf}");

    let input = read_vcf(&input_vcf)?;

    // validateInputHeader
    let version = input.meta_lines.iter().find_map(|l| l.strip_prefix("##sageVersion=")).unwrap_or("0");
    if version_numbers(version) < version_numbers(MIN_PRIOR_VERSION) {
        bail!("Sage VCF version({version}) older than required({MIN_PRIOR_VERSION})");
    }
    eprintln!("existing VCF samples: {}", input.samples.join(", "));
    for sample in &cfg.reference_ids {
        if input.samples.contains(sample) {
            bail!("config reference sample({sample}) already exits in input VCF");
        }
    }

    // the header as htsjdk holds it after reading
    let mut meta_lines: Vec<String> = input.meta_lines.iter().map(|l| repair_standard_header_line(l)).collect();
    // addGenotypeHeader, called again in case this version has new genotype fields, and addAppendHeader.
    // VCFHeader.addMetaDataLine does not add a line when the header already has one of the same type and
    // id, whatever its definition, so the fields of an older VCF keep their header lines.
    let mut new_lines: Vec<&str> = vcfout::header_field_lines()
        .filter(|l| l.starts_with("##FORMAT="))
        .filter(|l| !["GT", "AD", "DP"].iter().any(|id| l.starts_with(&format!("##FORMAT=<ID={id},"))))
        .collect();
    new_lines.push(LPS_APPEND_HEADER);
    for line in new_lines {
        let id = parse_field_line(line).map(|f| f.1);
        let exists = meta_lines.iter().any(|l| l.starts_with("##FORMAT=") && parse_field_line(l).map(|f| f.1) == id);
        if !exists {
            meta_lines.push(line.to_string());
        }
    }
    let mut defs = HeaderDefs::default();
    for line in &meta_lines {
        defs.add_line(line);
    }

    let mut existing_variants: Vec<VariantContext> = Vec::new();
    for line in &input.records {
        let variant = VariantContext::parse(line, &defs, &input.samples)?;
        if cfg.require_gene {
            // VariantImpactSerialiser: the gene is the first of the ten IMPACT values
            let impact: Vec<&str> = variant.info.get("IMPACT").map(|v| java_split(v, ',')).unwrap_or_default();
            if impact.len() != 10 || impact[0].is_empty() {
                continue;
            }
        }
        if !cfg.specific_variants.is_empty() {
            // SimpleVariant.matches(chromosome, position): the alleles are not compared
            if !cfg.specific_variants.iter().any(|(c, v)| *c == variant.chromosome && v.pos == variant.position) {
                continue;
            }
        } else if !cfg.specific_regions.is_empty() {
            let pos = variant.position;
            if !cfg.specific_regions.iter().any(|(c, s, e)| *c == variant.chromosome && pos >= *s && pos <= *e) {
                continue;
            }
        } else if !cfg.specific_chr.is_empty() && !cfg.specific_chr.contains(&variant.chromosome) {
            continue;
        }
        existing_variants.push(variant);
    }
    eprintln!("loaded {} variants", existing_variants.len());
    eprintln!("writing to file: {}", cfg.output_vcf);

    let mut samples = input.samples.clone();
    samples.extend(cfg.reference_ids.iter().cloned());
    let mut output = VcfOutput::create(&cfg.output_vcf, meta_lines, &samples)?;
    if existing_variants.is_empty() {
        eprintln!("writing empty output VCF");
        return output.close();
    }

    let mut main_worker = Worker::open(&cfg)?;
    if cfg.read_length == 0 {
        cfg.read_length = sample_read_length(&cfg, &main_worker.genome, &mut main_worker.bams[0])?;
        if cfg.read_length != DEFAULT_READ_LENGTH {
            eprintln!("max observed read length set({})", cfg.read_length);
        }
    }
    let cfg = cfg;

    let mut bqr = Vec::new();
    let mut msi = Vec::new();
    for sample in &cfg.reference_ids {
        bqr.push(if cfg.skip_bqr {
            BqrMap::default()
        } else {
            BqrMap::load(&format!("{}/{}.redux.bqr.tsv", cfg.jitter_bqr_dir, sample))?
        });
        let dir = if cfg.skip_msi_jitter { None } else { Some(cfg.jitter_bqr_dir.as_str()) };
        msi.push(MsiJitter::load(dir, sample, cfg.high_depth_mode)?);
    }
    if is_ultima() {
        crate::ultima::set_max_raw_qual(crate::ultima_max_raw_qual(&cfg, &bqr));
    }
    // append does not load the gene panel, hotspots or Ensembl data, so no position is in an HLA gene
    let ref_data = RefData::load(&crate::refdata::RefDataPaths {
        ensembl_dir: None,
        driver_gene_panel: None,
        panel_bed: None,
        hotspots: None,
        high_confidence_bed: None,
        v38: cfg.v38,
    })?;

    // every region task with variants, in the order the BAM's sequence dictionary lists the chromosomes
    let mut tasks: Vec<(String, BaseRegion, Vec<&VariantContext>)> = Vec::new();
    for chr in main_worker.bams[0].sequence_names() {
        if !cfg.process_chromosome(&chr) {
            continue;
        }
        let chromosome_variants: Vec<&VariantContext> = existing_variants.iter().filter(|v| v.chromosome == chr).collect();
        let chr_length = main_worker.genome.chromosome_length(&chr);
        for region in crate::partition(&cfg, &chr, chr_length) {
            let region_variants: Vec<&VariantContext> =
                chromosome_variants.iter().filter(|v| region.contains(v.position)).copied().collect();
            if !region_variants.is_empty() {
                tasks.push((chr.clone(), region, region_variants));
            }
        }
    }

    let shared = Shared { cfg: &cfg, ref_data: &ref_data, bqr: &bqr, msi: &msi, defs: &defs };
    let mut frag_out = crate::create_fragment_length_writer(&cfg)?;
    let next = AtomicUsize::new(0);
    let results: Vec<Mutex<Option<Result<(Vec<String>, Vec<String>)>>>> = tasks.iter().map(|_| Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..cfg.threads.max(1).min(tasks.len().max(1)) {
            scope.spawn(|| {
                let mut worker = Worker::open(&cfg);
                loop {
                    let index = next.fetch_add(1, Ordering::SeqCst);
                    if index >= tasks.len() {
                        return;
                    }
                    let (chr, region, variants) = &tasks[index];
                    let result = match worker.as_mut() {
                        Ok(w) => process_region(&shared, w, chr, *region, variants),
                        Err(e) => Err(anyhow!("failed to open inputs: {e}")),
                    };
                    *results[index].lock().unwrap() = Some(result);
                }
            });
        }
    });
    for result in results {
        let (lines, frag_lines) = result.into_inner().unwrap().ok_or_else(|| anyhow!("region task did not complete"))??;
        for line in lines {
            output.write(&line)?;
        }
        if let Some(w) = frag_out.as_mut() {
            for line in frag_lines {
                writeln!(w, "{line}")?;
            }
        }
    }
    output.close()?;
    if let Some(w) = frag_out {
        w.finish()?.flush()?;
    }
    eprintln!("SageAppend complete");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn defs() -> HeaderDefs {
        let mut defs = HeaderDefs::default();
        for line in vcfout::header_field_lines() {
            defs.add_line(line);
        }
        defs.add_line(LPS_APPEND_HEADER);
        defs
    }

    // values taken from a Java SAGE VCF before and after Java's append: floats written with "%.3f" by
    // SAGE (FSB, RSB) come back in htsjdk's own number format
    #[test]
    fn record_round_trip() {
        let line = "chr21\t5218922\t.\tTGTG\tT\t149\tmaxEdgeDistance;minMQF\tMED=0.257;MH=GTG;MQF=-5.513e+01;NEARBY_INDEL;RC_NM=1;TQP=1.221e-15\t\
                    GT:AD:AED:AF:DP:FSB:RC_JIT:RSB\t./.:123,3:0.223,0.067:0.024:126:0.545,0.000:0,0:1.000,0.005";
        let samples = vec!["S1".to_string()];
        let variant = VariantContext::parse(line, &defs(), &samples).unwrap();
        assert_eq!(
            variant.encode(&defs()),
            "chr21\t5218922\t.\tTGTG\tT\t149\tmaxEdgeDistance;minMQF\tMED=0.257;MH=GTG;MQF=-5.513e+01;NEARBY_INDEL;RC_NM=1;TQP=1.221e-15\t\
             GT:AD:AED:AF:DP:FSB:RC_JIT:RSB\t./.:123,3:0.223,0.067:0.024:126:0.545,0.00:0,0:1.00,5.000e-03"
        );
    }

    // a field which only some samples have is "." for the others (".,." where two values are declared),
    // and trailing missing or empty values are dropped
    #[test]
    fn missing_genotype_values() {
        let line = "1\t100\t.\tA\tC\t.\tPASS\t.\tGT:DP:LPSA:RSB\t./.:5:.:0.5,0.5\t./.:7:3=1/2";
        let samples = vec!["S1".to_string(), "S2".to_string()];
        let mut variant = VariantContext::parse(line, &defs(), &samples).unwrap();
        assert_eq!(variant.encode(&defs()), "1\t100\t.\tA\tC\t.\tPASS\t.\tGT:DP:LPSA:RSB\t./.:5:.:0.500,0.500\t./.:7:3=1/2");
        variant.genotypes[0].attributes.insert("ZZ".to_string(), String::new());
        variant.genotypes[1].attributes.remove("LPSA");
        assert_eq!(variant.encode(&defs()), "1\t100\t.\tA\tC\t.\tPASS\t.\tGT:DP:RSB:ZZ\t./.:5:0.500,0.500\t./.:7");
    }

    #[test]
    fn version_check() {
        assert!(version_numbers("5.1") >= version_numbers(MIN_PRIOR_VERSION));
        assert!(version_numbers("2.8") >= version_numbers(MIN_PRIOR_VERSION));
        assert!(version_numbers("2.7") < version_numbers(MIN_PRIOR_VERSION));
        assert!(version_numbers("0") < version_numbers(MIN_PRIOR_VERSION));
    }
}
