//! VCF files as htsjdk 2.23 reads and writes them: the header (VCFCodec with
//! VCFStandardHeaderLines.repairStandardHeaderLines, VCFHeader's sorted line order), records
//! (AbstractVCFCodec.decode and VariantContext.getAttributeAs*), and the writer (VCFEncoder, the BGZF
//! stream of BlockCompressedOutputStream and the tabix index TabixIndexCreator builds on the fly).
//!
//! A record keeps its fields as text; the encoder follows VCFEncoder: INFO keys sorted, flags without a
//! value, FILTER sorted, FORMAT keys sorted after GT, trailing missing genotype fields dropped.
//!
//! The record, header and encoder were esvee-rs's `caller::vcf`; the BGZF writer and the tabix index
//! amber-rs's `vcf`. The BGZF blocks hold 65498 bytes each and are deflated by zlib at level 5, as
//! htsjdk does, so the compressed bytes and therefore the index (whose chunks are virtual file offsets)
//! are those Java writes, as long as the zlib in use deflates as Java's does.

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::io::Write;

use crate::gz::LineReader;
use crate::java_text;

/// VariantContext.NO_LOG10_PERROR
pub const NO_LOG10_PERROR: f64 = 1.0;

// ---------------------------------------------------------------------------------------------
// header

#[derive(Clone, Debug, PartialEq)]
pub enum Count {
    Fixed(i32),
    A,
    R,
    G,
    Unbounded,
}

/// The Number and Type of an INFO or FORMAT key
#[derive(Clone, Debug, PartialEq)]
pub struct FieldDef {
    pub count: Count,
    pub ty: String,
}

#[derive(Clone, Debug, Default)]
pub struct VcfHeader {
    /// the `##` lines in input order and then those added, without the file format line
    pub meta_lines: Vec<String>,
    pub samples: Vec<String>,
    info: HashMap<String, FieldDef>,
    format: HashMap<String, FieldDef>,
    filters: std::collections::HashSet<String>,
}

/// (ID, Number, Type) of a `##INFO` or `##FORMAT` line
fn compound_fields(line: &str) -> Option<(String, String, String)> {
    let body = line.split_once("=<")?.1;
    let mut id = None;
    let mut number = None;
    let mut ty = None;
    // the three fields come before the description, which is the only quoted one
    let head = body.split("Description=").next().unwrap_or(body);
    for part in head.split(',') {
        if let Some((k, v)) = part.split_once('=') {
            match k {
                "ID" => id = Some(v.trim_end_matches('>').to_string()),
                "Number" => number = Some(v.trim_end_matches('>').to_string()),
                "Type" => ty = Some(v.trim_end_matches('>').to_string()),
                _ => {}
            }
        }
    }
    Some((id?, number?, ty?))
}

/// The ID of a `##FILTER` line
fn filter_id(line: &str) -> Option<String> {
    let body = line.strip_prefix("##FILTER=<")?;
    let head = body.split("Description=").next().unwrap_or(body);
    head.split(',').find_map(|p| p.strip_prefix("ID=")).map(|v| v.trim_end_matches('>').to_string())
}

fn parse_count(number: &str) -> Count {
    match number {
        "A" => Count::A,
        "R" => Count::R,
        "G" => Count::G,
        "." => Count::Unbounded,
        n => Count::Fixed(n.parse().unwrap_or(1)),
    }
}

/// VCFStandardHeaderLines: (ID, Number, Type, Description)
const STANDARD_FORMAT: &[(&str, &str, &str, &str)] = &[
    ("GT", "1", "String", "Genotype"),
    ("GQ", "1", "Integer", "Genotype Quality"),
    ("DP", "1", "Integer", "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"),
    ("PL", "G", "Integer", "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"),
    ("AD", "R", "Integer", "Allelic depths for the ref and alt alleles in the order listed"),
    ("FT", ".", "String", "Genotype-level filter"),
    ("PS", "1", "Integer", "Phasing set (typically the position of the first variant in the set)"),
    ("PQ", "1", "Float", "Read-backed phasing quality"),
];
const STANDARD_INFO: &[(&str, &str, &str, &str)] = &[
    ("END", "1", "Integer", "Stop position of the interval"),
    ("DB", "0", "Flag", "dbSNP Membership"),
    ("DP", "1", "Integer", "Approximate read depth; some reads may have been filtered"),
    ("SB", "1", "Float", "Strand Bias"),
    ("AF", "A", "Float", "Allele Frequency, for each ALT allele, in the same order as listed"),
    ("AC", "A", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed"),
    ("AN", "1", "Integer", "Total number of alleles in called genotypes"),
    ("MQ0", "1", "Integer", "Total Mapping Quality Zero Reads"),
    ("MQ", "1", "Float", "RMS Mapping Quality"),
    ("SOMATIC", "0", "Flag", "Somatic event"),
];

/// `##INFO=<...>` as VCFInfoHeaderLine.toString writes it
pub fn info_line(id: &str, number: &str, ty: &str, desc: &str) -> String {
    format!("##INFO=<ID={id},Number={number},Type={ty},Description=\"{desc}\">")
}
/// `##FORMAT=<...>`
pub fn format_line(id: &str, number: &str, ty: &str, desc: &str) -> String {
    format!("##FORMAT=<ID={id},Number={number},Type={ty},Description=\"{desc}\">")
}
/// `##FILTER=<...>`
pub fn filter_line(id: &str, desc: &str) -> String {
    format!("##FILTER=<ID={id},Description=\"{desc}\">")
}

/// VCFStandardHeaderLines.repairStandardHeaderLines: a standard key declared with another count or
/// type is replaced by the standard line
fn repair_standard_line(line: &str) -> String {
    let (kind, table): (&str, &[(&str, &str, &str, &str)]) = if line.starts_with("##INFO=<") {
        ("INFO", STANDARD_INFO)
    } else if line.starts_with("##FORMAT=<") {
        ("FORMAT", STANDARD_FORMAT)
    } else {
        return line.to_string();
    };
    let Some((id, number, ty)) = compound_fields(line) else { return line.to_string() };
    let Some(standard) = table.iter().find(|s| s.0 == id) else { return line.to_string() };
    if parse_count(&number) != parse_count(standard.1) || ty != standard.2 {
        if kind == "INFO" {
            info_line(standard.0, standard.1, standard.2, standard.3)
        } else {
            format_line(standard.0, standard.1, standard.2, standard.3)
        }
    } else {
        line.to_string()
    }
}

impl VcfHeader {
    pub fn from_lines(meta: Vec<String>, samples: Vec<String>) -> VcfHeader {
        let mut h = VcfHeader { samples, ..Default::default() };
        for l in meta.iter().filter(|l| !l.starts_with("##fileformat=")) {
            h.push_line(repair_standard_line(l));
        }
        h
    }

    fn push_line(&mut self, line: String) {
        // the definitions as the header's maps hold them: the last line of a key
        if line.starts_with("##INFO=<") || line.starts_with("##FORMAT=<") {
            if let Some((id, number, ty)) = compound_fields(&line) {
                let def = FieldDef { count: parse_count(&number), ty };
                if line.starts_with("##INFO=<") {
                    self.info.insert(id, def);
                } else {
                    self.format.insert(id, def);
                }
            }
        } else if let Some(id) = filter_id(&line) {
            self.filters.insert(id);
        }
        self.meta_lines.push(line);
    }

    /// VCFHeader.addMetaDataLine: a line equal to one present is not added again
    pub fn add_meta_line(&mut self, line: String) {
        if !self.meta_lines.contains(&line) {
            self.push_line(line);
        }
    }

    pub fn has_info_line(&self, id: &str) -> bool {
        self.info.contains_key(id)
    }
    pub fn has_format_line(&self, id: &str) -> bool {
        self.format.contains_key(id)
    }
    pub fn has_filter_line(&self, id: &str) -> bool {
        self.filters.contains(id)
    }
    pub fn info_def(&self, id: &str) -> Option<&FieldDef> {
        self.info.get(id)
    }
    pub fn format_def(&self, id: &str) -> Option<&FieldDef> {
        self.format.get(id)
    }

    /// The value of a `##key=value` line (VCFHeader.getMetaDataLine(key).getValue())
    pub fn meta_value(&self, key: &str) -> Option<&str> {
        let prefix = format!("##{key}=");
        self.meta_lines.iter().find_map(|l| l.strip_prefix(prefix.as_str()))
    }

    /// The header text as VCFWriter writes it: the file format line, the other lines sorted as strings
    /// without duplicates (the contig lines kept together in their own order), the column line.
    pub fn text(&self) -> String {
        let mut text = sorted_header_lines(&self.meta_lines).join("\n");
        text.push_str("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            text.push_str("\tFORMAT");
            for sample in &self.samples {
                text.push('\t');
                text.push_str(sample);
            }
        }
        text.push('\n');
        text
    }
}

/// Header lines in the order htsjdk writes them: the file format first, then the rest sorted as
/// strings without duplicates, with the contig lines kept together in their own order.
pub fn sorted_header_lines(lines: &[String]) -> Vec<String> {
    let mut contigs: Vec<&String> = Vec::new();
    let mut others: Vec<&String> = Vec::new();
    for line in lines {
        if line.starts_with("##fileformat=") {
            continue;
        } else if line.starts_with("##contig=") {
            if !contigs.contains(&line) {
                contigs.push(line);
            }
        } else {
            others.push(line);
        }
    }
    others.sort();
    others.dedup();
    let split = others.partition_point(|x| x.as_str() < "##contig=");
    let mut out: Vec<String> = vec!["##fileformat=VCFv4.2".to_string()];
    out.extend(others[..split].iter().map(|s| s.to_string()));
    out.extend(contigs.iter().map(|s| s.to_string()));
    out.extend(others[split..].iter().map(|s| s.to_string()));
    out
}

/// Reads the header of a VCF (plain, gzip or BGZF): the `##` lines and the samples of the column line.
/// The reader is left at the first record.
pub fn read_header(path: &str) -> Result<(VcfHeader, LineReader, Option<String>)> {
    let mut reader = LineReader::open(path).with_context(|| format!("opening VCF {path}"))?;
    let mut meta = Vec::new();
    let mut samples = Vec::new();
    let mut has_column_line = false;
    let mut first_record: Option<String> = None;
    while let Some(line) = reader.next_line()? {
        let text = String::from_utf8_lossy(line);
        if text.starts_with("##") {
            meta.push(text.into_owned());
        } else if text.starts_with('#') {
            samples = text.split('\t').skip(9).map(|s| s.to_string()).collect();
            has_column_line = true;
            break;
        } else {
            first_record = Some(text.into_owned());
            break;
        }
    }
    if !has_column_line {
        bail!("invalid VCF({path}): no header line");
    }
    Ok((VcfHeader::from_lines(meta, samples), reader, first_record))
}

// ---------------------------------------------------------------------------------------------
// records

/// An INFO value as htsjdk holds it: a flag (Boolean true) or text (a String, or a list of Strings
/// kept joined by commas)
#[derive(Clone, Debug, PartialEq)]
pub enum InfoValue {
    Flag,
    Text(String),
}

#[derive(Clone, Debug, Default)]
pub struct Genotype {
    pub sample: String,
    /// GT as written, None when the record has no GT key
    pub gt: Option<String>,
    /// extended attributes: every FORMAT key other than GT, GQ, AD, PL, DP and FT, without missing values
    pub attrs: Vec<(String, String)>,
    pub ad: Option<String>,
    pub dp: Option<String>,
    pub gq: Option<String>,
    pub pl: Option<String>,
    pub ft: Option<String>,
}

impl Genotype {
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }
    pub fn has_attr(&self, key: &str) -> bool {
        self.attr(key).is_some()
    }
    pub fn set_attr(&mut self, key: &str, value: String) {
        match self.attrs.iter_mut().find(|(k, _)| k == key) {
            Some(entry) => entry.1 = value,
            None => self.attrs.push((key.to_string(), value)),
        }
    }
    /// CommonVcfTags.getGenotypeAttributeAsInt
    pub fn attr_int(&self, key: &str, default: i32) -> i32 {
        match self.attr(key) {
            Some(v) => java_text::parse_int(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\"")),
            None => default,
        }
    }
    /// CommonVcfTags.getGenotypeAttributeAsDouble
    pub fn attr_double(&self, key: &str, default: f64) -> f64 {
        match self.attr(key) {
            Some(v) => java_text::parse_double(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\"")),
            None => default,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: i32,
    pub id: String,
    pub ref_allele: String,
    pub alts: Vec<String>,
    /// QUAL / -10, or NO_LOG10_PERROR
    pub log10_p_error: f64,
    /// the filters other than PASS
    pub filters: Vec<String>,
    /// whether filters were applied (FILTER not '.')
    pub filters_applied: bool,
    /// INFO in input order, a repeated key replacing the earlier value in place
    pub info: Vec<(String, InfoValue)>,
    pub genotypes: Vec<Genotype>,
}

/// Allele.create: bases are upper-cased unless the allele is symbolic
fn normalise_alt(alt: &str) -> String {
    let b = alt.as_bytes();
    let symbolic = b.len() > 1
        && (b[0] == b'<' || b[b.len() - 1] == b'>' || alt.contains('[') || alt.contains(']') || b[0] == b'.' || b[b.len() - 1] == b'.');
    if symbolic { alt.to_string() } else { alt.to_ascii_uppercase() }
}

/// The FILTER column of a line without parsing the rest (None for a short line)
pub fn filter_column(line: &str) -> Option<&str> {
    line.split('\t').nth(6)
}

impl VcfRecord {
    pub fn new(chrom: &str, pos: i32, id: &str, ref_allele: &str, alt: &str) -> VcfRecord {
        VcfRecord {
            chrom: chrom.to_string(),
            pos,
            id: id.to_string(),
            ref_allele: ref_allele.to_string(),
            alts: vec![alt.to_string()],
            log10_p_error: NO_LOG10_PERROR,
            ..Default::default()
        }
    }

    /// AbstractVCFCodec.decode of one tab-separated record line
    pub fn parse(line: &str, header: &VcfHeader) -> Result<VcfRecord> {
        let mut cols = line.split('\t');
        let mut next = |name: &str| cols.next().with_context(|| format!("VCF record without {name}: {line}"));
        let chrom = next("CHROM")?.to_string();
        let pos_str = next("POS")?;
        let pos: i32 = java_text::parse_int(pos_str).with_context(|| format!("invalid VCF position: {pos_str}"))?;
        let id = next("ID")?.to_string();
        let ref_str = next("REF")?;
        let alt_str = next("ALT")?;
        let qual_str = next("QUAL")?;
        let filter_str = next("FILTER")?;
        let info_str = next("INFO")?;

        let log10_p_error = if qual_str == "." {
            NO_LOG10_PERROR
        } else {
            let val = java_text::parse_double(qual_str).with_context(|| format!("invalid QUAL: {qual_str}"))?;
            if val < 0.0 && (val - -1.0).abs() < 0.1 { NO_LOG10_PERROR } else { val / -10.0 }
        };

        let filters_applied = filter_str != ".";
        let filters: Vec<String> = if filter_str == "." || filter_str == "PASS" {
            Vec::new()
        } else {
            let mut out: Vec<String> = Vec::new();
            for f in filter_str.split(';') {
                if !out.iter().any(|x| x == f) {
                    out.push(f.to_string());
                }
            }
            out
        };

        let mut rec = VcfRecord {
            chrom,
            pos,
            id,
            ref_allele: ref_str.to_ascii_uppercase(),
            alts: if alt_str == "." { Vec::new() } else { alt_str.split(',').map(normalise_alt).collect() },
            log10_p_error,
            filters,
            filters_applied,
            info: Vec::new(),
            genotypes: Vec::new(),
        };

        if info_str.is_empty() {
            bail!("'.' must be used to specify empty INFO field: {line}");
        }
        if info_str != "." {
            for field in info_str.split(';') {
                if field.is_empty() {
                    continue;
                }
                let (key, value) = match field.split_once('=') {
                    Some((key, value)) => {
                        let is_flag = header.info_def(key).is_some_and(|d| d.ty == "Flag");
                        if is_flag && value == "0" {
                            continue;
                        }
                        let value = if value.is_empty() { "." } else { value };
                        (key, InfoValue::Text(value.to_string()))
                    }
                    None => {
                        // a key without a value whose header type is not Flag reads as missing
                        if header.info_def(field).is_some_and(|d| d.ty != "Flag") {
                            (field, InfoValue::Text(".".to_string()))
                        } else {
                            (field, InfoValue::Flag)
                        }
                    }
                };
                rec.put_attr(key, value);
            }
        }

        if let Some(format) = cols.next() {
            let keys: Vec<&str> = format.split(':').collect();
            for (i, sample) in header.samples.iter().enumerate() {
                let Some(text) = cols.next() else {
                    bail!("VCF record has fewer genotype columns than samples ({} of {}): {}:{}", i, header.samples.len(), rec.chrom, rec.pos);
                };
                let mut g = Genotype { sample: sample.clone(), ..Default::default() };
                for (key, value) in keys.iter().zip(text.split(':')) {
                    match *key {
                        "GT" => g.gt = Some(value.to_string()),
                        _ if value == "." => {}
                        "FT" => g.ft = Some(value.to_string()),
                        "GQ" => g.gq = Some(value.to_string()),
                        "AD" => g.ad = Some(value.to_string()),
                        "PL" => g.pl = Some(value.to_string()),
                        "DP" => g.dp = Some(value.to_string()),
                        _ => g.attrs.push((key.to_string(), value.to_string())),
                    }
                }
                rec.genotypes.push(g);
            }
        }
        Ok(rec)
    }

    /// VariantContext.getEnd: the END attribute, else the last reference base
    pub fn end(&self) -> i32 {
        if let Some(InfoValue::Text(v)) = self.raw_attr("END") {
            if let Some(e) = java_text::parse_int(v) {
                return e;
            }
        }
        self.pos + self.ref_allele.len() as i32 - 1
    }

    /// VariantContext.getPhredScaledQual
    pub fn phred_scaled_qual(&self) -> f64 {
        (self.log10_p_error * -10.0) + 0.0
    }

    /// the first alternate allele ("" without one)
    pub fn alt(&self) -> &str {
        self.alts.first().map(String::as_str).unwrap_or("")
    }

    pub fn has_attr(&self, key: &str) -> bool {
        self.info.iter().any(|(k, _)| k == key)
    }

    pub fn raw_attr(&self, key: &str) -> Option<&InfoValue> {
        self.info.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// CommonInfo.putAttribute: replaces any value of the key
    pub fn put_attr(&mut self, key: &str, value: InfoValue) {
        match self.info.iter_mut().find(|(k, _)| k == key) {
            Some(entry) => entry.1 = value,
            None => self.info.push((key.to_string(), value)),
        }
    }
    pub fn put_text(&mut self, key: &str, value: String) {
        self.put_attr(key, InfoValue::Text(value));
    }
    pub fn put_flag(&mut self, key: &str) {
        self.put_attr(key, InfoValue::Flag);
    }

    /// getAttributeAsString: a list prints as `[a, b]`, a flag as `true`
    pub fn attr_string(&self, key: &str, default: &str) -> String {
        match self.raw_attr(key) {
            None => default.to_string(),
            Some(InfoValue::Flag) => "true".to_string(),
            Some(InfoValue::Text(v)) => {
                if v.contains(',') {
                    format!("[{}]", v.split(',').collect::<Vec<_>>().join(", "))
                } else {
                    v.clone()
                }
            }
        }
    }

    /// getAttributeAsStringList: the items of a list, or the value itself
    pub fn attr_string_list(&self, key: &str) -> Vec<String> {
        match self.raw_attr(key) {
            None => Vec::new(),
            Some(InfoValue::Flag) => vec!["true".to_string()],
            Some(InfoValue::Text(v)) => v.split(',').map(|s| s.to_string()).collect(),
        }
    }

    /// getAttributeAsInt
    pub fn attr_int(&self, key: &str, default: i32) -> i32 {
        match self.raw_attr(key) {
            None => default,
            Some(InfoValue::Text(v)) if v == "." => default,
            Some(InfoValue::Text(v)) => {
                java_text::parse_int(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\" ({key})"))
            }
            Some(InfoValue::Flag) => panic!("ClassCastException: flag {key} read as an integer"),
        }
    }

    /// getAttributeAsDouble
    pub fn attr_double(&self, key: &str, default: f64) -> f64 {
        match self.raw_attr(key) {
            None => default,
            Some(InfoValue::Text(v)) if v == "." => default,
            Some(InfoValue::Text(v)) => {
                java_text::parse_double(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\" ({key})"))
            }
            Some(InfoValue::Flag) => panic!("ClassCastException: flag {key} read as a decimal"),
        }
    }

    /// getAttributeAsBoolean
    pub fn attr_bool(&self, key: &str, default: bool) -> bool {
        match self.raw_attr(key) {
            None => default,
            Some(InfoValue::Flag) => true,
            Some(InfoValue::Text(v)) => v.eq_ignore_ascii_case("true"),
        }
    }

    /// getAttributeAsIntList
    pub fn attr_int_list(&self, key: &str, default: i32) -> Vec<i32> {
        match self.raw_attr(key) {
            None => Vec::new(),
            Some(InfoValue::Flag) => panic!("ClassCastException: flag {key} read as an integer list"),
            Some(InfoValue::Text(v)) => v
                .split(',')
                .map(|x| {
                    if x == "." {
                        default
                    } else {
                        java_text::parse_int(x).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{x}\" ({key})"))
                    }
                })
                .collect(),
        }
    }

    pub fn genotype(&self, sample: &str) -> Option<&Genotype> {
        self.genotypes.iter().find(|g| g.sample == sample)
    }
    pub fn genotype_mut(&mut self, sample: &str) -> Option<&mut Genotype> {
        self.genotypes.iter_mut().find(|g| g.sample == sample)
    }

    /// VCFEncoder.encode: the record line without its line end
    pub fn encode(&self, header: &VcfHeader) -> Result<String> {
        let mut s = String::with_capacity(1024);
        s.push_str(&self.chrom);
        s.push('\t');
        s.push_str(&self.pos.to_string());
        s.push('\t');
        s.push_str(if self.id.is_empty() { "." } else { &self.id });
        s.push('\t');
        s.push_str(&self.ref_allele);
        s.push('\t');
        if self.alts.is_empty() {
            s.push('.');
        } else {
            s.push_str(&self.alts.join(","));
        }
        s.push('\t');
        if self.log10_p_error == NO_LOG10_PERROR {
            s.push('.');
        } else {
            s.push_str(&format_qual((self.log10_p_error * -10.0) + 0.0));
        }
        s.push('\t');
        if !self.filters.is_empty() {
            for f in &self.filters {
                if !header.has_filter_line(f) {
                    bail!("Key {f} found in VariantContext field FILTER at {}:{} but this key isn't defined in the VCFHeader", self.chrom, self.pos);
                }
            }
            let mut filters: Vec<&str> = self.filters.iter().map(|f| f.as_str()).collect();
            filters.sort();
            s.push_str(&filters.join(";"));
        } else if self.filters_applied {
            s.push_str("PASS");
        } else {
            s.push('.');
        }
        s.push('\t');

        let mut keys: Vec<&(String, InfoValue)> = self.info.iter().collect();
        keys.sort_by(|a, b| a.0.cmp(&b.0));
        if keys.is_empty() {
            s.push('.');
        }
        for (i, (key, value)) in keys.iter().enumerate() {
            let Some(def) = header.info_def(key) else {
                bail!("Key {key} found in VariantContext field INFO at {}:{} but this key isn't defined in the VCFHeader", self.chrom, self.pos);
            };
            if i > 0 {
                s.push(';');
            }
            s.push_str(key);
            if let InfoValue::Text(v) = value {
                if !v.is_empty() && def.count != Count::Fixed(0) {
                    s.push('=');
                    s.push_str(v);
                }
            }
        }

        if !header.samples.is_empty() {
            self.encode_genotypes(header, &mut s)?;
        }
        Ok(s)
    }

    /// calcVCFGenotypeKeys and appendGenotypeData
    fn encode_genotypes(&self, header: &VcfHeader, s: &mut String) -> Result<()> {
        let genotypes = &self.genotypes;
        let mut format_keys: Vec<&str> = Vec::new();
        let (mut saw_gt, mut saw_gq, mut saw_dp, mut saw_ad, mut saw_pl, mut saw_ft) = (false, false, false, false, false, false);
        for g in genotypes {
            for (k, _) in &g.attrs {
                if !format_keys.contains(&k.as_str()) {
                    format_keys.push(k);
                }
            }
            saw_gt |= g.gt.is_some();
            saw_gq |= g.gq.is_some();
            saw_dp |= g.dp.is_some();
            saw_ad |= g.ad.is_some();
            saw_pl |= g.pl.is_some();
            saw_ft |= g.ft.as_deref().is_some_and(|f| f != "PASS");
        }
        for (saw, key) in [(saw_gq, "GQ"), (saw_dp, "DP"), (saw_ad, "AD"), (saw_pl, "PL"), (saw_ft, "FT")] {
            if saw && !format_keys.contains(&key) {
                format_keys.push(key);
            }
        }
        format_keys.sort();
        if saw_gt {
            format_keys.insert(0, "GT");
        } else if format_keys.is_empty() {
            format_keys.push("GT");
        }
        let has_gt = format_keys.contains(&"GT");
        for key in &format_keys {
            if !header.has_format_line(key) {
                bail!("Key {key} found in VariantContext field FORMAT at {}:{} but this key isn't defined in the VCFHeader", self.chrom, self.pos);
            }
        }

        s.push('\t');
        s.push_str(&format_keys.join(":"));

        let missing = Genotype::default();
        for sample in &header.samples {
            s.push('\t');
            let (g, absent) = match genotypes.iter().find(|g| &g.sample == sample) {
                Some(g) => (g, false),
                None => (&missing, true),
            };
            let mut values: Vec<String> = Vec::with_capacity(format_keys.len());
            for key in &format_keys {
                match *key {
                    "GT" => {}
                    "FT" => values.push(g.ft.clone().unwrap_or_else(|| "PASS".to_string())),
                    "DP" => values.push(g.dp.clone().unwrap_or_else(|| ".".to_string())),
                    "AD" => values.push(g.ad.clone().unwrap_or_else(|| ".".to_string())),
                    "PL" => values.push(g.pl.clone().unwrap_or_else(|| ".".to_string())),
                    "GQ" => values.push(match &g.gq {
                        // read with Math.round(Double.valueOf())
                        Some(v) => match java_text::parse_double(v) {
                            Some(d) => ((d + 0.5).floor() as i64).to_string(),
                            None => v.clone(),
                        },
                        None => ".".to_string(),
                    }),
                    _ => match g.attr(key) {
                        Some(v) => values.push(v.to_string()),
                        None => {
                            let n = match header.format_def(key).map(|d| &d.count) {
                                Some(Count::Fixed(n)) => *n,
                                Some(Count::A) => self.alts.len() as i32,
                                Some(Count::R) => self.alts.len() as i32 + 1,
                                Some(Count::G) => {
                                    let alleles = self.alts.len() as i32 + 1;
                                    alleles * (alleles + 1) / 2
                                }
                                _ => 1,
                            };
                            values.push(if n > 1 { vec!["."; n as usize].join(",") } else { ".".to_string() });
                        }
                    },
                }
            }
            while values.last().is_some_and(|v| is_missing_value(v)) {
                values.pop();
            }
            if has_gt {
                match (&g.gt, absent) {
                    (Some(gt), _) => s.push_str(gt),
                    (None, true) => s.push_str("./."),
                    (None, false) => bail!("GTs cannot be missing for some samples if they are available for others in the record"),
                }
            }
            for (i, v) in values.iter().enumerate() {
                if i > 0 || has_gt {
                    s.push(':');
                }
                s.push_str(v);
            }
        }
        Ok(())
    }
}

/// VCFEncoder.formatQualValue
pub fn format_qual(phred: f64) -> String {
    let s = crate::fmt::java_fixed(phred, 2);
    match s.strip_suffix(".00") {
        Some(t) => t.to_string(),
        None => s,
    }
}

/// VCFEncoder.isMissingValue
fn is_missing_value(s: &str) -> bool {
    s.bytes().all(|c| c == b'.' || c == b',')
}

// ---------------------------------------------------------------------------------------------
// BGZF and the tabix index

/// htsjdk BlockCompressedStreamConstants.DEFAULT_UNCOMPRESSED_BLOCK_SIZE (seen in its files: 65498 bytes)
const BLOCK_SIZE: usize = 65498;

/// The empty block that ends a BGZF file
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
];

/// A BGZF writer that reports htsjdk's virtual file pointers (BlockCompressedOutputStream).
pub struct BgzfWriter<W: Write> {
    out: W,
    buffer: Vec<u8>,
    /// compressed bytes written so far: the address of the block being filled
    block_address: u64,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(out: W) -> Self {
        BgzfWriter { out, buffer: Vec::with_capacity(BLOCK_SIZE), block_address: 0 }
    }
    /// BlockCompressedOutputStream.getFilePointer
    pub fn file_pointer(&self) -> u64 {
        (self.block_address << 16) | self.buffer.len() as u64
    }
    pub fn write_all(&mut self, mut data: &[u8]) -> Result<()> {
        while !data.is_empty() {
            let n = (BLOCK_SIZE - self.buffer.len()).min(data.len());
            self.buffer.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.buffer.len() == BLOCK_SIZE {
                self.deflate_block()?;
            }
        }
        Ok(())
    }
    /// Writes the bytes buffered as one block (none if empty)
    pub fn deflate_block(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let mut enc = flate2::write::DeflateEncoder::new(Vec::with_capacity(self.buffer.len()), flate2::Compression::new(5));
        enc.write_all(&self.buffer)?;
        let mut data = enc.finish()?;
        if data.len() + 26 > 65536 {
            // incompressible: store
            let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::none());
            enc.write_all(&self.buffer)?;
            data = enc.finish()?;
        }
        let mut crc = flate2::Crc::new();
        crc.update(&self.buffer);
        let bsize = (data.len() + 25) as u16;
        let mut block = Vec::with_capacity(data.len() + 26);
        block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00]);
        block.extend_from_slice(&bsize.to_le_bytes());
        block.extend_from_slice(&data);
        block.extend_from_slice(&crc.sum().to_le_bytes());
        block.extend_from_slice(&(self.buffer.len() as u32).to_le_bytes());
        self.out.write_all(&block)?;
        self.block_address += block.len() as u64;
        self.buffer.clear();
        Ok(())
    }
    pub fn finish(mut self) -> Result<W> {
        self.deflate_block()?;
        self.out.write_all(&BGZF_EOF)?;
        self.out.flush()?;
        Ok(self.out)
    }
}

/// GenomicIndexUtil.regionToBin (reg2bin), zero-based half-open
fn region_to_bin(beg: i32, end: i32) -> usize {
    let end = end - 1;
    let b = if beg >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (beg >> 26)
    } else {
        0
    };
    b as usize
}

/// LinearIndex.convertToLinearIndexOffset
fn linear_index_offset(pos: i32) -> usize {
    ((if pos <= 0 { 0 } else { pos - 1 }) >> 14) as usize
}

/// BinningIndexBuilder for one reference sequence
struct SequenceIndex {
    /// bin number to chunks (start, end); kept in bin order when written
    bins: std::collections::BTreeMap<usize, Vec<(u64, u64)>>,
    linear: Vec<u64>,
    largest_window: i64,
}

impl SequenceIndex {
    fn new() -> Self {
        SequenceIndex { bins: Default::default(), linear: Vec::new(), largest_window: -1 }
    }
    /// processFeature: a record at [start, end] (1-based) from file pointer chunk_start to chunk_end
    fn process(&mut self, start: i32, end: i32, chunk_start: u64, chunk_end: u64) {
        let bin = region_to_bin(start - 1, if end <= 0 { start } else { end });
        let chunks = self.bins.entry(bin).or_default();
        match chunks.last_mut() {
            // coalesce chunks in the same or adjacent blocks (BlockCompressedFilePointerUtil.areInSameOrAdjacentBlocks)
            Some(last) if (last.1 >> 16) == (chunk_start >> 16) || (last.1 >> 16) + 1 == (chunk_start >> 16) => last.1 = chunk_end,
            _ => chunks.push((chunk_start, chunk_end)),
        }
        let start_window = linear_index_offset(start);
        let end_window = if end == 0 { start_window } else { linear_index_offset(end) };
        if end_window as i64 > self.largest_window {
            self.largest_window = end_window as i64;
        }
        if self.linear.len() <= end_window {
            self.linear.resize(end_window + 1, 0);
        }
        for w in start_window..=end_window {
            if self.linear[w] == 0 || chunk_start < self.linear[w] {
                self.linear[w] = chunk_start;
            }
        }
    }
    /// generateIndexContent: the linear index up to the largest window, zeros after the first entry
    /// filled with the one before
    fn write(&mut self, out: &mut Vec<u8>) {
        out.extend_from_slice(&(self.bins.len() as i32).to_le_bytes());
        for (bin, chunks) in &self.bins {
            out.extend_from_slice(&(*bin as u32).to_le_bytes());
            out.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
            for (s, e) in chunks {
                out.extend_from_slice(&s.to_le_bytes());
                out.extend_from_slice(&e.to_le_bytes());
            }
        }
        let n = (self.largest_window + 1) as usize;
        self.linear.resize(n, 0);
        let mut last = 0u64;
        for v in self.linear.iter_mut() {
            if *v == 0 {
                *v = last;
            } else {
                last = *v;
            }
        }
        out.extend_from_slice(&(n as i32).to_le_bytes());
        for v in &self.linear {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
}

/// A bgzipped VCF written as htsjdk's VariantContextWriter with BLOCK_COMPRESSED_VCF and
/// INDEX_ON_THE_FLY writes it: the BGZF stream and the tabix index (`<path>.tbi`) on closing.
pub struct VcfFileWriter {
    path: String,
    out: BgzfWriter<std::io::BufWriter<std::fs::File>>,
    names: Vec<String>,
    indexes: Vec<SequenceIndex>,
    /// the record waiting for the next one's start as its chunk end: (start, end, file pointer)
    pending: Option<(i32, i32, u64)>,
}

impl VcfFileWriter {
    pub fn create(path: &str, header: &VcfHeader) -> Result<VcfFileWriter> {
        Self::create_with_header_text(path, &header.text())
    }

    /// The same, for a header held elsewhere (`variant_context`): `text` is the whole header, `#CHROM`
    /// line included
    pub fn create_with_header_text(path: &str, text: &str) -> Result<VcfFileWriter> {
        let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
        let mut out = BgzfWriter::new(std::io::BufWriter::with_capacity(1 << 20, file));
        out.write_all(text.as_bytes())?;
        Ok(VcfFileWriter { path: path.to_string(), out, names: Vec::new(), indexes: Vec::new(), pending: None })
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    /// Writes one encoded record (without its line end) spanning `start` to `end` of `chrom`
    pub fn write_line(&mut self, chrom: &str, start: i32, end: i32, line: &str) -> Result<()> {
        let pos = self.out.file_pointer();
        if self.names.last().map(|s| s.as_str()) != Some(chrom) {
            if self.names.iter().any(|n| n == chrom) {
                bail!("Sequence {chrom} added out sequence of order");
            }
            if let Some((s, e, at)) = self.pending.take() {
                self.indexes.last_mut().unwrap().process(s, e, at, pos);
            }
            self.names.push(chrom.to_string());
            self.indexes.push(SequenceIndex::new());
        } else if let Some((s, e, at)) = self.pending.take() {
            // TabixIndexCreator.addFeature: positions may not go back within a sequence
            if start < s {
                bail!("Features added out of order: previous start {s} > next start {start} on {chrom}");
            }
            self.indexes.last_mut().unwrap().process(s, e, at, pos);
        }
        self.out.write_all(line.as_bytes())?;
        self.out.write_all(b"\n")?;
        self.pending = Some((start, end, pos));
        Ok(())
    }

    /// Writes a record
    pub fn write(&mut self, rec: &VcfRecord, header: &VcfHeader) -> Result<()> {
        let line = rec.encode(header)?;
        self.write_line(&rec.chrom, rec.pos, rec.end(), &line)
    }

    /// Ends the BGZF stream and writes the tabix index
    pub fn close(mut self) -> Result<()> {
        // TabixIndexCreator.finalizeIndex takes the file pointer once the last block is written: the start
        // of the next block
        self.out.deflate_block()?;
        let end = self.out.file_pointer();
        if let Some((s, e, at)) = self.pending.take() {
            self.indexes.last_mut().unwrap().process(s, e, at, end);
        }
        self.out.finish()?;

        // TabixIndex.write, bgzipped
        let mut idx: Vec<u8> = Vec::new();
        idx.extend_from_slice(b"TBI\x01");
        for v in [self.names.len() as i32, 2, 1, 2, 0, b'#' as i32, 0] {
            idx.extend_from_slice(&v.to_le_bytes());
        }
        let name_block: usize = self.names.iter().map(|n| n.len() + 1).sum();
        idx.extend_from_slice(&(name_block as i32).to_le_bytes());
        for n in &self.names {
            idx.extend_from_slice(n.as_bytes());
            idx.push(0);
        }
        for s in self.indexes.iter_mut() {
            s.write(&mut idx);
        }
        let tbi = format!("{}.tbi", self.path);
        let file = std::fs::File::create(&tbi).with_context(|| format!("creating {tbi}"))?;
        let mut w = BgzfWriter::new(std::io::BufWriter::new(file));
        w.write_all(&idx)?;
        w.finish()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> VcfHeader {
        let meta = vec![
            "##fileformat=VCFv4.2".to_string(),
            filter_line("PASS", "All filters passed"),
            filter_line("LowQual", "Low quality"),
            info_line("DP", "1", "Float", "Depth"),
            info_line("FL", "0", "Flag", "A flag"),
            info_line("LST", ".", "String", "A list"),
            format_line("GT", "1", "String", "Genotype"),
            format_line("AD", "R", "Integer", "Allelic depths for the ref and alt alleles in the order listed"),
            format_line("AF", "1", "Float", "AF"),
            "##contig=<ID=chr2,length=10>".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "##source=x".to_string(),
        ];
        VcfHeader::from_lines(meta, vec!["N".to_string(), "T".to_string()])
    }

    #[test]
    fn standard_lines_repaired_and_sorted() {
        let h = header();
        // DP declared with another type is replaced by the standard line
        assert!(h.meta_lines.iter().any(|l| l.contains("Approximate read depth; some reads may have been filtered")));
        let text = h.text();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        // contigs kept in their order where "##contig" sorts
        let c2 = lines.iter().position(|l| l.contains("ID=chr2")).unwrap();
        let c1 = lines.iter().position(|l| l.contains("ID=chr1")).unwrap();
        assert!(c2 < c1);
        assert!(lines.iter().position(|l| l.starts_with("##source")).unwrap() > c1);
        assert!(lines.last().unwrap().ends_with("FORMAT\tN\tT"));
    }

    #[test]
    fn record_round_trip() {
        let h = header();
        let line = "chr1\t5\t.\tac\tG,T\t45.5\tPASS\tLST=a,b;FL;DP=7;FL2=0\tGT:AD:AF\t0/1:3,4:.\t./.:5,6:0.5";
        // FL2 is not in the header: a flag-typed check cannot skip it, so encoding fails
        let r = VcfRecord::parse(line, &h).unwrap();
        assert!(r.encode(&h).is_err());
        let line = "chr1\t5\t.\tac\tG,T\t45.5\tPASS\tLST=a,b;FL;DP=7\tGT:AD:AF\t0/1:3,4:.\t./.:5,6:0.5";
        let mut r = VcfRecord::parse(line, &h).unwrap();
        assert_eq!(r.ref_allele, "AC");
        assert_eq!(r.attr_string("LST", ""), "[a, b]");
        assert_eq!(r.attr_string_list("LST"), vec!["a", "b"]);
        assert_eq!(r.attr_int("DP", 0), 7);
        assert!(r.attr_bool("FL", false));
        assert_eq!(r.end(), 6);
        assert_eq!(r.encode(&h).unwrap(), "chr1\t5\t.\tAC\tG,T\t45.50\tPASS\tDP=7;FL;LST=a,b\tGT:AD:AF\t0/1:3,4\t./.:5,6:0.5");
        r.filters.push("LowQual".to_string());
        r.put_flag("FL");
        r.put_text("DP", "8".to_string());
        assert!(r.encode(&h).unwrap().contains("\tLowQual\tDP=8;FL;LST=a,b\t"));
        // a flag given as =0 is dropped; an empty value reads as missing
        let r = VcfRecord::parse("chr1\t5\t.\tA\tG\t.\t.\tFL=0;LST=\tGT\t0/1\t0/1", &h).unwrap();
        assert!(!r.has_attr("FL"));
        assert_eq!(r.attr_string("LST", ""), ".");
        assert_eq!(r.encode(&h).unwrap(), "chr1\t5\t.\tA\tG\t.\t.\tLST=.\tGT\t0/1\t0/1");
    }

    #[test]
    fn bins() {
        assert_eq!(region_to_bin(17592104, 17592105), 5754);
        assert_eq!(region_to_bin(0, 1), 4681);
        assert_eq!(linear_index_offset(17592105), 1073);
    }

    /// The index answers queries through htslib
    #[test]
    fn writer_and_index_read_back() {
        let dir = std::env::temp_dir().join(format!("hmf_common_vcf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.vcf.gz").to_str().unwrap().to_string();
        let h = VcfHeader::from_lines(
            vec![
                "##fileformat=VCFv4.2".to_string(),
                "##contig=<ID=chr1,length=100000>".to_string(),
                "##contig=<ID=chr2,length=100000>".to_string(),
            ],
            Vec::new(),
        );
        let mut w = VcfFileWriter::create(&path, &h).unwrap();
        for (chr, pos) in [("chr1", 100), ("chr1", 20000), ("chr1", 70000), ("chr2", 5)] {
            let r = VcfRecord::new(chr, pos, ".", "A", "G");
            w.write(&r, &h).unwrap();
        }
        assert!(w.write_line("chr1", 5, 5, "x").is_err());
        assert!(w.write_line("chr2", 4, 4, "x").is_err());
        w.close().unwrap();
        use rust_htslib::bcf::{IndexedReader, Read};
        let mut r = IndexedReader::from_path(&path).unwrap();
        let rid = r.header().name2rid(b"chr1").unwrap();
        r.fetch(rid, 19000, Some(80000)).unwrap();
        let positions: Vec<i64> = r.records().map(|x| x.unwrap().pos() + 1).collect();
        assert_eq!(positions, vec![20000, 70000]);
        std::fs::remove_dir_all(&dir).ok();
    }
}
