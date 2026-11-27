//! VCF records and headers as htsjdk 2.23 reads, holds and writes them (VCFCodec, VariantContext,
//! CommonInfo, Genotype, VCFEncoder, VCFHeader), for the tools that read a VCF, change some INFO or
//! genotype fields and write it again (PURPLE, and PAVE and LINX to come).
//!
//! What is reproduced:
//! - INFO values keep the type Java gives them: read values are strings (a value with commas is a
//!   List<String>), flags are Boolean.TRUE; values a tool sets keep their Java type (Double, Integer,
//!   Boolean, String, lists), so that the getAttributeAs* accessors return what Java returns and the
//!   writer formats them as VCFEncoder.formatVCFField does (doubles with formatVCFDouble).
//! - Genotypes stay unparsed until something asks for them (LazyGenotypesContext); a record whose
//!   genotypes were never decoded is written with its genotype text as read, one whose genotypes were
//!   decoded is encoded again (FORMAT keys sorted with GT first, trailing missing values dropped,
//!   missing values of a multi-valued field written as ".,.").
//! - INFO keys are written sorted, FILTER values sorted, QUAL with "%.2f" less a trailing ".00".
//! - Header lines are written sorted as strings without duplicates, with the contig lines kept in
//!   their order at the place where "##contig" sorts (VCFHeader's TreeSet of lines).
//! - The output is BGZF with a tabix index, both as htsjdk writes them (`crate::vcf::VcfFileWriter`).

use anyhow::{anyhow, bail, Context, Result};
use std::io::{BufRead, Write};

use crate::fmt::vcf_double;

/// VariantContext.NO_LOG10_PERROR
pub const NO_LOG10_PERROR: f64 = 1.0;

// ---------------------------------------------------------------------------------------------
// values

/// An INFO or genotype attribute value with the Java type it would have.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// read from a file: a String, or a List<String> when it holds commas
    Raw(String),
    /// Boolean (true for a flag read from a file)
    Bool(bool),
    Int(i32),
    Double(f64),
    /// a String set by the program (commas do not make it a list)
    Str(String),
    StrList(Vec<String>),
    IntList(Vec<i32>),
    DoubleList(Vec<f64>),
}

impl Value {
    pub fn flag() -> Value {
        Value::Bool(true)
    }

    /// VCFEncoder.formatVCFField; None where Java returns null (Boolean false: the key is left out)
    pub fn encode(&self) -> Option<String> {
        match self {
            Value::Raw(s) | Value::Str(s) => Some(s.clone()),
            Value::Bool(true) => Some(String::new()),
            Value::Bool(false) => None,
            Value::Int(i) => Some(i.to_string()),
            Value::Double(d) => Some(vcf_double(*d)),
            Value::StrList(v) => Some(if v.is_empty() { ".".into() } else { v.join(",") }),
            Value::IntList(v) => {
                Some(if v.is_empty() { ".".into() } else { v.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",") })
            }
            Value::DoubleList(v) => Some(if v.is_empty() { ".".into() } else { v.iter().map(|d| vcf_double(*d)).collect::<Vec<_>>().join(",") }),
        }
    }

    /// the list items as Java's getAttributeAsList sees them: a read value with commas is a list
    fn items(&self) -> Vec<Item<'_>> {
        match self {
            Value::Raw(s) => {
                if s.contains(',') {
                    s.split(',').map(Item::Str).collect()
                } else {
                    vec![Item::Str(s)]
                }
            }
            Value::Str(s) => vec![Item::Str(s)],
            Value::Bool(b) => vec![Item::Bool(*b)],
            Value::Int(i) => vec![Item::Int(*i)],
            Value::Double(d) => vec![Item::Double(*d)],
            Value::StrList(v) => v.iter().map(|s| Item::Str(s)).collect(),
            Value::IntList(v) => v.iter().map(|i| Item::Int(*i)).collect(),
            Value::DoubleList(v) => v.iter().map(|d| Item::Double(*d)).collect(),
        }
    }

    /// String.valueOf(value)
    pub fn java_string(&self) -> String {
        match self {
            Value::Raw(s) => {
                if s.contains(',') {
                    format!("[{}]", s.split(',').collect::<Vec<_>>().join(", "))
                } else {
                    s.clone()
                }
            }
            Value::Str(s) => s.clone(),
            Value::Bool(b) => b.to_string(),
            Value::Int(i) => i.to_string(),
            Value::Double(d) => crate::fmt::java_double_to_string(*d),
            Value::StrList(v) => format!("[{}]", v.join(", ")),
            Value::IntList(v) => format!("[{}]", v.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(", ")),
            Value::DoubleList(v) => {
                format!("[{}]", v.iter().map(|d| crate::fmt::java_double_to_string(*d)).collect::<Vec<_>>().join(", "))
            }
        }
    }
}

#[derive(Clone, Copy)]
enum Item<'a> {
    Str(&'a str),
    Bool(bool),
    Int(i32),
    Double(f64),
}

impl Item<'_> {
    fn java_string(&self) -> String {
        match self {
            Item::Str(s) => s.to_string(),
            Item::Bool(b) => b.to_string(),
            Item::Int(i) => i.to_string(),
            Item::Double(d) => crate::fmt::java_double_to_string(*d),
        }
    }
}

fn parse_int_or_panic(s: &str, key: &str) -> i32 {
    crate::java_text::parse_int(s).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{s}\" ({key})"))
}

fn parse_double_or_panic(s: &str, key: &str) -> f64 {
    crate::java_text::parse_double(s).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{s}\" ({key})"))
}

/// A set of attributes as CommonInfo holds them: insertion order is kept (the writer sorts the keys).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Attributes(pub Vec<(String, Value)>);

impl Attributes {
    pub fn get(&self, key: &str) -> Option<&Value> {
        self.0.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }
    pub fn has(&self, key: &str) -> bool {
        self.get(key).is_some()
    }
    /// putAttribute (allowOverwrites)
    pub fn put(&mut self, key: &str, value: Value) {
        match self.0.iter_mut().find(|(k, _)| k == key) {
            Some(entry) => entry.1 = value,
            None => self.0.push((key.to_string(), value)),
        }
    }
    pub fn remove(&mut self, key: &str) -> Option<Value> {
        let i = self.0.iter().position(|(k, _)| k == key)?;
        Some(self.0.remove(i).1)
    }

    /// getAttributeAsString(key, default)
    pub fn as_string(&self, key: &str, default: &str) -> String {
        match self.get(key) {
            None => default.to_string(),
            Some(v) => v.java_string(),
        }
    }
    /// getAttributeAsInt(key, default)
    pub fn as_int(&self, key: &str, default: i32) -> i32 {
        match self.get(key) {
            None => default,
            Some(Value::Int(i)) => *i,
            Some(Value::Raw(s)) | Some(Value::Str(s)) => {
                if s == "." {
                    default
                } else {
                    parse_int_or_panic(s, key)
                }
            }
            Some(v) => panic!("ClassCastException: {key} ({v:?}) read as an integer"),
        }
    }
    /// getAttributeAsDouble(key, default)
    pub fn as_double(&self, key: &str, default: f64) -> f64 {
        match self.get(key) {
            None => default,
            Some(Value::Double(d)) => *d,
            Some(Value::Int(i)) => *i as f64,
            Some(Value::Raw(s)) | Some(Value::Str(s)) => {
                if s == "." {
                    default
                } else {
                    parse_double_or_panic(s, key)
                }
            }
            Some(v) => panic!("ClassCastException: {key} ({v:?}) read as a decimal"),
        }
    }
    /// getAttributeAsBoolean(key, default)
    pub fn as_bool(&self, key: &str, default: bool) -> bool {
        match self.get(key) {
            None => default,
            Some(Value::Bool(b)) => *b,
            Some(Value::Int(i)) => *i != 0,
            Some(Value::Raw(s)) | Some(Value::Str(s)) => s.eq_ignore_ascii_case("true"),
            Some(v) => panic!("ClassCastException: {key} ({v:?}) read as a boolean"),
        }
    }
    /// getAttributeAsStringList(key, default)
    pub fn as_string_list(&self, key: &str, default: &str) -> Vec<String> {
        match self.get(key) {
            None => Vec::new(),
            Some(v) => v.items().iter().map(|i| if let Item::Str(".") = i { default.to_string() } else { i.java_string() }).collect(),
        }
    }
    /// getAttributeAsIntList(key, default)
    pub fn as_int_list(&self, key: &str, default: i32) -> Vec<i32> {
        match self.get(key) {
            None => Vec::new(),
            Some(v) => v
                .items()
                .iter()
                .map(|i| match i {
                    Item::Int(n) => *n,
                    Item::Str(".") => default,
                    Item::Str(s) => parse_int_or_panic(s, key),
                    other => parse_int_or_panic(&other.java_string(), key),
                })
                .collect(),
        }
    }
    /// getAttributeAsDoubleList(key, default)
    pub fn as_double_list(&self, key: &str, default: f64) -> Vec<f64> {
        match self.get(key) {
            None => Vec::new(),
            Some(v) => v
                .items()
                .iter()
                .map(|i| match i {
                    Item::Double(d) => *d,
                    Item::Int(n) => *n as f64,
                    Item::Str(".") => default,
                    Item::Str(s) => parse_double_or_panic(s, key),
                    other => parse_double_or_panic(&other.java_string(), key),
                })
                .collect(),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// genotypes

/// A decoded genotype (GenotypeBuilder's fields).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Genotype {
    pub sample: String,
    /// GT as allele indices (None for a no-call "."); empty when the record has no GT
    pub alleles: Vec<Option<u32>>,
    pub phased: bool,
    pub gq: Option<i32>,
    pub dp: Option<i32>,
    pub ad: Option<Vec<i32>>,
    pub pl: Option<Vec<i32>>,
    /// FT other than PASS
    pub filters: Option<String>,
    /// every other key with a value other than "."
    pub attrs: Attributes,
}

impl Genotype {
    pub fn has_ad(&self) -> bool {
        self.ad.is_some()
    }
    /// Genotype.getType(): None for UNAVAILABLE (no alleles)
    pub fn genotype_type(&self) -> GenotypeType {
        if self.alleles.is_empty() {
            return GenotypeType::Unavailable;
        }
        let calls: Vec<u32> = self.alleles.iter().flatten().copied().collect();
        if calls.len() < self.alleles.len() {
            return if calls.is_empty() { GenotypeType::NoCall } else { GenotypeType::Mixed };
        }
        let first = calls[0];
        if calls.iter().all(|&a| a == first) {
            if first == 0 { GenotypeType::HomRef } else { GenotypeType::HomVar }
        } else {
            GenotypeType::Het
        }
    }

    /// CommonVcfTags.getGenotypeAttributeAsInt
    pub fn attr_int(&self, key: &str, default: i32) -> i32 {
        match self.attrs.get(key) {
            None => default,
            Some(v) => parse_int_or_panic(&v.java_string(), key),
        }
    }
    /// CommonVcfTags.getGenotypeAttributeAsDouble
    pub fn attr_double(&self, key: &str, default: f64) -> f64 {
        match self.attrs.get(key) {
            None => default,
            Some(v) => parse_double_or_panic(&v.java_string(), key),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenotypeType {
    NoCall,
    HomRef,
    Het,
    HomVar,
    Unavailable,
    Mixed,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Genotypes {
    /// the FORMAT column and the sample columns as read (LazyGenotypesContext)
    Lazy(String),
    Decoded(Vec<Genotype>),
}

impl Default for Genotypes {
    fn default() -> Self {
        Genotypes::Decoded(Vec::new())
    }
}

fn decode_ints(s: &str) -> Vec<i32> {
    s.split(',').map(|x| if x == "." { -1 } else { crate::java_text::parse_int(x).unwrap_or(-1) }).collect()
}

/// VCFCodec.createGenotypeMap
fn decode_genotypes(text: &str, samples: &[String]) -> Result<Vec<Genotype>> {
    let mut cols = text.split('\t');
    let Some(format) = cols.next() else { return Ok(Vec::new()) };
    let keys: Vec<&str> = format.split(':').collect();
    let mut out = Vec::with_capacity(samples.len());
    for (i, sample) in samples.iter().enumerate() {
        let Some(sample_text) = cols.next() else {
            bail!("VCF record has fewer genotype columns than samples ({} of {})", i, samples.len());
        };
        let mut g = Genotype { sample: sample.clone(), ..Default::default() };
        for (key, value) in keys.iter().zip(sample_text.split(':')) {
            if *key == "GT" {
                g.phased = value.contains('|');
                g.alleles = value
                    .split(['/', '|'])
                    .map(|a| if a == "." { None } else { a.parse::<u32>().ok() })
                    .collect();
                continue;
            }
            if value == "." {
                continue;
            }
            match *key {
                "FT" => {
                    if value != "PASS" {
                        g.filters = Some(value.to_string());
                    }
                }
                "GQ" => g.gq = Some(crate::java_text::parse_double(value).map(|d| crate::java_math::java_round(d) as i32).unwrap_or(-1)),
                "AD" => g.ad = Some(decode_ints(value)),
                "PL" => g.pl = Some(decode_ints(value)),
                "DP" => g.dp = Some(crate::java_text::parse_int(value).unwrap_or(-1)),
                _ => g.attrs.0.push((key.to_string(), Value::Raw(value.to_string()))),
            }
        }
        out.push(g);
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// records

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: i32,
    pub id: String,
    pub ref_allele: String,
    pub alts: Vec<String>,
    /// QUAL / -10, or NO_LOG10_PERROR
    pub log10_p_error: f64,
    /// None: "." (filters not applied); Some(empty): PASS
    pub filters: Option<Vec<String>>,
    pub info: Attributes,
    pub genotypes: Genotypes,
}

/// Allele.create: bases are upper-cased unless the allele is symbolic
fn normalise_alt(alt: &str) -> String {
    let b = alt.as_bytes();
    let symbolic = b.len() > 1
        && (b[0] == b'<' || b[b.len() - 1] == b'>' || alt.contains('[') || alt.contains(']') || b[0] == b'.' || b[b.len() - 1] == b'.');
    if symbolic { alt.to_string() } else { alt.to_ascii_uppercase() }
}

/// VariantContext.Type
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

pub fn is_symbolic_allele(a: &str) -> bool {
    let b = a.as_bytes();
    b.len() > 1 && (b[0] == b'<' || b[b.len() - 1] == b'>' || a.contains('[') || a.contains(']') || b[0] == b'.' || b[b.len() - 1] == b'.')
        || a == "*"
}

impl VcfRecord {
    /// VariantContext.getEnd as AbstractVCFCodec sets the stop: the END attribute, else the last reference
    /// base (the tabix index takes it)
    pub fn end(&self) -> i32 {
        match self.info.get("END") {
            Some(Value::Int(e)) => *e,
            Some(Value::Raw(v)) | Some(Value::Str(v)) => crate::java_text::parse_int(v).unwrap_or(self.pos + self.ref_allele.len() as i32 - 1),
            _ => self.pos + self.ref_allele.len() as i32 - 1,
        }
    }

    /// AbstractVCFCodec.decode of one tab-separated line; genotypes are kept unparsed
    pub fn parse(line: &str, header: &VcfHeader) -> Result<VcfRecord> {
        let mut cols = line.splitn(9, '\t');
        let mut next = |name: &str| cols.next().with_context(|| format!("VCF record without {name}: {line}"));
        let chrom = next("CHROM")?.to_string();
        let pos_str = next("POS")?;
        let pos: i32 = pos_str.parse().with_context(|| format!("invalid VCF position: {pos_str}"))?;
        let id = next("ID")?.to_string();
        let ref_str = next("REF")?;
        let alt_str = next("ALT")?;
        let qual_str = next("QUAL")?;
        let filter_str = next("FILTER")?;
        let info_str = next("INFO")?;
        let rest = cols.next();

        let log10_p_error = if qual_str == "." {
            NO_LOG10_PERROR
        } else {
            let val = crate::java_text::parse_double(qual_str).with_context(|| format!("invalid QUAL: {qual_str}"))?;
            if val < 0.0 && (val - -1.0).abs() < 0.1 { NO_LOG10_PERROR } else { val / -10.0 }
        };

        let filters = if filter_str == "." {
            None
        } else if filter_str == "PASS" {
            Some(Vec::new())
        } else {
            let mut out: Vec<String> = Vec::new();
            for f in filter_str.split(';') {
                if !out.iter().any(|x| x == f) {
                    out.push(f.to_string());
                }
            }
            Some(out)
        };

        let mut info = Attributes::default();
        if info_str != "." && !info_str.is_empty() {
            for field in info_str.split(';') {
                if field.is_empty() {
                    continue;
                }
                match field.split_once('=') {
                    Some((key, value)) => {
                        if !value.contains(',') && value == "0" && header.info_type(key).as_deref() == Some("Flag") {
                            continue;
                        }
                        info.put(key, Value::Raw(value.to_string()));
                    }
                    None => {
                        let is_flag = header.info_type(field).is_none_or(|t| t == "Flag");
                        info.put(field, if is_flag { Value::Bool(true) } else { Value::Raw(".".into()) });
                    }
                }
            }
        }

        let alts = if alt_str == "." { Vec::new() } else { alt_str.split(',').map(normalise_alt).collect() };
        Ok(VcfRecord {
            chrom,
            pos,
            id,
            ref_allele: ref_str.to_ascii_uppercase(),
            alts,
            log10_p_error,
            filters,
            info,
            genotypes: match rest {
                Some(text) if !header.samples.is_empty() => Genotypes::Lazy(text.to_string()),
                _ => Genotypes::Decoded(Vec::new()),
            },
        })
    }

    pub fn alt(&self) -> &str {
        self.alts.first().map(String::as_str).unwrap_or("")
    }

    /// VariantContext.getPhredScaledQual
    pub fn phred_scaled_qual(&self) -> f64 {
        (self.log10_p_error * -10.0) + 0.0
    }

    pub fn is_filtered(&self) -> bool {
        self.filters.as_ref().is_some_and(|f| !f.is_empty())
    }

    /// VariantContextDecorator.displayFilter: the filters joined by ';' in their order, or PASS
    pub fn display_filter(&self) -> String {
        match &self.filters {
            Some(f) if !f.is_empty() => f.join(";"),
            _ => "PASS".to_string(),
        }
    }

    /// VariantContext.getType (determineType)
    pub fn variant_type(&self) -> VariantType {
        if self.alts.is_empty() {
            return VariantType::NoVariation;
        }
        let mut t: Option<VariantType> = None;
        for alt in &self.alts {
            let this = if is_symbolic_allele(alt) {
                VariantType::Symbolic
            } else if self.ref_allele.len() == alt.len() {
                if alt.len() == 1 { VariantType::Snp } else { VariantType::Mnp }
            } else {
                VariantType::Indel
            };
            match t {
                None => t = Some(this),
                Some(prev) if prev != this => return VariantType::Mixed,
                _ => {}
            }
        }
        t.unwrap()
    }

    /// decodes the genotypes if they are still unparsed (any genotype access in Java does)
    pub fn decode(&mut self, header: &VcfHeader) -> Result<()> {
        if let Genotypes::Lazy(text) = &self.genotypes {
            let g = decode_genotypes(text, &header.samples).with_context(|| format!("{}:{}", self.chrom, self.pos))?;
            self.genotypes = Genotypes::Decoded(g);
        }
        Ok(())
    }

    pub fn genotypes(&mut self, header: &VcfHeader) -> &mut Vec<Genotype> {
        self.decode(header).expect("invalid genotype data");
        match &mut self.genotypes {
            Genotypes::Decoded(g) => g,
            Genotypes::Lazy(_) => unreachable!(),
        }
    }

    /// VariantContext.getGenotype(sample) (decodes the genotypes)
    pub fn genotype(&mut self, header: &VcfHeader, sample: &str) -> Option<&mut Genotype> {
        self.genotypes(header).iter_mut().find(|g| g.sample == sample)
    }

    /// VariantContext.getGenotype(int) if the genotypes are decoded (no decoding)
    pub fn decoded_genotype_at(&self, index: usize) -> Option<&Genotype> {
        match &self.genotypes {
            Genotypes::Decoded(g) => g.get(index),
            Genotypes::Lazy(_) => None,
        }
    }

    /// the genotype of a sample if the genotypes are decoded (no decoding)
    pub fn decoded_genotype(&self, sample: &str) -> Option<&Genotype> {
        match &self.genotypes {
            Genotypes::Decoded(g) => g.iter().find(|g| g.sample == sample),
            Genotypes::Lazy(_) => None,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// header

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Count {
    Fixed(i32),
    A,
    R,
    G,
    Unbounded,
}

#[derive(Clone, Debug, Default)]
pub struct VcfHeader {
    /// the `##` lines in input order, without the file format line
    pub meta_lines: Vec<String>,
    pub samples: Vec<String>,
    info: std::collections::HashMap<String, (Count, String)>,
    format: std::collections::HashMap<String, (Count, String)>,
    filters: std::collections::HashSet<String>,
}

/// (ID, Number, Type) of a `##INFO` or `##FORMAT` line
fn compound_fields(line: &str) -> Option<(String, String, String)> {
    let body = line.split_once("=<")?.1;
    let mut id = None;
    let mut number = None;
    let mut ty = None;
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

fn line_id(line: &str) -> Option<String> {
    let body = line.split_once("=<")?.1;
    let id = body.strip_prefix("ID=")?;
    let end = id.find([',', '>']).unwrap_or(id.len());
    Some(id[..end].to_string())
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

/// VCFInfoHeaderLine(id, count, type, description).toString()
pub fn info_line(id: &str, number: &str, ty: &str, desc: &str) -> String {
    format!("##INFO=<ID={id},Number={number},Type={ty},Description=\"{desc}\">")
}
/// VCFFormatHeaderLine(...).toString()
pub fn format_line(id: &str, number: &str, ty: &str, desc: &str) -> String {
    format!("##FORMAT=<ID={id},Number={number},Type={ty},Description=\"{desc}\">")
}
/// VCFFilterHeaderLine(id, description).toString()
pub fn filter_line(id: &str, desc: &str) -> String {
    format!("##FILTER=<ID={id},Description=\"{desc}\">")
}
/// VCFHeaderLine(key, value).toString()
pub fn meta_line(key: &str, value: &str) -> String {
    format!("##{key}={value}")
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
    /// VCFHeader(Set<VCFHeaderLine>, samples) of the lines read: every line kept, the lookup maps
    /// holding the first line of each ID
    pub fn from_lines(meta: Vec<String>, samples: Vec<String>) -> VcfHeader {
        let mut header = VcfHeader { samples, ..Default::default() };
        // the lines kept so far, for the duplicate test (a linear search made this quadratic in the
        // thousands of ##contig lines of a GRCh38 header)
        let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        for line in meta {
            if line.starts_with("##fileformat=") {
                continue;
            }
            let line = repair_standard_line(&line);
            header.bind(&line);
            if seen.insert(line.clone()) {
                header.meta_lines.push(line);
            }
        }
        header
    }

    /// addMetadataLineLookupEntry: false where an INFO, FORMAT, FILTER or contig line of the ID is bound
    fn bind(&mut self, line: &str) -> bool {
        if line.starts_with("##INFO=<") {
            if let Some((id, number, ty)) = compound_fields(line) {
                if self.info.contains_key(&id) {
                    return false;
                }
                self.info.insert(id, (parse_count(&number), ty));
            }
        } else if line.starts_with("##FORMAT=<") {
            if let Some((id, number, ty)) = compound_fields(line) {
                if self.format.contains_key(&id) {
                    return false;
                }
                self.format.insert(id, (parse_count(&number), ty));
            }
        } else if line.starts_with("##FILTER=<") {
            if let Some(id) = line_id(line) {
                if !self.filters.insert(id) {
                    return false;
                }
            }
        } else if line.starts_with("##contig=<") {
            if let Some(id) = line_id(line) {
                if self.meta_lines.iter().any(|l| l.starts_with("##contig=<") && line_id(l).as_deref() == Some(id.as_str())) {
                    return false;
                }
            }
        }
        true
    }

    /// VCFHeader.addMetaDataLine: an INFO, FORMAT, FILTER or contig line whose ID is already declared is
    /// not added (the first is kept); other lines are added unless an equal line is present
    pub fn add_line(&mut self, line: String) {
        if self.bind(&line) && !self.meta_lines.contains(&line) {
            self.meta_lines.push(line);
        }
    }

    pub fn has_info(&self, id: &str) -> bool {
        self.info.contains_key(id)
    }
    pub fn info_count(&self, id: &str) -> Option<&Count> {
        self.info.get(id).map(|(c, _)| c)
    }
    pub fn info_type(&self, id: &str) -> Option<String> {
        self.info.get(id).map(|(_, t)| t.clone())
    }
    pub fn has_format(&self, id: &str) -> bool {
        self.format.contains_key(id)
    }
    pub fn has_filter(&self, id: &str) -> bool {
        self.filters.contains(id)
    }

    /// VCFHeader.getMetaDataLine(key).getValue(): the first `##key=value` line
    pub fn meta_value(&self, key: &str) -> Option<String> {
        let prefix = format!("##{key}=");
        self.meta_lines.iter().find_map(|l| l.strip_prefix(&prefix).map(|v| v.to_string()))
    }

    /// the contig lines' IDs in order
    pub fn contigs(&self) -> Vec<String> {
        self.meta_lines.iter().filter(|l| l.starts_with("##contig=<")).filter_map(|l| line_id(l)).collect()
    }

    /// Header lines in the order htsjdk writes them: the file format first, then the rest sorted as
    /// strings without duplicates, with the contig lines kept together in their own order.
    pub fn sorted_lines(&self) -> Vec<String> {
        let mut contigs: Vec<&String> = Vec::new();
        let mut others: Vec<&String> = Vec::new();
        for line in &self.meta_lines {
            if line.starts_with("##contig=") {
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

    pub fn text(&self) -> String {
        let mut text = self.sorted_lines().join("\n");
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

// ---------------------------------------------------------------------------------------------
// reading

/// A VCF file read record by record (VcfFileReader / AbstractFeatureReader with VCFCodec).
pub struct VcfReader {
    pub header: VcfHeader,
    reader: Box<dyn BufRead>,
    line: String,
    path: String,
}

impl VcfReader {
    pub fn open(path: &str) -> Result<VcfReader> {
        let raw: Box<dyn std::io::Read> = if path.ends_with(".gz") || path.ends_with(".bgz") {
            Box::new(rust_htslib::bgzf::Reader::from_path(path).map_err(|e| anyhow!("opening VCF {path}: {e}"))?)
        } else {
            Box::new(std::fs::File::open(path).with_context(|| format!("opening VCF {path}"))?)
        };
        let mut reader: Box<dyn BufRead> = Box::new(std::io::BufReader::with_capacity(1 << 17, raw));
        let mut meta = Vec::new();
        let mut samples = Vec::new();
        let mut line = String::new();
        let mut has_column_line = false;
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let text = line.trim_end_matches(['\n', '\r']);
            if text.starts_with("##") {
                meta.push(text.to_string());
            } else if text.starts_with('#') {
                samples = text.split('\t').skip(9).map(|s| s.to_string()).collect();
                has_column_line = true;
                break;
            } else {
                break;
            }
        }
        if !has_column_line {
            bail!("invalid VCF({path}): no header line");
        }
        Ok(VcfReader { header: VcfHeader::from_lines(meta, samples), reader, line: String::new(), path: path.to_string() })
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        loop {
            self.line.clear();
            if self.reader.read_line(&mut self.line).with_context(|| format!("reading {}", self.path))? == 0 {
                return Ok(None);
            }
            let text = self.line.trim_end_matches(['\n', '\r']);
            if text.is_empty() {
                continue;
            }
            return VcfRecord::parse(text, &self.header).map(Some);
        }
    }

    /// The next data line as read, without its line end and without parsing it (empty lines skipped), for
    /// a reader that looks at a few columns before parsing a record with `VcfRecord::parse`. (Added
    /// 2026-11-27 for SIGS, which skips filtered records unparsed.)
    pub fn next_line(&mut self) -> Result<Option<&str>> {
        loop {
            self.line.clear();
            if self.reader.read_line(&mut self.line).with_context(|| format!("reading {}", self.path))? == 0 {
                return Ok(None);
            }
            let len = self.line.trim_end_matches(['\n', '\r']).len();
            if len == 0 {
                continue;
            }
            return Ok(Some(&self.line[..len]));
        }
    }

    /// every record
    pub fn read_all(&mut self) -> Result<Vec<VcfRecord>> {
        let mut out = Vec::new();
        while let Some(r) = self.next_record()? {
            out.push(r);
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------------------------
// writing

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

/// VCFEncoder options
#[derive(Clone, Copy, Debug, Default)]
pub struct EncodeOptions {
    /// Options.ALLOW_MISSING_FIELDS_IN_HEADER
    pub allow_missing_fields: bool,
}

/// VCFEncoder.encode of one record, without the line end
pub fn encode_record(header: &VcfHeader, opts: EncodeOptions, rec: &VcfRecord) -> Result<String> {
    let mut s = String::with_capacity(512);
    s.push_str(&rec.chrom);
    s.push('\t');
    s.push_str(&rec.pos.to_string());
    s.push('\t');
    s.push_str(&rec.id);
    s.push('\t');
    s.push_str(&rec.ref_allele);
    s.push('\t');
    if rec.alts.is_empty() {
        s.push('.');
    } else {
        s.push_str(&rec.alts.join(","));
    }
    s.push('\t');
    if rec.log10_p_error == NO_LOG10_PERROR {
        s.push('.');
    } else {
        s.push_str(&format_qual(rec.phred_scaled_qual()));
    }
    s.push('\t');
    match &rec.filters {
        None => s.push('.'),
        Some(f) if f.is_empty() => s.push_str("PASS"),
        Some(f) => {
            if !opts.allow_missing_fields {
                for filter in f {
                    if !header.has_filter(filter) {
                        bail!(
                            "Key {filter} found in VariantContext field FILTER at {}:{} but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.",
                            rec.chrom,
                            rec.pos
                        );
                    }
                }
            }
            let mut sorted: Vec<&String> = f.iter().collect();
            sorted.sort();
            for (i, x) in sorted.iter().enumerate() {
                if i > 0 {
                    s.push(';');
                }
                s.push_str(x);
            }
        }
    }
    s.push('\t');

    // INFO: a TreeMap of key -> formatted value, values Java formats as null left out
    let mut fields: Vec<(&str, String)> = Vec::with_capacity(rec.info.0.len());
    for (key, value) in &rec.info.0 {
        if !opts.allow_missing_fields && !header.has_info(key) {
            bail!(
                "Key {key} found in VariantContext field INFO at {}:{} but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.",
                rec.chrom,
                rec.pos
            );
        }
        if let Some(v) = value.encode() {
            fields.push((key, v));
        }
    }
    fields.sort_by(|a, b| a.0.cmp(b.0));
    if fields.is_empty() {
        s.push('.');
    }
    for (i, (key, value)) in fields.iter().enumerate() {
        if i > 0 {
            s.push(';');
        }
        s.push_str(key);
        if !value.is_empty() && header.info_count(key) != Some(&Count::Fixed(0)) {
            s.push('=');
            s.push_str(value);
        }
    }

    if !header.samples.is_empty() {
        match &rec.genotypes {
            Genotypes::Lazy(text) => {
                s.push('\t');
                s.push_str(text);
            }
            Genotypes::Decoded(genotypes) => encode_genotypes(header, rec, genotypes, &mut s)?,
        }
    }
    Ok(s)
}

fn join_ints(v: &[i32]) -> String {
    v.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",")
}

/// VariantContext.calcVCFGenotypeKeys and VCFEncoder.appendGenotypeData
fn encode_genotypes(header: &VcfHeader, rec: &VcfRecord, genotypes: &[Genotype], s: &mut String) -> Result<()> {
    let mut keys: Vec<&str> = Vec::new();
    let (mut saw_gt, mut saw_gq, mut saw_dp, mut saw_ad, mut saw_pl, mut saw_ft) = (false, false, false, false, false, false);
    for g in genotypes {
        for (k, _) in &g.attrs.0 {
            if !keys.contains(&k.as_str()) {
                keys.push(k);
            }
        }
        saw_gt |= !g.alleles.is_empty();
        saw_gq |= g.gq.is_some();
        saw_dp |= g.dp.is_some();
        saw_ad |= g.ad.is_some();
        saw_pl |= g.pl.is_some();
        saw_ft |= g.filters.is_some();
    }
    for (saw, key) in [(saw_gq, "GQ"), (saw_dp, "DP"), (saw_ad, "AD"), (saw_pl, "PL"), (saw_ft, "FT")] {
        if saw && !keys.contains(&key) {
            keys.push(key);
        }
    }
    keys.sort();
    if saw_gt {
        keys.insert(0, "GT");
    }
    if keys.is_empty() {
        // all samples are no-calls
        keys.push("GT");
    }
    s.push('\t');
    s.push_str(&keys.join(":"));

    let missing = Genotype { alleles: vec![None, None], ..Default::default() };
    for sample in &header.samples {
        s.push('\t');
        let g = genotypes.iter().find(|g| &g.sample == sample).unwrap_or(&missing);
        let mut values: Vec<String> = Vec::with_capacity(keys.len());
        for key in &keys {
            match *key {
                "GT" => {}
                "FT" => values.push(g.filters.clone().unwrap_or_else(|| "PASS".to_string())),
                "DP" => values.push(g.dp.map(|v| v.to_string()).unwrap_or_else(|| ".".into())),
                "GQ" => values.push(g.gq.map(|v| v.to_string()).unwrap_or_else(|| ".".into())),
                "AD" => values.push(g.ad.as_ref().map(|v| join_ints(v)).unwrap_or_else(|| ".".into())),
                "PL" => values.push(g.pl.as_ref().map(|v| join_ints(v)).unwrap_or_else(|| ".".into())),
                _ => match g.attrs.get(key).and_then(|v| v.encode()) {
                    Some(v) => values.push(v),
                    None => {
                        let n = match header.format.get(*key).map(|(c, _)| c) {
                            Some(Count::Fixed(n)) => *n,
                            Some(Count::A) => rec.alts.len() as i32,
                            Some(Count::R) => rec.alts.len() as i32 + 1,
                            Some(Count::G) => {
                                let alleles = rec.alts.len() as i32 + 1;
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
        if saw_gt || keys == ["GT"] {
            if g.alleles.is_empty() {
                s.push('.');
            } else {
                let sep = if g.phased { "|" } else { "/" };
                for (i, a) in g.alleles.iter().enumerate() {
                    if i > 0 {
                        s.push_str(sep);
                    }
                    match a {
                        Some(n) => s.push_str(&n.to_string()),
                        None => s.push('.'),
                    }
                }
            }
        }
        let has_gt = keys.first() == Some(&"GT");
        for (i, v) in values.iter().enumerate() {
            if i > 0 || has_gt {
                s.push(':');
            }
            s.push_str(v);
        }
    }
    Ok(())
}

/// A VCF as htsjdk's VariantContextWriter writes it with BLOCK_COMPRESSED_VCF and INDEX_ON_THE_FLY: the BGZF
/// stream and the tabix index (`<path>.tbi`) of `crate::vcf::VcfFileWriter`, so the compressed bytes and the
/// index are Java's too; a plain VCF for a path without ".gz".
pub struct VcfWriter {
    out: Option<crate::vcf::VcfFileWriter>,
    plain: Option<std::io::BufWriter<std::fs::File>>,
    header: VcfHeader,
    opts: EncodeOptions,
}

impl VcfWriter {
    pub fn create(path: &str, header: &VcfHeader, opts: EncodeOptions) -> Result<VcfWriter> {
        let (out, plain) = if path.ends_with(".gz") {
            (Some(crate::vcf::VcfFileWriter::create_with_header_text(path, &header.text())?), None)
        } else {
            let mut p = std::io::BufWriter::with_capacity(1 << 16, std::fs::File::create(path).with_context(|| format!("creating {path}"))?);
            p.write_all(header.text().as_bytes())?;
            (None, Some(p))
        };
        Ok(VcfWriter { out, plain, header: header.clone(), opts })
    }

    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    pub fn write(&mut self, rec: &VcfRecord) -> Result<()> {
        let line = encode_record(&self.header, self.opts, rec)?;
        if let Some(o) = self.out.as_mut() {
            o.write_line(&rec.chrom, rec.pos, rec.end(), &line)?;
        } else if let Some(p) = self.plain.as_mut() {
            p.write_all(line.as_bytes())?;
            p.write_all(b"\n")?;
        }
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        if let Some(mut p) = self.plain.take() {
            p.flush()?;
        }
        if let Some(o) = self.out.take() {
            o.close()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> VcfHeader {
        VcfHeader::from_lines(
            vec![
                "##fileformat=VCFv4.2".into(),
                info_line("F", "0", "Flag", "a flag"),
                info_line("X", "1", "Float", "x"),
                info_line("L", ".", "Integer", "list"),
                format_line("GT", "1", "String", "Genotype"),
                format_line("AD", "R", "Integer", "Allelic depths for the ref and alt alleles in the order listed"),
                format_line("Q2", "2", "Integer", "two"),
                format_line("Z", "1", "String", "z"),
                filter_line("PASS", "All filters passed"),
                filter_line("lowQ", "low"),
                "##contig=<ID=chr1,length=100>".into(),
            ],
            vec!["N".into(), "T".into()],
        )
    }

    #[test]
    fn round_trip_lazy_and_decoded() {
        let h = header();
        let line = "chr1\t10\t.\tA\tc\t37.50\tPASS\tX=0.5;F;L=1,2\tGT:AD:Q2:Z\t0/0:10,0:.:.\t0/1:5,6:1,2:a";
        let mut rec = VcfRecord::parse(line, &h).unwrap();
        assert_eq!(rec.alts, vec!["C"]);
        assert_eq!(rec.info.as_int_list("L", 0), vec![1, 2]);
        assert!(rec.info.as_bool("F", false));
        assert_eq!(rec.info.as_string("L", ""), "[1, 2]");
        let enc = encode_record(&h, EncodeOptions::default(), &rec).unwrap();
        assert_eq!(enc, "chr1\t10\t.\tA\tC\t37.50\tPASS\tF;L=1,2;X=0.5\tGT:AD:Q2:Z\t0/0:10,0:.:.\t0/1:5,6:1,2:a");
        rec.genotype(&h, "T").unwrap().alleles = vec![Some(1), Some(1)];
        rec.info.put("X", Value::Double(0.123456));
        rec.info.put("G", Value::Bool(false));
        let enc = encode_record(&h, EncodeOptions { allow_missing_fields: true }, &rec).unwrap();
        assert_eq!(enc, "chr1\t10\t.\tA\tC\t37.50\tPASS\tF;L=1,2;X=0.123\tGT:AD:Q2:Z\t0/0:10,0\t1/1:5,6:1,2:a");
        assert!(encode_record(&h, EncodeOptions::default(), &{
            let mut r = rec.clone();
            r.info.put("NEW", Value::Int(1));
            r
        })
        .is_err());
    }

    #[test]
    fn values() {
        assert_eq!(Value::Double(0.0).encode().unwrap(), "0.00");
        assert_eq!(Value::Double(0.005).encode().unwrap(), "5.000e-03");
        assert_eq!(Value::Double(1.5).encode().unwrap(), "1.50");
        assert_eq!(Value::StrList(vec![]).encode().unwrap(), ".");
        assert_eq!(format_qual(200.0), "200");
        assert_eq!(format_qual(12.345), "12.35");
    }

    #[test]
    fn header_order() {
        let mut h = header();
        h.add_line(meta_line("purpleVersion", "4.5"));
        let lines = h.sorted_lines();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert!(lines.last().unwrap().starts_with("##purpleVersion"));
        let contig = lines.iter().position(|l| l.starts_with("##contig")).unwrap();
        assert!(lines[contig - 1].starts_with("##INFO"));
    }
}
