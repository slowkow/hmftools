//! VCF records and headers as htsjdk 2.23 reads and writes them, limited to what the caller needs.
//!
//! A record keeps its fields as text. The accessors follow VariantContext.getAttributeAs* and the
//! writer follows VCFEncoder: INFO keys sorted, flags without a value, FILTER sorted, FORMAT keys
//! sorted after GT, trailing missing genotype fields dropped.

use anyhow::{bail, Context, Result};
use std::io::{BufRead, Read, Write};

/// A gzip file as a stream of its decompressed bytes, member after member (as GZIPInputStream reads
/// it). Inflating is the slowest part of loading the reference data; zlib-rs does it about twice as
/// fast as flate2's default backend.
pub struct GzReader {
    file: std::fs::File,
    input: Vec<u8>,
    pos: usize,
    len: usize,
    eof: bool,
    inflater: Option<zlib_rs::Inflate>,
}

impl GzReader {
    pub fn new(file: std::fs::File) -> GzReader {
        GzReader { file, input: vec![0; 1 << 17], pos: 0, len: 0, eof: false, inflater: None }
    }

    fn fill(&mut self) -> std::io::Result<()> {
        if self.pos < self.len {
            self.input.copy_within(self.pos..self.len, 0);
        }
        self.len -= self.pos;
        self.pos = 0;
        while !self.eof && self.len < self.input.len() {
            match self.file.read(&mut self.input[self.len..]) {
                Ok(0) => self.eof = true,
                Ok(n) => self.len += n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
}

impl Read for GzReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            if self.len - self.pos < 1 << 15 && !self.eof {
                self.fill()?;
            }
            if self.inflater.is_none() {
                // the next member, if what follows is one
                let rest = &self.input[self.pos..self.len];
                if rest.len() < 2 || rest[0] != 0x1f || rest[1] != 0x8b {
                    return Ok(0);
                }
                self.inflater = Some(zlib_rs::Inflate::new(true, 31));
            }
            let inflater = self.inflater.as_mut().unwrap();
            let (before_in, before_out) = (inflater.total_in(), inflater.total_out());
            let status = inflater
                .decompress(&self.input[self.pos..self.len], out, zlib_rs::InflateFlush::NoFlush)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("invalid gzip data: {e:?}")))?;
            let consumed = (inflater.total_in() - before_in) as usize;
            let produced = (inflater.total_out() - before_out) as usize;
            self.pos += consumed;
            if status == zlib_rs::Status::StreamEnd {
                self.inflater = None;
            } else if consumed == 0 && produced == 0 && self.eof {
                return Err(std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Unexpected end of ZLIB input stream"));
            }
            if produced > 0 {
                return Ok(produced);
            }
        }
    }
}

/// A file as text: gzip and bgzip files (by extension, as in Java) are decompressed.
pub fn open_text(path: &str) -> std::io::Result<Box<dyn Read + Send>> {
    let file = std::fs::File::open(path)?;
    if path.ends_with(".gz") || path.ends_with(".bgz") { Ok(Box::new(GzReader::new(file))) } else { Ok(Box::new(file)) }
}

/// VariantContext.NO_LOG10_PERROR
pub const NO_LOG10_PERROR: f64 = 1.0;

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
            Some(v) => java_parse_int(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\"")),
            None => default,
        }
    }
    /// CommonVcfTags.getGenotypeAttributeAsDouble
    pub fn attr_double(&self, key: &str, default: f64) -> f64 {
        match self.attr(key) {
            Some(v) => java_parse_double(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\"")),
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
    /// INFO in input order; None is a flag
    pub info: Vec<(String, Option<String>)>,
    pub genotypes: Vec<Genotype>,
}

/// Integer.parseInt
pub fn java_parse_int(s: &str) -> Option<i32> {
    let b = s.as_bytes();
    if b.is_empty() {
        return None;
    }
    let digits = if b[0] == b'-' || b[0] == b'+' { &b[1..] } else { b };
    if digits.is_empty() || !digits.iter().all(|c| c.is_ascii_digit()) {
        return None;
    }
    s.strip_prefix('+').unwrap_or(s).parse::<i32>().ok()
}

/// Double.parseDouble for the decimal forms found in VCFs
pub fn java_parse_double(s: &str) -> Option<f64> {
    let t = s.trim();
    match t {
        "NaN" | "+NaN" | "-NaN" => return Some(f64::NAN),
        "Infinity" | "+Infinity" => return Some(f64::INFINITY),
        "-Infinity" => return Some(f64::NEG_INFINITY),
        _ => {}
    }
    let t = t.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(t);
    if t.is_empty() || !t.bytes().all(|c| c.is_ascii_digit() || matches!(c, b'.' | b'-' | b'+' | b'e' | b'E')) {
        return None;
    }
    t.parse::<f64>().ok()
}

impl VcfRecord {
    #[cfg(test)]
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

    /// AbstractVCFCodec: one tab-separated record line
    pub fn parse(line: &str, samples: &[String]) -> Result<VcfRecord> {
        let mut cols = line.split('\t');
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

        let log10_p_error = if qual_str == "." {
            NO_LOG10_PERROR
        } else {
            let val = java_parse_double(qual_str).with_context(|| format!("invalid QUAL: {qual_str}"))?;
            if val < 0.0 && (val - -1.0).abs() < 0.1 { NO_LOG10_PERROR } else { val / -10.0 }
        };

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

        let mut info: Vec<(String, Option<String>)> = Vec::new();
        if info_str != "." && !info_str.is_empty() {
            for field in info_str.split(';') {
                if field.is_empty() {
                    continue;
                }
                match field.split_once('=') {
                    Some((key, value)) => {
                        let value = if value.is_empty() { "." } else { value };
                        info.push((key.to_string(), Some(value.to_string())));
                    }
                    None => info.push((field.to_string(), None)),
                }
            }
        }

        let mut genotypes: Vec<Genotype> = Vec::new();
        if let Some(format) = cols.next() {
            let keys: Vec<&str> = format.split(':').collect();
            for (i, sample) in samples.iter().enumerate() {
                let Some(text) = cols.next() else {
                    bail!("VCF record has fewer genotype columns than samples ({} of {}): {chrom}:{pos}", i, samples.len());
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
                genotypes.push(g);
            }
        }

        Ok(VcfRecord {
            chrom,
            pos,
            id,
            ref_allele: ref_str.to_ascii_uppercase(),
            alts: alt_str.split(',').map(normalise_alt).collect(),
            log10_p_error,
            filters,
            info,
            genotypes,
        })
    }

    /// VariantContext.getPhredScaledQual
    pub fn phred_scaled_qual(&self) -> f64 {
        (self.log10_p_error * -10.0) + 0.0
    }

    pub fn alt(&self) -> &str {
        self.alts.first().map(String::as_str).unwrap_or("")
    }

    pub fn has_attr(&self, key: &str) -> bool {
        self.info.iter().any(|(k, _)| k == key)
    }

    fn raw_attr(&self, key: &str) -> Option<&Option<String>> {
        self.info.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    pub fn set_attr(&mut self, key: &str, value: Option<String>) {
        match self.info.iter_mut().find(|(k, _)| k == key) {
            Some(entry) => entry.1 = value,
            None => self.info.push((key.to_string(), value)),
        }
    }

    /// getAttributeAsString: a value with commas is a list and prints as `[a, b]`
    pub fn attr_string(&self, key: &str, default: &str) -> String {
        match self.raw_attr(key) {
            None => default.to_string(),
            Some(None) => "true".to_string(),
            Some(Some(v)) => {
                if v.contains(',') {
                    format!("[{}]", v.split(',').collect::<Vec<_>>().join(", "))
                } else {
                    v.clone()
                }
            }
        }
    }

    /// getAttributeAsString(key, null)
    pub fn attr_string_opt(&self, key: &str) -> Option<String> {
        if self.has_attr(key) { Some(self.attr_string(key, "")) } else { None }
    }

    /// getAttributeAsInt
    pub fn attr_int(&self, key: &str, default: i32) -> i32 {
        match self.raw_attr(key) {
            None => default,
            Some(Some(v)) if v == "." => default,
            Some(Some(v)) => java_parse_int(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\" ({key})")),
            Some(None) => panic!("ClassCastException: flag {key} read as an integer"),
        }
    }

    /// getAttributeAsDouble
    pub fn attr_double(&self, key: &str, default: f64) -> f64 {
        match self.raw_attr(key) {
            None => default,
            Some(Some(v)) if v == "." => default,
            Some(Some(v)) => java_parse_double(v).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{v}\" ({key})")),
            Some(None) => panic!("ClassCastException: flag {key} read as a decimal"),
        }
    }

    /// getAttributeAsIntList
    pub fn attr_int_list(&self, key: &str, default: i32) -> Vec<i32> {
        match self.raw_attr(key) {
            None => Vec::new(),
            Some(None) => panic!("ClassCastException: flag {key} read as an integer list"),
            Some(Some(v)) => v
                .split(',')
                .map(|x| {
                    if x == "." {
                        default
                    } else {
                        java_parse_int(x).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{x}\" ({key})"))
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
}

/// Allele.create: bases are upper-cased unless the allele is symbolic
fn normalise_alt(alt: &str) -> String {
    let b = alt.as_bytes();
    let symbolic = b.len() > 1
        && (b[0] == b'<' || b[b.len() - 1] == b'>' || alt.contains('[') || alt.contains(']') || b[0] == b'.' || b[b.len() - 1] == b'.');
    if symbolic { alt.to_string() } else { alt.to_ascii_uppercase() }
}

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

#[derive(Default)]
pub struct HeaderCounts {
    pub info: hashbrown::HashMap<String, Count>,
    pub format: hashbrown::HashMap<String, Count>,
}

#[derive(Clone, Debug)]
pub struct VcfHeader {
    /// the `##` lines in input order, without the file format line
    pub meta_lines: Vec<String>,
    pub samples: Vec<String>,
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

pub fn info_line(id: &str, number: &str, ty: &str, desc: &str) -> String {
    format!("##INFO=<ID={id},Number={number},Type={ty},Description=\"{desc}\">")
}
pub fn format_line(id: &str, number: &str, ty: &str, desc: &str) -> String {
    format!("##FORMAT=<ID={id},Number={number},Type={ty},Description=\"{desc}\">")
}
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
        let meta_lines = meta.iter().filter(|l| !l.starts_with("##fileformat=")).map(|l| repair_standard_line(l)).collect();
        VcfHeader { meta_lines, samples }
    }

    pub fn has_format_line(&self, id: &str) -> bool {
        self.compound_line("##FORMAT=<", id).is_some()
    }

    /// the last line declaring the key, as the header's map keeps the last one added
    fn compound_line(&self, prefix: &str, id: &str) -> Option<(Count, String)> {
        self.meta_lines.iter().rev().filter(|l| l.starts_with(prefix)).find_map(|l| {
            let (line_id, number, ty) = compound_fields(l)?;
            if line_id == id { Some((parse_count(&number), ty)) } else { None }
        })
    }

    /// the declared count of every INFO and FORMAT key
    pub fn counts(&self) -> HeaderCounts {
        let mut counts = HeaderCounts::default();
        for line in &self.meta_lines {
            let map = if line.starts_with("##INFO=<") {
                &mut counts.info
            } else if line.starts_with("##FORMAT=<") {
                &mut counts.format
            } else {
                continue;
            };
            if let Some((id, number, _)) = compound_fields(line) {
                map.insert(id, parse_count(&number));
            }
        }
        counts
    }

    pub fn add_meta_line(&mut self, line: String) {
        if !self.meta_lines.contains(&line) {
            self.meta_lines.push(line);
        }
    }
}

pub struct VcfInput {
    pub header: VcfHeader,
    reader: Box<dyn BufRead>,
    line: String,
}

impl VcfInput {
    pub fn open(path: &str) -> Result<VcfInput> {
        let file = open_text(path).with_context(|| format!("opening VCF {path}"))?;
        let mut reader: Box<dyn BufRead> = Box::new(std::io::BufReader::with_capacity(1 << 16, file));
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
        Ok(VcfInput { header: VcfHeader::from_lines(meta, samples), reader, line: String::new() })
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        loop {
            self.line.clear();
            if self.reader.read_line(&mut self.line)? == 0 {
                return Ok(None);
            }
            let text = self.line.trim_end_matches(['\n', '\r']);
            if text.is_empty() {
                continue;
            }
            return VcfRecord::parse(text, &self.header.samples).map(Some);
        }
    }
}

// ---------------------------------------------------------------------------------------------
// writing

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

pub fn header_text(header: &VcfHeader) -> String {
    let mut text = sorted_header_lines(&header.meta_lines).join("\n");
    text.push_str("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !header.samples.is_empty() {
        text.push_str("\tFORMAT");
        for sample in &header.samples {
            text.push('\t');
            text.push_str(sample);
        }
    }
    text.push('\n');
    text
}

/// VCFEncoder.formatQualValue
pub fn format_qual(phred: f64) -> String {
    let s = crate::common::fmt::java_fixed(phred, 2);
    match s.strip_suffix(".00") {
        Some(t) => t.to_string(),
        None => s,
    }
}

/// VCFEncoder.isMissingValue
fn is_missing_value(s: &str) -> bool {
    s.bytes().all(|c| c == b'.' || c == b',')
}

/// A record split into the columns before FILTER and those after it, so that the same record can be
/// written with several filter sets.
pub struct EncodedRecord {
    pub prefix: String,
    pub suffix: String,
}

/// VCFEncoder.encode without the FILTER column. `genotypes` are looked up by header sample name.
pub fn encode_record(header: &VcfHeader, counts: &HeaderCounts, rec: &VcfRecord, log10_p_error: f64, genotypes: &[Genotype]) -> Result<EncodedRecord> {
    let mut prefix = String::with_capacity(256);
    prefix.push_str(&rec.chrom);
    prefix.push('\t');
    prefix.push_str(&rec.pos.to_string());
    prefix.push('\t');
    prefix.push_str(&rec.id);
    prefix.push('\t');
    prefix.push_str(&rec.ref_allele);
    prefix.push('\t');
    if rec.alts.is_empty() {
        prefix.push('.');
    } else {
        prefix.push_str(&rec.alts.join(","));
    }
    prefix.push('\t');
    if log10_p_error == NO_LOG10_PERROR {
        prefix.push('.');
    } else {
        prefix.push_str(&format_qual((log10_p_error * -10.0) + 0.0));
    }

    let mut suffix = String::with_capacity(512);
    let mut keys: Vec<&(String, Option<String>)> = rec.info.iter().collect();
    keys.sort_by(|a, b| a.0.cmp(&b.0));
    if keys.is_empty() {
        suffix.push('.');
    }
    for (i, (key, value)) in keys.iter().enumerate() {
        let Some(count) = counts.info.get(key.as_str()) else {
            bail!("Key {key} found in VariantContext field INFO at {}:{} but this key isn't defined in the VCFHeader", rec.chrom, rec.pos);
        };
        if i > 0 {
            suffix.push(';');
        }
        suffix.push_str(key);
        if let Some(v) = value {
            if !v.is_empty() && *count != Count::Fixed(0) {
                suffix.push('=');
                suffix.push_str(v);
            }
        }
    }

    if !header.samples.is_empty() {
        // calcVCFGenotypeKeys
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

        suffix.push('\t');
        suffix.push_str(&format_keys.join(":"));

        let missing = Genotype::default();
        for sample in &header.samples {
            suffix.push('\t');
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
                        Some(v) => match java_parse_double(v) {
                            Some(d) => ((d + 0.5).floor() as i64).to_string(),
                            None => v.clone(),
                        },
                        None => ".".to_string(),
                    }),
                    _ => match g.attr(key) {
                        Some(v) => values.push(v.to_string()),
                        None => {
                            let n = match counts.format.get(*key) {
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
            if has_gt {
                match (&g.gt, absent) {
                    (Some(gt), _) => suffix.push_str(gt),
                    (None, true) => suffix.push_str("./."),
                    (None, false) => bail!("GTs cannot be missing for some samples if they are available for others in the record"),
                }
            }
            for (i, v) in values.iter().enumerate() {
                if i > 0 || has_gt {
                    suffix.push(':');
                }
                suffix.push_str(v);
            }
        }
    }

    Ok(EncodedRecord { prefix, suffix })
}

/// A bgzip-compressed VCF with a tabix index written when it is closed.
pub struct VcfWriter {
    path: String,
    out: rust_htslib::bgzf::Writer,
    buf: Vec<u8>,
}

impl VcfWriter {
    pub fn create(path: &str, header: &VcfHeader) -> Result<VcfWriter> {
        let out = rust_htslib::bgzf::Writer::from_path(path).with_context(|| format!("creating {path}"))?;
        let mut buf = Vec::with_capacity(1 << 16);
        buf.extend_from_slice(header_text(header).as_bytes());
        Ok(VcfWriter { path: path.to_string(), out, buf })
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    /// `filters` are the tags to write; none writes PASS
    pub fn write(&mut self, rec: &EncodedRecord, filters: &mut Vec<&str>) -> Result<()> {
        self.buf.extend_from_slice(rec.prefix.as_bytes());
        self.buf.push(b'\t');
        if filters.is_empty() {
            self.buf.extend_from_slice(b"PASS");
        } else {
            filters.sort();
            for (i, f) in filters.iter().enumerate() {
                if i > 0 {
                    self.buf.push(b';');
                }
                self.buf.extend_from_slice(f.as_bytes());
            }
        }
        self.buf.push(b'\t');
        self.buf.extend_from_slice(rec.suffix.as_bytes());
        self.buf.push(b'\n');
        if self.buf.len() >= 1 << 16 {
            self.out.write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        self.out.write_all(&self.buf)?;
        self.out.flush()?;
        let path = self.path.clone();
        drop(self.out);
        let c_path = std::ffi::CString::new(path.as_str())?;
        let ret = unsafe { rust_htslib::htslib::tbx_index_build(c_path.as_ptr(), 0, &rust_htslib::htslib::tbx_conf_vcf) };
        if ret != 0 {
            bail!("failed to build tabix index for {path}");
        }
        Ok(())
    }
}
