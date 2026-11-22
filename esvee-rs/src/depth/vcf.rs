//! VCF reading and writing that reproduces what htsjdk does to a record which is read with VCFCodec,
//! has attributes added, and is written again with VariantContextWriter (block-compressed, tabix index).
//!
//! The Java step never decodes values with the header's types, so INFO and FORMAT values pass through
//! as text. What does change on the way through htsjdk: QUAL is re-formatted, FILTER names are sorted,
//! INFO and FORMAT keys are sorted, AD/DP/PL/GQ are parsed as integers, trailing missing genotype values
//! are dropped, and the header lines are sorted.

use super::fmt::java_fixed;
use anyhow::{anyhow, bail, Context, Result};
use rust_htslib::bgzf;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{BufRead, Write};

const MISSING: &str = ".";

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
    if def.count != standard_def.count || def.ty != standard_def.ty { standard_line } else { line.to_string() }
}

impl HeaderDefs {
    pub fn add_line(&mut self, line: &str) {
        if let Some((is_info, id, def)) = parse_field_line(line) {
            let map = if is_info { &mut self.info } else { &mut self.format };
            map.entry(id).or_insert(def);
        }
    }
    pub fn has_info(&self, id: &str) -> bool {
        self.info.contains_key(id)
    }
    pub fn has_format(&self, id: &str) -> bool {
        self.format.contains_key(id)
    }
}

/// A VCF header: the `##` lines (standard lines repaired as htsjdk does on reading) and the sample names.
pub struct VcfHeader {
    pub meta_lines: Vec<String>,
    pub samples: Vec<String>,
    /// whether the column line has a FORMAT column
    pub has_genotypes: bool,
    pub defs: HeaderDefs,
}

impl VcfHeader {
    pub fn add_meta_line(&mut self, line: String) {
        self.defs.add_line(&line);
        self.meta_lines.push(line);
    }
    pub fn add_info_line(&mut self, id: &str, number: i32, ty: &str, description: &str) {
        self.add_meta_line(format!("##INFO=<ID={id},Number={number},Type={ty},Description=\"{description}\">"));
    }
    pub fn add_format_line(&mut self, id: &str, number: i32, ty: &str, description: &str) {
        self.add_meta_line(format!("##FORMAT=<ID={id},Number={number},Type={ty},Description=\"{description}\">"));
    }
}

/// A genotype as htsjdk holds it: GT, the fields with their own accessors (AD, DP, GQ, PL), the genotype
/// filter and the remaining ("extended") attributes. Values are stored as they will be written.
#[derive(Clone, Debug, Default)]
pub struct Genotype {
    pub gt: Option<String>,
    pub standard: BTreeMap<String, String>,
    pub filters: Option<String>,
    pub attributes: BTreeMap<String, String>,
}

impl Genotype {
    /// CommonVcfTags.getGenotypeAttributeAsInt: Integer.parseInt of the value, which fails on anything else
    pub fn attribute_as_int(&self, key: &str, default: i32) -> Result<i32> {
        match self.attributes.get(key) {
            None => Ok(default),
            Some(v) => parse_java_int(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\"")),
        }
    }
}

/// Integer.parseInt: an optional sign and decimal digits, within the int range
pub fn parse_java_int(s: &str) -> Option<i32> {
    let digits = s.strip_prefix(['+', '-']).unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse::<i32>().ok()
}

/// A VCF record as htsjdk holds it after reading; every value is kept in the form the encoder writes.
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

fn parse_ints(field: &str, value: &str) -> Result<String> {
    let ints: Result<Vec<String>> = value
        .split(',')
        .map(|x| x.parse::<i32>().map(|v| v.to_string()).map_err(|_| anyhow!("invalid {field} value {value}")))
        .collect();
    Ok(ints?.join(","))
}

/// AbstractVCFCodec.parseFilters; the names are held in a set
fn parse_filters(s: &str) -> Result<Option<Vec<String>>> {
    if s == MISSING {
        Ok(None)
    } else if s == "PASS" {
        Ok(Some(Vec::new()))
    } else if s == "0" {
        bail!("0 is an invalid filter name in vcf4")
    } else if s.is_empty() {
        bail!("The VCF specification requires a valid filter status: filter was {s}")
    } else {
        let mut names: Vec<String> = s.split(';').map(|x| x.to_string()).collect();
        names.sort();
        names.dedup();
        Ok(Some(names))
    }
}

/// VCFEncoder.isMissingValue
fn is_missing_value(s: &str) -> bool {
    s.bytes().all(|b| b == b'.' || b == b',')
}

/// Math.round
fn jround(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

const STANDARD_GENOTYPE_KEYS: [&str; 4] = ["AD", "DP", "GQ", "PL"];

/// Allele.wouldBeSymbolicAllele: `<...>`, breakpoints and single breakends keep their text; other alleles
/// are bases, which htsjdk holds in upper case
fn allele(a: &str) -> String {
    let symbolic = a.len() > 1 && (a.starts_with('<') || a.ends_with('>') || a.contains('[') || a.contains(']') || a.starts_with('.') || a.ends_with('.'));
    if symbolic { a.to_string() } else { a.to_ascii_uppercase() }
}

impl VariantContext {
    /// One data line, read as VCFCodec does (with the genotypes decoded).
    pub fn parse(line: &str, header: &VcfHeader) -> Result<VariantContext> {
        let defs = &header.defs;
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            bail!("invalid VCF record: {line}");
        }
        let position: i32 = cols[1].parse().with_context(|| format!("invalid VCF position: {}", cols[1]))?;
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
            if cols[7].is_empty() || cols[7].contains(' ') {
                bail!("The VCF specification does not allow for whitespace in the INFO field. Offending field value was \"{}\"", cols[7]);
            }
            for field in cols[7].split(';') {
                let info_def = defs.info.get(field.split('=').next().unwrap()).copied();
                let (key, value) = match field.split_once('=') {
                    Some((key, value_string)) => {
                        if !value_string.contains(',') && info_def.is_some_and(|d| d.ty == FieldType::Flag) && value_string == "0" {
                            // a flag written as KEY=0 is not set
                            continue;
                        }
                        // KEY= is read as missing; a list is written as it was read
                        (key, if value_string.is_empty() { MISSING } else { value_string })
                    }
                    // a key without a value is a flag unless the header says otherwise
                    None => (field, if info_def.is_some_and(|d| d.ty != FieldType::Flag) { MISSING } else { "" }),
                };
                info.insert(key.to_string(), value.to_string());
            }
        }

        let mut genotypes = Vec::new();
        if cols.len() > 9 {
            let keys: Vec<&str> = cols[8].split(':').collect();
            if cols.len() - 9 != header.samples.len() {
                bail!("there are {} genotypes while the header requires that {} genotypes be present for all records at {}:{}", cols.len() - 9, header.samples.len(), cols[0], position);
            }
            for (column, sample) in cols[9..].iter().zip(&header.samples) {
                let values: Vec<&str> = column.split(':').collect();
                if keys.len() < values.len() {
                    bail!("there are too many keys for the sample {sample} at {}:{}", cols[0], position);
                }
                let mut g = Genotype::default();
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
                            if let Some(filters) = parse_filters(value)?.filter(|f| !f.is_empty()) {
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
                            g.attributes.insert(key.to_string(), value.to_string());
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
            filters: parse_filters(cols[6])?,
            info,
            genotypes,
        })
    }

    /// getAttributeAsIntList: None when the attribute is absent; a value which is not an integer fails
    pub fn info_ints(&self, key: &str) -> Result<Option<Vec<i32>>> {
        let Some(value) = self.info.get(key) else { return Ok(None) };
        let ints: Option<Vec<i32>> = value.split(',').map(parse_java_int).collect();
        ints.map(Some).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string in {key}={value}"))
    }

    /// The record as htsjdk's VCFEncoder writes it.
    pub fn encode(&self, header: &VcfHeader) -> Result<String> {
        let defs = &header.defs;
        let mut line = String::with_capacity(1024);
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
            Some(f) => line.push_str(&f.join(";")),
        }
        line.push('\t');
        if self.info.is_empty() {
            line.push_str(MISSING);
        } else {
            let mut first = true;
            for (key, value) in &self.info {
                let Some(def) = defs.info.get(key) else {
                    bail!("Key {key} found in VariantContext field INFO at {}:{} but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.", self.chromosome, self.position);
                };
                if !first {
                    line.push(';');
                }
                first = false;
                line.push_str(key);
                if !value.is_empty() && def.count != FieldCount::Fixed(0) {
                    line.push('=');
                    line.push_str(value);
                }
            }
        }
        if !header.has_genotypes {
            return Ok(line);
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
        for key in &format {
            if !defs.format.contains_key(*key) {
                bail!("Key {key} found in VariantContext field FORMAT at {}:{} but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.", self.chromosome, self.position);
            }
        }
        line.push('\t');
        line.push_str(&format.join(":"));

        let has_gt = format[0] == "GT";
        let mut values: Vec<&str> = Vec::new();
        for g in &self.genotypes {
            line.push('\t');
            if has_gt {
                line.push_str(g.gt.as_deref().unwrap_or("./."));
            }
            values.clear();
            for key in format.iter().filter(|k| **k != "GT") {
                // a value the genotype does not have is written as a single ".", whatever the field's Number
                let value: &str = if *key == "FT" {
                    g.filters.as_deref().unwrap_or("PASS")
                } else if STANDARD_GENOTYPE_KEYS.contains(key) {
                    g.standard.get(*key).map(|s| s.as_str()).unwrap_or(MISSING)
                } else {
                    g.attributes.get(*key).map(|s| s.as_str()).unwrap_or(MISSING)
                };
                values.push(value);
            }
            // trailing missing values are dropped
            let mut n = values.len();
            while n > 0 && is_missing_value(values[n - 1]) {
                n -= 1;
            }
            for (i, value) in values[..n].iter().enumerate() {
                if i > 0 || has_gt {
                    line.push(':');
                }
                line.push_str(value);
            }
        }
        Ok(line)
    }
}

/// Header lines in the order htsjdk writes them: the file format first, then the rest sorted as
/// strings without duplicates, with the contig lines kept together in their own order.
pub fn sorted_header_lines(lines: Vec<String>) -> Vec<String> {
    let mut format: Option<String> = None;
    let mut contigs: Vec<String> = Vec::new();
    let mut others: Vec<String> = Vec::new();
    for line in lines {
        if line.starts_with("##fileformat=") {
            format.get_or_insert(line);
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
    let tail = others.split_off(split);
    let mut out: Vec<String> = format.into_iter().collect();
    out.extend(others);
    out.extend(contigs);
    out.extend(tail);
    out
}

/// The header and the data lines of a VCF (plain or gzip/bgzip compressed).
pub fn read_vcf(path: &str) -> Result<(VcfHeader, Vec<String>)> {
    let file = std::fs::File::open(path).with_context(|| format!("opening input VCF {path}"))?;
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") || path.ends_with(".bgz") {
        Box::new(std::io::BufReader::with_capacity(1 << 16, flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::with_capacity(1 << 16, file))
    };
    let mut header = VcfHeader { meta_lines: Vec::new(), samples: Vec::new(), has_genotypes: false, defs: HeaderDefs::default() };
    let mut records = Vec::new();
    let mut has_column_line = false;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("##") {
            header.add_meta_line(repair_standard_header_line(&line));
        } else if line.starts_with('#') {
            let cols: Vec<&str> = line.split('\t').collect();
            header.has_genotypes = cols.len() > 8;
            header.samples = cols.iter().skip(9).map(|s| s.to_string()).collect();
            has_column_line = true;
        } else if !line.is_empty() {
            records.push(line);
        }
    }
    if !has_column_line {
        bail!("Your input file has a malformed header: We never saw the required CHROM header line (starting with one #) for the input VCF file");
    }
    Ok((header, records))
}

/// A block-compressed VCF with a tabix index written on closing.
pub struct VcfWriter {
    path: String,
    out: bgzf::Writer,
}

impl VcfWriter {
    pub fn create(path: &str, header: &VcfHeader) -> Result<VcfWriter> {
        let mut out = bgzf::Writer::from_path(path).with_context(|| format!("creating {path}"))?;
        let mut lines = vec!["##fileformat=VCFv4.2".to_string()];
        lines.extend(header.meta_lines.iter().filter(|l| !l.starts_with("##fileformat=")).cloned());
        let mut text = sorted_header_lines(lines).join("\n");
        text.push_str("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if header.has_genotypes {
            text.push_str("\tFORMAT");
            for sample in &header.samples {
                text.push('\t');
                text.push_str(sample);
            }
        }
        text.push('\n');
        out.write_all(text.as_bytes())?;
        Ok(VcfWriter { path: path.to_string(), out })
    }

    pub fn write(&mut self, line: &str) -> Result<()> {
        self.out.write_all(line.as_bytes())?;
        self.out.write_all(b"\n")?;
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
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


#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> VcfHeader {
        let mut h = VcfHeader { meta_lines: Vec::new(), samples: vec!["S1".into(), "S2".into()], has_genotypes: true, defs: HeaderDefs::default() };
        for line in [
            "##INFO=<ID=AA,Number=1,Type=Integer,Description=\"a\">",
            "##INFO=<ID=ZZ,Number=1,Type=Float,Description=\"z\">",
            "##INFO=<ID=FLG,Number=0,Type=Flag,Description=\"f\">",
            "##INFO=<ID=LST,Number=2,Type=Integer,Description=\"l\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ad\">",
            "##FORMAT=<ID=VF,Number=1,Type=Integer,Description=\"vf\">",
            "##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"sb\">",
            "##FORMAT=<ID=PR,Number=2,Type=Integer,Description=\"pr\">",
        ] {
            h.add_meta_line(line.to_string());
        }
        h
    }

    #[test]
    fn round_trip_sorts_keys_and_sets_pass() {
        let h = header();
        let line = "1\t100\tid1\tA\tA[2:500[\t37.5\t.\tZZ=0.500;LST=1,2;FLG;AA=3\tGT:VF:AD:SB\t./.:5:0,5:0.250\t./.:7:0,7:0.5";
        let v = VariantContext::parse(line, &h).unwrap();
        assert_eq!(v.filters, None);
        assert_eq!(v.info_ints("LST").unwrap(), Some(vec![1, 2]));
        assert_eq!(v.genotypes[1].attribute_as_int("VF", 0).unwrap(), 7);
        let out = v.encode(&h).unwrap();
        assert_eq!(out, "1\t100\tid1\tA\tA[2:500[\t37.50\t.\tAA=3;FLG;LST=1,2;ZZ=0.500\tGT:AD:SB:VF\t./.:0,5:0.250:5\t./.:0,7:0.5:7");
    }

    #[test]
    fn missing_genotype_values() {
        let h = header();
        let line = "1\t100\t.\tn\t.A\t10\tb;a\t.\tGT:VF:PR\t./.:.:1,2\t./.";
        let mut v = VariantContext::parse(line, &h).unwrap();
        assert_eq!(v.ref_, "N");
        assert_eq!(v.alts, vec![".A".to_string()]);
        assert!(v.genotypes[1].attributes.is_empty());
        v.genotypes[0].attributes.insert("VF".into(), "4".into());
        // trailing missing values are dropped; a missing value in the middle is one "."
        assert_eq!(v.encode(&h).unwrap(), "1\t100\t.\tN\t.A\t10\ta;b\t.\tGT:PR:VF\t./.:1,2:4\t./.");
        v.genotypes[1].attributes.insert("VF".into(), "1".into());
        assert_eq!(v.encode(&h).unwrap(), "1\t100\t.\tN\t.A\t10\ta;b\t.\tGT:PR:VF\t./.:1,2:4\t./.:.:1");
    }

    #[test]
    fn undefined_key_is_an_error() {
        let h = header();
        let mut v = VariantContext::parse("1\t100\t.\tA\tT\t.\tPASS\tAA=1", &h).unwrap();
        v.info.insert("NEW".into(), "1".into());
        assert!(v.encode(&h).is_err());
    }

    #[test]
    fn header_order() {
        let lines = vec!["##b=2".to_string(), "##contig=<ID=2>".to_string(), "##fileformat=VCFv4.2".to_string(), "##contig=<ID=1>".to_string(), "##a=1".to_string(), "##zz=1".to_string()];
        assert_eq!(sorted_header_lines(lines), vec!["##fileformat=VCFv4.2", "##a=1", "##b=2", "##contig=<ID=2>", "##contig=<ID=1>", "##zz=1"]);
    }

    #[test]
    fn java_ints() {
        assert_eq!(parse_java_int("-12"), Some(-12));
        assert_eq!(parse_java_int("+7"), Some(7));
        assert_eq!(parse_java_int("1.0"), None);
        assert_eq!(parse_java_int(""), None);
    }
}
