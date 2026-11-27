//! DataItem and DataItem.Index (with its ordering), DataItemsIO's writers, and the gzip stream as
//! java.util.zip.GZIPOutputStream writes it (header with OS 255, zlib's default level), so that the
//! `.tsv.gz` files equal Java's byte for byte.

use std::cmp::Ordering;
use std::io::Write;

use anyhow::{Context, Result};
use hmf_common::info;

use super::ItemMatrix;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Source {
    Dna,
    Rna,
}

impl Source {
    pub fn name(self) -> &'static str {
        match self {
            Source::Dna => "DNA",
            Source::Rna => "RNA",
        }
    }
}

/// ItemType, in declaration order (which Index.compareTo sorts by first)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ItemType {
    Signature,
    Snv96,
    GenPos,
    TumorMutationalBurden,
    SampleTrait,
    SvCount,
    Virus,
    Fusion,
    Driver,
    Expression,
    AltSj,
}

impl ItemType {
    pub fn alias(self) -> &'static str {
        match self {
            ItemType::Signature => "sig",
            ItemType::Snv96 => "snv96",
            ItemType::GenPos => "gen_pos",
            ItemType::TumorMutationalBurden => "event.tmb",
            ItemType::SampleTrait => "event.trait",
            ItemType::SvCount => "event.sv",
            ItemType::Virus => "event.virus",
            ItemType::Fusion => "event.fusion",
            ItemType::Driver => "event.driver",
            ItemType::Expression => "gene_exp",
            ItemType::AltSj => "alt_sj",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Index {
    pub source: Source,
    pub item_type: ItemType,
    pub key: String,
}

/// `String.format("%2s", chrom).replace(' ', pad)`: digits padded with '0', others with 'a'
fn pad_chrom(chrom: &str) -> String {
    let pad = if chrom.chars().next().is_some_and(|c| c.is_ascii_digit()) { '0' } else { 'a' };
    let n = chrom.encode_utf16().count();
    if n >= 2 { chrom.replace(' ', &pad.to_string()) } else { format!("{pad}{chrom}").replace(' ', &pad.to_string()) }
}

fn java_parse_int(s: &str) -> i32 {
    // Integer.parseInt; Java throws for a malformed key, which the matrix keys never are
    hmf_common::java_text::parse_int(s).unwrap_or(0)
}

/// Index.compareGenomicLoci: chromosomes padded and compared as strings, then start (and end)
fn compare_loci(a: &str, b: &str, delim: char) -> Ordering {
    let s1: Vec<&str> = hmf_common::java_text::split_all(a, delim);
    let s2: Vec<&str> = hmf_common::java_text::split_all(b, delim);
    let c = hmf_common::java_text::compare(&pad_chrom(s1[0]), &pad_chrom(s2[0]));
    if c != Ordering::Equal {
        return c;
    }
    let d = java_parse_int(s1[1]).wrapping_sub(java_parse_int(s2[1]));
    if d != 0 {
        return d.cmp(&0);
    }
    if s1.len() == 3 {
        let e = java_parse_int(s1[2]).wrapping_sub(java_parse_int(s2[2]));
        if e != 0 {
            return e.cmp(&0);
        }
    }
    Ordering::Equal
}

impl Index {
    /// DataItem.Index.compareTo
    pub fn compare(&self, other: &Index) -> Ordering {
        let t = self.item_type.cmp(&other.item_type);
        if t != Ordering::Equal {
            return t;
        }
        match self.item_type {
            ItemType::GenPos => compare_loci(&self.key, &other.key, '_'),
            ItemType::AltSj => compare_loci(&self.key, &other.key, ';'),
            ItemType::Signature => {
                let n1 = java_parse_int(self.key.split('_').nth(1).unwrap_or(""));
                let n2 = java_parse_int(other.key.split('_').nth(1).unwrap_or(""));
                n1.wrapping_sub(n2).cmp(&0)
            }
            _ => hmf_common::java_text::compare(&self.key, &other.key),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DataItem {
    pub index: Index,
    pub value: String,
}

impl DataItem {
    pub fn int(source: Source, item_type: ItemType, key: &str, v: i64) -> DataItem {
        DataItem { index: Index { source, item_type, key: key.to_string() }, value: v.to_string() }
    }
    pub fn bool(source: Source, item_type: ItemType, key: &str, v: bool) -> DataItem {
        DataItem { index: Index { source, item_type, key: key.to_string() }, value: if v { "1" } else { "0" }.to_string() }
    }
    /// String.format(Locale.ENGLISH, "%.<places>f", v)
    pub fn fixed(source: Source, item_type: ItemType, key: &str, v: f64, places: usize) -> DataItem {
        DataItem { index: Index { source, item_type, key: key.to_string() }, value: hmf_common::fmt::java_fixed(v, places) }
    }
    /// String.format(Locale.ENGLISH, "%6.3e", v)
    pub fn sci(source: Source, item_type: ItemType, key: &str, v: f64) -> DataItem {
        let s = hmf_common::fmt::java_sci(v, 3);
        let s = if s.len() < 6 { format!("{s:>6}") } else { s };
        DataItem { index: Index { source, item_type, key: key.to_string() }, value: s }
    }
}

/// A file written as FileWriterUtils.createBufferedWriter(path, append) writes it: plain, or for a ".gz" name
/// one gzip member (appended to what the file holds when `append`)
pub struct JavaWriter {
    path: String,
    gz: bool,
    buf: Vec<u8>,
    append: bool,
}

impl JavaWriter {
    pub fn create(path: &str, append: bool) -> Result<JavaWriter> {
        // FileOutputStream(path, append) opens (and creates) the file at once
        std::fs::OpenOptions::new().create(true).append(append).write(true).truncate(!append).open(path).with_context(|| {
            format!("java.io.FileNotFoundException: {path} (No such file or directory)")
        })?;
        Ok(JavaWriter { path: path.to_string(), gz: path.ends_with(".gz"), buf: Vec::new(), append })
    }

    pub fn line(&mut self, s: &str) {
        self.buf.extend_from_slice(s.as_bytes());
        self.buf.push(b'\n');
    }

    pub fn close(self) -> Result<()> {
        let mut f = std::fs::OpenOptions::new().append(self.append).write(true).open(&self.path)?;
        if self.gz {
            f.write_all(&java_gzip(&self.buf)?)?;
        } else {
            f.write_all(&self.buf)?;
        }
        Ok(())
    }
}

/// One gzip member as GZIPOutputStream writes it: magic, deflate, no flags, mtime 0, XFL 0, OS 255; the
/// deflate stream of zlib at level 6 (Deflater.DEFAULT_COMPRESSION); CRC-32 and size
pub fn java_gzip(data: &[u8]) -> Result<Vec<u8>> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::new(6));
    enc.write_all(data)?;
    out.extend_from_slice(&enc.finish()?);
    let mut crc = flate2::Crc::new();
    crc.update(data);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    Ok(out)
}

/// DataItemsIO.writeDataItemList
pub fn write_item_list(items: &[DataItem], path: &str) -> Result<()> {
    info!("Writing data to: {path}");
    let mut w = JavaWriter::create(path, false).inspect_err(|e| hmf_common::error!("Failed to write features: {e:#}"))?;
    w.line("Source\tCategory\tKey\tValue");
    for d in items {
        w.line(&format!("{}\t{}\t{}\t{}", d.index.source.name(), d.index.item_type.alias(), d.index.key, d.value));
    }
    w.close()
}

/// DataItemsIO.writeDataItemMatrix (nulls written as "null", as StringJoiner adds them)
pub fn write_matrix(m: &ItemMatrix, path: &str, append: bool) -> Result<()> {
    info!("Writing data to: {path}");
    let mut w = JavaWriter::create(path, append).inspect_err(|e| hmf_common::error!("Failed to write multi-sample feature matrix: {e:#}"))?;
    if !append {
        let mut h = vec!["Source".to_string(), "Category".to_string(), "Key".to_string()];
        h.extend(m.sample_ids.iter().cloned());
        w.line(&h.join("\t"));
    }
    for idx in &m.indexes {
        let mut f = vec![idx.source.name().to_string(), idx.item_type.alias().to_string(), idx.key.clone()];
        for v in &m.values[idx] {
            f.push(v.clone().unwrap_or_else(|| "null".to_string()));
        }
        w.line(&f.join("\t"));
    }
    w.close()
}

#[cfg_attr(not(test), allow(dead_code))] // DataItemsIO's readers: used by the ported Java tests
impl ItemType {
    /// ItemType.fromAlias
    pub fn from_alias(alias: &str) -> Result<ItemType> {
        const ALL: [ItemType; 11] = [
            ItemType::Signature,
            ItemType::Snv96,
            ItemType::GenPos,
            ItemType::TumorMutationalBurden,
            ItemType::SampleTrait,
            ItemType::SvCount,
            ItemType::Virus,
            ItemType::Fusion,
            ItemType::Driver,
            ItemType::Expression,
            ItemType::AltSj,
        ];
        ALL.into_iter().find(|t| t.alias() == alias).ok_or_else(|| anyhow::anyhow!("java.lang.IllegalArgumentException: Invalid alias: {alias}"))
    }
}

#[cfg_attr(not(test), allow(dead_code))]
fn source_value_of(s: &str) -> Result<Source> {
    match s {
        "DNA" => Ok(Source::Dna),
        "RNA" => Ok(Source::Rna),
        _ => anyhow::bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.cup.prep.DataSource.{s}"),
    }
}

/// DataItemsIO.readDataItemList
#[cfg_attr(not(test), allow(dead_code))]
pub fn read_item_list(path: &str) -> Result<Vec<DataItem>> {
    let lines = hmf_common::java_text::read_lines(path)?;
    let mut out = Vec::new();
    for line in lines.iter().skip(1) {
        let v: Vec<&str> = hmf_common::java_text::split_all(line, '\t');
        out.push(DataItem {
            index: Index { source: source_value_of(v[0])?, item_type: ItemType::from_alias(v[1])?, key: v[2].to_string() },
            value: v[3].to_string(),
        });
    }
    Ok(out)
}

/// DataItemsIO.readDataItemMatrix ("null" read as null)
#[cfg_attr(not(test), allow(dead_code))]
pub fn read_matrix(path: &str) -> Result<ItemMatrix> {
    let lines = hmf_common::java_text::read_lines(path)?;
    let header: Vec<&str> = hmf_common::java_text::split(&lines[0], '\t');
    let sample_ids: Vec<String> = header[3..].iter().map(|s| s.to_string()).collect();
    let mut indexes = Vec::new();
    let mut values = std::collections::HashMap::new();
    for line in lines.iter().skip(1) {
        let v: Vec<&str> = hmf_common::java_text::split_all(line, '\t');
        let idx = Index { source: source_value_of(v[0])?, item_type: ItemType::from_alias(v[1])?, key: v[2].to_string() };
        let row: Vec<Option<String>> = v[3..header.len()].iter().map(|s| if *s == "null" { None } else { Some(s.to_string()) }).collect();
        indexes.push(idx.clone());
        values.insert(idx, row);
    }
    Ok(ItemMatrix { sample_ids, indexes, values })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn idx(t: ItemType, k: &str) -> Index {
        Index { source: Source::Dna, item_type: t, key: k.to_string() }
    }

    #[test]
    fn index_order() {
        // v38 chromosome names sort as strings, v37 names numerically (padded with '0'), X after digits
        assert_eq!(idx(ItemType::GenPos, "chr10_0").compare(&idx(ItemType::GenPos, "chr2_0")), Ordering::Less);
        assert_eq!(idx(ItemType::GenPos, "2_0").compare(&idx(ItemType::GenPos, "10_0")), Ordering::Less);
        assert_eq!(idx(ItemType::GenPos, "X_0").compare(&idx(ItemType::GenPos, "22_500000")), Ordering::Greater);
        assert_eq!(idx(ItemType::GenPos, "1_500000").compare(&idx(ItemType::GenPos, "1_0")), Ordering::Greater);
        assert_eq!(idx(ItemType::Signature, "SIG_10_POLE").compare(&idx(ItemType::Signature, "SIG_2_13_AID_APOBEC")), Ordering::Greater);
        assert_eq!(idx(ItemType::Snv96, "C>A_ACA").compare(&idx(ItemType::Signature, "SIG_1")), Ordering::Greater);
        assert_eq!(idx(ItemType::AltSj, "chr1;100;300").compare(&idx(ItemType::AltSj, "chr1;100;200")), Ordering::Greater);
    }

    #[test]
    fn formats() {
        assert_eq!(DataItem::sci(Source::Rna, ItemType::Expression, "G", 1.2345).value, "1.235e+00");
        assert_eq!(DataItem::fixed(Source::Dna, ItemType::Signature, "S", 2521.66, 1).value, "2521.7");
    }
}
