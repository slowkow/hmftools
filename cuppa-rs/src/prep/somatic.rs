//! SomaticVariantPrep (category SNV) with SomaticVariantsLoader, SomaticVariant, TrinucleotideCounts,
//! SomaticSigs and GenomicPositions.

use std::io::BufRead;

use anyhow::{anyhow, bail, Context, Result};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::sigs::{self, LeastSquaresFit, PositionCounts};
use hmf_common::variant_context::{VariantType as HtsType, VcfHeader, VcfRecord};
use hmf_common::{error, warn};

use super::config::{PrepConfig, GEN_POS_BUCKET_SIZE};
use super::items::{DataItem, ItemType, Source};
use crate::compat;

/// hmf-common VariantType as CUPPA uses it
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantType {
    Snp,
    Mnp,
    Indel,
    Sv,
    Undefined,
}

impl VariantType {
    fn value_of(s: &str) -> Result<VariantType> {
        Ok(match s {
            "SNP" => VariantType::Snp,
            "MNP" => VariantType::Mnp,
            "INDEL" => VariantType::Indel,
            "SV" => VariantType::Sv,
            "UNDEFINED" => VariantType::Undefined,
            _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.variant.VariantType.{s}"),
        })
    }
}

#[derive(Clone, Debug)]
pub struct SomaticVariant {
    pub chromosome: String,
    pub position: i32,
    pub ref_: String,
    pub alt: String,
    pub variant_type: VariantType,
    pub gene: String,
    pub trinucleotide_context: String,
    pub repeat_count: i32,
}

/// htsjdk 2.23 Allele.wouldBeSymbolicAllele: unlike hmf-common-rs's is_symbolic_allele, the spanning deletion
/// '*' is not symbolic (it is a one-base allele, so a C>* record is an SNP)
fn symbolic(a: &str) -> bool {
    let b = a.as_bytes();
    b.len() > 1 && (b[0] == b'<' || b[b.len() - 1] == b'>' || a.contains('[') || a.contains(']') || b[0] == b'.' || b[b.len() - 1] == b'.')
}

/// htsjdk VariantContext.getType (determinePolymorphicType), as VariantType.type maps it
fn hts_type(rec: &VcfRecord) -> VariantType {
    if rec.alts.is_empty() {
        return VariantType::Undefined;
    }
    let mut t: Option<HtsType> = None;
    for alt in &rec.alts {
        let this = if symbolic(alt) {
            HtsType::Symbolic
        } else if rec.ref_allele.len() == alt.len() {
            if alt.len() == 1 { HtsType::Snp } else { HtsType::Mnp }
        } else {
            HtsType::Indel
        };
        match t {
            None => t = Some(this),
            Some(prev) if prev != this => return VariantType::Undefined,
            _ => {}
        }
    }
    match t {
        Some(HtsType::Snp) => VariantType::Snp,
        Some(HtsType::Mnp) => VariantType::Mnp,
        Some(HtsType::Indel) => VariantType::Indel,
        _ => VariantType::Undefined,
    }
}

const CODING_EFFECTS: [&str; 6] = ["NONSENSE_OR_FRAMESHIFT", "SPLICE", "MISSENSE", "SYNONYMOUS", "NONE", "UNDEFINED"];

/// A passing record as read, its IMPACT values kept for SomaticVariant.fromContext to parse when the variant
/// is used (Java parses them only for the variants of the requested types)
#[derive(Clone)]
struct RawVariant {
    variant: SomaticVariant,
    impact: Option<Vec<String>>,
}

fn raw_from_context(rec: &VcfRecord) -> RawVariant {
    let mut alt = rec.alts.first().cloned().unwrap_or_else(|| rec.ref_allele.clone());
    if alt == "*" || alt == "N" {
        alt = rec.ref_allele.clone();
    }
    RawVariant {
        variant: SomaticVariant {
            chromosome: rec.chrom.clone(),
            position: rec.pos,
            ref_: rec.ref_allele.clone(),
            alt,
            variant_type: hts_type(rec),
            gene: String::new(),
            trinucleotide_context: rec.info.as_string("TNC", ""),
            repeat_count: rec.info.as_int("REP_C", 0),
        },
        impact: if rec.info.has("IMPACT") { Some(rec.info.as_string_list("IMPACT", "")) } else { None },
    }
}

/// SomaticVariant.fromContext's gene (VariantImpactSerialiser.fromAttributeValues: the first of exactly 10
/// values, whose coding effects and gene count must parse)
fn from_raw(raw: &RawVariant) -> Result<SomaticVariant> {
    let mut v = raw.variant.clone();
    if let Some(values) = &raw.impact {
        if values.len() == 10 {
            for i in [3, 8] {
                if !CODING_EFFECTS.contains(&values[i].as_str()) {
                    bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.variant.CodingEffect.{}", values[i]);
                }
            }
            if hmf_common::java_text::parse_int(&values[9]).is_none() {
                bail!("java.lang.NumberFormatException: For input string: \"{}\"", values[9]);
            }
            v.gene = values[0].clone();
        }
    }
    Ok(v)
}

/// In single-sample mode the somatic VCF is read once for the SNV and DRIVER categories (Java reads it
/// twice); a cohort reads each sample's file per category, holding none
pub static CACHE_VCF: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
static VCF_CACHE: std::sync::Mutex<Option<(String, std::sync::Arc<Vec<RawVariant>>)>> = std::sync::Mutex::new(None);

/// The header of a VCF, built from its INFO, FORMAT and FILTER lines (the record parsing needs no more;
/// VcfHeader::from_lines is quadratic in PURPLE's 3,366 contig lines)
fn read_header(reader: &mut dyn BufRead, path: &str) -> Result<VcfHeader> {
    let mut meta = Vec::new();
    let samples: Vec<String>;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            bail!("invalid VCF({path}): no header line");
        }
        let t = line.trim_end_matches(['\n', '\r']);
        if let Some(m) = t.strip_prefix("##") {
            if !m.starts_with("contig=") {
                meta.push(t.to_string());
            }
        } else if t.starts_with('#') {
            samples = t.split('\t').skip(9).map(String::from).collect();
            break;
        } else {
            bail!("invalid VCF({path}): no header line");
        }
    }
    Ok(VcfHeader::from_lines(meta, samples))
}

/// a plain or gzip / BGZF compressed text file; gzip is inflated with hmf-common's zlib-rs reader, on a
/// thread of its own in single-sample mode (two threads in all; a cohort's sample threads are the two)
pub fn open_text(path: &str) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let mut r = std::io::BufReader::with_capacity(1 << 16, f);
    let gz = r.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    if !gz {
        return Ok(Box::new(r));
    }
    if (path.ends_with(".gz") || path.ends_with(".bgz")) && CACHE_VCF.load(std::sync::atomic::Ordering::Relaxed) {
        return Ok(Box::new(std::io::BufReader::with_capacity(1 << 17, hmf_common::gz::open_threaded(path)?)));
    }
    Ok(Box::new(std::io::BufReader::with_capacity(1 << 17, hmf_common::gz::GzReader::new(r))))
}

/// the passing records (FILTER PASS or '.') of a VCF
fn read_passing(path: &str) -> Result<Vec<RawVariant>> {
    let mut reader = open_text(path).map_err(|e| {
        error!("Invalid somatic VCF file({path})");
        e
    })?;
    let header = read_header(&mut *reader, path)?;
    let mut out = Vec::new();
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let t = line.trim_end_matches(['\n', '\r']);
        if t.is_empty() {
            continue;
        }
        // htsjdk (AbstractVCFCodec.decodeLine): 9 columns with genotypes in the header, else exactly 8
        let parts = t.splitn(9, '\t').count();
        let expected = if header.samples.is_empty() { 8 } else { 9 };
        if parts != expected {
            bail!("htsjdk.tribble.TribbleException: there aren't enough columns for line {t} (we expected {expected} tokens, and saw {parts} )");
        }
        // a cheap look at FILTER before parsing (the 7th column)
        let filter = t.split('\t').nth(6).unwrap_or("");
        if filter != "PASS" && filter != "." {
            continue;
        }
        let rec = VcfRecord::parse(t, &header)?;
        if rec.is_filtered() {
            continue;
        }
        out.push(raw_from_context(&rec));
    }
    Ok(out)
}

/// SomaticVariantsLoader.loadFromVcf: the passing records of the given types
fn load_from_vcf(path: &str, types: Option<&[VariantType]>) -> Result<Vec<SomaticVariant>> {
    let raws: std::sync::Arc<Vec<RawVariant>> = if CACHE_VCF.load(std::sync::atomic::Ordering::Relaxed) {
        let mut cache = VCF_CACHE.lock().unwrap();
        match cache.as_ref() {
            Some((p, v)) if p == path => v.clone(),
            _ => {
                let v = std::sync::Arc::new(read_passing(path)?);
                *cache = Some((path.to_string(), v.clone()));
                v
            }
        }
    } else {
        std::sync::Arc::new(read_passing(path)?)
    };
    raws.iter().filter(|r| types.is_none_or(|ts| ts.contains(&r.variant.variant_type))).map(from_raw).collect()
}

/// SomaticVariantsLoader.loadFromGenericFile (comma-separated, String.split dropping trailing empties)
fn load_from_generic(path: &str, types: Option<&[VariantType]>) -> Result<Vec<SomaticVariant>> {
    let lines = hmf_common::java_text::read_lines(path).map_err(|e| {
        error!("failed to load generic variants file({path}): {e}");
        anyhow!("{e}")
    })?;
    let Some(header) = lines.first() else { bail!("java.lang.IndexOutOfBoundsException: empty file {path}") };
    let f = hmf_common::java_text::fields_index_map(header, ',');
    let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {n} in {path}"));
    let (ci, pi, ri, ai, ti, rci, tni, gi) =
        (col("Chromosome")?, col("Position")?, col("Ref")?, col("Alt")?, col("Type")?, col("RepeatCount")?, col("TriNucContext")?, col("Gene")?);
    let mut out = Vec::new();
    for line in &lines[1..] {
        let v: Vec<&str> = hmf_common::java_text::split(line, ',');
        let get = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds: {line}"));
        let int = |s: &str| hmf_common::java_text::parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""));
        let variant = SomaticVariant {
            chromosome: get(ci)?.to_string(),
            position: int(get(pi)?)?,
            ref_: get(ri)?.to_string(),
            alt: get(ai)?.to_string(),
            variant_type: VariantType::value_of(get(ti)?)?,
            gene: get(gi)?.to_string(),
            trinucleotide_context: get(tni)?.to_string(),
            repeat_count: int(get(rci)?)?,
        };
        // Java: `if(variantTypes == null && !variantTypes.contains(...))`, which never filters and throws a
        // NullPointerException when no types are given (CU3); corrected, the types filter as for the VCF
        match types {
            None if compat::java("CU3") => bail!("java.lang.NullPointerException: Cannot invoke \"java.util.List.contains(Object)\" because \"variantTypes\" is null"),
            Some(ts) if !compat::java("CU3") && !ts.contains(&variant.variant_type) => continue,
            _ => {}
        }
        out.push(variant);
    }
    Ok(out)
}

/// SomaticVariantsLoader.loadFromConfig: the generic file if it exists, else PURPLE's somatic VCF
pub fn load_variants(config: &PrepConfig, sample: &str, types: Option<&[VariantType]>) -> Result<Vec<SomaticVariant>> {
    let generic = config.somatic_variants_generic(sample);
    let vcf = config.purple_somatic_vcf(sample);
    let is_file = |p: &str| std::fs::metadata(p).is_ok_and(|m| m.is_file());
    if is_file(&generic) {
        if is_file(&vcf) {
            warn!("Both VCF and generic variants files exist for sample({sample}). Using generic variants file: {generic}");
        }
        return load_from_generic(&generic, types);
    }
    if is_file(&vcf) {
        return load_from_vcf(&vcf, types);
    }
    bail!("java.nio.file.NoSuchFileException: Invalid VCF file ('{vcf}') or generic variants file ('{generic}')")
}

const AID_APOBEC_CONTEXTS: [&str; 8] = ["C>T_TCA", "C>T_TCC", "C>T_TCG", "C>T_TCT", "C>G_TCA", "C>G_TCC", "C>G_TCG", "C>G_TCT"];

/// SomaticSigs.REPORTABLE_SIGS, a HashMap (iterated in its hash order)
fn reportable_sigs() -> JavaHashMap<&'static str, &'static str> {
    let mut m = JavaHashMap::new();
    for (k, v) in [
        ("Sig1", "SIG_1"),
        ("Sig2", "SIG_2_13_AID_APOBEC"),
        ("Sig13", "SIG_2_13_AID_APOBEC"),
        ("Sig4", "SIG_4_SMOKING"),
        ("Sig6", "SIG_6_MMR"),
        ("Sig7", "SIG_7_UV"),
        ("Sig10", "SIG_10_POLE"),
        ("Sig11", "SIG_11"),
        ("Sig17", "SIG_17"),
    ] {
        m.insert(k, v);
    }
    m
}

/// SomaticVariantPrep.extractSampleData
pub fn extract(config: &PrepConfig, sample: &str) -> Result<Vec<DataItem>> {
    let variants = load_variants(config, sample, Some(&[VariantType::Snp]))?;
    let mut items = Vec::new();

    // the 96 trinucleotide counts (TrinucleotideCounts)
    let names = sigs::bucket_names();
    let mut counts = vec![0.0f64; names.len()];
    for v in &variants {
        if v.variant_type != VariantType::Snp || v.alt.len() != 1 || v.trinucleotide_context.contains('N') {
            continue;
        }
        let bucket = sigs::variant_context(&v.ref_, &v.alt, &v.trinucleotide_context)?;
        match names.iter().position(|n| *n == bucket) {
            Some(i) => counts[i] += 1.0,
            None => error!("invalid bucketName({bucket}) from var({}>{}) context={})", v.ref_, v.alt, v.trinucleotide_context),
        }
    }
    let mut total: i64 = 0;
    for (n, c) in names.iter().zip(&counts) {
        let count = *c as i32;
        total += count as i64;
        items.push(DataItem::int(Source::Dna, ItemType::Snv96, n, count as i64));
    }
    // Java sums into an int
    let total = total as i32;
    items.push(DataItem::int(Source::Dna, ItemType::TumorMutationalBurden, "snv_count", total as i64));

    // signature allocations (SomaticSigs.fitSampleCounts)
    let (matrix, sig_names) = sigs::load_snv_signatures()?;
    let sample_total: f64 = counts.iter().sum();
    let allocations: Vec<f64> = if sample_total == 0.0 {
        vec![0.0; sig_names.len()]
    } else {
        let mut fit = LeastSquaresFit::new(matrix.rows, matrix.cols)?;
        fit.initialise(&matrix.data, &counts);
        fit.solve()?;
        fit.contribs().to_vec()
    };
    let alloc_of = |name: &str| -> Result<f64> {
        sig_names.iter().position(|n| n == name).map(|i| allocations[i]).ok_or_else(|| anyhow!("java.lang.NullPointerException: no signature {name}"))
    };
    for (sig, reported) in reportable_sigs().iter() {
        if sig.eq_ignore_ascii_case("Sig13") {
            continue;
        }
        let mut a = alloc_of(sig)?;
        if sig.eq_ignore_ascii_case("Sig2") {
            a += alloc_of("Sig13")?;
        }
        items.push(DataItem::fixed(Source::Dna, ItemType::Signature, reported, a, 1));
    }

    // SNV counts in 500 kb bins, without AID/APOBEC contexts or chromosome Y (GenomicPositions)
    let mut pos = PositionCounts::standard(config.v38, GEN_POS_BUCKET_SIZE);
    for v in &variants {
        if v.variant_type != VariantType::Snp {
            continue;
        }
        if v.chromosome == "Y" || v.chromosome == "chrY" {
            continue;
        }
        if v.trinucleotide_context.contains('N') {
            continue;
        }
        let bucket = sigs::variant_context(&v.ref_, &v.alt, &v.trinucleotide_context)?;
        if AID_APOBEC_CONTEXTS.contains(&bucket.as_str()) {
            continue;
        }
        if !pos.is_valid_chromosome(&v.chromosome) {
            warn!("variant chr({}) position({}) cannot map to genomic position", v.chromosome, v.position);
            continue;
        }
        pos.add_position(&v.chromosome, v.position);
    }
    let chr_y = sigs::versioned_chromosome(config.v38, "chrY");
    for b in 0..pos.counts.len() {
        let chromosome = sigs::chromosome_from_index(config.v38, &pos.index, b as i32)?;
        if chromosome == chr_y {
            continue;
        }
        let position = sigs::position_from_index(&pos.index, &chromosome, b as i32, GEN_POS_BUCKET_SIZE)?;
        items.push(DataItem::int(Source::Dna, ItemType::GenPos, &format!("{chromosome}_{position}"), pos.counts[b] as i64));
    }
    Ok(items)
}
