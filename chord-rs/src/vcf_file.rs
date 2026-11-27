//! VcfFile.loadVariants: the records of a VCF in file order, only those with FILTER PASS unless
//! -include_non_pass (a record with FILTER "." has no filters applied and is not taken as passing).

use anyhow::{bail, Result};

use hmf_common::debug;
use hmf_common::variant_context::{is_symbolic_allele, VcfHeader, VcfRecord};

/// The records of `path`; Err where Java throws (NoSuchFileException, a record htsjdk cannot parse).
///
/// The header is read here rather than by hmf-common's VcfReader: records are parsed against the INFO,
/// FORMAT and FILTER definitions only, and VcfHeader::from_lines checks each `##contig` line against all
/// earlier ones, which takes 0.3 s for the 3,366 contigs of PURPLE's SV VCF header on GRCh38.
pub fn load_variants(path: &str, include_non_pass: bool, log_prefix: &str) -> Result<Vec<VcfRecord>> {
    if !std::path::Path::new(path).is_file() {
        bail!("java.nio.file.NoSuchFileException: {path}");
    }
    let mut reader = hmf_common::gz::LineReader::open(path).map_err(|e| anyhow::anyhow!("reading {path}: {e}"))?;
    let mut meta: Vec<String> = Vec::new();
    let mut header: Option<VcfHeader> = None;
    let mut variants = Vec::new();
    while let Some(raw) = reader.next_line().map_err(|e| anyhow::anyhow!("reading {path}: {e}"))? {
        let line = String::from_utf8_lossy(raw);
        let Some(h) = &header else {
            if line.starts_with("##") {
                if line.starts_with("##INFO=") || line.starts_with("##FORMAT=") || line.starts_with("##FILTER=") {
                    meta.push(line.into_owned());
                }
                continue;
            }
            if line.starts_with('#') {
                let samples = line.split('\t').skip(9).map(str::to_string).collect();
                header = Some(VcfHeader::from_lines(std::mem::take(&mut meta), samples));
                continue;
            }
            bail!("htsjdk.tribble.TribbleException$MalformedFeatureFile: Unable to parse header with error: We never saw the required CHROM header line (starting with one #) for the input VCF file, for input source: {path}");
        };
        if line.is_empty() {
            continue;
        }
        let rec = VcfRecord::parse(&line, h)?;
        check_alleles(&rec)?;
        let is_explicit_pass = rec.filters.as_ref().is_some_and(|f| f.is_empty());
        if !include_non_pass && !is_explicit_pass {
            continue;
        }
        variants.push(rec);
    }
    if header.is_none() {
        bail!("htsjdk.tribble.TribbleException$MalformedFeatureFile: Unable to parse header with error: We never saw the required CHROM header line (starting with one #) for the input VCF file, for input source: {path}");
    }
    let filter_type = if include_non_pass { "" } else { "PASS " };
    debug!("{log_prefix}Loaded {} {filter_type}variants from: {path}", variants.len());
    Ok(variants)
}

/// htsjdk's Allele checks as the VCF codec applies them: bases A, C, G, T or N (either case), a
/// symbolic allele, or '*' for an ALT
fn check_alleles(rec: &VcfRecord) -> Result<()> {
    let valid_bases = |a: &str| !a.is_empty() && a.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'));
    if !valid_bases(&rec.ref_allele) {
        bail!(
            "htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number {}:{}: unparsable vcf record with allele {}",
            rec.chrom,
            rec.pos,
            rec.ref_allele
        );
    }
    for alt in &rec.alts {
        if !(is_symbolic_allele(alt) || valid_bases(alt)) {
            bail!(
                "htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line number {}:{}: unparsable vcf record with allele {alt}",
                rec.chrom,
                rec.pos
            );
        }
    }
    Ok(())
}

/// SmallVariant: what the SNV and INDEL preps take from a record
#[derive(Clone, Debug, PartialEq)]
pub struct SmallVariant {
    pub id: String,
    pub chromosome: String,
    pub position: i32,
    pub ref_bases: String,
    pub alt_bases: String,
}

impl SmallVariant {
    pub fn from(rec: &VcfRecord) -> SmallVariant {
        SmallVariant {
            id: rec.id.clone(),
            chromosome: rec.chrom.clone(),
            position: rec.pos,
            ref_bases: rec.ref_allele.clone(),
            alt_bases: rec.alt().to_string(),
        }
    }

    /// the @VisibleForTesting constructor
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn new(chromosome: &str, position: i32, ref_bases: &str, alt_bases: &str) -> SmallVariant {
        SmallVariant {
            id: String::new(),
            chromosome: chromosome.to_string(),
            position,
            ref_bases: ref_bases.to_string(),
            alt_bases: alt_bases.to_string(),
        }
    }
}

/// VariantContext.Type names as Java prints them
pub fn type_name(t: hmf_common::variant_context::VariantType) -> &'static str {
    use hmf_common::variant_context::VariantType::*;
    match t {
        NoVariation => "NO_VARIATION",
        Snp => "SNP",
        Mnp => "MNP",
        Indel => "INDEL",
        Symbolic => "SYMBOLIC",
        Mixed => "MIXED",
    }
}
