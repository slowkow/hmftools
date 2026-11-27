//! SigSnvLoader with the parts of SmallVariantFactory and VariantFilters it uses: the 96 trinucleotide
//! bucket counts of a sample's passing SNVs in a VCF.

use std::collections::HashMap;

use anyhow::{anyhow, bail, Context, Result};
use hmf_common::variant_context::{VariantType, VcfReader, VcfRecord};

use crate::compat;
use crate::config::{VariantFilters, TUMOR_ONLY_AF_THRESHOLD};
use crate::position_freq::PositionFrequencies;

pub const SNV_BUCKET_COUNT: usize = 96;

/// SnvSigUtils.populateBucketMap: "C>A_ACA" ... "T>G_TTT", in index order
pub fn bucket_names() -> Vec<String> {
    let ref_bases = ['C', 'T'];
    let bases = ['A', 'C', 'G', 'T'];
    let mut names = Vec::with_capacity(SNV_BUCKET_COUNT);
    for r in ref_bases {
        for alt in bases {
            if r == alt {
                continue;
            }
            for before in bases {
                for after in bases {
                    names.push(format!("{r}>{alt}_{before}{r}{after}"));
                }
            }
        }
    }
    names
}

/// Nucleotides.swapDnaBase
fn swap_dna_base(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// SnvSigUtils.variantContext: the base change on the C or T strand, with its context. Err where Java
/// throws StringIndexOutOfBoundsException (a context of fewer than 3 bases with an A or G reference).
pub fn variant_context(ref_: &str, alt: &str, trinuc: &str) -> Result<String> {
    let r = ref_.chars().next().ok_or_else(|| anyhow!("java.lang.StringIndexOutOfBoundsException: empty ref"))?;
    let a = alt.chars().next().ok_or_else(|| anyhow!("java.lang.StringIndexOutOfBoundsException: empty alt"))?;
    if r == 'A' || r == 'G' {
        let t: Vec<char> = trinuc.chars().collect();
        if t.len() < 3 {
            bail!("java.lang.StringIndexOutOfBoundsException: trinucleotide context '{trinuc}' of a {ref_}>{alt} SNV");
        }
        Ok(format!("{}>{}_{}{}{}", swap_dna_base(r), swap_dna_base(a), swap_dna_base(t[2]), swap_dna_base(t[1]), swap_dna_base(t[0])))
    } else {
        Ok(format!("{r}>{a}_{trinuc}"))
    }
}

/// SomaticLikelihood
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SomaticLikelihood {
    Unknown,
    High,
    Medium,
    Low,
}

impl SomaticLikelihood {
    fn value_of(s: &str) -> Result<SomaticLikelihood> {
        Ok(match s {
            "UNKNOWN" => SomaticLikelihood::Unknown,
            "HIGH" => SomaticLikelihood::High,
            "MEDIUM" => SomaticLikelihood::Medium,
            "LOW" => SomaticLikelihood::Low,
            _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.variant.SomaticLikelihood.{s}"),
        })
    }
}

/// The fields of SmallVariant (and its VariantContext) that SIGS reads
#[derive(Clone, Debug)]
pub struct SigVariant {
    pub chromosome: String,
    pub position: i32,
    pub ref_: String,
    /// the alternate alleles joined by ','
    pub alt: String,
    pub is_snp: bool,
    pub trinucleotide_context: String,
    pub subclonal_likelihood: f64,
    pub variant_copy_number: f64,
    pub somatic_likelihood: SomaticLikelihood,
    pub hotspot: bool,
    pub allele_read_count: i32,
    pub total_read_count: i32,
    /// getAttributeAsDouble(GND_FREQ, -1)
    pub gnomad_freq: f64,
}

impl SigVariant {
    pub fn allele_frequency(&self) -> f64 {
        if self.total_read_count > 0 { self.allele_read_count as f64 / self.total_read_count as f64 } else { 0.0 }
    }
}

/// HumanChromosome.contains: "chr" stripped, then 1 to 22, X or Y
pub fn human_chromosome_contains(contig: &str) -> bool {
    let trimmed = contig.strip_prefix("chr").unwrap_or(contig);
    if trimmed.bytes().all(|b| b.is_ascii_digit()) {
        return trimmed.parse::<i64>().is_ok_and(|v| (1..=22).contains(&v));
    }
    trimmed == "X" || trimmed == "Y"
}

/// SmallVariantFactory.createVariant(sample, context) for a record that is not filtered: None where Java
/// gives Optional.empty() (a non-human contig, an NT other than "ref", no genotype for the sample, no
/// allelic depths of the alt, or no reads).
pub fn create_variant(rec: &mut VcfRecord, header: &hmf_common::variant_context::VcfHeader, sample: &str) -> Result<Option<SigVariant>> {
    if !human_chromosome_contains(&rec.chrom) {
        return Ok(None);
    }
    if let Some(nt) = rec.info.get("NT") {
        let is_ref = matches!(nt, hmf_common::variant_context::Value::Raw(s) if s == "ref");
        if !is_ref {
            return Ok(None);
        }
    }
    rec.decode(header).with_context(|| format!("decoding the genotypes of {}:{}", rec.chrom, rec.pos))?;
    let Some(g) = rec.decoded_genotype(sample) else { return Ok(None) };
    let Some(ad) = &g.ad else { return Ok(None) };
    if ad.len() <= 1 {
        return Ok(None);
    }
    let allele_read_count = ad[1];
    let ad_sum: i32 = ad.iter().fold(0i32, |s, v| s.wrapping_add(*v));
    let total_read_count = match g.dp {
        Some(dp) => dp.max(ad_sum),
        None => ad_sum,
    };
    if total_read_count <= 0 {
        return Ok(None);
    }
    let somatic_likelihood = SomaticLikelihood::value_of(&rec.info.as_string("SOM_LH", "UNKNOWN"))?;
    Ok(Some(SigVariant {
        chromosome: rec.chrom.clone(),
        position: rec.pos,
        ref_: rec.ref_allele.clone(),
        alt: rec.alts.join(","),
        is_snp: rec.variant_type() == VariantType::Snp,
        trinucleotide_context: rec.info.as_string("TNC", ""),
        subclonal_likelihood: rec.info.as_double("SUBCL", 0.0),
        variant_copy_number: rec.info.as_double("PURPLE_VCN", 0.0),
        somatic_likelihood,
        hotspot: rec.info.as_bool("HOTSPOT", false),
        allele_read_count,
        total_read_count,
        gnomad_freq: rec.info.as_double("GND_FREQ", -1.0),
    }))
}

impl VariantFilters {
    /// VariantFilters.passesFilters
    pub fn passes(&self, v: &SigVariant) -> bool {
        if self.subclonal_min.is_some_and(|m| v.subclonal_likelihood < m) {
            return false;
        }
        if self.subclonal_max.is_some_and(|m| v.subclonal_likelihood > m) {
            return false;
        }
        if self.ploidy_min.is_some_and(|m| v.variant_copy_number < m) {
            return false;
        }
        if self.ploidy_max.is_some_and(|m| v.variant_copy_number > m) {
            return false;
        }
        if self.tumor_only {
            if v.somatic_likelihood == SomaticLikelihood::Low || v.hotspot {
                return false;
            }
            if v.gnomad_freq > 0.0 {
                return false;
            }
            if v.allele_frequency() < TUMOR_ONLY_AF_THRESHOLD {
                return false;
            }
        }
        true
    }
}

/// loadSomaticVariants: the SNVs of the sample among the VCF's unfiltered records. The VCF is read whole,
/// with or without an index, as htsjdk's reader does without requiring one.
pub fn load_somatic_variants(vcf_file: &str, sample: &str) -> Result<Vec<SigVariant>> {
    let mut reader = VcfReader::open(vcf_file)?;
    if !reader.header.samples.iter().any(|s| s == sample) && !compat::java("SI3") {
        bail!("sample {sample} is not in the VCF {vcf_file}");
    }
    let header = reader.header.clone();
    let mut out = Vec::new();
    while let Some(line) = reader.next_line()? {
        // filtered records and other contigs than 1-22, X and Y are skipped before the record is parsed
        let mut cols = line.splitn(8, '\t');
        let chrom = cols.next().unwrap_or("");
        let filter = cols.nth(5).unwrap_or("");
        if (filter != "PASS" && filter != ".") || !human_chromosome_contains(chrom) {
            continue;
        }
        let mut rec = VcfRecord::parse(line, &header)?;
        if rec.is_filtered() {
            continue;
        }
        let Some(v) = create_variant(&mut rec, &header, sample)? else { continue };
        if !v.is_snp {
            continue;
        }
        out.push(v);
    }
    Ok(out)
}

/// The bucket counts of one sample
pub struct SampleCounts {
    pub counts: Vec<f64>,
}

/// processSampleVariants: the counts per bucket, the positions added to `position_frequencies`. Err where
/// Java exits (a base change and context that make no bucket name, or a context too short to read).
pub fn count_variants(
    sample_id: &str,
    variants: &[SigVariant],
    filters: &VariantFilters,
    bucket_index: &HashMap<String, usize>,
    mut position_frequencies: Option<&mut PositionFrequencies>,
) -> Result<SampleCounts> {
    let mut counts = vec![0.0; SNV_BUCKET_COUNT];
    let mut filtered = 0;
    for v in variants {
        if !v.is_snp || v.alt.chars().count() != 1 {
            continue;
        }
        if v.trinucleotide_context.contains('N') {
            continue;
        }
        if !filters.passes(v) {
            filtered += 1;
            continue;
        }
        if let Some(pf) = position_frequencies.as_deref_mut() {
            pf.add_position(&v.chromosome, v.position);
        }
        let bucket = variant_context(&v.ref_, &v.alt, &v.trinucleotide_context)?;
        let Some(&i) = bucket_index.get(&bucket) else {
            hmf_common::error!(
                "sample({}) invalid bucketName({}) from var({}>{}) context={})",
                sample_id,
                bucket,
                v.ref_,
                v.alt,
                v.trinucleotide_context
            );
            bail!("invalid bucket name {bucket}");
        };
        counts[i] += 1.0;
    }
    if filtered > 0 {
        hmf_common::debug!("filtered {} variants", filtered);
    }
    Ok(SampleCounts { counts })
}

/// writeSampleCounts: "BucketName,<sample>", then each bucket's count as "%.0f"; an error writing is logged
/// and the run goes on, as in Java
pub fn write_sample_counts(filename: &str, sample_ids: &[String], names: &[String], counts: &[Vec<f64>]) {
    let mut out = String::with_capacity(2048);
    out.push_str("BucketName");
    for s in sample_ids {
        out.push(',');
        out.push_str(s);
    }
    out.push('\n');
    for (i, name) in names.iter().enumerate() {
        out.push_str(name);
        for c in counts {
            out.push(',');
            out.push_str(&hmf_common::fmt::java_fixed(c[i], 0));
        }
        out.push('\n');
    }
    if let Err(e) = std::fs::write(filename, out) {
        hmf_common::error!("error writing to outputFile: java.io.FileNotFoundException: {filename} ({e})");
    }
}
