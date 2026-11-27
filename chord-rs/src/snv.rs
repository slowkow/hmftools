//! SnvPrep and SnvDetails: SNVs counted by trinucleotide context (the 96 COSMIC-style bins, named as
//! the CHORD model expects them, e.g. A[C>A]A), with the reference read around each SNV.

use anyhow::{bail, Result};

use hmf_common::variant_context::{VariantType, VcfRecord};
use hmf_common::{debug, error, warn};

use crate::compat;
use crate::prep::MutContextCount;
use crate::prep::LazyRefGenome;
use crate::refgenome::RefGenomeSource;
use crate::vcf_file::{type_name, SmallVariant};

pub const SNV_DETAILS_FILE_SUFFIX: &str = ".chord.snv.details.tsv";

const VALID_NUCLEOTIDES: &[u8] = b"ACGTRYSWKMBDHVN";
const AMBIGUOUS_NUCLEOTIDES: &[u8] = b"RYSWKMBDHVN";

const MAX_AMBIGUOUS_TRI_NUCLEOTIDE_CONTEXT_FRACTION: f32 = 0.1;

/// SnvSigUtils.populateBucketMap: C>A_ACA, C>A_ACC, ... T>G_TTT in index order
pub fn bucket_names() -> Vec<String> {
    let ref_bases = ['C', 'T'];
    let bases = ['A', 'C', 'G', 'T'];
    let mut out = Vec::with_capacity(96);
    for r in ref_bases {
        for a in bases {
            if r == a {
                continue;
            }
            for before in bases {
                for after in bases {
                    out.push(format!("{r}>{a}_{before}{r}{after}"));
                }
            }
        }
    }
    out
}

/// Nucleotides.swapDnaBase
pub fn swap_dna_base(b: char) -> char {
    match b {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// SnvSigUtils.variantContext: the base change on the C or T strand, with its trinucleotide context
pub fn variant_context(ref_bases: &str, alt_bases: &str, tri_nuc: &str) -> String {
    let r = ref_bases.chars().next().unwrap_or('\0');
    let a = alt_bases.chars().next().unwrap_or('\0');
    let t: Vec<char> = tri_nuc.chars().collect();
    if r == 'A' || r == 'G' {
        format!("{}>{}_{}{}{}", swap_dna_base(r), swap_dna_base(a), swap_dna_base(t[2]), swap_dna_base(t[1]), swap_dna_base(t[0]))
    } else {
        format!("{r}>{a}_{tri_nuc}")
    }
}

/// SnvDetails.renameTriNucBin: "C>A_ACA" to "A[C>A]A"
pub fn rename_tri_nuc_bin(context: &str) -> String {
    let parts = hmf_common::java_text::split(context, '_');
    let sub = parts.first().copied().unwrap_or("");
    let tri: Vec<char> = parts.get(1).copied().unwrap_or("").chars().collect();
    format!("{}[{}]{}", tri.first().copied().unwrap_or('?'), sub, tri.get(2).copied().unwrap_or('?'))
}

#[derive(Clone, Debug)]
pub struct SnvDetails {
    pub chromosome: String,
    pub position: i32,
    pub ref_bases: String,
    pub alt_bases: String,
    pub tri_nuc_sequence: String,
    pub tri_nuc_context: String,
    pub tri_nuc_context_renamed: String,
}

impl SnvDetails {
    /// SnvDetails.from: Err where Java throws (the bases before the first or after the last of a contig)
    pub fn from(snv: &SmallVariant, ref_genome: &RefGenomeSource) -> Result<SnvDetails> {
        let tri_nuc_bytes = if compat::java("CH5") {
            match ref_genome.get_bases(&snv.chromosome, snv.position - 1, snv.position + 1)? {
                Some(b) => b,
                None => bail!("java.lang.NullPointerException: Cannot invoke \"java.lang.String.<init>(byte[])\" because bases are null"),
            }
        } else {
            ref_genome.get_bases_padded(&snv.chromosome, snv.position - 1, snv.position + 1)?
        };
        let tri_nuc_sequence = String::from_utf8_lossy(&tri_nuc_bytes).into_owned();
        let tri_nuc_context = variant_context(&snv.ref_bases, &snv.alt_bases, &tri_nuc_sequence);
        let renamed = rename_tri_nuc_bin(&tri_nuc_context);
        Ok(SnvDetails {
            chromosome: snv.chromosome.clone(),
            position: snv.position,
            ref_bases: snv.ref_bases.clone(),
            alt_bases: snv.alt_bases.clone(),
            tri_nuc_sequence,
            tri_nuc_context,
            tri_nuc_context_renamed: renamed,
        })
    }

    pub fn header() -> String {
        ["Chromosome", "Position", "RefBases", "AltBases", "TriNucSequence", "TriNucContext", "TriNucContextRenamed"].join("\t")
    }

    pub fn line(&self) -> String {
        [
            self.chromosome.as_str(),
            &self.position.to_string(),
            &self.ref_bases,
            &self.alt_bases,
            &self.tri_nuc_sequence,
            &self.tri_nuc_context,
            &self.tri_nuc_context_renamed,
        ]
        .join("\t")
    }
}

/// SnvPrep.isAmbiguousTriNucContext
fn is_ambiguous_tri_nuc_context(d: &SnvDetails) -> bool {
    let t = d.tri_nuc_sequence.as_bytes();
    let relevant = [d.alt_bases.as_bytes().first().copied().unwrap_or(0), t.first().copied().unwrap_or(0), t.get(1).copied().unwrap_or(0), t.get(2).copied().unwrap_or(0)];
    relevant.iter().all(|b| VALID_NUCLEOTIDES.contains(b)) && relevant.iter().any(|b| AMBIGUOUS_NUCLEOTIDES.contains(b))
}

/// SnvPrep.loadVariants: the SNP records (VariantContext.isSNP), with the skipped types counted
pub fn load_snvs(variants: &[VcfRecord], log_prefix: &str) -> Vec<SmallVariant> {
    let mut skipped: Vec<(VariantType, i32)> = Vec::new();
    let mut snvs = Vec::new();
    for v in variants {
        let t = v.variant_type();
        if t != VariantType::Snp {
            match skipped.iter_mut().find(|(k, _)| *k == t) {
                Some(e) => e.1 += 1,
                None => skipped.push((t, 1)),
            }
            continue;
        }
        snvs.push(SmallVariant::from(v));
    }
    debug!("{log_prefix}Found {} SNVs", snvs.len());
    if !skipped.is_empty() {
        // a HashMap keyed by an enum in Java (identity-hash order); the log line is not compared
        let text: Vec<String> = skipped.iter().map(|(k, n)| format!("{}={n}", type_name(*k))).collect();
        debug!("{log_prefix}Skipped variants: {{{}}}", text.join(", "));
    }
    snvs
}

/// SnvPrep.countMutationContexts after loading: Err where Java logs an error and returns null. With
/// `keep_details`, `write_details` gets the details of every SNV before the ambiguity check, where Java
/// writes the details file.
pub fn count_mutation_contexts(
    snvs: &[SmallVariant],
    ref_genome: &LazyRefGenome,
    keep_details: bool,
    log_prefix: &str,
    mut write_details: impl FnMut(&[SnvDetails]) -> Result<()>,
) -> Result<Vec<MutContextCount>> {
    let names = bucket_names();
    let mut counts = vec![0i32; names.len()];
    let index: std::collections::HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut details = Vec::new();
    let mut ambiguous_count = 0i32;
    for snv in snvs {
        let d = SnvDetails::from(snv, ref_genome.get()?)?;
        if let Some(&i) = index.get(d.tri_nuc_context.as_str()) {
            counts[i] += 1;
        } else if is_ambiguous_tri_nuc_context(&d) {
            warn!(
                "{log_prefix}Found ambiguous nucleotide in trinucleotide context '{}' for variant '{}:{}:{}:{}'",
                d.tri_nuc_context, d.chromosome, d.position, d.ref_bases, snv.alt_bases
            );
            ambiguous_count += 1;
        } else {
            error!(
                "{log_prefix}Found invalid trinucleotide context '{}' for variant '{}:{}:{}:{}'",
                d.tri_nuc_context, d.chromosome, d.position, d.ref_bases, snv.alt_bases
            );
            bail!("invalid trinucleotide context");
        }
        if keep_details {
            details.push(d);
        }
    }
    if keep_details {
        write_details(&details)?;
    }
    // float arithmetic as in Java: int > float * int
    if (ambiguous_count as f32) > MAX_AMBIGUOUS_TRI_NUCLEOTIDE_CONTEXT_FRACTION * (snvs.len() as i32 as f32) {
        error!(
            "{log_prefix}Fraction of trinucleotide contexts with ambiguous bases '{}' exceeds maximum '{}'",
            hmf_common::fmt::java_fixed((ambiguous_count as f32 / snvs.len() as f32) as f64, 2),
            hmf_common::fmt::java_fixed(MAX_AMBIGUOUS_TRI_NUCLEOTIDE_CONTEXT_FRACTION as f64, 2)
        );
        bail!("ambiguous trinucleotide contexts");
    }
    let out = names.iter().zip(counts).map(|(n, c)| MutContextCount { name: rename_tri_nuc_bin(n), count: c }).collect();
    debug!("{log_prefix}Completed SnvPrep");
    Ok(out)
}
