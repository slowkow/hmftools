//! hmf-common `sigs`: the SNV signature code SIGS and CUPPA share.
//!
//! - `lsq`: LeastSquaresFit (Lawson and Hanson's non-negative least squares, step for step as the Java
//!   translation, so the contributions equal Java's to the last bit);
//! - `matrix`: utils.Matrix and MatrixFile.loadMatrixDataFile;
//! - here: SnvSigUtils (the 96 bucket names, variantContext, the COSMIC signatures resource
//!   `sigs/snv_cosmic_signatures.csv`) and PositionFrequencies (genomic position buckets).
//!
//! Moved from `sigs-rs` (lsq.rs, matrix.rs, parts of position_freq.rs and snv_loader.rs) so that CUPPA uses
//! the same code; behaviour unchanged.

pub mod lsq;
pub mod matrix;

use std::collections::BTreeMap;

use anyhow::{anyhow, bail, Result};

pub use lsq::LeastSquaresFit;
pub use matrix::{load_matrix_file, load_matrix_lines, Matrix};

pub const SNV_TRINUCLEOTIDE_BUCKET_COUNT: usize = 96;

/// hmf-common's resource sigs/snv_cosmic_signatures.csv (30 COSMIC v2 signatures over the 96 buckets)
pub const SNV_COSMIC_SIGNATURES: &str = include_str!("snv_cosmic_signatures.csv");

/// SnvSigUtils.populateBucketMap: "C>A_ACA" ... "T>G_TTT", in index order
pub fn bucket_names() -> Vec<String> {
    let ref_bases = ['C', 'T'];
    let bases = ['A', 'C', 'G', 'T'];
    let mut names = Vec::with_capacity(SNV_TRINUCLEOTIDE_BUCKET_COUNT);
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
pub fn swap_dna_base(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// SnvSigUtils.variantContext: the base change on the C or T strand, with its context. Err where Java
/// throws StringIndexOutOfBoundsException (an empty allele, or a context of fewer than 3 bases with an A or
/// G reference).
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

/// SnvSigUtils.loadSnvSignatures(names): the resource's matrix (96 rows, one column per signature) and the
/// signature names in column order
pub fn load_snv_signatures() -> Result<(Matrix, Vec<String>)> {
    let lines: Vec<String> = SNV_COSMIC_SIGNATURES.lines().map(String::from).collect();
    let mut names = Vec::new();
    let m = load_matrix_lines(&lines, &mut names, Some(&[]), false)?;
    Ok((m, names))
}

/// HumanChromosome.values(): 1-22, X, Y
pub const HUMAN_CHROMOSOMES: [&str; 24] = crate::regions::HUMAN_CHROMOSOMES;

/// RefGenomeVersion.versionedChromosome
pub fn versioned_chromosome(v38: bool, chromosome: &str) -> String {
    if v38 {
        crate::regions::enforce_chr_prefix(chromosome)
    } else {
        crate::regions::strip_chr_prefix(chromosome).to_string()
    }
}

/// PositionFrequencies.buildStandardChromosomeLengths: GRCh38 lengths, or for GRCh37 the larger of the two
/// versions' lengths (as the reference data of version 016 and before was built), keyed by versioned name
pub fn standard_chromosome_lengths(v38: bool) -> Vec<(String, i32)> {
    HUMAN_CHROMOSOMES
        .iter()
        .map(|c| {
            let l38 = crate::ref_coords::chromosome_length(c, false).unwrap();
            let len = if v38 { l38 } else { crate::ref_coords::chromosome_length(c, true).unwrap().max(l38) };
            (versioned_chromosome(v38, c), len)
        })
        .collect()
}

/// PositionFrequencies.initialisePositionCache: the first bucket index of each chromosome present in
/// `lengths`, and the number of buckets.
///
/// Java starts a chromosome at the previous one's last index plus one only if that index is above 0, so a
/// first chromosome of a single bucket is followed by one starting at index 0 as well; `java_start` keeps
/// that rule (SIGS finding SI2), false starts each chromosome after the previous one.
pub fn initialise_position_cache(
    v38: bool,
    bucket_size: i32,
    lengths: &[(String, i32)],
    java_start: bool,
) -> (BTreeMap<String, i32>, i32) {
    let mut index = BTreeMap::new();
    let mut cache_size = 0;
    if bucket_size == 0 {
        return (index, cache_size);
    }
    let mut last_end: i64 = -1;
    for c in HUMAN_CHROMOSOMES {
        let chr = versioned_chromosome(v38, c);
        let Some((_, length)) = lengths.iter().find(|(n, _)| *n == chr) else { continue };
        let start = if java_start {
            if last_end > 0 { last_end + 1 } else { 0 }
        } else {
            last_end + 1
        };
        index.insert(chr, start as i32);
        let count = (*length as f64 / bucket_size as f64).ceil() as i32;
        cache_size += count;
        last_end = start + count as i64 - 1;
    }
    (index, cache_size)
}

/// PositionFrequencies.getBucketIndex: -1 for a chromosome not in the index
pub fn bucket_index(bucket_size: i32, index: &BTreeMap<String, i32>, chromosome: &str, position: i32) -> i64 {
    match index.get(chromosome) {
        None => -1,
        Some(start) => *start as i64 + (position as f64 / bucket_size as f64).floor() as i64,
    }
}

/// PositionFrequencies.getChromosomeFromIndex. Java unboxes a missing chromosome's start index (a
/// NullPointerException); Err then.
pub fn chromosome_from_index(v38: bool, index: &BTreeMap<String, i32>, bucket: i32) -> Result<String> {
    let mut last_start = -1;
    let mut last_chr = String::new();
    for c in HUMAN_CHROMOSOMES {
        let chr = versioned_chromosome(v38, c);
        let start = *index.get(&chr).ok_or_else(|| anyhow!("java.lang.NullPointerException: no position index of {chr}"))?;
        if last_start >= 0 && bucket >= last_start && bucket < start {
            return Ok(last_chr);
        }
        last_start = start;
        last_chr = chr;
    }
    Ok(last_chr)
}

/// PositionFrequencies.getPositionFromIndex
pub fn position_from_index(index: &BTreeMap<String, i32>, chromosome: &str, bucket: i32, bucket_size: i32) -> Result<i32> {
    let start = *index.get(chromosome).ok_or_else(|| anyhow!("java.lang.NullPointerException: no position index of {chromosome}"))?;
    Ok((bucket - start).wrapping_mul(bucket_size))
}

/// PositionFrequencies with buildMap false (the counts array only), as CUPPA builds it
pub struct PositionCounts {
    pub v38: bool,
    pub bucket_size: i32,
    pub index: BTreeMap<String, i32>,
    pub counts: Vec<i32>,
}

impl PositionCounts {
    /// new PositionFrequencies(version, bucketSize, buildStandardChromosomeLengths(version), false)
    pub fn standard(v38: bool, bucket_size: i32) -> PositionCounts {
        let lengths = standard_chromosome_lengths(v38);
        let (index, size) = initialise_position_cache(v38, bucket_size, &lengths, true);
        PositionCounts { v38, bucket_size, index, counts: vec![0; size.max(0) as usize] }
    }

    pub fn is_valid_chromosome(&self, chromosome: &str) -> bool {
        self.index.contains_key(chromosome)
    }

    /// addPosition: a bucket outside the counts is ignored
    pub fn add_position(&mut self, chromosome: &str, position: i32) {
        let b = bucket_index(self.bucket_size, &self.index, chromosome, position);
        if b >= 0 && (b as usize) < self.counts.len() {
            self.counts[b as usize] += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn buckets_and_contexts() {
        let names = bucket_names();
        assert_eq!(names.len(), 96);
        assert_eq!(names[0], "C>A_ACA");
        assert_eq!(names[95], "T>G_TTT");
        assert_eq!(variant_context("G", "A", "AGC").unwrap(), "C>T_GCT");
        assert_eq!(variant_context("C", "T", "TCA").unwrap(), "C>T_TCA");
        assert!(variant_context("A", "C", "AA").is_err());
    }

    #[test]
    fn signatures_resource() {
        let (m, names) = load_snv_signatures().unwrap();
        assert_eq!((m.rows, m.cols), (96, 30));
        assert_eq!(names[0], "Sig1");
        assert_eq!(names[29], "Sig30");
    }

    #[test]
    fn position_cache() {
        let pc = PositionCounts::standard(true, 500_000);
        // ceil(length / 500 kb) per chromosome, chr1 first
        assert_eq!(pc.index["chr1"], 0);
        assert_eq!(pc.index["chr2"], 498);
        assert_eq!(pc.counts.len(), 6188);
        assert_eq!(chromosome_from_index(true, &pc.index, 497).unwrap(), "chr1");
        assert_eq!(chromosome_from_index(true, &pc.index, 498).unwrap(), "chr2");
        assert_eq!(chromosome_from_index(true, &pc.index, 6187).unwrap(), "chrY");
        assert_eq!(position_from_index(&pc.index, "chr2", 499, 500_000).unwrap(), 500_000);
        let pc37 = PositionCounts::standard(false, 500_000);
        assert!(pc37.is_valid_chromosome("1") && !pc37.is_valid_chromosome("chr1"));
    }
}
