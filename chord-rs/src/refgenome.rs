//! The reference genome as CHORD reads it through hmf-common's RefGenomeSource and htsjdk 2.23's
//! IndexedFastaSequenceFile:
//! - `getBases(chr, start, end)` gives null for start < 1 or end < start, throws for an unknown contig or
//!   an end past the contig, and returns the bases as stored (soft-masked bases stay lower case);
//! - `getChromosomeLength(chr)` reads the sequence dictionary (`<name>.dict` or `<file>.dict`), which
//!   htsjdk does not build from the .fai: without one, it throws.
//!
//! The corrected behaviour (CH4, CH6) upper-cases the bases and takes the lengths from the .fai.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;

use crate::compat;

pub struct RefGenomeSource {
    genome: hmf_common::refgenome::RefGenome,
    /// the .fai lengths
    lengths: HashMap<String, i32>,
    /// the .dict lengths, if htsjdk finds a dictionary
    dict_lengths: Option<HashMap<String, i32>>,
}

/// htsjdk ReferenceSequenceFileFactory.FASTA_EXTENSIONS
const FASTA_EXTENSIONS: [&str; 8] = [".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz", ".txt", ".txt.gz"];

/// AbstractFastaSequenceFile.findSequenceDictionary: the FASTA extension replaced by .dict, else .dict
/// appended
pub fn find_dictionary(path: &str) -> Option<String> {
    let mut candidates = Vec::new();
    match FASTA_EXTENSIONS.iter().find(|e| path.ends_with(*e)) {
        Some(ext) => candidates.push(format!("{}.dict", &path[..path.len() - ext.len()])),
        None => candidates.push(format!("{path}.dict")),
    }
    candidates.push(format!("{path}.dict"));
    candidates.into_iter().find(|c| std::path::Path::new(c).exists())
}

/// the @SQ lines of a sequence dictionary: SN and LN
fn read_dictionary(path: &str) -> Result<HashMap<String, i32>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading sequence dictionary {path}"))?;
    let mut out = HashMap::new();
    for line in text.lines() {
        if !line.starts_with("@SQ") {
            continue;
        }
        let mut name = None;
        let mut len = None;
        for field in line.split('\t').skip(1) {
            if let Some(v) = field.strip_prefix("SN:") {
                name = Some(v.to_string());
            } else if let Some(v) = field.strip_prefix("LN:") {
                len = v.parse::<i32>().ok();
            }
        }
        if let (Some(n), Some(l)) = (name, len) {
            out.entry(n).or_insert(l);
        }
    }
    Ok(out)
}

impl RefGenomeSource {
    /// RefGenomeSource.loadRefGenome: Err where Java logs "Reference file loading failed" and returns null
    pub fn load(path: &str) -> Result<RefGenomeSource> {
        if !std::path::Path::new(path).is_file() {
            bail!("java.io.FileNotFoundException: {path}");
        }
        if !std::path::Path::new(&format!("{path}.fai")).is_file() && !path.ends_with(".gz") {
            bail!("htsjdk.samtools.SAMException: Unable to find index for {path}");
        }
        let genome = hmf_common::refgenome::RefGenome::open(path).with_context(|| format!("opening reference {path}"))?;
        let lengths = genome.sequence_names().iter().map(|n| (n.clone(), genome.chromosome_length(n))).collect();
        let dict_lengths = match find_dictionary(path) {
            Some(d) => Some(read_dictionary(&d)?),
            None => None,
        };
        Ok(RefGenomeSource { genome, lengths, dict_lengths })
    }

    /// RefGenomeSource.getBases: None where Java returns null (start < 1 or end < start)
    pub fn get_bases(&self, chromosome: &str, start: i32, end: i32) -> Result<Option<Vec<u8>>> {
        if start < 1 || end < start {
            return Ok(None);
        }
        let Some(&length) = self.lengths.get(chromosome) else {
            bail!("htsjdk.samtools.SAMException: Unable to find entry for contig: {chromosome}");
        };
        if end > length {
            bail!("htsjdk.samtools.SAMException: Query asks for data past end of contig. Query contig {chromosome} start:{start} stop:{end} contigLength:{length}");
        }
        let mut bases = self.genome.bases(chromosome, start, end)?;
        if !compat::java("CH4") {
            bases.make_ascii_uppercase();
        }
        Ok(Some(bases))
    }

    /// The bases from start to end, with the positions outside the contig as 'N' (the corrected reading
    /// of CH5); an unknown contig is still an error
    pub fn get_bases_padded(&self, chromosome: &str, start: i32, end: i32) -> Result<Vec<u8>> {
        let length = *self.lengths.get(chromosome).ok_or_else(|| anyhow!("htsjdk.samtools.SAMException: Unable to find entry for contig: {chromosome}"))?;
        let from = start.max(1);
        let to = end.min(length);
        let mut out = Vec::with_capacity((end - start + 1).max(0) as usize);
        for _ in start..from.min(end + 1) {
            out.push(b'N');
        }
        if from <= to {
            out.extend(self.get_bases(chromosome, from, to)?.unwrap_or_default());
        }
        for _ in (to + 1).max(start)..=end {
            out.push(b'N');
        }
        Ok(out)
    }

    /// RefGenomeSource.getChromosomeLength: from the sequence dictionary in Java (CH6), else the .fai
    pub fn chromosome_length(&self, chromosome: &str) -> Result<i32> {
        if compat::java("CH6") {
            let Some(dict) = &self.dict_lengths else {
                bail!(
                    "java.lang.NullPointerException: Cannot invoke \"htsjdk.samtools.SAMSequenceDictionary.getSequence(String)\" because the return value of \"htsjdk.samtools.reference.IndexedFastaSequenceFile.getSequenceDictionary()\" is null"
                );
            };
            return dict.get(chromosome).copied().ok_or_else(|| {
                anyhow!("java.lang.NullPointerException: Cannot invoke \"htsjdk.samtools.SAMSequenceRecord.getSequenceLength()\" because the return value of \"htsjdk.samtools.SAMSequenceDictionary.getSequence(String)\" is null")
            });
        }
        self.lengths.get(chromosome).copied().ok_or_else(|| anyhow!("htsjdk.samtools.SAMException: Unable to find entry for contig: {chromosome}"))
    }
}
