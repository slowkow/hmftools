//! ClinvarAnnotation and ClinvarChrCache: ClinVar's clinical significance of variants it lists.

use anyhow::{Context, Result};
use std::collections::HashMap;

use hmf_common::regions::{is_human_chromosome, strip_chr_prefix};
use hmf_common::vcf::{self, InfoValue, VcfRecord};

use super::NeededPositions;
use crate::constants::{CLNSIG, CLNSIGCONF};

#[derive(Clone, Debug)]
pub struct ClinvarEntry {
    pub ref_: String,
    pub alt: String,
    pub significance: String,
    pub conflict: String,
}

/// The entries at the needed positions, by chromosome without the chr prefix
#[derive(Default)]
pub struct ClinvarData {
    by_chr: HashMap<String, HashMap<i32, Vec<ClinvarEntry>>>,
    pub entry_count: u64,
}

/// ClinvarChrCache.stripBrackets
fn strip_brackets(s: &str) -> String {
    s.chars().filter(|c| *c != '[' && *c != ']' && *c != ' ').collect()
}

impl ClinvarData {
    /// whether ClinVar has entries on the chromosome (Java's chromosome cache exists)
    pub fn has_chromosome(&self, chromosome: &str) -> bool {
        self.by_chr.contains_key(strip_chr_prefix(chromosome))
    }

    /// findClinvarEntry: the first entry at the position with the ref and alt
    pub fn find(&self, chromosome: &str, position: i32, ref_: &str, alt: &str) -> Option<&ClinvarEntry> {
        self.by_chr.get(strip_chr_prefix(chromosome))?.get(&position)?.iter().find(|e| e.ref_ == ref_ && e.alt == alt)
    }

    /// ClinvarChrCache.addEntry
    #[cfg(test)]
    pub fn add_entry(&mut self, chromosome: &str, position: i32, ref_: &str, alt: &str, significance: &str, conflict: &str) {
        let e = ClinvarEntry { ref_: ref_.to_string(), alt: alt.to_string(), significance: strip_brackets(significance), conflict: strip_brackets(conflict) };
        self.by_chr.entry(strip_chr_prefix(chromosome).to_string()).or_default().entry(position).or_default().push(e);
    }

    /// ClinvarChrCache.annotateVariant
    pub fn annotate(&self, rec: &mut VcfRecord, chromosome: &str, position: i32, ref_: &str, alt: &str) -> Result<()> {
        if let Some(e) = self.find(chromosome, position, ref_, alt) {
            super::put_once(rec, CLNSIG, InfoValue::Text(e.significance.clone()))?;
            if !e.conflict.is_empty() {
                super::put_once(rec, CLNSIGCONF, InfoValue::Text(e.conflict.clone()))?;
            }
        }
        Ok(())
    }
}

/// ClinvarAnnotation.loadEntries: records with an alt on human chromosomes and a significance or a
/// conflict; a chromosome seen again later starts its entries again, as Java's cache does. Only the
/// needed positions are kept, but every record is counted.
pub fn load(path: &str, needed: &NeededPositions) -> Result<ClinvarData> {
    let (header, mut reader, first) = vcf::read_header(path).with_context(|| format!("failed to read Clinvar VCF file({path})"))?;
    let mut data = ClinvarData::default();
    let mut current_chr = String::new();
    let mut current_needed: Vec<i32> = Vec::new();
    let mut pending = first;
    loop {
        let line: String = match pending.take() {
            Some(l) => l,
            None => match reader.next_line()? {
                Some(l) => String::from_utf8_lossy(l).into_owned(),
                None => break,
            },
        };
        if line.is_empty() {
            continue;
        }
        let mut cols = line.splitn(9, '\t');
        let (Some(contig), Some(pos_str), Some(_id), Some(_ref), Some(alt)) = (cols.next(), cols.next(), cols.next(), cols.next(), cols.next())
        else {
            anyhow::bail!("failed to read Clinvar VCF file: short line: {line}");
        };
        if alt == "." {
            continue; // fewer than 2 alleles
        }
        if !is_human_chromosome(contig) {
            continue;
        }
        let chromosome = strip_chr_prefix(contig);
        if chromosome != current_chr {
            current_chr = chromosome.to_string();
            data.by_chr.insert(current_chr.clone(), HashMap::new());
            // the needed positions are by the VCF's chromosome names
            let mut n: Vec<i32> = Vec::new();
            for name in [chromosome.to_string(), format!("chr{chromosome}")] {
                if let Some(v) = needed.get(&name) {
                    n.extend_from_slice(v);
                }
            }
            n.sort_unstable();
            n.dedup();
            current_needed = n;
        }
        let Some(position) = hmf_common::java_text::parse_int(pos_str) else {
            anyhow::bail!("failed to read Clinvar VCF file: invalid position in line: {line}");
        };
        if current_needed.binary_search(&position).is_err() {
            // counted as Java does for entries with a significance or conflict
            if line.contains("CLNSIG") {
                data.entry_count += 1;
            }
            continue;
        }
        let rec = VcfRecord::parse(&line, &header)?;
        let significance = rec.attr_string(CLNSIG, "");
        let conflict = rec.attr_string(CLNSIGCONF, "");
        if significance.is_empty() && conflict.is_empty() {
            continue;
        }
        data.entry_count += 1;
        data.by_chr.get_mut(&current_chr).unwrap().entry(position).or_default().push(ClinvarEntry {
            ref_: rec.ref_allele.clone(),
            alt: rec.alt().to_string(),
            significance: strip_brackets(&significance),
            conflict: strip_brackets(&conflict),
        });
    }
    crate::info!("loaded {} Clinvar entries from file({})", data.entry_count, path);
    Ok(data)
}
