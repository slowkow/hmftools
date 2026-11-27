//! The other inputs: the drugs and haplotype function files (DrugInfoLoader, HaplotypeFunctionLoader, with
//! DrugInfoStore and HaplotypeFunctionStore) and the events of the VCF (HaplotypeEventLoader).

use std::collections::{BTreeMap, BTreeSet, HashSet};

use anyhow::{anyhow, bail, Context, Result};
use hmf_common::variant_context::{GenotypeType, VcfReader, VcfRecord};

use crate::algo::EventCounts;
use crate::compat;
use crate::panel::{Chromosome, VariantHaplotypeEvent};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DrugInfo {
    pub drug_name: String,
    pub gene_name: String,
    pub prescription_info_url: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HaplotypeFunction {
    pub gene_name: String,
    pub haplotype_name: String,
    pub function: String,
}

/// the fields of a line split at tabs with trailing empty fields dropped (String.split), by column name
fn tsv_rows(filename: &str, columns: &[&str]) -> Result<Vec<Vec<String>>> {
    let lines = hmf_common::java_text::read_lines(filename)?;
    let header = lines.first().ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
    let fields = hmf_common::java_text::fields_index_map(header, '\t');
    let mut out = Vec::new();
    for line in &lines[1..] {
        let values = hmf_common::java_text::split(line, '\t');
        let mut row = Vec::with_capacity(columns.len());
        for c in columns {
            let i = *fields.get(*c).ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {c}"))?;
            let v = values
                .get(i)
                .ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", values.len()))?;
            row.push(v.to_string());
        }
        out.push(row);
    }
    Ok(out)
}

/// DrugInfoLoader.loadDrugInfos
pub fn load_drug_infos(filename: &str) -> Result<Vec<DrugInfo>> {
    let rows = tsv_rows(filename, &["drugName", "gene", "urlPrescriptionInfo"])
        .with_context(|| format!("java.lang.RuntimeException: failed to load drugs TSV: {filename}"))?;
    let infos: Vec<DrugInfo> = rows
        .into_iter()
        .map(|mut r| DrugInfo { prescription_info_url: r.pop().unwrap(), gene_name: r.pop().unwrap(), drug_name: r.pop().unwrap() })
        .collect();
    hmf_common::info!("loaded {} lines of drug info from file ({})", infos.len(), filename);
    Ok(infos)
}

/// HaplotypeFunctionLoader.loadFunctions. Java's message for a file it cannot read names the drugs file
/// (PE2); corrected, the haplotype function file.
pub fn load_functions(filename: &str) -> Result<Vec<HaplotypeFunction>> {
    let what = if compat::java("PE2") { "drugs" } else { "haplotype functions" };
    let rows = tsv_rows(filename, &["gene", "haplotype", "function"])
        .with_context(|| format!("java.lang.RuntimeException: failed to load {what} TSV: {filename}"))?;
    let functions: Vec<HaplotypeFunction> = rows
        .into_iter()
        .map(|mut r| HaplotypeFunction { function: r.pop().unwrap(), haplotype_name: r.pop().unwrap(), gene_name: r.pop().unwrap() })
        .collect();
    hmf_common::info!("loaded {} lines of haplotype functions from file ({})", functions.len(), filename);
    Ok(functions)
}

/// DrugInfoStore
pub struct DrugInfoStore {
    pub drug_infos: Vec<DrugInfo>,
}

impl DrugInfoStore {
    /// getRelevantDrugNames, sorted (as the extractor uses them)
    pub fn sorted_drug_names(&self, gene: &str) -> Vec<String> {
        let set: BTreeSet<&str> = self.drug_infos.iter().filter(|d| d.gene_name == gene).map(|d| d.drug_name.as_str()).collect();
        let mut v: Vec<String> = set.into_iter().map(String::from).collect();
        hmf_common::java_sort::java_sort(&mut v, |a, b| hmf_common::java_text::compare(a, b));
        v
    }
    /// getPrescriptionInfoUrl: None without a match, an error for more than one
    pub fn prescription_info_url(&self, gene: &str, drug: &str) -> Result<Option<String>> {
        let m: Vec<&DrugInfo> = self.drug_infos.iter().filter(|d| d.gene_name == gene && d.drug_name == drug).collect();
        match m.len() {
            0 => Ok(None),
            1 => Ok(Some(m[0].prescription_info_url.clone())),
            _ => bail!("java.lang.IllegalStateException: Multiple urls configured for drug {drug} with gene {gene}"),
        }
    }
}

/// HaplotypeFunctionStore
pub struct HaplotypeFunctionStore {
    pub functions: Vec<HaplotypeFunction>,
}

impl HaplotypeFunctionStore {
    pub fn function(&self, gene: &str, haplotype: &str) -> Result<Option<String>> {
        let m: Vec<&HaplotypeFunction> = self.functions.iter().filter(|f| f.gene_name == gene && f.haplotype_name == haplotype).collect();
        match m.len() {
            0 => Ok(None),
            1 => Ok(Some(m[0].function.clone())),
            _ => bail!("java.lang.IllegalStateException: Multiple functions configured for haplotype {haplotype} for gene {gene}"),
        }
    }
}

/// HaplotypeEventLoader.getEventCount
fn event_count(t: GenotypeType, event_id: &str) -> Result<Option<i32>> {
    Ok(match t {
        GenotypeType::NoCall => None,
        GenotypeType::HomRef => Some(0),
        GenotypeType::Het => Some(1),
        GenotypeType::HomVar => Some(2),
        GenotypeType::Mixed | GenotypeType::Unavailable => {
            let name = if t == GenotypeType::Mixed { "MIXED" } else { "UNAVAILABLE" };
            bail!("java.lang.IllegalStateException: cannot get occurrence count for event with ID '{event_id}' and genotypeType '{name}'")
        }
    })
}

/// VariantHaplotypeEvent.fromVariantContext: the record's single alternate allele as an event
fn event_from_record(rec: &VcfRecord) -> Result<VariantHaplotypeEvent> {
    let chromosome = Chromosome::from_string(&rec.chrom)?;
    if rec.alts.len() > 1 {
        bail!(
            "java.lang.IllegalArgumentException: Cannot handle variant with multiple alts: '{}:{}{}>...'",
            chromosome.name(),
            rec.pos,
            rec.ref_allele
        );
    }
    let alt = rec
        .alts
        .first()
        .ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0 (no alternate allele at {}:{})", rec.chrom, rec.pos))?;
    Ok(VariantHaplotypeEvent::new(chromosome, rec.pos, &rec.ref_allele, alt))
}

/// HaplotypeEventLoader.loadRelevantVariantHaplotypeEvents: the count of each event of the VCF whose
/// reference bases cover a panel position, for records that are not filtered.
///
/// Java turns every record into an event and reads its genotype before it checks whether the record is
/// relevant, so a record anywhere in the VCF on another contig than 1-22, X and Y, with more than one
/// alternate allele, without one, or with a genotype of the sample that is mixed or unavailable ends the
/// run (PE1). Corrected, records on other contigs or that cover no panel position are skipped first.
pub fn load_relevant_events(vcf: &str, sample: &str, relevant: &BTreeMap<Chromosome, HashSet<i32>>) -> Result<EventCounts> {
    let run = || -> Result<EventCounts> {
        let mut reader = VcfReader::open(vcf)?;
        let header = reader.header.clone();
        if !header.samples.iter().any(|s| s == sample) && !compat::java("PE1") {
            bail!("sample {sample} is not in the VCF {vcf}");
        }
        let empty = HashSet::new();
        let mut counts = EventCounts::new();
        while let Some(line) = reader.next_line()? {
            let mut cols = line.splitn(6, '\t');
            let chrom = cols.next().unwrap_or("");
            if !compat::java("PE1") {
                // skipped unparsed when it cannot cover a panel position
                let Ok(c) = Chromosome::from_string(chrom) else { continue };
                let Some(positions) = relevant.get(&c) else { continue };
                let pos = cols.next().and_then(hmf_common::java_text::parse_int);
                let ref_len = cols.nth(1).map(|r| r.len() as i32);
                if let (Some(p), Some(l)) = (pos, ref_len) {
                    if !(p..p + l).any(|x| positions.contains(&x)) {
                        continue;
                    }
                }
            }
            let mut rec = VcfRecord::parse(line, &header)?;
            let positions = match Chromosome::from_string(&rec.chrom) {
                Ok(c) => relevant.get(&c).unwrap_or(&empty),
                Err(e) => return Err(e),
            };
            let is_candidate = !rec.is_filtered() && (rec.pos..rec.pos + rec.ref_allele.chars().count() as i32).any(|x| positions.contains(&x));
            if !compat::java("PE1") && !is_candidate {
                continue;
            }
            let event = event_from_record(&rec)?;
            let id = event.id();
            rec.decode(&header)?;
            let genotype = rec
                .decoded_genotype(sample)
                .ok_or_else(|| anyhow!("java.lang.NullPointerException: no genotype of {sample} at {}:{}", rec.chrom, rec.pos))?;
            let count = event_count(genotype.genotype_type(), &id)?;
            if is_candidate {
                if counts.contains_key(&id) {
                    bail!("java.lang.IllegalStateException: encountered event with ID '{id}' more than once in input VCF");
                }
                counts.insert(id, count);
            }
        }
        Ok(counts)
    };
    run()
}
