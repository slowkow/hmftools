//! The five output files (EventsFile, EventsPerGeneFile, AllHaplotypeCombinationsFile, PeachGenotypeFile
//! with PeachGenotypeExtractor, QcStatusFile). Every file is sorted, so no output depends on Java's HashMap
//! order.

use std::collections::BTreeMap;

use anyhow::Result;
use hmf_common::java_sort::java_sort;
use hmf_common::java_text::compare;

use crate::algo::{EventCounts, HaplotypeAnalysis, HaplotypeCombination, GERMLINE_TOTAL_COPY_NUMBER};
use crate::loaders::{DrugInfoStore, HaplotypeFunctionStore};

pub const UNKNOWN_ALLELE_STRING: &str = "Unresolved Haplotype";
pub const DRUG_SEPARATOR: &str = ";";
pub const COMBINATION_SEPARATOR: &str = ";";

/// PeachUtils.convertCountToString
fn count_string(c: Option<i32>) -> String {
    c.map_or("UNKNOWN".to_string(), |c| c.to_string())
}

fn sorted(mut v: Vec<String>) -> Vec<String> {
    java_sort(&mut v, |a, b| compare(a, b));
    v
}

/// EventsFile.toLines
pub fn events_lines(counts: &EventCounts) -> Vec<String> {
    let mut lines = vec!["event\tcount".to_string()];
    lines.extend(sorted(counts.iter().map(|(e, c)| format!("{e}\t{}", count_string(*c))).collect()));
    lines
}

/// EventsPerGeneFile.toLines
pub fn events_per_gene_lines(analyses: &BTreeMap<String, HaplotypeAnalysis>) -> Vec<String> {
    let mut lines = vec!["gene\tevent\tcount".to_string()];
    let body: Vec<String> = analyses
        .iter()
        .flat_map(|(g, a)| a.event_counts.iter().map(move |(e, c)| format!("{g}\t{e}\t{}", count_string(*c))))
        .collect();
    lines.extend(sorted(body));
    lines
}

/// getHaplotypeCombinationString: "(name,count)" entries sorted and joined by ';'
pub fn combination_string(c: &HaplotypeCombination) -> String {
    sorted(c.haplotype_name_to_count.iter().map(|(n, k)| format!("({n},{k})")).collect()).join(COMBINATION_SEPARATOR)
}

/// genes sorted by name (Map.Entry.comparingByKey)
fn genes_sorted(analyses: &BTreeMap<String, HaplotypeAnalysis>) -> Vec<(&String, &HaplotypeAnalysis)> {
    let mut v: Vec<(&String, &HaplotypeAnalysis)> = analyses.iter().collect();
    java_sort(&mut v, |a, b| compare(a.0, b.0));
    v
}

/// AllHaplotypeCombinationsFile.toLines: per gene, the combinations by distance from two haplotypes, then by
/// the number of haplotypes other than the wild type, then by their text
pub fn all_combinations_lines(analyses: &BTreeMap<String, HaplotypeAnalysis>) -> Vec<String> {
    let mut lines = vec!["gene\tcombination\tcount\tnonWildTypeCount".to_string()];
    for (gene, a) in genes_sorted(analyses) {
        let wt = a.wild_type_haplotype_name.as_str();
        let mut combos: Vec<(i32, i32, String, &HaplotypeCombination)> = a
            .haplotype_combinations
            .iter()
            .map(|c| ((c.haplotype_count() - GERMLINE_TOTAL_COPY_NUMBER).abs(), c.haplotype_count_without(wt), combination_string(c), c))
            .collect();
        java_sort(&mut combos, |x, y| x.0.cmp(&y.0).then(x.1.cmp(&y.1)).then_with(|| compare(&x.2, &y.2)));
        for (_, non_wild, text, c) in combos {
            lines.push(format!("{gene}\t{text}\t{}\t{non_wild}", c.haplotype_count()));
        }
    }
    lines
}

/// PeachGenotype
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PeachGenotype {
    pub gene: String,
    pub allele: String,
    pub allele_count: i32,
    pub function: String,
    pub linked_drugs: String,
    pub url_prescription_info: String,
}

fn genotype(gene: &str, haplotype: &str, count: i32, drugs: Option<&DrugInfoStore>, functions: Option<&HaplotypeFunctionStore>) -> Result<PeachGenotype> {
    let function = match functions {
        None => String::new(),
        Some(f) => f.function(gene, haplotype)?.unwrap_or_default(),
    };
    let (linked_drugs, url_prescription_info) = match drugs {
        None => (String::new(), String::new()),
        Some(d) => {
            let names = d.sorted_drug_names(gene);
            let urls_names = d.sorted_drug_names(gene);
            let mut urls = Vec::with_capacity(urls_names.len());
            for n in &urls_names {
                urls.push(d.prescription_info_url(gene, n)?.unwrap_or_else(|| "null".to_string()));
            }
            (names.join(DRUG_SEPARATOR), urls.join(DRUG_SEPARATOR))
        }
    };
    Ok(PeachGenotype {
        gene: gene.to_string(),
        allele: haplotype.to_string(),
        allele_count: count,
        function,
        linked_drugs,
        url_prescription_info,
    })
}

/// PeachGenotypeExtractor.extract: per gene in name order, the best combination's haplotypes in name order,
/// or "Unresolved Haplotype" twice
pub fn extract_genotypes(
    analyses: &BTreeMap<String, HaplotypeAnalysis>,
    drugs: Option<&DrugInfoStore>,
    functions: Option<&HaplotypeFunctionStore>,
) -> Result<Vec<PeachGenotype>> {
    let mut out = Vec::new();
    for (gene, a) in genes_sorted(analyses) {
        match &a.best_haplotype_combination {
            Some(best) => {
                let mut entries: Vec<(&String, &i32)> = best.haplotype_name_to_count.iter().collect();
                java_sort(&mut entries, |x, y| compare(x.0, y.0));
                for (name, count) in entries {
                    out.push(genotype(gene, name, *count, drugs, functions)?);
                }
            }
            None => out.push(genotype(gene, UNKNOWN_ALLELE_STRING, GERMLINE_TOTAL_COPY_NUMBER, drugs, functions)?),
        }
    }
    Ok(out)
}

/// PeachGenotypeFile.toLines
pub fn genotype_lines(genotypes: &[PeachGenotype]) -> Vec<String> {
    let mut lines = vec!["gene\thaplotype\tcount\tfunction\tlinkedDrugs\tprescriptionUrls".to_string()];
    for g in genotypes {
        lines.push(format!("{}\t{}\t{}\t{}\t{}\t{}", g.gene, g.allele, g.allele_count, g.function, g.linked_drugs, g.url_prescription_info));
    }
    lines
}

/// QcStatusFile.toLines
pub fn qc_lines(analyses: &BTreeMap<String, HaplotypeAnalysis>) -> Vec<String> {
    let mut lines = vec!["gene\tstatus".to_string()];
    for (gene, a) in genes_sorted(analyses) {
        lines.push(format!("{gene}\t{}", a.qc_status.name()));
    }
    lines
}

/// Files.write: each line ended by a newline
pub fn write_lines(path: &str, lines: &[String]) -> Result<()> {
    let mut s = String::with_capacity(lines.iter().map(|l| l.len() + 1).sum());
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    std::fs::write(path, s).map_err(|e| anyhow::anyhow!("java.nio.file.NoSuchFileException: {path} ({e})"))
}
