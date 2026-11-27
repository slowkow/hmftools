//! HaplotypeCaller, HaplotypeCombination, PeachAlgo and HaplotypeAnalysis: the combinations of haplotypes
//! that explain a gene's events, the best of them, and the QC status.

use std::collections::BTreeMap;

use anyhow::{bail, Result};

use crate::panel::{DefaultHaplotype, HaplotypePanel, NonDefaultHaplotype};

pub const GERMLINE_TOTAL_COPY_NUMBER: i32 = 2;

/// event id to its count in the sample: 0, 1 or 2, or None for a no-call (Java's null)
pub type EventCounts = BTreeMap<String, Option<i32>>;

/// HaplotypeCombination: haplotype name to how often it is called
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HaplotypeCombination {
    pub haplotype_name_to_count: BTreeMap<String, i32>,
}

impl HaplotypeCombination {
    pub fn new(haplotype_name_to_count: BTreeMap<String, i32>) -> Result<HaplotypeCombination> {
        if haplotype_name_to_count.values().any(|c| *c <= 0) {
            bail!("java.lang.RuntimeException: Invalid haplotype combination with non-positive count(s): {haplotype_name_to_count:?}");
        }
        if haplotype_name_to_count.values().sum::<i32>() < GERMLINE_TOTAL_COPY_NUMBER {
            bail!("java.lang.RuntimeException: Less than {GERMLINE_TOTAL_COPY_NUMBER} haplotype calls in combination is not allowed: {haplotype_name_to_count:?}");
        }
        Ok(HaplotypeCombination { haplotype_name_to_count })
    }

    /// for tests: a combination from (name, count) pairs
    #[cfg(test)]
    pub fn of(pairs: &[(&str, i32)]) -> HaplotypeCombination {
        HaplotypeCombination::new(pairs.iter().map(|(n, c)| (n.to_string(), *c)).collect()).unwrap()
    }

    pub fn haplotype_count(&self) -> i32 {
        self.haplotype_name_to_count.values().sum()
    }

    pub fn haplotype_count_without(&self, ignored: &str) -> i32 {
        self.haplotype_name_to_count.iter().filter(|(n, _)| n.as_str() != ignored).map(|(_, c)| *c).sum()
    }
}

/// HaplotypeCaller.findPossibleHaplotypeCombinations
pub fn find_possible_haplotype_combinations(
    event_counts: &EventCounts,
    non_default: &[NonDefaultHaplotype],
    default: &DefaultHaplotype,
) -> Result<Vec<HaplotypeCombination>> {
    // the counts are known here, or there is no combination
    if event_counts.values().any(|c| c.is_none()) {
        return Ok(Vec::new());
    }
    let counts: BTreeMap<String, i32> = event_counts.iter().map(|(k, v)| (k.clone(), v.unwrap())).collect();
    let mut out = Vec::new();
    for list in possible_non_default(&counts, non_default, 0)? {
        out.push(construct_combination(&list, non_default, default)?);
    }
    Ok(out)
}

/// getPossibleNonDefaultHaplotypes by recursive descent: the lists (as indices into `candidates`, each list
/// in the order Java builds it, last try first) of haplotypes from `from` on that explain every count
fn possible_non_default(unexplained: &BTreeMap<String, i32>, candidates: &[NonDefaultHaplotype], from: usize) -> Result<Vec<Vec<usize>>> {
    if unexplained.values().all(|c| *c == 0) {
        return Ok(vec![Vec::new()]);
    }
    if let Some((event, _)) = unexplained.iter().find(|(_, c)| **c < 0) {
        let names: Vec<&str> = candidates[from..].iter().map(|h| h.name.as_str()).collect();
        bail!("java.lang.IllegalStateException: Negative count encountered for event '{event}' for candidate haplotypes ({})", names.join(", "));
    }
    let mut out = Vec::new();
    for i in from..candidates.len() {
        let h = &candidates[i];
        let possible = h.events.iter().all(|e| {
            let id = e.id();
            unexplained.get(&id).copied().unwrap_or(0) >= h.matching_count(&id)
        });
        if possible {
            let after: BTreeMap<String, i32> = unexplained.iter().map(|(e, c)| (e.clone(), c - h.matching_count(e))).collect();
            for mut sub in possible_non_default(&after, candidates, i)? {
                sub.push(i);
                out.push(sub);
            }
        }
    }
    Ok(out)
}

fn construct_combination(list: &[usize], candidates: &[NonDefaultHaplotype], default: &DefaultHaplotype) -> Result<HaplotypeCombination> {
    let mut counts: BTreeMap<String, i32> = BTreeMap::new();
    for &i in list {
        *counts.entry(candidates[i].name.clone()).or_insert(0) += 1;
    }
    if (list.len() as i32) < GERMLINE_TOTAL_COPY_NUMBER {
        counts.insert(default.name.clone(), GERMLINE_TOTAL_COPY_NUMBER - list.len() as i32);
    }
    HaplotypeCombination::new(counts)
}

/// PeachQCStatus
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeachQcStatus {
    Pass,
    FailNoCombinationFound,
    FailNoUniqueBestCombinationFound,
    FailEventWithUnknownCount,
    WarnTooManyAllelesFound,
}

impl PeachQcStatus {
    pub fn name(self) -> &'static str {
        match self {
            PeachQcStatus::Pass => "PASS",
            PeachQcStatus::FailNoCombinationFound => "FAIL_NO_COMBINATION_FOUND",
            PeachQcStatus::FailNoUniqueBestCombinationFound => "FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND",
            PeachQcStatus::FailEventWithUnknownCount => "FAIL_EVENT_WITH_UNKNOWN_COUNT",
            PeachQcStatus::WarnTooManyAllelesFound => "WARN_TOO_MANY_ALLELES_FOUND",
        }
    }
    /// hasBestHaplotypeCombination
    pub fn has_best_haplotype_combination(self) -> bool {
        matches!(self, PeachQcStatus::Pass | PeachQcStatus::WarnTooManyAllelesFound)
    }
}

/// HaplotypeAnalysis
#[derive(Clone, Debug)]
pub struct HaplotypeAnalysis {
    pub event_counts: EventCounts,
    pub haplotype_combinations: Vec<HaplotypeCombination>,
    #[cfg_attr(not(test), allow(dead_code))]
    pub default_haplotype_name: String,
    pub wild_type_haplotype_name: String,
    pub qc_status: PeachQcStatus,
    pub best_haplotype_combination: Option<HaplotypeCombination>,
}

/// selectBestCombinations: those nearest to two haplotypes, then of those the ones with the fewest
/// haplotypes other than the wild type
fn select_best_combinations(combinations: &[HaplotypeCombination], wild_type: &str) -> Vec<HaplotypeCombination> {
    if combinations.is_empty() {
        return Vec::new();
    }
    let deviation = |c: &HaplotypeCombination| (c.haplotype_count() - GERMLINE_TOTAL_COPY_NUMBER).abs();
    let min_dev = combinations.iter().map(deviation).min().unwrap();
    let candidates: Vec<HaplotypeCombination> = combinations.iter().filter(|c| deviation(c) == min_dev).cloned().collect();
    if candidates.len() <= 1 {
        return candidates;
    }
    let min_non_wild = candidates.iter().map(|c| c.haplotype_count_without(wild_type)).min().unwrap();
    candidates.into_iter().filter(|c| c.haplotype_count_without(wild_type) == min_non_wild).collect()
}

fn determine_qc_status(events: &EventCounts, best: &[HaplotypeCombination]) -> PeachQcStatus {
    if events.values().any(|c| c.is_none()) {
        PeachQcStatus::FailEventWithUnknownCount
    } else if best.is_empty() {
        PeachQcStatus::FailNoCombinationFound
    } else if best.len() != 1 {
        PeachQcStatus::FailNoUniqueBestCombinationFound
    } else if best[0].haplotype_count() > 2 {
        PeachQcStatus::WarnTooManyAllelesFound
    } else {
        PeachQcStatus::Pass
    }
}

/// PeachAlgo.getGeneToHaplotypeAnalysis
pub fn gene_to_haplotype_analysis(panel: &HaplotypePanel, event_counts: &EventCounts) -> Result<BTreeMap<String, HaplotypeAnalysis>> {
    let mut out = BTreeMap::new();
    for (gene, gene_panel) in &panel.genes {
        hmf_common::info!("handle gene: {}", gene);
        let mut relevant = EventCounts::new();
        for (id, count) in event_counts {
            if panel.is_relevant_for_id(id, gene)? {
                relevant.insert(id.clone(), *count);
            }
        }
        let combinations = find_possible_haplotype_combinations(&relevant, &gene_panel.non_default_haplotypes, &gene_panel.default_haplotype)?;
        let best = select_best_combinations(&combinations, &gene_panel.wild_type_haplotype_name);
        let qc_status = determine_qc_status(&relevant, &best);
        let best_haplotype_combination = if qc_status.has_best_haplotype_combination() { Some(best[0].clone()) } else { None };
        out.insert(
            gene.clone(),
            HaplotypeAnalysis {
                event_counts: relevant,
                haplotype_combinations: combinations,
                default_haplotype_name: gene_panel.default_haplotype.name.clone(),
                wild_type_haplotype_name: gene_panel.wild_type_haplotype_name.clone(),
                qc_status,
                best_haplotype_combination,
            },
        );
    }
    Ok(out)
}
