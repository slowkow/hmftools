//! The smaller sections: LILAC (LilacInterpreter), immune escape (ImmuneEscapeInterpreter), Virus Interpreter
//! (VirusInterpreter), CHORD, CUPPA (CuppaDataFactory), PEACH, SIGS (SigsInterpreter) and the configured
//! primary tumor (OrangeAlgo.formConfiguredPrimaryTumorDoid).

use anyhow::{anyhow, bail, Result};
use hmf_common::hla;
use hmf_common::tool_results::{self as tr, CuppaPredictionEntry};

use crate::compat::java;
use crate::config::OrangeConfig;
use crate::model::*;

fn dir_exists(d: &Option<String>) -> Option<&str> {
    d.as_deref().filter(|d| std::path::Path::new(d).exists())
}

// ---------------------------------------------------------------------------------------------
// LILAC

/// OrangeAlgo.loadLilacData and LilacInterpreter.build
pub fn load_lilac(config: &OrangeConfig) -> Result<Option<LilacRecord>> {
    let Some(dir) = dir_exists(&config.lilac_dir) else {
        hmf_common::debug!("skipping loading Lilac results since input dir not provided");
        return Ok(None);
    };
    let qc = hla::read_lilac_qc(&hla::lilac_qc_filename(dir, &config.tumor_id))?;
    let alleles = hla::read_lilac_alleles(&hla::lilac_allele_filename(dir, &config.tumor_id))?;
    let mut out = Vec::new();
    for a in &alleles {
        let q = qc
            .iter()
            .find(|q| q.genes == a.genes)
            .ok_or_else(|| anyhow!("java.lang.NullPointerException: no LILAC QC line for genes {}", a.genes))?;
        out.push(LilacAllele {
            gene_class: a.genes.clone(),
            allele: a.allele.clone(),
            qc_status: q.status.clone(),
            tumor_copy_number: a.tumor_copy_number,
            ref_fragments: if config.has_reference() { Some(a.ref_fragments) } else { None },
            tumor_fragments: a.tumor_fragments,
            rna_fragments: if config.has_rna() { Some(a.rna_fragments) } else { None },
            somatic_missense: a.somatic_missense,
            somatic_nonsense_or_frameshift: a.somatic_nonsense_or_frameshift,
            somatic_splice: a.somatic_splice,
            somatic_synonymous: a.somatic_synonymous,
            somatic_inframe_indel: a.somatic_inframe_indel,
        });
    }
    Ok(Some(LilacRecord { alleles: out }))
}

// ---------------------------------------------------------------------------------------------
// immune escape

/// ImmuneEscapeInterpreter.interpret. OR9: Java's antigen presentation genes list "TABBP" for TAPBP, so that
/// gene is never checked; corrected, TAPBP.
pub fn immune_escape(purple: &PurpleRecord, linx: &LinxRecord) -> ImmuneEscapeRecord {
    const HLA: [&str; 3] = ["HLA-A", "HLA-B", "HLA-C"];
    let tapbp = if java("OR9") { "TABBP" } else { "TAPBP" };
    let antigen: [&str; 8] = ["B2M", "CALR", "TAP1", "TAP2", tapbp, "NLRC5", "CIITA", "RFX5"];
    const IFN: [&str; 7] = ["JAK1", "JAK2", "IRF2", "IFNGR1", "IFNGR2", "APLNR", "STAT1"];

    let has_loh = |genes: &[&str]| {
        genes.iter().any(|g| purple.somatic_drivers.iter().any(|d| d.gene == *g && ["LOH", "DEL", "HET_DEL"].contains(&d.driver_type.as_str())))
    };
    // the first canonical gain or deletion of the gene decides
    let first_canonical = |gene: &str| purple.somatic_gains_dels.iter().find(|g| g.driver.gene == gene && g.driver.is_canonical);
    let is_deleted = |gene: &str| first_canonical(gene).is_some_and(|g| g.driver.driver_type == "DEL");
    let is_amplified = |gene: &str| first_canonical(gene).is_some_and(|g| g.driver.driver_type == "AMP");
    let hom_disrupted = |gene: &str| linx.somatic_breakends.iter().any(|b| b.gene == gene && b.undisrupted_copy_number < 0.5);
    let inactivated = |genes: &[&str]| genes.iter().any(|g| is_deleted(g) || hom_disrupted(g));

    ImmuneEscapeRecord {
        has_hla_escape: has_loh(&HLA) || inactivated(&HLA),
        has_antigen_presentation_pathway_escape: inactivated(&antigen),
        has_ifn_gamma_pathway_escape: inactivated(&IFN),
        has_pdl1_overexpression_escape: is_amplified("CD274"),
        has_cd58_inactivation_escape: inactivated(&["CD58"]),
        has_epigenetic_setdb1_escape: is_amplified("SETDB1"),
    }
}

// ---------------------------------------------------------------------------------------------
// Virus Interpreter

fn convert_virus(v: &tr::AnnotatedVirus) -> VirusInterpreterEntry {
    VirusInterpreterEntry {
        name: v.name.clone(),
        qc_status: v.qc_status.clone(),
        integrations: v.integrations,
        interpretation: v.interpretation.clone(),
        percentage_covered: v.percentage_covered,
        mean_coverage: v.mean_coverage,
        expected_clonal_coverage: v.expected_clonal_coverage,
        reported: v.reported,
        driver_interpretation: v.driver_likelihood.clone(),
    }
}

/// VirusInterpreter.filterBlacklistedViruses: Err for a virus without a blacklist status
pub(crate) fn filter_blacklisted(viruses: &[&tr::AnnotatedVirus]) -> Result<Vec<VirusInterpreterEntry>> {
    if let Some(v) = viruses.iter().find(|v| v.blacklisted.is_none()) {
        bail!("java.lang.RuntimeException: Encountered virus '{}' with unknown blacklist status", v.name);
    }
    Ok(viruses.iter().filter(|v| v.blacklisted == Some(false)).map(|v| convert_virus(v)).collect())
}

/// OrangeAlgo.loadVirusInterpreterData and VirusInterpreter.interpret
pub fn load_virus(config: &OrangeConfig) -> Result<Option<VirusInterpreterData>> {
    let Some(dir) = dir_exists(&config.virus_dir) else { return Ok(None) };
    let all = tr::read_annotated_viruses(&tr::annotated_virus_filename(dir, &config.tumor_id))?;
    let all_refs: Vec<&tr::AnnotatedVirus> = all.iter().collect();
    let reportable: Vec<&tr::AnnotatedVirus> = all.iter().filter(|v| v.reported).collect();
    Ok(Some(VirusInterpreterData { all_viruses: filter_blacklisted(&all_refs)?, reportable_viruses: filter_blacklisted(&reportable)? }))
}

// ---------------------------------------------------------------------------------------------
// CHORD

pub fn load_chord(config: &OrangeConfig) -> Result<Option<ChordRecord>> {
    if !config.has_reference() || config.is_targeted() {
        return Ok(None);
    }
    let Some(dir) = dir_exists(&config.chord_dir) else { return Ok(None) };
    let c = tr::read_chord(&tr::chord_filename(dir, &config.tumor_id))?;
    Ok(Some(ChordRecord { brca1_value: c.brca1_value, brca2_value: c.brca2_value, hrd_value: c.hrd_value, hr_status: c.hr_status, hrd_type: c.hrd_type }))
}

// ---------------------------------------------------------------------------------------------
// CUPPA

/// CuppaDataFactory.getSvFeatureValue: the first FEAT_CONTRIB value of the feature, rounded
fn sv_feature_value(entries: &[CuppaPredictionEntry], feature: &str) -> Result<i32> {
    let e = entries
        .iter()
        .find(|e| e.data_type == "FEAT_CONTRIB" && e.feature_name == feature)
        .ok_or_else(|| anyhow!("java.lang.Exception: Input CuppaPredictions is empty"))?;
    Ok(hmf_common::java_math::java_round(e.feature_value) as i32)
}

/// CuppaPredictions.hasRnaPredictions: the main combined classifier
fn main_classifier(entries: &[CuppaPredictionEntry]) -> &'static str {
    let rna = entries.iter().any(|e| e.data_type == "PROB" && e.classifier_name == "RNA_COMBINED" && !e.data_value.is_nan());
    if rna { "COMBINED" } else { "DNA_COMBINED" }
}

/// CuppaDataFactory.create
pub fn cuppa_data(entries: &[CuppaPredictionEntry]) -> Result<CuppaData> {
    let mut probs: Vec<&CuppaPredictionEntry> = entries.iter().filter(|e| e.data_type == "PROB").collect();
    if probs.is_empty() {
        bail!("java.lang.IllegalStateException: `CuppaPredictions` must be instantiated with a non-empty list");
    }
    // sortByRank: by rank group, then rank (a stable sort)
    probs.sort_by_key(|e| (e.rank_group, e.rank));
    let probs_owned: Vec<CuppaPredictionEntry> = probs.iter().map(|e| (*e).clone()).collect();
    let main = main_classifier(&probs_owned);
    let mut cancer_types: Vec<&str> = Vec::new();
    for e in &probs {
        if !cancer_types.contains(&e.cancer_type.as_str()) {
            cancer_types.push(&e.cancer_type);
        }
    }
    let mut predictions = Vec::new();
    for ct in cancer_types {
        // a HashMap by classifier: the last value of a classifier wins
        let value = |name: &str| probs.iter().filter(|e| e.cancer_type == ct && e.classifier_name == name).last().map(|e| e.data_value);
        let likelihood = value(main).ok_or_else(|| anyhow!("java.lang.NullPointerException: no {main} probability for {ct}"))?;
        let p = CuppaPrediction {
            cancer_type: ct.to_string(),
            likelihood,
            snv_pairwise_classifier: value("SNV96"),
            genomic_position_classifier: value("GEN_POS"),
            feature_classifier: value("EVENT"),
            alt_sj_cohort_classifier: value("ALT_SJ"),
            expression_pairwise_classifier: value("GENE_EXP"),
        };
        if !p.likelihood.is_nan() {
            predictions.push(p);
        }
    }
    if predictions.is_empty() {
        bail!("java.lang.IllegalStateException: `CuppaPredictions` must contain a non-empty list of predictions");
    }
    let mode = match main_classifier(entries) {
        "COMBINED" => "WGTS",
        _ => "WGS",
    };
    // OR5: Java reads the count of simple duplications of 32-200 bases from the feature of deletions of
    // 20 kb to 1 Mb
    let dup_feature = if java("OR5") { "sv.SIMPLE_DEL_20KB_1MB" } else { "sv.SIMPLE_DUP_32B_200B" };
    Ok(CuppaData {
        mode: mode.to_string(),
        best_prediction: predictions[0].clone(),
        simple_dups_32_to_200b: sv_feature_value(entries, dup_feature)?,
        max_complex_size: sv_feature_value(entries, "sv.MAX_COMPLEX_SIZE")?,
        telomeric_sgls: sv_feature_value(entries, "sv.TELOMERIC_SGL")?,
        line_count: sv_feature_value(entries, "sv.LINE")?,
        predictions,
    })
}

pub fn load_cuppa(config: &OrangeConfig) -> Result<Option<CuppaData>> {
    if config.reference_id.is_none() || config.is_targeted() {
        return Ok(None);
    }
    let Some(dir) = dir_exists(&config.cuppa_dir) else { return Ok(None) };
    let entries = tr::read_cuppa_predictions(&tr::cuppa_vis_data_filename(dir, &config.tumor_id))?;
    Ok(Some(cuppa_data(&entries)?))
}

// ---------------------------------------------------------------------------------------------
// PEACH

/// OrangeAlgo.loadPeachData (UGT1A1 left out) and OrangeConversion; None gives an empty list in the record
pub fn load_peach(config: &OrangeConfig) -> Result<Vec<PeachGenotype>> {
    let Some(reference) = &config.reference_id else { return Ok(Vec::new()) };
    let Some(dir) = dir_exists(&config.peach_dir) else {
        hmf_common::debug!("skipping Peach loading as no input has been provided");
        return Ok(Vec::new());
    };
    let genotypes = tr::read_peach_genotypes(&tr::peach_genotype_filename(dir, reference))?;
    let mut out: Vec<PeachGenotype> = Vec::new();
    for g in genotypes.iter().filter(|g| g.gene != "UGT1A1") {
        let p = PeachGenotype {
            gene: g.gene.clone(),
            allele: g.allele.clone(),
            allele_count: g.allele_count,
            haplotype: g.haplotype()?,
            function: g.function.clone(),
            linked_drugs: g.linked_drugs.clone(),
            url_prescription_info: g.url_prescription_info.clone(),
            panel_version: g.panel_version.clone(),
            repo_version: g.repo_version.clone(),
        };
        // a Set in the data model (LinkedHashSet): equal genotypes once
        if !out.contains(&p) {
            out.push(p);
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// SIGS

pub fn load_sigs(config: &OrangeConfig) -> Result<Option<Vec<SignatureAllocation>>> {
    let Some(dir) = dir_exists(&config.sigs_dir) else {
        hmf_common::debug!("skipping Signature loading as no input has been provided");
        return Ok(None);
    };
    let allocations = tr::read_signature_allocations(&tr::signature_allocation_filename(dir, &config.tumor_id))?;
    let etiologies = tr::snv_signature_etiologies();
    let mut out = Vec::new();
    for a in &allocations {
        let etiology = match etiologies.get(&a.signature) {
            Some(e) => e.clone(),
            None if a.signature == "MISALLOC" => "-".to_string(),
            None => bail!("java.lang.IllegalArgumentException: Signature {} is not in signatures etiology tsv", a.signature),
        };
        out.push(SignatureAllocation { signature: a.signature.clone(), etiology, allocation: a.allocation, percent: a.percent });
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// primary tumor

/// OrangeAlgo.formConfiguredPrimaryTumorDoid. OR1: Java looks the DOIDs up in an empty set, so none is ever
/// resolved; corrected, they are looked up in the DOID JSON (`-doid_json`).
pub fn configured_primary_tumor(config: &OrangeConfig) -> Result<Vec<OrangeDoidNode>> {
    let mut nodes: Vec<OrangeDoidNode> = Vec::new();
    if !config.primary_tumor_doids.is_empty() {
        hmf_common::debug!("determining configured primary tumor");
        let known = if java("OR1") {
            Vec::new()
        } else {
            crate::doid::read_doid_nodes(config.doid_json_file.as_deref().unwrap_or_default())?
        };
        for doid in &config.primary_tumor_doids {
            match known.iter().find(|n| &n.doid == doid) {
                Some(n) => {
                    hmf_common::debug!(" Adding DOID {} ({:?}) as configured primary tumor", doid, n.doid_term);
                    if !nodes.contains(n) {
                        nodes.push(n.clone());
                    }
                }
                None => hmf_common::warn!("Could not resolve doid '{doid}'"),
            }
        }
    } else if let Some(location) = &config.primary_tumor_location {
        nodes.push(OrangeDoidNode { doid: String::new(), doid_term: Some(location.clone()) });
    }
    Ok(nodes)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(data_type: &str, clf: &str, feature: &str, value: f64, cancer: &str, data: f64, rank: i32, group: i32) -> CuppaPredictionEntry {
        CuppaPredictionEntry {
            sample_id: "S".into(),
            data_type: data_type.into(),
            classifier_group: "DNA".into(),
            classifier_name: clf.into(),
            feature_name: feature.into(),
            feature_value: value,
            cancer_type: cancer.into(),
            data_value: data,
            rank,
            rank_group: group,
        }
    }

    #[test]
    fn cuppa_predictions() {
        let mut e = vec![
            entry("PROB", "DNA_COMBINED", "NONE", f64::NAN, "B", 0.3, 2, 0),
            entry("PROB", "DNA_COMBINED", "NONE", f64::NAN, "A", 0.7, 1, 0),
            entry("PROB", "SNV96", "NONE", f64::NAN, "A", 0.6, 1, 1),
            entry("PROB", "DNA_COMBINED", "NONE", f64::NAN, "C", f64::NAN, 3, 0),
        ];
        for (f, v) in [("sv.SIMPLE_DEL_20KB_1MB", 5.0), ("sv.SIMPLE_DUP_32B_200B", 1.4), ("sv.MAX_COMPLEX_SIZE", 2.5), ("sv.TELOMERIC_SGL", 0.0), ("sv.LINE", 3.0)] {
            e.push(entry("FEAT_CONTRIB", "EVENT", f, v, "A", 0.1, 1, 2));
        }
        let d = crate::compat::with_java(false, || cuppa_data(&e).unwrap());
        assert_eq!(d.mode, "WGS");
        assert_eq!(d.predictions.len(), 2);
        assert_eq!(d.best_prediction.cancer_type, "A");
        assert_eq!(d.best_prediction.snv_pairwise_classifier, Some(0.6));
        assert_eq!(d.predictions[1].snv_pairwise_classifier, None);
        assert_eq!(d.simple_dups_32_to_200b, 1);
        assert_eq!(d.max_complex_size, 3);
        let j = crate::compat::with_java(true, || cuppa_data(&e).unwrap());
        assert_eq!(j.simple_dups_32_to_200b, 5);
    }
}
