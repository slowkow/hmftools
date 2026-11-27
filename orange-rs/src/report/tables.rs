//! The report's tables (orange/report/tables): FrontPageTables, SomaticVariantTable, GermlineVariantTable,
//! GainDeletionTable, DnaFusionTable, DisruptionTable, ViralPresenceTable, ChrArmCopyNumberTable,
//! SignatureAllocationTable, PharmacogeneticsTable, HLAAlleleTable, RnaStatisticsTable, ExpressionTable,
//! RnaFusionTable, NovelSpliceJunctionTable. Rows are sorted with Java's TimSort, as the Java comparators are
//! not all consistent.

use std::cmp::Ordering;

use anyhow::{bail, Result};
use hmf_common::java_hash::JavaHashSet;
use hmf_common::java_sort::java_sort;
use hmf_common::regions::chromosome_rank_java;
use hmf_plot::pdf::table::Table;

use super::common::*;
use super::DocCtx;
use crate::compat::java;
use crate::config::OrangeConfig;
use crate::model::*;

pub const DRIVER_LIKELIHOOD_LOW_THRESHOLD: f64 = 0.2;
pub const DRIVER_LIKELIHOOD_MEDIUM_THRESHOLD: f64 = 0.8;
const SIG_MIN_PERCENT: f64 = 0.08;
const SIG_MIN_ALLOCATION: f64 = 8.0;
const UNKNOWN_ALLELE_STRING: &str = "Unresolved Haplotype";

/// DriverInterpretation.interpret
pub fn interpret(likelihood: f64) -> &'static str {
    if likelihood >= DRIVER_LIKELIHOOD_MEDIUM_THRESHOLD {
        "HIGH"
    } else if likelihood >= DRIVER_LIKELIHOOD_LOW_THRESHOLD {
        "MEDIUM"
    } else {
        "LOW"
    }
}

/// OrangeConstants.isCandidateLikelihood
fn is_candidate(likelihood: f64) -> bool {
    likelihood <= DRIVER_LIKELIHOOD_LOW_THRESHOLD
}

/// Double.compare
fn double_compare(a: f64, b: f64) -> Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        _ => {
            if a < b {
                Ordering::Less
            } else if a > b {
                Ordering::Greater
            } else {
                // -0.0 before 0.0
                (a.is_sign_positive() as u8).cmp(&(b.is_sign_positive() as u8))
            }
        }
    }
}

/// PeachUtil.convertToZygosityString
fn zygosity(allele_count: i32) -> Result<&'static str> {
    match allele_count {
        1 => Ok("HET"),
        2 => Ok("HOM"),
        n => bail!("java.lang.IllegalArgumentException: Could not convert allele count {n} to a zygosity"),
    }
}

fn joined(items: impl IntoIterator<Item = String>) -> String {
    items.into_iter().collect::<Vec<_>>().join(", ")
}

/// A HashSet<String> of the items, in Java's iteration order
fn hash_set_order(items: impl IntoIterator<Item = String>) -> Vec<String> {
    let mut set: JavaHashSet<String> = JavaHashSet::new();
    for i in items {
        set.insert(i);
    }
    set.to_vec()
}

// ------------------------------------------------------------------------------------------------
// FrontPageTables

fn purity_string(fit: &PurpleFit) -> String {
    if fit.qc.status.iter().any(|s| s == "FAIL_NO_TUMOR") {
        return "NA".into();
    }
    format!("{} ({}-{})", percentage(fit.purity), percentage(fit.min_purity), percentage(fit.max_purity))
}

fn ploidy_string(fit: &PurpleFit) -> String {
    format!("{} ({}-{})", two_digit(fit.ploidy), two_digit(fit.min_ploidy), two_digit(fit.max_ploidy))
}

/// configuredPrimaryTumor: "<term> (DOID <doid>)" of each node, through a HashSet
fn configured_primary_tumor(nodes: &[OrangeDoidNode]) -> String {
    let items = nodes.iter().map(|n| {
        let mut s = n.doid_term.clone().unwrap_or_else(|| "null".into());
        if !n.doid.is_empty() {
            s += &format!(" (DOID {})", n.doid);
        }
        s
    });
    joined(hash_set_order(items))
}

pub fn is_fail_no_tumor(qc: &PurpleQC) -> bool {
    qc.status.iter().any(|s| s == "FAIL_NO_TUMOR")
}

pub fn has_tumor_contaminated(qc: &PurpleQC) -> bool {
    qc.status.iter().any(|s| s == "FAIL_CONTAMINATION" || s == "FAIL_TINC")
}

/// FrontPageTables.buildSampleSummary
pub fn sample_summary(ctx: &mut DocCtx, record: &OrangeRecord, config: &OrangeConfig, width: f32) -> Result<Table> {
    let widths = [2f32, 2.0, 2.0, 2.0, 2.0];
    let headers: Vec<String> = ["Primary Tumor", "Purity", "Ploidy", "Fit Method", "QC"].iter().map(|s| s.to_string()).collect();
    let cancer_type = match &config.primary_tumor_location {
        Some(l) => l.clone(),
        None => configured_primary_tumor(&record.configured_primary_tumor),
    };
    let mut table = create_content(ctx, width, &widths, &headers)?;
    let pcts = to_percentages(&widths);
    let fit = &record.purple.fit;
    let row = vec![
        if cancer_type.is_empty() { "NOT SPECIFIED".to_string() } else { cancer_type },
        purity_string(fit),
        ploidy_string(fit),
        fit.fitted_purity_method.clone(),
        joined(hash_set_order(fit.qc.status.iter().cloned())),
    ];
    let res = ctx.res;
    add_row(&res, &mut table, &pcts, &row)?;
    // addQCWarningInCaseOfFail
    let (no_tumor, contaminated) = (is_fail_no_tumor(&fit.qc), has_tumor_contaminated(&fit.qc));
    if no_tumor || contaminated {
        let reason = if no_tumor && contaminated {
            "no tumor and contamination"
        } else if no_tumor {
            "no tumor"
        } else {
            "contamination"
        };
        let warning = format!("The QC status of this sample is fail ({reason}): all presented data in this report should be interpreted with caution");
        table.create_row(12.0);
        let r = table.last_row_index();
        add_cell(&res, &mut table, r, 100.0, &warning, Style::Warning)?;
    }
    Ok(table)
}

/// FrontPageTables.buildTechnicalSummary
pub fn technical_summary(ctx: &mut DocCtx, record: &OrangeRecord, config: &OrangeConfig, width: f32) -> Result<Table> {
    let widths = [2f32; 6];
    let headers: Vec<String> =
        ["Pipeline Version", "Genome Version", "Sequencing Type", "Pipeline", "Samples", "Date Analysed"].iter().map(|s| s.to_string()).collect();
    let mut table = create_content(ctx, width, &widths, &headers)?;
    let pcts = to_percentages(&widths);
    let mode = if record.experiment_type == "WHOLE_GENOME" {
        "WHOLE GENOME".to_string()
    } else if let Some(p) = &config.panel_name {
        format!("PANEL: {p}")
    } else {
        "TARGETED PANEL".to_string()
    };
    let tumor_only = record.reference_id.is_none();
    let has_rna = record.isofox.is_some();
    let samples = if tumor_only && !has_rna {
        "TUMOR-ONLY".to_string()
    } else {
        let mut v = vec!["TUMOR"];
        if !tumor_only {
            v.push("NORMAL");
        }
        if has_rna {
            v.push("RNA");
        }
        v.join(" / ")
    };
    let row = vec![
        record.pipeline_version.clone().unwrap_or_default(),
        record.ref_genome_version.clone(),
        config.sequencing_type.clone(),
        mode,
        samples,
        record.sampling_date.clone(),
    ];
    let res = ctx.res;
    add_row(&res, &mut table, &pcts, &row)?;
    Ok(table)
}

fn key_value(ctx: &DocCtx, table: &mut Table, name: &str, value: &str) -> Result<()> {
    let res = ctx.res;
    table.create_row(COMPACT_ROW_HEIGHT);
    let r = table.last_row_index();
    add_cell(&res, table, r, 50.0, name, Style::Key)?;
    add_cell(&res, table, r, 50.0, value, Style::Value)?;
    Ok(())
}

fn titled_box_table(ctx: &mut DocCtx, width: f32, x: f32, title: &str) -> Result<Table> {
    let res = ctx.res;
    let mut table = ctx.create_table_at_x(width, x)?;
    table.create_row(15.0);
    let r = table.last_row_index();
    add_cell(&res, &mut table, r, 100.0, title, Style::Title)?;
    Ok(table)
}

/// reportableEventString: the count, and the high-likelihood genes sorted
fn reportable_event_string(count: usize, genes: Vec<String>) -> String {
    if genes.is_empty() {
        return count.to_string();
    }
    let mut g = genes;
    g.sort();
    g.dedup();
    format!("{count} ({})", g.join(", "))
}

fn driver_genes<'a>(items: impl Iterator<Item = (f64, &'a str)>) -> Vec<String> {
    items.filter(|(l, _)| interpret(*l) != "LOW").map(|(_, g)| g.to_string()).collect()
}

fn variants_string(variants: Option<&Vec<PurpleVariant>>) -> String {
    match variants {
        Some(v) if !v.is_empty() => reportable_event_string(v.len(), driver_genes(v.iter().map(|x| (x.driver_likelihood, x.gene.as_str())))),
        _ => "None".into(),
    }
}

fn gains_dels_string(items: Option<&Vec<PurpleGainDeletion>>) -> String {
    match items {
        Some(v) if !v.is_empty() => reportable_event_string(v.len(), driver_genes(v.iter().map(|x| (x.driver.driver_likelihood, x.driver.gene.as_str())))),
        _ => "None".into(),
    }
}

fn breakends_string(items: Option<&Vec<LinxBreakend>>) -> String {
    match items {
        Some(v) if !v.is_empty() => reportable_event_string(v.len(), driver_genes(v.iter().map(|x| (x.driver_likelihood, x.gene.as_str())))),
        _ => "None".into(),
    }
}

fn fusions_string(fusions: &[LinxFusion]) -> String {
    if fusions.is_empty() {
        return "None".into();
    }
    let genes = fusions.iter().filter(|f| f.driver_interpretation != "LOW").map(|f| format!("{}::{}", f.gene_up, f.gene_down)).collect();
    reportable_event_string(fusions.len(), genes)
}

fn viruses_string(viridian: &ViridianRecord) -> String {
    if viridian.reportable_viruses.is_empty() {
        return "None".into();
    }
    let names = viridian.reportable_viruses.iter().map(|v| v.interpretation.clone().unwrap_or_else(|| v.name.clone())).collect();
    reportable_event_string(viridian.reportable_viruses.len(), names)
}

fn peach_gene_status(peach: Option<&Vec<PeachGenotype>>, gene: &str) -> Result<String> {
    let Some(peach) = peach else { return Ok(NOT_AVAILABLE.into()) };
    let mut items = Vec::new();
    for g in peach.iter().filter(|g| g.gene == gene) {
        let haplotype = if g.allele == UNKNOWN_ALLELE_STRING { g.allele.clone() } else { format!("{} {}", g.allele, zygosity(g.allele_count)?) };
        items.push(format!("{haplotype} ({})", g.function));
    }
    let set = hash_set_order(items);
    Ok(if set.is_empty() { NOT_AVAILABLE.into() } else { joined(set) })
}

/// FrontPageTables.buildDriverSummary
pub fn driver_summary(ctx: &mut DocCtx, record: &OrangeRecord, width: f32, x: f32) -> Result<Table> {
    let include_germline = record.reference_id.is_some();
    let mut table = titled_box_table(ctx, width, x, "Driver Summary")?;
    let purple = &record.purple;
    let linx = &record.linx;
    key_value(ctx, &mut table, "Somatic variant:", &variants_string(Some(&purple.somatic_variants)))?;
    key_value(ctx, &mut table, "Somatic copy number:", &gains_dels_string(Some(&purple.somatic_gains_dels)))?;
    key_value(ctx, &mut table, "Somatic disruption:", &breakends_string(Some(&linx.somatic_breakends)))?;
    if include_germline {
        key_value(ctx, &mut table, "Germline variant:", &variants_string(purple.germline_variants.as_ref()))?;
        key_value(ctx, &mut table, "Germline copy number:", &gains_dels_string(purple.germline_gains_dels.as_ref()))?;
        key_value(ctx, &mut table, "Germline disruption:", &breakends_string(linx.germline_breakends.as_ref()))?;
    }
    key_value(ctx, &mut table, "Fusion drivers:", &fusions_string(&linx.fusions))?;
    if let Some(v) = &record.viridian {
        key_value(ctx, &mut table, "Viral presence:", &viruses_string(v))?;
    }
    let wgd = if purple.characteristics.whole_genome_duplication { "Yes" } else { "No" };
    key_value(ctx, &mut table, "Whole genome duplicated:", wgd)?;
    if include_germline {
        let s = peach_gene_status(record.peach.as_ref(), "DPYD")?;
        key_value(ctx, &mut table, "DPYD status:", &s)?;
    }
    Ok(table)
}

fn msi_display(status: &str) -> &str {
    match status {
        "MSI" => "Unstable",
        "MSS" => "Stable",
        "UNKNOWN" => "Unknown",
        s => s,
    }
}

fn tumor_mutational_display(status: &str) -> &str {
    match status {
        "HIGH" => "High",
        "LOW" => "Low",
        "UNKNOWN" => "Unknown",
        s => s,
    }
}

fn chord_display(status: &str) -> &str {
    match status {
        "CANNOT_BE_DETERMINED" => "Cannot be determined",
        "HR_PROFICIENT" => "Proficient",
        "HR_DEFICIENT" => "Deficient",
        "UNKNOWN" => "Unknown",
        s => s,
    }
}

/// hrDeficiencyString. OR10: Java fails on a run with a reference but without CHORD (targeted runs, or no
/// -chord_dir); corrected: "NA"
fn hr_deficiency_string(chord: Option<&ChordRecord>) -> Result<String> {
    let Some(chord) = chord else {
        if java("OR10") {
            bail!(
                "java.lang.NullPointerException: Cannot invoke \"com.hartwig.hmftools.datamodel.chord.ChordRecord.hrStatus()\" because \"chord\" is null (FrontPageTables.hrDeficiencyString)"
            );
        }
        return Ok(NOT_AVAILABLE.into());
    };
    if chord.hr_status == "CANNOT_BE_DETERMINED" {
        return Ok(chord_display(&chord.hr_status).into());
    }
    let mut addon = String::new();
    if chord.hr_status == "HR_DEFICIENT" && !chord.hrd_type.is_empty() {
        addon = if chord.hrd_type.contains("BRCA1") {
            format!(" - BRCA1 ({})", two_digit(chord.brca1_value))
        } else if chord.hrd_type.contains("BRCA2") {
            format!(" - BRCA2 ({})", two_digit(chord.brca2_value))
        } else if chord.hrd_type == "cannot_be_determined" {
            " - Undetermined".into()
        } else {
            format!(" - {}", chord.hrd_type)
        };
    }
    Ok(format!("{} ({}{addon})", single_digit(chord.hrd_value), chord_display(&chord.hr_status)))
}

/// FrontPageTables.buildGenomeWideFeaturesAtX
pub fn genome_wide_features(ctx: &mut DocCtx, record: &OrangeRecord, width: f32, x: f32) -> Result<Table> {
    let include_germline = record.reference_id.is_some();
    let mut table = titled_box_table(ctx, width, x, "Genome Wide Biomarkers")?;
    let c = &record.purple.characteristics;
    let msi = format!("{} ({})", single_digit(c.microsatellite_indels_per_mb), msi_display(&c.microsatellite_status));
    key_value(ctx, &mut table, "Microsatellite indels per Mb:", &msi)?;
    let tmb = format!("{} ({})", single_digit(c.tumor_mutational_burden_per_mb), tumor_mutational_display(&c.tumor_mutational_burden_status));
    key_value(ctx, &mut table, "Tumor mutations per Mb:", &tmb)?;
    let tml = format!("{} ({})", c.tumor_mutational_load, tumor_mutational_display(&c.tumor_mutational_load_status));
    key_value(ctx, &mut table, "Tumor mutational load:", &tml)?;
    if include_germline {
        let hrd = hr_deficiency_string(record.chord.as_ref())?;
        key_value(ctx, &mut table, "HR deficiency score:", &hrd)?;
    }
    key_value(ctx, &mut table, "LOH proportion:", &percentage(c.loh_percentage))?;
    key_value(ctx, &mut table, "Number of SVs:", &c.sv_tumor_mutational_burden.to_string())?;
    if record.reference_id.is_some() {
        if let Some(cuppa) = &record.cuppa {
            let best = &cuppa.best_prediction;
            key_value(ctx, &mut table, "CUPPA cancer type:", &format!("{} ({})", best.cancer_type, percentage(best.likelihood)))?;
        }
    }
    Ok(table)
}

// ------------------------------------------------------------------------------------------------
// Variants

/// SomaticVariantTable.sort: by driver likelihood (high first), then gene
fn sort_variants(variants: &[PurpleVariant]) -> Vec<PurpleVariant> {
    let mut v = variants.to_vec();
    java_sort(&mut v, |a, b| {
        let c = double_compare(b.driver_likelihood, a.driver_likelihood);
        if c != Ordering::Equal {
            return c;
        }
        a.gene.cmp(&b.gene)
    });
    v
}

fn location_display(v: &PurpleVariant) -> String {
    let phased = v.local_phase_sets.as_ref().is_some_and(|p| !p.is_empty());
    let loc = format!("{}:{}", v.chromosome, v.position);
    if phased { format!("{loc} *") } else { loc }
}

fn hgvs_display(i: &PurpleTranscriptImpact) -> String {
    if i.hgvs_protein_impact.is_empty() {
        i.hgvs_coding_impact.clone()
    } else {
        format!("{} [{}]", i.hgvs_coding_impact, i.hgvs_protein_impact)
    }
}

fn copy_number_display(v: &PurpleVariant) -> String {
    format!("{} of {}", fixed(v.variant_copy_number, 1), fixed(v.adjusted_copy_number, 1))
}

fn hotspot_display(v: &PurpleVariant) -> &'static str {
    match v.hotspot.as_str() {
        "HOTSPOT" => "Yes",
        "NEAR_HOTSPOT" => "Near",
        _ => "No",
    }
}

fn allele_frequency(d: &AllelicDepth) -> f64 {
    if d.total_read_count > 0 { d.allele_read_count as f64 / d.total_read_count as f64 } else { 0.0 }
}

fn rna_display(v: &PurpleVariant) -> String {
    match &v.rna_depth {
        None => NOT_AVAILABLE.into(),
        Some(d) => format!("{}/{}", d.allele_read_count, d.total_read_count),
    }
}

/// SomaticVariantTable.build
pub fn somatic_variants(ctx: &mut DocCtx, title: &str, width: f32, variants: &[PurpleVariant], tumor_only: bool, has_rna: bool) -> Result<Table> {
    if variants.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut headers = vec![
        HeaderData::new("Gene", 0.0, false),
        HeaderData::new("Position", 0.0, false),
        HeaderData::new("HGVS", 0.0, true),
        HeaderData::new("AF", 0.0, false),
        HeaderData::new("Depth", 0.0, false),
        HeaderData::new("Copies", 0.0, false),
        HeaderData::new("Hotspot", 0.0, false),
        HeaderData::new("Biall", 0.0, false),
        HeaderData::new("Clonal", 0.0, false),
    ];
    if tumor_only {
        headers.push(HeaderData::new("Somatic", 0.0, false));
    }
    if has_rna {
        headers.push(HeaderData::new("RNA", 0.0, false));
    }
    headers.push(HeaderData::new("Driver", 0.0, false));
    let mut rows = Vec::new();
    let mut shade = Vec::new();
    for v in sort_variants(variants) {
        let impacts = std::iter::once(&v.canonical_impact).chain(v.other_impacts.iter());
        for imp in impacts {
            let mut row = vec![
                v.gene.clone(),
                location_display(&v),
                hgvs_display(imp),
                two_digit(allele_frequency(&v.tumor_depth)),
                v.tumor_depth.total_read_count.to_string(),
                copy_number_display(&v),
                hotspot_display(&v).to_string(),
            ];
            // formatPercentageField(biallelicProbability()): Java unboxes the Double (NullPointerException when null)
            match v.biallelic_probability {
                Some(p) => row.push(percentage_field(p)),
                None => bail!("java.lang.NullPointerException: Cannot invoke \"java.lang.Double.doubleValue()\" because the return value of \"PurpleVariant.biallelicProbability()\" is null (SomaticVariantTable.build)"),
            }
            row.push(percentage_field(1.0 - v.subclonal_likelihood));
            if tumor_only {
                match &v.somatic_likelihood {
                    Some(s) => row.push(s.clone()),
                    None => bail!("java.lang.NullPointerException: Cannot invoke \"Object.toString()\" because the return value of \"PurpleVariant.somaticLikelihood()\" is null (SomaticVariantTable.build)"),
                }
            }
            if has_rna {
                row.push(rna_display(&v));
            }
            row.push(percentage_field(v.driver_likelihood));
            shade.push(is_candidate(v.driver_likelihood));
            rows.push(row);
        }
    }
    scaled_table(ctx, title, width, headers, rows, &shade)
}

fn genotype_display(status: &str) -> &str {
    match status {
        "HOM_REF" | "HOM_ALT" => "HOM",
        "HET" => "HET",
        s => s,
    }
}

/// GermlineVariantTable.build. OR13: a variant without a gnomAD frequency shows "  nu" (String.format("%4.2e")
/// of null); corrected: "NA"
pub fn germline_variants(ctx: &mut DocCtx, title: &str, width: f32, variants: &[PurpleVariant]) -> Result<Table> {
    if variants.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let headers = vec![
        HeaderData::new("Gene", 0.0, false),
        HeaderData::new("Variant", 0.0, false),
        HeaderData::new("HGVS", 0.0, true),
        HeaderData::new("Zygosity", 0.0, false),
        HeaderData::new("AF", 0.0, false),
        HeaderData::new("Depth", 0.0, false),
        HeaderData::new("Copies", 0.0, false),
        HeaderData::new("Hotspot", 0.0, false),
        HeaderData::new("Biall", 0.0, false),
        HeaderData::new("Gnomad", 0.0, false),
        HeaderData::new("Clinvar", 0.0, false),
        HeaderData::new("Driver", 0.0, false),
    ];
    let mut rows = Vec::new();
    let mut shade = Vec::new();
    for v in sort_variants(variants) {
        let impacts = std::iter::once(&v.canonical_impact).chain(v.other_impacts.iter());
        for imp in impacts {
            let gnomad = if v.gnomad_frequency.is_none() && !java("OR13") { NOT_AVAILABLE.to_string() } else { sci_4_2(v.gnomad_frequency) };
            rows.push(vec![
                v.gene.clone(),
                location_display(&v),
                hgvs_display(imp),
                genotype_display(&v.genotype_status).to_string(),
                two_digit(allele_frequency(&v.tumor_depth)),
                v.tumor_depth.total_read_count.to_string(),
                copy_number_display(&v),
                hotspot_display(&v).to_string(),
                if v.biallelic { "Yes".into() } else { "No".into() },
                gnomad,
                v.clinvar_pathogenicity.clone().unwrap_or_default(),
                interpret(v.driver_likelihood).to_string(),
            ]);
            shade.push(is_candidate(v.driver_likelihood));
        }
    }
    scaled_table(ctx, title, width, headers, rows, &shade)
}

// ------------------------------------------------------------------------------------------------
// Copy number

/// GainDeletionTable.build
pub fn gains_deletions(ctx: &mut DocCtx, title: &str, width: f32, items: &[PurpleGainDeletion], has_rna: bool) -> Result<Table> {
    if items.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut headers = vec![
        HeaderData::new("Location", 0.0, false),
        HeaderData::new("Gene", 0.0, false),
        HeaderData::new("Type", 0.0, false),
        HeaderData::new("Range", 0.0, false),
        HeaderData::new("Min CN", 0.0, false),
        HeaderData::new("Max CN", 0.0, false),
        HeaderData::new("Rel CN", 0.0, false),
        HeaderData::new("Arm CN", 0.0, false),
    ];
    let any_rna = has_rna && items.iter().any(|x| x.tpm.is_some());
    if any_rna {
        headers.push(HeaderData::new("TPM", 0.0, false));
        headers.push(HeaderData::new("Percentile", 0.0, false));
        headers.push(HeaderData::new("Fold Change", 0.6, false));
    }
    headers.push(HeaderData::new("Driver", 0.0, false));
    let mut sorted = items.to_vec();
    java_sort(&mut sorted, |a, b| {
        let c = double_compare(-a.driver.driver_likelihood, -b.driver.driver_likelihood);
        if c != Ordering::Equal {
            return c;
        }
        let l1 = zero_prefixed(&format!("{}{}", a.chromosome, a.chromosome_band));
        let l2 = zero_prefixed(&format!("{}{}", b.chromosome, b.chromosome_band));
        if l1 == l2 { a.driver.gene.cmp(&b.driver.gene) } else { l1.cmp(&l2) }
    });
    let mut rows = Vec::new();
    let mut shade = Vec::new();
    for g in &sorted {
        let gene = if g.driver.is_canonical { g.driver.gene.clone() } else { format!("{} (alt)", g.driver.gene) };
        let kind = match g.driver.driver_type.as_str() {
            "GERMLINE_AMP" => "AMP".to_string(),
            "GERMLINE_DELETION" => match &g.germline_amp_del_fields {
                Some(f) => (if f.germline_status == "HOM_DELETION" { "HOM" } else { "HET" }).to_string(),
                None => bail!("java.lang.NullPointerException: germlineAmpDelFields() is null (GainDeletionTable.typeDisplay)"),
            },
            t => t.to_string(),
        };
        let range = if g.exon_start.is_none() && g.exon_end.is_none() {
            g.gene_range.clone()
        } else {
            let f = |e: Option<i32>| e.map(|x| x.to_string()).unwrap_or_else(|| "null".into());
            format!("Exon {} - Exon {}", f(g.exon_start), f(g.exon_end))
        };
        let mut row = vec![
            format!("{}{}", g.chromosome, g.chromosome_band),
            gene,
            kind,
            range,
            single_digit(g.min_copy_number),
            single_digit(g.max_copy_number),
            single_digit(g.relative_copy_number),
            single_digit(g.arm_copy_number),
        ];
        if any_rna {
            match g.tpm {
                Some(tpm) => {
                    row.push(tpm_field(tpm));
                    // formatPercentage(tpmPercentile()) and formatFoldChangeField(tpmFoldChange()) unbox the Doubles
                    let (Some(p), Some(f)) = (g.tpm_percentile, g.tpm_fold_change) else {
                        bail!("java.lang.NullPointerException: a TPM percentile or fold change is null (GainDeletionTable.build)");
                    };
                    row.push(percentage(p));
                    row.push(fold_change_field(f));
                }
                None => {
                    row.push(NOT_AVAILABLE.into());
                    row.push(NOT_AVAILABLE.into());
                    row.push(NOT_AVAILABLE.into());
                }
            }
        }
        row.push(g.driver.driver_interpretation.clone());
        shade.push(g.driver.driver_interpretation == "LOW");
        rows.push(row);
    }
    scaled_table(ctx, title, width, headers, rows, &shade)
}

/// ChrArmCopyNumberTable.build
pub fn chr_arm_copy_numbers(ctx: &mut DocCtx, title: &str, width: f32, arms: &[PurpleChrArmCopyNumber]) -> Result<Table> {
    if arms.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut sorted = arms.to_vec();
    java_sort(&mut sorted, |a, b| {
        if a.driver_interpretation != b.driver_interpretation {
            return if a.driver_interpretation == "HIGH" { Ordering::Less } else { Ordering::Greater };
        }
        let (c1, c2) = (chromosome_rank_java(&a.chromosome), chromosome_rank_java(&b.chromosome));
        if c1 != c2 {
            return c1.cmp(&c2);
        }
        a.arm.cmp(&b.arm)
    });
    let rows: Vec<Vec<String>> = sorted
        .iter()
        .map(|a| {
            vec![
                a.chromosome.clone(),
                a.arm.clone(),
                a.cn_type.clone(),
                single_digit(a.copy_number),
                single_digit(a.relative_copy_number),
                a.driver_interpretation.clone(),
                String::new(),
            ]
        })
        .collect();
    let shade: Vec<bool> = sorted.iter().map(|a| a.driver_interpretation == "LOW").collect();
    fixed_table(ctx, title, width, &[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0], &["Chromosome", "Arm", "Type", "CN", "Rel CN", "Driver", ""], rows, &shade)
}

// ------------------------------------------------------------------------------------------------
// Structural variants

/// DnaFusionTable.build. OR11: once one fusion has RNA support Java reads it for every fusion and fails on one
/// without; corrected: "NA" for those
pub fn dna_fusions(ctx: &mut DocCtx, title: &str, width: f32, fusions: &[LinxFusion]) -> Result<Table> {
    if fusions.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut headers = vec![
        HeaderData::new("Fusion", 0.0, false),
        HeaderData::new("Junctions", 0.0, false),
        HeaderData::new("JCN", 0.0, false),
        HeaderData::new("Phasing", 0.0, false),
        HeaderData::new("Type", 0.0, false),
    ];
    let has_rna = fusions.iter().any(|f| f.rna_support.is_some());
    if has_rna {
        headers.push(HeaderData::new("RNA Frags", 0.0, false));
    }
    headers.push(HeaderData::new("Driver", 0.0, false));
    let mut sorted = fusions.to_vec();
    java_sort(&mut sorted, |a, b| {
        if a.driver_interpretation == b.driver_interpretation {
            if a.gene_up == b.gene_up { a.gene_down.cmp(&b.gene_down) } else { a.gene_up.cmp(&b.gene_up) }
        } else if a.driver_interpretation == "HIGH" {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    });
    let mut rows = Vec::new();
    for f in &sorted {
        let phasing = match f.phased.as_str() {
            "INFRAME" => "Inframe",
            "SKIPPED_EXONS" => "Skipped exons",
            "OUT_OF_FRAME" => "Out of frame",
            p => p,
        };
        let mut row = vec![
            format!("{}::{}", f.gene_up, f.gene_down),
            format!("{} ({}) - {} ({})", f.context_up, f.transcript_up, f.context_down, f.transcript_down),
            single_digit(f.junction_copy_number),
            phasing.to_string(),
            f.reported_type.clone(),
        ];
        if has_rna {
            match &f.rna_support {
                Some(d) => row.push(d.allele_read_count.to_string()),
                None if java("OR11") => bail!(
                    "java.lang.NullPointerException: Cannot invoke \"AllelicDepth.alleleReadCount()\" because the return value of \"LinxFusion.rnaSupport()\" is null (DnaFusionTable.build)"
                ),
                None => row.push(NOT_AVAILABLE.into()),
            }
        }
        row.push(f.driver_interpretation.clone());
        rows.push(row);
    }
    scaled_table(ctx, title, width, headers, rows, &[])
}

fn breakend_location(b: &LinxBreakend) -> String {
    format!("{}:{}", b.chromosome, b.position)
}

fn context_str(b: &LinxBreakend, include_orientation: bool) -> String {
    let mut exon_range: Option<String> = None;
    if b.exon_up > 0 {
        if b.exon_up == b.exon_down {
            exon_range = Some(format!("Exon {}", b.exon_up));
        } else if b.exon_down - b.exon_up == 1 {
            exon_range = Some(format!("Intron {}", b.exon_up));
        }
    } else if b.exon_up == 0 && (b.exon_down == 1 || b.exon_down == 2) {
        exon_range = Some("Promoter Region".into());
    }
    let er = exon_range.unwrap_or_else(|| "null".into());
    if include_orientation {
        let o = match b.gene_orientation.as_str() {
            "UPSTREAM" => "Upstream",
            "DOWNSTREAM" => "Downstream",
            o => o,
        };
        format!("{er} {o}")
    } else {
        er
    }
}

/// DisruptionTable.build: breakends of one SV on consecutive rows are shown as one
pub fn disruptions(ctx: &mut DocCtx, title: &str, width: f32, breakends: &[LinxBreakend]) -> Result<Table> {
    if breakends.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let headers = vec![
        HeaderData::new("Gene", 0.0, false),
        HeaderData::new("Position", 0.0, false),
        HeaderData::new("Zygosity", 0.0, false),
        HeaderData::new("Context", 0.0, false),
        HeaderData::new("Type", 0.0, false),
        HeaderData::new("JCN", 0.0, false),
        HeaderData::new("Undisrupted CN", 0.0, false),
        HeaderData::new("Driver", 0.0, false),
    ];
    let mut sorted = breakends.to_vec();
    java_sort(&mut sorted, |a, b| {
        if a.sv_id == b.sv_id {
            return a.id.cmp(&b.id);
        }
        let l1 = format!("{}{}", zero_prefixed(&a.chromosome), a.chromosome_band);
        let l2 = format!("{}{}", zero_prefixed(&b.chromosome), b.chromosome_band);
        let c = l1.cmp(&l2);
        if c != Ordering::Equal {
            return c;
        }
        let g = a.gene.cmp(&b.gene);
        if g != Ordering::Equal {
            return g;
        }
        // breakend1.exonUp() - breakend2.exonUp()
        (a.exon_up.wrapping_sub(b.exon_up)).cmp(&0)
    });
    let mut rows = Vec::new();
    let mut shade = Vec::new();
    let mut i = 0;
    while i < sorted.len() {
        let b = &sorted[i];
        let other = if i + 1 < sorted.len() && sorted[i + 1].sv_id == b.sv_id { Some(&sorted[i + 1]) } else { None };
        let location = match other {
            None => breakend_location(b),
            Some(o) => format!("{} - {}", breakend_location(b), breakend_location(o)),
        };
        let context = match other {
            None => context_str(b, true),
            Some(o) => format!("{} - {}", context_str(b, false), context_str(o, false)),
        };
        rows.push(vec![
            b.gene.clone(),
            location,
            (if b.undisrupted_copy_number < 0.5 { "HOM" } else { "HET" }).to_string(),
            context,
            b.breakend_type.clone(),
            single_digit(b.junction_copy_number),
            single_digit(b.undisrupted_copy_number),
            interpret(b.driver_likelihood).to_string(),
        ]);
        shade.push(is_candidate(b.driver_likelihood));
        if other.is_some() {
            i += 1;
        }
        i += 1;
    }
    scaled_table(ctx, title, width, headers, rows, &shade)
}

// ------------------------------------------------------------------------------------------------
// Viruses, signatures, pharmacogenetics, HLA

/// ViralPresenceTable.build
pub fn viral_presence(ctx: &mut DocCtx, title: &str, width: f32, viruses: &[VirusReport]) -> Result<Table> {
    if viruses.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let rows: Vec<Vec<String>> = viruses
        .iter()
        .map(|v| {
            vec![
                v.name.clone(),
                v.qc_status.clone(),
                v.interpretation.clone().unwrap_or_default(),
                v.integrations.to_string(),
                percentage(v.percentage_covered),
                single_digit(v.mean_coverage),
                v.copies_per_tumor_cell.map(single_digit).unwrap_or_default(),
            ]
        })
        .collect();
    fixed_table(
        ctx,
        title,
        width,
        &[2.5, 2.5, 1.0, 2.0, 1.5, 1.5, 2.0],
        &["Virus", "QC Status", "Type", "Integrations", "% Covered", "Mean Cov", "Copies/Cell"],
        rows,
        &[],
    )
}

/// SignatureAllocationTable.sort: by allocation, highest first, MISALLOC last
pub fn sort_signature_allocations(allocations: &[SignatureAllocation]) -> Vec<SignatureAllocation> {
    let mut sorted = allocations.to_vec();
    java_sort(&mut sorted, |a, b| {
        let (s1, s2) = (a.signature.to_lowercase(), b.signature.to_lowercase());
        if s1.eq_ignore_ascii_case("MISALLOC") {
            return Ordering::Greater;
        } else if s2.eq_ignore_ascii_case("MISALLOC") {
            return Ordering::Less;
        }
        if !s1.starts_with("sig") || !s2.starts_with("sig") {
            hmf_common::warn!("Signatures do not start with 'sig': {} & {}", s1, s2);
            return s1.cmp(&s2);
        }
        double_compare(b.allocation, a.allocation)
    });
    sorted
}

/// SignatureAllocationTable.build: the allocations of at least 8% and 8 variants
pub fn signature_allocations(ctx: &mut DocCtx, title: &str, width: f32, allocations: &[SignatureAllocation]) -> Result<Table> {
    if allocations.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let sorted = sort_signature_allocations(allocations);
    let rows: Vec<Vec<String>> = sorted
        .iter()
        .filter(|s| !(s.percent < SIG_MIN_PERCENT || s.allocation < SIG_MIN_ALLOCATION))
        .map(|s| vec![s.signature.clone(), s.etiology.clone(), fixed(s.allocation, 0), percentage(s.percent), String::new()])
        .collect();
    fixed_table(ctx, title, width, &[1.0, 2.0, 1.0, 1.0, 4.0], &["Signature", "Etiology", "Allocation", "Percent", ""], rows, &[])
}

/// PharmacogeneticsTable.build
pub fn pharmacogenetics(ctx: &mut DocCtx, title: &str, width: f32, genotypes: &[PeachGenotype]) -> Result<Table> {
    if genotypes.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut sorted = genotypes.to_vec();
    java_sort(&mut sorted, |a, b| {
        if a.gene != b.gene {
            a.gene.cmp(&b.gene)
        } else if a.allele != b.allele {
            a.allele.cmp(&b.allele)
        } else {
            a.allele_count.cmp(&b.allele_count)
        }
    });
    let mut rows = Vec::new();
    for g in &sorted {
        let genotype = if g.allele == UNKNOWN_ALLELE_STRING { String::new() } else { zygosity(g.allele_count)?.to_string() };
        let url = g.url_prescription_info.split(';').next().unwrap_or("");
        let source = if url.starts_with("https://www.pharmgkb.org") { "PHARMGKB" } else { "" };
        rows.push(vec![g.gene.clone(), g.allele.clone(), genotype, g.function.clone(), g.linked_drugs.clone(), source.to_string()]);
    }
    fixed_table(ctx, title, width, &[1.0, 1.0, 1.0, 1.0, 2.0, 1.0], &["Gene", "Haplotype", "Genotype", "Function", "Linked drugs", "Source"], rows, &[])
}

fn fragment_string(f: Option<i32>) -> String {
    f.map(|x| x.to_string()).unwrap_or_else(|| NOT_AVAILABLE.into())
}

/// Doubles.positive
fn positive(v: f64) -> bool {
    v - 0.0 > 1e-10
}

fn mutation_string(a: &LilacAllele) -> String {
    let mut parts = Vec::new();
    for (v, what) in [
        (a.somatic_missense, "missense"),
        (a.somatic_nonsense_or_frameshift, "nonsense or frameshift"),
        (a.somatic_splice, "splice"),
        (a.somatic_synonymous, "synonymous"),
        (a.somatic_inframe_indel, "inframe indel"),
    ] {
        if positive(v) {
            parts.push(format!("{} {what}", single_digit(v)));
        }
    }
    if parts.is_empty() { NONE.into() } else { parts.join(", ") }
}

/// HLAAlleleTable.build
pub fn hla_alleles(ctx: &mut DocCtx, title: &str, width: f32, alleles: &[LilacAllele], has_rna: bool) -> Result<Table> {
    if alleles.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let has_warnings = alleles.iter().any(|a| a.qc_status != "PASS");
    let mut widths = vec![1.0f32, if has_warnings { 3.0 } else { 1.0 }, 1.0, 1.0];
    let mut headers = vec!["Allele", "QC Status", "Ref Frags", "Tumor Frags"];
    if has_rna {
        widths.push(1.0);
        headers.push("RNA Frags");
    }
    widths.push(1.0);
    headers.push("Tumor CN");
    widths.push(3.0);
    headers.push("Somatic Mutations");
    let mut sorted = alleles.to_vec();
    java_sort(&mut sorted, |a, b| a.allele.cmp(&b.allele).then(a.ref_fragments.unwrap_or(0).cmp(&b.ref_fragments.unwrap_or(0))));
    let rows: Vec<Vec<String>> = sorted
        .iter()
        .map(|a| {
            let mut row = vec![a.allele.clone(), a.qc_status.clone(), fragment_string(a.ref_fragments), fragment_string(Some(a.tumor_fragments))];
            if has_rna {
                row.push(fragment_string(a.rna_fragments));
            }
            row.push(single_digit(a.tumor_copy_number));
            row.push(mutation_string(a));
            row
        })
        .collect();
    fixed_table(ctx, title, width, &widths, &headers, rows, &[])
}

// ------------------------------------------------------------------------------------------------
// RNA

/// RnaStatisticsTable.build
pub fn rna_statistics(ctx: &mut DocCtx, title: &str, width: f32, s: &RnaStatistics) -> Result<Table> {
    let duplicate_rate = s.duplicate_fragments as f64 / s.total_fragments as f64;
    let row = vec![
        s.qc_status.join(", "),
        s.total_fragments.to_string(),
        percentage(duplicate_rate),
        percentage(s.spliced_fragment_perc),
        percentage(s.unspliced_fragment_perc),
        percentage(s.alt_fragment_perc),
        percentage(s.chimeric_fragment_perc),
    ];
    fixed_table(
        ctx,
        title,
        width,
        &[1.0; 7],
        &["QC", "Total Fragments", "Duplicate Rate", "Spliced Rate", "Unspliced Rate", "Alt-sliced Rate", "Chimeric Rate"],
        vec![row],
        &[],
    )
}

/// Expressions.toFoldChangeDisplay
pub fn fold_change_display(tpm: f64, median: Option<f64>) -> String {
    match median {
        // Doubles.lessOrEqual(median, 0): not when NaN
        Some(m) if !(m - 0.0 < 1e-10) => fold_change_field(tpm / m),
        _ => NOT_AVAILABLE.into(),
    }
}

/// ExpressionTable.build (high expression: sorted by cohort percentile, highest first)
pub fn expressions(ctx: &mut DocCtx, title: &str, width: f32, expressions: &[GeneExpression], sort_ascending: bool) -> Result<Table> {
    if expressions.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut sorted = expressions.to_vec();
    java_sort(&mut sorted, |a, b| {
        if sort_ascending { double_compare(a.percentile_cohort, b.percentile_cohort) } else { double_compare(b.percentile_cohort, a.percentile_cohort) }
    });
    let rows: Vec<Vec<String>> = sorted
        .iter()
        .map(|e| {
            let (pct, fold) = match (e.percentile_cancer, e.median_tpm_cancer) {
                (Some(p), Some(_)) => (percentage(p), fold_change_display(e.tpm, e.median_tpm_cancer)),
                _ => (percentage(e.percentile_cohort), fold_change_display(e.tpm, Some(e.median_tpm_cohort))),
            };
            vec![e.gene.clone(), tpm_field(e.tpm), pct, fold, String::new()]
        })
        .collect();
    fixed_table(ctx, title, width, &[1.0, 1.0, 1.0, 1.0, 3.0], &["Gene", "TPM", "Percentile", "Fold Change", ""], rows, &[])
}

/// (int) Math.round(double)
fn round_half(v: f64) -> i32 {
    hmf_common::java_math::java_round(v) as i32
}

/// RnaFusionTable.build
pub fn rna_fusions(ctx: &mut DocCtx, title: &str, width: f32, fusions: &[RnaFusion]) -> Result<Table> {
    if fusions.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut sorted = fusions.to_vec();
    java_sort(&mut sorted, |a, b| {
        let (l1, l2) = (zero_prefixed(&a.chromosome_start), zero_prefixed(&b.chromosome_start));
        if l1 == l2 { a.position_start.cmp(&b.position_start) } else { l1.cmp(&l2) }
    });
    let rows: Vec<Vec<String>> = sorted
        .iter()
        .map(|f| {
            let avg = round_half(f.depth_start.wrapping_add(f.depth_end) as f64 * 0.5);
            vec![
                format!("{}::{}", f.gene_start.clone().unwrap_or_default(), f.gene_end.clone().unwrap_or_default()),
                f.known_type.clone(),
                f.sv_type.clone(),
                f.junction_type_start.clone(),
                f.junction_type_end.clone(),
                support_field(f.split_fragments, avg),
                f.cohort_frequency.to_string(),
            ]
        })
        .collect();
    fixed_table(
        ctx,
        title,
        width,
        &[2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        &["Fusion", "Type", "SV Type", "Junc Start", "Junc End", "Support", "Cohort Freq"],
        rows,
        &[],
    )
}

/// NovelSpliceJunctionTable.build
pub fn novel_splice_junctions(ctx: &mut DocCtx, title: &str, width: f32, junctions: &[NovelSpliceJunction]) -> Result<Table> {
    if junctions.is_empty() {
        return create_empty_table(ctx, title, width);
    }
    let mut sorted = junctions.to_vec();
    java_sort(&mut sorted, |a, b| {
        let (l1, l2) = (zero_prefixed(&a.chromosome), zero_prefixed(&b.chromosome));
        if l1 == l2 { a.junction_start.cmp(&b.junction_start) } else { l1.cmp(&l2) }
    });
    let rows: Vec<Vec<String>> = sorted
        .iter()
        .map(|j| {
            let dup = j.junction_type == "CIRCULAR";
            let (ps, pe) = if dup { (j.junction_end, j.junction_start) } else { (j.junction_start, j.junction_end) };
            let (es, ee) = if dup { (j.exon_end, j.exon_start) } else { (j.exon_start, j.exon_end) };
            let avg = round_half(j.depth_start.wrapping_add(j.depth_end) as f64 * 0.5);
            vec![
                j.gene.clone(),
                format!("Exon {es} ({}:{ps}) - Exon {ee} ({}:{pe})", j.chromosome, j.chromosome),
                j.junction_type.clone(),
                j.region_start.clone(),
                j.region_end.clone(),
                support_field(j.fragment_count, avg),
                j.cohort_frequency.to_string(),
            ]
        })
        .collect();
    fixed_table(
        ctx,
        title,
        width,
        &[1.0, 3.0, 2.0, 1.0, 1.0, 1.0, 1.0],
        &["Gene", "Junctions", "Type", "Junc Start", "Junc End", "Support", "Cohort Freq"],
        rows,
        &[],
    )
}
