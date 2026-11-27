//! PURPLE's part of the report: PurpleDataLoader, PurpleInterpreter, PurpleVariantFactory and
//! GermlineGainDeletionFactory (orange/algo/purple), with PurpleConversion.

use anyhow::{anyhow, bail, Result};
use hmf_common::java_hash::JavaHashSet;
use hmf_common::java_math::java_max;
use hmf_common::purple_files::{
    self, chr_arm_copy_number_filename, driver_catalog_germline_filename, driver_catalog_somatic_filename, gene_copy_number_filename,
    germline_amp_del_filename, purple_purity_filename, purple_qc_filename, ChrArmCopyNumber, DriverCatalog, GeneCopyNumber,
    GermlineAmpDel, PurplePurity, PurpleQc,
};
use hmf_common::small_variant::{self, effects_to_list, AltTranscriptReportableInfo, SmallVariant, VariantTranscriptImpact};

use crate::compat::java;
use crate::config::OrangeConfig;
use crate::isofox::IsofoxData;
use crate::linx::{find_reportable_linx_plot, LinxData};
use crate::model::*;

pub struct PurpleData {
    pub qc: PurpleQc,
    pub purity: PurplePurity,
    pub somatic_drivers: Vec<DriverCatalog>,
    pub germline_drivers: Option<Vec<DriverCatalog>>,
    pub somatic_variants: Vec<SmallVariant>,
    /// the PAVE_TI transcript impacts of each somatic variant (for OR7's correction)
    pub somatic_transcript_impacts: Vec<Vec<VariantTranscriptImpact>>,
    pub germline_variants: Option<Vec<SmallVariant>>,
    pub germline_transcript_impacts: Vec<Vec<VariantTranscriptImpact>>,
    pub somatic_gene_copy_numbers: Vec<GeneCopyNumber>,
    pub germline_amp_dels: Option<Vec<GermlineAmpDel>>,
    pub chr_arm_copy_numbers: Vec<ChrArmCopyNumber>,
    pub variant_plots: Vec<String>,
}

/// PurpleDataLoader.resolveVcfPath: an absent `.gz` VCF falls back to the uncompressed name when that exists
fn resolve_vcf_path(path: &str) -> String {
    if !std::path::Path::new(path).exists() && path.ends_with(".gz") {
        let unzipped = &path[..path.len() - 3];
        if std::path::Path::new(unzipped).exists() {
            return unzipped.to_string();
        }
    }
    path.to_string()
}

fn dir_sep(d: &str) -> String {
    if d.is_empty() || d.ends_with('/') { d.to_string() } else { format!("{d}/") }
}

/// The variants and, for the reported ones, their PAVE_TI impacts
fn read_reported_variants(
    path: &str,
    tumor: &str,
    reference: Option<&str>,
    rna: Option<&str>,
    pass_only: bool,
) -> Result<(Vec<SmallVariant>, Vec<Vec<VariantTranscriptImpact>>)> {
    let mut reader = hmf_common::variant_context::VcfReader::open(path)?;
    let header = reader.header.clone();
    let has = |s: &str| header.samples.iter().any(|x| x == s);
    for s in [Some(tumor), reference, rna].into_iter().flatten() {
        if !has(s) {
            bail!("java.lang.IllegalArgumentException: Sample {s} not found in vcf file {path}");
        }
    }
    if !header.has_format("AD") {
        bail!("java.lang.IllegalArgumentException: Allelic depths is a required format field in vcf file {path}");
    }
    let mut variants = Vec::new();
    let mut impacts = Vec::new();
    loop {
        let Some(line) = reader.next_line()? else { break };
        // a filtered record is skipped unparsed where only PASS records are read (the FILTER column: "PASS" or ".")
        if pass_only {
            let filter = line.splitn(8, '\t').nth(6).unwrap_or("");
            if filter != "PASS" && filter != "." {
                continue;
            }
        }
        let mut rec = hmf_common::variant_context::VcfRecord::parse(line, &header)?;
        if !small_variant::passes_factory_filters(&rec, pass_only) {
            continue;
        }
        if let Some(v) = small_variant::create_variant(&mut rec, &header, tumor, reference, rna)? {
            if v.reported {
                impacts.push(VariantTranscriptImpact::from_record(&rec)?);
                variants.push(v);
            }
        }
    }
    Ok((variants, impacts))
}

/// File.list() of a directory: the file system's order (Java), or sorted by name (OR3)
pub fn list_dir(dir: &str) -> Result<Vec<String>> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .map_err(|e| anyhow!("java.lang.NullPointerException: cannot list {dir}: {e}"))?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    if !java("OR3") {
        names.sort();
    }
    Ok(names)
}

/// PurpleDataLoader.load
pub fn load(config: &OrangeConfig) -> Result<PurpleData> {
    let tumor = &config.tumor_id;
    let purple_dir = &config.purple_data_directory;
    let qc_file = purple_qc_filename(purple_dir, tumor);
    let purity_tsv = purple_purity_filename(purple_dir, tumor);
    let somatic_driver_tsv = driver_catalog_somatic_filename(purple_dir, tumor);
    let somatic_vcf = match &config.purple_rna_append_vcf {
        Some(v) => v.clone(),
        None => resolve_vcf_path(&format!("{}{tumor}.purple.somatic.vcf.gz", dir_sep(purple_dir))),
    };
    let germline_driver_tsv = driver_catalog_germline_filename(purple_dir, tumor);
    let germline_vcf = resolve_vcf_path(&format!("{}{tumor}.purple.germline.vcf.gz", dir_sep(purple_dir)));
    let gene_cn_tsv = gene_copy_number_filename(purple_dir, tumor);
    let germline_del_tsv = germline_amp_del_filename(purple_dir, tumor);
    let arm_tsv = chr_arm_copy_number_filename(purple_dir, tumor);

    let (qc, purity) = purple_files::read_purity_context_with_qc(&qc_file, &purity_tsv)?;

    let somatic_drivers: Vec<DriverCatalog> =
        purple_files::read_driver_catalog(&somatic_driver_tsv)?.into_iter().filter(|d| d.reported_status != "NOT_REPORTED").collect();

    let samples = small_variant::vcf_samples(&somatic_vcf)?;
    let rna = config.rna_sample_id.as_deref().filter(|r| samples.iter().any(|s| s == r));
    let (somatic_variants, somatic_transcript_impacts) =
        read_reported_variants(&somatic_vcf, tumor, config.reference_id.as_deref(), rna, true)?;

    let somatic_gene_copy_numbers = purple_files::read_gene_copy_numbers(&gene_cn_tsv)?;
    let chr_arm_copy_numbers = purple_files::read_chr_arm_copy_numbers(&arm_tsv)?;

    let mut germline_drivers = None;
    let mut germline_variants = None;
    let mut germline_transcript_impacts = Vec::new();
    let mut germline_amp_dels = None;
    if config.has_reference() {
        germline_drivers = Some(
            purple_files::read_driver_catalog(&germline_driver_tsv)?.into_iter().filter(|d| d.reported_status != "NOT_REPORTED").collect(),
        );
        let (v, i) = read_reported_variants(&germline_vcf, tumor, config.reference_id.as_deref(), None, false)?;
        germline_variants = Some(v);
        germline_transcript_impacts = i;
        germline_amp_dels = Some(
            purple_files::read_germline_amp_dels(&germline_del_tsv)?
                .into_iter()
                .filter(|x| x.filter == "PASS")
                .filter(|x| x.reported == "REPORTED")
                .collect(),
        );
    }

    let mut variant_plots = Vec::new();
    if let Some(dir) = &config.sage_plot_directory {
        for file in list_dir(dir)? {
            if file.ends_with(".html.gz") {
                variant_plots.push(format!("{dir}{file}"));
            }
        }
        hmf_common::debug!(" loaded {} Sage-vis plots from {}", variant_plots.len(), dir);
    }

    Ok(PurpleData {
        qc,
        purity,
        somatic_drivers,
        germline_drivers,
        somatic_variants,
        somatic_transcript_impacts,
        germline_variants,
        germline_transcript_impacts,
        somatic_gene_copy_numbers,
        germline_amp_dels,
        chr_arm_copy_numbers,
        variant_plots,
    })
}

// ---------------------------------------------------------------------------------------------
// conversions

/// DriverInterpretation.interpret
pub fn driver_interpretation(likelihood: f64) -> &'static str {
    if likelihood >= 0.8 {
        "HIGH"
    } else if likelihood >= 0.2 {
        "MEDIUM"
    } else {
        "LOW"
    }
}

const PURPLE_DRIVER_TYPES: [&str; 10] =
    ["AMP", "PARTIAL_AMP", "DEL", "HET_DEL", "LOH", "MUTATION", "GERMLINE_MUTATION", "GERMLINE_DELETION", "GERMLINE_AMP", "UNKNOWN"];
const LIKELIHOOD_METHODS: [&str; 10] =
    ["BIALLELIC", "HOTSPOT", "INFRAME", "DNDS", "SPLICE_REGION", "AMP", "DEL", "DISRUPTION", "GERMLINE", "NONE"];

fn value_of(value: &str, names: &[&str], enum_name: &str) -> Result<String> {
    if names.contains(&value) {
        Ok(value.to_string())
    } else {
        bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.datamodel.{enum_name}.{value}")
    }
}

/// PurpleConversion.convert(DriverCatalog)
pub fn convert_driver(d: &DriverCatalog) -> Result<PurpleDriver> {
    Ok(PurpleDriver {
        gene: d.gene.clone(),
        transcript: d.transcript.clone(),
        driver_type: value_of(&d.driver, &PURPLE_DRIVER_TYPES, "purple.PurpleDriverType")?,
        driver_likelihood: d.driver_likelihood,
        likelihood_method: value_of(&d.likelihood_method, &LIKELIHOOD_METHODS, "purple.PurpleLikelihoodMethod")?,
        reported_status: value_of(&d.reported_status, &["NONE", "NOT_REPORTED", "REPORTED"], "driver.ReportedStatus")?,
        driver_interpretation: driver_interpretation(d.driver_likelihood).to_string(),
        is_canonical: d.is_canonical,
        category: match d.category.as_str() {
            "ONCO" | "TSG" => d.category.clone(),
            c => bail!("java.lang.IllegalArgumentException: No enum constant DriverCategory.{c}"),
        },
    })
}

/// PurpleConversion.convert(PurpleQC)
fn convert_qc(qc: &PurpleQc) -> Result<PurpleQC> {
    for a in &qc.germline_aberrations {
        value_of(a, &PURPLE_GERMLINE_ABERRATIONS, "purple.PurpleGermlineAberration")?;
    }
    Ok(PurpleQC {
        status: enum_set(&qc.status, &PURPLE_QC_STATUSES),
        germline_aberrations: enum_set(&qc.germline_aberrations, &PURPLE_GERMLINE_ABERRATIONS),
        amber_mean_depth: qc.amber_mean_depth,
        contamination: qc.contamination,
        total_copy_number_segments: qc.copy_number_segments,
        unsupported_copy_number_segments: qc.unsupported_copy_number_segments,
        deleted_genes: qc.deleted_genes,
    })
}

/// PurpleConversion.convert(GermlineStatus): PurpleGermlineStatus.valueOf
fn convert_germline_status(s: &str) -> Result<String> {
    value_of(s, &PURPLE_GERMLINE_STATUSES, "purple.PurpleGermlineStatus")
}

/// the effects of a transcript impact as the data model's EnumSet<PurpleVariantEffect>
fn convert_effects(effects: &str) -> Vec<String> {
    enum_set(&effects_to_list(effects), &PURPLE_VARIANT_EFFECTS)
}

// ---------------------------------------------------------------------------------------------
// variants

/// MiscToolFiles.generateSageVisFilePrefix
fn sage_vis_prefix(v: &SmallVariant) -> String {
    let chr = if v.chromosome.starts_with("chr") { v.chromosome.clone() } else { format!("chr{}", v.chromosome) };
    format!(".sage.{}_{}_{}_{}_{}", v.gene, chr, v.position, v.ref_bases, v.alt)
}

fn find_variant_plot(v: &SmallVariant, plots: Option<&[String]>) -> Option<String> {
    let plots = plots?;
    let prefix = sage_vis_prefix(v);
    if java("OR4") {
        return plots.iter().find(|p| p.contains(&prefix)).cloned();
    }
    // OR4: the variant must end at the prefix's end (C>A is not C>AT)
    plots
        .iter()
        .find(|p| p.match_indices(&prefix).any(|(i, _)| matches!(p[i + prefix.len()..].chars().next(), None | Some('_') | Some('.'))))
        .cloned()
}

/// PurpleVariantFactory.findPhasedVariantIds: the local phase set ids the variant shares with another
/// (reported) variant within 1000 bases, in HashSet<Integer> order
fn find_phased_variant_ids(index: usize, all: &[SmallVariant]) -> Option<Vec<i32>> {
    let variant = &all[index];
    let lps = variant.local_phase_sets.as_ref().filter(|l| !l.is_empty())?;
    let mut ids: JavaHashSet<i32> = JavaHashSet::new();
    for (i, other) in all.iter().enumerate() {
        if i == index {
            continue;
        }
        let Some(other_lps) = other.local_phase_sets.as_ref().filter(|l| !l.is_empty()) else { continue };
        if other.chromosome != variant.chromosome {
            continue;
        }
        if (variant.position - other.position).abs() <= 1000 {
            for id in lps {
                if other_lps.contains(id) {
                    ids.insert(*id);
                }
            }
        }
    }
    if ids.is_empty() { None } else { Some(ids.to_vec()) }
}

fn build_purple_variant(
    index: usize,
    all: &[SmallVariant],
    impacts: &[VariantTranscriptImpact],
    driver: Option<&DriverCatalog>,
    germline: bool,
    plot_filename: Option<String>,
) -> Result<PurpleVariant> {
    let v = &all[index];
    let alt_infos = AltTranscriptReportableInfo::parse_all(&v.other_reported_effects)?;
    let canonical = PurpleTranscriptImpact {
        transcript: v.canonical_transcript.clone(),
        hgvs_coding_impact: v.canonical_hgvs_coding_impact.clone(),
        hgvs_protein_impact: v.canonical_hgvs_protein_impact.clone(),
        affected_codon: Some(v.canonical_affected_codon),
        affected_exon: Some(v.canonical_affected_exon),
        in_splice_region: v.splice_region,
        effects: convert_effects(&v.canonical_effect),
        coding_effect: v.canonical_coding_effect.clone(),
        reported: v.reported,
    };
    let mut other_impacts = Vec::new();
    for info in &alt_infos {
        if !v.reportable_transcripts.contains(&info.trans_name) {
            continue;
        }
        // OR7: Java gives another transcript the canonical transcript's codon, exon, splice region and coding
        // effect; corrected, they are its own (PAVE_TI and the reportable effect's coding effect)
        let own = if java("OR7") { None } else { impacts.iter().find(|t| t.transcript == info.trans_name) };
        other_impacts.push(PurpleTranscriptImpact {
            transcript: info.trans_name.clone(),
            hgvs_coding_impact: info.hgvs_coding.clone(),
            hgvs_protein_impact: info.hgvs_protein.clone(),
            affected_codon: own.map(|t| t.affected_codon).or(canonical.affected_codon),
            affected_exon: own.map(|t| t.affected_exon).or(canonical.affected_exon),
            in_splice_region: own.map(|t| t.splice_region).unwrap_or(canonical.in_splice_region),
            effects: convert_effects(&info.effects),
            coding_effect: if java("OR7") { canonical.coding_effect.clone() } else { info.effect.clone() },
            reported: true,
        });
    }
    Ok(PurpleVariant {
        variant_type: v.variant_type.to_string(),
        gene: v.gene.clone(),
        chromosome: v.chromosome.clone(),
        position: v.position,
        ref_bases: v.ref_bases.clone(),
        alt: v.alt.clone(),
        driver_likelihood: driver.map(|d| d.driver_likelihood).unwrap_or(0.0),
        worst_coding_effect: v.worst_coding_effect.clone(),
        canonical_impact: canonical,
        other_impacts,
        hotspot: v.hotspot.to_string(),
        tumor_depth: AllelicDepth {
            total_read_count: v.allelic_depth.total_read_count,
            allele_read_count: v.allelic_depth.allele_read_count,
        },
        rna_depth: v.rna_depth.map(|d| AllelicDepth { total_read_count: d.total_read_count, allele_read_count: d.allele_read_count }),
        adjusted_copy_number: v.adjusted_copy_number,
        adjusted_vaf: v.adjusted_vaf,
        minor_allele_copy_number: v.minor_allele_copy_number,
        variant_copy_number: v.variant_copy_number,
        biallelic: v.biallelic,
        biallelic_probability: Some(v.biallelic_probability),
        genotype_status: v.genotype_status.to_string(),
        repeat_count: v.repeat_count,
        subclonal_likelihood: v.subclonal_likelihood,
        somatic_likelihood: Some(v.somatic_likelihood.clone()),
        local_phase_sets: find_phased_variant_ids(index, all),
        clinvar_pathogenicity: if germline { Some(v.pathogenicity.clone()) } else { None },
        gnomad_frequency: if germline { Some(v.gnomad_frequency) } else { None },
        plot_filename,
    })
}

fn build_purple_variants(
    variants: Option<&[SmallVariant]>,
    impacts: &[Vec<VariantTranscriptImpact>],
    drivers: Option<&[DriverCatalog]>,
    germline: bool,
    plots: Option<&[String]>,
) -> Result<Option<Vec<PurpleVariant>>> {
    let Some(variants) = variants else { return Ok(None) };
    let required = if germline { "GERMLINE_MUTATION" } else { "MUTATION" };
    let mut out = Vec::with_capacity(variants.len());
    for (i, v) in variants.iter().enumerate() {
        let drivers = drivers.ok_or_else(|| anyhow!("java.lang.NullPointerException: no drivers"))?;
        let driver = drivers.iter().find(|d| d.driver == required && d.gene == v.gene);
        let plot = find_variant_plot(v, plots);
        out.push(build_purple_variant(i, variants, &impacts[i], driver, germline, plot)?);
    }
    Ok(Some(out))
}

// ---------------------------------------------------------------------------------------------
// gains and deletions

/// PurpleInterpreter.findChrArmCopyNumber: the arm (P when the band starts with "p") of the gene's chromosome,
/// or one of copy number 0 when the file has none (short arms)
pub(crate) fn find_chr_arm_copy_number(arms: &[ChrArmCopyNumber], gcn: &GeneCopyNumber) -> Result<ChrArmCopyNumber> {
    let arm = if gcn.chromosome_band.starts_with('p') { "P" } else { "Q" };
    let chromosome = purple_files::human_chromosome(&gcn.chromosome)?;
    Ok(arms.iter().find(|a| a.chromosome == chromosome && a.arm == arm).cloned().unwrap_or(ChrArmCopyNumber {
        chromosome,
        arm: arm.to_string(),
        mean_copy_number: 0.0,
        median_copy_number: 0.0,
        min_copy_number: 0.0,
        max_copy_number: 0.0,
    }))
}

fn expression_fields(gene: &str, isofox: Option<&IsofoxData>) -> (Option<f64>, Option<f64>, Option<f64>) {
    let Some(e) = isofox.and_then(|i| i.gene_expressions.iter().find(|e| e.gene_name == gene)) else { return (None, None, None) };
    let tpm = e.tpm;
    let fold = if e.median_tpm_cohort > 0.0 { tpm / e.median_tpm_cohort } else { 0.0 };
    (Some(tpm), Some(e.percentile_cohort), Some(fold))
}

fn find_linx_driver_plot(driver: &DriverCatalog, linx: &LinxData) -> Option<String> {
    const AMP: [&str; 3] = ["GAIN", "GAIN_ARM", "GAIN_CHR"];
    const DEL: [&str; 6] = ["DEL", "LOH", "LOH_SV_CENTRO", "LOH_SV_TELO", "LOH_ARM", "LOH_CHR"];
    let mut cluster_id = None;
    for d in &linx.somatic_driver_data {
        if d.gene != driver.gene {
            continue;
        }
        if driver.driver == "AMP" || driver.driver == "PARTIAL_AMP" {
            if AMP.contains(&d.event_type.as_str()) {
                cluster_id = Some(d.cluster_id);
                break;
            }
        } else if driver.driver == "DEL" && DEL.contains(&d.event_type.as_str()) {
            cluster_id = Some(d.cluster_id);
            break;
        }
    }
    cluster_id.and_then(|c| find_reportable_linx_plot(&linx.reportable_event_plots, c))
}

fn somatic_gains_dels(data: &PurpleData, linx: &LinxData, isofox: Option<&IsofoxData>) -> Result<Vec<PurpleGainDeletion>> {
    const AMP_DEL_TYPES: [&str; 5] = ["AMP", "PARTIAL_AMP", "DEL", "HET_DEL", "LOH"];
    let ploidy = data.purity.ploidy;
    let mut out = Vec::new();
    for driver in &data.somatic_drivers {
        if !AMP_DEL_TYPES.contains(&driver.driver.as_str()) {
            continue;
        }
        let Some(gcn) = data.somatic_gene_copy_numbers.iter().find(|g| g.gene_name == driver.gene && driver.transcript == g.trans_name)
        else {
            continue;
        };
        let arm = find_chr_arm_copy_number(&data.chr_arm_copy_numbers, gcn)?;
        let (tpm, tpm_percentile, tpm_fold_change) = expression_fields(&driver.gene, isofox);
        let relative = if gcn.relative_min_copy_number == -1.0 { gcn.min_copy_number / ploidy } else { gcn.relative_min_copy_number };
        out.push(PurpleGainDeletion {
            driver: convert_driver(driver)?,
            germline_amp_del_fields: None,
            chromosome: driver.chromosome.clone(),
            chromosome_band: driver.chromosome_band.clone(),
            gene_range: if driver.driver == "AMP" { "FULL" } else { "PARTIAL" }.to_string(),
            exon_start: None,
            exon_end: None,
            min_copy_number: java_max(0.0, driver.min_copy_number),
            max_copy_number: java_max(0.0, driver.max_copy_number),
            relative_copy_number: relative,
            arm_copy_number: arm.mean_copy_number,
            min_minor_allele_copies: gcn.min_minor_allele_copy_number,
            tpm,
            tpm_percentile,
            tpm_fold_change,
            plot_filename: find_linx_driver_plot(driver, linx),
        });
    }
    Ok(out)
}

/// Math.min
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }
    if a == 0.0 && b == 0.0 {
        return if a.is_sign_negative() || b.is_sign_negative() { -0.0 } else { 0.0 };
    }
    if a <= b { a } else { b }
}

/// GermlineGainDeletionFactory.createGermlineGainDeletions
pub(crate) fn germline_gains_dels(
    amp_dels: &[GermlineAmpDel],
    drivers: &[PurpleDriver],
    gene_copy_numbers: &[GeneCopyNumber],
    arms: &[ChrArmCopyNumber],
    isofox: Option<&IsofoxData>,
) -> Result<Vec<PurpleGainDeletion>> {
    let mut out = Vec::new();
    for driver in drivers {
        if driver.reported_status != "REPORTED" {
            continue;
        }
        if driver.driver_type != "GERMLINE_DELETION" && driver.driver_type != "GERMLINE_AMP" {
            continue;
        }
        let matched: Vec<&GermlineAmpDel> = amp_dels.iter().filter(|x| x.gene_name == driver.gene).collect();
        let gcn = gene_copy_numbers.iter().find(|g| g.gene_name == driver.gene && driver.transcript == g.trans_name);
        let Some(gcn) = gcn else { continue };
        if matched.is_empty() {
            continue;
        }
        let arm = find_chr_arm_copy_number(arms, gcn)?;
        let (tpm, tpm_percentile, tpm_fold_change) = expression_fields(&driver.gene, isofox);
        let first = matched[0];
        let covers_transcript = matched.iter().all(|x| !x.is_partial);
        let mut min_copies = first.tumor_copy_number;
        let mut max_copies =
            if covers_transcript { first.tumor_copy_number } else { java_max(gcn.max_copy_number, first.tumor_copy_number) };
        let mut germline_min_copies = first.germline_copy_number;
        if matched.len() > 1 {
            for x in &matched {
                min_copies = java_min(x.tumor_copy_number, min_copies);
                max_copies = java_max(x.tumor_copy_number, max_copies);
                germline_min_copies = java_min(x.germline_copy_number, germline_min_copies);
            }
        }
        out.push(PurpleGainDeletion {
            driver: driver.clone(),
            germline_amp_del_fields: Some(GermlineAmpDelFields {
                germline_status: convert_germline_status(&first.normal_status)?,
                somatic_status: convert_germline_status(&first.tumor_status)?,
                germline_min_copy_number: germline_min_copies,
            }),
            chromosome: first.chromosome.clone(),
            chromosome_band: first.chromosome_band.clone(),
            gene_range: if first.is_partial { "PARTIAL" } else { "FULL" }.to_string(),
            exon_start: Some(first.exon_start),
            exon_end: Some(first.exon_end),
            min_copy_number: min_copies,
            max_copy_number: max_copies,
            relative_copy_number: gcn.relative_min_copy_number,
            arm_copy_number: arm.mean_copy_number,
            min_minor_allele_copies: gcn.min_minor_allele_copy_number,
            tpm,
            tpm_percentile,
            tpm_fold_change,
            plot_filename: None,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// fit, characteristics, arms

fn create_fit(data: &PurpleData) -> Result<PurpleFit> {
    let p = &data.purity;
    Ok(PurpleFit {
        qc: convert_qc(&data.qc)?,
        fitted_purity_method: p.fit_method.clone(),
        purity: p.purity,
        min_purity: p.min_purity,
        max_purity: p.max_purity,
        ploidy: p.ploidy,
        min_ploidy: p.min_ploidy,
        max_ploidy: p.max_ploidy,
    })
}

fn create_characteristics(data: &PurpleData) -> PurpleCharacteristics {
    let p = &data.purity;
    PurpleCharacteristics {
        whole_genome_duplication: p.whole_genome_duplication,
        microsatellite_indels_per_mb: p.ms_indels_per_mb,
        microsatellite_status: p.ms_status.clone(),
        tumor_mutational_burden_per_mb: p.tmb_per_mb,
        tumor_mutational_burden_status: p.tmb_status.clone(),
        tumor_mutational_load: p.tml,
        tumor_mutational_load_status: p.tml_status.clone(),
        sv_tumor_mutational_burden: p.sv_tumor_mutational_burden,
        loh_percentage: data.qc.loh_percent,
    }
}

/// ChrArmCopyNumber.driverInterpretation
fn arm_driver_interpretation(mean: f64, ploidy: f64) -> Option<&'static str> {
    let (high, low) = (0.4, 0.25);
    if mean > (1.0 + high) * ploidy || mean < (1.0 - high) * ploidy {
        return Some("HIGH");
    }
    if mean > (1.0 + low) * ploidy || mean < (1.0 - low) * ploidy {
        return Some("LOW");
    }
    None
}

fn build_arm_abberations(data: &PurpleData) -> Result<Vec<PurpleChrArmCopyNumber>> {
    let ploidy = data.purity.ploidy;
    let gender = data.purity.sex.as_str();
    let mut out = Vec::new();
    for arm in &data.chr_arm_copy_numbers {
        let mut expected = ploidy;
        if arm.chromosome == "Y" && gender == "FEMALE" {
            continue;
        }
        let allosome = arm.chromosome == "X" || arm.chromosome == "Y";
        if allosome && gender == "MALE" {
            expected *= 0.5;
        }
        let cn_type = if arm.mean_copy_number > (1.0 + 0.25) * expected {
            "GAIN"
        } else if arm.mean_copy_number < (1.0 - 0.25) * expected {
            "LOSS"
        } else {
            continue;
        };
        let interpretation = arm_driver_interpretation(arm.mean_copy_number, expected)
            .ok_or_else(|| anyhow!("java.lang.NullPointerException: arm driver interpretation"))?;
        out.push(PurpleChrArmCopyNumber {
            chromosome: format!("chr{}", arm.chromosome),
            arm: arm.arm.clone(),
            cn_type: cn_type.to_string(),
            copy_number: arm.median_copy_number,
            relative_copy_number: if expected > 0.0 { arm.mean_copy_number / expected } else { 0.0 },
            driver_interpretation: interpretation.to_string(),
        });
    }
    Ok(out)
}

/// PurpleInterpreter.interpret
pub fn interpret(data: &PurpleData, linx: &LinxData, isofox: Option<&IsofoxData>) -> Result<PurpleRecord> {
    let somatic_variants = build_purple_variants(
        Some(&data.somatic_variants),
        &data.somatic_transcript_impacts,
        Some(&data.somatic_drivers),
        false,
        Some(&data.variant_plots),
    )?
    .unwrap();
    let germline_variants = build_purple_variants(
        data.germline_variants.as_deref(),
        &data.germline_transcript_impacts,
        data.germline_drivers.as_deref(),
        true,
        None,
    )?;
    let germline_drivers: Option<Vec<PurpleDriver>> = match &data.germline_drivers {
        Some(d) => Some(d.iter().map(convert_driver).collect::<Result<_>>()?),
        None => None,
    };
    let somatic_gains_dels = somatic_gains_dels(data, linx, isofox)?;
    let arm_abberations = build_arm_abberations(data)?;
    let germline_gains_dels = match &data.germline_amp_dels {
        Some(amp_dels) => {
            let drivers = germline_drivers.as_deref().ok_or_else(|| anyhow!("java.lang.NullPointerException: germline drivers"))?;
            Some(germline_gains_dels(amp_dels, drivers, &data.somatic_gene_copy_numbers, &data.chr_arm_copy_numbers, isofox)?)
        }
        None => None,
    };
    Ok(PurpleRecord {
        fit: create_fit(data)?,
        characteristics: create_characteristics(data),
        somatic_drivers: data.somatic_drivers.iter().map(convert_driver).collect::<Result<_>>()?,
        germline_drivers,
        somatic_variants,
        germline_variants,
        somatic_gains_dels,
        germline_gains_dels,
        arm_copy_number_abberations: arm_abberations,
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    #[test]
    fn driver_interpretations() {
        // DriverInterpretationTest
        assert_eq!(driver_interpretation(0.0), "LOW");
        assert_eq!(driver_interpretation(0.2), "MEDIUM");
        assert_eq!(driver_interpretation(0.5), "MEDIUM");
        assert_eq!(driver_interpretation(0.8), "HIGH");
        assert_eq!(driver_interpretation(1.0), "HIGH");
    }

    #[test]
    fn java_min_semantics() {
        assert!(java_min(f64::NAN, 1.0).is_nan());
        assert!(java_min(0.0, -0.0).is_sign_negative());
        assert_eq!(java_min(1.0, 2.0), 1.0);
    }
}
