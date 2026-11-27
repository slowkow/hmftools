//! The germline package: GermlineVariants (genotype, hotspot, purity, low-VAF rescue, low tumor VCN,
//! pathogenicity and reported enrichment of the PAVE germline VCF, written to purple.germline.vcf.gz),
//! GermlineDrivers, GermlineAmpDelFinder and GermlineAmpDel.

use anyhow::Result;
use std::collections::{BTreeMap, BTreeSet, HashMap};

use hmf_common::commons_math::Poisson;
use hmf_common::driver_panel::{DriverGene, GermlineReporting as Reporting};
use hmf_common::ensembl::EnsemblDataCache;
use hmf_common::java_hash::JavaHashSet;
use hmf_common::variant_context::{self as vcf, EncodeOptions, Value, VcfHeader, VcfReader, VcfWriter};

use crate::constants::*;
use crate::copynumber::PurpleCopyNumber;
use crate::drivers::{category_of, has_transcript_coding_effect, reportable_transcript_list, worst_reportable_coding_effect, DriverCatalog};
use crate::fitting::PurityAdjuster;
use crate::gene::GeneCopyNumber;
use crate::reference::{GermlineAmpDelFrequencyCache, ReferenceData};
use crate::region::ObservedRegion;
use crate::somatic::hotspot_enrich;
use crate::somatic_stream::RegionLookup;
use crate::sv::StructuralVariant;
use crate::types::*;
use crate::variant::{driver_impact, pathogenic_summary, AllelicDepth, Variant};

pub struct GermlineVariant {
    pub v: Variant,
    /// GermlineVariant.mFilters: what is written; the record's filters are the context's, which the
    /// decorator (isPass of the reported enrichment) reads
    pub filters: BTreeSet<String>,
}

impl GermlineVariant {
    pub fn is_pass(&self) -> bool {
        self.filters.is_empty() || (self.filters.len() == 1 && self.filters.contains("PASS"))
    }
    fn add_context_filter(&mut self, f: &str) {
        let filters = self.v.rec.filters.get_or_insert_with(Vec::new);
        if !filters.iter().any(|x| x == f) {
            filters.push(f.to_string());
        }
        self.v.refresh();
    }
}

#[derive(Clone, Copy, PartialEq)]
enum GermlineGenotypeStatus {
    Het,
    HomAlt,
    LowVaf,
}

/// GermlineGenotypeEnrichment.status
fn genotype_status(d: AllelicDepth) -> GermlineGenotypeStatus {
    if d.allele_read_count == d.total_read_count {
        return GermlineGenotypeStatus::HomAlt;
    }
    if greater_than(d.allele_read_count as f64, 0.75 * d.total_read_count as f64) {
        let p = Poisson::new(d.total_read_count as f64 / 2.0).cumulative_probability(d.total_read_count - d.allele_read_count);
        return if less_than(p, 0.005) { GermlineGenotypeStatus::HomAlt } else { GermlineGenotypeStatus::Het };
    }
    if less_than(d.allele_read_count as f64, 0.3 * d.total_read_count as f64) {
        let p = Poisson::new(d.total_read_count as f64 / 2.0).cumulative_probability(d.allele_read_count);
        return if less_than(p, 0.002) { GermlineGenotypeStatus::LowVaf } else { GermlineGenotypeStatus::Het };
    }
    GermlineGenotypeStatus::Het
}

/// VariantContextDecorator.isGermlinePathogenic
pub fn is_germline_pathogenic(v: &Variant) -> bool {
    let (_, status) = pathogenic_summary(&v.rec);
    if status == Pathogenicity::BenignBlacklist {
        return false;
    }
    if v.hotspot_flag() || status.is_pathogenic() {
        return true;
    }
    status == Pathogenicity::Unknown && matches!(v.impact.canonical_coding_effect, CodingEffect::NonsenseOrFrameshift | CodingEffect::Splice)
}

fn germline_header(input: &VcfHeader, version: &str) -> VcfHeader {
    let mut h = input.clone();
    // GermlinePurityEnrichment
    h.add_line(vcf::meta_line("purpleVersion", version));
    h.add_line(vcf::info_line("PURPLE_AF", "1", "Float", "Purity adjusted variant allelic frequency"));
    h.add_line(vcf::info_line("PURPLE_CN", "1", "Float", "Purity adjusted copy number surrounding variant location"));
    h.add_line(vcf::info_line("PURPLE_VCN", "1", "Float", "Purity adjusted variant copy number"));
    h.add_line(vcf::info_line("PURPLE_MACN", "1", "Float", "Purity adjusted minor allele ploidy surrounding variant location"));
    h.add_line(vcf::info_line("BIALLELIC", "0", "Flag", "If the probability this variant is biallelic equals or is greater than 0.50"));
    h.add_line(vcf::info_line("REPORTED", "0", "Flag", "Variant is reported in the driver catalog"));
    h.add_line(vcf::info_line("REPORTABLE_TRANSCRIPTS", "1", "String", "List of reportable transcript when non-canonical are reportable"));
    // HotspotEnrichment
    h.add_line(vcf::info_line("HOTSPOT", "0", "Flag", "Site is at a known hotspot location"));
    h.add_line(vcf::info_line("NEAR_HOTSPOT", "0", "Flag", "Variant within 5 bases of hotspot"));
    // GermlineLowTumorVCNFilter
    h.add_line(vcf::filter_line("LOW_TUMOR_VCN", "Germline variant has very low tumor variant copy number"));
    // GermlineReportedEnrichment
    h.add_line(vcf::info_line("REPORTED", "0", "Flag", "Variant is reported in the driver catalog"));
    // GermlineGenotypeEnrichment
    h.add_line(vcf::filter_line("LOW_VAF", "Germline variant has very low allelic frequency"));
    // GermlinePathogenicEnrichment
    let names: Vec<&str> = Pathogenicity::ALL.iter().map(|p| p.name()).collect();
    h.add_line(vcf::info_line("PATH", "1", "String", &format!("Pathogenicity [{}]", names.join(","))));
    h
}

/// GermlineVariants.processAndWrite: the candidate reportable variants
#[allow(clippy::too_many_arguments)]
pub fn process_germline_variants(
    reference_id: &str,
    tumor_id: &str,
    path: &str,
    output: &str,
    version: &str,
    reference: &ReferenceData,
    adjuster: Option<&PurityAdjuster>,
    copy_numbers: &[PurpleCopyNumber],
    somatic_reported_genes: &std::collections::HashSet<String>,
) -> Result<Vec<GermlineVariant>> {
    if path.is_empty() {
        return Ok(Vec::new());
    }
    let mut reader = VcfReader::open(path)?;
    let input_header = reader.header.clone();
    let mut variants: Vec<GermlineVariant> = Vec::new();
    while let Some(rec) = reader.next_record()? {
        // GermlineVariant(context): the decorator; no genotype is read yet
        let v = Variant {
            vtype: crate::variant::variant_type(&rec),
            filter: rec.display_filter(),
            tier: crate::variant::variant_tier(&rec),
            impact: crate::variant::VariantImpact::from_record(&rec),
            tumor_depth: None,
            reference_depth: None,
            rec,
        };
        let filters: BTreeSet<String> = v.rec.filters.clone().unwrap_or_default().into_iter().collect();
        variants.push(GermlineVariant { v, filters });
    }
    crate::info!("load {} germline variants from {}", variants.len(), path);
    let header = germline_header(&input_header, version);
    let cn_lookup = RegionLookup::new(copy_numbers.iter().enumerate().map(|(i, c)| (c.chromosome.clone(), c.start, c.end, i)));
    for gv in variants.iter_mut() {
        genotype_enrich(gv, &input_header, reference_id, tumor_id);
        hotspot_enrich(&mut gv.v.rec, &reference.germline_hotspots);
        purity_enrich(gv, &input_header, reference_id, tumor_id, adjuster, copy_numbers, &cn_lookup);
        rescue_low_vaf(gv, &input_header, reference_id);
        low_tumor_vcn_filter(gv);
        let status = pathogenic_summary(&gv.v.rec).1;
        gv.v.rec.info.put("PATH", Value::Str(status.name().to_string()));
    }
    let candidates = reported_enrich(&mut variants, &reference.driver_genes.driver_gene_list, somatic_reported_genes);
    let mut writer = VcfWriter::create(output, &header, EncodeOptions { allow_missing_fields: true })?;
    for gv in &variants {
        let mut r = gv.v.rec.clone();
        // VariantContextBuilder.filters(variant.filters()): a set, written sorted, PASS when empty
        r.filters = Some(gv.filters.iter().cloned().collect());
        writer.write(&r)?;
    }
    writer.close()?;
    Ok(candidates.into_iter().map(|i| GermlineVariant { v: variants[i].v.clone(), filters: variants[i].filters.clone() }).collect())
}

fn genotype_enrich(gv: &mut GermlineVariant, header: &VcfHeader, reference_id: &str, tumor_id: &str) {
    let gs = gv.v.rec.genotypes(header);
    let germline = gs.iter().find(|g| g.sample == reference_id).unwrap_or_else(|| panic!("NullPointerException: no genotype {reference_id}")).clone();
    let tumor = gs.iter().find(|g| g.sample == tumor_id).unwrap_or_else(|| panic!("NullPointerException: no genotype {tumor_id}")).clone();
    let gd = AllelicDepth::from_genotype(&germline);
    let td = AllelicDepth::from_genotype(&tumor);
    let (gs_status, ts_status) = (genotype_status(gd), genotype_status(td));
    let status = if gs_status == GermlineGenotypeStatus::HomAlt {
        if ts_status == GermlineGenotypeStatus::HomAlt { GermlineGenotypeStatus::HomAlt } else { GermlineGenotypeStatus::Het }
    } else {
        gs_status
    };
    let mut updated: Vec<hmf_common::variant_context::Genotype> = gs.iter().filter(|g| g.sample != reference_id).cloned().collect();
    let mut g = germline;
    g.alleles = if status == GermlineGenotypeStatus::HomAlt { vec![Some(1), Some(1)] } else { vec![Some(0), Some(1)] };
    updated.push(g);
    *gs = updated;
    if status == GermlineGenotypeStatus::LowVaf {
        gv.add_context_filter("LOW_VAF");
        gv.filters.remove("PASS");
        gv.filters.insert("LOW_VAF".into());
    }
    gv.v.refresh();
}

fn purity_enrich(
    gv: &mut GermlineVariant,
    header: &VcfHeader,
    reference_id: &str,
    tumor_id: &str,
    adjuster: Option<&PurityAdjuster>,
    copy_numbers: &[PurpleCopyNumber],
    cn_lookup: &RegionLookup<usize>,
) {
    let rec = &mut gv.v.rec;
    let Some(adj) = adjuster else {
        rec.info.put("PURPLE_VCN", Value::Int(1));
        rec.info.put("PURPLE_CN", Value::Int(2));
        rec.info.put("PURPLE_AF", Value::Double(0.5));
        rec.info.put("PURPLE_MACN", Value::Int(1));
        rec.info.put("BIALLELIC", Value::Bool(false));
        return;
    };
    let contig = rec.chrom.clone();
    let gs = rec.genotypes(header);
    let tumor = gs.iter().find(|g| g.sample == tumor_id).cloned();
    let normal = gs.iter().find(|g| g.sample == reference_id).cloned();
    let (Some(tumor), Some(normal)) = (tumor, normal) else { return };
    if !tumor.has_ad() || !Chromosome::contains(&contig) {
        return;
    }
    let td = AllelicDepth::from_genotype(&tumor);
    let hom_alt = normal.genotype_type() == hmf_common::variant_context::GenotypeType::HomVar;
    let het = normal.genotype_type() == hmf_common::variant_context::GenotypeType::Het;
    let Some(&ci) = cn_lookup.get(&contig, rec.pos) else { return };
    let cn = &copy_numbers[ci];
    let copy_number = cn.average_tumor_copy_number;
    let vaf = if hom_alt {
        1.0
    } else if td.total_read_count == 0 || td.allele_read_count == 0 || less_or_equal(copy_number, 0.001) {
        0.0
    } else {
        let raw = td.allele_frequency();
        let constrained = jmax(0.001, copy_number);
        if het {
            adj.purity_adjusted_vaf_with_heterozygous_normal(&cn.chromosome, constrained, raw)
        } else {
            adj.purity_adjusted_vaf_with_homozygous_normal(&cn.chromosome, constrained, raw)
        }
    };
    let vcn = jmax(0.0, vaf * copy_number);
    let biallelic = less_or_equal(copy_number, 0.0) || greater_or_equal(vcn, copy_number - 0.5);
    rec.info.put("PURPLE_VCN", Value::Double(vcn));
    rec.info.put("PURPLE_CN", Value::Double(copy_number));
    rec.info.put("PURPLE_AF", Value::Str(hmf_common::fmt::java_fixed(vaf, 4)));
    rec.info.put("PURPLE_MACN", Value::Double(cn.minor_allele_copy_number()));
    rec.info.put("BIALLELIC", Value::Bool(biallelic));
}

fn rescue_low_vaf(gv: &mut GermlineVariant, header: &VcfHeader, reference_id: &str) {
    if gv.filters.len() == 1 && gv.filters.contains("LOW_VAF") {
        let g = gv.v.rec.genotypes(header).iter().find(|g| g.sample == reference_id).cloned().expect("germline genotype");
        let d = AllelicDepth::from_genotype(&g);
        if d.allele_read_count >= 3 && greater_or_equal(gv.v.copy_number(), 0.5) {
            gv.filters.clear();
            gv.filters.insert("PASS".into());
        }
    }
}

fn low_tumor_vcn_filter(gv: &mut GermlineVariant) {
    if less_than(gv.v.copy_number(), 0.5) {
        let min_qual = if gv.v.is_hotspot() { 45.0 } else { 110.0 };
        if gv.v.rec.phred_scaled_qual() < min_qual {
            gv.filters.insert("LOW_TUMOR_VCN".into());
            gv.add_context_filter("LOW_TUMOR_VCN");
        }
    }
}

/// GermlineReportedEnrichment.flush: REPORTED set; the candidate reportable variants' indices
fn reported_enrich(variants: &mut [GermlineVariant], driver_genes: &[DriverGene], somatic_reported: &std::collections::HashSet<String>) -> Vec<usize> {
    let gene_map: HashMap<&str, &DriverGene> = driver_genes.iter().map(|g| (g.gene.as_str(), g)).collect();
    let mut candidates = Vec::new();
    let mut by_gene: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, gv) in variants.iter().enumerate() {
        let gene = gv.v.gene();
        if !gene_map.contains_key(gene) {
            continue;
        }
        let v = &gv.v;
        let candidate = !gene.is_empty() && v.is_pass() && (is_germline_pathogenic(v) || v.impact.canonical_splice_region);
        if candidate {
            candidates.push(i);
            by_gene.entry(gene.to_string()).or_default().push(i);
        }
    }
    for (gene, indices) in &by_gene {
        let dg = gene_map[gene.as_str()];
        let mut multiple_unphased = false;
        if indices.len() > 1 {
            for &i in indices {
                let lps_i = local_phase_set(&variants[i].v);
                for &j in indices {
                    if i == j || variants[j].v.gene() != variants[i].v.gene() {
                        continue;
                    }
                    let lps_j = local_phase_set(&variants[j].v);
                    if lps_i.is_none() || lps_j.is_none() || lps_i != lps_j {
                        multiple_unphased = true;
                        break;
                    }
                }
            }
        }
        if dg.report_germline_variant != Reporting::None || dg.report_germline_hotspot != Reporting::None {
            let mut multi: std::collections::HashSet<String> = somatic_reported.clone();
            if multiple_unphased {
                multi.insert(gene.clone());
            }
            for &i in indices {
                let v = &variants[i].v;
                if is_reportable(v, dg.report_germline_hotspot, dg.report_germline_variant, &multi) {
                    let t = reportable_transcript_list(v.vtype, &v.impact);
                    let rec = &mut variants[i].v.rec;
                    rec.info.put("REPORTED", Value::Bool(true));
                    if let Some(t) = t {
                        rec.info.put("REPORTABLE_TRANSCRIPTS", Value::Str(t));
                    }
                }
            }
        }
    }
    candidates
}

fn local_phase_set(v: &Variant) -> Option<i32> {
    if !v.rec.info.has("LPS") {
        return None;
    }
    v.rec.info.as_int_list("LPS", 0).first().copied()
}

fn is_reportable(v: &Variant, hotspot_reporting: Reporting, variant_reporting: Reporting, multi: &std::collections::HashSet<String>) -> bool {
    if !is_germline_pathogenic(v) {
        return false;
    }
    let reporting = if v.is_hotspot() { hotspot_reporting } else { variant_reporting };
    match reporting {
        Reporting::None => false,
        Reporting::Any => true,
        _ => {
            if less_than(v.copy_number(), 0.5) {
                return false;
            }
            match reporting {
                Reporting::VariantNotLost => true,
                Reporting::WildtypeLost => v.biallelic() || multi.contains(v.gene()),
                _ => false,
            }
        }
    }
}

/// GermlineDrivers.findDrivers
pub fn germline_drivers(driver_genes: &[DriverGene], variants: &[GermlineVariant], gene_copy_numbers: &HashMap<String, Vec<GeneCopyNumber>>) -> Vec<DriverCatalog> {
    let gene_map: HashMap<&str, &DriverGene> = driver_genes.iter().map(|g| (g.gene.as_str(), g)).collect();
    // Collectors.toSet(): HashSet::new and add, iterated in Java's order
    let mut genes: JavaHashSet<String> = JavaHashSet::new();
    for v in variants {
        genes.insert(v.v.gene().to_string());
    }
    let mut out = Vec::new();
    for gene in genes.iter() {
        let Some(dg) = gene_map.get(gene.as_str()) else { continue };
        let gene_variants: Vec<&GermlineVariant> = variants.iter().filter(|x| x.v.gene() == gene).collect();
        let Some(gcns) = gene_copy_numbers.get(gene) else { continue };
        for gcn in gcns {
            if gcns.len() == 1 || gene_variants.iter().any(|x| has_transcript_coding_effect(&x.v.impact, x.v.vtype, &gcn.trans_name)) {
                out.push(germline_driver(dg, gene, &gene_variants, gcn));
            }
        }
    }
    out
}

fn germline_driver(dg: &DriverGene, gene: &str, variants: &[&GermlineVariant], gcn: &GeneCopyNumber) -> DriverCatalog {
    let mut counts: HashMap<DriverImpact, i32> = HashMap::new();
    let mut coding = false;
    for gv in variants {
        let impact = driver_impact(gv.v.vtype, worst_reportable_coding_effect(&gv.v.impact));
        coding |= gv.v.is_hotspot() || impact != DriverImpact::Unknown;
        *counts.entry(impact).or_insert(0) += 1;
    }
    let c = |i: DriverImpact| *counts.get(&i).unwrap_or(&0);
    let mut status = ReportedStatus::NotReported;
    if dg.report_germline_variant != Reporting::None || dg.report_germline_hotspot != Reporting::None {
        for gv in variants {
            if gv.v.rec.info.has("REPORTED") {
                status = ReportedStatus::Reported;
            }
        }
    }
    DriverCatalog {
        chromosome: variants[0].v.chromosome().to_string(),
        chromosome_band: gcn.chromosome_band.clone(),
        gene: gene.to_string(),
        transcript: gcn.trans_name.clone(),
        is_canonical: gcn.is_canonical,
        category: category_of(dg.likelihood_type),
        driver: DriverType::GermlineMutation,
        likelihood_method: if coding { LikelihoodMethod::Germline } else { LikelihoodMethod::SpliceRegion },
        reported_status: status,
        driver_likelihood: if coding { 1.0 } else { 0.0 },
        missense: c(DriverImpact::Missense),
        nonsense: c(DriverImpact::Nonsense),
        splice: c(DriverImpact::Splice),
        inframe: c(DriverImpact::Inframe),
        frameshift: c(DriverImpact::Frameshift),
        biallelic: variants.iter().any(|x| x.v.biallelic()),
        min_copy_number: gcn.min_copy_number,
        max_copy_number: gcn.max_copy_number,
    }
}

// ---------------------------------------------------------------------------------------------
// GermlineAmpDelFinder

#[derive(Clone, Debug)]
pub struct GermlineAmpDel {
    pub gene_name: String,
    pub transcript: String,
    pub chromosome: String,
    pub chromosome_band: String,
    pub region_start: i32,
    pub region_end: i32,
    pub depth_window_count: i32,
    pub exon_start: i32,
    pub exon_end: i32,
    pub is_partial: bool,
    pub normal_status: GermlineStatus,
    pub tumor_status: GermlineStatus,
    pub germline_copy_number: f64,
    pub tumor_copy_number: f64,
    pub filter: String,
    pub cohort_frequency: i32,
    pub reported: ReportedStatus,
}

/// GermlineAmpDel.write
pub fn write_germline_amp_dels(path: &str, events: &[GermlineAmpDel]) -> Result<()> {
    let mut out = String::from("gene\ttranscript\tchromosome\tchromosomeBand\tregionStart\tregionEnd\tdepthWindowCount\texonStart\texonEnd\tisPartial\tdetectionMethod\tgermlineStatus\ttumorStatus\tgermlineCopyNumber\ttumorCopyNumber\tfilter\tcohortFrequency\treportedStatus\n");
    for e in events {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tSEGMENT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            e.gene_name,
            e.transcript,
            e.chromosome,
            e.chromosome_band,
            e.region_start,
            e.region_end,
            e.depth_window_count,
            e.exon_start,
            e.exon_end,
            e.is_partial,
            e.normal_status,
            e.tumor_status,
            hmf_common::fmt::java_fixed(e.germline_copy_number, 2),
            hmf_common::fmt::java_fixed(e.tumor_copy_number, 2),
            e.filter,
            e.cohort_frequency,
            e.reported
        ));
    }
    std::fs::write(path, out)?;
    Ok(())
}

pub struct AmpDelFinder<'a> {
    pub driver_genes: HashMap<String, &'a DriverGene>,
    pub ensembl: &'a EnsemblDataCache,
    pub cohort: &'a GermlineAmpDelFrequencyCache,
    pub events: Vec<GermlineAmpDel>,
    pub drivers: Vec<DriverCatalog>,
}

fn report_event(dg: &DriverGene, status: GermlineStatus) -> bool {
    match status {
        GermlineStatus::Amplification => dg.report_germline_amplification,
        GermlineStatus::HetDeletion | GermlineStatus::HomDeletion => match dg.report_germline_deletion {
            Reporting::Any | Reporting::VariantNotLost => true,
            other => other == Reporting::WildtypeLost && status == GermlineStatus::HomDeletion,
        },
        _ => false,
    }
}

#[derive(Clone, Copy)]
struct MatchedSv<'b> {
    sv: &'b StructuralVariant,
    is_start: bool,
}

impl<'a> AmpDelFinder<'a> {
    pub fn new(driver_genes: &'a [DriverGene], ensembl: &'a EnsemblDataCache, cohort: &'a GermlineAmpDelFrequencyCache) -> AmpDelFinder<'a> {
        AmpDelFinder { driver_genes: driver_genes.iter().map(|g| (g.gene.clone(), g)).collect(), ensembl, cohort, events: Vec::new(), drivers: Vec::new() }
    }

    /// GermlineAmpDelFinder.findEvents
    pub fn find_events(&mut self, copy_numbers: &[PurpleCopyNumber], fitted: &[ObservedRegion], svs: &[StructuralVariant]) {
        let mut state = (0usize, 0usize, String::new());
        for i in 0..fitted.len() {
            let region = &fitted[i];
            let status = region.germline_status;
            if !matches!(status, GermlineStatus::HomDeletion | GermlineStatus::HetDeletion | GermlineStatus::Amplification) {
                continue;
            }
            let matched = find_matching_copy_number(region, copy_numbers, &mut state);
            let next = fitted.get(i + 1).filter(|n| n.chromosome == region.chromosome);
            let svs_matched = find_matching_svs(region, next, svs, status);
            self.find_overlapping_driver_gene(region, next, matched, svs_matched, status);
        }
    }

    fn find_overlapping_driver_gene(&mut self, region: &ObservedRegion, next: Option<&ObservedRegion>, cn: Option<&PurpleCopyNumber>, svs: [Option<MatchedSv>; 2], status: GermlineStatus) {
        let Some(genes) = self.ensembl.chr_gene_data_map().get(&region.chromosome) else { return };
        let pos = |m: &MatchedSv| m.sv.position(m.is_start).unwrap();
        let start = svs[0].as_ref().map(pos).unwrap_or(region.start);
        let end = svs[1].as_ref().map(pos).unwrap_or(region.end);
        let lower = if svs[0].is_some() { start } else { start - GERMLINE_AMP_DEL_GENE_BUFFER };
        let high = if svs[1].is_some() { end } else { end + GERMLINE_AMP_DEL_GENE_BUFFER };
        let tumor_cn = region.ref_normalised_copy_number;
        let mut tumor_status = status;
        if status == GermlineStatus::HetDeletion && tumor_cn < 0.5 {
            tumor_status = GermlineStatus::HomDeletion;
        } else if status == GermlineStatus::Amplification && tumor_cn >= 1.3 {
            tumor_status = GermlineStatus::Amplification;
        }
        let mut overlapping = Vec::new();
        let mut panel: Vec<&DriverGene> = Vec::new();
        for g in genes {
            if lower > g.gene_end {
                continue;
            }
            if high < g.gene_start {
                break;
            }
            if positions_overlap(g.gene_start, g.gene_end, lower, high) {
                overlapping.push(g);
                if let Some(dg) = self.driver_genes.get(&g.gene_name) {
                    panel.push(dg);
                }
            }
        }
        if overlapping.is_empty() {
            return;
        }
        let cohort = self.cohort.region_frequency(&region.chromosome, region.start, region.end, GERMLINE_AMP_DEL_REGION_MATCH_BUFFER, status == GermlineStatus::Amplification);
        // checkFilters
        let mut filters: Vec<&str> = Vec::new();
        let sv_bounded = region.support.is_sv() || next.is_some_and(|n| n.support.is_sv());
        if !sv_bounded && region.depth_window_count < GERMLINE_AMP_DEL_MIN_DEPTH_WINDOW_COUNT {
            filters.push("LOW_DEPTH_NO_SVS");
        }
        if let Some(c) = cn {
            if status != GermlineStatus::Amplification {
                let inconsistency = region.ref_normalised_copy_number - c.major_allele_copy_number();
                let limit = jmax(GERMLINE_AMP_DEL_CN_CONSISTENCY_MIN, c.major_allele_copy_number() * GERMLINE_AMP_DEL_CN_CONSISTENCY_MACN_PERC);
                if inconsistency > limit {
                    filters.push("INCONSISTENT_CN");
                }
            }
        }
        if cohort >= GERMLINE_AMP_DEL_COHORT_FREQ {
            filters.push("COHORT_FREQ");
        }
        let mut deleted = Vec::new();
        for g in &overlapping {
            let Some(t) = self.ensembl.canonical_transcript_data(&g.gene_id) else { continue };
            let exons: Vec<_> = t.exons.iter().filter(|x| positions_overlap(x.start, x.end, lower, high)).collect();
            if exons.is_empty() {
                continue;
            }
            let rmin = exons.iter().map(|x| x.rank).min().unwrap_or(0);
            let rmax = exons.iter().map(|x| x.rank).max().unwrap_or(0);
            deleted.push((*g, t, rmin, rmax));
        }
        if deleted.is_empty() {
            return;
        }
        let germline_cn = region.observed_normal_ratio * 2.0;
        let filter = if filters.is_empty() { "PASS".to_string() } else { filters.join(";") };
        for (g, t, rmin, rmax) in &deleted {
            let mut reported = ReportedStatus::None;
            if filters.is_empty() {
                for dg in &panel {
                    if dg.gene == g.gene_name {
                        reported = if report_event(dg, status) { ReportedStatus::Reported } else { ReportedStatus::NotReported };
                    }
                }
            }
            let partial = *rmin > 1 || (*rmax as usize) < t.exons.len();
            self.events.push(GermlineAmpDel {
                gene_name: g.gene_name.clone(),
                transcript: t.trans_name.clone(),
                chromosome: region.chromosome.clone(),
                chromosome_band: g.karyotype_band.clone(),
                region_start: start,
                region_end: end,
                depth_window_count: region.depth_window_count,
                exon_start: *rmin,
                exon_end: *rmax,
                is_partial: partial,
                normal_status: region.germline_status,
                tumor_status,
                germline_copy_number: germline_cn,
                tumor_copy_number: tumor_cn,
                filter: filter.clone(),
                cohort_frequency: cohort,
                reported,
            });
        }
        if !filters.is_empty() {
            return;
        }
        for (g, t, _, _) in &deleted {
            let Some(dg) = panel.iter().find(|x| x.gene == g.gene_name) else { continue };
            if self.drivers.iter().any(|x| x.gene == g.gene_name && x.transcript == t.trans_name) {
                continue;
            }
            let reported = if report_event(dg, status) { ReportedStatus::Reported } else { ReportedStatus::NotReported };
            // PU9: Java gives germline amplifications the driver type GERMLINE_DELETION too
            let driver = if status == GermlineStatus::Amplification && !crate::compat::java("PU9") { DriverType::GermlineAmp } else { DriverType::GermlineDeletion };
            self.drivers.push(DriverCatalog {
                chromosome: region.chromosome.clone(),
                chromosome_band: g.karyotype_band.clone(),
                gene: g.gene_name.clone(),
                transcript: t.trans_name.clone(),
                is_canonical: t.is_canonical,
                category: category_of(dg.likelihood_type),
                driver,
                likelihood_method: LikelihoodMethod::Germline,
                reported_status: reported,
                driver_likelihood: 1.0,
                missense: 0,
                nonsense: 0,
                splice: 0,
                inframe: 0,
                frameshift: 0,
                biallelic: region.ref_normalised_copy_number < 0.5,
                min_copy_number: region.ref_normalised_copy_number,
                max_copy_number: region.ref_normalised_copy_number,
            });
        }
    }
}

fn find_matching_copy_number<'b>(region: &ObservedRegion, copy_numbers: &'b [PurpleCopyNumber], state: &mut (usize, usize, String)) -> Option<&'b PurpleCopyNumber> {
    if copy_numbers.is_empty() {
        return None;
    }
    let chromosome = &region.chromosome;
    if &state.2 != chromosome {
        for index in state.0..copy_numbers.len() {
            let c = &copy_numbers[index];
            if &c.chromosome == chromosome && &state.2 != chromosome {
                state.2 = chromosome.clone();
                state.0 = index;
            } else if &state.2 == chromosome && &c.chromosome != chromosome {
                state.1 = index - 1;
                break;
            }
        }
    }
    if state.1 < state.0 {
        state.1 = copy_numbers.len() - 1;
    }
    for index in state.0..=state.1 {
        let c = &copy_numbers[index];
        if &c.chromosome != chromosome {
            break;
        }
        if positions_overlap(c.start, c.end, region.start, region.end) {
            return Some(c);
        }
    }
    None
}

fn find_matching_svs<'b>(region: &ObservedRegion, next: Option<&ObservedRegion>, svs: &'b [StructuralVariant], status: GermlineStatus) -> [Option<MatchedSv<'b>>; 2] {
    let mut matched: [Option<MatchedSv>; 2] = [None, None];
    let preferred = if status == GermlineStatus::Amplification { SvType::Dup } else { SvType::Del };
    for se in 0..2 {
        let mut orientation: i8 = if se == 0 { 1 } else { -1 };
        if status == GermlineStatus::Amplification {
            orientation = -orientation;
        }
        let (mut rs, mut re) = if se == 0 { (region.min_start, region.max_start) } else { (region.end, next.map(|n| n.max_start - 1).unwrap_or(region.end)) };
        rs -= crate::inputs::WINDOW_SIZE;
        re += crate::inputs::WINDOW_SIZE;
        let leg_matches = |leg: &crate::sv::SvLeg| leg.chromosome == region.chromosome && position_within(leg.position, rs, re) && leg.orientation == orientation;
        if let Some(m) = matched[0] {
            if let Some(end) = &m.sv.end {
                if leg_matches(end) {
                    matched[1] = Some(MatchedSv { sv: m.sv, is_start: false });
                    break;
                }
            }
        }
        for v in svs {
            for is_start in [true, false] {
                let leg = if is_start { Some(&v.start) } else { v.end.as_ref() };
                let Some(leg) = leg else { continue };
                if !leg_matches(leg) {
                    continue;
                }
                let replace = match matched[se] {
                    None => true,
                    Some(m) => m.sv.sv_type != preferred && v.sv_type == preferred,
                };
                if replace {
                    matched[se] = Some(MatchedSv { sv: v, is_start });
                }
            }
        }
    }
    matched
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{driver_gene, parse_record};
    use hmf_common::driver_panel::DriverCategory;
    use std::collections::HashSet;

    const GENE_1: &str = "BRCA1";
    const GENE_2: &str = "BRCA2";

    // GermlineReportedTest
    fn germline(gene: &str, biallelic: bool, hotspot: bool, clin_sig: &str, effect: &str, vcn: f64, lps: Option<i32>) -> GermlineVariant {
        let flag = if hotspot { "HOTSPOT" } else { "NEAR_HOTSPOT" };
        let impact = format!("{gene},ENST00000393562,UTR_variant,{effect},c.-275T>G,,false,,{effect},1");
        let lps = lps.map(|l| format!(";LPS={l}")).unwrap_or_default();
        let line = format!("11\t1000\tCOSM123;COSM456\tG\tA\t100\tPASS\tBIALLELIC={biallelic};{flag};IMPACT={impact};CLNSIG={clin_sig};PURPLE_VCN={vcn}{lps}\tGT:AD:DP\t0/1:73,17:91");
        let rec = parse_record(&line);
        let header = VcfHeader::from_lines(Vec::new(), vec!["SAMPLE".to_string()]);
        let v = Variant::new(rec, &header, None, None);
        let filters = v.rec.filters.clone().unwrap_or_default().into_iter().collect();
        GermlineVariant { v, filters }
    }

    fn gene(name: &str, variant: Reporting, hotspot: Reporting) -> DriverGene {
        let mut g = driver_gene(name, DriverCategory::Tsg);
        g.report_germline_variant = variant;
        g.report_germline_hotspot = hotspot;
        g.report_het_deletion = true;
        g.report_germline_disruption = Reporting::Any;
        g.report_germline_deletion = Reporting::Any;
        g
    }

    fn reported(variants: &mut [GermlineVariant], genes: &[DriverGene], somatic: &[&str]) -> Vec<bool> {
        let somatic: HashSet<String> = somatic.iter().map(|s| s.to_string()).collect();
        reported_enrich(variants, genes, &somatic);
        variants.iter().map(|v| v.v.reported()).collect()
    }

    #[test]
    fn report_hotspot() {
        let mut v = vec![germline(GENE_1, false, true, "UNKNOWN", "NONE", 0.5, None)];
        assert_eq!(reported(&mut v, &[gene(GENE_1, Reporting::None, Reporting::Any)], &[]), vec![true]);
    }

    #[test]
    fn report_hotspot_when_multiple_hits() {
        let g = gene(GENE_1, Reporting::WildtypeLost, Reporting::WildtypeLost);
        let mut v = vec![germline(GENE_1, false, true, "UNKNOWN", "NONE", 0.5, None), germline(GENE_1, false, false, "Pathogenic", "NONE", 0.5, None)];
        assert_eq!(reported(&mut v, &[g], &[]), vec![true, true]);
    }

    #[test]
    fn report_hotspot_when_multiple_hits_with_same_lps() {
        let g = gene(GENE_1, Reporting::WildtypeLost, Reporting::WildtypeLost);
        let mut v = vec![germline(GENE_1, false, true, "UNKNOWN", "NONE", 0.5, Some(4)), germline(GENE_1, false, false, "Pathogenic", "NONE", 0.5, Some(4))];
        assert_eq!(reported(&mut v, &[g], &[]), vec![false, false]);
    }

    #[test]
    fn report_unknown_frameshift() {
        let mut v = vec![germline(GENE_1, false, false, "UNKNOWN", "NONSENSE_OR_FRAMESHIFT", 0.5, None)];
        assert_eq!(reported(&mut v, &[gene(GENE_1, Reporting::Any, Reporting::None)], &[]), vec![true]);
    }

    #[test]
    fn report_pathogenic() {
        let mut v = vec![germline(GENE_1, false, false, "Pathogenic", "NONE", 0.5, None)];
        assert_eq!(reported(&mut v, &[gene(GENE_1, Reporting::Any, Reporting::None)], &[]), vec![true]);
    }

    #[test]
    fn report_biallelic_as_multiple_hit() {
        let mut v = vec![germline(GENE_1, true, false, "Pathogenic", "NONE", 0.5, None)];
        assert_eq!(reported(&mut v, &[gene(GENE_1, Reporting::Any, Reporting::WildtypeLost)], &[]), vec![true]);
    }

    #[test]
    fn ignore_single_pathogenic_when_no_multiple_hits() {
        let mut v = vec![germline(GENE_1, false, false, "Pathogenic", "NONE", 0.5, None)];
        assert_eq!(reported(&mut v, &[gene(GENE_1, Reporting::WildtypeLost, Reporting::None)], &[]), vec![false]);
    }

    #[test]
    fn somatic_count_towards_multiple_hits() {
        let mut v = vec![germline(GENE_1, false, false, "Pathogenic", "NONE", 0.5, None)];
        assert_eq!(reported(&mut v, &[gene(GENE_1, Reporting::WildtypeLost, Reporting::None)], &[GENE_1]), vec![true]);
    }

    #[test]
    fn variant_not_lost() {
        let genes = [gene(GENE_1, Reporting::VariantNotLost, Reporting::None), gene(GENE_2, Reporting::VariantNotLost, Reporting::None)];
        let mut v = vec![germline(GENE_1, false, false, "Pathogenic", "MISSENSE", 0.4, None), germline(GENE_2, false, false, "Pathogenic", "MISSENSE", 0.6, None)];
        assert_eq!(reported(&mut v, &genes, &[]), vec![false, true]);
    }

    // GermlineGenotypeEnrichment.status
    #[test]
    fn genotype_status_calls() {
        let d = |alt: i32, total: i32| AllelicDepth { total_read_count: total, allele_read_count: alt };
        assert!(genotype_status(d(30, 30)) == GermlineGenotypeStatus::HomAlt);
        assert!(genotype_status(d(15, 30)) == GermlineGenotypeStatus::Het);
        assert!(genotype_status(d(2, 60)) == GermlineGenotypeStatus::LowVaf);
        assert!(genotype_status(d(29, 30)) == GermlineGenotypeStatus::HomAlt);
    }
}
