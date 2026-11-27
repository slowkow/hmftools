//! ChromosomeTask: the variants of one chromosome, classified, annotated, filtered and encoded.

use anyhow::{Context, Result};

use hmf_common::vcf::{VcfHeader, VcfRecord};

use crate::annotation::blacklist::Blacklistings;
use crate::annotation::clinvar::ClinvarData;
use crate::annotation::gnomad::GnomadData;
use crate::annotation::mappability::MappabilityData;
use crate::annotation::pon::{PonData, PonFilters};
use crate::annotation::reportability::Reportability;
use crate::config::PaveConfig;
use crate::constants::*;
use crate::gene_cache::GeneDataCache;
use crate::impact::builder::create_variant_impact;
use crate::impact::phased::{process_phased_variants, PhasedVariantClassifier};
use crate::impact::utils::{create_right_aligned_variant, find_variant_impacts};
use crate::impact::RefGenome;
use crate::output;
use crate::variant::{MultiPonStatus, VariantData, VariantTier, NO_LOCAL_PHASE_SET};

/// What every chromosome task shares (ReferenceData and the config)
pub struct Shared<'a> {
    pub config: &'a PaveConfig,
    pub gene_cache: &'static GeneDataCache,
    pub ref_genome: &'a dyn RefGenome,
    pub input_header: &'a VcfHeader,
    pub output_header: &'a VcfHeader,
    /// the standard PON when it is enabled
    pub pon: Option<&'a PonData>,
    pub pon_filters: &'a PonFilters,
    pub artefact_pon: Option<&'a PonData>,
    pub gnomad: Option<&'a GnomadData>,
    pub gnomad_apply_filter: bool,
    pub mappability: Option<&'a MappabilityData>,
    pub clinvar: Option<&'a ClinvarData>,
    pub blacklist: &'a Blacklistings,
    pub reportability: &'a Reportability,
}

/// One output record: its position and end (for the index) and its line
pub struct OutputRecord {
    pub chrom: String,
    pub position: i32,
    pub end: i32,
    pub line: String,
}

#[derive(Default)]
pub struct ChromosomeOutput {
    pub records: Vec<OutputRecord>,
    pub transcript_rows: Vec<String>,
}

pub fn run_chromosome(shared: &Shared, chromosome: &str, lines: &[String]) -> Result<ChromosomeOutput> {
    let mut out = ChromosomeOutput::default();
    let mut indexing = shared.gene_cache.create_indexing(chromosome);
    let mut phased = PhasedVariantClassifier::default();

    for line in lines {
        let rec = VcfRecord::parse(line, shared.input_header)?;
        let mut variant = VariantData::from_context(rec);
        let classified = (|| -> Result<()> {
            let ra = create_right_aligned_variant(&variant, shared.ref_genome)?;
            variant.set_realigned_variant(ra);
            find_variant_impacts(&mut variant, shared.gene_cache, Some(&mut indexing), shared.ref_genome)?;
            Ok(())
        })();
        if let Err(e) = classified {
            // PA4: Java ends the run on an error classifying a variant (complex INDELs, whose ref and alt
            // are both several bases long, can index past their codons). Corrected: the variant is kept
            // without transcript impacts.
            if crate::compat::java("PA4") {
                return Err(e).with_context(|| format!("error processing var({variant})"));
            }
            crate::warn!("var({}) not classified: {:#}", variant, e);
            variant.impacts_mut().clear();
            variant.set_realigned_variant_none();
        }

        let lps = variant.local_phase_set();
        let pending = if variant.has_local_phase_set() {
            // all phased variants are held until their set is complete
            phased.add_variant(variant);
            None
        } else {
            Some(variant)
        };
        for v in process_phased_variants(&mut phased, lps, shared.ref_genome) {
            process_variant(shared, v, &mut out)?;
        }
        if let Some(v) = pending {
            process_variant(shared, v, &mut out)?;
        }
    }
    for v in process_phased_variants(&mut phased, NO_LOCAL_PHASE_SET, shared.ref_genome) {
        process_variant(shared, v, &mut out)?;
    }
    Ok(out)
}

/// processVariant(VariantData): the summary impact, annotations and filters, then the output record
fn process_variant(shared: &Shared, mut variant: VariantData, out: &mut ChromosomeOutput) -> Result<()> {
    // None if no transcript has an impact
    let vi = create_variant_impact(&variant, shared.gene_cache);

    annotate_and_filter(shared, &mut variant).with_context(|| format!("error processing var({variant})"))?;

    if shared.config.set_reportable {
        shared.reportability.set_reportability(&mut variant, vi.as_ref()).with_context(|| format!("error processing var({variant})"))?;
    }

    if shared.config.write_pass_only && !variant.filters.is_empty() {
        return Ok(());
    }

    let mut rec = variant.take_context().unwrap();
    // the transcript file shows the record's filters as read (VariantContextBuilder copies them)
    let input_filters = (rec.filters.clone(), rec.filters_applied);
    output::finish_record(&mut rec, &variant, vi.as_ref(), shared.config.write_detailed);
    let line = rec.encode(shared.output_header)?;
    out.records.push(OutputRecord { chrom: rec.chrom.clone(), position: rec.pos, end: rec.end(), line });

    if shared.config.write_transcript_file {
        // the record as the transcript file reads it
        let mut v = variant;
        rec.filters = input_filters.0;
        rec.filters_applied = input_filters.1;
        restore_context(&mut v, rec);
        for (gene_name, _) in v.impacts().iter() {
            output::transcript_rows(&v, gene_name, &shared.config.sample_id, &mut out.transcript_rows)?;
        }
    }
    Ok(())
}

fn restore_context(v: &mut VariantData, rec: VcfRecord) {
    v.set_context(rec);
}

fn annotate_and_filter(shared: &Shared, variant: &mut VariantData) -> Result<()> {
    let chr = variant.chromosome.clone();
    let (pos, ref_, alt) = (variant.position, variant.ref_.clone(), variant.alt.clone());

    if let Some(m) = shared.mappability.filter(|m| m.has_chromosome(&chr)) {
        if let Some(v) = m.value(&chr, pos) {
            if !variant.context().has_attr(MAPPABILITY) {
                variant.context_mut().put_text(MAPPABILITY, hmf_common::fmt::vcf_double(v));
            }
        }
    }

    if let Some(c) = shared.clinvar.filter(|c| c.has_chromosome(&chr)) {
        c.annotate(variant.context_mut(), &chr, pos, &ref_, &alt)?;
    }

    shared.blacklist.annotate(variant)?;

    match shared.gnomad.filter(|g| g.has_chromosome(&chr)) {
        Some(g) => {
            if let Some(f) = g.frequency(&chr, variant.is_mnv(), &ref_, &alt, pos) {
                variant.set_gnomad_frequency(Some(f));
            }
        }
        None => {
            // GnomadAnnotation.annotateFromContext
            if variant.context().has_attr(GNOMAD_FREQ) {
                let f = variant.context().attr_double(GNOMAD_FREQ, -1.0);
                variant.set_gnomad_frequency(Some(f));
            }
        }
    }

    match shared.pon.filter(|p| p.has_chromosome(&chr)) {
        Some(p) => {
            if let Some(d) = p.get(&chr, pos, &ref_, &alt) {
                let (s, m, mean, multi) = (d.samples, d.max_sample_reads, d.mean_read_count(), d.multi_pon);
                variant.set_pon_frequency(s, m, mean, multi);
            }
        }
        None => {
            // PonAnnotation.annotateFromContext
            let ctx = variant.context();
            if ctx.has_attr(PON_COUNT) {
                let status = MultiPonStatus::value_of(&ctx.attr_string(PON_MULTI_STATUS, "BASE")).unwrap_or_else(|e| panic!("{e}"));
                let (s, m, mean) = (ctx.attr_int(PON_COUNT, 0), ctx.attr_int(PON_MAX, 0), ctx.attr_int(PON_AVG_READS, 0));
                variant.set_pon_frequency(s, m, mean, status);
            }
        }
    }

    apply_filters(variant, shared);
    Ok(())
}

/// ChromosomeTask.applyFilters
pub fn apply_filters(variant: &mut VariantData, shared: &Shared) {
    let artefact = shared.artefact_pon.filter(|p| p.has_chromosome(&variant.chromosome));
    apply_filters_with(variant, &shared.config.sample_id, shared.pon_filters, artefact, shared.gnomad_apply_filter, shared.config.is_sbx(), shared.config.is_ultima())
}

fn below_pon_threshold(pon_sample_count: i32, repeat_count: i32) -> bool {
    pon_sample_count < PON_SAMPLE_COUNT_THRESHOLD && repeat_count < PON_REPEAT_COUNT_THRESHOLD
}

/// PathogenicSummaryFactory.fromContext(...).Status.isPathogenic()
pub fn is_clinvar_pathogenic(rec: &VcfRecord) -> bool {
    let cln_sig = rec.attr_string_list(CLNSIG).join(",");
    let cln_sig_conf = rec.attr_string_list(CLNSIGCONF).join(",");
    let blacklisted = rec.attr_bool(BLACKLIST_BED_FLAG, false) || rec.attr_bool(BLACKLIST_VCF_FLAG, false);
    if blacklisted {
        return false; // BENIGN_BLACKLIST
    }
    pathogenicity_from_clinvar(&cln_sig, &cln_sig_conf).is_pathogenic()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Pathogenicity {
    Pathogenic,
    LikelyPathogenic,
    Benign,
    LikelyBenign,
    Conflicting,
    Unknown,
}

impl Pathogenicity {
    pub fn is_pathogenic(self) -> bool {
        matches!(self, Pathogenicity::Pathogenic | Pathogenicity::LikelyPathogenic)
    }
}

/// Pathogenicity.fromClinvarAnnotation
pub fn pathogenicity_from_clinvar(cln_sig: &str, cln_sig_conf: &str) -> Pathogenicity {
    let has_pathogenic = |s: &str| s.contains("Pathogenic");
    let has_likely_pathogenic = |s: &str| s.contains("Likely_pathogenic");
    let has_benign = |s: &str| s.contains("Benign");
    let has_likely_benign = |s: &str| s.contains("Likely_benign");
    use Pathogenicity::*;
    if cln_sig.is_empty() {
        return Unknown;
    }
    if has_pathogenic(cln_sig) {
        return Pathogenic;
    }
    if has_likely_pathogenic(cln_sig) {
        return LikelyPathogenic;
    }
    if has_benign(cln_sig) {
        return Benign;
    }
    if has_likely_benign(cln_sig) {
        return LikelyBenign;
    }
    if cln_sig.contains("Conflicting") {
        let c = cln_sig_conf;
        if has_pathogenic(c) && !has_benign(c) && !has_likely_benign(c) {
            return Pathogenic;
        }
        if has_likely_pathogenic(c) && !has_pathogenic(c) && !has_benign(c) && !has_likely_benign(c) {
            return LikelyPathogenic;
        }
        if has_benign(c) && !has_pathogenic(c) && !has_likely_pathogenic(c) {
            return Benign;
        }
        if has_likely_benign(c) && !has_pathogenic(c) && !has_likely_pathogenic(c) && !has_benign(c) {
            return LikelyBenign;
        }
        if (has_pathogenic(c) || has_likely_pathogenic(c)) && (has_benign(c) || has_likely_benign(c)) {
            return Conflicting;
        }
    }
    Unknown
}

/// applyFilters(variant, sampleId, standardPon, artefactsPon, applyGnomadFilter)
pub fn apply_filters_with(
    variant: &mut VariantData, sample_id: &str, pon_filters: &PonFilters, artefact_pon: Option<&PonData>, apply_gnomad_filter: bool, is_sbx: bool,
    is_ultima: bool,
) {
    variant.filters.clear();

    let is_hotspot = variant.tier() == VariantTier::Hotspot;
    let clinvar_pathogenic = is_clinvar_pathogenic(variant.context());
    let hotspot_or_pathogenic = is_hotspot || clinvar_pathogenic;
    let repeat_count = variant.repeat_count();

    // the WGS PON
    let mut pon_filter = pon_filters.filter_on_tier_criteria(variant.tier(), variant.pon_sample_count, variant.pon_max_read_count);

    if hotspot_or_pathogenic {
        if variant.pon_sample_count == 0 || below_pon_threshold(variant.pon_sample_count, repeat_count) {
            pon_filter = false;
        } else {
            let variant_vaf = variant.sample_vaf(sample_id);
            let repeat_base_length = variant.repeat_sequence().len() as i32 * repeat_count;
            let vaf_limit = PON_VAF_THRESHOLD.max(0.01 * repeat_base_length as f64);
            if variant.pon_mean_read_count < PON_MEAN_READ_THRESHOLD && variant_vaf > vaf_limit {
                pon_filter = false;
            }
        }
    }

    if pon_filter && variant.is_indel() && repeat_count >= PON_INDEL_ARTEFACT_REPEAT_COUNT && (is_sbx || is_ultima) && variant.multi_pon_status == MultiPonStatus::Artefact {
        // the VAF-based INDEL filtering applies only to variants in the sequencing technology's PON alone,
        // MSI INDELs with a repeat count of 7 or more: kept if PON_MAX < min(60 * (AF - 0.1), 18) (Ultima)
        // or min(40 * (AF - 0.1), 8) (SBX)
        let variant_vaf = variant.sample_vaf(sample_id);
        let (min_threshold, factor) = if is_ultima {
            (PON_INDEL_ARTEFACT_ULTIMA_MIN_THRESHOLD, PON_INDEL_ARTEFACT_ULTIMA_FACTOR)
        } else {
            (PON_INDEL_ARTEFACT_SBX_MIN_THRESHOLD, PON_INDEL_ARTEFACT_SBX_FACTOR)
        };
        let pon_threshold = (factor as f64 * (variant_vaf - PON_INDEL_ARTEFACT_VAF_REDUCTION)).min(min_threshold as f64);
        if (variant.pon_max_read_count as f64) < pon_threshold {
            pon_filter = false;
        }
    }

    if pon_filter {
        variant.add_filter(PON_FILTER);
    }

    if let Some(artefacts) = artefact_pon {
        if let Some(d) = artefacts.get(&variant.chromosome, variant.position, &variant.ref_, &variant.alt) {
            if !hotspot_or_pathogenic || !below_pon_threshold(d.samples, repeat_count) {
                variant.add_filter(PON_ARTEFACT_FILTER);
            }
        }
    }

    if apply_gnomad_filter {
        let threshold = if hotspot_or_pathogenic { GNMOAD_FILTER_HOTSPOT_PATHOGENIC_THRESHOLD } else { GNMOAD_FILTER_THRESHOLD };
        if let Some(f) = variant.gnomad_frequency() {
            if f >= threshold {
                variant.add_filter(PON_GNOMAD_FILTER);
            }
        }
    }
}
