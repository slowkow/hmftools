//! A small variant read from a VCF with the accessors of VariantContextDecorator and SomaticVariant:
//! type, filter, tier, allelic depths, the PAVE impact, hotspot and ClinVar status.

use hmf_common::variant_context::{Genotype, VcfHeader, VcfRecord, VariantType as VcfType};

use crate::types::*;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AllelicDepth {
    pub total_read_count: i32,
    pub allele_read_count: i32,
}

impl AllelicDepth {
    pub fn allele_frequency(&self) -> f64 {
        if self.total_read_count > 0 { self.allele_read_count as f64 / self.total_read_count as f64 } else { 0.0 }
    }

    /// AllelicDepth.fromGenotype (Java requires AD)
    pub fn from_genotype(g: &Genotype) -> AllelicDepth {
        let ad = g.ad.as_ref().unwrap_or_else(|| panic!("IllegalArgumentException: genotype {} without AD", g.sample));
        let sum: i32 = ad.iter().sum();
        let total = match g.dp {
            Some(dp) => dp.max(sum),
            None => sum,
        };
        AllelicDepth { total_read_count: total, allele_read_count: ad[1] }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct VariantImpact {
    pub gene_name: String,
    pub canonical_transcript: String,
    pub canonical_effect: String,
    pub canonical_coding_effect: CodingEffect,
    pub canonical_hgvs_coding: String,
    pub canonical_hgvs_protein: String,
    pub canonical_splice_region: bool,
    pub other_reportable_effects: String,
    pub worst_coding_effect: CodingEffect,
    pub genes_affected: i32,
}

impl VariantImpact {
    /// VariantImpactSerialiser.fromVariantContext
    pub fn from_record(rec: &VcfRecord) -> VariantImpact {
        let v = rec.info.as_string_list("IMPACT", "");
        if v.len() != 10 {
            return VariantImpact {
                gene_name: String::new(),
                canonical_transcript: String::new(),
                canonical_effect: String::new(),
                canonical_coding_effect: CodingEffect::Undefined,
                canonical_hgvs_coding: String::new(),
                canonical_hgvs_protein: String::new(),
                canonical_splice_region: false,
                other_reportable_effects: String::new(),
                worst_coding_effect: CodingEffect::Undefined,
                genes_affected: 0,
            };
        }
        let effect = |s: &str| CodingEffect::from_name(s).unwrap_or_else(|| panic!("No enum constant CodingEffect.{s}"));
        VariantImpact {
            gene_name: v[0].clone(),
            canonical_transcript: v[1].clone(),
            canonical_effect: v[2].clone(),
            canonical_coding_effect: effect(&v[3]),
            canonical_splice_region: v[4].eq_ignore_ascii_case("true"),
            canonical_hgvs_coding: v[5].clone(),
            canonical_hgvs_protein: v[6].clone(),
            other_reportable_effects: v[7].clone(),
            worst_coding_effect: effect(&v[8]),
            genes_affected: hmf_common::java_text::parse_int(&v[9]).unwrap_or_else(|| panic!("NumberFormatException: {}", v[9])),
        }
    }
}

/// VariantType.type(context)
pub fn variant_type(rec: &VcfRecord) -> VariantType {
    match rec.variant_type() {
        VcfType::Mnp => VariantType::Mnp,
        VcfType::Snp => VariantType::Snp,
        VcfType::Indel => VariantType::Indel,
        _ => VariantType::Undefined,
    }
}

/// VariantTier.fromContext
pub fn variant_tier(rec: &VcfRecord) -> VariantTier {
    match rec.info.as_string("TIER", "UNKNOWN").as_str() {
        "HOTSPOT" => VariantTier::Hotspot,
        "PANEL" => VariantTier::Panel,
        "HIGH_CONFIDENCE" => VariantTier::HighConfidence,
        "LOW_CONFIDENCE" => VariantTier::LowConfidence,
        _ => VariantTier::Unknown,
    }
}

/// Pathogenicity.fromClinvarAnnotation
pub fn pathogenicity_from_clinvar(sig: &str, conf: &str) -> Pathogenicity {
    let p = |s: &str| s.contains("Pathogenic");
    let lp = |s: &str| s.contains("Likely_pathogenic");
    let b = |s: &str| s.contains("Benign");
    let lb = |s: &str| s.contains("Likely_benign");
    if sig.is_empty() {
        return Pathogenicity::Unknown;
    }
    if p(sig) {
        return Pathogenicity::Pathogenic;
    }
    if lp(sig) {
        return Pathogenicity::LikelyPathogenic;
    }
    if b(sig) {
        return Pathogenicity::Benign;
    }
    if lb(sig) {
        return Pathogenicity::LikelyBenign;
    }
    if sig.contains("Conflicting") {
        if p(conf) && !b(conf) && !lb(conf) {
            return Pathogenicity::Pathogenic;
        }
        if lp(conf) && !p(conf) && !b(conf) && !lb(conf) {
            return Pathogenicity::LikelyPathogenic;
        }
        if b(conf) && !p(conf) && !lp(conf) {
            return Pathogenicity::Benign;
        }
        if lb(conf) && !p(conf) && !lp(conf) && !b(conf) {
            return Pathogenicity::LikelyBenign;
        }
        if (p(conf) || lp(conf)) && (b(conf) || lb(conf)) {
            return Pathogenicity::Conflicting;
        }
    }
    Pathogenicity::Unknown
}

/// PathogenicSummaryFactory.fromContext: (clinvar info, status)
pub fn pathogenic_summary(rec: &VcfRecord) -> (String, Pathogenicity) {
    let sig = rec.info.as_string_list("CLNSIG", "").join(",");
    let conf = rec.info.as_string_list("CLNSIGCONF", "").join(",");
    let info = if !conf.is_empty() { conf.clone() } else { sig.clone() };
    let blacklisted = rec.info.as_bool("BLACKLIST_BED", false) || rec.info.as_bool("BLACKLIST_VCF", false);
    let status = if blacklisted { Pathogenicity::BenignBlacklist } else { pathogenicity_from_clinvar(&sig, &conf) };
    (info, status)
}

/// DriverImpact.select
pub fn driver_impact(t: VariantType, effect: CodingEffect) -> DriverImpact {
    let indel = t == VariantType::Indel;
    if indel && effect == CodingEffect::NonsenseOrFrameshift {
        DriverImpact::Frameshift
    } else if !indel && effect == CodingEffect::NonsenseOrFrameshift {
        DriverImpact::Nonsense
    } else if !indel && effect == CodingEffect::Missense {
        DriverImpact::Missense
    } else if effect == CodingEffect::Splice {
        DriverImpact::Splice
    } else if indel && effect == CodingEffect::Missense {
        DriverImpact::Inframe
    } else {
        DriverImpact::Unknown
    }
}

/// SomaticVariant / GermlineVariant: the record and what VariantContextDecorator caches
#[derive(Clone, Debug)]
pub struct Variant {
    pub rec: VcfRecord,
    pub vtype: VariantType,
    pub filter: String,
    pub tier: VariantTier,
    pub impact: VariantImpact,
    pub tumor_depth: Option<AllelicDepth>,
    pub reference_depth: Option<AllelicDepth>,
}

impl Variant {
    /// SomaticVariant(context, sampleId, referenceId): the genotypes are decoded for the allelic depths
    pub fn new(mut rec: VcfRecord, header: &VcfHeader, tumor: Option<&str>, reference: Option<&str>) -> Variant {
        let depth = |rec: &mut VcfRecord, s: Option<&str>| s.map(|s| rec.genotype(header, s).map(|g| AllelicDepth::from_genotype(g)).unwrap_or_default());
        let tumor_depth = depth(&mut rec, tumor);
        let reference_depth = depth(&mut rec, reference);
        Variant {
            vtype: variant_type(&rec),
            filter: rec.display_filter(),
            tier: variant_tier(&rec),
            impact: VariantImpact::from_record(&rec),
            tumor_depth,
            reference_depth,
            rec,
        }
    }

    /// recomputes the decorator after the record changed (SomaticVariant.setContext)
    pub fn refresh(&mut self) {
        self.vtype = variant_type(&self.rec);
        self.filter = self.rec.display_filter();
        self.tier = variant_tier(&self.rec);
    }

    pub fn chromosome(&self) -> &str {
        &self.rec.chrom
    }
    pub fn position(&self) -> i32 {
        self.rec.pos
    }
    pub fn is_pass(&self) -> bool {
        self.filter == "PASS"
    }
    pub fn is_filtered(&self) -> bool {
        !self.is_pass()
    }
    pub fn is_hotspot(&self) -> bool {
        self.rec.info.has("HOTSPOT")
    }
    /// HotspotType.fromVariant == HOTSPOT
    pub fn hotspot_flag(&self) -> bool {
        self.rec.info.as_bool("HOTSPOT", false)
    }
    pub fn is_hotspot_type(&self) -> bool {
        self.is_hotspot() || (self.tier == VariantTier::Panel && pathogenic_summary(&self.rec).1.is_pathogenic())
    }
    pub fn copy_number(&self) -> f64 {
        self.rec.info.as_double("PURPLE_VCN", 0.0)
    }
    pub fn biallelic(&self) -> bool {
        self.rec.info.as_bool("BIALLELIC", false)
    }
    pub fn gene(&self) -> &str {
        &self.impact.gene_name
    }
    pub fn allele_frequency(&self) -> f64 {
        self.tumor_depth.map(|d| d.allele_frequency()).unwrap_or(0.0)
    }
    pub fn total_read_count(&self) -> i32 {
        self.tumor_depth.map(|d| d.total_read_count).unwrap_or(0)
    }
    pub fn allele_read_count(&self) -> i32 {
        self.tumor_depth.map(|d| d.allele_read_count).unwrap_or(0)
    }
    pub fn reference_allele_read_count(&self) -> i32 {
        self.reference_depth.map(|d| d.allele_read_count).unwrap_or(0)
    }
    pub fn repeat_count(&self) -> i32 {
        self.rec.info.as_int("REP_C", 0)
    }
    pub fn ref_bases(&self) -> &str {
        &self.rec.ref_allele
    }
    /// VariantContextDecorator.alt: the alternates joined by ','
    pub fn alt(&self) -> String {
        self.rec.alts.join(",")
    }
    pub fn driver_impact(&self) -> DriverImpact {
        driver_impact(self.vtype, self.impact.canonical_coding_effect)
    }
    pub fn reported(&self) -> bool {
        self.rec.info.as_bool("REPORTED", false)
    }
    /// SomaticVariant.toString
    pub fn describe(&self) -> String {
        format!(
            "{} {}:{} {}>{} filter({}) tier({}) codingEffect({})",
            self.vtype,
            self.chromosome(),
            self.position(),
            self.ref_bases(),
            self.alt(),
            self.filter,
            self.tier,
            self.impact.canonical_coding_effect
        )
    }
}
