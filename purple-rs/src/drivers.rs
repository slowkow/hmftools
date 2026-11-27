//! The drivers package: DriverCatalog and its file, ReportablePredicate, the somatic variant drivers
//! (OncoDrivers, TsgDrivers with the DNDS likelihoods), AmplificationDrivers and DeletionDrivers.

use std::collections::HashMap;

use hmf_common::commons_math::Poisson;
use hmf_common::driver_panel::{DriverCategory as PanelCategory, DriverGene};
use hmf_common::java_hash::JavaHashMap;

use crate::gene::GeneCopyNumber;
use crate::reference::{DndsGeneLikelihood, DndsImpactLikelihood, DriverGeneResource};
use crate::types::*;
use crate::variant::{driver_impact, Variant, VariantImpact};

#[derive(Clone, Debug, PartialEq)]
pub struct DriverCatalog {
    pub chromosome: String,
    pub chromosome_band: String,
    pub gene: String,
    pub transcript: String,
    pub is_canonical: bool,
    pub category: DriverCategory,
    pub driver: DriverType,
    pub likelihood_method: LikelihoodMethod,
    pub reported_status: ReportedStatus,
    pub driver_likelihood: f64,
    pub missense: i32,
    pub nonsense: i32,
    pub splice: i32,
    pub inframe: i32,
    pub frameshift: i32,
    pub biallelic: bool,
    pub min_copy_number: f64,
    pub max_copy_number: f64,
}

/// DriverCatalogFile.write
pub fn write_driver_catalog(path: &str, catalog: &[DriverCatalog]) -> anyhow::Result<()> {
    let f = |v: f64| hmf_common::fmt::decimal_format_min(v, 4, 4);
    let mut out = String::from("chromosome\tchromosomeBand\tgene\ttranscript\tisCanonical\tdriver\tcategory\tlikelihoodMethod\treportedStatus\tdriverLikelihood\tmissense\tnonsense\tsplice\tinframe\tframeshift\tbiallelic\tminCopyNumber\tmaxCopyNumber\n");
    for d in catalog {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            d.chromosome,
            d.chromosome_band,
            d.gene,
            d.transcript,
            d.is_canonical,
            d.driver,
            d.category,
            d.likelihood_method,
            d.reported_status,
            f(d.driver_likelihood),
            d.missense,
            d.nonsense,
            d.splice,
            d.inframe,
            d.frameshift,
            d.biallelic,
            f(d.min_copy_number),
            f(d.max_copy_number)
        ));
    }
    std::fs::write(path, out).map_err(|e| anyhow::anyhow!("{path}: {e}"))
}

pub fn category_of(c: PanelCategory) -> DriverCategory {
    match c {
        PanelCategory::Onco => DriverCategory::Onco,
        PanelCategory::Tsg => DriverCategory::Tsg,
    }
}

/// a driver and what it came from, for the Circos driver files (DriverSourceData)
#[derive(Clone, Debug)]
pub struct DriverSource {
    pub driver: DriverCatalog,
    pub chromosome: String,
    /// the variant's position, or the gene copy number's start
    pub position: i32,
    /// the gene copy number's end (None for a variant source)
    pub gene_end: Option<i32>,
    pub gene: String,
    pub is_indel: bool,
    /// PURPLE_AF of a variant source
    pub adjusted_vaf: f64,
    /// SOM_LH is LOW
    pub likely_germline: bool,
}

impl DriverSource {
    pub fn from_variant(driver: DriverCatalog, v: &Variant) -> DriverSource {
        DriverSource {
            driver,
            chromosome: v.chromosome().to_string(),
            position: v.position(),
            gene_end: None,
            gene: v.gene().to_string(),
            is_indel: v.vtype == VariantType::Indel,
            adjusted_vaf: v.rec.info.as_double("PURPLE_AF", 0.0),
            likely_germline: v.rec.info.as_string("SOM_LH", "UNKNOWN") == "LOW",
        }
    }
    pub fn from_gene(driver: DriverCatalog, g: &GeneCopyNumber) -> DriverSource {
        DriverSource {
            driver,
            chromosome: g.chromosome.clone(),
            position: g.start,
            gene_end: Some(g.end),
            gene: g.gene_name.clone(),
            is_indel: false,
            adjusted_vaf: 0.0,
            likely_germline: false,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// VariantEffect and ReportablePredicate

/// CodingEffect.effect(VariantEffect.fromEffect(name))
fn coding_effect_of_effect(name: &str) -> CodingEffect {
    match name {
        "stop_gained" | "stop_lost" | "start_lost" | "frameshift_variant" => CodingEffect::NonsenseOrFrameshift,
        "splice_acceptor_variant" | "splice_donor_variant" => CodingEffect::Splice,
        "inframe_insertion" | "inframe_deletion" | "missense_variant" | "phased_inframe_insertion" | "phased_inframe_deletion" | "phased_missense" => {
            CodingEffect::Missense
        }
        "synonymous_variant" | "phased_synonymous" => CodingEffect::Synonymous,
        _ => CodingEffect::None,
    }
}

/// AltTranscriptReportableInfo.parse of each item of OtherReportableEffects
#[derive(Clone, Debug)]
pub struct AltTranscriptInfo {
    pub gene_name: String,
    pub trans_name: String,
    pub hgvs_coding: String,
    pub hgvs_protein: String,
    pub effects: String,
    pub effect: CodingEffect,
}

pub fn parse_alt_transcript_info(other: &str) -> Vec<AltTranscriptInfo> {
    let mut out = Vec::new();
    for item in other.split("--") {
        let v: Vec<&str> = item.split('|').collect();
        if v.len() < 5 {
            continue;
        }
        let mut i = 0;
        let mut gene = String::new();
        if v.len() == 6 {
            gene = v[0].to_string();
            i = 1;
        }
        let effect = CodingEffect::from_name(v[i + 4]).unwrap_or_else(|| panic!("No enum constant CodingEffect.{}", v[i + 4]));
        out.push(AltTranscriptInfo {
            gene_name: gene,
            trans_name: v[i].to_string(),
            hgvs_coding: v[i + 1].to_string(),
            hgvs_protein: v[i + 2].to_string(),
            effects: v[i + 3].to_string(),
            effect,
        });
    }
    out
}

pub struct ReportablePredicate<'a> {
    genes: HashMap<String, &'a DriverGene>,
}

impl<'a> ReportablePredicate<'a> {
    pub fn new(category: DriverCategory, genes: &'a [DriverGene]) -> ReportablePredicate<'a> {
        let mut map = HashMap::new();
        for g in genes.iter().filter(|g| category_of(g.likelihood_type) == category && g.report_somatic()) {
            if map.insert(g.gene.clone(), g).is_some() {
                panic!("IllegalStateException: Duplicate key {}", g.gene);
            }
        }
        ReportablePredicate { genes: map }
    }
    pub fn is_driver_gene(&self, gene: &str) -> bool {
        self.genes.contains_key(gene)
    }
    pub fn is_candidate_reportable(&self, impact: &VariantImpact, hotspot: bool) -> bool {
        if !self.genes.contains_key(&impact.gene_name) {
            return false;
        }
        hotspot || impact.canonical_coding_effect.has_protein_impact() || impact.canonical_splice_region
    }
    pub fn is_reportable_effects(&self, gene: &str, hotspot: bool, effects: &str) -> bool {
        let Some(g) = self.genes.get(gene) else { return false };
        if hotspot && g.report_somatic_hotspot {
            return true;
        }
        let effects: Vec<CodingEffect> = effects.split('&').map(coding_effect_of_effect).collect();
        if g.report_missense_and_inframe && effects.contains(&CodingEffect::Missense) {
            return true;
        }
        if g.report_splice && effects.contains(&CodingEffect::Splice) {
            return true;
        }
        if g.report_nonsense_and_frameshift && effects.contains(&CodingEffect::NonsenseOrFrameshift) {
            return true;
        }
        false
    }
    pub fn is_reportable(&self, impact: &VariantImpact, hotspot: bool) -> bool {
        if self.is_reportable_effects(&impact.gene_name, hotspot, &impact.canonical_effect) {
            return true;
        }
        if impact.other_reportable_effects.is_empty() {
            return false;
        }
        parse_alt_transcript_info(&impact.other_reportable_effects).iter().any(|x| self.is_reportable_effects(&impact.gene_name, hotspot, &x.effects))
    }
}

fn has_coding_effect(t: VariantType, e: CodingEffect) -> bool {
    driver_impact(t, e) != DriverImpact::Unknown
}

/// SomaticVariantDrivers.hasTranscriptCodingEffect
pub fn has_transcript_coding_effect(impact: &VariantImpact, t: VariantType, transcript: &str) -> bool {
    if impact.canonical_transcript == transcript {
        return has_coding_effect(t, impact.canonical_coding_effect);
    }
    if !impact.other_reportable_effects.is_empty() {
        return parse_alt_transcript_info(&impact.other_reportable_effects).iter().filter(|x| x.trans_name == transcript).any(|x| has_coding_effect(t, x.effect));
    }
    false
}

/// SomaticVariantDrivers.getWorstReportableCodingEffect
pub fn worst_reportable_coding_effect(impact: &VariantImpact) -> CodingEffect {
    let mut top = impact.canonical_coding_effect;
    if !impact.other_reportable_effects.is_empty() {
        for a in parse_alt_transcript_info(&impact.other_reportable_effects) {
            if a.effect.ordinal() < top.ordinal() {
                top = a.effect;
            }
        }
    }
    top
}

/// SomaticVariantDrivers.addReportableTranscriptList: the REPORTABLE_TRANSCRIPTS value, if any
pub fn reportable_transcript_list(t: VariantType, impact: &VariantImpact) -> Option<String> {
    if impact.other_reportable_effects.is_empty() {
        return None;
    }
    let mut out: Vec<String> = Vec::new();
    if has_coding_effect(t, impact.canonical_coding_effect) {
        out.push(impact.canonical_transcript.clone());
    }
    for a in parse_alt_transcript_info(&impact.other_reportable_effects) {
        if has_coding_effect(t, a.effect) {
            out.push(a.trans_name.clone());
        }
    }
    Some(out.join("|"))
}

// ---------------------------------------------------------------------------------------------
// DndsCalculator

fn probability_driver_same_impact(count: i32, sample_count: i64, l: &DndsImpactLikelihood) -> f64 {
    let lambda = sample_count as f64 * l.passengers_per_mutation;
    if less_or_equal(lambda, 0.0) {
        return 0.0;
    }
    let p_non_driver = 1.0 - Poisson::new(lambda).cumulative_probability(count);
    l.drivers_per_sample / (l.drivers_per_sample + p_non_driver * (1.0 - l.drivers_per_sample))
}

pub fn probability_driver_variant(sample_count: i64, l: &DndsImpactLikelihood) -> f64 {
    probability_driver_same_impact(0, sample_count, l)
}

pub fn probability_driver_variant2(first_count: i64, second_count: i64, first: &DndsImpactLikelihood, second: &DndsImpactLikelihood) -> f64 {
    if first == second {
        return probability_driver_same_impact(1, first_count, first);
    }
    let lambda1 = first_count as f64 * first.passengers_per_mutation;
    let lambda2 = second_count as f64 * second.passengers_per_mutation;
    if less_or_equal(lambda1, 0.0) || less_or_equal(lambda2, 0.0) {
        return jmax(probability_driver_variant(first_count, first), probability_driver_variant(second_count, second));
    }
    let p_driver = jmax(first.drivers_per_sample, second.drivers_per_sample);
    let p1 = 1.0 - Poisson::new(lambda1).cumulative_probability(0);
    let p2 = 1.0 - Poisson::new(lambda2).cumulative_probability(0);
    let p = p1 * p2;
    p_driver / (p_driver + p * (1.0 - p_driver))
}

// ---------------------------------------------------------------------------------------------
// somatic variant drivers

const NO_GENE_DNDS: DndsGeneLikelihood = DndsGeneLikelihood {
    gene: String::new(),
    missense: DndsImpactLikelihood { drivers_per_sample: 0.0, passengers_per_mutation: 0.0 },
    nonsense: DndsImpactLikelihood { drivers_per_sample: 0.0, passengers_per_mutation: 0.0 },
    splice: DndsImpactLikelihood { drivers_per_sample: 0.0, passengers_per_mutation: 0.0 },
    indel: DndsImpactLikelihood { drivers_per_sample: 0.0, passengers_per_mutation: 0.0 },
};

/// SomaticVariantDriverFinder (one per category); variants are indices into the somatic variant list
pub struct DriverFinder<'a> {
    pub category: DriverCategory,
    pub predicate: ReportablePredicate<'a>,
    likelihoods: &'a HashMap<String, DndsGeneLikelihood>,
    pub driver_variants: Vec<usize>,
}

#[derive(Default, Clone)]
pub struct VariantTypeCounts {
    pub snp: Option<i32>,
    pub indel: Option<i32>,
    pub mnp: Option<i32>,
    pub undefined: Option<i32>,
}

impl VariantTypeCounts {
    pub fn get(&self, t: VariantType) -> i32 {
        match t {
            VariantType::Snp => self.snp,
            VariantType::Indel => self.indel,
            VariantType::Mnp => self.mnp,
            VariantType::Undefined => self.undefined,
        }
        .unwrap_or(0)
    }
    pub fn increment(&mut self, t: VariantType) {
        let slot = match t {
            VariantType::Snp => &mut self.snp,
            VariantType::Indel => &mut self.indel,
            VariantType::Mnp => &mut self.mnp,
            VariantType::Undefined => &mut self.undefined,
        };
        *slot = Some(slot.unwrap_or(0) + 1);
    }
}

/// SomaticVariantDrivers.groupByImpact
fn group_by_impact(variants: &[&Variant]) -> HashMap<DriverImpact, i32> {
    let mut m = HashMap::new();
    for v in variants {
        *m.entry(driver_impact(v.vtype, v.impact.canonical_coding_effect)).or_insert(0) += 1;
    }
    m
}

/// TsgImpactComparator
fn tsg_compare(a: &Variant, b: &Variant) -> std::cmp::Ordering {
    let score = |v: &Variant| match driver_impact(v.vtype, v.impact.canonical_coding_effect) {
        DriverImpact::Frameshift => 0,
        DriverImpact::Nonsense => 1,
        DriverImpact::Splice => 2,
        DriverImpact::Missense => 3,
        DriverImpact::Inframe => 4,
        DriverImpact::Unknown => 5,
    };
    score(a).cmp(&score(b)).then(a.position().cmp(&b.position()))
}

impl<'a> DriverFinder<'a> {
    pub fn new(resource: &'a DriverGeneResource, category: DriverCategory) -> DriverFinder<'a> {
        DriverFinder {
            category,
            predicate: ReportablePredicate::new(category, &resource.driver_gene_list),
            likelihoods: if category == DriverCategory::Onco { &resource.onco_likelihood } else { &resource.tsg_likelihood },
            driver_variants: Vec::new(),
        }
    }

    pub fn check_variant(&mut self, index: usize, v: &Variant) -> bool {
        if self.predicate.is_candidate_reportable(&v.impact, v.is_hotspot()) {
            self.driver_variants.push(index);
            return self.predicate.is_reportable(&v.impact, v.is_hotspot());
        }
        false
    }

    pub fn find_drivers(
        &self,
        variants: &[Variant],
        gene_copy_numbers: &HashMap<String, Vec<GeneCopyNumber>>,
        counts: &VariantTypeCounts,
        biallelic_counts: &VariantTypeCounts,
        sources: &mut Vec<DriverSource>,
    ) -> Vec<DriverCatalog> {
        // Collectors.groupingBy(SomaticVariant::gene): a HashMap filled by computeIfAbsent, iterated in Java's order
        let mut groups: JavaHashMap<String, Vec<usize>> = JavaHashMap::new();
        for &i in &self.driver_variants {
            groups.compute_if_absent(variants[i].gene().to_string(), Vec::new).push(i);
        }
        let mut out = Vec::new();
        let no_gene = NO_GENE_DNDS;
        for (gene, indices) in groups.iter() {
            let dnds = self.likelihoods.get(gene).unwrap_or(&no_gene);
            let Some(gcns) = gene_copy_numbers.get(gene) else { continue };
            if !self.predicate.is_driver_gene(&gcns[0].gene_name) {
                continue;
            }
            let mut gene_variants: Vec<&Variant> = indices.iter().map(|&i| &variants[i]).collect();
            if self.category == DriverCategory::Tsg {
                // geneVariants.sort(TsgImpactComparator) happens inside each driver build; the first build sorts the list
                gene_variants.sort_by(|a, b| tsg_compare(a, b));
            }
            for gcn in gcns {
                if gcns.len() == 1 || gene_variants.iter().any(|x| has_transcript_coding_effect(&x.impact, x.vtype, &gcn.trans_name)) {
                    let coding = gene_variants.iter().any(|x| x.is_hotspot() || has_transcript_coding_effect(&x.impact, x.vtype, &gcn.trans_name));
                    let reportable = coding && gene_variants.iter().any(|x| self.predicate.is_reportable(&x.impact, x.is_hotspot()));
                    let status = if reportable { ReportedStatus::Reported } else { ReportedStatus::NotReported };
                    let method = if coding { LikelihoodMethod::Dnds } else { LikelihoodMethod::SpliceRegion };
                    let record = match self.category {
                        DriverCategory::Onco => onco_driver(&gene_variants, counts, gcn, dnds, method, status),
                        DriverCategory::Tsg => tsg_driver(&gene_variants, counts, biallelic_counts, gcn, dnds, method, status),
                    };
                    sources.push(DriverSource::from_variant(record.clone(), gene_variants[0]));
                    out.push(record);
                }
            }
        }
        out
    }
}

fn base_driver(variants: &[&Variant], gcn: &GeneCopyNumber, category: DriverCategory, method: LikelihoodMethod, status: ReportedStatus) -> DriverCatalog {
    let counts = group_by_impact(variants);
    let c = |i: DriverImpact| *counts.get(&i).unwrap_or(&0);
    DriverCatalog {
        chromosome: variants[0].chromosome().to_string(),
        chromosome_band: gcn.chromosome_band.clone(),
        gene: gcn.gene_name.clone(),
        transcript: gcn.trans_name.clone(),
        is_canonical: gcn.is_canonical,
        category,
        driver: DriverType::Mutation,
        likelihood_method: method,
        reported_status: status,
        driver_likelihood: if method != LikelihoodMethod::SpliceRegion { 1.0 } else { 0.0 },
        missense: c(DriverImpact::Missense),
        nonsense: c(DriverImpact::Nonsense),
        splice: c(DriverImpact::Splice),
        inframe: c(DriverImpact::Inframe),
        frameshift: c(DriverImpact::Frameshift),
        biallelic: variants.iter().any(|v| v.biallelic()),
        min_copy_number: gcn.min_copy_number,
        max_copy_number: gcn.max_copy_number,
    }
}

/// OncoDrivers.buildDriverCatalog
fn onco_driver(variants: &[&Variant], counts: &VariantTypeCounts, gcn: &GeneCopyNumber, dnds: &DndsGeneLikelihood, method: LikelihoodMethod, status: ReportedStatus) -> DriverCatalog {
    let mut d = base_driver(variants, gcn, DriverCategory::Onco, method, status);
    if method == LikelihoodMethod::SpliceRegion {
        return d;
    }
    if variants.iter().any(|v| v.is_hotspot()) {
        d.likelihood_method = LikelihoodMethod::Hotspot;
        return d;
    }
    if variants.iter().any(|v| v.vtype == VariantType::Indel && v.impact.canonical_coding_effect == CodingEffect::Missense && v.repeat_count() <= crate::constants::MAX_INDEL_DRIVER_REPEAT_COUNT) {
        d.likelihood_method = LikelihoodMethod::Inframe;
        return d;
    }
    let (snv, indel) = (counts.get(VariantType::Snp) as i64, counts.get(VariantType::Indel) as i64);
    let mut likelihood = 0.0;
    for v in variants {
        let effect = worst_reportable_coding_effect(&v.impact);
        if effect.has_protein_impact() {
            let impact = driver_impact(v.vtype, effect);
            let l = dnds.select(impact);
            let n = if matches!(impact, DriverImpact::Frameshift | DriverImpact::Inframe) { indel } else { snv };
            likelihood = jmax(likelihood, probability_driver_variant(n, &l));
        }
    }
    d.driver_likelihood = likelihood;
    d
}

/// TsgDrivers.buildDriverCatalog (the variants already sorted by TsgImpactComparator)
fn tsg_driver(
    variants: &[&Variant],
    counts: &VariantTypeCounts,
    biallelic: &VariantTypeCounts,
    gcn: &GeneCopyNumber,
    dnds: &DndsGeneLikelihood,
    method: LikelihoodMethod,
    status: ReportedStatus,
) -> DriverCatalog {
    let mut d = base_driver(variants, gcn, DriverCategory::Tsg, method, status);
    if method == LikelihoodMethod::SpliceRegion {
        return d;
    }
    if variants.iter().any(|v| v.is_hotspot()) {
        d.likelihood_method = LikelihoodMethod::Hotspot;
        return d;
    }
    let coding: Vec<&&Variant> = variants.iter().filter(|x| worst_reportable_coding_effect(&x.impact).has_protein_impact()).collect();
    if coding.iter().any(|x| x.biallelic() && !(x.vtype != VariantType::Indel && x.impact.canonical_coding_effect == CodingEffect::Missense)) {
        d.likelihood_method = LikelihoodMethod::Biallelic;
        return d;
    }
    let variant_count = |v: &Variant| {
        let map = if v.biallelic() { biallelic } else { counts };
        if v.vtype == VariantType::Indel { map.get(VariantType::Indel) } else { map.get(VariantType::Snp) }
    };
    let top = coding[0];
    let first_impact = driver_impact(top.vtype, worst_reportable_coding_effect(&top.impact));
    let first_l = dnds.select(first_impact);
    let first_count = variant_count(top) as i64;
    let missense_count = counts.get(VariantType::Snp) as i64;
    if coding.len() == 1 {
        let single = probability_driver_variant(first_count, &first_l);
        let substitute = if first_impact != DriverImpact::Missense { probability_driver_variant(missense_count, &dnds.missense) } else { single };
        d.driver_likelihood = jmax(single, substitute);
        return d;
    }
    let second = coding[1];
    let second_impact = driver_impact(second.vtype, worst_reportable_coding_effect(&second.impact));
    let second_l = dnds.select(second_impact);
    let second_count = variant_count(second) as i64;
    let multi = probability_driver_variant2(first_count, second_count, &first_l, &second_l);
    let sub_first = if first_impact == DriverImpact::Missense { multi } else { probability_driver_variant2(missense_count, second_count, &dnds.missense, &second_l) };
    let sub_second = if second_impact == DriverImpact::Missense { multi } else { probability_driver_variant2(first_count, missense_count, &first_l, &dnds.missense) };
    let sub_both = if first_impact == DriverImpact::Missense || second_impact == DriverImpact::Missense {
        multi
    } else {
        probability_driver_variant2(missense_count, missense_count, &dnds.missense, &dnds.missense)
    };
    d.driver_likelihood = jmax(jmax(sub_first, sub_second), sub_both);
    d
}

// ---------------------------------------------------------------------------------------------
// copy number drivers

pub fn copy_number_driver(category: DriverCategory, driver: DriverType, method: LikelihoodMethod, biallelic: bool, likelihood: f64, g: &GeneCopyNumber, status: ReportedStatus) -> DriverCatalog {
    DriverCatalog {
        chromosome: g.chromosome.clone(),
        chromosome_band: g.chromosome_band.clone(),
        gene: g.gene_name.clone(),
        transcript: g.trans_name.clone(),
        is_canonical: g.is_canonical,
        category,
        driver,
        likelihood_method: method,
        reported_status: status,
        driver_likelihood: likelihood,
        missense: 0,
        nonsense: 0,
        splice: 0,
        inframe: 0,
        frameshift: 0,
        biallelic,
        min_copy_number: g.min_copy_number,
        max_copy_number: g.max_copy_number,
    }
}

/// AmplificationDrivers.findAmplifications
pub fn find_amplifications(qc: &[QcStatus], gender: Gender, panel: &DriverGeneResource, ploidy: f64, genes: &mut [GeneCopyNumber], targeted: bool) -> Vec<DriverCatalog> {
    const PARTIAL_GENES: [&str; 7] = ["BRAF", "EGFR", "CTNNB1", "CBL", "MET", "ALK", "PDGFRA"];
    let high_noise = qc.contains(&QcStatus::WarnHighCopyNumberNoise);
    let mut out = Vec::new();
    for g in genes.iter_mut() {
        let Some(dg) = panel.gene(&g.gene_name) else { continue };
        if high_noise && !(g.min_region_start_support.is_sv() || g.min_region_end_support.is_sv()) {
            continue;
        }
        let amp_threshold = ploidy * dg.amplification_ratio;
        let is_x_male = gender == Gender::Male && Chromosome::from_name(&g.chromosome).is_some_and(is_x);
        let threshold = if is_x_male { amp_threshold * 0.5 } else { amp_threshold };
        let check_partials = !targeted || PARTIAL_GENES.contains(&g.gene_name.as_str());
        let driver_type = |t: f64| {
            if g.min_copy_number > t {
                DriverType::Amp
            } else if check_partials && g.max_copy_number > t {
                DriverType::PartialAmp
            } else {
                DriverType::Unknown
            }
        };
        let mut candidate = false;
        let mut dt = driver_type(threshold);
        if dt == DriverType::Unknown {
            dt = driver_type(threshold * 0.85);
            candidate = dt != DriverType::Unknown;
        }
        if dt == DriverType::Unknown {
            continue;
        }
        let status = if dg.report_amplification { ReportedStatus::Reported } else { ReportedStatus::NotReported };
        g.driver_type = dt;
        g.reported_status = status;
        out.push(copy_number_driver(DriverCategory::Onco, dt, LikelihoodMethod::Amp, false, if candidate { 0.0 } else { 1.0 }, g, status));
    }
    out
}

/// DeletionDrivers.deletedGenes
pub fn deleted_genes(genes: &[GeneCopyNumber]) -> i32 {
    genes.iter().filter(|g| !is_y(chr(&g.chromosome)) && less_than(g.min_copy_number, 0.5)).count() as i32
}

/// DeletionDrivers.findDeletions
pub fn find_deletions(qc: &[QcStatus], ploidy: f64, panel: &DriverGeneResource, genes: &mut [GeneCopyNumber], targeted: bool) -> Vec<DriverCatalog> {
    let check_qc = qc.contains(&QcStatus::WarnDeletedGenes) || qc.contains(&QcStatus::WarnHighCopyNumberNoise);
    let mere = |s: SegmentSupport| matches!(s, SegmentSupport::Centromere | SegmentSupport::Telomere);
    let mut out = Vec::new();
    for g in genes.iter_mut() {
        let Some(dg) = panel.gene(&g.gene_name) else { continue };
        let full_sv = g.min_region_start_support.is_sv() && g.min_region_end_support.is_sv();
        let is_short = g.min_region_bases() < 10_000_000;
        if targeted {
            if !full_sv {
                if check_qc {
                    continue;
                }
                if g.depth_window_count == 1 && (g.gc_content < 0.35 || g.gc_content > 0.6) {
                    continue;
                }
            }
        } else if check_qc {
            let one_sv_and_mere = (mere(g.min_region_start_support) && g.min_region_end_support.is_sv()) || (mere(g.min_region_end_support) && g.min_region_start_support.is_sv());
            if !full_sv && !(is_short && one_sv_and_mere) {
                continue;
            }
        }
        let mut dt = DriverType::Unknown;
        let mut status = ReportedStatus::None;
        if g.min_copy_number < 0.5 {
            dt = DriverType::Del;
            status = if dg.report_deletion { ReportedStatus::Reported } else { ReportedStatus::NotReported };
        } else if ploidy > 0.0 && is_autosome(chr(&g.chromosome)) {
            if g.min_copy_number / ploidy < dg.het_deletion_threshold {
                dt = DriverType::HetDel;
                status = if dg.report_het_deletion { ReportedStatus::Reported } else { ReportedStatus::NotReported };
            }
        }
        if dt == DriverType::Unknown && g.min_minor_allele_copy_number < 0.5 {
            dt = DriverType::Loh;
            status = if dg.report_loh { ReportedStatus::Reported } else { ReportedStatus::NotReported };
        }
        if dt == DriverType::Unknown {
            continue;
        }
        g.driver_type = dt;
        g.reported_status = status;
        let is_del = dt == DriverType::Del;
        out.push(copy_number_driver(category_of(dg.likelihood_type), dt, LikelihoodMethod::Del, is_del, if is_del { 1.0 } else { 0.0 }, g, status));
    }
    out
}

java_enum!(QcStatus {
    Pass = "PASS",
    WarnDeletedGenes = "WARN_DELETED_GENES",
    WarnHighCopyNumberNoise = "WARN_HIGH_COPY_NUMBER_NOISE",
    WarnGenderMismatch = "WARN_GENDER_MISMATCH",
    WarnLowPurity = "WARN_LOW_PURITY",
    WarnTinc = "WARN_TINC",
    WarnContamination = "WARN_CONTAMINATION",
    FailContamination = "FAIL_CONTAMINATION",
    FailNoTumor = "FAIL_NO_TUMOR",
    FailTinc = "FAIL_TINC",
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::DndsImpactLikelihood;
    use crate::test_support::{gene_copy_number, somatic, Hotspot};

    fn counts(snp: i32, indel: i32) -> VariantTypeCounts {
        VariantTypeCounts { snp: Some(snp), indel: Some(indel), mnp: None, undefined: None }
    }

    fn lh(drivers: f64, passengers: f64) -> DndsImpactLikelihood {
        DndsImpactLikelihood { drivers_per_sample: drivers, passengers_per_mutation: passengers }
    }

    // TsgDriversTest
    fn tsg(variants: &[&Variant], c: &VariantTypeCounts) -> DriverCatalog {
        let dnds = DndsGeneLikelihood {
            gene: "TP53".into(),
            missense: lh(0.0166790589295988, 9.89844130535209e-08),
            nonsense: lh(0.00730132326080717, 7.40370091523547e-09),
            splice: lh(0.00540540540540541, 1e-9),
            indel: lh(0.0137214137214137, 1e-9),
        };
        // TsgDrivers.buildDriverCatalog sorts the variants by TsgImpactComparator first
        let mut v = variants.to_vec();
        v.sort_by(|a, b| tsg_compare(a, b));
        tsg_driver(&v, c, c, &gene_copy_number("GENE_01"), &dnds, LikelihoodMethod::Dnds, ReportedStatus::Reported)
    }

    fn tsg_variants() -> [Variant; 5] {
        [
            somatic(VariantType::Snp, "MISSENSE", 0, Hotspot::None, 0.5),
            somatic(VariantType::Snp, "NONSENSE_OR_FRAMESHIFT", 0, Hotspot::None, 0.5),
            somatic(VariantType::Indel, "NONSENSE_OR_FRAMESHIFT", 0, Hotspot::None, 0.5),
            somatic(VariantType::Indel, "MISSENSE", 0, Hotspot::Hotspot, 0.5),
            somatic(VariantType::Mnp, "NONSENSE_OR_FRAMESHIFT", 0, Hotspot::None, 0.9),
        ]
    }

    fn approx(a: f64, b: f64, e: f64) {
        assert!((a - b).abs() <= e, "{a} != {b}");
    }

    #[test]
    fn tsg_drivers() {
        let [missense, nonsense, indel, hotspot, biallelic] = tsg_variants();
        let c = counts(5161, 10000);
        let d = tsg(&[&hotspot, &biallelic, &missense], &c);
        assert_eq!(d.likelihood_method, LikelihoodMethod::Hotspot);
        approx(d.driver_likelihood, 1.0, 0.01);
        let d = tsg(&[&biallelic, &missense], &c);
        assert_eq!(d.likelihood_method, LikelihoodMethod::Biallelic);
        approx(d.driver_likelihood, 1.0, 0.01);
        let c = counts(351610, 10000);
        for (vs, expected) in [(vec![&missense], 0.33), (vec![&missense, &missense], 0.97), (vec![&nonsense], 0.74), (vec![&missense, &nonsense], 1.0), (vec![&indel], 1.0)] {
            let d = tsg(&vs, &c);
            assert_eq!(d.likelihood_method, LikelihoodMethod::Dnds);
            approx(d.driver_likelihood, expected, 0.01);
        }
    }

    // OncoDriversTest
    const PASSENGERS: f64 = 9.26e-08;

    fn onco(variants: &[&Variant], c: &VariantTypeCounts) -> DriverCatalog {
        let l = |d| lh(d, PASSENGERS);
        let dnds = DndsGeneLikelihood { gene: "Gene".into(), missense: l(0.002), nonsense: l(0.003), splice: l(0.002), indel: l(0.001) };
        onco_driver(variants, c, &gene_copy_number("GENE_01"), &dnds, LikelihoodMethod::Dnds, ReportedStatus::Reported)
    }

    #[test]
    fn onco_drivers() {
        let fs_hotspot = somatic(VariantType::Indel, "NONSENSE_OR_FRAMESHIFT", 3, Hotspot::Hotspot, 0.5);
        let fs_near = somatic(VariantType::Indel, "NONSENSE_OR_FRAMESHIFT", 3, Hotspot::Near, 0.5);
        let inframe = somatic(VariantType::Indel, "MISSENSE", 3, Hotspot::None, 0.5);
        let unknown_inframe = somatic(VariantType::Indel, "MISSENSE", crate::constants::MAX_INDEL_DRIVER_REPEAT_COUNT + 1, Hotspot::None, 0.5);
        let frameshift = somatic(VariantType::Indel, "NONSENSE_OR_FRAMESHIFT", 3, Hotspot::None, 0.5);
        let missense = somatic(VariantType::Snp, "MISSENSE", 0, Hotspot::None, 0.5);
        let c = counts(10000, 0);
        let d = onco(&[&fs_hotspot, &fs_near, &inframe, &missense], &c);
        assert_eq!(d.likelihood_method, LikelihoodMethod::Hotspot);
        approx(d.driver_likelihood, 1.0, 0.01);
        for vs in [vec![&fs_near, &inframe, &missense], vec![&inframe, &missense]] {
            let d = onco(&vs, &c);
            assert_eq!(d.likelihood_method, LikelihoodMethod::Inframe);
            approx(d.driver_likelihood, 1.0, 0.01);
        }
        let c = counts(10000, 1000);
        let expected = |p: f64, n: i64| probability_driver_variant(n, &lh(p, PASSENGERS));
        for (vs, p, n) in [(vec![&frameshift, &missense], 0.001, 1000), (vec![&unknown_inframe], 0.001, 1000), (vec![&missense], 0.002, 10000), (vec![&missense, &missense], 0.002, 10000)] {
            let d = onco(&vs, &c);
            assert_eq!(d.likelihood_method, LikelihoodMethod::Dnds);
            approx(d.driver_likelihood, expected(p, n), 1e-10);
        }
    }

    // CopyNumberDriversTest
    fn panel() -> crate::reference::DriverGeneResource {
        let g = |name: &str| {
            let mut d = crate::test_support::driver_gene(name, hmf_common::driver_panel::DriverCategory::Onco);
            d.report_amplification = true;
            d.report_deletion = true;
            d.report_het_deletion = true;
            d
        };
        crate::reference::DriverGeneResource::new(vec![g("GENE1"), g("GENE2")])
    }

    fn gcn(chromosome: &str, gene: &str, min: f64, max: f64) -> GeneCopyNumber {
        let mut g = gene_copy_number(gene);
        g.chromosome = chromosome.into();
        g.chromosome_band = "p36.13".into();
        g.min_copy_number = min;
        g.max_copy_number = max;
        g.relative_min_copy_number = min / 2.0;
        g
    }

    fn gcn_support(gene: &str, min: f64, start: i32, end: i32, s: SegmentSupport, e: SegmentSupport) -> GeneCopyNumber {
        let mut g = gene_copy_number(gene);
        g.min_copy_number = min;
        g.min_region_start = start;
        g.min_region_end = end;
        g.min_region_start_support = s;
        g.min_region_end_support = e;
        g.relative_min_copy_number = min / 2.0;
        g
    }

    #[test]
    fn onco_amplification_without_dnds() {
        let p = panel();
        let pass = [QcStatus::Pass];
        let d = find_amplifications(&pass, Gender::Female, &p, 2.0, &mut [gcn("1", "GENE1", 7.0, 7.0)], false);
        assert_eq!((d[0].gene.as_str(), d[0].category), ("GENE1", DriverCategory::Onco));
        let d = find_amplifications(&pass, Gender::Male, &p, 2.0, &mut [gcn("X", "GENE1", 7.0, 7.0)], false);
        assert_eq!((d[0].gene.as_str(), d[0].category), ("GENE1", DriverCategory::Onco));
    }

    #[test]
    fn partial_amp() {
        let p = panel();
        let d = find_amplifications(&[QcStatus::Pass], Gender::Female, &p, 2.0, &mut [gcn("1", "GENE1", 0.1, 7.0), gcn("1", "GENE2", 7.0, 7.0)], false);
        assert_eq!(d.len(), 2);
        assert_eq!((d[0].gene.as_str(), d[0].driver), ("GENE1", DriverType::PartialAmp));
        assert_eq!((d[1].gene.as_str(), d[1].driver), ("GENE2", DriverType::Amp));
    }

    #[test]
    fn amp_with_warn() {
        let p = panel();
        let (n, b) = (SegmentSupport::None, SegmentSupport::Bnd);
        let amp = |s, e| gcn_support("GENE1", 100.0, 1, 10000, s, e);
        assert_eq!(find_amplifications(&[QcStatus::WarnDeletedGenes], Gender::Female, &p, 1.0, &mut [amp(n, n)], false).len(), 1);
        let noise = [QcStatus::WarnHighCopyNumberNoise];
        assert_eq!(find_amplifications(&noise, Gender::Female, &p, 1.0, &mut [amp(n, n)], false).len(), 0);
        for (s, e) in [(b, n), (n, b), (b, b)] {
            assert_eq!(find_amplifications(&noise, Gender::Female, &p, 1.0, &mut [amp(s, e)], false).len(), 1);
        }
    }

    #[test]
    fn del_with_warn() {
        let p = panel();
        let (n, b, t, c) = (SegmentSupport::None, SegmentSupport::Bnd, SegmentSupport::Telomere, SegmentSupport::Centromere);
        let del = |end: i32, s, e| gcn_support("GENE1", 0.01, 1, end, s, e);
        let long = [del(100_000_000, n, n), del(100_000_000, b, t), del(100_000_000, c, b), del(100_000_000, b, b)];
        let short = [del(100, n, n), del(100, b, t), del(100, c, b), del(100, b, b)];
        let count = |qc: &[QcStatus], g: &GeneCopyNumber, targeted: bool| find_deletions(qc, 2.0, &p, &mut [g.clone()], targeted).len();
        assert_eq!(count(&[], &long[0], false), 1);
        let warn = [QcStatus::WarnDeletedGenes];
        assert_eq!(long.iter().map(|g| count(&warn, g, false)).collect::<Vec<_>>(), vec![0, 0, 0, 1]);
        // (Java's test does not check the short deletion with start support here)
        assert_eq!([&short[0], &short[2], &short[3]].iter().map(|g| count(&warn, g, false)).collect::<Vec<_>>(), vec![0, 1, 1]);
        assert_eq!(long.iter().map(|g| count(&warn, g, true)).collect::<Vec<_>>(), vec![0, 0, 0, 1]);
        assert_eq!(short.iter().map(|g| count(&warn, g, true)).collect::<Vec<_>>(), vec![0, 0, 0, 1]);
        let noise = [QcStatus::WarnHighCopyNumberNoise];
        assert_eq!(long.iter().map(|g| count(&noise, g, false)).collect::<Vec<_>>(), vec![0, 0, 0, 1]);
        assert_eq!(short.iter().map(|g| count(&noise, g, false)).collect::<Vec<_>>(), vec![0, 1, 1, 1]);
    }

    #[test]
    fn heterozygous_deletion() {
        let p = panel();
        let mut g = [gcn_support("GENE1", 1.1, 1, 100_000_000, SegmentSupport::Bnd, SegmentSupport::Bnd)];
        let d = find_deletions(&[], 2.0, &p, &mut g, false);
        assert_eq!(d.len(), 1);
        assert_eq!(d[0].driver, DriverType::HetDel);
    }

    // TsgMultiHitTest
    fn jak1() -> DndsGeneLikelihood {
        DndsGeneLikelihood { gene: "JAK1".into(), missense: lh(0.0045, 3.6e-7), nonsense: lh(2.5e-4, 2.5e-8), splice: lh(2.1e-4, 9.2e-9), indel: lh(0.002, 2.5e-7) }
    }

    fn tsg_with(variants: &[&Variant], dnds: &DndsGeneLikelihood) -> DriverCatalog {
        let c = counts(83135, 241917);
        let mut v = variants.to_vec();
        v.sort_by(|a, b| tsg_compare(a, b));
        tsg_driver(&v, &c, &c, &gene_copy_number("GENE_01"), dnds, LikelihoodMethod::Dnds, ReportedStatus::Reported)
    }

    #[test]
    fn multi_hit_never_less_than_equivalent_missense() {
        let dnds = jak1();
        let missense = somatic(VariantType::Snp, "MISSENSE", 0, Hotspot::None, 0.5);
        let frameshift = somatic(VariantType::Indel, "NONSENSE_OR_FRAMESHIFT", 0, Hotspot::None, 0.5);
        let multi_fs = tsg_with(&[&frameshift, &frameshift], &dnds).driver_likelihood;
        let mixed = tsg_with(&[&missense, &frameshift], &dnds).driver_likelihood;
        let single_fs = tsg_with(&[&frameshift], &dnds).driver_likelihood;
        let multi_missense = tsg_with(&[&missense, &missense], &dnds).driver_likelihood;
        assert!(less_than(single_fs, mixed));
        assert!(equal(multi_missense, mixed));
        assert!(equal(multi_missense, multi_fs));
    }

    #[test]
    fn multi_hit_with_non_impacting_variant() {
        let dnds = jak1();
        let missense = somatic(VariantType::Snp, "MISSENSE", 0, Hotspot::None, 0.5);
        let splice_intron = somatic(VariantType::Snp, "UNDEFINED", 0, Hotspot::None, 0.5);
        approx(tsg_with(&[&missense, &splice_intron], &dnds).driver_likelihood, 0.13, 0.1);
    }
}
