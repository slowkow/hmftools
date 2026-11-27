//! VariantImpactBuilder and hmf-common VariantImpact: the variant's summary impact, from its canonical
//! transcript impacts and the worst impact, driver genes taking precedence over other genes.

use super::effects::CodingEffect;
use super::VariantTransImpact;
use crate::gene_cache::GeneDataCache;
use crate::variant::VariantData;

/// AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM
pub const VAR_IMPACT_OTHER_REPORT_DELIM: &str = "--";

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
    /// VariantImpactSerialiser.toVcfData
    pub fn to_vcf_data(&self) -> Vec<String> {
        vec![
            self.gene_name.clone(),
            self.canonical_transcript.clone(),
            self.canonical_effect.clone(),
            self.canonical_coding_effect.name().to_string(),
            self.canonical_splice_region.to_string(),
            self.canonical_hgvs_coding.clone(),
            self.canonical_hgvs_protein.clone(),
            self.other_reportable_effects.clone(),
            self.worst_coding_effect.name().to_string(),
            self.genes_affected.to_string(),
        ]
    }
}

/// AltTranscriptReportableInfo.serialise, eg CDKN2A|ENST00000579755|c.209_210delCCinsTT|p.Pro70Leu|missense_variant|MISSENSE
pub fn serialise_alt_transcript(gene: &str, trans: &str, hgvs_coding: &str, hgvs_protein: &str, effects: &str, ce: CodingEffect) -> String {
    format!("{gene}|{trans}|{hgvs_coding}|{hgvs_protein}|{effects}|{}", ce.name())
}

/// The effects of the other reportable transcripts, as parseAltTranscriptInfo reads them
pub struct AltTranscriptInfo {
    pub effects: String,
}

/// AltTranscriptReportableInfo.parseAltTranscriptInfo (only the effects are needed)
pub fn parse_alt_transcript_info(other: &str) -> Vec<AltTranscriptInfo> {
    let mut out = Vec::new();
    // String.split(delim, -1): trailing empty items kept
    for item in other.split(VAR_IMPACT_OTHER_REPORT_DELIM) {
        let values: Vec<&str> = item.split('|').collect();
        if values.len() < 5 {
            continue;
        }
        let mut idx = 0;
        if values.len() == 6 {
            idx += 1;
        }
        // transcript, hgvs coding, hgvs protein, effects
        out.push(AltTranscriptInfo { effects: values[idx + 3].to_string() });
    }
    out
}

pub fn determine_coding_effect(ti: &VariantTransImpact) -> CodingEffect {
    let simplified: Vec<CodingEffect> = ti.effects().iter().map(|e| CodingEffect::effect(*e)).collect();
    for ce in [CodingEffect::NonsenseOrFrameshift, CodingEffect::Splice, CodingEffect::Missense, CodingEffect::Synonymous] {
        if simplified.contains(&ce) {
            return ce;
        }
    }
    CodingEffect::None
}

/// createVariantImpact: None without impacts on any transcript
pub fn create_variant_impact(variant: &VariantData, gene_cache: &GeneDataCache) -> Option<VariantImpact> {
    if variant.impacts().is_empty() {
        return None;
    }

    // the coding impacts of the canonical transcript and of any other transcripts with worse impacts;
    // driver genes are favoured over others, otherwise the gene with the worst impact is taken
    let mut worst: Option<&VariantTransImpact> = None;
    let mut worst_rank = -1;
    let mut worst_canonical: Option<&VariantTransImpact> = None;
    let mut worst_canonical_rank = -1;
    let mut canonical_gene_name = String::new(); // can be the gene of other reportable transcripts if the canonical has no impact
    let mut canonical_splice_region = false;
    let mut worst_gene_name = String::new();
    let mut other_reportable_trans_data = String::new();
    let mut has_driver_gene = false;

    for (gene_name, gene_impacts) in variant.impacts().iter() {
        let is_driver_gene = gene_cache.is_driver_panel_gene(gene_name);
        if !has_driver_gene && is_driver_gene {
            // ignore any prior non-driver gene
            has_driver_gene = true;
            worst = None;
            worst_rank = -1;
            worst_canonical = None;
            worst_canonical_rank = -1;
            canonical_gene_name.clear();
            worst_gene_name.clear();
        } else if has_driver_gene && !is_driver_gene {
            continue;
        }

        for ti in gene_impacts {
            let rank = ti.top_rank();
            if rank > worst_rank || (rank == worst_rank && ti.trans.is_canonical) {
                worst_rank = rank;
                worst = Some(ti);
                worst_gene_name = gene_name.clone();
            }
            if ti.trans.is_canonical && rank > worst_canonical_rank {
                worst_canonical_rank = rank;
                worst_canonical = Some(ti);
                canonical_gene_name = gene_name.clone();
            }
        }

        if is_driver_gene {
            if let Some(other_trans) = gene_cache.other_reportable_transcripts(gene_name) {
                let mut items: Vec<String> = Vec::new();
                for ti in gene_impacts {
                    if other_trans.contains(&ti.trans.trans_name) {
                        let ce = determine_coding_effect(ti);
                        items.push(serialise_alt_transcript(gene_name, &ti.trans.trans_name, ti.hgvs_coding(), ti.hgvs_protein(), &ti.effects_str(), ce));
                        if canonical_gene_name.is_empty() {
                            canonical_gene_name = gene_name.clone();
                        }
                        canonical_splice_region |= ti.in_splice_region();
                    }
                }
                other_reportable_trans_data = items.join(VAR_IMPACT_OTHER_REPORT_DELIM);
            }
        }
    }

    let worst_coding_effect = worst.map(determine_coding_effect).unwrap_or(CodingEffect::None);
    let gene_name = if !canonical_gene_name.is_empty() { canonical_gene_name } else { worst_gene_name };

    let mut vi = VariantImpact {
        gene_name,
        canonical_transcript: String::new(),
        canonical_effect: String::new(),
        canonical_coding_effect: CodingEffect::None,
        canonical_hgvs_coding: String::new(),
        canonical_hgvs_protein: String::new(),
        canonical_splice_region,
        other_reportable_effects: other_reportable_trans_data,
        worst_coding_effect,
        genes_affected: variant.impacts().len() as i32,
    };
    if let Some(c) = worst_canonical {
        vi.canonical_effect = c.effects_str();
        vi.canonical_coding_effect = determine_coding_effect(c);
        vi.canonical_hgvs_coding = c.hgvs_coding().to_string();
        vi.canonical_hgvs_protein = c.hgvs_protein().to_string();
        vi.canonical_transcript = c.trans.trans_name.clone();
        vi.canonical_splice_region |= c.in_splice_region();
    }
    Some(vi)
}
