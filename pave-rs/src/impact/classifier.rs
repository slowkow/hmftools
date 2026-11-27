//! ImpactClassifier: the impact of a variant on one transcript.

use anyhow::Result;

use hmf_common::ensembl::TranscriptData;

use super::effects::{is_splice, VariantEffect};
use super::hgvs_protein::{self, report_protein_impact, HGVS_SPLICE_UNKNOWN};
use super::splice::{self, check_straddles_splice_region, is_insert_into_exon_start, is_within_splice_region};
use super::utils::within_trans_range;
use super::{coding, coding_end, coding_start, hgvs_coding, protein, ProteinContext, RefGenome, SpliceImpactType, VariantTransImpact};
use crate::variant::VariantData;

/// classifyVariant: None when the variant is outside the transcript's range
pub fn classify_variant(variant: &VariantData, trans: &'static TranscriptData, ref_genome: &dyn RefGenome) -> Result<Option<VariantTransImpact>> {
    // a conservative check by 1 base for an INDEL at the boundary
    if !within_trans_range(trans, variant.position, variant.end_position) {
        return Ok(None);
    }

    let mut ti = VariantTransImpact::new(trans);

    coding::determine_context(variant, trans, &mut ti.coding)?;
    hgvs_coding::set(variant, &mut ti.coding);

    if ti.coding.is_coding() {
        let pc = protein::determine_context(variant, &ti.coding, trans, ref_genome)?;
        ti.protein = Some(pc);
    }

    check_pre_post_coding_impact(variant, &mut ti)?;

    let mut in_splice_region = false;

    // each exon: whether the variant falls in a splice region and is intronic or exonic
    let exon_count = trans.exons.len();
    for i in 0..exon_count {
        let exon = &trans.exons[i];
        let next_exon = if i < exon_count - 1 { Some(&trans.exons[i + 1]) } else { None };

        if !in_splice_region && is_within_splice_region(variant, trans, exon) {
            in_splice_region = true;
            splice::classify_variant(variant, &mut ti, exon, ref_genome)?;
        } else if !in_splice_region && next_exon.is_some_and(|ne| is_within_splice_region(variant, trans, ne)) {
            in_splice_region = true;
            splice::classify_variant(variant, &mut ti, next_exon.unwrap(), ref_genome)?;
        }

        if variant.alt_positions_overlap(exon.start, exon.end) || is_insert_into_exon_start(variant, exon, trans.pos_strand()) {
            classify_exonic_position(variant, &mut ti);
            break;
        }

        if let Some(ne) = next_exon {
            if variant.alt_positions_within(exon.end + 1, ne.start - 1) {
                ti.add_effect(VariantEffect::Intronic);
                break;
            }
        }

        // a variant which fully straddles the splice region
        if !in_splice_region {
            if let Some(effect) = check_straddles_splice_region(variant, trans, exon) {
                ti.add_effect(effect);
                break;
            }
        }
    }

    if in_splice_region {
        ti.mark_splice_region();
    }

    check_stop_start_codons_impact(&mut ti);

    if ti.has_protein_context() {
        let protein_effects: Vec<VariantEffect> = ti.effects().iter().copied().filter(|e| report_protein_impact(*e)).collect();
        if !protein_effects.is_empty() {
            let hgvs = hgvs_protein::generate(ti.protein.as_ref().unwrap(), &protein_effects)?;
            ti.protein.as_mut().unwrap().hgvs = hgvs;
        }
    }

    // override for splice-effect variants
    if ti.top_effect().is_some_and(is_splice) && !trans.non_coding() {
        if !ti.has_protein_context() {
            ti.protein = Some(ProteinContext::default());
        }
        ti.protein.as_mut().unwrap().hgvs = HGVS_SPLICE_UNKNOWN.to_string();
    }

    Ok(Some(ti))
}

fn check_pre_post_coding_impact(variant: &VariantData, ti: &mut VariantTransImpact) -> Result<bool> {
    let trans = ti.trans;

    // the pre-gene region
    if (trans.pos_strand() && variant.alt_bases_below(trans.trans_start)) || (!trans.pos_strand() && variant.alt_bases_above(trans.trans_end)) {
        ti.add_effect(VariantEffect::UpstreamGene);
        return Ok(true);
    }

    if !ti.is_exonic() || trans.non_coding() {
        return Ok(false);
    }

    let (cs, ce) = (coding_start(trans)?, coding_end(trans)?);
    if trans.pos_strand() {
        // 5' and 3' UTR
        if variant.alt_bases_below(cs) {
            ti.add_effect(VariantEffect::FivePrimeUtr);
            return Ok(true);
        }
        if variant.alt_bases_above(ce) {
            ti.add_effect(VariantEffect::ThreePrimeUtr);
            return Ok(true);
        }
    } else {
        if variant.alt_bases_below(cs) {
            ti.add_effect(VariantEffect::ThreePrimeUtr);
            return Ok(true);
        }
        if variant.alt_bases_above(ce) {
            ti.add_effect(VariantEffect::FivePrimeUtr);
            return Ok(true);
        }
    }
    Ok(false)
}

fn classify_exonic_position(variant: &VariantData, ti: &mut VariantTransImpact) {
    if ti.trans.non_coding() {
        ti.add_effect(VariantEffect::NonCodingTranscript);
        return;
    }
    if ti.has_effect(VariantEffect::FivePrimeUtr) || ti.has_effect(VariantEffect::ThreePrimeUtr) {
        return;
    }
    if variant.is_indel() {
        if ti.coding.is_frameshift {
            ti.add_effect(VariantEffect::Frameshift);
        } else if variant.is_deletion() {
            ti.add_effect(VariantEffect::InframeDeletion);
        } else {
            ti.add_effect(VariantEffect::InframeInsertion);
        }
    } else if ti.has_amino_acids() && ti.protein.as_ref().unwrap().has_protein_change() {
        // Java means to leave out an alt of just the stop codon here, but compares the String with the char
        // 'X', which never matches; the missense it adds is removed again once the stop is gained
        ti.add_effect(VariantEffect::Missense);
    } else {
        ti.add_effect(VariantEffect::Synonymous);
    }
}

fn check_stop_start_codons_impact(ti: &mut VariantTransImpact) {
    if !ti.has_amino_acids() {
        return;
    }
    let pc = ti.protein.as_ref().unwrap();
    let Some(ss_effect) = check_stop_start_codons(pc.codon_index, &pc.ref_amino_acids, &pc.alt_amino_acids) else { return };
    // remove any superseded effects
    if ss_effect == VariantEffect::StopLost || ss_effect == VariantEffect::StartLost {
        ti.remove_effect(VariantEffect::Missense);
        ti.remove_effect(VariantEffect::Frameshift);
    }
    if ss_effect == VariantEffect::StopGained {
        ti.remove_effect(VariantEffect::Missense);
    }
    ti.add_effect(ss_effect);
}

/// checkStopStartCodons
pub fn check_stop_start_codons(amino_acid_start_pos: i32, ref_aas: &str, alt_aas: &str) -> Option<VariantEffect> {
    if ref_aas.is_empty() || alt_aas.is_empty() {
        return None;
    }
    let r = ref_aas.as_bytes();
    let a = alt_aas.as_bytes();
    if amino_acid_start_pos == 1 && r[0] == b'M' && a[0] != b'M' {
        return Some(VariantEffect::StartLost);
    }
    if r[r.len() - 1] != b'X' {
        // has the alt gained a stop codon anywhere
        if a.contains(&b'X') {
            return Some(VariantEffect::StopGained);
        }
    } else if !a.contains(&b'X') {
        return Some(VariantEffect::StopLost);
    }
    None
}

/// selectAlignedImpacts: whether to take the right-aligned variant's impact (marked realigned, with the
/// splice impact type of the original if it has none)
pub fn select_aligned_impacts(ti: &VariantTransImpact, ra: &mut VariantTransImpact) -> bool {
    // no right-alignment if a stop codon was created
    if ti.has_effect(VariantEffect::StopGained) && ti.trans.pos_strand() {
        return false;
    }

    // the least impact, favouring non-splice over splice
    let non_ra_has_splice = ti.effects().iter().any(|e| is_splice(*e));
    let ra_has_splice = ra.effects().iter().any(|e| is_splice(*e));

    let mut use_realigned = false;
    if ra_has_splice || non_ra_has_splice {
        if non_ra_has_splice && !ra_has_splice {
            use_realigned = true;
        }
    } else if ra.top_rank() != ti.top_rank() {
        use_realigned = ra.top_rank() < ti.top_rank();
    } else if let (Some(rp), Some(tp)) = (&ra.protein, &ti.protein) {
        if rp.is_duplication != tp.is_duplication {
            use_realigned = rp.is_duplication;
        }
    }

    if use_realigned {
        ra.mark_realigned();
        if ti.splice_impact_type() != SpliceImpactType::Unknown && ra.splice_impact_type() == SpliceImpactType::Unknown {
            ra.set_splice_impact_type(ti.splice_impact_type());
        }
    }
    use_realigned
}
