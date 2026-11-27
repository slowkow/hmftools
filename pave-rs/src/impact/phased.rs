//! PhasedVariantClassifier: variants of one local phase set are held until the set is complete, then
//! INDELs (and SNVs/MNVs between or overlapping them) that together keep the reading frame are given
//! their combined, phased impact.

use anyhow::Result;

use hmf_common::codon::{amino_acid_from_bases, reverse_complement_bases};

use super::classifier::check_stop_start_codons;
use super::effects::VariantEffect;
use super::protein::trim_amino_acids;
use super::{hgvs_protein, jsub, jsub_from, or_null, ProteinContext, RefGenome, VariantTransImpact, SE_END, SE_START};
use crate::variant::VariantData;

const CODON_LENGTH: i32 = 3;

/// PhasedVariants: the variants of one local phase set, in the order they were added
pub struct PhasedVariants {
    pub local_phase_id: i32,
    pub variants: Vec<VariantData>,
}

#[derive(Default)]
pub struct PhasedVariantClassifier {
    groups: Vec<PhasedVariants>,
}

impl PhasedVariantClassifier {
    /// checkAddVariant: the variant joins its phase set's group
    pub fn add_variant(&mut self, variant: VariantData) {
        debug_assert!(variant.has_local_phase_set());
        let lps = variant.local_phase_set();
        match self.groups.iter_mut().find(|g| g.local_phase_id == lps) {
            Some(g) => g.variants.push(variant),
            None => self.groups.push(PhasedVariants { local_phase_id: lps, variants: vec![variant] }),
        }
    }

    pub fn has_complete_variants(&self, current_local_phase_set: i32) -> bool {
        !self.groups.is_empty() && self.groups[0].local_phase_id != current_local_phase_set
    }

    /// the groups at the front which are not the current phase set
    pub fn pop_complete_phased_variants(&mut self, current_local_phase_set: i32) -> Vec<PhasedVariants> {
        let mut complete = Vec::new();
        while !self.groups.is_empty() && self.groups[0].local_phase_id != current_local_phase_set {
            complete.push(self.groups.remove(0));
        }
        complete
    }

    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }
}

/// ImpactClassifier.processPhasedVariants: the variants of the completed groups, reclassified
pub fn process_phased_variants(classifier: &mut PhasedVariantClassifier, current_local_phase_set: i32, ref_genome: &dyn RefGenome) -> Vec<VariantData> {
    if !classifier.has_complete_variants(current_local_phase_set) {
        return Vec::new();
    }
    let mut out = Vec::new();
    for mut group in classifier.pop_complete_phased_variants(current_local_phase_set) {
        reclassify_phased_variants(&mut group, ref_genome);
        out.extend(group.variants);
    }
    out
}

/// reclassifyPhasedVariants
pub fn reclassify_phased_variants(group: &mut PhasedVariants, ref_genome: &dyn RefGenome) {
    if group.variants.len() < 2 {
        return;
    }
    // groups without INDELs are left as they are
    if !group.variants.iter().any(|x| x.is_indel()) {
        return;
    }
    // variants are taken to be in ascending position order
    let v0 = &group.variants[0];
    if v0.impacts().is_empty() {
        return;
    }

    // the impacts of the first variant by gene and transcript, realigned ones ignored
    let mut candidates: Vec<(String, usize)> = Vec::new();
    for (gene, impacts) in v0.impacts().iter() {
        for (k, ti) in impacts.iter().enumerate() {
            if !ti.has_coding_bases() || !ti.protein.as_ref().unwrap().valid_ref_codon() || ti.realigned() {
                continue;
            }
            candidates.push((gene.clone(), k));
        }
    }

    for (gene, k) in candidates {
        let trans_name = group.variants[0].gene_impacts(&gene).unwrap()[k].trans.trans_name.clone();
        // (variant index, impact index in the gene's list)
        let mut members: Vec<(usize, usize)> = vec![(0, k)];
        for i in 1..group.variants.len() {
            if let Some(list) = group.variants[i].gene_impacts(&gene) {
                if let Some(j) = list.iter().position(|x| x.trans.trans_name == trans_name) {
                    let ni = &list[j];
                    if ni.has_coding_bases() && ni.protein.as_ref().unwrap().valid_ref_codon() {
                        members.push((i, j));
                    }
                }
            }
        }

        let result = {
            let variants: Vec<&VariantData> = members.iter().map(|(i, _)| &group.variants[*i]).collect();
            let impacts: Vec<&VariantTransImpact> = members.iter().map(|(i, j)| &group.variants[*i].gene_impacts(&gene).unwrap()[*j]).collect();
            reclassify_impacts(group.local_phase_id, &variants, &impacts, ref_genome).map_err(|e| (e, variants, impacts))
        };

        match result {
            Ok(Some((combined_pc, combined_effects, ignored))) => {
                for (m, (i, j)) in members.iter().enumerate() {
                    // SNVs/MNVs not overlapping the codons of an INDEL, unless between INDELs, are left alone
                    if ignored[m] {
                        continue;
                    }
                    let ti = &mut group.variants[*i].impacts_mut().get_mut(&gene).unwrap()[*j];
                    ti.mark_phased_frameshift();
                    ti.protein = Some(combined_pc.clone());
                    ti.coding.is_frameshift = false;
                    // splice effects are kept
                    let splice_effects: Vec<VariantEffect> =
                        ti.effects().iter().copied().filter(|e| *e == VariantEffect::SpliceAcceptor || *e == VariantEffect::SpliceDonor).collect();
                    ti.clear_effects();
                    for e in &combined_effects {
                        ti.add_effect(*e);
                    }
                    for e in splice_effects {
                        ti.add_effect(e);
                    }
                }
            }
            Ok(None) => {}
            Err((e, variants, impacts)) => {
                let info: Vec<String> = variants.iter().map(|v| v.to_string()).collect();
                crate::warn!("failed to phase variants({}): {}", info.join(";"), e);
                for ti in impacts {
                    crate::warn!("transImpact({}) coding({:?}) protein({:?})", ti.trans.trans_name, ti.coding, ti.protein);
                }
            }
        }
    }
}

/// the impacted ref coding positions (ImpactedRefCodingData)
fn impacted_ref_coding_base_position(ti: &VariantTransImpact) -> (i32, i32) {
    let inframe_indel = ti.has_effect(VariantEffect::InframeDeletion) || ti.has_effect(VariantEffect::InframeInsertion);
    let pc = ti.protein.as_ref().unwrap();
    let mut start = pc.ref_coding_base_position(SE_START);
    let mut end = pc.ref_coding_base_position(SE_END);
    if !inframe_indel {
        return (start, end);
    }
    let ref_codon_length_diff = (pc.ref_codon_bases.len() as i32 - pc.alt_codon_bases.len() as i32).abs();
    let net_codon_impact = pc.net_codon_index_range[SE_END] - pc.net_codon_index_range[SE_START] + 1;
    if net_codon_impact > 0 && ref_codon_length_diff == net_codon_impact * CODON_LENGTH {
        // strip the unaffected ref codon from the start and/or end
        if net_codon_impact > 1 {
            start += CODON_LENGTH;
            end -= CODON_LENGTH;
        } else if start == ti.coding.coding_position_range[SE_START] {
            end -= CODON_LENGTH;
        } else {
            start += CODON_LENGTH;
        }
    }
    (start, end)
}

type Reclassified = (ProteinContext, Vec<VariantEffect>, Vec<bool>);

/// reclassifyImpacts: the combined protein context and effects, and which members are ignored; None
/// where Java returns without changes
fn reclassify_impacts(
    local_phase_set: i32, variants: &[&VariantData], impacts: &[&VariantTransImpact], ref_genome: &dyn RefGenome,
) -> Result<Option<Reclassified>> {
    // at least 2 with coding impacts
    if impacts.len() < 2 {
        return Ok(None);
    }
    if impacts.iter().any(|x| x.has_effect(VariantEffect::StopLost) || x.has_effect(VariantEffect::StartLost)) {
        return Ok(None);
    }

    let n = impacts.len();
    let mut indel_base_total = 0;
    let mut frameshift_count = 0;
    let mut min_indel_index = n;
    let mut max_indel_index = 0;

    for i in 0..n {
        let v = variants[i];
        if !v.is_indel() {
            continue;
        }
        indel_base_total += if v.is_insert() { v.base_diff() } else { -impacts[i].coding.deleted_coding_bases };
        if impacts[i].coding.is_frameshift {
            frameshift_count += 1;
        }
        min_indel_index = min_indel_index.min(i);
        max_indel_index = max_indel_index.max(i);
    }

    // SNVs and MNVs neither between INDELs nor overlapping another variant are ignored
    let mut ignored = vec![false; n];
    let refs: Vec<(i32, i32)> = impacts.iter().map(|ti| impacted_ref_coding_base_position(ti)).collect();

    // and an overlap between an inframe INDEL and an SNV/MNV
    let mut has_inframe_snv_mnv_overlap = false;

    for i in 0..n {
        let v = variants[i];
        if i > 0 && variants[i - 1].position >= v.position {
            return Ok(None); // overlapping or out-of-order variants are not handled
        }
        if v.is_indel() {
            continue;
        }
        let mut overlaps_on_start = false;
        if i > 0 {
            overlaps_on_start = refs[i - 1].1 >= refs[i].0;
        }
        let mut overlaps_on_end = false;
        if !overlaps_on_start && i < n - 1 {
            overlaps_on_end = refs[i].1 >= refs[i + 1].0;
        }
        if (i < min_indel_index || i > max_indel_index) && !overlaps_on_start && !overlaps_on_end {
            // an SNV/MNV outside the INDEL range and not overlapping
            ignored[i] = true;
            continue;
        }
        // whether the overlap is with an inframe INDEL
        if overlaps_on_start && variants[i - 1].is_indel() && !impacts[i - 1].coding.is_frameshift {
            has_inframe_snv_mnv_overlap = true;
        }
        if overlaps_on_end && variants[i + 1].is_indel() && !impacts[i + 1].coding.is_frameshift {
            has_inframe_snv_mnv_overlap = true;
        }
    }

    // 2+ out-of-frame INDELs, or an INDEL and an overlapping SNV/MNV
    if !has_inframe_snv_mnv_overlap && frameshift_count < 2 {
        return Ok(None);
    }
    if indel_base_total % 3 != 0 {
        crate::trace!("lps({}) varCount({}) out-of-frame with indel base count({})", local_phase_set, variants.len(), indel_base_total);
        return Ok(None);
    }

    let chromosome = &variants[0].chromosome;
    let pos_strand = impacts[0].trans.pos_strand();

    let mut combined_ref = String::new();
    let mut combined_alt = String::new();
    let mut min_codon_index = -1;
    let mut last_impacted_ref_codon_end = 0;
    let mut last_ref_codon_end = 0;

    for i in 0..n {
        if ignored[i] {
            continue;
        }
        let v = variants[i];
        let pc = impacts[i].protein.as_ref().unwrap();

        // the recorded ref codon bases differ from those actually involved in or impacted by the variant;
        // overlaps are judged on the impacted ones
        let (impacted_start, _impacted_end) = refs[i];
        let impacted_end = refs[i].1;
        let ref_codon_start = pc.ref_coding_base_start();
        let ref_codon_end = pc.ref_coding_base_end();

        let overlaps_on_start = last_impacted_ref_codon_end > 0 && impacted_start <= last_impacted_ref_codon_end;

        if min_codon_index < 0 {
            min_codon_index = pc.codon_index;
        } else {
            min_codon_index = min_codon_index.min(pc.codon_index);
        }

        if overlaps_on_start {
            // the ref codons first
            let ref_overlap = last_ref_codon_end - ref_codon_start + 1;
            if ref_overlap <= pc.ref_codon_bases.len() as i32 {
                combined_ref += jsub_from(&pc.ref_codon_bases, ref_overlap)?;
            }

            // the alt codons, extending no further than the ref codons and ending at a codon boundary; the
            // alt bases come from the current variant since they hold the change
            let lower_non_alt_bases = v.position - ref_codon_start;
            let (current_alt_bases_trimmed, prev_alt_bases_trimmed) =
                if ref_overlap >= lower_non_alt_bases { (lower_non_alt_bases, ref_overlap - lower_non_alt_bases) } else { (ref_overlap, 0) };

            if combined_alt.len() as i32 - prev_alt_bases_trimmed > 0 {
                let mut previous_extra_alt_bases = String::new();
                if last_ref_codon_end > ref_codon_end {
                    previous_extra_alt_bases = jsub_from(&combined_alt, last_ref_codon_end - ref_codon_end + 1)?.to_string();
                }
                combined_alt = jsub(&combined_alt, 0, combined_alt.len() as i32 - prev_alt_bases_trimmed)?.to_string();
                combined_alt += jsub_from(&pc.alt_codon_bases, current_alt_bases_trimmed)?;
                // restore the trimmed ref bases the previous variants had
                combined_alt += &previous_extra_alt_bases;
            } else {
                crate::trace!(
                    "phasing variants LPS({}) var({}) skip adjusting combinedAltCodons({}) prevAltBasesTrimmed({})",
                    local_phase_set,
                    v,
                    combined_alt,
                    prev_alt_bases_trimmed
                );
            }
        } else if last_ref_codon_end == 0 {
            combined_ref += &pc.ref_codon_bases;
            combined_alt += &pc.alt_codon_bases;
        } else if ref_codon_start > last_ref_codon_end + 1 {
            // fill in any gaps between the codons
            let gap = or_null(ref_genome.base_string(chromosome, last_ref_codon_end + 1, ref_codon_start - 1)?);
            combined_ref += &gap;
            combined_alt += &gap;
            combined_ref += &pc.ref_codon_bases;
            combined_alt += &pc.alt_codon_bases;
        } else {
            combined_ref += jsub_from(&pc.ref_codon_bases, last_ref_codon_end - ref_codon_start + 1)?;
            combined_alt += jsub_from(&pc.alt_codon_bases, last_ref_codon_end - ref_codon_start + 1)?;
        }

        last_ref_codon_end = ref_codon_end.max(last_ref_codon_end);
        last_impacted_ref_codon_end = impacted_end.max(last_impacted_ref_codon_end);
    }

    if combined_ref.len() % 3 != 0 || combined_alt.len() % 3 != 0 {
        return Ok(None);
    }

    let mut pc = ProteinContext { codon_index: min_codon_index, ..Default::default() };
    if pos_strand {
        pc.ref_amino_acids = amino_acid_from_bases(combined_ref.as_bytes());
        pc.alt_amino_acids = amino_acid_from_bases(combined_alt.as_bytes());
    } else {
        pc.ref_amino_acids = amino_acid_from_bases(&reverse_complement_bases(combined_ref.as_bytes()));
        pc.alt_amino_acids = amino_acid_from_bases(&reverse_complement_bases(combined_alt.as_bytes()));
    }
    pc.ref_codon_bases = combined_ref;
    pc.alt_codon_bases = combined_alt;

    trim_amino_acids(&mut pc, pos_strand, true, true, false);

    // synonymous or missense from the combined bases
    let combined_effect = if indel_base_total == 0 {
        let e = if pc.ref_amino_acids == pc.alt_amino_acids { VariantEffect::PhasedSynonymous } else { VariantEffect::PhasedMissense };
        crate::trace!(
            "lps({}) varCount({}) combinedEffect({:?}) from aminoAcids({} -> {})",
            local_phase_set,
            variants.len(),
            e,
            pc.ref_amino_acids,
            pc.alt_amino_acids
        );
        e
    } else {
        crate::trace!("lps({}) varCount({}) indelBaseTotal({})", local_phase_set, variants.len(), indel_base_total);
        if indel_base_total > 0 { VariantEffect::PhasedInframeInsertion } else { VariantEffect::PhasedInframeDeletion }
    };

    let mut combined_effects = vec![combined_effect];
    if let Some(ss) = check_stop_start_codons(min_codon_index, &pc.ref_amino_acids, &pc.alt_amino_acids) {
        if ss == VariantEffect::StopLost || ss == VariantEffect::StartLost {
            // superseded if present
            combined_effects.retain(|e| *e != VariantEffect::PhasedMissense);
        }
        combined_effects.push(ss);
    }

    pc.is_phased = true;
    pc.hgvs = hgvs_protein::generate(&pc, &combined_effects)?;

    Ok(Some((pc, combined_effects, ignored)))
}
