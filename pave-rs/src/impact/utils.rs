//! PaveUtils: finding a variant's transcript impacts, its right-aligned form, and codon numbering.

use anyhow::Result;

use hmf_common::ensembl::TranscriptData;
use hmf_common::regions::positions_overlap;

use super::classifier::{classify_variant, select_aligned_impacts};
use super::RefGenome;
use crate::constants::GENE_UPSTREAM_DISTANCE;
use crate::gene_cache::{GeneCacheIndexing, GeneDataCache};
use crate::variant::{VariantData, NO_LOCAL_PHASE_SET};

/// findVariantImpacts: the impacts on every transcript of the genes near the variant, the right-aligned
/// variant's where it has one. Returns whether any transcript was classified.
pub fn find_variant_impacts(
    variant: &mut VariantData, gene_cache: &'static GeneDataCache, indexing: Option<&mut GeneCacheIndexing>, ref_genome: &dyn RefGenome,
) -> Result<bool> {
    let mut processed = false;
    let genes = gene_cache.find_genes(&variant.chromosome, variant.position, variant.end_position, indexing);

    // each of the genes and their transcripts
    for gene in genes {
        let trans_list = gene_cache.find_transcripts(&gene.gene_id, variant.position, variant.end_position);
        // non-coding transcripts are skipped for now
        if trans_list.is_empty() {
            continue;
        }
        for trans in trans_list {
            let ti = classify_variant(variant, trans, ref_genome)?;
            processed = true;

            let mut selected = ti;

            // the right-aligned variant if it has microhomology
            let ra = match variant.realigned_variant() {
                Some(rv) => classify_variant(rv, trans, ref_genome)?,
                None => None,
            };
            if let Some(mut ra_ti) = ra {
                let use_realigned = match &selected {
                    Some(t) => select_aligned_impacts(t, &mut ra_ti),
                    None => anyhow::bail!("NullPointerException: no impact for {variant} on {}", trans.trans_name),
                };
                variant.realigned_variant_mut().unwrap().add_impact(&gene.gene_name, ra_ti.clone());
                if use_realigned {
                    selected = Some(ra_ti);
                } else if trans.pos_strand() {
                    // the right-aligned coding string, whichever variant is prioritised, on the positive strand
                    selected.as_mut().unwrap().coding.hgvs = ra_ti.coding.hgvs.clone();
                }
            }

            if let Some(s) = selected {
                variant.add_impact(&gene.gene_name, s);
            }
        }
    }
    Ok(processed)
}

pub fn within_trans_range(trans: &TranscriptData, pos_start: i32, pos_end: i32) -> bool {
    if trans.pos_strand() {
        positions_overlap(pos_start, pos_end, trans.trans_start - GENE_UPSTREAM_DISTANCE, trans.trans_end)
    } else {
        positions_overlap(pos_start, pos_end, trans.trans_start, trans.trans_end + GENE_UPSTREAM_DISTANCE)
    }
}

/// createRightAlignedVariant: an INDEL with microhomology shifted to its most 3' position
pub fn create_right_aligned_variant(variant: &VariantData, ref_genome: &dyn RefGenome) -> Result<Option<VariantData>> {
    if variant.is_base_change() {
        return Ok(None);
    }
    if variant.microhomology().is_empty() || variant.microhomology() == "." {
        return Ok(None);
    }
    // no realigning if the ref and alt bases differ, eg T>AA
    if variant.ref_.as_bytes()[0] != variant.alt.as_bytes()[0] {
        return Ok(None);
    }

    // the repeat count is used only where the alt bases match the microhomology and the repeat sequence;
    // otherwise shift by the microhomology:
    // - MH = indel sequence = repeat sequence: to the end of the repeat (-1 repeat for a DEL)
    // - MH = indel sequence = N * repeat sequence: likewise, -N repeats for a DEL
    let alt_bases = if variant.is_deletion() { &variant.ref_[1..] } else { &variant.alt[1..] };
    let mh = variant.microhomology();
    let mc_length = mh.len() as i32;

    let shift_count = if mh == alt_bases && variant.repeat_count() > 0 {
        let mut repeat_count = 0;
        if variant.repeat_sequence() == alt_bases {
            repeat_count = variant.repeat_count();
        } else {
            // whether MH = N * repeat sequence
            let repeat_seq_len = variant.repeat_sequence().len() as i32;
            if repeat_seq_len == 0 {
                anyhow::bail!("ArithmeticException: / by zero (empty repeat sequence for {variant})");
            }
            let n_count = mc_length / repeat_seq_len;
            let test_mh = variant.repeat_sequence().repeat(n_count.max(0) as usize);
            if test_mh == mh {
                repeat_count = variant.repeat_count() / n_count;
            }
        }
        if repeat_count >= 1 {
            if variant.is_insert() { repeat_count } else { repeat_count - 1 }
        } else {
            1
        }
    } else {
        1
    };

    let shift_length = mc_length * shift_count;

    // eg pos 100 AG>A with G repeated 5 times
    let new_position = variant.position + shift_length;
    let Some(new_ref) = ref_genome.base_string(&variant.chromosome, new_position, new_position + variant.ref_.len() as i32 - 1)? else {
        anyhow::bail!("NullPointerException: no ref bases for the realigned {variant}");
    };
    let new_alt = if variant.is_insert() { format!("{}{}", &new_ref[..1], &variant.alt[1..]) } else { new_ref[..1].to_string() };

    let mut ra = VariantData::new(&variant.chromosome, new_position, &new_ref, &new_alt);
    // phasing realigned variants is not attempted
    ra.set_variant_details(NO_LOCAL_PHASE_SET, variant.microhomology(), variant.repeat_sequence(), variant.repeat_count());
    Ok(Some(ra))
}

/// the amino acid of a coding base, both starting at 1
pub fn codon_for_base(coding_base: i32) -> i32 {
    // Java's integer division rounds towards zero
    (coding_base - 1) / 3 + 1
}
