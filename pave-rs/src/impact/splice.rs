//! SpliceClassifier and hmf-common SpliceSites: splice region membership and splice donor and acceptor
//! effects.

use anyhow::Result;

use hmf_common::ensembl::{ExonData, TranscriptData};
use hmf_common::regions::position_within;

use super::effects::VariantEffect;
use super::{jsub, jsub_from, or_null, RefGenome, SpliceImpactType, VariantTransImpact, POS_STRAND};
use crate::constants::*;
use crate::variant::VariantData;

/// SpliceSites.getDonorPosition: D-1 is the last base of the exon, D+1 the first of the intron
pub fn get_donor_position(position: i32, exon_end: i32, strand: i8) -> i32 {
    if strand == POS_STRAND {
        if position > exon_end { position - exon_end } else { position - exon_end - 1 }
    } else if position < exon_end {
        exon_end - position
    } else {
        exon_end - position - 1
    }
}

/// SpliceSites.getAcceptorPosition: A-1 is the first base of the exon, A+1 the base before it
pub fn get_acceptor_position(position: i32, exon_end: i32, strand: i8) -> i32 {
    if strand == POS_STRAND {
        if position < exon_end { exon_end - position } else { exon_end - position - 1 }
    } else if position > exon_end {
        position - exon_end
    } else {
        position - exon_end - 1
    }
}

pub fn is_within_splice_region(variant: &VariantData, trans: &TranscriptData, exon: &ExonData) -> bool {
    let exon_start_region_start = exon.start - SPLICE_REGION_INTRON_RANGE;
    let exon_start_region_end = exon.start + SPLICE_REGION_EXON_RANGE - 1;
    let exon_end_region_start = exon.end - SPLICE_REGION_EXON_RANGE + 1;
    let exon_end_region_end = exon.end + SPLICE_REGION_INTRON_RANGE;

    let at_trans_start = exon.start == trans.trans_start;
    let at_trans_end = exon.end == trans.trans_end;

    if !at_trans_start && variant.alt_positions_overlap(exon_start_region_start, exon_start_region_end) {
        return true;
    }
    if !at_trans_end && variant.alt_positions_overlap(exon_end_region_start, exon_end_region_end) {
        return true;
    }
    false
}

pub fn check_straddles_splice_region(variant: &VariantData, trans: &TranscriptData, exon: &ExonData) -> Option<VariantEffect> {
    if !variant.is_deletion() {
        return None;
    }
    if exon.start != trans.trans_start && position_within(exon.start, variant.position, variant.end_position) {
        return Some(if trans.pos_strand() { VariantEffect::SpliceAcceptor } else { VariantEffect::SpliceDonor });
    } else if exon.end != trans.trans_end && position_within(exon.end, variant.position, variant.end_position) {
        return Some(if trans.pos_strand() { VariantEffect::SpliceDonor } else { VariantEffect::SpliceAcceptor });
    }
    None
}

/// SpliceClassifier.classifyVariant. The rules:
/// - a variant hitting splice donor sites (D-1, D+1, D+2, D+5) is a splice_donor_variant
/// - one hitting splice acceptor sites (A+1, A+2, A+3 if the alt is G, C on the negative strand) a
///   splice_acceptor_variant
/// - INDELs are judged on the splice bases they leave (determineIndelSpliceImpact)
pub fn classify_variant(variant: &VariantData, ti: &mut VariantTransImpact, exon: &ExonData, ref_genome: &dyn RefGenome) -> Result<()> {
    let trans = ti.trans;
    let pos_strand = trans.pos_strand();

    let donor_exon_pos = if pos_strand { exon.end } else { exon.start };
    let acceptor_exon_pos = if pos_strand { exon.start } else { exon.end };
    let is_donor_candidate = (variant.position - donor_exon_pos).abs() < (variant.position - acceptor_exon_pos).abs();

    let mut splice_effect: Option<VariantEffect> = None;
    let mut impact_type = SpliceImpactType::Unknown;
    let mut alt_positions: Vec<i32> = Vec::new();

    if variant.is_indel() {
        impact_type = determine_indel_splice_impact(variant, ref_genome, exon, pos_strand)?;
        if impact_type.is_disruptive() {
            splice_effect = Some(if is_donor_candidate { VariantEffect::SpliceDonor } else { VariantEffect::SpliceAcceptor });
        }
        // the affected bases, for annotation only
        if variant.is_deletion() {
            alt_positions.extend_from_slice(variant.alt_positions());
        } else {
            // the actual bases affected by the inserted bases are not captured
            alt_positions.push(variant.position);
        }
    } else {
        let alt = variant.alt.as_bytes();
        for i in 0..variant.ref_.len() {
            let position = variant.position + i as i32;
            let alt_base = alt[i];
            if is_donor_candidate {
                let donor_pos = get_donor_position(position, donor_exon_pos, trans.strand);
                if SPLICE_DONOR_POSITIONS.contains(&donor_pos) {
                    splice_effect = Some(VariantEffect::SpliceDonor);
                    alt_positions.push(position);
                }
            } else {
                let acceptor_pos = get_acceptor_position(position, acceptor_exon_pos, trans.strand);
                if SPLICE_ACCEPTOR_POSITIONS.contains(&acceptor_pos) {
                    if acceptor_pos == SPLICE_ACCEPTOR_END_RANGE {
                        let required_a3_base = if pos_strand { b'G' } else { b'C' };
                        if alt_base == required_a3_base {
                            splice_effect = Some(VariantEffect::SpliceAcceptor);
                            alt_positions.push(position);
                        }
                    } else {
                        splice_effect = Some(VariantEffect::SpliceAcceptor);
                        alt_positions.push(position);
                    }
                }
            }
        }
        if splice_effect.is_some() {
            impact_type = SpliceImpactType::BaseChange;
        }
    }

    if let Some(effect) = splice_effect {
        ti.add_effect(effect);
        ti.set_splice_impact_type(impact_type);

        // the splice bases affected
        let mut base_info: Vec<String> = Vec::new();
        for position in alt_positions {
            if is_donor_candidate {
                let donor_pos = get_donor_position(position, donor_exon_pos, trans.strand);
                if SPLICE_DONOR_POSITIONS.contains(&donor_pos) {
                    base_info.push(format!("D{donor_pos}"));
                }
            } else {
                let acceptor_pos = get_acceptor_position(position, acceptor_exon_pos, trans.strand);
                if SPLICE_ACCEPTOR_POSITIONS.contains(&acceptor_pos) {
                    base_info.push(format!("A{acceptor_pos}"));
                }
            }
        }
        ti.coding.splice_donor_acceptor_bases = base_info.join(";");
    }
    Ok(())
}

fn base_string(ref_genome: &dyn RefGenome, chr: &str, start: i32, end: i32) -> Result<String> {
    Ok(or_null(ref_genome.base_string(chr, start, end)?))
}

pub fn determine_indel_splice_impact(variant: &VariantData, ref_genome: &dyn RefGenome, exon: &ExonData, pos_strand: bool) -> Result<SpliceImpactType> {
    if variant.is_base_change() {
        return Ok(SpliceImpactType::BaseChange);
    }

    let donor_exon_pos = if pos_strand { exon.end } else { exon.start };
    let acceptor_exon_pos = if pos_strand { exon.start } else { exon.end };
    let is_donor_candidate = (variant.position - donor_exon_pos).abs() < (variant.position - acceptor_exon_pos).abs();
    let exon_boundary = if is_donor_candidate == pos_strand { exon.end } else { exon.start };

    // the offsets of A1-3 and D-1 to D5 relative to the exon boundary: 0 is the last base of the exon, -1
    // the first of the intron and 5 the 5th intronic base
    let (range_start, range_end) = if pos_strand {
        if is_donor_candidate { (0, 5) } else { (-3, -1) }
    } else if is_donor_candidate {
        (-5, 0)
    } else {
        (1, 3)
    };

    let pos_range_start = exon_boundary + range_start;
    let pos_range_end = exon_boundary + range_end;

    let mut ref_splice_bases = base_string(ref_genome, &variant.chromosome, pos_range_start, pos_range_end)?;

    // the alt bases for the same splice region
    let mut alt_splice_bases = String::new();
    let chr = &variant.chromosome;

    if variant.is_insert() {
        // checked already that the variant impacts the splice region, but check again
        if variant.position >= pos_range_end || variant.end_position <= pos_range_start {
            return Ok(SpliceImpactType::OutsideRange);
        }

        // inserting bases into the acceptor creating matching homology is left to the realigned variant
        if !is_donor_candidate {
            return Ok(SpliceImpactType::BaseShift);
        }

        let pre_insert_bases = variant.position - pos_range_start;
        if pre_insert_bases > 0 {
            alt_splice_bases = base_string(ref_genome, chr, pos_range_start, pos_range_start + pre_insert_bases - 1)?;
        }

        // the last exon base before the acceptor on the negative strand
        if !pos_strand && variant.position == exon_boundary && !is_donor_candidate {
            alt_splice_bases += jsub_from(&variant.alt, 1)?;
        } else {
            alt_splice_bases += &variant.alt;
        }

        // eg an insert whose ref is at 0 (the last exon base) of 2 bases needs the next 3 ref bases
        let next_ref_base = variant.position + variant.base_diff() + 1;
        let post_insert_bases = pos_range_end - next_ref_base;
        if post_insert_bases >= 0 {
            alt_splice_bases += &base_string(ref_genome, chr, variant.end_position, variant.end_position + post_insert_bases)?;
        }

        // too many bases: take the first or last ones to match the length of the ref
        if alt_splice_bases.len() > ref_splice_bases.len() {
            let n = ref_splice_bases.len() as i32;
            alt_splice_bases = if is_donor_candidate == pos_strand {
                jsub(&alt_splice_bases, 0, n)?.to_string()
            } else {
                jsub_from(&alt_splice_bases, alt_splice_bases.len() as i32 - n)?.to_string()
            };
        }
    } else {
        // no impact
        if variant.position >= pos_range_end || variant.end_position <= pos_range_start {
            return Ok(SpliceImpactType::OutsideRange);
        }
        // the entire region deleted
        if variant.position < pos_range_start && variant.end_position > pos_range_end {
            return Ok(SpliceImpactType::RegionDeleted);
        }
        // pulling upstream exonic bases to create matching homology is left to the realigned variant
        if variant.position < exon_boundary {
            return Ok(SpliceImpactType::BaseShift);
        }

        // otherwise the DEL partially overlaps the splice region and ref bases are pulled from one side
        if variant.position >= pos_range_start {
            // the start of the region is kept: the ref up to the first deleted base
            let pre_del_bases = variant.position - pos_range_start;
            if pre_del_bases >= 0 {
                alt_splice_bases = base_string(ref_genome, chr, variant.position - pre_del_bases, variant.position)?;
            }
            let post_del_bases = pos_range_end - variant.position;
            if post_del_bases > 0 {
                alt_splice_bases += &base_string(ref_genome, chr, variant.end_position, variant.end_position + post_del_bases - 1)?;
            }
        } else {
            // a realigned variant cannot then consider how exonic bases can create matching homology
            if variant.is_realigned_variant() {
                return Ok(SpliceImpactType::BaseChange);
            }
            // the DEL starts before the splice range and ends in it, so the alt bases come from lower positions
            let pre_del_bases = variant.end_position - pos_range_start;
            alt_splice_bases = base_string(ref_genome, chr, variant.position - (pre_del_bases - 1), variant.position)?;
            let post_del_bases = pos_range_end - variant.end_position;
            if post_del_bases >= 0 {
                alt_splice_bases += &base_string(ref_genome, chr, variant.end_position, variant.end_position + post_del_bases)?;
            }
        }
    }

    if alt_splice_bases.len() != ref_splice_bases.len() {
        crate::error!("splice base mismatch: ref({}) vs alt({})", ref_splice_bases, alt_splice_bases);
        return Ok(SpliceImpactType::Unknown);
    }

    if !is_donor_candidate {
        // the A3 base of the acceptor region only counts if it becomes a G (C on the negative strand)
        let r = ref_splice_bases.as_bytes();
        let a = alt_splice_bases.as_bytes();
        let ref_a3 = if pos_strand { r[0] } else { r[2] };
        let alt_a3 = if pos_strand { a[0] } else { a[2] };
        let required_a3 = if pos_strand { b'G' } else { b'C' };
        if ref_a3 == alt_a3 || alt_a3 != required_a3 {
            // this base is left out of the comparison
            ref_splice_bases = if pos_strand { jsub(&ref_splice_bases, 1, 3)?.to_string() } else { jsub(&ref_splice_bases, 0, 2)?.to_string() };
            alt_splice_bases = if pos_strand { jsub(&alt_splice_bases, 1, 3)?.to_string() } else { jsub(&alt_splice_bases, 0, 2)?.to_string() };
            if ref_splice_bases == alt_splice_bases {
                return Ok(SpliceImpactType::OutsideRange);
            }
        }
    }

    Ok(if ref_splice_bases == alt_splice_bases { SpliceImpactType::HomologyShift } else { SpliceImpactType::BaseShift })
}

/// whether the inserted bases are incorporated into the exon's coding sequence
pub fn is_insert_into_exon_start(variant: &VariantData, exon: &ExonData, pos_strand: bool) -> bool {
    if !variant.is_insert() {
        return false;
    }
    if pos_strand { variant.position == exon.start - 1 } else { variant.position == exon.end }
}
