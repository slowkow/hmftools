//! CodingUtils: a variant's coding context in a transcript (region, coding type, exon, coding base,
//! phase, distance to the nearest exon, frameshift).

use anyhow::Result;

use hmf_common::ensembl::{ExonData, TranscriptData};
use hmf_common::regions::{position_within, positions_overlap};

use super::effects::{CodingType, RegionType};
use super::protein::get_open_codon_bases;
use super::splice::is_insert_into_exon_start;
use super::{coding_end, coding_start, CodingContext, PHASE_0, PHASE_1, PHASE_NONE, POS_STRAND, SE_END, SE_START};
use crate::constants::{PROMOTOR_UPSTREAM_DISTANCE, PROMOTOR_UPSTREAM_GENE_IDS};
use crate::variant::VariantData;

/// TranscriptUtils.calcExonicCodingPhase
pub fn calc_exonic_coding_phase(exon: &ExonData, coding_start: i32, coding_end: i32, strand: i8, position: i32) -> i32 {
    if !position_within(position, coding_start, coding_end) {
        return PHASE_NONE;
    }
    // coding bases up to and including the position
    let coding_bases =
        if strand == POS_STRAND { position - exon.start.max(coding_start) + 1 } else { exon.end.min(coding_end) - position + 1 };
    if (strand == POS_STRAND && coding_start == exon.start) || (strand != POS_STRAND && coding_end == exon.end) {
        // coding starts on the first base of this exon: phasing starts on 1 if not specified
        let start_phase = if exon.phase_start == PHASE_NONE { PHASE_0 } else { exon.phase_start };
        (start_phase + coding_bases) % 3
    } else if (strand == POS_STRAND && coding_start < exon.start) || (strand != POS_STRAND && coding_end > exon.end) {
        // coding started before this exon
        (exon.phase_start + coding_bases) % 3
    } else {
        // coding starts within the exon
        coding_bases % 3
    }
}

pub fn determine_context(variant: &VariantData, trans: &TranscriptData, cc: &mut CodingContext) -> Result<()> {
    cc.strand = trans.strand;
    set_coding_context(variant, trans, cc)?;

    // pre-coding exonic bases are set separately since they don't fit the way exons are iterated through
    if matches!(cc.coding_type, CodingType::Utr5p | CodingType::Utr3p | CodingType::Enhancer) {
        set_utr_coding_exonic_distance(variant, trans, cc)?;
    }

    if variant.is_deletion() && cc.region_type == RegionType::Exonic {
        // trimmed deleted bases within exons
        let r = cc.coding_ref(variant)?;
        let a = cc.coding_alt(variant)?;
        cc.deleted_coding_bases = r.len() as i32 - a.len() as i32;
    }

    // mark frameshift INDELs now that the coding bases are clarified
    if variant.is_indel() && cc.coding_type == CodingType::Coding && cc.region_type == RegionType::Exonic {
        let adjusted_coding_bases = if variant.is_insert() {
            variant.base_diff()
        } else {
            if cc.nearest_exon_distance != 0 {
                // used to set the HGVS coding string
                if trans.pos_strand() && variant.alt_positions().contains(&coding_end(trans)?) {
                    cc.spans_coding_end = true;
                } else if !trans.pos_strand() && variant.alt_positions().contains(&coding_start(trans)?) {
                    cc.spans_coding_end = true;
                }
            }
            cc.deleted_coding_bases
        };
        if adjusted_coding_bases % 3 != 0 {
            cc.is_frameshift = true;
        }
    }
    Ok(())
}

fn set_coding_context(variant: &VariantData, trans: &TranscriptData, cc: &mut CodingContext) -> Result<()> {
    let mut in_coding = false;
    let is_coding_trans = !trans.non_coding();
    let coding_start = trans.coding_start.unwrap_or(0);
    let coding_end = trans.coding_end.unwrap_or(0);

    if trans.non_coding() {
        cc.coding_type = CodingType::NonCoding;
    }

    if variant.alt_bases_below(trans.trans_start) || variant.alt_bases_above(trans.trans_end) {
        cc.region_type = RegionType::Upstream;
        if !trans.non_coding() && PROMOTOR_UPSTREAM_GENE_IDS.contains(&trans.gene_id.as_str()) {
            cc.coding_type = CodingType::Enhancer;
        }
        return Ok(());
    }

    let pos_strand = trans.pos_strand();

    if is_coding_trans {
        let below_coding_start = variant.alt_bases_below(coding_start);
        let above_coding_end = variant.alt_bases_above(coding_end);
        if below_coding_start {
            cc.coding_type = if pos_strand { CodingType::Utr5p } else { CodingType::Utr3p };
        } else if above_coding_end {
            cc.coding_type = if pos_strand { CodingType::Utr3p } else { CodingType::Utr5p };
        } else {
            in_coding = true;
            cc.coding_type = CodingType::Coding;
        }
    }

    let mut pre_pos_exonic_base_count = 0; // for non-coding transcripts only
    let mut coding_base_count = 0; // accumulated since the start of coding before the exon closest to the variant

    let upstream_start_pos = if pos_strand { variant.position } else { variant.end_position };

    // INDELs don't affect the last base of an exon, so need other criteria to skip it
    let skip_exon_pos = if (!variant.is_indel()) == pos_strand { variant.position } else { variant.end_position };

    let exons = &trans.exons;
    let n = exons.len();

    if pos_strand {
        let mut i = 0;
        while i < n {
            let exon = &exons[i];
            if i == 0 && is_coding_trans && coding_start == exon.start {
                // a rare incomplete transcript whose coding start phase is not 1, usually at the first base of exon 1
                let start_phase = calc_exonic_coding_phase(exon, coding_start, coding_end, trans.strand, exon.start);
                if start_phase != PHASE_1 {
                    coding_base_count -= 3 - get_open_codon_bases(start_phase);
                }
            }
            let next_exon = if i < n - 1 { Some(&exons[i + 1]) } else { None };

            let within_exon = variant.alt_positions_within(exon.start, exon.end) || is_insert_into_exon_start(variant, exon, pos_strand);
            let overlaps_exon = variant.alt_positions_overlap(exon.start, exon.end);

            if skip_exon_pos > exon.end {
                // keep track of coding bases if this exon overlaps the coding region but is before the variant
                if in_coding && !within_exon && !overlaps_exon && positions_overlap(exon.start, exon.end, coding_start, coding_end) {
                    coding_base_count += coding_end.min(exon.end) - coding_start.max(exon.start) + 1;
                }
                if variant.position > exon.end {
                    pre_pos_exonic_base_count += exon.base_length();
                }
                if let Some(ne) = next_exon {
                    if upstream_start_pos >= ne.start {
                        i += 1;
                        continue;
                    }
                }
            }

            if within_exon || overlaps_exon {
                cc.region_type = RegionType::Exonic;
                cc.exon_rank = exon.rank;
                cc.spans_splice_junction = !within_exon;
                pre_pos_exonic_base_count += upstream_start_pos - exon.start + 1;

                if in_coding {
                    cc.upstream_phase = calc_exonic_coding_phase(exon, coding_start, coding_end, trans.strand, upstream_start_pos.max(exon.start));
                    cc.coding_position_range[SE_START] = exon.start.max(coding_start).max(variant.position);
                    cc.coding_position_range[SE_END] = if !variant.is_insert() {
                        exon.end.min(coding_end).min(variant.end_position)
                    } else {
                        cc.coding_position_range[SE_START]
                    };
                    cc.coding_base = coding_base_count + cc.coding_position_range[SE_START] - coding_start.max(exon.start) + 1;

                    // distance of a variant into the intron if applicable
                    if !variant.is_insert() {
                        let exon_coding_start = exon.start.max(coding_start);
                        let exon_coding_end = exon.end.min(coding_end);
                        if variant.end_position > exon_coding_end {
                            cc.nearest_exon_distance = variant.end_position - exon_coding_end;
                        } else if variant.position < exon_coding_start {
                            cc.nearest_exon_distance = variant.position - exon_coding_start; // -ve back into the previous intron
                        }
                        cc.spans_coding_start = variant.position < coding_start;
                    }
                } else {
                    // the exonic range, for variants spanning from exon to intron too
                    cc.coding_position_range[SE_START] = exon.start.max(variant.position);
                    cc.coding_position_range[SE_END] =
                        if !variant.is_insert() { exon.end.min(variant.end_position) } else { cc.coding_position_range[SE_START] };
                    if !variant.is_insert() {
                        if variant.end_position > exon.end {
                            cc.nearest_exon_distance = variant.end_position - exon.end;
                        } else if variant.position < exon.start {
                            cc.nearest_exon_distance = variant.position - exon.start;
                        }
                    }
                }
            } else if next_exon.is_some_and(|ne| variant.alt_positions_overlap(ne.start, ne.end))
                || insert_into_exon_start_or_npe(variant, next_exon, pos_strand)?
            {
                // the next exon is dealt with the next time around
                i += 1;
                continue;
            } else {
                cc.region_type = RegionType::Intronic;
                let ne = next_exon.ok_or_else(|| anyhow::anyhow!("NullPointerException: no next exon for {variant}"))?;
                let distance_to_prev = upstream_start_pos - exon.end;
                let distance_to_next = ne.start - upstream_start_pos;
                if distance_to_prev < distance_to_next {
                    cc.exon_rank = exon.rank;
                    cc.nearest_exon_distance = upstream_start_pos - exon.end;
                } else {
                    cc.exon_rank = ne.rank;
                    cc.nearest_exon_distance = upstream_start_pos - ne.start; // -ve for back into the previous intron
                }
                if in_coding {
                    cc.upstream_phase = exon.phase_end;
                    cc.coding_base = if distance_to_prev < distance_to_next { coding_base_count } else { coding_base_count + 1 };
                }
            }
            break;
        }
    } else {
        // negative strand transcripts
        let mut i = n as i64 - 1;
        while i >= 0 {
            let exon = &exons[i as usize];
            if exon.rank == 1 && is_coding_trans && coding_end == exon.end {
                let start_phase = calc_exonic_coding_phase(exon, coding_start, coding_end, trans.strand, exon.end);
                if start_phase != PHASE_1 {
                    coding_base_count -= 3 - get_open_codon_bases(start_phase);
                }
            }
            let next_exon = if i >= 1 { Some(&exons[i as usize - 1]) } else { None };

            let within_exon = variant.alt_positions_within(exon.start, exon.end) || is_insert_into_exon_start(variant, exon, pos_strand);
            let overlaps_exon = variant.alt_positions_overlap(exon.start, exon.end);

            if skip_exon_pos < exon.start {
                if in_coding && !within_exon && !overlaps_exon && positions_overlap(exon.start, exon.end, coding_start, coding_end) {
                    coding_base_count += coding_end.min(exon.end) - coding_start.max(exon.start) + 1;
                }
                if variant.end_position < exon.start {
                    pre_pos_exonic_base_count += exon.base_length();
                }
                if let Some(ne) = next_exon {
                    if upstream_start_pos <= ne.end {
                        i -= 1;
                        continue;
                    }
                }
            }

            if within_exon || overlaps_exon {
                cc.region_type = RegionType::Exonic;
                cc.exon_rank = exon.rank;
                cc.spans_splice_junction = !within_exon;
                pre_pos_exonic_base_count += exon.end - upstream_start_pos + 1;

                if in_coding {
                    cc.upstream_phase = calc_exonic_coding_phase(exon, coding_start, coding_end, trans.strand, upstream_start_pos.min(exon.end));
                    if !variant.is_insert() {
                        cc.coding_position_range[SE_START] = exon.start.max(coding_start).max(variant.position);
                        cc.coding_position_range[SE_END] = exon.end.min(coding_end).min(variant.end_position);
                    } else {
                        cc.coding_position_range[SE_END] = exon.end.min(coding_end).min(variant.end_position);
                        cc.coding_position_range[SE_START] = cc.coding_position_range[SE_END];
                    }
                    cc.coding_base = coding_base_count + coding_end.min(exon.end) - cc.coding_position_range[SE_END] + 1;

                    if !variant.is_insert() {
                        let exon_coding_start = exon.start.max(coding_start);
                        let exon_coding_end = exon.end.min(coding_end);
                        if variant.end_position > exon_coding_end {
                            cc.nearest_exon_distance = exon_coding_end - variant.end_position; // -ve back into the previous intron
                        } else if variant.position < exon_coding_start {
                            cc.nearest_exon_distance = exon_coding_start - variant.position;
                        }
                        cc.spans_coding_start = variant.end_position > coding_end;
                    }
                } else if !variant.is_insert() {
                    cc.coding_position_range[SE_START] = exon.start.max(variant.position);
                    cc.coding_position_range[SE_END] = exon.end.min(variant.end_position);
                    if variant.end_position > exon.end {
                        cc.nearest_exon_distance = exon.end - variant.end_position;
                    } else if variant.position < exon.start {
                        cc.nearest_exon_distance = exon.start - variant.position;
                    }
                } else {
                    cc.coding_position_range[SE_END] = exon.end.min(variant.end_position);
                    cc.coding_position_range[SE_START] = cc.coding_position_range[SE_END];
                }
            } else if next_exon.is_some_and(|ne| variant.alt_positions_overlap(ne.start, ne.end))
                || insert_into_exon_start_or_npe(variant, next_exon, pos_strand)?
            {
                i -= 1;
                continue;
            } else {
                cc.region_type = RegionType::Intronic;
                let ne = next_exon.ok_or_else(|| anyhow::anyhow!("NullPointerException: no next exon for {variant}"))?;
                let distance_to_prev = exon.start - upstream_start_pos;
                let distance_to_next = upstream_start_pos - ne.end;
                if distance_to_prev < distance_to_next {
                    cc.exon_rank = exon.rank;
                    cc.nearest_exon_distance = exon.start - upstream_start_pos;
                } else {
                    cc.exon_rank = ne.rank;
                    cc.nearest_exon_distance = ne.end - upstream_start_pos; // negative
                }
                if in_coding {
                    cc.upstream_phase = ne.phase_end;
                    cc.coding_base = if distance_to_prev < distance_to_next { coding_base_count } else { coding_base_count + 1 };
                }
            }
            break;
        }
    }

    if cc.coding_type == CodingType::NonCoding {
        cc.coding_base = pre_pos_exonic_base_count;
        if cc.nearest_exon_distance < 0 {
            cc.coding_base += 1; // to the next exonic base
        }
    }
    Ok(())
}

/// `(nextExon != null && overlaps) || isInsertIntoExonStart(variant, nextExon, ...)` in Java: the second
/// test runs on a null next exon too, which throws for an insert (a non-insert returns before touching it)
fn insert_into_exon_start_or_npe(variant: &VariantData, next_exon: Option<&ExonData>, pos_strand: bool) -> Result<bool> {
    match next_exon {
        Some(ne) => Ok(is_insert_into_exon_start(variant, ne, pos_strand)),
        None if variant.is_insert() => anyhow::bail!("NullPointerException: no next exon for {variant}"),
        None => Ok(false),
    }
}

/// For positions in the UTR, the exonic distance between the position and the start or end of the
/// coding region
pub fn set_utr_coding_exonic_distance(variant: &VariantData, trans: &TranscriptData, cc: &mut CodingContext) -> Result<()> {
    let pos_strand = trans.pos_strand();
    let pos_before_coding_start = pos_strand == (cc.coding_type == CodingType::Utr5p || cc.coding_type == CodingType::Enhancer);

    let coding_start = coding_start(trans)?;
    let coding_end = coding_end(trans)?;
    let position = if pos_strand { variant.position } else { variant.end_position };
    let mut total_coding_bases = 0;
    let mut coding_starts_ends_on_exon_boundary = false;

    for exon in &trans.exons {
        if pos_before_coding_start {
            if exon.start > coding_start {
                break;
            }
            if position <= exon.end {
                // any exonic bases between the position and coding
                if position_within(coding_start, exon.start, exon.end) {
                    cc.coding_base += coding_start - position.max(exon.start);
                } else {
                    cc.coding_base += exon.end - position.max(exon.start) + 1;
                }
            }
        } else {
            if exon.start > position {
                break;
            }
            if coding_end <= exon.end {
                if position_within(coding_end, exon.start, exon.end) {
                    cc.coding_base += position.min(exon.end) - coding_end;
                } else {
                    cc.coding_base += position.min(exon.end) - exon.start + 1;
                }
            }
        }
        if positions_overlap(coding_start, coding_end, exon.start, exon.end) {
            total_coding_bases += coding_end.min(exon.end) - coding_start.max(exon.start) + 1;
        }
        if coding_end == exon.end || coding_start == exon.start {
            coding_starts_ends_on_exon_boundary = true;
        }
    }

    // push the base by 1 if intronic and closest to the next exon
    if cc.region_type == RegionType::Intronic {
        if cc.coding_type == CodingType::Utr5p && cc.nearest_exon_distance > 0 {
            cc.coding_base += 1;
        } else if cc.coding_type == CodingType::Utr3p && cc.nearest_exon_distance < 0 {
            cc.coding_base += 1;
        }
    }

    if cc.coding_base == 0 {
        if coding_starts_ends_on_exon_boundary {
            cc.coding_ends_on_exon_boundary = true;
        }
        cc.coding_base = if cc.coding_type == CodingType::Utr5p { 1 } else { total_coding_bases };
    }

    if cc.coding_type == CodingType::Enhancer {
        cc.coding_base = if pos_strand { coding_start - position } else { position - coding_end };
        // excluded if too far upstream
        if cc.coding_base > PROMOTOR_UPSTREAM_DISTANCE {
            cc.coding_base = 0;
            cc.coding_type = CodingType::Unknown;
        }
    }
    Ok(())
}
