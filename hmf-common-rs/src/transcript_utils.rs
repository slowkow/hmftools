//! hmf-common `gene.TranscriptUtils` and `CodingBaseData`: coding bases and phases of a position in a
//! transcript.

use crate::ensembl::{ExonData, TranscriptData};

pub const PHASE_NONE: i32 = -1;
pub const PHASE_0: i32 = 0;
pub const PHASE_1: i32 = 1;
pub const PHASE_2: i32 = 2;

/// CodingBaseData
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CodingBaseData {
    pub total_coding_bases: i32,
    pub coding_bases: i32,
    pub phase: i32,
}

impl Default for CodingBaseData {
    fn default() -> Self {
        CodingBaseData { total_coding_bases: 0, coding_bases: 0, phase: PHASE_NONE }
    }
}

fn position_within(pos: i32, start: i32, end: i32) -> bool {
    pos >= start && pos <= end
}

/// TranscriptUtils.codingBasesToPhase
pub fn coding_bases_to_phase(mut coding_bases: i32, start_phase: i32) -> i32 {
    if start_phase == PHASE_2 {
        coding_bases += 1;
    } else if start_phase == PHASE_0 {
        coding_bases += 2;
    }
    if coding_bases >= 0 {
        return coding_bases % 3;
    }
    let m = coding_bases.abs() % 3;
    if m == 0 {
        return 0;
    }
    if m == 1 { 2 } else { 1 }
}

/// TranscriptUtils.tickPhaseForward(phase, increment)
pub fn tick_phase_forward(phase: i32, increment: i32) -> i32 {
    if increment != 0 { (phase + increment) % 3 } else { phase }
}

/// TranscriptUtils.calcExonicCodingPhase
pub fn calc_exonic_coding_phase(exon: &ExonData, coding_start: i32, coding_end: i32, strand: i8, position: i32) -> i32 {
    if !position_within(position, coding_start, coding_end) {
        return PHASE_NONE;
    }
    let coding_bases =
        if strand == 1 { position - exon.start.max(coding_start) + 1 } else { exon.end.min(coding_end) - position + 1 };
    if (strand == 1 && coding_start == exon.start) || (strand == -1 && coding_end == exon.end) {
        let start_phase = if exon.phase_start == PHASE_NONE { PHASE_0 } else { exon.phase_start };
        (start_phase + coding_bases) % 3
    } else if (strand == 1 && coding_start < exon.start) || (strand == -1 && coding_end > exon.end) {
        (exon.phase_start + coding_bases) % 3
    } else {
        coding_bases % 3
    }
}

/// TranscriptUtils.codingBaseLength
pub fn coding_base_length(trans: &TranscriptData) -> i32 {
    let (Some(cs), Some(ce)) = (trans.coding_start, trans.coding_end) else { return 0 };
    let mut bases = 0;
    for exon in &trans.exons {
        if cs > exon.end {
            continue;
        }
        if ce < exon.start {
            break;
        }
        bases += exon.end.min(ce) - exon.start.max(cs) + 1;
    }
    bases
}

/// TranscriptUtils.calcCodingBases: coding bases up to the position in the direction of transcription,
/// the transcript's total, and the position's phase
pub fn calc_coding_bases(trans: &TranscriptData, position: i32) -> CodingBaseData {
    let mut cb = CodingBaseData::default();
    let (Some(coding_start), Some(coding_end)) = (trans.coding_start, trans.coding_end) else { return cb };
    let mut in_coding_region = false;
    let mut coding_region_ended = false;
    let pos_in_coding_region = position_within(position, coding_start, coding_end);
    let mut is_exonic = false;
    let exons = &trans.exons;
    for i in 0..exons.len() {
        let exon = &exons[i];
        let next = exons.get(i + 1);
        let (exon_start, exon_end) = (exon.start, exon.end);
        if !is_exonic && position_within(position, exon_start, exon_end) {
            is_exonic = true;
            if pos_in_coding_region {
                cb.phase = calc_exonic_coding_phase(exon, coding_start, coding_end, trans.strand, position);
            }
        }
        if let Some(n) = next {
            if position > exon_end && position < n.start {
                cb.phase = if trans.strand == 1 { exon.phase_end } else { n.phase_end };
            }
        }
        if !in_coding_region {
            if exon_end >= coding_start {
                in_coding_region = true;
                cb.total_coding_bases += exon_end - coding_start + 1;
                if position >= coding_start {
                    if position < exon_end {
                        cb.coding_bases += position - coding_start + 1;
                    } else {
                        cb.coding_bases += exon_end - coding_start + 1;
                    }
                }
            }
        } else if !coding_region_ended {
            if exon_start > coding_end {
                coding_region_ended = true;
            } else if exon_end >= coding_end {
                coding_region_ended = true;
                cb.total_coding_bases += coding_end - exon_start + 1;
                if position >= exon_start {
                    if position < coding_end {
                        cb.coding_bases += position - exon_start + 1;
                    } else {
                        cb.coding_bases += coding_end - exon_start + 1;
                    }
                }
            } else {
                cb.total_coding_bases += exon_end - exon_start + 1;
                if position >= exon_start {
                    if position < exon_end {
                        cb.coding_bases += position - exon_start + 1;
                    } else {
                        cb.coding_bases += exon_end - exon_start + 1;
                    }
                }
            }
            if coding_region_ended {
                break;
            }
        }
    }
    if trans.strand == -1 {
        if pos_in_coding_region && is_exonic {
            cb.coding_bases = cb.total_coding_bases - cb.coding_bases + 1;
        } else {
            cb.coding_bases = cb.total_coding_bases - cb.coding_bases;
        }
    }
    if !pos_in_coding_region {
        cb.phase = PHASE_NONE;
    }
    cb
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trans(strand: i8, exons: &[(i32, i32, i32, i32)], cs: i32, ce: i32) -> TranscriptData {
        TranscriptData {
            trans_id: 1,
            trans_name: "T".into(),
            gene_id: "G".into(),
            is_canonical: true,
            strand,
            trans_start: exons[0].0,
            trans_end: exons[exons.len() - 1].1,
            coding_start: Some(cs),
            coding_end: Some(ce),
            bio_type: "protein_coding".into(),
            ref_seq_id: None,
            exons: exons
                .iter()
                .enumerate()
                .map(|(i, e)| ExonData { trans_id: 1, start: e.0, end: e.1, rank: if strand == 1 { i as i32 + 1 } else { (exons.len() - i) as i32 }, phase_start: e.2, phase_end: e.3 })
                .collect(),
        }
    }

    #[test]
    fn coding_bases_forward_and_reverse() {
        // exons 100-200, 300-400, 500-600; coding 150-550
        let t = trans(1, &[(100, 200, -1, 0), (300, 400, 0, 2), (500, 600, 2, -1)], 150, 550);
        let cb = calc_coding_bases(&t, 350);
        assert_eq!(cb.total_coding_bases, 51 + 101 + 51);
        assert_eq!(cb.coding_bases, 51 + 51);
        // intronic: the preceding exon's end phase
        let cb = calc_coding_bases(&t, 250);
        assert_eq!((cb.coding_bases, cb.phase), (51, 0));
        assert_eq!(calc_coding_bases(&t, 50).phase, PHASE_NONE);
        let r = trans(-1, &[(100, 200, -1, -1), (300, 400, 2, 0), (500, 600, -1, 2)], 150, 550);
        let cb = calc_coding_bases(&r, 350);
        assert_eq!(cb.total_coding_bases, 203);
        assert_eq!(cb.coding_bases, 203 - (51 + 51) + 1);
        assert_eq!(coding_base_length(&t), 203);
        assert_eq!(coding_bases_to_phase(-1, PHASE_NONE), 2);
        assert_eq!(tick_phase_forward(2, 2), 1);
    }
}
