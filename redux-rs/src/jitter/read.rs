//! One read against one microsatellite site (Java MicrosatelliteRead and ConsensusMarker).

use rust_htslib::bam::{self, record::Aux};

pub const MIN_FLANKING_BASE_MATCHES: i32 = 5;
pub const LOW_BASE_QUAL_FLANKING_BASES: i32 = 1;

/// BaseQualAdjustment.BASE_QUAL_MINIMUM
const BASE_QUAL_MINIMUM: i8 = 1;
/// SbxBamUtils.SBX_MEDIUM_QUAL_LOWER / UPPER
const SBX_MEDIUM_QUAL_LOWER: i8 = 27;
const SBX_MEDIUM_QUAL_UPPER: i8 = 29;

const CONSENSUS_TYPE_ATTRIBUTE: &[u8] = b"UT";
const SBX_DUPLEX_READ_INDEX_TAG: &[u8] = b"YX";

// BAM CIGAR operator codes
const OP_M: u32 = 0;
const OP_I: u32 = 1;
const OP_D: u32 = 2;
const OP_N: u32 = 3;
const OP_S: u32 = 4;
const OP_H: u32 = 5;
const OP_EQ: u32 = 7;
const OP_X: u32 = 8;

#[inline]
fn consumes_read(op: u32) -> bool {
    matches!(op, OP_M | OP_I | OP_S | OP_EQ | OP_X)
}
#[inline]
fn consumes_ref(op: u32) -> bool {
    matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X)
}

/// Reference bases covered by a raw BAM CIGAR.
#[inline]
pub fn reference_length(cigar: &[u32]) -> i32 {
    let mut len = 0;
    for &c in cigar {
        if consumes_ref(c & 0xf) {
            len += (c >> 4) as i32;
        }
    }
    len
}

/// common.bam.ConsensusType, in the enum's order.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(u8)]
pub enum ConsensusType {
    None = 0,
    Single = 1,
    Dual = 2,
}

impl ConsensusType {
    pub const ALL: [ConsensusType; 3] = [ConsensusType::None, ConsensusType::Single, ConsensusType::Dual];
    pub fn name(self) -> &'static str {
        match self {
            ConsensusType::None => "NONE",
            ConsensusType::Single => "SINGLE",
            ConsensusType::Dual => "DUAL",
        }
    }
}

/// The fields of MicrosatelliteSite that the read analysis uses.
#[derive(Clone, Copy)]
pub struct SiteRef<'a> {
    pub start: i32,
    pub end: i32,
    pub unit: &'a [u8],
    pub repeat_count: i32,
}

/// The result fields of MicrosatelliteRead.
#[derive(Clone, Copy, Debug, Default)]
pub struct MsRead {
    pub valid: bool,
    pub has_medium_quals: bool,
    pub repeat_length: i32,
    pub repeat_units: i32,
    pub jitter: i32,
}

/// MicrosatelliteRead.analyse without the consensus type. `read_start` and `read_end` are the 1-based alignment start
/// and end, `cigar` the raw BAM CIGAR.
#[inline]
pub fn analyse(site: &SiteRef, rec: &bam::Record, read_start: i32, read_end: i32, cigar: &[u32], apply_bq_filter: bool) -> MsRead {
    let mut out = MsRead::default();
    let ms_start = site.start;
    let ms_end = site.end;

    // the read must have aligned bases around the repeat including the flanks
    if read_start > ms_start - MIN_FLANKING_BASE_MATCHES || read_end < ms_end + MIN_FLANKING_BASE_MATCHES {
        return out;
    }

    let mut aligned_bases = 0;
    let mut inserted_bases = 0;
    let mut aligned_before = 0;
    let mut aligned_after = 0;
    let mut ms_read_index_start = -1;
    let mut read_index: i32 = 0;
    let mut ref_pos = read_start;

    for &c in cigar {
        let op = c & 0xf;
        let len = (c >> 4) as i32;
        let end_ref_pos = ref_pos + len - 1;

        // Java applies this to every operator, using the element length even when it does not consume the reference
        if ms_read_index_start < 0 && ref_pos <= ms_start && end_ref_pos >= ms_start {
            ms_read_index_start = read_index + ms_start - ref_pos;
        }

        if op == OP_D {
            // a delete covering the start or the end of the repeat, or just before its start
            if (ref_pos < ms_start && end_ref_pos >= ms_start - 1) || (ref_pos <= ms_end + 1 && end_ref_pos > ms_end) {
                return out;
            }
        } else if op == OP_I {
            if ref_pos >= ms_start && ref_pos <= ms_end + 1 {
                // the insert must be a repetition of the unit
                let seq = rec.seq();
                let unit_len = site.unit.len();
                for j in 0..len as usize {
                    if seq[read_index as usize + j] != site.unit[j % unit_len] {
                        return out;
                    }
                }
                inserted_bases += len;
            } else if ref_pos >= ms_start - MIN_FLANKING_BASE_MATCHES && ref_pos <= ms_end + MIN_FLANKING_BASE_MATCHES + 1 {
                // no insert within the flanks
                return out;
            }
        } else if op == OP_M {
            // Java counts only M here, not = or X
            let to_ms = (end_ref_pos + 1).min(ms_end + 1) - ref_pos.max(ms_start);
            if to_ms > 0 {
                aligned_bases += to_ms;
            }
            let to_flank_start = (end_ref_pos + 1).min(ms_start) - ref_pos.max(ms_start - MIN_FLANKING_BASE_MATCHES);
            if to_flank_start > 0 {
                aligned_before += to_flank_start;
            }
            let to_flank_end = (end_ref_pos + 1).min(ms_end + 1 + MIN_FLANKING_BASE_MATCHES) - ref_pos.max(ms_end + 1);
            if to_flank_end > 0 {
                aligned_after += to_flank_end;
            }
        }

        if consumes_read(op) {
            read_index += len;
        }
        if consumes_ref(op) {
            ref_pos += len;
        }
    }

    if aligned_before < MIN_FLANKING_BASE_MATCHES || aligned_after < MIN_FLANKING_BASE_MATCHES {
        return out;
    }
    if ms_read_index_start < 0 {
        return out;
    }
    let ms_read_index_end = ms_read_index_start + ms_end - ms_start;

    if apply_bq_filter && !has_valid_base_qualities(rec, ms_read_index_start, ms_read_index_end, &mut out.has_medium_quals) {
        return out;
    }

    out.valid = true;
    out.repeat_length = aligned_bases + inserted_bases;
    out.repeat_units = out.repeat_length / site.unit.len() as i32;
    out.jitter = out.repeat_units - site.repeat_count;
    out
}

/// MicrosatelliteRead.hasValidBaseQualities
#[inline]
fn has_valid_base_qualities(rec: &bam::Record, ms_read_index_start: i32, ms_read_index_end: i32, has_medium: &mut bool) -> bool {
    let flank_start = ms_read_index_start - LOW_BASE_QUAL_FLANKING_BASES;
    let flank_end = ms_read_index_end + LOW_BASE_QUAL_FLANKING_BASES;
    let quals = rec.qual();
    // htsjdk returns an empty array when the qualities are missing (first byte 0xFF)
    let qual_len = if quals.first() == Some(&0xff) { 0 } else { quals.len() as i32 };
    if flank_start < 0 || flank_end >= qual_len {
        return false;
    }
    for &q in &quals[flank_start as usize..=flank_end as usize] {
        // Java compares signed bytes
        let q = q as i8;
        if q <= BASE_QUAL_MINIMUM {
            return false;
        }
        *has_medium |= (SBX_MEDIUM_QUAL_LOWER..=SBX_MEDIUM_QUAL_UPPER).contains(&q);
    }
    true
}

/// SamRecordUtils.extractConsensusType. Java throws on a value other than the three names; that is read as NONE here.
#[inline]
pub fn extract_consensus_type(rec: &bam::Record) -> ConsensusType {
    match rec.aux(CONSENSUS_TYPE_ATTRIBUTE) {
        Ok(Aux::String(s)) => match s {
            "SINGLE" => ConsensusType::Single,
            "DUAL" => ConsensusType::Dual,
            _ => ConsensusType::None,
        },
        _ => ConsensusType::None,
    }
}

/// SbxBamUtils.extractDuplexBaseIndex
fn extract_duplex_base_index(rec: &bam::Record) -> i32 {
    match rec.aux(SBX_DUPLEX_READ_INDEX_TAG) {
        Ok(Aux::I8(v)) => v as i32,
        Ok(Aux::U8(v)) => v as i32,
        Ok(Aux::I16(v)) => v as i32,
        Ok(Aux::U16(v)) => v as i32,
        Ok(Aux::I32(v)) => v,
        Ok(Aux::U32(v)) => v as i32,
        _ => -1,
    }
}

/// SbxBamUtils.inDuplexRegion
fn in_duplex_region(pos_orientation_read: bool, duplex_base_index: i32, base_index: i32) -> bool {
    if duplex_base_index < 0 {
        return false;
    }
    if pos_orientation_read { base_index >= duplex_base_index } else { base_index <= duplex_base_index }
}

/// ConsensusMarker.getMicrosatelliteBoundaries: the read indices just outside the repeat.
pub fn microsatellite_boundaries(ms_start: i32, ms_end: i32, read_start: i32, cigar: &[u32]) -> (i32, i32) {
    const INVALID_INDEX: i32 = -1;
    let mut cigar_index = 0;
    let mut read_index: i32 = 0;
    if let Some(&first) = cigar.first() {
        if first & 0xf == OP_S {
            cigar_index += 1;
            read_index += (first >> 4) as i32;
        } else if first & 0xf == OP_H {
            cigar_index += 1;
        }
    }
    let mut ref_pos = read_start;
    let mut start_read_index = INVALID_INDEX;
    let mut end_read_index = INVALID_INDEX;

    for &c in &cigar[cigar_index.min(cigar.len())..] {
        let op = c & 0xf;
        let len = (c >> 4) as i32;
        if op == OP_S || op == OP_H {
            continue;
        }
        let is_read = consumes_read(op);
        let is_ref = consumes_ref(op);
        if is_read && is_ref {
            let end_ref_pos = ref_pos + len - 1;
            if start_read_index == INVALID_INDEX && end_ref_pos >= ms_start {
                start_read_index = if ref_pos >= ms_start { read_index } else { read_index + ms_start - ref_pos };
            }
            if end_ref_pos <= ms_end {
                end_read_index = read_index + len - 1;
            } else if ref_pos <= ms_end {
                end_read_index = read_index + ms_end - ref_pos;
            }
            read_index += len;
            ref_pos += len;
        } else if is_read {
            if start_read_index == INVALID_INDEX && ref_pos - 1 >= ms_start {
                start_read_index = read_index;
            }
            if ref_pos - 1 <= ms_end {
                end_read_index = read_index + len - 1;
            }
            read_index += len;
        } else if is_ref {
            ref_pos += len;
        }
        // Java throws IllegalStateException for an operator that consumes neither (P); it is skipped here
    }
    (start_read_index - 1, end_read_index + 1)
}

/// ConsensusMarker.SBXConsensusMarker.consensusType for a read whose UT attribute is `base`: a duplex read counts
/// as DUAL only when the repeat lies in its duplex part.
#[inline]
pub fn sbx_consensus_type(base: ConsensusType, site: &SiteRef, rec: &bam::Record, read_start: i32, cigar: &[u32]) -> ConsensusType {
    if base != ConsensusType::Dual {
        return base;
    }
    let duplex_base_index = extract_duplex_base_index(rec);
    let (left, right) = microsatellite_boundaries(site.start, site.end, read_start, cigar);
    let in_duplex =
        if rec.is_reverse() { in_duplex_region(false, duplex_base_index, right) } else { in_duplex_region(true, duplex_base_index, left) };
    if in_duplex { ConsensusType::Dual } else { ConsensusType::Single }
}
