//! ReadTranscriptUtils: a read's match against exon regions and transcripts.

use crate::common::*;
use crate::gene::RegionReadData;
use crate::read::{IntJavaMap, MappedCoords, RKey, Read, INVALID_INDEX};
use hmf_common::ensembl::TranscriptData;
use hmf_common::regions::BaseRegion;

const MIN_SC_BASE_MATCH: i32 = 2;
const MAX_SC_BASE_MATCH: i32 = 10;
const MAX_SC_WITHIN_EXON_LENGTH: i32 = 2;

fn rkey(regions: &[RegionReadData], idx: usize) -> RKey {
    RKey { idx: idx as u32, start: regions[idx].start }
}

/// Read.findOverlappingRegions: the regions any of the read's sections overlap, in list order
pub fn find_overlapping_regions(regions: &[RegionReadData], candidates: Option<&[usize]>, read: &Read) -> Vec<usize> {
    match candidates {
        Some(c) => c.iter().copied().filter(|i| read.mapped_coords.alignments_overlap(regions[*i].start, regions[*i].end)).collect(),
        None => (0..regions.len()).filter(|i| read.mapped_coords.alignments_overlap(regions[*i].start, regions[*i].end)).collect(),
    }
}

/// getRegionMatchType for one read section
pub fn region_match_type(mc: &MappedCoords, region: &RegionReadData, mapping_index: i32) -> RegionMatchType {
    if mapping_index == INVALID_INDEX || mapping_index as usize >= mc.alignment_count() {
        return RegionMatchType::None;
    }
    let s = mc.alignments[mapping_index as usize];
    if s.end < region.start || s.start > region.end {
        return RegionMatchType::None;
    }
    if s.start < region.start || s.end > region.end {
        return RegionMatchType::ExonIntron;
    }
    if s.start > region.start && s.end < region.end {
        return RegionMatchType::WithinExon;
    }
    RegionMatchType::ExonBoundary
}

/// processOverlappingRegions
pub fn process_overlapping_regions(read: &mut Read, regions: &[RegionReadData], overlapping: &[usize]) {
    let mut transcripts: IntJavaMap<()> = IntJavaMap::new();
    let has_soft_clipping = read.is_left_clipped() || read.is_right_clipped();

    for &ri in overlapping {
        let region = &regions[ri];
        for r in &region.trans_exon_refs {
            transcripts.insert(r.trans_id, ());
        }
        // setRegionMatchType, then the put of its result
        let mapping_index = read.mapped_coords.find_region_index(region.start, region.end);
        let match_type = if mapping_index == INVALID_INDEX {
            RegionMatchType::None
        } else {
            let m = region_match_type(&read.mapped_coords, region, mapping_index);
            read.mapped_regions.put(rkey(regions, ri), m);
            m
        };
        read.mapped_regions.put(rkey(regions, ri), match_type);

        let check_missed = match_type == RegionMatchType::ExonIntron || (has_soft_clipping && match_type.exon_boundary());
        if check_missed {
            check_missed_junctions(read, regions, ri);
        }
    }

    let mapped_region_count = read.mapped_region_count();

    for trans_id in transcripts.keys() {
        let mut trans_match = TransMatchType::Unknown;
        let mut trans_regions: Vec<usize> = overlapping.iter().copied().filter(|i| regions[*i].has_trans_id(trans_id)).collect();

        if trans_regions.len() == 1 && mapped_region_count == 1 {
            let mt = read.mapped_regions.get(trans_regions[0]).unwrap_or(RegionMatchType::None);
            if mt == RegionMatchType::None {
                trans_match = TransMatchType::Alt;
            } else if mt == RegionMatchType::ExonIntron {
                trans_match = TransMatchType::Unspliced;
            }
        } else if mapped_region_count > trans_regions.len() {
            trans_match = TransMatchType::Alt;
        } else {
            let mut min_rank = 0;
            let mut max_rank = 0;
            // Collections.sort by compareTo (start), stable
            trans_regions.sort_by_key(|i| regions[*i].start);
            let n = trans_regions.len();
            for (region_index, &ri) in trans_regions.iter().enumerate() {
                let region = &regions[ri];
                let rank = region.exon_rank(trans_id);
                max_rank = max_rank.max(rank);
                min_rank = if min_rank == 0 { rank } else { rank.min(min_rank) };

                let mapping_index = read.mapped_coords.find_region_index(region.start, region.end);
                let mut adjusted = mapping_index;
                if read.mapped_coords.inferred(SE_START) {
                    adjusted -= 1;
                }
                if adjusted < 0 || adjusted as usize != region_index {
                    trans_match = TransMatchType::Alt;
                    break;
                }
                let mt = read.mapped_regions.get(ri).unwrap_or(RegionMatchType::None);
                if mt == RegionMatchType::ExonIntron {
                    trans_match = TransMatchType::Unspliced;
                    break;
                } else {
                    let s = read.mapped_coords.alignments[mapping_index as usize];
                    let miss_start = s.start > region.start;
                    let miss_end = s.end < region.end;
                    if region_index == 0 {
                        if miss_end {
                            trans_match = TransMatchType::Alt;
                            break;
                        }
                    } else if region_index == n - 1 {
                        if miss_start {
                            trans_match = TransMatchType::Alt;
                            break;
                        }
                    } else if miss_start || miss_end {
                        trans_match = TransMatchType::Alt;
                        break;
                    }
                }
            }
            if trans_match == TransMatchType::Unknown {
                let expected = max_rank - min_rank + 1;
                if (trans_regions.len() as i32) < expected {
                    trans_match = TransMatchType::Alt;
                }
            }
        }

        if trans_match == TransMatchType::Unknown {
            trans_match = if trans_regions.len() > 1 { TransMatchType::SpliceJunction } else { TransMatchType::Exonic };
        }

        if valid_transcript_type(trans_match) && read.contains_soft_clipping() && !read.likely_adapter_soft_clipping() && crate::version::v201() {
            // ISOFOX 2.0.1: a soft clip not matched to the next exon makes the transcript ALT, left side first
            if read.is_left_clipped() && !read.mapped_coords.is_soft_clip_region_matched(SE_START) {
                trans_match = TransMatchType::Alt;
            } else if read.is_right_clipped() && !read.mapped_coords.is_soft_clip_region_matched(SE_END) {
                trans_match = TransMatchType::Alt;
            }
        } else if valid_transcript_type(trans_match) && read.contains_soft_clipping() && !read.likely_adapter_soft_clipping() {
            let side = if read.is_left_clipped() { SE_START } else { SE_END };
            if !read.mapped_coords.is_soft_clip_region_matched(side) && !short_clip_within_exon(read, side, regions, &trans_regions) {
                trans_match = TransMatchType::Alt;
            }
        }

        read.trans_class.insert(trans_id, trans_match);
    }
}

fn matches_other_region_bases(extra: &[u8], other: &RegionReadData, match_to_start: bool) -> bool {
    let len = other.length();
    if extra.len() > len {
        return false;
    }
    let other_bases = if match_to_start { &other.ref_bases[..extra.len()] } else { &other.ref_bases[len - extra.len()..] };
    other_bases == extra
}

fn check_missed_junctions(read: &mut Read, regions: &[RegionReadData], ri: usize) {
    if read.has_supp_alignment() {
        return;
    }
    let region = &regions[ri];

    // the start of the read
    let li = read.mapped_coords.lowest_index(false);
    let section = read.mapped_coords.alignments[li];
    let (read_start, read_end) = (section.start, section.end);
    let mut extra_len = 0;
    let mut sc_len = 0;
    let has_overhang = region.start > read_start && read_end > region.start && region.start - read_start <= MAX_SC_BASE_MATCH;
    if has_overhang {
        extra_len = region.start - read_start;
    }
    if read.is_left_clipped() && read_start <= region.start {
        sc_len = read.left_clip_length();
        extra_len += sc_len;
    }
    if (1..=MAX_SC_BASE_MATCH).contains(&extra_len) && sc_len <= MAX_SC_BASE_MATCH {
        let extra: Vec<u8> = read.read_bases()[..extra_len as usize].to_vec();
        let matched: Vec<usize> = region.pre_regions.iter().copied().filter(|p| matches_other_region_bases(&extra, &regions[*p], false)).collect();
        if !matched.is_empty() {
            read.mapped_coords.add_soft_clip_region_matched(true, matched.len() as i32);
            read.mapped_regions.put(rkey(regions, ri), RegionMatchType::ExonBoundary);
            if (matched.len() == 1 || extra_len < MIN_SC_BASE_MATCH) && !read.mapped_coords.inferred(SE_START) && has_overhang {
                // the section is the same object Java changes: its index is unchanged as no inferred
                // region has been added below it
                read.mapped_coords.alignments[li].start += region.start - read_start;
            }
            if matched.len() == 1 || extra_len >= MIN_SC_BASE_MATCH {
                for p in matched {
                    read.mapped_regions.put(rkey(regions, p), RegionMatchType::ExonBoundary);
                    read.mapped_coords.add_inferred_region(true, regions[p].end - extra_len + 1, regions[p].end);
                }
            }
        }
    }

    // the end of the read
    let hi = read.mapped_coords.highest_index(false);
    let section = read.mapped_coords.alignments[hi];
    let (read_start, read_end) = (section.start, section.end);
    let mut extra_len = 0;
    let mut sc_len = 0;
    let has_overhang = read_end > region.end && read_start < region.end && read_end - region.end <= MAX_SC_BASE_MATCH;
    if has_overhang {
        extra_len = read_end - region.end;
    }
    if read.is_right_clipped() && read_end >= region.end {
        sc_len = read.right_clip_length();
        extra_len += sc_len;
    }
    if (1..=MAX_SC_BASE_MATCH).contains(&extra_len) && sc_len <= MAX_SC_BASE_MATCH {
        let bases = read.read_bases();
        let n = bases.len();
        let extra: Vec<u8> = bases[n - extra_len as usize..n].to_vec();
        let matched: Vec<usize> = region.post_regions.iter().copied().filter(|p| matches_other_region_bases(&extra, &regions[*p], true)).collect();
        if !matched.is_empty() {
            read.mapped_coords.add_soft_clip_region_matched(false, matched.len() as i32);
            read.mapped_regions.put(rkey(regions, ri), RegionMatchType::ExonBoundary);
            if (matched.len() == 1 || extra_len < MIN_SC_BASE_MATCH) && !read.mapped_coords.inferred(SE_END) && has_overhang {
                read.mapped_coords.alignments[hi].end -= read_end - region.end;
            }
            if matched.len() == 1 || extra_len >= MIN_SC_BASE_MATCH {
                for p in matched {
                    read.mapped_regions.put(rkey(regions, p), RegionMatchType::ExonBoundary);
                    read.mapped_coords.add_inferred_region(false, regions[p].start, regions[p].start + extra_len - 1);
                }
            }
        }
    }
}

fn short_clip_within_exon(read: &Read, se: usize, regions: &[RegionReadData], trans_regions: &[usize]) -> bool {
    let clip = if se == SE_START { read.left_clip_length() } else { read.right_clip_length() };
    if clip > MAX_SC_WITHIN_EXON_LENGTH {
        return false;
    }
    let boundary = read.coords_boundary(se);
    let trans_boundary = if se == SE_START {
        trans_regions.iter().map(|i| regions[*i].start).min().unwrap_or(0)
    } else {
        trans_regions.iter().map(|i| regions[*i].end).max().unwrap_or(0)
    };
    if se == SE_START { boundary - clip >= trans_boundary } else { boundary + clip <= trans_boundary }
}

/// markRegionBases: adds depth to the region's bases from the read sections
pub fn mark_region_bases(coords: &[BaseRegion], region: &mut RegionReadData) {
    let (rs, re) = (region.start, region.end);
    let Some(depth) = region.ref_bases_matched.as_mut() else { return };
    let len = depth.len() as i32;
    for s in coords {
        if s.start > re || s.end < rs {
            continue;
        }
        let idx = if s.start > rs { s.start - rs } else { 0 };
        let overlap = s.end.min(re) - s.start.max(rs) + 1;
        if idx + overlap > len {
            hmf_common::error!("region({}-{}) read coords({} -> {}) regionBaseIndex({}) overlap({}) regionLength({})", rs, re, s.start, s.end, idx, overlap, len);
            return;
        }
        for d in &mut depth[idx as usize..(idx + overlap) as usize] {
            *d += 1;
        }
    }
}

/// calcFragmentLength: the fragment's length within a transcript, assuming it is spliced
pub fn calc_fragment_length(trans: &TranscriptData, min_pos: i32, max_pos: i32) -> i32 {
    let mut bases = 0;
    let mut start_found = false;
    for exon in &trans.exons {
        if !start_found {
            if min_pos < exon.start - MAX_SC_BASE_MATCH {
                break;
            }
            if min_pos > exon.end {
                continue;
            }
            if max_pos <= exon.end {
                return max_pos - min_pos + 1;
            }
            start_found = true;
            bases = exon.end - exon.start.max(min_pos) + 1;
        } else if max_pos > exon.end {
            bases += exon.end - exon.start + 1;
        } else if max_pos < exon.start {
            break;
        } else {
            bases += max_pos - exon.start + 1;
            break;
        }
    }
    bases
}
