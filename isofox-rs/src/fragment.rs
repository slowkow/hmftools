//! Fragment: a read, or a read and its mate (in the order given: the later read first, as the allocator pairs them).

use crate::common::*;
use crate::gene::RegionReadData;
use crate::read::{consensus_duplicate_count, AltAlignment, IntJavaMap, Read};
use hmf_common::regions::BaseRegion;

/// fragmentCount: one plus the first read's consensus duplicates
pub fn fragment_count(reads: &[&Read]) -> i32 {
    1 + consensus_duplicate_count(reads[0].consensus.as_deref())
}

pub fn min_num_loci(reads: &[&Read]) -> usize {
    reads.iter().map(|r| r.num_loci()).min().unwrap()
}

/// altLoci: those of the read with the fewest loci (the first of equals)
#[allow(dead_code)]
pub fn alt_loci<'a>(reads: &[&'a Read]) -> Option<&'a Vec<AltAlignment>> {
    let mut min = reads[0];
    for r in reads {
        if r.num_loci() < min.num_loci() {
            min = r;
        }
    }
    min.alt_loci.as_ref()
}

pub fn is_fully_intronic(reads: &[&Read]) -> bool {
    reads.iter().all(|r| r.mapped_regions.is_empty())
}

pub fn contains_split(reads: &[&Read]) -> bool {
    reads.iter().any(|r| r.has_split)
}

pub fn min_alignment_start(reads: &[&Read]) -> i32 {
    reads.iter().map(|r| r.pos_start).min().unwrap()
}

pub fn max_alignment_end(reads: &[&Read]) -> i32 {
    reads.iter().map(|r| r.pos_end).max().unwrap()
}

pub fn merged_mappings(reads: &[&Read]) -> Vec<BaseRegion> {
    if reads.len() == 1 {
        return reads[0].mapped_coords.alignments.clone();
    }
    derive_common_regions(&reads[0].mapped_coords.alignments, &reads[1].mapped_coords.alignments)
}

pub fn unique_valid_regions(reads: &[&Read]) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::new();
    for r in reads {
        for (reg, mt) in r.mapped_regions.entries() {
            if mt.valid_exon_match() && !out.contains(&reg) {
                out.push(reg);
            }
        }
    }
    out
}

pub fn spans_multiple_regions(reads: &[&Read], regions: &[RegionReadData], t: i32) -> bool {
    let mut found: Vec<usize> = Vec::new();
    for r in reads {
        for (reg, mt) in r.mapped_regions.entries() {
            if regions[reg].has_trans_id(t) && mt != RegionMatchType::ExonIntron && !found.contains(&reg) {
                found.push(reg);
                if found.len() > 1 {
                    return true;
                }
            }
        }
    }
    false
}

/// validTypeTranscripts: the first read's transcripts that every read supports
pub fn valid_type_transcripts(reads: &[&Read]) -> Vec<i32> {
    reads[0].trans_class.keys().into_iter().filter(|t| reads.iter().all(|r| valid_transcript_type(r.trans_classification(*t)))).collect()
}

/// invalidTranscripts: a HashSet<Integer>
pub fn invalid_transcripts(reads: &[&Read], valid: &[i32]) -> IntJavaMap<()> {
    let mut out: IntJavaMap<()> = IntJavaMap::new();
    for r in reads {
        for t in r.trans_class.keys() {
            if !valid.contains(&t) {
                out.insert(t, ());
            }
        }
    }
    out
}

/// setOtherTranscripts
pub fn set_other_transcripts(reads: &mut [&mut Read], valid: &[i32]) {
    for r in reads.iter_mut() {
        for (t, m) in r.trans_class.values_mut_unordered() {
            if valid_transcript_type(*m) && !valid.contains(&t) {
                *m = TransMatchType::OtherTrans;
            }
        }
    }
}

pub fn has_transcript_classification(reads: &[&Read], t: i32, m: TransMatchType) -> bool {
    reads.iter().any(|r| r.trans_classification(t) == m)
}

pub fn reads_in_different_exons(reads: &[&Read], regions: &[RegionReadData], t: i32) -> bool {
    if reads.len() < 2 {
        return false;
    }
    let k2 = reads[1].mapped_regions.keys();
    for r1 in reads[0].mapped_regions.keys() {
        let rank = regions[r1].exon_rank(t);
        if rank == NO_EXON {
            continue;
        }
        if k2.iter().any(|r2| regions[*r2].exon_rank(t) == rank) {
            return false;
        }
    }
    true
}

/// orientation: the first-of-pair read's, None if both reads face the same way
pub fn orientation(reads: &[&Read]) -> Option<i8> {
    let r1 = reads[0];
    if reads.len() == 1 {
        return Some(r1.orient());
    }
    let r2 = reads[1];
    if r1.orient() == r2.orient() {
        return None;
    }
    Some(if r1.is_first_of_pair() { r1.orient() } else { r2.orient() })
}
