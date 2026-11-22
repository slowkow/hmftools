//! RemoteRegionFinder: collects the remote regions of an assembly's mates, discordant reads and supplementaries.

use crate::assembly::common_utils::{below_min_qual, positions_overlap};
use crate::assembly::constants::*;
use crate::assembly::globals::{has_paired_reads, is_illumina};
use crate::assembly::ptr::P;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::sam::*;
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::remote_region::{RemoteRegion, RemoteRegionP};
use crate::assembly::types::RemoteReadType;
use crate::assembly::common_utils::human_chromosome_contains as is_human_chromosome;

const MAX_REMOTE_REGION_LENGTH: i32 = 10000;

/// prep ReadFilters.filterLowQualRead(read): also replaces non-standard bases by N as it goes
pub fn filter_low_qual_read(read: &mut SamRecord) -> bool {
    filter_low_qual_read_perc(read, 0.5)
}

pub fn filter_low_qual_read_perc(read: &mut SamRecord, max_permitted_low_percent: f64) -> bool {
    if !is_illumina() {
        return false;
    }
    let base_length = read.read_bases.len();
    let qual_count_threshold = (base_length as f64 * max_permitted_low_percent).floor() as i32 + 1;
    let mut low_qual_count = 0;
    for i in 0..base_length {
        if !matches!(read.read_bases[i], b'A' | b'C' | b'G' | b'T' | b'N') {
            read.read_bases[i] = b'N';
        }
        if below_min_qual(read.base_qualities[i]) {
            low_qual_count += 1;
            if low_qual_count >= qual_count_threshold {
                return true;
            }
        } else {
            let high_qual_count = i as i32 + 1 - low_qual_count;
            if high_qual_count >= qual_count_threshold {
                return false;
            }
        }
    }
    false
}

pub fn find_remote_regions(assembly: AsmP, discordant_reads: &[ReadP], remote_junction_mates: &[ReadP], supp_junction_reads: &[ReadP]) {
    if remote_junction_mates.is_empty() && discordant_reads.is_empty() && supp_junction_reads.is_empty() {
        return;
    }
    let mut remote_regions: Vec<RemoteRegionP> = Vec::new();
    for x in discordant_reads {
        add_or_create_mate_remote_region(&mut remote_regions, x, false);
    }
    for read in remote_junction_mates {
        add_or_create_mate_remote_region(&mut remote_regions, read, true);
    }
    for read in supp_junction_reads {
        let sc_length = if assembly.is_forward_junction() { read.right_clip_length() } else { read.left_clip_length() };
        add_or_create_supplementary_remote_region(&mut remote_regions, read, sc_length);
    }
    RemoteRegion::merge_regions(&mut remote_regions);
    if has_paired_reads() {
        RemoteRegion::purge_weak_supplementary_regions(&mut remote_regions);
        RemoteRegion::purge_low_qual_discordant_only_regions(&mut remote_regions);
    }
    let mut index = 0;
    while index < remote_regions.len() {
        let remote_region = remote_regions[index];
        if assembly_overlaps_remote_region(&assembly, &remote_region) && !remote_region.has_discordant_reads() {
            remote_regions.remove(index);
        } else if remote_region.end() - remote_region.start() + 1 > MAX_REMOTE_REGION_LENGTH {
            crate::warn!("assembly({}) excluding long remote region({})", *assembly, *remote_region);
            remote_regions.remove(index);
        } else {
            index += 1;
        }
    }
    assembly.m().add_remote_regions(&remote_regions);
}

fn assembly_overlaps_remote_region(assembly: &JunctionAssembly, remote_region: &RemoteRegion) -> bool {
    if assembly.junction().chromosome != remote_region.chromosome {
        return false;
    }
    if assembly.is_forward_junction() {
        positions_overlap(assembly.ref_base_position(), assembly.junction().position, remote_region.start(), remote_region.end())
    } else {
        positions_overlap(assembly.junction().position, assembly.ref_base_position(), remote_region.start(), remote_region.end())
    }
}

pub fn add_or_create_mate_remote_region(remote_regions: &mut Vec<RemoteRegionP>, read: &Read, is_junction_read: bool) {
    if read.is_mate_unmapped() || read.is_supplementary() {
        return;
    }
    let Some(mate_chr) = read.mate_chromosome() else { return };
    if !is_human_chromosome(mate_chr) {
        return;
    }
    add_or_create_remote_region(
        remote_regions,
        read,
        if is_junction_read { RemoteReadType::JunctionMate } else { RemoteReadType::Discordant },
        mate_chr,
        read.mate_alignment_start(),
        read.mate_alignment_end(),
    );
}

fn add_or_create_remote_region(remote_regions: &mut Vec<RemoteRegionP>, read: &Read, read_type: RemoteReadType, remote_chr: &str, remote_pos_start: i32, remote_pos_end: i32) -> RemoteRegionP {
    let remote_region = match remote_regions.iter().find(|x| x.overlaps(remote_chr, remote_pos_start, remote_pos_end)) {
        None => {
            let r = P::new(RemoteRegion::new(intern(remote_chr), remote_pos_start, remote_pos_end, read.id(), read_type));
            remote_regions.push(r);
            r
        }
        Some(r) => {
            r.m().add_read_details(read.id(), remote_pos_start, remote_pos_end, read_type);
            *r
        }
    };
    if read_type == RemoteReadType::Discordant && !remote_region.has_high_qual_discordant_read() && is_high_quality_discordant_read(read) {
        remote_region.m().set_has_high_qual_discordant_read();
    }
    remote_region
}

fn is_high_quality_discordant_read(read: &Read) -> bool {
    let alignment_score = read.bam_record().get_integer_attribute(ALIGNMENT_SCORE_ATTRIBUTE).unwrap_or(0);
    if alignment_score < REMOTE_REGION_DISC_READ_BASE_MIN_AS {
        return false;
    }
    // the Java filter may replace non-standard bases in the record; reads reaching here had that done when loaded
    let mut record = read.bam_record().clone();
    !filter_low_qual_read_perc(&mut record, REMOTE_REGION_DISC_READ_BASE_MIN_QUAL_PERC)
}

fn add_or_create_supplementary_remote_region(remote_regions: &mut Vec<RemoteRegionP>, read: &Read, read_sc_length: i32) {
    let supp_data_list = SupplementaryReadData::extract_alignments(read.bam_record().get_string_attribute(SUPPLEMENTARY_ATTRIBUTE).as_deref())
        .expect("supplementary alignments for a read with supplementary data");
    for supp_data in supp_data_list {
        if !is_human_chromosome(supp_data.chromosome) {
            continue;
        }
        let remote_pos_end = get_alignment_end_from_cigar(supp_data.position, &supp_data.cigar);
        if read.chromosome() == supp_data.chromosome && positions_overlap(read.alignment_start(), read.alignment_end(), supp_data.position, remote_pos_end) {
            continue;
        }
        let region = add_or_create_remote_region(remote_regions, read, RemoteReadType::Supplementary, supp_data.chromosome, supp_data.position, remote_pos_end);
        region.m().add_soft_clip_map_qual(read_sc_length, supp_data.map_quality);
    }
}
