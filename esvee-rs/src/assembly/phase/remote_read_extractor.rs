//! RemoteReadExtractor: reads the mates of an assembly's reads from their remote regions in the BAMs.

use std::collections::HashSet;

use crate::assembly::config::AssemblyConfig;
use crate::assembly::read::bam_reader::BamReader;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::remote_region::{RemoteRegion, RemoteRegionP};
use crate::assembly::types::support_read::SupP;
use crate::common::java_sort::java_sort;

const MAX_MATCHED_READ_CHECK: i32 = 100;
const MAX_REGIONS: usize = 100;
const MAX_TOTAL_READS: usize = 5000;
const MAX_REGION_READS: usize = 1000;

pub struct RemoteReadExtractor<'a> {
    bam_reader: Option<BamReader>,
    config: Option<&'a AssemblyConfig>,
    total_remote_reads_search: i32,
    total_remote_reads_matched: i32,
    total_remote_read_slices: i32,
    remote_reads_matched: i32,
    remote_read_slices: i32,
}

pub fn collect_candidate_remote_regions(assembly: AsmP, phased_assemblies: &[AsmP], check_assembly_matches: bool) -> Vec<RemoteRegionP> {
    let mut remote_regions: Vec<RemoteRegionP> = Vec::new();
    for remote_region in assembly.remote_regions().iter().copied() {
        if remote_region.is_supp_only_region() {
            continue;
        }
        if !check_assembly_matches {
            remote_regions.push(remote_region);
            continue;
        }
        let mut matches_assembly = false;
        for other_assembly in phased_assemblies.iter().copied() {
            if other_assembly == assembly {
                continue;
            }
            if other_assembly.discordant_only() {
                continue;
            }
            if !remote_region.overlaps_assembly(&other_assembly) {
                continue;
            }
            let mut reads_checked = 0;
            for read_id in remote_region.read_ids().iter() {
                if other_assembly.support().iter().filter(|x| x.type_().is_split_support()).any(|x| x.id() == read_id) {
                    matches_assembly = true;
                    break;
                }
                if other_assembly.candidate_support().iter().any(|x| x.id() == read_id) {
                    matches_assembly = true;
                    break;
                }
                reads_checked += 1;
                if reads_checked > MAX_MATCHED_READ_CHECK {
                    break;
                }
            }
            if matches_assembly {
                break;
            }
        }
        if !matches_assembly {
            remote_regions.push(remote_region);
        }
    }
    remote_regions
}

impl<'a> RemoteReadExtractor<'a> {
    pub fn new(bam_reader: Option<BamReader>, config: Option<&'a AssemblyConfig>) -> Self {
        RemoteReadExtractor { bam_reader, config, total_remote_reads_search: 0, total_remote_reads_matched: 0, total_remote_read_slices: 0, remote_reads_matched: 0, remote_read_slices: 0 }
    }
    pub fn remote_read_slices(&self) -> i32 {
        self.total_remote_read_slices
    }
    pub fn remote_reads_search(&self) -> i32 {
        self.total_remote_reads_search
    }
    pub fn remote_reads_matched(&self) -> i32 {
        self.total_remote_reads_matched
    }
    pub fn reset_counts(&mut self) {
        self.total_remote_reads_search = 0;
        self.total_remote_reads_matched = 0;
        self.total_remote_read_slices = 0;
    }

    /// The list of regions is merged and sorted in place, as in Java.
    pub fn extract_remote_region_reads(&mut self, phase_group_id: i32, remote_regions: &mut Vec<RemoteRegionP>, mut apply_thresholds: bool) -> Vec<ReadP> {
        RemoteRegion::merge_regions(remote_regions);
        java_sort(remote_regions, |a, b| (-a.read_count()).cmp(&(-b.read_count())));
        let mut min_remote_region_read_count = 1;
        let mut remote_region_slice_count = 0usize;
        let mut remote_region_read_count = 0usize;
        self.remote_reads_matched = 0;
        self.remote_read_slices = 0;
        let total_remote_reads: i32 = remote_regions.iter().map(|x| x.read_count()).sum();
        apply_thresholds |= remote_regions.len() > MAX_REGIONS;
        if apply_thresholds && total_remote_reads as usize >= MAX_REGION_READS {
            min_remote_region_read_count = (total_remote_reads as f64).log10().floor() as i32;
        }
        let mut remote_region_reads: Vec<ReadP> = Vec::new();
        for remote_region in remote_regions.iter() {
            if min_remote_region_read_count > 1 && remote_region.read_count() < min_remote_region_read_count {
                break;
            }
            let remote_reads = self.extract_remote_reads(phase_group_id, remote_region, MAX_REGION_READS);
            remote_region_read_count += remote_reads.len();
            remote_region_reads.extend(remote_reads);
            remote_region_slice_count += 1;
            if remote_region_slice_count >= MAX_REGIONS || remote_region_read_count >= MAX_TOTAL_READS {
                break;
            }
        }
        remote_region_reads
    }

    fn extract_remote_reads(&mut self, phase_group_id: i32, remote_region: &RemoteRegion, max_read_count: usize) -> Vec<ReadP> {
        let mut source_read_ids: Vec<&String> = Vec::new();
        if remote_region.read_count() as usize <= max_read_count {
            source_read_ids.extend(remote_region.read_ids().iter());
        } else {
            for read_id in remote_region.read_ids().iter() {
                source_read_ids.push(read_id);
                if source_read_ids.len() >= max_read_count {
                    break;
                }
            }
        }
        self.total_remote_reads_search += source_read_ids.len() as i32;
        let mut matched_remote_reads: Vec<ReadP> = Vec::new();
        if let (Some(bam_reader), Some(config)) = (self.bam_reader.as_mut(), self.config) {
            // read ids are trimmed, so a record matches when its name contains one of them;
            // ids of one length are looked up as suffix-free substrings through a set per length
            let exact: HashSet<&str> = source_read_ids.iter().map(|s| s.as_str()).collect();
            let mut lengths: Vec<usize> = exact.iter().map(|s| s.len()).collect();
            lengths.sort();
            lengths.dedup();
            bam_reader.slice_bam(config, remote_region.chromosome, remote_region.start(), remote_region.end(), |record, is_reference| {
                let name = record.read_name.as_str();
                let contained = lengths.iter().any(|&len| {
                    if len > name.len() {
                        return false;
                    }
                    if len == 0 {
                        return true;
                    }
                    (0..=name.len() - len).any(|i| name.is_char_boundary(i) && name.is_char_boundary(i + len) && exact.contains(&name[i..i + len]))
                });
                if !contained {
                    return;
                }
                let remote_read = Read::new(record);
                if is_reference {
                    remote_read.m().mark_reference();
                }
                matched_remote_reads.push(remote_read);
            });
            self.total_remote_reads_matched += matched_remote_reads.len() as i32;
            self.remote_reads_matched += matched_remote_reads.len() as i32;
        }
        self.total_remote_read_slices += 1;
        self.remote_read_slices += 1;
        if self.remote_read_slices % 100 == 0 && self.remote_reads_matched > 10000 {
            crate::debug!(
                "pgId({}) remote region read extraction: slices({}) matched({}) last region({}:{}-{})",
                phase_group_id, self.remote_read_slices, self.remote_reads_matched, remote_region.chromosome, remote_region.start(), remote_region.end()
            );
        }
        matched_remote_reads.into_iter().filter(|x| !x.is_supplementary()).collect()
    }
}

pub fn purge_supplementary_reads(assembly: &JunctionAssembly, remote_reads: &mut Vec<ReadP>) {
    let junction_supplementaries: Vec<SupP> = assembly.support().iter().filter(|x| x.supplementary_data().is_some()).cloned().collect();
    if junction_supplementaries.is_empty() {
        return;
    }
    if junction_supplementaries.len() > MAX_REGION_READS || remote_reads.len() > MAX_REGION_READS {
        remote_reads.retain(|read| read.supplementary_data().is_none());
        return;
    }
    remote_reads.retain(|read| !junction_supplementaries.iter().any(|x| x.matches_fragment_read(read, false) && x.first_in_pair() == read.first_in_pair()));
}
