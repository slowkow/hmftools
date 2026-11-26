//! The split of the genome into partitions (PartitionThread.splitRegionsIntoPartitions). Chromosomes,
//! regions and the -specific_regions / -specific_chr options are in `hmf_common::regions`.

pub use hmf_common::regions::*;

use crate::refgenome::RefGenome;

/// `ref_genome` None limits the regions to the human chromosomes (-std_chr_only).
pub fn split_regions_into_partitions(
    specific: &SpecificRegions, thread_count: usize, v38: bool, genome: &RefGenome, std_chr_only: bool,
) -> Vec<Vec<ChrBaseRegion>> {
    let mut partitions: Vec<Vec<ChrBaseRegion>> = Vec::new();
    if !specific.regions.is_empty() {
        let n = specific.regions.len();
        if thread_count == n && n > 1 {
            for r in &specific.regions {
                partitions.push(vec![r.clone()]);
            }
        } else if n == 1 && thread_count > 1 {
            let r = &specific.regions[0];
            let interval = (r.base_length() as f64 / thread_count as f64).ceil() as i32;
            let mut start = r.start;
            for _ in 0..thread_count {
                let end = (start + interval - 1).min(r.end);
                partitions.push(vec![ChrBaseRegion::new(&r.chromosome, start, end)]);
                start = end + 1;
            }
        } else {
            partitions.push(specific.regions.clone());
        }
        return partitions;
    }

    let mut genome_regions: Vec<ChrBaseRegion> = Vec::new();
    if std_chr_only {
        // the Java code takes these lengths from its built-in coordinates; they equal the reference's
        for c in HUMAN_CHROMOSOMES {
            let name = versioned_chromosome(v38, c);
            if !specific.include_chromosome(&name) {
                continue;
            }
            genome_regions.push(ChrBaseRegion::new(&name, 1, genome.chromosome_length(&name)));
        }
    } else {
        for name in genome.sequence_names() {
            if !specific.include_chromosome(name) {
                continue;
            }
            genome_regions.push(ChrBaseRegion::new(name, 1, genome.chromosome_length(name)));
        }
    }

    let total: i64 = genome_regions.iter().map(|r| if is_alt_region_contig(&r.chromosome) { 0 } else { r.base_length() as i64 }).sum();
    // Java: (int)ceil(totalLength / (double)threadCount) held in a long
    let interval: i64 = ((total as f64 / thread_count as f64).ceil() as i32) as i64;
    let chr_end_buffer = (interval as f64 * 0.05).round() as i32;
    let mut current_length: i64 = 0;
    partitions.push(Vec::new());
    for r in &genome_regions {
        if is_alt_region_contig(&r.chromosome) {
            partitions.last_mut().unwrap().push(r.clone());
            continue;
        }
        let mut next_start = 1;
        let region_length = r.base_length();
        let mut remaining = region_length;
        while current_length + remaining as i64 >= interval {
            let remaining_interval = (interval - current_length) as i32;
            let mut region_end = next_start + remaining_interval - 1;
            if region_length - region_end < chr_end_buffer {
                region_end = region_length;
                remaining = 0;
            }
            partitions.last_mut().unwrap().push(ChrBaseRegion::new(&r.chromosome, next_start, region_end));
            partitions.push(Vec::new());
            current_length = 0;
            if remaining == 0 {
                break;
            }
            next_start = region_end + 1;
            remaining = region_length - next_start + 1;
        }
        if remaining <= 0 {
            continue;
        }
        current_length += remaining as i64;
        partitions.last_mut().unwrap().push(ChrBaseRegion::new(&r.chromosome, next_start, region_length));
    }
    partitions.into_iter().filter(|p| !p.is_empty()).collect()
}
