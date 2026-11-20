//! The first pass over the unmapping regions (Java RegionUnmapper): reads that qualify are unmapped and moved to
//! their mate's position, or to the unmapped reads at the end of the BAM when the mate is unmapped too.

use crate::bamin::BamReader;
use crate::config::Config;
use crate::out::RawBamOut;
use crate::partition::{fully_unmapped, should_filter_read};
use crate::rec::*;
use crate::regions::{compare_chromosomes, ChrBaseRegion};
use crate::unmap::*;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const REGION_PROXIMITY_BUFFER: i32 = 50;

/// UnmapRegion: unmapping regions that lie close together, sliced as one
pub struct UnmapRegion {
    pub chromosome: String,
    pub tid: i32,
    pub start: i32,
    pub end: i32,
    pub regions: Vec<UnmappingRegion>,
}

/// RegionUnmapper.createThreadTasks: the slices, ordered by chromosome and start
pub fn create_unmap_regions(config: &Config, unmapper: &ReadUnmapper) -> Vec<UnmapRegion> {
    let mut out: Vec<UnmapRegion> = Vec::new();
    for tid in 0..unmapper.chromosome_count() {
        let chromosome = unmapper.chromosome(tid);
        let regions = unmapper.regions(tid);
        if regions.is_empty() || !config.specific.include_chromosome(chromosome) {
            continue;
        }
        let mut current: Option<usize> = None;
        for region in regions {
            if !config.specific.regions.is_empty() && !config.specific.regions.iter().any(|r| r.start <= region.end && r.end >= region.start) {
                continue;
            }
            let region_start = (region.start - UNMAP_MAX_NON_OVERLAPPING_BASES * 2).max(1);
            match current {
                Some(i) if out[i].end >= region_start - REGION_PROXIMITY_BUFFER => {
                    out[i].end = region.end;
                    out[i].regions.push(*region);
                }
                _ => {
                    out.push(UnmapRegion { chromosome: chromosome.to_string(), tid: tid as i32, start: region_start, end: region.end, regions: vec![*region] });
                    current = Some(out.len() - 1);
                }
            }
        }
    }
    out.sort_by(|a, b| if a.chromosome == b.chromosome { a.start.cmp(&b.start) } else { compare_chromosomes(&a.chromosome, &b.chromosome) });
    out
}

/// A range of one worker's file of fully unmapped reads, for one region.
pub struct FullUnmappedChunk {
    pub region_index: usize,
    pub path: String,
    pub start: u64,
    pub end: u64,
    pub reads: u64,
}

pub struct Stage1Result {
    pub stats: UnmapStats,
    pub processed_reads: u64,
    /// in region order
    pub full_unmapped_chunks: Vec<FullUnmappedChunk>,
    pub temp_files: Vec<String>,
}

/// The tag that carries a relocated read's order of discovery to the second pass; removed there.
const ORDER_TAG: &[u8; 3] = b"zoZ";
const ORDER_TAG_LEN: usize = 3 + 16 + 1;

fn append_order_tag(rec: &mut bam::Record, region_index: usize, seq: u64) {
    let mut aux = aux_bytes(rec).to_vec();
    aux.extend_from_slice(ORDER_TAG);
    aux.extend_from_slice(format!("{:06x}{:010x}", region_index, seq).as_bytes());
    aux.push(0);
    set_aux_bytes(rec, &aux);
}

/// Removes the order tag and returns its value.
pub fn take_order_tag(rec: &mut bam::Record) -> u64 {
    let aux = aux_bytes(rec);
    let n = aux.len();
    assert!(n >= ORDER_TAG_LEN && &aux[n - ORDER_TAG_LEN..n - ORDER_TAG_LEN + 3] == ORDER_TAG);
    let value = u64::from_str_radix(std::str::from_utf8(&aux[n - 17..n - 1]).unwrap(), 16).unwrap();
    let trimmed = aux[..n - ORDER_TAG_LEN].to_vec();
    set_aux_bytes(rec, &trimmed);
    value
}

/// Finds the partition whose regions hold a position.
pub struct PartitionLookup {
    /// (tid, start, end, partition index), sorted
    spans: Vec<(i32, i32, i32, usize)>,
}

impl PartitionLookup {
    pub fn new(partitions: &[Vec<ChrBaseRegion>], names: &[String]) -> Self {
        let mut spans = Vec::new();
        for (i, p) in partitions.iter().enumerate() {
            for r in p {
                if let Some(tid) = names.iter().position(|n| *n == r.chromosome) {
                    spans.push((tid as i32, r.start, r.end, i));
                }
            }
        }
        spans.sort();
        PartitionLookup { spans }
    }
    pub fn find(&self, tid: i32, pos: i32) -> Option<usize> {
        let i = self.spans.partition_point(|s| (s.0, s.1) <= (tid, pos));
        if i == 0 {
            return None;
        }
        let s = &self.spans[i - 1];
        if s.0 == tid && pos >= s.1 && pos <= s.2 { Some(s.3) } else { None }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    config: &Config, unmapper: &ReadUnmapper, regions: &[UnmapRegion], expected_valid_tid: &[bool], header: &bam::HeaderView,
    lookup: &PartitionLookup, buckets: &[Mutex<RawBamOut>], temp_prefix: &str, extract_fully_unmapped: bool,
) -> Result<Stage1Result> {
    // one extra task extracts the reads that are already unmapped and unplaced
    let task_count = regions.len() + usize::from(extract_fully_unmapped);
    let next = AtomicUsize::new(0);
    let worker_count = config.threads.min(task_count).max(1);
    let results: Vec<Result<Stage1Result>> = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for w in 0..worker_count {
            let next = &next;
            handles.push(scope.spawn(move || -> Result<Stage1Result> {
                let mut reader = BamReader::open(&config.bam_files, &config.ref_genome_file, 1)?;
                let path = format!("{temp_prefix}.full_unmapped.{w}.tmp");
                let mut full_out = RawBamOut::create(&path, config.compress_level)?;
                let mut result = Stage1Result { stats: UnmapStats::default(), processed_reads: 0, full_unmapped_chunks: Vec::new(), temp_files: vec![path.clone()] };
                let _ = header;
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= task_count {
                        break;
                    }
                    let chunk_start = full_out.compressed_offset();
                    let written_before = full_out.written;
                    if i == regions.len() {
                        // RegionUnmapper.processFullyUnmappedReads
                        for bam_file in &config.bam_files {
                            let mut r = bam::IndexedReader::from_path(bam_file)?;
                            r.set_reference(&config.ref_genome_file)?;
                            r.fetch(bam::FetchDefinition::Unmapped)?;
                            let mut rec = bam::Record::new();
                            while let Some(x) = r.read(&mut rec) {
                                x?;
                                full_out.write(&rec)?;
                                result.stats.existing_unmapped += 1;
                            }
                        }
                    } else {
                        let region = &regions[i];
                        let mut state = UnmapRegionState {
                            chromosome: region.chromosome.clone(),
                            slice_start: region.start,
                            slice_end: region.end,
                            regions: region.regions.clone(),
                            last_matched_index: Some(0),
                        };
                        let mut seq: u64 = 0;
                        let mut err: Result<()> = Ok(());
                        reader.slice_region(&region.chromosome, region.start, region.end, |rec| {
                            if err.is_err() || should_filter_read(rec, expected_valid_tid) {
                                return;
                            }
                            if !has(rec, FLAG_UNMAPPED) {
                                let (s, e) = (start(rec), end(rec));
                                if !region.regions.iter().any(|x| overlaps_unmap_region(x, s, e)) {
                                    return;
                                }
                            }
                            result.processed_reads += 1;
                            set_flag(rec, FLAG_DUPLICATE, false);
                            let already_unmapped = has(rec, FLAG_UNMAPPED);
                            unmapper.check_transform_read(rec, &mut state, &mut result.stats);
                            let internally_unmapped = !already_unmapped && has(rec, FLAG_UNMAPPED);
                            let fully = fully_unmapped(rec);
                            if internally_unmapped || fully {
                                if has(rec, FLAG_SUPPLEMENTARY | FLAG_SECONDARY) {
                                    return;
                                }
                                if fully {
                                    err = full_out.write(rec);
                                } else if let Some(p) = lookup.find(rec.tid(), start(rec)) {
                                    append_order_tag(rec, i, seq);
                                    seq += 1;
                                    err = buckets[p].lock().unwrap().write(rec);
                                }
                            }
                        });
                        err?;
                    }
                    if full_out.written > written_before {
                        full_out.flush()?;
                        result.full_unmapped_chunks.push(FullUnmappedChunk { region_index: i, path: path.clone(), start: chunk_start, end: full_out.compressed_offset(), reads: full_out.written - written_before });
                    }
                }
                full_out.close()?;
                Ok(result)
            }));
        }
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut total = Stage1Result { stats: UnmapStats::default(), processed_reads: 0, full_unmapped_chunks: Vec::new(), temp_files: Vec::new() };
    for r in results {
        let r = r?;
        total.stats.merge(&r.stats);
        total.processed_reads += r.processed_reads;
        total.full_unmapped_chunks.extend(r.full_unmapped_chunks);
        total.temp_files.extend(r.temp_files);
    }
    total.full_unmapped_chunks.sort_by_key(|c| c.region_index);
    Ok(total)
}

/// Sorts one partition's relocated reads into a new BAM, in the order the Java run gives them by sorting its
/// interim BAM: by reference, position and strand, then by order of discovery.
pub fn sort_bucket(path: &str, sorted_path: &str, header: &bam::HeaderView) -> Result<()> {
    let mut reader = bam::Reader::from_path(path)?;
    let mut keyed: Vec<((u32, i64, bool, u64), bam::Record)> = Vec::new();
    let mut rec = bam::Record::new();
    while let Some(x) = reader.read(&mut rec) {
        x?;
        let mut r = std::mem::replace(&mut rec, bam::Record::new());
        let order = take_order_tag(&mut r);
        keyed.push(((r.tid() as u32, r.pos(), has(&r, FLAG_REVERSE), order), r));
    }
    keyed.sort_by_key(|k| k.0);
    let mut out = RawBamOut::create(sorted_path, 1)?;
    out.write_header(header)?;
    for (_, r) in &keyed {
        out.write(r)?;
    }
    out.close()
}
