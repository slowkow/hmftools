//! FragmentSizeDistribution: fragment lengths sampled from the first BAM, and the bounds taken from them.

use super::config::*;
use super::filters::mate_cigar_is_aligned;
use super::rec::*;
use super::types::inferred_insert_size_abs;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::collections::BTreeMap;

pub const FRAG_LENGTH_DIST_SAMPLE_SIZE: i32 = 100000;
pub const FRAG_LENGTH_DIST_MIN_QUAL: u8 = 60;
pub const FRAG_LENGTH_DIST_MAX_LENGTH: i32 = 1500;
pub const FRAG_LENGTH_DIST_PERCENTILE: f64 = 0.9975;
pub const FRAG_LENGTH_1_STD_DEV_PERCENTILE: f64 = 0.16;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FragmentLengthBounds {
    pub lower_bound: i32,
    pub upper_bound: i32,
    pub median: i32,
    pub std_deviation: f64,
}

impl FragmentLengthBounds {
    pub fn is_valid(&self) -> bool {
        self.lower_bound > 0 && self.upper_bound > self.lower_bound
    }
}

/// calculateFragmentLengthBounds over (length, frequency) in ascending length
pub fn calculate_fragment_length_bounds(length_frequencies: &[(i32, i32)]) -> FragmentLengthBounds {
    if length_frequencies.is_empty() {
        return FragmentLengthBounds { lower_bound: 0, upper_bound: 0, median: 0, std_deviation: 0.0 };
    }
    let total: i32 = length_frequencies.iter().map(|x| x.1).sum();
    let mut cumulative = 0;
    let required_min_total = (total as f64 * (1.0 - FRAG_LENGTH_DIST_PERCENTILE)).floor() as i32;
    let required_std_dev_total = (total as f64 * FRAG_LENGTH_1_STD_DEV_PERCENTILE).floor() as i32;
    let required_max_total = (total as f64 * FRAG_LENGTH_DIST_PERCENTILE).floor() as i32;
    let median_fragment = total / 2;
    let (mut lower_bound, mut upper_bound, mut std_dev_length, mut median) = (0, 0, 0, 0);
    for &(length, frequency) in length_frequencies {
        if lower_bound == 0 && frequency + cumulative >= required_min_total {
            lower_bound = length;
        }
        if std_dev_length == 0 && frequency + cumulative >= required_std_dev_total {
            std_dev_length = length;
        }
        if median == 0 && frequency + cumulative >= median_fragment {
            median = length;
        }
        if upper_bound == 0 && frequency + cumulative >= required_max_total {
            upper_bound = length;
        } else {
            cumulative += frequency;
        }
    }
    FragmentLengthBounds { lower_bound, upper_bound, median, std_deviation: (median - std_dev_length) as f64 }
}

#[derive(Default)]
struct ChromosomeSample {
    lengths: BTreeMap<i32, i32>,
    max_read_length: i32,
    has_paired_reads: bool,
}

fn is_candidate_record(rec: &bam::Record) -> bool {
    let flags = rec.flags();
    if flags & FLAG_DUPLICATE != 0 {
        return false;
    }
    let mut fragment_length = inferred_insert_size_abs(rec);
    if fragment_length > FRAG_LENGTH_DIST_MAX_LENGTH {
        return false;
    }
    let read_length = rec.seq_len() as i32;
    let paired = flags & FLAG_PAIRED != 0;
    if !paired && fragment_length == 0 {
        fragment_length = read_length;
    }
    if fragment_length < read_length {
        return false;
    }
    // fully aligned reads only
    let cigar = rec.raw_cigar();
    if !(cigar.len() == 1 && op(cigar[0]) == OP_M && oplen(cigar[0]) == read_length) {
        return false;
    }
    if aux_int(rec, b"NM").is_some_and(|nm| nm as i32 > 2) {
        return false;
    }
    if paired {
        if flags & FLAG_SECOND != 0 {
            return false;
        }
        // no translocations or inversions
        if rec.mtid() != rec.tid() {
            return false;
        }
        if (flags & FLAG_MATE_REVERSE != 0) == (flags & FLAG_REVERSE != 0) {
            return false;
        }
        if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            return false;
        }
        if aux_str(rec, b"MC").is_some_and(|mc| !mate_cigar_is_aligned(mc, read_length as usize)) {
            return false;
        }
    }
    true
}

/// round to the nearest unit up to 1000, then to 10s up to 3000, then to 100s
fn length_bucket(fragment_length: i32) -> i32 {
    if fragment_length < 1000 {
        fragment_length
    } else if fragment_length < 3000 {
        10 * (fragment_length as f64 / 10.0 + 0.5).floor() as i32
    } else {
        100 * (fragment_length as f64 / 100.0 + 0.5).floor() as i32
    }
}

fn sample_chromosome(config: &PrepConfig, chromosome: &str) -> Result<ChromosomeSample> {
    let mut sample = ChromosomeSample::default();
    let (start, end) = if !config.specific.regions.is_empty() {
        match config.specific.regions.iter().find(|r| r.0 == chromosome) {
            Some(r) => (r.1, r.2),
            None => return Ok(sample),
        }
    } else {
        (1_000_000, 10_000_000)
    };
    let mut reader = crate::prep::open_indexed_reader(&config.bam_files[0], &config.ref_genome_file)?;
    let Some(tid) = reader.header().tid(chromosome.as_bytes()) else {
        crate::warn!("cannot find sequence index for chromosome {chromosome} in bam header");
        return Ok(sample);
    };
    reader.fetch((tid, (start - 1) as i64, end as i64))?;
    let mut rec = bam::Record::new();
    let mut processed = 0;
    while let Some(r) = reader.read(&mut rec) {
        r?;
        // the slicer's filter: mapping quality, no unmapped, secondary, supplementary or duplicate reads
        let flags = rec.flags();
        if rec.mapq() < FRAG_LENGTH_DIST_MIN_QUAL || flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_DUPLICATE) != 0 {
            continue;
        }
        if !is_candidate_record(&rec) {
            continue;
        }
        sample.max_read_length = sample.max_read_length.max(rec.seq_len() as i32);
        sample.has_paired_reads |= flags & FLAG_PAIRED != 0;
        processed += 1;
        if processed >= FRAG_LENGTH_DIST_SAMPLE_SIZE {
            crate::debug!("chr({chromosome}) reached max fragment count");
            break;
        }
        let bucket = length_bucket(inferred_insert_size_abs(&rec));
        if bucket > 0 {
            *sample.lengths.entry(bucket).or_insert(0) += 1;
        }
    }
    Ok(sample)
}

pub struct FragmentSizeDistribution {
    pub length_frequencies: Vec<(i32, i32)>,
    pub max_read_length: i32,
    pub has_paired_reads: bool,
}

impl FragmentSizeDistribution {
    /// Samples the chromosomes (1 to 10, or the specific ones) on the configured threads and merges the counts.
    pub fn run(config: &PrepConfig) -> Result<FragmentSizeDistribution> {
        crate::info!("calculating fragment size distribution");
        let chromosomes: Vec<String> = if config.specific.has_filters() {
            config.specific.chromosomes.clone()
        } else {
            (1..=10).map(|i| versioned_chromosome(config.v38, &i.to_string())).collect()
        };
        let next = std::sync::atomic::AtomicUsize::new(0);
        let results: Vec<std::sync::Mutex<Option<Result<ChromosomeSample>>>> = chromosomes.iter().map(|_| std::sync::Mutex::new(None)).collect();
        std::thread::scope(|scope| {
            for _ in 0..config.threads.min(chromosomes.len()).max(1) {
                scope.spawn(|| loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= chromosomes.len() {
                        break;
                    }
                    *results[i].lock().unwrap() = Some(sample_chromosome(config, &chromosomes[i]));
                });
            }
        });
        let mut merged: BTreeMap<i32, i32> = BTreeMap::new();
        let mut dist = FragmentSizeDistribution { length_frequencies: Vec::new(), max_read_length: 0, has_paired_reads: false };
        for r in results {
            let Some(sample) = r.into_inner().unwrap() else { continue };
            let sample = sample?;
            for (length, count) in sample.lengths {
                *merged.entry(length).or_insert(0) += count;
            }
            dist.max_read_length = dist.max_read_length.max(sample.max_read_length);
            dist.has_paired_reads |= sample.has_paired_reads;
        }
        dist.length_frequencies = merged.into_iter().collect();
        crate::info!("maxReadLength({}) {}", dist.max_read_length, if !dist.has_paired_reads { "unpaired reads" } else { "" });
        match (dist.length_frequencies.first(), dist.length_frequencies.last()) {
            (Some(min), Some(max)) => crate::debug!("fragment size distribution complete: min({}) max({})", min.0, max.0),
            _ => crate::debug!("no fragment lengths recorded"),
        }
        if config.has_write_type(WriteType::FragmentLengthDist) {
            dist.write(config)?;
        }
        Ok(dist)
    }

    fn write(&self, config: &PrepConfig) -> Result<()> {
        let mut text = String::from("FragmentLength\tCount\n");
        for &(length, frequency) in &self.length_frequencies {
            // entries below the read length come from trimmed fragments
            if self.has_paired_reads && self.max_read_length > 0 && length < self.max_read_length {
                continue;
            }
            text.push_str(&format!("{length}\t{frequency}\n"));
        }
        std::fs::write(config.form_filename(WriteType::FragmentLengthDist, config.sample_id()), text)?;
        Ok(())
    }

    pub fn bounds(&self) -> FragmentLengthBounds {
        calculate_fragment_length_bounds(&self.length_frequencies)
    }
}
