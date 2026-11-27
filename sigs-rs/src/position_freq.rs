//! hmf-common sigs.PositionFrequencies and SIGS's PositionFreqBuilder.writeSampleCounts: SNV counts per
//! genomic position bucket (`-position_bucket_size`), written to `<sample>.sig.pos_freq_counts.csv`.

use std::collections::BTreeMap;

use anyhow::Result;

use crate::compat;

pub const HUMAN_CHROMOSOMES: [&str; 24] = hmf_common::regions::HUMAN_CHROMOSOMES;

/// RefGenomeVersion.versionedChromosome
pub fn versioned_chromosome(v38: bool, chromosome: &str) -> String {
    if v38 {
        hmf_common::regions::enforce_chr_prefix(chromosome)
    } else {
        hmf_common::regions::strip_chr_prefix(chromosome).to_string()
    }
}

/// PositionFrequencies.buildStandardChromosomeLengths: GRCh38 lengths, or for GRCh37 the larger of the
/// two versions' lengths (as the reference data of version 016 and before was built)
pub fn standard_chromosome_lengths(v38: bool) -> Vec<(String, i32)> {
    HUMAN_CHROMOSOMES
        .iter()
        .map(|c| {
            let l38 = hmf_common::ref_coords::chromosome_length(c, false).unwrap();
            let len = if v38 { l38 } else { hmf_common::ref_coords::chromosome_length(c, true).unwrap().max(l38) };
            (versioned_chromosome(v38, c), len)
        })
        .collect()
}

/// PositionFrequencies.initialisePositionCache: the first bucket index of each chromosome, and the number
/// of buckets.
///
/// Java starts a chromosome at the previous one's last index plus one only if that index is above 0, so a
/// first chromosome of a single bucket (a bucket size of at least its length) is followed by chromosomes
/// starting at index 0 as well (SI2); corrected, each starts after the previous one.
pub fn initialise_position_cache(v38: bool, bucket_size: i32, lengths: &[(String, i32)]) -> (BTreeMap<String, i32>, i32) {
    let mut index = BTreeMap::new();
    let mut cache_size = 0;
    if bucket_size == 0 {
        return (index, cache_size);
    }
    let mut last_end: i64 = -1;
    for c in HUMAN_CHROMOSOMES {
        let chr = versioned_chromosome(v38, c);
        let Some((_, length)) = lengths.iter().find(|(n, _)| *n == chr) else { continue };
        let start = if compat::java("SI2") {
            if last_end > 0 { last_end + 1 } else { 0 }
        } else {
            last_end + 1
        };
        index.insert(chr, start as i32);
        let count = (*length as f64 / bucket_size as f64).ceil() as i32;
        cache_size += count;
        last_end = start + count as i64 - 1;
    }
    (index, cache_size)
}

/// PositionFrequencies.getBucketIndex: -1 for a chromosome not in the index
pub fn bucket_index(bucket_size: i32, index: &BTreeMap<String, i32>, chromosome: &str, position: i32) -> i64 {
    match index.get(chromosome) {
        None => -1,
        Some(start) => *start as i64 + (position as f64 / bucket_size as f64).floor() as i64,
    }
}

/// PositionFrequencies.getChromosomeFromIndex
pub fn chromosome_from_index(v38: bool, index: &BTreeMap<String, i32>, bucket: i32) -> String {
    let mut last_start = -1;
    let mut last_chr = String::new();
    for c in HUMAN_CHROMOSOMES {
        let chr = versioned_chromosome(v38, c);
        let start = *index.get(&chr).unwrap_or(&0);
        if last_start >= 0 && bucket >= last_start && bucket < start {
            return last_chr;
        }
        last_start = start;
        last_chr = chr;
    }
    last_chr
}

/// The position frequencies SigSnvLoader keeps for a sample.
///
/// Java builds them as `new PositionFrequencies(V37, bucketSize)`, which keeps no per-bucket map
/// (buildMap false), so the counts the writer reads are always empty and every bucket is written as 0
/// (SI1). Corrected, the map is kept, and the writer matches its chromosome names to those of
/// `-ref_genome_version` with or without the "chr" prefix.
pub struct PositionFrequencies {
    pub bucket_size: i32,
    build_map: bool,
    /// chromosome as in the VCF -> position bucket (position / size * size) -> count
    pub chr_pos_bucket_frequencies: BTreeMap<String, BTreeMap<i32, i32>>,
}

impl PositionFrequencies {
    pub fn new(bucket_size: i32) -> PositionFrequencies {
        PositionFrequencies { bucket_size, build_map: !compat::java("SI1"), chr_pos_bucket_frequencies: BTreeMap::new() }
    }

    pub fn add_position(&mut self, chromosome: &str, position: i32) {
        if !self.build_map {
            return;
        }
        let bucket = position / self.bucket_size * self.bucket_size;
        *self.chr_pos_bucket_frequencies.entry(chromosome.to_string()).or_default().entry(bucket).or_insert(0) += 1;
    }
}

/// PositionFreqBuilder(refGenomeVersion, outputDir, outputId, bucketSize).writeSampleCounts
pub fn write_sample_counts(
    v38: bool,
    output_dir: &str,
    bucket_size: i32,
    sample_id: &str,
    frequencies: &PositionFrequencies,
    write_coords: bool,
) -> Result<()> {
    let lengths = standard_chromosome_lengths(v38);
    let (index, cache_size) = initialise_position_cache(v38, bucket_size, &lengths);
    hmf_common::info!("position cache size({}) from position bucket position({})", cache_size, bucket_size);
    hmf_common::info!("sample({}) writing position frequencies data", sample_id);

    let filename = format!("{output_dir}{sample_id}.sig.pos_freq_counts.csv");
    let mut buckets = vec![0i32; cache_size.max(0) as usize];
    let mut skipped = 0;
    for (chromosome, positions) in &frequencies.chr_pos_bucket_frequencies {
        // corrected (SI1), the VCF's chromosome names are matched to the version's with or without "chr"
        let chromosome = if crate::compat::java("SI1") { chromosome.clone() } else { versioned_chromosome(v38, chromosome) };
        for (position, frequency) in positions {
            let b = bucket_index(bucket_size, &index, &chromosome, *position);
            if b < 0 || b >= buckets.len() as i64 {
                // Java fails here (ArrayIndexOutOfBoundsException); only reached with SI1 corrected
                skipped += frequency;
                continue;
            }
            // Java assigns (a chromosome's buckets are its own); with names matched, "chr1" and "1" add up
            if crate::compat::java("SI1") {
                buckets[b as usize] = *frequency;
            } else {
                buckets[b as usize] += *frequency;
            }
        }
    }
    if skipped > 0 {
        hmf_common::warn!(
            "{} variants on chromosomes or positions outside the position buckets of ref genome version {} were not counted",
            skipped,
            if v38 { "V38" } else { "V37" }
        );
    }

    let mut out = String::with_capacity(buckets.len() * 8);
    if write_coords {
        out.push_str("Chromosome,Position,Frequency\n");
    } else {
        out.push_str(sample_id);
        out.push('\n');
    }
    for (b, frequency) in buckets.iter().enumerate() {
        if write_coords {
            let chromosome = chromosome_from_index(v38, &index, b as i32);
            // an int in Java: with SI2 the products overflow and wrap, as they do there
            let position = (b as i32).wrapping_sub(index.get(&chromosome).copied().unwrap_or(0)).wrapping_mul(bucket_size);
            out.push_str(&format!("{chromosome},{position},{frequency}\n"));
        } else {
            out.push_str(&format!("{frequency}\n"));
        }
    }
    if let Err(e) = std::fs::write(&filename, out) {
        hmf_common::error!("failed to write sample pos data output: java.io.FileNotFoundException: {filename} ({e})");
    }
    Ok(())
}
