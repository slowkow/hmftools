//! GeneCoverage and ExonCoverage: the read depth over the coding exons (with splice margins) of the
//! canonical transcripts of the driver panel's genes, and the gene and exon coverage files.

use anyhow::{bail, Result};

use hmf_common::driver_panel::{self, GeneRegion};
use hmf_common::ensembl::EnsemblDataCache;
use hmf_common::java_hash::JavaHashMap;
use hmf_common::metrics::{DepthRange, ExonDepth, GeneDepth, GeneDepthFile};
use hmf_common::regions::positions_overlap;
use hmf_common::segmentation::HumanChromosome;

use crate::config::MetricsConfig;

pub const MAX_DEPTH_BUCKET: i32 = 10000;
pub const GENE_COVERAGE_MIN_MAP_QUALITY: u8 = 10;

/// The depth buckets: 1 up to 30, 10 up to 100, 50 up to 500, 100 up to 2000, 1000 up to 10000
pub fn depth_buckets() -> Vec<i32> {
    let mut b = Vec::new();
    let mut d = 0;
    while d < 30 {
        b.push(d);
        d += 1;
    }
    while d < 100 {
        b.push(d);
        d += 10;
    }
    while d < 500 {
        b.push(d);
        d += 50;
    }
    while d < 2000 {
        b.push(d);
        d += 100;
    }
    while d <= 10000 {
        b.push(d);
        d += 1000;
    }
    b
}

/// DEPTH_RANGES: [bucket, next bucket - 1], the last one [10000, 10000]
pub fn depth_ranges() -> Vec<DepthRange> {
    let b = depth_buckets();
    let mut r: Vec<DepthRange> = b.windows(2).map(|w| DepthRange { lower: w[0], upper: w[1] - 1 }).collect();
    let last = *b.last().unwrap();
    r.push(DepthRange { lower: last, upper: last });
    r
}

/// GeneCoverage.bucket: the index of the bucket a depth falls in
#[cfg(test)]
pub fn bucket(depth: i32) -> usize {
    bucket_of(depth, &depth_buckets())
}

fn bucket_of(depth: i32, buckets: &[i32]) -> usize {
    for (i, &b) in buckets.iter().enumerate() {
        if depth <= b {
            return if depth < b { i.wrapping_sub(1) } else { i };
        }
    }
    buckets.len() - 1
}

/// A depth to bucket table for depths up to MAX_DEPTH_BUCKET (above it, the last bucket)
fn bucket_table() -> Vec<u8> {
    let b = depth_buckets();
    (0..=MAX_DEPTH_BUCKET).map(|d| bucket_of(d, &b) as u8).collect()
}

/// An exon region of the panel
#[derive(Clone, Debug)]
pub struct Exon {
    pub region: GeneRegion,
    /// index into the per-exon depth arrays
    pub index: usize,
}

/// GeneCoverage: the panel's exon regions by chromosome, iterated as Java's HashMap
pub struct GeneCoverage {
    pub chr_exons: JavaHashMap<String, Vec<Exon>>,
    pub exon_count: usize,
    /// per chromosome: the exons sorted by start (their list order) and the running maximum of their ends
    max_ends: JavaHashMap<String, Vec<i32>>,
}

impl GeneCoverage {
    pub fn empty() -> Self {
        GeneCoverage { chr_exons: JavaHashMap::new(), exon_count: 0, max_ends: JavaHashMap::new() }
    }

    /// new GeneCoverage(configBuilder): nothing unless both the Ensembl data and the panel are given
    pub fn load(config: &MetricsConfig) -> Result<Self> {
        let (Some(ensembl_dir), Some(_)) = (&config.ensembl_data_dir, &config.driver_gene_panel) else {
            return Ok(Self::empty());
        };
        let mut cache = EnsemblDataCache::new(ensembl_dir, config.gene_ref_gen_version);
        cache.set_required_data(true, false, false, true);
        // Java ignores load()'s result; a missing file leaves the cache empty, then every panel gene is invalid
        if let Err(e) = cache.load(false) {
            crate::warn!("Ensembl data cache: {e:#}");
        }
        cache.create_gene_name_id_map();
        let driver_genes = driver_panel::load_driver_genes(config.driver_gene_panel.as_deref())?;
        let invalid = driver_panel::find_invalid_driver_genes(&driver_genes, &cache);
        if !invalid.is_empty() {
            let names: Vec<&str> = invalid.iter().map(|g| g.gene.as_str()).collect();
            crate::error!("Invalid non-Ensembl gene names: {}", names.join(";"));
            bail!("invalid driver genes");
        }
        let genes: Vec<String> = driver_genes.iter().map(|g| g.gene.clone()).collect();
        let regions = driver_panel::build_driver_gene_regions(&cache, &genes)?;
        Ok(Self::from_regions(&regions))
    }

    /// the GeneCoverage map of the panel regions, in their iteration order
    pub fn from_regions(regions: &JavaHashMap<String, Vec<GeneRegion>>) -> Self {
        let mut chr_exons: JavaHashMap<String, Vec<Exon>> = JavaHashMap::new();
        let mut max_ends: JavaHashMap<String, Vec<i32>> = JavaHashMap::new();
        let mut index = 0;
        for (chr, list) in regions.iter() {
            let exons: Vec<Exon> = list
                .iter()
                .map(|r| {
                    let e = Exon { region: r.clone(), index };
                    index += 1;
                    e
                })
                .collect();
            let mut running = i32::MIN;
            let ends: Vec<i32> = exons
                .iter()
                .map(|e| {
                    running = running.max(e.region.end);
                    running
                })
                .collect();
            max_ends.insert(chr.clone(), ends);
            chr_exons.insert(chr.clone(), exons);
        }
        GeneCoverage { chr_exons, exon_count: index, max_ends }
    }

    /// the exons of a chromosome, sorted by start, with the running maximum of their ends
    pub fn chromosome(&self, chromosome: &str) -> Option<(&[Exon], &[i32])> {
        let key = chromosome.to_string();
        let exons = self.chr_exons.get(&key)?;
        Some((exons, self.max_ends.get(&key).unwrap()))
    }

    /// GeneCoverage.getGeneRegions(region): the indices (into the chromosome's list) of the exons
    /// overlapping the partition, ChrBaseRegion.overlaps comparing the chromosomes as HumanChromosome
    pub fn partition_exons(&self, chromosome: &str, start: i32, end: i32) -> Vec<usize> {
        let Some((exons, _)) = self.chromosome(chromosome) else { return Vec::new() };
        let region_chr = HumanChromosome::from_name(chromosome);
        exons
            .iter()
            .enumerate()
            .filter(|(_, e)| HumanChromosome::from_name(&e.region.chromosome) == region_chr && positions_overlap(start, end, e.region.start, e.region.end))
            .map(|(i, _)| i)
            .collect()
    }

    /// GeneCoverage.writeResults: the exon and gene coverage files from the summed depths
    pub fn write_results(&self, output_dir: &str, sample_id: &str, depths: &[Vec<u32>]) -> Result<()> {
        if self.chr_exons.is_empty() {
            return Ok(());
        }
        let ranges = depth_ranges();
        let table = bucket_table();
        let last = ranges.len() - 1;
        let counts_of = |d: &[u32], counts: &mut [i32]| {
            for &v in d {
                let b = if v as i64 > MAX_DEPTH_BUCKET as i64 { last } else { table[v as usize] as usize };
                counts[b] += 1;
            }
        };
        // `depths` holds each exon's depth array (an empty one for an exon no read reached)
        let depth_of = |e: &Exon| -> std::borrow::Cow<'_, [u32]> {
            match depths.get(e.index) {
                Some(d) if !d.is_empty() => std::borrow::Cow::Borrowed(&d[..]),
                _ => std::borrow::Cow::Owned(vec![0u32; e.region.base_length().max(1) as usize]),
            }
        };
        // the exons in the map's iteration order
        let mut exon_depths: Vec<ExonDepth> = Vec::new();
        // genes keyed by name, their exons in that order (only the first one's chromosome is written)
        let mut gene_exons: JavaHashMap<String, Vec<(String, i32, i32, Vec<i32>)>> = JavaHashMap::new();
        for (_, exons) in self.chr_exons.iter() {
            for e in exons {
                let d = depth_of(e);
                let mut counts = vec![0i32; ranges.len()];
                counts_of(&d, &mut counts);
                let mut sorted: Vec<u32> = d.to_vec();
                sorted.sort_unstable();
                let median = sorted[sorted.len() / 2] as f64;
                let mut total = 0.0f64;
                for &v in d.iter() {
                    total += v as f64;
                }
                let mean = total / d.len() as f64;
                exon_depths.push(ExonDepth {
                    depth: GeneDepth {
                        gene: e.region.gene_name.clone(),
                        chromosome: e.region.chromosome.clone(),
                        pos_start: e.region.start,
                        pos_end: e.region.end,
                        depth_counts: ranges.iter().copied().zip(counts.iter().copied()).collect(),
                    },
                    exon_rank: e.region.exon_rank,
                    median_depth: median,
                    mean_depth: mean,
                });
                let entry = (e.region.chromosome.clone(), e.region.start, e.region.end, counts);
                match gene_exons.get_mut(&e.region.gene_name) {
                    Some(l) => l.push(entry),
                    None => {
                        gene_exons.insert(e.region.gene_name.clone(), vec![entry]);
                    }
                }
            }
        }
        GeneDepthFile::write_exon_coverage(&GeneDepthFile::generate_exon_coverage_filename(output_dir, sample_id), &mut exon_depths)?;
        let mut gene_depths: Vec<GeneDepth> = Vec::new();
        for (gene, exons) in gene_exons.iter() {
            let mut counts = vec![0i32; ranges.len()];
            for (_, _, _, c) in exons {
                for (a, b) in counts.iter_mut().zip(c) {
                    *a += *b;
                }
            }
            let min = exons.iter().map(|e| e.1).min().unwrap_or(0);
            let max = exons.iter().map(|e| e.2).max().unwrap_or(0);
            gene_depths.push(GeneDepth {
                gene: gene.clone(),
                chromosome: exons[0].0.clone(),
                pos_start: min,
                pos_end: max,
                depth_counts: ranges.iter().copied().zip(counts).collect(),
            });
        }
        GeneDepthFile::write_gene_coverage(&GeneDepthFile::generate_gene_coverage_filename(output_dir, sample_id), &mut gene_depths)?;
        Ok(())
    }
}

/// ExonCoverage.processRead: one more read over the exon's bases within [read_start, read_end]
#[inline]
pub fn add_read(depth: &mut [u32], exon_start: i32, exon_end: i32, read_start: i32, read_end: i32) {
    let a = exon_start.max(read_start);
    let b = exon_end.min(read_end);
    if a > b {
        return;
    }
    for d in &mut depth[(a - exon_start) as usize..=(b - exon_start) as usize] {
        *d += 1;
    }
}

#[cfg(test)]
mod tests {
    //! GeneCoverageTest
    use super::*;

    fn assert_bucket(min: i32, max: i32, expected: usize) {
        for d in min..max {
            assert_eq!(bucket(d), expected, "depth {d}");
        }
    }

    #[test]
    fn test_bucket() {
        for d in 0..30 {
            assert_eq!(bucket(d), d as usize);
        }
        assert_bucket(30, 39, 30);
        assert_bucket(40, 49, 31);
        assert_bucket(50, 59, 32);
        assert_bucket(60, 69, 33);
        assert_bucket(70, 79, 34);
        assert_bucket(80, 89, 35);
        assert_bucket(90, 99, 36);
        assert_bucket(100, 149, 37);
        assert_bucket(150, 200, 38);
        assert_bucket(500, 600, 45);
        assert_bucket(1900, 2000, 59);
        assert_bucket(9000, 10000, 67);
        assert_eq!(bucket(10000), 68);
        assert_eq!(bucket(MAX_DEPTH_BUCKET), 68);
        assert_eq!(bucket(20000), 68);
        let table = bucket_table();
        for d in [0, 29, 30, 99, 100, 451, 1999, 2000, 9999, 10000] {
            assert_eq!(table[d as usize] as usize, bucket(d));
        }
        assert_eq!(depth_ranges().len(), 69);
    }

    fn exon() -> Vec<u32> {
        vec![0; 5] // 100-104
    }

    #[test]
    fn test_alignments() {
        let mut v = exon();
        add_read(&mut v, 100, 104, 30, 99);
        assert_eq!(v, vec![0, 0, 0, 0, 0]);
        let mut v = exon();
        add_read(&mut v, 100, 104, 90, 100);
        assert_eq!(v, vec![1, 0, 0, 0, 0]);
        add_read(&mut v, 100, 104, 91, 102);
        assert_eq!(v, vec![2, 1, 1, 0, 0]);
        let mut v = exon();
        add_read(&mut v, 100, 104, 103, 110);
        assert_eq!(v, vec![0, 0, 0, 1, 1]);
        add_read(&mut v, 100, 104, 104, 110);
        assert_eq!(v, vec![0, 0, 0, 1, 2]);
        let mut v = exon();
        add_read(&mut v, 100, 104, 90, 110);
        add_read(&mut v, 100, 104, 90, 110);
        assert_eq!(v, vec![2, 2, 2, 2, 2]);
        let mut v = exon();
        add_read(&mut v, 100, 104, 101, 103);
        assert_eq!(v, vec![0, 1, 1, 1, 0]);
        add_read(&mut v, 100, 104, 101, 105);
        assert_eq!(v, vec![0, 2, 2, 2, 1]);
        add_read(&mut v, 100, 104, 99, 101);
        assert_eq!(v, vec![1, 3, 2, 2, 1]);
        let mut v = exon();
        add_read(&mut v, 100, 104, 106, 200);
        assert_eq!(v, vec![0, 0, 0, 0, 0]);
    }
}
