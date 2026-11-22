//! Java UnmappedRegionDepth: after all counts are in, a variant in a high-depth unmap region whose
//! allele fraction is 0.9 or more has each sample's reference depth raised to that sample's median.

use super::task::{ALLELE_FRACTION, REF_DEPTH, REF_DEPTH_PAIR, TOTAL_FRAGS};
use super::fmt::vcf_double;
use super::vcf::{parse_java_int, VariantContext};
use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;

const UNMAPPED_VAF_THRESHOLD: f64 = 0.9;

struct UnmappingRegion {
    start: i32,
    end: i32,
    max_depth: i32,
}

pub struct UnmappedRegionDepth {
    /// None when Java's loader returns null (regions that overlap or are out of order)
    regions: Option<HashMap<String, Vec<UnmappingRegion>>>,
}

fn column(names: &[&str], candidates: &[&str]) -> Option<usize> {
    candidates.iter().find_map(|c| names.iter().position(|n| n == c))
}

impl UnmappedRegionDepth {
    /// UnmappedRegions.loadUnmapRegions
    pub fn load(filename: &str) -> Result<UnmappedRegionDepth> {
        if !std::path::Path::new(filename).exists() {
            return Ok(UnmappedRegionDepth { regions: Some(HashMap::new()) });
        }
        let text = std::fs::read_to_string(filename).with_context(|| format!("failed to read high-depth regions file {filename}"))?;
        let delim = if filename.ends_with(".csv") || filename.ends_with(".csv.gz") { ',' } else { '\t' };
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: empty unmap regions file"))?.split(delim).collect();
        let missing = || anyhow!("java.lang.ArrayIndexOutOfBoundsException: unmap regions file lacks a required column");
        let chr_index = column(&header, &["Chromosome", "chromosome"]).ok_or_else(missing)?;
        let start_index = column(&header, &["PositionStart", "positionstart", "PosStart", "posstart"]).ok_or_else(missing)?;
        let end_index = column(&header, &["PositionEnd", "positionend", "PosEnd", "posend"]).ok_or_else(missing)?;
        let depth_index = column(&header, &["MaxDepth"]).ok_or_else(missing)?;

        let mut map: HashMap<String, Vec<UnmappingRegion>> = HashMap::new();
        for line in lines {
            let values: Vec<&str> = line.split(delim).collect();
            let field = |i: usize| values.get(i).copied().ok_or_else(missing);
            let int = |i: usize| -> Result<i32> {
                let v = field(i)?;
                parse_java_int(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
            };
            let region = UnmappingRegion { start: int(start_index)?, end: int(end_index)?, max_depth: int(depth_index)? };
            let regions = map.entry(field(chr_index)?.to_string()).or_default();
            if let Some(last) = regions.last() {
                if !(last.start > region.end || last.end < region.start) {
                    error!("unmap regions overlap: current({}-{} depth({})) next({}-{} depth({}))", last.start, last.end, last.max_depth, region.start, region.end, region.max_depth);
                    return Ok(UnmappedRegionDepth { regions: None });
                } else if region.end < last.start {
                    error!("unmap regions not sorted: current({}-{} depth({})) next({}-{} depth({}))", last.start, last.end, last.max_depth, region.start, region.end, region.max_depth);
                    return Ok(UnmappedRegionDepth { regions: None });
                }
            }
            regions.push(region);
        }
        Ok(UnmappedRegionDepth { regions: Some(map) })
    }

    fn in_high_depth_unmapped_region(&self, chromosome: &str, position: i32) -> Result<bool> {
        let Some(map) = &self.regions else { bail!("java.lang.NullPointerException: unmap regions were not loaded") };
        Ok(map.get(chromosome).is_some_and(|regions| regions.iter().any(|r| r.max_depth > 0 && position >= r.start && position <= r.end)))
    }

    /// setUnmappedRegionsDepth over the variants of all tasks, in task order. `sample_count` is the
    /// number of configured samples, while the genotypes gone over are all of those in the VCF; Java
    /// fails with an index error where the two disagree, and so does this.
    pub fn set_unmapped_regions_depth(&self, sample_count: usize, variants: &mut [&mut VariantContext]) -> Result<()> {
        let mut sample_variant_depths: Vec<Vec<i32>> = vec![Vec::new(); sample_count];
        // a variant is listed once for each sample that qualifies
        let mut low_ref_depth_variants: Vec<usize> = Vec::new();
        let index_error = |index: usize| anyhow!("java.lang.IndexOutOfBoundsException: Index {index} out of bounds for length {sample_count}");

        for (variant_index, variant) in variants.iter().enumerate() {
            for (sample_index, genotype) in variant.genotypes.iter().enumerate() {
                let ref_depth = genotype.attribute_as_int(REF_DEPTH, 0)?;
                let ref_pair_depth = genotype.attribute_as_int(REF_DEPTH_PAIR, 0)?;
                let tumor_frags = genotype.attribute_as_int(TOTAL_FRAGS, 0)?;
                // 0 / 0 gives NaN, which is not below the threshold
                let af = tumor_frags as f64 / ref_depth.wrapping_add(ref_pair_depth).wrapping_add(tumor_frags) as f64;

                if af < UNMAPPED_VAF_THRESHOLD {
                    sample_variant_depths.get_mut(sample_index).ok_or_else(|| index_error(sample_index))?.push(ref_depth);
                    continue;
                }
                if self.in_high_depth_unmapped_region(&variant.chromosome, variant.position)? {
                    low_ref_depth_variants.push(variant_index);
                }
            }
        }

        if low_ref_depth_variants.is_empty() {
            return Ok(());
        }

        let mut sample_median_depths: Vec<i32> = Vec::new();
        for depth_values in sample_variant_depths.iter_mut() {
            if depth_values.len() < 2 {
                sample_median_depths.push(0);
                continue;
            }
            depth_values.sort_unstable();
            let median_depth = depth_values[depth_values.len() / 2];
            info!("sample median depth({}) from {} variants", median_depth, depth_values.len());
            sample_median_depths.push(median_depth);
        }

        info!("setting ref depth to median for {} variants", low_ref_depth_variants.len());

        for variant_index in low_ref_depth_variants {
            let variant = &mut *variants[variant_index];
            let mut total_ref_depth = 0i32;

            for (sample_index, genotype) in variant.genotypes.iter_mut().enumerate() {
                let ref_depth = genotype.attribute_as_int(REF_DEPTH, 0)?;
                let median_depth = *sample_median_depths.get(sample_index).ok_or_else(|| index_error(sample_index))?;

                if ref_depth < median_depth {
                    total_ref_depth = total_ref_depth.wrapping_add(median_depth);
                    if genotype.attributes.is_empty() {
                        // htsjdk's shared empty attribute map cannot be written to
                        bail!("java.lang.UnsupportedOperationException: genotype without attributes at {}:{}", variant.chromosome, variant.position);
                    }
                    genotype.attributes.insert(REF_DEPTH.to_string(), median_depth.to_string());

                    let variant_frags = genotype.attribute_as_int(TOTAL_FRAGS, 0)?;
                    let ref_pair_depth = genotype.attribute_as_int(REF_DEPTH_PAIR, 0)?;
                    let total = variant_frags.wrapping_add(median_depth).wrapping_add(ref_pair_depth) as f64;
                    let af = if total > 0.0 { variant_frags as f64 / total } else { 0.0 };
                    genotype.attributes.insert(ALLELE_FRACTION.to_string(), vcf_double(af));
                } else {
                    total_ref_depth = total_ref_depth.wrapping_add(ref_depth);
                }
            }

            variant.info.insert(REF_DEPTH.to_string(), total_ref_depth.to_string());
        }
        Ok(())
    }
}
