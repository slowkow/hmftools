//! Java UnmappedRegionDepth: after all counts are in, a variant in a high-depth unmap region whose
//! allele fraction is 0.9 or more has each sample's reference depth raised to that sample's median.

use super::task::{ALLELE_FRACTION, REF_DEPTH, REF_DEPTH_PAIR, TOTAL_FRAGS};
use super::fmt::vcf_double;
use super::vcf::{parse_java_int, VariantContext};
use crate::common::compat;
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

    /// setUnmappedRegionsDepth over the variants of all tasks, in task order. `genotype_ids` are the VCF
    /// genotype columns of the configured samples, in sample order; `ref_tag` and `ref_pair_tag` are the
    /// tags the counts were written under (with the -vcf_tag_prefix, if any). Returns the number of
    /// variants gone over for raising, as logged.
    pub fn set_unmapped_regions_depth(&self, genotype_ids: &[usize], ref_tag: &str, ref_pair_tag: &str, variants: &mut [&mut VariantContext]) -> Result<usize> {
        let sample_count = genotype_ids.len();
        // D8: Java goes over every genotype of the record, taking its column as the sample index, while
        // its lists have one entry per configured sample: an index error when the VCF has more samples
        // than were configured. Corrected: the configured samples, each through its own column.
        let all_genotypes = compat::java("D8");
        // D9: Java writes into the attribute map of a genotype that has none, which is htsjdk's shared
        // unmodifiable map, and fails. Corrected: such a genotype is left as it is, as when the counts
        // are set.
        let fail_on_empty_genotype = compat::java("D9");
        // D16: Java reads and writes REF and REFPAIR without the -vcf_tag_prefix, and lists a variant
        // once for each sample that qualifies.
        let java_tags = compat::java("D16");
        let (ref_tag, ref_pair_tag) = if java_tags { (REF_DEPTH, REF_DEPTH_PAIR) } else { (ref_tag, ref_pair_tag) };

        // (sample index, genotype column) pairs to go over for a record with `genotype_count` genotypes
        let columns = |genotype_count: usize| -> Vec<(usize, usize)> {
            if all_genotypes {
                (0..genotype_count).map(|i| (i, i)).collect()
            } else {
                genotype_ids.iter().copied().enumerate().filter(|(_, column)| *column < genotype_count).collect()
            }
        };

        let mut sample_variant_depths: Vec<Vec<i32>> = vec![Vec::new(); sample_count];
        let mut low_ref_depth_variants: Vec<usize> = Vec::new();
        let index_error = |index: usize| anyhow!("java.lang.IndexOutOfBoundsException: Index {index} out of bounds for length {sample_count}");

        for (variant_index, variant) in variants.iter().enumerate() {
            for (sample_index, column) in columns(variant.genotypes.len()) {
                let genotype = &variant.genotypes[column];
                let ref_depth = genotype.attribute_as_int(ref_tag, 0)?;
                let ref_pair_depth = genotype.attribute_as_int(ref_pair_tag, 0)?;
                let tumor_frags = genotype.attribute_as_int(TOTAL_FRAGS, 0)?;
                // 0 / 0 gives NaN, which is not below the threshold
                let af = tumor_frags as f64 / ref_depth.wrapping_add(ref_pair_depth).wrapping_add(tumor_frags) as f64;

                if af < UNMAPPED_VAF_THRESHOLD {
                    sample_variant_depths.get_mut(sample_index).ok_or_else(|| index_error(sample_index))?.push(ref_depth);
                    continue;
                }
                if self.in_high_depth_unmapped_region(&variant.chromosome, variant.position)? && (java_tags || low_ref_depth_variants.last() != Some(&variant_index)) {
                    low_ref_depth_variants.push(variant_index);
                }
            }
        }

        if low_ref_depth_variants.is_empty() {
            return Ok(0);
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

        let listed = low_ref_depth_variants.len();
        info!("setting ref depth to median for {} variants", listed);

        for variant_index in low_ref_depth_variants {
            let variant = &mut *variants[variant_index];
            let mut total_ref_depth = 0i32;

            for (sample_index, column) in columns(variant.genotypes.len()) {
                let genotype = &mut variant.genotypes[column];
                let ref_depth = genotype.attribute_as_int(ref_tag, 0)?;
                let median_depth = *sample_median_depths.get(sample_index).ok_or_else(|| index_error(sample_index))?;

                if ref_depth < median_depth {
                    if genotype.attributes.is_empty() {
                        if fail_on_empty_genotype {
                            bail!("java.lang.UnsupportedOperationException: genotype without attributes at {}:{}", variant.chromosome, variant.position);
                        }
                        continue;
                    }
                    total_ref_depth = total_ref_depth.wrapping_add(median_depth);
                    genotype.attributes.insert(ref_tag.to_string(), median_depth.to_string());

                    let variant_frags = genotype.attribute_as_int(TOTAL_FRAGS, 0)?;
                    let ref_pair_depth = genotype.attribute_as_int(ref_pair_tag, 0)?;
                    let total = variant_frags.wrapping_add(median_depth).wrapping_add(ref_pair_depth) as f64;
                    let af = if total > 0.0 { variant_frags as f64 / total } else { 0.0 };
                    genotype.attributes.insert(ALLELE_FRACTION.to_string(), vcf_double(af));
                } else {
                    total_ref_depth = total_ref_depth.wrapping_add(ref_depth);
                }
            }

            variant.info.insert(ref_tag.to_string(), total_ref_depth.to_string());
        }
        Ok(listed)
    }
}
