//! RegionsFilter: drops the positions inside any of a set of regions (here the HLA, IG and TR regions),
//! walking each chromosome's regions alongside positions that are in order within a chromosome.

use hmf_common::regions::ChrBaseRegion;
use std::collections::HashMap;

pub struct RegionsFilter {
    regions: HashMap<String, Vec<ChrBaseRegion>>,
}

impl RegionsFilter {
    pub fn new(regions_to_exclude: &[ChrBaseRegion]) -> Self {
        let mut regions: HashMap<String, Vec<ChrBaseRegion>> = HashMap::new();
        for region in regions_to_exclude {
            regions.entry(region.chromosome.clone()).or_default().push(region.clone());
        }
        // Collections.sort: ChrBaseRegion.compareTo within a chromosome is by start only; the sort is stable
        for list in regions.values_mut() {
            list.sort_by_key(|r| r.start);
        }
        RegionsFilter { regions }
    }

    /// filter: the indices of the positions kept, in order. As in Java, the region index only moves
    /// forward within a run of one chromosome, so positions out of order in a run may be missed.
    pub fn filter<'p>(&self, positions: impl IntoIterator<Item = (&'p str, i32)>) -> Vec<u32> {
        let mut result = Vec::new();
        let mut current_chromosome: Option<&str> = None;
        let mut current_regions: Option<&Vec<ChrBaseRegion>> = None;
        let mut region_index = 0usize;
        for (i, (chromosome, position)) in positions.into_iter().enumerate() {
            if current_chromosome != Some(chromosome) {
                current_chromosome = Some(chromosome);
                current_regions = self.regions.get(chromosome);
                region_index = 0;
            }
            let Some(regions) = current_regions else {
                result.push(i as u32);
                continue;
            };
            while region_index < regions.len() && regions[region_index].end < position {
                region_index += 1;
            }
            if region_index >= regions.len() || !regions[region_index].contains_position(position) {
                result.push(i as u32);
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // RegionsFilterTest.multipleRegionsTest
    #[test]
    fn multiple_regions() {
        let r = |c: &str, s, e| ChrBaseRegion::new(c, s, e);
        let regions = vec![
            r("3", 40_001, 50_000),
            r("2", 40_001, 50_000),
            r("1", 40_001, 50_000),
            r("2", 10_001, 20_000),
            r("3", 10_001, 20_000),
            r("1", 10_001, 20_000),
            r("2", 1001, 2000),
            r("3", 1001, 2000),
            r("1", 1001, 2000),
        ];
        let filter = RegionsFilter::new(&regions);
        let positions: Vec<(&str, i32)> = vec![
            ("1", 100),
            ("1", 200),
            ("1", 900),
            ("1", 1000),
            ("1", 1001),
            ("1", 1200),
            ("1", 1300),
            ("1", 1900),
            ("1", 2000),
            ("1", 2001),
            ("1", 20_001),
            ("1", 30_001),
            ("1", 40_001),
            ("1", 50_001),
            ("1", 60_001),
            ("3", 500),
            ("3", 1500),
            ("3", 2500),
            ("3", 9500),
            ("3", 10_500),
            ("3", 15_500),
            ("3", 20_500),
            ("3", 30_500),
            ("3", 40_500),
            ("3", 50_500),
            ("5", 2500),
        ];
        let kept = filter.filter(positions.iter().copied());
        assert_eq!(kept.len(), 16);
        let expected_pairs = [(0, 0), (1, 1), (3, 3), (9, 4), (10, 5), (11, 6), (13, 7), (14, 8), (15, 9), (17, 10), (18, 11), (21, 12), (22, 13), (24, 14), (25, 15)];
        for (input, output) in expected_pairs {
            assert_eq!(kept[output], input);
        }
    }
}
