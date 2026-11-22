//! Known fusion hotspots (HotspotCache, KnownHotspot).

use super::types::*;
use anyhow::{bail, Result};
use std::collections::HashMap;
use std::sync::Arc;

pub struct KnownHotspot {
    pub region_start: Region,
    pub orient_start: Orientation,
    pub region_end: Region,
    pub orient_end: Orientation,
    #[allow(dead_code)]
    pub gene_info: String,
}

impl KnownHotspot {
    pub fn same_orientation(&self) -> bool {
        self.orient_start == self.orient_end
    }
}

#[derive(Default)]
pub struct HotspotCache {
    /// by chromosome of either end
    regions: HashMap<i32, Vec<Arc<KnownHotspot>>>,
    pub count: usize,
}

impl HotspotCache {
    pub fn empty() -> Self {
        Self::default()
    }

    /// The BEDPE file: chromosome, start (0-based), end of each side, name, score, the two strands.
    pub fn load(filename: &str, contigs: &Contigs) -> Result<Self> {
        let text = std::fs::read_to_string(filename).map_err(|e| anyhow::anyhow!("failed to load hotspot data file({filename}): {e}"))?;
        let mut cache = HotspotCache::default();
        for line in text.lines() {
            let v: Vec<&str> = line.split('\t').collect();
            if v.len() < 10 {
                bail!("invalid hotspot entry: {line}");
            }
            let num = |s: &str| -> Result<i32> { s.parse::<i32>().map_err(|_| anyhow::anyhow!("invalid hotspot entry: {line}")) };
            let orient = |s: &str| -> Result<Orientation> {
                s.bytes().next().and_then(Orientation::from_char).ok_or_else(|| anyhow::anyhow!("invalid hotspot entry: {line}"))
            };
            let (chr_start, chr_end) = (contigs.id(v[0]), contigs.id(v[3]));
            let hotspot = Arc::new(KnownHotspot {
                region_start: Region::new(chr_start, num(v[1])? + 1, num(v[2])?),
                orient_start: orient(v[8])?,
                region_end: Region::new(chr_end, num(v[4])? + 1, num(v[5])?),
                orient_end: orient(v[9])?,
                gene_info: v[6].to_string(),
            });
            cache.regions.entry(chr_start).or_default().push(hotspot.clone());
            if v[0] != v[3] {
                cache.regions.entry(chr_end).or_default().push(hotspot);
            }
            cache.count += 1;
        }
        Ok(cache)
    }

    pub fn find_region_hotspots(&self, region: &Region) -> Vec<Arc<KnownHotspot>> {
        match self.regions.get(&region.chr) {
            Some(list) => list.iter().filter(|x| x.region_start.overlaps_region(region) || x.region_end.overlaps_region(region)).cloned().collect(),
            None => Vec::new(),
        }
    }
}
