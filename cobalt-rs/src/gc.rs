//! The GC profile (GCProfileFactory.loadGCContent), the excluded regions (ExcludedRegionsFile,
//! SuppliedExcludedRegions), the tumor-only diploid regions (DiploidRegionLoader) and WindowStatuses, which
//! combine them into the windows to mask.

use anyhow::{bail, Context, Result};
use std::io::BufRead;

use crate::chromosome::{Chr, PerChr};
use crate::constants::*;

/// One line of the GC profile: a 1,000-base window from `start` (1-based).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GcProfile {
    pub start: i32,
    pub gc_content: f64,
    pub mappable_percentage: f64,
}

impl GcProfile {
    pub fn end(&self) -> i32 {
        self.start + WINDOW_SIZE - 1
    }
    /// Doubles.greaterOrEqual(mappablePercentage, 0.85)
    pub fn is_mappable(&self) -> bool {
        doubles::greater_or_equal(self.mappable_percentage, MIN_MAPPABLE_PERCENTAGE)
    }
}

/// GCProfileFactory.loadGCContent: the lines of the human chromosomes, each chromosome's in file order.
pub fn load_gc_profile(path: &str) -> Result<PerChr<Vec<GcProfile>>> {
    let file = std::fs::File::open(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut result: PerChr<Vec<GcProfile>> = PerChr::default();
    let mut current: Option<(String, Option<Chr>)> = None;
    for line in reader.split(b'\n') {
        let line = line?;
        let line = std::str::from_utf8(&line).context("GC profile is not UTF-8")?;
        // String.split drops trailing empty strings: an empty last line has no fields
        let mut fields = line.split('\t');
        let chromosome = fields.next().unwrap_or("").trim();
        let start: i32 = field(fields.next(), line)?.parse().with_context(|| format!("java.lang.NumberFormatException in GC profile line: {line}"))?;
        let gc_content = parse_double(field(fields.next(), line)?, line)?;
        let _non_n = parse_double(field(fields.next(), line)?, line)?;
        let mappable_percentage = parse_double(field(fields.next(), line)?, line)?;
        let chr = match &current {
            Some((name, chr)) if name == chromosome => *chr,
            _ => {
                let chr = Chr::from_contig(chromosome);
                current = Some((chromosome.to_string(), chr));
                chr
            }
        };
        if let Some(chr) = chr {
            result.get_mut(chr).push(GcProfile { start: start + 1, gc_content, mappable_percentage });
        }
    }
    Ok(result)
}

fn field<'a>(f: Option<&'a str>, line: &str) -> Result<&'a str> {
    f.map(|s| s.trim()).ok_or_else(|| anyhow::anyhow!("java.lang.ArrayIndexOutOfBoundsException in GC profile line: {line}"))
}

fn parse_double(s: &str, line: &str) -> Result<f64> {
    crate::config::java_parse_double(s).ok_or_else(|| anyhow::anyhow!("java.lang.NumberFormatException in GC profile line: {line}"))
}

/// An excluded region of the resource file (chromosome as written, start, end)
#[derive(Clone, Debug)]
pub struct Region {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

/// ExcludedRegionsFile on the resource of the genome version: a header, then chromosome, start, end.
pub fn excluded_regions(v38: bool) -> Result<Vec<Region>> {
    let text = if v38 { include_str!("../resources/excluded_regions_v38.tsv") } else { include_str!("../resources/excluded_regions_v37.tsv") };
    let mut regions = Vec::new();
    for line in text.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        regions.push(Region { chromosome: f[0].to_string(), start: f[1].parse()?, end: f[2].parse()? });
    }
    Ok(regions)
}

/// SuppliedExcludedRegions.findIntersections: for each chromosome, the indices of its GC profile windows that
/// the regions overlap, searched as Java does (the regions in order, each from after the last window found).
pub fn excluded_windows(regions: &[Region], profiles: &PerChr<Vec<GcProfile>>) -> Result<PerChr<Vec<bool>>> {
    let mut by_chr: PerChr<Vec<(i32, i32)>> = PerChr::default();
    for r in regions {
        let Some(chr) = Chr::from_name(&r.chromosome) else { bail!("java.lang.IllegalArgumentException: No enum constant for {}", r.chromosome) };
        by_chr.get_mut(chr).push((r.start, r.end));
    }
    let mut result: PerChr<Vec<bool>> = PerChr::default();
    for chr in Chr::all() {
        let list = profiles.get(chr);
        let mut excluded = vec![false; list.len()];
        let mut filters = by_chr.get(chr).clone();
        // ChrBaseRegion.compareTo: start, then end
        filters.sort();
        let mut last: i64 = -1;
        for &(start, end) in &filters {
            let mut keep_searching = true;
            let mut found = false;
            let mut j = (last + 1) as usize;
            while j < list.len() && keep_searching {
                let p = &list[j];
                if start <= p.end() && end >= p.start {
                    excluded[j] = true;
                    last = j as i64;
                    found = true;
                } else {
                    keep_searching = !found;
                }
                j += 1;
            }
        }
        *result.get_mut(chr) = excluded;
    }
    Ok(result)
}

/// DiploidRegionLoader: for each chromosome, one status per 1,000 bases from position 1, false (not
/// diploid) up to each BED region and true across it, as Java builds them (a chromosome met again later in
/// the file starts again at position 1, appended to its list).
pub fn load_diploid_regions(path: &str) -> Result<PerChr<Vec<bool>>> {
    crate::info!("Reading diploid regions from {path}");
    let file = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") || path.ends_with(".bgz") {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };
    let mut regions: PerChr<Vec<bool>> = PerChr::default();
    let mut current: Option<(String, Chr)> = None;
    let mut position = 0i32;
    for line in reader.lines() {
        let line = line?;
        // BEDCodec.decode: blank, comment, track and browser lines are skipped; fields split at a tab or spaces
        if line.trim().is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let tokens: Vec<&str> = line.split(['\t', ' ']).filter(|t| !t.is_empty()).collect();
        if tokens.len() < 2 {
            continue;
        }
        let contig = tokens[0];
        let start: i32 = tokens[1].parse::<i32>().with_context(|| format!("invalid BED line: {line}"))? + 1;
        let end: i32 = if tokens.len() > 2 { tokens[2].parse().with_context(|| format!("invalid BED line: {line}"))? } else { start };
        if current.as_ref().is_none_or(|(name, _)| name != contig) {
            let Some(chr) = Chr::from_name(contig) else { bail!("java.lang.IllegalArgumentException: No enum constant for {contig}") };
            current = Some((contig.to_string(), chr));
            position = 1;
        }
        let list = regions.get_mut(current.as_ref().unwrap().1);
        let mut p = position;
        while p < start {
            list.push(false);
            p += WINDOW_SIZE;
        }
        let mut p = start;
        while p < end {
            list.push(true);
            p += WINDOW_SIZE;
        }
        position = end + 1;
    }
    Ok(regions)
}

/// WindowStatuses: per chromosome, for each GC profile window, whether it is masked (excluded region,
/// unmappable or, with diploid regions, not diploid), and its GC content for windows of low depth.
pub struct WindowStatuses {
    pub masked: PerChr<Vec<bool>>,
    pub gc_content: PerChr<Vec<f64>>,
}

impl WindowStatuses {
    pub fn new(profiles: PerChr<Vec<GcProfile>>, excluded: &PerChr<Vec<bool>>, diploid: Option<&PerChr<Vec<bool>>>) -> WindowStatuses {
        let mut masked: PerChr<Vec<bool>> = PerChr::default();
        let mut gc_content: PerChr<Vec<f64>> = PerChr::default();
        for chr in Chr::all() {
            let list = profiles.get(chr);
            let ex = excluded.get(chr);
            let statuses: Vec<bool> = list
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    let non_diploid = match diploid {
                        Some(d) => {
                            let index = (p.start / WINDOW_SIZE) as usize;
                            d.get(chr).get(index).is_none_or(|is_diploid| !is_diploid)
                        }
                        None => false,
                    };
                    ex[i] || !p.is_mappable() || non_diploid
                })
                .collect();
            *masked.get_mut(chr) = statuses;
            *gc_content.get_mut(chr) = list.iter().map(|p| p.gc_content).collect();
        }
        WindowStatuses { masked, gc_content }
    }

    /// exclude: the status of the window at `position` (index position / 1000). Java throws when the GC
    /// profile has no such window.
    pub fn exclude(&self, chr: Chr, position: i32) -> Result<bool> {
        match self.masked.get(chr).get((position / WINDOW_SIZE) as usize) {
            Some(m) => Ok(*m),
            None => bail!("java.lang.IndexOutOfBoundsException: no GC profile window for {}:{position}", chr.short_name()),
        }
    }

    /// referenceGcValueForWindow; None where Java throws (and the caller catches it)
    pub fn reference_gc(&self, chr: Chr, position: i32) -> Option<f64> {
        self.gc_content.get(chr).get((position / WINDOW_SIZE) as usize).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn profile(n: usize, gc: f64, mappable: f64) -> Vec<GcProfile> {
        (0..n).map(|i| GcProfile { start: i as i32 * 1000 + 1, gc_content: gc, mappable_percentage: mappable }).collect()
    }

    /// SuppliedExcludedRegionsTest
    #[test]
    fn excluded_regions_intersect_windows() {
        let mut profiles: PerChr<Vec<GcProfile>> = PerChr::default();
        *profiles.get_mut(Chr::new(0)) = profile(20, 0.5, 1.0);
        *profiles.get_mut(Chr::new(1)) = profile(20, 0.5, 1.0);
        let regions = vec![
            Region { chromosome: "1".into(), start: 1500, end: 3200 },
            Region { chromosome: "chr1".into(), start: 10001, end: 10001 },
            Region { chromosome: "2".into(), start: 19999, end: 25000 },
        ];
        let ex = excluded_windows(&regions, &profiles).unwrap();
        let idx = |c: usize| -> Vec<usize> { ex.get(Chr::new(c)).iter().enumerate().filter(|(_, e)| **e).map(|(i, _)| i).collect() };
        assert_eq!(idx(0), vec![1, 2, 3, 10]);
        assert_eq!(idx(1), vec![19]);
    }

    #[test]
    fn excluded_resources_load() {
        let r38 = excluded_regions(true).unwrap();
        let r37 = excluded_regions(false).unwrap();
        assert_eq!((r38.len(), r37.len()), (21, 21));
        assert!(r38.iter().all(|r| r.chromosome == "15" || r.chromosome == "16"));
    }

    #[test]
    fn mappability_and_diploid_status() {
        let mut profiles: PerChr<Vec<GcProfile>> = PerChr::default();
        let mut p = profile(5, 0.5, 1.0);
        p[1].mappable_percentage = 0.85 - 5e-11;
        p[2].mappable_percentage = 0.84;
        *profiles.get_mut(Chr::new(0)) = p;
        let excluded = excluded_windows(&[], &profiles).unwrap();
        let mut diploid: PerChr<Vec<bool>> = PerChr::default();
        *diploid.get_mut(Chr::new(0)) = vec![true, true, true, false];
        let s = WindowStatuses::new(profiles.clone(), &excluded, None);
        assert_eq!(s.masked.get(Chr::new(0)), &vec![false, false, true, false, false]);
        let s = WindowStatuses::new(profiles, &excluded, Some(&diploid));
        assert_eq!(s.masked.get(Chr::new(0)), &vec![false, false, true, true, true]);
    }
}
