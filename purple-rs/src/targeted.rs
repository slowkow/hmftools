//! Targeted mode: TargetRegionDataBuilder and purple.target_region_cn.tsv (each COBALT window over a
//! target region with its copy number, split where a copy number starts inside the window).

use anyhow::Result;
use std::collections::HashMap;

use hmf_common::fmt::decimal_format_min;

use crate::copynumber::PurpleCopyNumber;
use crate::inputs::CobaltRatio;
use crate::reference::TargetRegionsData;
use crate::region::ObservedRegion;
use crate::types::*;

pub struct TargetRegionCopyNumber<'a> {
    pub ratio: CobaltRatio,
    pub regions: Vec<&'a (i32, i32, String)>,
    pub copy_number: &'a PurpleCopyNumber,
    pub germline_status: GermlineStatus,
}

fn window(r: &CobaltRatio) -> (i32, i32) {
    (r.position, r.position + 1000 - 1)
}

fn overlaps(w: (i32, i32), start: i32, end: i32) -> bool {
    positions_overlap(w.0, w.1, start, end)
}

/// the find*Index helpers: the next element from `current` overlapping the window, or None
fn find_index<T>(items: &[T], current: usize, w: (i32, i32), bounds: impl Fn(&T) -> (i32, i32)) -> Option<usize> {
    if bounds(&items[current]).0 > w.1 {
        return None;
    }
    let mut i = current + 1;
    while i < items.len() {
        let (s, e) = bounds(&items[i]);
        if overlaps(w, s, e) {
            return Some(i);
        }
        if s > w.1 {
            break;
        }
        i += 1;
    }
    None
}

/// TargetRegionDataBuilder.buildTargetRegionData
pub fn build<'a>(
    v38: bool,
    targets: &'a TargetRegionsData,
    segments: &'a [ObservedRegion],
    ratios: &'a std::collections::BTreeMap<Chromosome, Vec<CobaltRatio>>,
    copy_numbers: &'a [PurpleCopyNumber],
) -> Vec<TargetRegionCopyNumber<'a>> {
    let mut cn_map: HashMap<&str, Vec<&PurpleCopyNumber>> = HashMap::new();
    for c in copy_numbers {
        cn_map.entry(c.chromosome.as_str()).or_default().push(c);
    }
    let mut seg_map: HashMap<&str, Vec<&ObservedRegion>> = HashMap::new();
    for s in segments {
        seg_map.entry(s.chromosome.as_str()).or_default().push(s);
    }
    let mut out = Vec::new();
    for c in Chromosome::ALL.iter().copied() {
        let chr_str = if v38 { format!("chr{}", c.short_name()) } else { c.short_name().to_string() };
        let (Some(regions), Some(cobalt)) = (targets.tagged.get(&chr_str), ratios.get(&c)) else { continue };
        if regions.is_empty() || cobalt.is_empty() {
            continue;
        }
        let (Some(cns), Some(segs)) = (cn_map.get(chr_str.as_str()), seg_map.get(chr_str.as_str())) else { continue };
        if cns.is_empty() || segs.is_empty() {
            continue;
        }
        let (mut ti, mut ci, mut si) = (0usize, 0usize, 0usize);
        let mut current_segment: Option<&ObservedRegion> = Some(segs[0]);
        // the current target region and copy number are never null once set (Java starts at the first)
        let mut target_ok = true;
        let mut cn_ok = true;
        for ratio in cobalt {
            let w = window(ratio);
            if !target_ok || !overlaps(w, regions[ti].0, regions[ti].1) {
                match find_index(regions, ti, w, |r| (r.0, r.1)) {
                    None => continue,
                    Some(i) => {
                        ti = i;
                        target_ok = true;
                    }
                }
            }
            if !cn_ok || !overlaps(w, cns[ci].start, cns[ci].end) {
                match find_index(cns, ci, w, |c| (c.start, c.end)) {
                    None => continue,
                    Some(i) => {
                        ci = i;
                        cn_ok = true;
                    }
                }
            }
            let segment_overlaps = current_segment.is_some_and(|s| overlaps(w, s.start, s.end));
            if !segment_overlaps {
                match find_index(segs, si, w, |s| (s.start, s.end)) {
                    None => current_segment = None,
                    Some(i) => {
                        si = i;
                        current_segment = Some(segs[si]);
                    }
                }
            }
            // findOverlappingTargetRegions: the current region, then those after, then those before
            let mut overlapping = vec![&regions[ti]];
            let mut i = ti + 1;
            while i < regions.len() && overlaps(w, regions[i].0, regions[i].1) {
                overlapping.push(&regions[i]);
                i += 1;
            }
            let mut i = ti as i64 - 1;
            while i >= 0 && overlaps(w, regions[i as usize].0, regions[i as usize].1) {
                overlapping.push(&regions[i as usize]);
                i -= 1;
            }
            // findInternalCopyNumbers
            let mut internal = vec![cns[ci]];
            for c in cns.iter().skip(ci + 1) {
                if overlaps(w, c.start, c.end) {
                    internal.push(c);
                } else {
                    break;
                }
            }
            let status = current_segment.map(|s| s.germline_status).unwrap_or(GermlineStatus::Unknown);
            for (k, c) in internal.iter().enumerate() {
                let mut r = ratio.clone();
                if k > 0 {
                    r.position = c.start;
                }
                out.push(TargetRegionCopyNumber { ratio: r, regions: overlapping.clone(), copy_number: c, germline_status: status });
            }
        }
    }
    out
}

/// TargetRegionsCopyNumberFile.write
pub fn write(path: &str, data: &[TargetRegionCopyNumber]) -> Result<()> {
    let f = |v: f64| decimal_format_min(v, 4, 4);
    let mut text = String::from("chromosome\twindowStart\twindowEnd\tbedRegions\tmasked\taverageDepth\twindowGCContent\twindowTumorRatio\tregionStart\tregionEnd\tcopyNumber\tminorAlleleCopyNumber\tdepthWindowCount\tbafCount\tgermlineStatus\tcopyNumberMethod\n");
    for d in data {
        let (ws, we) = window(&d.ratio);
        let regions: Vec<String> = d.regions.iter().map(|(s, e, t)| if t.is_empty() { format!("{s}-{e}") } else { format!("{t}:{s}-{e}") }).collect();
        let c = d.copy_number;
        text.push_str(&format!(
            "{}\t{ws}\t{we}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            d.ratio.chromosome,
            regions.join(";"),
            d.ratio.tumor_gc_ratio < 0.0,
            f(d.ratio.tumor_read_depth),
            f(d.ratio.tumor_gc_content),
            f(d.ratio.tumor_gc_ratio),
            c.start,
            c.end,
            f(c.average_tumor_copy_number),
            f(c.minor_allele_copy_number()),
            c.depth_window_count,
            c.baf_count,
            d.germline_status.name(),
            c.method.name()
        ));
    }
    std::fs::write(path, text)?;
    Ok(())
}
