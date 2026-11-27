//! TelbamReader: the read groups of a telbam, and the excluded regions BreakEndParams loads.

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;

use crate::prep::FxBuild;
use crate::read::{Contigs, Read, ReadGroup};

/// A region of the excluded (or included) list, 1-based inclusive
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GenomeRegion {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

/// The read groups of a telbam. Java keeps them in a HashMap; nothing TEAL writes depends on their
/// order, so they are kept in the order their names first appear.
pub struct TelbamGroups {
    pub groups: Vec<Option<ReadGroup>>,
}

impl TelbamGroups {
    pub fn iter(&self) -> impl Iterator<Item = &ReadGroup> {
        self.groups.iter().flatten()
    }
}

/// ContigComparator.compare == 0: the same chromosome rank (all non-human contigs share one)
fn same_rank(a: &str, b: &str) -> bool {
    hmf_common::regions::chromosome_rank(a) == hmf_common::regions::chromosome_rank(b)
}

fn should_exclude(rg: &ReadGroup, excluded: &[GenomeRegion], included: Option<&[GenomeRegion]>) -> bool {
    let overlaps = |r: &Read, g: &GenomeRegion| r.start <= g.end && r.end() >= g.start && same_rank(&r.chrom, &g.chromosome);
    match included {
        // none included: excluded
        Some(inc) => !rg.all_reads().any(|r| inc.iter().any(|g| overlaps(r, g))),
        // all excluded: excluded
        None => rg.all_reads().all(|r| excluded.iter().any(|g| overlaps(r, g))),
    }
}

/// TelbamReader.read
pub fn read_telbam(path: &str, excluded: &[GenomeRegion], included: Option<&[GenomeRegion]>) -> Result<TelbamGroups> {
    crate::info!("processing telbam: {path}");
    let mut reader = bam::Reader::from_path(path).with_context(|| format!("opening {path}"))?;
    let contigs = Contigs::from_header(reader.header());
    let mut index: HashMap<Box<[u8]>, usize, FxBuild> = HashMap::default();
    let mut groups: Vec<Option<ReadGroup>> = Vec::new();
    let mut rec = bam::Record::new();
    while let Some(r) = reader.read(&mut rec) {
        r.with_context(|| format!("reading {path}"))?;
        let read = Read::from_record(&rec, &contigs, true, false);
        let gi = match index.get(rec.qname()) {
            Some(&i) => i,
            None => {
                index.insert(Box::from(rec.qname()), groups.len());
                groups.push(Some(ReadGroup::new(rec.qname())));
                groups.len() - 1
            }
        };
        let rg = groups[gi].as_mut().unwrap();
        rg.accept_read(read);
        if rg.is_complete()? && should_exclude(rg, excluded, included) {
            groups[gi] = None;
            index.remove(rec.qname());
        }
    }
    Ok(TelbamGroups { groups })
}

/// BreakEndParams.loadExcludedRegionBed: the regions of blacklistedTelomereRegions.37.bed for V37, none for
/// V38. Java copies the resource to a temporary file without the .bed extension, so the BED's 0-based
/// starts are taken as 1-based (TE7); corrected, they are converted.
pub fn excluded_regions(v37: bool) -> Vec<GenomeRegion> {
    if !v37 {
        return Vec::new();
    }
    let bed_start_offset = if crate::compat::java("TE7") { 0 } else { 1 };
    let text = include_str!("../resources/blacklistedTelomereRegions.37.bed");
    // loadChrBaseRegions: per chromosome, sorted and merged where they overlap
    let mut by_chr: Vec<(String, Vec<(i32, i32)>)> = Vec::new();
    for line in text.lines().filter(|l| !l.is_empty()) {
        let v: Vec<&str> = line.split('\t').collect();
        let start = v[1].parse::<i32>().unwrap() + bed_start_offset;
        let end = v[2].parse::<i32>().unwrap();
        match by_chr.iter_mut().find(|(c, _)| c == v[0]) {
            Some((_, l)) => l.push((start, end)),
            None => by_chr.push((v[0].to_string(), vec![(start, end)])),
        }
    }
    let mut out = Vec::new();
    for (chr, mut l) in by_chr {
        l.sort();
        let mut merged: Vec<(i32, i32)> = Vec::new();
        for (s, e) in l {
            match merged.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        out.extend(merged.into_iter().map(|(start, end)| GenomeRegion { chromosome: chr.clone(), start, end }));
    }
    out
}
