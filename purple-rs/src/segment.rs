//! The segment package: SV and PCF breakpoints clustered per chromosome (ClusterFactory), the support
//! segments between them with the centromere split (PurpleSupportSegmentFactory), the excluded-region
//! refining (SegmentRefiner, ExcludedRegionsRefiner) and Segmentation.createObservedRegions.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use hmf_common::segmentation::pcf::{PcfPosition, PcfSource};

use crate::inputs::{AmberData, CobaltData, CobaltRatio, WINDOW_SIZE};
use crate::region::{ObservedRegion, ObservedRegionFactory};
use crate::types::*;

/// SvPosition: a breakend's copy-number position
#[derive(Clone, Debug, PartialEq)]
pub struct SvPosition {
    pub chromosome: String,
    pub position: i32,
    pub sv_type: SvType,
}

#[derive(Clone, Debug)]
pub enum Pos {
    Sv(SvPosition),
    Pcf(PcfPosition),
}

impl Pos {
    fn chromosome(&self) -> &str {
        match self {
            Pos::Sv(s) => &s.chromosome,
            Pos::Pcf(p) => &p.chromosome,
        }
    }
    fn position(&self) -> i32 {
        match self {
            Pos::Sv(s) => s.position,
            Pos::Pcf(p) => p.position,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Cluster {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub pcf_positions: Vec<PcfPosition>,
    pub variants: Vec<SvPosition>,
}

impl Cluster {
    fn has_ratios(&self) -> bool {
        self.pcf_positions.iter().any(|p| p.source != PcfSource::TumorBaf)
    }
}

/// ClusterFactory.findFirstValidCobaltRatio
fn find_first_valid_cobalt_ratio(position: i32, index: usize, ratios: &[CobaltRatio]) -> i32 {
    let window = Window { size: WINDOW_SIZE };
    let min = window.start(position) - WINDOW_SIZE + 1;
    if !ratios.is_empty() {
        let mut i = index as i64;
        while i >= 0 {
            let r = &ratios[i as usize];
            if r.position <= min && greater_than(r.tumor_gc_ratio, -1.0) {
                return r.position + 1;
            }
            i -= 1;
        }
    }
    min
}

/// ClusterFactory.buildChromosomeClusters
pub fn build_chromosome_clusters(sv_positions: &[SvPosition], pcf_positions: &[PcfPosition], cobalt_ratios: &[CobaltRatio]) -> Vec<Cluster> {
    let mut all: Vec<Pos> = Vec::with_capacity(sv_positions.len() + pcf_positions.len());
    all.extend(sv_positions.iter().cloned().map(Pos::Sv));
    all.extend(pcf_positions.iter().cloned().map(Pos::Pcf));
    all.sort_by(|a, b| position_compare(a.chromosome(), a.position(), b.chromosome(), b.position()));

    let mut clusters: Vec<Cluster> = Vec::new();
    let mut cobalt_index: usize = 0;
    let mut have_segment = false;
    let mut last_position_was_amber_pcf_end = false;
    let n = all.len();
    for i in 0..n {
        let position = &all[i];
        if position.position() == 1 {
            continue;
        }
        while cobalt_index + 1 < cobalt_ratios.len() && cobalt_ratios[cobalt_index].position < position.position() {
            cobalt_index += 1;
        }
        let mut can_segment = true;
        let pcf = match position {
            Pos::Pcf(p) => Some(p),
            _ => None,
        };
        let earliest = find_first_valid_cobalt_ratio(position.position(), cobalt_index, cobalt_ratios);
        if have_segment {
            let segment_end = clusters.last().unwrap().end;
            if earliest > segment_end {
                if let Some(p) = pcf {
                    if p.source == PcfSource::TumorBaf {
                        let use_amber = last_position_was_amber_pcf_end && p.is_segment_start() && i + 2 < n;
                        if !use_amber {
                            continue;
                        }
                    }
                }
            } else {
                can_segment = false;
            }
        } else {
            can_segment = pcf.is_none_or(|p| p.source != PcfSource::TumorBaf);
        }
        if !have_segment || can_segment {
            clusters.push(Cluster {
                chromosome: position.chromosome().to_string(),
                start: earliest,
                end: 0,
                pcf_positions: Vec::new(),
                variants: Vec::new(),
            });
            have_segment = true;
        }
        let segment = clusters.last_mut().unwrap();
        segment.end = position.position();
        last_position_was_amber_pcf_end = false;
        match position {
            Pos::Sv(s) => segment.variants.push(s.clone()),
            Pos::Pcf(p) => {
                segment.pcf_positions.push(p.clone());
                if p.source == PcfSource::TumorBaf && p.is_segment_end() {
                    last_position_was_amber_pcf_end = true;
                }
            }
        }
    }
    clusters
}

#[derive(Clone, Debug, PartialEq)]
pub struct SupportSegment {
    pub chromosome: String,
    pub ratio_support: bool,
    pub support: SegmentSupport,
    pub sv_cluster: bool,
    pub start: i32,
    pub end: i32,
    pub min_start: i32,
    pub max_start: i32,
}

impl SupportSegment {
    pub fn new(chromosome: &str, start: i32, end: i32, ratio_support: bool, support: SegmentSupport, sv_cluster: bool, min_start: i32, max_start: i32) -> SupportSegment {
        SupportSegment { chromosome: chromosome.to_string(), ratio_support, support, sv_cluster, start, end, min_start, max_start }
    }

    /// PurpleSupportSegment.split
    pub fn split(&self, position: i32) -> Vec<SupportSegment> {
        if position <= self.start || position > self.end {
            return vec![self.clone()];
        }
        let end_of_left = position - 1;
        let new_max_start = self.max_start.min(end_of_left);
        vec![
            SupportSegment::new(&self.chromosome, self.start, end_of_left, self.ratio_support, self.support, self.sv_cluster, self.min_start, new_max_start),
            SupportSegment::new(&self.chromosome, position, self.end, self.ratio_support, SegmentSupport::Excl, self.sv_cluster, position, position),
        ]
    }

    /// PurpleSupportSegment.splitBy
    pub fn split_by(&self, chromosome: &str, start: i32, end: i32) -> Vec<SupportSegment> {
        if self.chromosome != chromosome {
            return vec![self.clone()];
        }
        self.split(start).iter().flat_map(|s| s.split(end + 1)).collect()
    }

    fn contains(&self, chromosome: &str, position: i32) -> bool {
        self.chromosome == chromosome && self.start <= position && self.end >= position
    }
}

fn create_pcf_segment(chromosome: &str, start: i32, pcf_positions: &[PcfPosition]) -> SupportSegment {
    let tumor: Vec<&PcfPosition> = pcf_positions.iter().filter(|p| p.source == PcfSource::TumorRatio).collect();
    let min_start = tumor.iter().map(|p| p.min_position).min().unwrap_or(start);
    let max_start = tumor.iter().map(|p| p.max_position).max().unwrap_or(start);
    SupportSegment::new(chromosome, start, 0, true, SegmentSupport::None, false, min_start, max_start)
}

/// PurpleSupportSegmentFactory.createClusterSegments
fn create_cluster_segments(chromosome: &str, length: i32, clusters: &[Cluster]) -> Vec<SupportSegment> {
    let mut result = Vec::new();
    let mut segment = SupportSegment::new(chromosome, 1, 0, true, SegmentSupport::Telomere, false, 1, 1);
    for cluster in clusters {
        let ratio_support = cluster.has_ratios();
        if !cluster.variants.is_empty() {
            for variant in &cluster.variants {
                if variant.position != segment.start {
                    segment.end = variant.position - 1;
                    result.push(segment);
                    segment = SupportSegment::new(
                        &cluster.chromosome,
                        variant.position,
                        0,
                        ratio_support,
                        SegmentSupport::from_variant(variant.sv_type),
                        true,
                        variant.position,
                        variant.position,
                    );
                } else {
                    segment.support = SegmentSupport::Multiple;
                }
            }
            segment.sv_cluster = false;
        } else {
            let first = &cluster.pcf_positions[0];
            segment.end = first.position - 1;
            result.push(segment);
            segment = create_pcf_segment(&first.chromosome, first.position, &cluster.pcf_positions);
        }
    }
    segment.end = length;
    result.push(segment);
    result
}

/// PurpleSupportSegmentFactory.addCentromere
fn add_centromere(centromere: i32, chromosome: &str, segments: Vec<SupportSegment>) -> Vec<SupportSegment> {
    let mut result = Vec::with_capacity(segments.len() + 1);
    for segment in segments {
        if segment.contains(chromosome, centromere) {
            if segment.start == centromere {
                let mut start = segment.clone();
                start.support = SegmentSupport::Centromere;
                result.push(start);
            } else {
                let mut start = segment.clone();
                start.end = centromere - 1;
                start.max_start = start.max_start.min(start.end);
                let mut end = segment.clone();
                end.support = SegmentSupport::Centromere;
                end.start = centromere;
                end.min_start = centromere;
                end.max_start = centromere;
                result.push(start);
                result.push(end);
            }
        } else {
            result.push(segment);
        }
    }
    result
}

/// PurpleSupportSegmentFactory.createChromosomeSegments
pub fn create_chromosome_segments(chromosome: &str, centromere: i32, length: i32, clusters: &[Cluster]) -> Vec<SupportSegment> {
    let mut segments = create_cluster_segments(chromosome, length, clusters);
    for i in 1..segments.len() {
        let prev_end = segments[i - 1].end;
        segments[i].min_start = segments[i].min_start.max(prev_end + 1);
    }
    add_centromere(centromere, chromosome, segments)
}

/// GenomePositions.union of two sorted lists: equal positions are kept once (the first list's)
fn union_positions(first: &[PcfPosition], second: &[PcfPosition]) -> Vec<PcfPosition> {
    let cmp = |a: &PcfPosition, b: &PcfPosition| position_compare(&a.chromosome, a.position, &b.chromosome, b.position);
    let mut out = Vec::with_capacity(first.len() + second.len());
    let (mut i, mut j) = (0, 0);
    while i < first.len() || j < second.len() {
        if i >= first.len() || (j < second.len() && cmp(&second[j], &first[i]) == Ordering::Less) {
            out.push(second[j].clone());
            j += 1;
        } else if j >= second.len() || cmp(&first[i], &second[j]) == Ordering::Less {
            out.push(first[i].clone());
            i += 1;
        } else {
            out.push(first[i].clone());
            i += 1;
            j += 1;
        }
    }
    out
}

fn union_of_maps(a: &BTreeMap<Chromosome, Vec<PcfPosition>>, b: &BTreeMap<Chromosome, Vec<PcfPosition>>) -> BTreeMap<Chromosome, Vec<PcfPosition>> {
    let mut out = BTreeMap::new();
    for c in a.keys().chain(b.keys()) {
        if out.contains_key(c) {
            continue;
        }
        let v = match (a.get(c), b.get(c)) {
            (Some(x), Some(y)) => union_positions(x, y),
            (Some(x), None) => x.clone(),
            (None, Some(y)) => y.clone(),
            (None, None) => unreachable!(),
        };
        out.insert(*c, v);
    }
    out
}

/// PCFPositionsSupplier.createPositions
pub fn pcf_positions(amber: &AmberData, cobalt: &CobaltData) -> BTreeMap<Chromosome, Vec<PcfPosition>> {
    union_of_maps(&union_of_maps(&cobalt.tumor_segments, &cobalt.reference_segments), &amber.tumor_segments)
}

/// SVSegmentFactory.create plus Multimaps.fromPositions: the SV breakends (not INS) by chromosome, sorted
pub fn sv_positions(variants: &[SvPosition]) -> BTreeMap<Chromosome, Vec<SvPosition>> {
    let mut positions: Vec<SvPosition> = variants.iter().filter(|v| v.sv_type != SvType::Ins).cloned().collect();
    positions.sort_by(|a, b| position_compare(&a.chromosome, a.position, &b.chromosome, b.position));
    let mut map: BTreeMap<Chromosome, Vec<SvPosition>> = BTreeMap::new();
    for p in positions {
        if Chromosome::contains(&p.chromosome) {
            map.entry(chr(&p.chromosome)).or_default().push(p);
        }
    }
    map
}

pub struct ReferenceCoords {
    pub v38: bool,
}

impl ReferenceCoords {
    pub fn length(&self, c: Chromosome) -> i32 {
        hmf_common::ref_coords::chromosome_length(c.short_name(), !self.v38).unwrap_or(0)
    }
    pub fn centromere(&self, c: Chromosome) -> i32 {
        hmf_common::segmentation::centromere(c.short_name(), self.v38).unwrap_or(0)
    }
    pub fn chromosome_name(&self, c: Chromosome) -> String {
        c.versioned(self.v38)
    }
}

/// PurpleSupportSegmentFactory.createSegments: `sv_breakends` are the start and end legs of every SV
pub fn create_segments(
    coords: &ReferenceCoords,
    sv_breakends: &[SvPosition],
    pcf: &BTreeMap<Chromosome, Vec<PcfPosition>>,
    ratios: &BTreeMap<Chromosome, Vec<CobaltRatio>>,
) -> Vec<SupportSegment> {
    let svs = sv_positions(sv_breakends);
    let empty_r: Vec<CobaltRatio> = Vec::new();
    let empty_s: Vec<SvPosition> = Vec::new();
    let mut all = Vec::new();
    for (c, positions) in pcf {
        let clusters = build_chromosome_clusters(svs.get(c).unwrap_or(&empty_s), positions, ratios.get(c).unwrap_or(&empty_r));
        // ClusterFactory keys each chromosome's clusters in a multimap: a chromosome without clusters has none
        if clusters.is_empty() {
            continue;
        }
        let name = coords.chromosome_name(*c);
        all.extend(create_chromosome_segments(&name, coords.centromere(*c), coords.length(*c), &clusters));
    }
    all.sort_by(|a, b| region_compare(&a.chromosome, a.start, &b.chromosome, b.start));
    all
}

/// PurpleSupportSegmentFactory.validateSegments
pub fn validate_segments(segments: &[SupportSegment]) -> bool {
    let mut valid = true;
    for i in 1..segments.len() {
        let s = &segments[i];
        if !positions_within(s.min_start, s.max_start, s.start, s.end) {
            crate::error!("purple segment({}:{}-{}) has invalid min/maxStart({}-{})", s.chromosome, s.start, s.end, s.min_start, s.max_start);
            valid = false;
        }
        let prev = &segments[i - 1];
        if s.chromosome != prev.chromosome {
            continue;
        }
        if s.start <= prev.end {
            crate::error!("purple segment({}:{}-{}) overlaps previous({}-{})", s.chromosome, s.start, s.end, prev.start, prev.end);
            valid = false;
        }
    }
    valid
}

/// ExcludedRegionsRefiner over the excluded (immune) regions: segments overlapping a region are split at
/// its start and after its end. The regions of a chromosome are kept as a TreeMultimap value set, so
/// two with the same start are one.
pub fn refine_excluded(segments: Vec<SupportSegment>, excluded: &[ChrRegion]) -> Vec<SupportSegment> {
    let mut by_chr: BTreeMap<String, Vec<(i32, i32)>> = BTreeMap::new();
    for r in excluded {
        let list = by_chr.entry(r.chromosome.clone()).or_default();
        if !list.iter().any(|(s, _)| *s == r.start) {
            list.push((r.start, r.end));
        }
    }
    for list in by_chr.values_mut() {
        list.sort_by_key(|(s, _)| *s);
    }
    let mut result = Vec::with_capacity(segments.len() + 16);
    let mut i = 0;
    while i < segments.len() {
        let chromosome = segments[i].chromosome.clone();
        let mut j = i;
        while j < segments.len() && segments[j].chromosome == chromosome {
            j += 1;
        }
        let mut current: Vec<SupportSegment> = segments[i..j].to_vec();
        if let Some(regions) = by_chr.get(&chromosome) {
            for (start, end) in regions {
                current = current.iter().flat_map(|s| s.split_by(&chromosome, *start, *end)).collect();
            }
        }
        result.extend(current);
        i = j;
    }
    result
}

/// Segmentation.createObservedRegions
pub fn create_observed_regions(coords: &ReferenceCoords, sv_breakends: &[SvPosition], amber: &AmberData, cobalt: &CobaltData, factory: &ObservedRegionFactory) -> Vec<ObservedRegion> {
    let pcf = pcf_positions(amber, cobalt);
    let segments = create_segments(coords, sv_breakends, &pcf, &cobalt.ratios);
    if !validate_segments(&segments) {
        return Vec::new();
    }
    factory.form_observed_regions(segments, &amber.chromosome_bafs, &cobalt.ratios)
}

/// Segmentation.validateObservedRegions
pub fn validate_observed_regions(regions: &[ObservedRegion]) -> bool {
    let mut valid = true;
    for i in 1..regions.len() {
        let r = &regions[i];
        if r.start > r.end {
            crate::error!("observed region({}:{}-{}) is invalid", r.chromosome, r.start, r.end);
            valid = false;
        }
        let prev = &regions[i - 1];
        if r.chromosome != prev.chromosome {
            continue;
        }
        if r.start <= prev.end {
            crate::error!("observed region({}:{}-{}) overlaps previous({}-{})", r.chromosome, r.start, r.end, prev.start, prev.end);
            valid = false;
        }
    }
    valid
}
