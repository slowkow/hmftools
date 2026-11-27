//! Packages `visualiser.data` and the data parts of `visualiser.circos`: segment tracks and frames,
//! connectors, copy number limits, highlights, position scaling (ScalePosition, ScaleContig), downsampling,
//! colours and line thickness, genes and exons.

use std::collections::{BTreeMap, HashMap};

use hmf_common::java_hash::{JavaHashMap, JavaHashSet};

use super::files::*;
use super::genome::*;

// ---------------------------------------------------------------------------------------------------------
// RefGenomeCoordinates

#[derive(Clone, Copy, Debug)]
pub struct RefCoords {
    pub v37: bool,
}

impl RefCoords {
    /// length(chromosome): 0 for a non-human chromosome
    pub fn length(&self, chromosome: &str) -> i32 {
        human_chromosome_name(chromosome).and_then(|n| hmf_common::ref_coords::chromosome_length(n, self.v37)).unwrap_or(0)
    }
    pub fn centromere(&self, chromosome: &str) -> i32 {
        human_chromosome_name(chromosome).and_then(|n| hmf_common::ref_coords::centromere(n, self.v37)).unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------------------------------------
// VisLinks

pub fn find_link<'a>(chr: &str, pos: i32, links: &'a [VisSvData]) -> Option<&'a VisSvData> {
    links.iter().find(|l| (l.chr_start == chr && l.pos_start == pos) || (l.chr_end == chr && l.pos_end == pos))
}

pub fn find_start_link<'a>(chr: &str, pos: i32, links: &'a [VisSvData]) -> Option<&'a VisSvData> {
    links.iter().find(|l| l.chr_start == chr && l.pos_start == pos)
}

pub fn find_end_link<'a>(chr: &str, pos: i32, links: &'a [VisSvData]) -> Option<&'a VisSvData> {
    links.iter().find(|l| l.chr_end == chr && l.pos_end == pos)
}

/// VisLinks.allPositions
pub fn link_positions(links: &[VisSvData]) -> Vec<GenomePosition> {
    let mut out = Vec::new();
    for l in links {
        if l.is_valid_start() && l.pos_start != -1 {
            out.push(GenomePosition::new(&l.chr_start, l.pos_start));
        }
        if l.is_valid_end() && l.pos_end != -1 {
            out.push(GenomePosition::new(&l.chr_end, l.pos_end));
        }
    }
    sort_positions(&mut out);
    out
}

fn connected(link: &VisSvData, s: &VisSegment) -> bool {
    let at_start = s.chromosome == link.chr_start && (s.start == link.pos_start || s.end == link.pos_start);
    let at_end = s.chromosome == link.chr_end && (s.start == link.pos_end || s.end == link.pos_end);
    at_start || at_end
}

/// VisLinks.addFrame: each link's frame is one more than the lowest frame of the segments it touches
pub fn add_frame(segments: &[VisSegment], links: &[VisSvData]) -> Vec<VisSvData> {
    links
        .iter()
        .map(|l| {
            let min = segments.iter().filter(|s| connected(l, s)).map(|s| s.frame).min().unwrap_or(0);
            let mut n = l.clone();
            n.frame = min + 1;
            n
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------------------
// VisSegments

pub fn segment_ploidy_before(track: i32, chr: &str, pos: i32, segments: &[VisSegment]) -> f64 {
    let mut sum = 0.0;
    for s in segments.iter().filter(|s| s.track < track) {
        if s.chromosome == chr && (s.start == pos || s.end == pos) {
            sum += s.link_ploidy;
        }
    }
    sum
}

/// VisSegments.entireChromosome
pub fn entire_chromosome(chromosome: &str, coords: &RefCoords) -> VisSegment {
    let len = coords.length(chromosome).to_string();
    let mut s = VisSegment::new(-1, -1, chromosome, "1", &len, 0.0, false).expect("numeric positions");
    s.terminal_start = SegmentTerminal::TELOMERE;
    s.terminal_end = SegmentTerminal::TELOMERE;
    s
}

fn centromere_segment(chromosome: &str, coords: &RefCoords) -> VisSegment {
    let pos = coords.centromere(chromosome).to_string();
    let mut s = VisSegment::new(-1, -1, chromosome, &pos, &pos, 0.0, false).expect("numeric positions");
    s.terminal_start = SegmentTerminal::CENTROMERE;
    s.terminal_end = SegmentTerminal::CENTROMERE;
    s
}

/// VisSegments.extendTerminals
pub fn extend_terminals(
    terminal_distance: i32, segments: &[VisSegment], links: &[VisSvData], all_positions: &[GenomePosition], show_simple_sv_segments: bool,
    coords: &RefCoords,
) -> anyhow::Result<Vec<VisSegment>> {
    let mm = min_max_per_chromosome(all_positions);
    let mut results = Vec::new();
    for s in segments {
        let mut seg = s.clone();
        let centromere = coords.centromere(&seg.chromosome);
        let length = coords.length(&seg.chromosome);
        if seg.terminal_start != SegmentTerminal::NONE {
            let min = mm.get(&seg.chromosome).map(|m| m.0).ok_or_else(|| anyhow::anyhow!("NullPointerException: no positions on {}", seg.chromosome))?;
            seg.start = if seg.terminal_start == SegmentTerminal::CENTROMERE && min < centromere {
                centromere
            } else {
                i32::max(1, min - terminal_distance)
            };
        }
        if seg.terminal_end != SegmentTerminal::NONE {
            let max = mm.get(&seg.chromosome).map(|m| m.1).ok_or_else(|| anyhow::anyhow!("NullPointerException: no positions on {}", seg.chromosome))?;
            seg.end = if seg.terminal_end == SegmentTerminal::CENTROMERE && max > centromere {
                centromere
            } else {
                i32::min(length, max + terminal_distance)
            };
        }
        results.push(seg);
    }
    Ok(increment_on_chromosome(add_centromeres(results, coords), links, show_simple_sv_segments))
}

/// VisSegments.addCentromeres: a centromere segment for each chromosome whose segments span its centromere
/// and that has none yet (in the iteration order of Java's HashSet of chromosome names)
pub fn add_centromeres(segments: Vec<VisSegment>, coords: &RefCoords) -> Vec<VisSegment> {
    if segments.is_empty() {
        return segments;
    }
    let existing: Vec<&str> = segments
        .iter()
        .filter(|s| s.terminal_start == SegmentTerminal::CENTROMERE || s.terminal_end == SegmentTerminal::CENTROMERE)
        .map(|s| s.chromosome.as_str())
        .collect();
    let mut required: JavaHashSet<String> = JavaHashSet::new();
    for r in span_regions(&segments) {
        let c = coords.centromere(&r.chromosome);
        if r.start < c && r.end > c {
            required.insert(r.chromosome.clone());
        }
    }
    required.retain(|c| !existing.contains(&c.as_str()));
    let mut result = segments.clone();
    for chr in required.iter() {
        result.push(centromere_segment(chr, coords));
    }
    result
}

/// VisSegments.incrementOnChromosome: tracks per chromosome and frames for the shown segments
pub fn increment_on_chromosome(segments: Vec<VisSegment>, links: &[VisSvData], show_simple_sv_segments: bool) -> Vec<VisSegment> {
    let mut track_map: HashMap<String, i32> = HashMap::new();
    let mut result = Vec::new();
    let mut frame = -1;
    for s in segments {
        if s.cluster_id == -1 {
            let mut n = s.clone();
            n.track = 0;
            result.push(n);
        } else if show_segment(show_simple_sv_segments, &s, links) {
            let current = match track_map.get(&s.chromosome) {
                None => 1,
                Some(t) => t + 1,
            };
            track_map.insert(s.chromosome.clone(), current);
            let mut n = s.clone();
            n.track = current;
            frame += 2;
            n.frame = frame;
            result.push(n);
        }
    }
    result
}

fn show_segment(show_simple: bool, s: &VisSegment, links: &[VisSvData]) -> bool {
    if s.terminal_start == SegmentTerminal::NONE && s.terminal_end == SegmentTerminal::NONE {
        return true;
    }
    let shown = |pos: i32| find_link(&s.chromosome, pos, links).is_some_and(|l| !l.connectors_only(show_simple));
    if s.terminal_start == SegmentTerminal::NONE && shown(s.start) {
        return true;
    }
    s.terminal_end == SegmentTerminal::NONE && shown(s.end)
}

// ---------------------------------------------------------------------------------------------------------
// VisCopyNumbers

/// VisCopyNumbers.copyNumbers: the copy number regions within the span (plus 10% up to 1000 bases each side),
/// a base wider where the next or previous region has a higher copy number, clipped to the span
pub fn limit_copy_numbers(copy_numbers: &[VisCopyNumber], span: &[Region]) -> Vec<VisCopyNumber> {
    let mut results = Vec::new();
    for (i, cn) in copy_numbers.iter().enumerate() {
        let segs: Vec<&Region> = span.iter().filter(|r| r.chromosome == cn.chromosome).collect();
        if segs.is_empty() {
            continue;
        }
        let mut min_pos = segs.iter().map(|r| r.start).min().unwrap_or(0);
        let mut max_pos = segs.iter().map(|r| r.end).max().unwrap_or(0);
        let distance = max_pos - min_pos;
        let additional = i64::max(1, i64::min(1000, java_round(0.1 * distance as f64))) as i32;
        min_pos -= additional;
        max_pos += additional;
        if cn.end >= min_pos && cn.start <= max_pos {
            let start_decreasing = i > 0 && cn_less_than(cn, &copy_numbers[i - 1]);
            let start = if start_decreasing { cn.start - 1 } else { cn.start };
            let end_increasing = i + 1 < copy_numbers.len() && cn_less_than(cn, &copy_numbers[i + 1]);
            let end = if end_increasing { cn.end + 1 } else { cn.end };
            let mut n = cn.clone();
            n.start = i32::max(min_pos, start);
            n.end = i32::min(max_pos, end);
            n.truncated = min_pos > start || max_pos < end;
            results.push(n);
        }
    }
    results
}

fn cn_less_than(a: &VisCopyNumber, b: &VisCopyNumber) -> bool {
    a.chromosome == b.chromosome && doubles_less_than(a.copy_number, b.copy_number)
}

/// Doubles.lessThan (epsilon 1e-10)
pub fn doubles_less_than(a: f64, b: f64) -> bool {
    a - b < -1e-10
}

/// Math.round(double)
pub fn java_round(x: f64) -> i64 {
    if x.is_nan() {
        return 0;
    }
    (x + 0.5).floor() as i64
}

// ---------------------------------------------------------------------------------------------------------
// Connectors

#[derive(Clone, Debug)]
pub struct Connector {
    pub chromosome: String,
    pub position: i32,
    pub frame: i32,
    pub cluster_id: i32,
    pub chain_id: i32,
    pub track: i32,
    pub ploidy: f64,
}

/// Connectors.createConnectors
pub fn create_connectors(show_simple: bool, segments: &[VisSegment], links: &[VisSvData]) -> Vec<Connector> {
    let mut result = Vec::new();
    for s in segments {
        let base = |position: i32, ploidy: f64, frame: i32| Connector {
            chromosome: s.chromosome.clone(),
            position,
            frame,
            cluster_id: s.cluster_id,
            chain_id: s.chain_id,
            track: s.track,
            ploidy,
        };
        if let Some(l) = find_link(&s.chromosome, s.start, links) {
            let before = segment_ploidy_before(s.track, &s.chromosome, s.start, segments);
            if l.jcn > 0.0 {
                result.push(base(s.start, f64::max(0.0, l.jcn - before), l.frame));
            }
        }
        if let Some(l) = find_link(&s.chromosome, s.end, links) {
            if l.jcn > 0.0 {
                let before = segment_ploidy_before(s.track, &s.chromosome, s.end, segments);
                result.push(base(s.end, f64::max(0.0, l.jcn - before), l.frame));
            }
        }
    }
    for l in links {
        if l.connectors_only(show_simple) {
            let c = |chr: &str, pos: i32| Connector {
                chromosome: chr.to_string(),
                position: pos,
                frame: 0,
                cluster_id: l.cluster_id,
                chain_id: l.chain_id,
                track: 0,
                ploidy: l.jcn,
            };
            if l.is_valid_start() {
                result.push(c(&l.chr_start, l.pos_start));
            }
            if l.is_valid_end() {
                result.push(c(&l.chr_end, l.pos_end));
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------------------------------------
// Highlights

pub struct KnownSites {
    pub centromeres: Vec<Region>,
    pub fragile_sites: Vec<Region>,
    pub line_elements: Vec<Region>,
}

fn load_config_regions(lines: &[&str], v37: bool, delim: char) -> Vec<Region> {
    let mut out = Vec::new();
    for line in lines {
        if line.contains("Chromosome") {
            continue;
        }
        let items: Vec<&str> = line.split(delim).collect();
        if items.len() < 3 {
            continue;
        }
        let chr = hmf_common::regions::versioned_chromosome(!v37, items[0]);
        let (Ok(s), Ok(e)) = (items[1].parse::<i32>(), items[2].parse::<i32>()) else { continue };
        out.push(Region::new(&chr, s, e));
    }
    out
}

impl KnownSites {
    /// Highlights.populateKnownSites
    pub fn load(v37: bool) -> KnownSites {
        let gaps: Vec<String> = hmf_common::ref_coords::centromere_gap_lines(v37).iter().map(|l| l.strip_prefix("hs").unwrap_or(l).to_string()).collect();
        let gaps: Vec<&str> = gaps.iter().map(|s| s.as_str()).collect();
        let fragile = crate::annotators::fragile_sites_resource(v37);
        let line = crate::annotators::line_elements_resource(v37);
        KnownSites {
            centromeres: load_config_regions(&gaps, v37, '\t'),
            fragile_sites: load_config_regions(&fragile.lines().collect::<Vec<_>>(), v37, ','),
            line_elements: load_config_regions(&line.lines().collect::<Vec<_>>(), v37, ','),
        }
    }
}

/// Highlights.limitHighlightsToRegions
pub fn limit_highlights(highlights: &[Region], segments: &[Region]) -> Vec<Region> {
    let mut out = Vec::new();
    for h in highlights {
        let segs: Vec<&Region> = segments.iter().filter(|s| s.chromosome == h.chromosome).collect();
        if segs.is_empty() {
            continue;
        }
        let min = segs.iter().map(|s| s.start).min().unwrap_or(0);
        let max = segs.iter().map(|s| s.end).max().unwrap_or(0);
        if h.end >= min && h.start <= max {
            out.push(Region::new(&h.chromosome, i32::max(min, h.start), i32::min(max, h.end)));
        }
    }
    out
}

// ---------------------------------------------------------------------------------------------------------
// ScaleContig, ScalePosition

#[derive(Clone, Debug)]
pub struct ScaleContig {
    pub contig: String,
    max_adjusted: i32,
    map: BTreeMap<i32, i32>,
}

/// ScaleContig.logDistance
pub fn log_distance(distance: i32) -> i32 {
    (distance as f64).log10().powf(3.0).floor() as i32 + 10
}

impl ScaleContig {
    pub fn new(contig: &str, positions: &[i32]) -> ScaleContig {
        let mut sorted: Vec<i32> = positions.to_vec();
        sorted.sort();
        sorted.dedup();
        let mut map = BTreeMap::new();
        let mut max_adjusted = 0;
        if !sorted.is_empty() {
            let mut current = 1;
            let mut prev = sorted[0];
            map.insert(prev, current);
            for &p in &sorted[1..] {
                current += log_distance(p - prev);
                map.insert(p, current);
                prev = p;
            }
            max_adjusted = current;
        }
        ScaleContig { contig: contig.to_string(), max_adjusted, map }
    }

    #[cfg(test)]
    pub fn from_map(contig: &str, map: BTreeMap<i32, i32>) -> ScaleContig {
        ScaleContig { contig: contig.to_string(), max_adjusted: 0, map }
    }

    pub fn length(&self) -> i32 {
        self.max_adjusted
    }

    pub fn expand(&mut self, factor: f64) {
        for v in self.map.values_mut() {
            if *v > 1 {
                *v = java_round(factor * *v as f64) as i32;
            }
        }
        self.max_adjusted = java_round(factor * self.max_adjusted as f64) as i32;
    }

    /// scale: the exact position (Java throws for any other)
    pub fn scale(&self, position: i32) -> anyhow::Result<i32> {
        self.map.get(&position).copied().ok_or_else(|| anyhow::anyhow!("Invalid position"))
    }

    pub fn interpolate(&self, value: i32) -> i32 {
        if let Some(&v) = self.map.get(&value) {
            return v;
        }
        let min_value = self.map.keys().next().copied().unwrap_or(0);
        let max_value = self.map.keys().next_back().copied().unwrap_or(0);
        let start = self.map.range(..value).next_back().map(|(k, _)| *k).unwrap_or(min_value);
        let end = self.map.range(value + 1..).next().map(|(k, _)| *k).unwrap_or(max_value);
        if start == end {
            return self.map[&start];
        }
        let proportion = (value as i64 - start as i64).abs() as f64 / ((end as i64 - start as i64).abs() as f64);
        let a = self.map[&start];
        let b = self.map[&end];
        a + (proportion * (b - a).abs() as f64).floor() as i32
    }
}

pub struct ScalePosition {
    pub total_length: i32,
    /// (contig, length), one per contig
    pub contig_lengths: Vec<GenomePosition>,
    pub contigs: HashMap<String, ScaleContig>,
}

const MIN_CONTIG_PERCENTAGE: f64 = 0.01;
const GENE_NAME_DISTANCE: f64 = 1.5 / 360.0;

impl ScalePosition {
    pub fn new(positions: &[GenomePosition]) -> ScalePosition {
        let mut contigs: HashMap<String, ScaleContig> = HashMap::new();
        let mut names: Vec<String> = Vec::new();
        let mut total = 0;
        for p in positions {
            if !names.contains(&p.chromosome) {
                names.push(p.chromosome.clone());
            }
        }
        for n in &names {
            let ps: Vec<i32> = positions.iter().filter(|p| p.chromosome == *n).map(|p| p.position).collect();
            let sc = ScaleContig::new(n, &ps);
            total += sc.length();
            contigs.insert(n.clone(), sc);
        }
        let min_contig_distance = java_round(MIN_CONTIG_PERCENTAGE * total as f64);
        let mut total = 0;
        let mut lengths = Vec::new();
        for n in &names {
            let sc = contigs.get_mut(n).unwrap();
            if (sc.length() as i64) < min_contig_distance {
                let factor = min_contig_distance as f64 / sc.length() as f64;
                sc.expand(factor);
            }
            total += sc.length();
            lengths.push(GenomePosition::new(n, sc.length()));
        }
        ScalePosition { total_length: total, contig_lengths: lengths, contigs }
    }

    pub fn contig(&self, chr: &str) -> anyhow::Result<&ScaleContig> {
        self.contigs.get(chr).ok_or_else(|| anyhow::anyhow!("NullPointerException: no scale for {chr}"))
    }
}

// ---------------------------------------------------------------------------------------------------------
// Downsampler

const DOWNSAMPLE_TARGET: usize = 5000;
const MIN_PRESERVED: i32 = 5;

pub fn downsample_list<T: Clone>(list: &[T], target: i32) -> Vec<T> {
    if list.len() as i64 <= target as i64 {
        return list.to_vec();
    }
    let step = list.len() as f64 / target as f64;
    let mut out = Vec::new();
    for i in 0..target.max(0) {
        let idx = (i as f64 * step).floor() as usize;
        if idx >= list.len() {
            break;
        }
        out.push(list[idx].clone());
    }
    out
}

pub fn group_downsample_targets(sizes: &[i32], target: i32, min_preserved: i32) -> Vec<i32> {
    let n = sizes.len();
    let mut reserved = vec![0i32; n];
    let mut extra = vec![0i32; n];
    for i in 0..n {
        if sizes[i] > 0 {
            reserved[i] = min_preserved;
            extra[i] = sizes[i] - min_preserved;
        }
    }
    let reserved_sum: i32 = reserved.iter().sum();
    let extra_sum: i32 = extra.iter().sum();
    let distribute = target - reserved_sum;
    (0..n)
        .map(|i| {
            let fraction = extra[i] as f64 / extra_sum as f64;
            reserved[i] + java_round(fraction * distribute as f64) as i32
        })
        .collect()
}

pub fn regions_from_positions(positions: &[GenomePosition]) -> Vec<Region> {
    let mut chrs: Vec<&str> = Vec::new();
    for p in positions {
        if !chrs.contains(&p.chromosome.as_str()) {
            chrs.push(&p.chromosome);
        }
    }
    let mut regions = Vec::new();
    for c in chrs {
        let mut ps: Vec<GenomePosition> = positions.iter().filter(|p| p.chromosome == c).cloned().collect();
        ps.push(GenomePosition::new(c, 1));
        ps.push(GenomePosition::new(c, i32::MAX));
        let mut ps = distinct_positions(ps);
        sort_positions(&mut ps);
        for i in 1..ps.len() {
            let mut start = ps[i - 1].position;
            if start == 1 {
                start = 0;
            }
            regions.push(Region::new(c, start, ps[i].position));
        }
    }
    regions
}

/// Downsampler.downsampleWithMinimumPerRegion for (chromosome, position, value) points
pub fn downsample_points(points: &[(String, i32, f64)], positions: &[GenomePosition]) -> Vec<(String, i32, f64)> {
    if points.len() < DOWNSAMPLE_TARGET {
        return points.to_vec();
    }
    let regions = regions_from_positions(positions);
    let groups: Vec<Vec<(String, i32, f64)>> = regions
        .iter()
        .map(|r| points.iter().filter(|p| p.0 == r.chromosome && p.1 > r.start && p.1 <= r.end).cloned().collect())
        .collect();
    let sizes: Vec<i32> = groups.iter().map(|g| g.len() as i32).collect();
    let targets = group_downsample_targets(&sizes, DOWNSAMPLE_TARGET as i32, MIN_PRESERVED);
    let mut out = Vec::new();
    for (g, t) in groups.iter().zip(targets) {
        out.extend(downsample_list(g, t));
    }
    out
}

/// CircosData.centerCobaltPositions: a COBALT window's position moved to its middle, or to the middle
/// breakpoint (or middle segment) of the plot positions it contains
pub fn center_cobalt_positions(ratios: &[(String, i32, f64)], positions: &[GenomePosition]) -> Vec<(String, i32, f64)> {
    const WINDOW_SIZE: i32 = 1000;
    ratios
        .iter()
        .map(|(chr, pos, v)| {
            let (s, e) = (*pos, *pos + WINDOW_SIZE);
            let overlapping: Vec<i32> = positions.iter().filter(|p| p.chromosome == *chr && s <= p.position && e >= p.position).map(|p| p.position).collect();
            let mid = if overlapping.is_empty() {
                pos + WINDOW_SIZE / 2
            } else {
                let m = overlapping.len() / 2;
                if !overlapping.len().is_multiple_of(2) { overlapping[m] } else { (overlapping[m - 1] + overlapping[m]) / 2 }
            };
            (chr.clone(), mid, *v)
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------------------
// ColorPicker, Thickness

const BLACK: (u8, u8, u8) = (5, 5, 5);
const COLOURS: [(u8, u8, u8); 8] = [(106, 61, 154), (140, 81, 10), (1, 102, 94), (255, 127, 0), (212, 193, 23), (31, 120, 180), (51, 160, 44), (152, 51, 160)];
const DEL: (u8, u8, u8) = (251, 154, 153);
const DUP: (u8, u8, u8) = (178, 223, 138);
const INS: (u8, u8, u8) = (255, 255, 153);
const LINE: (u8, u8, u8) = (97, 171, 227);
const DOUBLE_MINUTE: (u8, u8, u8) = (227, 26, 28);

fn colour_str(c: (u8, u8, u8)) -> String {
    format!("color=({},{},{})", c.0, c.1, c.2)
}

fn simple_sv_colour(sv_type: &str) -> String {
    match sv_type {
        "DEL" => colour_str(DEL),
        "DUP" => colour_str(DUP),
        _ => colour_str(INS),
    }
}

/// ColorPicker.contigColour as "#RRGGBB" (hexContigColor)
pub fn hex_contig_colour(contig: &str) -> String {
    const MAP: [(&str, (u8, u8, u8)); 24] = [
        ("1", (128, 125, 186)), ("2", (145, 142, 179)), ("3", (161, 159, 173)), ("4", (179, 176, 166)), ("5", (196, 193, 160)), ("6", (213, 210, 153)),
        ("7", (230, 228, 147)), ("8", (202, 218, 138)), ("9", (175, 209, 129)), ("10", (147, 199, 120)), ("11", (120, 190, 111)), ("12", (92, 180, 102)),
        ("13", (65, 171, 93)), ("14", (65, 166, 110)), ("15", (65, 162, 128)), ("16", (65, 158, 145)), ("17", (65, 154, 163)), ("18", (65, 150, 180)),
        ("19", (66, 146, 198)), ("20", (76, 142, 196)), ("21", (86, 139, 194)), ("22", (97, 135, 192)), ("X", (107, 132, 190)), ("Y", (117, 128, 188)),
    ];
    let key = contig.replace("chr", "");
    let c = MAP.iter().find(|(k, _)| *k == key).map(|(_, c)| *c).unwrap_or((0, 0, 0));
    format!("#{:02X}{:02X}{:02X}", c.0, c.1, c.2)
}

/// Collectors.toMap(x -> x, x -> 1L, Math::addExact) over integer keys, iterated: the keys with their counts
/// in Java's HashMap order. HashMap.merge resizes at the start of a call once the size exceeds the threshold
/// (so a map that has just reached 13 keys keeps 16 buckets until the next call), puts a new key at the head
/// of its bucket, and a resize keeps each bucket's order.
pub fn merge_count_order(keys: &[i32]) -> Vec<(i32, i64)> {
    let mut cap: usize = 16;
    let mut order: Vec<(i32, i64)> = Vec::new(); // insertion order
    for &k in keys {
        if order.len() > cap * 3 / 4 {
            cap *= 2;
        }
        match order.iter_mut().find(|e| e.0 == k) {
            Some(e) => e.1 += 1,
            None => order.push((k, 1)),
        }
    }
    let spread = |h: i32| ((h ^ ((h as u32) >> 16) as i32) as u32 as usize) & (cap - 1);
    let mut idx: Vec<usize> = (0..order.len()).collect();
    // by bucket, then the later inserted first
    idx.sort_by(|&a, &b| spread(order[a].0).cmp(&spread(order[b].0)).then(b.cmp(&a)));
    idx.into_iter().map(|i| order[i]).collect()
}

pub struct ColorPicker {
    cluster_mode: bool,
    map: HashMap<i32, String>,
    transparency: f64,
}

fn connector_transparency(links: usize) -> f64 {
    if links < 10 {
        1.0
    } else if links < 50 {
        0.8
    } else if links < 100 {
        0.6
    } else if links < 200 {
        0.4
    } else if links < 400 {
        0.2
    } else {
        0.1
    }
}

impl ColorPicker {
    /// ColorPicker.clusterColors: the largest clusters get the colours in turn (ties in the iteration order
    /// of Java's HashMap of cluster ids), simple SVs and LINE clusters their own
    pub fn cluster_colours(links: &[VisSvData]) -> ColorPicker {
        let mut sizes: Vec<(i32, i64)> = merge_count_order(&links.iter().map(|l| l.cluster_id).collect::<Vec<_>>());
        sizes.sort_by_key(|s| std::cmp::Reverse(s.1));
        let mut map = HashMap::new();
        for (i, (id, _)) in sizes.iter().enumerate() {
            map.insert(*id, if i < COLOURS.len() { colour_str(COLOURS[i]) } else { "color=black".to_string() });
        }
        for l in links {
            if l.is_simple_sv() {
                map.insert(l.cluster_id, simple_sv_colour(l.sv_type));
            } else if l.is_line_element() {
                map.insert(l.cluster_id, colour_str(LINE));
            }
        }
        ColorPicker { cluster_mode: true, map, transparency: connector_transparency(links.len()) }
    }

    /// ColorPicker.chainColors
    pub fn chain_colours(links: &[VisSvData]) -> ColorPicker {
        let mut map = HashMap::new();
        if let Some(first) = links.first() {
            if first.is_simple_sv() && !first.in_double_minute {
                let c = simple_sv_colour(first.sv_type);
                for l in links {
                    map.insert(l.chain_id, c.clone());
                }
            } else if first.is_line_element() {
                for l in links {
                    map.insert(l.chain_id, colour_str(LINE));
                }
            } else {
                let mut dm: Vec<i32> = Vec::new();
                let mut chains: Vec<i32> = Vec::new();
                for l in links {
                    if l.in_double_minute && !dm.contains(&l.chain_id) {
                        dm.push(l.chain_id);
                    }
                    if !chains.contains(&l.chain_id) {
                        chains.push(l.chain_id);
                    }
                }
                for (i, id) in chains.iter().enumerate() {
                    let c = if dm.contains(id) {
                        colour_str(DOUBLE_MINUTE)
                    } else if i < COLOURS.len() {
                        colour_str(COLOURS[i])
                    } else {
                        colour_str(BLACK)
                    };
                    map.insert(*id, c);
                }
            }
        }
        ColorPicker { cluster_mode: false, map, transparency: connector_transparency(links.len()) }
    }

    /// color(clusterId, chainId); "null" where Java's map has no entry
    pub fn colour(&self, cluster_id: i32, chain_id: i32) -> String {
        if cluster_id == -1 || chain_id == -1 {
            return colour_str(BLACK);
        }
        let key = if self.cluster_mode { cluster_id } else { chain_id };
        self.map.get(&key).cloned().unwrap_or_else(|| "null".to_string())
    }

    pub fn transparent_colour(&self, cluster_id: i32, chain_id: i32) -> String {
        let c = self.colour(cluster_id, chain_id);
        c.replace(')', &format!(",{})", hmf_common::fmt::java_double_to_string(self.transparency)))
    }
}

pub struct Thickness {
    min: f64,
    max: f64,
    gradient: f64,
}

const MAX_TOTAL_PLOIDY: f64 = 1000.0;

fn bound(v: f64, min: f64, max: f64) -> f64 {
    f64::min(f64::max(v, min), max)
}

impl Thickness {
    pub fn new(min: f64, max: f64, connectors: &[Connector]) -> Thickness {
        let max_ploidy = bound(connectors.iter().map(|c| c.ploidy).fold(None, |m: Option<f64>, v| Some(m.map_or(v, |m| m.max(v)))).unwrap_or(0.0), 6.0, 60.0);
        let unadjusted = (max - min) / (max_ploidy - 1.0);
        let mut t = Thickness { min, max, gradient: unadjusted };
        let total: f64 = connectors.iter().map(|c| t.pixels_with(c.ploidy, unadjusted)).sum();
        t.gradient = if connectors.len() as f64 - MAX_TOTAL_PLOIDY >= -1e-10 {
            0.0
        } else if total - MAX_TOTAL_PLOIDY > 1e-10 {
            (MAX_TOTAL_PLOIDY - connectors.len() as f64) / (total - connectors.len() as f64) * unadjusted
        } else {
            unadjusted
        };
        t
    }

    fn pixels_with(&self, ploidy: f64, gradient: f64) -> f64 {
        java_round(bound(gradient * ploidy + (self.min - gradient), self.min, self.max)) as f64
    }

    pub fn pixels(&self, ploidy: f64) -> f64 {
        self.pixels_with(ploidy, self.gradient)
    }
}

// ---------------------------------------------------------------------------------------------------------
// Genes and exons

#[derive(Clone, Debug)]
pub struct Gene {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub name: String,
    pub annotation_type: AnnotationType,
    pub transcript: String,
    pub name_position: i64,
    pub strand: i32,
}

impl GenomeRegion for Gene {
    fn chromosome(&self) -> &str {
        &self.chromosome
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn end(&self) -> i32 {
        self.end
    }
}

/// GeneUtils.genes: one gene per transcript (not exon-lost), from its first and last exon by position
pub fn genes(exons: &[VisGeneExon]) -> Vec<Gene> {
    let mut transcripts: JavaHashSet<String> = JavaHashSet::new();
    for e in exons.iter().filter(|e| e.annotation_type != AnnotationType::EXON_LOST) {
        transcripts.insert(e.transcript.clone());
    }
    let mut result = Vec::new();
    for t in transcripts.iter() {
        let mut te: Vec<&VisGeneExon> = exons.iter().filter(|e| e.transcript == *t).collect();
        te.sort_by(|a, b| region_cmp(*a, *b));
        let first = te[0];
        let last = te[te.len() - 1];
        let name_position = if first.exon_rank <= last.exon_rank { first.exon_start } else { last.exon_end } as i64;
        result.push(Gene {
            chromosome: first.chromosome.clone(),
            start: first.exon_start,
            end: last.exon_end,
            name: first.gene.clone(),
            annotation_type: first.annotation_type,
            transcript: first.transcript.clone(),
            name_position,
            strand: if first.exon_rank <= last.exon_rank { 1 } else { -1 },
        });
    }
    result
}

/// GeneUtils.unique: one gene per name (the longer, or the later of equal length), in Java's HashMap order
pub fn unique_genes(genes: Vec<Gene>) -> Vec<Gene> {
    let mut m: JavaHashMap<String, Gene> = JavaHashMap::new();
    for g in genes {
        if let Some(prior) = m.get(&g.name) {
            if g.bases() < prior.bases() {
                continue;
            }
        }
        m.insert(g.name.clone(), g);
    }
    m.values().cloned().collect()
}

/// VisExons.transcriptExons: the exons of a transcript, the first of each rank
fn transcript_exons(transcript: &str, exons: &[VisGeneExon]) -> Vec<VisGeneExon> {
    let mut ranks: Vec<i32> = Vec::new();
    let mut out = Vec::new();
    for e in exons.iter().filter(|e| e.transcript == transcript) {
        if !ranks.contains(&e.exon_rank) {
            ranks.push(e.exon_rank);
            out.push(e.clone());
        }
    }
    out
}

/// VisExons.geneExons
pub fn gene_exons(genes: &[Gene], exons: &[VisGeneExon]) -> Vec<VisGeneExon> {
    genes.iter().flat_map(|g| transcript_exons(&g.transcript, exons)).collect()
}

/// VisExons.fusionExons
pub fn fusion_exons(f: &VisFusion, exons: &[VisGeneExon]) -> Vec<VisGeneExon> {
    let mut out = transcript_exons(&f.transcript_up, exons);
    if f.transcript_up != f.transcript_down {
        out.extend(transcript_exons(&f.transcript_down, exons));
    }
    out
}

fn sorted_by_rank(mut v: Vec<VisGeneExon>) -> Vec<VisGeneExon> {
    v.sort_by_key(|e| e.exon_rank);
    v
}

pub fn sorted_upstream_exons(f: &VisFusion, exons: &[VisGeneExon]) -> Vec<VisGeneExon> {
    sorted_by_rank(exons.iter().filter(|e| e.gene == f.gene_name_up && e.transcript == f.transcript_up).cloned().collect())
}

pub fn sorted_downstream_exons(f: &VisFusion, exons: &[VisGeneExon]) -> Vec<VisGeneExon> {
    sorted_by_rank(exons.iter().filter(|e| e.gene == f.gene_name_down && e.transcript == f.transcript_down).cloned().collect())
}

/// DisruptedExons.disruptedGeneRegions
pub fn disrupted_gene_regions(fusions: &[VisFusion], exons: &[VisGeneExon]) -> Vec<Region> {
    let mut out: Vec<Region> =
        exons.iter().filter(|e| e.annotation_type == AnnotationType::EXON_LOST).map(|e| Region::new(&e.chromosome, e.exon_start, e.exon_end)).collect();
    for f in fusions {
        let up: Vec<VisGeneExon> = sorted_upstream_exons(f, exons).into_iter().filter(|e| e.exon_rank >= f.fused_exon_up).collect();
        let down: Vec<VisGeneExon> = sorted_downstream_exons(f, exons).into_iter().filter(|e| e.exon_rank <= f.fused_exon_down).collect();
        if up.is_empty() || down.is_empty() {
            continue;
        }
        let (inc_up, exc_up) = (&up[0], &up[up.len() - 1]);
        let up_region = if f.strand_up < 0 {
            Region::new(&exc_up.chromosome, exc_up.exon_start, i32::max(f.pos_up, inc_up.exon_start))
        } else {
            Region::new(&exc_up.chromosome, i32::min(f.pos_up, inc_up.exon_end), exc_up.exon_end)
        };
        let (exc_down, inc_down) = (&down[0], &down[down.len() - 1]);
        let down_region = if f.strand_down < 0 {
            Region::new(&exc_down.chromosome, i32::min(f.pos_down, inc_down.exon_end), exc_down.exon_end)
        } else {
            Region::new(&exc_down.chromosome, exc_down.exon_start, i32::max(f.pos_down, inc_down.exon_start))
        };
        out.push(up_region);
        out.push(down_region);
    }
    out
}
