//! GenomePosition, GenomeRegion and HumanChromosome as the visualiser uses them, and `circos.Span`.

use std::cmp::Ordering;

pub fn strip_chr(chromosome: &str) -> &str {
    chromosome.strip_prefix("chr").unwrap_or(chromosome)
}

/// HumanChromosome.contains: 1-22, X and Y, with or without "chr"
pub fn human_chromosome_contains(chromosome: &str) -> bool {
    let t = strip_chr(chromosome);
    if t.chars().all(|c| c.is_ascii_digit()) {
        return match t.parse::<i64>() {
            Ok(n) => (1..=22).contains(&n),
            Err(_) => false,
        };
    }
    t == "X" || t == "Y"
}

/// HumanChromosome.chromosomeRank
pub fn chromosome_rank(chromosome: &str) -> i32 {
    let t = strip_chr(chromosome);
    if t.eq_ignore_ascii_case("X") {
        23
    } else if t.eq_ignore_ascii_case("Y") {
        24
    } else if t.eq_ignore_ascii_case("MT") || t.eq_ignore_ascii_case("M") {
        25
    } else {
        t.parse::<i32>().unwrap_or(26)
    }
}

/// HumanChromosome.fromString(chromosome).toString(): "1".."22", "X", "Y" (a "chr" prefix of any case is
/// removed); None for any other name (Java throws).
pub fn human_chromosome_name(chromosome: &str) -> Option<&str> {
    let t = if chromosome.len() >= 3 && chromosome[..3].eq_ignore_ascii_case("chr") { &chromosome[3..] } else { chromosome };
    const NAMES: [&str; 24] = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"];
    NAMES.iter().find(|n| **n == t).copied()
}

/// CircosDataWriter.circosContig: "hs" + the chromosome's short name
pub fn circos_contig(chromosome: &str) -> String {
    format!("hs{}", human_chromosome_name(chromosome).unwrap_or(chromosome))
}

/// ContigComparator
pub fn contig_cmp(a: &str, b: &str) -> Ordering {
    chromosome_rank(a).cmp(&chromosome_rank(b))
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct GenomePosition {
    pub chromosome: String,
    pub position: i32,
}

impl GenomePosition {
    pub fn new(chromosome: &str, position: i32) -> GenomePosition {
        GenomePosition { chromosome: chromosome.to_string(), position }
    }
}

/// GenomePosition.compare
pub fn position_cmp(a: &GenomePosition, b: &GenomePosition) -> Ordering {
    match contig_cmp(&a.chromosome, &b.chromosome) {
        Ordering::Equal => a.position.cmp(&b.position),
        o => o,
    }
}

pub trait GenomeRegion {
    fn chromosome(&self) -> &str;
    fn start(&self) -> i32;
    fn end(&self) -> i32;
    fn bases(&self) -> i32 {
        1 + self.end() - self.start()
    }
    fn overlaps_region(&self, other: &dyn GenomeRegion) -> bool {
        other.chromosome() == self.chromosome() && other.end() > self.start() && other.start() < self.end()
    }
    fn contains_position(&self, p: &GenomePosition) -> bool {
        self.chromosome() == p.chromosome && self.start() <= p.position && self.end() >= p.position
    }
}

/// GenomeRegion.compareTo: by chromosome rank (when the names differ), then start only
pub fn region_cmp<A: GenomeRegion + ?Sized, B: GenomeRegion + ?Sized>(a: &A, b: &B) -> Ordering {
    if a.chromosome() == b.chromosome() {
        a.start().cmp(&b.start())
    } else {
        contig_cmp(a.chromosome(), b.chromosome())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

impl Region {
    pub fn new(chromosome: &str, start: i32, end: i32) -> Region {
        Region { chromosome: chromosome.to_string(), start, end }
    }
}

impl GenomeRegion for Region {
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

/// Collections.sort of positions
pub fn sort_positions(p: &mut [GenomePosition]) {
    p.sort_by(position_cmp);
}

/// stream().distinct(): the first of equal positions, in order
pub fn distinct_positions(p: Vec<GenomePosition>) -> Vec<GenomePosition> {
    let mut seen = std::collections::HashSet::new();
    p.into_iter().filter(|x| seen.insert(x.clone())).collect()
}

/// Span.spanPositions: per chromosome (in order of first appearance) the min and max position, sorted
pub fn span_positions(positions: &[GenomePosition]) -> Vec<Region> {
    let mut chrs: Vec<&str> = Vec::new();
    for p in positions {
        if !chrs.contains(&p.chromosome.as_str()) {
            chrs.push(&p.chromosome);
        }
    }
    let mut out: Vec<Region> = chrs
        .iter()
        .map(|c| {
            let ps = positions.iter().filter(|p| p.chromosome == *c).map(|p| p.position);
            let (mut mn, mut mx) = (i32::MAX, i32::MIN);
            for v in ps {
                mn = mn.min(v);
                mx = mx.max(v);
            }
            Region::new(c, mn, mx)
        })
        .collect();
    out.sort_by(region_cmp);
    out
}

/// Span.spanRegions: per chromosome the min start and max end, sorted. Java collects the chromosomes into a
/// HashSet first, but the stable sort by rank and start gives one order for distinct chromosomes.
pub fn span_regions<R: GenomeRegion>(regions: &[R]) -> Vec<Region> {
    let mut chrs: Vec<&str> = Vec::new();
    for r in regions {
        if !chrs.contains(&r.chromosome()) {
            chrs.push(r.chromosome());
        }
    }
    let mut out: Vec<Region> = chrs
        .iter()
        .map(|c| {
            let rs = regions.iter().filter(|r| r.chromosome() == *c);
            let (mut mn, mut mx) = (i32::MAX, i32::MIN);
            for r in rs {
                mn = mn.min(r.start());
                mx = mx.max(r.end());
            }
            Region::new(c, mn, mx)
        })
        .collect();
    out.sort_by(region_cmp);
    out
}

/// Span.minMaxPositions
pub fn min_max_positions<R: GenomeRegion>(regions: &[R]) -> Vec<GenomePosition> {
    let mut out = Vec::new();
    for r in span_regions(regions) {
        out.push(GenomePosition::new(&r.chromosome, r.start));
        out.push(GenomePosition::new(&r.chromosome, r.end));
    }
    out
}

/// Span.allPositions: the start and end of every region on a human chromosome, sorted
pub fn all_positions<R: GenomeRegion>(regions: &[R]) -> Vec<GenomePosition> {
    let mut out = Vec::new();
    for r in regions {
        if human_chromosome_contains(r.chromosome()) {
            out.push(GenomePosition::new(r.chromosome(), r.start()));
            out.push(GenomePosition::new(r.chromosome(), r.end()));
        }
    }
    sort_positions(&mut out);
    out
}

/// Span.minPositionPerChromosome / maxPositionPerChromosome
pub fn min_max_per_chromosome(positions: &[GenomePosition]) -> std::collections::HashMap<String, (i32, i32)> {
    let mut m: std::collections::HashMap<String, (i32, i32)> = std::collections::HashMap::new();
    for p in positions {
        let e = m.entry(p.chromosome.clone()).or_insert((p.position, p.position));
        e.0 = e.0.min(p.position);
        e.1 = e.1.max(p.position);
    }
    m
}
