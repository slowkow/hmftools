//! Chromosome helpers, regions and the -specific_regions / -specific_chr options.

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChrBaseRegion {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

impl ChrBaseRegion {
    pub fn new(chromosome: &str, start: i32, end: i32) -> Self {
        ChrBaseRegion { chromosome: chromosome.to_string(), start, end }
    }
    pub fn base_length(&self) -> i32 {
        self.end - self.start + 1
    }
    pub fn contains_position(&self, pos: i32) -> bool {
        pos >= self.start && pos <= self.end
    }
    pub fn overlaps(&self, chr: &str, start: i32, end: i32) -> bool {
        self.chromosome == chr && self.start <= end && self.end >= start
    }
}

impl std::fmt::Display for ChrBaseRegion {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chromosome, self.start, self.end)
    }
}

pub const HUMAN_CHROMOSOMES: [&str; 24] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y",
];

fn strip_chr(contig: &str) -> &str {
    if contig.len() > 3 && contig[..3].eq_ignore_ascii_case("chr") { &contig[3..] } else { contig }
}

/// HumanChromosome.contains
pub fn is_human_chromosome(contig: &str) -> bool {
    HUMAN_CHROMOSOMES.contains(&strip_chr(contig))
}

/// MitochondrialChromosome.contains
pub fn is_mitochondrial(contig: &str) -> bool {
    contig == "MT" || contig == "chrM" || contig == "M" || contig == "chrMT"
}

/// Chromosome.isAltRegionContig
pub fn is_alt_region_contig(contig: &str) -> bool {
    !is_human_chromosome(contig) && !is_mitochondrial(contig)
}

pub fn versioned_chromosome(v38: bool, chr: &str) -> String {
    let bare = strip_chr(chr);
    if v38 { format!("chr{bare}") } else { bare.to_string() }
}

#[derive(Clone, Debug, Default)]
pub struct SpecificRegions {
    pub regions: Vec<ChrBaseRegion>,
    pub chromosomes: Vec<String>,
}

impl SpecificRegions {
    /// `-specific_regions chr:start-end;...` (or chr:start:end) and `-specific_chr a;b`
    pub fn from(regions: Option<&str>, chromosomes: Option<&str>) -> anyhow::Result<Self> {
        let mut s = SpecificRegions::default();
        if let Some(r) = regions.filter(|r| !r.is_empty()) {
            for item in r.split(';') {
                let (chr, range) = item.split_once(':').ok_or_else(|| anyhow::anyhow!("invalid region: {item}"))?;
                let (a, b) = range.split_once('-').or_else(|| range.split_once(':')).ok_or_else(|| anyhow::anyhow!("invalid region: {item}"))?;
                s.regions.push(ChrBaseRegion::new(chr, a.parse()?, b.parse()?));
            }
            // ChrBaseRegion.compareTo: chromosome rank, then start
            s.regions.sort_by(|a, b| if a.chromosome == b.chromosome { a.start.cmp(&b.start) } else { compare_chromosomes(&a.chromosome, &b.chromosome) });
            for r in &s.regions {
                if !s.chromosomes.contains(&r.chromosome) {
                    s.chromosomes.push(r.chromosome.clone());
                }
            }
        } else if let Some(c) = chromosomes.filter(|c| !c.is_empty()) {
            s.chromosomes = c.split(';').map(|x| x.to_string()).collect();
        }
        Ok(s)
    }
    pub fn has_filters(&self) -> bool {
        !self.regions.is_empty() || !self.chromosomes.is_empty()
    }
    pub fn include_chromosome(&self, chr: &str) -> bool {
        self.chromosomes.is_empty() || self.chromosomes.iter().any(|c| c == chr)
    }
}

/// HumanChromosome.chromosomeRank: 1-22, X 23, Y 24, MT 25, any other contig 26
pub fn chromosome_rank(chr: &str) -> i32 {
    let bare = chr.strip_prefix("chr").unwrap_or(chr);
    if bare.eq_ignore_ascii_case("X") {
        23
    } else if bare.eq_ignore_ascii_case("Y") {
        24
    } else if bare.eq_ignore_ascii_case("MT") || bare.eq_ignore_ascii_case("M") {
        25
    } else {
        bare.parse::<i32>().unwrap_or(26)
    }
}

/// ChrBaseRegion.compareChromosomes
pub fn compare_chromosomes(a: &str, b: &str) -> std::cmp::Ordering {
    let (ra, rb) = (chromosome_rank(a), chromosome_rank(b));
    if ra == rb { a.cmp(b) } else { ra.cmp(&rb) }
}

