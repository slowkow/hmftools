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


// ---- additions for BamTools and later ports (2026-11-26): Java-exact helpers of hmf-common `region`,
// `genome.refgenome` and `genome.chromosome` ----

use crate::java_hash::JavaHashMap;
use crate::java_text;
use crate::segmentation::HumanChromosome;
use anyhow::{anyhow, bail, Context};

/// RefGenomeVersion
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RefGenomeVersion {
    V37,
    V38,
}

/// RefGenomeCoordinates lengths (hmf-common resources lengths.37.tsv and lengths.38.tsv), 1-22, X, Y
const LENGTHS_37: [i32; 24] = [
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566,
];
const LENGTHS_38: [i32; 24] = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415,
];

impl RefGenomeVersion {
    /// RefGenomeVersion.from(String): "V37", "37" or "HG37"; "V38", "38" or "HG38"
    pub fn from_java(s: &str) -> anyhow::Result<Self> {
        match s {
            "V37" | "37" | "HG37" => Ok(RefGenomeVersion::V37),
            "V38" | "38" | "HG38" => Ok(RefGenomeVersion::V38),
            _ => bail!("Cannot resolve ref genome version: {s}"),
        }
    }
    pub fn is_37(self) -> bool {
        self == RefGenomeVersion::V37
    }
    pub fn is_38(self) -> bool {
        self == RefGenomeVersion::V38
    }
    /// identifier(): "37" or "38", as in resource file names
    pub fn identifier(self) -> &'static str {
        if self.is_37() { "37" } else { "38" }
    }
    /// the enum name, as Java logs it
    pub fn name(self) -> &'static str {
        if self.is_37() { "V37" } else { "V38" }
    }
    /// versionedChromosome: "chr" added for 38 where missing, a leading "chr" removed for 37
    pub fn versioned_chromosome(self, chromosome: &str) -> String {
        if self.is_38() { enforce_chr_prefix(chromosome) } else { strip_chr_prefix(chromosome).to_string() }
    }
    /// RefGenomeCoordinates.length(chromosome): 0 for a chromosome HumanChromosome.contains rejects
    pub fn chromosome_length(self, chromosome: &str) -> i32 {
        if !HumanChromosome::contains(chromosome) {
            return 0;
        }
        let Some(c) = HumanChromosome::from_name(chromosome) else { return 0 };
        if self.is_37() { LENGTHS_37[c.ordinal()] } else { LENGTHS_38[c.ordinal()] }
    }
}

/// RefGenomeFunctions.stripChrPrefix: a leading lower-case "chr" removed
pub fn strip_chr_prefix(chromosome: &str) -> &str {
    chromosome.strip_prefix("chr").unwrap_or(chromosome)
}

/// RefGenomeFunctions.enforceChrPrefix
pub fn enforce_chr_prefix(chromosome: &str) -> String {
    if chromosome.starts_with("chr") { chromosome.to_string() } else { format!("chr{chromosome}") }
}

/// BaseRegion: 1-based inclusive positions without a chromosome
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub struct BaseRegion {
    pub start: i32,
    pub end: i32,
}

impl BaseRegion {
    pub fn new(start: i32, end: i32) -> Self {
        BaseRegion { start, end }
    }
    pub fn base_length(&self) -> i32 {
        self.end - self.start + 1
    }
    pub fn contains_position(&self, pos: i32) -> bool {
        position_within(pos, self.start, self.end)
    }
    pub fn overlaps(&self, start: i32, end: i32) -> bool {
        positions_overlap(self.start, self.end, start, end)
    }
}

/// BaseRegion.positionsOverlap
#[inline]
pub fn positions_overlap(start1: i32, end1: i32, start2: i32, end2: i32) -> bool {
    !(start1 > end2 || end1 < start2)
}

/// BaseRegion.positionWithin
#[inline]
pub fn position_within(pos: i32, start: i32, end: i32) -> bool {
    pos >= start && pos <= end
}

/// ChrBaseRegion.compareTo: chromosome rank (then name), then start only (ties are equal, so a stable sort
/// keeps their order)
pub fn compare_chr_regions(a: &ChrBaseRegion, b: &ChrBaseRegion) -> std::cmp::Ordering {
    if a.chromosome == b.chromosome { a.start.cmp(&b.start) } else { compare_chromosomes_java(&a.chromosome, &b.chromosome) }
}

/// ChrBaseRegion.compareChromosomes with HumanChromosome.chromosomeRank as Java has it: a lower-case
/// "chr" prefix stripped, then X 23, Y 24, MT or M 25 (any case), a number its value, anything else 26;
/// equal ranks compare the names (String.compareTo)
pub fn compare_chromosomes_java(a: &str, b: &str) -> std::cmp::Ordering {
    let (ra, rb) = (chromosome_rank_java(a), chromosome_rank_java(b));
    if ra == rb { java_text::compare(a, b) } else { ra.cmp(&rb) }
}

/// HumanChromosome.chromosomeRank
pub fn chromosome_rank_java(chromosome: &str) -> i32 {
    let t = strip_chr_prefix(chromosome);
    if t.eq_ignore_ascii_case("X") {
        23
    } else if t.eq_ignore_ascii_case("Y") {
        24
    } else if t.eq_ignore_ascii_case("MT") || t.eq_ignore_ascii_case("M") {
        25
    } else if !t.is_empty() && t.bytes().all(|b| b.is_ascii_digit()) {
        t.parse::<i32>().unwrap_or(26)
    } else {
        26
    }
}

/// BaseRegion.checkMergeOverlaps(regions, true): a stable sort by start, then each region absorbs the
/// next one that starts within 2 bases of its end. Java sets the merged end to the next region's end even
/// when that is smaller, so a region lying inside the one before cuts it short; `keep_larger_end` keeps
/// the larger end instead (ChrBaseRegion.checkMergeOverlaps does that).
pub fn merge_base_regions(regions: &mut Vec<BaseRegion>, keep_larger_end: bool) {
    regions.sort_by_key(|r| r.start);
    let mut i = 0;
    while i + 1 < regions.len() {
        let next = regions[i + 1];
        if regions[i].end >= next.start - 2 {
            if !keep_larger_end || next.end > regions[i].end {
                regions[i].end = next.end;
            }
            regions.remove(i + 1);
        } else {
            i += 1;
        }
    }
}

/// ChrBaseRegion.checkMergeOverlaps(regions, true): sorted by compareTo (stable), then regions of the
/// same chromosome that overlap or lie within 2 bases are merged, keeping the larger end
pub fn merge_chr_regions(regions: &mut Vec<ChrBaseRegion>) {
    regions.sort_by(compare_chr_regions);
    let mut i = 0;
    while i + 1 < regions.len() {
        if regions[i].chromosome == regions[i + 1].chromosome && regions[i].end >= regions[i + 1].start - 2 {
            let next = regions.remove(i + 1);
            if next.end > regions[i].end {
                regions[i].end = next.end;
            }
        } else {
            i += 1;
        }
    }
}

/// ChrBaseRegion.loadChrBaseRegions(filename): a BED file (name ending .bed or .bed.gz; start + 1) or a
/// TSV/CSV with a header naming Chromosome (or chromosome), PositionStart/PosStart and PositionEnd/PosEnd
/// (or lower case), or without a header (columns 1-3). The header is recognised by the text
/// "Chromosome" in the first line; the delimiter is a comma for a .csv or .csv.gz name, else a tab.
/// Each chromosome's regions go through BaseRegion.checkMergeOverlaps (see `merge_base_regions`). The map
/// iterates as Java's HashMap does.
pub fn load_chr_base_regions(filename: &str, keep_larger_end: bool) -> anyhow::Result<JavaHashMap<String, Vec<BaseRegion>>> {
    let is_bed = filename.ends_with(".bed") || filename.ends_with(".bed.gz");
    let delim = if filename.ends_with(".csv") || filename.ends_with(".csv.gz") { ',' } else { '\t' };
    let lines = java_text::read_lines(filename).with_context(|| format!("reading {filename}"))?;
    let mut map: JavaHashMap<String, Vec<BaseRegion>> = JavaHashMap::new();
    let Some(first) = lines.first() else {
        // Java: NullPointerException on the first readLine of an empty file
        bail!("empty regions file: {filename}");
    };
    let (mut chr_i, mut start_i, mut end_i) = (0usize, 1usize, 2usize);
    let mut rest = &lines[..];
    if first.contains("Chromosome") {
        let f = java_text::fields_index_map(first, delim);
        let get = |names: &[&str]| -> Option<usize> { names.iter().find_map(|n| f.get(*n).copied()) };
        // INVALID_FIELD (-1) in Java: the read of values[-1] throws
        chr_i = get(&["Chromosome", "chromosome"]).ok_or_else(|| anyhow!("{filename}: no chromosome column"))?;
        start_i = get(&["PositionStart", "positionstart", "PosStart", "posstart"]).ok_or_else(|| anyhow!("{filename}: no start column"))?;
        end_i = get(&["PositionEnd", "positionend", "PosEnd", "posend"]).ok_or_else(|| anyhow!("{filename}: no end column"))?;
        rest = &lines[1..];
    }
    for line in rest {
        let v = java_text::split_all(line, delim);
        let field = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("{filename}: missing column {} in line: {line}", i + 1));
        let chromosome = field(chr_i)?;
        let mut start = java_text::parse_int(field(start_i)?).ok_or_else(|| anyhow!("{filename}: invalid start in line: {line}"))?;
        if is_bed {
            start += 1;
        }
        let end = java_text::parse_int(field(end_i)?).ok_or_else(|| anyhow!("{filename}: invalid end in line: {line}"))?;
        let key = chromosome.to_string();
        if let Some(list) = map.get_mut(&key) {
            list.push(BaseRegion::new(start, end));
        } else {
            map.insert(key, vec![BaseRegion::new(start, end)]);
        }
    }
    let keys: Vec<String> = map.iter_unordered().map(|(k, _)| k.clone()).collect();
    for k in keys {
        merge_base_regions(map.get_mut(&k).unwrap(), keep_larger_end);
    }
    Ok(map)
}

/// BedFileReader.loadBedFile(lines) with BedLine.factory(): lines containing "Chromosome" are skipped;
/// each other line is tab-separated with at least 3 fields, start + 1 as BED has it, the chromosome
/// as HumanChromosome names it ("1" for "chr1": BedLine.toRegion takes HumanChromosome.toString); then
/// ChrBaseRegion.checkMergeOverlaps. A non-human contig makes Java add a null region and fail; that is
/// an error here.
pub fn load_bed_lines(lines: &[String]) -> anyhow::Result<Vec<ChrBaseRegion>> {
    let mut regions = Vec::new();
    for line in lines {
        if line.contains("Chromosome") {
            continue;
        }
        let d = java_text::split_all(line, '\t');
        if d.len() < 3 {
            bail!("invalid slice BED entry: {line}");
        }
        if !HumanChromosome::contains(d[0]) {
            bail!("non-human contig in BED entry: {line}");
        }
        let c = HumanChromosome::from_name(d[0]).ok_or_else(|| anyhow!("No enum constant HumanChromosome for chromosome '{}'", d[0]))?;
        let start = java_text::parse_int(d[1]).ok_or_else(|| anyhow!("invalid BED start: {line}"))? + 1;
        let end = java_text::parse_int(d[2]).ok_or_else(|| anyhow!("invalid BED end: {line}"))?;
        regions.push(ChrBaseRegion::new(c.short_name(), start, end));
    }
    merge_chr_regions(&mut regions);
    Ok(regions)
}

/// PartitionUtils.buildPartitions: [1, size-1], [size, 2*size-1], ... up to `max_position`
pub fn build_partitions(chromosome: &str, max_position: i32, partition_size: i32) -> Vec<ChrBaseRegion> {
    let mut partitions = Vec::new();
    let mut i: i64 = 0;
    loop {
        let start = (i * partition_size as i64) as i32;
        let end = (start as i64 + partition_size as i64 - 1).min(max_position as i64) as i32;
        partitions.push(ChrBaseRegion::new(chromosome, start.max(1), end));
        if end >= max_position {
            break;
        }
        i += 1;
    }
    partitions
}

/// PartitionUtils.partitionSpecificRegions: each region cut at multiples of the partition size
pub fn partition_specific_regions(regions: &[ChrBaseRegion], partition_size: i32) -> Vec<ChrBaseRegion> {
    let mut partitions = Vec::new();
    for r in regions {
        let (first, last) = (r.start / partition_size, r.end / partition_size);
        if first == last {
            partitions.push(r.clone());
        } else {
            for i in first..=last {
                let p_start = i * partition_size;
                let p_end = p_start + partition_size - 1;
                partitions.push(ChrBaseRegion::new(&r.chromosome, p_start.max(r.start).max(1), p_end.min(r.end)));
            }
        }
    }
    partitions
}

/// PartitionUtils.partitionChromosome(chromosome, version, specificRegions, partitionSize): the specific
/// regions on this chromosome cut into partitions, or the whole chromosome (its RefGenomeCoordinates
/// length) when there are none
pub fn partition_chromosome(
    chromosome: &str, version: RefGenomeVersion, specific_regions: &[ChrBaseRegion], partition_size: i32,
) -> Vec<ChrBaseRegion> {
    if !specific_regions.is_empty() {
        let on_chr: Vec<ChrBaseRegion> = specific_regions.iter().filter(|r| r.chromosome == chromosome).cloned().collect();
        return partition_specific_regions(&on_chr, partition_size);
    }
    let length = version.chromosome_length(strip_chr_prefix(chromosome));
    build_partitions(&version.versioned_chromosome(chromosome), length, partition_size)
}

impl SpecificRegions {
    /// SpecificRegions.from(configBuilder) as Java parses it: regions separated by ';' as chr:start-end or
    /// chr:start:end, each with start > 0 and end >= start, sorted by ChrBaseRegion.compareTo, their
    /// chromosomes in that order; else the chromosomes separated by ';'. An error where Java logs "invalid
    /// specific config" and returns null.
    pub fn from_java(regions: Option<&str>, chromosomes: Option<&str>) -> anyhow::Result<Self> {
        let mut s = SpecificRegions::default();
        if let Some(r) = regions.filter(|r| !r.is_empty()) {
            for item in java_text::split_all(r, ';') {
                let parts = java_text::split(item, ':');
                let bad = || anyhow!("invalid specific region: {item}");
                let (start, end) = if parts.len() == 3 {
                    (java_text::parse_int(parts[1]).ok_or_else(bad)?, java_text::parse_int(parts[2]).ok_or_else(bad)?)
                } else if parts.len() == 2 && parts[1].contains('-') {
                    let (a, b) = parts[1].split_once('-').unwrap();
                    (java_text::parse_int(a).ok_or_else(bad)?, java_text::parse_int(b).ok_or_else(bad)?)
                } else {
                    return Err(bad());
                };
                if !(start > 0 && end >= start) {
                    return Err(bad());
                }
                s.regions.push(ChrBaseRegion::new(parts[0], start, end));
            }
            s.regions.sort_by(compare_chr_regions);
            for r in &s.regions {
                if !s.chromosomes.contains(&r.chromosome) {
                    s.chromosomes.push(r.chromosome.clone());
                }
            }
        } else if let Some(c) = chromosomes.filter(|c| !c.is_empty()) {
            s.chromosomes = java_text::split(c, ';').into_iter().map(|x| x.to_string()).collect();
        }
        Ok(s)
    }
    /// excludeChromosome
    pub fn exclude_chromosome(&self, chr: &str) -> bool {
        !self.include_chromosome(chr)
    }
    /// includeRegion(posStart, posEnd): no regions, or one overlapping (the chromosome is not checked)
    pub fn include_positions(&self, start: i32, end: i32) -> bool {
        self.regions.is_empty() || self.regions.iter().any(|r| positions_overlap(r.start, r.end, start, end))
    }
}

#[cfg(test)]
mod java_region_tests {
    use super::*;

    #[test]
    fn versions_and_lengths() {
        assert_eq!(RefGenomeVersion::from_java("38").unwrap(), RefGenomeVersion::V38);
        assert!(RefGenomeVersion::from_java("hg38").is_err());
        assert_eq!(RefGenomeVersion::V38.versioned_chromosome("1"), "chr1");
        assert_eq!(RefGenomeVersion::V37.versioned_chromosome("chrX"), "X");
        assert_eq!(RefGenomeVersion::V38.chromosome_length("chr21"), 46709983);
        assert_eq!(RefGenomeVersion::V37.chromosome_length("Y"), 59373566);
        assert_eq!(RefGenomeVersion::V38.chromosome_length("chrM"), 0);
    }

    #[test]
    fn merging() {
        let mut r = vec![BaseRegion::new(100, 500), BaseRegion::new(200, 300), BaseRegion::new(502, 510), BaseRegion::new(600, 700)];
        let mut fixed = r.clone();
        merge_base_regions(&mut r, false);
        // Java: the contained region cuts the first short, then 502 is not within 2 bases of 300
        assert_eq!(r, vec![BaseRegion::new(100, 300), BaseRegion::new(502, 510), BaseRegion::new(600, 700)]);
        merge_base_regions(&mut fixed, true);
        assert_eq!(fixed, vec![BaseRegion::new(100, 510), BaseRegion::new(600, 700)]);
        let mut c = vec![ChrBaseRegion::new("chr2", 10, 20), ChrBaseRegion::new("chr1", 50, 60), ChrBaseRegion::new("chr1", 5, 100), ChrBaseRegion::new("chr1", 102, 103)];
        merge_chr_regions(&mut c);
        assert_eq!(c, vec![ChrBaseRegion::new("chr1", 5, 103), ChrBaseRegion::new("chr2", 10, 20)]);
    }

    #[test]
    fn partitions() {
        let p = build_partitions("chr1", 2_500_000, 1_000_000);
        assert_eq!(p, vec![ChrBaseRegion::new("chr1", 1, 999_999), ChrBaseRegion::new("chr1", 1_000_000, 1_999_999), ChrBaseRegion::new("chr1", 2_000_000, 2_500_000)]);
        assert_eq!(build_partitions("chr1", 999_999, 1_000_000).len(), 1);
        let s = partition_specific_regions(&[ChrBaseRegion::new("chr1", 999_000, 2_000_500), ChrBaseRegion::new("chr1", 5, 10)], 1_000_000);
        assert_eq!(
            s,
            vec![
                ChrBaseRegion::new("chr1", 999_000, 999_999),
                ChrBaseRegion::new("chr1", 1_000_000, 1_999_999),
                ChrBaseRegion::new("chr1", 2_000_000, 2_000_500),
                ChrBaseRegion::new("chr1", 5, 10)
            ]
        );
        assert_eq!(partition_chromosome("chr21", RefGenomeVersion::V38, &[], 10_000_000).last().unwrap().end, 46709983);
    }

    #[test]
    fn specific_regions_java() {
        let s = SpecificRegions::from_java(Some("chr2:5-10;chr1:3:4"), None).unwrap();
        assert_eq!(s.regions, vec![ChrBaseRegion::new("chr1", 3, 4), ChrBaseRegion::new("chr2", 5, 10)]);
        assert_eq!(s.chromosomes, vec!["chr1", "chr2"]);
        assert!(SpecificRegions::from_java(Some("chr1:0-10"), None).is_err());
        assert!(SpecificRegions::from_java(Some("chr1:10-5"), None).is_err());
        assert!(SpecificRegions::from_java(Some("chr1:1-5;"), None).is_err());
        assert!(SpecificRegions::from_java(Some("chr1"), None).is_err());
        let c = SpecificRegions::from_java(None, Some("chr1;chrX;")).unwrap();
        assert_eq!(c.chromosomes, vec!["chr1", "chrX"]);
        assert!(c.exclude_chromosome("chr2") && c.include_positions(1, 2));
    }

    #[test]
    fn region_files() {
        let dir = std::env::temp_dir().join(format!("hmf_common_regions_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bed = dir.join("t.bed");
        std::fs::write(&bed, "chr1\t99\t200\nchr1\t150\t180\nchr2\t9\t20\n").unwrap();
        let m = load_chr_base_regions(bed.to_str().unwrap(), false).unwrap();
        assert_eq!(m.get(&"chr1".to_string()).unwrap(), &vec![BaseRegion::new(100, 180)]);
        let m = load_chr_base_regions(bed.to_str().unwrap(), true).unwrap();
        assert_eq!(m.get(&"chr1".to_string()).unwrap(), &vec![BaseRegion::new(100, 200)]);
        let tsv = dir.join("t.tsv");
        std::fs::write(&tsv, "Chromosome\tPositionStart\tPositionEnd\nchr3\t10\t20\n").unwrap();
        let m = load_chr_base_regions(tsv.to_str().unwrap(), true).unwrap();
        assert_eq!(m.get(&"chr3".to_string()).unwrap(), &vec![BaseRegion::new(10, 20)]);
        let lines: Vec<String> = ["Chromosome\tStart\tEnd", "chr1\t0\t10000", "chr1\t10001\t10010", "chr2\t5\t6"].iter().map(|s| s.to_string()).collect();
        assert_eq!(load_bed_lines(&lines).unwrap(), vec![ChrBaseRegion::new("1", 1, 10010), ChrBaseRegion::new("2", 6, 6)]);
        std::fs::remove_dir_all(&dir).ok();
    }
}
