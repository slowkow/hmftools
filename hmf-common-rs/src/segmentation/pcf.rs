//! hmf-common `utils.pcf`: the .pcf file names of AMBER and COBALT, PURPLE's reader of segment boundaries
//! (`PCFFile.readPositions`, as PCFPosition), and the older reader and writer of the R copynumber layout
//! (`sampleID chrom arm start.pos end.pos n.probes mean`), which readPositions still accepts.

use super::HumanChromosome;
use anyhow::{anyhow, bail, Context, Result};
use std::collections::BTreeMap;
use std::io::Write;

const HEADER_PREFIX: &str = "sampleID";
pub const RATIO_EXTENSION: &str = ".cobalt.ratio.pcf";
pub const BAF_EXTENSION: &str = ".amber.baf.pcf";
const COL_CHROMOSOME: usize = 1;
const COL_POS_START: usize = 3;
const COL_POS_END: usize = 4;

/// PcfSegment: a region (1-based, inclusive) and the mean of the values in it.
#[derive(Clone, Debug, PartialEq)]
pub struct PcfSegment {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub mean_ratio: f64,
}

/// PCFSource
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum PcfSource {
    ReferenceRatio,
    TumorRatio,
    TumorBaf,
}

/// PCFPosition: a segment boundary and the range [min, max] in which PURPLE may place it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PcfPosition {
    pub source: PcfSource,
    pub chromosome: String,
    pub position: i32,
    pub min_position: i32,
    pub max_position: i32,
}

impl PcfPosition {
    pub fn new(source: PcfSource, chromosome: &str, position: i32) -> Self {
        PcfPosition { source, chromosome: chromosome.to_string(), position, min_position: position, max_position: position }
    }
    pub fn is_segment_end(&self) -> bool {
        self.position == self.min_position
    }
    pub fn is_segment_start(&self) -> bool {
        !self.is_segment_end()
    }
}

impl std::fmt::Display for PcfPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let source = match self.source {
            PcfSource::ReferenceRatio => "REFERENCE_RATIO",
            PcfSource::TumorRatio => "TUMOR_RATIO",
            PcfSource::TumorBaf => "TUMOR_BAF",
        };
        write!(
            f,
            "src({source}) pos({}:{}) range({} - {}) segment {}",
            self.chromosome,
            self.position,
            self.min_position,
            self.max_position,
            if self.position == self.min_position { "end" } else { "start" }
        )
    }
}

/// FileWriterUtils.checkAddDirSeparator
fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

/// PCFFile.generateRatioFilename: `<dir>/<sample>.cobalt.ratio.pcf`
pub fn generate_ratio_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{RATIO_EXTENSION}", check_add_dir_separator(base_path))
}

/// PCFFile.generateBAFFilename: `<dir>/<sample>.amber.baf.pcf`
pub fn generate_baf_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{BAF_EXTENSION}", check_add_dir_separator(base_path))
}

/// Integer.parseInt
fn parse_int(s: &str) -> Result<i32> {
    let ok = !s.is_empty() && s.strip_prefix(['+', '-']).unwrap_or(s).bytes().all(|b| b.is_ascii_digit()) && s.len() > usize::from(s.starts_with(['+', '-']));
    if !ok {
        bail!("For input string: \"{s}\"");
    }
    s.parse::<i32>().map_err(|_| anyhow!("For input string: \"{s}\""))
}

/// Double.parseDouble for decimal text: surrounding whitespace and an f/d suffix allowed, "NaN" and
/// "Infinity" spelled as Java spells them (Rust's "inf" and "nan" are refused).
pub fn parse_double(s: &str) -> Result<f64> {
    let t = s.trim_matches(|c: char| c <= ' ');
    let body = t.strip_prefix(['+', '-']).unwrap_or(t);
    let negative = t.starts_with('-');
    if body == "NaN" {
        return Ok(f64::NAN);
    }
    if body == "Infinity" {
        return Ok(if negative { f64::NEG_INFINITY } else { f64::INFINITY });
    }
    let numeric = body.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(body);
    let valid = !numeric.is_empty()
        && numeric.bytes().all(|b| b.is_ascii_digit() || matches!(b, b'.' | b'e' | b'E' | b'+' | b'-'))
        && numeric.bytes().any(|b| b.is_ascii_digit());
    let suffixed = numeric.len() != body.len();
    let v = if valid { t[..t.len() - usize::from(suffixed)].parse::<f64>().ok() } else { None };
    match v {
        // a float suffix makes Java parse as float and widen
        Some(x) if suffixed && body.ends_with(['f', 'F']) => Ok(x as f32 as f64),
        Some(x) => Ok(x),
        None => bail!("For input string: \"{s}\""),
    }
}

/// Window.start: the first base of the window holding `position`.
fn window_start(window_size: i32, position: i32) -> i32 {
    (position - 1) / window_size * window_size + 1
}

/// PCFFile.readPositions: each segment gives a start and an end boundary (end = segment end + 1); the gap
/// before a segment bounds its start from below, the next segment bounds the previous end from above, and
/// boundaries at the same position are merged. Reads the current layout (Chromosome Start End MeanRatio)
/// and the older R layout (header starting with "sampleID", positions as windows of `window_size`).
/// Rows of other contigs are skipped.
///
/// Java returns a HashMap keyed by HumanChromosome, iterated in identity-hash order (PURPLE's
/// ChimerismDetection loops over its keys); this map is in chromosome order. When a chromosome comes back
/// after another, Java starts a new list for it and drops the earlier one; so does this.
pub fn read_positions(window_size: i32, source: PcfSource, filename: &str) -> Result<BTreeMap<HumanChromosome, Vec<PcfPosition>>> {
    let text = std::fs::read_to_string(filename).with_context(|| format!("failed to read {filename}"))?;
    read_positions_from_text(window_size, source, &text)
}

pub fn read_positions_from_text(window_size: i32, source: PcfSource, text: &str) -> Result<BTreeMap<HumanChromosome, Vec<PcfPosition>>> {
    let lines: Vec<&str> = java_lines(text);
    let Some(header) = lines.first() else { bail!("Index 0 out of bounds for length 0") };
    let old_format = header.starts_with(HEADER_PREFIX);
    let (chr_index, start_index, end_index) = if old_format { (COL_CHROMOSOME, COL_POS_START, COL_POS_END) } else { (0, 1, 2) };

    let mut map: BTreeMap<HumanChromosome, Vec<PcfPosition>> = BTreeMap::new();
    let mut current_chromosome = String::new();
    let mut current: Option<HumanChromosome> = None;
    let mut min_position = 1;
    let mut have_last = false;
    for line in &lines[1..] {
        let values = java_split_tab(line);
        let chr = *values.get(chr_index).ok_or_else(|| anyhow!("Index {chr_index} out of bounds for length {}", values.len()))?;
        let bare = chr.strip_prefix("chr").unwrap_or(chr);
        if bare.is_empty() {
            // Java: HumanChromosome.contains parses the empty name as an integer
            bail!("For input string: \"\"");
        }
        if !HumanChromosome::contains(chr) {
            continue;
        }
        if chr != current_chromosome {
            current_chromosome = chr.to_string();
            let c = HumanChromosome::parse(chr)?;
            map.insert(c, Vec::new());
            current = Some(c);
            min_position = 1;
            have_last = false;
        }
        let field = |i: usize| values.get(i).copied().ok_or_else(|| anyhow!("Index {i} out of bounds for length {}", values.len()));
        let raw_start = parse_int(field(start_index)?)?;
        let raw_end = parse_int(field(end_index)?)?;
        let start = if old_format { window_start(window_size, raw_start) } else { raw_start };
        let end = if old_format { window_start(window_size, raw_end) + window_size } else { raw_end + 1 };
        let positions = map.get_mut(&current.unwrap()).unwrap();
        if have_last {
            positions.last_mut().unwrap().max_position = start;
        }
        let mut pcf_start = PcfPosition::new(source, chr, start);
        pcf_start.min_position = min_position;
        pcf_start.max_position = start;
        positions.push(pcf_start);
        min_position = end;
        positions.push(PcfPosition::new(source, chr, end));
        have_last = true;
    }
    for positions in map.values_mut() {
        merge_positions(positions);
    }
    Ok(map)
}

/// PCFFile.mergePositions: sorted by chromosome name and position (a stable sort), equal positions merged
/// into the first, neighbours bounding each other.
fn merge_positions(positions: &mut Vec<PcfPosition>) {
    positions.sort_by(|a, b| a.chromosome.cmp(&b.chromosome).then(a.position.cmp(&b.position)));
    let mut i = 0;
    while i + 1 < positions.len() {
        let (same_chr, same_pos) = (positions[i].chromosome == positions[i + 1].chromosome, positions[i].position == positions[i + 1].position);
        if same_chr && same_pos {
            let next = positions.remove(i + 1);
            let current = &mut positions[i];
            current.min_position = current.min_position.max(next.min_position);
            current.max_position = current.max_position.min(next.max_position);
        } else if same_chr {
            let next_position = positions[i + 1].position;
            let current_position = positions[i].position;
            positions[i].max_position = positions[i].max_position.min(next_position);
            positions[i + 1].min_position = positions[i + 1].min_position.max(current_position);
            i += 1;
        } else {
            i += 1;
        }
    }
}

/// PCFFile.readPcfFile: the segments of the current layout by chromosome, in file order. Java's
/// ListMultimap is keyed by HumanChromosome (it throws for other contigs).
pub fn read_pcf_file(path: &str) -> Result<BTreeMap<HumanChromosome, Vec<PcfSegment>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("failed to read {path}"))?;
    let mut result: BTreeMap<HumanChromosome, Vec<PcfSegment>> = BTreeMap::new();
    let mut lines = text.split('\n');
    lines.next(); // header
    for line in lines {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if line.is_empty() {
            continue; // BufferedReader gives no line after the final newline
        }
        let values: Vec<&str> = line.split('\t').collect();
        let field = |i: usize| values.get(i).copied().ok_or_else(|| anyhow!("Index {i} out of bounds for length {}", values.len()));
        let chr = field(0)?;
        let c = HumanChromosome::parse(chr)?;
        let segment = PcfSegment { chromosome: chr.to_string(), start: parse_int(field(1)?)?, end: parse_int(field(2)?)?, mean_ratio: parse_double(field(3)?)? };
        result.entry(c).or_default().push(segment);
    }
    Ok(result)
}

/// A region of the older layout (GenomeRegion): chromosome, start and end, 1-based and inclusive.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

/// PCFFile.read: the regions of the older layout, the end extended to the end of its window; lines starting
/// with "sampleID" are skipped. Java's ArrayListMultimap is keyed by the chromosome name.
pub fn read_regions(window_size: i32, filename: &str) -> Result<BTreeMap<String, Vec<Region>>> {
    let text = std::fs::read_to_string(filename).with_context(|| format!("failed to read {filename}"))?;
    let mut result: BTreeMap<String, Vec<Region>> = BTreeMap::new();
    for line in java_lines(&text) {
        if line.starts_with(HEADER_PREFIX) {
            continue;
        }
        let values = java_split_tab(line);
        let field = |i: usize| values.get(i).copied().ok_or_else(|| anyhow!("Index {i} out of bounds for length {}", values.len()));
        let chromosome = field(COL_CHROMOSOME)?.to_string();
        let start = parse_int(field(COL_POS_START)?)?;
        let end = parse_int(field(COL_POS_END)?)? + window_size - 1;
        result.entry(chromosome.clone()).or_default().push(Region { chromosome, start, end });
    }
    Ok(result)
}

/// PCFIntervals: the regions of one chromosome, each after the previous one.
#[derive(Clone, Debug)]
pub struct PcfIntervals {
    pub chromosome: HumanChromosome,
    pub intervals: Vec<Region>,
}

impl PcfIntervals {
    pub fn new(chromosome: HumanChromosome, intervals: Vec<Region>) -> Result<Self> {
        let mut previous_end: Option<i32> = None;
        for r in &intervals {
            if HumanChromosome::from_name(&r.chromosome) != Some(chromosome) {
                bail!("PCFIntervals: region {} is not on chromosome {}", r.chromosome, chromosome.short_name());
            }
            if previous_end.is_some_and(|e| r.start <= e) {
                bail!("PCFIntervals: regions overlap or are out of order");
            }
            previous_end = Some(r.end);
        }
        Ok(PcfIntervals { chromosome, intervals })
    }
}

/// GenomeIntervals: PCFIntervals of increasing chromosome.
#[derive(Clone, Debug)]
pub struct GenomeIntervals {
    pub intervals: Vec<PcfIntervals>,
}

impl GenomeIntervals {
    pub fn new(intervals: Vec<PcfIntervals>) -> Result<Self> {
        if intervals.windows(2).any(|w| w[0].chromosome >= w[1].chromosome) {
            bail!("GenomeIntervals: chromosomes out of order");
        }
        Ok(GenomeIntervals { intervals })
    }
    pub fn regions_list(&self) -> Vec<Region> {
        self.intervals.iter().flat_map(|i| i.intervals.iter().cloned()).collect()
    }
}

/// PCFFile.write: the regions in the older layout, the unknown fields "unused" and -1.
pub fn write_intervals(filename: &str, data: &GenomeIntervals) -> Result<()> {
    let file = std::fs::File::create(filename).with_context(|| format!("failed to create {filename}"))?;
    let mut w = std::io::BufWriter::new(file);
    w.write_all(b"sampleID\tchrom\tarm\tstart.pos\tend.pos\tn.probes\tmean\n")?;
    for r in data.regions_list() {
        writeln!(w, "unused\t{}\tunused\t{}\t{}\t-1\t-1", r.chromosome, r.start, r.end)?;
    }
    w.flush()?;
    Ok(())
}

/// Files.readAllLines: lines split at \n, \r or \r\n, none after a final terminator.
fn java_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let bytes = text.as_bytes();
    let mut start = 0;
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'\n' => {
                lines.push(&text[start..i]);
                start = i + 1;
            }
            b'\r' => {
                lines.push(&text[start..i]);
                if i + 1 < bytes.len() && bytes[i + 1] == b'\n' {
                    i += 1;
                }
                start = i + 1;
            }
            _ => {}
        }
        i += 1;
    }
    if start < bytes.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// String.split("\t"): trailing empty fields removed (all of them for a line of tabs only).
fn java_split_tab(line: &str) -> Vec<&str> {
    let mut v: Vec<&str> = line.split('\t').collect();
    if line.is_empty() {
        return v; // [""]
    }
    while v.last() == Some(&"") {
        v.pop();
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(c: &str, start: i32, end: i32) -> Region {
        Region { chromosome: format!("chr{c}"), start, end }
    }

    fn chr(c: &str) -> HumanChromosome {
        HumanChromosome::parse(c).unwrap()
    }

    fn assert_position(position: i32, min: i32, max: i32, p: &PcfPosition) {
        assert_eq!((position, min, max), (p.position, p.min_position, p.max_position), "{p}");
    }

    fn three_chromosomes() -> (GenomeIntervals, Vec<Region>, Vec<Region>) {
        let r1 = vec![region("1", 101, 200), region("1", 201, 300), region("1", 401, 500)];
        let r2 = vec![region("2", 1101, 1200), region("2", 1201, 1300), region("2", 1401, 1500)];
        let data = GenomeIntervals::new(vec![
            PcfIntervals::new(chr("1"), r1.clone()).unwrap(),
            PcfIntervals::new(chr("2"), r2.clone()).unwrap(),
            PcfIntervals::new(chr("3"), vec![]).unwrap(),
        ])
        .unwrap();
        (data, r1, r2)
    }

    // ---- GenomeIntervalsTest
    #[test]
    fn genome_intervals() {
        assert!(GenomeIntervals::new(vec![]).unwrap().regions_list().is_empty());
        let r1 = vec![region("1", 101, 200)];
        let one = GenomeIntervals::new(vec![PcfIntervals::new(chr("1"), r1.clone()).unwrap()]).unwrap();
        assert_eq!(r1, one.regions_list());
        let (data, r1, r2) = three_chromosomes();
        assert_eq!([r1, r2].concat(), data.regions_list());
        assert!(PcfIntervals::new(chr("1"), vec![region("1", 101, 200), region("1", 200, 300)]).is_err());
        assert!(PcfIntervals::new(chr("1"), vec![region("2", 101, 200)]).is_err());
        let a = PcfIntervals::new(chr("2"), vec![]).unwrap();
        let b = PcfIntervals::new(chr("1"), vec![]).unwrap();
        assert!(GenomeIntervals::new(vec![a, b]).is_err());
    }

    // ---- PCFFileTest
    #[test]
    fn write_test() {
        let (data, _, _) = three_chromosomes();
        let dir = std::env::temp_dir().join(format!("hmf_common_pcf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("whatever.pcf");
        let file = file.to_str().unwrap();
        write_intervals(file, &data).unwrap();
        let read = read_positions(100, PcfSource::TumorRatio, file).unwrap();
        assert_eq!(2, read.len());
        let p = &read[&chr("1")];
        assert_eq!(5, p.len());
        assert_position(101, 1, 101, &p[0]);
        assert_position(201, 201, 201, &p[1]);
        assert_position(301, 301, 401, &p[2]);
        assert_position(401, 301, 401, &p[3]);
        assert_position(501, 501, 501, &p[4]);
        let regions = read_regions(100, file).unwrap();
        assert_eq!(regions["chr1"], vec![region("1", 101, 299), region("1", 201, 399), region("1", 401, 599)]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn baf_file() {
        let Some(dir) = super::super::tests::resources_dir() else { return };
        let path = dir.join("pcf/baf.pcf");
        let map = read_positions(1000, PcfSource::TumorBaf, path.to_str().unwrap()).unwrap();
        for c in ["1", "2"] {
            let p = &map[&chr(c)];
            assert_eq!(4, p.len());
            let positions: Vec<i32> = p.iter().map(|x| x.position).collect();
            assert_eq!(positions, [93548001, 193800001, 193803001, 193804001]);
        }
    }

    #[test]
    fn ratio_file() {
        let Some(dir) = super::super::tests::resources_dir() else { return };
        let path = dir.join("pcf/ratio.pcf");
        let map = read_positions(1000, PcfSource::TumorBaf, path.to_str().unwrap()).unwrap();
        let p = &map[&chr("1")];
        assert_eq!(5, p.len());
        assert_position(835001, 1, 835001, &p[0]);
        assert_position(2583001, 2583001, 2583001, &p[1]);
        assert_position(2584001, 2584001, 2695001, &p[2]);
        assert_position(2695001, 2584001, 2695001, &p[3]);
        assert_position(4363001, 4363001, 4363001, &p[4]);
        let p = &map[&chr("3")];
        assert_eq!(5, p.len());
        assert_position(90449001, 1, 90449001, &p[0]);
        assert_position(90452001, 90452001, 90452001, &p[1]);
        assert_position(90455001, 90455001, 90457001, &p[2]);
        assert_position(90457001, 90455001, 90457001, &p[3]);
        assert_position(90458001, 90458001, 90458001, &p[4]);
    }

    /// The current layout, as SegmentsFile writes it
    #[test]
    fn current_layout() {
        let text = "Chromosome\tStart\tEnd\tMeanRatio\nchr1\t1000\t5000\t0.5\nchr1\t5100\t9000\t.75\nchrM\t1\t2\t0.5\nchr2\t100\t100\t1\n";
        let map = read_positions_from_text(1000, PcfSource::TumorBaf, text).unwrap();
        let p = &map[&chr("1")];
        assert_eq!(p.len(), 4);
        assert_position(1000, 1, 1000, &p[0]);
        assert_position(5001, 5001, 5100, &p[1]);
        assert_position(5100, 5001, 5100, &p[2]);
        assert_position(9001, 9001, 9001, &p[3]);
        assert!(p[1].is_segment_end() && p[2].is_segment_start());
        // a one-base segment: start and end boundaries 1 apart
        let p = &map[&chr("2")];
        assert_position(100, 1, 100, &p[0]);
        assert_position(101, 101, 101, &p[1]);
        assert!(read_positions_from_text(1000, PcfSource::TumorBaf, "").is_err());
        assert!(read_positions_from_text(1000, PcfSource::TumorBaf, "h\n\n").is_err());
    }

    #[test]
    fn segments_file_round_trip() {
        let dir = std::env::temp_dir().join(format!("hmf_common_pcf_rt_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.pcf");
        std::fs::write(&path, "Chromosome\tStart\tEnd\tMeanRatio\nchr1\t1\t10\t0.5432\nchr1\t11\t20\t1\nchrX\t5\t9\tNaN\n").unwrap();
        let m = read_pcf_file(path.to_str().unwrap()).unwrap();
        assert_eq!(m[&chr("1")].len(), 2);
        assert_eq!(m[&chr("1")][0], PcfSegment { chromosome: "chr1".into(), start: 1, end: 10, mean_ratio: 0.5432 });
        assert!(m[&chr("X")][0].mean_ratio.is_nan());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn names_and_numbers() {
        assert_eq!(generate_baf_filename("/out", "S"), "/out/S.amber.baf.pcf");
        assert_eq!(generate_ratio_filename("/out/", "S"), "/out/S.cobalt.ratio.pcf");
        assert_eq!(parse_double(" .5 ").unwrap(), 0.5);
        assert_eq!(parse_double("1e3d").unwrap(), 1000.0);
        assert_eq!(parse_double("0.1f").unwrap(), 0.1f32 as f64);
        assert_eq!(parse_double("-Infinity").unwrap(), f64::NEG_INFINITY);
        assert!(parse_double("inf").is_err() && parse_double("nan").is_err() && parse_double("").is_err() && parse_double("∞").is_err());
        assert_eq!(parse_int("+12").unwrap(), 12);
        assert!(parse_int("1.0").is_err() && parse_int("").is_err() && parse_int("-").is_err());
        assert_eq!(java_split_tab("a\tb\t\t"), ["a", "b"]);
        assert_eq!(java_lines("a\r\nb\rc\n"), ["a", "b", "c"]);
    }
}
