//! Reference data other than the PON: known hotspots (HotspotCache), target regions (TargetRegions)
//! and the repeat mask (RepeatMaskAnnotator with hmf-common's RepeatMaskAnnotations).

use super::types::*;
use super::variant::*;
use anyhow::{bail, Context, Result};
use hashbrown::HashMap;
use std::io::Read;

// ---------------------------------------------------------------------------------------------
// hotspots

#[derive(Clone, Debug)]
pub struct KnownHotspot {
    pub region_start: ChrBaseRegion,
    pub orient_start: Orient,
    pub region_end: ChrBaseRegion,
    pub orient_end: Orient,
}

impl KnownHotspot {
    pub fn matches(&self, chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, orient_start: Orient, orient_end: Orient) -> bool {
        if self.region_start.contains_position(chr_start, pos_start)
            && self.orient_start == orient_start
            && self.region_end.contains_position(chr_end, pos_end)
            && self.orient_end == orient_end
        {
            return true;
        }
        self.region_start.contains_position(chr_end, pos_end)
            && self.orient_start == orient_end
            && self.region_end.contains_position(chr_start, pos_start)
            && self.orient_end == orient_start
    }
}

#[derive(Default)]
pub struct HotspotCache {
    hotspots: Vec<KnownHotspot>,
    /// indices into `hotspots` by chromosome, of either end
    by_chromosome: HashMap<String, Vec<usize>>,
}

const POLY_A_HOMOLOGY: &str = "AAAAAAA";
const POLY_T_HOMOLOGY: &str = "TTTTTTT";

impl HotspotCache {
    /// HotspotCache.loadFile: a BEDPE file without a header. Returns the number of records.
    pub fn load(&mut self, filename: &str) -> Result<usize> {
        let text = std::fs::read_to_string(filename)?;
        // Java (C13): a line with too few fields is an uncaught ArrayIndexOutOfBoundsException. Corrected:
        // any line that cannot be read is reported with its number and what is wrong with it.
        let java_c13 = crate::common::compat::java("C13");
        let mut count = 0;
        // BufferedReader.readLine: a final line without a terminator is still a line
        for (index, line) in text.lines().enumerate() {
            let items: Vec<&str> = line.split('\t').collect();
            let parsed = (|| -> Result<KnownHotspot> {
                if items.len() < 10 {
                    if java_c13 {
                        bail!("ArrayIndexOutOfBoundsException: hotspot line with {} fields", items.len());
                    }
                    bail!("{} tab-separated fields where 10 are expected", items.len());
                }
                let int = |s: &str| {
                    super::vcf::java_parse_int(s).with_context(|| {
                        if java_c13 { format!("NumberFormatException: For input string: \"{s}\"") } else { format!("'{s}' is not a whole number") }
                    })
                };
                let orient = |s: &str| match s.chars().next() {
                    Some('+') => Ok(Orient::Forward),
                    Some('-') => Ok(Orient::Reverse),
                    other if java_c13 => bail!("Invalid orientation: {other:?}"),
                    _ => bail!("orientation '{s}' is neither + nor -"),
                };
                Ok(KnownHotspot {
                    region_start: ChrBaseRegion::new(items[0], int(items[1])? + 1, int(items[2])?),
                    orient_start: orient(items[8])?,
                    region_end: ChrBaseRegion::new(items[3], int(items[4])? + 1, int(items[5])?),
                    orient_end: orient(items[9])?,
                })
            })();
            match parsed {
                Ok(hotspot) => self.add_hotspot(hotspot),
                Err(e) if java_c13 => return Err(e),
                Err(e) => return Err(e.context(format!("line {}", index + 1))),
            }
            count += 1;
        }
        Ok(count)
    }

    pub fn add_hotspot(&mut self, hotspot: KnownHotspot) {
        let index = self.hotspots.len();
        self.by_chromosome.entry(hotspot.region_start.chromosome.clone()).or_default().push(index);
        if hotspot.region_start.chromosome != hotspot.region_end.chromosome {
            self.by_chromosome.entry(hotspot.region_end.chromosome.clone()).or_default().push(index);
        }
        self.hotspots.push(hotspot);
    }

    pub fn is_hotspot_variant(&self, var: &Variant) -> bool {
        if var.is_sgl() {
            return self.sgl_matches_hotspot(var);
        }
        let is_poly_at = |b: &Breakend| {
            let homology = b.context.attr_string(HOMSEQ, "");
            homology.contains(POLY_A_HOMOLOGY) || homology.contains(POLY_T_HOMOLOGY)
        };
        if is_poly_at(var.start()) || var.end().is_some_and(is_poly_at) {
            return false;
        }
        let Some(regions) = self.by_chromosome.get(var.chromosome_start()) else { return false };
        regions.iter().any(|&i| {
            self.hotspots[i].matches(var.chromosome_start(), var.chromosome_end(), var.pos_start(), var.pos_end(), var.orient_start(), var.orient_end())
        })
    }

    /// a single breakend is a hotspot when one of its insert sequence alignments completes a known pair
    fn sgl_matches_hotspot(&self, sgl: &Variant) -> bool {
        let alignments_str = sgl.start().context.attr_string(INSALN, "");
        if alignments_str.is_empty() {
            return false;
        }
        let alignments = AlternativeAlignment::from_vcf_tag(&alignments_str);
        let Some(regions) = self.by_chromosome.get(sgl.chromosome_start()) else { return false };
        let chr_start = sgl.chromosome_start();
        let pos_start = sgl.pos_start();
        let orient_start = sgl.orient_start();
        let java_c18 = crate::common::compat::java("C18");
        for alt in &alignments {
            // Java (C18): an INSALN item is "chr:start|strand|cigar|score", the leftmost aligned base and
            // the strand of the insert sequence, and Java takes them as the position and orientation of
            // the breakend at the other end of the insert. That breakend's orientation is the opposite
            // of the strand at a forward single breakend (the same at a reverse one), and when it faces
            // forward its position is the alignment's last base, as Linx reads the tag
            // (`SglMapping.from`). So at a forward single breakend Java matches a hotspot only through
            // inverted orientations, and the position can be off by the alignment's length.
            // Corrected: the partner breakend is formed as Linx forms it.
            let (alt_position, alt_orient) = if java_c18 {
                (alt.position, alt.orient)
            } else {
                let orient = if alt.orient == orient_start { Orient::Reverse } else { Orient::Forward };
                let position = if orient.is_forward() { alt.position + calc_cigar_aligned_length(&alt.cigar) - 1 } else { alt.position };
                (position, orient)
            };
            for &i in regions {
                let region = &self.hotspots[i];
                if region.matches(chr_start, &alt.chromosome, pos_start, alt_position, orient_start, alt_orient)
                    || region.matches(&alt.chromosome, chr_start, alt_position, pos_start, alt_orient, orient_start)
                {
                    return true;
                }
            }
        }
        false
    }
}

// ---------------------------------------------------------------------------------------------
// target regions

#[derive(Default)]
pub struct TargetRegions {
    regions: HashMap<String, Vec<(i32, i32)>>,
}

fn read_text(filename: &str) -> Result<String> {
    let mut text = String::new();
    super::vcf::open_text(filename)?.read_to_string(&mut text)?;
    Ok(text)
}

impl TargetRegions {
    /// TargetRegions.loadTargetRegionsBed with BedFileReader.loadBedFileChrMap. Returns the number of
    /// regions kept.
    pub fn load(&mut self, filename: &str, version: RefGenomeVersion) -> Result<usize> {
        let text = read_text(filename).with_context(|| format!("failed to load BED file({filename})"))?;
        let mut by_chromosome: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
        let mut current: Option<String> = None;
        let java_c11 = crate::common::compat::java("C11");
        for line in text.lines() {
            if line.contains("Chromosome") {
                continue;
            }
            let data: Vec<&str> = line.split('\t').collect();
            if data.len() < 3 {
                bail!("IllegalArgumentException: invalid slice BED entry: {line}");
            }
            if !is_human_chromosome(data[0]) {
                continue;
            }
            let chromosome = human_chromosome_name(data[0]);
            let int = |s: &str| super::vcf::java_parse_int(s).with_context(|| format!("NumberFormatException: For input string: \"{s}\""));
            let region = (int(data[1])? + 1, int(data[2])?);
            if current.as_deref() != Some(chromosome.as_str()) {
                if java_c11 {
                    // Java (C11): a chromosome that comes back after another one starts a new list, which
                    // replaces its first
                    by_chromosome.insert(chromosome.clone(), Vec::new());
                }
                current = Some(chromosome.clone());
            }
            by_chromosome.entry(chromosome).or_default().push(region);
        }
        for (chromosome, regions) in by_chromosome {
            self.regions.insert(version.versioned_chromosome(&chromosome), regions);
        }
        Ok(self.regions.values().map(|r| r.len()).sum())
    }

    #[cfg(test)]
    pub fn add_region(&mut self, chromosome: &str, start: i32, end: i32) {
        self.regions.entry(chromosome.to_string()).or_default().push((start, end));
    }

    pub fn has_target_regions(&self) -> bool {
        !self.regions.is_empty()
    }

    pub fn in_target_regions(&self, chromosome: &str, position: i32, buffer: i32) -> bool {
        match self.regions.get(chromosome) {
            Some(regions) => regions.iter().any(|&(start, end)| position >= start - buffer && position <= end + buffer),
            None => false,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// repeat mask

#[derive(Clone, Copy, Debug)]
struct RepeatMaskEntry {
    start: i32,
    end: i32,
    /// index into RepeatMask::names of the repeat and of its class
    repeat: u32,
    class_type: u32,
}

#[derive(Default)]
pub struct RepeatMask {
    by_chromosome: HashMap<String, Vec<RepeatMaskEntry>>,
    names: Vec<String>,
}

const MIN_OVERLAP: i32 = 20;
const MIN_COVERAGE_PERC: f64 = 0.1;
const POLY_A_T_PERC: f64 = 0.9;
const SMALL_LINEAR_SEARCH_SIZE: usize = 50;

/// Java String.trim then split on whitespace, keeping trailing empty fields out as the trim does
fn split_whitespace(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(|c| matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)).filter(|f| !f.is_empty())
}

fn parse_int_bytes(field: &[u8]) -> Option<i32> {
    std::str::from_utf8(field).ok().and_then(super::vcf::java_parse_int)
}

const PARSE_THREADS: usize = 4;
const BLOCK_SIZE: usize = 1 << 21;

enum BlockError {
    ShortLine,
    InvalidEntry(String),
}

/// The entries of a block of lines: runs of entries of one chromosome, with the names of repeats
/// and classes numbered within the block. Lines after a failing line are not parsed.
struct ParsedBlock {
    segments: Vec<(String, Vec<RepeatMaskEntry>)>,
    names: Vec<Vec<u8>>,
    error: Option<BlockError>,
}

fn parse_block(block: &[u8], version: RefGenomeVersion) -> ParsedBlock {
    let mut parsed = ParsedBlock { segments: Vec::new(), names: Vec::new(), error: None };
    let mut name_index: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut current_raw: Vec<u8> = Vec::new();
    let mut current_is_human = false;
    // a block ends with a line end, apart from the last one of a file without a final line end
    let lines = block.strip_suffix(b"\n").unwrap_or(block);

    for line in lines.split(|&c| c == b'\n') {
        let mut fields = split_whitespace(line);
        let sw_score = fields.next();
        let (Some(_), Some(_), Some(_), Some(_), Some(chr_raw)) = (sw_score, fields.next(), fields.next(), fields.next(), fields.next()) else {
            parsed.error = Some(BlockError::ShortLine);
            break;
        };
        if chr_raw != current_raw.as_slice() {
            // the name is only worked out when it changes; other contigs do not end a chromosome's run
            let name = String::from_utf8_lossy(chr_raw).into_owned();
            let chromosome = version.versioned_chromosome(&name);
            current_is_human = is_human_chromosome(&chromosome);
            if current_is_human && parsed.segments.last().is_none_or(|(c, _)| *c != chromosome) {
                parsed.segments.push((chromosome, Vec::new()));
            }
            current_raw.clear();
            current_raw.extend_from_slice(chr_raw);
        }
        if !current_is_human {
            continue;
        }
        // begin, end, (left), strand, repeat, class, begin, end, (left), id
        let mut rest: [Option<&[u8]>; 10] = [None; 10];
        for slot in rest.iter_mut() {
            *slot = fields.next();
        }
        let mut name_id = |name: &[u8]| -> u32 {
            if let Some(i) = name_index.get(name) {
                return *i;
            }
            let i = parsed.names.len() as u32;
            parsed.names.push(name.to_vec());
            name_index.insert(name.to_vec(), i);
            i
        };
        let entry = (|| {
            let start = parse_int_bytes(rest[0]?)?.wrapping_add(1);
            let end = parse_int_bytes(rest[1]?)?;
            parse_int_bytes(rest[9]?)?;
            parse_int_bytes(sw_score?)?;
            rest[3]?.first()?;
            let class_type = name_id(rest[5]?);
            let repeat = name_id(rest[4]?);
            Some(RepeatMaskEntry { start, end, repeat, class_type })
        })();
        match entry {
            Some(entry) => parsed.segments.last_mut().unwrap().1.push(entry),
            None => {
                parsed.error = Some(BlockError::InvalidEntry(String::from_utf8_lossy(line).into_owned()));
                break;
            }
        }
    }
    parsed
}

struct NameTable {
    index: HashMap<Vec<u8>, u32>,
    names: Vec<String>,
}

impl NameTable {
    fn id(&mut self, name: &[u8]) -> u32 {
        if let Some(i) = self.index.get(name) {
            return *i;
        }
        let i = self.names.len() as u32;
        self.names.push(String::from_utf8_lossy(name).into_owned());
        self.index.insert(name.to_vec(), i);
        i
    }
}

impl RepeatMask {
    pub fn has_data(&self) -> bool {
        !self.by_chromosome.is_empty()
    }

    pub fn entry_count(&self) -> usize {
        self.by_chromosome.values().map(|v| v.len()).sum()
    }

    /// RepeatMaskAnnotations.load: RepeatMasker output with three header lines. An entry that cannot
    /// be parsed is an error, reported with the line.
    ///
    /// The file is the largest input of the step, so it is decompressed on one thread and its lines
    /// are parsed in blocks on `PARSE_THREADS` others; the blocks are then joined in file order.
    ///
    /// `java_c11` is `compat::java("C11")`, read by the caller because this may run on another thread:
    /// Java lets a chromosome that comes back after another one replace its earlier entries. Otherwise
    /// all are kept, and the list of such a chromosome is put back in position order for the lookup.
    pub fn load(&mut self, filename: &str, version: RefGenomeVersion, java_c11: bool) -> Result<()> {
        // checked here so that a missing file is reported as such; the reader itself is opened on its thread
        drop(std::fs::File::open(filename)?);

        let (block_tx, block_rx) = std::sync::mpsc::sync_channel::<(usize, Vec<u8>)>(PARSE_THREADS * 2);
        let block_rx = std::sync::Arc::new(std::sync::Mutex::new(block_rx));
        let (parsed_tx, parsed_rx) = std::sync::mpsc::sync_channel::<(usize, ParsedBlock)>(PARSE_THREADS * 2);

        let read_error = std::thread::scope(|scope| -> Result<Option<std::io::Error>> {
            // blocks of whole lines, without the three header lines
            let reader_thread = scope.spawn(move || -> Option<std::io::Error> {
                let mut reader = match super::vcf::open_text(filename) {
                    Ok(r) => r,
                    Err(e) => return Some(e),
                };
                let mut carry: Vec<u8> = Vec::new();
                let mut header_lines = 3;
                let mut seq = 0;
                loop {
                    let mut block = std::mem::take(&mut carry);
                    let filled = block.len();
                    block.resize(filled + BLOCK_SIZE, 0);
                    let mut n = filled;
                    let mut eof = false;
                    while n < block.len() {
                        match reader.read(&mut block[n..]) {
                            Ok(0) => {
                                eof = true;
                                break;
                            }
                            Ok(read) => n += read,
                            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                            Err(e) => return Some(e),
                        }
                    }
                    block.truncate(n);
                    if !eof {
                        match memchr::memrchr(b'\n', &block) {
                            Some(i) => carry = block.split_off(i + 1),
                            None => {
                                // no line end yet: keep reading into the same block
                                carry = block;
                                continue;
                            }
                        }
                    }
                    let mut start = 0;
                    while header_lines > 0 && start < block.len() {
                        start = memchr::memchr(b'\n', &block[start..]).map(|i| start + i + 1).unwrap_or(block.len());
                        header_lines -= 1;
                    }
                    if start > 0 {
                        block.drain(..start);
                    }
                    if !block.is_empty() && block_tx.send((seq, block)).is_err() {
                        return None;
                    }
                    seq += 1;
                    if eof {
                        return None;
                    }
                }
            });
            for _ in 0..PARSE_THREADS {
                let block_rx = block_rx.clone();
                let parsed_tx = parsed_tx.clone();
                scope.spawn(move || loop {
                    let next = block_rx.lock().unwrap().recv();
                    let Ok((seq, block)) = next else { break };
                    if parsed_tx.send((seq, parse_block(&block, version))).is_err() {
                        break;
                    }
                });
            }
            drop(parsed_tx);
            drop(block_rx);

            // join the blocks in file order
            let mut names = NameTable { index: HashMap::new(), names: Vec::new() };
            let mut pending: HashMap<usize, ParsedBlock> = HashMap::new();
            let mut next_seq = 0;
            let mut current_chr = String::new();
            let mut returned: Vec<String> = Vec::new();
            let mut index = 0usize;
            let mut failure: Option<anyhow::Error> = None;
            for (seq, block) in parsed_rx.iter() {
                pending.insert(seq, block);
                while let Some(block) = pending.remove(&next_seq) {
                    next_seq += 1;
                    if failure.is_some() {
                        continue;
                    }
                    let remap: Vec<u32> = block.names.iter().map(|n| names.id(n)).collect();
                    for (chromosome, mut entries) in block.segments {
                        if chromosome != current_chr {
                            current_chr = chromosome;
                            if java_c11 {
                                // Java (C11): a chromosome that comes back later starts again from nothing
                                self.by_chromosome.insert(current_chr.clone(), Vec::new());
                            } else if self.by_chromosome.contains_key(&current_chr) {
                                if !returned.contains(&current_chr) {
                                    returned.push(current_chr.clone());
                                }
                            } else {
                                self.by_chromosome.insert(current_chr.clone(), Vec::new());
                            }
                        }
                        for entry in entries.iter_mut() {
                            entry.repeat = remap[entry.repeat as usize];
                            entry.class_type = remap[entry.class_type as usize];
                        }
                        index += entries.len();
                        self.by_chromosome.get_mut(&current_chr).unwrap().append(&mut entries);
                    }
                    match block.error {
                        Some(BlockError::ShortLine) => failure = Some(anyhow::anyhow!("ArrayIndexOutOfBoundsException: repeat-mask line has fewer than 5 fields")),
                        Some(BlockError::InvalidEntry(line)) => failure = Some(anyhow::anyhow!("invalid RM file entry: index({line}) line({index})")),
                        None => {}
                    }
                }
            }
            let read_error = reader_thread.join().expect("repeat mask reader thread");
            self.names = names.names;
            // the lookup bisects, so the runs of a chromosome that came back are merged by position
            for chromosome in &returned {
                if let Some(list) = self.by_chromosome.get_mut(chromosome) {
                    list.sort_by_key(|e| e.start);
                }
            }
            match failure {
                Some(e) => Err(e),
                None => Ok(read_error),
            }
        })?;
        if let Some(e) = read_error {
            return Err(e.into());
        }
        for list in self.by_chromosome.values_mut() {
            list.shrink_to_fit();
        }
        Ok(())
    }

    #[cfg(test)]
    pub fn add_entry(&mut self, chromosome: &str, start: i32, end: i32, repeat: &str, class_type: &str) {
        let mut id = |name: &str| match self.names.iter().position(|n| n == name) {
            Some(i) => i as u32,
            None => {
                self.names.push(name.to_string());
                self.names.len() as u32 - 1
            }
        };
        let entry = RepeatMaskEntry { start, end, repeat: id(repeat), class_type: id(class_type) };
        self.by_chromosome.entry(chromosome.to_string()).or_default().push(entry);
    }

    /// RepeatMaskAnnotations.findMatches: a bisection to one overlapping entry, then its neighbours
    fn find_matches(&self, chromosome: &str, region_start: i32, region_end: i32) -> Vec<RepeatMaskEntry> {
        let Some(regions) = self.by_chromosome.get(chromosome) else { return Vec::new() };
        let overlaps = |e: &RepeatMaskEntry| positions_overlap(e.start, e.end, region_start, region_end);
        if regions.len() < SMALL_LINEAR_SEARCH_SIZE {
            return regions.iter().filter(|e| overlaps(e)).copied().collect();
        }
        let mut current = regions.len() / 2;
        let mut lower = 0usize;
        let mut upper = regions.len() - 1;
        let mut matched = Vec::new();
        loop {
            let entry = &regions[current];
            if region_end < entry.start {
                if lower + 1 == current {
                    break;
                }
                upper = current;
                current = (lower + upper) / 2;
            } else if region_start > entry.end {
                if current + 1 == upper {
                    break;
                }
                lower = current;
                current = (lower + upper) / 2;
            } else if overlaps(entry) {
                matched.push(*entry);
                break;
            }
        }
        for search_up in [true, false] {
            let mut index = if search_up { current as i64 + 1 } else { current as i64 - 1 };
            while index >= 0 && (index as usize) < regions.len() {
                let entry = &regions[index as usize];
                if !overlaps(entry) {
                    break;
                }
                matched.push(*entry);
                index += if search_up { 1 } else { -1 };
            }
        }
        matched
    }

    /// RepeatMaskAnnotator.annotateVariants. Returns the number of variants annotated.
    pub fn annotate_variants(&self, variants: &mut [Variant]) -> usize {
        let mut annotated = 0;
        for var in variants.iter_mut() {
            if var.insert_sequence.is_empty() {
                continue;
            }
            let alignments = var.start().context.attr_string(INSALN, "");
            if alignments.is_empty() {
                continue;
            }
            if let Some(annotation) = self.annotate(&var.insert_sequence, &alignments) {
                var.rm_annotation = Some(annotation);
                annotated += 1;
            }
        }
        annotated
    }

    /// RepeatMaskAnnotator.annotate
    pub fn annotate(&self, insert_sequence: &str, alignments_str: &str) -> Option<RepeatMaskAnnotation> {
        let alignments = AlternativeAlignment::from_vcf_tag(alignments_str);
        if alignments.is_empty() {
            return None;
        }
        let ins_seq_length = insert_sequence.len() as f64;

        // first check for a poly-A or T
        for alignment in &alignments {
            let matched_bases = extract_matched_bases(insert_sequence, alignment);
            let at_count = matched_bases.bytes().filter(|&c| c == b'T' || c == b'A').count();
            let poly_at_perc = at_count as f64 / matched_bases.len() as f64;
            if poly_at_perc >= POLY_A_T_PERC {
                let coverage = matched_bases.len() as f64 / ins_seq_length;
                let a_count = matched_bases.bytes().filter(|&c| c == b'A').count();
                let repeat = if a_count > matched_bases.len() / 2 { "A(n)" } else { "T(n)" };
                return Some(RepeatMaskAnnotation { class_type: "Simple_repeat".to_string(), repeat: repeat.to_string(), coverage });
            }
        }

        let mut top: Option<RepeatMaskEntry> = None;
        let mut top_coverage = 0.0;
        for alignment in &alignments {
            let aligned_length = calc_cigar_aligned_length(&alignment.cigar);
            let alignment_start = alignment.position;
            let alignment_end = alignment.position + aligned_length - 1;
            for entry in self.find_matches(&alignment.chromosome, alignment_start, alignment_end) {
                let overlap = alignment_end.min(entry.end) - alignment_start.max(entry.start) + 1;
                if overlap < MIN_OVERLAP {
                    continue;
                }
                let coverage = overlap as f64 / ins_seq_length;
                if coverage < MIN_COVERAGE_PERC {
                    continue;
                }
                if coverage > top_coverage {
                    top_coverage = coverage;
                    top = Some(entry);
                }
            }
        }
        top.map(|entry| RepeatMaskAnnotation {
            class_type: self.names[entry.class_type as usize].clone(),
            repeat: self.names[entry.repeat as usize].clone(),
            coverage: top_coverage,
        })
    }
}

/// RepeatMaskAnnotator.extractMatchedBases: the insert sequence bases that the alignment matched
fn extract_matched_bases<'a>(insert_sequence: &'a str, alignment: &AlternativeAlignment) -> &'a str {
    let elements = cigar_elements(&alignment.cigar);
    let Some(first) = elements.first() else {
        panic!("IndexOutOfBoundsException: alignment without a cigar");
    };
    let match_start_pos = if first.1 == 'S' { first.0 } else { 0 };
    let match_bases: i32 = elements.iter().filter(|e| e.1 == 'M' || e.1 == 'I').map(|e| e.0).sum();
    let ins_seq_length = insert_sequence.len() as i32;
    let start_pos = match_start_pos.min(ins_seq_length);
    let end_pos = (match_start_pos + match_bases).min(ins_seq_length);
    if start_pos < 0 || start_pos > end_pos {
        panic!("StringIndexOutOfBoundsException: begin {start_pos}, end {end_pos}, length {ins_seq_length}");
    }
    &insert_sequence[start_pos as usize..end_pos as usize]
}
