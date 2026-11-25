//! Read unmapping in problematic regions (Java ReadUnmapper and UnmappedRegions).

use crate::rec::*;
use anyhow::{bail, Result};
use rust_htslib::bam;
use std::collections::HashMap;

pub const UNMAP_MAX_NON_OVERLAPPING_BASES: i32 = 10;
pub const UNMAP_MIN_SOFT_CLIP: i32 = 20;
pub const UNMAP_MIN_HIGH_DEPTH: i32 = 1000;
pub const UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX: i64 = 1000;
const READ_END_APPROX_BUFFER: i32 = 100;

#[derive(Clone, Copy, Debug)]
pub struct UnmappingRegion {
    pub start: i32,
    pub end: i32,
    pub max_depth: i32,
}

#[derive(Default, Clone, Debug)]
pub struct UnmapStats {
    pub read_count: i64,
    pub mate_count: i64,
    pub supp_alignment_count: i64,
    pub supplementary_count: i64,
    pub secondary_count: i64,
    pub fully_unmapped_count: i64,
    pub existing_unmapped: i64,
    pub high_depth_count: i64,
    pub long_soft_clip_count: i64,
    pub chimeric_count: i64,
}

impl UnmapStats {
    pub fn merge(&mut self, o: &UnmapStats) {
        self.read_count += o.read_count;
        self.mate_count += o.mate_count;
        self.supp_alignment_count += o.supp_alignment_count;
        self.supplementary_count += o.supplementary_count;
        self.secondary_count += o.secondary_count;
        self.fully_unmapped_count += o.fully_unmapped_count;
        self.existing_unmapped += o.existing_unmapped;
        self.high_depth_count += o.high_depth_count;
        self.long_soft_clip_count += o.long_soft_clip_count;
        self.chimeric_count += o.chimeric_count;
    }
    pub fn paired_string(&self) -> String {
        format!(
            "reads({}) mates({}) both({} existing={}) supps(read={} align={}) seconds({}) reasons(depth={} softClip={} chimeric={})",
            self.read_count, self.mate_count, self.fully_unmapped_count, self.existing_unmapped, self.supplementary_count,
            self.supp_alignment_count, self.secondary_count, self.high_depth_count, self.long_soft_clip_count, self.chimeric_count
        )
    }
    pub fn unpaired_string(&self) -> String {
        format!(
            "reads({}) existingUnmapped({}) supps(read={} align={}) seconds({}) reasons(depth={} softClip={})",
            self.fully_unmapped_count, self.existing_unmapped, self.supplementary_count, self.supp_alignment_count,
            self.secondary_count, self.high_depth_count, self.long_soft_clip_count
        )
    }
}

/// UnmapRegionState: the regions of the current slice and the index of the last one matched.
pub struct UnmapRegionState {
    pub chromosome: String,
    pub slice_start: i32,
    pub slice_end: i32,
    pub regions: Vec<UnmappingRegion>,
    pub last_matched_index: Option<usize>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum RegionMatch {
    None,
    HighDepth,
    Other,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum UnmapReason {
    HighDepth,
    SoftClip,
    Chimeric,
    None,
}

/// One supplementary alignment from the SA attribute.
pub struct SuppData<'a> {
    pub chromosome: &'a [u8],
    pub position: i32,
    pub strand: u8,
    pub cigar: &'a [u8],
}

fn supp_from_alignment(alignment: &[u8]) -> Option<SuppData<'_>> {
    // Java String.split(",") drops trailing empty items; six fields are required
    let mut items: Vec<&[u8]> = alignment.split(|&c| c == b',').collect();
    while items.last().is_some_and(|i| i.is_empty()) {
        items.pop();
    }
    if items.len() != 6 {
        return None;
    }
    let position = std::str::from_utf8(items[1]).ok()?.parse().ok()?;
    Some(SuppData { chromosome: items[0], position, strand: *items[2].first()?, cigar: items[3] })
}

/// SupplementaryReadData.extractAlignments: None when absent or when any alignment is malformed
pub fn extract_supp_alignments(sa: Option<&[u8]>) -> Option<Vec<SuppData<'_>>> {
    let sa = sa?;
    if sa.is_empty() {
        return None;
    }
    let mut parts: Vec<&[u8]> = sa.split(|&c| c == b';').collect();
    if sa.last() == Some(&b';') {
        parts.pop();
    }
    let mut out = Vec::with_capacity(parts.len());
    for p in parts {
        out.push(supp_from_alignment(p)?);
    }
    Some(out)
}

/// SupplementaryReadData.extractAlignment: the first alignment
pub fn extract_supp_alignment(sa: Option<&[u8]>) -> Option<SuppData<'_>> {
    let sa = sa?;
    if sa.is_empty() {
        return None;
    }
    let first = sa.split(|&c| c == b';').next()?;
    supp_from_alignment(first)
}

pub struct ReadUnmapper {
    /// by reference index of the BAM header
    regions_by_tid: Vec<Vec<UnmappingRegion>>,
    by_name: HashMap<Vec<u8>, usize>,
    names: Vec<String>,
    enabled: bool,
}

/// Collections.binarySearch on region starts: the last region starting at or before `read_start`, else 0
pub fn binary_search(read_start: i32, regions: &[UnmappingRegion]) -> usize {
    match regions.binary_search_by(|r| r.start.cmp(&read_start)) {
        Ok(i) => i,
        Err(0) => 0,
        Err(i) => i - 1,
    }
}

#[inline]
fn positions_overlap(s1: i32, e1: i32, s2: i32, e2: i32) -> bool {
    !(s1 > e2 || e1 < s2)
}
#[inline]
fn positions_within(inner_s: i32, inner_e: i32, outer_s: i32, outer_e: i32) -> bool {
    inner_s <= inner_e && inner_s >= outer_s && inner_e <= outer_e
}

/// ReadUnmapper.overlapsUnmapRegion
pub fn overlaps_unmap_region(region: &UnmappingRegion, read_start: i32, read_end: i32) -> bool {
    if !positions_overlap(read_start, read_end, region.start, region.end) {
        return false;
    }
    if positions_within(read_start, read_end, region.start, region.end) {
        return true;
    }
    let overlap = region.end.min(read_end) - region.start.max(read_start) + 1;
    let non_overlapping = (read_end - read_start + 1) - overlap;
    non_overlapping < UNMAP_MAX_NON_OVERLAPPING_BASES
}

/// ReadUnmapper.getClipLengthFromCigarStr: soft and hard clips
pub fn clip_length_from_cigar_str(cigar: Option<&[u8]>) -> i32 {
    let Some(cigar) = cigar else { return 0 };
    let mut total = 0;
    let mut len: i32 = 0;
    for &c in cigar {
        if c == b'S' || c == b'H' {
            total += len;
            len = 0;
            continue;
        }
        let digit = c as i32 - b'0' as i32;
        if (0..=9).contains(&digit) {
            len = len * 10 + digit;
        } else {
            len = 0;
        }
    }
    total
}

fn read_end_from_cigar_str(read_start: i32, cigar: &[u8]) -> i32 {
    let mut current = read_start;
    let mut len: i32 = 0;
    for &c in cigar {
        if matches!(c, b'D' | b'M' | b'N') {
            current += len;
            len = 0;
            continue;
        }
        let digit = c as i32 - b'0' as i32;
        if (0..=9).contains(&digit) {
            len = len * 10 + digit;
        } else {
            len = 0;
        }
    }
    current - 1
}

fn exceeds_soft_clip_length(cigar: Option<&[u8]>) -> bool {
    clip_length_from_cigar_str(cigar) >= UNMAP_MIN_SOFT_CLIP
}

impl ReadUnmapper {
    pub fn disabled(header: &bam::HeaderView) -> Self {
        Self::build(header, HashMap::new())
    }

    pub(crate) fn build(header: &bam::HeaderView, mut by_chr: HashMap<String, Vec<UnmappingRegion>>) -> Self {
        let enabled = !by_chr.is_empty();
        let mut names = Vec::new();
        let mut by_name = HashMap::new();
        let mut regions_by_tid = Vec::new();
        for tid in 0..header.target_count() {
            let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
            regions_by_tid.push(by_chr.remove(&name).unwrap_or_default());
            by_name.insert(name.as_bytes().to_vec(), tid as usize);
            names.push(name);
        }
        ReadUnmapper { regions_by_tid, by_name, names, enabled }
    }

    /// UnmappedRegions.loadUnmapRegions: a TSV with Chromosome, PosStart, PosEnd, MaxDepth
    pub fn load(filename: &str, header: &bam::HeaderView) -> Result<Self> {
        let text = std::fs::read_to_string(filename)?;
        let delim = if filename.ends_with(".csv") { ',' } else { '\t' };
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap_or("").split(delim).collect();
        let col = |names: &[&str]| head.iter().position(|h| names.contains(h));
        let (Some(ci), Some(si), Some(ei), Some(di)) =
            (col(&["Chromosome", "chromosome"]), col(&["PosStart", "posStart"]), col(&["PosEnd", "posEnd"]), col(&["MaxDepth"]))
        else {
            bail!("unmap regions file {filename}: missing columns");
        };
        let mut by_chr: HashMap<String, Vec<UnmappingRegion>> = HashMap::new();
        let mut count = 0;
        for line in lines {
            let v: Vec<&str> = line.split(delim).collect();
            let region = UnmappingRegion { start: v[si].parse()?, end: v[ei].parse()?, max_depth: v[di].parse()? };
            let list = by_chr.entry(v[ci].to_string()).or_default();
            if let Some(last) = list.last() {
                if positions_overlap(last.start, last.end, region.start, region.end) {
                    bail!("unmap regions overlap: current({}-{}) next({}-{})", last.start, last.end, region.start, region.end);
                } else if region.end < last.start {
                    bail!("unmap regions not sorted: current({}-{}) next({}-{})", last.start, last.end, region.start, region.end);
                }
            }
            list.push(region);
            count += 1;
        }
        crate::info!("loaded {count} unmapping regions from {filename}");
        Ok(Self::build(header, by_chr))
    }

    pub fn enabled(&self) -> bool {
        self.enabled
    }
    pub fn regions(&self, tid: usize) -> &[UnmappingRegion] {
        &self.regions_by_tid[tid]
    }
    pub fn chromosome(&self, tid: usize) -> &str {
        &self.names[tid]
    }
    pub fn chromosome_count(&self) -> usize {
        self.names.len()
    }

    /// addMitochondrialRegion / addNonStandardContigs
    pub fn add_region(&mut self, tid: usize, region: UnmappingRegion) {
        self.regions_by_tid[tid].push(region);
        self.enabled = true;
    }

    fn regions_by_name(&self, chr: &[u8]) -> Option<&[UnmappingRegion]> {
        let tid = *self.by_name.get(chr)?;
        let r = &self.regions_by_tid[tid];
        // Java: a chromosome without regions has no map entry
        if r.is_empty() { None } else { Some(r) }
    }

    /// ReadUnmapper.checkTransformRead: true when the read or its mate or supplementary data was changed
    pub fn check_transform_read(&self, rec: &mut bam::Record, state: &mut UnmapRegionState, stats: &mut UnmapStats) -> bool {
        if !self.enabled {
            return false;
        }
        if has(rec, FLAG_SECONDARY) {
            self.check_secondary_read(rec, state, stats)
        } else if has(rec, FLAG_PAIRED) {
            self.check_transform_paired_read(rec, state, stats)
        } else {
            self.check_transform_unpaired_read(rec, state, stats)
        }
    }

    fn count_reason(reason: UnmapReason, stats: &mut UnmapStats) {
        match reason {
            UnmapReason::HighDepth => stats.high_depth_count += 1,
            UnmapReason::SoftClip => stats.long_soft_clip_count += 1,
            UnmapReason::Chimeric => stats.chimeric_count += 1,
            UnmapReason::None => {}
        }
    }

    fn check_transform_paired_read(&self, rec: &mut bam::Record, state: &mut UnmapRegionState, stats: &mut UnmapStats) -> bool {
        let read_unmapped = has(rec, FLAG_UNMAPPED);
        let mate_unmapped = has(rec, FLAG_MATE_UNMAPPED);
        if read_unmapped && mate_unmapped {
            return false;
        }
        let is_supplementary = has(rec, FLAG_SUPPLEMENTARY);
        let mut unmap_read = false;
        let mut unmap_mate = false;
        let mut read_region_type = RegionMatch::None;
        let mut unmap_reason = UnmapReason::None;
        if !read_unmapped {
            read_region_type = Self::find_max_depth_region_overlap_state(start(rec), end(rec), state, true);
            unmap_reason = self.check_unmap_read(rec, read_region_type);
            unmap_read = unmap_reason != UnmapReason::None;
        }
        let mut mate_region_type = RegionMatch::None;
        if !mate_unmapped {
            mate_region_type = self.mate_max_depth_region_overlap(rec, state);
            unmap_mate = if is_supplementary && !crate::compat::java("R43") {
                self.check_unmap_mate_of_supplementary(rec, mate_region_type)
            } else {
                Self::check_unmap_mate(rec, mate_region_type)
            };
        }
        // Java compares against a null region type for an unmapped read or mate; null != NONE is true there
        let read_type_not_none = read_unmapped || read_region_type != RegionMatch::None;
        let mate_type_not_none = mate_unmapped || mate_region_type != RegionMatch::None;
        if !read_unmapped && unmap_reason == UnmapReason::None && read_type_not_none && unmap_mate {
            unmap_reason = UnmapReason::Chimeric;
            unmap_read = true;
        }
        if !mate_unmapped && !unmap_mate && mate_type_not_none && unmap_read {
            unmap_mate = true;
        }
        if !read_unmapped {
            if unmap_reason != UnmapReason::None {
                Self::count_reason(unmap_reason, stats);
            } else if is_supplementary && !unmap_read {
                unmap_read = self.check_unmap_supplementary_read(rec);
            }
            if unmap_read && is_supplementary {
                set_unmapped_attributes(rec, self);
                stats.supplementary_count += 1;
                return true;
            }
        }
        if unmap_read {
            unmap_read_alignment(rec, mate_unmapped, unmap_mate, self);
        }
        if unmap_mate {
            unmap_mate_alignment(rec, read_unmapped, unmap_read);
        }
        if (read_unmapped || unmap_read) && (mate_unmapped || unmap_mate) {
            stats.fully_unmapped_count += 1;
        } else if unmap_read {
            stats.read_count += 1;
        } else if unmap_mate {
            stats.mate_count += 1;
        }
        if is_supplementary {
            return unmap_mate;
        }
        let mut unmap_supp_alignment = false;
        if !unmap_read && self.check_unmap_supplementary_alignments(rec) {
            set_attribute(rec, b"SA", None);
            stats.supp_alignment_count += 1;
            unmap_supp_alignment = true;
        }
        unmap_read || unmap_mate || unmap_supp_alignment
    }

    fn check_transform_unpaired_read(&self, rec: &mut bam::Record, state: &mut UnmapRegionState, stats: &mut UnmapStats) -> bool {
        if has(rec, FLAG_UNMAPPED) {
            return false;
        }
        let is_supplementary = has(rec, FLAG_SUPPLEMENTARY);
        let mut unmap_read = false;
        let match_type = Self::find_max_depth_region_overlap_state(start(rec), end(rec), state, true);
        let reason = self.check_unmap_read(rec, match_type);
        if reason != UnmapReason::None {
            unmap_read = true;
            // the Java switch counts depth and soft clip only
            if reason != UnmapReason::Chimeric {
                Self::count_reason(reason, stats);
            }
        } else if is_supplementary {
            unmap_read = self.check_unmap_supplementary_read(rec);
        }
        if unmap_read && is_supplementary {
            set_unmapped_attributes(rec, self);
            stats.supplementary_count += 1;
            return true;
        }
        if unmap_read {
            unmap_read_alignment(rec, false, true, self);
            stats.fully_unmapped_count += 1;
        }
        if is_supplementary {
            return false;
        }
        let mut unmap_supp_alignment = false;
        if !unmap_read && self.check_unmap_supplementary_alignments(rec) {
            set_attribute(rec, b"SA", None);
            stats.supp_alignment_count += 1;
            unmap_supp_alignment = true;
        }
        unmap_read || unmap_supp_alignment
    }

    fn check_secondary_read(&self, rec: &mut bam::Record, state: &mut UnmapRegionState, stats: &mut UnmapStats) -> bool {
        if has(rec, FLAG_UNMAPPED) {
            return false;
        }
        let match_type = Self::find_max_depth_region_overlap_state(start(rec), end(rec), state, true);
        let reason = self.check_unmap_read(rec, match_type);
        if reason == UnmapReason::None {
            return false;
        }
        if reason != UnmapReason::Chimeric {
            Self::count_reason(reason, stats);
        }
        set_flag(rec, FLAG_UNMAPPED, true);
        stats.secondary_count += 1;
        true
    }

    fn check_unmap_read(&self, rec: &bam::Record, match_type: RegionMatch) -> UnmapReason {
        match match_type {
            RegionMatch::None => return UnmapReason::None,
            RegionMatch::HighDepth => return UnmapReason::HighDepth,
            RegionMatch::Other => {}
        }
        let cigar = rec.raw_cigar();
        // clips at either end, soft or hard, as the cigar string parse counts them
        let clip: i32 = cigar.iter().filter(|c| matches!(op(**c), OP_S | OP_H)).map(|c| oplen(*c)).sum();
        if clip >= UNMAP_MIN_SOFT_CLIP {
            return UnmapReason::SoftClip;
        }
        if has(rec, FLAG_PAIRED) && self.is_chimeric_read(rec, false) {
            return UnmapReason::Chimeric;
        }
        UnmapReason::None
    }

    fn check_unmap_mate(rec: &bam::Record, match_type: RegionMatch) -> bool {
        match match_type {
            RegionMatch::None => return false,
            RegionMatch::HighDepth => return true,
            RegionMatch::Other => {}
        }
        if Self::is_chimeric(rec, true) {
            return true;
        }
        exceeds_soft_clip_length(aux_str(rec, b"MC"))
    }

    /// R43: Java decides the mate of a supplementary with the supplementary's own fragment length, chromosome and strand
    /// (isChimericRead(read, true)), so the supplementary can mark its mate unmapped while the mate, which checks
    /// itself against the primary, stays mapped (or the other way round), against "must be identical to how the mate
    /// checks itself". Corrected: the pair is taken as the mate sees it, from the primary's alignment (the first SA
    /// entry) and the mate's fields, with bwa's fragment length; Java's check where the SA tag or mate cigar is missing.
    fn check_unmap_mate_of_supplementary(&self, rec: &bam::Record, match_type: RegionMatch) -> bool {
        match match_type {
            RegionMatch::None => return false,
            RegionMatch::HighDepth => return true,
            RegionMatch::Other => {}
        }
        let chimeric = match self.pair_is_chimeric_for_mate(rec) {
            Some(c) => c,
            None => Self::is_chimeric(rec, true),
        };
        chimeric || exceeds_soft_clip_length(aux_str(rec, b"MC"))
    }

    /// isChimericRead(mate, false) for the mate of a supplementary's primary, or None without the data for it
    fn pair_is_chimeric_for_mate(&self, rec: &bam::Record) -> Option<bool> {
        if !has(rec, FLAG_PAIRED) || rec.mtid() < 0 {
            return None;
        }
        let primary = extract_supp_alignment(aux_str(rec, b"SA"))?;
        let mate_cigar = aux_str(rec, b"MC")?;
        let mate_chr = self.names.get(rec.mtid() as usize)?;
        if primary.chromosome != mate_chr.as_bytes() {
            return Some(true);
        }
        let primary_reverse = primary.strand == b'-';
        let mate_reverse = has(rec, FLAG_MATE_REVERSE);
        if primary_reverse == mate_reverse {
            return Some(true);
        }
        // bwa's fragment length: between the 5' ends of the aligned parts, plus one
        let five_prime = |pos: i32, cigar: &[u8], reverse: bool| -> i64 {
            if reverse { read_boundary_position(pos, cigar, false, false) as i64 } else { pos as i64 }
        };
        let p0 = five_prime(primary.position, primary.cigar, primary_reverse);
        let p1 = five_prime(mate_start(rec), mate_cigar, mate_reverse);
        let fragment_length = if p0 == p1 { 0 } else { (p0 - p1).abs() + 1 };
        Some(fragment_length > UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX)
    }

    fn is_chimeric_read(&self, rec: &bam::Record, check_for_mate: bool) -> bool {
        Self::is_chimeric(rec, check_for_mate)
    }

    /// ReadUnmapper.isChimericRead
    fn is_chimeric(rec: &bam::Record, check_for_mate: bool) -> bool {
        let read_unmapped = has(rec, FLAG_UNMAPPED);
        let mate_unmapped = has(rec, FLAG_MATE_UNMAPPED);
        if rec.insert_size().abs() > UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX {
            return true;
        }
        if !check_for_mate && mate_unmapped {
            return true;
        }
        if check_for_mate && read_unmapped {
            return true;
        }
        if has(rec, FLAG_PAIRED) {
            if rec.tid() != rec.mtid() {
                return true;
            }
            if has(rec, FLAG_REVERSE) == has(rec, FLAG_MATE_REVERSE) {
                return true;
            }
        }
        false
    }

    fn check_unmap_supplementary_read(&self, rec: &bam::Record) -> bool {
        let Some(alignments) = extract_supp_alignments(aux_str(rec, b"SA")) else { return false };
        for supp in &alignments {
            match self.supplementary_max_depth_region_overlap(supp) {
                RegionMatch::None => continue,
                RegionMatch::HighDepth => return true,
                RegionMatch::Other => {}
            }
            if exceeds_soft_clip_length(Some(supp.cigar)) {
                return true;
            }
            if self.is_supplementary_chimeric_read(rec, supp) {
                return true;
            }
        }
        false
    }

    fn is_supplementary_chimeric_read(&self, rec: &bam::Record, supp: &SuppData) -> bool {
        if has(rec, FLAG_PAIRED) {
            if has(rec, FLAG_MATE_UNMAPPED) {
                return true;
            }
            let mate_chr: &[u8] = if rec.mtid() >= 0 { self.names[rec.mtid() as usize].as_bytes() } else { b"*" };
            if supp.chromosome != mate_chr {
                return true;
            }
            let primary_negative = supp.strand == b'-';
            if primary_negative == has(rec, FLAG_MATE_REVERSE) {
                return true;
            }
        }
        false
    }

    fn check_unmap_supplementary_alignments(&self, rec: &bam::Record) -> bool {
        let Some(alignments) = extract_supp_alignments(aux_str(rec, b"SA")) else { return false };
        for supp in &alignments {
            match self.supplementary_max_depth_region_overlap(supp) {
                RegionMatch::None => continue,
                RegionMatch::HighDepth => return true,
                RegionMatch::Other => {}
            }
            if exceeds_soft_clip_length(Some(supp.cigar)) {
                return true;
            }
        }
        false
    }

    /// ReadUnmapper.mateInUnmapRegion
    pub fn mate_in_unmap_region(&self, rec: &bam::Record) -> bool {
        self.mate_max_depth_region_overlap_opt(rec, None) != RegionMatch::None
    }

    fn mate_max_depth_region_overlap(&self, rec: &bam::Record, state: &mut UnmapRegionState) -> RegionMatch {
        self.mate_max_depth_region_overlap_opt(rec, Some(state))
    }

    fn mate_max_depth_region_overlap_opt(&self, rec: &bam::Record, state: Option<&mut UnmapRegionState>) -> RegionMatch {
        let mate_start = mate_start(rec);
        let mate_chr: &str = if rec.mtid() >= 0 { &self.names[rec.mtid() as usize] } else { "*" };
        let mate_end = match aux_str(rec, b"MC") {
            Some(mc) => read_end_from_cigar_str(mate_start, mc),
            None => mate_start + rec.seq_len() as i32 + READ_END_APPROX_BUFFER,
        };
        if let Some(state) = state {
            if mate_chr == state.chromosome && mate_start >= state.slice_start && mate_start <= state.slice_end {
                return Self::find_max_depth_region_overlap_state(mate_start, mate_end, state, false);
            }
        }
        let Some(regions) = self.regions_by_name(mate_chr.as_bytes()) else { return RegionMatch::None };
        if !Self::is_within_region_range(mate_start, regions) {
            return RegionMatch::None;
        }
        Self::find_max_depth_region_overlap(mate_start, mate_end, regions, binary_search(mate_start, regions))
    }

    fn supplementary_max_depth_region_overlap(&self, supp: &SuppData) -> RegionMatch {
        let Some(regions) = self.regions_by_name(supp.chromosome) else { return RegionMatch::None };
        if !Self::is_within_region_range(supp.position, regions) {
            return RegionMatch::None;
        }
        let read_end = read_end_from_cigar_str(supp.position, supp.cigar);
        Self::find_max_depth_region_overlap(supp.position, read_end, regions, binary_search(supp.position, regions))
    }

    fn is_within_region_range(read_start: i32, regions: &[UnmappingRegion]) -> bool {
        let start_index = binary_search(read_start, regions);
        for region in &regions[start_index..] {
            if positions_overlap(read_start, read_start + READ_END_APPROX_BUFFER, region.start, region.end) {
                return true;
            }
            if region.start > read_start + READ_END_APPROX_BUFFER {
                break;
            }
        }
        false
    }

    /// checkRegionStateMatch: the index to start from when the last matched region still applies
    fn check_region_state_match(read_start: i32, read_end: i32, state: &UnmapRegionState) -> Option<usize> {
        let last = state.last_matched_index?;
        if state.regions.is_empty() {
            return None;
        }
        if !positions_within(read_start, read_end, state.slice_start, state.slice_end) {
            return None;
        }
        let region = &state.regions[last];
        if read_start < region.start {
            if last == 0 {
                return Some(0);
            }
            let prev = &state.regions[last - 1];
            return if read_start >= prev.start { Some(last - 1) } else { None };
        }
        if read_start >= region.start && read_start <= region.end {
            return Some(last);
        }
        if last + 1 >= state.regions.len() {
            return None;
        }
        let next = &state.regions[last + 1];
        if read_start < next.start { Some(last) } else { None }
    }

    /// findMaxDepthRegionOverlap against the slice's regions, with the cached index
    fn find_max_depth_region_overlap_state(read_start: i32, read_end: i32, state: &mut UnmapRegionState, update: bool) -> RegionMatch {
        if state.regions.is_empty() {
            return RegionMatch::None;
        }
        let start_index = match Self::check_region_state_match(read_start, read_end, state) {
            Some(i) => i,
            None => {
                let i = binary_search(read_start, &state.regions);
                if update {
                    state.last_matched_index = Some(i);
                }
                i
            }
        };
        Self::find_max_depth_region_overlap(read_start, read_end, &state.regions, start_index)
    }

    fn find_max_depth_region_overlap(read_start: i32, read_end: i32, regions: &[UnmappingRegion], start_index: usize) -> RegionMatch {
        let mut match_type = RegionMatch::None;
        for region in &regions[start_index..] {
            if region.start > read_end {
                break;
            }
            if !overlaps_unmap_region(region, read_start, read_end) {
                continue;
            }
            if positions_within(read_start, read_end, region.start, region.end) {
                return if region.max_depth >= UNMAP_MIN_HIGH_DEPTH { RegionMatch::HighDepth } else { RegionMatch::Other };
            }
            // overlapsUnmapRegion already established fewer than the maximum non-overlapping bases
            if region.max_depth >= UNMAP_MIN_HIGH_DEPTH {
                return RegionMatch::HighDepth;
            }
            match_type = RegionMatch::Other;
        }
        match_type
    }
}

/// setUnmappedAttributes: unmapped flag, map quality 0 and, for a primary read, the original location in UM
fn set_unmapped_attributes(rec: &mut bam::Record, unmapper: &ReadUnmapper) {
    let was_secondary_or_supp = has(rec, FLAG_SECONDARY | FLAG_SUPPLEMENTARY);
    set_flag(rec, FLAG_UNMAPPED, true);
    rec.set_mapq(0);
    if !was_secondary_or_supp {
        let chr = if rec.tid() >= 0 { unmapper.chromosome(rec.tid() as usize) } else { "*" };
        // contig names lose any ':' since it is the delimiter
        let value = format!("{}:{}", chr.replace(':', ""), start(rec));
        set_attribute(rec, b"UM", Some(&encode_z(value.as_bytes())));
    }
}

/// ReadUnmapper.unmapReadAlignment
pub fn unmap_read_alignment(rec: &mut bam::Record, mate_unmapped: bool, unmap_mate: bool, unmapper: &ReadUnmapper) {
    set_unmapped_attributes(rec, unmapper);
    set_flag(rec, FLAG_PROPER_PAIR, false);
    rec.set_insert_size(0);
    if mate_unmapped || unmap_mate || !has(rec, FLAG_PAIRED) {
        rec.set_pos(-1);
        rec.set_tid(-1);
    } else {
        rec.set_pos(rec.mpos());
        rec.set_tid(rec.mtid());
    }
    if mate_unmapped {
        rec.set_mpos(-1);
        rec.set_mtid(-1);
    }
    // setCigarString makes htsjdk write the attributes in their decoded form
    clear_cigar(rec);
    decode_attributes(rec);
    set_attribute(rec, b"SA", None);
    fix_bin(rec);
}

/// ReadUnmapper.unmapMateAlignment
pub fn unmap_mate_alignment(rec: &mut bam::Record, read_unmapped: bool, unmap_read: bool) {
    set_flag(rec, FLAG_MATE_UNMAPPED, true);
    set_flag(rec, FLAG_PROPER_PAIR, false);
    rec.set_insert_size(0);
    if read_unmapped || unmap_read {
        rec.set_mpos(-1);
        rec.set_mtid(-1);
        if read_unmapped {
            rec.set_pos(-1);
            rec.set_tid(-1);
            fix_bin(rec);
        }
    } else {
        rec.set_mpos(rec.pos());
        rec.set_mtid(rec.tid());
    }
    set_attribute(rec, b"MC", None);
}
