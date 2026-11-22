//! The per-chromosome work: Java DepthTask with VariantInfo, RefSupportCounts, SliceRegionState and
//! ReadGroup.
//!
//! Variants of a chromosome are grouped by proximity, each group's region is read from every sample's
//! BAM, and each fragment is counted as reference support (a read over the breakend) or reference-pair
//! support (a read either side of it) for the variants it reaches. The control flow follows Java
//! statement by statement, including the parts that look unintended (see the notes on
//! `process_read` and `process_single_read`), because the counts depend on it once a variant reaches
//! its cap.

use super::config::DepthConfig;
use super::fmt::vcf_double;
use super::vcf::{parse_java_int, VariantContext};
use anyhow::{anyhow, bail, Result};
use hashbrown::HashMap;
use std::time::Instant;

pub const DEFAULT_MAX_FRAGMENT_LENGTH: i32 = 1100;
const READ_POSITION_MARGIN: i32 = 100;
const SHORT_INDEL_LENGTH: i32 = 1000;

pub const REF_DEPTH: &str = "REF";
pub const REF_DEPTH_PAIR: &str = "REFPAIR";
pub const ALLELE_FRACTION: &str = "AF";
pub const TOTAL_FRAGS: &str = "VF";
const CIPOS: &str = "CIPOS";
const SV_TYPE: &str = "SVTYPE";
const SV_TYPES: [&str; 7] = ["BND", "DEL", "DUP", "INF", "INS", "INV", "SGL"];

pub const FLAG_PAIRED: u16 = 0x1;
#[cfg(test)]
pub const FLAG_PROPER_PAIR: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
#[cfg(test)]
pub const FLAG_FIRST_OF_PAIR: u16 = 0x40;
#[cfg(test)]
pub const FLAG_SECOND_OF_PAIR: u16 = 0x80;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// What the depth step uses of an alignment record.
#[derive(Clone, Copy, Debug, Default)]
pub struct ReadInfo {
    /// 1-based alignment start
    pub start: i32,
    /// 1-based inclusive alignment end, as htsjdk's getAlignmentEnd (start - 1 when the cigar covers no
    /// reference bases)
    pub end: i32,
    pub flags: u16,
    pub ref_index: i32,
    pub mate_ref_index: i32,
    /// 1-based, 0 when not set
    pub mate_start: i32,
    pub insert_size: i32,
    /// number of bases stored in the record
    pub read_length: i32,
    pub left_soft_clipped: bool,
    pub right_soft_clipped: bool,
}

impl ReadInfo {
    fn paired(&self) -> bool {
        self.flags & FLAG_PAIRED != 0
    }
    fn mate_unmapped(&self) -> bool {
        self.flags & FLAG_MATE_UNMAPPED != 0
    }
    fn supplementary(&self) -> bool {
        self.flags & FLAG_SUPPLEMENTARY != 0
    }
    fn negative_strand(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }
    /// Java compares `read.getMateReferenceIndex()` and `read.getReferenceIndex()` with `==` / `!=`.
    /// Both are `Integer` objects, so the comparison is by identity: equal values are the same object
    /// only inside the Integer cache (-128 to 127). On a sequence with a higher index a mate on the same
    /// sequence therefore counts as being elsewhere.
    fn mate_on_same_reference(&self) -> bool {
        self.ref_index == self.mate_ref_index && (-128..=127).contains(&self.ref_index)
    }
    /// SamRecordUtils.inferredInsertSizeAbs
    fn inferred_insert_size_abs(&self) -> i32 {
        if self.paired() { self.insert_size.wrapping_abs() } else { (self.end - self.start + 1).wrapping_abs() }
    }
}

/// An alignment record as the task sees it.
pub trait SamRead {
    fn name(&self) -> &[u8];
    fn info(&self) -> ReadInfo;
    /// the SA attribute
    fn supplementary_data(&self) -> Option<&[u8]>;
}

/// String.split(delim) without a limit: trailing empty strings are dropped
fn java_split(s: &[u8], delim: u8) -> Vec<&[u8]> {
    let mut parts: Vec<&[u8]> = s.split(|b| *b == delim).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

fn java_int(s: &[u8]) -> Result<i32> {
    std::str::from_utf8(s)
        .ok()
        .and_then(parse_java_int)
        .ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{}\"", String::from_utf8_lossy(s)))
}

/// SupplementaryReadData.extractAlignment: the chromosome and position of the first alignment in an SA
/// attribute. None when the first alignment does not have six fields; an error where Java throws.
pub fn extract_supplementary_alignment(sa: &[u8]) -> Result<Option<(&[u8], i32)>> {
    if sa.is_empty() {
        return Ok(None);
    }
    let items = if sa.contains(&b',') {
        let first = sa.split(|b| *b == b';').next().unwrap_or(&[]);
        java_split(first, b',')
    } else {
        // the alternative delimiter, which Java only has for tests
        java_split(sa, b';')
    };
    if items.len() != 6 {
        return Ok(None);
    }
    let position = java_int(items[1])?;
    if items[2].is_empty() {
        bail!("java.lang.StringIndexOutOfBoundsException: empty strand in supplementary alignment");
    }
    java_int(items[4])?;
    java_int(items[5])?;
    Ok(Some((items[0], position)))
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Orientation {
    Forward,
    Reverse,
}

/// The groups of `^(.*)([\[\]])(.+)[\[\]](.*)$` (VariantAltInsertCoords.BREAKEND_REGEX) as the Java
/// matcher finds them: the first group is as long as it can be, then the third.
fn match_breakend(alt: &str) -> Option<(&str, &str, &str, &str)> {
    let brackets: Vec<usize> = alt.bytes().enumerate().filter(|(_, b)| *b == b'[' || *b == b']').map(|(i, _)| i).collect();
    let last = *brackets.last()?;
    for &first in brackets.iter().rev() {
        if last >= first + 2 {
            return Some((&alt[..first], &alt[first..first + 1], &alt[first + 1..last], &alt[last + 1..]));
        }
    }
    None
}

#[derive(Clone, Debug)]
pub struct RefSupportCounts {
    pub ref_support: i32,
    pub ref_pair_support: i32,
    pub vaf_cap: i32,
}

impl RefSupportCounts {
    fn new(vaf_cap: i32) -> Self {
        RefSupportCounts { ref_support: 0, ref_pair_support: 0, vaf_cap }
    }
    pub fn exceeds_max_depth(&self) -> bool {
        self.vaf_cap > 0 && self.ref_support.wrapping_add(self.ref_pair_support) >= self.vaf_cap
    }
}

#[derive(Clone, Debug)]
pub struct VariantInfo {
    pub position: i32,
    pub position_min: i32,
    pub position_max: i32,
    pub is_sgl: bool,
    pub is_short_indel: bool,
    pub orient: Orientation,
    pub sample_support_counts: Vec<RefSupportCounts>,
}

impl VariantInfo {
    /// `genotype_ids` are the VCF genotype columns of the configured samples, in sample order.
    pub fn new(variant: &VariantContext, genotype_ids: &[usize], vaf_cap: f64) -> Result<VariantInfo> {
        let position = variant.position;
        let Some(alt) = variant.alts.first() else {
            bail!("java.lang.IndexOutOfBoundsException: variant at {}:{} has no alternate allele", variant.chromosome, position);
        };
        // StructuralVariantFactory.isSingleBreakend, parseSingleOrientation and parseSvOrientation
        let is_sgl = alt.starts_with('.') || alt.ends_with('.');
        let breakend = if is_sgl { None } else { match_breakend(alt) };
        let orient = if is_sgl {
            if alt.starts_with('.') { Orientation::Reverse } else { Orientation::Forward }
        } else {
            match breakend {
                Some((before, _, _, _)) => if !before.is_empty() { Orientation::Forward } else { Orientation::Reverse },
                None => bail!("java.lang.IllegalArgumentException: Invalid orientation: 0"),
            }
        };

        let mut is_short_indel = false;
        if !is_sgl {
            if let Some(sv_type) = variant.info.get(SV_TYPE) {
                if !SV_TYPES.contains(&sv_type.as_str()) {
                    bail!("java.lang.IllegalArgumentException: No enum constant StructuralVariantType.{sv_type}");
                }
                if matches!(sv_type.as_str(), "DEL" | "DUP" | "INS") {
                    // VariantAltInsertCoords.fromRefAlt, for the position of the other breakend
                    let (before, _, location, after) = breakend.unwrap();
                    let anchor = if !before.is_empty() { before } else { after };
                    if !anchor.is_empty() && anchor.len() < variant.ref_.len() {
                        bail!("java.lang.StringIndexOutOfBoundsException: breakend bases shorter than the reference allele at {}:{}", variant.chromosome, position);
                    }
                    let chr_pos = java_split(location.as_bytes(), b':');
                    if chr_pos.len() < 2 {
                        bail!("java.lang.ArrayIndexOutOfBoundsException: Index 1 out of bounds for length {}", chr_pos.len());
                    }
                    let other_position = java_int(chr_pos[1])?;
                    let sv_length = position.wrapping_sub(other_position).wrapping_abs();
                    is_short_indel = sv_length <= SHORT_INDEL_LENGTH;
                }
            }
        }

        let mut homology = [0, 0];
        if let Some(values) = variant.info_ints(CIPOS)? {
            if values.len() < 2 {
                bail!("java.lang.IndexOutOfBoundsException: Index 1 out of bounds for length {}", values.len());
            }
            homology = [values[0], values[1]];
        }

        let mut sample_support_counts = Vec::with_capacity(genotype_ids.len());
        for &genotype_index in genotype_ids {
            let sample_fragments = match variant.genotypes.get(genotype_index) {
                Some(g) => g.attribute_as_int(TOTAL_FRAGS, 0)?,
                None => 0,
            };
            // (int) of a double: saturating, as the Rust cast
            let ref_frags_cap = if vaf_cap > 0.0 { (sample_fragments.max(1) as f64 / vaf_cap) as i32 } else { 0 };
            sample_support_counts.push(RefSupportCounts::new(ref_frags_cap));
        }

        Ok(VariantInfo {
            position,
            position_min: position.wrapping_add(homology[0]),
            position_max: position.wrapping_add(homology[1]),
            is_sgl,
            is_short_indel,
            orient,
            sample_support_counts,
        })
    }

    /// (ref, ref-pair) over the samples
    pub fn total_support(&self) -> (i32, i32) {
        let mut total = (0i32, 0i32);
        for c in &self.sample_support_counts {
            total.0 = total.0.wrapping_add(c.ref_support);
            total.1 = total.1.wrapping_add(c.ref_pair_support);
        }
        total
    }
}

/// The variants being sliced together. Variants are held as indexes into the task's list.
#[derive(Default)]
pub struct SliceRegionState {
    pub position_max: i32,
    pub position_min: i32,
    pub variants: Vec<usize>,
    pub uncapped_variants: Vec<usize>,
    pub min_position_index: usize,
}

impl SliceRegionState {
    pub fn reset(&mut self) {
        self.variants.clear();
        self.uncapped_variants.clear();
        self.position_min = 0;
        self.position_max = 0;
        self.min_position_index = 0;
    }
    pub fn reset_uncapped_variants(&mut self) {
        self.uncapped_variants.clear();
        self.uncapped_variants.extend_from_slice(&self.variants);
        self.min_position_index = 0;
    }
    pub fn add_variant(&mut self, index: usize, variant: &VariantInfo) {
        if self.variants.is_empty() {
            self.position_min = variant.position_min;
        } else {
            self.position_min = variant.position_min.min(self.position_min);
        }
        self.variants.push(index);
        self.position_max = variant.position_max.max(self.position_max);
    }
    pub fn variant_count(&self) -> usize {
        self.variants.len()
    }
}

pub struct ReadGroup {
    pub reads: Vec<ReadInfo>,
    pub wait_for_all: bool,
}

/// PerformanceCounter, as far as the step uses it
#[derive(Default)]
pub struct PerfCounter {
    started: Option<Instant>,
    pub last_time: f64,
    pub max_time: f64,
    pub total_time: f64,
    pub interval_count: i32,
}

impl PerfCounter {
    fn start(&mut self) {
        self.started = Some(Instant::now());
    }
    fn stop(&mut self) {
        let Some(started) = self.started.take() else { return };
        let seconds = started.elapsed().as_secs_f64();
        self.last_time = seconds;
        self.max_time = self.max_time.max(seconds);
        self.total_time += seconds;
        self.interval_count += 1;
    }
    /// one more interval of `seconds`
    pub fn add(&mut self, seconds: f64) {
        self.last_time = seconds;
        self.max_time = self.max_time.max(seconds);
        self.total_time += seconds;
        self.interval_count += 1;
    }
    pub fn merge(&mut self, other: &PerfCounter) {
        self.total_time += other.total_time;
        self.max_time = self.max_time.max(other.max_time);
        self.interval_count += other.interval_count;
    }
}

/// Where a task's reads come from: one source per sample, in sample order.
pub trait ReadSource {
    /// Passes the reads overlapping `start..=end` (1-based) of `chromosome` that pass the slicer's
    /// filters to `consumer`, in file order, until it returns false.
    fn slice(&mut self, chromosome: &str, start: i32, end: i32, consumer: &mut dyn FnMut(&dyn SamRead) -> Result<bool>) -> Result<()>;
}

pub struct DepthTask {
    chromosome: String,
    sample_ids: Vec<String>,
    trace: bool,
    /// VCF genotype column of each configured sample
    genotype_ids: Vec<usize>,
    vaf_cap: f64,
    proximity_distance: i32,
    perf_log_time: f64,
    ref_tag: String,
    ref_pair_tag: String,

    pub variants: Vec<VariantContext>,
    infos: Vec<VariantInfo>,
    read_groups: HashMap<Vec<u8>, ReadGroup>,
    state: SliceRegionState,
    current_sample_index: usize,
    total_read_count: i32,
    halt: bool,
    pub perf: PerfCounter,
}

impl DepthTask {
    pub fn new(chromosome: &str, config: &DepthConfig, genotype_ids: Vec<usize>) -> DepthTask {
        DepthTask {
            chromosome: chromosome.to_string(),
            sample_ids: config.sample_ids.clone(),
            trace: super::log_enabled(super::LEVEL_TRACE),
            genotype_ids,
            vaf_cap: config.vaf_cap,
            proximity_distance: config.proximity_distance,
            perf_log_time: config.perf_log_time,
            ref_tag: config.vcf_tag(REF_DEPTH),
            ref_pair_tag: config.vcf_tag(REF_DEPTH_PAIR),
            variants: Vec::new(),
            infos: Vec::new(),
            read_groups: HashMap::new(),
            state: SliceRegionState::default(),
            current_sample_index: 0,
            total_read_count: 0,
            halt: false,
            perf: PerfCounter::default(),
        }
    }

    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    pub fn add_variants(&mut self, variants: Vec<VariantContext>) -> Result<()> {
        for variant in variants {
            self.infos.push(VariantInfo::new(&variant, &self.genotype_ids, self.vaf_cap)?);
            self.variants.push(variant);
        }
        Ok(())
    }

    pub fn variant_count(&self) -> usize {
        self.infos.len()
    }

    /// The index ranges of the variants which are sliced together: the grouping loop of DepthTask.call.
    /// A variant joins the group while it is within the proximity distance of the one before it.
    pub fn group_ranges(&self) -> Vec<std::ops::Range<usize>> {
        let mut ranges = Vec::new();
        let mut index = 0usize;
        while index < self.infos.len() {
            let mut pos_end = self.infos[index].position;
            let mut next_index = index + 1;
            while next_index < self.infos.len() {
                let next_position = self.infos[next_index].position;
                if next_position.wrapping_sub(pos_end) > self.proximity_distance {
                    break;
                }
                pos_end = next_position;
                next_index += 1;
            }
            ranges.push(index..next_index);
            index = next_index;
        }
        ranges
    }

    /// The per-variant state the slicing works on, for handing out a group at a time.
    pub fn infos_mut(&mut self) -> &mut [VariantInfo] {
        &mut self.infos
    }

    /// Slices one group of a chromosome's variants from every sample and adds to the counts of `group`.
    /// Returns the number of reads gone through and the time taken in seconds. Groups do not share any
    /// state in Java either (the read cache and the slice state are reset for each), so they can be
    /// taken in any order and by any worker.
    pub fn slice_group(&mut self, chromosome: &str, group: &mut [VariantInfo], sources: &mut [Box<dyn ReadSource + '_>]) -> Result<(i32, f64)> {
        self.chromosome.clear();
        self.chromosome.push_str(chromosome);
        self.infos.clear();
        self.infos.extend_from_slice(group);
        self.state.reset();
        for (index, variant) in self.infos.iter().enumerate() {
            self.state.add_variant(index, variant);
        }
        self.total_read_count = 0;

        self.slice_sample_bams(sources)?;

        for (to, from) in group.iter_mut().zip(&self.infos) {
            to.sample_support_counts.clone_from(&from.sample_support_counts);
        }
        self.read_groups.clear();
        Ok((self.total_read_count, self.perf.last_time))
    }

    /// The counts go into each sample's genotype and, summed, into INFO.
    pub fn set_depth_values(&mut self) -> Result<()> {
        for (variant, info) in self.variants.iter_mut().zip(&self.infos) {
            // single breakends and short indels do not take ref-pair support into the allele fraction
            let skip_ref_pair_support = info.is_sgl || info.is_short_indel;

            for (counts, &genotype_index) in info.sample_support_counts.iter().zip(&self.genotype_ids) {
                let Some(genotype) = variant.genotypes.get_mut(genotype_index) else {
                    bail!("java.lang.IndexOutOfBoundsException: no genotype {genotype_index} at {}:{}", variant.chromosome, variant.position);
                };
                if genotype.attributes.is_empty() {
                    continue;
                }
                genotype.attributes.insert(self.ref_tag.clone(), counts.ref_support.to_string());
                let variant_frags = genotype.attribute_as_int(TOTAL_FRAGS, 0)?;
                let ref_pair_support = if skip_ref_pair_support { 0 } else { counts.ref_pair_support };
                let total = variant_frags.wrapping_add(counts.ref_support).wrapping_add(ref_pair_support) as f64;
                genotype.attributes.insert(self.ref_pair_tag.clone(), ref_pair_support.to_string());
                let af = if total > 0.0 { variant_frags as f64 / total } else { 0.0 };
                genotype.attributes.insert(ALLELE_FRACTION.to_string(), vcf_double(af));
            }

            let (total_ref, total_ref_pair) = info.total_support();
            variant.info.insert(self.ref_tag.clone(), total_ref.to_string());
            variant.info.insert(self.ref_pair_tag.clone(), (if skip_ref_pair_support { 0 } else { total_ref_pair }).to_string());
        }
        Ok(())
    }

    fn slice_sample_bams(&mut self, sources: &mut [Box<dyn ReadSource + '_>]) -> Result<()> {
        self.perf.start();

        // reads may start and end before the first variant, for the fragments which span the variants
        let region_start = self.state.position_min.wrapping_sub(DEFAULT_MAX_FRAGMENT_LENGTH);
        let region_end = self.state.position_max.wrapping_add(DEFAULT_MAX_FRAGMENT_LENGTH);

        let mut times: Vec<f64> = Vec::new();
        let mut read_counts: Vec<i32> = Vec::new();
        let chromosome = self.chromosome.clone();

        for (i, source) in sources.iter_mut().enumerate() {
            self.current_sample_index = i;
            self.read_groups.clear();
            self.state.reset_uncapped_variants();

            let start_time = Instant::now();
            let read_count = self.total_read_count;
            self.halt = false;
            if self.trace {
                trace!("sample({}) slice for {} variants", self.sample_ids[i], self.state.variant_count());
            }
            source.slice(&chromosome, region_start, region_end, &mut |read| {
                self.process_read(read)?;
                Ok(!self.halt)
            })?;
            times.push(start_time.elapsed().as_secs_f64());
            read_counts.push(self.total_read_count.wrapping_sub(read_count));

            // the fragments still waiting for a read; the order they are taken in does not change the counts
            let groups = std::mem::take(&mut self.read_groups);
            for (id, group) in &groups {
                self.process_read_group(&group.reads, id);
            }
            self.read_groups = groups;
        }

        self.perf.stop();

        if self.perf_log_time > 0.0 && self.perf.last_time > self.perf_log_time {
            let times: Vec<String> = times.iter().map(|x| super::fmt::java_fixed(*x, 3)).collect();
            let counts: Vec<String> = read_counts.iter().map(|x| x.to_string()).collect();
            debug!(
                "chr({}) slice(variants({}) uncapped({}) range({} - {}) minPosIndex({})) high depth retrieval time({}) totalFrags(0) times({}) readCounts({})",
                self.chromosome,
                self.state.variants.len(),
                self.state.uncapped_variants.len(),
                self.state.position_min,
                self.state.position_max,
                self.state.min_position_index,
                super::fmt::java_fixed(self.perf.last_time, 3),
                times.join(","),
                counts.join(",")
            );
        }
        Ok(())
    }

    /// Whether the read's SA attribute points into the current slice on this chromosome.
    fn supplementary_in_slice(read: &dyn SamRead, chromosome: &str, position_min: i32, max_slice_position: i32) -> Result<bool> {
        let Some(sa) = read.supplementary_data() else { return Ok(false) };
        Ok(match extract_supplementary_alignment(sa)? {
            Some((supp_chromosome, supp_position)) => {
                supp_chromosome == chromosome.as_bytes() && supp_position >= position_min && supp_position <= max_slice_position
            }
            None => false,
        })
    }

    /// DepthTask.processRead
    pub fn process_read(&mut self, read: &dyn SamRead) -> Result<()> {
        self.total_read_count = self.total_read_count.wrapping_add(1);

        let max_slice_position = self.state.position_max.wrapping_add(DEFAULT_MAX_FRAGMENT_LENGTH);
        let info = read.info();

        if let Some(group) = self.read_groups.get_mut(read.name()) {
            group.reads.push(info);

            if !group.wait_for_all {
                if info.supplementary() {
                    group.wait_for_all = true;
                } else if Self::supplementary_in_slice(read, &self.chromosome, self.state.position_min, max_slice_position)? {
                    group.wait_for_all = true;
                }
            }

            if !group.wait_for_all {
                let group = self.read_groups.remove(read.name()).unwrap();
                self.process_read_group(&group.reads, read.name());
            }
            return Ok(());
        }

        // Look for an uncapped variant this read is relevant to. Capped variants are dropped from the list
        // only as this scan meets them, so one which lies after the first relevant variant stays in the
        // list and keeps being counted.
        let mut is_relevant = false;
        let mut any_uncapped = false;
        let sample = self.current_sample_index;

        let mut index = 0;
        while index < self.state.uncapped_variants.len() {
            let variant = &self.infos[self.state.uncapped_variants[index]];

            if variant.sample_support_counts[sample].exceeds_max_depth() {
                self.state.uncapped_variants.remove(index);
                self.state.min_position_index = self.state.min_position_index.saturating_sub(1);
                continue;
            }

            any_uncapped = true;

            if is_relevant_read(&info, variant) {
                is_relevant = true;
                break;
            }
            index += 1;
        }

        // stop slicing this region for this sample once every variant has reached its cap
        if !any_uncapped {
            debug!(
                "chr({}) variants({}) range({} - {}) sampleIndex({}) VAF cap exceeded for all variants",
                self.chromosome,
                self.state.variant_count(),
                self.state.position_min,
                self.state.position_max,
                sample
            );
            self.halt = true;
            return Ok(());
        }

        if !is_relevant {
            return Ok(());
        }

        // whether the mate or a supplementary is expected within this slice
        let mut expect_supplementaries = info.supplementary();
        let expect_mate = info.paired() && !info.mate_unmapped() && info.mate_on_same_reference() && info.mate_start <= max_slice_position;

        if Self::supplementary_in_slice(read, &self.chromosome, self.state.position_min, max_slice_position)? {
            expect_supplementaries = true;
        }

        if !expect_mate && !expect_supplementaries {
            self.process_single_read(&info, read.name());
        } else {
            self.read_groups.insert(read.name().to_vec(), ReadGroup { reads: vec![info], wait_for_all: expect_supplementaries });
        }
        Ok(())
    }

    /// DepthTask.processSingleRead. With ten or more uncapped variants the scan starts at the index
    /// remembered from the previous single read rather than at the start of the list.
    fn process_single_read(&mut self, read: &ReadInfo, name: &[u8]) {
        if self.state.uncapped_variants.len() < 10 {
            self.process_read_group(std::slice::from_ref(read), name);
            return;
        }
        let trace = if self.trace { Some((self.chromosome.as_str(), name)) } else { None };

        let mut new_min_position_index = self.state.min_position_index;
        let sample = self.current_sample_index;

        for i in self.state.min_position_index..self.state.uncapped_variants.len() {
            let variant = &mut self.infos[self.state.uncapped_variants[i]];
            if read.start > variant.position_max.wrapping_add(READ_POSITION_MARGIN) {
                new_min_position_index = i + 1;
                continue;
            } else if read.end < variant.position_min.wrapping_sub(READ_POSITION_MARGIN) {
                break;
            }
            check_read_group_support(variant, sample, std::slice::from_ref(read), trace);
        }

        self.state.min_position_index = new_min_position_index;
    }

    /// DepthTask.processReadGroup
    fn process_read_group(&mut self, reads: &[ReadInfo], id: &[u8]) {
        let trace = if self.trace { Some((self.chromosome.as_str(), id)) } else { None };
        let group_min_position = reads.iter().map(|x| x.start).min().unwrap_or(0);
        let group_max_position = reads.iter().map(|x| x.end).max().unwrap_or(0);
        let sample = self.current_sample_index;

        for &variant_index in &self.state.uncapped_variants {
            let variant = &mut self.infos[variant_index];
            if group_min_position > variant.position_max.wrapping_add(READ_POSITION_MARGIN) {
                continue;
            } else if group_max_position < variant.position_min.wrapping_sub(READ_POSITION_MARGIN) {
                break;
            }
            check_read_group_support(variant, sample, reads, trace);
        }
    }

    // the accessors and entry points the Java unit tests use

    #[cfg(test)]
    pub fn reset(&mut self) {
        self.infos.clear();
        self.variants.clear();
        self.read_groups.clear();
        self.current_sample_index = 0;
        self.total_read_count = 0;
        self.state.reset();
    }

    #[cfg(test)]
    pub fn read_group_count(&self) -> usize {
        self.read_groups.len()
    }

    #[cfg(test)]
    pub fn variant_infos(&self) -> &[VariantInfo] {
        &self.infos
    }

    /// addSliceVariants(variantInfos()): all of the task's variants become the current slice
    #[cfg(test)]
    pub fn add_slice_variants(&mut self) {
        self.state.reset();
        for (index, variant) in self.infos.iter().enumerate() {
            self.state.add_variant(index, variant);
        }
        self.state.reset_uncapped_variants();
    }

    #[cfg(test)]
    pub fn slice_region_state(&self) -> &SliceRegionState {
        &self.state
    }

}

/// DepthTask.isRelevantRead: false for a read which cannot span the variant
fn is_relevant_read(read: &ReadInfo, variant: &VariantInfo) -> bool {
    if read.end < variant.position_min {
        if !read.paired() || read.mate_unmapped() || !read.mate_on_same_reference() {
            return false;
        }
        // both reads before the variant
        if read.mate_start.wrapping_add(read.read_length) < variant.position_min {
            return false;
        }
    } else if read.start > variant.position_max {
        return false;
    }
    true
}

/// DepthTask.checkReadGroupSupport: counts the fragment for the variant unless one of its reads is
/// soft-clipped at the breakend. Counting does not stop at the cap. `trace` is the chromosome and the
/// fragment's read name when trace logging is on.
fn check_read_group_support(variant: &mut VariantInfo, sample: usize, reads: &[ReadInfo], trace: Option<(&str, &[u8])>) {
    let mut read_supports_ref = false;
    let mut has_lower_pos_read = false;
    let mut has_upper_pos_read = false;
    let mut strand_count = 0;

    for read in reads {
        let read_start = read.start;
        let read_end = read.end;
        let is_supplementary = read.supplementary();

        // a soft clip exactly at the breakend
        if (variant.orient == Orientation::Reverse && read_start >= variant.position_min && read_start <= variant.position_max && read.left_soft_clipped)
            || (variant.orient == Orientation::Forward && read_end >= variant.position_min && read_end <= variant.position_max && read.right_soft_clipped)
        {
            if let Some((_, id)) = trace {
                trace!(
                    "var({}) pos({}-{}) read({}-{}) id({}) at junction",
                    variant.position,
                    variant.position_min,
                    variant.position_max,
                    read_start,
                    read_end,
                    String::from_utf8_lossy(id)
                );
            }
            return;
        }

        if !is_supplementary && reads.len() > 1 {
            let forward = !read.negative_strand();
            let fragment_size = read.inferred_insert_size_abs();

            if forward && read_end <= variant.position.max(variant.position_max) && !has_lower_pos_read && fragment_size < DEFAULT_MAX_FRAGMENT_LENGTH {
                has_lower_pos_read = true;
                strand_count += 1;
            } else if !forward && read_start >= variant.position.min(variant.position_min) && !has_upper_pos_read && fragment_size < DEFAULT_MAX_FRAGMENT_LENGTH {
                has_upper_pos_read = true;
                strand_count -= 1;
            }
        }

        if !(variant.position_min > read_end || variant.position_max < read_start) {
            read_supports_ref = true;
        }
    }

    let counts = &mut variant.sample_support_counts[sample];
    if read_supports_ref {
        counts.ref_support = counts.ref_support.wrapping_add(1);
    } else if has_lower_pos_read && has_upper_pos_read && strand_count == 0 {
        counts.ref_pair_support = counts.ref_pair_support.wrapping_add(1);
        if let Some((_, id)) = trace {
            let primaries = || reads.iter().filter(|x| !x.supplementary());
            trace!(
                "var({}) pos({}-{}) fragment(id={} {}-{}) has ref-pair support",
                variant.position,
                variant.position_min,
                variant.position_max,
                String::from_utf8_lossy(id),
                primaries().map(|x| x.start).min().unwrap_or(0),
                primaries().map(|x| x.end).max().unwrap_or(0).max(0)
            );
        }
    }

    if let Some((chromosome, _)) = trace {
        let counts = &variant.sample_support_counts[sample];
        if counts.exceeds_max_depth() {
            trace!(
                "var({}:{}) sampleIndex({}) ref limit({}) reached with support(ref={} pair={})",
                chromosome,
                variant.position,
                sample,
                counts.vaf_cap,
                counts.ref_support,
                counts.ref_pair_support
            );
        }
    }
}

#[cfg(test)]
mod helper_tests {
    use super::*;

    #[test]
    fn breakend_groups() {
        assert_eq!(match_breakend("A[2:500["), Some(("A", "[", "2:500", "")));
        assert_eq!(match_breakend("]2:500]ACG"), Some(("", "]", "2:500", "ACG")));
        assert_eq!(match_breakend("AC]x]y]G"), Some(("AC]x", "]", "y", "G")));
        assert_eq!(match_breakend("A[]"), None);
        assert_eq!(match_breakend("A[x[]"), Some(("A", "[", "x[", "")));
        assert_eq!(match_breakend("ACGT"), None);
    }

    #[test]
    fn supplementary_alignment() {
        assert_eq!(extract_supplementary_alignment(b"2,3000,+,94S57M,5,0;").unwrap(), Some((&b"2"[..], 3000)));
        assert_eq!(extract_supplementary_alignment(b"7,100,-,23S32M,0,6;11,200,+,115S32M,0,0;").unwrap(), Some((&b"7"[..], 100)));
        assert_eq!(extract_supplementary_alignment(b"2,3000,+,94S57M,5;").unwrap(), None);
        assert_eq!(extract_supplementary_alignment(b"2;3000;+;94S57M;5;0").unwrap(), Some((&b"2"[..], 3000)));
        assert_eq!(extract_supplementary_alignment(b"").unwrap(), None);
        assert!(extract_supplementary_alignment(b"2,x,+,94S57M,5,0").is_err());
    }

    #[test]
    fn integer_identity() {
        let mut read = ReadInfo { ref_index: 20, mate_ref_index: 20, ..Default::default() };
        assert!(read.mate_on_same_reference());
        read.mate_ref_index = 21;
        assert!(!read.mate_on_same_reference());
        read.ref_index = 200;
        read.mate_ref_index = 200;
        assert!(!read.mate_on_same_reference());
    }
}
