//! BamReader: the reads of one partition (region) counted into the partition's and the worker's
//! statistics.
//!
//! As in Java, each read is counted in the partition where it starts, and a partition's depth covers only
//! its own bases. BT1: Java drops the bases a read has beyond the end of its partition, so the first read
//! length of every partition is under-counted. Corrected: the reads that start in an earlier partition
//! and reach into this one add their bases here (`carry_over`), with the read-level counts left to the
//! partition they start in.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use anyhow::{anyhow, Result};
use hmf_common::regions::{BaseRegion, ChrBaseRegion, RefGenomeVersion};
use hmf_common::segmentation::HumanChromosome;

use crate::compat;
use crate::config::MetricsConfig;
use crate::coverage::{BaseCoverage, CoverageMetrics, CoverageParams};
use crate::genes::{self, GeneCoverage, GENE_COVERAGE_MIN_MAP_QUALITY};
use crate::record::{ConsensusType, Rec, OP_D, OP_EQ, OP_M, OP_N, OP_X};
use crate::stats::{FlagStats, FragmentLengths, PartitionStats, ReadCounts};
use crate::targets::{OffTargetAnalyser, TargetRegionStats};

/// FxHash, for the read names of the read groups
#[derive(Default)]
pub struct FxHasher {
    hash: u64,
}

const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut h = self.hash;
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            h = (h.rotate_left(5) ^ u64::from_le_bytes(c.try_into().unwrap())).wrapping_mul(SEED);
        }
        for &b in chunks.remainder() {
            h = (h.rotate_left(5) ^ b as u64).wrapping_mul(SEED);
        }
        self.hash = h;
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

type FxMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// ReadGroup: the aligned coordinates of a fragment's reads seen so far, kept until the read starting
/// at `max_read_start` (the mate or a supplementary within the read) arrives
struct ReadGroup {
    max_read_start: i32,
    #[allow(dead_code)]
    is_consensus: bool,
    coords: Vec<(i32, i32)>,
}

/// A partition to read
#[derive(Clone, Debug)]
pub struct Task {
    pub id: usize,
    pub region: ChrBaseRegion,
    /// with -only_target: the region before it on the same chromosome (BT18)
    pub previous: Option<(i32, i32)>,
    /// the reads without a position, for the flag counts' total (Java reads them after the partitions)
    pub unmapped: bool,
}

/// The rows of the per-partition files for one partition: (task id, partition stats row, target region
/// rows, enriched region rows)
pub type PartitionRows = (usize, Option<String>, Vec<String>, Vec<String>);

/// What a worker collects over its partitions (CombinedStats, and the rows of the per-partition files)
pub struct WorkerStats {
    pub coverage: CoverageMetrics,
    pub fragment_lengths: FragmentLengths,
    pub read_counts: ReadCounts,
    pub flag_stats: FlagStats,
    pub off_target_overlaps: HashMap<i32, i64>,
    /// per exon (GeneCoverage index), empty until a read reaches it
    pub exon_depths: Vec<Vec<u32>>,
    pub rows: Vec<PartitionRows>,
    pub reads: u64,
}

impl WorkerStats {
    pub fn new(config: &MetricsConfig, exon_count: usize) -> Self {
        WorkerStats {
            coverage: CoverageMetrics::new(config.max_coverage),
            fragment_lengths: FragmentLengths::default(),
            read_counts: ReadCounts::default(),
            flag_stats: FlagStats::default(),
            off_target_overlaps: HashMap::new(),
            exon_depths: vec![Vec::new(); exon_count],
            rows: Vec::new(),
            reads: 0,
        }
    }

    pub fn merge(&mut self, o: WorkerStats) {
        self.coverage.merge(&o.coverage);
        self.fragment_lengths.merge(&o.fragment_lengths);
        self.read_counts.merge(&o.read_counts);
        self.flag_stats.merge(&o.flag_stats);
        for (k, v) in o.off_target_overlaps {
            *self.off_target_overlaps.entry(k).or_insert(0) += v;
        }
        for (a, b) in self.exon_depths.iter_mut().zip(o.exon_depths) {
            if b.is_empty() {
                continue;
            }
            if a.is_empty() {
                *a = b;
            } else {
                for (x, y) in a.iter_mut().zip(b) {
                    *x += y;
                }
            }
        }
        self.rows.extend(o.rows);
        self.reads += o.reads;
    }
}

/// ExcludedRegions.getPolyGRegion
fn poly_g_region(version: RefGenomeVersion) -> ChrBaseRegion {
    if version.is_37() { ChrBaseRegion::new("2", 33141260, 33141700) } else { ChrBaseRegion::new("chr2", 32916190, 32916630) }
}

/// ChrBaseRegion.overlaps: the chromosomes compared as HumanChromosome (Java throws for a non-human
/// contig: BT16), then the positions
fn human_overlaps(a: &ChrBaseRegion, b: &ChrBaseRegion) -> Result<bool> {
    let (ca, cb) = (HumanChromosome::from_name(&a.chromosome), HumanChromosome::from_name(&b.chromosome));
    if ca.is_none() || cb.is_none() {
        if compat::java("BT16") {
            let bad = if ca.is_none() { &a.chromosome } else { &b.chromosome };
            return Err(anyhow!("java.lang.IllegalArgumentException: No enum constant HumanChromosome._{}", bad.strip_prefix("chr").unwrap_or(bad)));
        }
        return Ok(a.chromosome == b.chromosome && a.start <= b.end && a.end >= b.start);
    }
    Ok(ca == cb && a.start <= b.end && a.end >= b.start)
}

/// The exons a partition's reads are matched against: indices into the chromosome's list, sorted by
/// start, and the running maximum of their ends (so that exons ending before a read are skipped)
struct ExonCursor {
    idx: Vec<usize>,
    max_end: Vec<i32>,
}

/// The per-partition state of BamReader
pub struct PartitionReader<'c> {
    config: &'c MetricsConfig,
    genes: &'c GeneCoverage,
    coverage: BaseCoverage,
    groups: FxMap<Vec<u8>, ReadGroup>,
    spare_coords: Vec<Vec<(i32, i32)>>,
    // the current partition
    tid: i32,
    chromosome: String,
    start: i32,
    end: i32,
    previous: Option<(i32, i32)>,
    targets: Vec<TargetRegionStats>,
    record_off_target: bool,
    off_target: Option<OffTargetAnalyser>,
    excluded: Option<(i32, i32)>,
    exons: ExonCursor,
    read_counts: ReadCounts,
    flag_stats: FlagStats,
    fragment_lengths: FragmentLengths,
    partition_stats: PartitionStats,
    last_purge_check: i64,
    log_read_ids: bool,
    // compat switches, read once
    java_bt1: bool,
    java_bt5: bool,
    java_bt18: bool,
}

const FRAGMENT_MAP_SIZE_READ_COUNT_CHECK: i64 = 10000;
const FRAGMENT_MAP_SIZE_THRESHOLD: usize = 10000;

impl<'c> PartitionReader<'c> {
    pub fn new(config: &'c MetricsConfig, genes: &'c GeneCoverage) -> Self {
        let params = CoverageParams {
            map_qual_threshold: config.map_quality_threshold,
            base_qual_threshold: config.base_quality_threshold,
            max_coverage: config.max_coverage,
            expect_duplicates: config.expect_duplicates(),
        };
        PartitionReader {
            config,
            genes,
            coverage: BaseCoverage::new(params),
            groups: FxMap::default(),
            spare_coords: Vec::new(),
            tid: -1,
            chromosome: String::new(),
            start: 0,
            end: 0,
            previous: None,
            targets: Vec::new(),
            record_off_target: !config.target_regions.is_empty(),
            off_target: None,
            excluded: None,
            exons: ExonCursor { idx: Vec::new(), max_end: Vec::new() },
            read_counts: ReadCounts::default(),
            flag_stats: FlagStats::default(),
            fragment_lengths: FragmentLengths::default(),
            partition_stats: PartitionStats::default(),
            last_purge_check: 0,
            log_read_ids: !config.log_read_ids.is_empty(),
            java_bt1: compat::java("BT1"),
            java_bt5: compat::java("BT5"),
            java_bt18: compat::java("BT18"),
        }
    }

    /// new BamReader(region, ...): the partition's target regions, unmappable regions, excluded region and exons
    pub fn begin(&mut self, task: &Task, tid: i32) -> Result<()> {
        let config = self.config;
        let region = &task.region;
        self.tid = tid;
        self.chromosome.clone_from(&region.chromosome);
        self.start = region.start;
        self.end = region.end;
        self.previous = task.previous;
        self.targets.clear();
        if let Some(list) = config.target_regions_of(&region.chromosome) {
            for t in list {
                // BaseRegion.overlaps(ChrBaseRegion): positions only
                if t.start <= region.end && t.end >= region.start {
                    self.targets.push(TargetRegionStats { chromosome: region.chromosome.clone(), start: t.start, end: t.end, fragment_counts: ReadCounts::default() });
                }
            }
        }
        self.off_target = if !config.target_regions.is_empty() && config.write_off_target {
            Some(OffTargetAnalyser::new(&region.chromosome, config.high_fragment_overlap_threshold))
        } else {
            None
        };
        let mut unmappable = Vec::new();
        for u in &config.unmappable_regions {
            if human_overlaps(u, region)? {
                unmappable.push(BaseRegion::new(u.start, u.end));
            }
        }
        let poly_g = poly_g_region(config.ref_gen_version);
        self.excluded = if human_overlaps(&poly_g, region)? { Some((poly_g.start, poly_g.end)) } else { None };
        self.coverage.reset(region.start, region.end, unmappable);
        self.read_counts = ReadCounts::default();
        self.flag_stats = FlagStats::default();
        self.fragment_lengths = FragmentLengths::default();
        self.partition_stats = PartitionStats::default();
        self.last_purge_check = 0;
        self.groups.clear();
        // BT2: Java matches reads only against the exons overlapping the partition, so a read reaching
        // into an exon that starts after the partition's end is not counted there. Corrected: every exon
        // of the chromosome.
        self.exons.idx.clear();
        self.exons.max_end.clear();
        if let Some((list, max_ends)) = self.genes.chromosome(&region.chromosome) {
            if compat::java("BT2") {
                let mut running = i32::MIN;
                for i in self.genes.partition_exons(&region.chromosome, region.start, region.end) {
                    running = running.max(list[i].region.end);
                    self.exons.idx.push(i);
                    self.exons.max_end.push(running);
                }
            } else {
                self.exons.idx.extend(0..list.len());
                self.exons.max_end.extend_from_slice(max_ends);
            }
        }
        Ok(())
    }

    /// BamReader.processSamRecord
    #[inline]
    pub fn process(&mut self, rec: &Rec, out: &mut WorkerStats) -> Result<()> {
        let read_start = rec.start;
        let only_target = self.config.only_target_regions;
        if (read_start < self.start || read_start > self.end) && !only_target {
            if !self.java_bt1 && read_start < self.start {
                self.carry_over(rec)?;
            }
            return Ok(());
        }
        // BT18: with -only_target Java counts a read in every target region it overlaps. Corrected: the
        // read-level counts in the first one; its bases count in each region they fall in.
        if only_target && !self.java_bt18 {
            if let Some((ps, pe)) = self.previous {
                let read_end = rec.end.max(rec.start);
                if rec.start <= pe && read_end >= ps {
                    // this region's fragment counts (per region, as Java), then its bases
                    if self.overlaps_target(rec) || !self.java_bt5 {
                        let is_consensus = rec.has_cr;
                        let is_dual = is_consensus && rec.consensus_type().map_err(|e| anyhow!(e))? == ConsensusType::Dual;
                        self.check_target_regions(rec, is_consensus, is_dual);
                    }
                    return self.carry_over(rec);
                }
            }
        }
        if self.log_read_ids && self.config.log_read_ids.iter().any(|id| id.as_bytes() == rec.name) {
            crate::debug!("specific readId({}) unmapped({})", String::from_utf8_lossy(rec.name), rec.unmapped());
        }
        out.reads += 1;
        let is_consensus = rec.has_cr;
        let expect_duplicates = self.config.expect_duplicates();
        let partition = (self.tid, self.start, self.end);
        if self.record_off_target && !self.overlaps_target(rec) {
            self.read_counts.off_target += 1;
            self.partition_stats.process_read(rec, partition, is_consensus, expect_duplicates);
            // BT5: Java returns here, before checkTargetRegions, its only caller of OffTargetAnalyser.addRead,
            // which then sees only reads that overlap a target: the off-target analysis is never fed, and a
            // fragment whose first read lies outside the targets is not counted in them. Corrected: these
            // reads go through checkTargetRegions too.
            if !self.java_bt5 {
                let is_dual = is_consensus && rec.consensus_type().map_err(|e| anyhow!(e))? == ConsensusType::Dual;
                self.check_target_regions(rec, is_consensus, is_dual);
            }
            return Ok(());
        }
        let mut is_dual = false;
        if is_consensus {
            is_dual = rec.consensus_type().map_err(|e| anyhow!(e))? == ConsensusType::Dual;
            if expect_duplicates {
                self.read_counts.duplicates -= 1;
            } else {
                self.read_counts.total += 1;
            }
            if is_dual {
                self.read_counts.dual_strand += 1;
            }
        } else {
            self.read_counts.total += 1;
            self.check_log_partition_read_count();
            self.purge_read_groups(read_start);
            if rec.duplicate() {
                self.read_counts.duplicates += 1;
            }
        }
        self.flag_stats.process_read(rec, is_consensus, expect_duplicates);
        self.fragment_lengths.process_read(rec);
        self.partition_stats.process_read(rec, partition, is_consensus, expect_duplicates);
        self.check_target_regions(rec, is_consensus, is_dual);
        if rec.secondary() {
            return Ok(());
        }
        if let Some((es, ee)) = self.excluded {
            if !(es > rec.end || ee < read_start) {
                return Ok(());
            }
        }
        self.update_base_coverage(rec, is_consensus, true)?;
        self.check_gene_coverage(rec, out);
        Ok(())
    }

    /// a read already counted in an earlier partition (BT1) or target region (BT18): only its bases in
    /// this one, with its fragment's overlaps removed as for any read
    fn carry_over(&mut self, rec: &Rec) -> Result<()> {
        if rec.secondary() {
            return Ok(());
        }
        if let Some((es, ee)) = self.excluded {
            if !(es > rec.end || ee < rec.start) {
                return Ok(());
            }
        }
        self.update_base_coverage(rec, rec.has_cr, false)
    }

    fn overlaps_target(&self, rec: &Rec) -> bool {
        self.targets.iter().any(|t| !(rec.start > t.end || rec.end < t.start))
    }

    /// BamReader.checkTargetRegions: fragment counts of the target regions; off-target fragments to the analyser
    fn check_target_regions(&mut self, rec: &Rec, is_consensus: bool, is_dual: bool) {
        if self.targets.is_empty() && self.off_target.is_none() {
            return;
        }
        if rec.supplementary() {
            return;
        }
        let (s, e);
        if rec.paired() {
            let local = !rec.is_discordant();
            let mate_end = rec.mate_alignment_end();
            if local {
                // counted once per fragment, through its lower read
                if rec.start > rec.mate_start || (rec.start == rec.mate_start && rec.second_of_pair()) {
                    // BT5, with -only_target: a lower read lying wholly before the region is never read for
                    // it, so Java loses the fragment. Corrected: the upper read counts it then.
                    let lower_unseen = self.config.only_target_regions && !self.java_bt5 && mate_end > 0 && mate_end < self.start;
                    if !lower_unseen {
                        return;
                    }
                }
            }
            s = if local { rec.start.min(rec.mate_start) } else { rec.start };
            e = if local { rec.end.max(mate_end) } else { rec.end };
        } else {
            s = rec.start;
            e = rec.end;
        }
        let mut overlaps = false;
        let dup = rec.duplicate();
        for t in &mut self.targets {
            if s > t.end || e < t.start {
                continue;
            }
            overlaps = true;
            if is_consensus {
                t.fragment_counts.duplicates -= 1;
            } else {
                t.fragment_counts.total += 1;
            }
            if dup {
                t.fragment_counts.duplicates += 1;
            } else if is_dual {
                t.fragment_counts.dual_strand += 1;
            }
        }
        if !dup && !overlaps {
            if let Some(ot) = &mut self.off_target {
                ot.add_read(rec.name, rec.unmapped(), s, e);
            }
        }
    }

    /// the group-start test of findOverlappingFragmentStart: Java requires the partition to contain it
    #[inline]
    fn group_start_ok(&self, max_start: i32) -> bool {
        max_start > 0 && max_start <= self.end && (!self.java_bt1 || max_start >= self.start)
    }

    /// BamReader.findOverlappingFragmentStart: the start of the mate or of the first supplementary when
    /// it lies within the read (the later one), 0 for none
    fn overlapping_fragment_start(&self, rec: &Rec) -> Result<i32> {
        let mut max = 0;
        if rec.paired() && !rec.mate_unmapped() && rec.tid == rec.mtid {
            let mate = rec.mate_start;
            if mate >= rec.start && mate <= rec.end {
                max = mate;
            }
        }
        if let Some(sa) = rec.first_supplementary() {
            let (chr, pos) = sa.map_err(|e| anyhow!(e))?;
            if chr == self.chromosome.as_bytes() && pos >= rec.start && pos <= rec.end {
                max = max.max(pos);
            }
        }
        Ok(max)
    }

    /// BamReader.updateBaseCoverage with the read group look-up before it
    fn update_base_coverage(&mut self, rec: &Rec, is_consensus: bool, count_read: bool) -> Result<()> {
        let mut found = false;
        let mut removed: Option<ReadGroup> = None;
        if !self.groups.is_empty() {
            if let Some(g) = self.groups.get(rec.name) {
                found = true;
                if rec.start >= g.max_read_start {
                    removed = self.groups.remove(rec.name);
                }
            }
        }
        let eq_x = !compat::java("BT11");
        if let Some(g) = removed {
            self.coverage.process_read(rec, Some(&g.coords), is_consensus, count_read).map_err(|e| anyhow!(e))?;
            let mut coords = g.coords;
            coords.clear();
            self.spare_coords.push(coords);
            return Ok(());
        }
        if found {
            let g = self.groups.get_mut(rec.name).unwrap();
            self.coverage.process_read(rec, Some(&g.coords), is_consensus, count_read).map_err(|e| anyhow!(e))?;
            add_aligned_coords(rec, &mut g.coords, eq_x);
            return Ok(());
        }
        self.coverage.process_read(rec, None, is_consensus, count_read).map_err(|e| anyhow!(e))?;
        let max_start = self.overlapping_fragment_start(rec)?;
        if !self.group_start_ok(max_start) {
            return Ok(());
        }
        let mut coords = self.spare_coords.pop().unwrap_or_default();
        add_aligned_coords(rec, &mut coords, eq_x);
        self.groups.insert(rec.name.to_vec(), ReadGroup { max_read_start: max_start, is_consensus, coords });
        Ok(())
    }

    /// BamReader.checkGeneCoverage
    fn check_gene_coverage(&self, rec: &Rec, out: &mut WorkerStats) {
        if self.exons.idx.is_empty() {
            return;
        }
        if rec.supplementary() || rec.duplicate() || rec.unmapped() {
            return;
        }
        if rec.mapq < GENE_COVERAGE_MIN_MAP_QUALITY {
            return;
        }
        let (rs, re) = (rec.start, rec.end);
        let Some((list, _)) = self.genes.chromosome(&self.chromosome) else { return };
        // exons before the first whose running maximum end reaches the read cannot overlap it
        let first = self.exons.max_end.partition_point(|&m| m < rs);
        for &i in &self.exons.idx[first..] {
            let e = &list[i];
            if e.region.start > re {
                break;
            }
            if rs <= e.region.end && re >= e.region.start {
                let depth = &mut out.exon_depths[e.index];
                if depth.is_empty() {
                    depth.resize(e.region.base_length() as usize, 0);
                }
                genes::add_read(depth, e.region.start, e.region.end, rs, re);
            }
        }
    }

    fn check_log_partition_read_count(&self) {
        let check = self.config.partition_read_count_check as i64;
        if check == 0 {
            return;
        }
        let total = self.partition_stats.total_reads;
        if total > 0 && total % check == 0 {
            crate::debug!(
                "partition({}:{}-{}) reads(total={} chimeric={}, interpartition={}) cachedFragments({})",
                self.chromosome,
                self.start,
                self.end,
                total,
                self.partition_stats.chimeric_reads,
                self.partition_stats.inter_partition,
                self.groups.len()
            );
        }
    }

    /// BamReader.purgeReadGroups: a fall-back clearing of groups whose expected read was not seen
    fn purge_read_groups(&mut self, current_start: i32) {
        if self.groups.len() < FRAGMENT_MAP_SIZE_THRESHOLD {
            return;
        }
        if self.partition_stats.total_reads < self.last_purge_check + FRAGMENT_MAP_SIZE_READ_COUNT_CHECK {
            return;
        }
        self.last_purge_check = self.partition_stats.total_reads;
        let spare = &mut self.spare_coords;
        self.groups.retain(|_, g| {
            if g.max_read_start < current_start {
                spare.push(std::mem::take(&mut g.coords));
                false
            } else {
                true
            }
        });
    }

    /// BamReader.postSliceProcess, and what PartitionThread writes after each partition
    pub fn end(&mut self, out: &mut WorkerStats, task_id: usize, process_time: f64) {
        for (_, g) in self.groups.drain() {
            let mut c = g.coords;
            c.clear();
            self.spare_coords.push(c);
        }
        let config = self.config;
        let target_regions: Option<Vec<BaseRegion>> =
            if !config.target_regions.is_empty() { Some(self.targets.iter().map(|t| BaseRegion::new(t.start, t.end)).collect()) } else { None };
        self.coverage.add_metrics(target_regions.as_deref(), config.exclude_zero_coverage, &mut out.coverage);
        out.fragment_lengths.merge(&self.fragment_lengths);
        out.read_counts.merge(&self.read_counts);
        out.flag_stats.merge(&self.flag_stats);
        let mut enriched = Vec::new();
        if let Some(ot) = &mut self.off_target {
            ot.finish();
            for (k, v) in &ot.overlap_counts {
                *out.off_target_overlaps.entry(*k).or_insert(0) += *v;
            }
            enriched = ot.take_enriched();
        }
        let ps = &self.partition_stats;
        let partition_row = if !config.only_target_regions {
            Some(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.chromosome,
                self.start,
                self.end,
                ps.total_reads,
                ps.duplicate_reads,
                ps.chimeric_reads,
                ps.inter_partition,
                ps.unmapped_reads,
                hmf_common::fmt::java_fixed(process_time, 3)
            ))
        } else {
            None
        };
        let target_rows: Vec<String> = self.targets.iter().map(|t| t.row()).collect();
        out.rows.push((task_id, partition_row, target_rows, enriched));
    }

    #[cfg(test)]
    pub fn coverage(&self) -> &BaseCoverage {
        &self.coverage
    }
    #[cfg(test)]
    pub fn groups_empty(&self) -> bool {
        self.groups.is_empty()
    }
    #[cfg(test)]
    pub fn read_counts(&self) -> ReadCounts {
        self.read_counts
    }
}

/// BamReader.addAlignedCoords: the read's M blocks (and = and X, BT11), each once
fn add_aligned_coords(rec: &Rec, coords: &mut Vec<(i32, i32)>, eq_x_as_match: bool) {
    let mut position = rec.start;
    for &c in rec.cigar {
        let op = c & 0xf;
        let len = (c >> 4) as i32;
        if op == OP_M || (eq_x_as_match && (op == OP_EQ || op == OP_X)) {
            let block = (position, position + len - 1);
            if !coords.contains(&block) {
                coords.push(block);
            }
            position += len;
        } else if op == OP_D || op == OP_N {
            position += len;
        }
    }
}

#[cfg(test)]
mod tests {
    //! CoverageTest and MetricsTest
    use super::*;
    use crate::coverage::*;
    use crate::record::test_records::*;
    use hmf_common::metrics::FlagStatType;

    fn reader<'a>(config: &'a MetricsConfig, genes: &'a GeneCoverage, start: i32, end: i32) -> PartitionReader<'a> {
        let mut r = PartitionReader::new(config, genes);
        r.begin(&Task { id: 0, region: ChrBaseRegion::new("1", start, end), previous: None, unmapped: false }, 0).unwrap();
        r
    }

    fn feed(r: &mut PartitionReader, out: &mut WorkerStats, rec: &rust_htslib::bam::Record) {
        r.process(&Rec::of(rec), out).unwrap();
    }

    /// CoverageTest.testReadFilterTypes
    #[test]
    fn read_filter_types() {
        let config = MetricsConfig::for_tests(10);
        let params = CoverageParams { map_qual_threshold: 20, base_qual_threshold: 10, max_coverage: 10, expect_duplicates: true };
        let mut bc = BaseCoverage::new(params);
        bc.reset(1, 200, Vec::new());
        let base = Spec { name: "R1", start: 20, mate_start: 100, ..Default::default() };
        let read = make(&base);
        bc.process_read(&Rec::of(&read), None, false, true).unwrap();
        let mut low_mq = make(&Spec { name: "R2", ..base });
        low_mq.set_mapq((config.map_quality_threshold - 1) as u8);
        bc.process_read(&Rec::of(&low_mq), None, false, true).unwrap();
        let unmapped_mate = make(&Spec { name: "R3", mtid: -1, ..base });
        bc.process_read(&Rec::of(&unmapped_mate), None, false, true).unwrap();
        let mut dup = make(&Spec { name: "R4", ..base });
        dup.set_flags(dup.flags() | crate::record::FLAG_DUPLICATE);
        bc.process_read(&Rec::of(&dup), None, false, true).unwrap();
        bc.process_read(&Rec::of(&dup), None, false, true).unwrap();
        let mut consensus = make(&Spec { name: "R5", ..base });
        add_consensus(&mut consensus, 2, 1, "SINGLE");
        bc.process_read(&Rec::of(&consensus), None, true, true).unwrap();
        // low base qual
        let low_bq = {
            let mut r = rust_htslib::bam::Record::new();
            let cigar = rust_htslib::bam::record::CigarString::try_from("10M".as_bytes()).unwrap();
            r.set(b"R6", Some(&cigar), &[b'A'; 10], &[(config.base_quality_threshold - 1) as u8; 10]);
            r.set_tid(0);
            r.set_pos(19);
            r.set_mapq(60);
            r.set_mtid(0);
            r.set_mpos(99);
            r.set_flags(crate::record::FLAG_PAIRED | crate::record::FLAG_FIRST | crate::record::FLAG_PROPER_PAIR);
            r
        };
        bc.process_read(&Rec::of(&low_bq), None, false, true).unwrap();
        let mut m = CoverageMetrics::new(10);
        bc.add_metrics(None, false, &mut m);
        assert_eq!(m.filter_counts[UNFILTERED], 20);
        assert_eq!(m.filter_counts[LOW_MAP_QUAL], 10);
        assert_eq!(m.filter_counts[DUPLICATE], 10);
        assert_eq!(m.filter_counts[LOW_BASE_QUAL], 10);

        // exceeding max coverage: bases 30-39 reach 30x
        bc.reset(1, 200, Vec::new());
        let r1 = make(&Spec { name: "A", start: 10, bases_len: 30, cigar: "5S20M5S", mate_start: 100, ..Default::default() });
        let r2 = make(&Spec { name: "B", start: 15, bases_len: 30, cigar: "5S20M5S", mate_start: 100, ..Default::default() });
        let r3 = make(&Spec { name: "C", start: 20, bases_len: 30, cigar: "5S20M5S", mate_start: 100, ..Default::default() });
        for _ in 0..10 {
            for r in [&r1, &r2, &r3] {
                bc.process_read(&Rec::of(r), None, false, true).unwrap();
            }
        }
        let mut m = CoverageMetrics::new(10);
        bc.add_metrics(None, false, &mut m);
        assert!(bc.depth().iter().all(|d| *d <= 10));
        assert_eq!(m.filter_counts[UNFILTERED], 300);
        assert_eq!(m.filter_counts[MAX_COVERAGE], 300);
    }

    /// CoverageTest.testOverlappingReadCoverage
    #[test]
    fn overlapping_read_coverage() {
        for java in [true, false] {
            compat::with_java(java, || {
                let config = MetricsConfig::for_tests(10);
                let genes = GeneCoverage::empty();
                let mut out = WorkerStats::new(&config, 0);
                let mut r = reader(&config, &genes, 1, 500);
                let read1 = make(&Spec { name: "F1", start: 20, mate_start: 25, ..Default::default() });
                feed(&mut r, &mut out, &read1);
                let mate1 = make(&Spec { name: "F1", start: 25, mate_start: 20, is_reversed: true, ..Default::default() });
                feed(&mut r, &mut out, &mate1);
                let mut m = CoverageMetrics::new(10);
                r.coverage().add_metrics(None, false, &mut m);
                assert_eq!(m.filter_counts[UNFILTERED], 15);
                assert!(r.groups_empty());
                for (i, d) in r.coverage().depth().iter().enumerate() {
                    assert_eq!(*d, if (19..=33).contains(&i) { 1 } else { 0 }, "base {i}");
                }

                // overlapping alignments with a supplementary
                r.begin(&Task { id: 1, region: ChrBaseRegion::new("1", 1, 500), previous: None, unmapped: false }, 0).unwrap();
                let read_cigar = "5S10M10D10M5I10M5S";
                let mate_cigar = "5S10M5I10M10N10M5S";
                let supp_cigar = "5S10M10N10M5S";
                let read1 = make(&Spec { name: "F2", start: 10, bases_len: 45, cigar: read_cigar, mate_start: 30, supp: Some(("1", 40, supp_cigar)), ..Default::default() });
                feed(&mut r, &mut out, &read1);
                let mate1 = make(&Spec { name: "F2", start: 30, bases_len: 45, cigar: mate_cigar, mate_start: 10, is_reversed: true, ..Default::default() });
                feed(&mut r, &mut out, &mate1);
                let supp1 = make(&Spec {
                    name: "F2",
                    start: 40,
                    bases_len: 45,
                    cigar: supp_cigar,
                    mate_start: 30,
                    supplementary: true,
                    supp: Some(("1", 10, read_cigar)),
                    ..Default::default()
                });
                feed(&mut r, &mut out, &supp1);
                for (i, d) in r.coverage().depth().iter().enumerate() {
                    let expect = (9..=18).contains(&i) || (29..=48).contains(&i) || (59..=68).contains(&i);
                    assert_eq!(*d, expect as u32, "base {i}");
                }
                assert!(r.groups_empty());

                // the supplementary overlaps but not the mate
                r.begin(&Task { id: 2, region: ChrBaseRegion::new("1", 1, 500), previous: None, unmapped: false }, 0).unwrap();
                let read1 = make(&Spec { name: "F3", start: 10, bases_len: 20, cigar: "20M", mate_start: 30, supp: Some(("1", 20, "20M")), ..Default::default() });
                feed(&mut r, &mut out, &read1);
                let supp1 = make(&Spec { name: "F3", start: 20, bases_len: 20, cigar: "20M", mate_start: 10, supplementary: true, supp: Some(("1", 10, "20M")), ..Default::default() });
                feed(&mut r, &mut out, &supp1);
                assert!(r.groups_empty());
                let mate1 = make(&Spec { name: "F3", start: 50, bases_len: 20, cigar: "20M", mate_start: 10, is_reversed: true, ..Default::default() });
                feed(&mut r, &mut out, &mate1);
                for (i, d) in r.coverage().depth().iter().enumerate() {
                    let expect = (9..=38).contains(&i) || (49..=68).contains(&i);
                    assert_eq!(*d, expect as u32, "base {i}");
                }
                assert!(r.groups_empty());
            });
        }
    }

    /// MetricsTest.testReadCounts
    #[test]
    fn read_counts() {
        let config = MetricsConfig::for_tests(10);
        let genes = GeneCoverage::empty();
        let mut out = WorkerStats::new(&config, 0);
        let mut r = reader(&config, &genes, 1, 1000);
        let base = Spec { start: 20, mate_start: 100, ..Default::default() };
        feed(&mut r, &mut out, &make(&Spec { name: "A", ..base }));
        for name in ["B", "C"] {
            let mut d = make(&Spec { name, ..base });
            d.set_flags(d.flags() | crate::record::FLAG_DUPLICATE);
            feed(&mut r, &mut out, &d);
        }
        let mut c = make(&Spec { name: "D", ..base });
        add_consensus(&mut c, 2, 1, "SINGLE");
        feed(&mut r, &mut out, &c);
        r.end(&mut out, 0, 0.0);
        assert_eq!(out.read_counts.total, 3);
        assert_eq!(out.read_counts.duplicates, 1);
        assert_eq!(out.flag_stats.pass_count(FlagStatType::Primary), 3);
        assert_eq!(out.flag_stats.pass_count(FlagStatType::Duplicate), 1);
        assert_eq!(out.flag_stats.pass_count(FlagStatType::PrimaryDuplicate), 1);
    }

    /// MetricsTest.testDualStrandSupplementaries
    #[test]
    fn dual_strand_supplementaries() {
        for java in [true, false] {
            compat::with_java(java, || {
                let config = MetricsConfig::for_tests(10);
                let genes = GeneCoverage::empty();
                let mut out = WorkerStats::new(&config, 0);
                let mut r = reader(&config, &genes, 1, 1000);
                let base = Spec { start: 20, mate_start: 100, supplementary: true, ..Default::default() };
                let mut d = make(&Spec { name: "A", ..base });
                d.set_flags(d.flags() | crate::record::FLAG_DUPLICATE);
                feed(&mut r, &mut out, &d);
                let mut d = make(&Spec { name: "B", is_reversed: true, ..base });
                d.set_flags(d.flags() | crate::record::FLAG_DUPLICATE);
                feed(&mut r, &mut out, &d);
                let mut c = make(&Spec { name: "C", ..base });
                add_consensus(&mut c, 2, 1, "DUAL");
                feed(&mut r, &mut out, &c);
                r.end(&mut out, 0, 0.0);
                assert_eq!(out.read_counts.total, 2);
                assert_eq!(out.read_counts.duplicates, 1);
                assert_eq!(out.read_counts.dual_strand, 1);
                assert_eq!(out.flag_stats.pass_count(FlagStatType::Primary), 0);
                assert_eq!(out.flag_stats.pass_count(FlagStatType::Duplicate), 1);
                assert_eq!(out.flag_stats.pass_count(FlagStatType::PrimaryDuplicate), 0);
                assert_eq!(out.flag_stats.pass_count(FlagStatType::Supplementary), 2);
            });
        }
    }

    /// BT1: a read crossing the partition start adds its bases in the second partition in the corrected mode
    #[test]
    fn partition_boundary_bases() {
        for java in [true, false] {
            compat::with_java(java, || {
                let config = MetricsConfig::for_tests(10);
                let genes = GeneCoverage::empty();
                let mut out = WorkerStats::new(&config, 0);
                let mut r = reader(&config, &genes, 101, 200);
                // a pair: the first read 95-104 (starts in the earlier partition), its mate 100-109
                let read1 = make(&Spec { name: "X", start: 95, mate_start: 102, ..Default::default() });
                feed(&mut r, &mut out, &read1);
                let mate = make(&Spec { name: "X", start: 102, mate_start: 95, is_reversed: true, ..Default::default() });
                feed(&mut r, &mut out, &mate);
                let d = r.coverage().depth();
                if java {
                    // only the mate's bases
                    assert_eq!(&d[..12], &[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]);
                } else {
                    // 101-104 from the first read, 102-104 overlapped, 105-111 from the mate
                    assert_eq!(&d[..12], &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]);
                    assert_eq!(r.coverage().filter_counts[OVERLAPPED], 3);
                }
                assert_eq!(r.read_counts().total, 1);
            });
        }
    }
}
