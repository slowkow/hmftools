//! Per-partition processing (Java PartitionReader and BamWriter): filter, unmap, group duplicates, form consensus
//! reads and write in coordinate order.

use crate::bamin::BamReader;
use crate::capture::ReadCapture;
use crate::config::{Config, FilterReadsType};
use crate::consensus::{calc_base_qual_average, ConsensusReads};
use crate::dup::*;
use crate::out::RawBamOut;
use crate::rec::*;
use crate::regions::ChrBaseRegion;
use crate::sorted::SortedBamWriter;
use crate::unmap::{ReadUnmapper, UnmapRegionState, UnmapStats};
use anyhow::Result;
use rust_htslib::bam;

pub const SUPP_ALIGNMENT_SCORE_MIN: i64 = 30;

/// PartitionReader.shouldFilterRead: earlier consensus reads, weak supplementaries, and supplementaries
/// without supplementary data on contigs other than the chromosomes, MT and EBV
pub fn should_filter_read(rec: &bam::Record, expected_valid_tid: &[bool]) -> bool {
    if has_aux(rec, b"CR") {
        return true;
    }
    if has(rec, FLAG_SUPPLEMENTARY) {
        if let Some(score) = aux_int(rec, b"AS") {
            if score < SUPP_ALIGNMENT_SCORE_MIN {
                return true;
            }
        }
        if !has_aux(rec, b"SA") && !(rec.tid() >= 0 && expected_valid_tid[rec.tid() as usize]) {
            return true;
        }
    }
    false
}

/// PartitionReader.fullyUnmapped
pub fn fully_unmapped(rec: &bam::Record) -> bool {
    has(rec, FLAG_UNMAPPED) && (!has(rec, FLAG_PAIRED) || has(rec, FLAG_MATE_UNMAPPED))
}

pub struct PartitionOutput {
    pub stats: Statistics,
    pub consensus_stats: crate::consensus::ConsensusStatistics,
    pub unmap_stats: UnmapStats,
    pub unsorted_write_count: u64,
    pub processed_reads: u64,
    /// reads written to the partition's file of fully unmapped reads
    pub full_unmapped_written: u64,
}

pub struct PartitionReader<'a> {
    config: &'a Config,
    names: &'a [String],
    expected_valid_tid: &'a [bool],
    unmapper: &'a ReadUnmapper,
    unmap_paired_reads: bool,
    read_cache: ReadCache,
    /// SBX with a duplicate distance above zero (ReadCache.mSbxDuplicateCollapser)
    sbx_collapser: Option<crate::sbx_dup::SbxDuplicateCollapser>,
    /// the unmapper called a setter on the current read, so htsjdk writes it from its decoded fields
    read_decoded: bool,
    consensus: ConsensusReads,
    sorted: Option<SortedBamWriter>,
    out: Option<RawBamOut>,
    /// unpaired sequencing types: primary reads unmapped here have no mate to move to and go to the end of the BAM
    full_unmapped_out: Option<RawBamOut>,
    capture: ReadCapture<'a>,
    current_region: ChrBaseRegion,
    current_tid: i32,
    last_write_lower_position: i32,
    unmap_state: UnmapRegionState,
    pub stats: Statistics,
    pub unmap_stats: UnmapStats,
    unsorted_write_count: u64,
    processed_reads: u64,
    /// Test seam: without an output the records are kept here in the order written (Java TestBamWriter)
    #[cfg(test)]
    pub test_written: Vec<bam::Record>,
    ultima: crate::ultima::UltimaRoutines,
}

impl<'a> PartitionReader<'a> {
    /// See dup::CacheSync
    pub fn set_cache_sync(&mut self, sync: std::sync::Arc<crate::dup::CacheSync>, index: usize) {
        self.read_cache.set_sync(sync, index);
    }

    pub fn publish_cache_sync(&mut self) {
        self.read_cache.publish_sync();
    }

    pub fn new(
        config: &'a Config, names: &'a [String], expected_valid_tid: &'a [bool], unmapper: &'a ReadUnmapper, out: Option<RawBamOut>,
        full_unmapped_out: Option<RawBamOut>, capture: ReadCapture<'a>,
    ) -> Result<Self> {
        let mut consensus = ConsensusReads::new(&config.ref_genome_file, config.sequencing_type)?;
        consensus.run_checks = config.args.has_flag("run_checks");
        Ok(PartitionReader {
            config,
            names,
            expected_valid_tid,
            unmapper,
            unmap_paired_reads: unmapper.enabled() && config.is_illumina(),
            // non-Illumina reads are long and unpaired: wider position groups and soft-clip allowance
            // (SBX_READ_CACHE_GROUP_SIZE, SBX_READ_CACHE_MAX_SOFT_CLIP, SBX_READ_CACHE_LOG_READ_COUNT_THRESHOLD)
            read_cache: if config.is_illumina() {
                ReadCache::illumina(config.umi_enabled)
            } else {
                ReadCache::new(
                    crate::sbx::SBX_READ_CACHE_GROUP_SIZE, crate::sbx::SBX_READ_CACHE_MAX_SOFT_CLIP, false, DEFAULT_POP_DISTANCE_CHECK,
                    DEFAULT_DYNAMIC_READ_COUNT_THRESHOLD, crate::sbx::SBX_READ_CACHE_LOG_READ_COUNT_THRESHOLD,
                )
            },
            sbx_collapser: if config.is_sbx() && config.sbx_max_dup_dist > 0 {
                Some(crate::sbx_dup::SbxDuplicateCollapser::new(config.sbx_max_dup_dist))
            } else {
                None
            },
            read_decoded: false,
            consensus,
            sorted: if config.skip_duplicate_marking { None } else { Some(SortedBamWriter::new()) },
            out,
            full_unmapped_out,
            capture,
            current_region: ChrBaseRegion::new("", 0, 0),
            current_tid: -1,
            last_write_lower_position: 0,
            unmap_state: UnmapRegionState { chromosome: String::new(), slice_start: 0, slice_end: 0, regions: Vec::new(), last_matched_index: None },
            stats: Statistics::default(),
            unmap_stats: UnmapStats::default(),
            unsorted_write_count: 0,
            processed_reads: 0,
            #[cfg(test)]
            test_written: Vec::new(),
            ultima: crate::ultima::UltimaRoutines::new(config.ultima_drop_qual_tags),
        })
    }

    #[cfg(test)]
    pub fn consensus_mut(&mut self) -> &mut ConsensusReads {
        &mut self.consensus
    }
    /// PartitionReader.flushReadPositions
    #[cfg(test)]
    pub fn flush_read_positions(&mut self) -> Result<()> {
        let evicted = self.read_cache.evict_all();
        self.process_read_groups(evicted)
    }

    pub fn process_partition(&mut self, regions: &[ChrBaseRegion], reader: &mut BamReader) -> Result<()> {
        for region in regions {
            self.setup_region(region);
            let (rs, re) = (region.start, region.end);
            let mut result = Ok(());
            reader.slice_region(&region.chromosome, rs, re, |rec| {
                if result.is_ok() {
                    result = self.process_sam_record(rec);
                }
            });
            result?;
            self.post_process_region()?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<PartitionOutput> {
        if let Some(out) = self.out.take() {
            out.close()?;
        }
        let mut full_unmapped_written = 0;
        if let Some(out) = self.full_unmapped_out.take() {
            full_unmapped_written = out.written;
            out.close()?;
        }
        Ok(PartitionOutput {
            stats: self.stats,
            consensus_stats: self.consensus.stats.clone(),
            unmap_stats: self.unmap_stats,
            unsorted_write_count: self.unsorted_write_count,
            processed_reads: self.processed_reads,
            full_unmapped_written,
        })
    }

    /// PartitionReader.setupRegion
    pub(crate) fn setup_region(&mut self, region: &ChrBaseRegion) {
        self.current_region = region.clone();
        self.current_tid = self.names.iter().position(|n| *n == region.chromosome).map(|i| i as i32).unwrap_or(-1);
        let chromosome_length = self.consensus.refs().chromosome_length(&region.chromosome);
        self.consensus.set_chromosome_length(chromosome_length);
        self.last_write_lower_position = 0;
        let regions = if self.current_tid >= 0 {
            self.unmapper.regions(self.current_tid as usize).iter().filter(|r| r.start <= region.end && r.end >= region.start).copied().collect()
        } else {
            Vec::new()
        };
        self.unmap_state = UnmapRegionState {
            chromosome: region.chromosome.clone(),
            slice_start: region.start,
            slice_end: region.end,
            regions,
            last_matched_index: None,
        };
        self.capture.bqr.initialise(region);
        if let Some(s) = &mut self.sorted {
            s.initialise_start_position(self.current_tid, region.start);
        }
    }

    /// PartitionReader.postProcessRegion
    pub(crate) fn post_process_region(&mut self) -> Result<()> {
        let evicted = self.read_cache.evict_all();
        self.process_read_groups(evicted)?;
        self.capture.bqr.on_region_complete();
        if let (Some(s), Some(out)) = (&mut self.sorted, &mut self.out) {
            s.flush(out)?;
        }
        Ok(())
    }

    /// PartitionReader.processSamRecord
    pub(crate) fn process_sam_record(&mut self, rec: &mut bam::Record) -> Result<()> {
        let read_start = start(rec);
        if !self.current_region.contains_position(read_start) {
            return Ok(());
        }
        self.processed_reads += 1;
        if self.config.bqr_jitter_msi_only {
            self.capture.capture(rec, &self.current_region.chromosome);
            return Ok(());
        }
        if should_filter_read(rec, self.expected_valid_tid) {
            return Ok(());
        }
        self.check_missing_mate_cigar(rec)?;
        self.stats.total_reads += 1;
        if has(rec, FLAG_DUPLICATE) {
            if self.config.is_ultima() {
                return Ok(());
            }
            set_flag(rec, FLAG_DUPLICATE, false);
        }
        if self.config.specific_region_filter != FilterReadsType::None && self.read_outside_specified_regions(rec) {
            return Ok(());
        }
        if self.check_read_unmapping(rec)? {
            return Ok(());
        }
        self.pre_process_read(rec)?;
        let fresh = self.sorted.as_mut().and_then(|s| s.pool.pop()).unwrap_or_else(bam::Record::new);
        let read = std::mem::replace(rec, fresh);
        if has(&read, FLAG_SECONDARY) || self.config.skip_duplicate_marking {
            let mut read = read;
            if self.config.skip_duplicate_marking {
                self.post_process_read(&mut read);
            }
            self.set_boundary_position(read_start, false);
            self.write_read(read, false)?;
            return Ok(());
        }
        self.set_boundary_position(read_start, false);
        self.read_cache.process_read(read, self.names);
        let popped = self.read_cache.pop_reads();
        self.process_read_groups(popped)
    }

    /// PartitionReader.checkMissingMateCigar: with -fail_supp_no_mate_cigar a paired primary read with a mapped mate
    /// must have the mate cigar attribute; Java logs the read and exits
    fn check_missing_mate_cigar(&self, rec: &bam::Record) -> Result<()> {
        if !self.config.fail_on_missing_supp_mate_cigar {
            return Ok(());
        }
        if has(rec, FLAG_SUPPLEMENTARY | FLAG_SECONDARY) || !has(rec, FLAG_PAIRED) || has(rec, FLAG_MATE_UNMAPPED) {
            return Ok(());
        }
        if !has_aux(rec, b"MC") {
            anyhow::bail!("read({}) missing mate CIGAR", String::from_utf8_lossy(rec.qname()));
        }
        Ok(())
    }

    /// FilterReadsType.readOutsideSpecifiedRegions
    fn read_outside_specified_regions(&self, rec: &bam::Record) -> bool {
        let filter = self.config.specific_region_filter;
        let chromosomes = &self.config.specific.chromosomes;
        let regions = &self.config.specific.regions;
        let name = |tid: i32| -> &str { if tid >= 0 { &self.names[tid as usize] } else { "*" } };
        let filter_mates = matches!(filter, FilterReadsType::Mate | FilterReadsType::MateAndSupp);
        let paired = has(rec, FLAG_PAIRED);
        if !chromosomes.is_empty() && filter_mates && paired && !chromosomes.iter().any(|c| c == name(rec.mtid())) {
            return true;
        }
        if !regions.is_empty() {
            if filter == FilterReadsType::Read && !regions.iter().any(|r| r.chromosome == name(rec.tid()) && r.contains_position(start(rec))) {
                return true;
            }
            if filter_mates && paired && !regions.iter().any(|r| r.chromosome == name(rec.mtid()) && r.contains_position(mate_start(rec))) {
                return true;
            }
        }
        if filter == FilterReadsType::MateAndSupp {
            if let Some(sd) = crate::unmap::extract_supp_alignment(aux_str(rec, b"SA")) {
                let chr = String::from_utf8_lossy(sd.chromosome);
                if !regions.is_empty() && !regions.iter().any(|r| r.chromosome == chr && r.contains_position(sd.position)) {
                    return true;
                }
                if !chromosomes.is_empty() && !chromosomes.iter().any(|c| *c == chr) {
                    return true;
                }
            }
        }
        false
    }

    /// PartitionReader.checkReadUnmapping: true when the read is dropped from further processing
    fn check_read_unmapping(&mut self, rec: &mut bam::Record) -> Result<bool> {
        if !self.unmapper.enabled() {
            return Ok(false);
        }
        if self.unmap_paired_reads {
            // a read unmapped by the first pass arrives again at its mate's position, marked with UM
            let is_unmapped = has(rec, FLAG_UNMAPPED);
            let internally_unmapped = is_unmapped && has_aux(rec, b"UM");
            if !internally_unmapped {
                self.unmapper.check_transform_read(rec, &mut self.unmap_state, &mut self.unmap_stats);
                let fully = fully_unmapped(rec);
                let unmapped = !is_unmapped && has(rec, FLAG_UNMAPPED);
                if unmapped || fully {
                    return Ok(true);
                }
            }
            return Ok(false);
        }
        // unpaired sequencing types have no first pass: reads are unmapped here
        self.read_decoded = false;
        if self.unmapper.check_transform_read(rec, &mut self.unmap_state, &mut self.unmap_stats) {
            if !has(rec, FLAG_UNMAPPED) {
                self.read_decoded = true;
                return Ok(false);
            }
            if has(rec, FLAG_SUPPLEMENTARY | FLAG_SECONDARY) {
                return Ok(true);
            }
            if !self.config.skip_fully_unmapped_reads {
                self.strip_unmapped_read_attributes(rec)?;
                if let Some(out) = &mut self.full_unmapped_out {
                    out.write(rec)?;
                }
            }
            return Ok(true);
        }
        Ok(false)
    }

    /// PartitionReader.preProcessRead: sequencing-type specific changes to a read before duplicate marking.
    /// Illumina: none. Ultima: UltimaRoutines.preProcessRead. SBX: SbxRoutines.prepProcessRead.
    fn pre_process_read(&mut self, rec: &mut bam::Record) -> Result<()> {
        if self.config.is_ultima() {
            self.ultima.pre_process_read(rec, &mut self.stats.ultima)?;
        } else if self.config.is_sbx() {
            crate::sbx::prep_process_read(rec, self.read_decoded);
        }
        Ok(())
    }

    /// PartitionReader.postProcessRead: sequencing-type specific changes to a read that is written as a
    /// non-duplicate or as the primary or consensus read of a group. Illumina: none. SBX: SbxRoutines.finaliseRead.
    fn post_process_read(&mut self, rec: &mut bam::Record) {
        if self.config.is_sbx() {
            crate::sbx::finalise_read(self.consensus.refs(), self.names, rec);
        }
    }

    /// Ultima: UltimaRoutines.stripAttributes on a read that goes to the unmapped reads.
    fn strip_unmapped_read_attributes(&mut self, rec: &mut bam::Record) -> Result<()> {
        if self.config.is_ultima() {
            self.ultima.strip_attributes(rec, Some(&mut self.stats.ultima))?;
        }
        Ok(())
    }

    fn set_boundary_position(&mut self, position: i32, is_lower: bool) {
        if let Some(s) = &mut self.sorted {
            if is_lower {
                s.set_upper_writable_position(position);
            } else {
                s.set_upper_bound_position(position);
            }
        }
    }

    /// BamWriter.processRecord: capture for jitter and BQR, then write
    fn process_record(&mut self, mut rec: bam::Record) -> Result<()> {
        if self.config.is_sbx() {
            crate::sbx::restore_stored_quals(&mut rec);
        }
        if rec.tid() >= 0 {
            self.capture.capture(&rec, &self.names[rec.tid() as usize]);
        }
        #[cfg(test)]
        if self.out.is_none() {
            self.test_written.push(rec);
            return Ok(());
        }
        let Some(out) = &mut self.out else { return Ok(()) };
        match &mut self.sorted {
            Some(s) => {
                if s.can_write_record(&rec) {
                    s.add_record(rec, out)?;
                } else {
                    self.unsorted_write_count += 1;
                }
            }
            None => out.write(&rec)?,
        }
        Ok(())
    }

    /// BamWriter.writeRead
    fn write_read(&mut self, mut rec: bam::Record, duplicate: bool) -> Result<()> {
        if duplicate {
            if self.config.drop_duplicates {
                return Ok(());
            }
            set_flag(&mut rec, FLAG_DUPLICATE, true);
        }
        self.process_record(rec)
    }

    /// PartitionReader.processReadGroups
    fn process_read_groups(&mut self, popped: Option<FragmentCoordReads>) -> Result<()> {
        let Some(mut fcr) = popped else { return Ok(()) };
        // ReadCache.collectFragments
        if let Some(collapser) = &self.sbx_collapser {
            fcr = collapser.collapse_groups(fcr, self.names);
        }
        let min_cached_read_position = self.read_cache.min_cached_read_start();
        let min_popped_reads_position = fcr.min_read_position_start();
        let read_flush_position = if min_cached_read_position > 0 { min_cached_read_position - 1 } else { min_popped_reads_position - 1 };
        if self.last_write_lower_position > 0 && min_popped_reads_position < self.last_write_lower_position {
            crate::warn!(
                "position({}:{})) processing earlier popped read min start({}) vs last({}) minCachedPos({})",
                self.current_region.chromosome, self.read_cache.current_read_min_position(), min_popped_reads_position,
                self.last_write_lower_position, min_cached_read_position
            );
        }
        self.last_write_lower_position = read_flush_position;

        // DuplicateGroupBuilder.processDuplicateGroups
        if self.config.umi.enabled {
            crate::umi::process_duplicate_groups(&self.config.umi, self.config.is_illumina(), &mut fcr, &mut self.stats)?;
        } else {
            for group in fcr.duplicate_groups.iter_mut() {
                self.stats.add_duplicate_group(group.total_read_count());
                if !self.config.form_consensus {
                    set_primary_read(group);
                }
            }
        }

        for mut group in fcr.duplicate_groups.drain(..) {
            if self.config.form_consensus && !group.poly_g_unmapped {
                // the header is only needed for the read group of new records, which is copied from the template
                self.consensus.form_consensus_read_plain(&mut group, self.names);
                if let Some(read_id) = self.consensus.undetermined_base_read.take() {
                    return Err(crate::ultima::undetermined_consensus_base_error(&read_id));
                }
                let pos = start(group.consensus_read.as_ref().unwrap());
                self.set_boundary_position(pos, false);
            }
            // PartitionReader.postProcessPrimaryRead
            if let Some(consensus) = group.consensus_read.as_mut() {
                self.post_process_read(consensus);
                if self.config.is_sbx() && self.config.form_consensus && !group.poly_g_unmapped {
                    if let Some(i) = self.consensus.built_from_read {
                        crate::sbx::share_finalised_bases(consensus, &mut group.reads[i]);
                    }
                }
            } else if let Some(i) = group.primary_read {
                self.post_process_read(&mut group.reads[i]);
            }
            // BamWriter.writeDuplicateGroup
            if let Some(consensus) = group.consensus_read.take() {
                self.process_record(consensus)?;
            }
            let primary = group.primary_read;
            let poly_g = group.poly_g_unmapped;
            // with UMIs every read of the group gets the group's UMI as its UI attribute
            let umi_attribute = if self.config.umi.enabled { group.umi.as_ref().map(|u| encode_z(u.as_bytes())) } else { None };
            for (i, mut read) in group.reads.drain(..).enumerate() {
                if let Some(value) = &umi_attribute {
                    set_attribute(&mut read, b"UI", Some(value));
                }
                let is_primary = !poly_g && primary == Some(i);
                self.write_read(read, !is_primary)?;
            }
            for mut read in group.non_consensus_reads.drain(..) {
                if let Some(value) = &umi_attribute {
                    set_attribute(&mut read, b"UI", Some(value));
                }
                self.write_read(read, true)?;
            }
        }
        let single_count = fcr.single_reads.len();
        for info in fcr.single_reads.iter_mut() {
            self.post_process_read(&mut info.read);
        }
        for info in fcr.single_reads.drain(..) {
            self.write_read(info.read, false)?;
        }
        self.stats.add_non_duplicate_counts(single_count);
        self.set_boundary_position(read_flush_position, true);
        Ok(())
    }
}

#[cfg(test)]
#[path = "umi_partition_test.rs"]
mod umi_partition_test;

/// DuplicateGroupBuilder.setPrimaryRead: the first read with the highest average base quality
pub(crate) fn set_primary_read(group: &mut DuplicateGroup) {
    let mut max_read = None;
    let mut max_qual = 0.0;
    for (i, read) in group.reads.iter().enumerate() {
        let avg = calc_base_qual_average(read.qual());
        if avg > max_qual {
            max_qual = avg;
            max_read = Some(i);
        }
    }
    group.primary_read = max_read;
}
