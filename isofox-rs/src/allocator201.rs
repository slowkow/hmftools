//! FragmentAllocator as ISOFOX 2.0.1 has it, for STAR-aligned BAMs duplicate-marked by GATK MarkDuplicates.
//!
//! What differs from 2.1 (allocator.rs):
//! - duplicates are read (unless `-drop_dups`, or only fusions or novel junctions are run) and counted in
//!   expression; they add no depth, no chimeric junction data and no novel junctions or retained introns;
//! - secondary alignments are read when transcript counts run, each named `<read name>_<HI>` so that a
//!   secondary pair forms its own fragment; they are counted in expression (weighted by map quality) but
//!   not in the fragment type counts or strandedness, nor as chimeric fragments;
//! - a fragment whose lower map quality is 3 or less is multi-mapped (STAR's convention), and its weight in
//!   expression is 0.5 (map quality 3), 0.33 (2), 0.2 (1) or 0.1 (0); such fragments are LowMapQual;
//! - supplementary alignments pair with whichever read of the same name comes next;
//! - a read's exon regions and transcripts are worked out when it is read, before adapter trimming;
//! - reads in an excluded region (`-excluded_regions`, and each enriched gene with 1 kb) are skipped, and
//!   the BAM is sliced around those regions;
//! - fragments in an enriched gene's exons (RN7SL1-5, RN7SK) are counted without further classification.

use crate::allocator::FragmentAllocator;
use crate::bam::{aux_string, contig_name, rec_end, rec_mate_start, rec_start, to_read, SliceFilter};
use crate::chimeric::is_realigned_fragment_candidate_201;
use crate::common::*;
use crate::config::{IsofoxConfig, IsofoxFunction, MULTI_MAP_QUALITY_THRESHOLD};
use crate::fragment;
use crate::gene::GeneCollection;
use crate::read::{Read, SuppData};
use crate::read_trans::{calc_fragment_length, find_overlapping_regions, mark_region_bases, process_overlapping_regions};
use crate::refgen::RefSource;
use anyhow::{bail, Result};
use hmf_common::regions::BaseRegion;
use rust_htslib::bam;

const NON_GENIC_BASE_DEPTH_WIDTH: i32 = 250_000;

/// the BamSlicer of 2.0.1's FragmentAllocator
pub fn slice_filter_201(config: &IsofoxConfig) -> SliceFilter {
    let keep_duplicates =
        (config.run_function(IsofoxFunction::TranscriptCounts) || config.run_function(IsofoxFunction::Statistics)) && !config.drop_duplicates;
    let keep_supplementaries = config.run_function(IsofoxFunction::Fusions) || config.run_function(IsofoxFunction::AltSpliceJunctions);
    let keep_secondaries = config.run_function(IsofoxFunction::TranscriptCounts);
    let min_map_quality = if keep_secondaries { 0 } else { config.single_map_quality };
    SliceFilter { min_map_quality, keep_duplicates, keep_supplementaries, keep_secondaries }
}

/// the expression weight of a fragment by its lower map quality (ExpressionReadTracker.addGcCounts)
pub fn map_quality_weight(min_map_quality: u8) -> f64 {
    match min_map_quality {
        3 => 0.5,
        2 => 0.33,
        1 => 0.2,
        q if q <= MULTI_MAP_QUALITY_THRESHOLD => 0.1,
        _ => 1.0,
    }
}

impl<'a> FragmentAllocator<'a> {
    /// processBam
    pub(crate) fn process_bam_201(&mut self, gc: &mut GeneCollection, region: [i32; 2], refg: &mut RefSource) -> Result<()> {
        self.clear_cache();
        self.gene_read_count = 0;
        self.enriched_gene_fragments = 0;
        let depth_range = [region[0].max(gc.region_bounds[0] - NON_GENIC_BASE_DEPTH_WIDTH), region[1].min(gc.region_bounds[1] + NON_GENIC_BASE_DEPTH_WIDTH)];
        self.base_depth.initialise(depth_range);
        self.chimeric.initialise(gc);
        if self.expression_enabled {
            self.set_expression_gene_data();
        }
        self.valid_read_start = region;

        let chromosome = gc.chromosome;
        let excluded = self.config.filters.find_excluded_regions(chromosome, region[0], region[1]);
        let mut slices: Vec<[i32; 2]> = Vec::new();
        if !excluded.is_empty() {
            // a gene collection wholly within an excluded region is skipped, without chimeric post-processing
            if excluded.iter().any(|x| x.start <= region[0] && x.end >= region[1]) {
                return Ok(());
            }
            let mut start = region[0];
            for x in &excluded {
                if start < x.start {
                    slices.push([start, x.start - 1]);
                }
                start = x.end + 1;
            }
            if start < region[1] {
                slices.push([start, region[1]]);
            }
            // ChrBaseRegion.isValid
            slices.retain(|s| s[0] > 0 && s[1] >= s[0]);
        } else {
            slices.push(region);
        }

        let Some(mut reader) = self.reader.take() else { return Ok(()) };
        let filter = self.filter;
        let mut res = Ok(());
        for s in slices {
            // BamSlicer.slice resets the halt for each region
            self.halted = false;
            res = reader.slice(chromosome, s[0], s[1], &filter, |rec: &bam::Record| {
                self.process_sam_record_201(rec, gc, refg);
                !self.halted
            });
            if res.is_err() || self.failure.is_some() {
                break;
            }
        }
        self.reader = Some(reader);
        res?;
        if let Some(f) = self.failure.take() {
            bail!("{f}");
        }

        if self.enriched_gene_fragments > 0 {
            self.process_enriched_gene_fragments(gc, refg);
        }
        if self.chimeric.enabled {
            let lone: Vec<usize> = self.tracker.values().copied().collect();
            let immune = std::mem::take(&mut self.immune);
            self.chimeric.post_process_chimeric_reads(&mut self.base_depth, &lone, &mut self.reads, gc, &immune, self.known_splice_sites);
            self.immune = immune;
            if let Some(f) = self.chimeric.failure.take() {
                bail!("{f}");
            }
            self.process_chimeric_novel_junctions(gc, refg);
        }
        Ok(())
    }

    /// GeneRegionFilters.skipRead(record, true) in 2.0.1
    fn skip_record_201(&mut self, rec: &bam::Record) -> bool {
        let f = &self.config.filters;
        if f.skip_read(contig_name(rec.tid()), rec_start(rec), rec_end(rec), false) {
            return true;
        }
        if !rec.is_paired() {
            // htsjdk's getMateUnmappedFlag throws for an unpaired read, which ends Java's run (IF15)
            self.failure = Some(format!(
                "read({}) is not paired: ISOFOX 2.0.1 fails on unpaired reads (htsjdk: Inappropriate call if not paired read)",
                crate::bam::record_name(rec)
            ));
            self.halted = true;
            return true;
        }
        if !rec.is_mate_unmapped() {
            let ms = rec_mate_start(rec);
            if f.skip_read(contig_name(rec.mtid()), ms, ms, true) {
                return true;
            }
        }
        if let Some(sa) = aux_string(rec, b"SA") {
            if let Some(s) = SuppData::extract(&sa) {
                if f.skip_read(&s.chromosome, s.position, s.position, true) {
                    return true;
                }
            }
        }
        false
    }

    fn process_sam_record_201(&mut self, rec: &bam::Record, gc: &mut GeneCollection, refg: &mut RefSource) {
        let start = rec_start(rec);
        if !position_within(start, self.valid_read_start[0], self.valid_read_start[1]) {
            return;
        }
        if self.skip_record_201(rec) {
            if self.failure.is_none() {
                self.chimeric.stats.excluded += 1;
            }
            return;
        }
        // trackFragmentCounts: each fragment once, by its first read, without supplementaries and secondaries
        self.total_bam_read_count += 1;
        self.gene_read_count += 1;
        if !rec.is_supplementary() && !rec.is_secondary() && (!rec.is_paired() || rec.is_first_in_template()) {
            gc.add_count(FragmentType::Total, 1);
            if rec.is_duplicate() {
                gc.add_count(FragmentType::Duplicate, 1);
            }
        }
        if gc.in_enriched_region(start, rec_end(rec)) {
            self.process_enriched_region_read(rec, gc);
            return;
        }
        let read = to_read(rec);
        if self.config.log_read_ids.iter().any(|x| **x == *read.id) {
            hmf_common::debug!("specific read: {} range({}: {} -> {}) flags({}) cigar({})", read.id, read.chromosome, read.pos_start, read.pos_end, read.flags, read.cigar_str());
        }
        self.reads.push(read);
        let ri = self.reads.len() - 1;
        self.process_read_201(ri, gc, refg);
    }

    fn reached_gene_read_limit(&mut self, gc: &GeneCollection) -> bool {
        if self.config.gene_read_limit == 0 || self.gene_read_count < self.config.gene_read_limit {
            return false;
        }
        self.halted = true;
        hmf_common::warn!("chr({}) genes({}) readCount({}) exceeds max read count", gc.chromosome, gc.gene_names(10), self.gene_read_count);
        true
    }

    /// processEnrichedRegionRead
    fn process_enriched_region_read(&mut self, rec: &bam::Record, gc: &mut GeneCollection) {
        if self.reached_gene_read_limit(gc) {
            return;
        }
        if !rec.is_paired() || rec.is_secondary() || rec.is_supplementary() {
            return;
        }
        if rec.is_reverse() == rec.is_mate_reverse() {
            return;
        }
        if rec.tid() != rec.mtid() {
            return;
        }
        if !rec.is_first_in_template() {
            return;
        }
        self.enriched_gene_fragments += 1;
        gc.add_count(FragmentType::TransSupporting, 1);
    }

    /// ExpressionReadTracker.processEnrichedGeneFragments
    fn process_enriched_gene_fragments(&mut self, gc: &mut GeneCollection, refg: &mut RefSource) {
        if !self.expression_enabled {
            return;
        }
        let region = gc.enriched_region.unwrap_or([0, 0]);
        let genes: Vec<String> = gc.find_genes_covering_range(region[0], region[1], true).iter().map(|g| gc.genes[*g].gene.gene_id.clone()).collect();
        let trans: Vec<i32> = gc.enriched_transcripts.as_ref().map(|l| l.iter().map(|t| t.trans_id).collect()).unwrap_or_default();
        let ci = self.category_counts(&trans, &genes);
        let n = self.enriched_gene_fragments as f64;
        let mut idx = [-1, -1];
        let mut counts = [0.0, 0.0];
        if let Some(g) = &self.gc_ratio_counts {
            let ratio = refg.gc_ratio_of_regions(gc.chromosome, &[BaseRegion::new(region[0], region[1])]);
            g.determine_ratio_data(ratio, &mut idx, &mut counts);
            counts[0] *= n;
            counts[1] *= n;
        }
        self.add_gc_counts_raw_201(ci, &idx, &counts, n);
    }

    fn process_read_201(&mut self, ri: usize, gc: &mut GeneCollection, refg: &mut RefSource) {
        if self.reached_gene_read_limit(gc) {
            return;
        }
        // the read's regions and transcripts, before the pair's adapter trimming
        let overlapping = find_overlapping_regions(&gc.exon_regions, None, &self.reads[ri]);
        if !overlapping.is_empty() {
            process_overlapping_regions(&mut self.reads[ri], &gc.exon_regions, &overlapping);
        }
        let vr = self.valid_read_start;
        gc.set_read_gene_collections(&mut self.reads[ri], vr);
        // checkFragmentRead: paired with the read of the same name, supplementaries included
        let id = self.reads[ri].id.to_string();
        match self.tracker.remove(&id) {
            Some(other) => self.process_fragment_reads_201(ri, other, gc, refg),
            None => {
                self.tracker.insert(id, ri);
            }
        }
    }

    fn process_fragment_reads_201(&mut self, r1: usize, r2: usize, gc: &mut GeneCollection, refg: &mut RefSource) {
        {
            let (a, b) = self.two_reads(r1, r2);
            let (o, s, e) = (b.orient(), b.unclipped_start, b.unclipped_end);
            a.trim_adapter_soft_clip_bases_201(o, s, e);
            let (o, s, e) = (a.orient(), a.unclipped_start, a.unclipped_end);
            b.trim_adapter_soft_clip_bases_201(o, s, e);
        }
        let frag = [r1, r2];
        let is_duplicate = self.reads[r1].is_duplicate() || self.reads[r2].is_duplicate();
        let min_mq = self.reads[r1].mapq.min(self.reads[r2].mapq);
        let multi_mapped = min_mq <= MULTI_MAP_QUALITY_THRESHOLD;
        let is_chimeric = self.chimeric.is_chimeric_201(&frag, &mut self.reads, gc, is_duplicate, multi_mapped);
        let secondary = self.reads[r1].is_secondary() || self.reads[r2].is_secondary();
        let fully_intronic = fragment::is_fully_intronic(&self.frag_reads(&frag));

        if self.stats_only {
            if !is_duplicate {
                if is_chimeric {
                    gc.add_count(FragmentType::Chimeric, 1);
                } else if fully_intronic {
                    gc.add_count(FragmentType::Unspliced, 1);
                } else {
                    gc.add_count(FragmentType::TransSupporting, 1);
                }
            }
            return;
        }

        let common = derive_common_regions(&self.reads[r1].mapped_coords.alignments, &self.reads[r2].mapped_coords.alignments);
        if !is_duplicate {
            self.base_depth.process_read(&common);
        }

        if is_chimeric {
            if !multi_mapped && !secondary {
                if self.chimeric.enabled {
                    let immune = std::mem::take(&mut self.immune);
                    self.chimeric.add_chimeric_read_pair_201(&frag, &mut self.reads, gc, &immune);
                    self.immune = immune;
                } else {
                    gc.add_count(FragmentType::Chimeric, 1);
                }
            }
            return;
        }

        if self.run_fusions {
            if frag.iter().any(|r| is_realigned_fragment_candidate_201(&self.reads[*r])) {
                self.chimeric.add_realignment_candidates(&frag, &self.reads);
            }
            if self.fusions_only {
                return;
            }
        }

        let min_pos = self.reads[r1].pos_start.min(self.reads[r2].pos_start);
        let max_pos = self.reads[r1].pos_end.max(self.reads[r2].pos_end);
        let overlap_genes = gc.find_genes_covering_range(min_pos, max_pos, true);
        let contains_split = self.reads[r1].has_split || self.reads[r2].has_split;
        if fully_intronic {
            // processIntronicReads. IF17: Java counts a secondary pair here (only here) as a fragment too
            let count_it = !secondary || crate::compat::java("IF17");
            if contains_split {
                if count_it {
                    gc.add_count(FragmentType::Alt, 1);
                }
                if self.alt_sj.enabled {
                    let reads: Vec<&Read> = frag.iter().map(|r| &self.reads[*r]).collect();
                    self.alt_sj.evaluate_fragment_reads(gc, refg, &overlap_genes, &reads, &[]);
                }
                return;
            }
            if self.expression_enabled && !overlap_genes.is_empty() {
                let genes: Vec<String> = overlap_genes.iter().map(|g| gc.genes[*g].gene.gene_id.clone()).collect();
                let ci = self.category_counts(&[], &genes);
                self.add_gc_counts_201(ci, &common, min_mq, gc, refg);
            }
            if count_it {
                gc.add_count(FragmentType::Unspliced, 1);
            }
            return;
        }

        let valid_regions = fragment::unique_valid_regions(&self.frag_reads(&frag));
        let mut valid_trans: Vec<i32> = Vec::new();
        for t in fragment::valid_type_transcripts(&self.frag_reads(&frag)) {
            if let Some(td) = gc.transcript(t) {
                let len = calc_fragment_length(td, min_pos, max_pos);
                if len > 0 && len <= self.config.max_fragment_length {
                    valid_trans.push(t);
                }
            }
        }
        let invalid = fragment::invalid_transcripts(&self.frag_reads(&frag), &valid_trans);

        let mut fragment_type = FragmentType::Unspliced;
        if valid_trans.is_empty() {
            if contains_split {
                fragment_type = FragmentType::Alt;
                if self.alt_sj.enabled {
                    let reads: Vec<&Read> = frag.iter().map(|r| &self.reads[*r]).collect();
                    let related = invalid.keys();
                    self.alt_sj.evaluate_fragment_reads(gc, refg, &overlap_genes, &reads, &related);
                }
            } else {
                for t in invalid.keys() {
                    if fragment::spans_multiple_regions(&self.frag_reads(&frag), &gc.exon_regions, t) {
                        fragment_type = FragmentType::Alt;
                        break;
                    }
                }
            }
            if self.retained_introns.enabled {
                let reads: Vec<&Read> = frag.iter().map(|r| &self.reads[*r]).collect();
                self.retained_introns.evaluate_fragment_reads(gc, &reads);
            }
            if fragment_type == FragmentType::Unspliced && self.expression_enabled && !overlap_genes.is_empty() {
                let genes: Vec<String> = overlap_genes.iter().map(|g| gc.genes[*g].gene.gene_id.clone()).collect();
                let ci = self.category_counts(&valid_trans, &genes);
                self.add_gc_counts_201(ci, &common, min_mq, gc, refg);
            }
        } else {
            fragment_type = FragmentType::TransSupporting;
            {
                let mut rs = self.frag_reads_mut(&frag);
                fragment::set_other_transcripts(&mut rs, &valid_trans);
            }
            for &reg in &valid_regions {
                mark_region_bases(&common, &mut gc.exon_regions[reg]);
            }
            let is_unique = valid_trans.len() == 1;
            let mut gene_is_forward: Option<bool> = None;
            let mut combo = FragmentMatchType::Short;
            for &t in &valid_trans {
                let region_count = valid_regions.iter().filter(|r| gc.exon_regions[**r].has_trans_id(t)).count();
                if gene_is_forward.is_none() {
                    for &ri in &frag {
                        gene_is_forward = self.find_gene_strand(ri, gc, &valid_trans);
                        if gene_is_forward.is_some() {
                            break;
                        }
                    }
                }
                let mt = if fragment::has_transcript_classification(&self.frag_reads(&frag), t, TransMatchType::SpliceJunction) {
                    combo = FragmentMatchType::Spliced;
                    FragmentMatchType::Spliced
                } else if region_count > 1 {
                    if combo != FragmentMatchType::Spliced {
                        combo = FragmentMatchType::Long;
                    }
                    FragmentMatchType::Long
                } else {
                    FragmentMatchType::Short
                };
                gc.add_transcript_read_match(t, is_unique, mt);
                if !contains_split && fragment::reads_in_different_exons(&self.frag_reads(&frag), &gc.exon_regions, t) {
                    gc.add_transcript_read_match_type(t, FragmentMatchType::Discordant);
                }
                self.process_valid_transcript(t, &frag, gc, is_unique);
            }
            if self.expression_enabled {
                let genes: Vec<String> = if combo == FragmentMatchType::Short { overlap_genes.iter().map(|g| gc.genes[*g].gene.gene_id.clone()).collect() } else { Vec::new() };
                let ci = self.category_counts(&valid_trans, &genes);
                self.add_gc_counts_201(ci, &common, min_mq, gc, refg);
            }
            if !secondary {
                if let Some(fwd) = gene_is_forward {
                    // the first-of-pair read's strand against the gene's
                    let (a, b) = (&self.reads[r1], &self.reads[r2]);
                    let first_fwd = if a.is_first_of_pair() { !a.is_read_reversed() } else { !b.is_read_reversed() };
                    let second_fwd = if !a.is_first_of_pair() { !a.is_read_reversed() } else { !b.is_read_reversed() };
                    if first_fwd != second_fwd {
                        if first_fwd == fwd {
                            gc.add_count(FragmentType::ForwardStrand, 1);
                        } else {
                            gc.add_count(FragmentType::ReverseStrand, 1);
                        }
                    }
                }
            }
        }
        if !secondary {
            gc.add_count(fragment_type, 1);
        }
    }

    /// ExpressionReadTracker.addGcCounts(catCounts, readRegions, minMapQuality): the fragment's weight by its
    /// lower map quality; a fragment of map quality 3 or less is LowMapQual
    fn add_gc_counts_201(&mut self, ci: usize, regions: &[BaseRegion], min_map_quality: u8, gc: &mut GeneCollection, refg: &mut RefSource) {
        let mut idx = [-1, -1];
        let mut counts = [0.0, 0.0];
        if let Some(g) = &self.gc_ratio_counts {
            let ratio = refg.gc_ratio_of_regions(gc.chromosome, regions);
            g.determine_ratio_data(ratio, &mut idx, &mut counts);
        }
        let weight = map_quality_weight(min_map_quality);
        if min_map_quality <= MULTI_MAP_QUALITY_THRESHOLD {
            gc.add_count(FragmentType::MultiMapped, 1);
        }
        self.add_gc_counts_raw_201(ci, &idx, &counts, weight);
    }

    /// addGcCounts(catCounts, gcRatioIndices, gcRatioCounts, count): CategoryCountsData counts a fragment below
    /// a weight of 1 as LowMapQual
    fn add_gc_counts_raw_201(&mut self, ci: usize, idx: &[i32; 2], counts: &[f64; 2], count: f64) {
        let low = count < 1.0;
        if self.gc_ratio_counts.is_some() {
            for i in 0..2 {
                if idx[i] >= 0 {
                    self.gc_ratio_counts.as_mut().unwrap().add_gc_ratio_count(idx[i], counts[i]);
                    self.gene_gc_ratio_counts.as_mut().unwrap().add_gc_ratio_count(idx[i], counts[i]);
                }
            }
            if self.config.apply_gc_bias_adjust() {
                self.trans_combo_data[ci].add_gc_ratio_counts(count, low, idx, counts);
            } else {
                self.trans_combo_data[ci].add_counts(count, low);
            }
        } else {
            self.trans_combo_data[ci].add_counts(count, low);
        }
    }

    /// processReadRecords (tests) in 2.0.1: the reads of a gene collection without a BAM
    #[cfg(test)]
    pub fn process_read_records_201(&mut self, gc: &mut GeneCollection, reads: Vec<Read>, refg: &mut RefSource) {
        self.base_depth.initialise(gc.region_bounds);
        self.valid_read_start = [if gc.non_genic_positions[0] >= 0 { gc.non_genic_positions[0] } else { gc.region_bounds[0] }, gc.region_bounds[1]];
        self.set_expression_gene_data();
        self.alt_sj.clear();
        self.retained_introns.introns.clear();
        self.chimeric.initialise(gc);
        let n = reads.len();
        for (i, mut r) in reads.into_iter().enumerate() {
            if n == 2 && i == 0 {
                r.flags |= crate::read::FLAG_FIRST;
            }
            self.reads.push(r);
            let ri = self.reads.len() - 1;
            self.process_read_201(ri, gc, refg);
        }
    }
}
