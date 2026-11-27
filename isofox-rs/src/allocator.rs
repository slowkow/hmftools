//! FragmentAllocator and ExpressionReadTracker: the reads of a gene collection, paired into fragments and
//! assigned to transcripts, unspliced genes, novel junctions, retained introns and chimeric candidates.

use crate::adjusts::GcRatioCounts;
use crate::bam::{aux_string, rec_mate_start, rec_start, to_read, BamReader, SliceFilter};
use crate::chimeric::{is_realigned_fragment_candidate, ChimericReadTracker};
use crate::common::*;
use crate::config::{IsofoxConfig, IsofoxFunction};
use crate::expression::CategoryCountsData;
use crate::gene::GeneCollection;
use crate::novel::{AltSjCohortCache, AltSpliceJunctionFinder, RetainedIntronFinder};
use crate::fragment;
use crate::read::{consensus_duplicate_count, trim_adapter_bases, Read};
use crate::read_trans::{calc_fragment_length, find_overlapping_regions, mark_region_bases, process_overlapping_regions};
use crate::refgen::RefSource;
use anyhow::Result;
use hmf_common::java_seq_hash::SeqHashMap;
use hmf_common::regions::{BaseRegion, ChrBaseRegion};
use rust_htslib::bam;
use std::collections::{HashMap, HashSet};

const NON_GENIC_BASE_DEPTH_WIDTH: i32 = 250_000;

pub struct FragmentAllocator<'a> {
    config: &'a IsofoxConfig,
    reader: Option<BamReader>,
    filter: SliceFilter,
    pub reads: Vec<Read>,
    tracker: SeqHashMap<String, usize>,
    gene_read_count: i32,
    pub total_bam_read_count: i32,
    valid_read_start: [i32; 2],
    pub base_depth: BaseDepth,
    run_fusions: bool,
    fusions_only: bool,
    stats_only: bool,
    expression_enabled: bool,
    pub gc_ratio_counts: Option<GcRatioCounts>,
    pub gene_gc_ratio_counts: Option<GcRatioCounts>,
    pub trans_combo_data: Vec<CategoryCountsData>,
    pub chimeric: ChimericReadTracker,
    pub alt_sj: AltSpliceJunctionFinder<'a>,
    pub retained_introns: RetainedIntronFinder,
    immune: Vec<ChrBaseRegion>,
    known_splice_sites: &'a HashMap<String, HashSet<i32>>,
    halted: bool,
}

impl<'a> FragmentAllocator<'a> {
    pub fn new(
        config: &'a IsofoxConfig,
        alt_sj_cohort: Option<&'a AltSjCohortCache>,
        known_splice_sites: &'a HashMap<String, HashSet<i32>>,
    ) -> Result<FragmentAllocator<'a>> {
        let reader = match &config.bam_file {
            Some(b) => Some(BamReader::open(b, config.ref_genome_file.as_deref())?),
            None => None,
        };
        let keep_supps = config.run_function(IsofoxFunction::Fusions) || config.run_function(IsofoxFunction::AltSpliceJunctions);
        let gc = config.require_gc_ratio_calcs();
        Ok(FragmentAllocator {
            config,
            reader,
            filter: SliceFilter { min_map_quality: 0, keep_duplicates: false, keep_supplementaries: keep_supps, keep_secondaries: false },
            reads: Vec::new(),
            tracker: SeqHashMap::new(),
            gene_read_count: 0,
            total_bam_read_count: 0,
            valid_read_start: [0, 0],
            base_depth: BaseDepth::new(),
            run_fusions: config.run_function(IsofoxFunction::Fusions),
            fusions_only: config.run_fusions_only(),
            stats_only: config.run_statistics_only(),
            expression_enabled: config.run_function(IsofoxFunction::TranscriptCounts),
            gc_ratio_counts: if gc { Some(GcRatioCounts::new()) } else { None },
            gene_gc_ratio_counts: if gc { Some(GcRatioCounts::new()) } else { None },
            trans_combo_data: Vec::new(),
            chimeric: ChimericReadTracker::new(config),
            alt_sj: AltSpliceJunctionFinder::new(config.run_function(IsofoxFunction::AltSpliceJunctions), alt_sj_cohort, &config.fusions.known_fusions),
            retained_introns: RetainedIntronFinder { enabled: config.run_function(IsofoxFunction::RetainedIntrons), introns: Vec::new() },
            immune: config.filters.immune_gene_regions.clone(),
            known_splice_sites,
            halted: false,
        })
    }

    pub fn clear_cache(&mut self) {
        self.tracker.clear();
        self.chimeric.clear_data();
        self.set_expression_gene_data();
        self.alt_sj.clear();
        self.retained_introns.introns.clear();
        self.reads.clear();
    }

    fn set_expression_gene_data(&mut self) {
        if let Some(g) = &mut self.gene_gc_ratio_counts {
            g.clear_counts();
        }
        self.trans_combo_data.clear();
    }

    /// processBam
    pub fn process_bam(&mut self, gc: &mut GeneCollection, region: [i32; 2], refg: &mut RefSource) -> Result<()> {
        self.clear_cache();
        self.gene_read_count = 0;
        let depth_range = [region[0].max(gc.region_bounds[0] - NON_GENIC_BASE_DEPTH_WIDTH), region[1].min(gc.region_bounds[1] + NON_GENIC_BASE_DEPTH_WIDTH)];
        self.base_depth.initialise(depth_range);
        self.chimeric.initialise(gc);
        if self.expression_enabled {
            self.set_expression_gene_data();
        }
        self.valid_read_start = region;
        self.halted = false;

        let Some(mut reader) = self.reader.take() else { return Ok(()) };
        let filter = self.filter;
        let chromosome = gc.chromosome;
        let res = reader.slice(chromosome, region[0], region[1], &filter, |rec: &bam::Record| {
            self.process_sam_record(rec, gc, refg);
            !self.halted
        });
        self.reader = Some(reader);
        res?;

        if self.chimeric.enabled {
            self.process_incomplete_reads(gc, refg);
        }
        Ok(())
    }

    fn process_sam_record(&mut self, rec: &bam::Record, gc: &mut GeneCollection, refg: &mut RefSource) {
        let start = rec_start(rec);
        if !position_within(start, self.valid_read_start[0], self.valid_read_start[1]) {
            return;
        }
        if self.skip_record(rec) {
            self.chimeric.stats.excluded += 1;
            return;
        }
        // trackFragmentCounts
        self.total_bam_read_count += 1;
        self.gene_read_count += 1;
        let cr = aux_string(rec, b"CR");
        if !rec.is_supplementary() && (!rec.is_paired() || rec.is_first_in_template()) {
            gc.add_count(FragmentType::Total, 1);
            if cr.is_some() {
                let d = consensus_duplicate_count(cr.as_deref()) as i64;
                gc.add_count(FragmentType::Total, d);
                gc.add_count(FragmentType::Duplicate, d);
            }
        }
        let read = to_read(rec);
        if self.config.log_read_ids.iter().any(|x| **x == *read.id) {
            hmf_common::debug!("specific read: {} range({}: {} -> {}) flags({}) cigar({})", read.id, read.chromosome, read.pos_start, read.pos_end, read.flags, read.cigar_str());
        }
        self.reads.push(read);
        let ri = self.reads.len() - 1;
        self.process_read(ri, gc, refg);
    }

    /// GeneRegionFilters.skipRead(record, true)
    fn skip_record(&self, rec: &bam::Record) -> bool {
        let f = &self.config.filters;
        let chr = crate::bam::contig_name(rec.tid());
        if f.skip_read(chr, rec_start(rec), crate::bam::rec_end(rec), false) {
            return true;
        }
        if rec.is_paired() && !rec.is_mate_unmapped() {
            let ms = rec_mate_start(rec);
            if f.skip_read(crate::bam::contig_name(rec.mtid()), ms, ms, true) {
                return true;
            }
        }
        if let Some(sa) = aux_string(rec, b"SA") {
            if let Some(s) = crate::read::SuppData::extract(&sa) {
                if f.skip_read(&s.chromosome, s.position, s.position, true) {
                    return true;
                }
            }
        }
        false
    }

    fn process_read(&mut self, ri: usize, gc: &mut GeneCollection, refg: &mut RefSource) {
        if self.config.gene_read_limit > 0 && self.gene_read_count >= self.config.gene_read_limit {
            self.halted = true;
            hmf_common::warn!("chr({}) genes({}) readCount({}) exceeds max read count", gc.chromosome, gc.gene_names(10), self.gene_read_count);
            return;
        }
        let vr = self.valid_read_start;
        gc.set_read_gene_collections(&mut self.reads[ri], vr);
        if self.reads[ri].is_supplementary() {
            self.handle_supplementary_read(ri, gc);
        } else if !self.reads[ri].is_read_paired() {
            self.process_fragment_reads(&[ri], gc, refg);
        } else {
            let id = self.reads[ri].id.to_string();
            match self.tracker.remove(&id) {
                Some(other) => self.process_fragment_reads(&[ri, other], gc, refg),
                None => {
                    self.tracker.insert(id, ri);
                }
            }
        }
    }

    fn mark_gene_data_regions(&mut self, ri: usize, gc: &GeneCollection) {
        let overlapping = find_overlapping_regions(&gc.exon_regions, None, &self.reads[ri]);
        if !overlapping.is_empty() {
            process_overlapping_regions(&mut self.reads[ri], &gc.exon_regions, &overlapping);
        }
    }

    fn handle_supplementary_read(&mut self, ri: usize, gc: &GeneCollection) {
        if self.reads[ri].is_mate_unmapped() {
            return;
        }
        let coords = self.reads[ri].mapped_coords.alignments.clone();
        self.base_depth.process_read(&coords);
        self.mark_gene_data_regions(ri, gc);
        let immune = std::mem::take(&mut self.immune);
        self.chimeric.add_supplementary_read(ri, &mut self.reads, gc, &immune);
        self.immune = immune;
    }

    fn frag_reads(&self, frag: &[usize]) -> Vec<&Read> {
        frag.iter().map(|r| &self.reads[*r]).collect()
    }

    fn frag_reads_mut(&mut self, frag: &[usize]) -> Vec<&mut Read> {
        if frag.len() == 1 {
            return vec![&mut self.reads[frag[0]]];
        }
        let (a, b) = self.two_reads(frag[0], frag[1]);
        vec![a, b]
    }

    fn two_reads(&mut self, a: usize, b: usize) -> (&mut Read, &mut Read) {
        assert!(a != b);
        if a < b {
            let (x, y) = self.reads.split_at_mut(b);
            (&mut x[a], &mut y[0])
        } else {
            let (x, y) = self.reads.split_at_mut(a);
            (&mut y[0], &mut x[b])
        }
    }


    fn process_fragment_reads(&mut self, frag: &[usize], gc: &mut GeneCollection, refg: &mut RefSource) {
        if frag.len() == 2 {
            let (r1, r2) = self.two_reads(frag[0], frag[1]);
            trim_adapter_bases(r1, r2);
        }
        for &ri in frag {
            self.mark_gene_data_regions(ri, gc);
        }
        let fragment_count = fragment::fragment_count(&self.frag_reads(frag));
        let num_loci = fragment::min_num_loci(&self.frag_reads(frag));
        let multi_mapped = num_loci > 1;
        let is_chimeric = self.chimeric.is_chimeric(frag, &mut self.reads, gc, multi_mapped);

        if self.stats_only {
            if is_chimeric {
                gc.add_count(FragmentType::Chimeric, 1);
            } else if fragment::is_fully_intronic(&self.frag_reads(frag)) {
                gc.add_count(FragmentType::Unspliced, fragment_count as i64);
            } else {
                gc.add_count(FragmentType::TransSupporting, fragment_count as i64);
            }
            return;
        }

        let common = fragment::merged_mappings(&self.frag_reads(frag));
        self.base_depth.process_read(&common);

        if is_chimeric {
            if !multi_mapped {
                if self.chimeric.enabled {
                    let immune = std::mem::take(&mut self.immune);
                    self.chimeric.add_chimeric_fragment(frag, &mut self.reads, gc, &immune);
                    self.immune = immune;
                } else {
                    gc.add_count(FragmentType::Chimeric, 1);
                }
            }
            return;
        }

        if self.run_fusions {
            if frag.iter().any(|r| is_realigned_fragment_candidate(&self.reads[*r])) {
                self.chimeric.add_realignment_candidates(frag, &self.reads);
            }
            if self.fusions_only {
                return;
            }
        }

        let min_pos = fragment::min_alignment_start(&self.frag_reads(frag));
        let max_pos = fragment::max_alignment_end(&self.frag_reads(frag));
        let overlap_genes = gc.find_genes_covering_range(min_pos, max_pos, true);
        let contains_split = fragment::contains_split(&self.frag_reads(frag));

        if fragment::is_fully_intronic(&self.frag_reads(frag)) {
            // processIntronicReads
            if contains_split {
                gc.add_count(FragmentType::Alt, 1);
                if self.alt_sj.enabled {
                    let reads: Vec<&Read> = frag.iter().map(|r| &self.reads[*r]).collect();
                    self.alt_sj.evaluate_fragment_reads(gc, refg, &overlap_genes, &reads, &[]);
                }
                return;
            }
            if self.expression_enabled && !overlap_genes.is_empty() {
                let genes: Vec<String> = overlap_genes.iter().map(|g| gc.genes[*g].gene.gene_id.clone()).collect();
                let ci = self.category_counts(&[], &genes);
                self.add_gc_counts(ci, &common, fragment_count, multi_mapped, gc, refg);
            }
            gc.add_count(FragmentType::Unspliced, fragment_count as i64);
            return;
        }

        let valid_regions = fragment::unique_valid_regions(&self.frag_reads(frag));
        // validTypeTranscripts, then the fragment length check
        let mut valid_trans: Vec<i32> = Vec::new();
        for t in fragment::valid_type_transcripts(&self.frag_reads(frag)) {
            if let Some(td) = gc.transcript(t) {
                let len = calc_fragment_length(td, min_pos, max_pos);
                if len > 0 && len <= self.config.max_fragment_length {
                    valid_trans.push(t);
                }
            }
        }
        let invalid = fragment::invalid_transcripts(&self.frag_reads(frag), &valid_trans);

        let mut fragment_type = FragmentType::Unspliced;
        if valid_trans.is_empty() {
            let check_retained;
            if contains_split {
                fragment_type = FragmentType::Alt;
                if self.alt_sj.enabled {
                    let reads: Vec<&Read> = frag.iter().map(|r| &self.reads[*r]).collect();
                    let related = invalid.keys();
                    self.alt_sj.evaluate_fragment_reads(gc, refg, &overlap_genes, &reads, &related);
                }
                check_retained = true;
            } else {
                for t in invalid.keys() {
                    if fragment::spans_multiple_regions(&self.frag_reads(frag), &gc.exon_regions, t) {
                        fragment_type = FragmentType::Alt;
                        break;
                    }
                }
                check_retained = true;
            }
            if check_retained && self.retained_introns.enabled {
                let reads: Vec<&Read> = frag.iter().map(|r| &self.reads[*r]).collect();
                self.retained_introns.evaluate_fragment_reads(gc, &reads);
            }
            if fragment_type == FragmentType::Unspliced && self.expression_enabled && !overlap_genes.is_empty() {
                let genes: Vec<String> = overlap_genes.iter().map(|g| gc.genes[*g].gene.gene_id.clone()).collect();
                let ci = self.category_counts(&valid_trans, &genes);
                self.add_gc_counts(ci, &common, fragment_count, multi_mapped, gc, refg);
            }
        } else {
            fragment_type = FragmentType::TransSupporting;
            {
                let mut rs = self.frag_reads_mut(frag);
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
                    for &ri in frag {
                        gene_is_forward = self.find_gene_strand(ri, gc, &valid_trans);
                        if gene_is_forward.is_some() {
                            break;
                        }
                    }
                }
                let mt = if fragment::has_transcript_classification(&self.frag_reads(frag), t, TransMatchType::SpliceJunction) {
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
                if !contains_split && fragment::reads_in_different_exons(&self.frag_reads(frag), &gc.exon_regions, t) {
                    gc.add_transcript_read_match_type(t, FragmentMatchType::Discordant);
                }
                self.process_valid_transcript(t, frag, gc, is_unique);
            }
            if self.expression_enabled {
                let genes: Vec<String> = if combo == FragmentMatchType::Short { overlap_genes.iter().map(|g| gc.genes[*g].gene.gene_id.clone()).collect() } else { Vec::new() };
                let ci = self.category_counts(&valid_trans, &genes);
                self.add_gc_counts(ci, &common, fragment_count, multi_mapped, gc, refg);
            }
            if let Some(fwd) = gene_is_forward {
                if let Some(o) = fragment::orientation(&self.frag_reads(frag)) {
                    if (o == ORIENT_FWD) == fwd {
                        gc.add_count(FragmentType::ForwardStrand, fragment_count as i64);
                    } else {
                        gc.add_count(FragmentType::ReverseStrand, fragment_count as i64);
                    }
                }
            }
        }
        gc.add_count(fragment_type, fragment_count as i64);
    }


    fn find_gene_strand(&self, ri: usize, gc: &GeneCollection, transcripts: &[i32]) -> Option<bool> {
        let keys = self.reads[ri].mapped_regions.keys();
        for &t in transcripts {
            for &reg in &keys {
                if let Some(te) = gc.exon_regions[reg].trans_exon_refs.iter().find(|x| x.trans_id == t) {
                    if let Some(g) = gc.genes.iter().find(|g| g.gene.gene_id == te.gene_id) {
                        return Some(g.gene.strand == ORIENT_FWD);
                    }
                }
            }
        }
        None
    }



    /// ExpressionReadTracker.processValidTranscript
    fn process_valid_transcript(&self, t: i32, frag: &[usize], gc: &mut GeneCollection, unique: bool) {
        let mut processed: Vec<usize> = Vec::new();
        for &ri in frag {
            let read = &self.reads[ri];
            let regions: Vec<usize> = read
                .mapped_regions
                .entries()
                .into_iter()
                .filter(|(reg, mt)| gc.exon_regions[*reg].has_trans_id(t) && mt.valid_exon_match())
                .map(|(reg, _)| reg)
                .collect();
            for &reg in &regions {
                if !processed.contains(&reg) {
                    gc.exon_regions[reg].add_transcript_read_match(t, unique);
                }
            }
            if regions.len() > 1 && read.trans_classification(t) == TransMatchType::SpliceJunction {
                for a in 0..regions.len() - 1 {
                    let r1 = regions[a];
                    for &r2 in &regions[a + 1..] {
                        if processed.contains(&r1) && processed.contains(&r2) {
                            continue;
                        }
                        if gc.exon_regions[r1].post_regions.contains(&r2) {
                            gc.exon_regions[r1].add_transcript_junction_match(t, SE_END, unique);
                            gc.exon_regions[r2].add_transcript_junction_match(t, SE_START, unique);
                        } else if gc.exon_regions[r1].pre_regions.contains(&r2) {
                            gc.exon_regions[r1].add_transcript_junction_match(t, SE_START, unique);
                            gc.exon_regions[r2].add_transcript_junction_match(t, SE_END, unique);
                        }
                    }
                }
            }
            processed.extend(regions);
        }
    }

    /// ExpressionReadTracker.getCategoryCountsData: the index of the category
    fn category_counts(&mut self, transcripts: &[i32], genes: &[String]) -> usize {
        let gene_refs: Vec<&str> = genes.iter().map(|s| s.as_str()).collect();
        if let Some(i) = self.trans_combo_data.iter().position(|c| c.matches(transcripts, &gene_refs)) {
            return i;
        }
        let mut c = CategoryCountsData::new(transcripts.to_vec(), genes.to_vec());
        if let Some(g) = &self.gc_ratio_counts {
            if self.config.apply_gc_bias_adjust() {
                c.counts_by_gc = Some(vec![0.0; g.counts.len()]);
            }
        }
        self.trans_combo_data.push(c);
        self.trans_combo_data.len() - 1
    }

    /// ExpressionReadTracker.addGcCounts
    fn add_gc_counts(&mut self, ci: usize, regions: &[BaseRegion], count: i32, multi_mapped: bool, gc: &mut GeneCollection, refg: &mut RefSource) {
        let mut idx = [-1, -1];
        let mut counts = [0.0, 0.0];
        if let Some(g) = &self.gc_ratio_counts {
            let ratio = refg.gc_ratio_of_regions(gc.chromosome, regions);
            g.determine_ratio_data(ratio, &mut idx, &mut counts);
        }
        if multi_mapped {
            gc.add_count(FragmentType::MultiMapped, 1);
        }
        if self.gc_ratio_counts.is_some() {
            for i in 0..2 {
                if idx[i] >= 0 {
                    self.gc_ratio_counts.as_mut().unwrap().add_gc_ratio_count(idx[i], counts[i]);
                    self.gene_gc_ratio_counts.as_mut().unwrap().add_gc_ratio_count(idx[i], counts[i]);
                }
            }
            if self.config.apply_gc_bias_adjust() {
                self.trans_combo_data[ci].add_gc_ratio_counts(count as f64, multi_mapped, &idx, &counts);
            } else {
                self.trans_combo_data[ci].add_counts(count as f64, multi_mapped);
            }
        } else {
            self.trans_combo_data[ci].add_counts(count as f64, multi_mapped);
        }
    }

    /// processIncompleteReads: the reads whose mate was not found in this gene collection
    fn process_incomplete_reads(&mut self, gc: &mut GeneCollection, refg: &mut RefSource) {
        let lone: Vec<usize> = self.tracker.values().copied().collect();
        for &ri in &lone {
            self.mark_gene_data_regions(ri, gc);
            let coords = self.reads[ri].mapped_coords.alignments.clone();
            self.base_depth.process_read(&coords);
        }
        let immune = std::mem::take(&mut self.immune);
        self.chimeric.post_process_chimeric_reads(&mut self.base_depth, &lone, &mut self.reads, gc, &immune, self.known_splice_sites);
        self.immune = immune;
        self.process_chimeric_novel_junctions(gc, refg);
    }

    fn process_chimeric_novel_junctions(&mut self, gc: &mut GeneCollection, refg: &mut RefSource) {
        if !self.alt_sj.enabled || self.chimeric.local_chimeric_reads.is_empty() {
            return;
        }
        let groups = self.chimeric.local_chimeric_reads.clone();
        for reads in groups {
            let (r1, r2) = if reads.len() == 2 {
                (Some(reads[0]), Some(reads[1]))
            } else if reads.len() == 3 {
                let mut a = None;
                let mut b = None;
                for r in &reads {
                    if self.reads[*r].has_supp_alignment() {
                        if a.is_none() {
                            a = Some(*r);
                        } else {
                            b = Some(*r);
                            break;
                        }
                    }
                }
                (a, b)
            } else {
                (None, None)
            };
            let (Some(r1), Some(r2)) = (r1, r2) else { continue };
            let (a, b) = (&self.reads[r1], &self.reads[r2]);
            let min_pos = a.pos_start.min(b.pos_start);
            let max_pos = a.pos_end.max(b.pos_end);
            let genes = gc.find_genes_covering_range(min_pos, max_pos, false);
            let rs = [a, b];
            self.alt_sj.evaluate_fragment_reads(gc, refg, &genes, &rs, &[]);
        }
    }

    /// annotateNovelLocations
    pub fn annotate_novel_locations(&mut self, gc: &GeneCollection, refg: &mut RefSource, alt_sj_unfiltered: &mut Vec<String>, alt_sj_pass: &mut Vec<String>, retained: &mut Vec<String>) {
        if !(self.alt_sj.junctions.is_empty() && self.retained_introns.introns.is_empty()) {
            if self.alt_sj.enabled {
                self.alt_sj.set_position_depth(&self.base_depth);
            }
            if self.retained_introns.enabled {
                self.retained_introns.set_position_depth(&self.base_depth);
            }
        }
        if self.alt_sj.enabled {
            self.alt_sj.prioritise_genes(gc, refg);
            self.alt_sj.finalise(gc);
            self.alt_sj.write(alt_sj_unfiltered, alt_sj_pass);
        }
        if self.retained_introns.enabled {
            self.retained_introns.write(gc, retained);
        }
    }

    /// processReadRecords (tests): the reads of a gene collection without a BAM
    #[cfg(test)]
    pub fn process_read_records(&mut self, gc: &mut GeneCollection, reads: Vec<Read>, refg: &mut RefSource) {
        self.base_depth.initialise(gc.region_bounds);
        self.valid_read_start = [if gc.non_genic_positions[0] >= 0 { gc.non_genic_positions[0] } else { gc.region_bounds[0] }, gc.region_bounds[1]];
        self.set_expression_gene_data();
        self.alt_sj.clear();
        self.retained_introns.introns.clear();
        self.chimeric.initialise(gc);
        for r in reads {
            self.reads.push(r);
            let ri = self.reads.len() - 1;
            self.process_read(ri, gc, refg);
        }
    }

}
