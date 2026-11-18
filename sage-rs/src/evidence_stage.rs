//! ReadContextEvidence: gathers evidence for a list of candidates from one sample's BAM.

use crate::bamio::BamSlicer;
use crate::candidate::{calc_adjusted_num_mutations, Candidate};
use crate::config::Config;
use crate::evidence::*;
use crate::quality::*;
use crate::refdata::RefData;
use crate::refseq::{RefGenome, RefSequence};
use crate::sync::{FragmentSync, SyncOutput};
use crate::types::*;
use anyhow::Result;

const SLICE_SOFT_CLIP_BUFFER: i32 = 30;

pub struct Sample<'a> {
    pub bam: &'a mut BamSlicer,
    pub bqr: &'a BqrMap,
    pub msi: &'a MsiJitter,
    pub is_reference: bool,
}

/// Receives the counters supporting / contradicting each read, for phasing.
pub trait PhaseSink {
    fn register(&mut self, counters: &[ReadContextCounter], positive: &[usize], negative: &[usize]);
}

/// buildCandidateRegions
fn build_candidate_regions(cfg: &Config, candidates: &[Candidate]) -> Vec<BaseRegion> {
    let first = candidates[0].position();
    let last = candidates[candidates.len() - 1].position();
    let min_gap = (cfg.read_length * 2).max(100);
    let average_gap = (last - first) / candidates.len() as i32;
    if average_gap <= min_gap || cfg.max_partition_slices == 1 {
        return vec![BaseRegion::new(first, last)];
    }
    let mut gaps: Vec<i32> = candidates.windows(2).map(|w| w[1].position() - w[0].position()).collect();
    gaps.sort_by(|a, b| b.cmp(a));
    let nth = (cfg.max_partition_slices as usize).min(gaps.len());
    let nth_gap = gaps[nth - 1];
    let mut regions = vec![BaseRegion::new(first, first)];
    for c in &candidates[1..] {
        let cur = regions.last_mut().unwrap();
        if c.position() - cur.end <= min_gap.max(nth_gap) {
            cur.end = c.position();
        } else {
            regions.push(BaseRegion::new(c.position(), c.position()));
        }
    }
    regions
}

fn consider_read(rc: &ReadContextCounter, unclipped_start: i32, unclipped_end: i32) -> bool {
    if rc.exceeds_max_coverage() {
        return false;
    }
    if rc.variant().is_delete() {
        return positions_overlap(rc.position(), rc.max_position_vs_read_start, unclipped_start, unclipped_end);
    }
    position_within(rc.position(), unclipped_start, unclipped_end)
}

struct Collector<'a> {
    counters: Vec<ReadContextCounter>,
    last_index: i64,
    selected: Vec<usize>,
    ref_seq: &'a RefSequence,
    qc: &'a QualityCalculator<'a>,
}

impl<'a> Collector<'a> {
    fn process(&mut self, read: &Read, fragment: Option<&FragmentData>, phaser: &mut Option<&mut dyn PhaseSink>) {
        self.selected.clear();
        let n = self.counters.len() as i64;
        let mut read_start = read.start - read.left_soft_clip();
        let mut read_end = read.end + read.right_soft_clip();
        let mut next_index = self.last_index + 1;
        let mut prev_index = self.last_index;
        while prev_index >= 0 && prev_index < n {
            let rc = &self.counters[prev_index as usize];
            if consider_read(rc, read_start, read_end) {
                self.last_index = prev_index;
                self.selected.insert(0, prev_index as usize);
                read_start -= rc.max_candidate_delete_length;
                read_end += rc.max_candidate_delete_length;
            } else if rc.position() < read_start {
                break;
            }
            prev_index -= 1;
        }
        while next_index < n {
            let rc = &self.counters[next_index as usize];
            if consider_read(rc, read_start, read_end) {
                self.selected.push(next_index as usize);
            } else if rc.position() < read_start {
                self.last_index = next_index;
            } else if rc.position() > read_end {
                break;
            }
            next_index += 1;
        }
        if self.selected.is_empty() {
            return;
        }
        let number_of_events = calc_adjusted_num_mutations(read, self.ref_seq);
        let mut positive: Vec<usize> = Vec::new();
        let mut negative: Vec<usize> = Vec::new();
        for i in 0..self.selected.len() {
            let idx = self.selected[i];
            let m = self.counters[idx].process_read(read, number_of_events, fragment, self.qc);
            match m {
                ReadMatchType::AltSupportExact => positive.push(idx),
                ReadMatchType::RefSupport | ReadMatchType::NoAltRefMatch => negative.push(idx),
                _ => {}
            }
        }
        if let Some(p) = phaser.as_mut() {
            p.register(&self.counters, &positive, &negative);
        }
    }
}

/// ReadContextEvidence.collectEvidence
#[allow(clippy::too_many_arguments)]
pub fn collect_evidence(
    cfg: &Config, chr: &str, candidates: &[Candidate], genome: &RefGenome, ref_data: &RefData, sample: &mut Sample,
    sync: &mut FragmentSync, mut phaser: Option<&mut dyn PhaseSink>,
) -> Result<Vec<ReadContextCounter>> {
    if candidates.is_empty() {
        return Ok(Vec::new());
    }
    let mut slices = build_candidate_regions(cfg, candidates);
    for s in slices.iter_mut() {
        s.start = (s.start - SLICE_SOFT_CLIP_BUFFER).max(1);
        s.end += SLICE_SOFT_CLIP_BUFFER;
    }
    let bounds = BaseRegion::new(slices[0].start, slices[slices.len() - 1].end);
    let ref_seq = RefSequence::new(genome, chr, bounds.start, bounds.end)?;
    let qc = QualityCalculator { cfg, bqr: sample.bqr, ref_seq: &ref_seq, msi: sample.msi, ref_data, chr };

    let mut counters: Vec<ReadContextCounter> = candidates
        .iter()
        .map(|c| {
            ReadContextCounter::new(c.ctx.clone(), c.tier, max_coverage(cfg, c.tier), c.min_events, &qc, sample.is_reference)
        })
        .collect();

    let deletes: Vec<(i32, i32)> =
        candidates.iter().filter(|c| c.variant().is_delete()).map(|c| (c.position(), c.variant().ref_len() - 1)).collect();
    for rc in counters.iter_mut() {
        let max_close_del =
            deletes.iter().filter(|(p, _)| (p - rc.position()).abs() < 50).map(|(_, l)| *l).max().unwrap_or(0);
        if max_close_del >= 5 {
            rc.max_candidate_delete_length = max_close_del;
        }
    }

    let mut collector = Collector { counters, last_index: 0, selected: Vec::new(), ref_seq: &ref_seq, qc: &qc };
    let mut out: SyncOutput = Vec::new();
    let sync_fragments = cfg.sync_fragments;
    sample.bam.slice(chr, &slices, false, |read| {
        if sync_fragments {
            match sync.handle(read, &mut out) {
                Some(read) => collector.process(&read, None, &mut phaser),
                None => {
                    for (r, f) in out.drain(..) {
                        collector.process(&r, f.as_ref(), &mut phaser);
                    }
                }
            }
        } else {
            collector.process(&read, None, &mut phaser);
        }
    })?;
    for read in sync.empty_cached_reads() {
        collector.process(&read, None, &mut phaser);
    }

    let mut counters = collector.counters;
    for rc in counters.iter_mut() {
        set_jitter_qual_filter_state(sample.msi, rc);
    }
    Ok(counters)
}

/// VariantFilters.passesHardFilters
pub fn passes_hard_filters(cfg: &Config, rc: &ReadContextCounter) -> bool {
    if rc.tier == Tier::Hotspot {
        return true;
    }
    if rc.alt_support() < cfg.hard_min_tumor_raw_alt_support {
        return false;
    }
    if rc.tumor_quality() < cfg.hard_min_tumor_qual as f64 {
        return false;
    }
    if rc.vaf() < cfg.hard_min_tumor_vaf {
        return false;
    }
    if rc.jitter.hard_filter_on_noise {
        return false;
    }
    true
}
