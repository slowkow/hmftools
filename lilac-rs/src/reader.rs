//! BamRecordReader and GeneCodingRegions: the fragments of each HLA gene's coding region.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use anyhow::Result;

use crate::bam::{filter_low_qual_read, Bam, Read, ReadCtx, Rec};
use crate::constants::MAX_DISTANCE;
use crate::fragment::{merge_fragments, read_id_map_order, Factory, Fragment};
use crate::hla::HlaGene;
use crate::indel::Indel;
use crate::refdata::GeneCache;

/// GeneCodingRegions
#[derive(Clone, Debug)]
pub struct GeneRegions {
    pub gene: HlaGene,
    pub chromosome: String,
    pub strand: i8,
    pub coding_start: i32,
    pub coding_end: i32,
    pub coding_regions: Vec<(i32, i32)>,
}

impl GeneRegions {
    pub fn within_coding_bounds(&self, pos: i32, buffer: i32) -> bool {
        pos >= self.coding_start - buffer && pos <= self.coding_end + buffer
    }
}

/// The state of a BamRecordReader: counters and the known stop-loss fragments
pub struct GeneReader {
    pub regions: Vec<GeneRegions>,
    pub filtered_count: i32,
    /// fragment uids with the known stop-loss indel (mKnownStopLossFragments)
    pub stop_loss_frags: Vec<u64>,
    pub has_stop_loss: bool,
    unmatched_indels: Vec<(Indel, i32)>,
    unmatched_pon: Vec<(Indel, i32)>,
    discard_ids: HashSet<String>,
}

fn increment(v: &mut Vec<(Indel, i32)>, i: &Indel) {
    match v.iter_mut().find(|x| x.0 == *i) {
        Some(x) => x.1 += 1,
        None => v.push((i.clone(), 1)),
    }
}

impl GeneReader {
    pub fn new(gc: &GeneCache, hla_chr: &str) -> GeneReader {
        let regions = gc
            .gene_names
            .iter()
            .zip(gc.transcripts.iter())
            .map(|(g, t)| GeneRegions {
                gene: *g,
                chromosome: hla_chr.to_string(),
                strand: t.strand,
                coding_start: t.coding_start,
                coding_end: t.coding_end,
                coding_regions: t
                    .exons
                    .iter()
                    .filter(|e| t.coding_start <= e.end && t.coding_end >= e.start)
                    .map(|e| (t.coding_start.max(e.start), t.coding_end.min(e.end)))
                    .collect(),
            })
            .collect();
        GeneReader {
            regions,
            filtered_count: 0,
            stop_loss_frags: Vec::new(),
            has_stop_loss: false,
            unmatched_indels: Vec::new(),
            unmatched_pon: Vec::new(),
            discard_ids: HashSet::new(),
        }
    }

    /// bothEndsInRangeOfCodingTranscripts
    pub fn both_ends_in_range(&self, r: &Rec, hla_chr: &str) -> bool {
        let read_in = self.regions.iter().any(|x| x.within_coding_bounds(r.start, MAX_DISTANCE));
        if !r.paired() {
            return read_in;
        }
        if r.mate_contig != hla_chr {
            return false;
        }
        let mate_in = self.regions.iter().any(|x| x.within_coding_bounds(r.mate_start, MAX_DISTANCE));
        read_in && mate_in
    }

    /// findGeneFragments: the fragments of every gene, gene by gene
    pub fn find_gene_fragments(&mut self, bam: &mut Bam, factory: &Factory, ctx: &ReadCtx, stop_loss: Option<&Indel>, pon: &[Indel]) -> Result<Vec<Fragment>> {
        let mut out = Vec::new();
        for gi in 0..self.regions.len() {
            let g = self.regions[gi].clone();
            let records = bam.slice(&g.chromosome, g.coding_start, g.coding_end)?;
            let mut reads: Vec<Read> = Vec::new();
            for rec in records {
                if filter_low_qual_read(&rec, ctx.params) {
                    continue;
                }
                let overlaps: Vec<(i32, i32)> = g.coding_regions.iter().copied().filter(|x| x.0 <= rec.end && x.1 >= rec.start).collect();
                if overlaps.is_empty() {
                    continue;
                }
                if !self.both_ends_in_range(&rec, &ctx.params.hla_chr) {
                    crate::debug!("filter read: id({}) coords({}-{})", rec.name, rec.start, rec.end);
                    self.filtered_count += 1;
                } else {
                    for o in overlaps {
                        reads.push(Read::create(ctx, o, &rec, true, true)?);
                    }
                }
            }
            let frags = self.create_fragments(g.gene, g.strand, reads, factory, stop_loss, pon)?;
            out.extend(frags.into_iter().filter(|f| !self.discard_ids.contains(f.id())));
        }
        Ok(out)
    }

    /// createFragments: one fragment per read id, in the HashMap's order of ids
    fn create_fragments(&mut self, gene: HlaGene, strand: i8, reads: Vec<Read>, factory: &Factory, stop_loss: Option<&Indel>, pon: &[Indel]) -> Result<Vec<Fragment>> {
        let mut frags: Vec<Fragment> = Vec::new();
        let mut index: HashMap<Arc<str>, usize> = HashMap::new();
        for read in reads {
            let fragment = match factory.create_fragment(&read, gene, strand) {
                Ok(f) => f,
                Err(e) => {
                    crate::error!("read({} {}) failed to build fragment: {e}", read.id, read.read_info());
                    None
                }
            };
            if let Some(f) = fragment {
                let pos = match index.get(&read.id) {
                    None => {
                        index.insert(read.id.clone(), frags.len());
                        frags.push(f);
                        frags.len() - 1
                    }
                    Some(p) => {
                        let p = *p;
                        merge_fragments(&mut frags[p], &f);
                        p
                    }
                };
                if stop_loss.is_some_and(|s| read.valid_indels.contains(s)) {
                    self.has_stop_loss = true;
                    let uid = frags[pos].uid;
                    if !self.stop_loss_frags.contains(&uid) {
                        self.stop_loss_frags.push(uid);
                    }
                }
                continue;
            }
            if read.contains_valid_indel() {
                if stop_loss.is_some_and(|s| read.valid_indels.contains(s)) {
                    crate::debug!("missing known indel fragment: {} {}", read.id, read.read_info());
                }
                self.discard_ids.insert(read.id.to_string());
            }
            for indel in &read.valid_indels {
                if pon.contains(indel) {
                    increment(&mut self.unmatched_pon, indel);
                } else {
                    increment(&mut self.unmatched_indels, indel);
                }
            }
        }
        let ids: Vec<&str> = frags.iter().map(|f| f.id()).collect();
        let order = read_id_map_order(&ids);
        let mut slots: Vec<Option<Fragment>> = frags.into_iter().map(Some).collect();
        Ok(order.into_iter().map(|i| slots[i].take().unwrap()).collect())
    }

    /// unmatchedIndels(minCount), without the known stop-loss indel
    pub fn unmatched_indels(&self, min: i32, stop_loss: Option<&Indel>) -> Vec<(Indel, i32)> {
        self.unmatched_indels.iter().filter(|x| x.1 >= min && stop_loss.is_none_or(|s| x.0 != *s)).cloned().collect()
    }
    pub fn unmatched_pon_indels(&self, min: i32) -> Vec<(Indel, i32)> {
        self.unmatched_pon.iter().filter(|x| x.1 >= min).cloned().collect()
    }
}
