//! GeneDataCache: the Ensembl genes and transcripts and the driver gene panel, with the per-chromosome
//! gene search the chromosome tasks use (GeneCacheIndexing).

use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};

use hmf_common::driver_panel::{self, DriverGene};
use hmf_common::ensembl::{EnsemblDataCache, GeneData, TranscriptData};
use hmf_common::regions::{positions_overlap, RefGenomeVersion};

use crate::compat;
use crate::constants::GENE_UPSTREAM_DISTANCE;
use crate::impact::utils::within_trans_range;
use crate::impact::POS_STRAND;

pub struct GeneDataCache {
    ensembl: EnsemblDataCache,
    driver_genes: Vec<DriverGene>,
    driver_gene_names: HashSet<String>,
    other_reportable_transcripts: HashMap<String, Vec<String>>,
}

/// GeneCacheIndexing: a chromosome's genes and those near the current position
pub struct GeneCacheIndexing {
    pub chromosome_genes: Option<&'static Vec<GeneData>>,
    pub current_genes: Vec<&'static GeneData>,
    pub current_pos_strand_gene_index: usize,
    pub current_neg_strand_gene_index: usize,
}

impl GeneDataCache {
    /// new GeneDataCache(...) and loadCache(canonicalOnly, false)
    pub fn load(ensembl_dir: &str, version: RefGenomeVersion, driver_gene_file: Option<&str>, canonical_only: bool) -> Result<GeneDataCache> {
        let mut cache = GeneDataCache {
            ensembl: EnsemblDataCache::new(ensembl_dir, version),
            driver_genes: Vec::new(),
            driver_gene_names: HashSet::new(),
            other_reportable_transcripts: HashMap::new(),
        };
        if let Some(f) = driver_gene_file.filter(|f| !f.is_empty()) {
            cache.driver_genes = driver_panel::read_driver_genes(f).with_context(|| format!("failed to load driver gene panel file({f})"))?;
            for g in &cache.driver_genes {
                cache.driver_gene_names.insert(g.gene.clone());
                if !g.additional_reported_transcripts.is_empty() {
                    cache.other_reportable_transcripts.insert(g.gene.clone(), g.additional_reported_transcripts.clone());
                }
            }
        }
        cache.ensembl.set_required_data(true, false, false, canonical_only);
        cache.ensembl.set_require_non_ensembl_transcripts();
        if let Err(e) = cache.ensembl.load(false) {
            bail!("invalid Ensembl data cache: {e:#}");
        }
        cache.ensembl.create_gene_id_data_map();
        Ok(cache)
    }

    #[allow(dead_code)]
    pub fn ensembl(&self) -> &EnsemblDataCache {
        &self.ensembl
    }

    /// new GeneDataCache("", version, null) without loading, for tests to fill
    #[cfg(test)]
    pub fn empty(version: RefGenomeVersion) -> GeneDataCache {
        GeneDataCache {
            ensembl: EnsemblDataCache::new("", version),
            driver_genes: Vec::new(),
            driver_gene_names: HashSet::new(),
            other_reportable_transcripts: HashMap::new(),
        }
    }

    #[cfg(test)]
    pub fn ensembl_mut(&mut self) -> &mut EnsemblDataCache {
        &mut self.ensembl
    }

    #[cfg(test)]
    #[allow(dead_code)]
    pub fn add_driver_gene(&mut self, gene: DriverGene) {
        self.driver_gene_names.insert(gene.gene.clone());
        if !gene.additional_reported_transcripts.is_empty() {
            self.other_reportable_transcripts.insert(gene.gene.clone(), gene.additional_reported_transcripts.clone());
        }
        self.driver_genes.push(gene);
    }
    pub fn is_driver_panel_gene(&self, gene_name: &str) -> bool {
        self.driver_gene_names.contains(gene_name)
    }
    pub fn driver_panel(&self) -> &[DriverGene] {
        &self.driver_genes
    }
    pub fn other_reportable_transcripts(&self, gene_name: &str) -> Option<&Vec<String>> {
        self.other_reportable_transcripts.get(gene_name)
    }

    /// findTranscripts: the gene's transcripts whose range (with the upstream distance) overlaps
    pub fn find_transcripts(&'static self, gene_id: &str, start: i32, end: i32) -> Vec<&'static TranscriptData> {
        match self.ensembl.transcripts(gene_id) {
            None => Vec::new(),
            Some(list) => list.iter().filter(|t| within_trans_range(t, start, end)).collect(),
        }
    }

    pub fn create_indexing(&'static self, chromosome: &str) -> GeneCacheIndexing {
        GeneCacheIndexing {
            chromosome_genes: self.ensembl.chr_gene_data_map().get(&chromosome.to_string()),
            current_genes: Vec::new(),
            current_pos_strand_gene_index: 0,
            current_neg_strand_gene_index: 0,
        }
    }

    /// findGenes: the genes whose range (with the upstream distance) the variant overlaps. With an
    /// indexing, a forward search from the previous variant's position that keeps the genes still near.
    pub fn find_genes(&'static self, chromosome: &str, start: i32, end: i32, indexing: Option<&mut GeneCacheIndexing>) -> Vec<&'static GeneData> {
        let Some(ix) = indexing else {
            let mut genes = Vec::new();
            let Some(list) = self.ensembl.chr_gene_data_map().get(&chromosome.to_string()) else { return genes };
            for g in list {
                if is_within_gene_range(g, start, end) {
                    genes.push(g);
                } else if g.gene_start > end + GENE_UPSTREAM_DISTANCE {
                    break;
                }
            }
            return genes;
        };

        let Some(chromosome_genes) = ix.chromosome_genes else { return Vec::new() };

        // PA1: Java drops a current gene once its end (plus the upstream distance on the negative strand) is
        // before the variant's end position. A variant reaching past the end of a gene it overlaps, or any
        // variant after a long deletion which did, then misses the gene, since the search goes on from
        // genes not yet seen. Corrected: a gene is dropped once it ends before the variant's start, which
        // no later variant precedes.
        let purge_pos = if compat::java("PA1") { end } else { start };
        ix.current_genes.retain(|g| {
            let is_past_gene = (g.forward_strand() && g.gene_end < purge_pos) || (!g.forward_strand() && g.gene_end + GENE_UPSTREAM_DISTANCE < purge_pos);
            !is_past_gene
        });

        // search forward and add any further genes overlapping the position; one index per strand since
        // the range checks can stop each at another gene
        for i in 0..=1 {
            let current_strand: i8 = if i == 0 { POS_STRAND } else { -1 };
            let mut gene_index = if i == 0 { ix.current_pos_strand_gene_index } else { ix.current_neg_strand_gene_index };
            while gene_index < chromosome_genes.len() {
                let g = &chromosome_genes[gene_index];
                if g.strand != current_strand {
                    gene_index += 1;
                    continue;
                }
                if is_within_gene_range(g, start, end) {
                    if ix.current_genes.iter().any(|x| std::ptr::eq(*x, g)) {
                        crate::error!("adding current gene({}:{}) index({}) twice", g.gene_id, g.gene_name, gene_index);
                        break;
                    }
                    ix.current_genes.push(g);
                } else if start > g.gene_end {
                    gene_index += 1;
                    continue;
                } else if g.gene_start > start {
                    break;
                }
                gene_index += 1;
            }
            if i == 0 {
                ix.current_pos_strand_gene_index = gene_index;
            } else {
                ix.current_neg_strand_gene_index = gene_index;
            }
        }
        ix.current_genes.clone()
    }
}

fn is_within_gene_range(g: &GeneData, start: i32, end: i32) -> bool {
    let mut range_start = g.gene_start;
    let mut range_end = g.gene_end;
    if g.strand == POS_STRAND {
        range_start -= GENE_UPSTREAM_DISTANCE;
    } else {
        range_end += GENE_UPSTREAM_DISTANCE;
    }
    positions_overlap(start, end, range_start, range_end)
}
