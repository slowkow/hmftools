//! RegionReadData, GeneReadData and GeneCollection.

use crate::common::*;
use crate::read::{Read, RKey};
use hmf_common::ensembl::{EnsemblDataCache, GeneData, TranscriptData};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::regions::BaseRegion;
use std::collections::HashMap;

/// RegionReadData: a unique exon (start and end) of one or more transcripts
#[derive(Clone, Debug)]
pub struct RegionReadData {
    pub start: i32,
    pub end: i32,
    pub trans_exon_refs: Vec<TransExonRef>,
    pub ref_bases: Vec<u8>,
    /// set with the reference bases
    pub ref_bases_matched: Option<Vec<i32>>,
    pub unique_ref_bases: Vec<bool>,
    pub trans_read_counts: HashMap<i32, [i32; 2]>,
    pub trans_junction_counts: JavaHashMap<i32, [[i32; 2]; 2]>,
    pub pre_regions: Vec<usize>,
    pub post_regions: Vec<usize>,
}

impl RegionReadData {
    pub fn new(start: i32, end: i32) -> Self {
        RegionReadData {
            start,
            end,
            trans_exon_refs: Vec::new(),
            ref_bases: Vec::new(),
            ref_bases_matched: None,
            unique_ref_bases: Vec::new(),
            trans_read_counts: HashMap::new(),
            trans_junction_counts: JavaHashMap::new(),
            pre_regions: Vec::new(),
            post_regions: Vec::new(),
        }
    }
    pub fn exon_rank(&self, trans_id: i32) -> i32 {
        self.trans_exon_refs.iter().find(|r| r.trans_id == trans_id).map(|r| r.exon_rank).unwrap_or(NO_EXON)
    }
    pub fn has_trans_id(&self, trans_id: i32) -> bool {
        self.trans_exon_refs.iter().any(|r| r.trans_id == trans_id)
    }
    fn add_exon_ref(&mut self, r: TransExonRef) {
        if self.trans_exon_refs.iter().any(|x| x.trans_id == r.trans_id && x.exon_rank == r.exon_rank) {
            return;
        }
        self.trans_exon_refs.push(r);
    }
    pub fn set_ref_bases(&mut self, bases: Vec<u8>) {
        let n = bases.len();
        self.ref_bases = bases;
        self.ref_bases_matched = Some(vec![0; n]);
        self.unique_ref_bases = vec![true; n];
    }
    /// length(): the length of the reference bases
    pub fn length(&self) -> usize {
        self.ref_bases.len()
    }
    #[allow(dead_code)]
    pub fn transcript_read_count(&self, trans_id: i32) -> [i32; 2] {
        self.trans_read_counts.get(&trans_id).copied().unwrap_or([0, 0])
    }
    pub fn add_transcript_read_match(&mut self, trans_id: i32, is_unique: bool) {
        let c = self.trans_read_counts.entry(trans_id).or_insert([0, 0]);
        if is_unique {
            c[UNIQUE_TRANS_COUNT] += 1;
        }
        c[TRANS_COUNT] += 1;
    }
    pub fn transcript_junction_match_count(&self, trans_id: i32, se: usize) -> [i32; 2] {
        self.trans_junction_counts.get(&trans_id).map(|c| c[se]).unwrap_or([0, 0])
    }
    pub fn add_transcript_junction_match(&mut self, trans_id: i32, se: usize, is_unique: bool) {
        if !self.trans_junction_counts.contains_key(&trans_id) {
            self.trans_junction_counts.insert(trans_id, [[0, 0], [0, 0]]);
        }
        let c = self.trans_junction_counts.get_mut(&trans_id).unwrap();
        if is_unique {
            c[se][UNIQUE_TRANS_COUNT] += 1;
        }
        c[se][TRANS_COUNT] += 1;
    }
    pub fn base_coverage(&self, min_count: i32) -> i32 {
        match &self.ref_bases_matched {
            None => 0,
            Some(m) => m.iter().filter(|x| **x >= min_count).count() as i32,
        }
    }
}

/// RegionReadData.findUniqueBases
pub fn find_unique_bases(regions: &mut [RegionReadData]) {
    let n = regions.len();
    for i in 0..n.saturating_sub(1) {
        for j in i + 1..n {
            let (a, b) = regions.split_at_mut(j);
            let r1 = &mut a[i];
            let r2 = &mut b[0];
            if r1.start > r2.end || r1.end < r2.start {
                continue;
            }
            let overlap = (r1.end.min(r2.end) - r1.start.max(r2.start)) + 1;
            let this_off = (r2.start - r1.start).max(0);
            let other_off = (r1.start - r2.start).max(0);
            for k in 0..overlap {
                r1.unique_ref_bases[(k + this_off) as usize] = false;
                r2.unique_ref_bases[(k + other_off) as usize] = false;
            }
        }
    }
}

pub fn find_exon_region(regions: &[RegionReadData], indices: Option<&[usize]>, start: i32, end: i32) -> Option<usize> {
    match indices {
        Some(idx) => idx.iter().copied().find(|i| regions[*i].start == start && regions[*i].end == end),
        None => regions.iter().position(|r| r.start == start && r.end == end),
    }
}

/// GeneReadData: a gene, its transcripts and its exon regions
#[derive(Clone, Debug)]
pub struct GeneReadData {
    pub gene: &'static GeneData,
    /// indices into the gene collection's exon regions
    pub exon_regions: Vec<usize>,
    pub transcripts: Vec<&'static TranscriptData>,
    pub has_unspliced_regions: bool,
}

impl GeneReadData {
    pub fn name(&self) -> &str {
        &self.gene.gene_name
    }
}

/// createGeneReadData
pub fn create_gene_read_data(genes: &[&'static GeneData], cache: &'static EnsemblDataCache) -> Vec<GeneReadData> {
    let mut out = Vec::new();
    for g in genes {
        let trans = cache.transcripts(&g.gene_id);
        match trans {
            Some(t) if !t.is_empty() => out.push(GeneReadData { gene: g, exon_regions: Vec::new(), transcripts: t.iter().collect(), has_unspliced_regions: false }),
            _ => {
                hmf_common::warn!("no transcripts found for gene({}:{} chr={})", g.gene_id, g.gene_name, g.chromosome);
            }
        }
    }
    out
}

/// GeneReadData.generateCommonExonicRegions
pub fn generate_common_exonic_regions(regions: &[RegionReadData], indices: &[usize]) -> Vec<[i32; 2]> {
    if indices.is_empty() {
        return Vec::new();
    }
    let first = &regions[indices[0]];
    let mut common = vec![BaseRegion::new(first.start, first.end)];
    for i in &indices[1..] {
        let r = &regions[*i];
        common = derive_common_regions(&common, &[BaseRegion::new(r.start, r.end)]);
    }
    common.iter().map(|x| [x.start, x.end]).collect()
}

/// GeneCollection: overlapping genes processed together
pub struct GeneCollection {
    pub id: i32,
    pub genes: Vec<GeneReadData>,
    pub gene_ids: Vec<&'static str>,
    pub chromosome: &'static str,
    pub region_bounds: [i32; 2],
    pub non_genic_positions: [i32; 2],
    pub end_of_chromosome: bool,
    trans_ids_gene: HashMap<i32, usize>,
    pub exon_regions: Vec<RegionReadData>,
    pub common_exonic_regions: Vec<[i32; 2]>,
    pub transcripts: Vec<&'static TranscriptData>,
    pub transcript_read_counts: HashMap<i32, [[i32; 2]; 5]>,
    pub fragment_counts: FragmentTypeCounts,
}

impl GeneCollection {
    pub fn new(id: i32, genes: Vec<GeneReadData>) -> GeneCollection {
        let chromosome: &'static str = &genes[0].gene.chromosome;
        let gene_ids = genes.iter().map(|g| g.gene.gene_id.as_str()).collect();
        let mut gc = GeneCollection {
            id,
            genes,
            gene_ids,
            chromosome,
            region_bounds: [0, 0],
            non_genic_positions: [0, 0],
            end_of_chromosome: false,
            trans_ids_gene: HashMap::new(),
            exon_regions: Vec::new(),
            common_exonic_regions: Vec::new(),
            transcripts: Vec::new(),
            transcript_read_counts: HashMap::new(),
            fragment_counts: FragmentTypeCounts::default(),
        };
        gc.build_cache();
        gc
    }

    fn build_cache(&mut self) {
        for gi in 0..self.genes.len() {
            let gene_trans = self.genes[gi].transcripts.clone();
            for t in &gene_trans {
                self.trans_ids_gene.insert(t.trans_id, gi);
                self.transcripts.push(t);
                self.region_bounds[0] = if self.region_bounds[0] == 0 { t.trans_start } else { self.region_bounds[0].min(t.trans_start) };
                self.region_bounds[1] = self.region_bounds[1].max(t.trans_end);
            }
            // generateExonicRegions
            for t in &gene_trans {
                let mut prev: Option<usize> = None;
                for exon in &t.exons {
                    let idx = match find_exon_region(&self.exon_regions, None, exon.start, exon.end) {
                        Some(i) => i,
                        None => {
                            self.exon_regions.push(RegionReadData::new(exon.start, exon.end));
                            self.exon_regions.len() - 1
                        }
                    };
                    self.exon_regions[idx].add_exon_ref(TransExonRef::new(t, exon.rank));
                    if let Some(p) = prev {
                        if !self.exon_regions[p].post_regions.contains(&idx) {
                            self.exon_regions[p].post_regions.push(idx);
                        }
                        if !self.exon_regions[idx].pre_regions.contains(&p) {
                            self.exon_regions[idx].pre_regions.push(p);
                        }
                    }
                    prev = Some(idx);
                }
            }
            for t in &gene_trans {
                for exon in &t.exons {
                    let Some(idx) = find_exon_region(&self.exon_regions, None, exon.start, exon.end) else {
                        hmf_common::error!("genes({}) failed to create exonic regions", self.gene_names(10));
                        return;
                    };
                    let g = &mut self.genes[gi];
                    let exists = g.exon_regions.iter().any(|i| self.exon_regions[*i].start == self.exon_regions[idx].start && self.exon_regions[*i].end == self.exon_regions[idx].end);
                    if !exists {
                        g.exon_regions.push(idx);
                    }
                }
            }
            self.set_has_unspliced_regions(gi);
        }
        let all: Vec<usize> = (0..self.exon_regions.len()).collect();
        self.common_exonic_regions = generate_common_exonic_regions(&self.exon_regions, &all);
    }

    fn set_has_unspliced_regions(&mut self, gi: usize) {
        let regions = &self.exon_regions;
        let g = &mut self.genes[gi];
        g.has_unspliced_regions = false;
        if g.exon_regions.is_empty() {
            return;
        }
        let min = g.exon_regions.iter().map(|i| regions[*i].start).min().unwrap_or(0);
        let max = g.exon_regions.iter().map(|i| regions[*i].end).max().unwrap_or(0);
        for i in &g.exon_regions {
            let r = &regions[*i];
            if r.start > min && !g.exon_regions.iter().any(|x| position_within(r.start - 1, regions[*x].start, regions[*x].end)) {
                g.has_unspliced_regions = true;
                return;
            }
            if r.end < max && !g.exon_regions.iter().any(|x| position_within(r.end + 1, regions[*x].start, regions[*x].end)) {
                g.has_unspliced_regions = true;
                return;
            }
        }
    }

    pub fn chr_id(&self) -> String {
        format!("{}_{}", self.chromosome, self.id)
    }
    pub fn gene_names(&self, max_count: usize) -> String {
        if self.genes.len() == 1 {
            return self.genes[0].name().to_string();
        }
        self.genes.iter().take(max_count).map(|g| g.name()).collect::<Vec<_>>().join(";")
    }
    pub fn strand(&self, trans_id: i32) -> i8 {
        self.trans_ids_gene.get(&trans_id).map(|g| self.genes[*g].gene.strand).unwrap_or(0)
    }
    #[allow(dead_code)]
    pub fn rkey(&self, idx: usize) -> RKey {
        RKey { idx: idx as u32, start: self.exon_regions[idx].start }
    }

    /// setReadGeneCollections
    pub fn set_read_gene_collections(&self, read: &mut Read, non_genic: [i32; 2]) {
        if positions_within(read.pos_start, read.pos_end, self.region_bounds[0], self.region_bounds[1]) {
            read.set_gene_collection(SE_START, self.id, true);
            read.set_gene_collection(SE_END, self.id, true);
        } else {
            for se in [SE_START, SE_END] {
                let b = read.coords_boundary(se);
                if position_within(b, self.region_bounds[0], self.region_bounds[1]) {
                    read.set_gene_collection(se, self.id, true);
                } else if position_within(b, non_genic[0], non_genic[1]) {
                    read.set_gene_collection(se, self.id, false);
                }
            }
        }
    }

    /// findGenesCoveringRange: gene indices
    pub fn find_genes_covering_range(&self, start: i32, end: i32, check_unspliced: bool) -> Vec<usize> {
        if check_unspliced {
            (0..self.genes.len())
                .filter(|i| self.genes[*i].has_unspliced_regions && positions_within(start, end, self.genes[*i].gene.gene_start, self.genes[*i].gene.gene_end))
                .collect()
        } else {
            (0..self.genes.len()).filter(|i| positions_overlap(self.genes[*i].gene.gene_start, self.genes[*i].gene.gene_end, start, end)).collect()
        }
    }

    pub fn transcript_read_count(&self, trans_id: i32) -> [[i32; 2]; 5] {
        self.transcript_read_counts.get(&trans_id).copied().unwrap_or([[0; 2]; 5])
    }
    pub fn add_transcript_read_match(&mut self, trans_id: i32, is_unique: bool, t: FragmentMatchType) {
        let c = self.transcript_read_counts.entry(trans_id).or_insert([[0; 2]; 5]);
        if is_unique {
            c[t as usize][UNIQUE_TRANS_COUNT] += 1;
        }
        c[t as usize][TRANS_COUNT] += 1;
    }
    pub fn add_transcript_read_match_type(&mut self, trans_id: i32, t: FragmentMatchType) {
        let c = self.transcript_read_counts.entry(trans_id).or_insert([[0; 2]; 5]);
        c[t as usize][TRANS_COUNT] += 1;
    }
    pub fn add_count(&mut self, t: FragmentType, count: i64) {
        self.fragment_counts.add_count(t, count);
    }
    pub fn transcript(&self, trans_id: i32) -> Option<&'static TranscriptData> {
        self.transcripts.iter().copied().find(|t| t.trans_id == trans_id)
    }
    /// GeneReadData.calcExonicRegionLength: the sum of (end - start) of its merged exons
    pub fn exonic_region_length(&self, gi: usize) -> i64 {
        generate_common_exonic_regions(&self.exon_regions, &self.genes[gi].exon_regions).iter().map(|x| (x[1] - x[0]) as i64).sum()
    }
}

/// ChromosomeTaskExecutor.findNextOverlappingGenes
pub fn find_next_overlapping_genes(genes: &[&'static GeneData], mut index: usize, out: &mut Vec<&'static GeneData>) -> usize {
    out.clear();
    while index < genes.len() {
        let g = genes[index];
        if out.is_empty() || out.iter().any(|x| positions_overlap(g.gene_start, g.gene_end, x.gene_start, x.gene_end)) {
            out.push(g);
            index += 1;
        } else {
            break;
        }
    }
    index
}
