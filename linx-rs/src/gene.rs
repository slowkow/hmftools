//! Package `gene`: the genes and transcripts a breakend falls in or before (BreakendGeneData,
//! BreakendTransData) and their annotation (BreakendGenePrep). Genes and transcripts are references into
//! the Ensembl data cache, which lives for the whole run.

use std::cell::{Cell, RefCell};

use hmf_common::ensembl::{EnsemblDataCache, ExonData, GeneData, TranscriptData};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::java_id_hash::{IdHashMap, IdKey};
use hmf_common::sv::{StructuralVariantData, StructuralVariantType};
use hmf_common::transcript_utils::{calc_coding_bases, tick_phase_forward, PHASE_NONE};

use crate::types::*;

pub type Gene = &'static GeneData;
pub type Trans = &'static TranscriptData;

/// A GeneData key compared by identity (GeneData has no equals or hashCode)
#[derive(Clone, Copy)]
pub struct GeneKey(pub Gene);
impl IdKey for GeneKey {
    fn id_addr(&self) -> usize {
        self.0 as *const GeneData as usize
    }
}
impl PartialEq for GeneKey {
    fn eq(&self, o: &Self) -> bool {
        std::ptr::eq(self.0, o.0)
    }
}

/// TranscriptCodingType
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TranscriptCodingType {
    UNKNOWN,
    CODING,
    UTR_5P,
    UTR_3P,
    NON_CODING,
    ENHANCER,
}
impl std::fmt::Display for TranscriptCodingType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// TranscriptRegionType
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TranscriptRegionType {
    UNKNOWN,
    UPSTREAM,
    EXONIC,
    INTRONIC,
    IG,
    DOWNSTREAM,
}
impl std::fmt::Display for TranscriptRegionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// The ensembl data cache plus what LINX adds to it at run time: downstream gene annotations and
/// alternative gene data (from the known fusions), which Java keeps inside its EnsemblDataCache
pub struct GeneCache {
    pub ensembl: &'static EnsemblDataCache,
    /// mDownstreamGeneAnnotations: a HashMap<GeneData,Integer> iterated in identity-hash order
    pub downstream_gene_annotations: RefCell<IdHashMap<GeneKey, i32>>,
    pub alternative_gene_data: RefCell<Vec<Gene>>,
}

impl GeneCache {
    pub fn new(ensembl: &'static EnsemblDataCache) -> GeneCache {
        GeneCache { ensembl, downstream_gene_annotations: RefCell::new(IdHashMap::new()), alternative_gene_data: RefCell::new(Vec::new()) }
    }
    pub fn add_downstream_gene_annotations(&self, gene: Gene, distance: i32) {
        self.downstream_gene_annotations.borrow_mut().insert(GeneKey(gene), distance);
    }
    pub fn has_downstream_gene_annotation(&self, gene: Gene) -> bool {
        self.downstream_gene_annotations.borrow().contains_key(&GeneKey(gene))
    }
    pub fn gene_data_by_name(&self, name: &str) -> Option<Gene> {
        self.ensembl.gene_data_by_name(name)
    }
    pub fn transcripts(&self, gene_id: &str) -> Option<&'static Vec<TranscriptData>> {
        self.ensembl.transcripts(gene_id)
    }
    pub fn transcript_data(&self, gene_id: &str, trans_name: &str) -> Option<Trans> {
        self.ensembl.transcript_data(gene_id, trans_name)
    }

    /// EnsemblDataCache.findGeneRegions, with the downstream gene annotations
    pub fn find_gene_regions(&self, chromosome: &str, position: i32, upstream_distance: i32) -> Vec<Gene> {
        let mut matched: Vec<Gene> = self.ensembl.find_gene_regions(chromosome, position, upstream_distance);
        for (key, dist) in self.downstream_gene_annotations.borrow().iter() {
            let g = key.0;
            if matched.iter().any(|m| std::ptr::eq(*m, g)) || g.chromosome != chromosome {
                continue;
            }
            if (g.strand == 1 && position >= g.gene_end && position <= g.gene_end + dist)
                || (g.strand == -1 && position <= g.gene_start && position >= g.gene_start - dist)
            {
                matched.push(g);
            }
        }
        matched
    }
}

/// BreakendGeneData: a gene a breakend touches
pub struct BreakendGeneData {
    pub gene_data: Gene,
    sv_id: i32,
    sv_vcf_id: RefCell<String>,
    is_start: bool,
    upstream: Cell<bool>,
    transcripts: RefCell<Vec<BreakendTrans>>,
    sv_type: Cell<Option<StructuralVariantType>>,
    chromosome: RefCell<String>,
    orientation: Cell<i8>,
    position: Cell<i32>,
    junction_copy_number: Cell<f64>,
    insert_sequence: RefCell<String>,
}

handle!(BreakendGene, BreakendGeneData);

impl BreakendGeneData {
    pub fn gene_name(&self) -> &'static str {
        &self.gene_data.gene_name
    }
    pub fn gene_id(&self) -> &'static str {
        &self.gene_data.gene_id
    }
    pub fn strand(&self) -> i8 {
        self.gene_data.strand
    }
    pub fn set_positional_data(&self, chromosome: &str, position: i32, orientation: i8) {
        *self.chromosome.borrow_mut() = chromosome.to_string();
        self.position.set(position);
        self.orientation.set(orientation);
        self.upstream.set(self.strand() as i32 * orientation as i32 > 0);
    }
    pub fn set_sv_data(&self, var: &StructuralVariantData, jcn: f64) {
        *self.sv_vcf_id.borrow_mut() = if self.is_start { var.vcf_id_start.clone() } else { var.vcf_id_end.clone() };
        self.orientation.set(if self.is_start { var.start_orientation } else { var.end_orientation });
        self.junction_copy_number.set(jcn);
        self.position.set(if self.is_start { var.start_position } else { var.end_position });
        *self.chromosome.borrow_mut() = if self.is_start { var.start_chromosome.clone() } else { var.end_chromosome.clone() };
        self.sv_type.set(Some(var.sv_type()));
        *self.insert_sequence.borrow_mut() = var.insert_sequence.clone();
        self.upstream.set(self.strand() as i32 * self.orientation.get() as i32 > 0);
    }
    pub fn set_type(&self, t: StructuralVariantType) {
        self.sv_type.set(Some(t))
    }
    pub fn var_id(&self) -> i32 {
        self.sv_id
    }
    pub fn vcf_id(&self) -> String {
        self.sv_vcf_id.borrow().clone()
    }
    pub fn orientation(&self) -> i8 {
        self.orientation.get()
    }
    pub fn position(&self) -> i32 {
        self.position.get()
    }
    pub fn sv_type(&self) -> Option<StructuralVariantType> {
        self.sv_type.get()
    }
    pub fn chromosome(&self) -> String {
        self.chromosome.borrow().clone()
    }
    pub fn jcn(&self) -> f64 {
        self.junction_copy_number.get()
    }
    pub fn insert_sequence(&self) -> String {
        self.insert_sequence.borrow().clone()
    }
    pub fn is_upstream(&self) -> bool {
        self.upstream.get()
    }
    pub fn is_start(&self) -> bool {
        self.is_start
    }
    pub fn is_end(&self) -> bool {
        !self.is_start
    }
    pub fn transcripts(&self) -> Vec<BreakendTrans> {
        self.transcripts.borrow().clone()
    }
    pub fn transcript_count(&self) -> usize {
        self.transcripts.borrow().len()
    }
    pub fn add_transcript(&self, t: BreakendTrans) {
        self.transcripts.borrow_mut().push(t);
    }
    pub fn retain_transcripts(&self, f: impl FnMut(&BreakendTrans) -> bool) {
        self.transcripts.borrow_mut().retain(f);
    }
    pub fn canonical(&self) -> Option<BreakendTrans> {
        self.transcripts.borrow().iter().find(|t| t.is_canonical()).cloned()
    }

    /// breakendWithinGene(preGeneDistance)
    pub fn breakend_within_gene(&self, pre_gene_distance: i32) -> bool {
        let g = self.gene_data;
        let pos = self.position.get();
        if pos >= g.gene_start && pos <= g.gene_end {
            return true;
        }
        if pre_gene_distance <= 0 {
            return false;
        }
        for trans in self.transcripts.borrow().iter() {
            if !trans.is_upstream() && trans.has_negative_prev_splice_acceptor_distance() {
                continue;
            }
            let distance = if g.strand == 1 { trans.trans_data.trans_start - pos } else { pos - trans.trans_data.trans_end };
            if distance > 0 && distance <= pre_gene_distance {
                return true;
            }
        }
        false
    }

    pub fn coords_str(&self) -> String {
        coords_str(&self.chromosome(), self.position(), self.orientation())
    }
}

impl BreakendGene {
    pub fn create(sv_id: i32, is_start: bool, gene: Gene) -> BreakendGene {
        BreakendGene::new(BreakendGeneData {
            gene_data: gene,
            sv_id,
            sv_vcf_id: RefCell::new(String::new()),
            is_start,
            upstream: Cell::new(false),
            transcripts: RefCell::new(Vec::new()),
            sv_type: Cell::new(None),
            chromosome: RefCell::new(String::new()),
            orientation: Cell::new(0),
            position: Cell::new(-1),
            junction_copy_number: Cell::new(0.0),
            insert_sequence: RefCell::new(String::new()),
        })
    }
}

impl std::fmt::Display for BreakendGeneData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "gene({}:{} strand={}) breakend(sv={} coords={}) trans({})",
            self.gene_id(),
            self.gene_name(),
            self.strand(),
            self.var_id(),
            self.coords_str(),
            self.transcript_count()
        )
    }
}

pub const POST_CODING_PHASE: i32 = -2;
const STOP_CODON_LENGTH: i32 = 3;
pub const NO_NEXT_SPLICE_ACCEPTOR: i32 = -1;

/// BreakendTransData: a transcript a breakend touches
pub struct BreakendTransData {
    pub trans_data: Trans,
    pub exon_upstream: i32,
    pub exon_downstream: i32,
    gene: BreakendGene,
    pub coding_bases: i32,
    pub total_coding_bases: i32,
    pub phase: i32,
    pub exonic_base_phase: i32,
    coding_type: Cell<TranscriptCodingType>,
    region_type: Cell<TranscriptRegionType>,
    prev_splice_acceptor_distance: Cell<Option<i32>>,
    next_splice_acceptor_distance: Cell<Option<i32>>,
    alternative_phasing: RefCell<JavaHashMap<i32, i32>>,
    is_disruptive: Cell<bool>,
    reportable_disruption: Cell<bool>,
    undisrupted_copy_number: Cell<f64>,
    undisrupted_copy_number_set: Cell<bool>,
    protein_features_kept: RefCell<String>,
    protein_features_lost: RefCell<String>,
}

handle!(BreakendTrans, BreakendTransData);

impl BreakendTrans {
    #[allow(clippy::too_many_arguments)]
    pub fn create(
        gene: &BreakendGene, trans_data: Trans, exon_upstream: i32, exon_downstream: i32, phase: i32, mut exonic_base_phase: i32,
        mut coding_bases: i32, total_coding_bases: i32,
    ) -> BreakendTrans {
        let total = if total_coding_bases > STOP_CODON_LENGTH { total_coding_bases - STOP_CODON_LENGTH } else { 0 };
        let ins = gene.insert_sequence();
        if gene.is_upstream() && !ins.is_empty() && exon_upstream == exon_downstream && coding_bases > 0 && coding_bases < total_coding_bases {
            let len = ins.chars().count() as i32;
            coding_bases += len;
            exonic_base_phase = tick_phase_forward(phase, len);
        }
        let cb = if coding_bases > total { total } else { coding_bases };
        let mut t = BreakendTransData {
            trans_data,
            exon_upstream,
            exon_downstream,
            gene: gene.clone(),
            coding_bases: cb,
            total_coding_bases: total,
            phase: PHASE_NONE,
            exonic_base_phase,
            coding_type: Cell::new(TranscriptCodingType::UNKNOWN),
            region_type: Cell::new(TranscriptRegionType::UNKNOWN),
            prev_splice_acceptor_distance: Cell::new(None),
            next_splice_acceptor_distance: Cell::new(None),
            alternative_phasing: RefCell::new(JavaHashMap::new()),
            is_disruptive: Cell::new(false),
            reportable_disruption: Cell::new(false),
            undisrupted_copy_number: Cell::new(0.0),
            undisrupted_copy_number_set: Cell::new(false),
            protein_features_kept: RefCell::new(String::new()),
            protein_features_lost: RefCell::new(String::new()),
        };
        t.init_coding_type();
        t.init_region_type();
        t.phase = if t.coding_type.get() == TranscriptCodingType::UTR_3P { POST_CODING_PHASE } else { phase };
        let rt = t.region_type.get();
        t.is_disruptive.set(rt == TranscriptRegionType::EXONIC || rt == TranscriptRegionType::INTRONIC);
        BreakendTrans::new(t)
    }
}

impl BreakendTransData {
    pub fn gene(&self) -> &BreakendGene {
        &self.gene
    }
    pub fn sv_position(&self) -> i32 {
        self.gene.position()
    }
    pub fn gene_name(&self) -> &'static str {
        self.gene.gene_name()
    }
    pub fn is_upstream(&self) -> bool {
        self.gene.is_upstream()
    }
    pub fn trans_id(&self) -> i32 {
        self.trans_data.trans_id
    }
    pub fn trans_name(&self) -> &'static str {
        &self.trans_data.trans_name
    }
    pub fn trans_start(&self) -> i32 {
        self.trans_data.trans_start
    }
    pub fn trans_end(&self) -> i32 {
        self.trans_data.trans_end
    }
    pub fn bio_type(&self) -> &'static str {
        &self.trans_data.bio_type
    }
    pub fn coding_start(&self) -> i32 {
        self.trans_data.coding_start.unwrap_or(0)
    }
    pub fn coding_end(&self) -> i32 {
        self.trans_data.coding_end.unwrap_or(0)
    }
    pub fn exon_count(&self) -> i32 {
        self.trans_data.exons.len() as i32
    }
    pub fn is_exonic(&self) -> bool {
        self.exon_upstream > 0 && self.exon_upstream == self.exon_downstream
    }
    pub fn is_promoter(&self) -> bool {
        self.exon_upstream == 0 && (self.exon_downstream == 1 || self.exon_downstream == 2)
    }
    pub fn is_intronic(&self) -> bool {
        self.exon_upstream > 0 && (self.exon_downstream - self.exon_upstream) == 1
    }
    pub fn is_post_transcript(&self) -> bool {
        if self.gene.strand() == 1 { self.sv_position() > self.trans_data.trans_end } else { self.sv_position() < self.trans_data.trans_start }
    }
    pub fn distance_upstream(&self) -> i32 {
        if !self.is_promoter() {
            return 0;
        }
        if self.gene.strand() == 1 { self.trans_data.trans_start - self.sv_position() } else { self.sv_position() - self.trans_data.trans_end }
    }
    pub fn distance_downstream(&self) -> i32 {
        if self.coding_type.get() != TranscriptCodingType::UTR_3P {
            return 0;
        }
        if self.gene.strand() == 1 { self.sv_position() - self.trans_data.trans_end } else { self.trans_data.trans_start - self.sv_position() }
    }
    pub fn coding_type(&self) -> TranscriptCodingType {
        self.coding_type.get()
    }
    pub fn region_type(&self) -> TranscriptRegionType {
        self.region_type.get()
    }
    pub fn is_coding(&self) -> bool {
        self.coding_type.get() == TranscriptCodingType::CODING
    }
    pub fn pre_coding(&self) -> bool {
        self.coding_type.get() == TranscriptCodingType::UTR_5P
    }
    pub fn post_coding(&self) -> bool {
        self.coding_type.get() == TranscriptCodingType::UTR_3P
    }
    pub fn non_coding(&self) -> bool {
        self.coding_type.get() == TranscriptCodingType::NON_CODING
    }

    fn init_region_type(&mut self) {
        let t = if self.is_promoter() {
            TranscriptRegionType::UPSTREAM
        } else if self.is_post_transcript() {
            TranscriptRegionType::DOWNSTREAM
        } else if self.is_intronic() {
            TranscriptRegionType::INTRONIC
        } else if self.is_exonic() {
            TranscriptRegionType::EXONIC
        } else {
            TranscriptRegionType::UNKNOWN
        };
        self.region_type.set(t);
    }
    pub fn set_region_type(&self, t: TranscriptRegionType) {
        self.region_type.set(t)
    }

    pub fn next_splice_exon_rank(&self) -> i32 {
        if self.is_exonic() {
            if self.is_upstream() { self.exon_upstream - 1 } else { self.exon_downstream + 1 }
        } else if self.is_upstream() {
            self.exon_upstream
        } else {
            self.exon_downstream
        }
    }

    fn init_coding_type(&mut self) {
        let (Some(cs), Some(ce)) = (self.trans_data.coding_start, self.trans_data.coding_end) else {
            self.coding_type.set(TranscriptCodingType::NON_CODING);
            return;
        };
        let position = self.gene.position();
        let t = if position >= cs && position <= ce {
            TranscriptCodingType::CODING
        } else if self.trans_data.strand == 1 {
            if position > ce { TranscriptCodingType::UTR_3P } else { TranscriptCodingType::UTR_5P }
        } else if position < cs {
            TranscriptCodingType::UTR_3P
        } else {
            TranscriptCodingType::UTR_5P
        };
        self.coding_type.set(t);
    }
    pub fn set_coding_type(&self, t: TranscriptCodingType) {
        self.coding_type.set(t)
    }

    pub fn is_canonical(&self) -> bool {
        self.trans_data.is_canonical
    }
    pub fn alternative_phasing(&self) -> std::cell::Ref<'_, JavaHashMap<i32, i32>> {
        self.alternative_phasing.borrow()
    }
    pub fn alternative_phasing_mut(&self) -> std::cell::RefMut<'_, JavaHashMap<i32, i32>> {
        self.alternative_phasing.borrow_mut()
    }
    pub fn set_splice_acceptor_distance(&self, is_previous: bool, distance: Option<i32>) {
        if is_previous {
            self.prev_splice_acceptor_distance.set(distance);
        } else {
            self.next_splice_acceptor_distance.set(distance);
        }
    }
    pub fn has_prev_splice_acceptor_distance(&self) -> bool {
        self.prev_splice_acceptor_distance.get().is_some()
    }
    pub fn has_negative_prev_splice_acceptor_distance(&self) -> bool {
        self.prev_splice_acceptor_distance.get().is_some_and(|d| d < 0)
    }
    pub fn prev_splice_acceptor_distance(&self) -> i32 {
        self.prev_splice_acceptor_distance.get().unwrap_or(NO_NEXT_SPLICE_ACCEPTOR)
    }
    pub fn next_splice_acceptor_distance(&self) -> i32 {
        self.next_splice_acceptor_distance.get().unwrap_or(NO_NEXT_SPLICE_ACCEPTOR)
    }
    pub fn is_disruptive(&self) -> bool {
        self.is_disruptive.get()
    }
    pub fn set_is_disruptive(&self, t: bool) {
        self.is_disruptive.set(t)
    }
    pub fn reportable_disruption(&self) -> bool {
        self.reportable_disruption.get()
    }
    pub fn set_reportable_disruption(&self, t: bool) {
        self.reportable_disruption.set(t)
    }
    pub fn undisrupted_copy_number(&self) -> f64 {
        self.undisrupted_copy_number.get()
    }
    pub fn undisrupted_copy_number_set(&self) -> bool {
        self.undisrupted_copy_number_set.get()
    }
    pub fn set_undisrupted_copy_number(&self, cn: f64) {
        self.undisrupted_copy_number.set(cn);
        self.undisrupted_copy_number_set.set(true);
    }
    pub fn add_protein_feature(&self, feature: &str, is_preserved: bool) {
        let mut s = if is_preserved { self.protein_features_kept.borrow_mut() } else { self.protein_features_lost.borrow_mut() };
        if !s.is_empty() {
            s.push(';');
        }
        s.push_str(feature);
    }
    pub fn protein_features_kept(&self) -> String {
        self.protein_features_kept.borrow().clone()
    }
    pub fn protein_features_lost(&self) -> String {
        self.protein_features_lost.borrow().clone()
    }
}

impl std::fmt::Display for BreakendTransData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", self.gene_name(), self.trans_name())
    }
}

/// BreakendGenePrep.setSvGeneData
pub fn set_sv_gene_data(sv_list: &[Sv], cache: &GeneCache, default_pre_gene_distance: i32, specific_pre_gene_distances: &JavaHashMap<String, i32>) {
    let mut specific_ranges: Vec<hmf_common::regions::ChrBaseRegion> = Vec::new();
    let mut specific_distance = default_pre_gene_distance;
    for (gene_name, d) in specific_pre_gene_distances.iter() {
        specific_distance = *d;
        let Some(g) = cache.gene_data_by_name(gene_name) else { continue };
        let start = if g.strand == 1 { g.gene_start - specific_distance } else { g.gene_end };
        let end = if g.strand == 1 { g.gene_start } else { g.gene_end + specific_distance };
        specific_ranges.push(hmf_common::regions::ChrBaseRegion::new(&g.chromosome, start, end));
    }
    let pre_gene_distance = |chr: &str, pos: i32| -> i32 {
        if specific_ranges.iter().any(|r| r.chromosome == chr && r.contains_position(pos)) { specific_distance } else { default_pre_gene_distance }
    };
    for var in sv_list {
        for be in [SE_START, SE_END] {
            let is_start = is_start(be);
            if be == SE_END && var.is_sgl_breakend() {
                for mapping in var.sgl_mappings() {
                    let d = pre_gene_distance(&mapping.chromosome, mapping.position);
                    let genes = find_gene_annotations_by_sv(cache, var.id(), is_start, &mapping.chromosome, mapping.position, mapping.orientation, d);
                    for g in &genes {
                        g.set_type(var.sv_type());
                    }
                    var.genes_mut()[se(is_start)].extend(genes);
                }
            } else {
                let d = pre_gene_distance(var.chromosome(is_start), var.position(is_start));
                let genes = find_gene_annotations_by_sv(
                    cache,
                    var.id(),
                    is_start,
                    var.chromosome(is_start),
                    var.position(is_start),
                    var.orientation(is_start),
                    d,
                );
                var.genes_mut()[se(is_start)].extend(genes);
                for g in var.genes_list(is_start) {
                    g.set_sv_data(var.sv_data(), var.jcn());
                }
            }
        }
    }
}

/// BreakendGenePrep.findGeneAnnotationsBySv
pub fn find_gene_annotations_by_sv(
    cache: &GeneCache, sv_id: i32, is_start: bool, chromosome: &str, position: i32, orientation: i8, upstream_distance: i32,
) -> Vec<BreakendGene> {
    let mut out = Vec::new();
    for gene in cache.find_gene_regions(chromosome, position, upstream_distance) {
        let Some(trans_list) = cache.transcripts(&gene.gene_id) else { continue };
        if trans_list.is_empty() {
            continue;
        }
        let current = BreakendGene::create(sv_id, is_start, gene);
        current.set_positional_data(chromosome, position, orientation);
        for trans_data in trans_list.iter() {
            if let Some(t) = create_breakend_transcript_data(trans_data, position, &current) {
                current.add_transcript(t.clone());
                set_alternative_transcript_phasings(&t, &trans_data.exons, position, orientation);
                if !t.has_prev_splice_acceptor_distance() {
                    set_preceding_gene_distance(cache, &t, position);
                }
            }
        }
        if current.transcript_count() == 0 && cache.has_downstream_gene_annotation(gene) {
            if let Some(trans_data) = trans_list.iter().find(|x| x.is_canonical) {
                let cb = calc_coding_bases(trans_data, position);
                let t = BreakendTrans::create(&current, trans_data, 1, 1, PHASE_NONE, PHASE_NONE, cb.coding_bases, cb.total_coding_bases);
                current.add_transcript(t);
            }
        }
        out.push(current);
    }
    let alt: Vec<Gene> = cache
        .alternative_gene_data
        .borrow()
        .iter()
        .copied()
        .filter(|g| g.chromosome == chromosome && position >= g.gene_start && position <= g.gene_end)
        .collect();
    for alt_gene in alt {
        let Some(trans_list) = cache.transcripts(&alt_gene.gene_id) else { continue };
        if trans_list.is_empty() {
            continue;
        }
        let trans = trans_list.iter().find(|x| x.is_canonical);
        let ga = BreakendGene::create(sv_id, is_start, alt_gene);
        ga.set_positional_data(chromosome, position, orientation);
        if let Some(trans) = trans {
            let cb = calc_coding_bases(trans, position);
            let t = BreakendTrans::create(&ga, trans, 1, 1, PHASE_NONE, PHASE_NONE, cb.coding_bases, cb.total_coding_bases);
            ga.add_transcript(t);
        }
        out.push(ga);
    }
    out
}

/// BreakendGenePrep.setAlternativeTranscriptPhasings
pub fn set_alternative_transcript_phasings(transcript: &BreakendTrans, exons: &[ExonData], position: i32, orientation: i8) {
    let is_upstream = (transcript.gene().strand() as i32 * orientation as i32) > 0;
    let forward_strand = transcript.trans_data.pos_strand();
    let mut alt = transcript.alternative_phasing_mut();
    alt.clear();
    let trans_phase = transcript.phase;
    let trans_rank = if is_upstream { transcript.exon_upstream } else { transcript.exon_downstream };
    for exon in exons {
        if is_upstream == forward_strand {
            if exon.start > position || trans_rank == exon.rank {
                break;
            }
        } else if position > exon.end || trans_rank == exon.rank {
            continue;
        }
        let mut exon_phase = if is_upstream { exon.phase_end } else { exon.phase_start };
        if exon_phase == PHASE_NONE && !transcript.trans_data.non_coding() {
            if forward_strand && exon.end > transcript.trans_data.coding_end.unwrap() {
                exon_phase = POST_CODING_PHASE;
            } else if !forward_strand && exon.start < transcript.trans_data.coding_start.unwrap() {
                exon_phase = POST_CODING_PHASE;
            }
        }
        let exons_skipped = if is_upstream { (trans_rank - exon.rank).max(0) } else { (exon.rank - trans_rank).max(0) };
        if exon_phase != trans_phase {
            if is_upstream == forward_strand {
                alt.insert(exon_phase, exons_skipped);
            } else if !alt.contains_key(&exon_phase) {
                alt.insert(exon_phase, exons_skipped);
            }
        }
    }
}

/// BreakendGenePrep.createBreakendTranscriptData
pub fn create_breakend_transcript_data(trans_data: Trans, position: i32, gene: &BreakendGene) -> Option<BreakendTrans> {
    let exons = &trans_data.exons;
    if exons.is_empty() {
        return None;
    }
    let is_forward_strand = gene.strand() == 1;
    let is_upstream = gene.is_upstream();
    let mut up_exon_rank = -1;
    let mut down_exon_rank = -1;
    let mut next_up_distance = -1;
    let mut next_down_distance = -1;
    let mut is_coding_type_override = false;
    let mut phase = PHASE_NONE;
    let first = &exons[0];
    let last = &exons[exons.len() - 1];
    if position < first.start {
        if is_forward_strand {
            let first_spa = if exons.len() > 1 { &exons[1] } else { first };
            down_exon_rank = first_spa.rank;
            next_down_distance = first_spa.start - position;
            is_coding_type_override = trans_data.coding_start.is_some_and(|cs| first_spa.start > cs);
            if let Some(cs) = trans_data.coding_start {
                if first_spa.start > cs {
                    is_coding_type_override = true;
                }
                phase = if first_spa.start == cs { PHASE_NONE } else { first_spa.phase_start };
            }
            up_exon_rank = 0;
        } else {
            return None;
        }
    } else if position > last.end {
        if !is_forward_strand {
            let first_spa = if exons.len() > 1 { &exons[exons.len() - 2] } else { last };
            down_exon_rank = first_spa.rank;
            next_down_distance = position - last.end;
            if let Some(ce) = trans_data.coding_end {
                if first_spa.end < ce {
                    is_coding_type_override = true;
                }
                phase = if first_spa.end == ce { PHASE_NONE } else { first_spa.phase_start };
            }
            up_exon_rank = 0;
        } else {
            return None;
        }
    } else {
        for index in 0..exons.len() {
            let exon = &exons[index];
            if position >= exon.start && position <= exon.end {
                up_exon_rank = exon.rank;
                down_exon_rank = exon.rank;
                if is_forward_strand {
                    next_up_distance = position - exon.start;
                    if index < exons.len() - 1 {
                        next_down_distance = exons[index + 1].start - position;
                    }
                } else {
                    next_up_distance = exon.end - position;
                    if index > 1 {
                        next_down_distance = position - exons[index - 1].end;
                    }
                }
                phase = if is_upstream { exon.phase_start } else { exon.phase_end };
                break;
            } else if position < exon.start {
                let prev = &exons[index - 1];
                if is_forward_strand {
                    up_exon_rank = prev.rank;
                    down_exon_rank = exon.rank;
                    next_down_distance = exon.start - position;
                    next_up_distance = position - prev.end;
                } else {
                    up_exon_rank = exon.rank;
                    down_exon_rank = prev.rank;
                    next_up_distance = exon.start - position;
                    next_down_distance = position - prev.end;
                }
                if is_upstream {
                    phase = if is_forward_strand { prev.phase_end } else { exon.phase_end };
                } else if is_forward_strand {
                    phase = if trans_data.coding_start == Some(exon.start) { PHASE_NONE } else { exon.phase_start };
                } else {
                    phase = if trans_data.coding_end == Some(prev.end) { PHASE_NONE } else { prev.phase_start };
                }
                break;
            }
        }
    }
    let cb = calc_coding_bases(trans_data, position);
    let t = BreakendTrans::create(gene, trans_data, up_exon_rank, down_exon_rank, phase, cb.phase, cb.coding_bases, cb.total_coding_bases);
    t.set_splice_acceptor_distance(true, if next_up_distance >= 0 { Some(next_up_distance) } else { None });
    t.set_splice_acceptor_distance(false, if next_down_distance >= 0 { Some(next_down_distance) } else { None });
    if is_coding_type_override {
        t.set_coding_type(TranscriptCodingType::CODING);
    }
    Some(t)
}

fn set_preceding_gene_distance(cache: &GeneCache, transcript: &BreakendTrans, position: i32) {
    let sa_pos = cache.ensembl.find_preceding_gene_splice_acceptor_position(transcript.trans_id());
    if sa_pos >= 0 {
        let d = if transcript.gene().strand() == 1 { position - sa_pos } else { sa_pos - position };
        transcript.set_splice_acceptor_distance(true, Some(d));
    }
}
