//! Package `fusion`: gene fusions (FusionFinder, SpecialFusions, FusionReportability), disruptions
//! (DisruptionFinder), the annotation driver (FusionDisruptionAnalyser) and the neo-epitope candidates
//! (NeoEpitopeWriter). The RNA fusion matching and the cohort writers are left out.

pub mod analyser;
pub mod disruption;
pub mod finder;
pub mod neo;

use std::cell::{Cell, RefCell};

use hmf_common::ensembl::BIOTYPE_LINC_RNA;
use hmf_common::ensembl::{BIOTYPE_NONSENSE_MED_DECAY, BIOTYPE_PROCESSED_TRANS, BIOTYPE_PROTEIN_CODING, BIOTYPE_RETAINED_INTRON};
use hmf_common::known_fusion::KnownFusionType::{self, *};
use hmf_common::known_fusion::{fs_index, FS_DOWN, FS_UP, PROMISCUOUS_BOTH};
use hmf_common::sv::StructuralVariantType::SGL;

use crate::gene::BreakendTrans;
use crate::types::handle;

// FusionConstants
pub const MAX_UPSTREAM_DISTANCE_KNOWN: i32 = 100000;
/// 250000 until upstream 7a4490caa6 (IG fusions, AUS-509)
pub const MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN: i32 = 1_000_000;
pub const MAX_UPSTREAM_DISTANCE_OTHER: i32 = 10000;
pub const MAX_ENHANCER_PARTNER_GENE_DISTANCE_NO_ORIENT: i32 = 50000;
pub const MIN_ENHANCER_PARTNER_GENE_DISTANCE: i32 = 1_000_000;
pub const SHORT_UNPHASED_DISTANCE_KNOWN: i32 = 1000000;
pub const REQUIRED_BIOTYPES: [&str; 5] =
    [BIOTYPE_PROCESSED_TRANS, BIOTYPE_PROTEIN_CODING, BIOTYPE_NONSENSE_MED_DECAY, BIOTYPE_RETAINED_INTRON, BIOTYPE_LINC_RNA];
pub const FUSION_MAX_CHAIN_LENGTH: i64 = 150000;
pub const FUSION_MAX_CHAIN_LINKS: i32 = 4;
pub const ENHANCER_PROMISCUOUS_MIN_DISTANCE: i32 = 100000;
/// the name prefix of the transcripts LINX makes for enhancer regions ("@IGH")
pub const ENHANCER_TRANSCRIPT_PREFIX: &str = "@";

/// first entries for V37, second for V38 in Ensembl
pub const PROTEINS_REQUIRED_LOST: [&str; 2] = ["Raf-like Ras-binding", "RBD"];
pub const PROTEINS_REQUIRED_KEPT: [&str; 15] = [
    "Ets domain",
    "ETS_DOMAIN_3",
    "Protein kinase domain",
    "PROTEIN_KINASE_DOM",
    "Epidermal growth factor-like domain",
    "EGF_3",
    "Ankyrin repeat-containing domain",
    "ANK_REPEAT",
    "ANK_REP_REGION",
    "Basic-leucine zipper domain",
    "BZIP",
    "High mobility group box domain",
    "HMG_BOX_2",
    "Bromodomain",
    "BROMODOMAIN_2",
];

/// FusionConfig (the options LINX's annotation reads; the logging-only ones are kept for completeness). Upstream
/// 7a4490caa6 removed AllowExonSkipping, RequirePhaseMatch (`-skip_unphased_fusions`) and RequireUpstreamBiotypes:
/// exon skipping is always allowed, unphased fusions are always kept and upstream biotypes always required.
#[derive(Clone, Debug, Default)]
pub struct FusionConfig {
    pub log_reportable_only: bool,
    pub log_all_potentials: bool,
    pub write_all_vis_fusions: bool,
    /// fusion_debug_sv_genes: gene and SV id pairs to limit fusions to
    pub debug_sv_genes: Vec<(String, i32)>,
}

impl FusionConfig {
    pub fn process_sv_gene_debug(&self, gene: &str, sv_id: i32) -> bool {
        self.debug_sv_genes.is_empty() || self.debug_sv_genes.iter().any(|(g, id)| g == gene && *id == sv_id)
    }
}

/// FusionReportableReason
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FusionReportableReason {
    OK,
    NOT_KNOWN,
    UNPHASED_NOT_KNOWN,
    UNPHASED_5P_UTR,
    UNPHASED_SHORT,
    SGL_NOT_KNOWN,
    PRE_GENE_DISTANCE,
    NMD,
    NEG_SPLICE_ACC_DISTANCE,
    EXON_SKIPPING,
    CHAIN_TERMINATED,
    NON_DISRUPTIVE_CHAIN,
    INVALID_TRAVERSAL,
    CHAIN_LINKS,
    PROTEIN_DOMAINS,
}
impl std::fmt::Display for FusionReportableReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}
use FusionReportableReason as FR;

/// FusionPhasedType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FusionPhasedType {
    INFRAME,
    SKIPPED_EXONS,
    OUT_OF_FRAME,
}

/// FusionLikelihoodType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FusionLikelihoodType {
    HIGH,
    LOW,
    NA,
}

/// FusionChainInfo
#[derive(Clone, Debug)]
pub struct FusionChainInfo {
    pub chain_id: i32,
    pub chain_links: i32,
    pub chain_length: i32,
    pub valid_traversal: bool,
    pub traversal_assembled: bool,
    pub non_disruptive: bool,
}

/// FusionAnnotations
#[derive(Clone, Debug)]
pub struct FusionAnnotations {
    pub cluster_id: i32,
    pub cluster_count: i32,
    pub resolved_type: String,
    pub chain_info: Option<FusionChainInfo>,
    pub terminated_up: bool,
    pub terminated_down: bool,
}

/// GeneFusion: an upstream and a downstream transcript fused by an SV or a chain
pub struct GeneFusionData {
    id: Cell<i32>,
    transcripts: [BreakendTrans; 2],
    is_reportable: Cell<bool>,
    reportable_reasons: RefCell<Vec<FusionReportableReason>>,
    phase_matched: bool,
    exons_skipped: Cell<[i32; 2]>,
    known_type: Cell<KnownFusionType>,
    is_promiscuous: Cell<[bool; 2]>,
    known_exons: Cell<bool>,
    high_impact_promiscuous: Cell<bool>,
    protein_features_set: Cell<bool>,
    annotations: RefCell<Option<FusionAnnotations>>,
    priority: Cell<f64>,
}

handle!(GeneFusion, GeneFusionData);

impl GeneFusion {
    pub fn create(up: &BreakendTrans, down: &BreakendTrans, phase_matched: bool) -> GeneFusion {
        GeneFusion::new(GeneFusionData {
            id: Cell::new(0),
            transcripts: [up.clone(), down.clone()],
            is_reportable: Cell::new(false),
            reportable_reasons: RefCell::new(Vec::new()),
            phase_matched,
            exons_skipped: Cell::new([0, 0]),
            known_type: Cell::new(NONE),
            is_promiscuous: Cell::new([false, false]),
            known_exons: Cell::new(false),
            high_impact_promiscuous: Cell::new(false),
            protein_features_set: Cell::new(false),
            annotations: RefCell::new(None),
            priority: Cell::new(0.0),
        })
    }
}

impl GeneFusionData {
    pub fn set_id(&self, id: i32) {
        self.id.set(id)
    }
    pub fn id(&self) -> i32 {
        self.id.get()
    }
    pub fn name(&self) -> String {
        if self.is_enhancer() {
            format!("{}_{}", self.transcripts[FS_UP].trans_name(), self.transcripts[FS_DOWN].gene_name())
        } else {
            format!("{}_{}", self.transcripts[FS_UP].gene_name(), self.transcripts[FS_DOWN].gene_name())
        }
    }
    pub fn sv_id(&self, is_upstream: bool) -> i32 {
        self.transcripts[fs_index(is_upstream)].gene().var_id()
    }
    pub fn transcripts(&self) -> &[BreakendTrans; 2] {
        &self.transcripts
    }
    pub fn upstream_trans(&self) -> &BreakendTrans {
        &self.transcripts[FS_UP]
    }
    pub fn downstream_trans(&self) -> &BreakendTrans {
        &self.transcripts[FS_DOWN]
    }
    pub fn gene_name(&self, fs: usize) -> String {
        if fs == FS_UP && self.is_enhancer() {
            self.transcripts[FS_UP].trans_name().to_string()
        } else {
            self.transcripts[fs].gene_name().to_string()
        }
    }
    pub fn reportable(&self) -> bool {
        self.is_reportable.get()
    }
    pub fn set_reportable(&self, t: bool) {
        self.is_reportable.set(t);
        self.reportable_reasons.borrow_mut().clear();
    }
    pub fn set_reportable_reasons(&self, reasons: &[FusionReportableReason]) {
        self.reportable_reasons.borrow_mut().extend_from_slice(reasons);
    }
    pub fn add_reportable_reason(&self, r: FusionReportableReason) {
        self.reportable_reasons.borrow_mut().push(r);
    }
    pub fn reportable_reasons(&self) -> Vec<FusionReportableReason> {
        let r = self.reportable_reasons.borrow();
        if r.is_empty() {
            return vec![FR::OK];
        }
        if self.known_type.get() == NONE {
            return vec![r[0]];
        }
        r.clone()
    }
    pub fn protein_features_set(&self) -> bool {
        self.protein_features_set.get()
    }
    pub fn set_protein_features_set(&self) {
        self.protein_features_set.set(true)
    }
    pub fn known_type_str(&self) -> String {
        let t = self.known_type.get();
        if t == PROMISCUOUS_5 || t == PROMISCUOUS_3 {
            let p = self.is_promiscuous.get();
            if p[FS_UP] && p[FS_DOWN] {
                return PROMISCUOUS_BOTH.to_string();
            }
        }
        t.name()
    }
    pub fn known_type(&self) -> KnownFusionType {
        self.known_type.get()
    }
    pub fn set_known_type(&self, t: KnownFusionType) {
        self.known_type.set(t)
    }
    pub fn is_promiscuous(&self) -> [bool; 2] {
        self.is_promiscuous.get()
    }
    pub fn set_promiscuous(&self, fs: usize) {
        let mut p = self.is_promiscuous.get();
        p[fs] = true;
        self.is_promiscuous.set(p);
    }
    pub fn is_enhancer(&self) -> bool {
        let t = self.known_type.get();
        t == ENHANCER_KNOWN_PAIR || t == ENHANCER_PROMISCUOUS
    }
    pub fn phase_matched(&self) -> bool {
        self.phase_matched
    }
    pub fn phase_type(&self) -> FusionPhasedType {
        if !self.phase_matched {
            return FusionPhasedType::OUT_OF_FRAME;
        }
        let e = self.exons_skipped.get();
        if e[FS_UP] > 0 || e[FS_DOWN] > 0 { FusionPhasedType::SKIPPED_EXONS } else { FusionPhasedType::INFRAME }
    }
    pub fn likelihood_type(&self) -> FusionLikelihoodType {
        if !self.is_reportable.get() {
            return FusionLikelihoodType::NA;
        }
        let t = self.known_type.get();
        if t == KNOWN_PAIR || t == EXON_DEL_DUP || t == ENHANCER_KNOWN_PAIR || self.known_exons.get() {
            FusionLikelihoodType::HIGH
        } else {
            FusionLikelihoodType::LOW
        }
    }
    pub fn set_known_exons(&self) {
        self.known_exons.set(true)
    }
    pub fn known_exons(&self) -> bool {
        self.known_exons.get()
    }
    pub fn set_high_impact_promiscuous(&self) {
        self.high_impact_promiscuous.set(true)
    }
    pub fn is_high_impact_promiscuous(&self) -> bool {
        self.high_impact_promiscuous.get()
    }
    pub fn same_chromosome(&self) -> bool {
        self.transcripts[FS_UP].gene().chromosome() == self.transcripts[FS_DOWN].gene().chromosome()
    }
    pub fn distance(&self) -> i32 {
        if !self.same_chromosome() {
            return 0;
        }
        (self.transcripts[FS_UP].gene().position() - self.transcripts[FS_DOWN].gene().position()).abs()
    }
    pub fn set_exons_skipped(&self, up: i32, down: i32) {
        self.exons_skipped.set([up, down]);
    }
    pub fn exons_skipped(&self, is_upstream: bool) -> i32 {
        self.exons_skipped.get()[fs_index(is_upstream)]
    }
    pub fn fused_exon(&self, is_upstream: bool) -> i32 {
        let up = &self.transcripts[FS_UP];
        let down = &self.transcripts[FS_DOWN];
        if self.is_exonic() {
            return if is_upstream { up.exon_upstream } else { down.exon_downstream };
        }
        let e = self.exons_skipped.get();
        if is_upstream {
            (up.next_splice_exon_rank() - e[FS_UP]).max(1)
        } else {
            (down.next_splice_exon_rank() + e[FS_DOWN]).min(down.exon_count())
        }
    }
    pub fn breakend_exon(&self, is_upstream: bool) -> i32 {
        let up = &self.transcripts[FS_UP];
        let down = &self.transcripts[FS_DOWN];
        if self.is_exonic() {
            return if is_upstream { up.exon_upstream } else { down.exon_downstream };
        }
        if is_upstream { up.next_splice_exon_rank() } else { down.next_splice_exon_rank() }
    }
    pub fn is_exonic(&self) -> bool {
        self.transcripts[FS_UP].is_exonic() && self.transcripts[FS_DOWN].is_exonic()
    }
    pub fn annotations(&self) -> Option<FusionAnnotations> {
        self.annotations.borrow().clone()
    }
    pub fn has_annotations(&self) -> bool {
        self.annotations.borrow().is_some()
    }
    pub fn set_annotations(&self, a: FusionAnnotations) {
        *self.annotations.borrow_mut() = Some(a);
    }
    pub fn set_priority(&self, p: f64) {
        self.priority.set(p)
    }
    pub fn priority(&self) -> f64 {
        self.priority.get()
    }
    pub fn is_terminated(&self) -> bool {
        self.annotations.borrow().as_ref().is_some_and(|a| a.terminated_up || a.terminated_down)
    }
    pub fn chain_length(&self) -> i32 {
        self.annotations.borrow().as_ref().and_then(|a| a.chain_info.as_ref()).map(|c| c.chain_length).unwrap_or(0)
    }
    pub fn chain_links(&self) -> i32 {
        self.annotations.borrow().as_ref().and_then(|a| a.chain_info.as_ref()).map(|c| c.chain_links).unwrap_or(0)
    }
    pub fn valid_chain_traversal(&self) -> bool {
        self.annotations.borrow().as_ref().and_then(|a| a.chain_info.as_ref()).map(|c| c.valid_traversal).unwrap_or(true)
    }
    pub fn non_disruptive_chain(&self) -> bool {
        self.annotations.borrow().as_ref().and_then(|a| a.chain_info.as_ref()).map(|c| c.non_disruptive).unwrap_or(false)
    }
    pub fn sv_id_pair(&self) -> String {
        format!("{}_{}", self.transcripts[FS_UP].gene().var_id(), self.transcripts[FS_DOWN].gene().var_id())
    }
}

impl std::fmt::Display for GeneFusionData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let up = self.transcripts[FS_UP].gene();
        let down = self.transcripts[FS_DOWN].gene();
        let reasons: Vec<String> = self.reportable_reasons.borrow().iter().map(|r| r.to_string()).collect();
        write!(
            f,
            "{}: {} type({}) reportable({} reason=[{}]) phased({}) SVs({} & {}) up({}:{}:{} exon={}) down({}:{}:{} exon={})",
            self.id(),
            self.name(),
            self.known_type_str(),
            self.reportable(),
            reasons.join(", "),
            self.phase_matched,
            up.var_id(),
            down.var_id(),
            up.chromosome(),
            up.orientation(),
            up.position(),
            self.fused_exon(true),
            down.chromosome(),
            down.orientation(),
            down.position(),
            self.fused_exon(false)
        )
    }
}

// FusionReportability

/// determineReportability: the reasons a fusion is not reportable (none if it is)
pub fn determine_reportability(fusion: &GeneFusionData) -> Vec<FusionReportableReason> {
    let kt = fusion.known_type();
    if kt == NONE || kt == ENHANCER_PROMISCUOUS {
        return vec![FR::NOT_KNOWN];
    }
    let mut reasons = Vec::new();
    let up = fusion.upstream_trans();
    let down = fusion.downstream_trans();
    if !fusion.phase_matched() {
        if kt == EXON_DEL_DUP && fusion.is_exonic() {
            // allowed for specific exon-to-exon fusions
        } else {
            if kt != KNOWN_PAIR {
                reasons.push(FR::UNPHASED_NOT_KNOWN);
            }
            if !up.non_coding() && !up.pre_coding() {
                reasons.push(FR::UNPHASED_5P_UTR);
            }
            if fusion.same_chromosome() && fusion.distance() <= SHORT_UNPHASED_DISTANCE_KNOWN {
                reasons.push(FR::UNPHASED_SHORT);
            }
        }
    }
    if (up.gene().sv_type() == Some(SGL) || down.gene().sv_type() == Some(SGL)) && kt != KNOWN_PAIR && kt != ENHANCER_KNOWN_PAIR {
        reasons.push(FR::SGL_NOT_KNOWN);
    }
    let max_upstream = max_upstream_distance(fusion);
    if up.distance_upstream() > max_upstream {
        reasons.push(FR::PRE_GENE_DISTANCE);
    }
    if !fusion.is_high_impact_promiscuous() && down.distance_upstream() > max_upstream {
        reasons.push(FR::PRE_GENE_DISTANCE);
    }
    if down.bio_type() == BIOTYPE_NONSENSE_MED_DECAY {
        reasons.push(FR::NMD);
    }
    if !fusion.is_enhancer() && down.has_negative_prev_splice_acceptor_distance() {
        reasons.push(FR::NEG_SPLICE_ACC_DISTANCE);
    }
    if !permitted_exon_skipping(fusion) {
        reasons.push(FR::EXON_SKIPPING);
    }
    if fusion.is_terminated() && !allow_suspect_chains(kt) && !fusion.is_high_impact_promiscuous() {
        reasons.push(FR::CHAIN_TERMINATED);
    }
    if fusion.non_disruptive_chain() {
        reasons.push(FR::NON_DISRUPTIVE_CHAIN);
    }
    if !fusion.valid_chain_traversal() && !allow_suspect_chains(kt) {
        reasons.push(FR::INVALID_TRAVERSAL);
    }
    if fusion.chain_links() > FUSION_MAX_CHAIN_LINKS {
        reasons.push(FR::CHAIN_LINKS);
    }
    reasons
}

pub fn is_reportable(fusion: &GeneFusionData) -> bool {
    determine_reportability(fusion).is_empty()
}

fn max_upstream_distance(fusion: &GeneFusionData) -> i32 {
    let kt = fusion.known_type();
    if kt == ENHANCER_KNOWN_PAIR {
        MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN
    } else if kt == KNOWN_PAIR || fusion.is_high_impact_promiscuous() {
        MAX_UPSTREAM_DISTANCE_KNOWN
    } else if kt == ENHANCER_PROMISCUOUS {
        MAX_ENHANCER_PARTNER_GENE_DISTANCE_NO_ORIENT
    } else {
        MAX_UPSTREAM_DISTANCE_OTHER
    }
}

pub fn allow_suspect_chains(t: KnownFusionType) -> bool {
    t == KNOWN_PAIR || t == EXON_DEL_DUP || t == ENHANCER_KNOWN_PAIR
}

/// findTopPriorityFusion: the highest priority score, the first of equal scores
pub fn find_top_priority_fusion(fusions: &[GeneFusion]) -> Option<GeneFusion> {
    let mut top: Option<GeneFusion> = None;
    let mut highest = 0.0;
    for fusion in fusions {
        let score = calc_fusion_priority(fusion);
        fusion.set_priority(score);
        if score > highest {
            top = Some(fusion.clone());
            highest = score;
        }
    }
    top
}

fn calc_fusion_priority(fusion: &GeneFusionData) -> f64 {
    let up = fusion.upstream_trans();
    let down = fusion.downstream_trans();
    let mut score = 0.0;
    let mut factor = 1000000.0;
    if fusion.reportable() {
        score += factor;
    }
    factor /= 10.0;
    if is_high_priority_type(fusion.known_type()) {
        score += factor;
    }
    factor /= 10.0;
    if fusion.phase_matched() {
        score += factor;
    }
    factor /= 10.0;
    if !fusion.is_terminated() && fusion.valid_chain_traversal() {
        score += factor;
    }
    factor /= 10.0;
    if down.bio_type() == BIOTYPE_PROTEIN_CODING {
        score += factor;
    }
    factor /= 10.0;
    if fusion.exons_skipped(true) == 0 && fusion.exons_skipped(false) == 0 {
        score += factor;
    }
    factor /= 10.0;
    if down.is_canonical() {
        score += factor;
    }
    factor /= 10.0;
    if down.bio_type() != BIOTYPE_NONSENSE_MED_DECAY {
        score += factor;
    }
    factor /= 100.0;
    // round(length/10) of an int: the integer quotient
    let mut length = if down.is_coding() { down.coding_bases } else { down.exon_count() };
    length = (length / 10).min(99);
    score += length as f64 * factor;
    factor /= 10.0;
    if up.is_canonical() {
        score += factor;
    }
    factor /= 10.0;
    if up.bio_type() == BIOTYPE_PROTEIN_CODING {
        score += factor;
    }
    factor /= 100.0;
    length = if up.is_coding() { up.coding_bases } else { up.exon_count() };
    length = (length / 10).min(99);
    score += length as f64 * factor;
    factor /= 10.0;
    score += up.gene().jcn() * factor;
    score
}

pub fn is_high_priority_type(t: KnownFusionType) -> bool {
    t == KNOWN_PAIR || t == ENHANCER_KNOWN_PAIR || t == EXON_DEL_DUP || t == PROMISCUOUS_ENHANCER_TARGET
}

pub fn check_protein_domains(t: KnownFusionType) -> bool {
    t == PROMISCUOUS_5 || t == PROMISCUOUS_3 || t == KNOWN_PAIR
}

pub fn valid_protein_domains(fusion: &GeneFusionData) -> bool {
    let down = fusion.downstream_trans();
    let lost = down.protein_features_lost();
    let kept = down.protein_features_kept();
    let kept_but_lost = PROTEINS_REQUIRED_KEPT.iter().filter(|f| lost.contains(*f)).count();
    let lost_but_kept = PROTEINS_REQUIRED_LOST.iter().filter(|f| kept.contains(*f)).count();
    kept_but_lost == 0 && lost_but_kept == 0
}

fn permitted_exon_skipping(fusion: &GeneFusionData) -> bool {
    let kt = fusion.known_type();
    if kt == KNOWN_PAIR || kt == ENHANCER_KNOWN_PAIR || kt == ENHANCER_PROMISCUOUS || kt == EXON_DEL_DUP {
        return true;
    }
    if fusion.is_high_impact_promiscuous() || fusion.known_exons() {
        return true;
    }
    fusion.exons_skipped(true) == 0 && fusion.exons_skipped(false) == 0
}

/// SvDisruptionData: a disrupted transcript of a driver gene
pub struct SvDisruptionData {
    pub var: crate::types::Sv,
    pub is_start: bool,
    pub gene: crate::gene::Gene,
    pub transcript: crate::gene::Trans,
    pub exons: [i32; 2],
    pub coding_type: crate::gene::TranscriptCodingType,
    pub region_type: crate::gene::TranscriptRegionType,
    pub undisrupted_copy_number: f64,
    pub max_copy_number: f64,
    reportable: Cell<bool>,
    pseudogene_deletion: Cell<bool>,
}

impl SvDisruptionData {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        var: &crate::types::Sv, is_start: bool, gene: crate::gene::Gene, transcript: crate::gene::Trans, exons: [i32; 2],
        coding_type: crate::gene::TranscriptCodingType, region_type: crate::gene::TranscriptRegionType, undisrupted_copy_number: f64,
        max_copy_number: f64,
    ) -> SvDisruptionData {
        SvDisruptionData {
            var: var.clone(),
            is_start,
            gene,
            transcript,
            exons,
            coding_type,
            region_type,
            undisrupted_copy_number,
            max_copy_number,
            reportable: Cell::new(true),
            pseudogene_deletion: Cell::new(false),
        }
    }
    pub fn reportable(&self) -> bool {
        self.reportable.get()
    }
    pub fn set_reportable(&self, t: bool) {
        self.reportable.set(t)
    }
    pub fn is_pseudogene_deletion(&self) -> bool {
        self.pseudogene_deletion.get()
    }
    pub fn mark_pseudogene_deletion(&self) {
        self.pseudogene_deletion.set(true)
    }
}
