//! LINX's Java unit tests (hmftools/linx/src/test), ported with their synthetic fixtures. The test
//! utilities (LinxTester, SvTestUtils, GeneTestUtils and hmf-common's GeneTestUtils) are here.
//!
//! Left out: the visualiser tests (not ported), RnaFusionTest (RNA matching is not ported) and
//! ChainingActualTest, whose fixtures are real sample data (the port's fixtures are synthetic).

mod analyser;
mod chaining;
mod clustering;
mod drivers;
mod fusion;
mod germline;
mod misc;
mod upstream;

use std::rc::Rc;

use hmf_common::ensembl::{EnsemblDataCache, ExonData, GeneData, TranscriptData};
use hmf_common::regions::RefGenomeVersion;
use hmf_common::sv::StructuralVariantData;
use hmf_common::sv::StructuralVariantType::{self, *};
use hmf_common::transcript_utils::{coding_bases_to_phase, tick_phase_forward, PHASE_NONE};

use crate::analysis::clustering::{link_sgl_mapped_inferreds, populate_chromosome_breakend_map};
use crate::analysis::ClusterAnalyser;
use crate::annotators::LineElementAnnotator;
use crate::cn::{set_sv_copy_number_data, CnDataLoader, CnSegmentBuilder, Loh, LohEvent};
use crate::gene::{BreakendGene, BreakendTrans, GeneCache};
use crate::types::*;
use crate::utils::get_chromosomal_arm;
use crate::vis::VisSampleData;

pub const PHASE_0: i32 = 0;
pub const PHASE_1: i32 = 1;
pub const PHASE_2: i32 = 2;

// SvTestUtils

pub fn set_assembled_link_info(var1: &Sv, use_start1: bool, var2: &Sv, use_start2: bool) {
    let id2 = if use_start2 { var2.sv_data().vcf_id_start.clone() } else { var2.sv_data().vcf_id_end.clone() };
    let id1 = if use_start1 { var1.sv_data().vcf_id_start.clone() } else { var1.sv_data().vcf_id_end.clone() };
    var1.add_assembly_info(use_start1, &id2);
    var2.add_assembly_info(use_start2, &id1);
}

#[allow(clippy::too_many_arguments)]
pub fn create_sv(id: i32, chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, o_start: i32, o_end: i32, t: StructuralVariantType, ins: &str) -> Sv {
    create_test_sv_full(id, chr_start, chr_end, pos_start, pos_end, o_start, o_end, t, 2.0, 2.0, 1.0, 1.0, 1.0, ins, "PASS", "", "")
}
pub fn create_del(id: i32, chr: &str, s: i32, e: i32) -> Sv {
    create_test_sv_full(id, chr, chr, s, e, 1, -1, DEL, 2.0, 2.0, 1.0, 1.0, 1.0, "", "PASS", "", "")
}
pub fn create_ins(id: i32, chr: &str, s: i32, e: i32) -> Sv {
    create_test_sv_full(id, chr, chr, s, e, 1, -1, INS, 2.0, 2.0, 1.0, 1.0, 1.0, "", "PASS", "", "")
}
pub fn create_dup(id: i32, chr: &str, s: i32, e: i32) -> Sv {
    create_test_sv_full(id, chr, chr, s, e, -1, 1, DUP, 3.0, 3.0, 1.0, 1.0, 1.0, "", "PASS", "", "")
}
pub fn create_inv(id: i32, chr: &str, s: i32, e: i32, o: i32) -> Sv {
    let (cs, ce) = if o == 1 { (4.0, 3.0) } else { (3.0, 4.0) };
    create_test_sv_full(id, chr, chr, s, e, o, o, INV, cs, ce, 1.0, 1.0, 1.0, "", "PASS", "", "")
}
pub fn create_sgl(id: i32, chr: &str, pos: i32, o: i32) -> Sv {
    create_test_sv_full(id, chr, "0", pos, -1, o, -1, SGL, 3.0, 0.0, 1.0, 0.0, 1.0, "", "PASS", "", "")
}
pub fn create_inf(id: i32, chr: &str, pos: i32, o: i32) -> Sv {
    create_test_sv_full(id, chr, "0", pos, -1, o, -1, INF, 3.0, 0.0, 1.0, 0.0, 1.0, "", "PASS", "", "")
}
pub fn create_bnd(id: i32, chr_start: &str, pos_start: i32, o_start: i32, chr_end: &str, pos_end: i32, o_end: i32) -> Sv {
    create_test_sv_full(id, chr_start, chr_end, pos_start, pos_end, o_start, o_end, BND, 3.0, 3.0, 1.0, 1.0, 1.0, "", "PASS", "", "")
}
/// createTestSv(..., jcn): the copy number data from the test copy number routine
#[allow(clippy::too_many_arguments)]
pub fn create_test_sv(id: i32, chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, o_start: i32, o_end: i32, t: StructuralVariantType, jcn: f64) -> Sv {
    create_test_sv_full(id, chr_start, chr_end, pos_start, pos_end, o_start, o_end, t, 0.0, 0.0, jcn, jcn, jcn, "", "PASS", "", "")
}
#[allow(clippy::too_many_arguments)]
pub fn create_test_sv_cn(
    id: i32, chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, o_start: i32, o_end: i32, t: StructuralVariantType, cn_start: f64,
    cn_end: f64, cn_chg_start: f64, cn_chg_end: f64, jcn: f64, ins: &str,
) -> Sv {
    create_test_sv_full(id, chr_start, chr_end, pos_start, pos_end, o_start, o_end, t, cn_start, cn_end, cn_chg_start, cn_chg_end, jcn, ins, "PASS", "", "")
}

#[allow(clippy::too_many_arguments)]
pub fn create_test_sv_full(
    id: i32, chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, o_start: i32, o_end: i32, t: StructuralVariantType, cn_start: f64,
    cn_end: f64, cn_chg_start: f64, cn_chg_end: f64, jcn: f64, ins: &str, filter: &str, repeat_class: &str, repeat_type: &str,
) -> Sv {
    let d = StructuralVariantData {
        id,
        vcf_id_start: format!("{id}_start"),
        vcf_id_end: format!("{id}_end"),
        start_chromosome: chr_start.to_string(),
        end_chromosome: chr_end.to_string(),
        start_position: pos_start,
        end_position: pos_end,
        start_orientation: o_start as i8,
        end_orientation: o_end as i8,
        start_af: 1.0,
        end_af: 1.0,
        junction_copy_number: jcn,
        adjusted_start_af: 1.0,
        adjusted_end_af: 1.0,
        adjusted_start_copy_number: cn_start,
        adjusted_end_copy_number: cn_end,
        adjusted_start_copy_number_change: cn_chg_start,
        adjusted_end_copy_number_change: cn_chg_end,
        insert_sequence: ins.to_string(),
        sv_type: Some(t),
        filter: filter.to_string(),
        start_tumor_variant_fragment_count: 10,
        start_tumor_reference_fragment_count: 10,
        start_normal_variant_fragment_count: 10,
        start_normal_reference_fragment_count: 10,
        end_tumor_variant_fragment_count: 10,
        end_tumor_reference_fragment_count: 10,
        end_normal_variant_fragment_count: 10,
        end_normal_reference_fragment_count: 10,
        insert_sequence_repeat_class: repeat_class.to_string(),
        insert_sequence_repeat_type: repeat_type.to_string(),
        ..Default::default()
    };
    let var = Sv::from_data(Rc::new(d));
    initialise_sv(&var);
    var
}

pub fn initialise_sv(var: &Sv) {
    let start_arm = get_chromosomal_arm(var.chromosome(true), var.position(true));
    let end_arm = if !var.is_sgl_breakend() { get_chromosomal_arm(var.chromosome(false), var.position(false)) } else { Arm::P };
    var.set_chromosomal_arms(start_arm, end_arm);
    let jcn = var.sv_data().junction_copy_number;
    var.set_jcn_recalc_data(jcn, jcn);
}

// LinxTester

pub struct LinxTester {
    pub sample_id: String,
    pub all_variants: Vec<Sv>,
    pub analyser: ClusterAnalyser,
    pub cn_loader: CnDataLoader,
    pub line_annotator: LineElementAnnotator,
    pub vis: VisSampleData,
    pub is_germline: bool,
    cn_builder: CnSegmentBuilder,
    next_var_id: i32,
}

impl LinxTester {
    pub fn new() -> Self {
        Self::with_germline(false)
    }
    pub fn with_germline(is_germline: bool) -> Self {
        crate::utils::set_ref_genome_v37(true);
        let mut analyser = ClusterAnalyser::new(DEFAULT_PROXIMITY_DISTANCE, is_germline, 0, true, false);
        analyser.run_validation_checks = true;
        LinxTester {
            sample_id: "TEST".into(),
            all_variants: Vec::new(),
            analyser,
            cn_loader: CnDataLoader::new(""),
            line_annotator: LineElementAnnotator::new(DEFAULT_PROXIMITY_DISTANCE),
            vis: VisSampleData::default(),
            is_germline,
            cn_builder: CnSegmentBuilder::new(),
            next_var_id: 0,
        }
    }
    pub fn next_var_id(&mut self) -> i32 {
        let id = self.next_var_id;
        self.next_var_id += 1;
        id
    }
    pub fn add(&mut self, vars: &[&Sv]) {
        self.all_variants.extend(vars.iter().map(|v| (*v).clone()));
    }
    pub fn add_and_cluster(&mut self, var1: &Sv, var2: &Sv) {
        self.clear_clusters_and_svs();
        self.add(&[var1, var2]);
        self.pre_clustering_init();
        self.cluster_and_analyse();
    }
    pub fn pre_clustering_init(&mut self) {
        self.pre_clustering_init_ploidy(false);
    }
    pub fn pre_clustering_init_ploidy(&mut self, include_ploidy_calcs: bool) {
        self.analyser.state.reset();
        link_sgl_mapped_inferreds(&mut self.all_variants);
        for v in self.all_variants.iter().filter(|x| x.linked_svs().is_some() && !x.is_sgl_breakend()) {
            initialise_sv(v);
        }
        populate_chromosome_breakend_map(&self.all_variants, &mut self.analyser.state);
        if !self.is_germline {
            self.populate_copy_number_data(include_ploidy_calcs);
        }
        let all = self.all_variants.clone();
        self.analyser.set_sample_data(&self.sample_id, &all, &self.cn_loader);
        self.analyser.pre_clustering_preparation();
    }
    pub fn cluster_and_analyse(&mut self) -> bool {
        let tc = self.cn_loader.chr_ends_cn_map.clone();
        self.analyser.cluster_and_analyse(&self.line_annotator, None, &tc)
    }
    pub fn add_loh_event(&mut self, b1: &Be, b2: &Be) {
        if b1.reverse_orient() || b2.forward_orient() || b1.chromosome() != b2.chromosome() {
            return;
        }
        let loh = LohEvent::new(b1.chromosome(), b1.position(), b2.position(), b1.sv().type_str(), b2.sv().type_str(), 1, b1.sv().id(), b2.sv().id());
        self.cn_loader.loh_event_data.push(Loh::new(loh));
    }
    pub fn clear_clusters_and_svs(&mut self) {
        for v in &self.all_variants {
            v.set_cluster(None);
        }
        self.all_variants.clear();
        self.analyser.state.reset();
        self.analyser.clusters.clear();
    }
    pub fn clusters(&self) -> &Vec<Cluster> {
        &self.analyser.clusters
    }
    pub fn all_clusters(&self) -> Vec<Cluster> {
        self.analyser.all_clusters()
    }
    pub fn has_cluster_with_svs(&self, svs: &[&Sv]) -> bool {
        self.find_cluster_with_svs(svs).is_some()
    }
    pub fn find_cluster_with_svs(&self, svs: &[&Sv]) -> Option<Cluster> {
        self.all_clusters().into_iter().find(|c| c.sv_count() == svs.len() && svs.iter().all(|v| contains(&c.svs(), v)))
    }
    pub fn find_chain_with_svs(&self, cluster: &Cluster, svs: &[&Sv]) -> Option<Chain> {
        cluster.chains().into_iter().find(|c| c.sv_list().len() == svs.len() && svs.iter().all(|v| contains(&c.sv_list(), v)))
    }
    pub fn set_non_cluster_allele_ploidies(&mut self, other: f64, undisrupted: f64) {
        self.cn_builder.set_allele_ploidies(other, undisrupted);
    }
    pub fn populate_copy_number_data(&mut self, include_ploidy_calcs: bool) {
        self.cn_builder.create_copy_number_data(&mut self.cn_loader, &self.analyser.state.chr_breakend_map);
        if include_ploidy_calcs {
            let s = self.sample_id.clone();
            self.cn_loader.calculate_adjusted_jcn(&s);
        }
        CnSegmentBuilder::set_sample_purity(&mut self.cn_loader, 1.0, 2.0, "MALE");
        self.cn_loader.create_chr_copy_number_map();
        set_sv_copy_number_data(&self.all_variants, Some(&self.cn_loader.sv_jcn_calc_map), &self.cn_loader.sv_id_cn_data_map, &self.cn_loader.chr_cn_data_map);
    }
}

// hmf-common GeneTestUtils

pub const CHR_1: &str = "1";
pub const CHR_2: &str = "2";
pub const CHR_3: &str = "3";
pub const POS_STRAND: i8 = 1;
pub const NEG_STRAND: i8 = -1;

pub fn create_gene_data_cache() -> EnsemblDataCache {
    EnsemblDataCache::new("", RefGenomeVersion::V37)
}
pub fn create_ensembl_gene_data(gene_id: &str, gene_name: &str, chromosome: &str, strand: i8, start: i32, end: i32) -> GeneData {
    GeneData {
        gene_id: gene_id.into(),
        gene_name: gene_name.into(),
        chromosome: chromosome.into(),
        strand,
        gene_start: start,
        gene_end: end,
        karyotype_band: String::new(),
        synonyms: String::new(),
    }
}
pub fn add_trans_exon_data(cache: &mut EnsemblDataCache, gene_id: &str, list: Vec<TranscriptData>) {
    cache.transcript_data_map_mut().insert(gene_id.to_string(), list);
}
pub fn add_gene_data(cache: &mut EnsemblDataCache, chromosome: &str, list: Vec<GeneData>) {
    cache.chr_gene_data_map_mut().insert(chromosome.to_string(), list);
}
pub fn generate_exon_starts(start: i32, count: i32, exon_length: i32, intron_length: i32) -> Vec<i32> {
    (0..count).map(|i| start + i * (exon_length + intron_length)).collect()
}
pub fn generate_trans_name(trans_id: i32) -> String {
    format!("TRAN{trans_id:04}")
}

/// createTransExons, with its phases (and the negative-strand loop's `finishedCoding = false`)
#[allow(clippy::too_many_arguments)]
pub fn create_trans_exons(
    gene_id: &str, trans_id: i32, strand: i8, exon_starts: &[i32], exon_length: i32, coding_start: Option<i32>, coding_end: Option<i32>,
    canonical: bool, biotype: &str,
) -> TranscriptData {
    let n = exon_starts.len();
    let trans_start = exon_starts[0];
    let trans_end = exon_starts[n - 1] + exon_length;
    let mut exons: Vec<ExonData> = Vec::new();
    let has_coding = coding_start.is_some() && coding_end.is_some();
    let (cs, ce) = (coding_start.unwrap_or(0), coding_end.unwrap_or(0));
    let mut in_coding = false;
    let mut finished = false;
    let mut last_end_phase = PHASE_NONE;
    let exon = |s: i32, e: i32, r: i32, ps: i32, pe: i32| ExonData { trans_id, start: s, end: e, rank: r, phase_start: ps, phase_end: pe };
    if strand == POS_STRAND {
        for (i, &exon_start) in exon_starts.iter().enumerate() {
            let exon_end = exon_start + exon_length;
            let rank = i as i32 + 1;
            let mut phase = PHASE_0;
            let mut start_phase = PHASE_NONE;
            let mut end_phase = PHASE_NONE;
            let mut exon_coding_start = 0;
            if has_coding && !finished {
                if !in_coding {
                    if cs <= exon_end {
                        in_coding = true;
                        phase = PHASE_1;
                        start_phase = if cs == exon_start { PHASE_0 } else { PHASE_NONE };
                        exon_coding_start = cs;
                    }
                } else {
                    phase = tick_phase_forward(last_end_phase, 1);
                    start_phase = last_end_phase;
                    exon_coding_start = exon_start;
                }
                if in_coding {
                    let bases = exon_end.min(ce) - exon_coding_start + 1;
                    end_phase = coding_bases_to_phase(bases, phase);
                    last_end_phase = end_phase;
                    if ce <= exon_end {
                        finished = true;
                        in_coding = false;
                        end_phase = PHASE_NONE;
                    }
                }
            }
            exons.push(exon(exon_start, exon_end, rank, start_phase, end_phase));
        }
    } else {
        for i in (0..n).rev() {
            let exon_start = exon_starts[i];
            let exon_end = exon_start + exon_length;
            let rank = (n - i) as i32;
            let mut phase = PHASE_0;
            let mut start_phase = PHASE_NONE;
            let mut end_phase = PHASE_NONE;
            let mut exon_coding_end = 0;
            if has_coding && !finished {
                if !in_coding {
                    if ce >= exon_start {
                        in_coding = true;
                        phase = PHASE_1;
                        start_phase = if ce == exon_end { PHASE_0 } else { PHASE_NONE };
                        exon_coding_end = ce;
                    }
                } else {
                    phase = tick_phase_forward(last_end_phase, 1);
                    start_phase = last_end_phase;
                    exon_coding_end = exon_end;
                }
                if in_coding {
                    let bases = exon_coding_end - exon_start.max(cs) + 1;
                    end_phase = coding_bases_to_phase(bases, phase);
                    last_end_phase = end_phase;
                    if cs >= exon_start {
                        finished = false;
                        in_coding = false;
                        end_phase = PHASE_NONE;
                    }
                }
            }
            exons.insert(0, exon(exon_start, exon_end, rank, start_phase, end_phase));
        }
    }
    TranscriptData {
        trans_id,
        trans_name: generate_trans_name(trans_id),
        gene_id: gene_id.into(),
        is_canonical: canonical,
        strand,
        trans_start,
        trans_end,
        coding_start,
        coding_end,
        bio_type: biotype.into(),
        ref_seq_id: None,
        exons,
    }
}

/// a leaked gene cache over a finished Ensembl cache (the analyser holds references for the run)
pub fn leak_cache(cache: EnsemblDataCache) -> &'static GeneCache {
    let e: &'static EnsemblDataCache = Box::leak(Box::new(cache));
    Box::leak(Box::new(GeneCache::new(e)))
}

// LINX GeneTestUtils

/// createGeneAnnotation: a breakend gene with its own gene data (no range)
#[allow(clippy::too_many_arguments)]
pub fn create_gene_annotation(sv_id: i32, is_start: bool, gene_name: &str, gene_id: &str, strand: i8, chromosome: &str, position: i32, orientation: i32) -> BreakendGene {
    let g: &'static GeneData = Box::leak(Box::new(create_ensembl_gene_data(gene_id, gene_name, chromosome, strand, 0, 0)));
    let gene = BreakendGene::create(sv_id, is_start, g);
    gene.set_positional_data(chromosome, position, orientation as i8);
    gene
}

/// createTranscript
#[allow(clippy::too_many_arguments)]
pub fn create_transcript(
    gene: &BreakendGene, trans_id: i32, canonical: bool, trans_start: i32, trans_end: i32, coding_start: Option<i32>, coding_end: Option<i32>,
    bio_type: &str, exon_up: i32, exon_down: i32, phase: i32, coding_bases: i32, total_coding_bases: i32,
) -> BreakendTrans {
    create_transcript_ex(gene, trans_id, canonical, trans_start, trans_end, coding_start, coding_end, bio_type, exon_up, exon_down, phase, PHASE_NONE, coding_bases, total_coding_bases)
}

#[allow(clippy::too_many_arguments)]
pub fn create_transcript_ex(
    gene: &BreakendGene, trans_id: i32, canonical: bool, trans_start: i32, trans_end: i32, coding_start: Option<i32>, coding_end: Option<i32>,
    bio_type: &str, exon_up: i32, exon_down: i32, phase: i32, exonic_phase: i32, coding_bases: i32, total_coding_bases: i32,
) -> BreakendTrans {
    let t: &'static TranscriptData = Box::leak(Box::new(TranscriptData {
        trans_id,
        trans_name: format!("TRANS{trans_id:03}"),
        gene_id: gene.gene_id().to_string(),
        is_canonical: canonical,
        strand: gene.strand(),
        trans_start,
        trans_end,
        coding_start,
        coding_end,
        bio_type: bio_type.into(),
        ref_seq_id: None,
        exons: Vec::new(),
    }));
    BreakendTrans::create(gene, t, exon_up, exon_down, phase, exonic_phase, coding_bases, total_coding_bases)
}

// LINX GeneTestUtils: six genes on chromosomes 1 and 2

pub const GENE_NAME_1: &str = "GENE1";
pub const GENE_ID_1: &str = "ENSG0001";
pub const GENE_START_1: i32 = 1000;
pub const GENE_NAME_2: &str = "GENE2";
pub const GENE_ID_2: &str = "ENSG0002";
pub const GENE_START_2: i32 = 10000;
pub const GENE_NAME_3: &str = "GENE3";
pub const GENE_ID_3: &str = "ENSG0003";
pub const GENE_START_3: i32 = 20000;
pub const GENE_NAME_4: &str = "GENE4";
pub const GENE_ID_4: &str = "ENSG0004";
pub const GENE_START_4: i32 = 1000;
pub const GENE_NAME_5: &str = "GENE5";
pub const GENE_ID_5: &str = "ENSG0005";
pub const GENE_START_5: i32 = 10000;
pub const GENE_NAME_6: &str = "GENE6";
pub const GENE_ID_6: &str = "ENSG0006";
pub const GENE_START_6: i32 = 10400;
pub const EXON_LENGTH: i32 = 100;
pub const TRANS_1: i32 = 1;
pub const TRANS_2: i32 = 2;
pub const TRANS_3: i32 = 3;
pub const TRANS_4: i32 = 4;
pub const TRANS_5: i32 = 5;
/// as in Java, the same id as TRANS_5
pub const TRANS_6: i32 = 5;

pub fn add_test_genes(cache: &mut EnsemblDataCache) {
    add_gene_data(
        cache,
        CHR_1,
        vec![
            create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, GENE_START_1, GENE_START_1 + 1000),
            create_ensembl_gene_data(GENE_ID_2, GENE_NAME_2, CHR_1, POS_STRAND, GENE_START_2, GENE_START_2 + 1000),
            create_ensembl_gene_data(GENE_ID_3, GENE_NAME_3, CHR_1, NEG_STRAND, GENE_START_3, GENE_START_3 + 1000),
        ],
    );
    add_gene_data(
        cache,
        CHR_2,
        vec![
            create_ensembl_gene_data(GENE_ID_4, GENE_NAME_4, CHR_2, POS_STRAND, GENE_START_4, GENE_START_4 + 1000),
            create_ensembl_gene_data(GENE_ID_5, GENE_NAME_5, CHR_2, NEG_STRAND, GENE_START_5, GENE_START_5 + 1500),
            create_ensembl_gene_data(GENE_ID_6, GENE_NAME_6, CHR_2, NEG_STRAND, GENE_START_6, GENE_START_6 + 1500),
        ],
    );
}

pub fn add_test_transcripts(cache: &mut EnsemblDataCache) {
    let bio = hmf_common::ensembl::BIOTYPE_PROTEIN_CODING;
    let add = |cache: &mut EnsemblDataCache, gene_id: &str, trans_id: i32, strand: i8, start: i32, count: i32| {
        let t = create_trans_exons(
            gene_id,
            trans_id,
            strand,
            &generate_exon_starts(start, count, EXON_LENGTH, 100),
            EXON_LENGTH,
            Some(start + 250),
            Some(start + if gene_id == GENE_ID_1 { 1450 } else { 850 }),
            true,
            bio,
        );
        add_trans_exon_data(cache, gene_id, vec![t]);
    };
    add(cache, GENE_ID_1, TRANS_1, POS_STRAND, GENE_START_1, 7);
    add(cache, GENE_ID_2, TRANS_2, POS_STRAND, GENE_START_2, 5);
    add(cache, GENE_ID_3, TRANS_3, NEG_STRAND, GENE_START_3, 3);
    add(cache, GENE_ID_4, TRANS_4, POS_STRAND, GENE_START_4, 3);
    add(cache, GENE_ID_5, TRANS_5, NEG_STRAND, GENE_START_5, 5);
    add(cache, GENE_ID_6, TRANS_6, NEG_STRAND, GENE_START_6, 6);
}
