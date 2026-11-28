//! fusion: FusionTest, FusionRulesTest, FusionTranscriptTest, SpecialFusionsTest, DisruptionTest and
//! GeneCollectionTest (RnaFusionTest is left out: RNA matching is not ported).
//!
//! Java's tests add genes and known fusion data to the caches the analyser already holds; here each
//! stage takes a snapshot of the caches (FusionEnv::analyser), which the analyser then holds.

use hmf_common::driver_panel::DriverGene;
use hmf_common::ensembl::{BIOTYPE_PROCESSED_TRANS, BIOTYPE_PROTEIN_CODING};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::known_fusion::KnownFusionType::*;
use hmf_common::known_fusion::{KnownFusionCache, KnownFusionData};

use super::*;
use crate::fusion::analyser::FusionDisruptionAnalyser;
use crate::fusion::disruption::DisruptionFinder;
use crate::fusion::finder::FusionFinder;
use crate::fusion::*;
use crate::gene::{find_gene_annotations_by_sv, set_sv_gene_data};

pub struct FusionEnv {
    pub ensembl: EnsemblDataCache,
    pub kfc: KnownFusionCache,
    pub driver_genes: JavaHashMap<String, DriverGene>,
    pub config: FusionConfig,
    pub is_germline: bool,
}

impl FusionEnv {
    pub fn new() -> Self {
        FusionEnv {
            ensembl: create_gene_data_cache(),
            kfc: KnownFusionCache::new(),
            driver_genes: JavaHashMap::new(),
            config: FusionConfig::default(),
            is_germline: false,
        }
    }
    /// initialiseFusions over a snapshot of the caches
    pub fn analyser(&self) -> (&'static GeneCache, FusionDisruptionAnalyser) {
        let gc = leak_cache(self.ensembl.clone());
        let kfc: &'static KnownFusionCache = Box::leak(Box::new(self.kfc.clone()));
        let dg: &'static JavaHashMap<String, DriverGene> = Box::leak(Box::new(self.driver_genes.clone()));
        let fa = FusionDisruptionAnalyser::new(
            FusionFinder::new(self.config.clone(), gc, kfc),
            DisruptionFinder::new(self.is_germline, gc, dg),
            "",
            true,
            self.is_germline,
            false,
            false,
            dg,
        );
        (gc, fa)
    }
}

/// setSvGeneData(svList, cache, applyPreGeneDistance) of SampleDataLoader
pub fn set_sv_gene_data_test(svs: &[Sv], gc: &GeneCache, pre_gene_distance: i32) {
    set_sv_gene_data(svs, gc, pre_gene_distance, &JavaHashMap::new());
}

/// annotate the SVs' genes and transcripts, then run the fusion analyser
pub fn run_fusions(t: &mut LinxTester, gc: &GeneCache, fa: &mut FusionDisruptionAnalyser, pre_gene_distance: i32) {
    set_sv_gene_data_test(&t.all_variants, gc, pre_gene_distance);
    fa.annotate_transcripts(&t.all_variants, true);
    let clusters = t.clusters().clone();
    fa.run(&t.sample_id.clone(), &t.all_variants.clone(), &clusters, &mut t.vis).unwrap();
}

fn validate_fusion_annotations(f: &GeneFusion, valid_ends: bool, valid_traversal: bool) -> bool {
    let Some(a) = f.annotations() else { return false };
    let Some(ci) = &a.chain_info else { return false };
    if valid_traversal != ci.valid_traversal {
        return false;
    }
    if valid_ends {
        !f.is_terminated()
    } else {
        a.terminated_up || a.terminated_down
    }
}

fn exon(trans_id: i32, s: i32, e: i32, rank: i32, ps: i32, pe: i32) -> ExonData {
    ExonData { trans_id, start: s, end: e, rank, phase_start: ps, phase_end: pe }
}

#[allow(clippy::too_many_arguments)]
fn transcript(
    trans_id: i32, name: &str, gene_id: &str, canonical: bool, strand: i8, s: i32, e: i32, cs: Option<i32>, ce: Option<i32>, bio: &str,
    exons: Vec<ExonData>,
) -> TranscriptData {
    TranscriptData {
        trans_id,
        trans_name: name.into(),
        gene_id: gene_id.into(),
        is_canonical: canonical,
        strand,
        trans_start: s,
        trans_end: e,
        coding_start: cs,
        coding_end: ce,
        bio_type: bio.into(),
        ref_seq_id: None,
        exons,
    }
}

// FusionTest

#[test]
fn test_reportable_fusion_comparison() {
    let mut env = FusionEnv::new();
    let gene_id1 = "ENSG0001";
    let mut genes = vec![create_ensembl_gene_data(gene_id1, "GENE1", CHR_1, POS_STRAND, 100, 1000)];
    let t1 = create_trans_exons(gene_id1, 1, POS_STRAND, &[100, 300, 500, 700, 900], 99, Some(300), Some(998), true, BIOTYPE_PROTEIN_CODING);
    add_trans_exon_data(&mut env.ensembl, gene_id1, vec![t1]);
    let gene_id2 = "ENSG0002";
    genes.push(create_ensembl_gene_data(gene_id2, "GENE2", CHR_1, 1, 10000, 12000));
    add_gene_data(&mut env.ensembl, CHR_1, genes);
    // the Java test gives the second transcript the first gene's id
    let t2 = create_trans_exons(gene_id1, 2, POS_STRAND, &[10100, 10300, 10500, 10700, 10900], 99, Some(10100), Some(10998), true, "");
    add_trans_exon_data(&mut env.ensembl, gene_id2, vec![t2]);
    let pre_gene = 100;
    let (gc, mut fa) = env.analyser();

    let mut up = Vec::new();
    for (i, pos) in [250, 450, 650, 850].iter().enumerate() {
        up.extend(find_gene_annotations_by_sv(gc, i as i32, true, CHR_1, *pos, 1, pre_gene));
    }
    for (g, pos) in up.iter().zip([250, 450, 650, 850]) {
        g.set_positional_data(CHR_1, pos, 1);
    }
    let mut down = Vec::new();
    for (i, pos) in [10250, 10450, 10650, 10850].iter().enumerate() {
        down.extend(find_gene_annotations_by_sv(gc, i as i32, false, CHR_1, *pos, -1, pre_gene));
    }
    for (g, pos) in down.iter().zip([10250, 10450, 10650, 10850]) {
        g.set_positional_data(CHR_1, pos, -1);
    }
    let fusions = fa.finder.find_fusions(&up, &down);
    for f in &fusions {
        f.set_known_type(KNOWN_PAIR);
    }
    fa.finder.find_top_reportable_fusion(&fusions);
    assert_eq!(12, fusions.len());
    let reported = fusions.iter().find(|x| x.reportable()).unwrap().clone();
    assert_eq!(450, reported.upstream_trans().gene().position());
    assert_eq!(10850, reported.downstream_trans().gene().position());
    assert_eq!(0, reported.exons_skipped(true));
    assert_eq!(0, reported.exons_skipped(false));
    for f in &fusions {
        if *f != reported {
            assert!(!f.reportable());
        }
    }
}

#[test]
fn test_fusion_prioritisation() {
    let gene_id1 = "ENSG0001";
    let up_gene = create_gene_annotation(0, true, gene_id1, gene_id1, POS_STRAND, CHR_1, 450, 1);
    let up1 = create_transcript(&up_gene, 1, true, 100, 1000, Some(350), Some(950), BIOTYPE_PROTEIN_CODING, 2, 3, PHASE_1, 50, 300);
    let down_gene = create_gene_annotation(0, true, gene_id1, gene_id1, POS_STRAND, CHR_1, 450, 1);
    let (cs, ce) = (Some(10350), Some(10950));
    let down1 = create_transcript(&down_gene, 11, true, 10000, 11000, cs, ce, BIOTYPE_PROTEIN_CODING, 2, 3, PHASE_1, 50, 300);

    let fusion1 = GeneFusion::create(&up1, &down1, true);
    let fusion2 = GeneFusion::create(&up1, &down1, false);
    assert!(find_top_priority_fusion(&[fusion1.clone(), fusion2.clone()]) == Some(fusion1.clone()));

    let fusion2 = GeneFusion::create(&up1, &down1, true);
    fusion1.set_annotations(FusionAnnotations {
        cluster_id: 1,
        cluster_count: 10,
        resolved_type: String::new(),
        chain_info: None,
        terminated_up: true,
        terminated_down: false,
    });
    assert!(find_top_priority_fusion(&[fusion1.clone(), fusion2.clone()]) == Some(fusion2.clone()));

    let down2 = create_transcript(&down_gene, 12, true, 10000, 11000, cs, ce, BIOTYPE_PROCESSED_TRANS, 2, 3, PHASE_1, 50, 300);
    // setAnnotations(null): a fresh fusion with the same fields
    let fusion1 = GeneFusion::create(&up1, &down1, true);
    let fusion2 = GeneFusion::create(&up1, &down2, true);
    assert!(find_top_priority_fusion(&[fusion1.clone(), fusion2.clone()]) == Some(fusion1.clone()));

    let fusion2 = GeneFusion::create(&up1, &down1, true);
    fusion1.set_exons_skipped(2, 0);
    assert!(find_top_priority_fusion(&[fusion1.clone(), fusion2.clone()]) == Some(fusion2.clone()));

    let down3 = create_transcript(&down_gene, 13, false, 10000, 11000, cs, ce, BIOTYPE_PROTEIN_CODING, 2, 3, PHASE_1, 50, 300);
    let fusion1 = GeneFusion::create(&up1, &down3, true);
    assert!(find_top_priority_fusion(&[fusion1.clone(), fusion2.clone()]) == Some(fusion2.clone()));

    let up2 = create_transcript(&up_gene, 2, false, 100, 1000, cs, ce, BIOTYPE_PROTEIN_CODING, 2, 3, PHASE_1, 50, 300);
    let fusion1 = GeneFusion::create(&up2, &down1, true);
    assert!(find_top_priority_fusion(&[fusion1.clone(), fusion2.clone()]) == Some(fusion2.clone()));

    let up3 = create_transcript(&up_gene, 3, true, 100, 1000, cs, ce, BIOTYPE_PROTEIN_CODING, 2, 3, PHASE_1, 40, 200);
    let fusion1 = GeneFusion::create(&up3, &down1, true);
    assert!(find_top_priority_fusion(&[fusion1, fusion2.clone()]) == Some(fusion2));
}

fn chained_fusion_genes(env: &mut FusionEnv, gene_id1: &str, gene_name1: &str) {
    let mut genes = vec![create_ensembl_gene_data(gene_id1, gene_name1, CHR_1, 1, 1000, 2000)];
    let t1 = transcript(
        1,
        "ENST0001",
        gene_id1,
        true,
        POS_STRAND,
        1000,
        2000,
        Some(1401),
        Some(1900),
        BIOTYPE_PROTEIN_CODING,
        vec![exon(1, 1000, 1100, 1, -1, -1), exon(1, 1300, 1500, 2, -1, 1), exon(1, 1600, 1700, 3, 1, 0), exon(1, 1800, 1900, 4, 0, -1)],
    );
    add_trans_exon_data(&mut env.ensembl, gene_id1, vec![t1]);
    genes.push(create_ensembl_gene_data("ENSG0002", "GENE2", CHR_1, POS_STRAND, 10000, 12000));
    add_gene_data(&mut env.ensembl, CHR_1, genes);
    let t2 = transcript(
        2,
        "ENST0002",
        "ENSG0002",
        true,
        POS_STRAND,
        11000,
        12000,
        Some(11049),
        Some(11980),
        BIOTYPE_PROTEIN_CODING,
        vec![exon(1, 11000, 11100, 1, -1, 1), exon(1, 11300, 11501, 2, 1, 2), exon(1, 11600, 11699, 3, 2, 0), exon(1, 11950, 12000, 4, 0, -1)],
    );
    add_trans_exon_data(&mut env.ensembl, "ENSG0002", vec![t2]);
    env.kfc.add_data(KnownFusionData::new(KNOWN_PAIR, gene_name1, "GENE2", "", ""));
}

fn cluster_svs(t: &mut LinxTester, vars: &[&Sv]) -> Cluster {
    t.clear_clusters_and_svs();
    t.add(vars);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(1, t.clusters().len());
    t.clusters()[0].clone()
}

#[test]
fn test_chained_fusions() {
    let mut t = LinxTester::new();
    let mut env = FusionEnv::new();
    chained_fusion_genes(&mut env, "ENSG0001", "GENE1");
    let pre_gene = 200;
    let (gc, mut fa) = env.analyser();

    let var1 = create_del(0, CHR_1, 300, 400);
    let var2 = create_del(1, CHR_1, 500, 600);
    let var3 = create_del(2, CHR_1, 1550, 11200);
    let var4 = create_del(3, CHR_1, 15000, 16000);
    let c = cluster_svs(&mut t, &[&var1, &var2, &var3, &var4]);
    assert_eq!(1, c.chains().len());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(1, fa.fusions.len());
    let f = fa.fusions[0].clone();
    assert_eq!(var3.id(), f.upstream_trans().gene().var_id());
    assert_eq!(var3.id(), f.downstream_trans().gene().var_id());
    assert!(validate_fusion_annotations(&f, true, true));

    let var1 = create_del(0, CHR_1, 1550, 50000);
    let var2 = create_del(1, CHR_1, 50100, 51000);
    let var3 = create_del(2, CHR_1, 55000, 56000);
    let var4 = create_dup(3, CHR_1, 10900, 56500);
    let c = cluster_svs(&mut t, &[&var1, &var2, &var3, &var4]);
    assert_eq!(1, c.chains().len());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(1, fa.fusions.len());
    let f = fa.fusions[0].clone();
    assert_eq!(var1.id(), f.upstream_trans().gene().var_id());
    assert_eq!(var4.id(), f.downstream_trans().gene().var_id());
    assert!(validate_fusion_annotations(&f, true, true));

    let var1 = create_del(0, CHR_1, 1150, 1250);
    let var2 = create_dup(1, CHR_1, 1350, 1575);
    let var3 = create_del(2, CHR_1, 1750, 11525);
    let var4 = create_del(3, CHR_1, 11575, 111800);
    let c = cluster_svs(&mut t, &[&var1, &var2, &var3, &var4]);
    assert_eq!(1, c.chains().len());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(3, fa.fusions.len());
    let f = fa.fusions.iter().find(|x| x.reportable()).unwrap().clone();
    assert_eq!(var3.id(), f.upstream_trans().gene().var_id());
    assert_eq!(var3.id(), f.downstream_trans().gene().var_id());
    assert!(validate_fusion_annotations(&f, false, true));

    // a third gene on a remote chromosome
    let remote = "4";
    add_gene_data(&mut env.ensembl, remote, vec![create_ensembl_gene_data("ENSG0003", "GENE3", remote, 1, 1000, 3000)]);
    let t3 = create_trans_exons("ENSG0003", 3, POS_STRAND, &[1000, 2000, 3000], 98, Some(1003), Some(3098), true, "");
    add_trans_exon_data(&mut env.ensembl, "ENSG0003", vec![t3]);
    let (gc, mut fa) = env.analyser();
    let var1 = create_bnd(1, "1", 51000, 1, "2", 1000, -1);
    let var2 = create_inv(2, CHR_1, 500, 50000, -1);
    let var3 = create_bnd(3, CHR_1, 1750, 1, "2", 10000, -1);
    let var4 = create_bnd(4, CHR_1, 61000, 1, "2", 11000, 1);
    let var5 = create_bnd(5, CHR_1, 60000, -1, remote, 1200, -1);
    let var6 = create_bnd(6, remote, 2500, 1, "5", 1000, -1);
    let var7 = create_bnd(7, CHR_1, 11525, -1, "5", 2000, 1);
    let c = cluster_svs(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6, &var7]);
    assert_eq!(1, c.chains().len());
    assert_eq!(7, c.chains()[0].sv_count());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(2, fa.fusions.len());
    assert_eq!(1, fa.fusions.iter().filter(|x| x.valid_chain_traversal()).count());
    let f = fa.fusions.iter().find(|x| x.valid_chain_traversal()).unwrap().clone();
    assert!(f.phase_matched());
    assert_eq!(var3.id(), f.upstream_trans().gene().var_id());
    assert_eq!(var5.id(), f.downstream_trans().gene().var_id());
    assert!(!validate_fusion_annotations(&f, false, false));

    let clusters = t.clusters().clone();
    fa.run(&t.sample_id.clone(), &t.all_variants.clone(), &clusters, &mut t.vis).unwrap();
    assert_eq!(2, fa.fusions.len());
    let f = fa
        .fusions
        .iter()
        .find(|x| x.upstream_trans().gene().var_id() == var3.id() && x.downstream_trans().gene().var_id() == var7.id())
        .unwrap()
        .clone();
    assert_eq!(KNOWN_PAIR, f.known_type());
    assert!(!f.valid_chain_traversal());
}

#[test]
fn test_chained_fusion_disruptions() {
    let mut t = LinxTester::new();
    let mut env = FusionEnv::new();
    chained_fusion_genes(&mut env, "ENSG0001", "GENE1");
    let pre_gene = 30000;
    let (gc, mut fa) = env.analyser();

    let var1 = create_bnd(1, "1", 1510, 1, "2", 1000, -1);
    let var2 = create_bnd(2, "1", 1580, -1, "2", 1100, 1);
    let c = cluster_svs(&mut t, &[&var1, &var2]);
    assert_eq!(1, c.chains().len());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert!(fa.fusions.is_empty());
    assert_eq!(1, fa.invalid_fusions.len());
    assert!(fa.invalid_fusions[0].0.non_disruptive_chain());

    let var1 = create_del(0, CHR_1, 1150, 1250);
    let var2 = create_dup(1, CHR_1, 1525, 1575);
    let var3 = create_del(2, CHR_1, 1750, 11525);
    let var4 = create_bnd(3, CHR_1, 11750, 1, CHR_2, 50000, -1);
    let var5 = create_bnd(4, CHR_1, 11850, -1, CHR_2, 50100, 1);
    let c = cluster_svs(&mut t, &[&var1, &var2, &var3, &var4, &var5]);
    assert_eq!(1, c.chains().len());
    assert_eq!(4, c.chains()[0].link_count());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(1, fa.unique_fusions.len());
    let f = fa.unique_fusions.iter().find(|x| x.reportable()).unwrap().clone();
    assert_eq!(var3.id(), f.upstream_trans().gene().var_id());
    assert_eq!(var3.id(), f.downstream_trans().gene().var_id());
    assert!(validate_fusion_annotations(&f, true, true));
}

// FusionRulesTest

use crate::fusion::finder::check_fusion_logic;
use crate::gene::{TranscriptRegionType, POST_CODING_PHASE};
use hmf_common::known_fusion::{FS_DOWN, FS_UP};

fn coding_bases(s: Option<i32>, e: Option<i32>) -> i32 {
    match (s, e) {
        (Some(s), Some(e)) => e - s + 1,
        _ => 0,
    }
}
fn rules_config() -> FusionConfig {
    FusionConfig::default()
}

#[test]
fn test_fusion_types() {
    let gene1 = create_gene_annotation(0, true, "GENE1", "ENSG0001", POS_STRAND, CHR_1, 150, 1);
    let gene2 = create_gene_annotation(0, false, "GENE2", "ENSG0003", NEG_STRAND, CHR_1, 150, 1);
    let params = rules_config();
    let (none, some) = (None, None::<i32>);
    let tr = |g: &BreakendGene, id: i32, canon: bool, cs: Option<i32>, ce: Option<i32>, eu: i32, ed: i32, phase: i32, cb: i32| {
        create_transcript(g, id, canon, 50, 250, cs, ce, "", eu, ed, phase, cb, coding_bases(cs, ce))
    };
    // the upstream transcripts of the Java test are protein coding since upstream 7a4490caa6
    let tr_pc = |g: &BreakendGene, id: i32, canon: bool, cs: Option<i32>, ce: Option<i32>, eu: i32, ed: i32, phase: i32, cb: i32| {
        create_transcript(g, id, canon, 50, 250, cs, ce, BIOTYPE_PROTEIN_CODING, eu, ed, phase, cb, coding_bases(cs, ce))
    };
    let trans1 = tr_pc(&gene1, 1, false, none, some, 2, 3, PHASE_NONE, 0);
    let (cs, ce) = (Some(100), Some(200));
    let trans2 = tr(&gene2, 2, false, cs, ce, 2, 3, PHASE_NONE, 10);
    assert!(trans1.non_coding());
    assert!(trans2.is_coding());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());

    let trans2 = tr(&gene2, 2, false, None, None, 2, 3, PHASE_NONE, 10);
    let trans1 = tr(&gene2, 2, false, cs, ce, 2, 3, PHASE_NONE, 10);
    assert!(trans1.is_coding());
    assert!(trans2.non_coding());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());

    let trans2 = tr(&gene2, 2, false, cs, ce, 2, 3, PHASE_0, 10);
    gene1.set_positional_data(CHR_1, 260, 1);
    let trans1 = tr(&gene1, 1, true, cs, ce, 2, 3, PHASE_NONE, coding_bases(cs, ce));
    assert!(trans1.post_coding());
    assert!(trans2.is_coding());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());

    let trans1 = tr(&gene1, 1, true, cs, ce, 0, 1, PHASE_NONE, 0);
    assert!(trans1.is_promoter());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());

    gene1.set_positional_data(CHR_1, 40, 1);
    let trans1 = tr_pc(&gene1, 1, true, cs, ce, 2, 3, PHASE_NONE, 0);
    let trans2 = tr(&gene2, 2, true, cs, ce, 2, 3, PHASE_NONE, 10);
    assert!(trans1.pre_coding());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());

    trans1.set_is_disruptive(false);
    let trans2 = tr(&gene2, 2, true, cs, ce, 2, 3, PHASE_NONE, 10);
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());
    trans1.set_is_disruptive(true);

    // the Java test calls this "up coding exonic to not down coding exonic", but at position 40 the transcript is
    // pre-coding: it is the pre-coding exonic rule that rejects it, before and after upstream 7a4490caa6
    let trans1 = tr_pc(&gene1, 1, true, cs, ce, 2, 2, PHASE_NONE, 0);
    let trans2 = tr(&gene2, 2, false, cs, ce, 2, 3, PHASE_NONE, 0);
    assert!(trans1.is_exonic());
    assert!(trans1.pre_coding());
    assert!(trans2.is_intronic());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());

    // "unmatched phasing intronic" left the Java test with RequirePhaseMatch (upstream 7a4490caa6)

    gene1.set_positional_data(CHR_1, 40, 1);
    let trans1 = create_transcript(&gene1, 1, true, 50, 450, cs, ce, BIOTYPE_PROTEIN_CODING, 1, 2, PHASE_NONE, 0, coding_bases(cs, ce));
    let trans2 = create_transcript(&gene1, 1, true, 50, 550, cs, ce, "", 4, 5, PHASE_NONE, 0, coding_bases(cs, ce));
    assert!(trans1.pre_coding() && trans2.pre_coding());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_none());

    gene1.set_positional_data(CHR_1, 110, 1);
    let trans1 =
        create_transcript_ex(&gene1, 1, true, 50, 250, cs, ce, BIOTYPE_PROTEIN_CODING, 2, 2, PHASE_1, PHASE_2, 10, coding_bases(cs, ce));
    gene2.set_positional_data(CHR_1, 190, 1);
    let trans2 = create_transcript_ex(&gene2, 2, true, 50, 250, cs, ce, "", 2, 2, PHASE_0, PHASE_2, 10, coding_bases(cs, ce));
    assert!(trans1.is_exonic() && trans2.is_exonic());
    // unphased fusions are kept since upstream 7a4490caa6
    let f = check_fusion_logic(&trans1, &trans2, &params).unwrap();
    assert!(!f.phase_matched());

    let trans1 = tr_pc(&gene1, 1, true, cs, ce, 2, 2, PHASE_2, 11);
    let f = check_fusion_logic(&trans1, &trans2, &params).unwrap();
    assert!(!f.phase_matched());
    let trans2 = tr(&gene2, 2, true, cs, ce, 2, 2, PHASE_0, 12);
    let f = check_fusion_logic(&trans1, &trans2, &params).unwrap();
    assert!(f.phase_matched());

    gene1.set_positional_data(CHR_1, 90, 1);
    let trans1 = tr_pc(&gene1, 1, true, cs, ce, 2, 2, PHASE_NONE, 0);
    gene2.set_positional_data(CHR_1, 210, 1);
    let trans2 = tr(&gene2, 2, true, cs, ce, 2, 2, PHASE_NONE, 0);
    assert!(trans1.is_exonic() && trans2.is_exonic());
    assert!(trans1.pre_coding() && trans2.pre_coding());
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_some());

    let trans1 = tr_pc(&gene1, 1, true, cs, ce, 2, 3, PHASE_1, 10);
    let trans2 = tr(&gene2, 2, true, cs, ce, 2, 3, PHASE_1, 10);
    assert!(check_fusion_logic(&trans1, &trans2, &params).is_some());

    let mut kfc = KnownFusionCache::new();
    kfc.add_data(KnownFusionData::new(KNOWN_PAIR, gene1.gene_name(), gene2.gene_name(), "", ""));
    let kfc: &'static KnownFusionCache = Box::leak(Box::new(kfc));
    let finder = FusionFinder::new(rules_config(), leak_cache(create_gene_data_cache()), kfc);
    gene1.add_transcript(trans1);
    gene2.add_transcript(trans2);
    let fusions = finder.find_fusions(&[gene1], &[gene2]);
    assert_eq!(1, fusions.len());
    assert!(fusions[0].phase_matched());
}

#[test]
fn test_alternate_phasings() {
    let gene1 = create_gene_annotation(0, true, "GENE1", "ENSG0001", POS_STRAND, CHR_1, 150, 1);
    let gene2 = create_gene_annotation(0, false, "GENE2", "ENSG0003", NEG_STRAND, CHR_1, 150, 1);
    let (cs, ce) = (Some(100), Some(200));
    let cb = coding_bases(cs, ce);
    let up = create_transcript(&gene1, 1, true, 50, 250, cs, ce, BIOTYPE_PROTEIN_CODING, 2, 3, PHASE_1, 10, cb);
    let down = create_transcript(&gene2, 2, false, 50, 250, cs, ce, "", 2, 3, PHASE_0, 10, cb);
    let params = rules_config();
    assert!(up.is_coding() && down.is_coding());
    let f = check_fusion_logic(&up, &down, &params).unwrap();
    assert!(!f.phase_matched());

    up.alternative_phasing_mut().insert(0, 1);
    let f = check_fusion_logic(&up, &down, &params).unwrap();
    assert!(f.phase_matched());
    assert_eq!(1, f.exons_skipped(true));
    assert_eq!(0, f.exons_skipped(false));

    up.alternative_phasing_mut().clear();
    down.alternative_phasing_mut().insert(1, 1);
    let f = check_fusion_logic(&up, &down, &params).unwrap();
    assert!(f.phase_matched());
    assert_eq!(0, f.exons_skipped(true));
    assert_eq!(1, f.exons_skipped(false));

    gene1.set_positional_data(CHR_1, 210, 1);
    let up = create_transcript(&gene1, 1, true, 50, 250, cs, ce, BIOTYPE_PROTEIN_CODING, 6, 7, PHASE_NONE, 100, 100);
    assert!(up.post_coding());
    assert_eq!(POST_CODING_PHASE, up.phase);
    up.alternative_phasing_mut().clear();
    up.alternative_phasing_mut().insert(0, 3);
    let f = check_fusion_logic(&up, &down, &params).unwrap();
    assert!(f.phase_matched());
    assert_eq!(3, f.exons_skipped(true));
    assert_eq!(0, f.exons_skipped(false));

    let up = create_transcript(&gene1, 1, true, 50, 250, cs, ce, BIOTYPE_PROTEIN_CODING, 2, 3, PHASE_1, 10, cb);
    let down = create_transcript(&gene2, 2, true, 50, 250, cs, ce, "", 2, 3, PHASE_0, 10, cb);
    let f = check_fusion_logic(&up, &down, &params).unwrap();
    assert!(!f.phase_matched());
    up.alternative_phasing_mut().insert(2, 2);
    down.alternative_phasing_mut().insert(2, 3);
    let f = check_fusion_logic(&up, &down, &params).unwrap();
    assert!(f.phase_matched());
    assert_eq!(2, f.exons_skipped(true));
    assert_eq!(3, f.exons_skipped(false));
}

#[test]
fn test_exonic_fusions() {
    let mut env = FusionEnv::new();
    let mut genes = vec![create_ensembl_gene_data("ENSG0001", "GENE1", CHR_1, POS_STRAND, 100, 1000)];
    let t1 = create_trans_exons("ENSG0001", 1, POS_STRAND, &[100, 300, 500, 700, 900], 100, Some(150), Some(950), true, BIOTYPE_PROTEIN_CODING);
    add_trans_exon_data(&mut env.ensembl, "ENSG0001", vec![t1]);
    genes.push(create_ensembl_gene_data("ENSG0002", "GENE2", CHR_1, POS_STRAND, 10100, 11000));
    add_gene_data(&mut env.ensembl, CHR_1, genes);
    let t2 =
        create_trans_exons("ENSG0002", 2, POS_STRAND, &[10100, 10300, 10500, 10700, 10900], 100, Some(10150), Some(10950), true, BIOTYPE_PROTEIN_CODING);
    add_trans_exon_data(&mut env.ensembl, "ENSG0002", vec![t2]);
    let pre_gene = 100000;
    let (gc, fa) = env.analyser();

    let up = find_gene_annotations_by_sv(gc, 0, true, CHR_1, 701, 1, pre_gene);
    up[0].set_positional_data(CHR_1, 701, 1);
    let down = find_gene_annotations_by_sv(gc, 0, false, CHR_1, 10650, -1, pre_gene);
    down[0].set_positional_data(CHR_1, 10650, -1);
    let fusions = fa.finder.find_fusions(&up, &down);
    assert_eq!(1, fusions.len());
    let f = &fusions[0];
    assert!(f.phase_matched());
    assert_eq!(0, f.transcripts()[FS_UP].exonic_base_phase);
    assert_eq!(1, f.transcripts()[FS_UP].phase);
    assert_eq!(TranscriptRegionType::EXONIC, f.transcripts()[FS_UP].region_type());
    assert_eq!(4, f.transcripts()[FS_UP].exon_upstream);
    assert_eq!(3, f.transcripts()[FS_DOWN].exon_upstream);
    assert_eq!(3, f.fused_exon(true));
    assert_eq!(0, f.exons_skipped(true));
    assert_eq!(0, f.exons_skipped(false));

    let up = find_gene_annotations_by_sv(gc, 0, true, CHR_1, 701, 1, pre_gene);
    up[0].set_positional_data(CHR_1, 701, 1);
    let down = find_gene_annotations_by_sv(gc, 0, false, CHR_1, 10705, -1, pre_gene);
    down[0].set_positional_data(CHR_1, 10705, -1);
    let fusions = fa.finder.find_fusions(&up, &down);
    assert_eq!(1, fusions.len());
    let f = &fusions[0];
    assert!(f.phase_matched());
    assert_eq!(0, f.transcripts()[FS_UP].exonic_base_phase);
    assert_eq!(1, f.transcripts()[FS_DOWN].exonic_base_phase);
    assert_eq!(TranscriptRegionType::EXONIC, f.transcripts()[FS_UP].region_type());
    assert_eq!(TranscriptRegionType::EXONIC, f.transcripts()[FS_DOWN].region_type());
    assert_eq!(4, f.transcripts()[FS_UP].exon_upstream);
    assert_eq!(4, f.transcripts()[FS_DOWN].exon_upstream);
    assert_eq!(4, f.fused_exon(true));
    assert_eq!(4, f.fused_exon(false));
}

// FusionTranscriptTest

use crate::gene::{create_breakend_transcript_data, set_alternative_transcript_phasings, TranscriptCodingType};

fn leak_trans(t: TranscriptData) -> &'static TranscriptData {
    Box::leak(Box::new(t))
}
fn alt(t: &BreakendTrans, phase: i32) -> Option<i32> {
    t.alternative_phasing().get(&phase).copied()
}

#[test]
fn test_alternative_breakend_transcript_phasing() {
    let gene_pos = create_gene_annotation(0, true, "GENE1", "ENSG0001", POS_STRAND, CHR_1, 0, 1);
    let exon_starts = [100, 300, 500, 700, 900];
    let td = leak_trans(create_trans_exons("ENSG0001", 1, POS_STRAND, &exon_starts, 100, Some(349), Some(950), false, ""));

    gene_pos.set_positional_data(CHR_1, 250, 1);
    let t = create_breakend_transcript_data(td, 250, &gene_pos).unwrap();
    assert_eq!(5, t.exon_count());
    assert_eq!(1, t.exon_upstream);
    assert_eq!(2, t.exon_downstream);
    assert_eq!(PHASE_NONE, t.phase);
    set_alternative_transcript_phasings(&t, &td.exons, 250, 1);
    assert_eq!(0, t.alternative_phasing().len());
    set_alternative_transcript_phasings(&t, &td.exons, 250, -1);
    assert_eq!(3, t.alternative_phasing().len());
    assert_eq!(Some(1), alt(&t, PHASE_1));
    assert_eq!(Some(3), alt(&t, PHASE_2));
    assert_eq!(Some(2), alt(&t, PHASE_0));

    gene_pos.set_positional_data(CHR_1, 450, 1);
    let t = create_breakend_transcript_data(td, 450, &gene_pos).unwrap();
    assert_eq!(2, t.exon_upstream);
    assert_eq!(3, t.exon_downstream);
    assert_eq!(PHASE_1, t.phase);
    set_alternative_transcript_phasings(&t, &td.exons, 450, 1);
    assert_eq!(1, t.alternative_phasing().len());
    assert_eq!(Some(1), alt(&t, -1));
    set_alternative_transcript_phasings(&t, &td.exons, 450, -1);
    assert_eq!(2, t.alternative_phasing().len());
    assert_eq!(Some(2), alt(&t, PHASE_2));
    assert_eq!(Some(1), alt(&t, PHASE_0));

    gene_pos.set_positional_data(CHR_1, 650, 1);
    let t = create_breakend_transcript_data(td, 650, &gene_pos).unwrap();
    assert_eq!(3, t.exon_upstream);
    assert_eq!(4, t.exon_downstream);
    assert_eq!(PHASE_0, t.phase);
    set_alternative_transcript_phasings(&t, &td.exons, 650, 1);
    assert_eq!(2, t.alternative_phasing().len());
    assert_eq!(Some(2), alt(&t, -1));
    assert_eq!(Some(1), alt(&t, 1));
    set_alternative_transcript_phasings(&t, &td.exons, 650, -1);
    assert_eq!(1, t.alternative_phasing().len());
    assert_eq!(Some(1), alt(&t, PHASE_2));

    let gene_neg = create_gene_annotation(0, true, "GENE2", "ENSG0002", NEG_STRAND, CHR_1, 0, 1);
    let td = leak_trans(create_trans_exons("ENSG0002", 2, NEG_STRAND, &exon_starts, 100, Some(349), Some(950), false, ""));
    gene_neg.set_positional_data(CHR_1, 850, 1);
    let t = create_breakend_transcript_data(td, 850, &gene_neg).unwrap();
    assert_eq!(5, t.exon_count());
    assert_eq!(1, t.exon_upstream);
    assert_eq!(2, t.exon_downstream);
    assert_eq!(PHASE_0, t.phase);
    set_alternative_transcript_phasings(&t, &td.exons, 850, -1);
    assert_eq!(0, t.alternative_phasing().len());
    set_alternative_transcript_phasings(&t, &td.exons, 850, 1);
    assert_eq!(3, t.alternative_phasing().len());
    assert_eq!(Some(2), alt(&t, PHASE_1));
    assert_eq!(Some(3), alt(&t, POST_CODING_PHASE));
    assert_eq!(Some(1), alt(&t, PHASE_2));

    gene_neg.set_positional_data(CHR_1, 250, 1);
    let t = create_breakend_transcript_data(td, 250, &gene_neg).unwrap();
    assert_eq!(4, t.exon_upstream);
    assert_eq!(5, t.exon_downstream);
    assert_eq!(POST_CODING_PHASE, t.phase);
    set_alternative_transcript_phasings(&t, &td.exons, 250, -1);
    assert_eq!(3, t.alternative_phasing().len());
    assert_eq!(Some(1), alt(&t, PHASE_1));
    assert_eq!(Some(3), alt(&t, PHASE_0));
    assert_eq!(Some(2), alt(&t, PHASE_2));
    set_alternative_transcript_phasings(&t, &td.exons, 250, 1);
    assert_eq!(0, t.alternative_phasing().len());
}

#[test]
fn test_breakend_transcript_coding() {
    use TranscriptCodingType::*;
    let gene_pos = create_gene_annotation(0, true, "GENE_1", "GENE001", POS_STRAND, CHR_1, 0, 1);
    let gene_neg = create_gene_annotation(0, true, "GENE_2", "GENE002", NEG_STRAND, CHR_1, 0, 1);
    let starts = [100, 200, 300, 400, 500];
    let mut id = 1;
    let mut td = |gene: &str, strand: i8, cs: i32, ce: i32| {
        let t = leak_trans(create_trans_exons(gene, id, strand, &starts, 10, Some(cs), Some(ce), true, BIOTYPE_PROTEIN_CODING));
        id += 1;
        t
    };
    let check = |g: &BreakendGene, td: &'static TranscriptData, pos: i32, o: i8, up: i32, down: i32, ct: TranscriptCodingType, phase: i32| {
        g.set_positional_data(CHR_1, pos, o);
        let t = create_breakend_transcript_data(td, pos, g).unwrap();
        assert_eq!((up, down, ct, phase), (t.exon_upstream, t.exon_downstream, t.coding_type(), t.phase));
        t
    };
    let t = check(&gene_pos, td("GENE001", POS_STRAND, 100, 110), 150, 1, 1, 2, UTR_3P, -2);
    assert_eq!(5, t.exon_count());
    assert_eq!((8, 8), (t.coding_bases, t.total_coding_bases));
    let t = check(&gene_pos, td("GENE001", POS_STRAND, 105, 405), 350, 1, 3, 4, CODING, PHASE_1);
    assert_eq!((28, 31), (t.coding_bases, t.total_coding_bases));
    let t = check(&gene_neg, td("GENE002", NEG_STRAND, 500, 510), 450, 1, 1, 2, UTR_3P, -2);
    assert_eq!(5, t.exon_count());
    assert_eq!((8, 8), (t.coding_bases, t.total_coding_bases));
    let t = check(&gene_neg, td("GENE002", NEG_STRAND, 205, 505), 250, 1, 3, 4, CODING, PHASE_1);
    assert_eq!((28, 31), (t.coding_bases, t.total_coding_bases));
    check(&gene_pos, td("GENE001", POS_STRAND, 200, 405), 50, -1, 0, 2, UTR_5P, PHASE_NONE);
    check(&gene_pos, td("GENE001", POS_STRAND, 105, 405), 50, -1, 0, 2, CODING, PHASE_0);
    check(&gene_neg, td("GENE002", NEG_STRAND, 105, 310), 350, 1, 2, 3, UTR_5P, PHASE_NONE);
    check(&gene_neg, td("GENE002", NEG_STRAND, 105, 405), 350, 1, 2, 3, CODING, PHASE_0);
}

// GeneCollectionTest

use hmf_common::ensembl::{get_exon_rankings, get_protein_domain_positions, TranscriptProteinData};

#[test]
fn test_exon_data_extraction() {
    let mut cache = create_gene_data_cache();
    add_gene_data(&mut cache, "1", vec![create_ensembl_gene_data("ENSG0001", "GENE1", "1", 1, 10000, 20000)]);
    let t1 = create_trans_exons("ENSG0001", 1, 1, &[10500, 11500, 12500, 13500], 99, Some(11500), Some(13598), true, "");
    let t2 = create_trans_exons("ENSG0001", 2, 1, &[12500, 13500, 14500], 100, None, None, true, "");
    let name2 = t2.trans_name.clone();
    add_trans_exon_data(&mut cache, "ENSG0001", vec![t1, t2]);
    let t = cache.transcript_data("ENSG0001", "").unwrap();
    assert_eq!(1, t.trans_id);
    assert_eq!(4, t.exons.len());
    let t = cache.transcript_data("ENSG0001", &name2).unwrap();
    assert_eq!(2, t.trans_id);
    assert_eq!(3, t.exons.len());
    // getExonRankings(geneId, position): the canonical transcript's
    let canonical = cache.transcript_data("ENSG0001", "").unwrap();
    let r = |p: i32| get_exon_rankings(canonical.strand, &canonical.exons, p);
    assert_eq!([1, 2], r(11400));
    assert_eq!([0, 1], r(9000));
    assert_eq!([4, -1], r(16000));
    assert_eq!([3, 3], r(12500));
}

#[test]
fn test_protein_domain_positions() {
    let pd = |s: i32, e: i32| TranscriptProteinData {
        transcript_id: 2,
        translation_id: 0,
        protein_feature_id: 0,
        seq_start: s,
        seq_end: e,
        hit_description: "hd".into(),
    };
    let t = create_trans_exons("G0001", 1, POS_STRAND, &[100, 300, 500], 100, Some(150), Some(550), true, "");
    assert_eq!([Some(165), Some(515)], get_protein_domain_positions(&pd(5, 55), &t));
    assert_eq!([Some(515), Some(545)], get_protein_domain_positions(&pd(55, 65), &t));
    let t = create_trans_exons("G0001", 2, NEG_STRAND, &[100, 300, 500], 100, Some(350), Some(550), true, "");
    assert_eq!([Some(185), Some(535)], get_protein_domain_positions(&pd(5, 55), &t));
    assert_eq!([Some(155), Some(185)], get_protein_domain_positions(&pd(55, 65), &t));
}

// SpecialFusionsTest

fn test_genes_env() -> FusionEnv {
    let mut env = FusionEnv::new();
    add_test_genes(&mut env.ensembl);
    add_test_transcripts(&mut env.ensembl);
    env
}
fn breakend_genes(gc: &GeneCache, sv_id: i32, is_start: bool, chr: &str, pos: i32, o: i8, pre_gene: i32) -> Vec<BreakendGene> {
    let genes = find_gene_annotations_by_sv(gc, sv_id, is_start, chr, pos, o, pre_gene);
    genes[0].set_positional_data(chr, pos, o);
    genes
}
const DEFAULT_PRE_GENE: i32 = 100000;

#[test]
fn test_same_gene_fusions() {
    let mut env = test_genes_env();
    env.kfc.add_data(KnownFusionData::new(PROMISCUOUS_3, "", GENE_NAME_1, "", ""));
    let (gc, mut fa) = env.analyser();
    let up = breakend_genes(gc, 0, false, CHR_1, 1950, 1, DEFAULT_PRE_GENE);
    let down = breakend_genes(gc, 0, true, CHR_1, 1350, -1, DEFAULT_PRE_GENE);
    let fusions = fa.finder.find_fusions(&up, &down);
    assert_eq!(1, fusions.len());
    let f = &fusions[0];
    assert_eq!(1950, f.upstream_trans().gene().position());
    assert_eq!(1350, f.downstream_trans().gene().position());
    assert_eq!(0, f.exons_skipped(false));
    assert!(!f.reportable());
}

#[test]
fn test_exon_del_dup_fusion() {
    let mut env = test_genes_env();
    let mut kf = KnownFusionData::new(EXON_DEL_DUP, GENE_NAME_1, GENE_NAME_1, "", "");
    kf.set_known_exon_data(&generate_trans_name(TRANS_1), "2;3", "5;6").unwrap();
    env.kfc.add_data(kf);
    // a second gene for the second stage (on its own, it touches nothing else)
    let gene_id = "ENSG0010";
    let gene_start = 30000;
    let mut chr1 = env.ensembl.chr_gene_data_map().get(&CHR_1.to_string()).unwrap().clone();
    chr1.push(create_ensembl_gene_data(gene_id, gene_id, CHR_1, POS_STRAND, gene_start, gene_start + 15 * 200));
    add_gene_data(&mut env.ensembl, CHR_1, chr1);
    let t = create_trans_exons(
        gene_id,
        10,
        POS_STRAND,
        &generate_exon_starts(30000, 15, 100, 100),
        100,
        Some(gene_start + 250),
        Some(gene_start + 2850),
        true,
        BIOTYPE_PROTEIN_CODING,
    );
    add_trans_exon_data(&mut env.ensembl, gene_id, vec![t]);
    let mut kf = KnownFusionData::new(EXON_DEL_DUP, gene_id, gene_id, "", "");
    kf.set_known_exon_data(&generate_trans_name(10), "2;5", "10;14").unwrap();
    env.kfc.add_data(kf);
    let (gc, mut fa) = env.analyser();

    let up = breakend_genes(gc, 0, true, CHR_1, 1550, 1, DEFAULT_PRE_GENE);
    let down = breakend_genes(gc, 0, false, CHR_1, 2150, -1, DEFAULT_PRE_GENE);
    let mut fusions = fa.finder.find_fusions(&up, &down);
    let up = breakend_genes(gc, 1, true, CHR_1, 1350, 1, DEFAULT_PRE_GENE);
    let down = breakend_genes(gc, 1, false, CHR_1, 1950, -1, DEFAULT_PRE_GENE);
    fusions.extend(fa.finder.find_fusions(&up, &down));
    assert_eq!(1, fusions.len());
    let f = fa.finder.find_top_reportable_fusion(&fusions).unwrap();
    assert_eq!(EXON_DEL_DUP, f.known_type());
    assert_eq!(1350, f.upstream_trans().gene().position());
    assert_eq!(1950, f.downstream_trans().gene().position());
    assert_eq!(0, f.exons_skipped(true));
    assert_eq!(0, f.exons_skipped(false));
    assert!(f.reportable());

    let up = breakend_genes(gc, 1, true, CHR_1, gene_start + 350, 1, DEFAULT_PRE_GENE);
    let down = breakend_genes(gc, 1, false, CHR_1, gene_start + 1950, -1, DEFAULT_PRE_GENE);
    let fusions = fa.finder.find_fusions(&up, &down);
    assert_eq!(1, fusions.len());
    let f = fa.finder.find_top_reportable_fusion(&fusions).unwrap();
    assert_eq!(EXON_DEL_DUP, f.known_type());
    assert_eq!(gene_start + 350, f.upstream_trans().gene().position());
    assert_eq!(gene_start + 1950, f.downstream_trans().gene().position());
    assert_eq!(0, f.exons_skipped(true));
    assert_eq!(1, f.exons_skipped(false));
    assert_eq!(2, f.fused_exon(true));
    assert_eq!(12, f.fused_exon(false));
    assert!(f.reportable());
}

#[test]
fn test_enhancer_region_fusion() {
    let mut env = FusionEnv::new();
    let chr = CHR_1;
    add_gene_data(
        &mut env.ensembl,
        chr,
        vec![
            create_ensembl_gene_data("ENSG0001", "IGH", chr, POS_STRAND, 100, 1500),
            create_ensembl_gene_data("ENSG0002", "GENE2", chr, POS_STRAND, 10000, 11500),
            create_ensembl_gene_data("ENSG0003", "GENE3", chr, POS_STRAND, 10020000, 10021500),
        ],
    );
    let starts = generate_exon_starts(100, 7, 98, 100);
    add_trans_exon_data(&mut env.ensembl, "ENSG0001", vec![create_trans_exons("ENSG0001", 1, POS_STRAND, &starts, 98, Some(300), Some(1398), true, "")]);
    let starts = generate_exon_starts(10000, 7, 98, 100);
    add_trans_exon_data(&mut env.ensembl, "ENSG0002", vec![create_trans_exons("ENSG0002", 2, POS_STRAND, &starts, 98, Some(10300), Some(11390), true, "")]);
    let starts = generate_exon_starts(10020000, 7, 98, 100);
    add_trans_exon_data(
        &mut env.ensembl,
        "ENSG0003",
        vec![create_trans_exons("ENSG0003", 3, POS_STRAND, &starts, 98, Some(10020300), Some(10021390), true, "")],
    );
    let ig_region = format!("ENHANCER_RANGE={};{};{};{}", NEG_STRAND, chr, 50, 2000);
    let mut kf = KnownFusionData::new(ENHANCER_KNOWN_PAIR, "IGH", "GENE2", "", "");
    kf.apply_overrides(&ig_region).unwrap();
    env.kfc.add_data(kf);
    let mut kf = KnownFusionData::new(ENHANCER_PROMISCUOUS, "IGH", "", "", "");
    kf.apply_overrides(&ig_region).unwrap();
    env.kfc.add_data(kf);
    let (gc, mut fa) = env.analyser();

    let up = breakend_genes(gc, 0, true, chr, 200, -1, 1000);
    let down = breakend_genes(gc, 0, false, chr, 9500, -1, 1000);
    let mut fusions = fa.finder.find_fusions(&up, &down);
    let up = find_gene_annotations_by_sv(gc, 1, true, chr, 200, -1, 1000);
    up[0].set_positional_data(chr, 500, -1);
    let down = find_gene_annotations_by_sv(gc, 1, false, chr, 10019500, -1, 1000);
    down[0].set_positional_data(chr, 10020100, -1);
    fusions.extend(fa.finder.find_fusions(&up, &down));
    assert_eq!(2, fusions.len());
    let f = fa.finder.find_top_reportable_fusion(&fusions).unwrap();
    assert_eq!(ENHANCER_KNOWN_PAIR, f.known_type());
    assert_eq!(200, f.upstream_trans().gene().position());
    assert_eq!(9500, f.downstream_trans().gene().position());
    assert!(f.reportable());
    let f = fusions.iter().find(|x| x.known_type() == ENHANCER_PROMISCUOUS).unwrap();
    assert_eq!(500, f.upstream_trans().gene().position());
    assert_eq!(10020100, f.downstream_trans().gene().position());
    assert!(!f.reportable());
}

fn sgl_mapping(chr: &str, pos: i32, o: i8) -> SglMapping {
    SglMapping { chromosome: chr.into(), position: pos, orientation: o, cigar: String::new(), qual_score: 1 }
}

#[test]
fn test_single_breakend_fusions() {
    let mut t = LinxTester::new();
    let mut env = test_genes_env();
    env.kfc.add_data(KnownFusionData::new(KNOWN_PAIR, GENE_NAME_1, GENE_NAME_2, "", ""));
    let (gc, mut fa) = env.analyser();
    let pre_gene = 200;
    let mut var_id = 1;
    let mut next = || {
        var_id += 1;
        var_id - 1
    };

    let sgl1 = create_sgl(next(), CHR_1, 1150, 1);
    sgl1.add_sgl_mapping(SglMapping::from(&format!("{CHR_1}:10150|+|10M|19"), 1).unwrap());
    t.clear_clusters_and_svs();
    t.add(&[&sgl1]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    set_sv_gene_data_test(&t.all_variants, gc, 0);
    fa.annotate_transcripts(&t.all_variants, false);
    let clusters = t.clusters().clone();
    fa.run(&t.sample_id.clone(), &t.all_variants.clone(), &clusters, &mut t.vis).unwrap();
    assert_eq!(1, fa.fusions.len());
    let f = fa.fusions[0].clone();
    assert_eq!(sgl1.id(), f.upstream_trans().gene().var_id());
    assert_eq!(sgl1.id(), f.downstream_trans().gene().var_id());
    assert_eq!(10150, f.downstream_trans().gene().position());

    let sgl1 = create_sgl(next(), CHR_1, 10150, -1);
    sgl1.add_sgl_mapping(sgl_mapping(CHR_1, 1150, 1));
    let var1 = create_del(next(), CHR_1, 14000, 15000);
    let c = cluster_svs(&mut t, &[&sgl1, &var1]);
    assert_eq!(1, c.chains().len());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(1, fa.fusions.len());
    let f = fa.fusions[0].clone();
    assert_eq!(sgl1.id(), f.upstream_trans().gene().var_id());
    assert_eq!(sgl1.id(), f.downstream_trans().gene().var_id());

    let sgl1 = create_sgl(next(), CHR_2, 100, -1);
    sgl1.add_sgl_mapping(sgl_mapping(CHR_1, 1150, 1));
    let var1 = create_bnd(next(), CHR_1, 10150, -1, CHR_2, 200, 1);
    let var2 = create_bnd(next(), CHR_1, 500, 1, CHR_2, 500, -1);
    let c = cluster_svs(&mut t, &[&sgl1, &var1, &var2]);
    assert_eq!(1, c.chains().len());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(1, fa.fusions.len());
    let f = fa.fusions[0].clone();
    assert_eq!(sgl1.id(), f.upstream_trans().gene().var_id());
    assert_eq!(var1.id(), f.downstream_trans().gene().var_id());

    let sgl1 = create_sgl(next(), CHR_2, 100, -1);
    sgl1.add_sgl_mapping(sgl_mapping(CHR_1, 1150, 1));
    let var1 = create_bnd(next(), CHR_2, 200, 1, "3", 200, 1);
    let sgl2 = create_sgl(next(), "3", 100, -1);
    sgl2.add_sgl_mapping(sgl_mapping(CHR_1, 10150, -1));
    let c = cluster_svs(&mut t, &[&sgl1, &sgl2, &var1]);
    assert_eq!(1, c.chains().len());
    assert_eq!(3, c.chains()[0].sv_count());
    run_fusions(&mut t, gc, &mut fa, pre_gene);
    assert_eq!(1, fa.fusions.len());
    let f = fa.fusions[0].clone();
    assert_eq!(sgl1.id(), f.upstream_trans().gene().var_id());
    assert_eq!(sgl2.id(), f.downstream_trans().gene().var_id());
}

#[test]
fn test_sgl_mapped_inf_fusions() {
    let mut t = LinxTester::new();
    let mut env = test_genes_env();
    env.kfc.add_data(KnownFusionData::new(KNOWN_PAIR, GENE_NAME_1, GENE_NAME_2, "", ""));
    let (gc, mut fa) = env.analyser();
    let sgl = create_sgl(1, CHR_1, 1150, 1);
    let inf = create_inf(2, CHR_2, 100050, -1);
    sgl.add_sgl_mapping(SglMapping::from(&format!("{CHR_2}:100000|+|10M|19"), 1).unwrap());
    let bnd = create_bnd(3, CHR_1, 10020, -1, CHR_2, 100200, 1);
    t.add(&[&sgl, &inf, &bnd]);
    t.pre_clustering_init();
    assert_eq!(4, t.all_variants.len());
    let sgl_inf = t.all_variants[3].clone();
    t.cluster_and_analyse();
    set_sv_gene_data_test(&t.all_variants, gc, 0);
    fa.annotate_transcripts(&t.all_variants, false);
    let clusters = t.clusters().clone();
    fa.run(&t.sample_id.clone(), &t.all_variants.clone(), &clusters, &mut t.vis).unwrap();
    assert_eq!(1, fa.fusions.len());
    let f = fa.fusions[0].clone();
    assert_eq!(sgl_inf.id(), f.upstream_trans().gene().var_id());
    assert_eq!(bnd.id(), f.downstream_trans().gene().var_id());
}

#[test]
fn test_known_unmappable_fusions() {
    let mut t = LinxTester::new();
    let mut env = test_genes_env();
    let chr4 = "4";
    add_gene_data(&mut env.ensembl, chr4, vec![create_ensembl_gene_data("ENSG0IGH", "IGH", chr4, POS_STRAND, 100, 1500)]);
    let starts = generate_exon_starts(100, 7, 100, 100);
    add_trans_exon_data(&mut env.ensembl, "ENSG0IGH", vec![create_trans_exons("ENSG0IGH", 1, POS_STRAND, &starts, 100, Some(100), Some(799), true, "")]);
    let ig_region = format!("ENHANCER_RANGE={};{};{};{}", POS_STRAND, chr4, 50, 2000);
    let alts = "ALTS=ALT;1;20000;25000;ALT;GS;50000;60000";
    let mut kf = KnownFusionData::new(KNOWN_PAIR, GENE_NAME_1, GENE_NAME_3, "", "");
    kf.apply_overrides(alts).unwrap();
    env.kfc.add_data(kf);
    let mut kf = KnownFusionData::new(ENHANCER_KNOWN_PAIR, "IGH", GENE_NAME_3, "", "");
    kf.apply_overrides(&format!("{ig_region} {alts}")).unwrap();
    env.kfc.add_data(kf);
    let (gc, mut fa) = env.analyser();
    fa.cache_special_fusion_genes();

    let run_one = |t: &mut LinxTester, fa: &mut FusionDisruptionAnalyser, sgl: &Sv| {
        t.clear_clusters_and_svs();
        t.add(&[sgl]);
        t.pre_clustering_init();
        t.cluster_and_analyse();
        set_sv_gene_data_test(&t.all_variants, gc, 0);
        fa.annotate_transcripts(&t.all_variants, false);
        let clusters = t.clusters().clone();
        fa.run(&t.sample_id.clone(), &t.all_variants.clone(), &clusters, &mut t.vis).unwrap();
        assert_eq!(1, fa.fusions.len());
        fa.fusions[0].clone()
    };
    let sgl1 = create_sgl(1, CHR_1, 1150, 1);
    sgl1.add_sgl_mapping(sgl_mapping(CHR_1, 24000, 1));
    let f = run_one(&mut t, &mut fa, &sgl1);
    assert_eq!(sgl1.id(), f.upstream_trans().gene().var_id());
    assert_eq!(sgl1.id(), f.downstream_trans().gene().var_id());
    assert!(f.reportable());
    assert_eq!(KNOWN_PAIR, f.known_type());

    let sgl1 = create_sgl(2, CHR_1, 1150, 1);
    sgl1.add_sgl_mapping(sgl_mapping("GS", 55000, 1));
    let f = run_one(&mut t, &mut fa, &sgl1);
    assert_eq!(sgl1.id(), f.upstream_trans().gene().var_id());
    assert_eq!(sgl1.id(), f.downstream_trans().gene().var_id());
    assert!(f.reportable());
    assert_eq!(KNOWN_PAIR, f.known_type());

    let sgl1 = create_sgl(3, chr4, 1000, 1);
    sgl1.add_sgl_mapping(sgl_mapping("GS", 55000, 1));
    let f = run_one(&mut t, &mut fa, &sgl1);
    assert!(f.reportable());
    assert_eq!(ENHANCER_KNOWN_PAIR, f.known_type());
    assert_eq!("IGH", f.upstream_trans().gene_name());
    assert_eq!(GENE_NAME_3, f.downstream_trans().gene_name());
}

// DisruptionTest

fn add_disruption_test_gene(cache: &mut EnsemblDataCache, chromosome: &str, gene_id: &str) {
    add_gene_data(cache, chromosome, vec![create_ensembl_gene_data(gene_id, gene_id, chromosome, 1, 10000, 120000)]);
    let starts = [10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000];
    add_trans_exon_data(cache, gene_id, vec![create_trans_exons(gene_id, 1, 1, &starts, 100, Some(20002), Some(110098), true, "")]);
}
fn disruptive(v: &Sv, is_start: bool) -> bool {
    v.genes_list(is_start)[0].transcripts()[0].is_disruptive()
}
fn prepare_gene_annotations(t: &mut LinxTester, gc: &GeneCache, fa: &FusionDisruptionAnalyser, vars: &[&Sv], pre_gene: i32) {
    t.clear_clusters_and_svs();
    t.add(vars);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    set_sv_gene_data_test(&t.all_variants, gc, pre_gene);
    fa.annotate_transcripts(&t.all_variants, true);
}

#[test]
fn test_simple_disruptions() {
    let mut t = LinxTester::new();
    let mut env = FusionEnv::new();
    let chr = CHR_1;
    add_disruption_test_gene(&mut env.ensembl, chr, GENE_ID_1);
    let (gc, mut fa) = env.analyser();
    fa.disruption_finder.add_disruption_gene(GENE_ID_1);

    let var1 = create_del(t.next_var_id(), chr, 11000, 12000);
    let var2 = create_dup(t.next_var_id(), chr, 21000, 22000);
    let var3 = create_ins(t.next_var_id(), chr, 31000, 32000);
    let var4 = create_inv(t.next_var_id(), chr, 41000, 42000, -1);
    let var5 = create_bnd(t.next_var_id(), chr, 51000, -1, "2", 1000, -1);
    let var6 = create_sgl(t.next_var_id(), chr, 61000, -1);
    let var7 = create_del(t.next_var_id(), chr, 68000, 72000);
    let var8 = create_dup(t.next_var_id(), chr, 78000, 96000);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8], DEFAULT_PRE_GENE);
    for v in &t.all_variants {
        assert_eq!(1, v.genes_list(true).len());
    }
    for v in [&var1, &var2, &var3] {
        assert!(!disruptive(v, true));
    }
    for v in [&var4, &var5, &var6, &var7, &var8] {
        assert!(disruptive(v, true));
    }
    let clusters = t.clusters().clone();
    fa.disruption_finder.find_reportable_disruptions(&t.all_variants, &clusters, None);
    assert_eq!(8, fa.disruption_finder.disruptions.len());

    for (s, e, empty_start) in [(10050, 10500, false), (9500, 10500, true), (200, 10500, true), (10500, 15000, false)] {
        let v = create_dup(t.next_var_id(), chr, s, e);
        prepare_gene_annotations(&mut t, gc, &fa, &[&v], 0);
        if empty_start && s < 10000 {
            if s == 9500 {
                assert!(v.genes_list(true).is_empty());
            }
            assert!(!disruptive(&v, false));
        } else {
            assert!(!disruptive(&v, true));
            if s == 10050 {
                assert!(!disruptive(&v, false));
            }
        }
        let clusters = t.clusters().clone();
        fa.disruption_finder.find_reportable_disruptions(&t.all_variants, &clusters, None);
        assert!(fa.disruption_finder.disruptions.is_empty());
    }
}

#[test]
fn test_chained_disruptions() {
    let mut t = LinxTester::new();
    let mut env = FusionEnv::new();
    let (chr, chr2) = (CHR_1, CHR_2);
    add_disruption_test_gene(&mut env.ensembl, chr, GENE_ID_1);
    add_disruption_test_gene(&mut env.ensembl, chr2, GENE_ID_2);
    let (gc, fa) = env.analyser();

    let var1 = create_bnd(t.next_var_id(), chr, 11000, 1, "3", 100, -1);
    let var2 = create_bnd(t.next_var_id(), chr, 12000, -1, "3", 200, 1);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2], 0);
    assert_eq!(1, t.clusters().len());
    for v in &t.all_variants {
        assert_eq!(1, v.genes_list(true).len());
    }
    assert!(!disruptive(&var1, true));
    assert!(!disruptive(&var2, true));

    let var1 = create_inv(t.next_var_id(), chr, 11000, 14000, 1);
    let var2 = create_inv(t.next_var_id(), chr, 13900, 18000, -1);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2], 0);
    for v in [&var1, &var2] {
        assert!(!disruptive(v, true));
        assert!(!disruptive(v, false));
    }

    let var1 = create_del(t.next_var_id(), chr, 11000, 12000);
    let var2 = create_bnd(t.next_var_id(), chr, 13000, 1, "3", 100, -1);
    let var3 = create_bnd(t.next_var_id(), chr, 14000, -1, "3", 200, 1);
    let var4 = create_del(t.next_var_id(), chr, 15000, 16000);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3, &var4], 0);
    for v in [&var1, &var2, &var3, &var4] {
        assert!(!disruptive(v, true));
    }

    let var1 = create_del(t.next_var_id(), chr, 11000, 12000);
    let var2 = create_bnd(t.next_var_id(), chr, 13000, 1, "3", 100, -1);
    let var3 = create_bnd(t.next_var_id(), chr, 14000, 1, "3", 200, 1);
    let var4 = create_dup(t.next_var_id(), chr, 13500, 16000);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3, &var4], 0);
    assert!(!disruptive(&var1, true));
    assert!(!disruptive(&var1, false));
    assert!(disruptive(&var2, true));
    assert!(disruptive(&var3, true));
    assert!(!disruptive(&var4, true));

    let var1 = create_del(t.next_var_id(), chr, 11000, 12000);
    let var2 = create_bnd(t.next_var_id(), chr, 13000, 1, chr2, 18000, -1);
    let var3 = create_bnd(t.next_var_id(), chr2, 22000, 1, "3", 1000, -1);
    let var4 = create_bnd(t.next_var_id(), chr, 14000, -1, "3", 3000, 1);
    let var5 = create_del(t.next_var_id(), chr, 15000, 16000);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3, &var4, &var5], 0);
    assert_eq!(1, t.clusters().len());
    assert!(!disruptive(&var1, true));
    assert!(disruptive(&var2, true));
    assert!(disruptive(&var2, false));
    assert!(disruptive(&var3, true));
    assert!(disruptive(&var4, true));
    assert!(!disruptive(&var5, true));
}

#[test]
fn test_chained_disruptions2() {
    let mut t = LinxTester::new();
    let mut env = FusionEnv::new();
    let (chr, chr2) = (CHR_1, CHR_2);
    add_disruption_test_gene(&mut env.ensembl, chr, GENE_ID_1);
    add_disruption_test_gene(&mut env.ensembl, chr2, GENE_ID_2);
    let (gc, fa) = env.analyser();

    let var1 = create_inv(t.next_var_id(), chr, 1000, 15000, 1);
    let var2 = create_inv(t.next_var_id(), chr, 2000, 14000, -1);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2], 0);
    assert_eq!(1, t.clusters().len());
    assert!(!disruptive(&var1, false));
    assert!(!disruptive(&var2, false));

    let var1 = create_inv(t.next_var_id(), chr, 1000, 15000, 1);
    let var2 = create_bnd(t.next_var_id(), chr, 14000, -1, chr2, 25000, -1);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2], 0);
    assert_eq!(1, t.clusters().len());
    assert!(var1.genes_list(true).is_empty());
    assert!(disruptive(&var1, false));
    assert!(disruptive(&var2, true));
    assert!(disruptive(&var2, false));

    let remote = "3";
    let var1 = create_inv(t.next_var_id(), chr, 12000, 14000, 1);
    let var2 = create_bnd(t.next_var_id(), chr, 13000, -1, remote, 100, -1);
    let var3 = create_bnd(t.next_var_id(), chr, 17000, 1, remote, 200, 1);
    let var4 = create_inv(t.next_var_id(), chr, 16000, 18000, -1);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3, &var4], 0);
    assert_eq!(1, t.clusters().len());
    assert!(!disruptive(&var1, true));
    assert!(!disruptive(&var1, false));
    assert!(!disruptive(&var2, true));
    assert!(!disruptive(&var3, true));
    assert!(!disruptive(&var4, true));

    let var1 = create_del(t.next_var_id(), chr, 11800, 14500);
    let var2 = create_bnd(t.next_var_id(), chr, 12000, -1, remote, 100, 1);
    let var3 = create_inv(t.next_var_id(), chr, 12100, 14100, 1);
    set_assembled_link_info(&var2, true, &var3, true);
    let var4 = create_bnd(t.next_var_id(), chr, 14000, -1, remote, 200, -1);
    set_assembled_link_info(&var3, false, &var4, true);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3, &var4], 0);
    assert_eq!(1, t.clusters().len());
    assert!(!disruptive(&var1, true));
    assert!(!disruptive(&var1, false));
    assert!(!disruptive(&var3, true));
    assert!(!disruptive(&var3, false));
    assert!(!disruptive(&var4, true));

    let var1 = create_del(t.next_var_id(), chr, 12000, 25500);
    let var2 = create_bnd(t.next_var_id(), chr, 12100, -1, remote, 100, 1);
    let var3 = create_bnd(t.next_var_id(), chr, 12500, 1, remote, 200, -1);
    set_assembled_link_info(&var2, true, &var3, true);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3], 0);
    assert_eq!(1, t.clusters().len());
    assert!(disruptive(&var1, true));
    assert!(disruptive(&var1, false));
    assert!(disruptive(&var2, true));
    assert!(disruptive(&var3, true));

    let var1 = create_del(t.next_var_id(), chr, 12000, 25500);
    let var2 = create_bnd(t.next_var_id(), chr, 12100, -1, remote, 100, 1);
    let var3 = create_inv(t.next_var_id(), chr, 12500, 25400, 1);
    let var4 = create_bnd(t.next_var_id(), chr, 25000, -1, remote, 200, -1);
    prepare_gene_annotations(&mut t, gc, &fa, &[&var1, &var2, &var3, &var4], 0);
    assert_eq!(1, t.clusters().len());
    assert!(disruptive(&var1, true));
    assert!(disruptive(&var1, false));
    assert!(disruptive(&var2, true));
    assert!(disruptive(&var3, true));
    assert!(disruptive(&var3, false));
    assert!(disruptive(&var4, true));
}
