//! germline: GermlineDisruptionTest

use hmf_common::java_hash::JavaHashMap;
use hmf_common::purple_files::DriverCatalog;

use super::drivers::create_driver_gene;
use super::fusion::{set_sv_gene_data_test, FusionEnv};
use super::*;
use crate::fusion::analyser::FusionDisruptionAnalyser;
use crate::output::LinxBreakend;

struct GermlineTest {
    t: LinxTester,
    gc: &'static GeneCache,
    fa: FusionDisruptionAnalyser,
}

impl GermlineTest {
    fn new() -> Self {
        let mut env = FusionEnv::new();
        env.is_germline = true;
        let mut dg = JavaHashMap::new();
        dg.insert(GENE_NAME_1.to_string(), create_driver_gene(GENE_NAME_1));
        env.driver_genes = dg;
        add_gene_data(&mut env.ensembl, CHR_1, vec![create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, 10000, 120000)]);
        let td = create_trans_exons(GENE_ID_1, 1, POS_STRAND, &[10000, 11000, 12000, 13000, 14000], 200, Some(11100), Some(14100), true, "");
        add_trans_exon_data(&mut env.ensembl, GENE_ID_1, vec![td]);
        let (gc, mut fa) = env.analyser();
        fa.disruption_finder.add_disruption_gene(GENE_ID_1);
        GermlineTest { t: LinxTester::with_germline(true), gc, fa }
    }
    fn cluster(&mut self, vars: &[&Sv]) {
        self.t.clear_clusters_and_svs();
        self.t.add(vars);
        self.t.pre_clustering_init();
        self.t.cluster_and_analyse();
    }
    fn find_disruptions(&mut self) {
        let clusters = self.t.clusters().clone();
        self.fa.disruption_finder.find_reportable_disruptions(&self.t.all_variants, &clusters, None);
    }
    /// populateGermlineResults: (germline SVs, drivers)
    fn results(&self) -> (usize, usize) {
        let mut svs: Vec<String> = Vec::new();
        let mut breakends: Vec<LinxBreakend> = Vec::new();
        let mut drivers: Vec<DriverCatalog> = Vec::new();
        let df = &self.fa.disruption_finder;
        df.germline.as_ref().unwrap().populate_germline_svs(&df.disruptions, &mut svs, &mut breakends, &mut drivers);
        (svs.len(), drivers.len())
    }
}

fn sgl_mapped(t: &mut LinxTester, pos: i32, map_pos: i32) -> Sv {
    let v = create_sgl(t.next_var_id(), CHR_1, pos, 1);
    v.add_sgl_mapping(SglMapping { chromosome: CHR_1.into(), position: map_pos, orientation: -1, cigar: String::new(), qual_score: 1 });
    v
}

#[test]
fn test_whole_gene_deletions() {
    let mut g = GermlineTest::new();
    let var1 = create_del(g.t.next_var_id(), CHR_1, 5000, 150000);
    g.cluster(&[&var1]);
    g.find_disruptions();
    assert_eq!((1, 0), g.results());
}

#[test]
fn test_sgl_mapped_del() {
    let mut g = GermlineTest::new();
    for (pos, map_pos, expected) in [(5000, 150000, 1), (5000, 12100, 1), (11500, 13500, 1), (11500, 11800, 0)] {
        let v = sgl_mapped(&mut g.t, pos, map_pos);
        g.cluster(&[&v]);
        g.find_disruptions();
        let (svs, drivers) = g.results();
        assert_eq!(expected, svs, "SGL at {pos} mapped to {map_pos}");
        assert_eq!(0, drivers);
    }
    let v = sgl_mapped(&mut g.t, 11500, 13500);
    let var2 = create_sgl(g.t.next_var_id(), CHR_1, 15000, 1);
    g.cluster(&[&v, &var2]);
    assert_eq!(1, g.t.clusters().len());
    g.find_disruptions();
    assert_eq!((1, 0), g.results());
}

#[test]
fn test_sgl_mapped_dup() {
    let mut g = GermlineTest::new();
    let v = sgl_mapped(&mut g.t, 13500, 11500);
    g.cluster(&[&v]);
    g.find_disruptions();
    assert_eq!((1, 0), g.results());
    let v = sgl_mapped(&mut g.t, 13500, 9000);
    g.cluster(&[&v]);
    g.find_disruptions();
    assert_eq!(0, g.results().0);
}

#[test]
fn test_pseudogene_deletions() {
    let mut g = GermlineTest::new();
    let var1 = create_del(g.t.next_var_id(), CHR_1, 11200, 12000);
    g.cluster(&[&var1]);
    set_sv_gene_data_test(&g.t.all_variants, g.gc, 100000);
    g.fa.annotate_transcripts(&g.t.all_variants, true);
    g.find_disruptions();
    assert_eq!((1, 0), g.results());
}
