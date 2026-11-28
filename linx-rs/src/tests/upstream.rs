//! The port's own tests of what upstream 7a4490caa6 (IG fusions, AUS-509) changed: coding exonic upstream breakends fuse
//! with any valid downstream transcript, the enhancer 3' gene upstream distance (1 Mb) and its scope (LX6), and the
//! visualiser's enhancer exon regions. Synthetic fixtures.

use hmf_common::ensembl::BIOTYPE_PROTEIN_CODING;
use hmf_common::java_hash::JavaHashMap;
use hmf_common::known_fusion::KnownFusionData;
use hmf_common::known_fusion::KnownFusionType::*;

use super::fusion::FusionEnv;
use super::*;
use crate::compat::with_java;
use crate::fusion::finder::specific_pre_gene_distances;
use crate::fusion::MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN;
use crate::gene::{find_gene_annotations_by_sv, set_sv_gene_data, TranscriptRegionType};
use crate::vis::{VisDataWriter, VisGeneAnnotationType, VisGeneData};

const POS: i8 = 1;

/// GENE1 (5', exons 100-1000, coding 150-950) and GENE2 (3', exons 10100-11000, coding from exon 1), both positive strand
fn two_gene_env() -> FusionEnv {
    let mut env = FusionEnv::new();
    let genes = vec![
        create_ensembl_gene_data("ENSG0001", "GENE1", CHR_1, POS, 100, 1000),
        create_ensembl_gene_data("ENSG0002", "GENE2", CHR_1, POS, 10100, 11000),
    ];
    add_gene_data(&mut env.ensembl, CHR_1, genes);
    let t1 = create_trans_exons("ENSG0001", 1, POS, &[100, 300, 500, 700, 900], 100, Some(150), Some(950), true, BIOTYPE_PROTEIN_CODING);
    add_trans_exon_data(&mut env.ensembl, "ENSG0001", vec![t1]);
    let t2 =
        create_trans_exons("ENSG0002", 2, POS, &[10100, 10300, 10500, 10700, 10900], 100, Some(10150), Some(10950), true, BIOTYPE_PROTEIN_CODING);
    add_trans_exon_data(&mut env.ensembl, "ENSG0002", vec![t2]);
    env
}

#[test]
fn test_coding_exonic_upstream_to_pre_gene() {
    // a coding exonic 5' breakend and a 3' breakend before GENE2's start (the first splice acceptor, exon 2, is coding):
    // rejected until 7a4490caa6 ("coding exonic to non-exonic"), now a fusion through exon 2
    let env = two_gene_env();
    let (gc, fa) = env.analyser();
    let up = find_gene_annotations_by_sv(gc, 0, true, CHR_1, 701, 1, 100000);
    up[0].set_positional_data(CHR_1, 701, 1);
    let down = find_gene_annotations_by_sv(gc, 0, false, CHR_1, 10050, -1, 100000);
    down[0].set_positional_data(CHR_1, 10050, -1);
    let d = &down[0].transcripts()[0];
    assert_eq!(TranscriptRegionType::UPSTREAM, d.region_type());
    assert!(d.is_coding());
    let fusions = fa.finder.find_fusions(&up, &down);
    assert_eq!(1, fusions.len());
    let f = &fusions[0];
    assert_eq!(TranscriptRegionType::EXONIC, f.upstream_trans().region_type());
    assert!(f.upstream_trans().is_coding());
    // the exonic breakend in exon 4 has exon 4's start phase (1); GENE2 enters at exon 2 (phase 0): matched through the
    // alternative phasing of exon 1's donor (end phase 0), skipping exons 2-4
    assert_eq!(1, f.upstream_trans().phase);
    assert_eq!(0, f.downstream_trans().phase);
    assert!(f.phase_matched());
    assert_eq!(3, f.exons_skipped(true));
    assert_eq!(0, f.exons_skipped(false));
}

#[test]
fn test_coding_exonic_upstream_other_sv() {
    // the 5' breakend (SV 0) and the 3' breakend (SV 1) on different SVs, as in a chain: rejected until 7a4490caa6
    // ("up coding exonic diff SVs"), now as for one SV
    let env = two_gene_env();
    let (gc, fa) = env.analyser();
    let up = find_gene_annotations_by_sv(gc, 0, true, CHR_1, 701, 1, 100000);
    up[0].set_positional_data(CHR_1, 701, 1);
    let down_same = find_gene_annotations_by_sv(gc, 0, false, CHR_1, 10650, -1, 100000);
    down_same[0].set_positional_data(CHR_1, 10650, -1);
    let down_other = find_gene_annotations_by_sv(gc, 1, false, CHR_1, 10650, -1, 100000);
    down_other[0].set_positional_data(CHR_1, 10650, -1);
    let same = fa.finder.find_fusions(&up, &down_same);
    let other = fa.finder.find_fusions(&up, &down_other);
    assert_eq!(1, same.len());
    assert_eq!(1, other.len());
    assert_eq!(same[0].phase_matched(), other[0].phase_matched());
    assert!(other[0].phase_matched());
    assert_eq!(1, other[0].downstream_trans().gene().var_id());
}

/// an enhancer known pair (IGX enhancer region on chromosome 2, 3' gene TGT) and genes on chromosome 1: TGT starting at
/// 2,000,000; UNR, unrelated, 200 kb downstream of the breakend at 1,400,000 (in TGT's 1 Mb pre-gene window); NEAR, 50 kb
fn window_env() -> FusionEnv {
    let mut env = FusionEnv::new();
    let mut genes = Vec::new();
    for (id, name, start) in
        [("ENSG0011", "NEAR", 1_450_000), ("ENSG0012", "UNR", 1_600_000), ("ENSG0013", "TGT", 2_000_000), ("ENSG0014", "FAR", 2_600_000)]
    {
        genes.push(create_ensembl_gene_data(id, name, CHR_1, POS, start, start + 5000));
        let exons = [start, start + 2000, start + 4000];
        let t = create_trans_exons(id, start / 1000, POS, &exons, 100, Some(start + 50), Some(start + 4050), true, BIOTYPE_PROTEIN_CODING);
        add_trans_exon_data(&mut env.ensembl, id, vec![t]);
    }
    add_gene_data(&mut env.ensembl, CHR_1, genes);
    let mut kf = KnownFusionData::new(ENHANCER_KNOWN_PAIR, "IGX", "TGT", "", "");
    kf.apply_overrides(&format!("ENHANCER_RANGE=1;{CHR_2};100000;200000")).unwrap();
    env.kfc.add_data(kf);
    env
}

fn window_genes(java: bool) -> Vec<String> {
    window_genes_with(java, 100000)
}

fn window_genes_with(java: bool, default_distance: i32) -> Vec<String> {
    let env = window_env();
    let (gc, _fa) = env.analyser();
    let dist = specific_pre_gene_distances(&env.kfc);
    assert_eq!(Some(&MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN), dist.get(&"TGT".to_string()));
    let var = create_bnd(0, CHR_1, 1_400_000, -1, CHR_2, 150_000, 1);
    with_java(java, || set_sv_gene_data(std::slice::from_ref(&var), gc, default_distance, &dist));
    var.genes_list(true).iter().map(|g| g.gene_name().to_string()).collect()
}

#[test]
fn test_enhancer_window_pre_gene_distance() {
    assert_eq!(1_000_000, MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN);
    // Java: every gene within 1 Mb of a breakend in the enhancer 3' gene's window
    assert_eq!(vec!["NEAR", "UNR", "TGT"], window_genes(true));
    // corrected (LX6): the 1 Mb only for TGT, the default 100 kb for the others
    assert_eq!(vec!["NEAR", "TGT"], window_genes(false));
    // a default distance above 1 Mb: Java searches the window with the smaller 1 Mb, so FAR (1.2 Mb) is missed; corrected,
    // every gene within the default distance
    assert_eq!(vec!["NEAR", "UNR", "TGT"], window_genes_with(true, 2_000_000));
    assert_eq!(vec!["NEAR", "UNR", "TGT", "FAR"], window_genes_with(false, 2_000_000));
}

#[test]
fn test_enhancer_window_outside() {
    // a breakend outside every window: the default distance in both modes
    let env = window_env();
    let (gc, _fa) = env.analyser();
    let dist = specific_pre_gene_distances(&env.kfc);
    for java in [true, false] {
        let var = create_bnd(0, CHR_1, 900_000, -1, CHR_2, 150_000, 1);
        with_java(java, || set_sv_gene_data(std::slice::from_ref(&var), gc, 100000, &dist));
        assert!(var.genes_list(true).is_empty());
    }
    let _ = JavaHashMap::<String, i32>::new();
}

#[test]
fn test_enhancer_vis_exon_region() {
    // an enhancer fusion's upstream side is drawn over the region it was made from; without one, nothing (as Java)
    let env = FusionEnv::new();
    let gc = leak_cache(env.ensembl.clone());
    let dir = std::env::temp_dir().join(format!("linx_rs_vis_enhancer_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let writer = VisDataWriter::new(dir.to_str().unwrap(), Some(gc), true, false);
    let mut sample = VisSampleData::default();
    sample.set_sample_id("S1");
    sample.add_gene_exon_data(
        VisGeneData::new(3, "ENSG0021", "IGXV1", "@IGX", 0, "chr2", VisGeneAnnotationType::FUSION).with_region_override(Some((100000, 200000))),
    );
    sample.add_gene_exon_data(VisGeneData::new(4, "ENSG0022", "TRXV1", "@TRX", 0, "chr7", VisGeneAnnotationType::FUSION));
    writer.write_gene_data(&sample).unwrap();
    let text = std::fs::read_to_string(dir.join("S1.linx.vis_gene_exon.tsv")).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(2, lines.len());
    assert_eq!("3\t@IGX\t@IGX\tchr2\tFUSION\t1\t100000\t200000", lines[1]);
    std::fs::remove_dir_all(&dir).ok();
}
