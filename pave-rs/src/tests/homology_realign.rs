//! HomologyRealignTest

use super::*;
use hmf_common::ensembl::GeneData;
use hmf_common::regions::RefGenomeVersion;

use crate::gene_cache::GeneDataCache;
use crate::impact::effects::VariantEffect;
use crate::impact::utils::{create_right_aligned_variant, find_variant_impacts};
use crate::variant::NO_LOCAL_PHASE_SET;

/// GeneTestUtils.createEnsemblGeneData
pub fn create_ensembl_gene_data(gene_id: &str, gene_name: &str, chromosome: &str, strand: i8, start: i32, end: i32) -> GeneData {
    GeneData {
        gene_id: gene_id.to_string(),
        gene_name: gene_name.to_string(),
        chromosome: chromosome.to_string(),
        strand,
        gene_start: start,
        gene_end: end,
        karyotype_band: String::new(),
        synonyms: String::new(),
    }
}

/// a gene cache holding one gene and its transcript
pub fn gene_cache_with(gene: GeneData, trans: &TranscriptData) -> &'static GeneDataCache {
    let mut cache = GeneDataCache::empty(RefGenomeVersion::V37);
    cache.ensembl_mut().chr_gene_data_map_mut().insert(gene.chromosome.clone(), vec![gene.clone()]);
    cache.ensembl_mut().transcript_data_map_mut().insert(gene.gene_id.clone(), vec![trans.clone()]);
    Box::leak(Box::new(cache))
}

fn realign_and_find(var: &mut VariantData, ref_genome: &MockRefGenome, cache: &'static GeneDataCache) {
    let ra = create_right_aligned_variant(var, ref_genome).unwrap();
    var.set_realigned_variant(ra);
    assert!(var.realigned_variant().is_some());
    find_variant_impacts(var, cache, None, ref_genome).unwrap();
}

fn gene_impact(var: &VariantData) -> &crate::impact::VariantTransImpact {
    let impacts = var.gene_impacts(GENE_NAME_1).expect("no impacts for the gene");
    assert_eq!(1, impacts.len());
    &impacts[0]
}

#[test]
fn test_realignment() {
    let ref_genome = create_mock_genome(201);
    let rb = ref_genome.bases().to_string();
    let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[20, 50, 80], 10, Some(25), Some(65), false, "");
    let cache = gene_cache_with(create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, 20, 90), trans_data_pos);

    // intronic to splice
    let pos = 40; // intronic DEL
    let r = sub(&rb, pos, pos + 2);
    let a = sub(&rb, pos, pos + 1);
    let mut var = VariantData::new(CHR_1, pos, &r, &a);
    let alt_bases = r[1..].to_string();
    var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 9);
    realign_and_find(&mut var, &ref_genome, cache);
    let impact = gene_impact(&var);
    assert!(!impact.realigned());
    assert_eq!(Some(VariantEffect::Intronic), impact.top_effect());

    // realigning from the previous donor into the intron
    let pos = 33;
    let r = sub(&rb, pos, pos + 2);
    let a = sub(&rb, pos, pos + 1);
    let mut var = VariantData::new(CHR_1, pos, &r, &a);
    let alt_bases = r[1..].to_string();
    var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 9);
    realign_and_find(&mut var, &ref_genome, cache);
    let impact = gene_impact(&var);
    assert!(impact.realigned());
    assert_eq!(Some(VariantEffect::Intronic), impact.top_effect());

    // splice to an exonic inframe INS
    let pos = 48;
    let r = sub(&rb, pos, pos + 1);
    let alt_bases = "AAA";
    let a = cat!(r, alt_bases);
    let mut var = VariantData::new(CHR_1, pos, &r, &a);
    var.set_variant_details(NO_LOCAL_PHASE_SET, alt_bases, alt_bases, 2);
    realign_and_find(&mut var, &ref_genome, cache);
    let impact = gene_impact(&var);
    assert!(impact.realigned());
    assert_eq!(Some(VariantEffect::InframeInsertion), impact.top_effect());
}

#[test]
fn test_splice_realignment() {
    let mut rb = generate_test_bases(10);
    rb += &generate_test_bases(10);
    rb += &rc("TAA");
    rb += "ATACCTGCT"; // Y-R-S going up, X (20-22), Y (23-25), R (26-28), S (29-31)
    rb += "CTATAGAGCG"; // intron 32-41
    rb += "CTTCTCCCT"; // K (42-44), E (45-47), R (48-50)
    rb += &rc("ATG"); // M (51-52)
    rb += &generate_test_bases(10);

    let neg_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[10, 42], 21, Some(20), Some(52), false, "");
    let mut ref_genome = MockRefGenome::default();
    ref_genome.map.insert(CHR_1.to_string(), rb.clone());

    // exonic to splice
    let pos = 26;
    let r = sub(&rb, pos, pos + 6);
    let a = sub(&rb, pos, pos + 1);
    let mut var = VariantData::new(CHR_1, pos, &r, &a);
    var.set_variant_details(NO_LOCAL_PHASE_SET, "CT", "C", 2);
    let cache = gene_cache_with(create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, NEG_STRAND, neg_trans.trans_start, neg_trans.trans_end), neg_trans);
    realign_and_find(&mut var, &ref_genome, cache);
    let impact = gene_impact(&var);
    assert_eq!(Some(VariantEffect::Frameshift), impact.top_effect());
    assert!(!impact.realigned());
    let ra_impact = var.realigned_impact(GENE_NAME_1, impact).unwrap();
    assert!(ra_impact.has_effect(VariantEffect::SpliceAcceptor));
}

#[test]
fn test_splice_vs_frameshift_realignment() {
    let mut rb = generate_test_bases(10);
    rb += &generate_test_bases(10);
    let mut coding = get_amino_acid_codon(b'M', 0);
    coding += &get_amino_acids_codons("CFR", false);
    coding += &get_amino_acid_codon(b'G', 2); // force to 'GGG'
    coding += &get_amino_acids_codons("D", false);
    coding += "TAA";
    rb += &coding[..10];
    rb += &(generate_test_bases(17) + "CAG"); // intron ending in a splice acceptor motif with 'G'
    rb += &coding[10..];
    rb += &generate_test_bases(10);
    rb += &generate_test_bases(10);

    let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 50], 19, Some(20), Some(60), false, "");
    let mut ref_genome = MockRefGenome::default();
    ref_genome.map.insert(CHR_1.to_string(), rb.clone());

    // the initial variant is a splice acceptor, the realigned one a frameshift, which realignment's
    // priority rules select
    let pos = 48;
    let r = sub(&rb, pos, pos + 2);
    let a = sub(&rb, pos, pos + 1);
    let mut var = VariantData::new(CHR_1, pos, &r, &a);
    var.set_variant_details(NO_LOCAL_PHASE_SET, "G", "G", 7);
    let cache = gene_cache_with(create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, pos_trans.trans_start, pos_trans.trans_end), pos_trans);
    realign_and_find(&mut var, &ref_genome, cache);
    let impact = gene_impact(&var);
    assert_eq!(Some(VariantEffect::Frameshift), impact.top_effect());
    assert!(impact.realigned());

    // an insertion of a base at A2 which keeps the acceptor motif but with a spare base before the exon
    let r = sub(&rb, pos, pos + 1);
    let mut var = VariantData::new(CHR_1, pos, &r, "AG");
    var.set_variant_details(NO_LOCAL_PHASE_SET, "G", "G", 7);
    realign_and_find(&mut var, &ref_genome, cache);
    let impact = gene_impact(&var);
    assert_eq!(Some(VariantEffect::Frameshift), impact.top_effect());
    assert!(impact.realigned());
}
