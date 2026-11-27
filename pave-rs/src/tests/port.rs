//! Tests of the port: the corrected Java faults (both behaviours), the resource readers, the output.

use super::*;
use hmf_common::vcf::VcfRecord;

use crate::annotation::gnomad::GnomadData;
use crate::annotation::NeededPositions;
use crate::compat::with_java;
use crate::variant::VariantData;

fn temp_dir(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("pave_rs_{name}_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn needed(chr: &str, positions: &[i32]) -> NeededPositions {
    let mut n = NeededPositions::default();
    n.by_chr.insert(chr.to_string(), positions.to_vec());
    n.finish();
    n
}

/// PA1: a deletion reaching past the end of a gene no longer drops the gene for itself and later variants
#[test]
fn gene_search_keeps_genes_a_deletion_overlaps() {
    use super::homology_realign::create_ensembl_gene_data;
    let mut cache = crate::gene_cache::GeneDataCache::empty(hmf_common::regions::RefGenomeVersion::V38);
    let genes = vec![
        create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, "chr1", POS_STRAND, 1100, 1500),
        create_ensembl_gene_data(GENE_ID_2, GENE_NAME_2, "chr1", POS_STRAND, 5000, 6000),
    ];
    cache.ensembl_mut().chr_gene_data_map_mut().insert("chr1".to_string(), genes);
    let cache: &'static crate::gene_cache::GeneDataCache = Box::leak(Box::new(cache));
    let names = |v: Vec<&hmf_common::ensembl::GeneData>| v.iter().map(|g| g.gene_name.clone()).collect::<Vec<_>>();
    for java in [true, false] {
        with_java(java, || {
            let mut ix = cache.create_indexing("chr1");
            assert_eq!(names(cache.find_genes("chr1", 1200, 1200, Some(&mut ix))), vec![GENE_NAME_1]);
            // a DEL from 1450 to 1600 overlaps the gene's end
            let found = names(cache.find_genes("chr1", 1450, 1600, Some(&mut ix)));
            let later = names(cache.find_genes("chr1", 1480, 1480, Some(&mut ix)));
            if java {
                assert!(found.is_empty() && later.is_empty());
            } else {
                assert_eq!(found, vec![GENE_NAME_1]);
                assert_eq!(later, vec![GENE_NAME_1]);
            }
            // without indexing every overlapping gene is found
            assert_eq!(names(cache.find_genes("chr1", 1450, 1600, None)), vec![GENE_NAME_1]);
        });
    }
}

/// PA2: the mappability file's positions are 1-based; a variant at a listed position takes its value
#[test]
fn mappability_positions() {
    let dir = temp_dir("map");
    let f = dir.join("map.bed");
    std::fs::write(&f, "chr1\t1\t1\t0.0000\nchr1\t100\t100\t0.5000\nchr1\t200\t200\t0.9000\nchr2\t1\t1\t1.0000\n").unwrap();
    let n = needed("chr1", &[50, 100, 150, 200, 250]);
    for java in [true, false] {
        with_java(java, || {
            let m = crate::annotation::mappability::load(f.to_str().unwrap(), &n, java).unwrap();
            assert!(m.has_chromosome("chr1") && m.has_chromosome("chr2") && !m.has_chromosome("chr3"));
            assert_eq!(m.value("chr1", 50), Some(0.0));
            assert_eq!(m.value("chr1", 150), Some(0.5));
            assert_eq!(m.value("chr1", 250), None); // past the last entry
            if java {
                assert_eq!(m.value("chr1", 100), Some(0.0));
                assert_eq!(m.value("chr1", 200), Some(0.5));
            } else {
                assert_eq!(m.value("chr1", 100), Some(0.5));
                assert_eq!(m.value("chr1", 200), Some(0.9));
            }
        });
    }
    // a line of another field count ends the run, as Java exits
    std::fs::write(&f, "chr1\t1\t1\n").unwrap();
    assert!(crate::annotation::mappability::load(f.to_str().unwrap(), &n, false).is_err());
    std::fs::remove_dir_all(&dir).ok();
}

/// PA6: an MNV's gnomAD frequency is the lowest of its bases'
#[test]
fn gnomad_mnv_frequency() {
    let mut g = GnomadData::default();
    g.add_entry("chr1", 100, "A", "C", 0.0);
    g.add_entry("chr1", 101, "G", "T", 0.5);
    g.add_entry("chr1", 102, "T", "A", 0.2);
    assert_eq!(with_java(true, || g.frequency("chr1", true, "AGT", "CTA", 100)), Some(0.2));
    assert_eq!(with_java(false, || g.frequency("chr1", true, "AGT", "CTA", 100)), Some(0.0));
    assert_eq!(g.frequency("chr1", true, "AG", "CA", 100), None);
    assert_eq!(g.frequency("chr1", false, "G", "T", 101), Some(0.5));
}

/// the PON reader keeps the entries at the needed positions, the first of repeated ones
#[test]
fn pon_reader() {
    let dir = temp_dir("pon");
    let f = dir.join("pon.tsv");
    std::fs::write(
        &f,
        "Chromosome\tPosition\tRef\tAlt\tSampleCount\tMaxSampleReadCount\tTotalReadCount\tMultiPonStatus\n\
         chr1\t10\tA\tC\t3\t4\t9\tBASE\nchr1\t20\tA\tC\t5\t6\t12\tARTEFACT\nchr1\t20\tA\tC\t9\t9\t9\tBASE\nchr2\t5\tG\tT\t1\t1\t1\tMULTI\n",
    )
    .unwrap();
    let n = needed("chr1", &[20]);
    let p = crate::annotation::pon::load_pon(f.to_str().unwrap(), Some(&n), true).unwrap();
    assert!(p.get("chr1", 10, "A", "C").is_none());
    let e = p.get("chr1", 20, "A", "C").unwrap();
    assert_eq!((e.samples, e.max_sample_reads, e.mean_read_count()), (5, 6, 2));
    assert_eq!(e.multi_pon, crate::variant::MultiPonStatus::Artefact);
    assert!(p.has_chromosome("chr2") && p.get("chr2", 5, "G", "T").is_none());
    std::fs::remove_dir_all(&dir).ok();
}

/// buildVariant: PAVE's filters replace PASS, sorted with the record's own; PASS when there are none
#[test]
fn output_filters_and_annotations() {
    let mut v = VariantData::new("chr1", 100, "A", "C");
    let mut rec = VcfRecord::new("chr1", 100, ".", "A", "C");
    rec.filters = vec!["minTumorQual".to_string()];
    v.set_context(rec.clone());
    v.add_filter("PON");
    v.set_gnomad_frequency(Some(0.0012));
    v.set_pon_frequency(3, 7, 2, crate::variant::MultiPonStatus::Base);
    let mut out = rec.clone();
    crate::output::finish_record(&mut out, &v, None, false);
    assert_eq!(out.filters, vec!["minTumorQual".to_string(), "PON".to_string()]);
    assert_eq!(out.attr_string("GND_FREQ", ""), "1.200e-03");
    assert_eq!(out.attr_int("PON_COUNT", 0), 3);
    assert_eq!(out.attr_int("PON_MAX", 0), 7);
    // without filters: PASS, and no impacts written
    let mut v2 = VariantData::new("chr1", 100, "A", "C");
    v2.set_context(VcfRecord::new("chr1", 100, ".", "A", "C"));
    let mut out2 = v2.context().clone();
    crate::output::finish_record(&mut out2, &v2, None, false);
    assert_eq!(out2.filters, vec!["PASS".to_string()]);
    assert!(!out2.has_attr("PAVE_TI") && !out2.has_attr("IMPACT"));
}

/// Pathogenicity.fromClinvarAnnotation
#[test]
fn clinvar_pathogenicity() {
    use crate::task::{pathogenicity_from_clinvar as p, Pathogenicity::*};
    assert_eq!(p("", ""), Unknown);
    assert_eq!(p("Pathogenic/Likely_pathogenic", ""), Pathogenic);
    assert_eq!(p("Likely_pathogenic", ""), LikelyPathogenic);
    assert_eq!(p("Benign", ""), Benign);
    assert_eq!(p("Conflicting_classifications_of_pathogenicity", "Pathogenic(1)|Uncertain_significance(2)"), Pathogenic);
    assert_eq!(p("Conflicting_classifications_of_pathogenicity", "Pathogenic(1)|Benign(1)"), Conflicting);
    assert_eq!(p("Uncertain_significance", ""), Unknown);
}
