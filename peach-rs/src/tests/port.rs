//! Tests of the port: event ids, the chromosome names, and the corrected findings in both modes.

use std::collections::{BTreeMap, HashSet};

use crate::compat::with_java;
use crate::loaders::load_relevant_events;
use crate::panel::*;

fn temp(name: &str, text: &str) -> String {
    let dir = std::env::temp_dir().join(format!("peach_port_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_str().unwrap().to_string()
}

#[test]
fn event_ids() {
    let e = VariantHaplotypeEvent::from_id("VAR_chr1_97450058_C_T").unwrap();
    assert_eq!((e.chromosome.name(), e.position, e.ref_.as_str(), e.alt.as_str()), ("1", 97450058, "C", "T"));
    assert_eq!(e.id(), "VAR_chr1_97450058_C_T");
    assert!(VariantHaplotypeEvent::from_id("VAR_chr1_097450058_C_T").is_err());
    assert!(VariantHaplotypeEvent::from_id("VAR_chr1_97450058_C").is_err());
    assert!(VariantHaplotypeEvent::from_id("SV_chr1_97450058_C_T").is_err());
    assert!(VariantHaplotypeEvent::from_id("VAR_chrM_100_C_T").is_err());
    assert_eq!(Chromosome::from_string("CHRX").unwrap().name(), "X");
    assert!(Chromosome::from_string("chrx").is_err());
    assert_eq!(Chromosome::from_string("22").unwrap().name(), "22");
}

const HEADER: &str = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

fn positions() -> BTreeMap<Chromosome, HashSet<i32>> {
    let mut m = BTreeMap::new();
    m.insert(Chromosome::from_string("1").unwrap(), [1000].into_iter().collect());
    m
}

#[test]
fn pe1_records_that_cover_no_panel_position() {
    // a relevant SNV, then records Java cannot turn into events, none of them at a panel position
    let vcf = temp(
        "pe1.vcf",
        &format!(
            "{HEADER}1\t1000\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n1\t2000\t.\tA\tC,G\t.\tPASS\t.\tGT\t1/2\nchrM\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\n1\t3000\t.\tA\tG\t.\tPASS\t.\tGT\t./1\n1\t4000\t.\tA\t.\t.\tPASS\t.\tGT\t0/0\n"
        ),
    );
    let corrected = with_java(false, || load_relevant_events(&vcf, "S1", &positions())).unwrap();
    assert_eq!(corrected.len(), 1);
    assert_eq!(corrected["VAR_chr1_1000_C_T"], Some(1));
    assert!(with_java(true, || load_relevant_events(&vcf, "S1", &positions())).is_err());
    // each kind of record alone ends the Java run
    for rec in ["1\t2000\t.\tA\tC,G\t.\tPASS\t.\tGT\t1/2", "chrM\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1/1", "1\t3000\t.\tA\tG\t.\tPASS\t.\tGT\t./1", "1\t4000\t.\tA\t.\t.\tPASS\t.\tGT\t0/0"] {
        let v = temp("pe1_one.vcf", &format!("{HEADER}{rec}\n"));
        assert!(with_java(true, || load_relevant_events(&v, "S1", &positions())).is_err(), "{rec}");
        assert!(with_java(false, || load_relevant_events(&v, "S1", &positions())).unwrap().is_empty(), "{rec}");
    }
    // a sample that is not in the VCF
    assert!(with_java(false, || load_relevant_events(&vcf, "S2", &positions())).is_err());
}

#[test]
fn relevant_records_keep_java_checks() {
    // a multi-allelic record at a panel position ends the run in both modes
    let vcf = temp("rel.vcf", &format!("{HEADER}1\t1000\t.\tC\tT,G\t.\tPASS\t.\tGT\t1/2\n"));
    assert!(with_java(false, || load_relevant_events(&vcf, "S1", &positions())).is_err());
    // a deletion that covers the panel position from before it is relevant; a filtered record is not
    let vcf = temp("del.vcf", &format!("{HEADER}1\t998\t.\tTACG\tT\t.\tPASS\t.\tGT\t1/1\n1\t1000\t.\tC\tA\t.\tLOW_VAF\t.\tGT\t0/1\n"));
    for java in [true, false] {
        let e = with_java(java, || load_relevant_events(&vcf, "S1", &positions())).unwrap();
        assert_eq!(e.len(), 1);
        assert_eq!(e["VAR_chr1_998_TACG_T"], Some(2));
    }
}
