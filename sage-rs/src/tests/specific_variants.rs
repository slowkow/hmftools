//! Not from Java: the helpers behind -specific_var_file and -write_frag_lengths.

use crate::config::{add_specific_region, load_simple_variants, Config};
use crate::evidence::add_fragment_length;
use crate::phase::java_int_hash_order;
use crate::types::SimpleVariant;
use std::io::Write;

fn temp_file(name: &str, content: &str) -> String {
    let path = std::env::temp_dir().join(format!("sage_rs_test_{}_{name}", std::process::id()));
    std::fs::File::create(&path).unwrap().write_all(content.as_bytes()).unwrap();
    path.to_string_lossy().to_string()
}

#[test]
fn test_add_specific_region() {
    let mut chromosomes = Vec::new();
    let mut regions = Vec::new();
    add_specific_region(&mut chromosomes, &mut regions, "1", 100, 500);
    add_specific_region(&mut chromosomes, &mut regions, "1", 1000, 1400);
    add_specific_region(&mut chromosomes, &mut regions, "2", 300, 700);
    assert_eq!(chromosomes, vec!["1", "2"]);
    assert_eq!(regions.len(), 3);

    // an overlap widens the first region it meets
    add_specific_region(&mut chromosomes, &mut regions, "1", 400, 800);
    assert_eq!(regions[0], ("1".to_string(), 100, 800));
    assert_eq!(regions.len(), 3);

    // as in SpecificRegions.addRegion, regions which come to overlap each other stay separate
    add_specific_region(&mut chromosomes, &mut regions, "1", 700, 1100);
    assert_eq!(regions[0], ("1".to_string(), 100, 1100));
    assert_eq!(regions[1], ("1".to_string(), 1000, 1400));

    // the same range on another chromosome is not an overlap
    add_specific_region(&mut chromosomes, &mut regions, "3", 100, 500);
    assert_eq!(regions.len(), 4);
    assert_eq!(chromosomes, vec!["1", "2", "3"]);
}

#[test]
fn test_load_simple_variants() {
    // tab-delimited, columns in any order
    let tsv = temp_file("vars.tsv", "Gene\tPosition\tChromosome\tAlt\tRef\nTP53\t150\tchr1\tT\tA\n\t200\tchr2\tG\tGAT\n");
    let variants = load_simple_variants(&tsv).unwrap();
    assert_eq!(variants.len(), 2);
    assert_eq!(variants[0], ("chr1".to_string(), SimpleVariant::new(150, b"A", b"T")));
    assert_eq!(variants[1], ("chr2".to_string(), SimpleVariant::new(200, b"GAT", b"G")));

    // a header with a comma is read as comma-delimited
    let csv = temp_file("vars.csv", "Chromosome,Position,Ref,Alt\n7,55,C,CA\n");
    let variants = load_simple_variants(&csv).unwrap();
    assert_eq!(variants, vec![("7".to_string(), SimpleVariant::new(55, b"C", b"CA"))]);

    let bad = temp_file("bad.tsv", "Chromosome\tPosition\tRef\nchr1\t1\tA\n");
    assert!(load_simple_variants(&bad).is_err());
    for f in [&tsv, &csv, &bad] {
        let _ = std::fs::remove_file(f);
    }
}

fn config_with(extra: &[&str]) -> Config {
    let mut args: Vec<String> =
        ["-tumor", "TUMOR_ID", "-ref_genome", "none.fasta", "-ref_genome_version", "37", "-output_vcf", "none.vcf.gz"]
            .iter()
            .map(|x| x.to_string())
            .collect();
    args.extend(extra.iter().map(|x| x.to_string()));
    Config::from_args(&args).unwrap()
}

#[test]
fn test_specific_var_file_config() {
    let file = temp_file("cfg_vars.tsv", "Chromosome\tPosition\tRef\tAlt\nchr1\t1000\tA\tT\nchr1\t1300\tC\tG\nchr2\t50\tG\tA\n");
    // each variant adds 200 bases either side; the second overlaps the first
    let cfg = config_with(&["-specific_var_file", &file]);
    assert_eq!(cfg.specific_variants.len(), 3);
    assert_eq!(cfg.specific_chr, vec!["chr1", "chr2"]);
    assert_eq!(cfg.specific_regions, vec![("chr1".to_string(), 800, 1500), ("chr2".to_string(), -150, 250)]);
    assert!(cfg.process_chromosome("chr1"));
    assert!(!cfg.process_chromosome("chr3"));

    // SageConfig ignores the file when a region or a chromosome is given
    let cfg = config_with(&["-specific_var_file", &file, "-specific_chr", "chr1"]);
    assert!(cfg.specific_variants.is_empty());
    assert!(cfg.specific_regions.is_empty());
    let cfg = config_with(&["-specific_var_file", &file, "-specific_regions", "chr2:1-100"]);
    assert!(cfg.specific_variants.is_empty());
    assert_eq!(cfg.specific_regions, vec![("chr2".to_string(), 1, 100)]);
    let _ = std::fs::remove_file(&file);
}

#[test]
fn test_fragment_length_counts() {
    let mut lengths = Vec::new();
    for (length, is_alt) in [(300, false), (150, true), (300, true), (0, true), (501, false), (-20, false), (500, false), (1, true), (300, false)] {
        add_fragment_length(&mut lengths, length, is_alt);
    }
    // lengths outside 1 to 500 are dropped; counts are [ref, alt]
    assert_eq!(lengths, vec![(300, [2, 1]), (150, [0, 1]), (500, [1, 0]), (1, [0, 1])]);

    // written in the order of Java's HashMap<Integer, int[]>: by bucket, then by insertion
    let keys: Vec<i32> = lengths.iter().map(|(l, _)| *l).collect();
    assert_eq!(java_int_hash_order(&keys), vec![1, 500, 150, 300]);

    // enough keys to grow the table to 64 buckets: 70 and 6 share a bucket and keep insertion order
    let keys: Vec<i32> = (60..=95).rev().chain([6]).collect();
    let order = java_int_hash_order(&keys);
    assert_eq!(&order[..4], &[64, 65, 66, 67]);
    let at = |k: i32| order.iter().position(|x| *x == k).unwrap();
    assert_eq!(at(6), at(70) + 1);
    assert_eq!(*order.last().unwrap(), 63);
}
