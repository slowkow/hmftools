//! Tests of the port: the summary reader, the output format, and the corrected findings.

use std::collections::HashSet;

use crate::files::*;

fn temp(name: &str, text: &str) -> String {
    let dir = std::env::temp_dir().join(format!("vi_port_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_str().unwrap().to_string()
}

const SUMMARY_HEADER: &str = "taxid_genus\tname_genus\treads_genus_tree\ttaxid_species\tname_species\treads_species_tree\ttaxid_assigned\tname_assigned\treads_assigned_tree\treads_assigned_direct\treference\treference_taxid\treference_kmer_count\talternate_kmer_count\t#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tintegrations\tQCStatus";

#[test]
fn summary_lines_with_and_without_a_qc_status() {
    let row = "10375\tLymphocryptovirus\t100\t10376\tHuman gammaherpesvirus 4\t100\t10376\tHuman gammaherpesvirus 4\t100\t100\tadjusted_taxid|10376|NC_007605.1\t10376\t10\t0\tadjusted_taxid|10376|NC_007605.1\t1\t171823\t1000\t160000\t93.1\t40.5\t35.2\t50.1\t2";
    let path = temp("summary.tsv", &format!("{SUMMARY_HEADER}\n{row}\t\n{row}\tLOW_VIRAL_COVERAGE\n"));
    let b = read_virus_breakends(&path).unwrap();
    assert_eq!(b.len(), 2);
    assert_eq!(b[0].qc_status, "NO_ABNORMALITIES");
    assert_eq!(b[1].qc_status, "LOW_VIRAL_COVERAGE");
    assert_eq!((b[0].taxid_genus, b[0].taxid_species, b[0].reference_taxid, b[0].integrations), (10375, 10376, 10376, 2));
    assert_eq!((b[0].coverage, b[0].mean_depth), (93.1, 40.5));
    // an empty file holds no viruses (VIRUSBreakend writes one when it finds none)
    assert!(read_virus_breakends(&temp("empty.tsv", "")).unwrap().is_empty());
    // a short line ends the run, as Java's ArrayIndexOutOfBoundsException
    assert!(read_virus_breakends(&temp("short.tsv", &format!("{SUMMARY_HEADER}\n1\t2\n"))).is_err());
    let bad = row.replacen("\t2", "\tX", 1);
    assert!(read_virus_breakends(&temp("bad.tsv", &format!("{SUMMARY_HEADER}\n{bad}\n"))).is_err());
}

#[test]
fn taxonomy_keeps_the_wanted_names_and_the_last_of_a_taxid() {
    let path = temp("tax.tsv", "1\t root \n2\tBacteria\n2\tBacteria again\nbad line\n");
    let wanted: HashSet<i32> = [2].into_iter().collect();
    let db = load_taxonomy_db(&path, Some(&wanted)).unwrap();
    assert_eq!(db.lookup_name(2), "Bacteria again");
    assert!(!db.taxid_exists(1));
    let all = load_taxonomy_db(&path, None).unwrap();
    assert_eq!(all.count, 2);
    assert_eq!(all.lookup_name(1), "root");
    // String.split drops trailing empty fields: a name followed by tabs is still two fields
    let t = load_taxonomy_db(&temp("tax_tabs.tsv", "5\tname\t\t\n6\t\t\n7\ta\tb\n"), None).unwrap();
    assert_eq!(t.count, 1);
    assert_eq!(t.lookup_name(5), "name");
    assert!(load_taxonomy_db(&temp("tax_utf8.tsv", "5\tna\u{e9}me\n"), None).is_ok());
    // Files.readAllLines fails on bytes that are not UTF-8
    let bin = temp("tax_bin.tsv", "");
    std::fs::write(&bin, [b'5', b'\t', 0xff, b'\n']).unwrap();
    assert!(load_taxonomy_db(&bin, None).is_err());
    assert!(load_taxonomy_db(&temp("tax_bad.tsv", "x\tname\n"), None).is_err());
}

#[test]
fn annotated_virus_format() {
    let v = AnnotatedVirus {
        taxid: 10376,
        name: "Human gammaherpesvirus 4".into(),
        qc_status: "NO_ABNORMALITIES",
        integrations: 2,
        interpretation: Some("EBV"),
        percentage_covered: 93.0,
        mean_coverage: 40.5,
        expected_clonal_coverage: Some(31.935483870967744),
        reported: true,
        blacklisted: Some(false),
        virus_driver_likelihood_type: VirusLikelihoodType::High,
    };
    let mut none = v.clone();
    none.interpretation = None;
    none.expected_clonal_coverage = None;
    none.percentage_covered = 0.0001;
    assert_eq!(
        annotated_virus_text(&[v, none]),
        format!(
            "{ANNOTATED_VIRUS_HEADER}\n10376\tHuman gammaherpesvirus 4\tNO_ABNORMALITIES\t2\tEBV\t93.0\t40.5\t31.935483870967744\ttrue\tfalse\tHIGH\n10376\tHuman gammaherpesvirus 4\tNO_ABNORMALITIES\t2\tnull\t1.0E-4\t40.5\tnull\ttrue\tfalse\tHIGH\n"
        )
    );
}
