//! Tests of the port: the corrected findings in both modes, the forest and the prediction file.

use super::*;
use crate::compat::with_java;
use crate::predict::forest::{Forest, Node};
use crate::predict::{predict, predictions_text, ContextTable};
use crate::refgenome::RefGenomeSource;

/// A synthetic genome in a fresh directory: (path of the FASTA, its directory)
fn write_genome(name: &str, contigs: &[(&str, &str)], dict: bool) -> String {
    let dir = tmp_output_dir(name);
    let fasta = format!("{dir}genome.fasta");
    let mut text = String::new();
    let mut fai = String::new();
    let mut dict_text = String::from("@HD\tVN:1.6\n");
    for (name, seq) in contigs {
        text.push_str(&format!(">{name}\n"));
        let offset = text.len();
        for chunk in seq.as_bytes().chunks(60) {
            text.push_str(std::str::from_utf8(chunk).unwrap());
            text.push('\n');
        }
        fai.push_str(&format!("{name}\t{}\t{offset}\t60\t61\n", seq.len()));
        dict_text.push_str(&format!("@SQ\tSN:{name}\tLN:{}\n", seq.len()));
    }
    std::fs::write(&fasta, text).unwrap();
    std::fs::write(format!("{fasta}.fai"), fai).unwrap();
    if dict {
        std::fs::write(format!("{dir}genome.dict"), dict_text).unwrap();
    }
    fasta
}

fn write_vcf(dir: &str, name: &str, records: &[&str]) -> String {
    let path = format!("{dir}{name}");
    let mut text = String::from("##fileformat=VCFv4.2\n##INFO=<ID=MATEID,Number=1,Type=String,Description=\"Mate\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for r in records {
        text.push_str(r);
        text.push('\n');
    }
    std::fs::write(&path, text).unwrap();
    path
}

const SEQ: &str = "ACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTTAGCTAGGCTACGATCGATCGATGCATCGATCGTAGCTAGCATGCTAGCTAGCTAGCTGACTGATCGTAGCTAGCTAGCTGATCG";

fn snv_context_of(fasta: &str, pos: i32, r: &str, a: &str) -> anyhow::Result<String> {
    let g = RefGenomeSource::load(fasta).unwrap();
    let d = crate::snv::SnvDetails::from(&crate::vcf_file::SmallVariant::new("chr1", pos, r, a), &g)?;
    Ok(d.tri_nuc_context)
}

/// CH4: soft-masked (lower case) bases make Java's trinucleotide context invalid
#[test]
fn ch4_soft_masked_reference() {
    let lower: String = SEQ.to_ascii_lowercase();
    let fasta = write_genome("ch4", &[("chr1", &lower)], true);
    with_java(true, || assert_eq!(snv_context_of(&fasta, 10, "C", "A").unwrap(), "C>A_acg"));
    with_java(false, || assert_eq!(snv_context_of(&fasta, 10, "C", "A").unwrap(), "C>A_ACG"));
    // INDELs: the flanks are compared with the upper-case allele; without homology or repeats nothing changes
    let indel = crate::indel::IndelVariant::new(crate::vcf_file::SmallVariant::new("chr1", 22, "GGAT", "G"));
    let g = RefGenomeSource::load(&fasta).unwrap();
    with_java(true, || {
        let d = crate::indel::IndelDetails::from(&indel, &g).unwrap();
        assert_eq!((d.repeat_units_count, d.context_name().as_str()), (1, "del.none.len.3"));
    });
    with_java(false, || {
        let d = crate::indel::IndelDetails::from(&indel, &g).unwrap();
        assert_eq!((d.repeat_units_count, d.context_name().as_str()), (1, "del.none.len.3"));
    });
    // a deletion in a repeat (CTAG CTAG) is counted as one without context in Java
    let repeat = crate::indel::IndelVariant::new(crate::vcf_file::SmallVariant::new("chr1", 14, "GCTAG", "G"));
    with_java(true, || assert_eq!(crate::indel::IndelDetails::from(&repeat, &g).unwrap().context_name(), "del.none.len.4"));
    with_java(false, || assert_eq!(crate::indel::IndelDetails::from(&repeat, &g).unwrap().context_name(), "del.rep.len.4"));
}

/// CH5: an SNV at the first or last base of a contig, INDELs whose flanks reach past either end
#[test]
fn ch5_contig_edges() {
    let fasta = write_genome("ch5", &[("chr1", SEQ)], true);
    let last = SEQ.len() as i32;
    with_java(true, || {
        assert!(snv_context_of(&fasta, 1, "A", "C").is_err());
        assert!(snv_context_of(&fasta, last, "G", "C").is_err());
    });
    with_java(false, || {
        assert_eq!(snv_context_of(&fasta, 1, "A", "C").unwrap(), "T>G_GTN");
        assert_eq!(snv_context_of(&fasta, last, "G", "C").unwrap(), "C>G_NCG");
    });
    let g = RefGenomeSource::load(&fasta).unwrap();
    let start = crate::indel::IndelVariant::new(crate::vcf_file::SmallVariant::new("chr1", 2, "CGTT", "C"));
    let end = crate::indel::IndelVariant::new(crate::vcf_file::SmallVariant::new("chr1", last, "G", "GAAA"));
    with_java(true, || {
        assert!(crate::indel::IndelDetails::from(&start, &g).is_err());
        assert!(crate::indel::IndelDetails::from(&end, &g).is_err());
    });
    with_java(false, || {
        let d = crate::indel::IndelDetails::from(&start, &g).unwrap();
        assert_eq!(d.left_flank_bases, "AC");
        let d = crate::indel::IndelDetails::from(&end, &g).unwrap();
        assert_eq!((d.right_flank_bases.as_str(), d.context_name().as_str()), ("", "ins.none.len.3"));
    });
}

/// CH6: without a sequence dictionary Java has no chromosome lengths
#[test]
fn ch6_no_dictionary() {
    let fasta = write_genome("ch6", &[("chr1", SEQ)], false);
    let g = RefGenomeSource::load(&fasta).unwrap();
    with_java(true, || assert!(g.chromosome_length("chr1").is_err()));
    with_java(false, || assert_eq!(g.chromosome_length("chr1").unwrap(), SEQ.len() as i32));
    let with_dict = write_genome("ch6b", &[("chr1", SEQ)], true);
    let g = RefGenomeSource::load(&with_dict).unwrap();
    with_java(true, || assert_eq!(g.chromosome_length("chr1").unwrap(), SEQ.len() as i32));
}

/// CH8: a breakend on a contig whose name has ':' in it
#[test]
fn ch8_contig_names_with_colons() {
    let mk = || crate::sv::StructuralVariant::new("x", "HLA-A*01:01:01:01", 100, "A", "A[HLA-A*01:01:01:01:2100[").unwrap();
    with_java(true, || {
        let s = mk();
        assert_eq!((s.sv_type, s.alt_position), (hmf_common::sv::StructuralVariantType::BND, 1));
    });
    with_java(false, || {
        let s = mk();
        assert_eq!((s.sv_type, s.length), (hmf_common::sv::StructuralVariantType::DEL, Some(2000)));
    });
}

/// A forest of two stumps, one voting BRCA1 and one BRCA2: every sample ties at 0.5 each
fn tie_forest() -> Forest {
    let stump = |class| vec![Node { status: -1, var: 0, left: 0, right: 0, pred: class, split: 0.0 }];
    Forest { classes: vec!["BRCA1".into(), "BRCA2".into(), "none".into()], variables: vec!["C.A".into()], trees: vec![stump(1), stump(2)] }
}

fn table(ids: &[&str], indels: f64, svs: f64) -> ContextTable {
    let mut columns: Vec<String> = crate::snv::bucket_names().iter().map(|b| crate::snv::rename_tri_nuc_bin(b)).collect();
    columns.extend(crate::indel::context_names());
    columns.extend(crate::sv::context_names());
    let mut values = vec![0.0; columns.len()];
    values[0] = 10.0;
    values[96 + 20] = indels;
    values[96 + 30] = svs;
    ContextTable { columns, rows: ids.iter().map(|id| (id.to_string(), values.clone())).collect() }
}

/// CH1: ties of the BRCA1 and BRCA2 probabilities
#[test]
fn ch1_ties() {
    let f = tie_forest();
    let t = table(&["A", "B", "C", "D"], 150.0, 40.0);
    // R under set.seed(1): runif 0.266, 0.372, 0.573, 0.908: BRCA2, BRCA2, BRCA1, BRCA1
    with_java(true, || {
        let p = predict(&t, &f).unwrap();
        let types: Vec<&str> = p.iter().map(|x| x.hrd_type.as_str()).collect();
        assert_eq!(types, ["BRCA2_type", "BRCA2_type", "BRCA1_type", "BRCA1_type"]);
    });
    with_java(false, || {
        let p = predict(&t, &f).unwrap();
        assert!(p.iter().all(|x| x.hrd_type == "BRCA1_type" && x.hr_status == "HR_deficient"));
    });
}

/// CH3: the sample ids as R's read.delim converts them
#[test]
fn ch3_sample_ids() {
    let f = tie_forest();
    with_java(true, || {
        let p = predict(&table(&["0123"], 150.0, 40.0), &f).unwrap();
        assert_eq!(p[0].sample, "123");
        assert!(predict(&table(&["NA"], 150.0, 40.0), &f).is_err());
        assert!(predict(&table(&["S1", "S1"], 150.0, 40.0), &f).is_err());
    });
    with_java(false, || {
        assert_eq!(predict(&table(&["0123"], 150.0, 40.0), &f).unwrap()[0].sample, "0123");
        assert_eq!(predict(&table(&["NA"], 150.0, 40.0), &f).unwrap()[0].sample, "NA");
        assert_eq!(predict(&table(&["S1", "S1"], 150.0, 40.0), &f).unwrap().len(), 2);
    });
}

/// The QC rules and remarks
#[test]
fn qc_rules() {
    let f = tie_forest();
    with_java(false, || {
        let line = |indels: f64, svs: f64| predictions_text(&predict(&table(&["S"], indels, svs), &f).unwrap()).lines().nth(1).unwrap().to_string();
        assert_eq!(line(150.0, 40.0), "S\t0.5\t0.5\t1\tHR_deficient\tBRCA1_type\t\t");
        assert_eq!(line(150.0, 29.0), "S\t0.5\t0.5\t1\tHR_deficient\tBRCA1_type\t\t<30 SVs");
        assert_eq!(line(99.0, 29.0), "S\t0.5\t0.5\t1\tcannot_be_determined\tnone\t<100 indels\t<30 SVs");
        assert_eq!(line(0.0, 0.0), "S\t0.5\t0.5\t1\tcannot_be_determined\tnone\t<100 indels\t<30 SVs");
    });
    // MSI: more than 14,000 INDELs in repeats (column 0 of the INDELs is del.rep.len.1)
    let mut t = table(&["S"], 0.0, 40.0);
    t.rows[0].1[96] = 14001.0;
    with_java(false, || {
        let p = predict(&t, &f).unwrap();
        assert_eq!((p[0].hr_status.as_str(), p[0].remarks_hr_status.as_str()), ("cannot_be_determined", "Has MSI (>14000 indel.rep)"));
    });
}

/// The column groups must account for every column (groupColumns)
#[test]
fn unexpected_columns_fail() {
    let mut t = table(&["S"], 150.0, 40.0);
    t.columns[0] = "odd".into();
    assert!(predict(&t, &tie_forest()).is_err());
    assert!(crate::predict::r_grepl(r"^\w\[\w>\w\]\w$", "A[C>A]T"));
    assert!(!crate::predict::r_grepl(r"^\w\[\w>\w\]\w$", "A[C>A]TT"));
    assert!(crate::predict::r_grepl("del.mh.bimh.[2,3,4,5]", "del.mh.bimh.4"));
    assert!(!crate::predict::r_grepl("del.mh.bimh.[2,3,4,5]", "del.mh.bimh.1"));
    assert!(crate::predict::r_grepl("^[A-Z]{3}", "DUP_1e03_1e04_bp") && !crate::predict::r_grepl("^[A-Z]{3}", "del.rep.len.1"));
}

/// The embedded forest is CHORD.rds's (read here with the RDS reader, where the Java sources are at hand)
#[test]
fn embedded_forest_equals_chord_rds() {
    let embedded = Forest::embedded().unwrap();
    assert_eq!(embedded.trees.len(), 500);
    assert_eq!(embedded.classes, ["BRCA1", "BRCA2", "none"]);
    let candidates = [
        format!("{}/../chord/src/main/R/CHORD.rds", env!("CARGO_MANIFEST_DIR")),
        format!("{}/../../../../hmftools/chord/src/main/R/CHORD.rds", env!("CARGO_MANIFEST_DIR")),
    ];
    let Some(rds) = candidates.iter().find(|p| std::path::Path::new(p).exists()) else {
        eprintln!("CHORD.rds not found; skipped");
        return;
    };
    let from_rds = crate::predict::rds::forest_from_object(&crate::predict::rds::read_rds(rds).unwrap()).unwrap();
    assert_eq!(from_rds, embedded);
}

/// HG008 (public GIAB tumor): the counts of the pipeline run, and R's prediction from them
#[test]
fn hg008_prediction() {
    let counts = "107 128 10 88 157 90 13 218 109 98 13 71 227 116 19 266 49 25 7 44 37 27 9 47 22 43 5 24 364 86 17 444 109 57 188 90 176 455 147 381 96 89 211 115 873 254 143 460 69 41 40 123 104 98 118 155 21 22 20 23 74 34 36 63 97 41 59 66 46 55 46 82 61 43 30 58 63 41 33 82 39 17 29 46 22 19 28 31 14 11 18 22 44 17 38 64 295 41 25 11 19 442 40 14 13 14 67 46 29 14 16 3 3 8 0 5 85 43 26 11 68 25 14 11 2 2 43 15 8 1 3 1 9 1 10 39 7 0 4 7 1 1 0 5 26";
    let mut t = table(&["HG008-T"], 0.0, 0.0);
    t.rows[0].1 = counts.split(' ').map(|v| v.parse().unwrap()).collect();
    let text = with_java(true, || predictions_text(&predict(&t, &Forest::embedded().unwrap()).unwrap()));
    assert_eq!(text.lines().nth(1).unwrap(), "HG008-T\t0.004\t0.002\t0.006\tHR_proficient\tnone\t\t");
}

/// The text form of a forest is checked as it is read: here node 1 points to a node 3 that is missing
#[test]
fn forest_text_is_checked() {
    let bad = "type\tclassification\nclasses\tA\tB\nvariables\tx\nncat\t1\nntree\t1\ntree\t1\t2\n1\t1\t2\t3\t0\t0x1p-1\n-1\t0\t0\t0\t1\t0x0p+0\n";
    assert!(Forest::parse(bad).is_err());
    let good = bad.replace("1\t1\t2\t3\t0\t0x1p-1\n-1", "1\t1\t2\t3\t0\t0x1p-1\n-1\t0\t0\t0\t1\t0x0p+0\n-1").replace("tree\t1\t2", "tree\t1\t3");
    let f = Forest::parse(&good).unwrap();
    assert_eq!(f.votes(&[0.4]), [1, 0]);
    assert_eq!(f.votes(&[0.6]), [1, 0]);
}

/// The mutation contexts file of a multi-sample run: the header of the first sample, a line per sample
#[test]
fn mutation_contexts_text() {
    let rows = vec![("A".to_string(), vec![mc("x", 1), mc("y", 2)]), ("B".to_string(), vec![mc("x", 0), mc("y", 5)])];
    assert_eq!(crate::prep::mutation_contexts_text(&rows), "sample_id\tx\ty\nA\t1\t2\nB\t0\t5\n");
}

/// The sample wildcard, with the sample id as a Java replacement string
#[test]
fn wildcard_sample_paths() {
    use crate::config::convert_wildcard_sample_path as w;
    assert_eq!(w("/d/*.purple.sv.vcf.gz", "S1").unwrap(), "/d/S1.purple.sv.vcf.gz");
    assert_eq!(w("/d/x.vcf.gz", "S$1").unwrap(), "/d/x.vcf.gz");
    assert!(w("/d/*.vcf.gz", "S$1").is_err());
    assert_eq!(w("/d/*.vcf.gz", "S\\$1").unwrap(), "/d/S$1.vcf.gz");
    assert_eq!(w("/d/*.vcf.gz", "S$0").unwrap(), "/d/S*.vcf.gz");
}

/// The sample id file: a SampleId column, or a single column
#[test]
fn sample_id_file() {
    let dir = tmp_output_dir("ids");
    let p = format!("{dir}ids.csv");
    std::fs::write(&p, "SampleId,Other\nA,1\n#B,2\n\nC,3\n").unwrap();
    assert_eq!(crate::config::load_sample_ids_file(&p).unwrap(), ["A", "C"]);
    std::fs::write(&p, "A\nB\n").unwrap();
    assert_eq!(crate::config::load_sample_ids_file(&p).unwrap(), ["A", "B"]);
    std::fs::write(&p, "A,1\nB,2\n").unwrap();
    assert!(crate::config::load_sample_ids_file(&p).unwrap().is_empty());
}

/// PASS only, unless -include_non_pass; FILTER "." is not a pass
#[test]
fn pass_filter() {
    let dir = tmp_output_dir("pass");
    let p = write_vcf(&dir, "x.vcf", &["chr1\t5\t.\tA\tC\t.\tPASS\t.", "chr1\t6\t.\tA\tC\t.\tPON\t.", "chr1\t7\t.\tA\tC\t.\t.\t."]);
    assert_eq!(crate::vcf_file::load_variants(&p, false, "").unwrap().len(), 1);
    assert_eq!(crate::vcf_file::load_variants(&p, true, "").unwrap().len(), 3);
    let bad = write_vcf(&dir, "y.vcf", &["chr1\t5\t.\tR\tC\t.\tPASS\t."]);
    assert!(crate::vcf_file::load_variants(&bad, false, "").is_err());
}
