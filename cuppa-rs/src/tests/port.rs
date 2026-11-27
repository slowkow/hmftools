//! Tests of the port: numpy's arithmetic against numbers numpy wrote (testdata/numpy_fixture.tsv, from
//! scripts/numpy_fixture.py), the model resource, HG008's prediction in both modes, and the findings.

use std::sync::OnceLock;

use super::{temp_dir, testdata};
use crate::compat;
use crate::predict::model::{Model, Resource, EMBEDDED_MODEL};
use crate::predict::{numpy, pyfmt};

fn hex(s: &str) -> f64 {
    // float.hex(): [-]0x1.<hex>p<exp>, 0x0.0p+0, inf, nan
    let (neg, t) = match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s),
    };
    let v = match t {
        "inf" => f64::INFINITY,
        "nan" => f64::NAN,
        _ => {
            let t = t.strip_prefix("0x").unwrap();
            let (mant, exp) = t.split_once('p').unwrap();
            let (int, frac) = mant.split_once('.').unwrap_or((mant, ""));
            let mut bits: u64 = u64::from_str_radix(int, 16).unwrap();
            let mut e: i32 = exp.parse().unwrap();
            for c in frac.chars() {
                bits = (bits << 4) | c.to_digit(16).unwrap() as u64;
                e -= 4;
            }
            if bits == 0 {
                0.0
            } else {
                // exact: bits has at most 53 significant bits; scale by powers of two step by step
                let mut v = bits as f64;
                while e > 0 {
                    let s = e.min(1000);
                    v *= 2f64.powi(s);
                    e -= s;
                }
                while e < 0 {
                    let s = (-e).min(1000);
                    v /= 2f64.powi(s);
                    e += s;
                }
                v
            }
        }
    };
    if neg { -v } else { v }
}

fn list(s: &str) -> Vec<f64> {
    if s.is_empty() { Vec::new() } else { s.split(',').map(hex).collect() }
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

#[test]
fn numpy_fixture() {
    let text = include_str!("../../testdata/numpy_fixture.tsv");
    let mut n = 0;
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "sum" => assert!(same(numpy::pairwise_sum(&list(f[1])), hex(f[2])), "{line:.80}"),
            "dot" => assert!(same(numpy::dot(&list(f[1]), &list(f[2])), hex(f[3])), "dot {}", list(f[1]).len()),
            "dot_n" => assert!(same(numpy::dot_n(&list(f[1]), &list(f[2])), hex(f[3])), "dot_n {}", list(f[1]).len()),
            "gemm" => assert!(same(numpy::gemm_dot(&list(f[1]), &list(f[2])), hex(f[3])), "gemm {}", list(f[1]).len()),
            "interp" => assert!(same(numpy::interp(hex(f[1]), &list(f[2]), &list(f[3])), hex(f[4])), "{line}"),
            "argsort" => {
                let want: Vec<usize> = f[2].split(',').map(|i| i.parse().unwrap()).collect();
                assert_eq!(numpy::argsort(&list(f[1])), want, "{line}");
            }
            "repr" => assert_eq!(pyfmt::repr(hex(f[1])), f[2]),
            "round1" => assert!(same(numpy::round(hex(f[1]), 1), hex(f[2])), "{line}"),
            "signif2" => {
                let x = hex(f[1]);
                let m = 10f64.powf(1.0 - x.abs().log10().floor());
                assert!(same(numpy::round(x * m, 0) / m, hex(f[2])), "{line}");
            }
            k => panic!("unknown fixture kind {k}"),
        }
        n += 1;
    }
    assert!(n > 900);
}

pub fn model() -> &'static Model {
    static M: OnceLock<Model> = OnceLock::new();
    M.get_or_init(|| Model::load(&Resource::parse(EMBEDDED_MODEL).unwrap()).unwrap())
}

#[test]
fn embedded_model() {
    let m = model();
    assert_eq!(m.classes.len(), 40);
    assert_eq!(m.req_gen_pos.len(), 6188);
    assert_eq!((m.req_snv96.len(), m.req_event.len(), m.req_sig.len()), (96, 336, 5));
    assert_eq!((m.req_gene_exp.len(), m.req_alt_sj.len()), (5927, 68533));
    assert_eq!(m.rna_combined.lr.classes.len(), 33);
    assert!(m.alt_sj.lr.f32 && !m.gene_exp.lr.f32 && !m.dna_combined.lr.f32);
    assert_eq!(m.prob_floor, 0.01);
    assert_eq!(m.source_sha256.len(), 64);
    // the resource is the one recorded in REPORT.md
    assert_eq!(m.source_sha256, "3e0f9490a9f37ac252b23b38b6b991245b72bcdf8f963b2a91ade82547c9fe70");
}

fn predict_hg008(java: bool, group: &str) -> (String, String, String) {
    let out = temp_dir(&format!("hg008-{java}-{group}"));
    let args: Vec<String> = [
        "--sample_id",
        "HG008-T",
        "--features_path",
        &testdata("hg008/HG008-T.cuppa_data.tsv.gz"),
        "--classifier_path",
        &format!("{}/resources/cuppa_model.38.bin.gz", env!("CARGO_MANIFEST_DIR")),
        "--output_dir",
        &out,
        "--clf_group",
        group,
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(compat::with_java(java, || crate::predict::main(&args)), 0);
    let r = (
        std::fs::read_to_string(format!("{out}HG008-T.cuppa.pred_summ.tsv")).unwrap(),
        std::fs::read_to_string(format!("{out}HG008-T.cuppa.vis_data.tsv")).unwrap_or_default(),
        out.clone(),
    );
    let _ = std::fs::remove_dir_all(&out);
    r
}

/// HG008's prediction equals pycuppa's (macOS arm64 reference), in both modes (no finding applies)
#[test]
fn hg008_prediction_equals_pycuppa() {
    let want_summ = std::fs::read_to_string(testdata("hg008/HG008-T.cuppa.pred_summ.tsv")).unwrap();
    let want_vis = std::fs::read_to_string(testdata("hg008/HG008-T.cuppa.vis_data.tsv")).unwrap();
    for java in [true, false] {
        let (summ, vis, _) = predict_hg008(java, "DNA");
        assert_eq!(summ, want_summ, "java mode {java}");
        assert_eq!(vis, want_vis, "java mode {java}");
    }
}

/// CU2: without the sex feature, pycuppa treats the sample as male (female classes zeroed)
#[test]
fn cu2_missing_sex_feature() {
    use crate::predict::features::{fill_missing, Features};
    let m = model();
    let names: Vec<String> = m.req_snv96.clone();
    let f = Features { samples: vec!["S".into()], names: names.clone(), values: vec![vec![10.0; names.len()]] };
    let x = fill_missing(&f, m);
    assert!(x.sex_missing[0]);
    let ovary = m.dna_combined.lr.classes.iter().position(|c| c.contains("Ovary")).unwrap();
    let prostate = m.dna_combined.lr.classes.iter().position(|c| c == "Prostate").unwrap();
    let java = compat::with_java(true, || crate::predict::classifier::predict(m, &x).unwrap());
    assert_eq!(java[0].dna_combined[ovary], 0.0, "Java mode: treated as male");
    let fixed = compat::with_java(false, || crate::predict::classifier::predict(m, &x).unwrap());
    assert!(fixed[0].dna_combined[ovary] > 0.0 && fixed[0].dna_combined[prostate] > 0.0, "corrected: unknown sex");
}

/// CU7: a sample's numbers in a cohort equal its numbers alone (corrected), or pycuppa's cohort arithmetic
#[test]
fn cu7_cohort_independence() {
    use crate::predict::features::{fill_missing, load};
    let m = model();
    let f = load(&testdata("hg008/HG008-T.cuppa_data.tsv.gz"), Some("HG008-T")).unwrap();
    let mut two = load(&testdata("hg008/HG008-T.cuppa_data.tsv.gz"), Some("HG008-T")).unwrap();
    two.samples.push("COPY".into());
    let v = two.values[0].iter().map(|x| x * 2.0).collect();
    two.values.push(v);
    let alone = fill_missing(&f, m);
    let both = fill_missing(&two, m);
    let a = compat::with_java(false, || crate::predict::classifier::predict(m, &alone).unwrap());
    let b = compat::with_java(false, || crate::predict::classifier::predict(m, &both).unwrap());
    assert_eq!(a[0].dna_combined, b[0].dna_combined);
    assert_eq!(a[0].gen_pos, b[0].gen_pos);
    let j = compat::with_java(true, || crate::predict::classifier::predict(m, &both).unwrap());
    assert_ne!(a[0].gen_pos, j[0].gen_pos, "pycuppa's cohort arithmetic differs in the last bits");
}

/// CU6: a GERMLINE_AMP row of the driver catalog
#[test]
fn cu6_germline_amp_is_germline() {
    let dir = temp_dir("cu6");
    let data = testdata("java/pipeline_output/TUMOR_SAMPLE");
    for f in std::fs::read_dir(&data).unwrap() {
        let p = f.unwrap().path();
        std::fs::copy(&p, format!("{dir}{}", p.file_name().unwrap().to_string_lossy())).unwrap();
    }
    let cat = format!("{dir}TUMOR_SAMPLE.linx.driver.catalog.tsv");
    let text = std::fs::read_to_string(&cat).unwrap();
    let mut lines: Vec<String> = text.lines().map(String::from).collect();
    let h: Vec<&str> = lines[0].split('\t').collect();
    let mut row: Vec<String> = lines[1].split('\t').map(String::from).collect();
    row[h.iter().position(|c| *c == "gene").unwrap()] = "ERBB2".into();
    row[h.iter().position(|c| *c == "driver").unwrap()] = "GERMLINE_AMP".into();
    lines.push(row.join("\t"));
    std::fs::write(&cat, lines.join("\n") + "\n").unwrap();
    let c = super::with_sample_data_dir(super::test_config(&["TUMOR_SAMPLE"], None), &dir);
    let java = compat::with_java(true, || crate::prep::drivers::extract(&c, "TUMOR_SAMPLE").unwrap());
    let fixed = compat::with_java(false, || crate::prep::drivers::extract(&c, "TUMOR_SAMPLE").unwrap());
    assert!(java.iter().any(|d| d.index.key == "ERBB2.mut"));
    assert!(!fixed.iter().any(|d| d.index.key.starts_with("ERBB2")));
    let _ = std::fs::remove_dir_all(&dir);
}

/// CU3: the generic somatic variants file with no type filter (the indel drivers of an MSS sample)
#[test]
fn cu3_generic_file_without_types() {
    let c = super::test_config(&["TUMOR_SAMPLE"], None);
    let java = compat::with_java(true, || crate::prep::somatic::load_variants(&c, "TUMOR_SAMPLE", None));
    assert!(java.is_err(), "Java: NullPointerException");
    let fixed = compat::with_java(false, || crate::prep::somatic::load_variants(&c, "TUMOR_SAMPLE", None)).unwrap();
    let snvs = compat::with_java(false, || crate::prep::somatic::load_variants(&c, "TUMOR_SAMPLE", Some(&[crate::prep::somatic::VariantType::Snp]))).unwrap();
    assert!(fixed.len() >= snvs.len());
}

/// CU4: a single sample whose RNA file is missing
#[test]
fn cu4_single_sample_missing_rna() {
    let out = temp_dir("cu4");
    let mut c = super::with_sample_data_dir(super::test_config(&["MINIMAL_SAMPLE_NO_RNA"], None), &format!("{}*", testdata("java/pipeline_output/")));
    c.categories = vec![crate::prep::config::CategoryType::Snv, crate::prep::config::CategoryType::GeneExp];
    c.output_dir = out.clone();
    assert!(compat::with_java(true, || crate::prep::run(&c)).is_err());
    compat::with_java(false, || crate::prep::run(&c)).unwrap();
    let items = crate::prep::items::read_item_list(&crate::prep::output_path(&c, None)).unwrap();
    assert!(items.iter().all(|d| d.index.source == crate::prep::items::Source::Dna));
    let _ = std::fs::remove_dir_all(&out);
}

/// CU8: alt-SJ counts are log-transformed before the classifier (corrected) or taken as they are
#[test]
fn cu8_alt_sj_counts() {
    use crate::predict::features::{fill_missing, Features};
    let m = model();
    let mut names = m.req_alt_sj[..3].to_vec();
    names.push(m.req_gene_exp[0].clone());
    let f = Features { samples: vec!["S".into()], names, values: vec![vec![7.0, 0.0, 1.0, 2.5]] };
    let java = compat::with_java(true, || fill_missing(&f, m));
    let fixed = compat::with_java(false, || fill_missing(&f, m));
    let sj = |x: &crate::predict::features::Prepared| x.alt_sj(0).unwrap()[..3].to_vec();
    assert_eq!(sj(&java), vec![7.0, 0.0, 1.0]);
    assert_eq!(sj(&fixed), vec![7f64.ln_1p(), 0.0, 2f64.ln()]);
    assert_eq!(java.gene_exp(0).unwrap()[0], fixed.gene_exp(0).unwrap()[0]);
}

#[test]
fn sha256_of_the_embedded_resource_is_stable() {
    // the digest function on a known input; the pickle itself (CUPPA_PICKLE=<cuppa_classifier.38.pickle.gz>) is
    // checked when given, and recognised by --classifier_path
    assert_eq!(crate::sha256::hex_digest(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    if let Ok(pickle) = std::env::var("CUPPA_PICKLE") {
        let bytes = std::fs::read(&pickle).unwrap();
        assert_eq!(crate::sha256::hex_digest(&bytes), model().source_sha256);
        assert!(crate::predict::model::load_for_path(&pickle).is_ok());
    }
}
