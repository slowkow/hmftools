//! NEO's TpmCalcsTest, ported, and tests of the scorer's own: the findings (NE6 to NE11), the binding
//! classes, the writers' formats and the readers.

use super::*;
use crate::bind::*;
use crate::compat;
use crate::score::*;
use hmf_common::neo::NeoEpitopeType;

const TEST_COHORT_CANCER_TPM: [f64; 2] = [20.0, 20.0];
const TEST_COHORT_PAN_CANCER_TPM: [f64; 2] = [30.0, 30.0];

const TRANS_NAME_1: &str = "TRANS_01";
const TRANS_NAME_2: &str = "TRANS_02";
const TRANS_NAME_3: &str = "TRANS_03";
const TRANS_NAME_4: &str = "TRANS_04";

/// TpmCalcsTest.createNeoepitope
fn create_neoepitope(id: i32, variant_info: &str, ty: NeoEpitopeType, rna_fragments: i32, tpms: [f64; 2], up_aas: &str, down_aas: &str) -> NeoEpitopeData {
    let mut n = NeoEpitopeData::new(
        id,
        ty,
        variant_info,
        GENE_NAME_1,
        up_aas,
        "",
        down_aas,
        Vec::new(),
        Vec::new(),
        [0, 0],
        1.0,
        1.0,
        0.0,
        [0, 0],
        0,
        [0, 0],
    );
    n.rna.set_coverage(rna_fragments, rna_fragments * 2, rna_fragments * 2);
    n.rna.set_expression(Some(tpms));
    n.rna.set_cohort_values(TEST_COHORT_CANCER_TPM, TEST_COHORT_PAN_CANCER_TPM);
    n
}

fn get_peptide<'a>(peptide: &str, peptides: &'a [PeptideScoreData]) -> &'a PeptideScoreData {
    peptides.iter().find(|p| p.peptide == peptide).unwrap_or_else(|| panic!("peptide {peptide}"))
}

#[test]
fn test_effective_tpm_calcs() {
    let mut calc = TpmCalculator::new(0, [8, 8]);
    let mut neos = vec![create_neoepitope(0, "2:1000:A:T", NeoEpitopeType::Missense, 5, [10.0, 10.0], "ABCDEFGH", "IJKLMNOP")];
    let mut peptides = Vec::new();
    calc.compute("", &mut neos, 2.0, &mut peptides);
    assert_eq!(7, neos[0].peptides.len());
}

#[test]
fn test_complex_peptide_tpm_calcs() {
    let mut calc = TpmCalculator::new(0, [4, 6]);
    let var_info = "1:1000-2:2000";
    let mut n1 = create_neoepitope(0, var_info, NeoEpitopeType::InframeFusion, 5, [10.0, 50.0], "ABC", "DEF");
    n1.transcripts[0].push(TRANS_NAME_1.into());
    n1.transcripts[1].push(TRANS_NAME_1.into());
    let mut n2 = create_neoepitope(1, var_info, NeoEpitopeType::InframeFusion, 5, [20.0, 60.0], "XBC", "DEF");
    n2.transcripts[0].push(TRANS_NAME_2.into());
    n2.transcripts[1].push(TRANS_NAME_3.into());
    let mut n3 = create_neoepitope(2, var_info, NeoEpitopeType::InframeFusion, 5, [30.0, 70.0], "ZBC", "DEG");
    n3.transcripts[0].push(TRANS_NAME_3.into());
    n3.transcripts[1].push(TRANS_NAME_4.into());
    let mut n4 = create_neoepitope(3, var_info, NeoEpitopeType::InframeFusion, 5, [10.0, 80.0], "ABC", "DEH");
    n4.transcripts[0].push(TRANS_NAME_1.into()); // matches 1
    n4.transcripts[1].push(TRANS_NAME_2.into());
    let mut neos = vec![n1, n2, n3, n4];
    let mut peptides = Vec::new();
    calc.compute("", &mut neos, 1.0, &mut peptides);

    assert_eq!(6, neos[0].peptides.len());
    let close = |p: &str, v: f64| {
        let t = get_peptide(p, &peptides).effective_tpm;
        assert!((t - v).abs() < 0.1, "{p}: {t} vs {v}");
    };
    close("ABCD", 10.0);
    close("ABCDEF", 3.8);
    close("XBCD", 20.0);
    close("BCDE", 60.0);
    close("BCDEF", 23.8);
    close("BCDEH", 6.2);
}

#[test]
fn test_poisson_range_calcs() {
    let near = |a: f64, b: f64| assert!((a - b).abs() < 0.1, "{a} vs {b}");
    near(calc_poisson_mean(10, LOW_PROBABILITY), 13.0);
    near(calc_poisson_mean(10, HIGH_PROBABILITY), 8.6);
    near(calc_poisson_mean(0, LOW_PROBABILITY), 1.39);
    near(calc_poisson_mean(0, HIGH_PROBABILITY), 0.0);
    near(calc_poisson_mean(1, LOW_PROBABILITY), 2.7);
    near(calc_poisson_mean(1, HIGH_PROBABILITY), 0.876);
    near(calc_poisson_mean(2, LOW_PROBABILITY), 3.9);
    near(calc_poisson_mean(2, HIGH_PROBABILITY), 1.7);
}

// ---------------------------------------------------------------------------------------------
// findings

fn medians_cache(dir: &TempDir, cancer_type: &str) -> TpmMediansCache {
    let f = dir.path("medians.csv");
    std::fs::write(&f, "GeneId,GeneName,TransName,All,Lung,Skin\nG1,N1,T1,1.5,2.5,9\nG1,N1,T2,0.25,0.5,9\nG2,N2,T3,4,8,9\n").unwrap();
    TpmMediansCache::new(Some(&f), &[cancer_type.to_string()]).unwrap()
}

fn sample(cancer_type: &str) -> SampleData {
    SampleData { tumor_id: "S1".into(), rna_sample_id: None, cancer_type: cancer_type.into(), has_rna: false }
}

/// NE6: the cohort TPM sums start at -1 in Java
#[test]
fn ne6_cohort_tpm_sums() {
    let dir = TempDir::new("ne6");
    let cache = medians_cache(&dir, "Lung");
    assert!(cache.has_cancer_type("Lung") && cache.has_cancer_type("All") && !cache.has_cancer_type("Other"));
    let mut n = create_neoepitope(0, "1:100:A:T", NeoEpitopeType::Missense, 0, [0.0, 0.0], "AAA", "CCC");
    n.rna = NeoRnaData::default();
    n.transcripts = [vec!["T1".into(), "T2".into()], vec!["T9".into()]];
    let tpms = std::collections::HashMap::new();
    for java in [true, false] {
        let mut n = n.clone();
        compat::with_java(java, || n.set_expression_data(&sample("Lung"), &tpms, None, &cache));
        let off = if java { -1.0 } else { 0.0 };
        assert_eq!(n.rna.tpm_pan_cancer[0], 1.75 + off);
        assert_eq!(n.rna.tpm_cancer[0], 3.0 + off);
        // no transcript with a value: the sentinel in both modes
        assert_eq!(n.rna.tpm_pan_cancer[1], NO_TPM_VALUE);
        assert_eq!(n.rna.tpm_cancer[1], NO_TPM_VALUE);
        // the TPM falls back to the cancer type's, then the pan-cancer one
        assert_eq!(n.rna.tpm(0), 3.0 + off);
        assert!(!n.rna.has_expression());
    }
    // a sample without a cancer type takes the pan-cancer medians
    let mut n2 = n.clone();
    compat::with_java(false, || n2.set_expression_data(&sample(""), &tpms, None, &cache));
    assert_eq!(n2.rna.tpm_cancer[0], NO_TPM_VALUE);
    assert_eq!(n2.rna.tpm(0), 1.75);
}

/// NE7: the skip of neo-epitopes already in a shared upstream group never applies in Java; no effect on the
/// TPMs, as a later group is a subset of an earlier one
#[test]
fn ne7_shared_upstream_groups() {
    let make = || {
        let mut v = Vec::new();
        for (i, (up, down)) in [("ABC", "DEF"), ("ABC", "DEG"), ("ABC", "DEH"), ("XBC", "DEF")].iter().enumerate() {
            let mut n = create_neoepitope(i as i32, "1:1000-2:2000", NeoEpitopeType::InframeFusion, 3, [10.0 + i as f64, 20.0 * (i + 1) as f64], up, down);
            n.transcripts[0].push(if i < 3 { TRANS_NAME_1.into() } else { TRANS_NAME_2.into() });
            n.transcripts[1].push(format!("D{i}"));
            v.push(n);
        }
        v
    };
    let run = |java: bool| {
        compat::with_java(java, || {
            let mut neos = make();
            let mut peptides = Vec::new();
            TpmCalculator::new(3, [4, 6]).compute("", &mut neos, 2.0, &mut peptides);
            peptides.iter().map(|p| (p.peptide.clone(), p.effective_tpm.to_bits(), p.ne_ids_str())).collect::<Vec<_>>()
        })
    };
    let (j, d) = (run(true), run(false));
    assert!(!j.is_empty());
    assert_eq!(j, d);
}

#[test]
fn ne8_score_file_with_dir() {
    let dir = TempDir::new("ne8");
    let f = dir.path("my_recognition.csv");
    std::fs::write(&f, "x").unwrap();
    let mut cb = hmf_common::config_builder::ConfigBuilder::new("Neo", "1.4");
    register_score_config(&mut cb);
    let args: Vec<String> = ["-score_file_dir", &dir.dir(), "-recognition_file", &f].iter().map(|s| s.to_string()).collect();
    assert!(matches!(cb.parse(&args), Ok(hmf_common::config_builder::Parsed::Run)));
    let d = dir.dir();
    let java = compat::with_java(true, || get_scoring_filename(&cb, Some(&d), None, FILE_ID_RECOGNITION));
    let fixed = compat::with_java(false, || get_scoring_filename(&cb, Some(&d), None, FILE_ID_RECOGNITION));
    assert_eq!(java, Some(format!("{d}{f}")));
    assert_eq!(fixed, Some(f.clone()));
    // without a file option, the file of that type in the directory
    assert_eq!(get_scoring_filename(&cb, Some(&d), Some("v2"), FILE_ID_POS_WEIGHT), Some(format!("{d}neo_train_pos_weight_v2.csv")));
    assert_eq!(get_scoring_filename(&cb, None, None, FILE_ID_POS_WEIGHT), None);
}

fn scorer_config(args: &[&str], java: bool) -> anyhow::Result<NeoScorerConfig> {
    let mut cb = hmf_common::config_builder::ConfigBuilder::new("Neo", "1.4");
    register(&mut cb);
    let args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
    assert!(matches!(cb.parse(&args), Ok(hmf_common::config_builder::Parsed::Run)));
    compat::with_java(java, || NeoScorerConfig::new(&cb))
}

/// NE9: Java registers -peptide_lengths and ignores it
#[test]
fn ne9_peptide_lengths() {
    let c = scorer_config(&["-sample", "S1", "-peptide_lengths", "9-10"], true).unwrap();
    assert_eq!(c.peptide_length_range, DEFAULT_PEPTIDE_LENGTH_RANGE);
    let c = scorer_config(&["-sample", "S1", "-peptide_lengths", "9-10"], false).unwrap();
    assert_eq!(c.peptide_length_range, [9, 10]);
    assert!(scorer_config(&["-sample", "S1", "-peptide_lengths", "10"], false).is_err());
    assert!(scorer_config(&["-sample", "S1", "-peptide_lengths", "11-9"], false).is_err());
    let c = scorer_config(&["-sample", "S1"], false).unwrap();
    assert_eq!(c.peptide_length_range, [8, 12]);
}

/// NE11: Java's RNA somatic VCF defaults to the sample data directory
#[test]
fn ne11_rna_vcf_default() {
    let dir = TempDir::new("ne11");
    let d = dir.dir();
    let c = scorer_config(&["-sample", "S1", "-sample_data_dir", &d], true).unwrap();
    assert_eq!(c.rna_somatic_vcf, d);
    assert!(!c.samples[0].has_rna);
    let c = scorer_config(&["-sample", "S1", "-sample_data_dir", &d], false).unwrap();
    assert_eq!(c.rna_somatic_vcf, "");
    // the directories default to the sample data directory in both modes
    assert_eq!(c.neo_dir, d);
    assert_eq!(c.lilac_dir, d);
    assert_eq!(c.output_dir, d);
}

#[test]
fn config_and_file_names() {
    let c = scorer_config(&["-sample", "S1", "-output_dir", "/tmp/x", "-output_id", "v1"], false).unwrap();
    assert_eq!(c.form_filename("neoepitope"), "/tmp/x/S1.neo.v1.neoepitope.tsv");
    assert_eq!(c.write_types, vec![OutputType::Neoepitope, OutputType::AllelePeptide]);
    assert_eq!(c.similarity_threshold, 10.0);
    assert_eq!(c.likelihood_threshold, 0.0);
    let c = scorer_config(&["-sample", "S1", "-write_types", "PEPTIDE;NEOEPITOPE"], false).unwrap();
    assert_eq!(c.write_types, vec![OutputType::Peptide, OutputType::Neoepitope]);
    assert!(scorer_config(&["-sample", "S1", "-write_types", "PEPTIDES"], false).is_err());
    // RNA needs both the ISOFOX directory and the RNA somatic VCF
    let dir = TempDir::new("cfg");
    let d = dir.dir();
    let c = scorer_config(&["-sample", "S1", "-isofox_dir", &d, "-rna_sample", "S1_RNA"], false).unwrap();
    assert!(!c.samples[0].has_rna);
    assert_eq!(c.samples[0].rna_sample_id.as_deref(), Some("S1_RNA"));
    let vcf = dir.path("v.vcf");
    std::fs::write(&vcf, "").unwrap();
    let c = scorer_config(&["-sample", "S1", "-isofox_dir", &d, "-rna_somatic_vcf", &vcf], false).unwrap();
    assert!(c.samples[0].has_rna);
    // several samples from a sample id file: neo_cohort files, no RNA sample id
    let f = dir.path("samples.csv");
    std::fs::write(&f, "SampleId,CancerType,HasRna\nS1,Lung,true\n#S2,,\nS3,,false\n").unwrap();
    let c = scorer_config(&["-sample_id_file", &f, "-output_dir", &d], false).unwrap();
    assert_eq!(c.samples.len(), 2);
    assert_eq!((c.samples[0].cancer_type.as_str(), c.samples[0].has_rna, c.samples[0].rna_sample_id.is_none()), ("Lung", true, true));
    assert_eq!(c.form_filename("peptide_scores"), format!("{d}neo_cohort.peptide_scores.tsv"));
}

// ---------------------------------------------------------------------------------------------
// binding classes

#[test]
fn rank_in_distributions() {
    let asc = [
        ScoreDistributionData { score_bucket: 0.0, score: 1.0 },
        ScoreDistributionData { score_bucket: 0.5, score: 2.0 },
        ScoreDistributionData { score_bucket: 1.0, score: 4.0 },
    ];
    assert_eq!(get_rank(&asc, 0.5), 0.0);
    assert_eq!(get_rank(&asc, 5.0), 1.0);
    assert_eq!(get_rank(&asc, 2.0), 0.5);
    assert!((get_rank(&asc, 3.0) - 0.75).abs() < 1e-12);
    // descending, as the score distributions are (the best score first)
    let desc = [
        ScoreDistributionData { score_bucket: 0.0, score: 10.0 },
        ScoreDistributionData { score_bucket: 0.1, score: 6.0 },
        ScoreDistributionData { score_bucket: 0.2, score: 2.0 },
    ];
    assert_eq!(get_rank(&desc, 11.0), 0.0);
    assert_eq!(get_rank(&desc, 1.0), 1.0);
    assert!((get_rank(&desc, 8.0) - 0.05).abs() < 1e-12);
    assert_eq!(get_rank(&desc, 6.0 + 1e-11), 0.1);
}

#[test]
fn binding_likelihoods() {
    let dir = TempDir::new("lik");
    let f = dir.path("lik.csv");
    let row = |a: &str, l: i32, v: f64| {
        let vals: Vec<String> = (0..15).map(|i| format!("{}", v / (1 << i) as f64)).collect();
        format!("{a},{l},{}\n", vals.join(","))
    };
    let mut text = String::from("Allele,PeptideLength,P0\n");
    for l in 8..=12 {
        text += &row("A0201", l, 0.4);
    }
    for l in 8..=12 {
        text += &row("B07", l, 0.2);
    }
    std::fs::write(&f, text).unwrap();
    let mut b = BindingLikelihood::default();
    assert!(b.load_likelihoods(&f).unwrap());
    // below the first bucket: its value; between buckets: interpolated
    assert_eq!(b.binding_likelihood("A0201", "AAAAAAAAA", 0.00001), 0.4);
    let mid = b.binding_likelihood("A0201", "AAAAAAAAA", 0.000075);
    assert!((mid - 0.3).abs() < 1e-12, "{mid}");
    // the last bucket and above: 0
    assert_eq!(b.binding_likelihood("A0201", "AAAAAAAAA", 0.9), 0.0);
    // the two-digit type, then the gene, stand in for an allele without data
    assert_eq!(b.binding_likelihood("B0702", "AAAAAAAAA", 0.00001), 0.2);
    assert_eq!(b.binding_likelihood("C0702", "AAAAAAAAA", 0.00001), -1.0);
    // lengths outside 8 to 12
    assert_eq!(b.binding_likelihood("A0201", "AAAAAAA", 0.00001), -1.0);
    // NE13: an invalid rank (no score distribution) gives Java the highest likelihood
    assert_eq!(compat::with_java(true, || b.binding_likelihood("A0201", "AAAAAAAAA", INVALID_SCORE)), 0.4);
    assert_eq!(compat::with_java(false, || b.binding_likelihood("A0201", "AAAAAAAAA", INVALID_SCORE)), -1.0);
    assert!(!b.load_likelihoods(&dir.path("missing.csv")).unwrap());
}

#[test]
fn expression_likelihoods() {
    let dir = TempDir::new("exp");
    let f = dir.path("exp.csv");
    std::fs::write(&f, "DataType,Bucket,BindingRate,TotalCount\nDecile,1,0.9,1\nTpmLevel,1,0.1,10\nTpmLevel,3,0.3,10\nTpmLevel,7,0.5,10\n").unwrap();
    let mut e = ExpressionLikelihood::default();
    assert!(!e.has_data());
    assert!(e.load_tpm_rates(&f).unwrap());
    assert_eq!(e.calc_likelihood(-1.0), 0.1);
    assert_eq!(e.calc_likelihood(0.0), 0.1);
    assert_eq!(e.calc_likelihood(10.0), 0.5);
    assert_eq!(e.calc_likelihood(3.0), 0.3);
    assert!((e.calc_likelihood(2.0) - 0.2).abs() < 1e-12);
    assert!((e.calc_likelihood(6.0) - 0.45).abs() < 1e-12);
}

#[test]
fn matrix_and_flank_scores() {
    let mut m = BindScoreMatrix::new("A0201", 3);
    for aa in 0..20 {
        for pos in 0..3 {
            m.set_weight(aa, pos, (aa * 10 + pos) as f64 * 0.01);
        }
    }
    // A=0, C=1, D=2
    assert!((m.calc_score("ACD") - (0.0 + 0.11 + 0.22)).abs() < 1e-12);
    assert_eq!(m.calc_score("ACDE"), -1000.0);
    assert_eq!(m.calc_score("ACU"), -1000.0);

    let dir = TempDir::new("flank");
    let f = dir.path("flank.csv");
    let mut text = String::from("DataType,AminoAcid,U3,U2,U1,D1,D2,D3\n");
    for (i, aa) in "ACDEFGHIKLMNPQRSTVWYmX".chars().enumerate() {
        let v: Vec<String> = (0..6).map(|p| format!("{}", (i * 6 + p) as f64)).collect();
        text += &format!("PosWeights,{aa},{}\n", v.join(","));
    }
    text += "BindCounts,A,9,9,9,9,9,9\n";
    std::fs::write(&f, text).unwrap();
    let mut fs = FlankScores::default();
    assert_eq!(fs.calc_score("A", "A"), 0.0);
    assert!(fs.load_pos_weights(&f).unwrap());
    // A row 0: 0..5; C row 1: 6..11; m row 20; X row 21
    assert_eq!(fs.calc_score("", ""), 0.0);
    assert_eq!(fs.calc_score("AAC", "CAA"), 0.0 + 1.0 + 8.0 + 9.0 + 4.0 + 5.0);
    assert_eq!(fs.calc_score("AC", "C"), 1.0 + 8.0 + 9.0);
    // an empty flank takes the start / stop placeholder; an unknown letter weighs 0
    assert_eq!(fs.calc_score("", "U"), (20 * 6 + 2) as f64);
    assert_eq!(fs.calc_score("C", ""), 8.0 + (21 * 6 + 3) as f64);
}

#[test]
fn recognition_similarity() {
    let dir = TempDir::new("recog");
    let f = dir.path("recog.csv");
    std::fs::write(
        &f,
        "Peptide,Allele,Immunogenic,Source\nAAAAAAAA,HLA-A*02:01,TRUE,x\nAAAAAAAC,A0201,1,x\nCCCCCCCC,A0201,0,x\nAAAAAAAW,B0702,true,x\nAAAAAAA,A0201,TRUE,x\n",
    )
    .unwrap();
    let mut r = RecognitionSimilarity::default();
    assert!(r.load_data(&f).unwrap());
    // the identical peptide: 0; one substitution A->C: (4 + 9) / 2 - 0 = 6.5
    assert_eq!(r.calc_similarity("A0201", "AAAAAAAA"), 0.0);
    assert_eq!(r.calc_similarity("A0201", "AAAAAAAD"), 7.0);
    assert_eq!(r.calc_similarity("C0101", "AAAAAAAA"), -1.0);
    assert_eq!(r.calc_similarity("A0201", "AAAAAAAAA"), -1.0);
    // the other alleles: B0702's AAAAAAAW, A->W (4 + 11) / 2 + 3 = 10.5
    assert_eq!(r.calc_other_allele_similarity("A0201", "AAAAAAAA"), 10.5);
    assert_eq!(r.calc_other_allele_similarity("B0702", "AAAAAAAA"), 0.0);
    // with the self check, the same peptide does not count
    r.check_self_similarity = true;
    assert_eq!(r.calc_similarity("A0201", "AAAAAAAA"), 6.5);
    r.check_self_similarity = false;
    // the scorer's per-peptide similarities give the same values for every allele
    let mut s = BindScorer::new(ScoreConfig::default());
    s.recognition_similarity = r;
    let sims = s.peptide_similarities("AAAAAAAD").unwrap();
    for a in ["A0201", "B0702", "C0101"] {
        assert_eq!(sims.similarity(a), s.recognition_similarity.calc_similarity(a, "AAAAAAAD"));
        assert_eq!(sims.other_allele_similarity(a), s.recognition_similarity.calc_other_allele_similarity(a, "AAAAAAAD"));
    }
}

#[test]
fn score_data_without_a_matrix() {
    let s = BindScorer::new(ScoreConfig::default());
    let mut b = BindData::new("DPB10401", "AAAAAAAA", "", "");
    s.calc_score_data(&mut b, None);
    assert_eq!((b.score, b.likelihood_rank, b.recognition_similarity), (INVALID_SCORE, INVALID_SCORE, INVALID_SCORE));
    // such rows are not written: no value passes the thresholds
    assert!(!passes_threshold(b.likelihood_rank, 0.0));
    assert!(passes_threshold(0.5, 0.0));
    assert!(passes_threshold(0.01, 0.02) && !passes_threshold(0.03, 0.02));
    assert!(!passes_threshold(0.01, -1.0));
    // setTPM tests the current TPM, which starts at 0: a negative TPM is kept
    b.set_tpm(-1.5);
    assert_eq!(b.tpm(), -1.5);
    b.set_tpm(2.0);
    assert_eq!(b.tpm(), -1.5);
}

// ---------------------------------------------------------------------------------------------
// peptides, files

#[test]
fn neo_id_order() {
    // Sets.newHashSet(neId): a table of 2 buckets, doubling
    let mut p = PeptideScoreData::new("AAAAAAAA", "", "", 20);
    assert_eq!(p.ne_ids_str(), "20");
    p.ne_ids.insert(5);
    assert_eq!(p.ne_ids_str(), "20;5");
    p.ne_ids.insert(3);
    assert_eq!(p.ne_ids_str(), "20;5;3");
    let mut q = PeptideScoreData::new("AAAAAAAA", "", "", 1);
    q.ne_ids.insert(2);
    q.ne_ids.insert(17);
    assert_eq!(q.ne_ids_str(), "1;17;2");
}

#[test]
fn readers() {
    let dir = TempDir::new("readers");
    // neo-epitopes: the gene names joined when they differ, transcripts split
    let header = hmf_common::neo::neo_epitope_header();
    let cols: Vec<&str> = header.split('\t').collect();
    let line = |vals: &[(&str, &str)]| {
        cols.iter().map(|c| vals.iter().find(|(k, _)| k == c).map(|(_, v)| *v).unwrap_or("0")).collect::<Vec<_>>().join("\t")
    };
    let l1 = line(&[
        ("VariantType", "INFRAME_FUSION"),
        ("VariantInfo", "1:100:1;2:200:-1"),
        ("GeneNameUp", "G1"),
        ("GeneNameDown", "G2"),
        ("UpstreamAA", "ABC"),
        ("DownstreamAA", "DEF"),
        ("NovelAA", ""),
        ("UpTranscripts", "T1;T2"),
        ("DownTranscripts", "T3"),
        ("VariantCopyNumber", "1.5"),
        ("NeId", "4"),
    ]);
    std::fs::write(dir.path("S1.neo.neo_data.tsv"), format!("{header}\n{l1}\n")).unwrap();
    let neos = load_neo_epitopes("S1", &dir.dir()).unwrap().unwrap();
    assert_eq!(neos.len(), 1);
    assert_eq!((neos[0].id, neos[0].gene_name.as_str(), neos[0].variant_copy_number), (4, "G1_G2", 1.5));
    assert_eq!(neos[0].transcripts, [vec!["T1".to_string(), "T2".to_string()], vec!["T3".to_string()]]);
    assert!(load_neo_epitopes("S2", &dir.dir()).unwrap().is_none());

    // LILAC alleles: '*' and ':' removed, the variant counts summed
    let lh = "Genes\tAllele\tRefTotal\tRefUnique\tRefShared\tRefWild\tTumorTotal\tTumorUnique\tTumorShared\tTumorWild\tRnaTotal\tRnaUnique\tRnaShared\tRnaWild\tTumorCopyNumber\tSomaticMissense\tSomaticNonsenseOrFrameshift\tSomaticSplice\tSomaticSynonymous\tSomaticInframeIndel";
    std::fs::write(
        dir.path("S1.lilac.tsv"),
        format!("{lh}\nMHC_CLASS_1\tA*02:01\t1\t1\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0.98\t0.0\t1.0\t0.0\t5.0\t1.0\nMHC_CLASS_1\tB*07:02\t1\t1\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0.4\t0\t0\t0\t0\t0\n"),
    )
    .unwrap();
    let a = load_allele_coverage("S1", &dir.dir()).unwrap().unwrap();
    assert_eq!((a[0].allele.as_str(), a[0].copy_number, a[0].variant_count, a[0].is_lost()), ("A0201", 0.98, 2.0, true));
    assert_eq!((a[1].allele.as_str(), a[1].is_lost()), ("B0702", true));
    assert!(load_allele_coverage("S2", &dir.dir()).unwrap().is_none());

    // ISOFOX: absent files give an empty list without RNA, none with it
    let s = SampleData { tumor_id: "S1".into(), rna_sample_id: Some("S1_RNA".into()), cancer_type: String::new(), has_rna: false };
    assert_eq!(load_rna_neo_data(&s, &dir.dir()).unwrap().map(|v| v.len()), Some(0));
    let s_rna = SampleData { has_rna: true, ..s.clone() };
    assert!(load_rna_neo_data(&s_rna, &dir.dir()).unwrap().is_none());
    std::fs::write(
        dir.path("S1_RNA.isf.neoepitope.tsv"),
        "NeId\tVariantType\tVariantInfo\tFragmentCount\tBaseDepthUp\tBaseDepthDown\n4\tINFRAME_FUSION\t1:100:1;2:200:-1\t7\t20\t30\n1\tMISSENSE\t1:1:A:T\t3\t4\t4\n",
    )
    .unwrap();
    let r = load_rna_neo_data(&s, &dir.dir()).unwrap().unwrap();
    assert_eq!(r.len(), 1);
    let mut n = neos[0].clone();
    n.set_fusion_rna_support(&r);
    assert_eq!((n.rna.fragment_support(), n.rna.average_base_depth(), n.rna.has_coverage()), (7, 25.0, true));
    std::fs::write(dir.path("S1_RNA.isf.transcript_data.tsv"), "GeneId\tGeneName\tTransName\tAdjTPM\nG\tN\tT1\t2.5\nG\tN\tT1\t3\n").unwrap();
    assert!(load_transcript_expression(&dir.dir(), "S1_RNA").is_err());
}

#[test]
fn rna_somatic_variants() {
    let dir = TempDir::new("rnavcf");
    let vcf = dir.path("S1.rna.vcf");
    let text = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
        ##FORMAT=<ID=RC_CNT,Number=7,Type=Integer,Description=\"r\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS1_RNA\n\
        1\t100\t.\tA\tT\t.\tPASS\t.\tGT:DP:RC_CNT\t0/1:50:5,1,2,3,0,40,50\t0/1:30:4,0,1,1,0,25,30\n\
        1\t200\t.\tC\tG\t.\tLOW\t.\tGT:DP:RC_CNT\t0/1:50:5,1,2,3,0,40,50\t0/1:30:4,0,1,1,0,25,30\n\
        1\t300\t.\tC\tG\t.\tPASS\t.\tGT:DP:RC_CNT\t0/1:50:5,1,2,3,0,40,50\t0/1:30:4,0,1,1,0,25,30\n";
    std::fs::write(&vcf, text).unwrap();
    let mut n1 = create_neoepitope(0, "1:100:A:T", NeoEpitopeType::Missense, 0, [0.0, 0.0], "AAA", "CCC");
    let n2 = create_neoepitope(1, "1:200:C:G", NeoEpitopeType::Missense, 0, [0.0, 0.0], "AAA", "CCC");
    let s = SampleData { tumor_id: "S1".into(), rna_sample_id: Some("S1_RNA".into()), cancer_type: String::new(), has_rna: true };
    let v = load_somatic_variants(&s, &dir.path("*.rna.vcf"), &[&n1, &n2]).unwrap().unwrap();
    // the filtered variant and the one without a neo-epitope are left out
    assert_eq!(v.len(), 1);
    assert_eq!((v[0].tumor_fragments, v[0].tumor_depth, v[0].rna_fragments, v[0].rna_depth), (11, 50, 6, 30));
    n1.rna = NeoRnaData::default();
    n1.set_mutation_rna_support(&v);
    assert_eq!((n1.rna.fragment_support(), n1.rna.average_base_depth()), (6, 30.0));
    // no RNA genotype: nothing matched
    let s2 = SampleData { rna_sample_id: Some("OTHER".into()), ..s.clone() };
    assert_eq!(load_somatic_variants(&s2, &vcf, &[&n1]).unwrap().map(|v| v.len()), Some(0));
    // an absent VCF: none with RNA
    assert!(load_somatic_variants(&s, &dir.path("missing.vcf"), &[&n1]).unwrap().is_none());
}

#[test]
fn tpm_normalisation_and_effective_tpm() {
    // the missense neo-epitopes with coverage and expression: summed TPM over summed depth
    let mut a = create_neoepitope(0, "1:1:A:T", NeoEpitopeType::Missense, 5, [10.0, 10.0], "AAA", "CCC");
    let b = create_neoepitope(1, "1:2:A:T", NeoEpitopeType::Missense, 10, [30.0, 30.0], "AAA", "CCC");
    let c = create_neoepitope(2, "1:3:A:TT", NeoEpitopeType::Frameshift, 10, [99.0, 99.0], "AAA", "CCC");
    let f = calculate_tpm_normalisation(&[a.clone(), b.clone(), c]);
    assert!((f - 40.0 / 30.0).abs() < 1e-12);
    assert_eq!(calculate_tpm_normalisation(&[]), NO_TPM_VALUE);
    // effective TPM: 0.8 raw + 0.2 expected clamped to the Poisson range
    assert!((calc_effective_tpm(2.0, 10.0, 100.0, 7.0, 3.0) - (8.0 + 0.2 * 14.0)).abs() < 1e-12);
    assert!((calc_effective_tpm(2.0, 10.0, 1.0, 7.0, 3.0) - (8.0 + 0.2 * 6.0)).abs() < 1e-12);
    // the copy number factor: subclonal likelihood 0 for a point mutation keeps it at 1
    a.set_calculated_tpm_values(f, 2.0, 7.0, 3.0);
    assert_eq!(a.copy_number_factor, 1.0);
    assert_eq!(a.expected_tpm, 10.0);
    assert!((a.raw_effective_tpm - 5.0 * f).abs() < 1e-12);
    // fusions count as clonal: variant CN 1 over the segment CN 2 halves the expected TPM
    let mut fu = create_neoepitope(3, "1:1:1;2:2:-1", NeoEpitopeType::InframeFusion, 0, [10.0, 10.0], "AAA", "CCC");
    fu.copy_number = 2.0;
    fu.set_calculated_tpm_values(NO_TPM_VALUE, 2.0, 1.39, 0.0);
    assert_eq!(fu.copy_number_factor, 0.5);
    assert_eq!((fu.expected_tpm, fu.effective_tpm, fu.raw_effective_tpm), (5.0, 5.0, NO_TPM_VALUE));
}
