//! ChordApplicationTest: ignored in Java, which needs R for it; the port predicts without R, so the tests
//! run here.

use super::*;
use crate::output;

fn run_application(config: &crate::config::ChordConfig) {
    let rows = crate::prep::run(config).unwrap();
    let table = crate::predict::ContextTable::from_counts(&rows);
    crate::predict::predict_to_file(&table, config.chord_model_file.as_deref(), &output::predictions_file(config)).unwrap();
}

#[test]
fn can_prep_and_predict_one_sample() {
    let mut c = config();
    c.sample_ids = vec![MINIMAL_SAMPLE.into()];
    c.ref_genome_file = Some(dummy_genome_fasta());
    c.purple_dir = Some(input_vcf_dir());
    c.output_dir = Some(tmp_output_dir("one_sample"));
    run_application(&c);
    assert!(std::path::Path::new(&output::mutation_contexts_file(&c)).exists());
    assert!(std::path::Path::new(&output::predictions_file(&c)).exists());
}

#[test]
fn can_prep_and_predict_multiple_samples() {
    let mut c = config();
    c.sample_ids = vec![MINIMAL_SAMPLE.into(), EMPTY_SAMPLE.into()];
    c.ref_genome_file = Some(dummy_genome_fasta());
    c.purple_dir = Some(input_vcf_dir());
    c.output_dir = Some(tmp_output_dir("multiple_samples"));
    c.threads = 2;
    run_application(&c);
    assert!(output::mutation_contexts_file(&c).ends_with("/cohort.chord.mutation_contexts.tsv"));
    assert!(std::path::Path::new(&output::mutation_contexts_file(&c)).exists());
    let text = std::fs::read_to_string(output::predictions_file(&c)).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 3);
    assert!(lines[1].starts_with("MINIMAL_SAMPLE\t") && lines[2].starts_with("EMPTY_SAMPLE\t"));
    // no mutations at all: the forest still predicts, and the QC rules call it undetermined
    assert!(lines[2].ends_with("\tcannot_be_determined\tnone\t<100 indels\t<30 SVs"), "{}", lines[2]);
}
