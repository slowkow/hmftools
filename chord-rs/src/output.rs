//! ChordOutput: the output file names, `<dir>/<sample or "cohort">.chord.<kind>[.<output id>].tsv`.

use crate::config::ChordConfig;

pub const COHORT_FILE_PREFIX: &str = "cohort";
pub const CHORD_MUTATION_CONTEXTS: &str = ".chord.mutation_contexts";
pub const CHORD_PREDICTION: &str = ".chord.prediction";

/// formOutputFile: the directory (with its separator, as parseOutputDir leaves it) and "/" are joined
/// as Java joins them
fn form_output_file(output_dir: &str, prefix: &str, suffix: &str, output_id: Option<&str>) -> String {
    let mut f = format!("{output_dir}/{prefix}{suffix}");
    if let Some(id) = output_id.filter(|id| !id.is_empty()) {
        f.push('.');
        f.push_str(id);
    }
    f.push_str(".tsv");
    f
}

fn prefix(config: &ChordConfig) -> &str {
    if config.is_single_sample() { config.sample_ids.first().map(String::as_str).unwrap_or("") } else { COHORT_FILE_PREFIX }
}

pub fn mutation_contexts_file(config: &ChordConfig) -> String {
    form_output_file(&config.output_dir_str(), prefix(config), CHORD_MUTATION_CONTEXTS, config.output_id.as_deref())
}

pub fn predictions_file(config: &ChordConfig) -> String {
    form_output_file(&config.output_dir_str(), prefix(config), CHORD_PREDICTION, config.output_id.as_deref())
}
