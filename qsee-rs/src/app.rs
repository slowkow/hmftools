//! QseeApplication.run and FeaturePrep: the features of each sample, the vis data and QC status files,
//! the gene limit of detection, and the plot.

use anyhow::{anyhow, Result};

use hmf_common::config_builder::ConfigBuilder;

use crate::cohort::{prep_cohort_category, read_cohort_percentiles, FeatureMatrix};
use crate::common::SampleType;
use crate::config::{log_level_name, PrepCategory, QseeConfig, QseePlotConfig};
use crate::feature::SampleFeatures;
use crate::prep::{run_category_task, PrepError};
use crate::status::write_qc_status_file;
use crate::{lod, vis};

/// An error already logged that ends the run with exit status 1 (Java's System.exit(1))
#[derive(Debug)]
pub struct Exit;

impl std::fmt::Display for Exit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "exit")
    }
}
impl std::error::Error for Exit {}

pub fn prep_error(e: PrepError) -> anyhow::Error {
    match e {
        PrepError::Exit => anyhow::Error::new(Exit),
        PrepError::Fatal(e) => e,
        PrepError::Missing(m) => anyhow!("missing input: {m}"),
    }
}

/// FeaturePrep.createCategoryPreps: the selected categories in preparation order, the others logged
pub fn selected_categories(config: &QseeConfig) -> Vec<PrepCategory> {
    let mut out = Vec::new();
    for c in PrepCategory::PREP_ORDER {
        if !config.categories.contains(&c) {
            hmf_common::info!("Skipping unselected category({})", c.name());
            continue;
        }
        out.push(c);
    }
    out
}

/// FeaturePrep.prepSample: the categories in order, their features one after the other
pub fn prep_sample(config: &QseeConfig, sample_type: SampleType, sample_id: &str) -> Result<SampleFeatures> {
    hmf_common::info!("Extracting sample data - sampleType({}) sample({sample_id})", sample_type.name());
    let mut features = Vec::new();
    for c in selected_categories(config) {
        hmf_common::debug!("Extracting category({})", c.name());
        features.extend(run_category_task(config, c, sample_id, sample_type, None).map_err(prep_error)?);
    }
    Ok(SampleFeatures { sample_id: sample_id.to_string(), sample_type, features })
}

/// FeaturePrep.prepMultiSample: each category for all samples, the features of each sample in key order
pub fn prep_multi_sample(config: &QseeConfig, sample_type: SampleType) -> Result<Vec<SampleFeatures>> {
    hmf_common::info!("Extracting multi-sample data - sampleType({})", sample_type.name());
    let ids = config.sample_ids(sample_type).to_vec();
    let mut m = FeatureMatrix::with_samples(&ids);
    for c in selected_categories(config) {
        hmf_common::info!("Extracting category({})", c.name());
        prep_cohort_category(config, c, sample_type, &mut m).map_err(prep_error)?;
    }
    m.sort_feature_keys();
    ids.iter().map(|id| Ok(SampleFeatures { sample_id: id.clone(), sample_type, features: m.row_of(id)? })).collect()
}

fn prep_sample_features(config: &QseeConfig, sample_type: SampleType) -> Result<Vec<SampleFeatures>> {
    let ids = config.sample_ids(sample_type);
    match ids.len() {
        0 => Ok(Vec::new()),
        1 => Ok(vec![prep_sample(config, sample_type, &ids[0])?]),
        _ => prep_multi_sample(config, sample_type),
    }
}

/// QseeApplication.run
pub fn run(config: &QseeConfig, cb: &ConfigBuilder) -> Result<()> {
    hmf_common::info!("Running Qsee");
    let mut samples = prep_sample_features(config, SampleType::Tumor)?;
    samples.extend(prep_sample_features(config, SampleType::Normal)?);
    if !config.is_single_patient() {
        samples.sort_by(|a, b| hmf_common::java_text::compare(&a.sample_id, &b.sample_id));
    }
    if let Some(f) = &config.cohort_percentiles_file {
        let cohort = read_cohort_percentiles(f)?;
        for s in &samples {
            for feature in &s.features {
                cohort.warn_if_missing(s.sample_type, &feature.key)?;
            }
        }
    }
    vis::write_vis_data(&vis::vis_data_filename_for(config), &samples)?;
    write_qc_status_file(&vis::qc_status_filename_for(config), &samples)?;
    lod::run(config)?;
    if config.skip_plot {
        hmf_common::info!("Skipping plotting");
        return Ok(());
    }
    let pc = QseePlotConfig::from(cb)?;
    vis::run_plot(&pc, &log_level_name(cb))
}
