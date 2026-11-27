//! CUPPA's feature extraction (com.hartwig.hmftools.cup.prep.CuppaDataPrep), and the jar's main class
//! (com.hartwig.hmftools.cup.cli.PredictionRunner: the extraction, then the prediction).
//!
//! Per sample and category: SNV (the 96 trinucleotide counts, the SNV count, 8 signature allocations by
//! least squares, SNV counts in 500 kb bins), SV (counts of six SV types from PURPLE's SV VCF and LINX's
//! clusters), SAMPLE_TRAIT (sex, MS indels per Mb, WGD from PURPLE's purity file), DRIVER (drivers from
//! the LINX or PURPLE catalog, indel drivers of three genes, LINX fusions, reported viruses), GENE_EXP and
//! ALT_SJ (Isofox's gene expression and alternative splice junctions of the RNA sample).

pub mod config;
pub mod drivers;
pub mod items;
pub mod rna;
pub mod somatic;
pub mod sv;
pub mod traits;

use std::collections::HashMap;

use anyhow::{bail, Result};
use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::{error, info, warn};

use config::{CategoryType, PrepConfig, APP_NAME, VERSION};
use items::DataItem;

/// log4j's standard levels (ConfigUtils.setLogLevel's Level.valueOf rejects other names)
const LOG4J_LEVELS: [&str; 8] = ["OFF", "FATAL", "ERROR", "WARN", "INFO", "DEBUG", "TRACE", "ALL"];

fn check_log_level(cb: &ConfigBuilder) -> Result<()> {
    if cb.has_value("log_level") {
        let level = cb.value("log_level").unwrap_or_default().trim().to_ascii_uppercase();
        if !LOG4J_LEVELS.contains(&level.as_str()) {
            bail!("java.lang.IllegalArgumentException: Unknown level constant [{level}].");
        }
    }
    Ok(())
}

/// CuppaDataPrep.main: the exit status
pub fn main(argv: &[String]) -> i32 {
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    config::register(&mut cb, false);
    match cb.parse(argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    }
    if let Err(e) = check_log_level(&cb) {
        error!("{e:#}");
        return 1;
    }
    let config = match PrepConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            error!("{e:#}");
            return 1;
        }
    };
    match run(&config) {
        Ok(()) => 0,
        Err(e) => {
            error!("{e:#}");
            1
        }
    }
}

/// CuppaDataPrep.getOutputPath
pub fn output_path(config: &PrepConfig, category: Option<CategoryType>) -> String {
    let mut path = format!("{}/", config.output_dir);
    if config.is_single_sample() {
        path.push_str(&format!("{}.cuppa_data.tsv.gz", config.sample_ids[0]));
        return path;
    }
    path.push_str("cuppa_data.cohort");
    if let Some(c) = category {
        path.push('.');
        path.push_str(&c.name().to_lowercase());
    }
    if let Some(id) = &config.output_id {
        path.push('.');
        path.push_str(id);
    }
    path.push_str(".tsv.gz");
    path
}

/// CuppaDataPrep.run
pub fn run(config: &PrepConfig) -> Result<()> {
    info!("Starting Cuppa feature extraction");
    if config.sample_ids.is_empty() {
        error!("No sample ID(s) loaded");
        bail!("no sample ids");
    }
    let start = std::time::Instant::now();
    if config.is_single_sample() {
        somatic::CACHE_VCF.store(true, std::sync::atomic::Ordering::Relaxed);
        let mut data_items: Vec<DataItem> = Vec::new();
        for category in &config.categories {
            info!("Extracting category({})", category.name());
            match extract(config, *category, 0)? {
                Some(items) => data_items.extend(items),
                None => {
                    if crate::compat::java("CU4") {
                        // SamplePrepTask leaves the items null and dataItems.addAll(null) throws
                        bail!("java.lang.NullPointerException: no {} data for sample({})", category.name(), config.sample_ids[0]);
                    }
                    warn!("sample({}) has no {} data, its features are left out", config.sample_ids[0], category.name());
                }
            }
        }
        items::write_item_list(&data_items, &output_path(config, None))?;
    } else {
        info!("Extracting CUPPA features in multi sample mode: {} samples, {} threads", config.sample_ids.len(), config.threads);
        for (i, category) in config.categories.iter().enumerate() {
            let matrix = extract_multi(config, *category)?;
            if config.write_by_category {
                items::write_matrix(&matrix, &output_path(config, Some(*category)), false)?;
            } else {
                items::write_matrix(&matrix, &output_path(config, None), i != 0)?;
            }
        }
    }
    let mins = start.elapsed().as_secs_f64() / 60.0;
    info!("Cuppa data extraction complete, mins({})", hmf_common::fmt::java_fixed(mins, 3));
    Ok(())
}

/// The items of one category for one sample (SamplePrepTask.run): None where the category prep returns
/// null (RNA data missing, or the sample marked NO_RNA)
pub fn extract(config: &PrepConfig, category: CategoryType, sample_index: usize) -> Result<Option<Vec<DataItem>>> {
    let sample_id = if category.is_rna() { &config.rna_sample_ids[sample_index] } else { &config.sample_ids[sample_index] };
    if category.is_rna() && sample_id == config::NO_RNA_SAMPLE_ID {
        return Ok(None);
    }
    let result = match category {
        CategoryType::Snv => somatic::extract(config, sample_id).map(Some),
        CategoryType::Sv => sv::extract(config, sample_id).map(Some),
        CategoryType::SampleTrait => traits::extract(config, sample_id).map(Some),
        CategoryType::Driver => drivers::extract(config, sample_id).map(Some),
        CategoryType::GeneExp => Ok(rna::extract_gene_expression(config, sample_id)),
        CategoryType::AltSj => rna::extract_alt_sj(config, sample_id),
    };
    result.map_err(|e| {
        error!("sample({sample_id}) failed to extract category({}): {e:#}", category.name());
        e
    })
}

/// The cohort matrix of one category (extractMultiSampleOneCategory): the items' indexes sorted, a value
/// per sample (None where a sample has no such item)
pub struct ItemMatrix {
    pub sample_ids: Vec<String>,
    pub indexes: Vec<items::Index>,
    pub values: HashMap<items::Index, Vec<Option<String>>>,
}

fn extract_multi(config: &PrepConfig, category: CategoryType) -> Result<ItemMatrix> {
    info!("Extracting category({})", category.name());
    let n = config.sample_ids.len();
    let results: Vec<Result<Option<Vec<DataItem>>>> = if config.threads > 1 {
        let next = std::sync::atomic::AtomicUsize::new(0);
        let slots: Vec<std::sync::Mutex<Option<Result<Option<Vec<DataItem>>>>>> = (0..n).map(|_| std::sync::Mutex::new(None)).collect();
        std::thread::scope(|s| {
            for _ in 0..config.threads.min(n) {
                s.spawn(|| loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= n {
                        break;
                    }
                    *slots[i].lock().unwrap() = Some(sample_task(config, category, i));
                });
            }
        });
        slots.into_iter().map(|m| m.into_inner().unwrap().unwrap()).collect()
    } else {
        (0..n).map(|i| sample_task(config, category, i)).collect()
    };
    let mut values: HashMap<items::Index, Vec<Option<String>>> = HashMap::new();
    let mut indexes: Vec<items::Index> = Vec::new();
    for (i, r) in results.into_iter().enumerate() {
        match r? {
            None => warn!("multi-sample feature matrix will contain nulls for sample({}) category({})", config.sample_ids[i], category.name()),
            Some(items) => {
                for item in items {
                    let slot = values.entry(item.index.clone()).or_insert_with(|| {
                        indexes.push(item.index.clone());
                        vec![None; n]
                    });
                    slot[i] = Some(item.value);
                }
            }
        }
    }
    // Collections.sort (stable) of the ConcurrentHashMap's keys; the comparator ties only equal keys
    indexes.sort_by(|a, b| a.compare(b));
    Ok(ItemMatrix { sample_ids: config.sample_ids.clone(), indexes, values })
}

fn sample_task(config: &PrepConfig, category: CategoryType, i: usize) -> Result<Option<Vec<DataItem>>> {
    const PROGRESS_INTERVAL: usize = 100;
    if i < PROGRESS_INTERVAL || i % PROGRESS_INTERVAL == 0 {
        let id = if category.is_rna() { &config.rna_sample_ids[i] } else { &config.sample_ids[i] };
        hmf_common::debug!("{}/{}: sample({id})", i + 1, config.sample_ids.len());
    }
    extract(config, category, i)
}

/// PredictionRunner.main (the jar's main class): the features (unless -features_path is given), then the
/// prediction in this process, with the arguments Java passes to `python -m cuppa.predict`
pub fn runner_main(argv: &[String]) -> i32 {
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    config::register(&mut cb, true);
    match cb.parse(argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    }
    if let Err(e) = check_log_level(&cb) {
        error!("{e:#}");
        return 1;
    }
    let config = match PrepConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            error!("{e:#}");
            return 1;
        }
    };
    let python_path = cb.value("python_path");
    // PredictionRunner.createOutputDirIfNotExist (File.mkdir: one level)
    let out = std::path::Path::new(&config.output_dir);
    if !out.exists() {
        info!("Creating output dir: {}", config.output_dir.trim_end_matches('/'));
        let _ = std::fs::create_dir(out);
    }
    let features_path = match cb.value("features_path") {
        Some(p) => {
            info!("Using pre-extracted features at: {p}");
            p
        }
        None => {
            if let Err(e) = run(&config) {
                error!("{e:#}");
                return 1;
            }
            output_path(&config, None)
        }
    };
    // PythonInterpreter.checkValid: Path.getFileName().startsWith("python") compares whole names, so only an
    // interpreter file named exactly "python" passes (CU5); the port needs no interpreter
    if crate::compat::java("CU5") {
        let p = python_path.clone().unwrap_or_default();
        let path = std::path::Path::new(&p);
        let executable = std::fs::metadata(path).map(|m| {
            use std::os::unix::fs::PermissionsExt;
            m.is_file() && m.permissions().mode() & 0o111 != 0
        });
        let named_python = path.file_name().is_some_and(|n| n == "python");
        if !executable.unwrap_or(false) || !named_python {
            error!("Invalid python interpreter: {p}");
            return 1;
        }
    }
    let mut args = vec![
        "--classifier_path".to_string(),
        cb.value("classifier_path").unwrap_or_default(),
        "--features_path".to_string(),
        features_path,
        "--output_dir".to_string(),
        config.output_dir.clone(),
    ];
    if let Some(s) = cb.value("sample") {
        args.push("--sample_id".into());
        args.push(s);
    }
    if let Some(g) = cb.value("clf_group") {
        args.push("--clf_group".into());
        args.push(g);
    }
    if let Some(cv) = cb.value("cv_predictions_path") {
        args.push("--cv_predictions_path".into());
        args.push(cv);
    }
    args.push("--log_format".into());
    args.push("%(asctime)s.%(msecs)03d [python] %(levelname)s %(name)s.%(funcName)s | %(message)s".into());
    info!("Running command: [python -m cuppa.predict {}] (in process)", args.join(" "));
    let status = crate::predict::main(&args);
    if status != 0 {
        error!("Failed to run command: python -m cuppa.predict");
        return 1;
    }
    info!("Command complete: [python -m cuppa.predict]");
    0
}
