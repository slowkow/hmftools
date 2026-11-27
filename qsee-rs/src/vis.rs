//! qsee `plot`: VisDataFile (`<tumor>.qsee.vis.data.tsv.gz`, what the R script draws from), and
//! QseePlot / QseePlotTask: the R script `plot_qc.R` run as Java runs it (RExecutor), unless -skip_plot.

use anyhow::{anyhow, bail, Result};
use std::io::Write;

use hmf_common::delim_file::JavaWriter;

use crate::common::{decimal_8, generate_filename, SampleType, COL_FEATURE_NAME, COL_FEATURE_TYPE, COL_FEATURE_VALUE, COL_SAMPLE_ID, COL_SAMPLE_TYPE, COL_SOURCE_TOOL, MULTISAMPLE_SAMPLE_ID};
use crate::config::{QseeConfig, QseePlotConfig};
use crate::feature::SampleFeatures;

/// plot_qc.R as the jar carries it (identical to qsee/src/main/R/plot_qc.R)
pub const PLOT_SCRIPT: &str = include_str!("../resources/plot_qc.R");

pub fn vis_data_filename(base: &str, sample_id: &str, output_id: Option<&str>) -> String {
    generate_filename(base, sample_id, "vis.data", output_id, "tsv.gz")
}

/// VisDataFile.generateFilename(config)
pub fn vis_data_filename_for(config: &QseeConfig) -> String {
    let sample = if config.is_single_patient() { config.sample_ids(SampleType::Tumor)[0].as_str() } else { MULTISAMPLE_SAMPLE_ID };
    vis_data_filename(&config.output_dir, sample, config.output_id.as_deref())
}

/// QcStatusFile.generateFilename(config)
pub fn qc_status_filename_for(config: &QseeConfig) -> String {
    let sample = if config.is_single_patient() { config.sample_ids(SampleType::Tumor)[0].as_str() } else { MULTISAMPLE_SAMPLE_ID };
    generate_filename(&config.output_dir, sample, "status", config.output_id.as_deref(), "tsv.gz")
}

/// VisDataFile.write: a line per feature with a value (NaN values are left out). A line that cannot be
/// formed ends the writing; the file is closed with the lines before it, as Java's try-with-resources.
pub fn write_vis_data(path: &str, samples: &[SampleFeatures]) -> Result<()> {
    let mut w = JavaWriter::create(path).map_err(|e| anyhow!("Failed to write to file: {path}: {e:#}"))?;
    hmf_common::info!("Writing sample vis data to: {path}");
    let header = [COL_SAMPLE_ID, COL_SAMPLE_TYPE, COL_SOURCE_TOOL, COL_FEATURE_TYPE, COL_FEATURE_NAME, COL_FEATURE_VALUE, "PlotMetadata"].join("\t");
    w.write_all(header.as_bytes())?;
    w.write_all(b"\n")?;
    let mut result = Ok(());
    'outer: for s in samples {
        for f in &s.features {
            if f.value.is_nan() {
                continue;
            }
            let metadata = match f.metadata.display_string() {
                Ok(t) => t,
                Err(e) => {
                    result = Err(e);
                    break 'outer;
                }
            };
            let line = [
                s.sample_id.as_str(),
                s.sample_type.name(),
                f.key.source_tool.map(|t| t.name()).unwrap_or("null"),
                f.key.feature_type.map(|t| t.name()).unwrap_or("null"),
                f.key.name.as_str(),
                &decimal_8(f.value),
                &metadata,
            ]
            .join("\t");
            w.write_all(line.as_bytes())?;
            w.write_all(b"\n")?;
        }
    }
    w.finish()?;
    result
}

/// QseePlotTask.run: `Rscript plot_qc.R <tumor> <reference|NA> <vis data> <cohort percentiles|NA> <pdf>
/// <show warnings> <single patient> <log level>`, the script written to a temporary file first
pub fn run_plot_task(tumor: &str, reference: Option<&str>, vis_data: &str, cohort: Option<&str>, plot_path: &str, show_warnings: bool, single_patient: bool, log_level: &str) -> Result<()> {
    let dir = std::env::temp_dir();
    let script = dir.join(format!("script{}{}.R", std::process::id(), unique_suffix()));
    std::fs::write(&script, PLOT_SCRIPT).map_err(|e| anyhow!("Failed to run QC plot script: {e}"))?;
    let args = [tumor, reference.unwrap_or("NA"), vis_data, cohort.unwrap_or("NA"), plot_path, if show_warnings { "true" } else { "false" }, if single_patient { "true" } else { "false" }, log_level];
    hmf_common::debug!("executing R script via command: Rscript {} {}", script.display(), args.join(" "));
    let status = std::process::Command::new("Rscript").arg(&script).args(args).status();
    let _ = std::fs::remove_file(&script);
    match status {
        Err(e) => bail!("Failed to run QC plot script: {e}"),
        Ok(s) if !s.success() => {
            hmf_common::error!("error executing R script");
            bail!("java.lang.RuntimeException: R script exited with {s}")
        }
        Ok(_) => Ok(()),
    }
}

fn unique_suffix() -> u128 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0)
}

/// QseePlot.run: one report for a single patient, else one per tumor sample (on -threads threads); the
/// PDFs of several patients are not merged by the port (Java merges them with PDFBox for -merge_plots)
pub fn run_plot(pc: &QseePlotConfig, log_level: &str) -> Result<()> {
    hmf_common::info!("Plotting QC metrics");
    let single = pc.samples.is_single_patient();
    let vis_path = |tumor: &str| -> String {
        match &pc.vis_data_file {
            Some(v) => {
                hmf_common::debug!("Using existing vis data file: {v}");
                v.clone()
            }
            None => vis_data_filename(&pc.output_dir, if single { tumor } else { MULTISAMPLE_SAMPLE_ID }, pc.output_id.as_deref()),
        }
    };
    let plot_path = |tumor: &str| generate_filename(&pc.output_dir, tumor, "vis.report", pc.output_id.as_deref(), "pdf");
    let reference = |i: usize| -> Option<&str> { if pc.samples.reference_ids.is_empty() { None } else { pc.samples.reference_ids.get(i).map(|s| s.as_str()) } };
    if single {
        let t = &pc.samples.tumor_ids[0];
        run_plot_task(t, reference(0), &vis_path(t), pc.cohort_percentiles_file.as_deref(), &plot_path(t), pc.show_plot_warnings, true, log_level)?;
    } else {
        for (i, t) in pc.samples.tumor_ids.iter().enumerate() {
            run_plot_task(t, reference(i), &vis_path(t), pc.cohort_percentiles_file.as_deref(), &plot_path(t), pc.show_plot_warnings, false, log_level)?;
        }
        if pc.merge_plots {
            hmf_common::warn!("-merge_plots: the port does not merge PDFs; the per-sample reports are kept");
        }
    }
    hmf_common::info!("Completed plotting QC metrics");
    Ok(())
}
