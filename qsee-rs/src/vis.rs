//! qsee `plot`: VisDataFile (`<tumor>.qsee.vis.data.tsv.gz`, what the R script draws from), and
//! QseePlot / QseePlotTask: the report the R script `plot_qc.R` draws, drawn natively (`crate::plot`), unless
//! -skip_plot.

use anyhow::{anyhow, Result};
use std::io::Write;

use hmf_common::delim_file::JavaWriter;

use crate::common::{decimal_8, generate_filename, SampleType, COL_FEATURE_NAME, COL_FEATURE_TYPE, COL_FEATURE_VALUE, COL_SAMPLE_ID, COL_SAMPLE_TYPE, COL_SOURCE_TOOL, MULTISAMPLE_SAMPLE_ID};
use crate::config::{QseeConfig, QseePlotConfig};
use crate::feature::SampleFeatures;


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

/// VisDataFile.KEEP_NAN_FEATURE_TYPES (upstream 883ab11bba): the BQR features keep a NaN value (a group with too
/// few counts), which the report hides
fn keeps_nan(t: Option<crate::feature::FeatureType>) -> bool {
    matches!(t, Some(crate::feature::FeatureType::BqrPerOrigQual) | Some(crate::feature::FeatureType::BqrPerSnv96Context))
}

/// VisDataFile.write: a line per feature with a value (NaN values are left out, except the BQR features'). A
/// line that cannot be formed ends the writing; the file is closed with the lines before it, as Java's
/// try-with-resources.
pub fn write_vis_data(path: &str, samples: &[SampleFeatures]) -> Result<()> {
    let mut w = JavaWriter::create(path).map_err(|e| anyhow!("Failed to write to file: {path}: {e:#}"))?;
    hmf_common::info!("Writing sample vis data to: {path}");
    let header = [COL_SAMPLE_ID, COL_SAMPLE_TYPE, COL_SOURCE_TOOL, COL_FEATURE_TYPE, COL_FEATURE_NAME, COL_FEATURE_VALUE, "PlotMetadata"].join("\t");
    w.write_all(header.as_bytes())?;
    w.write_all(b"\n")?;
    let mut result = Ok(());
    'outer: for s in samples {
        for f in &s.features {
            if f.value.is_nan() && !keeps_nan(f.key.feature_type) {
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

/// QseePlotTask.run: the report plot_qc.R draws (`Rscript plot_qc.R <tumor> <reference|NA> <vis data> <cohort
/// percentiles|NA> <pdf> <show warnings> <single patient> <log level>`), drawn natively (`crate::plot`): the PDF (two
/// pages since upstream 883ab11bba), and for a single patient a PNG per page beside it. Returns the PDF's pages.
/// (show warnings only sets R's warning output.)
pub fn run_plot_task(tumor: &str, reference: Option<&str>, vis_data: &str, cohort: Option<&str>, plot_path: &str, show_warnings: bool, single_patient: bool, log_level: &str) -> Result<Vec<hmf_plot::display::Figure>> {
    let _ = (show_warnings, log_level);
    hmf_common::debug!("drawing the QC report natively: {tumor} {} {vis_data} {} {plot_path}", reference.unwrap_or("NA"), cohort.unwrap_or("NA"));
    crate::plot::report(tumor, reference, vis_data, cohort, plot_path, single_patient).map_err(|e| {
        hmf_common::error!("error drawing the QC report: {e:#}");
        anyhow!("java.lang.RuntimeException: {e:#}")
    })
}

/// QseePlot.run: one report for a single patient, else one per tumor sample; with -merge_plots the reports'
/// pages go into `multisample.qsee.vis.report.pdf` and the per-sample PDFs are deleted (as Java with PDFBox).
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
        let mut pages = Vec::new();
        for (i, t) in pc.samples.tumor_ids.iter().enumerate() {
            pages.extend(run_plot_task(t, reference(i), &vis_path(t), pc.cohort_percentiles_file.as_deref(), &plot_path(t), pc.show_plot_warnings, false, log_level)?);
        }
        if pc.merge_plots {
            let dest = plot_path(MULTISAMPLE_SAMPLE_ID);
            hmf_common::info!("Merging plots to file: {dest}");
            let merged = crate::plot::write_merged(&dest, &pages);
            for t in &pc.samples.tumor_ids {
                let p = plot_path(t);
                if std::path::Path::new(&p).exists() {
                    hmf_common::debug!("Deleting plot file: {p}");
                    let _ = std::fs::remove_file(&p);
                }
            }
            merged.map_err(|e| {
                hmf_common::error!("Failed to merge plots to file: {dest}: {e:#}");
                e
            })?;
        }
    }
    hmf_common::info!("Completed plotting QC metrics");
    Ok(())
}
