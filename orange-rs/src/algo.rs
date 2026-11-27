//! OrangeAlgo.run: every tool's results loaded and interpreted into the OrangeRecord, the plots copied.

use anyhow::Result;
use hmf_common::cytobands::CytoBands;
use hmf_common::tool_results;

use crate::config::OrangeConfig;
use crate::model::*;
use crate::plots::PlotManager;
use crate::{isofox, linx, others, purple};

/// OrangeAlgo.determinePipelineVersion: the file's major.minor version, else the default pipeline output
/// structure's name ("OA_V3_1")
fn pipeline_version(config: &OrangeConfig) -> Result<String> {
    const DEFAULT: &str = "OA_V3_1";
    let Some(file) = &config.pipeline_version_file else {
        hmf_common::info!("no pipeline version file, defaulting to {DEFAULT}");
        return Ok(DEFAULT.to_string());
    };
    match tool_results::pipeline_major_dot_minor_version(file)? {
        Some(v) => {
            hmf_common::info!("determined pipeline version to be 'v{v}'");
            Ok(v)
        }
        None => {
            hmf_common::warn!("no pipeline version determined from {file}");
            Ok(DEFAULT.to_string())
        }
    }
}

fn purple_plot(config: &OrangeConfig, suffix: &str) -> String {
    format!("{}{}{suffix}", crate::config::check_add_dir_separator(Some(&config.purple_plot_directory)).unwrap(), config.tumor_id)
}

/// OrangeAlgo.buildPlots
fn build_plots(config: &OrangeConfig, purple: &PurpleRecord, linx_data: &linx::LinxData) -> Result<OrangePlots> {
    let pm = PlotManager::new(&config.output_dir);
    pm.create_plot_directory()?;

    // PurpleDataLoader.addPurplePlots: a missing plot ends the run (Java: "failed to find required plots", exit 1)
    let required = |suffix: &str| -> Result<String> {
        match pm.process_plot_file(Some(&purple_plot(config, suffix)))? {
            Some(p) => Ok(p),
            None => {
                hmf_common::warn!("failed to find required plots: java.lang.NullPointerException: purple plot {suffix}");
                std::process::exit(1);
            }
        }
    };
    let final_circos = required(".circos.png")?;
    let input_circos = required(".input.png")?;
    let copy_number = required(".copynumber.png")?;
    let clonality = required(".somatic.clonality.png")?;
    let purity_range = required(".purity.range.png")?;
    let minor_allele = required(".map.png")?;
    let somatic_cn = required(".somatic.png")?;
    let rainfall = required(".somatic.rainfall.png")?;

    for v in &purple.somatic_variants {
        if let Some(p) = &v.plot_filename {
            pm.process_plot_file(Some(p))?;
        }
    }

    let mut linx_driver_plots = Vec::new();
    for (plot, title) in linx_data.reportable_event_plots.iter().zip(&linx_data.reportable_event_plot_titles) {
        let filename = pm.process_plot_file(Some(plot))?.ok_or_else(|| anyhow::anyhow!("java.lang.NullPointerException: plot filename"))?;
        linx_driver_plots.push(OrangePlot { filename, title: title.clone() });
    }

    // MiscToolFiles.generateQseeVisPlot: checkAddDirSeparator(null) + id gives "null<id>..." without a directory
    let qsee_source = format!("{}{}.qsee.vis.report.png", config.qsee_directory.as_deref().unwrap_or("null"), config.tumor_id);
    let q_see_plot = if std::path::Path::new(&qsee_source).exists() { pm.process_plot_file(Some(&qsee_source))? } else { None };

    let cuppa_source = config.reference_id.as_ref().map(|_| {
        tool_results::cuppa_vis_plot_filename(config.cuppa_dir.as_deref().unwrap_or("null"), &config.tumor_id)
    });
    // CuppaPredictions.generateVisPlotFilename(null dir): "null/<id>.cuppa.vis.png" (checkAddDirSeparator of null is null)
    let cuppa_source = cuppa_source.map(|s| if config.cuppa_dir.is_none() { format!("null{}.cuppa.vis.png", config.tumor_id) } else { s });
    let cuppa_summary_plot = pm.process_plot_file(cuppa_source.as_deref())?;

    Ok(OrangePlots {
        purple_input_circos_plot: input_circos,
        purple_final_circos_plot: final_circos,
        purple_clonality_plot: clonality,
        purple_copy_number_plot: copy_number,
        purple_minor_allele_map_plot: minor_allele,
        purple_variant_copy_number_plot: somatic_cn,
        purple_purity_range_plot: purity_range,
        purple_rainfall_plot: rainfall,
        linx_driver_plots,
        cuppa_summary_plot,
        q_see_plot,
    })
}

/// OrangeAlgo.run
pub fn run(config: &OrangeConfig) -> Result<OrangeRecord> {
    hmf_common::info!("building report data");
    let t0 = std::time::Instant::now();
    let lap = |what: &str| hmf_common::debug!("time: {what} done at {:.3} s", t0.elapsed().as_secs_f64());
    let primary_tumor = others::configured_primary_tumor(config)?;
    let pipeline_version = pipeline_version(config)?;

    let purple_data = purple::load(config)?;
    lap("PURPLE data");
    hmf_common::debug!(" loaded {} somatic driver catalog entries", purple_data.somatic_drivers.len());
    hmf_common::debug!(" loaded {} somatic variants", purple_data.somatic_variants.len());
    let linx_data = linx::load(config)?;
    lap("LINX data");
    let chord = others::load_chord(config)?;
    let lilac = others::load_lilac(config)?;
    let virus = others::load_virus(config)?;
    let cuppa = others::load_cuppa(config)?;
    let peach = others::load_peach(config)?;
    let sigs = others::load_sigs(config)?;
    let isofox_data = match &config.isofox_dir {
        None => None,
        Some(dir) => match isofox::check_rna_vs_tumor_id(config.rna_sample_id.as_deref(), &config.tumor_id, dir) {
            None => {
                hmf_common::warn!(" invalid RNA sample ID");
                None
            }
            Some(id) => Some(isofox::load(&id, dir)?),
        },
    };

    lap("other tools' data");
    let cyto_bands = CytoBands::new(config.ref_genome_v37);
    let linx = linx::interpret(&linx_data, isofox_data.as_ref(), &cyto_bands)?;
    let purple = purple::interpret(&purple_data, &linx_data, isofox_data.as_ref())?;
    let immune_escape = others::immune_escape(&purple, &linx);
    let isofox = match &isofox_data {
        Some(d) => Some(isofox::interpret(d, &linx.fusions)?),
        None => None,
    };
    lap("interpretation");
    let plots = build_plots(config, &purple, &linx_data)?;
    lap("plots");

    Ok(OrangeRecord {
        sample_id: config.tumor_id.clone(),
        reference_id: config.reference_id.clone(),
        // OR6: Java never sets the RNA sample id; corrected, the configured one
        rna_id: if crate::compat::java("OR6") { None } else { config.rna_sample_id.clone() },
        sampling_date: config.sampling_date.clone(),
        experiment_type: config.run_type.clone(),
        configured_primary_tumor: primary_tumor,
        ref_genome_version: if config.ref_genome_v37 { "V37" } else { "V38" }.to_string(),
        pipeline_version: Some(pipeline_version),
        panel_name: config.panel_name.clone(),
        sequencing_type: config.sequencing_type.clone(),
        purple,
        linx,
        isofox,
        lilac,
        immune_escape,
        viridian: virus,
        chord,
        cuppa,
        peach: Some(peach),
        sig_allocations: sigs,
        plots,
    })
}
