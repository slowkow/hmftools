//! LINX's part of the report: LinxDataLoader, LinxInterpreter and LinxBreakendInterpreter (orange/algo/linx).

use anyhow::{anyhow, bail, Result};
use hmf_common::cytobands::CytoBands;
use hmf_common::java_hash::JavaHashSet;
use hmf_common::linx_files::{self as lf, LinxBreakend, LinxCluster, LinxDriver, LinxFusion, LinxSvAnnotation};
use hmf_common::purple_files::{self, DriverCatalog};

use crate::compat::java;
use crate::config::OrangeConfig;
use crate::isofox::IsofoxData;
use crate::model;
use crate::purple::list_dir;

pub struct LinxData {
    /// the SV annotations of the reportable breakends
    pub somatic_sv_annotations: Vec<LinxSvAnnotation>,
    pub somatic_driver_data: Vec<LinxDriver>,
    /// the LINX driver catalog's disruptions
    pub somatic_drivers: Vec<DriverCatalog>,
    pub fusions: Vec<LinxFusion>,
    pub somatic_breakends: Vec<LinxBreakend>,
    pub germline_drivers: Option<Vec<DriverCatalog>>,
    pub germline_breakends: Option<Vec<LinxBreakend>>,
    pub germline_sv_annotations: Option<Vec<LinxSvAnnotation>>,
    pub reportable_event_plots: Vec<String>,
    pub reportable_event_plot_titles: Vec<String>,
}

/// LinxDataLoader.loadBreakends: reported breakends, and those of a (reported) fusion
fn load_breakends(path: &str, fusions: &[LinxFusion]) -> Result<Vec<LinxBreakend>> {
    Ok(lf::read_breakends(path)?
        .into_iter()
        .filter(|b| {
            b.reported_status == "REPORTED"
                || fusions.iter().any(|f| f.five_prime_breakend_id == b.id || f.three_prime_breakend_id == b.id)
        })
        .collect())
}

fn restrict_to_matching_breakends(annotations: &mut Vec<LinxSvAnnotation>, breakends: &[LinxBreakend]) {
    annotations.retain(|a| breakends.iter().any(|b| b.sv_id == a.sv_id));
}

/// LinxDataLoader.formPlotTitle: "Cluster(<id>: <resolved type>) SVs(<count>) Genes(<sorted genes>)", or "" with
/// a warning where the file name or cluster does not resolve
fn form_plot_title(file: &str, svs: &[LinxSvAnnotation], breakends: &[LinxBreakend], clusters: &[LinxCluster]) -> String {
    let parse = || -> Option<i32> {
        let idx = file.find(lf::SV_VIS_CLUSTER_PREFIX)?;
        let remainder = &file[idx + lf::SV_VIS_CLUSTER_PREFIX.len()..];
        let end = remainder.find('.')?;
        if end < 1 {
            return None;
        }
        hmf_common::java_text::parse_int(&remainder[1..end])
    };
    let Some(cluster_id) = parse() else {
        hmf_common::warn!("Linx plot filename({file}) failed to extract cluster ID");
        return String::new();
    };
    let Some(cluster) = clusters.iter().find(|c| c.cluster_id == cluster_id) else {
        hmf_common::warn!("Linx plot filename({file}) cluster({cluster_id}) not found");
        return String::new();
    };
    let mut genes: JavaHashSet<String> = JavaHashSet::new();
    let mut sv_count = 0;
    for sv in svs.iter().filter(|s| s.cluster_id == cluster_id) {
        sv_count += 1;
        for b in breakends.iter().filter(|b| b.sv_id == sv.sv_id) {
            genes.insert(b.gene.clone());
        }
    }
    let mut list: Vec<String> = genes.to_vec().into_iter().filter(|g| !g.is_empty()).collect();
    list.sort_by(|a, b| hmf_common::java_text::compare(a, b));
    format!("Cluster({cluster_id}: {}) SVs({sv_count}) Genes({})", cluster.resolved_type, list.join(", "))
}

/// LinxDataLoader.load
pub fn load(config: &OrangeConfig) -> Result<LinxData> {
    let dir = &config.linx_somatic_data_directory;
    let tumor = &config.tumor_id;
    let somatic_driver_data = lf::read_drivers(&lf::driver_filename(dir, tumor))?;
    let fusions: Vec<LinxFusion> = lf::read_fusions(&lf::fusion_filename(dir, tumor))?.into_iter().filter(|f| f.reported).collect();
    let somatic_breakends = load_breakends(&lf::breakend_filename(dir, tumor, false), &fusions)?;
    let somatic_clusters = lf::read_clusters(&lf::cluster_filename(dir, tumor, false))?;
    let somatic_drivers: Vec<DriverCatalog> = purple_files::read_driver_catalog(&lf::driver_catalog_filename(dir, tumor, true))?
        .into_iter()
        .filter(|d| ["HOM_DUP_DISRUPTION", "HOM_DEL_DISRUPTION", "DISRUPTION"].contains(&d.driver.as_str()))
        .collect();
    let all_somatic_svs = lf::read_sv_annotations(&lf::sv_annotation_filename(dir, tumor, false))?;
    let mut somatic_sv_annotations = all_somatic_svs.clone();
    restrict_to_matching_breakends(&mut somatic_sv_annotations, &somatic_breakends);

    let mut germline_drivers = None;
    let mut germline_breakends = None;
    let mut germline_sv_annotations = None;
    if let Some(gdir) = &config.linx_germline_data_directory {
        germline_drivers = Some(
            purple_files::read_driver_catalog(&lf::driver_catalog_filename(gdir, tumor, false))?
                .into_iter()
                .filter(|d| d.driver == "GERMLINE_DISRUPTION")
                .collect::<Vec<_>>(),
        );
        let mut svs = lf::read_sv_annotations(&lf::sv_annotation_filename(gdir, tumor, true))?;
        let breakends = load_breakends(&lf::breakend_filename(gdir, tumor, true), &[])?;
        restrict_to_matching_breakends(&mut svs, &breakends);
        // read as Java does (the reportable disruptions are not part of the report)
        lf::read_germline_disruptions(&lf::germline_disruption_filename(gdir, tumor))?;
        germline_sv_annotations = Some(svs);
        germline_breakends = Some(breakends);
    }

    let mut reportable_event_plots = Vec::new();
    let mut reportable_event_plot_titles = Vec::new();
    if let Some(plot_dir) = &config.linx_plot_directory {
        for file in list_dir(plot_dir)? {
            if file.contains(lf::SV_VIS_CLUSTER_PREFIX) && file.ends_with(".png") {
                reportable_event_plots.push(format!("{plot_dir}{file}"));
                reportable_event_plot_titles.push(form_plot_title(&file, &all_somatic_svs, &somatic_breakends, &somatic_clusters));
            }
        }
        hmf_common::debug!(" loaded {} Linx plots from {}", reportable_event_plots.len(), plot_dir);
    }

    Ok(LinxData {
        somatic_sv_annotations,
        somatic_driver_data,
        somatic_drivers,
        fusions,
        somatic_breakends,
        germline_drivers,
        germline_breakends,
        germline_sv_annotations,
        reportable_event_plots,
        reportable_event_plot_titles,
    })
}

/// LinxInterpreter.findReportableLinxPlot: the first plot whose name contains ".cluster-<id>" (Java, OR4), or
/// whose name has exactly that cluster id (".cluster-<id>." or ".cluster-<id>_")
pub fn find_reportable_linx_plot(plots: &[String], cluster_id: i32) -> Option<String> {
    let prefix = lf::vis_plot_cluster_prefix(&cluster_id.to_string());
    if java("OR4") {
        return plots.iter().find(|p| p.contains(&prefix)).cloned();
    }
    plots
        .iter()
        .find(|p| {
            let name = p.rsplit('/').next().unwrap_or(p);
            name.match_indices(&prefix).any(|(i, _)| {
                let next = name[i + prefix.len()..].chars().next();
                !matches!(next, Some(c) if c.is_ascii_digit())
            })
        })
        .cloned()
}

/// LinxBreakendInterpreter.junctionCopyNumber
pub(crate) fn junction_copy_number(sv: &LinxSvAnnotation) -> f64 {
    0.5 * (sv.junction_copy_number_min + sv.junction_copy_number_max)
}

/// Collectors.toMap(svId): a repeated id is an IllegalStateException
fn sv_map(svs: &[LinxSvAnnotation]) -> Result<std::collections::HashMap<i32, &LinxSvAnnotation>> {
    let mut m = std::collections::HashMap::new();
    for s in svs {
        if m.insert(s.sv_id, s).is_some() {
            bail!("java.lang.IllegalStateException: Duplicate key {}", s.sv_id);
        }
    }
    Ok(m)
}

fn build_breakend(
    b: &LinxBreakend,
    svs: &std::collections::HashMap<i32, &LinxSvAnnotation>,
    drivers: &[DriverCatalog],
    cyto_bands: &CytoBands,
    plots: &[String],
) -> Result<model::LinxBreakend> {
    let sv = svs.get(&b.sv_id).ok_or_else(|| anyhow!("java.lang.NullPointerException: no SV annotation for svId {}", b.sv_id))?;
    let coords = if b.is_start { &sv.coords_start } else { &sv.coords_end };
    let chromosome = lf::chromosome_from_coords(coords);
    let position = lf::position_from_coords(coords)?;
    let orientation = lf::orientation_from_coords(coords)?;
    let band = cyto_bands.cyto_band_name(&chromosome, position).map(String::from);
    let mut driver_type = "DISRUPTION".to_string();
    let mut driver_likelihood = 0.0;
    if let Some(d) = drivers.iter().find(|d| d.gene == b.gene) {
        if d.driver != "GERMLINE_DISRUPTION" {
            if !["DISRUPTION", "HOM_DEL_DISRUPTION", "HOM_DUP_DISRUPTION"].contains(&d.driver.as_str()) {
                bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.datamodel.linx.LinxDriverType.{}", d.driver);
            }
            driver_type = d.driver.clone();
        }
        driver_likelihood = d.driver_likelihood;
    }
    Ok(model::LinxBreakend {
        id: b.id,
        sv_id: b.sv_id,
        cluster_id: sv.cluster_id,
        gene: b.gene.clone(),
        position,
        orientation: orientation as i32,
        chromosome,
        chromosome_band: band.ok_or_else(|| anyhow!("java.lang.NullPointerException: chromosomeBand"))?,
        transcript: b.transcript_id.clone(),
        is_canonical: b.canonical,
        gene_orientation: if b.gene_orientation == lf::BREAKEND_ORIENTATION_UPSTREAM { "UPSTREAM" } else { "DOWNSTREAM" }.to_string(),
        disruptive: b.disruptive,
        reported_status: b.reported_status.clone(),
        undisrupted_copy_number: b.undisrupted_copy_number,
        breakend_type: sv.sv_type.clone(),
        region_type: b.region_type.clone(),
        coding_type: b.coding_type.clone(),
        next_splice_exon_rank: b.next_splice_exon_rank,
        exon_up: b.exon_up,
        exon_down: b.exon_down,
        junction_copy_number: junction_copy_number(sv),
        driver_type,
        driver_likelihood,
        plot_filename: find_reportable_linx_plot(plots, sv.cluster_id),
    })
}

/// LinxInterpreter.buildContextStr
fn context_str(region_type: &str, known_type: &str, fused_exon: i32) -> Result<String> {
    Ok(match region_type {
        "UPSTREAM" => "Promoter Region".to_string(),
        "DOWNSTREAM" => "Post-coding".to_string(),
        "IG" => "IG".to_string(),
        "EXONIC" | "INTRONIC" => format!("Exon {fused_exon}"),
        "UNKNOWN" => {
            if known_type == "PROMISCUOUS_ENHANCER_TARGET" {
                "Unknown".to_string()
            } else {
                "ERROR: UNKNOWN".to_string()
            }
        }
        other => bail!("java.lang.IllegalStateException: TranscriptRegionType not supported in determination of fusion context: {other}"),
    })
}

const LINX_FUSION_TYPES: [&str; 9] = [
    "NONE",
    "KNOWN_PAIR",
    "PROMISCUOUS_5",
    "PROMISCUOUS_3",
    "PROMISCUOUS_BOTH",
    "ENHANCER_KNOWN_PAIR",
    "ENHANCER_PROMISCUOUS",
    "EXON_DEL_DUP",
    "PROMISCUOUS_ENHANCER_TARGET",
];

fn find_rna_support(fusion: &LinxFusion, isofox: Option<&IsofoxData>) -> Option<model::AllelicDepth> {
    let isofox = isofox?;
    let (mut coverage, mut fragments, mut matched) = (0i32, 0i32, false);
    for r in isofox.fusions.iter().filter(|r| r.name == fusion.name) {
        matched = true;
        let frags = r.split_fragments + r.discordant_frags + r.realigned_frags;
        let average_depth = (hmf_common::java_math::java_round(r.depth_up.wrapping_add(r.depth_down) as f64 * 0.5) as i32).min(frags);
        coverage += average_depth;
        fragments += frags;
    }
    if !matched {
        return None;
    }
    Some(model::AllelicDepth { total_read_count: coverage, allele_read_count: fragments })
}

fn build_fusions(linx: &LinxData, isofox: Option<&IsofoxData>) -> Result<Vec<model::LinxFusion>> {
    let mut out = Vec::new();
    for f in &linx.fusions {
        let up = linx.somatic_breakends.iter().find(|b| b.id == f.five_prime_breakend_id);
        let down = linx.somatic_breakends.iter().find(|b| b.id == f.three_prime_breakend_id);
        let (Some(up), Some(down)) = (up, down) else {
            hmf_common::error!("fusion({}) missing corresponding breakends", f.name);
            continue;
        };
        let Some(sv) = linx.somatic_sv_annotations.iter().find(|s| s.sv_id == up.sv_id) else {
            hmf_common::error!("fusion({}) missing corresponding svAnnotation", f.name);
            continue;
        };
        let plot_filename = find_reportable_linx_plot(&linx.reportable_event_plots, sv.cluster_id);
        if !LINX_FUSION_TYPES.contains(&f.reported_type.as_str()) {
            bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.datamodel.linx.LinxFusionType.{}", f.reported_type);
        }
        let fusion_type = f.reported_type.as_str();
        let mut gene_up = up.gene.clone();
        let gene_down = down.gene.clone();
        if fusion_type == "ENHANCER_KNOWN_PAIR" || fusion_type == "ENHANCER_PROMISCUOUS" {
            gene_up = up.transcript_id.clone();
        }
        let context_up = context_str(&up.region_type, fusion_type, f.fused_exon_up)?;
        let context_down = context_str(&down.region_type, fusion_type, f.fused_exon_down)?;
        // OR8: Java averages the breakends' undisrupted copy numbers; corrected, the junction copy number of
        // the fusion's SV (as for the breakends)
        let jcn = if java("OR8") { (up.undisrupted_copy_number + down.undisrupted_copy_number) * 0.5 } else { junction_copy_number(sv) };
        out.push(model::LinxFusion {
            gene_up,
            cluster_id: sv.cluster_id,
            context_up,
            transcript_up: up.transcript_id.clone(),
            gene_down,
            context_down,
            transcript_down: down.transcript_id.clone(),
            reported_type: fusion_type.to_string(),
            phased: f.phased.clone(),
            driver_interpretation: match f.likelihood.as_str() {
                "HIGH" => "HIGH",
                "LOW" => "LOW",
                _ => "UNKNOWN",
            }
            .to_string(),
            fused_exon_up: f.fused_exon_up,
            fused_exon_down: f.fused_exon_down,
            chain_links: f.chain_links,
            chain_terminated: f.chain_terminated,
            domains_kept: f.domains_kept.clone(),
            domains_lost: f.domains_lost.clone(),
            junction_copy_number: jcn,
            rna_support: find_rna_support(f, isofox),
            plot_filename,
        });
    }
    Ok(out)
}

/// LinxInterpreter.interpret
pub fn interpret(linx: &LinxData, isofox: Option<&IsofoxData>, cyto_bands: &CytoBands) -> Result<model::LinxRecord> {
    let svs = sv_map(&linx.somatic_sv_annotations)?;
    let mut somatic_breakends = Vec::new();
    for b in linx.somatic_breakends.iter().filter(|b| b.reported_status == "REPORTED") {
        somatic_breakends.push(build_breakend(b, &svs, &linx.somatic_drivers, cyto_bands, &linx.reportable_event_plots)?);
    }
    let germline_breakends = match &linx.germline_sv_annotations {
        None => None,
        Some(gsvs) => {
            let map = sv_map(gsvs)?;
            let drivers = linx.germline_drivers.as_deref().unwrap_or(&[]);
            let mut v = Vec::new();
            for b in linx.germline_breakends.as_deref().unwrap_or(&[]).iter().filter(|b| b.reported_status == "REPORTED") {
                v.push(build_breakend(b, &map, drivers, cyto_bands, &[])?);
            }
            Some(v)
        }
    };
    Ok(model::LinxRecord { fusions: build_fusions(linx, isofox)?, somatic_breakends, germline_breakends })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plot_lookup() {
        let plots = vec!["d/S.cluster-40.001.png".to_string(), "d/S.cluster-4.003.png".to_string()];
        crate::compat::with_java(true, || {
            assert_eq!(find_reportable_linx_plot(&plots, 4).as_deref(), Some("d/S.cluster-40.001.png"));
        });
        crate::compat::with_java(false, || {
            assert_eq!(find_reportable_linx_plot(&plots, 4).as_deref(), Some("d/S.cluster-4.003.png"));
            assert_eq!(find_reportable_linx_plot(&plots, 40).as_deref(), Some("d/S.cluster-40.001.png"));
            assert_eq!(find_reportable_linx_plot(&plots, 5), None);
        });
    }

    #[test]
    fn contexts() {
        assert_eq!(context_str("UPSTREAM", "KNOWN_PAIR", 1).unwrap(), "Promoter Region");
        assert_eq!(context_str("EXONIC", "KNOWN_PAIR", 3).unwrap(), "Exon 3");
        assert_eq!(context_str("UNKNOWN", "PROMISCUOUS_ENHANCER_TARGET", 0).unwrap(), "Unknown");
        assert_eq!(context_str("UNKNOWN", "KNOWN_PAIR", 0).unwrap(), "ERROR: UNKNOWN");
    }
}
