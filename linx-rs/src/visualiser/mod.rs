//! LINX's visualiser (`com.hartwig.hmftools.linx.visualiser.SvVisualiser`): per cluster, chromosome and
//! driver gene, the Circos data and configuration files (as Java writes them) and the plot, drawn natively:
//! the Circos image by hmf-plot-rs's Circos engine, the chromosome range strip (R's Gviz in Java) and the
//! fusion panels (R's ggplot2 in Java), composed as Java's R scripts compose them.
//!
//! Findings (corrected by default; `-java_compat` keeps Java's behaviour):
//! - LV1: VisProteinDomain.name() tests the 5' UTR twice and never the 3' UTR.
//! - LV2: with a fusion panel, the chromosome range script scales the image (Circos plus fusions) into a
//!   square as tall as the image is wide, which shrinks and narrows it.
//! - LV3: R draws the fusion plot's protein domain legend on a black row, hiding the domain names.
//! - LV4: links (and their connectors) whose frame is past the last are drawn faint grey in the final image:
//!   all links of plots without segment tracks, and those of the last segment.

pub mod circos_data;
pub mod data;
pub mod files;
pub mod fusion;
pub mod genome;
pub mod plot;
#[cfg(test)]
mod tests;

use std::collections::HashSet;
use std::sync::Mutex;

use anyhow::{bail, Context, Result};
use hmf_common::args::Args;
use hmf_common::ensembl::EnsemblDataCache;
use hmf_common::regions::RefGenomeVersion;
use hmf_common::{debug, error, info, warn};

use circos_data::*;
use data::*;
use files::*;
use genome::*;

pub const MAIN_CLASS: &str = "com.hartwig.hmftools.linx.visualiser.SvVisualiser";

const FLAGS: [&str; 21] = [
    "load_cohort_files",
    "germline",
    "restrict_cluster_by_gene",
    "include_fragile_sites",
    "include_line_elements",
    "plot_reportable",
    "plot_cluster_genes",
    "debug",
    "interpolate_cna_positions",
    "interpolate_exon_positions",
    "exact_position",
    "show_sv_id",
    "step",
    "log_debug",
    "no_plots",
    "no_data_files",
    "perf_debug",
    "version",
    "help",
    "h",
    // `-circos` given without a path (the pipeline's `-circos $(which circos)` where Circos is not installed)
    "circos_flag",
];

const ITEMS: [&str; 38] = [
    "sample",
    "vis_file_dir",
    "amber_dir",
    "cobalt_dir",
    "purple_dir",
    "plot_out",
    "data_out",
    "circos",
    "ref_genome_version",
    "ensembl_data_dir",
    "gene",
    "clusterId",
    "chainId",
    "chromosome",
    "specific_regions",
    "threads",
    "log_level",
    "outer_radius",
    "inner_radius",
    "gap_radius",
    "exon_rank_radius",
    "max_plot_svs",
    "fusion_height",
    "fusion_legend_rows",
    "fusion_legend_height_per_row",
    "chr_range_height",
    "chr_range_columns",
    "gene_relative_size",
    "segment_relative_size",
    "cna_relative_size",
    "min_label_size",
    "max_label_size",
    "max_distance_labels",
    "max_position_labels",
    "max_gene_characters",
    "min_line_size",
    "max_line_size",
    "glyph_size",
];

/// CircosConfig and the drawing switches of VisualiserConfig
#[derive(Clone, Debug)]
pub struct VisOptions {
    pub fusion_legend_rows: i32,
    pub fusion_legend_height_per_row: i32,
    pub fusion_height: i32,
    pub chromosome_range_height: i32,
    pub chromosome_range_columns: i32,
    pub max_plot_sv_count: i32,
    pub min_label_size: i32,
    pub max_label_size: i32,
    pub max_gene_characters: i32,
    pub max_position_labels: i32,
    pub max_distance_labels: i32,
    pub gene_relative_size: f64,
    pub segment_relative_size: f64,
    pub copy_number_relative_size: f64,
    pub min_line_size: i32,
    pub max_line_size: i32,
    pub glyph_size: i32,
    pub outer_radius: f64,
    pub inner_radius: f64,
    pub gap_radius: f64,
    pub exon_rank_radius: f64,
    pub interpolate_copy_number_positions: bool,
    pub interpolate_exon_positions: bool,
    pub step: bool,
    pub exact_position: bool,
    pub show_sv_id: bool,
}

impl VisOptions {
    fn from_args(a: &Args) -> Result<VisOptions> {
        Ok(VisOptions {
            fusion_legend_rows: a.integer("fusion_legend_rows", 1)?,
            fusion_legend_height_per_row: a.integer("fusion_legend_height_per_row", 35)?,
            fusion_height: a.integer("fusion_height", 250)?,
            chromosome_range_height: a.integer("chr_range_height", 150)?,
            chromosome_range_columns: a.integer("chr_range_columns", 6)?,
            max_plot_sv_count: a.integer("max_plot_svs", 2000)?,
            min_label_size: a.integer("min_label_size", 35)?,
            max_label_size: a.integer("max_label_size", 40)?,
            max_gene_characters: a.integer("max_gene_characters", 5)?,
            max_position_labels: a.integer("max_position_labels", 60)?,
            max_distance_labels: a.integer("max_distance_labels", 100)?,
            gene_relative_size: a.decimal("gene_relative_size", 0.3)?,
            segment_relative_size: a.decimal("segment_relative_size", 1.0)?,
            copy_number_relative_size: a.decimal("cna_relative_size", 2.0)?,
            min_line_size: a.integer("min_line_size", 1)?,
            max_line_size: a.integer("max_line_size", 12)?,
            glyph_size: a.integer("glyph_size", 20)?,
            outer_radius: a.decimal("outer_radius", 0.88)?,
            inner_radius: a.decimal("inner_radius", 0.20)?,
            gap_radius: a.decimal("gap_radius", 0.025)?,
            exon_rank_radius: a.decimal("exon_rank_radius", 0.025)?,
            interpolate_copy_number_positions: a.has_flag("interpolate_cna_positions"),
            interpolate_exon_positions: a.has_flag("interpolate_exon_positions"),
            step: a.has_flag("step"),
            exact_position: a.has_flag("exact_position"),
            show_sv_id: a.has_flag("show_sv_id"),
        })
    }

    /// CircosConfig.labelSize
    pub fn label_size(&self, count: i64) -> i64 {
        if count > self.max_distance_labels as i64 {
            return self.min_label_size as i64;
        }
        java_round(self.max_label_size as f64 - 1.0 * count as f64 * (self.max_label_size - self.min_label_size) as f64 / self.max_distance_labels as f64)
    }

    pub fn exceeds_max_plot_sv_count(&self, n: usize) -> bool {
        self.max_plot_sv_count > 0 && n as i64 > self.max_plot_sv_count as i64
    }

    /// CircosConfig.isValid
    fn is_valid(&self) -> bool {
        let checks: [(bool, &str); 4] = [
            (self.min_label_size <= 0, "Parameter min_label_size should be > 0"),
            (self.max_label_size < self.min_label_size, "Parameter max_label_size should be > min_label_size"),
            (self.min_line_size <= 0, "Parameter min_line_size should be > 0"),
            (self.max_line_size < self.min_line_size, "Parameter max_line_size should be > min_line_size"),
        ];
        for (bad, msg) in checks {
            if bad {
                error!("{msg}");
                return false;
            }
        }
        for (name, v) in [("outer_radius", self.outer_radius), ("inner_radius", self.inner_radius), ("gap_radius", self.gap_radius), ("exon_rank_radius", self.exon_rank_radius)] {
            if v - 0.0 < -1e-10 || v - 1.0 > 1e-10 {
                error!("radial parameter({name}) value({v}), must be betwee 0 - 1");
                return false;
            }
        }
        for (name, v) in [("gene_relative_size", self.gene_relative_size), ("segment_relative_size", self.segment_relative_size), ("cna_relative_size", self.copy_number_relative_size)] {
            if v - 0.0 <= -1e-10 || v.is_nan() {
                error!("relative parameter({name}) value({v}), must be greater than 0");
                return false;
            }
        }
        true
    }
}

/// VisualiserConfig
pub struct VisConfig {
    pub sample: String,
    pub sample_data_dir: String,
    pub amber_dir: Option<String>,
    pub cobalt_dir: Option<String>,
    pub purple_dir: Option<String>,
    pub is_germline: bool,
    pub output_conf_path: String,
    pub output_plot_path: String,
    pub ensembl_data_dir: Option<String>,
    pub v37: bool,
    pub threads: usize,
    pub debug: bool,
    pub include_fragile_sites: bool,
    pub include_line_elements: bool,
    pub cluster_ids: Vec<i32>,
    pub chain_ids: Vec<i32>,
    pub chromosomes: Vec<String>,
    pub genes: Vec<String>,
    pub plot_reportable: bool,
    pub plot_cluster_genes: bool,
    pub restrict_cluster_by_gene: bool,
    pub write_plots: bool,
    pub write_data: bool,
    pub perf_debug: bool,
}

impl VisConfig {
    /// the genes in the iteration order of Java's HashSet
    fn genes_java_order(&self) -> Vec<String> {
        let mut set: hmf_common::java_hash::JavaHashSet<String> = hmf_common::java_hash::JavaHashSet::new();
        for g in &self.genes {
            set.insert(g.clone());
        }
        set.iter().cloned().collect()
    }
}

fn dir_sep(d: &str) -> String {
    if d.ends_with('/') { d.to_string() } else { format!("{d}/") }
}

fn parse_ints(s: &str, what: &str) -> Result<Vec<i32>> {
    s.split(',').map(|x| x.parse::<i32>().with_context(|| format!("{what}: For input string: \"{x}\""))).collect()
}

/// The chromosome names of HumanChromosome.values(), as toString gives them
const ALL_CHROMOSOMES: [&str; 24] = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"];

impl VisConfig {
    fn new(a: &Args) -> Result<VisConfig> {
        let sample = a.value("sample").context("-sample is required")?.to_string();
        let sample_data_dir = dir_sep(a.value("vis_file_dir").context("-vis_file_dir is required")?);
        if a.has_flag("load_cohort_files") {
            bail!("-load_cohort_files (cohort vis files) is not ported");
        }
        if a.has_value("specific_regions") {
            bail!("-specific_regions is not ported");
        }
        let v37 = match a.value_or("ref_genome_version", "V37") {
            "38" | "V38" | "HG38" | "GRCh38" => false,
            "37" | "V37" | "HG19" | "GRCh37" => true,
            v => bail!("invalid ref genome version: {v}"),
        };
        let chromosomes = match a.value("chromosome") {
            Some("All") => ALL_CHROMOSOMES.iter().map(|s| s.to_string()).collect(),
            Some(v) => v.split(',').map(|s| s.to_string()).collect(),
            None => Vec::new(),
        };
        let mut genes: Vec<String> = Vec::new();
        if let Some(g) = a.value("gene") {
            let delim = if g.contains(';') { ';' } else { ',' };
            for x in g.split(delim) {
                if !genes.iter().any(|y| y == x) {
                    genes.push(x.to_string());
                }
            }
        }
        let output_plot_path = dir_sep(a.value("plot_out").map(|s| s.to_string()).unwrap_or(format!("{sample_data_dir}plot/")).as_str());
        let output_conf_path = dir_sep(a.value("data_out").map(|s| s.to_string()).unwrap_or(format!("{sample_data_dir}data/")).as_str());
        std::fs::create_dir_all(&output_plot_path).with_context(|| format!("Unable to write to plot directory {output_plot_path}"))?;
        std::fs::create_dir_all(&output_conf_path).with_context(|| format!("Unable to write to data directory {output_conf_path}"))?;
        Ok(VisConfig {
            sample,
            sample_data_dir,
            amber_dir: a.value("amber_dir").map(dir_sep),
            cobalt_dir: a.value("cobalt_dir").map(dir_sep),
            purple_dir: a.value("purple_dir").map(dir_sep),
            is_germline: a.has_flag("germline"),
            output_conf_path,
            output_plot_path,
            ensembl_data_dir: a.value("ensembl_data_dir").map(|s| s.to_string()),
            v37,
            threads: a.integer("threads", 1)?.max(1) as usize,
            debug: a.has_flag("debug"),
            include_fragile_sites: a.has_flag("include_fragile_sites"),
            include_line_elements: a.has_flag("include_line_elements"),
            cluster_ids: a.value("clusterId").map(|v| parse_ints(v, "clusterId")).transpose()?.unwrap_or_default(),
            chain_ids: a.value("chainId").map(|v| parse_ints(v, "chainId")).transpose()?.unwrap_or_default(),
            chromosomes,
            genes,
            plot_reportable: a.has_flag("plot_reportable"),
            plot_cluster_genes: a.has_flag("plot_cluster_genes"),
            restrict_cluster_by_gene: a.has_flag("restrict_cluster_by_gene"),
            write_plots: !a.has_flag("no_plots"),
            write_data: !a.has_flag("no_data_files"),
            perf_debug: a.has_flag("perf_debug"),
        })
    }
}

/// SampleData
pub struct SampleData {
    pub sv_data: Vec<VisSvData>,
    pub segments: Vec<VisSegment>,
    pub copy_numbers: Vec<VisCopyNumber>,
    pub protein_domains: Vec<VisProteinDomain>,
    pub fusions: Vec<VisFusion>,
    pub exons: Vec<VisGeneExon>,
    pub amber: Vec<ChrPoints>,
    pub cobalt: Vec<ChrPoints>,
    pub purple_segments: Vec<PurpleSegment>,
    pub purple_gene_cnvs: Vec<(String, String)>,
}

fn vis_file(dir: &str, sample: &str, germline: bool, kind: &str) -> String {
    format!("{dir}{sample}{}{kind}", if germline { ".linx.germline." } else { ".linx." })
}

/// What SampleData does besides loading: None to exit with the code given (0: nothing to plot)
fn load_sample_data(config: &mut VisConfig) -> Result<Result<SampleData, i32>> {
    let s = config.sample.clone();
    let d = config.sample_data_dir.clone();
    info!("loading LINX visualiser data - sample({})", s);
    let mut germline = config.is_germline;
    if !germline && !std::path::Path::new(&vis_file(&d, &s, false, "vis_sv_data.tsv")).exists() && std::path::Path::new(&vis_file(&d, &s, true, "vis_sv_data.tsv")).exists() {
        germline = true;
    }
    let mut sv_data = read_vis_sv_data(&vis_file(&d, &s, germline, "vis_sv_data.tsv"))?;
    let mut sv_filtered = false;
    if !config.cluster_ids.is_empty() {
        sv_data.retain(|x| config.cluster_ids.contains(&x.cluster_id));
        sv_filtered = true;
        if !config.chain_ids.is_empty() {
            sv_data.retain(|x| config.chain_ids.contains(&x.chain_id));
        }
    }
    debug!("loaded {} LINX vis SV data entries", sv_data.len());
    let fusion_file = vis_file(&d, &s, germline, "vis_fusion.tsv");
    let fusions = if std::path::Path::new(&fusion_file).exists() { read_vis_fusions(&fusion_file)? } else { Vec::new() };
    let exons = read_vis_gene_exons(&vis_file(&d, &s, germline, "vis_gene_exon.tsv"))?;
    let segments = read_vis_segments(&vis_file(&d, &s, germline, "vis_segments.tsv"))?;
    let copy_numbers = read_vis_copy_numbers(&vis_file(&d, &s, germline, "vis_copy_number.tsv"))?;
    // VisProteinDomains.readProteinDomains: those of the fusions' transcripts
    let mut transcripts = HashSet::new();
    for f in &fusions {
        transcripts.insert(f.transcript_up.clone());
        transcripts.insert(f.transcript_down.clone());
    }
    let protein_domains: Vec<VisProteinDomain> =
        read_vis_protein_domains(&vis_file(&d, &s, germline, "vis_protein_domain.tsv"))?.into_iter().filter(|p| transcripts.contains(&p.transcript)).collect();

    if config.amber_dir.is_some() || config.cobalt_dir.is_some() || config.purple_dir.is_some() {
        info!("loading CNV data - sample({})", s);
    }
    let amber = match &config.amber_dir {
        Some(a) => read_amber_bafs(&format!("{a}{s}.amber.baf.tsv.gz"))?,
        None => Vec::new(),
    };
    let cobalt = match &config.cobalt_dir {
        Some(c) => read_cobalt_ratios(&format!("{c}{s}.cobalt.ratio.tsv.gz"))?,
        None => Vec::new(),
    };
    let (purple_segments, purple_gene_cnvs) = match &config.purple_dir {
        Some(p) => (read_purple_segments(&format!("{p}{s}.purple.segment.tsv"))?, read_purple_gene_cnvs(&format!("{p}{s}.purple.driver.catalog.somatic.tsv"))?),
        None => (Vec::new(), Vec::new()),
    };
    let clusters_given = !config.chain_ids.is_empty() || !config.cluster_ids.is_empty();
    let upstream = amber.iter().any(|c| !c.points.is_empty()) || cobalt.iter().any(|c| !c.points.is_empty()) || !purple_segments.is_empty();
    let mut missing = Vec::new();
    if copy_numbers.is_empty() {
        missing.push("copy numbers");
    }
    if segments.is_empty() {
        missing.push("segments");
    }
    if sv_data.is_empty() {
        missing.push("SV data");
    }
    if clusters_given && !missing.is_empty() {
        error!("sample({}) Cannot plot user specified cluster/chain IDs because LINX vis data was empty", s);
        return Ok(Err(1));
    }
    if !missing.is_empty() {
        if upstream {
            info!("proceeding to plotting CNV data only - the following LINX vis data are empty: {}", missing.join(", "));
        } else {
            info!("no plots to generate - sample({}) LINX vis data and AMBER/COBALT/PURPLE CNV data are all empty", s);
            return Ok(Err(0));
        }
    }

    // the SV annotations: genes of the given clusters, or the clusters of the given genes
    if config.plot_cluster_genes || config.restrict_cluster_by_gene {
        let f = hmf_common::linx_files::sv_annotation_filename(&d, &s, false);
        let fg = hmf_common::linx_files::sv_annotation_filename(&d, &s, true);
        let mut anns = hmf_common::linx_files::read_sv_annotations(if std::path::Path::new(&f).exists() { &f } else { &fg })?;
        if sv_filtered {
            anns.retain(|a| sv_data.iter().any(|y| y.sv_id == a.sv_id));
        }
        if config.plot_cluster_genes && !config.cluster_ids.is_empty() {
            for a in &anns {
                if config.cluster_ids.contains(&a.cluster_id) {
                    for g in [&a.gene_start, &a.gene_end] {
                        if !g.is_empty() {
                            for x in g.split(';') {
                                if !config.genes.iter().any(|y| y == x) {
                                    config.genes.push(x.to_string());
                                }
                            }
                        }
                    }
                }
            }
        } else if !config.genes.is_empty() && config.cluster_ids.is_empty() && config.restrict_cluster_by_gene {
            for a in &anns {
                if config.genes.iter().any(|x| *x == a.gene_start || *x == a.gene_end) && !config.cluster_ids.contains(&a.cluster_id) {
                    config.cluster_ids.push(a.cluster_id);
                }
            }
            for dr in hmf_common::linx_files::read_drivers(&hmf_common::linx_files::driver_filename(&d, &s))? {
                if config.genes.contains(&dr.gene) && !config.cluster_ids.contains(&dr.cluster_id) {
                    config.cluster_ids.push(dr.cluster_id);
                }
            }
        }
    }

    Ok(Ok(SampleData { sv_data, segments, copy_numbers, protein_domains, fusions, exons, amber, cobalt, purple_segments, purple_gene_cnvs }))
}

impl SampleData {
    /// findReportableClusters: clusters with a fusion, a reported disruption or a driver
    fn reportable_clusters(&self, config: &VisConfig) -> Vec<i32> {
        let mut ids: Vec<i32> = Vec::new();
        let add = |id: i32, ids: &mut Vec<i32>| {
            if !ids.contains(&id) {
                ids.push(id);
            }
        };
        for f in &self.fusions {
            add(f.cluster_id, &mut ids);
        }
        let d = &config.sample_data_dir;
        let s = &config.sample;
        let r = (|| -> Result<()> {
            let breakends = hmf_common::linx_files::read_breakends(&hmf_common::linx_files::breakend_filename(d, s, false))?;
            for b in breakends.iter().filter(|b| b.reported_status == "REPORTED") {
                if let Some(sv) = self.sv_data.iter().find(|x| x.sv_id == b.sv_id) {
                    add(sv.cluster_id, &mut ids);
                }
            }
            for dr in hmf_common::linx_files::read_drivers(&hmf_common::linx_files::driver_filename(d, s))? {
                if dr.cluster_id >= 0 {
                    add(dr.cluster_id, &mut ids);
                }
            }
            Ok(())
        })();
        if let Err(e) = r {
            error!("sample({}) could not read breakends or drivers: {}", s, e);
        }
        ids
    }
}

/// One plot to make (the arguments of submitFiltered)
pub struct PlotSpec {
    pub file_id: String,
    pub links: Vec<VisSvData>,
    pub segments: Vec<VisSegment>,
    pub exons: Vec<VisGeneExon>,
    pub protein_domains: Vec<VisProteinDomain>,
    pub fusions: Vec<VisFusion>,
    pub show_simple_sv_segments: bool,
    pub cluster_colours: bool,
}

struct Visualiser<'a> {
    config: &'a VisConfig,
    opts: &'a VisOptions,
    data: &'a SampleData,
    ensembl: Option<EnsemblDataCache>,
    coords: RefCoords,
    sites: KnownSites,
}

impl Visualiser<'_> {
    fn gene_chromosome(&self, gene: &str) -> Result<String> {
        let e = self.ensembl.as_ref().context("NullPointerException: no Ensembl data (-ensembl_data_dir)")?;
        Ok(e.gene_data_by_name(gene).with_context(|| format!("NullPointerException: gene({gene}) not found in Ensembl data"))?.chromosome.clone())
    }

    fn submit_chromosome(&self, chromosomes: &[String], gene: Option<(&str, Option<&str>)>, specs: &mut Vec<PlotSpec>) -> Result<()> {
        if !all_human(chromosomes) {
            warn!("invalid chromosomes: [{}]", chromosomes.join(", "));
            return Ok(());
        }
        let chr_str = if chromosomes.len() >= ALL_CHROMOSOMES.len() {
            "All".to_string()
        } else {
            chromosomes.iter().map(|c| strip_chr(c).to_string()).collect::<Vec<_>>().join("-")
        };
        let mut parts = vec![self.config.sample.clone(), format!("chr{chr_str}")];
        if let Some((g, driver)) = gene {
            parts.push(format!("gene-{g}"));
            if let Some(dt) = driver {
                parts.push(format!("driver-{dt}"));
            }
        }
        if self.config.debug {
            parts.push("debug".into());
        }
        let file_id = parts.join(".");
        let include_line = self.config.include_line_elements;
        let cluster_ids: HashSet<i32> = self
            .data
            .sv_data
            .iter()
            .filter(|x| (chromosomes.contains(&x.chr_start) || chromosomes.contains(&x.chr_end)) && (!x.is_line_element() || include_line))
            .map(|x| x.cluster_id)
            .collect();
        let links: Vec<VisSvData> = self.data.sv_data.iter().filter(|x| cluster_ids.contains(&x.cluster_id)).cloned().collect();
        if self.opts.exceeds_max_plot_sv_count(links.len()) {
            warn!("plot({}) - chromosomes({}) svCount({}) exceeds limit({})", file_id, chr_str, links.len(), self.opts.max_plot_sv_count);
            return Ok(());
        }
        let mut segments: Vec<VisSegment> = self.data.segments.iter().filter(|x| cluster_ids.contains(&x.cluster_id)).cloned().collect();
        for c in chromosomes {
            segments.push(entire_chromosome(c, &self.coords));
        }
        let mut of_interest: HashSet<String> = chromosomes.iter().cloned().collect();
        for l in &links {
            of_interest.insert(l.chr_start.clone());
            of_interest.insert(l.chr_end.clone());
        }
        for s in &segments {
            of_interest.insert(s.chromosome.clone());
        }
        let exons = match gene {
            None => {
                let mut e: Vec<VisGeneExon> = self.data.exons.iter().filter(|x| of_interest.contains(&x.chromosome)).cloned().collect();
                if !self.config.genes.is_empty() {
                    let extra = self.exons_from_cache(&self.config.genes_java_order(), &e, &file_id)?;
                    e.extend(extra);
                }
                e
            }
            Some((g, _)) => self.exons_from_cache(&[g.to_string()], &[], &file_id)?,
        };
        let protein_domains: Vec<VisProteinDomain> = self.data.protein_domains.iter().filter(|x| of_interest.contains(&x.chromosome)).cloned().collect();
        specs.push(PlotSpec {
            file_id,
            links,
            segments,
            exons,
            protein_domains,
            fusions: Vec::new(),
            show_simple_sv_segments: false,
            cluster_colours: true,
        });
        Ok(())
    }

    /// getExonDataFromCache: the canonical transcript's exons of each gene, once per given cluster (or for
    /// cluster 0)
    fn exons_from_cache(&self, genes: &[String], current: &[VisGeneExon], file_id: &str) -> Result<Vec<VisGeneExon>> {
        let mut out = Vec::new();
        let cluster_ids = &self.config.cluster_ids;
        let all: Vec<i32> = if cluster_ids.is_empty() { vec![0] } else { cluster_ids.clone() };
        for g in genes {
            if current.iter().any(|x| x.gene == *g && cluster_ids.contains(&x.cluster_id)) {
                continue;
            }
            debug!("plot({}) - loading exon data for gene({})", file_id, g);
            let e = self.ensembl.as_ref().context("NullPointerException: no Ensembl data (-ensembl_data_dir)")?;
            let gd = e.gene_data_by_name(g);
            let td = gd.and_then(|gd| e.canonical_transcript_data(&gd.gene_id));
            let (Some(gd), Some(td)) = (gd, td) else {
                warn!("data not found for specified gene({})", g);
                continue;
            };
            for &c in &all {
                for x in &td.exons {
                    out.push(VisGeneExon {
                        cluster_id: c,
                        gene: gd.gene_name.clone(),
                        transcript: td.trans_name.clone(),
                        chromosome: gd.chromosome.clone(),
                        annotation_type: AnnotationType::DRIVER,
                        exon_rank: x.rank,
                        exon_start: x.start,
                        exon_end: x.end,
                    });
                }
            }
        }
        Ok(out)
    }

    fn submit_cluster(&self, cluster_id: i32, chain_ids: &[i32], skip_singles: bool, reportable: bool, specs: &mut Vec<PlotSpec>) {
        let svs: Vec<VisSvData> =
            self.data.sv_data.iter().filter(|x| x.cluster_id == cluster_id && (chain_ids.is_empty() || chain_ids.contains(&x.chain_id))).cloned().collect();
        if svs.is_empty() {
            warn!("cluster({}) not present in svAnnotations file", cluster_id);
            return;
        }
        if self.opts.exceeds_max_plot_sv_count(svs.len()) {
            warn!("clusterId({}) svCount({}) exceeds limit({})", cluster_id, svs.len(), self.opts.max_plot_sv_count);
            return;
        }
        let segments: Vec<VisSegment> =
            self.data.segments.iter().filter(|x| x.cluster_id == cluster_id && (chain_ids.is_empty() || chain_ids.contains(&x.chain_id))).cloned().collect();
        let exons: Vec<VisGeneExon> = self.data.exons.iter().filter(|x| x.cluster_id == cluster_id).cloned().collect();
        if svs.len() == 1 && skip_singles && exons.is_empty() {
            debug!("skipping simple cluster({})", cluster_id);
            return;
        }
        let link_chains: HashSet<i32> = svs.iter().map(|x| x.chain_id).collect();
        let mut extra: Vec<i32> = segments.iter().map(|x| x.chain_id).filter(|c| !link_chains.contains(c)).collect();
        extra.sort();
        extra.dedup();
        if !extra.is_empty() {
            warn!("cluster({}) contains chain ids({:?}) not found in the links", cluster_id, extra);
            return;
        }
        let mut file_id = self.config.sample.clone();
        if !self.config.genes.is_empty() && self.config.restrict_cluster_by_gene {
            file_id += if self.config.genes.len() > 1 { ".genes-" } else { ".gene-" };
            file_id += &self.config.genes_java_order().join("-");
        } else {
            file_id += &hmf_common::linx_files::vis_plot_cluster_prefix(&cluster_id.to_string());
            if !reportable {
                file_id += &format!(".resolved_type-{}", svs[0].resolved_type);
            }
        }
        if !reportable {
            file_id += &format!(".sv_count-{}", svs.len());
        }
        if self.config.debug {
            file_id += ".debug";
        }
        let protein_domains: Vec<VisProteinDomain> = self.data.protein_domains.iter().filter(|x| x.cluster_id == cluster_id).cloned().collect();
        let fusions: Vec<VisFusion> = self.data.fusions.iter().filter(|x| x.cluster_id == cluster_id).cloned().collect();
        specs.push(PlotSpec { file_id, links: svs, segments, exons, protein_domains, fusions, show_simple_sv_segments: true, cluster_colours: false });
    }

    fn submit_gene_cnvs(&self, specs: &mut Vec<PlotSpec>) -> Result<()> {
        for (gene, driver) in &self.data.purple_gene_cnvs {
            let chr = self.gene_chromosome(gene)?;
            self.submit_chromosome(&[chr], Some((gene, Some(driver))), specs)?;
        }
        Ok(())
    }

    fn plan(&self) -> Result<Vec<PlotSpec>> {
        let mut specs = Vec::new();
        let c = self.config;
        if c.plot_reportable {
            for id in self.data.reportable_clusters(c) {
                self.submit_cluster(id, &[], true, true, &mut specs);
            }
            self.submit_gene_cnvs(&mut specs)?;
        } else if !c.cluster_ids.is_empty() || !c.chromosomes.is_empty() || !c.genes.is_empty() {
            for &id in &c.cluster_ids {
                self.submit_cluster(id, &c.chain_ids, false, false, &mut specs);
            }
            if !c.chromosomes.is_empty() {
                self.submit_chromosome(&c.chromosomes, None, &mut specs)?;
            }
            for g in &c.genes {
                let chr = self.gene_chromosome(g)?;
                let driver = self.data.purple_gene_cnvs.iter().find(|x| x.0 == *g).map(|x| x.1.as_str());
                if let Some(dt) = driver {
                    debug!("found CNV entry({}) in Purple driver catalog for manually selected gene({})", dt, g);
                }
                self.submit_chromosome(&[chr], Some((g, driver)), &mut specs)?;
            }
        } else {
            let mut ids: Vec<i32> = self.data.sv_data.iter().map(|x| x.cluster_id).collect();
            ids.sort();
            ids.dedup();
            for id in ids {
                self.submit_cluster(id, &[], true, false, &mut specs);
            }
            let mut chrs: Vec<String> = Vec::new();
            for x in &self.data.sv_data {
                for chr in [&x.chr_start, &x.chr_end] {
                    if human_chromosome_contains(chr) && !chrs.contains(chr) {
                        chrs.push(chr.clone());
                    }
                }
            }
            for chr in chrs {
                self.submit_chromosome(&[chr], None, &mut specs)?;
            }
            self.submit_gene_cnvs(&mut specs)?;
        }
        Ok(specs)
    }

    /// submitFiltered and what its callables do: the data files, the configuration and the image
    fn make_plot(&self, spec: &PlotSpec) -> Result<()> {
        let t0 = std::time::Instant::now();
        let mut positions: Vec<GenomePosition> = Vec::new();
        positions.extend(link_positions(&spec.links));
        positions.extend(all_positions(&spec.segments));
        positions.extend(all_positions(&spec.exons));
        let copy_numbers = limit_copy_numbers(&self.data.copy_numbers, &span_positions(&positions));
        positions.extend(all_positions(&copy_numbers));
        let mut positions = distinct_positions(positions);
        sort_positions(&mut positions);
        let regions = span_positions(&positions);
        let mut amber = Vec::new();
        let mut cobalt = Vec::new();
        let mut purple = Vec::new();
        for r in &regions {
            if self.config.amber_dir.is_some() {
                if let Some(c) = self.data.amber.iter().find(|c| c.chromosome == r.chromosome) {
                    amber.extend(c.in_range(r.start, r.end).into_iter().map(|(p, v)| (r.chromosome.clone(), p, v)));
                }
            }
            if self.config.cobalt_dir.is_some() {
                if let Some(c) = self.data.cobalt.iter().find(|c| c.chromosome == r.chromosome) {
                    cobalt.extend(c.in_range(r.start, r.end).into_iter().map(|(p, v)| (r.chromosome.clone(), p, v)));
                }
            }
            if self.config.purple_dir.is_some() {
                for p in &self.data.purple_segments {
                    if p.chromosome == r.chromosome && p.pos_end > r.start && p.pos_start < r.end {
                        let mut n = p.clone();
                        n.pos_start = i32::max(p.pos_start, r.start);
                        n.pos_end = i32::min(p.pos_end, r.end);
                        purple.push(n);
                    }
                }
            }
        }
        let segments = extend_terminals(0, &spec.segments, &spec.links, &positions, spec.show_simple_sv_segments, &self.coords)?;
        let links = add_frame(&segments, &spec.links);
        let colours = if spec.cluster_colours { ColorPicker::cluster_colours(&links) } else { ColorPicker::chain_colours(&links) };
        let data = CircosData::new(
            self.opts,
            &self.sites,
            CircosInputs {
                segments,
                links,
                copy_numbers,
                exons: &spec.exons,
                fusions: &spec.fusions,
                amber,
                cobalt,
                purple,
                show_simple_sv_segments: spec.show_simple_sv_segments,
                show_fragile_sites: self.config.include_fragile_sites,
                show_line_elements: self.config.include_line_elements,
            },
        )?;
        let cw = ConfigWriter::new(&data, self.opts, self.config.amber_dir.is_some(), self.config.cobalt_dir.is_some());
        let files = write_circos_data(&data, &colours, &cw, self.opts);
        let fusion = FusionData::new(&spec.fusions, &spec.exons, &spec.protein_domains)?;
        let prefix = format!("{}{}", self.config.output_conf_path, spec.file_id);
        if self.config.write_data {
            for (suffix, lines) in &files.files {
                write_lines(&format!("{prefix}.{suffix}"), lines)?;
            }
            if !fusion.exons.is_empty() {
                write_lines(&format!("{prefix}.fusions.tsv"), &fusion.exon_lines())?;
                write_lines(&format!("{prefix}.protein_domains.tsv"), &fusion.domain_lines())?;
            }
        }
        let min_frame = if self.opts.step { 0 } else { data.max_frame };
        let t_data = t0.elapsed().as_secs_f64();
        for frame in min_frame..=data.max_frame {
            let conf = cw.config(&data, self.opts, &spec.file_id, frame);
            let conf_name = format!("{}.circos.{:03}.conf", spec.file_id, frame);
            if self.config.write_data {
                std::fs::write(format!("{}{}", self.config.output_conf_path, conf_name), conf.as_bytes())?;
            }
            let plot_chromosome = !self.config.debug;
            if plot_chromosome && self.config.write_data {
                std::fs::write(format!("{prefix}.cytoBand.txt"), plot::cytoband_resource(self.config.v37))?;
            }
            if self.config.write_plots {
                let image = format!("{}{}.{:03}.png", self.config.output_plot_path, spec.file_id, frame);
                let fusion_panel = !fusion.exons.is_empty() && !self.config.debug;
                let ranges = files.get("chromosome.circos").cloned().unwrap_or_default();
                plot::draw(&plot::PlotInputs {
                    file_id: &spec.file_id,
                    conf: &conf,
                    files: &files,
                    fusion: if fusion_panel { Some(&fusion) } else { None },
                    chromosome_ranges: if plot_chromosome { Some(&ranges) } else { None },
                    label_size: 1.2 * data.label_size,
                    opts: self.opts,
                    v37: self.config.v37,
                    out: &image,
                    perf: self.config.perf_debug,
                    final_frame: frame == data.max_frame,
                })?;
            }
        }
        if self.config.perf_debug {
            info!("perf: plot({}) data {:.3} s, total {:.3} s", spec.file_id, t_data, t0.elapsed().as_secs_f64());
        }
        Ok(())
    }
}

/// Files.write(path, lines): every line followed by a newline
fn write_lines(path: &str, lines: &[String]) -> Result<()> {
    let mut s = String::with_capacity(lines.iter().map(|l| l.len() + 1).sum());
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    std::fs::write(path, s).with_context(|| format!("writing {path}"))
}

/// Runs the visualiser with its command line (the main class name already removed).
pub fn run(argv: &[String]) -> Result<i32> {
    // the pipeline passes `-circos $(which circos)`, which is empty where Circos is not installed
    let mut args: Vec<String> = Vec::new();
    let mut i = 0;
    while i < argv.len() {
        if argv[i] == "-circos" && (i + 1 >= argv.len() || argv[i + 1].starts_with('-')) {
            args.push("-circos_flag".into());
            i += 1;
            continue;
        }
        args.push(argv[i].clone());
        i += 1;
    }
    let a = Args::parse(&args, &FLAGS, &ITEMS)?;
    let level = a.value_or("log_level", "").to_ascii_uppercase();
    if a.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    let t0 = std::time::Instant::now();
    let opts = VisOptions::from_args(&a)?;
    let mut config = VisConfig::new(&a)?;
    let coords = RefCoords { v37: config.v37 };
    let sites = KnownSites::load(config.v37);
    let sample = match load_sample_data(&mut config)? {
        Ok(s) => s,
        Err(code) => return Ok(code),
    };
    if !opts.is_valid() {
        return Ok(1);
    }
    // the Ensembl data: gene names for the driver gene plots and the canonical transcripts of those genes
    let mut ensembl = None;
    if let Some(dir) = &config.ensembl_data_dir {
        debug!("loading Ensembl data");
        let version = if config.v37 { RefGenomeVersion::V37 } else { RefGenomeVersion::V38 };
        let mut e = EnsemblDataCache::new(dir, version);
        e.set_required_data(true, false, false, true);
        e.load(true)?;
        let mut names: Vec<String> = config.genes.clone();
        names.extend(sample.purple_gene_cnvs.iter().map(|g| g.0.clone()));
        let ids: Vec<String> = names.iter().filter_map(|n| e.gene_data_by_name(n).map(|g| g.gene_id.clone())).collect();
        if !ids.is_empty() {
            e.load_transcript_data(&ids, &[])?;
        }
        ensembl = Some(e);
    }
    let vis = Visualiser { config: &config, opts: &opts, data: &sample, ensembl, coords, sites };
    let specs = vis.plan()?;
    if config.perf_debug {
        info!("perf: loading and planning {:.3} s, {} plots", t0.elapsed().as_secs_f64(), specs.len());
    }
    let next = Mutex::new(0usize);
    let failed: Mutex<Option<anyhow::Error>> = Mutex::new(None);
    std::thread::scope(|scope| {
        for _ in 0..config.threads.min(specs.len().max(1)) {
            scope.spawn(|| loop {
                let k = {
                    let mut n = next.lock().unwrap();
                    let k = *n;
                    *n += 1;
                    k
                };
                if k >= specs.len() || failed.lock().unwrap().is_some() {
                    break;
                }
                if let Err(e) = vis.make_plot(&specs[k]) {
                    let mut f = failed.lock().unwrap();
                    if f.is_none() {
                        *f = Some(e.context(format!("plot({})", specs[k].file_id)));
                    }
                }
            });
        }
    });
    if let Some(e) = failed.into_inner().unwrap() {
        return Err(e);
    }
    info!("Linx Visualiser complete");
    Ok(0)
}
