//! LINX in Rust: a port of hmftools LINX 2.4 (structural variant clustering, chaining and annotation;
//! see ../../hmftools/linx).
//!
//! Usage: `linx-rs <the Java options>` for one sample, somatic or (`-germline`) germline. JVM options
//! (-Xmx..., -XX:..., -D...=...) are dropped, as the bioconda wrapper passes them to java. `-java_compat`
//! reproduces the Java faults the port corrects (REPORT.md).

#![allow(dead_code)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]
#![allow(clippy::too_many_arguments)]
// the Java control flow is kept as written (if / else if chains with the same body, nested ifs), and the
// handle types hash by identity, which clippy takes for mutable keys
#![allow(clippy::if_same_then_else, clippy::collapsible_if, clippy::collapsible_else_if, clippy::mutable_key_type)]

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod analysis;
mod annotators;
mod chaining;
mod cn;
mod compat;
mod drivers;
mod fusion;
mod gene;
mod idmap;
mod loader;
mod output;
mod types;
mod utils;
mod vis;
#[cfg(test)]
mod tests;

use std::rc::Rc;

use anyhow::{bail, Context, Result};
use hmf_common::args::Args;
use hmf_common::driver_panel::{find_invalid_driver_genes, load_driver_genes, DriverGene};
use hmf_common::ensembl::EnsemblDataCache;
use hmf_common::java_hash::JavaHashMap;
use hmf_common::known_fusion::{KnownFusionCache, DEFAULT_PRE_GENE_PROMOTOR_DISTANCE};
use hmf_common::purple_files::{driver_catalog_somatic_filename, read_driver_catalog, write_driver_catalog, write_lines};
use hmf_common::regions::RefGenomeVersion;
use hmf_common::sv::StructuralVariantData;
use hmf_common::{error, info};

use crate::analysis::clustering::link_sgl_mapped_inferreds;
use crate::analysis::ClusterAnalyser;
use crate::annotators::{FragileSiteAnnotator, LineElementAnnotator, PseudoGeneFinder};
use crate::cn::{set_sv_copy_number_data, CnDataLoader, CnSegmentBuilder};
use crate::drivers::{linx_driver_lines, DriverGeneAnnotator};
use crate::fusion::analyser::FusionDisruptionAnalyser;
use crate::fusion::disruption::DisruptionFinder;
use crate::fusion::finder::{specific_pre_gene_distances, FusionFinder};
use crate::fusion::FusionConfig;
use crate::gene::{set_sv_gene_data, GeneCache};
use crate::output::{sample_file, write_clusters, write_links, write_svs};
use crate::types::*;
use crate::utils::get_chromosomal_arm;
use crate::vis::{VisDataWriter, VisSampleData};

pub const VERSION: &str = "2.4";
const APP_NAME: &str = "Linx";

const FLAGS: [&str; 22] = [
    "germline",
    "write_all",
    "write_sv_data",
    "write_links",
    "write_cluster_history",
    "write_single_sv_clusters",
    "write_vis_data",
    "write_cohort",
    "no_vis_files",
    "fail_on_missing",
    "perf_debug",
    "log_verbose",
    "skip_unphased_fusions",
    "write_neo_epitopes",
    "log_reportable_fusions",
    "log_all_potential_fusions",
    "log_invalid_fusions",
    "write_all_vis_fusions",
    "log_debug",
    "version",
    "help",
    "h",
];

const ITEMS: [&str; 25] = [
    "sample",
    "sample_id_file",
    "purple_dir",
    "sample_data_dir",
    "ref_genome_version",
    "sv_vcf",
    "driver_gene_panel",
    "line_element_file",
    "fragile_site_file",
    "proximity_distance",
    "chaining_sv_limit",
    "annotations",
    "gene_id_file",
    "write_type",
    "log_chain_size",
    "output_dir",
    "threads",
    "known_fusion_file",
    "fusion_gene_distance",
    "rna_fusions_file",
    "fusion_debug_sv_genes",
    "ensembl_data_dir",
    "log_level",
    "rna_file_source",
    "sample_rna_file",
];

/// -perf_debug: the time of each phase, logged as it ends
struct PhaseTimer {
    on: bool,
    last: std::time::Instant,
}
impl PhaseTimer {
    fn new(on: bool) -> Self {
        PhaseTimer { on, last: std::time::Instant::now() }
    }
    fn lap(&mut self, phase: &str) {
        if self.on {
            info!("perf: {} {:.3} s", phase, self.last.elapsed().as_secs_f64());
        }
        self.last = std::time::Instant::now();
    }
}

fn is_jvm_option(a: &str) -> bool {
    a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-Xss") || a.starts_with("-XX:") || (a.starts_with("-D") && a.contains('='))
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).filter(|a| !is_jvm_option(a) && !a.is_empty()).collect();
    std::process::exit(match run(&argv) {
        Ok(code) => code,
        Err(e) => {
            error!("{e:#}");
            1
        }
    });
}

/// LinxConfig with LinxOutput and FusionConfig
struct LinxConfig {
    sample_id: String,
    purple_data_path: String,
    output_dir: String,
    sv_vcf_file: String,
    is_germline: bool,
    run_fusions: bool,
    write_vis: bool,
    write_neo_epitopes: bool,
    proximity_distance: i32,
    chaining_sv_limit: i32,
    fragile_site_file: Option<String>,
    line_element_file: Option<String>,
    log_dm_candidates: bool,
    show_dm_candidates: bool,
    driver_genes: &'static JavaHashMap<String, DriverGene>,
    restricted_gene_ids: Vec<String>,
    ensembl_dir: Option<String>,
    known_fusion_file: Option<String>,
    pre_gene_distance: i32,
    fusion: FusionConfig,
    ref_genome_version: RefGenomeVersion,
    perf_debug: bool,
}

fn dir_with_sep(d: &str) -> String {
    if d.ends_with('/') { d.to_string() } else { format!("{d}/") }
}

impl LinxConfig {
    fn new(a: &Args) -> Result<LinxConfig> {
        let Some(sample) = a.value("sample") else { bail!("-sample is required (sample id files are not supported)") };
        if sample.contains(',') {
            bail!("multiple samples are not supported");
        }
        let is_germline = a.has_flag("germline");
        let (purple_data_path, output_dir) = if let Some(sd) = a.value("sample_data_dir") {
            let p = dir_with_sep(sd);
            (p.clone(), a.value("output_dir").map(|d| d.to_string()).unwrap_or(p))
        } else {
            (a.value_or("purple_dir", "").to_string(), a.value("output_dir").context("-output_dir is required")?.to_string())
        };
        let output_dir = output_dir.trim_end_matches('/').to_string();
        let mut sv_vcf_file = a.value_or("sv_vcf", "").to_string();
        if sv_vcf_file.is_empty() {
            let p = dir_with_sep(&purple_data_path);
            sv_vcf_file =
                if is_germline { format!("{p}{sample}.purple.sv.germline.vcf.gz") } else { format!("{p}{sample}.purple.sv.vcf.gz") };
        }
        let run_fusions = a.has_value("known_fusion_file");
        let write_type = a.value_or("write_type", "STANDARD");
        let single_sample_somatic = !is_germline;
        let mut write_vis = write_type == "ALL" || write_type.split(';').any(|t| t == "VIS_DATA");
        if write_type != "write_all" && single_sample_somatic && !a.has_flag("no_vis_files") {
            write_vis = true;
        }
        let rgv = match a.value_or("ref_genome_version", "V37") {
            "38" | "V38" | "HG38" | "GRCh38" => RefGenomeVersion::V38,
            "37" | "V37" | "HG19" | "GRCh37" => RefGenomeVersion::V37,
            v => bail!("invalid ref genome version: {v}"),
        };
        let annotations = a.value_or("annotations", "");
        let exts: Vec<&str> = if annotations.is_empty() { Vec::new() } else { annotations.split(';').collect() };
        for e in &exts {
            if !matches!(*e, "DOUBLE_MINUTES" | "CANDIDATE_VIS_DOUBLE_MINUTES" | "LINE_CHAINS" | "UNDER_CLUSTERING") {
                bail!("No enum constant com.hartwig.hmftools.linx.analysis.AnnotationExtension.{e}");
            }
        }
        let genes = load_driver_genes(a.value("driver_gene_panel"))?;
        let mut driver_genes: JavaHashMap<String, DriverGene> = JavaHashMap::new();
        for g in genes {
            driver_genes.insert(g.gene.clone(), g);
        }
        let mut restricted_gene_ids = Vec::new();
        if let Some(f) = a.value("gene_id_file") {
            for (i, line) in std::fs::read_to_string(f).with_context(|| format!("reading {f}"))?.lines().enumerate() {
                if i == 0 && line.contains("GeneId") {
                    continue;
                }
                if !line.is_empty() {
                    restricted_gene_ids.push(line.to_string());
                }
            }
            info!("file({}) loaded {} restricted genes", f, restricted_gene_ids.len());
        }
        let mut debug_sv_genes = Vec::new();
        if let Some(v) = a.value("fusion_debug_sv_genes") {
            for p in v.split(';') {
                if let Some((g, id)) = p.split_once(':') {
                    debug_sv_genes.push((g.to_string(), id.parse()?));
                }
            }
        }
        let fusion = FusionConfig {
            log_reportable_only: a.has_flag("log_reportable_fusions"),
            log_all_potentials: a.has_flag("log_all_potential_fusions"),
            write_all_vis_fusions: a.has_flag("write_all_vis_fusions"),
            allow_exon_skipping: true,
            require_phase_match: a.has_flag("skip_unphased_fusions"),
            require_upstream_biotypes: true,
            debug_sv_genes,
        };
        if a.has_value("rna_fusions_file") {
            bail!("-rna_fusions_file (RNA fusion matching) is not ported");
        }
        Ok(LinxConfig {
            sample_id: sample.to_string(),
            purple_data_path,
            output_dir,
            sv_vcf_file,
            is_germline,
            run_fusions,
            write_vis,
            write_neo_epitopes: a.has_flag("write_neo_epitopes"),
            proximity_distance: a.integer("proximity_distance", DEFAULT_PROXIMITY_DISTANCE)?,
            chaining_sv_limit: a.integer("chaining_sv_limit", 0)?,
            fragile_site_file: a.value("fragile_site_file").map(|s| s.to_string()),
            line_element_file: a.value("line_element_file").map(|s| s.to_string()),
            log_dm_candidates: exts.contains(&"DOUBLE_MINUTES"),
            show_dm_candidates: exts.contains(&"CANDIDATE_VIS_DOUBLE_MINUTES"),
            driver_genes: Box::leak(Box::new(driver_genes)),
            restricted_gene_ids,
            ensembl_dir: a.value("ensembl_data_dir").map(|s| s.to_string()),
            known_fusion_file: a.value("known_fusion_file").map(|s| s.to_string()),
            pre_gene_distance: a.integer("fusion_gene_distance", DEFAULT_PRE_GENE_PROMOTOR_DISTANCE)?,
            fusion,
            ref_genome_version: rgv,
            perf_debug: a.has_flag("perf_debug"),
        })
    }
    fn run_drivers(&self) -> bool {
        !self.driver_genes.is_empty()
    }
}

fn run(argv: &[String]) -> Result<i32> {
    let argv = compat::take_options(argv)?;
    let a = Args::parse(&argv, &FLAGS, &ITEMS)?;
    if a.has_flag("version") {
        println!("{APP_NAME} version {VERSION}");
        return Ok(0);
    }
    let level = a.value_or("log_level", "").to_ascii_uppercase();
    if a.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("{APP_NAME} version {VERSION}");
    let config = LinxConfig::new(&a)?;
    utils::set_ref_genome_v37(config.ref_genome_version == RefGenomeVersion::V37);
    std::fs::create_dir_all(&config.output_dir).with_context(|| format!("failed to create output directory({})", config.output_dir))?;
    if !a.has_value("sv_vcf") && !a.has_value("sample_data_dir") && !a.has_value("purple_dir") {
        error!("missing SV VCF file or sample data directory");
        return Ok(1);
    }
    if !config.is_germline && config.purple_data_path.is_empty() {
        error!("missing purple directory");
        return Ok(1);
    }
    info!("running {} analysis for {}", if config.is_germline { "germline SV" } else { "SV" }, config.sample_id);

    let mut timer = PhaseTimer::new(config.perf_debug);
    let mut kfc = KnownFusionCache::new();
    if let Some(f) = &config.known_fusion_file {
        kfc.load_from_file(f);
    }
    if config.run_fusions && !kfc.has_valid_data() {
        info!("invalid fusion config, exiting");
        return Ok(1);
    }
    let kfc: &'static KnownFusionCache = Box::leak(Box::new(kfc));
    timer.lap("known fusion data");

    let gene_cache: Option<&'static GeneCache> = match &config.ensembl_dir {
        Some(dir) => {
            let mut e = EnsemblDataCache::new(dir, config.ref_genome_version);
            e.set_required_data(true, config.run_fusions, config.run_fusions, config.is_germline);
            if !config.restricted_gene_ids.is_empty() {
                e.set_restricted_gene_id_list(&config.restricted_gene_ids);
            }
            if let Err(err) = e.load(false) {
                error!("Ensembl data cache load failed, exiting: {err:#}");
                return Ok(1);
            }
            if config.run_drivers() {
                let genes: Vec<DriverGene> = config.driver_genes.values().cloned().collect();
                let invalid = find_invalid_driver_genes(&genes, &e);
                if !invalid.is_empty() {
                    let names: Vec<&str> = invalid.iter().map(|g| g.gene.as_str()).collect();
                    error!("Invalid non-Ensembl gene names: {}", names.join(";"));
                    return Ok(1);
                }
                e.create_gene_name_id_map();
            }
            let e: &'static EnsemblDataCache = Box::leak(Box::new(e));
            Some(Box::leak(Box::new(GeneCache::new(e))))
        }
        // Java builds the cache from an empty path, whose load fails
        None => {
            error!("Ensembl data cache load failed, exiting");
            return Ok(1);
        }
    };

    timer.lap("Ensembl data");
    let mut fragile = FragileSiteAnnotator::default();
    fragile.load_fragile_sites_file(config.fragile_site_file.as_deref());
    let mut line = LineElementAnnotator::new(config.proximity_distance);
    line.load_line_elements_file(config.line_element_file.as_deref());
    let pseudo = PseudoGeneFinder::new(gene_cache);

    let svs = process_sample(&config, gene_cache, kfc, &fragile, &line, &pseudo, &mut timer)?;
    output::write_version(&config.output_dir)?;
    let _ = svs;
    info!("Linx complete for {}", config.sample_id);
    Ok(0)
}

/// SampleAnalyser.processSample
fn process_sample(
    config: &LinxConfig, gene_cache: Option<&'static GeneCache>, kfc: &'static KnownFusionCache, fragile: &FragileSiteAnnotator,
    line: &LineElementAnnotator, pseudo: &PseudoGeneFinder, timer: &mut PhaseTimer,
) -> Result<usize> {
    let sample = config.sample_id.as_str();
    let mut vis = VisSampleData::default();
    vis.set_sample_id(sample);
    let records: Vec<StructuralVariantData> = if config.is_germline {
        loader::load_germline_variants(&config.sv_vcf_file, sample)?
    } else {
        loader::load_somatic_variants(&config.sv_vcf_file)?
    };
    let records: Vec<Rc<StructuralVariantData>> = records.into_iter().map(Rc::new).collect();
    let mut all: Vec<Sv> = records.iter().filter(|r| loader::keep_sv(&r.filter, config.is_germline)).map(|r| Sv::from_data(r.clone())).collect();
    timer.lap("SV VCF");
    let mut cn_loader = CnDataLoader::new(&config.purple_data_path);
    if !config.is_germline {
        cn_loader.load_sample_data(sample, &records);
    }
    timer.lap("copy number data");
    let vis_writer = VisDataWriter::new(&config.output_dir, gene_cache, config.write_vis, config.is_germline);
    if all.is_empty() {
        info!("sample({}) has no passing SVs", sample);
        write_sample_with_no_svs(config, &vis_writer, &vis, &cn_loader)?;
        return Ok(0);
    }
    let sv_count = all.len();

    // setSampleSVs
    link_sgl_mapped_inferreds(&mut all);
    if !config.is_germline {
        set_sv_copy_number_data(&all, Some(&cn_loader.sv_jcn_calc_map), &cn_loader.sv_id_cn_data_map, &cn_loader.chr_cn_data_map);
    }

    let fusion_finder_cache = gene_cache;
    let mut fusion_analyser = gene_cache.map(|gc| {
        FusionDisruptionAnalyser::new(
            FusionFinder::new(config.fusion.clone(), gc, kfc),
            DisruptionFinder::new(config.is_germline, gc, config.driver_genes),
            &config.output_dir,
            config.run_fusions,
            config.is_germline,
            config.write_vis,
            config.write_neo_epitopes,
            config.driver_genes,
        )
    });
    if config.run_fusions {
        if let Some(fa) = &fusion_analyser {
            fa.cache_special_fusion_genes();
        }
    }
    if let Some(gc) = fusion_finder_cache {
        let specific = specific_pre_gene_distances(kfc);
        let pre_gene = if config.run_fusions { config.pre_gene_distance } else { 0 };
        set_sv_gene_data(&all, gc, pre_gene, &specific);
    }
    timer.lap("breakend genes");

    // analyse
    for var in &all {
        var.set_fragile_sites(fragile.is_fragile_site(var, true), fragile.is_fragile_site(var, false));
        line.set_known_line_elements(var);
        let start_arm = get_chromosomal_arm(var.chromosome(true), var.position(true));
        let end_arm = if !var.is_sgl_breakend() { get_chromosomal_arm(var.chromosome(false), var.position(false)) } else { Arm::P };
        var.set_chromosomal_arms(start_arm, end_arm);
    }
    let mut analyser = ClusterAnalyser::new(
        config.proximity_distance,
        config.is_germline,
        config.chaining_sv_limit,
        config.log_dm_candidates,
        config.show_dm_candidates,
    );
    analyser.set_sample_data(sample, &all, &cn_loader);
    analyser.pre_clustering_preparation();
    if config.is_germline {
        let mut b = CnSegmentBuilder::new();
        b.set_allele_ploidies(1.0, 1.0);
        b.create_germline_copy_number_data(&mut cn_loader, &analyser.state.chr_breakend_map);
        CnSegmentBuilder::set_sample_purity(&mut cn_loader, 1.0, 2.0, "MALE");
        cn_loader.create_chr_copy_number_map();
        set_sv_copy_number_data(&all, Some(&cn_loader.sv_jcn_calc_map), &cn_loader.sv_id_cn_data_map, &cn_loader.chr_cn_data_map);
        for var in &all {
            var.set_jcn_recalc_data(1.0, 1.0);
        }
    }
    let valid = analyser.cluster_and_analyse(line, Some(pseudo), &cn_loader.chr_ends_cn_map);
    if !valid {
        info!("exiting after sample({}), in invalid state", sample);
        return Ok(sv_count);
    }
    timer.lap("clustering and chaining");
    if let Some(fa) = &fusion_analyser {
        fa.annotate_transcripts(&all, true);
    }

    // annotate
    analyser.annotate_clusters();
    pseudo.check_pseudo_gene_annotations(&analyser.clusters, &mut vis);

    let mut driver_annotator = match gene_cache {
        Some(gc) if config.run_drivers() => Some(DriverGeneAnnotator::new(gc, config.driver_genes, &config.purple_data_path)),
        _ => None,
    };
    if let Some(d) = driver_annotator.as_mut() {
        d.annotate_svs(sample, &analyser.state.chr_breakend_map, &cn_loader, &mut vis);
    }
    timer.lap("transcripts, cluster annotation and drivers");
    if config.run_fusions || config.is_germline {
        if let Some(fa) = fusion_analyser.as_mut() {
            fa.run(sample, &all, &analyser.clusters, &mut vis)?;
        }
    }

    timer.lap("fusions and disruptions");
    // writeOutput
    let all_clusters = analyser.all_clusters();
    let dir = &config.output_dir;
    vis_writer.write_output(&vis, &analyser.clusters, &all, &cn_loader.chr_cn_data_map)?;
    let germline = config.is_germline;
    write_svs(&sample_file(dir, sample, "svs", germline), &all)?;
    write_clusters(&sample_file(dir, sample, "clusters", germline), &all_clusters)?;
    write_links(&sample_file(dir, sample, "links", germline), &analyser.clusters)?;
    if !germline {
        let mut catalog = Vec::new();
        let mut linx_drivers = Vec::new();
        if let Some(d) = driver_annotator.as_mut() {
            catalog.append(&mut d.driver_catalog);
            linx_drivers.append(&mut d.driver_output);
        }
        if let Some(fa) = &fusion_analyser {
            fa.disruption_finder.add_reportable_disruptions(&mut catalog);
        }
        write_driver_catalog(&format!("{dir}/{sample}.linx.driver.catalog.tsv"), &catalog)?;
        write_lines(&format!("{dir}/{sample}.linx.drivers.tsv"), &linx_driver_lines(&linx_drivers))?;
    }
    if let Some(fa) = fusion_analyser.as_mut() {
        fa.close()?;
    }
    timer.lap("output files");
    info!("sample({}) procesed {} SVs", sample, sv_count);
    Ok(sv_count)
}

/// writeSampleWithNoSVs
fn write_sample_with_no_svs(config: &LinxConfig, vis_writer: &VisDataWriter, vis: &VisSampleData, cn: &CnDataLoader) -> Result<()> {
    let dir = &config.output_dir;
    let sample = &config.sample_id;
    let germline = config.is_germline;
    output::write_empty_sample_files(dir, sample, germline)?;
    let mut purple_drivers = Vec::new();
    if !germline {
        purple_drivers = read_driver_catalog(&driver_catalog_somatic_filename(config.purple_data_path.trim_end_matches('/'), sample))?;
    }
    let catalog_file =
        if germline { format!("{dir}/{sample}.linx.germline.driver.catalog.tsv") } else { format!("{dir}/{sample}.linx.driver.catalog.tsv") };
    write_driver_catalog(&catalog_file, &purple_drivers)?;
    if germline {
        write_lines(&format!("{dir}/{sample}.linx.germline.disruption.tsv"), &[output::germline_disruption_header()])?;
    } else {
        vis_writer.write_copy_number_data(vis, &cn.chr_cn_data_map)?;
        output::write_empty_fusions(&format!("{dir}/{sample}.linx.fusion.tsv"))?;
        write_lines(&format!("{dir}/{sample}.linx.drivers.tsv"), &linx_driver_lines(&[]))?;
        vis_writer.write_empty(sample)?;
    }
    Ok(())
}
