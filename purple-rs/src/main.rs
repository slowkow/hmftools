//! PURPLE in Rust: a port of hmftools PURPLE 4.5 (see ../purple).
//!
//! Usage: `purple-rs <the Java options>`: reads the AMBER, COBALT, ESVEE and PAVE outputs and writes the
//! same files as the Java PurpleApplication (purity and ploidy fit, copy numbers, drivers, the annotated
//! VCFs, QC, Circos input files). `-java_compat` reproduces the Java faults the port corrects (REPORT.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod compat;
mod charts;
mod config;
mod constants;
mod copynumber;
mod drivers;
mod fitting;
mod fittingsnv;
mod gene;
mod germline;
mod germline_sv;
mod inputs;
mod parallel;
mod purity_fitter;
mod reference;
mod region;
mod reseg;
mod segment;
mod somatic;
mod somatic_stream;
mod sv;
mod summary;
mod targeted;
#[cfg(test)]
mod test_support;
mod sv_cache;
mod types;
mod variant;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};
#[allow(unused_imports)]
use hmf_common::{fmt, log};

use anyhow::Result;
use std::collections::HashMap;
use hmf_common::config_builder::{ConfigBuilder, Parsed};

use crate::config::{PurpleConfig, APP_NAME, VERSION};
use crate::types::*;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    // oncoanalyser passes JVM options to the wrapper; they are ignored
    let argv: Vec<String> = argv.into_iter().filter(|a| !(a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-XX:"))).collect();
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    config::register(&mut cb);
    match cb.parse(&argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return,
        Err(_) => std::process::exit(1),
    }
    match run(&cb) {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            error!("failed processing sample: {e:#}");
            std::process::exit(1);
        }
    }
}

fn write_version(dir: &str) -> Result<()> {
    let text = format!("version={}\nbuild.date={}", VERSION, env!("PURPLE_BUILD_DATE"));
    std::fs::write(format!("{dir}purple.version"), text)?;
    Ok(())
}

fn run(cb: &ConfigBuilder) -> Result<i32> {
    let started = std::time::Instant::now();
    let (config, valid) = PurpleConfig::from(cb)?;
    if !valid {
        error!("initialisation error, exiting");
        return Ok(1);
    }
    info!(
        "reference({}) tumor({}) {}",
        config.reference_id.as_deref().unwrap_or("NONE"),
        config.tumor_id,
        if config.target_regions_mode { "running on target-regions only" } else { "" }
    );
    info!("output directory: {}", config.output_dir);
    let (reference, ref_valid) = reference::ReferenceData::load(cb, &config)?;
    if !ref_valid {
        error!("invalid reference data, exiting");
        return Ok(1);
    }
    if !config.fitting.has_valid_values() {
        error!("invalid configured fitting values");
        return Ok(1);
    }
    if !has_valid_sample_names(&config)? {
        return Ok(1);
    }
    let files = &config.sample_files;
    let tumor_id = config.tumor_id.clone();
    let amber_sample = if config.germline_mode() { config.reference_id.clone().unwrap() } else { tumor_id.clone() };
    let load_amber = || {
        inputs::AmberData::load(
            &amber_sample,
            files.amber_directory.as_deref().unwrap_or("null"),
            config.germline_mode(),
            reference.v38,
            reference.target_regions.chromosome_x_region_percentage(),
        )
    };
    let load_cobalt = |amber_gender: Gender| {
        inputs::CobaltData::load(
            config.reference_id.as_deref(),
            Some(&tumor_id),
            files.cobalt_directory.as_deref().unwrap_or("null"),
            amber_gender,
            config.tumor_only_mode(),
            config.germline_mode(),
        )
    };
    let load_variants = || -> Result<(sv_cache::SomaticSvCache, somatic::SomaticVariantCache)> {
        let sv_output = format!("{}{}.purple.sv.vcf.gz", config.output_dir, tumor_id);
        let sv_cache = match files.somatic_sv_vcf_file.as_deref() {
            Some(f) if !f.is_empty() => sv_cache::SomaticSvCache::load(VERSION, f, &sv_output, &config)?,
            _ => sv_cache::SomaticSvCache::empty(),
        };
        let mut somatic_cache = somatic::SomaticVariantCache::empty();
        if config.run_tumor() {
            // HlaCommon.populateGeneData in tumor-only mode: the HLA genes of the Ensembl cache
            let hla = config.tumor_only_mode().then(|| {
                let chromosome = if reference.v38 { "chr6" } else { "6" };
                let names = ["HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"];
                let genes = reference.ensembl.chr_gene_data_map().get(&chromosome.to_string()).map(|l| {
                    l.iter().filter(|g| names.contains(&g.gene_name.as_str())).map(|g| (g.chromosome.clone(), g.gene_start, g.gene_end)).collect()
                });
                somatic::HlaGenes { genes: genes.unwrap_or_default() }
            });
            somatic_cache.load(files.somatic_vcf_file.as_deref().unwrap_or(""), &config, &reference.somatic_hotspots, hla.as_ref())?;
        }
        Ok((sv_cache, somatic_cache))
    };
    // COBALT needs AMBER's gender only in tumor-only mode; otherwise it loads beside the rest
    let (amber, cobalt, mut sv_cache, mut somatic_cache) = if config.threads >= 2 && !config.tumor_only_mode() {
        std::thread::scope(|scope| -> Result<_> {
            let cobalt_thread = scope.spawn(|| load_cobalt(Gender::Female));
            let amber = load_amber()?;
            let (sv, somatic) = load_variants()?;
            let cobalt = cobalt_thread.join().map_err(|_| anyhow::anyhow!("COBALT loading failed"))??;
            Ok((amber, cobalt, sv, somatic))
        })?
    } else {
        let amber = load_amber()?;
        let cobalt = load_cobalt(amber.patient_gender)?;
        let (sv, somatic) = load_variants()?;
        (amber, cobalt, sv, somatic)
    };

    let msi_prediction = load_msi_prediction(&config);

    let amber_gender = amber.patient_gender;
    let cobalt_gender = cobalt.gender();
    let gender = if config.tumor_only_mode() { amber_gender } else { cobalt_gender };
    if cobalt_gender == amber_gender {
        info!("sample gender is {}", cobalt_gender.name().to_lowercase());
    } else {
        warn!("Cobalt gender {} does not match Amber gender {}", cobalt_gender, amber_gender);
    }
    info!("output directory: {}", config.output_dir);
    write_version(&config.output_dir)?;

    let tinc_level = if config.run_tumor() { somatic_cache.tinc_level } else { 0.0 };
    let mut germline_sv_cache = germline_sv::GermlineSvCache::empty();
    if config.run_germline() {
        germline_sv_cache = germline_sv::GermlineSvCache::load(VERSION, files.germline_sv_vcf_file.as_deref().unwrap_or(""), config.reference_id.as_deref(), &tumor_id)?;
        if config.run_tumor() && tinc_level > 0.0 {
            germline_sv_cache.check_tinc_variants(&mut sv_cache, tinc_level);
        }
    }

    info!("applying segmentation");
    let coords = segment::ReferenceCoords { v38: reference.v38 };
    let breakends: Vec<segment::SvPosition> = sv_cache
        .somatic_variants()
        .iter()
        .flat_map(|v| {
            let mut out = vec![segment::SvPosition { chromosome: v.start.chromosome.clone(), position: v.start.cna_position(), sv_type: v.sv_type }];
            if let Some(e) = &v.end {
                out.push(segment::SvPosition { chromosome: e.chromosome.clone(), position: e.cna_position(), sv_type: v.sv_type });
            }
            out
        })
        .collect();
    let factory = region::ObservedRegionFactory { cobalt: &cobalt.cobalt_chromosomes, regions: &reference.special_regions };
    let mut observed = segment::create_observed_regions(&coords, &breakends, &amber, &cobalt, &factory);
    if observed.is_empty() || !segment::validate_observed_regions(&observed) {
        warn!("no valid observed regions created, exiting");
        write_empty_result_files(&config, &reference, &mut sv_cache, &mut somatic_cache, msi_prediction, &observed)?;
        return Ok(0);
    }
    if !config.target_regions_mode {
        region::write_segments(&format!("{}{}.purple.pre_resegmented_regions.tsv", config.output_dir, tumor_id), &observed)?;
        observed = reseg::run(&observed, config.threads.max(1) as usize);
    }
    let ambiguous_baf = fitting::find_ambiguous_baf(&observed, amber.average_tumor_depth);
    let calculator = fitting::RegionFitCalculator::new(&cobalt.cobalt_chromosomes, &config.fitting, ambiguous_baf);

    let mut gene_copy_numbers: Vec<gene::GeneCopyNumber> = Vec::new();
    let mut copy_numbers: Vec<copynumber::PurpleCopyNumber> = Vec::new();
    let fitted_regions: Vec<region::ObservedRegion>;
    let best_fit: purity_fitter::BestFit;
    let mut purity_adjuster: Option<fitting::PurityAdjuster> = None;
    let mut stream: Option<somatic_stream::SomaticStream> = None;
    let mut stream_result: Option<somatic_stream::SomaticStreamResult> = None;
    let out = &config.output_dir;

    if config.run_tumor() {
        info!("fitting purity");
        let svs = sv_cache.current_variants().to_vec();
        let fitter = purity_fitter::PurityPloidyFitter {
            config: &config,
            reference: &reference,
            amber: &amber,
            cobalt: &cobalt,
            somatics: &somatic_cache.variants,
            svs: &svs,
            calculator: &calculator,
            observed: &observed,
            gender,
        };
        let Some(result) = fitter.run() else {
            warn!("insufficient input Amber and Cobalt data, exiting");
            write_empty_result_files(&config, &reference, &mut sv_cache, &mut somatic_cache, msi_prediction, &observed)?;
            return Ok(0);
        };
        fitted_regions = result.fitted_regions;
        best_fit = result.best_fit;
        copy_numbers = result.copy_numbers;
        sv_cache.infer_missing_variants(&copy_numbers);
        gene_copy_numbers = gene::create_gene_copy_numbers(reference.v38, &reference.ensembl, &copy_numbers, best_fit.fit.ploidy);
        somatic_stream::purity_enrich(&mut somatic_cache.variants, &result.purity_adjuster, &copy_numbers, &fitted_regions);
        info!("enriching somatic variants");
        let mut s = somatic_stream::SomaticStream::new(&config, &reference);
        let r = s.process_and_write(
            &mut somatic_cache,
            Some(&result.purity_adjuster),
            msi_prediction,
            &format!("{out}{tumor_id}.purple.somatic.vcf.gz"),
            &format!("{out}{tumor_id}.purple.somatic.hist.tsv"),
            None,
        )?;
        write_purity_range(&format!("{out}/{tumor_id}.purple.purity.range.tsv"), &best_fit.all_fits)?;
        copynumber::write_copy_numbers(&format!("{out}/{tumor_id}.purple.cnv.somatic.tsv"), &copy_numbers)?;
        fittingsnv::write_peak_model(&format!("{out}/{tumor_id}.purple.somatic.clonality.tsv"), &r.peak_model)?;
        purity_adjuster = Some(result.purity_adjuster);
        stream = Some(s);
        stream_result = Some(r);
    } else {
        best_fit = purity_fitter::BestFit::germline();
        fitted_regions = calculator.fit_regions(best_fit.fit.purity, best_fit.fit.norm_factor, &observed);
    }

    let max_deleted = if config.target_regions_mode { summary::TARGET_REGIONS_MAX_DELETED_GENES } else { summary::MAX_DELETED_GENES };
    let qc = summary::create_qc(
        amber.contamination,
        &best_fit,
        amber_gender,
        cobalt_gender,
        &copy_numbers,
        &gene_copy_numbers,
        &cobalt.cobalt_chromosomes.aberrations,
        amber.average_tumor_depth,
        max_deleted,
        tinc_level,
    );
    let sv_tmb = sv_cache.passing_bnd();
    let purity_context = summary::create_purity(&best_fit, gender, config.run_mode(), config.tumor_only_mode(), config.target_regions_mode, qc, &copy_numbers, stream_result.as_ref(), sv_tmb);
    summary::write_purity_context(out, &tumor_id, &purity_context)?;
    region::write_segments(&format!("{out}{tumor_id}.purple.segment.tsv"), &fitted_regions)?;

    let mut germline_candidates: Vec<germline::GermlineVariant> = Vec::new();
    let mut amp_del: Option<germline::AmpDelFinder> = None;
    if config.run_germline() {
        let reported = stream_result.as_ref().map(|r| r.reported_genes.clone()).unwrap_or_default();
        germline_candidates = germline::process_germline_variants(
            config.reference_id.as_deref().unwrap_or(""),
            &tumor_id,
            files.germline_vcf_file.as_deref().unwrap_or(""),
            &format!("{out}{tumor_id}.purple.germline.vcf.gz"),
            VERSION,
            &reference,
            purity_adjuster.as_ref(),
            &copy_numbers,
            &reported,
        )?;
        germline_sv_cache.annotate_copy_number_info(&fitted_regions, &copy_numbers, purity_context.fit.purity);
        germline_sv_cache.write(&format!("{out}{tumor_id}.purple.sv.germline.vcf.gz"))?;
        let mut finder = germline::AmpDelFinder::new(&reference.driver_genes.driver_gene_list, &reference.ensembl, &reference.germline_amp_del_frequencies);
        let germline_svs = germline_sv_cache.germline_variants().to_vec();
        finder.find_events(&copy_numbers, &fitted_regions, &germline_svs);
        germline::write_germline_amp_dels(&format!("{out}{tumor_id}.purple.germline_amp_del.tsv"), &finder.events)?;
        amp_del = Some(finder);
    }

    if config.run_tumor() {
        sv_cache.write(purity_adjuster.as_ref(), &copy_numbers, config.tumor_only_mode(), amber_gender)?;
    }

    let mut driver_sources: Vec<drivers::DriverSource> = Vec::new();
    if config.target_regions_mode {
        let data = targeted::build(reference.v38, &reference.target_regions, &fitted_regions, &cobalt.ratios, &copy_numbers);
        targeted::write(&format!("{out}{tumor_id}.purple.target_region_cn.tsv"), &data)?;
    }
    let arms = copynumber::chromosome_arm_copy_numbers(&copy_numbers, reference.v38);
    copynumber::write_chromosome_arms(&format!("{out}{tumor_id}.purple.chromosome_arm.tsv"), &arms)?;

    if config.run_drivers {
        info!("generating drivers");
        let mut gene_map: HashMap<String, Vec<gene::GeneCopyNumber>> = HashMap::new();
        for g in &gene_copy_numbers {
            gene_map.entry(g.gene_name.clone()).or_default().push(g.clone());
        }
        if config.run_tumor() {
            let s = stream.as_ref().unwrap();
            let r = stream_result.as_ref().unwrap();
            let mut catalog = s.build_drivers(&somatic_cache.variants, r, &gene_map, &mut driver_sources);
            let qc_status = &purity_context.qc.status;
            let amps = drivers::find_amplifications(qc_status, purity_context.qc.amber_gender, &reference.driver_genes, purity_context.fit.ploidy, &mut gene_copy_numbers, config.target_regions_mode);
            add_gene_driver_sources(&amps, &gene_copy_numbers, &mut driver_sources);
            catalog.extend(amps);
            let dels = drivers::find_deletions(qc_status, purity_context.fit.ploidy, &reference.driver_genes, &mut gene_copy_numbers, config.target_regions_mode);
            add_gene_driver_sources(&dels, &gene_copy_numbers, &mut driver_sources);
            catalog.extend(dels);
            drivers::write_driver_catalog(&format!("{out}{tumor_id}.purple.driver.catalog.somatic.tsv"), &catalog)?;
        }
        if config.run_germline() {
            let mut catalog = germline::germline_drivers(&reference.driver_genes.driver_gene_list, &germline_candidates, &gene_map);
            if let Some(f) = &amp_del {
                catalog.extend(f.drivers.iter().cloned());
            }
            drivers::write_driver_catalog(&format!("{out}{tumor_id}.purple.driver.catalog.germline.tsv"), &catalog)?;
        }
    }
    gene::write_gene_copy_numbers(&format!("{out}{tumor_id}.purple.cnv.gene.tsv"), &gene_copy_numbers)?;

    if !config.germline_mode() && !config.charting.disabled {
        info!("generating charts");
        let mut charting = || -> Result<()> {
            if config.charting.circos_binary.is_some() {
                let plotting: Vec<&variant::Variant> = stream_result.as_ref().map(|r| r.plotting_variants.iter().map(|&i| &somatic_cache.variants[i]).collect()).unwrap_or_default();
                let svs = sv_cache.somatic_variants().to_vec();
                let inputs = charts::CircosInputs {
                    reference_id: config.reference_id.as_deref(),
                    tumor_id: &tumor_id,
                    gender,
                    v38: reference.v38,
                    copy_numbers: &copy_numbers,
                    somatics: plotting,
                    svs: &svs,
                    regions: fitted_regions.iter().filter(|r| r.germline_status != GermlineStatus::Excluded).collect(),
                    amber: &amber,
                    immune_regions: &reference.special_regions.excluded_immune,
                    drivers: if config.skip_circos_genes { &[] } else { &driver_sources },
                };
                charts::write_circos(&config, &inputs)?;
                charts::run_circos(&config, &tumor_id)?;
            }
            // the R plots (copy number, somatic) are not ported
            Ok(())
        };
        if let Err(e) = charting() {
            if config.ignore_plot_errors {
                warn!("charting error: {e}");
            } else {
                error!("charting error: {e}");
                return Ok(1);
            }
        }
    }
    info!("Purple complete, mins({})", fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(0)
}

/// SampleDataFiles.hasValidSampleNames: the VCFs hold the tumor and reference genotypes the run needs.
/// PU4: Java checks the somatic VCFs only when the germline is run; corrected: when the tumor is run.
fn has_valid_sample_names(config: &PurpleConfig) -> Result<bool> {
    let check = |file: &Option<String>| -> Result<bool> {
        let Some(path) = file else { anyhow::bail!("NullPointerException: a VCF file named on the command line is missing") };
        if path.is_empty() || !std::path::Path::new(path).exists() {
            return Ok(true);
        }
        let samples = hmf_common::variant_context::VcfReader::open(path)?.header.samples;
        if config.run_tumor() && !samples.contains(&config.tumor_id) {
            error!("vcf({path}) missing tumorId({}) genotype", config.tumor_id);
            return Ok(false);
        }
        if let Some(reference) = config.reference_id.as_ref().filter(|_| config.run_germline()) {
            if !samples.contains(reference) {
                error!("vcf({path}) missing referenceId({reference}) genotype");
                return Ok(false);
            }
        }
        Ok(true)
    };
    let f = &config.sample_files;
    let germline_valid = !config.run_germline() || (check(&f.germline_sv_vcf_file)? && check(&f.germline_vcf_file)?);
    let check_tumor = if compat::java("PU4") { config.run_germline() } else { config.run_tumor() };
    let tumor_valid = !check_tumor || (check(&f.somatic_sv_vcf_file)? && check(&f.somatic_vcf_file)?);
    Ok(tumor_valid && germline_valid)
}

/// PurpleApplication.addGeneDriverSourceData: the first gene copy number of each driver's gene
fn add_gene_driver_sources(catalog: &[drivers::DriverCatalog], genes: &[gene::GeneCopyNumber], sources: &mut Vec<drivers::DriverSource>) {
    for d in catalog {
        if let Some(g) = genes.iter().find(|g| g.gene_name == d.gene) {
            sources.push(drivers::DriverSource::from_gene(d.clone(), g));
        }
    }
}

/// MsiModelPrediction.read of `<redux_tumor_dir>/<tumor>.redux.msi_prediction.tsv`: -1 when unreadable
fn load_msi_prediction(config: &PurpleConfig) -> Option<f64> {
    let Some(dir) = config.sample_files.redux_tumor_directory.as_deref() else {
        if config.target_regions_mode {
            warn!("target mode running without Redux MSI prediction");
        }
        return None;
    };
    let path = format!("{dir}/{}.redux.msi_prediction.tsv", config.tumor_id);
    let read = || -> Result<f64> {
        let text = std::fs::read_to_string(&path)?;
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let values: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let i = header.iter().position(|h| *h == "PredictedValue").ok_or_else(|| anyhow::anyhow!("no PredictedValue column"))?;
        hmf_common::java_text::parse_double(values.get(i).copied().unwrap_or("")).ok_or_else(|| anyhow::anyhow!("invalid PredictedValue"))
    };
    match read() {
        Ok(v) => Some(v),
        Err(e) => {
            error!("failed to read MSI prediction file({path}): {e}");
            Some(-1.0)
        }
    }
}

/// PurpleApplication.writeEmptyResultFiles
fn write_empty_result_files(
    config: &PurpleConfig,
    reference: &reference::ReferenceData,
    sv_cache: &mut sv_cache::SomaticSvCache,
    somatic_cache: &mut somatic::SomaticVariantCache,
    msi_prediction: Option<f64>,
    regions: &[region::ObservedRegion],
) -> Result<()> {
    let out = &config.output_dir;
    let tumor_id = &config.tumor_id;
    if config.run_tumor() {
        let mut s = somatic_stream::SomaticStream::new(config, reference);
        s.process_and_write(somatic_cache, None, msi_prediction, &format!("{out}{tumor_id}.purple.somatic.vcf.gz"), &format!("{out}{tumor_id}.purple.somatic.hist.tsv"), None)?;
        sv_cache.write(None, &[], config.tumor_only_mode(), Gender::Female)?;
        write_purity_range(&format!("{out}/{tumor_id}.purple.purity.range.tsv"), &[])?;
        copynumber::write_copy_numbers(&format!("{out}/{tumor_id}.purple.cnv.somatic.tsv"), &[])?;
        gene::write_gene_copy_numbers(&format!("{out}{tumor_id}.purple.cnv.gene.tsv"), &[])?;
        fittingsnv::write_peak_model(&format!("{out}/{tumor_id}.purple.somatic.clonality.tsv"), &[])?;
        let context = summary::empty_purity_context(config.run_mode(), config.target_regions_mode);
        summary::write_purity_context(out, tumor_id, &context)?;
        region::write_segments(&format!("{out}{tumor_id}.purple.segment.tsv"), regions)?;
        drivers::write_driver_catalog(&format!("{out}{tumor_id}.purple.driver.catalog.somatic.tsv"), &[])?;
        copynumber::write_chromosome_arms(&format!("{out}{tumor_id}.purple.chromosome_arm.tsv"), &[])?;
    }
    if config.run_germline() {
        if let Some(f) = config.sample_files.germline_sv_vcf_file.as_deref().filter(|f| !f.is_empty()) {
            let c = germline_sv::GermlineSvCache::load(VERSION, f, config.reference_id.as_deref(), tumor_id)?;
            c.write(&format!("{out}{tumor_id}.purple.sv.germline.vcf.gz"))?;
        }
        germline::write_germline_amp_dels(&format!("{out}{tumor_id}.purple.germline_amp_del.tsv"), &[])?;
        drivers::write_driver_catalog(&format!("{out}{tumor_id}.purple.driver.catalog.germline.tsv"), &[])?;
    }
    Ok(())
}

/// FittedPurityRangeFile.write
fn write_purity_range(path: &str, fits: &[fitting::FittedPurity]) -> Result<()> {
    let f = |v: f64| fmt::decimal_format_min(v, 4, 4);
    let mut out = String::from("purity\tnormFactor\tscore\tdiploidProportion\tploidy\tsomaticPenalty\n");
    for p in fits {
        out.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\n", f(p.purity), f(p.norm_factor), f(p.score), f(p.diploid_proportion), f(p.ploidy), f(p.somatic_penalty)));
    }
    std::fs::write(path, out)?;
    Ok(())
}
