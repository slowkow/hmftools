//! CIDER in Rust: a port of hmftools CIDER 1.2 (see ../../hmftools/cider), which finds the IG/TCR CDR3
//! sequences of a sample from the reads at the IG/TCR loci of its BAM.
//!
//! Usage: `cider-rs [com.hartwig.hmftools.cider.CiderApplication] <the Java options>`: reads the same files
//! and writes the same outputs as the Java CiderApplication. JVM options (-Xmx..., -XX:..., -D...=...) are
//! dropped, as oncoanalyser passes them to the Java tool's wrapper. `-java_compat` reproduces the Java
//! faults the port corrects (REPORT.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod alignment;
mod annotate;
mod bam;
mod blosum;
mod compat;
mod config;
mod constants;
mod genes;
mod layout;
mod layout_builder;
mod output;
mod primer;
mod read;
mod screener;
mod vdj;
#[cfg(test)]
mod vdj_tests;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};
// the paths the shared code's macros use (crate::log)
#[allow(unused_imports)]
use hmf_common::{fmt, log};

use std::collections::HashMap;
use std::sync::Arc;

use anyhow::{bail, Result};

use crate::alignment::AlignmentAnnotator;
use crate::blosum::AnchorBlosumSearcher;
use crate::config::{CiderParams, ConfigBuilder, Parsed};
use crate::constants::*;
use crate::genes::{GeneDatastore, VjGeneType};
use crate::layout::LayoutRead;
use crate::layout_builder::{LayoutBuildResult, VjReadLayoutBuilder};
use crate::screener::{CiderReadScreener, VjReadCandidate};
use crate::vdj::{parallel_map, VdjSequenceBuilder};

const MAIN_CLASS: &str = "com.hartwig.hmftools.cider.CiderApplication";

fn is_jvm_option(a: &str) -> bool {
    a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-Xss") || a.starts_with("-XX:") || (a.starts_with("-D") && a.contains('='))
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).filter(|a| !is_jvm_option(a)).collect();
    // the bioconda wrapper runs `java -cp cider.jar <class> ...` when the first argument is a class name
    let rest: Vec<String> = match argv.first().map(|s| s.as_str()) {
        Some(MAIN_CLASS) => argv[1..].to_vec(),
        Some(c) if c.starts_with("com.hartwig.") => {
            eprintln!("cider (Rust port): main class {c} is not ported; only {MAIN_CLASS}");
            std::process::exit(2);
        }
        _ => argv,
    };
    let rest = match compat::take_options(&rest) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    let cb = match ConfigBuilder::parse(&rest) {
        Ok(Parsed::Run(cb)) => cb,
        Ok(Parsed::Exit) => return,
        Err(_) => std::process::exit(1),
    };
    // ConfigUtils.setLogLevel
    let level = cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    if cb.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("{} version {}", config::APP_NAME, config::VERSION);
    let params = match CiderParams::from(&cb) {
        Ok(p) => p,
        Err(e) => {
            error!("{e:#}");
            std::process::exit(1);
        }
    };
    match run(&params) {
        Ok(()) => std::process::exit(0),
        Err(e) => {
            error!("{e:#}");
            std::process::exit(1);
        }
    }
}

/// CiderApplication.run
fn run(params: &CiderParams) -> Result<()> {
    let started = std::time::Instant::now();
    std::fs::create_dir_all(&params.output_dir).ok();
    // CI7: Executors.newFixedThreadPool(threads) throws for -threads 0 or less (after cider.bam is written),
    // though the option says 0 is single-threaded. Corrected: taken as 1.
    let java_bad_threads = params.thread_count < 1 && compat::java("CI7");
    let threads = params.thread_count.max(1) as usize;

    let ds = GeneDatastore::load(params.ref_genome_version)?;
    let candidate_searcher = AnchorBlosumSearcher::new(&ds, CANDIDATE_MIN_PARTIAL_ANCHOR_AA_LENGTH);
    let screener = CiderReadScreener::new(&ds, &candidate_searcher, MIN_READ_CDR3_OVERLAP, params.approx_max_fragment_length);

    // readBamFile
    let phase = std::time::Instant::now();
    let regions = bam::query_regions(&ds, params.approx_max_fragment_length)?;
    debug!("Processing {} potential sites in bam {}", regions.len(), params.bam_path);
    let screened = bam::read_bam(&params.bam_path, params.ref_genome_path.as_deref(), &regions, &screener, threads)?;
    debug!("PerfStats: reading the BAM {:.2}s", phase.elapsed().as_secs_f64());

    if params.write_filtered_bam {
        let phase = std::time::Instant::now();
        bam::write_cider_bam(&params.file(".cider.bam"), &screened.header_text, &screened.targets, &screened.all_matched_reads)?;
        debug!("PerfStats: writing cider.bam {:.2}s", phase.elapsed().as_secs_f64());
    }
    if java_bad_threads {
        bail!("Exception in thread \"main\" java.lang.IllegalArgumentException (Executors.newFixedThreadPool({}))", params.thread_count);
    }

    // buildLayouts: per gene type, in gene type order, then numbered
    let phase = std::time::Instant::now();
    let builder = VjReadLayoutBuilder::new(params.num_bases_to_trim, params.min_base_quality);
    let candidates: Vec<Arc<VjReadCandidate>> = screened.vj_read_candidates.into_iter().map(Arc::new).collect();
    let per_type: Vec<(VjGeneType, Vec<Arc<VjReadCandidate>>)> =
        VjGeneType::ALL.iter().map(|&t| (t, candidates.iter().filter(|c| c.gene_type == t).cloned().collect())).collect();
    let built = parallel_map(per_type, threads, |(t, c)| {
        builder.build_layouts(t, &c, LAYOUT_MIN_READ_OVERLAP_BASES, params.max_read_count_per_gene, params.max_low_qual_base_fraction).map(|r| (t, r))
    })?;
    let mut results: Vec<(VjGeneType, LayoutBuildResult)> = built;
    let mut next_id = 1;
    for (_, r) in results.iter_mut() {
        for l in r.layouts.iter_mut() {
            Arc::get_mut(l).expect("unshared layout").id = next_id.to_string();
            next_id += 1;
        }
    }
    output::write_layouts(&params.file(".cider.layout.gz"), &results)?;
    debug!("PerfStats: layouts {:.2}s", phase.elapsed().as_secs_f64());

    // the mates' layout reads, by read name
    let mut mate_layout_reads: HashMap<String, Vec<LayoutRead>> = HashMap::new();
    for c in screened.vj_mate_candidates {
        let name = c.read.name.clone();
        let c = Arc::new(c);
        let entry = mate_layout_reads.entry(name).or_default();
        if let Some(r) = builder.candidate_to_layout_read(&c) {
            entry.push(r);
        }
    }

    let phase = std::time::Instant::now();
    let vdj_searcher = AnchorBlosumSearcher::new(&ds, VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH);
    let vdj_builder =
        VdjSequenceBuilder { adaptor: &builder, ds: &ds, searcher: &vdj_searcher, min_base_quality: params.min_base_quality as u8, min_overlapped_bases: MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES };
    let layouts_by_type: Vec<(VjGeneType, Vec<Arc<layout::ReadLayout>>)> = results.iter().map(|(t, r)| (*t, r.layouts.clone())).collect();
    let vdjs = vdj_builder.build_vdj_sequences(&layouts_by_type, &mate_layout_reads, threads)?;
    debug!("PerfStats: VDJ sequences {:.2}s", phase.elapsed().as_secs_f64());

    let mut primer_matches = Vec::new();
    if let Some(path) = &params.primer_csv {
        let primers = primer::load_primers(path)?;
        primer_matches = primer::match_vdj_primers(&vdjs, &primers, params.primer_mismatch_max);
        output::write_primer_matches(&params.file(".cider.primer_match.tsv"), &vdjs, &primer_matches)?;
    }

    let phase = std::time::Instant::now();
    let alignments = match (&params.ref_genome_path, &params.bwa_index_image_path) {
        (Some(ref_path), Some(index_path)) => {
            // sequences that match the reference are not aligned
            let filtered: Vec<(usize, &vdj::VdjSequence)> = vdjs.iter().enumerate().filter(|(_, v)| !annotate::vdj_matches_ref(v)).collect();
            let annotator = AlignmentAnnotator::new(params.ref_genome_version, &format!("{ref_path}.dict"), index_path)?;
            let a = annotator.run_annotate(&filtered, threads)?;
            output::write_alignment_matches(&params.file(".cider.alignment_match.tsv.gz"), &vdjs, &a)?;
            a
        }
        _ => Vec::new(),
    };
    debug!("PerfStats: alignment {:.2}s", phase.elapsed().as_secs_f64());

    let counts = primer::primer_counts(vdjs.len(), &primer_matches);
    let annotations = annotate::sort_and_annotate(&vdjs, &alignments, &counts)?;
    let annotations: Vec<annotate::VdjAnnotation> = annotations
        .into_iter()
        .filter(|a| annotate::passes_hard_filters(a, a.alignment.map(|k| &alignments[k]), params.report_match_ref_seq))
        .collect();
    output::write_vdj_sequences(&params.file(".cider.vdj.tsv.gz"), &vdjs, &annotations, &alignments)?;
    output::write_locus_stats(&params.file(".cider.locus_stats.tsv"), &results, &vdjs, &annotations)?;
    info!("Cider complete, mins({})", fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(())
}
