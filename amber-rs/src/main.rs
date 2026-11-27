//! AMBER in Rust: a port of hmftools AMBER 4.4 (see ../amber).
//!
//! Usage: `amber-rs <the Java options>`: reads the same files and writes the same outputs as the Java
//! AmberApplication, in its three modes: tumor/normal, tumor-only and germline-only. `-java_compat`
//! reproduces the Java faults the port corrects (REPORT.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod bam;
mod compat;
mod config;
mod constants;
mod contamination;
mod evidence;
mod germline;
mod loh;
mod output;
mod poisson;
mod purity;
mod roh;
mod sites;
mod tumor;
#[allow(dead_code)]
mod types;
mod vcf;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};
// the paths the shared code uses (crate::fmt and the log macros' crate::log)
#[allow(unused_imports)]
use hmf_common::{fmt, log};

use anyhow::Result;

use crate::config::{AmberConfig, ConfigBuilder, Parsed};
use crate::constants::{doubles, QUAL_FILTERED_THRESHOLD};
use crate::evidence::Evidence;
use crate::germline::GermlineAnalysis;
use crate::sites::{human_chromosome_name, ByChr, Locus};
use crate::tumor::{above_qual_filter, TumorAnalysis};
use crate::types::{compare_positions, AmberBaf};

pub const VERSION: &str = "4.4";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    let cb = match ConfigBuilder::parse(&argv) {
        Ok(Parsed::Run(cb)) => cb,
        Ok(Parsed::Exit) => return,
        Err(_) => std::process::exit(1),
    };
    // ConfigUtils.setLogLevel
    let level = cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    if cb.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("Amber version {VERSION}");
    if let Some(kb) = std::env::var("AMBER_FILE_BLOCK_KB").ok().and_then(|v| v.parse::<i32>().ok()) {
        bam::FILE_BLOCK_KB.store(kb.max(0), std::sync::atomic::Ordering::Relaxed);
    }
    let config = match AmberConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    match run(&config) {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            error!("{e:#}");
            std::process::exit(1);
        }
    }
}

/// AmberApplication.run
fn run(config: &AmberConfig) -> Result<i32> {
    let started = std::time::Instant::now();
    let sites = sites::load_amber_sites(config)?;
    if !config.is_valid() {
        error!(" invalid config, exiting");
        // AM6: Java's main ignores the return value of run(), so an invalid configuration ends with exit
        // status 0 and no output. Corrected: status 1.
        return Ok(if compat::java("AM6") { 0 } else { 1 });
    }
    if config.is_tumor_only() {
        run_tumor_only(config, &sites)?;
    } else if config.is_germline_only() {
        run_germline_only(config, &sites)?;
    } else {
        run_normal_mode(config, &sites)?;
    }
    info!("Amber complete, mins({})", fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(0)
}

/// AmberUtils.fromTumorBaf
fn from_tumor_baf(b: &tumor::TumorBaf, loh_recovered: bool) -> AmberBaf {
    let t = &b.tumor.c;
    AmberBaf {
        chromosome: b.tumor.chr.to_string(),
        position: b.tumor.pos,
        tumor_baf: t.alt_support as f64 / (t.alt_support + t.ref_support) as f64,
        tumor_depth: t.read_depth,
        normal_baf: b.normal_alt_support as f64 / (b.normal_alt_support + b.normal_ref_support) as f64,
        normal_depth: b.normal_read_depth,
        loh_recovered,
    }
}

/// AmberUtils.isValid
fn is_valid_baf(b: &AmberBaf) -> bool {
    b.tumor_baf.is_finite() && b.normal_baf.is_finite()
}

/// Collections.sort of GenomePosition items: chromosome rank, then position (stable)
fn sort_by_position<T>(items: &mut [T], key: impl Fn(&T) -> (&str, i32)) {
    items.sort_by(|a, b| {
        let (ca, pa) = key(a);
        let (cb, pb) = key(b);
        compare_positions(ca, pa, cb, pb)
    });
}

/// ResultsWriter.persistBAF: the BAF file, then the segmentation of a tumor's BAFs
fn persist_baf(config: &AmberConfig, bafs: &[AmberBaf]) -> Result<()> {
    output::write_baf(&output::baf_filename(&config.output_dir, config.sample_id()), bafs)?;
    if let Some(tumor_id) = &config.tumor_id {
        if !config.skip_baf_segmentation {
            info!("applying PCF segmentation");
            // AM7: Executors.newFixedThreadPool(threads) throws for -threads 0 (which the option's
            // description allows, and the BAM reading takes as 1), so the run ends here without the
            // segments and the files written after the BAFs. Corrected: 0 is taken as 1.
            if config.threads < 1 && compat::java("AM7") {
                anyhow::bail!("java.lang.IllegalArgumentException (Executors.newFixedThreadPool({}))", config.threads);
            }
            let points: Vec<hmf_common::segmentation::SegmentPoint> = bafs
                .iter()
                .map(|b| hmf_common::segmentation::SegmentPoint {
                    chromosome: b.chromosome.clone(),
                    position: b.position,
                    value: b.tumor_modified_baf(),
                })
                .collect();
            let path = format!("{}{tumor_id}.amber.baf.pcf", config.output_dir);
            hmf_common::segmentation::write_absolute_segments(&points, config.v38, 100.0, config.reader_threads(), &path)?;
        }
    }
    Ok(())
}

fn persist_qc(config: &AmberConfig, consanguinity: f64, contamination: f64, upd: Option<usize>) -> Result<()> {
    let path = format!("{}{}.amber.qc", config.output_dir, config.sample_id());
    output::write_qc(&path, contamination, consanguinity, upd.map(human_chromosome_name))
}

fn persist_snp_check(config: &AmberConfig, snp: &ByChr<Evidence>) -> Result<()> {
    if snp.is_empty() {
        return Ok(());
    }
    let path = format!("{}{}.amber.snp.vcf.gz", config.output_dir, config.primary_reference());
    info!("writing {} germline snp records to {path}", snp.len());
    let mut list: Vec<&Evidence> = snp.values().collect();
    sort_by_position(&mut list, |e| (e.chr, e.pos));
    vcf::write_base_depths(&path, &list, config.primary_reference())
}

fn persist_homozygous_regions(config: &AmberConfig, germline: &GermlineAnalysis) -> Result<()> {
    let path = format!("{}/{}.amber.homozygousregion.tsv", config.output_dir, config.primary_reference());
    output::write_homozygous_regions(&path, &germline.regions_of_homozygosity)
}

fn run_germline_only(config: &AmberConfig, sites: &ByChr<Locus>) -> Result<()> {
    let germline = GermlineAnalysis::run(config, sites)?;
    debug!("BAM reading complete");
    let mut bafs: Vec<AmberBaf> = Vec::new();
    for e in germline.heterozygous.values() {
        // AmberUtils.fromBaseDepth
        let b = AmberBaf {
            chromosome: e.chr.to_string(),
            position: e.pos,
            tumor_baf: -1.0,
            tumor_depth: -1,
            normal_baf: e.c.alt_support as f64 / (e.c.alt_support + e.c.ref_support) as f64,
            normal_depth: e.c.read_depth,
            loh_recovered: false,
        };
        if config.write_unfiltered_germline || is_valid_baf(&b) {
            bafs.push(b);
        }
    }
    sort_by_position(&mut bafs, |b| (&b.chromosome, b.position));
    persist_qc(config, germline.consanguinity_proportion, 0.0, germline.uniparental_disomy)?;
    output::write_version(&config.output_dir)?;
    persist_snp_check(config, &germline.snp_checked)?;
    persist_baf(config, &bafs)?;
    persist_homozygous_regions(config, &germline)
}

fn run_normal_mode(config: &AmberConfig, sites: &ByChr<Locus>) -> Result<()> {
    let germline = GermlineAnalysis::run(config, sites)?;
    let tumor = TumorAnalysis::run(config, &germline.heterozygous, &germline.homozygous)?;
    debug!("BAM reading complete");
    let mut tumor_bafs: Vec<&tumor::TumorBaf> = tumor
        .bafs
        .values()
        .filter(|b| b.tumor.c.read_depth >= config.tumor_min_depth && above_qual_filter(&b.tumor.c, QUAL_FILTERED_THRESHOLD))
        .collect();
    sort_by_position(&mut tumor_bafs, |b| (b.tumor.chr, b.tumor.pos));
    let bafs: Vec<AmberBaf> = tumor_bafs.iter().map(|b| from_tumor_baf(b, false)).filter(is_valid_baf).collect();
    if config.write_tumor_data {
        let rows: Vec<crate::types::PositionEvidence> = tumor_bafs.iter().map(|b| b.tumor.to_position_evidence()).collect();
        let refs: Vec<&crate::types::PositionEvidence> = rows.iter().collect();
        output::write_tumor_data(&output::tumor_data_filename(&config.output_dir, config.tumor_id.as_deref().unwrap()), &refs)?;
    }
    let mut contamination_sites = tumor.contamination;
    let contamination = contamination::calc_contamination(&contamination_sites, bafs.len() as i64);
    persist_qc(config, germline.consanguinity_proportion, contamination, germline.uniparental_disomy)?;
    output::write_version(&config.output_dir)?;
    sort_by_position(&mut contamination_sites, |s| (s.chr, s.pos));
    let path = format!("{}/{}.amber.contamination.tsv", config.output_dir, config.tumor_id.as_deref().unwrap());
    output::write_contamination(&path, &contamination_sites)?;
    persist_snp_check(config, &germline.snp_checked)?;
    persist_baf(config, &bafs)?;
    persist_homozygous_regions(config, &germline)
}

fn run_tumor_only(config: &AmberConfig, sites: &ByChr<Locus>) -> Result<()> {
    let het_sites = sites::tumor_only_sites(config, sites);
    let het = het_sites.map(|s| Evidence::new(s.chr, s.pos, s.ref_, s.alt));
    let tumor = TumorAnalysis::run(config, &het, &ByChr::new())?;
    drop(het);
    debug!("BAM reading complete");

    // the LOH identification works on the tumor BAFs of each autosome, in the germline order
    let mut bafs: ByChr<crate::types::TumorBaf> = tumor.bafs.map(|b| b.to_types());
    drop(tumor);
    let loh = loh::LohIdentifier::new(config.v38, config.tumor_min_depth);
    for (ci, list) in bafs.lists.iter_mut().enumerate() {
        if ci == 0 || ci >= 23 || list.is_empty() {
            continue;
        }
        loh.find_chromosome_loh_segments(human_chromosome_name(ci), list);
    }
    debug!("LOH identification complete");

    let mut filtered: Vec<&crate::types::TumorBaf> = bafs
        .values()
        .filter(|b| {
            let t = &b.tumor;
            let c = evidence::Counts {
                read_depth: t.read_depth,
                base_qual_filtered: t.base_qual_filtered,
                map_qual_filtered: t.map_qual_filtered,
                seq_tech_filtered: t.seq_tech_filtered,
                ..Default::default()
            };
            t.read_depth >= config.tumor_min_depth
                && above_qual_filter(&c, QUAL_FILTERED_THRESHOLD)
                && ((t.ref_support >= config.tumor_only_min_support && t.alt_support >= config.tumor_only_min_support) || t.loh_recovered())
        })
        .collect();
    filtered.sort_by(|a, b| a.tumor.compare(&b.tumor));
    if config.write_tumor_data {
        let rows: Vec<&crate::types::PositionEvidence> = filtered.iter().map(|b| &b.tumor).collect();
        output::write_tumor_data(&output::tumor_data_filename(&config.output_dir, config.tumor_id.as_deref().unwrap()), &rows)?;
    }
    let evidence: Vec<&crate::types::PositionEvidence> = filtered.iter().filter(|b| !b.tumor.loh_recovered()).map(|b| &b.tumor).collect();
    let amber_sites: Vec<crate::types::AmberSite> = sites.values().map(|s| s.to_amber_site()).collect();
    let purity_config = purity::PurityAnalysisConfig {
        tumor_id: config.tumor_id.clone().unwrap(),
        v38: config.v38,
        output_dir: config.output_dir.clone(),
        threads: config.threads.max(0) as usize,
    };
    let analysis = purity::TumorOnlyPurityAnalysis::from_refs(&evidence, &amber_sites, &purity_config);
    drop(evidence);
    drop(amber_sites);
    let noise_floor = analysis.cutoff();
    // Stream.max(Double::compare): NaN above every number, the first of equal values kept
    let contamination = analysis
        .contamination_peaks()
        .iter()
        .map(|p| p.vaf())
        .fold(None, |m: Option<f64>, v| match m {
            Some(m) if m.total_cmp(&v) != std::cmp::Ordering::Less || (m.is_nan() && v.is_nan()) => Some(m),
            _ => Some(v),
        })
        .unwrap_or(0.0);
    debug!("contamination level({}) noiseFloor({})", fmt::java_fixed(contamination, 3), fmt::java_fixed(noise_floor, 3));

    let mut amber_bafs: Vec<AmberBaf> = Vec::new();
    for b in &filtered {
        // AmberUtils.aboveNoiseFloor
        let t = &b.tumor;
        let above = t.loh_recovered() || {
            let ref_freq = t.ref_support as f64 / t.read_depth as f64;
            let alt_freq = t.alt_support as f64 / t.read_depth as f64;
            ref_freq.is_finite()
                && doubles::greater_or_equal(ref_freq, noise_floor)
                && alt_freq.is_finite()
                && doubles::greater_or_equal(alt_freq, noise_floor)
        };
        if !above {
            continue;
        }
        let baf = AmberBaf {
            chromosome: t.chromosome.clone(),
            position: t.position,
            tumor_baf: t.alt_support as f64 / (t.alt_support + t.ref_support) as f64,
            tumor_depth: t.read_depth,
            normal_baf: b.normal_alt_support as f64 / (b.normal_alt_support + b.normal_ref_support) as f64,
            normal_depth: b.normal_read_depth,
            loh_recovered: t.loh_recovered(),
        };
        if baf.tumor_baf.is_finite() {
            amber_bafs.push(baf);
        }
    }
    persist_qc(config, 0.0, contamination, None)?;
    output::write_version(&config.output_dir)?;
    persist_baf(config, &amber_bafs)
}


#[cfg(test)]
mod tests {
    //! End-to-end runs on synthetic BAMs, and the e2e scenario of the Java tests (RunAmberTest).
    use super::*;
    use rust_htslib::bam::{self, header::HeaderRecord, record::CigarString, Format, Header};
    use std::io::Read as _;

    fn strings(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    /// Writes an indexed BAM of single-end reads of 21 bases centred on each site: (position, base, count)
    fn write_bam(path: &str, reads: &[(i32, u8, usize)]) {
        let mut header = Header::new();
        header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
        header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr1").push_tag(b"LN", 1_000_000));
        let mut w = bam::Writer::from_path(path, &header, Format::Bam).unwrap();
        let mut sorted: Vec<(i32, u8, usize)> = reads.to_vec();
        sorted.sort_by_key(|r| r.0);
        let cigar = CigarString::try_from("21M".as_bytes()).unwrap();
        let mut n = 0;
        for (pos, base, count) in sorted {
            for _ in 0..count {
                let mut bases = vec![b'A'; 21];
                bases[10] = base;
                let mut r = bam::Record::new();
                r.set(format!("r{n}").as_bytes(), Some(&cigar), &bases, &[37; 21]);
                r.set_tid(0);
                r.set_pos((pos - 11) as i64);
                r.set_mapq(60);
                r.set_mtid(-1);
                r.set_mpos(-1);
                r.set_flags(0);
                w.write(&r).unwrap();
                n += 1;
            }
        }
        drop(w);
        bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
    }

    fn read_text(path: &str) -> String {
        let mut s = String::new();
        if path.ends_with(".gz") {
            flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap()).read_to_string(&mut s).unwrap();
        } else {
            std::fs::File::open(path).unwrap().read_to_string(&mut s).unwrap();
        }
        s
    }

    /// Tumor/normal: two germline het sites with BAFs, one hom site with tumor alt reads, a SNP-check site
    #[test]
    fn tumor_normal_end_to_end() {
        let dir = std::env::temp_dir().join(format!("amber_rs_e2e_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let d = |f: &str| dir.join(f).to_str().unwrap().to_string();
        std::fs::write(
            d("sites.tsv"),
            "Chromosome\tPosition\tRef\tAlt\tSnpCheck\tFrequency\nchr1\t1000\tA\tG\tfalse\t0.5\nchr1\t5000\tC\tT\tfalse\t0.5\nchr1\t20000\tG\tA\ttrue\t0.5\n",
        )
        .unwrap();
        write_bam(&d("n.bam"), &[(1000, b'A', 10), (1000, b'G', 10), (5000, b'C', 20), (20000, b'G', 9), (20000, b'A', 11)]);
        write_bam(&d("t.bam"), &[(1000, b'A', 5), (1000, b'G', 15), (5000, b'C', 18), (5000, b'T', 2), (20000, b'G', 10), (20000, b'A', 10)]);
        let args = strings(&[
            "-tumor", "T", "-tumor_bam", &d("t.bam"), "-reference", "N", "-reference_bam", &d("n.bam"), "-loci", &d("sites.tsv"),
            "-ref_genome_version", "38", "-output_dir", &d("out"), "-threads", "2", "-skip_baf_segmentation", "-write_tumor_data",
        ]);
        let Parsed::Run(cb) = ConfigBuilder::parse(&args).unwrap() else { panic!() };
        let config = AmberConfig::from(&cb).unwrap();
        assert_eq!(run(&config).unwrap(), 0);
        let out = |f: &str| d(&format!("out/{f}"));
        assert_eq!(
            read_text(&out("T.amber.baf.tsv.gz")),
            "chromosome\tposition\ttumorBAF\ttumorModifiedBAF\ttumorDepth\tnormalBAF\tnormalModifiedBAF\tnormalDepth\tlohRecovered\n\
             chr1\t1000\t0.7500\t0.7500\t20\t0.5000\t0.5000\t20\tfalse\n\
             chr1\t20000\t0.5000\t0.5000\t20\t0.5500\t0.5500\t20\tfalse\n"
        );
        assert_eq!(
            read_text(&out("T.amber.contamination.tsv")),
            "chromosome\tposition\tref\talt\tnormalDepth\tnormalRefSupport\tnormalAltSupport\ttumorDepth\ttumorRefSupport\ttumorAltSupport\n\
             chr1\t5000\tC\tT\t20\t20\t0\t20\t18\t2\n"
        );
        assert_eq!(read_text(&out("T.amber.qc")), "QCStatus\tPASS\nContamination\t0.0000\nConsanguinityProportion\t0.0000\nUniparentalDisomy\tNONE\n");
        assert_eq!(
            read_text(&out("N.amber.homozygousregion.tsv")),
            "chromosome\tstartPosition\tendPosition\tsnpCount\thomCount\thetCount\tfilter\n"
        );
        let vcf = read_text(&out("N.amber.snp.vcf.gz"));
        assert!(vcf.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\nchr1\t20000\t.\tG\tA\t.\t.\t.\tGT:AD:DP\t0/1:9,11:20\n"));
        assert!(std::path::Path::new(&out("N.amber.snp.vcf.gz.tbi")).exists());
        let raw = read_text(&out("T.amber.tumor.raw.tsv.gz"));
        assert_eq!(raw.lines().nth(1), Some("chr1\t1000\tA\tG\t20\t0\t5\t15\t0\t0\t0"));
        assert!(read_text(&out("amber.version")).starts_with("version=4.4\nbuild.date="));
        std::fs::remove_dir_all(&dir).ok();
    }

    /// Germline-only on the same normal: the het sites with tumor BAF -1
    #[test]
    fn germline_only_end_to_end() {
        let dir = std::env::temp_dir().join(format!("amber_rs_e2e_go_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let d = |f: &str| dir.join(f).to_str().unwrap().to_string();
        std::fs::write(d("sites.tsv"), "Chromosome\tPosition\tRef\tAlt\tSnpCheck\nchr1\t1000\tA\tG\tfalse\nchr1\t5000\tC\tT\tfalse\n").unwrap();
        write_bam(&d("n.bam"), &[(1000, b'A', 10), (1000, b'G', 10), (5000, b'C', 20)]);
        let args = strings(&["-reference", "N", "-reference_bam", &d("n.bam"), "-loci", &d("sites.tsv"), "-ref_genome_version", "38", "-output_dir", &d("out")]);
        let Parsed::Run(cb) = ConfigBuilder::parse(&args).unwrap() else { panic!() };
        assert_eq!(run(&AmberConfig::from(&cb).unwrap()).unwrap(), 0);
        let baf = read_text(&d("out/N.amber.baf.tsv.gz"));
        assert_eq!(baf.lines().nth(1), Some("chr1\t1000\t-1.0000\t2.0000\t-1\t0.5000\t0.5000\t20\tfalse"));
        assert_eq!(baf.lines().count(), 2);
        assert!(!std::path::Path::new(&d("out/N.amber.snp.vcf.gz")).exists()); // no SNP-check sites
        std::fs::remove_dir_all(&dir).ok();
    }

    /// An invalid configuration: no output, exit status 1 (0 in Java, AM6)
    #[test]
    fn invalid_config() {
        let dir = std::env::temp_dir().join(format!("amber_rs_e2e_bad_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let d = |f: &str| dir.join(f).to_str().unwrap().to_string();
        std::fs::write(d("sites.tsv"), "Chromosome\tPosition\tRef\tAlt\tSnpCheck\nchr1\t1000\tA\tG\tfalse\n").unwrap();
        let args = strings(&["-reference", "A,B", "-reference_bam", &d("sites.tsv"), "-loci", &d("sites.tsv"), "-output_dir", &d("out")]);
        let Parsed::Run(cb) = ConfigBuilder::parse(&args).unwrap() else { panic!() };
        let config = AmberConfig::from(&cb).unwrap();
        assert_eq!(run(&config).unwrap(), 1);
        assert_eq!(compat::with_java(true, || run(&config).unwrap()), 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    /// Runs a RunAmberTest scenario (tumor-only on the Java test's BAM) and returns the output directory
    fn run_java_scenario(name: &str, target_bed: Option<&str>) -> String {
        let res = concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/e2e");
        let dir = std::env::temp_dir().join(format!("amber_rs_e2e_java_{name}_{}_{}", target_bed.unwrap_or("wgs"), std::process::id()));
        let out = dir.to_str().unwrap().to_string();
        let mut args = strings(&[
            "-loci", &format!("{res}/{name}.tsv"), "-output_dir", &out, "-ref_genome_version", "38", "-tumor",
            &format!("{name}_tumor"), "-tumor_bam", &format!("{res}/{name}.tumor.bam"),
        ]);
        if let Some(bed) = target_bed {
            args.extend(strings(&["-target_regions_bed", &format!("{res}/{bed}")]));
        }
        let Parsed::Run(cb) = ConfigBuilder::parse(&args).unwrap() else { panic!() };
        assert_eq!(run(&AmberConfig::from(&cb).unwrap()).unwrap(), 0);
        out
    }

    /// The sites of TwoChromosomes.tumor.reads.tsv with their (ref, alt, depth) counts
    fn two_chromosome_expectations() -> Vec<(String, i32, i32, i32, i32)> {
        let res = concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/e2e");
        let sites: Vec<(String, i32)> = read_text(&format!("{res}/TwoChromosomes.tsv"))
            .lines()
            .skip(1)
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[0].to_string(), f[1].parse().unwrap())
            })
            .collect();
        let mut expected = Vec::new();
        for l in read_text(&format!("{res}/TwoChromosomes.tumor.reads.tsv")).lines().filter(|l| !l.is_empty()) {
            let (label, counts) = l.split_once(' ').unwrap();
            let (chr, index) = label.split_once('_').unwrap();
            let c: Vec<i32> = counts.split(',').map(|x| x.parse().unwrap()).collect();
            let chr_sites: Vec<&(String, i32)> = sites.iter().filter(|s| s.0 == chr).collect();
            let s = chr_sites[index.parse::<usize>().unwrap()];
            expected.push((s.0.clone(), s.1, c[0], c[1], c[0] + c[1] + c[2]));
        }
        expected
    }

    fn check_bafs(out: &str, name: &str, expected: &[(String, i32, i32, i32, i32)]) {
        let baf = read_text(&format!("{out}/{name}_tumor.amber.baf.tsv.gz"));
        let rows: Vec<Vec<String>> = baf.lines().skip(1).map(|l| l.split('\t').map(|s| s.to_string()).collect()).collect();
        assert_eq!(rows.len(), expected.len());
        for (chr, pos, r, a, depth) in expected {
            let row = rows.iter().find(|x| &x[0] == chr && x[1] == pos.to_string()).unwrap();
            assert_eq!(row[2], output::decimal_format_4(*a as f64 / (a + r) as f64));
            assert_eq!(row[4], depth.to_string());
        }
    }

    /// RunAmberTest.twoChromosomes and newSegmenter: every site's BAF as the read counts give it, and the
    /// segments.
    #[test]
    fn java_e2e_two_chromosomes() {
        let out = run_java_scenario("TwoChromosomes", None);
        check_bafs(&out, "TwoChromosomes", &two_chromosome_expectations());
        let pcf = read_text(&format!("{out}/TwoChromosomes_tumor.amber.baf.pcf"));
        let segments: Vec<Vec<&str>> = pcf.lines().skip(1).map(|l| l.split('\t').collect()).collect();
        let chr1: Vec<&Vec<&str>> = segments.iter().filter(|s| s[0] == "chr1").collect();
        assert!(segments.iter().all(|s| s[0] == "chr1" || s[0] == "chr2") && segments.iter().any(|s| s[0] == "chr2"));
        assert_eq!(chr1.len(), 3);
        assert_eq!(chr1[0][1], "49554");
        std::fs::remove_dir_all(&out).ok();
    }

    /// RunAmberTest.targetedMode: only the sites within 300 bases of the panel's regions
    #[test]
    fn java_e2e_targeted_mode() {
        let out = run_java_scenario("TwoChromosomes", Some("PanelA.bed"));
        let regions = [("chr1", 52000 - 300, 60000 + 300), ("chr2", 11000 - 300, 12000 + 300)];
        let expected: Vec<_> = two_chromosome_expectations()
            .into_iter()
            .filter(|e| regions.iter().any(|(c, s, en)| e.0 == *c && e.1 >= *s && e.1 <= *en))
            .collect();
        check_bafs(&out, "TwoChromosomes", &expected);
        std::fs::remove_dir_all(&out).ok();
    }

    /// RunAmberTest.handleEmptyData: no BAFs, a segments file with no segments
    #[test]
    fn java_e2e_no_bafs() {
        let out = run_java_scenario("NoBafs", None);
        assert_eq!(read_text(&format!("{out}/NoBafs_tumor.amber.baf.tsv.gz")).lines().count(), 1);
        assert_eq!(read_text(&format!("{out}/NoBafs_tumor.amber.baf.pcf")).lines().count(), 1);
        std::fs::remove_dir_all(&out).ok();
    }
}
