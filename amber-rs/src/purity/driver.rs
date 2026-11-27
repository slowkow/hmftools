//! The Rust side of the comparison with Java (REPORT.md): reads a -write_tumor_data TSV and the loci file and
//! runs LohIdentifier and TumorOnlyPurityAnalysis as AmberApplication.runTumorOnly does after the BAM reading,
//! writing the same three files as the Java harness (work/amber/tumoronly/harness/TumorOnlyHarness.java).
//!
//! Run as an ignored test, with the inputs in the environment:
//!   AMBER_TO_TSV, AMBER_TO_LOCI, AMBER_TO_OUT (output prefix), AMBER_TO_MIN_DEPTH (25),
//!   AMBER_TO_MIN_SUPPORT (2), AMBER_TO_THREADS (2), AMBER_TO_VERSION (V38), AMBER_TO_SKIP_PURITY=1 (stop before
//!   the purity analysis, to measure the rest), and AMBER_JAVA_COMPAT=1 for Java behaviour:  <test binary> --ignored --exact purity::driver::tumor_only_driver

use crate::loh::{above_qual_filter, LohIdentifier};
use crate::purity::{PeakClass, PurityAnalysisConfig, TumorOnlyPurityAnalysis};
use crate::types::{AmberBase, AmberSite, PositionEvidence, TumorBaf};
use std::io::{BufRead, BufReader, BufWriter, Write};

fn reader(path: &str) -> BufReader<Box<dyn std::io::Read>> {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let inner: Box<dyn std::io::Read> = if path.ends_with(".gz") { Box::new(flate2::read::MultiGzDecoder::new(file)) } else { Box::new(file) };
    BufReader::with_capacity(1 << 16, inner)
}

fn hex(d: f64) -> String {
    format!("{:x}", d.to_bits())
}

fn env(name: &str, default: &str) -> String {
    std::env::var(name).unwrap_or_else(|_| default.to_string())
}

/// AmberSitesFile.loadFile, for the columns the purity analysis uses (human chromosomes only)
pub fn load_sites(path: &str) -> Vec<AmberSite> {
    let mut lines = reader(path).lines();
    let header = lines.next().unwrap().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let at = |name: &str| cols.iter().position(|c| *c == name);
    let (c, p, r, a, s) = (at("Chromosome").unwrap(), at("Position").unwrap(), at("Ref").unwrap(), at("Alt").unwrap(), at("SnpCheck").unwrap());
    let f = at("Frequency");
    let mut sites = Vec::new();
    for line in lines {
        let line = line.unwrap();
        let v: Vec<&str> = line.split('\t').collect();
        if !hmf_common::regions::is_human_chromosome(v[c]) {
            continue;
        }
        let frequency = f.map(|i| v[i].parse::<f64>().unwrap()).unwrap_or(0.0);
        sites.push(AmberSite::with_frequency(v[c], v[p].parse().unwrap(), v[r], v[a], v[s].eq_ignore_ascii_case("true"), frequency));
    }
    sites
}

fn read_tumor(path: &str) -> Vec<TumorBaf> {
    let mut lines = reader(path).lines();
    let header = lines.next().unwrap().unwrap();
    let cols: Vec<String> = header.split('\t').map(String::from).collect();
    let at = |name: &str| cols.iter().position(|c| c == name).unwrap();
    let idx = [
        at("Chromosome"),
        at("Position"),
        at("Ref"),
        at("Alt"),
        at("ReadDepth"),
        at("IndelCount"),
        at("RefCount"),
        at("AltCount"),
        at("BaseQualFiltered"),
        at("MapQualFiltered"),
        at("SeqTechFiltered"),
    ];
    let mut bafs = Vec::new();
    let mut last_chromosome = String::new();
    for line in lines {
        let line = line.unwrap();
        let v: Vec<&str> = line.split('\t').collect();
        let int = |k: usize| v[idx[k]].parse::<i32>().unwrap();
        if v[idx[0]] != last_chromosome {
            last_chromosome = v[idx[0]].to_string();
        }
        let mut e = PositionEvidence::new(&last_chromosome, int(1), AmberBase::from_str(v[idx[2]]).unwrap(), AmberBase::from_str(v[idx[3]]).unwrap());
        e.read_depth = int(4);
        e.indel_count = int(5);
        e.ref_support = int(6);
        e.alt_support = int(7);
        e.base_qual_filtered = int(8);
        e.map_qual_filtered = int(9);
        e.seq_tech_filtered = int(10);
        bafs.push(TumorBaf { normal_read_depth: 0, normal_ref_support: 0, normal_alt_support: 0, tumor: e });
    }
    bafs
}

#[test]
#[ignore]
fn tumor_only_driver() {
    let _ = crate::compat::take_options(&[]);
    if env("AMBER_DEBUG", "0") == "1" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    let tsv = env("AMBER_TO_TSV", "");
    let loci = env("AMBER_TO_LOCI", "");
    let out = env("AMBER_TO_OUT", "");
    let min_depth: i32 = env("AMBER_TO_MIN_DEPTH", "25").parse().unwrap();
    let min_support: i32 = env("AMBER_TO_MIN_SUPPORT", "2").parse().unwrap();
    let threads: usize = env("AMBER_TO_THREADS", "2").parse().unwrap();
    let v38 = env("AMBER_TO_VERSION", "V38") == "V38";
    println!("java compat: {}", crate::compat::java_all());

    let t0 = std::time::Instant::now();
    let sites = load_sites(&loci);
    let t_sites = t0.elapsed();
    let mut bafs = read_tumor(&tsv);
    let t_read = t0.elapsed();

    // findLohSegments: per chromosome (Java's multimap keys), the sites in file order
    let loh = LohIdentifier::new(v38, min_depth);
    let mut by_chromosome: Vec<(usize, usize)> = Vec::new();
    let mut seen = std::collections::HashMap::new();
    for i in 0..bafs.len() {
        let key = crate::purity::genome::chromosome_ordinal(&bafs[i].tumor.chromosome).expect("human chromosome");
        if i == 0 || bafs[i].tumor.chromosome != bafs[i - 1].tumor.chromosome {
            assert!(seen.insert(key, ()).is_none(), "the TSV has the sites of one chromosome together");
            by_chromosome.push((i, i + 1));
        } else {
            by_chromosome.last_mut().unwrap().1 = i + 1;
        }
    }
    for &(s, e) in &by_chromosome {
        let chromosome = bafs[s].tumor.chromosome.clone();
        loh.find_chromosome_loh_segments(&chromosome, &mut bafs[s..e]);
    }
    let t_loh = t0.elapsed();

    let mut w = BufWriter::with_capacity(1 << 16, std::fs::File::create(format!("{out}.loh.tsv")).unwrap());
    writeln!(w, "Chromosome\tPosition\tRef\tAlt\tValidBafHigh\tValidBafMed\tValidBafLow\tValidRateHigh\tValidRateMed\tValidRateLow\tLohCandidate\tLohRecovered").unwrap();
    for b in &bafs {
        let e = &b.tumor;
        write!(w, "{}\t{}\t{}\t{}", e.chromosome, e.position, e.ref_.as_str(), e.alt.as_str()).unwrap();
        match &e.loh_data {
            None => writeln!(w, "\t-").unwrap(),
            Some(d) => writeln!(
                w,
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                d.valid_baf_high,
                d.valid_baf_medium,
                d.valid_baf_low,
                hex(d.valid_rate_high),
                hex(d.valid_rate_medium),
                hex(d.valid_rate_low),
                d.loh_candidate,
                d.loh_recovered
            )
            .unwrap(),
        }
    }
    w.flush().unwrap();

    // AmberApplication.runTumorOnly: the depth, quality and support filters (the TSV is in genome order)
    let raw: Vec<&TumorBaf> = bafs
        .iter()
        .filter(|b| b.tumor.read_depth >= min_depth)
        .filter(|b| above_qual_filter(&b.tumor, 0.15))
        .filter(|b| (b.tumor.ref_support >= min_support && b.tumor.alt_support >= min_support) || b.tumor.loh_recovered())
        .collect();
    let evidence: Vec<&PositionEvidence> = raw.iter().map(|b| &b.tumor).filter(|e| !e.loh_recovered()).collect();
    let t_filter = t0.elapsed();
    if env("AMBER_TO_SKIP_PURITY", "0") == "1" {
        // to measure the memory of everything but the purity analysis
        println!("purity analysis skipped; evidence {}", evidence.len());
        return;
    }
    let config = PurityAnalysisConfig { tumor_id: "TUMOR".into(), v38, output_dir: ".".into(), threads };
    let analysis = TumorOnlyPurityAnalysis::from_refs(&evidence, &sites, &config);
    let t_purity = t0.elapsed();

    let report = analysis.report();
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.purity.tsv")).unwrap());
    writeln!(w, "evidence\t{}", report.evidence).unwrap();
    writeln!(w, "filtered\t{}", report.filtered).unwrap();
    writeln!(w, "hetVariants\t{}", report.het_variants).unwrap();
    writeln!(w, "baselineFrequency\t{}", hex(report.baseline_frequency)).unwrap();
    for (g, l) in report.levels.iter().enumerate() {
        // Java's mPointsTested holds the testable points only once there are enough of them
        let tested = if l.testable < crate::purity::peak::MINIMUM_CAPTURED_POINTS { 0 } else { l.testable };
        writeln!(w, "level\t{g}\t{}\t{}\t{tested}\t{}\t{}\t{}\t{}\t{}", hex(l.level), hex(l.step), l.homozygous, l.heterozygous, l.captured, hex(l.score), l.message).unwrap();
    }
    let mut any_overflow = false;
    for (peak, class) in &report.peaks {
        let s = peak.stats().unwrap();
        let g = &s.gnomad;
        write!(
            w,
            "peak\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            hex(peak.vaf()),
            hex(s.score),
            s.homozygous_count,
            s.heterozygous_count,
            s.captured,
            g.low_count,
            hex(g.low_mean),
            g.high_count,
            hex(g.high_mean),
            g.pass
        )
        .unwrap();
        if *class == PeakClass::FailedGnomad {
            writeln!(w, "\t-\t-\t-\t-\tFailedGnomad").unwrap();
        } else {
            any_overflow |= s.auc_overflow;
            writeln!(w, "\t{}\t{}\t{}\t{}\t{:?}", hex(s.homozygous_proportion), hex(s.chr_arm_auc), hex(s.mutation_auc), hex(s.chr_arm_lower_bound), class).unwrap();
        }
    }
    writeln!(w, "cutoff\t{}", hex(analysis.cutoff())).unwrap();
    writeln!(w, "contamination\t{}", hex(analysis.contamination())).unwrap();
    let list = |peaks: &[crate::purity::CandidatePeak]| peaks.iter().map(|p| hex(p.vaf())).collect::<Vec<_>>().join(",");
    writeln!(w, "copyNumberPeaks\t{}", list(analysis.copy_number_peaks())).unwrap();
    writeln!(w, "contaminationPeaks\t{}", list(analysis.contamination_peaks())).unwrap();
    w.flush().unwrap();

    // the BAF file's sites: AmberUtils.aboveNoiseFloor, fromTumorBaf, a finite tumor BAF
    let noise_floor = analysis.cutoff();
    let mut w = BufWriter::with_capacity(1 << 16, std::fs::File::create(format!("{out}.bafs.tsv")).unwrap());
    for b in &raw {
        let e = &b.tumor;
        let above = e.loh_recovered() || {
            let (rf, af) = (b.ref_frequency(), b.alt_frequency());
            rf.is_finite() && rf - noise_floor > -1e-10 && af.is_finite() && af - noise_floor > -1e-10
        };
        if !above {
            continue;
        }
        let tumor_baf = e.alt_support as f64 / (e.alt_support + e.ref_support) as f64;
        if tumor_baf.is_finite() {
            writeln!(w, "{}\t{}\t{}\t{}\t{}", e.chromosome, e.position, hex(tumor_baf), e.read_depth, e.loh_recovered()).unwrap();
        }
    }
    w.flush().unwrap();
    let t_end = t0.elapsed();
    println!(
        "times(ms): sites {}, read {}, loh {}, filter {}, purity {}, write {}; evidence {}, int overflow in an AUC (AM53): {}",
        t_sites.as_millis(),
        (t_read - t_sites).as_millis(),
        (t_loh - t_read).as_millis(),
        (t_filter - t_loh).as_millis(),
        (t_purity - t_filter).as_millis(),
        (t_end - t_purity).as_millis(),
        evidence.len(),
        any_overflow
    );
}
