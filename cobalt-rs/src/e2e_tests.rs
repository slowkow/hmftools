//! The scenarios of the Java ProcessBamTest (cobalt/src/test/java/.../e2e), run in-process on the Java
//! tests' synthetic BAMs (cobalt/src/test/resources/bam, found beside this crate in the hmftools repository
//! or where COBALT_TEST_RESOURCES points; the tests are skipped without them). The GC profile, panel and
//! diploid files are written as the Java test helpers write them. Java's own test class is @Ignore'd; each
//! scenario here checks the values the Java test asserts, in Java mode.
//!
//! With COBALT_E2E_KEEP=<dir>, each scenario's inputs and command line (`args`, one item per line) are kept
//! in <dir>/<scenario>/ so that scripts/e2e_compare.sh can run Java and the binary on them and compare
//! every output file.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn resources() -> Option<PathBuf> {
    let candidates = [
        std::env::var("COBALT_TEST_RESOURCES").ok().map(PathBuf::from),
        Some(Path::new(env!("CARGO_MANIFEST_DIR")).join("../cobalt/src/test/resources")),
    ];
    let found = candidates.into_iter().flatten().find(|d| d.join("bam").is_dir());
    if found.is_none() {
        eprintln!("COBALT test resources not found: test skipped");
    }
    found
}

/// A GC profile section (GcFileSection): lines from start to stop (0-based) every 1,000 bases
struct Gc {
    chr: &'static str,
    start: i32,
    stop: i32,
    gc: f64,
    mappable: f64,
}

fn gc(chr: &'static str, start: i32, stop: i32, gc: f64) -> Gc {
    Gc { chr, start, stop, gc, mappable: 1.0 }
}
fn mappable(chr: &'static str, start: i32, stop: i32, mappable: f64) -> Gc {
    Gc { chr, start, stop, gc: 0.5, mappable }
}

/// A panel section (PanelFileSection): lines at start + 1 to stop + 1 every 1,000 bases
struct Panel {
    chr: &'static str,
    start: i32,
    stop: i32,
    enrichment: f64,
}

fn chr_rank(chr: &str) -> usize {
    crate::chromosome::Chr::from_name(chr).unwrap().index()
}

fn section_positions(start: i32, stop: i32) -> Vec<i32> {
    let steps = (stop - start) / 1000 + 1;
    (0..steps).map(|i| start + i * 1000).collect()
}

/// SectionalData.write: the sections in chromosome order (stable), lines joined by "\n"
fn write_sections(path: &Path, header: Option<&str>, mut sections: Vec<(usize, Vec<String>)>) {
    sections.sort_by_key(|s| s.0);
    let mut lines: Vec<String> = header.map(|h| vec![h.to_string()]).unwrap_or_default();
    for (_, l) in sections {
        lines.extend(l);
    }
    std::fs::write(path, lines.join("\n")).unwrap();
}

struct Scenario {
    name: &'static str,
    tumor_bam: Option<&'static str>,
    reference_bam: Option<&'static str>,
    gc: Vec<Gc>,
    panel: Option<Vec<Panel>>,
    /// diploid BED sections: chromosome, start, stop
    diploid: Option<Vec<(&'static str, i32, i32)>>,
}

/// The ratio file's rows by chromosome (all nine columns, the chromosome first)
type Rows = BTreeMap<String, Vec<Vec<f64>>>;

struct Output {
    dir: PathBuf,
    rows: Rows,
}

impl Output {
    fn chr(&self, chr: &str) -> &Vec<Vec<f64>> {
        &self.rows[chr]
    }
    /// tumorGCRatio of chromosome `chr` window `i`
    fn tumor_ratio(&self, chr: &str, i: usize) -> f64 {
        self.chr(chr)[i][4]
    }
    fn exists(&self, name: &str) -> bool {
        self.dir.join(name).exists()
    }
}

const SAMPLE: &str = "SAMPLE";
const REFERENCE: &str = "SAMPLE_reference";

fn run(s: Scenario) -> Option<Output> {
    let res = resources()?;
    let base = match std::env::var("COBALT_E2E_KEEP") {
        Ok(d) => PathBuf::from(d).join(s.name),
        Err(_) => std::env::temp_dir().join(format!("cobalt_rs_e2e_{}_{}", std::process::id(), s.name)),
    };
    let _ = std::fs::remove_dir_all(&base);
    let out = base.join("output");
    std::fs::create_dir_all(&out).unwrap();
    // GcProfilesUtilities: chr<n>, position, GC (%.2f), 1, mappable (%.2f)
    let gc_path = base.join("GC_profile.1000bp.38.cnp");
    let gc_sections: Vec<(usize, Vec<String>)> = s
        .gc
        .iter()
        .map(|g| {
            let lines = section_positions(g.start, g.stop)
                .into_iter()
                .map(|p| format!("chr{}\t{p}\t{}\t1\t{}", g.chr, hmf_common::fmt::java_fixed(g.gc, 2), hmf_common::fmt::java_fixed(g.mappable, 2)))
                .collect();
            (chr_rank(g.chr), lines)
        })
        .collect();
    write_sections(&gc_path, None, gc_sections);
    let mut args: Vec<String> = vec![
        "-gc_profile".into(),
        gc_path.to_string_lossy().into(),
        "-pcf_gamma".into(),
        "50".into(),
        "-output_dir".into(),
        out.to_string_lossy().into(),
        "-ref_genome_version".into(),
        if s.name.ends_with("_v37") { "37".into() } else { "38".into() },
    ];
    let bam = |name: &str| res.join("bam").join(format!("{name}.sorted.bam")).to_string_lossy().to_string();
    if let Some(t) = s.tumor_bam {
        args.extend(["-tumor".into(), SAMPLE.into(), "-tumor_bam".into(), bam(t)]);
    }
    if let Some(r) = s.reference_bam {
        args.extend(["-reference".into(), REFERENCE.into(), "-reference_bam".into(), bam(r)]);
    }
    if let Some(panel) = &s.panel {
        let path = base.join("ThePanel.tsv");
        let sections: Vec<(usize, Vec<String>)> = panel
            .iter()
            .map(|p| (chr_rank(p.chr), section_positions(p.start, p.stop).into_iter().map(|x| format!("chr{}\t{}\t{:.9}", p.chr, x + 1, p.enrichment)).collect()))
            .collect();
        write_sections(&path, Some("chromosome\tposition\trelativeEnrichment"), sections);
        args.extend(["-target_region_norm_file".into(), path.to_string_lossy().into()]);
    }
    if let Some(d) = &s.diploid {
        let path = base.join("diploid.bed");
        let sections: Vec<(usize, Vec<String>)> = d.iter().map(|(chr, a, b)| (chr_rank(chr), vec![format!("chr{chr}\t{a}\t{b}")])).collect();
        write_sections(&path, None, sections);
        args.extend(["-tumor_only_diploid_bed".into(), path.to_string_lossy().into()]);
    }
    std::fs::write(base.join("args"), args.join("\n") + "\n").unwrap();
    let status = crate::compat::with_java(true, || crate::cobalt_main(&args));
    assert_eq!(status, 0, "scenario {}", s.name);
    let sample = if s.tumor_bam.is_some() { SAMPLE } else { REFERENCE };
    // the ratio file's chromosome names: chr1 for 38, 1 for 37 (keys below are without the prefix)
    let path = out.join(format!("{sample}.cobalt.ratio.tsv.gz"));
    let text = {
        use std::io::Read;
        let mut t = String::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap()).read_to_string(&mut t).unwrap();
        t
    };
    let mut rows: Rows = BTreeMap::new();
    for line in text.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        let values: Vec<f64> = f[1..].iter().map(|v| v.parse::<f64>().unwrap()).collect();
        rows.entry(f[0].trim_start_matches("chr").to_string()).or_default().push(values);
    }
    Some(Output { dir: out, rows })
}

fn close(expected: f64, actual: f64) {
    assert!((expected - actual).abs() <= 0.01, "expected {expected}, got {actual}");
}

fn one_window(name: &'static str, tumor: bool, reference: bool, targeted: bool) -> Scenario {
    Scenario {
        name,
        tumor_bam: tumor.then_some("one_window"),
        reference_bam: reference.then_some("one_window"),
        gc: vec![gc("1", 0, 3000, 0.5)],
        panel: targeted.then(|| vec![Panel { chr: "1", start: 0, stop: 3000, enrichment: 1.000001 }]),
        diploid: None,
    }
}

/// columns: position, referenceReadDepth, tumorReadDepth, referenceGCRatio, tumorGCRatio,
/// referenceGCDiploidRatio, referenceGCContent, tumorGCContent
fn check_row(row: &[f64], expected: [f64; 8]) {
    for (g, e) in row.iter().zip(expected) {
        close(e, *g);
    }
}

#[test]
fn e2e_tumor_only() {
    let Some(o) = run(one_window("tumor_only", true, false, true)) else { return };
    let r = o.chr("1");
    assert_eq!(r.len(), 3);
    check_row(&r[0], [1.0, -1.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.5]);
    check_row(&r[1], [1001.0, -1.0, 100.0, -1.0, 1.0, -1.0, -1.0, 0.5]);
    check_row(&r[2], [2001.0, -1.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.5]);
}

#[test]
fn e2e_germline_only() {
    let Some(o) = run(one_window("germline_only", false, true, true)) else { return };
    let r = o.chr("1");
    // Java's test expects referenceGCDiploidRatio 0 for the windows without reads; Java writes -1 (a
    // targeted reference's diploid ratios are set for positive ratios only), as this port does
    check_row(&r[0], [1.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.5, -1.0]);
    check_row(&r[1], [1001.0, 100.0, -1.0, 1.0, -1.0, 1.0, 0.5, -1.0]);
    check_row(&r[2], [2001.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.5, -1.0]);
    let medians = std::fs::read_to_string(o.dir.join(format!("{REFERENCE}.cobalt.ratio.median.tsv"))).unwrap();
    assert!(!o.exists(&format!("{SAMPLE}.cobalt.ratio.median.tsv")));
    // the panel's median ratios: every human chromosome, chr1's median of its positive ratios, 3 windows
    assert!(medians.lines().nth(1).unwrap().starts_with("chr1\t1.0000\t3"), "{medians}");
}

#[test]
fn e2e_tumor_and_germline() {
    let Some(o) = run(one_window("tumor_and_germline", true, true, true)) else { return };
    let r = o.chr("1");
    check_row(&r[1], [1001.0, 100.0, 100.0, 1.0, 1.0, 1.0, 0.5, 0.5]);
    // as in e2e_germline_only: -1 where Java's test expects 0
    check_row(&r[2], [2001.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.5, 0.5]);
}

#[test]
fn e2e_segmentation_by_mode() {
    let Some(o) = run(one_window("segmentation_tumor_germline", true, true, false)) else { return };
    assert!(o.exists(&format!("{SAMPLE}.cobalt.ratio.pcf")) && o.exists(&format!("{REFERENCE}.cobalt.ratio.pcf")));
    let Some(o) = run(one_window("segmentation_tumor_only", true, false, false)) else { return };
    assert!(o.exists(&format!("{SAMPLE}.cobalt.ratio.pcf")) && !o.exists(&format!("{REFERENCE}.cobalt.ratio.pcf")));
    let Some(o) = run(one_window("segmentation_germline_only", false, true, false)) else { return };
    assert!(!o.exists(&format!("{SAMPLE}.cobalt.ratio.pcf")) && o.exists(&format!("{REFERENCE}.cobalt.ratio.pcf")));
}

#[test]
fn e2e_apply_diploid_region_filtering() {
    let s = Scenario {
        name: "diploid_regions",
        tumor_bam: Some("chrs_1_and_2_10_windows_each"),
        reference_bam: None,
        gc: vec![gc("1", 0, 10_000, 0.5), gc("2", 0, 10_000, 0.5)],
        panel: Some(vec![Panel { chr: "1", start: 0, stop: 10_000, enrichment: 1.0001 }, Panel { chr: "2", start: 0, stop: 10_000, enrichment: 1.0001 }]),
        diploid: Some(vec![("1", 5000, 10000), ("2", 0, 5000)]),
    };
    let Some(o) = run(s) else { return };
    for i in 0..5 {
        close(-1.0, o.tumor_ratio("1", i));
        close(1.0, o.tumor_ratio("2", i));
    }
    for i in 5..10 {
        close(1.0, o.tumor_ratio("1", i));
        close(-1.0, o.tumor_ratio("2", i));
    }
}

#[test]
fn e2e_use_supplied_gc_if_coverage_is_low() {
    let s = Scenario {
        name: "supplied_gc_low_coverage",
        tumor_bam: Some("low_depth_regions"),
        reference_bam: None,
        gc: vec![gc("1", 1, 1001, 0.6), gc("1", 1001, 2001, 0.6), gc("1", 2001, 3001, 0.52), gc("1", 3001, 4001, 0.52), gc("1", 4001, 5001, 0.52)],
        panel: None,
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    let r = o.chr("1");
    assert_eq!(r.len(), 5);
    for (i, e) in [0.6, 0.6, 0.44, 0.46, 0.48].iter().enumerate() {
        close(*e, r[i][7]);
    }
}

fn three_windows(name: &'static str, gc_sections: Vec<Gc>, panel: Vec<Panel>) -> Scenario {
    Scenario { name, tumor_bam: Some("three_windows"), reference_bam: None, gc: gc_sections, panel: Some(panel), diploid: None }
}

#[test]
fn e2e_panel_normalisation_is_applied() {
    let panel = vec![
        Panel { chr: "1", start: 0, stop: 0, enrichment: 0.0 },
        Panel { chr: "1", start: 1000, stop: 1000, enrichment: 2.0 },
        Panel { chr: "1", start: 2000, stop: 2000, enrichment: 0.5 },
        Panel { chr: "1", start: 3000, stop: 3000, enrichment: 1.0 },
        Panel { chr: "1", start: 4000, stop: 4000, enrichment: 2.0 },
    ];
    let Some(o) = run(three_windows("panel_normalisation", vec![gc("1", 0, 4000, 0.5)], panel)) else { return };
    let seventh = 1.0 / 7.0;
    close(3.0 * seventh, o.tumor_ratio("1", 1));
    close(12.0 * seventh, o.tumor_ratio("1", 2));
    close(6.0 * seventh, o.tumor_ratio("1", 3));
}

#[test]
fn e2e_borderline_mappable_windows() {
    let m = crate::constants::MIN_MAPPABLE_PERCENTAGE;
    let sections = vec![mappable("1", 0, 1000, m + 0.01), mappable("1", 2000, 2000, m), mappable("1", 3000, 4000, m - 0.01)];
    let Some(o) = run(three_windows("borderline_mappable", sections, vec![Panel { chr: "1", start: 0, stop: 4000, enrichment: 1.0000001 }])) else { return };
    for (i, e) in [0.0, 1.0, 1.0, -1.0, -1.0].iter().enumerate() {
        close(*e, o.tumor_ratio("1", i));
    }
}

#[test]
fn e2e_filter_out_unmappable_windows() {
    let sections = vec![mappable("1", 0, 0, 0.99), mappable("1", 1000, 1000, 0.05), mappable("1", 2000, 2000, 0.99), mappable("1", 3000, 3000, 0.99), mappable("1", 4000, 4000, 0.99)];
    let Some(o) = run(three_windows("unmappable_windows", sections, vec![Panel { chr: "1", start: 0, stop: 4000, enrichment: 1.0000001 }])) else { return };
    for (i, e) in [0.0, -1.0, 1.0, 1.0, 0.0].iter().enumerate() {
        close(*e, o.tumor_ratio("1", i));
    }
}

#[test]
fn e2e_reads_are_apportioned_to_windows() {
    let s = Scenario {
        name: "reads_apportioned",
        tumor_bam: Some("three_windows_with_offset"),
        reference_bam: None,
        gc: vec![gc("1", 0, 6000, 0.5)],
        panel: Some(vec![Panel { chr: "1", start: 0, stop: 6000, enrichment: 1.0001 }]),
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    let r = o.chr("1");
    for (i, e) in [(1, 8.0), (2, 82.0), (3, 60.0), (4, 10.0)] {
        close(e, r[i][2]);
    }
}

#[test]
fn e2e_extreme_gc_filtered_targeted() {
    let s = Scenario {
        name: "extreme_gc_targeted",
        tumor_bam: Some("increasing_gc_per_window_depth_1"),
        reference_bam: None,
        gc: vec![gc("1", 0, 101_000, 0.5)],
        panel: Some(vec![Panel { chr: "1", start: 0, stop: 101_000, enrichment: 1.0001 }]),
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    let r = o.chr("1");
    for (i, row) in r.iter().enumerate().skip(1) {
        let expected = if (25..68).contains(&i) { 1.0 } else { -1.0 };
        close(expected, row[4]);
        if i < 25 {
            close(0.01 * i as f64, row[7]);
        }
    }
}

fn depth0_windows(name: &'static str, targeted: bool) -> Scenario {
    Scenario {
        name,
        tumor_bam: Some("depth_0_windows"),
        reference_bam: None,
        gc: vec![gc("1", 0, 8000, 0.5)],
        panel: targeted.then(|| vec![Panel { chr: "1", start: 0, stop: 8000, enrichment: 1.000001 }]),
        diploid: None,
    }
}

#[test]
fn e2e_depth0_windows() {
    for targeted in [true, false] {
        let Some(o) = run(depth0_windows(if targeted { "depth0_targeted" } else { "depth0_whole_genome" }, targeted)) else { return };
        let r = o.chr("1");
        assert_eq!(r[4][0], 4001.0);
        close(100.0, r[4][2]);
        close(0.52, r[4][7]);
        assert!(r[4][4] > 0.0);
        assert_eq!(r[5][0], 5001.0);
        close(0.0, r[5][4]);
        close(0.0, r[5][2]);
        close(0.5, r[5][7]);
    }
}

#[test]
fn e2e_low_coverage() {
    let s = Scenario {
        name: "low_coverage",
        tumor_bam: Some("chr1_depth_2_chr2_depth_4"),
        reference_bam: None,
        gc: vec![gc("1", 0, 100_000, 0.5), gc("2", 0, 100_000, 0.5)],
        panel: None,
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    assert_eq!((o.chr("1").len(), o.chr("2").len()), (100, 100));
    for i in 0..100 {
        let (e1, e2) = if [14, 44, 74, 95].contains(&i) { (2.0, 4.0) } else { (-1.0, -1.0) };
        close(e1, o.tumor_ratio("1", i));
        close(e2, o.tumor_ratio("2", i));
    }
}

#[test]
fn e2e_low_coverage_targeted() {
    let s = Scenario {
        name: "low_coverage_targeted",
        tumor_bam: Some("chr1_depth_2_chr2_depth_4"),
        reference_bam: None,
        gc: vec![gc("1", 0, 100_000, 0.5), gc("2", 0, 100_000, 0.5)],
        panel: Some(vec![Panel { chr: "1", start: 0, stop: 100_000, enrichment: 0.501 }, Panel { chr: "2", start: 0, stop: 100_000, enrichment: 1.001 }]),
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    for i in 0..100 {
        close(1.0, o.tumor_ratio("1", i));
        close(1.0, o.tumor_ratio("2", i));
    }
}

#[test]
fn e2e_gc_normalisation() {
    let s = Scenario {
        name: "gc_normalisation",
        tumor_bam: Some("gc_40_to_50"),
        reference_bam: None,
        gc: vec![gc("1", 1000, 61_000, 0.5)],
        panel: Some(vec![Panel { chr: "1", start: 1000, stop: 61_000, enrichment: 1.000001 }]),
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    let mean = 1.850231745;
    let check = |e: f64, idx: &[usize]| {
        for &i in idx {
            close(e, o.tumor_ratio("1", i));
        }
    };
    check(-1.0, &[0]);
    check(10.0 / 7.0 / mean, &[1, 2, 3]);
    check(20.0 / 7.0 / mean, &[4]);
    check(50.0 / 7.0 / mean, &[5]);
    check(1.0 / mean, &[6, 7, 8, 11, 12, 13, 46, 47, 48]);
    check(21.0 / 11.0 / mean, &[9]);
    check(59.0 / 19.0 / mean, &[50]);
    check(20.0 / 13.0 / mean, &[51, 52, 53]);
    check(60.0 / 13.0 / mean, &[55]);
    check(0.0, &[56, 57, 58, 59]);
}

#[test]
fn e2e_segmentation() {
    for targeted in [true, false] {
        let s = Scenario {
            name: if targeted { "segmentation_targeted" } else { "segmentation_whole_genome" },
            tumor_bam: Some("multiple_segments"),
            reference_bam: None,
            gc: vec![gc("1", 0, 100_000, 0.5), gc("2", 0, 100_000, 0.5)],
            panel: targeted.then(|| vec![Panel { chr: "1", start: 0, stop: 100_000, enrichment: 1.0001 }, Panel { chr: "2", start: 0, stop: 100_000, enrichment: 1.0001 }]),
            diploid: None,
        };
        let Some(o) = run(s) else { return };
        let pcf = std::fs::read_to_string(o.dir.join(format!("{SAMPLE}.cobalt.ratio.pcf"))).unwrap();
        let chr1: Vec<(i32, i32)> = pcf
            .lines()
            .skip(1)
            .filter(|l| l.starts_with("chr1\t"))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[1].parse().unwrap(), f[2].parse().unwrap())
            })
            .collect();
        assert_eq!(&chr1[..3], &[(1, 41000), (41001, 71000), (71001, 101000)], "{pcf}");
    }
}

fn whole_genome(name: &'static str, bam: &'static str, chrs: &[&'static str], length: i32) -> Scenario {
    Scenario { name, tumor_bam: Some(bam), reference_bam: None, gc: chrs.iter().map(|c| gc(c, 0, length, 0.5)).collect(), panel: None, diploid: None }
}

#[test]
fn e2e_whole_genome_female_and_male() {
    let Some(o) = run(whole_genome("female", "female", &["1", "2", "X"], 3000)) else { return };
    let f = 105.0 / 3.0;
    for (chr, values) in [("1", [100.0 / f, 110.0 / f, 0.0]), ("2", [110.0 / f, 100.0 / f, 0.0]), ("X", [120.0 / f, 120.0 / f, 0.0])] {
        for (i, e) in values.iter().enumerate() {
            close(*e, o.tumor_ratio(chr, i));
        }
    }
    let Some(o) = run(whole_genome("male", "male", &["1", "2", "X", "Y"], 3000)) else { return };
    let f = 104.0 / 3.0 / (104.0 / ((100.0 + 110.0 + 104.0) / 3.0));
    for (chr, values) in [("1", [100.0, 110.0, 104.0]), ("2", [110.0, 100.0, 104.0]), ("X", [60.0, 60.0, 50.0]), ("Y", [40.0, 40.0, 50.0])] {
        for (i, e) in values.iter().enumerate() {
            close(e / f, o.tumor_ratio(chr, i));
        }
    }
}

#[test]
fn e2e_whole_genome_depth0_region() {
    let s = Scenario {
        name: "depth0_region",
        tumor_bam: Some("depth_0_region_T"),
        reference_bam: Some("depth_0_region_R"),
        gc: ["1", "2", "X"].iter().map(|c| gc(c, 0, 10_000, 0.5)).collect(),
        panel: None,
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    let ref_gc_mean = 6.0 / 3.0 * (34.0 / 36.0);
    let tumor_gc_mean = 11.0 / 3.0;
    let r = o.chr("1");
    assert_eq!(r.len(), 10);
    check_row(&r[0], [1.0, 5.0, 10.0, 5.0 / ref_gc_mean, 10.0 / tumor_gc_mean, 5.0 / ref_gc_mean, 0.5, 0.5]);
    check_row(&r[3], [3001.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5]);
    check_row(&r[7], [7001.0, 6.0, 12.0, 6.0 / ref_gc_mean, 12.0 / tumor_gc_mean, 6.0 / ref_gc_mean, 0.5, 0.5]);
}

#[test]
fn e2e_pseudo_ig_regions() {
    // chr15 and chr16 at their 38 lengths; the windows that the excluded regions touch are masked
    let s = Scenario {
        name: "pseudo_ig_regions",
        tumor_bam: None,
        reference_bam: Some("pseudo_ig_regions"),
        gc: vec![gc("15", 0, 101_991_000, 0.5), gc("16", 0, 90_338_000, 0.5)],
        panel: None,
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    let r15 = o.chr("15");
    let masked = |rows: &Vec<Vec<f64>>, pos: i32| rows[(pos as usize - 1) / 1000][3] == -1.0;
    for pos in [19964001, 19972001, 19987001, 22160001, 22178001, 22184001, 22185001, 22194001, 22195001] {
        assert!(masked(r15, pos), "chr15:{pos}");
    }
    assert!(!masked(r15, 19965001) && !masked(r15, 22161001));
    let r16 = o.chr("16");
    for pos in [31962001, 32052001, 33949001, 33950001, 34013001] {
        assert!(masked(r16, pos), "chr16:{pos}");
    }
}

#[test]
fn e2e_extreme_gc_filtered_whole_genome() {
    let s = Scenario {
        name: "extreme_gc_whole_genome",
        tumor_bam: Some("increasing_gc_per_window_depth_10"),
        reference_bam: None,
        gc: vec![gc("1", 1000, 102_000, 0.5)],
        panel: None,
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    for i in 0..100 {
        let expected = if (25..68).contains(&i) { 1.0 } else { -1.0 };
        close(expected, o.tumor_ratio("1", i));
    }
}

#[test]
fn e2e_ref_genome_version_37() {
    // the 38 test BAM with ref genome version 37: chromosome names without "chr", the 37 exclusions
    let s = Scenario {
        name: "tumor_normal_v37",
        tumor_bam: Some("depth_0_region_T"),
        reference_bam: Some("depth_0_region_R"),
        gc: ["1", "2", "X"].iter().map(|c| gc(c, 0, 10_000, 0.5)).collect(),
        panel: None,
        diploid: None,
    };
    let Some(o) = run(s) else { return };
    let text = {
        use std::io::Read;
        let mut t = String::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(o.dir.join(format!("{SAMPLE}.cobalt.ratio.tsv.gz"))).unwrap()).read_to_string(&mut t).unwrap();
        t
    };
    assert!(text.lines().nth(1).unwrap().starts_with("1\t1\t"), "{}", text.lines().nth(1).unwrap());
    let pcf = std::fs::read_to_string(o.dir.join(format!("{SAMPLE}.cobalt.ratio.pcf"))).unwrap();
    assert!(pcf.lines().nth(1).unwrap().starts_with("1\t"));
    assert_eq!(o.chr("1").len(), 10);
}
