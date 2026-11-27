//! The Circos input files of circos/ (CircosCharts, CircosFileWriter, CircosSnvWriter, CircosIndelWriter,
//! CircosLinkWriter, CircosDriverWriter) and the Circos runs (CircosExecution). The R plots (RCharts) are
//! not ported: REPORT.md lists them.

use anyhow::{bail, Result};

use hmf_common::fmt::{java_double_to_string, java_fixed};
use hmf_common::java_id_hash::{IdHashMap, IdKey};
use hmf_common::java_strictmath;

use crate::config::PurpleConfig;
use crate::copynumber::PurpleCopyNumber;
use crate::drivers::DriverSource;
use crate::fitting::downsample;
use crate::inputs::{AmberBaf, AmberData};
use crate::region::ObservedRegion;
use crate::sv::StructuralVariant;
use crate::types::*;
use crate::variant::Variant;

const MAX_PLOT_POINTS: usize = 25000;
const CIRCOS_TEMPLATE: &str = include_str!("../resources/circos.template");
const INPUT_TEMPLATE: &str = include_str!("../resources/input.template");
const GAPS_37: &str = include_str!("../resources/centromere_gaps.37.tsv");
const GAPS_38: &str = include_str!("../resources/centromere_gaps.38.tsv");

/// CircosFileWriter.circosContig
fn contig(chromosome: &str) -> String {
    format!("hs{}", chr(chromosome).short_name())
}

/// Files.write(path, lines): every line ends in a newline
fn write_lines(path: &str, lines: &[String]) -> Result<()> {
    let mut text = String::with_capacity(lines.iter().map(|l| l.len() + 1).sum());
    for l in lines {
        text.push_str(l);
        text.push('\n');
    }
    std::fs::write(path, text).map_err(|e| anyhow::anyhow!("{path}: {e}"))
}

const HEADER: &str = "#chromosome\tstart\tend\tvalue";

fn region_line(chromosome: &str, start: i32, end: i32, value: f64, glyph: Option<f64>) -> String {
    let mut s = format!("{}\t{}\t{}\t{}", contig(chromosome), start, end, java_fixed(value, 4));
    if let Some(g) = glyph {
        s.push_str(&format!("\tglyph_size={}", java_fixed(g, 2)));
    }
    s
}

fn swap_base(b: &str) -> &str {
    match b {
        "A" => "T",
        "T" => "A",
        "C" => "G",
        "G" => "C",
        other => other,
    }
}

/// CircosSnvWriter.color
fn snv_colour(v: &Variant) -> &'static str {
    let r = v.ref_bases();
    let a = v.alt();
    let sig = |sr: &str, sa: &str| (r == sr && a == sa) || (r == swap_base(sr) && a == swap_base(sa));
    if sig("C", "A") {
        "20,176,239"
    } else if sig("C", "G") {
        "6,8,9"
    } else if sig("C", "T") {
        "224,7,20"
    } else if sig("T", "A") {
        "191,190,191"
    } else if sig("T", "C") {
        "144,202,75"
    } else if sig("T", "G") {
        "233,187,184"
    } else {
        "purple"
    }
}

/// CircosIndelWriter.color
fn indel_colour(v: &Variant) -> &'static str {
    if v.alt().len() > v.ref_bases().len() {
        "green"
    } else if v.rec.info.as_int("REP_C", 0) >= 4 {
        "orange"
    } else if v.rec.info.as_string("MH", "").len() > 2 {
        "vdred"
    } else {
        "vdyellow"
    }
}

fn sv_colour(t: SvType) -> &'static str {
    match t {
        SvType::Dup => "green",
        SvType::Del => "red",
        SvType::Bnd => "blue",
        SvType::Ins => "vdyellow",
        SvType::Inv => "black",
        _ => "purple",
    }
}

/// CircosLinkWriter.bezierRadius
fn bezier_radius(v: &StructuralVariant) -> Option<String> {
    let (s, e) = (v.chromosome(true)?, v.chromosome(false)?);
    if s != e {
        return None;
    }
    let length = (v.position(false)? - v.position(true)? + 1).max(1);
    let r = jmin(java_strictmath::log10(length as f64), 7.0);
    let r = 0.35 - r / 7.0 * (0.35 - 0.1);
    Some(format!("{}r", java_fixed(r, 3)))
}

/// CircosDriverWriter.filterAndPrioritise: reported somatic drivers of medium or high likelihood, one per
/// gene, copy-number events first
fn prioritise_drivers(sources: &[DriverSource]) -> Vec<&DriverSource> {
    let mut list: Vec<&DriverSource> = sources
        .iter()
        .filter(|s| s.driver.reported_status == ReportedStatus::Reported)
        .filter(|s| DriverInterpretation::interpret(s.driver.driver_likelihood) != DriverInterpretation::Low)
        .filter(|s| matches!(s.driver.driver, DriverType::Amp | DriverType::PartialAmp | DriverType::Del | DriverType::HetDel | DriverType::Loh | DriverType::Mutation))
        .collect();
    let rank = |t: DriverType| match t {
        DriverType::Amp | DriverType::PartialAmp | DriverType::Del => 0,
        DriverType::Loh | DriverType::HetDel => 1,
        DriverType::Mutation => 2,
        _ => 3,
    };
    // Collections.sort: stable
    list.sort_by(|a, b| {
        let (l1, l2) = (a.driver.driver_likelihood, b.driver.driver_likelihood);
        if l1 != l2 {
            return if l1 > l2 { std::cmp::Ordering::Less } else { std::cmp::Ordering::Greater };
        }
        rank(a.driver.driver).cmp(&rank(b.driver.driver))
    });
    let mut seen = std::collections::HashSet::new();
    list.retain(|s| seen.insert(s.driver.gene.clone()));
    list
}

fn driver_text(s: &DriverSource) -> String {
    match s.gene_end {
        Some(end) => {
            let mid = ((s.position as i64 + end as i64) / 2) as i32;
            let colour = if s.driver.driver == DriverType::Del { "vdred" } else { "vdgreen" };
            format!("{}\t{mid}\t{mid}\t{}\tcolor={colour}", contig(&s.chromosome), s.gene)
        }
        None => {
            let colour = if s.likely_germline { "grey" } else { "vdblue" };
            format!("{}\t{}\t{}\t{}\tcolor={colour}", contig(&s.chromosome), s.position, s.position, s.gene)
        }
    }
}

fn driver_pointer(s: &DriverSource) -> String {
    match s.gene_end {
        Some(end) => {
            let r0 = match s.driver.driver {
                DriverType::Amp | DriverType::PartialAmp => 0.7,
                DriverType::Del => 0.55,
                _ => 0.5,
            };
            let mid = ((s.position as i64 + end as i64) / 2) as i32;
            format!("{}\t{mid}\t{mid}\t{}\tr0={}r", contig(&s.chromosome), s.gene, java_fixed(r0, 2))
        }
        None => {
            let r0 = if s.is_indel { 0.775 } else { 0.775 + s.adjusted_vaf.clamp(0.0, 1.0) * (0.975 - 0.775) };
            format!("{}\t{}\t{}\t{}\tr0={}r", contig(&s.chromosome), s.position, s.position, s.gene, java_fixed(r0, 2))
        }
    }
}

/// a HumanChromosome key of the BAF multimap (an enum: identity hash)
#[derive(Clone)]
struct ChromosomeKey(Chromosome);

impl IdKey for ChromosomeKey {
    fn id_addr(&self) -> usize {
        self.0.ordinal()
    }
}

/// `Lists.newArrayList(amberData.ChromosomeBafs.values())`: the BAF points chromosome by chromosome.
/// PU2: Java's order is that of an ArrayListMultimap keyed by HumanChromosome, whose hash is the identity
/// hash, so it changes from run to run. Corrected: genome order. Java mode: the order of a JVM run with
/// -XX:hashCode=2 (every identity hash equal), the defined Java reference: ArrayListMultimap.create()
/// is a HashMap with a table of 16, filled in the order the chromosomes first appear.
fn bafs_in_java_order(amber: &AmberData) -> Vec<&AmberBaf> {
    let mut out = Vec::new();
    if crate::compat::java("PU2") {
        let mut map: IdHashMap<ChromosomeKey, ()> = IdHashMap::with_capacity(16);
        for c in &amber.chromosome_order {
            map.insert(ChromosomeKey(*c), ());
        }
        for k in map.keys() {
            out.extend(amber.chromosome_bafs[&k.0].iter());
        }
    } else {
        for bafs in amber.chromosome_bafs.values() {
            out.extend(bafs.iter());
        }
    }
    out
}

pub struct CircosInputs<'a> {
    pub reference_id: Option<&'a str>,
    pub tumor_id: &'a str,
    pub gender: Gender,
    pub v38: bool,
    pub copy_numbers: &'a [PurpleCopyNumber],
    /// SomaticStream.plottingVariants
    pub somatics: Vec<&'a Variant>,
    pub svs: &'a [StructuralVariant],
    /// the fitted regions other than EXCLUDED
    pub regions: Vec<&'a ObservedRegion>,
    pub amber: &'a AmberData,
    pub immune_regions: &'a [ChrRegion],
    pub drivers: &'a [DriverSource],
}

/// CircosCharts.write
pub fn write_circos(config: &PurpleConfig, d: &CircosInputs) -> Result<()> {
    let dir = &config.charting.circos_directory;
    let base = format!("{dir}{}", d.tumor_id);
    let reference_base = d.reference_id.map(|r| format!("{dir}{r}"));

    // writeConfig
    for (kind, template) in [("circos", CIRCOS_TEMPLATE), ("input", INPUT_TEMPLATE)] {
        let content = template
            .replace("SAMPLE", d.tumor_id)
            .replace("REFERENCE", d.reference_id.unwrap_or(d.tumor_id))
            .replace("COBALT_REF_COLOUR", if d.reference_id.is_some() { "green" } else { "cobalt" })
            .replace("EXCLUDE", if d.gender == Gender::Female { "hsY" } else { "hsZ" })
            .replace("KARYOTYPE", if d.v38 { "data/karyotype/karyotype.human.hg38.txt" } else { "data/karyotype/karyotype.human.hg19.txt" });
        std::fs::write(format!("{base}.{kind}.conf"), content)?;
    }
    let gaps: Vec<String> = (if d.v38 { GAPS_38 } else { GAPS_37 }).lines().map(|l| l.to_string()).collect();
    write_lines(&format!("{dir}gaps.txt"), &gaps)?;

    // writeCopyNumbers (the BAF file is written again below from the AMBER points)
    let with_header = |lines: Vec<String>| std::iter::once(HEADER.to_string()).chain(lines).collect::<Vec<_>>();
    let cn = d.copy_numbers;
    write_lines(
        &format!("{base}.map.circos"),
        &with_header(cn.iter().filter(|c| !c.minor_allele_copy_number().is_nan()).map(|c| region_line(&c.chromosome, c.start, c.end, c.minor_allele_copy_number() - 1.0, None)).collect()),
    )?;
    write_lines(
        &format!("{base}.cnv.circos"),
        &with_header(cn.iter().filter(|c| !c.average_tumor_copy_number.is_nan()).map(|c| region_line(&c.chromosome, c.start, c.end, c.average_tumor_copy_number - 2.0, None)).collect()),
    )?;

    // writeEnrichedSomatics
    let somatics: Vec<&Variant> = d.somatics.iter().copied().filter(|v| Chromosome::contains(v.chromosome())).collect();
    let snps: Vec<&Variant> = somatics.iter().copied().filter(|v| v.vtype == VariantType::Snp).collect();
    let indels: Vec<&Variant> = somatics.iter().copied().filter(|v| v.vtype == VariantType::Indel).collect();
    let snp_lines: Vec<String> = downsample(MAX_PLOT_POINTS, &snps)
        .iter()
        .map(|v| {
            format!(
                "{}\t{}\t{}\t{}\tcolor={},glyph=circle,glyph_size=10",
                contig(v.chromosome()),
                v.position(),
                v.position(),
                java_double_to_string(v.rec.info.as_double("PURPLE_AF", 0.0)),
                snv_colour(v)
            )
        })
        .collect();
    write_lines(&format!("{base}.snp.circos"), &with_header(snp_lines))?;
    let indel_lines: Vec<String> =
        downsample(MAX_PLOT_POINTS, &indels).iter().map(|v| format!("{}\t{}\t{}\t1\tcolor={}", contig(v.chromosome()), v.position(), v.position(), indel_colour(v))).collect();
    write_lines(&format!("{base}.indel.circos"), &with_header(indel_lines))?;

    // writeStructuralVariants (no header)
    let links: Vec<String> = d
        .svs
        .iter()
        .filter(|v| v.end.is_some())
        .map(|v| {
            let mut s = format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                contig(v.chromosome(true).unwrap()),
                v.position(true).unwrap(),
                v.position(true).unwrap(),
                contig(v.chromosome(false).unwrap()),
                v.position(false).unwrap(),
                v.position(false).unwrap()
            );
            match bezier_radius(v) {
                None => s.push_str(&format!("\tcolor={}", sv_colour(v.sv_type))),
                Some(r) => s.push_str(&format!("\tcolor={},bezier_radius={r}", sv_colour(v.sv_type))),
            }
            s
        })
        .collect();
    write_lines(&format!("{base}.link.circos"), &links)?;

    // writeObservedRegions: the diploid regions of the downsampled fitted regions
    let regions: Vec<&ObservedRegion> = downsample(MAX_PLOT_POINTS, &d.regions).into_iter().filter(|r| r.germline_status == GermlineStatus::Diploid).collect();
    let glyph = |r: &ObservedRegion| 4.0 * (java_strictmath::log10(r.depth_window_count.max(1) as f64) + 1.0);
    if let Some(rb) = &reference_base {
        let lines = regions
            .iter()
            .filter(|r| positive(r.unnormalised_observed_normal_ratio))
            .map(|r| region_line(&r.chromosome, r.start, r.end, r.unnormalised_observed_normal_ratio, Some(glyph(r))))
            .collect();
        write_lines(&format!("{rb}.ratio.circos"), &with_header(lines))?;
    }
    let lines =
        regions.iter().filter(|r| positive(r.observed_tumor_ratio)).map(|r| region_line(&r.chromosome, r.start, r.end, r.observed_tumor_ratio, Some(glyph(r)))).collect();
    write_lines(&format!("{base}.ratio.circos"), &with_header(lines))?;

    // writeBafs
    let bafs = bafs_in_java_order(d.amber);
    let lines = downsample(MAX_PLOT_POINTS, &bafs)
        .into_iter()
        .filter(|b| !d.immune_regions.iter().any(|r| r.contains_position(&b.chromosome, b.position)))
        .map(|b| format!("{}\t{}\t{}\t{}", contig(&b.chromosome), b.position, b.position, java_double_to_string(b.tumor_baf)))
        .collect();
    write_lines(&format!("{base}.baf.circos"), &with_header(lines))?;

    // writeDrivers
    let drivers = prioritise_drivers(d.drivers);
    write_lines(&format!("{base}.driver_text.circos"), &drivers.iter().map(|s| driver_text(s)).collect::<Vec<_>>())?;
    write_lines(&format!("{base}.driver_pointer.circos"), &drivers.iter().map(|s| driver_pointer(s)).collect::<Vec<_>>())?;
    Ok(())
}

/// CircosExecution.generateCircos for the "input" and "circos" plots
pub fn run_circos(config: &PurpleConfig, tumor_id: &str) -> Result<()> {
    let Some(binary) = config.charting.circos_binary.as_deref() else { return Ok(()) };
    let plot_dir = &config.charting.plot_directory;
    let abs = |p: &str| std::path::absolute(p).map(|p| p.to_string_lossy().into_owned()).unwrap_or_else(|_| p.to_string());
    for kind in ["input", "circos"] {
        let output_file = format!("{tumor_id}.{kind}.png");
        let plot_path = format!("{plot_dir}/{output_file}");
        let _ = std::fs::remove_file(&plot_path);
        let conf = format!("{}{tumor_id}.{kind}.conf", config.charting.circos_directory);
        // circos writes a temporary circos.colorlist: each run gets its own TMPDIR
        let tmp = std::env::temp_dir().join(format!("circos{}_{kind}", std::process::id()));
        std::fs::create_dir_all(&tmp)?;
        let out = std::process::Command::new(binary)
            .args(["-nosvg", "-conf", &abs(&conf), "-outputdir", &abs(plot_dir), "-outputfile", &output_file])
            .env("TMPDIR", &tmp)
            .output();
        let _ = std::fs::remove_dir_all(&tmp);
        let out = out.map_err(|e| anyhow::anyhow!("failed to run circos({binary}): {e}"))?;
        if !out.status.success() {
            eprint!("{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
            bail!("failed to create circos plot");
        }
        if !std::path::Path::new(&plot_path).exists() {
            eprint!("{}", String::from_utf8_lossy(&out.stdout));
            bail!("failed to create file {plot_path}");
        }
    }
    Ok(())
}
