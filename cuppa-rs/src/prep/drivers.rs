//! DriverPrep (category DRIVER): drivers of the LINX (or PURPLE) catalog, indel drivers of ALB, SFTPB and
//! SLC34A2, LINX's reported fusions and the reported viruses, de-duplicated by maximum likelihood; with
//! ViralInsertionType. Viruses come from Viridian's report (the jar, `-viridian_dir`) or from Virus
//! Interpreter's annotated file (`-virus_dir`, as CUPPA 2.5.1 and the pipeline).

use anyhow::{anyhow, bail, Result};
use hmf_common::error;

use super::config::PrepConfig;
use super::items::{DataItem, Index, ItemType, Source};
use super::somatic::{self, VariantType};
use super::traits;
use crate::compat;

/// ViralInsertionType.fromVirusName
pub fn viral_insertion_type(name: &str) -> &'static str {
    if name.contains("papillo") {
        "HPV"
    } else if name.contains("HBV") || name.contains("Hepatitis B") || name.contains("hepatitis B") {
        "HBV"
    } else if name.contains("Merkel") {
        "MERKEL"
    } else if name.contains("adenovirus") || name.contains("Adeno") {
        "AAV"
    } else if name.contains("herpes") {
        "HERPES"
    } else {
        "OTHER"
    }
}

/// the items in insertion order, a later item of the same index replacing an earlier one only if its
/// value is greater as a float (DriverPrep.addDataItem, a LinkedHashMap)
struct Items {
    items: Vec<DataItem>,
}

impl Items {
    fn add(&mut self, d: DataItem) {
        let parse = |s: &str| s.parse::<f32>().unwrap_or(f32::NAN);
        match self.items.iter_mut().find(|x| x.index == d.index) {
            None => self.items.push(d),
            Some(existing) => {
                if parse(&d.value) > parse(&existing.value) {
                    *existing = d;
                }
            }
        }
    }
}

fn is_file(p: &str) -> bool {
    std::fs::metadata(p).is_ok_and(|m| m.is_file())
}

/// the java.util.String.split of a line by tabs, with a field getter that fails as Java's array access
fn field<'a>(v: &[&'a str], i: usize) -> Result<&'a str> {
    v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len()))
}

fn int(s: &str) -> Result<i32> {
    hmf_common::java_text::parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

fn dbl(s: &str) -> Result<f64> {
    hmf_common::java_text::parse_double(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

fn enum_check(v: &str, allowed: &[&str], class: &str) -> Result<()> {
    if allowed.contains(&v) { Ok(()) } else { bail!("java.lang.IllegalArgumentException: No enum constant {class}.{v}") }
}

struct Fusion {
    name: String,
    reported: bool,
    reported_type: String,
    likelihood: String,
}

/// LinxFusion.read (String.split dropping trailing empty fields; the fields Java parses checked)
fn read_fusions(path: &str) -> Result<Vec<Fusion>> {
    let lines = hmf_common::java_text::read_lines(path).map_err(|_| anyhow!("java.nio.file.NoSuchFileException: {path}"))?;
    let Some(header) = lines.first() else { bail!("java.lang.IndexOutOfBoundsException: empty file {path}") };
    let f = hmf_common::java_text::fields_index_map(header, '\t');
    let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {n} in {path}"));
    let mut out = Vec::new();
    for line in &lines[1..] {
        let v: Vec<&str> = hmf_common::java_text::split(line, '\t');
        let g = |n: &str| -> Result<&str> { field(&v, col(n)?) };
        int(g("fivePrimeBreakendId")?)?;
        int(g("threePrimeBreakendId")?)?;
        let name = g("name")?.to_string();
        let reported = hmf_common::java_text::parse_bool(g("reported")?);
        let reported_type = g("reportedType")?.to_string();
        if let Some(i) = f.get("reportableReasons") {
            field(&v, *i)?;
        }
        enum_check(g("phased")?, &["INFRAME", "SKIPPED_EXONS", "OUT_OF_FRAME"], "com.hartwig.hmftools.common.linx.FusionPhasedType")?;
        let likelihood = g("likelihood")?.to_string();
        enum_check(&likelihood, &["HIGH", "LOW", "NA"], "com.hartwig.hmftools.common.linx.FusionLikelihoodType")?;
        for opt in ["fivePrimeVcfId", "threePrimeVcfId", "fivePrimeCoords", "threePrimeCoords"] {
            if let Some(i) = f.get(opt) {
                field(&v, *i)?;
            }
        }
        int(g("chainLength")?)?;
        int(g("chainLinks")?)?;
        g("chainTerminated")?;
        g("domainsKept")?;
        g("domainsLost")?;
        for n in ["skippedExonsUp", "skippedExonsDown", "fusedExonUp", "fusedExonDown"] {
            int(g(n)?)?;
        }
        out.push(Fusion { name, reported, reported_type, likelihood });
    }
    Ok(out)
}

/// the reported viruses: (name, likelihood) from Viridian's report (likelihood 1) or from Virus
/// Interpreter's annotated file (1 for a HIGH or UNKNOWN driver likelihood, else 0.5)
fn reported_viruses(config: &PrepConfig, sample: &str) -> Result<Vec<(String, f64)>> {
    if let Some(dir) = &config.virus_dir {
        // AnnotatedVirusFile.read
        let path = config.virus_annotated(dir, sample);
        let lines = hmf_common::java_text::read_lines(&path).map_err(|_| anyhow!("java.nio.file.NoSuchFileException: {path}"))?;
        let Some(header) = lines.first() else { bail!("java.lang.IndexOutOfBoundsException: empty file {path}") };
        let f = hmf_common::java_text::fields_index_map(header, '\t');
        let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {n} in {path}"));
        let mut out = Vec::new();
        for line in &lines[1..] {
            let v: Vec<&str> = hmf_common::java_text::split_all(line, '\t');
            let g = |n: &str| -> Result<&str> { field(&v, col(n)?) };
            if let Some(i) = f.get("expectedClonalCoverage") {
                let s = field(&v, *i)?;
                if s != "null" {
                    dbl(s)?;
                }
            }
            int(g("taxid")?)?;
            let name = g("name")?.to_string();
            enum_check(
                g("qcStatus")?,
                &["NO_ABNORMALITIES", "LOW_VIRAL_COVERAGE", "EXCESSIVE_VIRAL_COVERAGE", "ASSEMBLY_DOWNSAMPLED", "CHILD_TAXID_REFERENCE", "UNCLEAR_TAXID_ASSIGNMENT"],
                "com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus",
            )?;
            int(g("integrations")?)?;
            let interp = g("interpretation")?;
            if interp != "null" && !["MCV", "EBV", "HPV", "HBV", "HHV-8"].contains(&interp) {
                bail!("java.lang.IllegalStateException: Cannot resolve virus name: {interp}");
            }
            for opt in ["percentageCovered", "meanCoverage"] {
                if let Some(i) = f.get(opt) {
                    dbl(field(&v, *i)?)?;
                }
            }
            let reported = hmf_common::java_text::parse_bool(g("reported")?);
            if let Some(i) = f.get("blacklisted") {
                field(&v, *i)?;
            }
            let likelihood = match f.get("driverLikelihood") {
                Some(i) => {
                    let s = field(&v, *i)?;
                    enum_check(s, &["HIGH", "LOW", "UNKNOWN"], "com.hartwig.hmftools.common.virusinterpreter.VirusLikelihoodType")?;
                    s.to_string()
                }
                None => "UNKNOWN".to_string(),
            };
            if reported {
                out.push((name, if likelihood == "HIGH" || likelihood == "UNKNOWN" { 1.0 } else { 0.5 }));
            }
        }
        return Ok(out);
    }
    // VirusReportFile.read (DelimFileReader: "null" / "NULL" are nulls)
    let path = config.viridian_report(sample);
    let lines = hmf_common::java_text::read_lines(&path).map_err(|e| anyhow!("java.io.UncheckedIOException: java.nio.file.NoSuchFileException: {path} ({e})"))?;
    let Some(header) = lines.first() else { return Ok(Vec::new()) };
    let cols: Vec<&str> = hmf_common::java_text::split_all(header, '\t');
    let col = |n: &str| cols.iter().position(|c| *c == n).ok_or_else(|| anyhow!("java.util.NoSuchElementException: column: {n} not found"));
    let mut out = Vec::new();
    for line in &lines[1..] {
        let v: Vec<&str> = hmf_common::java_text::split_all(line, '\t');
        let raw = |n: &str| -> Result<Option<&str>> {
            let s = field(&v, col(n)?)?;
            Ok(if s == "null" || s == "NULL" { None } else { Some(s) })
        };
        let req = |n: &str| -> Result<&str> { raw(n)?.ok_or_else(|| anyhow!("java.util.NoSuchElementException: column {} is null", n)) };
        let name = req("name")?.to_string();
        req("virusGroup")?;
        enum_check(req("qcStatus")?, &["NO_ABNORMALITIES", "AMBIGUOUS_VIRUS"], "com.hartwig.hmftools.common.viridian.ReportQcStatus")?;
        if let Some(t) = raw("interpretation")? {
            if !["MCV", "EBV", "HPV", "HBV", "HHV-8"].contains(&t) {
                bail!("java.lang.IllegalStateException: Cannot resolve virus name: {t}");
            }
        }
        int(req("integrations")?)?;
        dbl(req("percentageCovered")?)?;
        dbl(req("meanCoverage")?)?;
        if let Some(c) = raw("copiesPerTumorCell")? {
            dbl(c)?;
        }
        enum_check(
            req("clonalityStatus")?,
            &["CLONAL", "NOT_CLONAL", "NOT_PROFILED", "PURPLE_FAIL_NO_TUMOR", "PURPLE_FAIL_CONTAMINATION"],
            "com.hartwig.hmftools.common.viridian.ClonalityStatus",
        )?;
        let reported = match req("reported")?.to_lowercase().as_str() {
            "true" => true,
            "false" => false,
            s => bail!("java.lang.IllegalArgumentException: Invalid boolean string: {s}"),
        };
        if reported {
            out.push((name, 1.0));
        }
    }
    Ok(out)
}

/// DriverType.isGermline; Java leaves out GERMLINE_AMP (CU6)
fn is_germline(driver: &str) -> bool {
    matches!(driver, "GERMLINE_DELETION" | "GERMLINE_DISRUPTION" | "GERMLINE_MUTATION") || (!compat::java("CU6") && driver == "GERMLINE_AMP")
}

pub fn extract(config: &PrepConfig, sample: &str) -> Result<Vec<DataItem>> {
    let mut items = Items { items: Vec::new() };
    let item = |t: ItemType, key: String, v: f64| DataItem {
        index: Index { source: Source::Dna, item_type: t, key },
        value: hmf_common::fmt::java_fixed(v, 4),
    };

    // the driver catalog: LINX's, else PURPLE's
    let linx = config.linx_driver_catalog(sample);
    let catalog = if is_file(&linx) {
        linx
    } else {
        let purple = config.purple_driver_catalog(sample);
        if !is_file(&purple) {
            error!("sample({sample}) has no LINX or PURPLE driver catalog file");
            bail!("java.nio.file.NoSuchFileException: {purple}");
        }
        purple
    };
    for d in hmf_common::purple_files::read_driver_catalog(&catalog)? {
        if is_germline(&d.driver) || d.reported_status != "REPORTED" {
            continue;
        }
        let suffix = if d.driver == "AMP" || d.driver == "PARTIAL_AMP" { ".amp" } else { ".mut" };
        items.add(item(ItemType::Driver, format!("{}{suffix}", d.gene), d.driver_likelihood));
    }

    // indel drivers of microsatellite-stable samples
    let purity = traits::read_purity(&config.purple_purity(sample))?;
    if purity.ms_status == "MSS" {
        for v in somatic::load_variants(config, sample, None)? {
            let known = v.gene == "ALB" || v.gene == "SFTPB" || v.gene == "SLC34A2";
            if known && v.variant_type == VariantType::Indel && v.repeat_count <= 6 {
                items.add(item(ItemType::Driver, format!("{}.indel", v.gene), 1.0));
            }
        }

        // fusions
        fusions(config, sample, &mut items, &item)?;
        viruses(config, sample, &mut items, &item)?;
    } else {
        fusions(config, sample, &mut items, &item)?;
        viruses(config, sample, &mut items, &item)?;
    }
    Ok(items.items)
}

fn fusions(config: &PrepConfig, sample: &str, items: &mut Items, item: &dyn Fn(ItemType, String, f64) -> DataItem) -> Result<()> {
    for f in read_fusions(&config.linx_fusions(sample))? {
        if !f.reported {
            continue;
        }
        let p5 = f.reported_type == "PROMISCUOUS_5";
        let p3 = f.reported_type == "PROMISCUOUS_3";
        let name = if p5 || p3 {
            if f.likelihood != "HIGH" {
                continue;
            }
            let genes: Vec<&str> = hmf_common::java_text::split(&f.name, '_');
            if p5 {
                format!("{}_PROM5", field(&genes, 0)?)
            } else {
                format!("{}_PROM3", field(&genes, 1)?)
            }
        } else {
            f.name.clone()
        };
        items.add(item(ItemType::Fusion, name, 1.0));
    }
    Ok(())
}

fn viruses(config: &PrepConfig, sample: &str, items: &mut Items, item: &dyn Fn(ItemType, String, f64) -> DataItem) -> Result<()> {
    for (name, likelihood) in reported_viruses(config, sample)? {
        items.add(item(ItemType::Virus, viral_insertion_type(&name).to_string(), likelihood));
    }
    Ok(())
}
