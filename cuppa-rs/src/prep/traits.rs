//! SampleTraitPrep (category SAMPLE_TRAIT) and PurplePurity.read.

use anyhow::{anyhow, bail, Result};

use super::config::PrepConfig;
use super::items::{DataItem, ItemType, Source};

/// The fields of PurplePurity that CUPPA uses; read checks every field as PurplePurity.read parses it
pub struct PurplePurity {
    pub male: bool,
    pub whole_genome_duplication: bool,
    pub ms_indels_per_mb: f64,
    /// MicrosatelliteStatus: MSI, MSS or UNKNOWN
    pub ms_status: String,
}

fn enum_value(v: &str, allowed: &[&str], class: &str) -> Result<String> {
    if allowed.contains(&v) {
        Ok(v.to_string())
    } else {
        bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.purple.{class}.{v}")
    }
}

pub fn read_purity(path: &str) -> Result<PurplePurity> {
    let lines = hmf_common::java_text::read_lines(path).map_err(|_| anyhow!("java.nio.file.NoSuchFileException: {path}"))?;
    if lines.len() < 2 {
        bail!("java.lang.IndexOutOfBoundsException: Index 1 out of bounds for length {}", lines.len());
    }
    let f = hmf_common::java_text::fields_index_map(&lines[0], '\t');
    let v: Vec<&str> = hmf_common::java_text::split_all(&lines[1], '\t');
    let get = |n: &str| -> Result<&str> {
        let i = f.get(n).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {n} in {path}"))?;
        v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len()))
    };
    let dbl = |n: &str| -> Result<f64> {
        let s = get(n)?;
        hmf_common::java_text::parse_double(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
    };
    let int = |n: &str| -> Result<i32> {
        let s = get(n)?;
        hmf_common::java_text::parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
    };
    for n in ["purity", "normFactor", "score", "diploidProportion", "ploidy"] {
        dbl(n)?;
    }
    let gender = enum_value(get("gender")?, &["MALE", "FEMALE"], "Gender")?;
    enum_value(get("status")?, &["NORMAL", "HIGHLY_DIPLOID", "SOMATIC", "NO_TUMOR"], "FittedPurityMethod")?;
    for n in ["polyclonalProportion", "minPurity", "maxPurity", "minPloidy", "maxPloidy", "minDiploidProportion", "maxDiploidProportion", "somaticPenalty"] {
        dbl(n)?;
    }
    let wgd = hmf_common::java_text::parse_bool(get("wholeGenomeDuplication")?);
    let ms_indels = dbl("msIndelsPerMb")?;
    let ms_status = enum_value(get("msStatus")?, &["MSI", "MSS", "UNKNOWN"], "MicrosatelliteStatus")?;
    int("tml")?;
    enum_value(get("tmlStatus")?, &["HIGH", "LOW", "UNKNOWN"], "TumorMutationalStatus")?;
    dbl("tmbPerMb")?;
    enum_value(get("tmbStatus")?, &["HIGH", "LOW", "UNKNOWN"], "TumorMutationalStatus")?;
    int("svTumorMutationalBurden")?;
    enum_value(get("runMode")?, &["TUMOR_GERMLINE", "TUMOR", "GERMLINE"], "RunMode")?;
    get("targeted")?;
    Ok(PurplePurity { male: gender == "MALE", whole_genome_duplication: wgd, ms_indels_per_mb: ms_indels, ms_status })
}

pub fn extract(config: &PrepConfig, sample: &str) -> Result<Vec<DataItem>> {
    let p = read_purity(&config.purple_purity(sample))?;
    Ok(vec![
        DataItem::bool(Source::Dna, ItemType::SampleTrait, "is_male", p.male),
        DataItem::fixed(Source::Dna, ItemType::TumorMutationalBurden, "indels_per_mb", p.ms_indels_per_mb, 4),
        DataItem::bool(Source::Dna, ItemType::SampleTrait, "whole_genome_duplication", p.whole_genome_duplication),
    ])
}
