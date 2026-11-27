//! StructuralVariantPrep (category SV): PURPLE's SV VCF (StructuralVariantFileLoader with AlwaysPassFilter,
//! EnrichedStructuralVariantFactory, convertSvData) and LINX's clusters.

use anyhow::{anyhow, bail, Result};
use hmf_common::sv::{always_pass, read_vcf, StructuralVariantData, StructuralVariantFactory, StructuralVariantType, PASS_FILTER};

use super::config::PrepConfig;
use super::items::{DataItem, ItemType, Source};

/// SvDataType, in declaration order
const SV_DATA_TYPES: [&str; 6] = ["LINE", "SIMPLE_DEL_20KB_1MB", "SIMPLE_DUP_32B_200B", "SIMPLE_DUP_100KB_5MB", "MAX_COMPLEX_SIZE", "TELOMERIC_SGL"];

pub struct LinxCluster {
    pub resolved_type: String,
    pub cluster_count: i32,
}

/// LinxCluster.read (String.split dropping trailing empty fields; the fields Java parses checked)
pub fn read_clusters(path: &str) -> Result<Vec<LinxCluster>> {
    let lines = hmf_common::java_text::read_lines(path).map_err(|_| anyhow!("java.nio.file.NoSuchFileException: {path}"))?;
    let Some(header) = lines.first() else { bail!("java.lang.IndexOutOfBoundsException: empty file {path}") };
    let f = hmf_common::java_text::fields_index_map(header, '\t');
    let col = |n: &str| f.get(n).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {n} in {path}"));
    let (id, cat, syn, rt, cc, desc) = (col("clusterId")?, col("category")?, col("synthetic")?, col("resolvedType")?, col("clusterCount")?, col("clusterDesc")?);
    let mut out = Vec::new();
    for line in &lines[1..] {
        let v: Vec<&str> = hmf_common::java_text::split(line, '\t');
        let get = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len()));
        let int = |s: &str| hmf_common::java_text::parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""));
        int(get(id)?)?;
        get(cat)?;
        get(syn)?;
        get(desc)?;
        out.push(LinxCluster { resolved_type: get(rt)?.to_string(), cluster_count: int(get(cc)?)? });
    }
    Ok(out)
}

/// StructuralVariantPrep.extractSvCounts
pub fn sv_counts(svs: &[StructuralVariantData], clusters: &[LinxCluster]) -> [i32; 6] {
    let line_count: i32 = clusters.iter().filter(|c| c.resolved_type == "LINE").map(|c| c.cluster_count).fold(0i32, |a, b| a.wrapping_add(b));
    let pass: Vec<&StructuralVariantData> = svs.iter().filter(|s| s.filter.is_empty() || s.filter == PASS_FILTER).collect();
    let telomeric = pass
        .iter()
        .filter(|s| s.sv_type() == StructuralVariantType::SGL)
        .filter(|s| s.insert_sequence_repeat_class == "Simple_repeat")
        .filter(|s| s.insert_sequence_repeat_type == "(TTAGGG)n" || s.insert_sequence_repeat_type == "(CCCTAA)n")
        .count() as i32;
    let lengths = |t: StructuralVariantType| -> Vec<i32> {
        pass.iter().filter(|s| s.sv_type() == t).map(|s| s.end_position.wrapping_sub(s.start_position)).collect()
    };
    let dels = lengths(StructuralVariantType::DEL);
    let dups = lengths(StructuralVariantType::DUP);
    let short_dels = dels.iter().filter(|x| **x as f64 >= 2e4 && **x as f64 <= 1e6).count() as i32;
    let short_dups = dups.iter().filter(|x| **x >= 32 && **x <= 200).count() as i32;
    let long_dups = dups.iter().filter(|x| **x as f64 >= 1e5 && **x as f64 <= 5e6).count() as i32;
    let max_event = clusters.iter().filter(|c| c.resolved_type != "LINE").map(|c| c.cluster_count).max().unwrap_or(0);
    [line_count, short_dels, short_dups, long_dups, max_event, telomeric]
}

pub fn extract(config: &PrepConfig, sample: &str) -> Result<Vec<DataItem>> {
    let sv_file = config.purple_sv_vcf(sample);
    if !std::fs::metadata(&sv_file).is_ok_and(|m| m.is_file()) {
        bail!("java.nio.file.NoSuchFileException: {sv_file}");
    }
    let vcf = read_vcf(&sv_file)?;
    let mut factory = StructuralVariantFactory::build(always_pass);
    factory.set_enrich(true);
    for r in vcf.records {
        factory.add_variant_context(r)?;
    }
    let svs: Vec<StructuralVariantData> = factory.results().iter().enumerate().map(|(i, v)| StructuralVariantData::convert(v, i as i32)).collect();
    let cluster_file = config.linx_clusters(sample);
    if !std::fs::metadata(&cluster_file).is_ok_and(|m| m.is_file()) {
        bail!("java.nio.file.NoSuchFileException: {cluster_file}");
    }
    let clusters = read_clusters(&cluster_file)?;
    let counts = sv_counts(&svs, &clusters);
    Ok(SV_DATA_TYPES.iter().zip(counts).map(|(t, c)| DataItem::int(Source::Dna, ItemType::SvCount, t, c as i64)).collect())
}
