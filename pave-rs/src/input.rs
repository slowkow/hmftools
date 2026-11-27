//! Reading the input VCF: its header, and the records each chromosome task processes.
//!
//! Java gives each chromosome task its own reader, which queries the chromosome through the index. The
//! port reads the file once: a record is kept for the task of its chromosome if it passes the filter
//! type and the specific regions, as the task's checks would keep it. The positions of the kept
//! variants are collected for the annotation resources.

use anyhow::{Context, Result};
use std::collections::HashMap;

use hmf_common::vcf::{self, VcfHeader};

use crate::annotation::NeededPositions;
use crate::config::{FilterType, PaveConfig};

pub struct InputVariants {
    pub header: VcfHeader,
    /// the record lines of each chromosome the run processes, by versioned chromosome name
    pub by_chr: HashMap<String, std::sync::Arc<Vec<String>>>,
    /// the records read per chromosome, kept or not (the task's log count)
    pub read_counts: HashMap<String, usize>,
    pub needed: NeededPositions,
}

/// The value of an INFO key in the INFO column text, without parsing the rest
fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|f| f.strip_prefix(key).and_then(|r| r.strip_prefix('=')))
}

pub fn read_header(path: &str) -> Result<VcfHeader> {
    let (header, _, _) = vcf::read_header(path).with_context(|| format!("invalid somatic VCF file({path})"))?;
    Ok(header)
}

/// Reads the records of the chromosomes in `chromosomes` (versioned names)
pub fn read_variants(config: &PaveConfig, chromosomes: &[String]) -> Result<InputVariants> {
    let (header, mut reader, first) = vcf::read_header(&config.vcf_file).with_context(|| format!("invalid somatic VCF file({})", config.vcf_file))?;
    let mut by_chr: HashMap<String, Vec<String>> = HashMap::new();
    // PA3: a VCF read without its index (not ending in .vcf.gz) is read whole by every chromosome task in
    // Java, which processes each human-chromosome record in every task; the writer then fails on the
    // records' order. Corrected: each task takes its chromosome's records, as through the index.
    let every_task_reads_all = !config.require_index() && crate::compat::java("PA3");
    let mut all: Vec<String> = Vec::new();
    let mut read_counts: HashMap<String, usize> = HashMap::new();
    let mut needed = NeededPositions::default();
    let regions = &config.specific_regions.regions;

    let mut pending = first;
    let mut current_chr = String::new();
    let mut current_included = false;
    loop {
        let line: String = match pending.take() {
            Some(l) => l,
            None => match reader.next_line()? {
                Some(l) => String::from_utf8_lossy(l).into_owned(),
                None => break,
            },
        };
        if line.is_empty() {
            continue;
        }
        let mut cols = line.splitn(9, '\t');
        let chrom = cols.next().unwrap_or("");
        if chrom != current_chr {
            current_chr = chrom.to_string();
            current_included = if every_task_reads_all {
                hmf_common::regions::is_human_chromosome(chrom)
            } else {
                chromosomes.iter().any(|c| c == chrom)
            };
        }
        if !current_included {
            continue;
        }
        let (Some(pos_s), Some(_id), Some(ref_), Some(alt), Some(_qual), Some(filter), Some(info)) =
            (cols.next(), cols.next(), cols.next(), cols.next(), cols.next(), cols.next(), cols.next())
        else {
            anyhow::bail!("invalid VCF record: {line}");
        };
        let Some(position) = hmf_common::java_text::parse_int(pos_s) else { anyhow::bail!("invalid VCF position: {line}") };

        // the task's region filter: a variant's start within one of the specific regions
        if config.specific_regions.has_filters() && !regions.is_empty() && !regions.iter().any(|r| r.chromosome == chrom && r.contains_position(position)) {
            continue;
        }
        *read_counts.entry(current_chr.clone()).or_default() += 1;

        let is_pass = filter == "." || filter == "PASS" || filter.split(';').any(|f| f == "PASS");
        match config.filter {
            FilterType::Pass if !is_pass => continue,
            FilterType::Panel if !is_pass => {
                // anything in the panel or passing
                let tier = info_value(info, "TIER").unwrap_or("UNKNOWN");
                if tier != "HOTSPOT" && tier != "PANEL" {
                    continue;
                }
            }
            _ => {}
        }

        // the positions the annotations look up: the variant's, and each base of an MNV
        let ref_len = ref_.len() as i32;
        let mut alt_first = alt.split(',').next().unwrap_or("");
        if alt_first == "*" || alt_first == "N" || alt_first == "." {
            alt_first = ref_; // VariantData.fromContext takes the ref
        }
        let entry = needed.by_chr.entry(current_chr.clone()).or_default();
        entry.push(position);
        if ref_len > 1 && alt_first.len() as i32 == ref_len {
            for i in 1..ref_len {
                entry.push(position + i);
            }
        }
        if every_task_reads_all {
            all.push(line);
        } else {
            by_chr.entry(current_chr.clone()).or_default().push(line);
        }
    }
    needed.finish();
    let by_chr: HashMap<String, std::sync::Arc<Vec<String>>> = if every_task_reads_all {
        let all = std::sync::Arc::new(all);
        let total: usize = read_counts.values().sum();
        read_counts = chromosomes.iter().map(|c| (c.clone(), total)).collect();
        chromosomes.iter().map(|c| (c.clone(), all.clone())).collect()
    } else {
        by_chr.into_iter().map(|(k, v)| (k, std::sync::Arc::new(v))).collect()
    };
    Ok(InputVariants { header, by_chr, read_counts, needed })
}
