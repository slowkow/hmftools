//! PonCache (hmf-common) and PonAnnotation: the panel of normals, its tier filters and the artefact PON.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;

use hmf_common::gz::LineReader;
use hmf_common::java_text;

use super::{parse_int, text, NeededPositions};
use crate::variant::{MultiPonStatus, VariantTier};

#[derive(Clone, Debug)]
pub struct PonVariantData {
    pub ref_: String,
    pub alt: String,
    pub samples: i32,
    pub max_sample_reads: i32,
    pub total_sample_reads: i32,
    pub multi_pon: MultiPonStatus,
}

impl PonVariantData {
    /// meanReadCount: the rounded mean of the samples' reads
    pub fn mean_read_count(&self) -> i32 {
        if self.samples > 0 { hmf_common::java_math::java_round(self.total_sample_reads as f64 / self.samples as f64) as i32 } else { 0 }
    }
}

/// The PON entries at the needed positions, by the file's chromosome name
#[derive(Default)]
pub struct PonData {
    by_chr: HashMap<String, HashMap<i32, Vec<PonVariantData>>>,
    pub entries_read: u64,
}

impl PonData {
    /// PonChrCache.getPonData: the first entry at the position with the ref and alt
    pub fn get(&self, chromosome: &str, position: i32, ref_: &str, alt: &str) -> Option<&PonVariantData> {
        self.by_chr.get(chromosome)?.get(&position)?.iter().find(|x| x.ref_ == ref_ && x.alt == alt)
    }
    pub fn has_chromosome(&self, chromosome: &str) -> bool {
        self.by_chr.contains_key(chromosome)
    }

    /// PonChrCache.addEntry
    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    pub fn add_entry(&mut self, chr: &str, position: i32, ref_: &str, alt: &str, samples: i32, max: i32, total: i32, multi: MultiPonStatus) {
        let e = PonVariantData { ref_: ref_.to_string(), alt: alt.to_string(), samples, max_sample_reads: max, total_sample_reads: total, multi_pon: multi };
        self.by_chr.entry(chr.to_string()).or_default().entry(position).or_default().push(e);
    }

    #[cfg(test)]
    pub fn clear(&mut self) {
        for m in self.by_chr.values_mut() {
            m.clear();
        }
    }
}

/// The byte ranges of the first `n` tab-separated fields of a line (fewer if the line has fewer)
fn field_ranges(line: &[u8], n: usize, out: &mut Vec<(usize, usize)>) {
    out.clear();
    let mut start = 0;
    for i in memchr::memchr_iter(b'\t', line) {
        out.push((start, i));
        start = i + 1;
        if out.len() == n {
            return;
        }
    }
    out.push((start, line.len()));
}

/// Reads a PON file (PonCache.initialiseFile and loadPonEntries), keeping the entries at the needed
/// positions only (all of them with `needed` None). Every chromosome block is read; a chromosome seen
/// again later adds to its entries, as Java's cache does.
pub fn load_pon(path: &str, needed: Option<&NeededPositions>, threaded: bool) -> Result<PonData> {
    let mut reader = if threaded { LineReader::open_threaded(path) } else { LineReader::open(path) }.with_context(|| format!("failed to load PON file({path})"))?;
    let Some(header) = reader.next_line()? else { bail!("empty PON file({path})") };
    let header = text(header);
    let mut data = PonData::default();
    let columns = header.split('\t').count();
    if columns < 4 {
        crate::error!("pon file({}) has insufficient column count({})", path, columns);
        return Ok(data);
    }
    let f = java_text::fields_index_map(&header, '\t');
    let col = |name: &str| f.get(name).copied();
    let req = |name: &str| col(name).ok_or_else(|| anyhow!("NullPointerException: PON file({path}) has no column {name}"));
    let chr_i = req("Chromosome")?;
    let pos_i = req("Position")?;
    let ref_i = req("Ref")?;
    let alt_i = req("Alt")?;
    // support for the v37 PON's column names
    let samples_i = col("SampleCount").or_else(|| col("SamplesCount"));
    let max_i = col("MaxSampleReadCount").or_else(|| col("MaxSampleReads"));
    let total_i = col("TotalReadCount").or_else(|| col("TotalReads"));
    let multi_i = col("MultiPonStatus");
    let max_col = [Some(chr_i), Some(pos_i), Some(ref_i), Some(alt_i), samples_i, max_i, total_i, multi_i].iter().flatten().copied().max().unwrap();
    let key_col = chr_i.max(pos_i);

    let mut ranges: Vec<(usize, usize)> = Vec::with_capacity(max_col + 2);
    let mut current_chr: Vec<u8> = Vec::new();
    let mut current_needed: Option<&[i32]> = None;
    let mut have_current = false;
    let mut cursor = 0usize;

    while let Some(line) = reader.next_line()? {
        data.entries_read += 1;
        field_ranges(line, key_col + 1, &mut ranges);
        if ranges.len() <= key_col {
            bail!("PON file({path}): short line: {}", text(line));
        }
        let chr = &line[ranges[chr_i].0..ranges[chr_i].1];
        if !have_current || chr != current_chr.as_slice() {
            current_chr = chr.to_vec();
            have_current = true;
            // the chromosome's cache exists once a line of it is read
            data.by_chr.entry(text(chr)).or_default();
            current_needed = needed.map(|n| n.get(&text(chr)).unwrap_or(&[]));
            cursor = 0;
        }
        let position = parse_int(&line[ranges[pos_i].0..ranges[pos_i].1])
            .ok_or_else(|| anyhow!("NumberFormatException: PON file({path}) position in line: {}", text(line)))?;
        if let Some(n) = current_needed {
            // the file is in position order within a chromosome, so the cursor moves forward; a step back is
            // searched for
            if cursor < n.len() && n[cursor] < position {
                cursor += n[cursor..].partition_point(|p| *p < position);
            } else if cursor > 0 && n[cursor - 1] >= position {
                cursor = n.partition_point(|p| *p < position);
            }
            if cursor >= n.len() || n[cursor] != position {
                continue;
            }
        }
        field_ranges(line, max_col + 1, &mut ranges);
        let field = |i: usize| -> Result<&[u8]> {
            ranges.get(i).map(|r| &line[r.0..r.1]).ok_or_else(|| anyhow!("PON file({path}): short line: {}", text(line)))
        };
        let int_field = |i: Option<usize>| -> Result<i32> {
            match i {
                None => Ok(0),
                Some(i) => parse_int(field(i)?).ok_or_else(|| anyhow!("NumberFormatException: PON file({path}) line: {}", text(line))),
            }
        };
        let entry = PonVariantData {
            ref_: text(field(ref_i)?),
            alt: text(field(alt_i)?),
            samples: int_field(samples_i)?,
            max_sample_reads: int_field(max_i)?,
            total_sample_reads: int_field(total_i)?,
            multi_pon: match multi_i {
                Some(i) => MultiPonStatus::value_of(&text(field(i)?))?,
                None => MultiPonStatus::Base,
            },
        };
        data.by_chr.entry(text(&current_chr)).or_default().entry(position).or_default().push(entry);
    }
    Ok(data)
}

/// The PON filter of each tier: (sample count, max read count) needed to filter
pub struct PonFilters {
    filters: HashMap<VariantTier, (i32, i32)>,
}

impl PonFilters {
    /// PonCache.loadFilters: "TIER:SAMPLES:MAX_READS" items separated by ';'
    pub fn parse(config: &str) -> Result<PonFilters> {
        let mut filters = HashMap::new();
        for item in java_text::split_all(config, ';') {
            let values = java_text::split_all(item, ':');
            let tier = VariantTier::from_string(values.first().copied().unwrap_or(""));
            let bad = || anyhow!("invalid PON filters({config})");
            let samples = values.get(1).and_then(|v| java_text::parse_int(v)).ok_or_else(bad)?;
            let max_reads = values.get(2).and_then(|v| java_text::parse_int(v)).ok_or_else(bad)?;
            filters.insert(tier, (samples, max_reads));
            crate::info!("loaded PON filter: tier({}) sampleCount({}) maxReads({})", crate::variant::tier_name(tier), samples, max_reads);
        }
        Ok(PonFilters { filters })
    }

    /// filterOnTierCriteria: the tier's filter, else that of UNKNOWN
    pub fn filter_on_tier_criteria(&self, tier: VariantTier, sample_count: i32, max_sample_reads: i32) -> bool {
        let f = self.filters.get(&tier).or_else(|| self.filters.get(&VariantTier::Unknown));
        match f {
            None => false,
            Some((samples, max_reads)) => sample_count >= *samples && max_sample_reads >= *max_reads,
        }
    }
}
