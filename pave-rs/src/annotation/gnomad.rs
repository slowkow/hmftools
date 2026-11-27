//! GnomadCache (hmf-common) and GnomadAnnotation: gnomAD population frequencies, from one file or a
//! directory of per-chromosome files.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::{HashMap, HashSet};

use hmf_common::gz::LineReader;
use hmf_common::regions::{RefGenomeVersion, HUMAN_CHROMOSOMES};

use super::{parse_double, parse_int, split_fields, text};

pub const GNOMAD_FILE_ID: &str = "gnomad_variants";

pub enum GnomadSource {
    /// gnomad_freq_file: all chromosomes in one file, the chromosome in the first column
    File(String),
    /// gnomad_freq_dir: a file per chromosome, by versioned chromosome name
    Dir(HashMap<String, String>),
}

/// loadAllFrequencyFiles: a file for every chromosome, named `gnomad_variants_chr<N>.csv` or containing
/// `gnomad_variants_chr<N>_`
pub fn find_chromosome_files(dir: &str, version: RefGenomeVersion) -> Result<HashMap<String, String>> {
    let mut files: Vec<String> = Vec::new();
    for entry in std::fs::read_dir(dir).with_context(|| format!("failed to find Gnomad chromosome files in dir({dir})"))? {
        let entry = entry?;
        let path = entry.path();
        if std::fs::metadata(&path).map(|m| m.is_dir()).unwrap_or(false) {
            continue;
        }
        let s = path.to_string_lossy().into_owned();
        if s.contains(GNOMAD_FILE_ID) {
            files.push(s);
        }
    }
    // Java takes the first match in directory listing order; sorted here, so that the choice is defined
    files.sort();
    let mut map = HashMap::new();
    for chr in HUMAN_CHROMOSOMES {
        let no_id = format!("{dir}{GNOMAD_FILE_ID}_chr{chr}.csv");
        let with_id = format!("{GNOMAD_FILE_ID}_chr{chr}_");
        let chr_str = version.versioned_chromosome(chr);
        let Some(f) = files.iter().find(|x| x.ends_with(&no_id) || x.contains(&with_id)) else {
            bail!("missing Gnomad chromosome({chr_str}) file");
        };
        map.insert(chr_str, f.clone());
    }
    Ok(map)
}

/// The gnomAD entries at the needed positions, by chromosome
#[derive(Default)]
pub struct GnomadData {
    by_chr: HashMap<String, HashMap<i32, Vec<(String, String, f64)>>>,
    /// the chromosomes with a cache (any entry read), as Java has them
    chromosomes: HashSet<String>,
}

impl GnomadData {
    /// adds another chromosome's data (later loads replace a chromosome's entries)
    pub fn merge(&mut self, other: GnomadData) {
        self.by_chr.extend(other.by_chr);
        self.chromosomes.extend(other.chromosomes);
    }

    pub fn has_chromosome(&self, chromosome: &str) -> bool {
        self.chromosomes.contains(chromosome)
    }

    /// GnomadChrCache.getFrequency(position, ref, alt)
    fn frequency_at(&self, chromosome: &str, position: i32, ref_: &str, alt: &str) -> Option<f64> {
        self.by_chr.get(chromosome)?.get(&position)?.iter().find(|x| x.0 == ref_ && x.1 == alt).map(|x| x.2)
    }

    /// GnomadChrCache.getFrequency(isMnv, ...): an MNV takes the lowest frequency of its bases, if each has one
    pub fn frequency(&self, chromosome: &str, is_mnv: bool, ref_: &str, alt: &str, position: i32) -> Option<f64> {
        if !is_mnv {
            return self.frequency_at(chromosome, position, ref_, alt);
        }
        // PA6: Java keeps the lowest frequency of the bases with 0 meaning none yet, so a base of frequency
        // 0 is passed over by the next (no entry of the resource files is 0). Corrected: the lowest.
        let java = crate::compat::java("PA6");
        let mut min_freq: Option<f64> = None;
        for i in 0..ref_.len() {
            let freq = self.frequency_at(chromosome, position + i as i32, &ref_[i..i + 1], &alt[i..i + 1])?;
            let replace = match min_freq {
                None => true,
                Some(m) if java => m == 0.0 || freq < m,
                Some(m) => freq < m,
            };
            if replace {
                min_freq = Some(freq);
            }
        }
        Some(min_freq.unwrap_or(0.0))
    }

    #[cfg(test)]
    pub fn add_entry(&mut self, chromosome: &str, position: i32, ref_: &str, alt: &str, freq: f64) {
        self.chromosomes.insert(chromosome.to_string());
        self.by_chr.entry(chromosome.to_string()).or_default().entry(position).or_default().push((ref_.to_string(), alt.to_string(), freq));
    }

    /// loadChromosomeEntries: one file, with the chromosome given (per-chromosome files) or in the first
    /// column (a new chromosome block replacing any earlier one). Only the needed positions are kept.
    pub fn load(&mut self, filename: &str, file_chromosome: Option<&str>, needed: &dyn Fn(&str) -> Option<Vec<i32>>) -> Result<u64> {
        let mut reader = LineReader::open(filename).with_context(|| format!("failed to load Gnomad frequency file({filename})"))?;
        let _ = reader.next_line()?; // the header
        let (chr_index, pos_index) = if file_chromosome.is_some() { (None, 0) } else { (Some(0), 1) };
        let (ref_index, alt_index, freq_index) = (pos_index + 1, pos_index + 2, pos_index + 3);
        let mut current_chr: Option<String> = file_chromosome.map(|s| s.to_string());
        let mut current_needed: Vec<i32> = Vec::new();
        let mut started = false;
        let mut count = 0u64;
        let mut cursor = 0usize;
        // the first fields of a line without splitting it all: only a needed position's line is split
        let first_fields = |line: &[u8], n: usize| -> Vec<(usize, usize)> {
            let mut out = Vec::with_capacity(n);
            let mut start = 0;
            for i in memchr::memchr_iter(b',', line) {
                out.push((start, i));
                start = i + 1;
                if out.len() == n {
                    return out;
                }
            }
            out.push((start, line.len()));
            out
        };
        while let Some(line) = reader.next_line()? {
            count += 1;
            let bad = || anyhow!("invalid Gnomad frequency line in file({filename}): {}", text(line));
            let head = first_fields(line, pos_index + 1);
            if let Some(ci) = chr_index {
                let (s0, e0) = *head.get(ci).ok_or_else(bad)?;
                let chr = &line[s0..e0];
                if !started || current_chr.as_deref().map(|c| c.as_bytes()) != Some(chr) {
                    let c = text(chr);
                    // a new cache for each chromosome change
                    self.by_chr.insert(c.clone(), HashMap::new());
                    self.chromosomes.insert(c.clone());
                    current_needed = needed(&c).unwrap_or_default();
                    current_chr = Some(c);
                    cursor = 0;
                    started = true;
                }
            } else if !started {
                let c = current_chr.clone().unwrap();
                self.by_chr.insert(c.clone(), HashMap::new());
                self.chromosomes.insert(c.clone());
                current_needed = needed(&c).unwrap_or_default();
                started = true;
            }
            let (ps, pe) = *head.get(pos_index).ok_or_else(bad)?;
            let position = parse_int(&line[ps..pe]).ok_or_else(bad)?;
            let n = &current_needed;
            if cursor < n.len() && n[cursor] < position {
                cursor += n[cursor..].partition_point(|p| *p < position);
            } else if cursor > 0 && n[cursor - 1] >= position {
                cursor = n.partition_point(|p| *p < position);
            }
            if cursor >= n.len() || n[cursor] != position {
                continue;
            }
            let values = split_fields(line, b',');
            let freq = values.get(freq_index).and_then(|f| parse_double(f)).ok_or_else(bad)?;
            let r = text(values.get(ref_index).ok_or_else(bad)?);
            let a = text(values.get(alt_index).ok_or_else(bad)?);
            let chr = current_chr.as_ref().unwrap();
            self.by_chr.get_mut(chr).unwrap().entry(position).or_default().push((r, a, freq));
        }
        Ok(count)
    }
}
