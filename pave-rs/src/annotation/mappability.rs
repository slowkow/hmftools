//! Mappability and MappabilityChrCache: the GEM mappability of each variant's position.
//!
//! The file lists positions where the mappability changes (`chr  pos  pos  value`, 1-based: PAVE's
//! MappabilityBuilder writes the start plus 1). A variant takes the value of the entry at its position,
//! else that of the entry before it.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::{HashMap, HashSet};

use hmf_common::gz::LineReader;

use super::{parse_double, parse_int, text, NeededPositions};
use crate::compat;

#[derive(Default)]
pub struct MappabilityData {
    values: HashMap<String, HashMap<i32, f64>>,
    chromosomes: HashSet<String>,
}

impl MappabilityData {
    /// whether the file has the chromosome (Java's chromosome cache exists)
    pub fn has_chromosome(&self, chromosome: &str) -> bool {
        self.chromosomes.contains(chromosome)
    }
    /// MappabilityChrCache.annotateVariant: the value for a variant at the position
    pub fn value(&self, chromosome: &str, position: i32) -> Option<f64> {
        self.values.get(chromosome)?.get(&position).copied()
    }
}

/// Reads the file (Mappability.loadEntries) and works out the value of each needed position.
pub fn load(path: &str, needed: &NeededPositions, threaded: bool) -> Result<MappabilityData> {
    let mut reader = if threaded { LineReader::open_threaded(path) } else { LineReader::open(path) }.with_context(|| format!("failed to load mappability file({path})"))?;
    let mut data = MappabilityData::default();

    // PA2: Java adds 1 to each start as for a 0-based BED, but the file is 1-based already (its single-base
    // entries become start = end + 1), so a variant at a listed position gets the previous entry's value.
    // Corrected: the positions are read as they are.
    let start_offset = if compat::java("PA2") { 1 } else { 0 };

    let mut current_chr: Vec<u8> = Vec::new();
    let mut started = false;
    let mut current: &[i32] = &[];
    let mut cursor = 0usize;
    let mut prev_value: Option<f64> = None;
    let mut prev_end = i32::MIN;
    let mut unsorted_warned = false;
    let mut values: HashMap<i32, f64> = HashMap::new();

    let finish_chromosome = |chr: &[u8], values: &mut HashMap<i32, f64>, data: &mut MappabilityData| {
        if !chr.is_empty() || !values.is_empty() {
            data.values.entry(text(chr)).or_default().extend(values.drain());
        }
    };

    // the previous entry's value field, parsed only when a position takes it
    let mut prev_value_text: Vec<u8> = Vec::new();
    while let Some(line) = reader.next_line()? {
        // eg: chr1  10001  10001  0.0023; any other number of fields ends the run in Java
        let mut tabs = [0usize; 3];
        let mut n = 0;
        for i in memchr::memchr_iter(b'\t', line) {
            if n < 3 {
                tabs[n] = i;
            }
            n += 1;
        }
        if n != 3 {
            bail!("invalid mappability entry({})", text(line));
        }
        let chr = &line[..tabs[0]];
        if !started || chr != current_chr.as_slice() {
            if started {
                finish_chromosome(&current_chr, &mut values, &mut data);
            }
            current_chr = chr.to_vec();
            started = true;
            let c = text(chr);
            data.chromosomes.insert(c.clone());
            current = needed.get(&c).unwrap_or(&[]);
            cursor = 0;
            prev_value = None;
            prev_value_text.clear();
            prev_end = i32::MIN;
        }
        let bad = || anyhow!("invalid mappability entry({})", text(line));
        let end = parse_int(&line[tabs[1] + 1..tabs[2]]).ok_or_else(bad)?;
        let value_field = &line[tabs[2] + 1..];
        if end < prev_end && !unsorted_warned {
            crate::warn!("mappability file({}) is not in position order at {}:{}", path, text(chr), end);
            unsorted_warned = true;
        }
        prev_end = end;

        // each needed position up to this entry's end has this entry as the first ending at or after it
        if cursor < current.len() && current[cursor] <= end {
            let start = parse_int(&line[tabs[0] + 1..tabs[1]]).ok_or_else(bad)? + start_offset;
            let value = parse_double(value_field).ok_or_else(bad)?;
            if prev_value.is_none() && !prev_value_text.is_empty() {
                prev_value = parse_double(&prev_value_text);
            }
            while cursor < current.len() && current[cursor] <= end {
                let p = current[cursor];
                if start <= p {
                    values.insert(p, value);
                } else if let Some(v) = prev_value {
                    // the previous entry's if the next is past the variant
                    values.insert(p, v);
                }
                cursor += 1;
            }
        }
        prev_value = None;
        prev_value_text.clear();
        prev_value_text.extend_from_slice(value_field);
    }
    if started {
        finish_chromosome(&current_chr, &mut values, &mut data);
    }
    Ok(data)
}
