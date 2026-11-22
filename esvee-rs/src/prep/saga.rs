//! Matching junctions to SAGA resource variants by location (SagaResource, SagaMatcherFactory,
//! SagaLocationMatcher). The resource is a FASTA whose record names describe an indel:
//! `<id>|<chr>:<pos>:<1 or -1>|<chr>:<pos>:<1 or -1>|<junction offset>[|<junction offset 2>]`.

use super::types::Orientation;
use anyhow::{bail, Result};
use std::collections::HashMap;
use std::sync::Arc;

pub const SAGA_LOCATION_MATCH_DISTANCE: i32 = 50;

#[derive(Debug)]
pub struct SagaVariant {
    pub id: String,
    pub chromosome: String,
    pub position1: i32,
    pub position2: i32,
    #[allow(dead_code)]
    pub insert_sequence: String,
}

impl SagaVariant {
    /// SagaVariant.toString: "<id> <chr>:<pos>:1-<chr>:<pos>:-1" (the first breakend is forward, the second reverse)
    pub fn display(&self) -> String {
        format!("{} {}:{}:1-{}:{}:-1", self.id, self.chromosome, self.position1, self.chromosome, self.position2)
    }
}

#[derive(Clone)]
pub struct IndexedBreakend {
    pub position: i32,
    pub orientation: Orientation,
    pub variant: Arc<SagaVariant>,
}

fn parse_breakend(s: &str) -> Result<(String, i32, Orientation)> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 {
        bail!("Invalid breakend string: {s}");
    }
    let orientation = match parts[2].parse::<i8>() {
        Ok(1) => Orientation::Forward,
        Ok(-1) => Orientation::Reverse,
        _ => bail!("Invalid breakend string: {s}"),
    };
    Ok((parts[0].to_string(), parts[1].parse::<i32>().map_err(|_| anyhow::anyhow!("Invalid breakend string: {s}"))?, orientation))
}

/// SagaAssembly.fromFastaRecord with the checks of SagaVariant and SagaAssembly
pub fn variant_from_fasta_record(label: &str, sequence: &str) -> Result<SagaVariant> {
    let mut parts: Vec<&str> = label.split('|').collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    if !(parts.len() == 4 || parts.len() == 5) {
        bail!("Invalid fasta label: {label}");
    }
    let (chr1, pos1, orient1) = parse_breakend(parts[1])?;
    let (chr2, pos2, orient2) = parse_breakend(parts[2])?;
    let offset1: usize = parts[3].parse().map_err(|_| anyhow::anyhow!("Invalid fasta label: {label}"))?;
    let offset2: Option<usize> = if parts.len() >= 5 { Some(parts[4].parse().map_err(|_| anyhow::anyhow!("Invalid fasta label: {label}"))?) } else { None };
    let insert_end = offset2.unwrap_or(offset1);
    if insert_end < offset1 || insert_end > sequence.len() {
        bail!("Invalid fasta label: {label}");
    }
    let insert_sequence = sequence[offset1..insert_end].to_string();
    if chr1 != chr2 || !(orient1.is_forward() && orient2.is_reverse()) || pos1 >= pos2 {
        bail!("Invalid SAGA variant: {label}");
    }
    let is_insert_by_position = pos1 + 1 == pos2;
    if is_insert_by_position && insert_sequence.is_empty() {
        bail!("Invalid SAGA variant: {label}");
    }
    if !insert_sequence.is_empty() && offset2.is_none() {
        bail!("Insert should have 2 junction offsets");
    }
    for &offset in [Some(offset1), offset2].iter().flatten() {
        if !(offset >= 1 && offset < sequence.len()) {
            bail!("Junction offsets out of bounds");
        }
    }
    if offset2.is_some_and(|o| offset1 >= o) {
        bail!("Junction offsets not in ascending order");
    }
    Ok(SagaVariant { id: parts[0].to_string(), chromosome: chr1, position1: pos1, position2: pos2, insert_sequence })
}

/// SagaMatcherFactory: the breakends of the resource by chromosome, in position order
pub struct SagaMatcherFactory {
    breakends: HashMap<String, Arc<Vec<IndexedBreakend>>>,
    pub variant_count: usize,
}

impl SagaMatcherFactory {
    pub fn from_variants(variants: Vec<SagaVariant>) -> Result<Self> {
        let mut ids = std::collections::HashSet::new();
        let mut by_chr: HashMap<String, Vec<IndexedBreakend>> = HashMap::new();
        let variant_count = variants.len();
        for v in variants {
            if !ids.insert(v.id.clone()) {
                bail!("Duplicate variant IDs in SAGA resource");
            }
            let v = Arc::new(v);
            let list = by_chr.entry(v.chromosome.clone()).or_default();
            list.push(IndexedBreakend { position: v.position1, orientation: Orientation::Forward, variant: v.clone() });
            list.push(IndexedBreakend { position: v.position2, orientation: Orientation::Reverse, variant: v });
        }
        let breakends = by_chr
            .into_iter()
            .map(|(chr, mut list)| {
                list.sort_by_key(|b| b.position);
                (chr, Arc::new(list))
            })
            .collect();
        Ok(SagaMatcherFactory { breakends, variant_count })
    }

    /// Reads the FASTA: each record name up to the first white space, and its bases.
    pub fn load(fasta_path: &str) -> Result<Self> {
        let text = std::fs::read_to_string(fasta_path).map_err(|e| anyhow::anyhow!("Failed to load SAGA resource: {e}"))?;
        let mut variants = Vec::new();
        let mut label: Option<String> = None;
        let mut sequence = String::new();
        for line in text.lines().chain(std::iter::once(">")) {
            if let Some(name) = line.strip_prefix('>') {
                if let Some(l) = label.take() {
                    variants.push(variant_from_fasta_record(&l, &sequence)?);
                }
                sequence.clear();
                let name = name.split_whitespace().next().unwrap_or("");
                if !name.is_empty() {
                    label = Some(name.to_string());
                }
            } else {
                sequence.push_str(line.trim_end());
            }
        }
        let factory = Self::from_variants(variants)?;
        crate::debug!("loaded {} SAGA variant assemblies from file({})", factory.variant_count, fasta_path);
        Ok(factory)
    }

    /// createLocationMatcher: a matcher over the breakends within the region and the match distance
    pub fn create_location_matcher(&self, chromosome: &str, start: i32, end: i32) -> SagaLocationMatcher {
        let (start, end) = (start - SAGA_LOCATION_MATCH_DISTANCE, end + SAGA_LOCATION_MATCH_DISTANCE);
        match self.breakends.get(chromosome) {
            Some(list) if !list.is_empty() => {
                let from = list.partition_point(|b| b.position < start);
                let to = list.partition_point(|b| b.position <= end);
                SagaLocationMatcher { chromosome: chromosome.to_string(), breakends: Some(list.clone()), from, to }
            }
            _ => SagaLocationMatcher { chromosome: chromosome.to_string(), breakends: None, from: 0, to: 0 },
        }
    }
}

pub struct SagaLocationMatcher {
    chromosome: String,
    breakends: Option<Arc<Vec<IndexedBreakend>>>,
    from: usize,
    to: usize,
}

impl SagaLocationMatcher {
    /// SagaLocationMatcher.match: the nearest breakend within the distance, of the junction's orientation when
    /// several are as near, then by variant id. Returns the variant as the junction file writes it.
    pub fn match_junction(&self, chromosome: &str, position: i32, orientation: Orientation) -> Option<String> {
        self.match_variant(chromosome, position, orientation).map(|v| v.display())
    }

    pub fn match_variant(&self, chromosome: &str, position: i32, orientation: Orientation) -> Option<Arc<SagaVariant>> {
        if chromosome != self.chromosome {
            return None;
        }
        let list = &self.breakends.as_ref()?[self.from..self.to];
        if list.is_empty() {
            return None;
        }
        let lo = list.partition_point(|b| b.position < position);
        let hi = list.partition_point(|b| b.position <= position);
        let candidates: Vec<(&IndexedBreakend, i32)> = if lo == hi {
            // not found exactly: the breakends either side
            let (s, e) = (lo.saturating_sub(1), (lo + 1).min(list.len()));
            let all: Vec<(&IndexedBreakend, i32)> =
                list[s..e].iter().map(|b| (b, (b.position - position).abs())).filter(|x| x.1 <= SAGA_LOCATION_MATCH_DISTANCE).collect();
            let min_distance = all.iter().map(|x| x.1).min()?;
            all.into_iter().filter(|x| x.1 == min_distance).collect()
        } else {
            list[lo..hi].iter().map(|b| (b, 0)).collect()
        };
        // Stream.min keeps the first of equal elements
        let mut best: Option<&IndexedBreakend> = None;
        for (b, _) in candidates {
            let key = |x: &IndexedBreakend| (x.orientation != orientation, x.variant.id.clone());
            if best.is_none_or(|cur| key(b) < key(cur)) {
                best = Some(b);
            }
        }
        best.map(|b| b.variant.clone())
    }
}
