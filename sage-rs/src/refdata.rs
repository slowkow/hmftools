//! Reference data: Ensembl genes/transcripts, driver gene panel, hotspots, high-confidence BED.
//! Ports ReferenceData plus the hmf-common loaders it uses.

use crate::types::*;
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read as IoRead};

pub fn open_text(path: &str) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| format!("opening {path}"))?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

pub fn read_lines(path: &str) -> Result<Vec<String>> {
    let mut s = String::new();
    open_text(path)?.read_to_string(&mut s)?;
    Ok(s.lines().map(|x| x.to_string()).collect())
}

fn header_index(header: &str, delim: char) -> HashMap<String, usize> {
    header.split(delim).enumerate().map(|(i, n)| (n.to_string(), i)).collect()
}

#[derive(Clone, Debug)]
pub struct GeneData {
    pub gene_id: String,
    pub name: String,
    pub start: i32,
    pub end: i32,
}

#[derive(Clone, Debug)]
pub struct ExonData {
    pub start: i32,
    pub end: i32,
    pub rank: i32,
    pub phase_start: i32,
}

#[derive(Clone, Debug)]
pub struct TranscriptData {
    pub strand: i8,
    pub trans_start: i32,
    pub trans_end: i32,
    pub coding_start: Option<i32>,
    pub coding_end: Option<i32>,
    pub exons: Vec<ExonData>,
}

pub struct RefData {
    /// versioned chromosome -> regions
    pub panel: HashMap<String, Vec<BaseRegion>>,
    pub hotspots: HashMap<String, Vec<SimpleVariant>>,
    pub high_confidence: HashMap<String, Vec<BaseRegion>>,
    /// canonical transcripts per chromosome, in gene file order
    pub transcripts: HashMap<String, Vec<TranscriptData>>,
    pub hla_genes: Vec<BaseRegion>,
    pub hla_chromosome: String,
}

fn versioned_chr(chr: &str, v38: bool) -> String {
    let stripped = chr.strip_prefix("chr").unwrap_or(chr);
    if v38 {
        format!("chr{stripped}")
    } else {
        stripped.to_string()
    }
}

const HLA_GENES: [&str; 6] = ["HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"];
const SPLICE_SIZE: i32 = 10;

fn load_genes(dir: &str, v38: bool) -> Result<(Vec<String>, HashMap<String, Vec<GeneData>>)> {
    let path = format!("{}/ensembl_gene_data.csv", dir.trim_end_matches('/'));
    let lines = read_lines(&path)?;
    let idx = header_index(&lines[0], ',');
    let (gi, ni, ci, si, ei) =
        (idx["GeneId"], idx["GeneName"], idx["Chromosome"], idx["GeneStart"], idx["GeneEnd"]);
    let mut order = Vec::new();
    let mut map: HashMap<String, Vec<GeneData>> = HashMap::new();
    for line in &lines[1..] {
        let items: Vec<&str> = line.split(',').collect();
        let chr = versioned_chr(items[ci], v38);
        if !map.contains_key(&chr) {
            order.push(chr.clone());
        }
        map.entry(chr).or_default().push(GeneData {
            gene_id: items[gi].to_string(),
            name: items[ni].to_string(),
            start: items[si].parse()?,
            end: items[ei].parse()?,
        });
    }
    Ok((order, map))
}

/// EnsemblDataLoader.loadTranscriptData with canonicalOnly=true, cacheExons=true, includeNonEnsembl=false.
/// Returns the canonical transcript per gene id.
fn load_canonical_transcripts(dir: &str) -> Result<HashMap<String, TranscriptData>> {
    let path = format!("{}/ensembl_trans_exon_data.csv", dir.trim_end_matches('/'));
    let mut reader = open_text(&path)?;
    let mut header = String::new();
    reader.read_line(&mut header)?;
    let idx = header_index(header.trim_end(), ',');
    let gi = idx["GeneId"];
    let ci = idx["CanonicalTranscriptId"];
    let sti = idx["Strand"];
    let ti = idx["TransId"];
    let tni = idx["TransName"];
    let tsi = idx["TransStart"];
    let tei = idx["TransEnd"];
    let eri = idx["ExonRank"];
    let esi = idx["ExonStart"];
    let eei = idx["ExonEnd"];
    let epi = idx["ExonPhase"];
    let csi = idx["CodingStart"];
    let cei = idx["CodingEnd"];

    let mut result: HashMap<String, TranscriptData> = HashMap::new();
    let mut current_gene = String::new();
    let mut current_trans_id: Option<i64> = None;
    let mut current: Option<TranscriptData> = None;
    let mut line = String::new();
    let flush = |gene: &str, t: Option<TranscriptData>, result: &mut HashMap<String, TranscriptData>| {
        if let Some(t) = t {
            // getCanonicalTranscriptData returns the first canonical transcript in the gene's list
            result.entry(gene.to_string()).or_insert(t);
        }
    };
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let items: Vec<&str> = line.trim_end_matches(['\n', '\r']).split(',').collect();
        let gene_id = items[gi];
        let trans_id: i64 = items[ti].parse()?;
        if current_trans_id != Some(trans_id) || current.is_none() {
            // note: Java keeps the previous transcript as current when the new one is skipped, and
            // skipped transcripts never reach the exon-add code since their id never matches
            let canonical: i64 = items[ci].parse()?;
            if canonical != trans_id || !items[tni].starts_with("ENST") {
                continue;
            }
            let prev = current.take();
            flush(&current_gene, prev, &mut result);
            current_gene = gene_id.to_string();
            let parse_opt = |s: &str| -> Option<i32> {
                if s.eq_ignore_ascii_case("NULL") {
                    None
                } else {
                    s.parse().ok()
                }
            };
            current = Some(TranscriptData {
                strand: items[sti].parse()?,
                trans_start: items[tsi].parse()?,
                trans_end: items[tei].parse()?,
                coding_start: parse_opt(items[csi]),
                coding_end: parse_opt(items[cei]),
                exons: Vec::new(),
            });
            current_trans_id = Some(trans_id);
        }
        if let Some(t) = current.as_mut() {
            t.exons.push(ExonData {
                start: items[esi].parse()?,
                end: items[eei].parse()?,
                rank: items[eri].parse()?,
                phase_start: items[epi].parse()?,
            });
        }
    }
    flush(&current_gene, current.take(), &mut result);
    Ok(result)
}

/// DriverGeneRegions.getTranscriptRegions(includeUTR=false, includeSplice=true)
fn transcript_regions(t: &TranscriptData) -> Vec<BaseRegion> {
    let (start, end) = match (t.coding_start, t.coding_end) {
        (Some(s), Some(e)) => (s, e),
        _ => (t.trans_start, t.trans_end),
    };
    let n = t.exons.len();
    let mut regions = Vec::new();
    for (i, exon) in t.exons.iter().enumerate() {
        let mut es = exon.start;
        let mut ee = exon.end;
        if i > 0 {
            es -= SPLICE_SIZE;
        }
        if i + 1 < n {
            ee += SPLICE_SIZE;
        }
        if positions_overlap(start, end, es, ee) {
            regions.push(BaseRegion::new(start.max(es), end.min(ee)));
        }
    }
    regions
}

fn load_driver_gene_names(path: &str) -> Result<Vec<String>> {
    let lines = read_lines(path)?;
    let idx = header_index(&lines[0], '\t');
    let gi = *idx.get("gene").context("driver gene panel missing 'gene' column")?;
    Ok(lines[1..].iter().map(|l| l.split('\t').nth(gi).unwrap_or("").to_string()).collect())
}

/// VariantHotspotFile.loadHotspotVcf: CHROM, POS, REF and first ALT, human chromosomes only.
fn load_hotspots(path: &str) -> Result<HashMap<String, Vec<SimpleVariant>>> {
    let mut map: HashMap<String, Vec<SimpleVariant>> = HashMap::new();
    for line in read_lines(path)? {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 || !is_human_chromosome(f[0]) {
            continue;
        }
        let alt = f[4].split(',').next().unwrap();
        map.entry(f[0].to_string())
            .or_default()
            .push(SimpleVariant::new(f[1].parse()?, f[3].as_bytes(), alt.as_bytes()));
    }
    Ok(map)
}

/// BedFileReader.loadBedFileChrMap, including its quirk of starting a fresh list whenever the
/// chromosome changes (a non-contiguous BED loses the earlier block).
fn load_bed_chr_map(path: &str, v38: bool) -> Result<HashMap<String, Vec<BaseRegion>>> {
    let mut map: HashMap<String, Vec<BaseRegion>> = HashMap::new();
    let mut current: Option<String> = None;
    for line in read_lines(path)? {
        if line.contains("Chromosome") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            bail!("invalid slice BED entry: {line}");
        }
        if !is_human_chromosome(f[0]) {
            continue;
        }
        let chr = versioned_chr(f[0], v38);
        if current.as_deref() != Some(chr.as_str()) {
            current = Some(chr.clone());
            map.insert(chr.clone(), Vec::new());
        }
        map.get_mut(&chr).unwrap().push(BaseRegion::new(f[1].parse::<i32>()? + 1, f[2].parse()?));
    }
    Ok(map)
}

/// ReferenceData.loadHotspots: add each hotspot to the panel if not already covered.
fn merge_hotspots_into_panel(panel: &mut Vec<BaseRegion>, hotspots: &[SimpleVariant]) {
    for hotspot in hotspots {
        let pos = hotspot.pos;
        let mut covered = false;
        let mut index = 0;
        while index < panel.len() {
            let region = panel[index];
            if region.contains(pos) {
                covered = true;
                break;
            } else if pos == region.start - 1 {
                covered = true;
                panel[index].start = pos;
                break;
            } else if pos == region.end + 1 {
                covered = true;
                panel[index].end = pos;
                if index + 1 < panel.len() {
                    let next = panel[index + 1];
                    if next.start <= panel[index].end + 1 {
                        panel[index].end = next.end;
                        panel.remove(index + 1);
                    }
                }
                break;
            } else if pos < region.start {
                break;
            }
            index += 1;
        }
        if !covered {
            panel.insert(index, BaseRegion::new(pos, pos));
        }
    }
}

pub struct RefDataPaths<'a> {
    pub ensembl_dir: Option<&'a str>,
    pub driver_gene_panel: Option<&'a str>,
    pub panel_bed: Option<&'a str>,
    pub hotspots: Option<&'a str>,
    pub high_confidence_bed: Option<&'a str>,
    pub v38: bool,
}

impl RefData {
    pub fn load(p: &RefDataPaths) -> Result<RefData> {
        let v38 = p.v38;
        let hla_chromosome = versioned_chr("6", v38);
        let mut panel: HashMap<String, Vec<BaseRegion>> = HashMap::new();
        let mut transcripts: HashMap<String, Vec<TranscriptData>> = HashMap::new();
        let mut hla_genes = Vec::new();

        let mut genes: Option<(Vec<String>, HashMap<String, Vec<GeneData>>)> = None;
        let mut canonical: HashMap<String, TranscriptData> = HashMap::new();
        if let Some(dir) = p.ensembl_dir {
            let g = load_genes(dir, v38)?;
            canonical = load_canonical_transcripts(dir)?;
            for chr in &g.0 {
                let list: Vec<TranscriptData> =
                    g.1[chr].iter().filter_map(|gene| canonical.get(&gene.gene_id).cloned()).collect();
                transcripts.insert(chr.clone(), list);
            }
            if let Some(list) = g.1.get(&hla_chromosome) {
                for gene in list {
                    if HLA_GENES.contains(&gene.name.as_str()) {
                        hla_genes.push(BaseRegion::new(gene.start, gene.end));
                    }
                }
            }
            genes = Some(g);
        }

        if let Some(path) = p.driver_gene_panel {
            let names: HashSet<String> = load_driver_gene_names(path)?.into_iter().collect();
            let (_, gene_map) = genes.as_ref().context("driver_gene_panel requires ensembl_data_dir")?;
            let known: HashSet<&str> = gene_map.values().flatten().map(|g| g.name.as_str()).collect();
            let invalid: Vec<&String> = names.iter().filter(|n| !known.contains(n.as_str())).collect();
            if !invalid.is_empty() {
                bail!("Invalid non-Ensembl gene names: {invalid:?}");
            }
            // DriverGeneRegions.buildDriverGeneBaseRegions over the human chromosomes
            let mut chrs: Vec<String> = (1..=22).map(|i| i.to_string()).collect();
            chrs.push("X".into());
            chrs.push("Y".into());
            for c in chrs {
                let chr = versioned_chr(&c, v38);
                let mut regions: Vec<BaseRegion> = Vec::new();
                if let Some(list) = gene_map.get(&chr) {
                    for gene in list {
                        if !names.contains(&gene.name) {
                            continue;
                        }
                        if let Some(t) = canonical.get(&gene.gene_id) {
                            regions.extend(transcript_regions(t));
                        }
                    }
                }
                regions.sort_by_key(|r| r.start); // stable, as Collections.sort
                let mut index = 0;
                while index + 1 < regions.len() {
                    let next_index = index + 1;
                    while next_index < regions.len() {
                        let next = regions[next_index];
                        if regions[index].end >= next.start {
                            if next.end > regions[index].end {
                                regions[index].end = next.end;
                            }
                            regions.remove(next_index);
                        } else {
                            break;
                        }
                    }
                    index += 1;
                }
                panel.insert(chr, regions);
            }
        } else if let Some(path) = p.panel_bed {
            panel = load_bed_chr_map(path, v38)?;
        }

        let mut hotspots = HashMap::new();
        if let Some(path) = p.hotspots {
            hotspots = load_hotspots(path)?;
            for (chr, list) in &hotspots {
                merge_hotspots_into_panel(panel.entry(chr.clone()).or_default(), list);
            }
        }

        let high_confidence = match p.high_confidence_bed {
            Some(path) => load_bed_chr_map(path, v38)?,
            None => HashMap::new(),
        };

        Ok(RefData { panel, hotspots, high_confidence, transcripts, hla_genes, hla_chromosome })
    }

    /// HlaCommon.overlaps
    pub fn hla_overlaps(&self, chr: &str, start: i32, end: i32) -> bool {
        chr == self.hla_chromosome && self.hla_genes.iter().any(|g| positions_overlap(start, end, g.start, g.end))
    }
    /// HlaCommon.containsPosition
    pub fn hla_contains(&self, chr: &str, pos: i32) -> bool {
        chr == self.hla_chromosome && self.hla_genes.iter().any(|g| g.contains(pos))
    }
}
