//! hmf-common `ensemblcache`: the Ensembl data cache (EnsemblDataCache and EnsemblDataLoader), as far as
//! genes, transcripts and exons go. The files are `ensembl_gene_data.<37|38>.csv` and
//! `ensembl_trans_exon_data.<37|38>.csv` in the data directory, or the same names without the version
//! (`checkOldFilename`).
//!
//! The maps iterate as Java's HashMaps do (`JavaHashMap`), since callers such as
//! `createGeneNameIdMap` (the last gene of a name wins) and the BamTools gene coverage depend on that
//! order. Protein domains, splice-acceptor positions and amino acids are not ported yet.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashSet;

use crate::java_hash::JavaHashMap;
use crate::java_text::{self, parse_byte, parse_int};
use crate::regions::RefGenomeVersion;

pub const ENSEMBL_DATA_DIR: &str = "ensembl_data_dir";
pub const ENSEMBL_DATA_DIR_CFG: &str = "Ensembl data file directory";
pub const ENSEMBL_TRANSCRIPT_PREFIX: &str = "ENST";

/// GeneData
#[derive(Clone, Debug, PartialEq)]
pub struct GeneData {
    pub gene_id: String,
    pub gene_name: String,
    pub chromosome: String,
    pub strand: i8,
    pub gene_start: i32,
    pub gene_end: i32,
    pub karyotype_band: String,
    /// set only when gene synonyms are required; "" otherwise
    pub synonyms: String,
}

impl GeneData {
    pub fn forward_strand(&self) -> bool {
        self.strand == 1
    }
    /// hasSynonym: a substring test on the synonyms text, as Java does it
    pub fn has_synonym(&self, name: &str) -> bool {
        self.synonyms.contains(name)
    }
    pub fn length(&self) -> i32 {
        self.gene_end - self.gene_start
    }
}

/// ExonData
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExonData {
    pub trans_id: i32,
    pub start: i32,
    pub end: i32,
    pub rank: i32,
    pub phase_start: i32,
    pub phase_end: i32,
}

impl ExonData {
    pub fn length(&self) -> i32 {
        self.end - self.start
    }
    pub fn base_length(&self) -> i32 {
        self.length() + 1
    }
}

/// TranscriptData
#[derive(Clone, Debug, PartialEq)]
pub struct TranscriptData {
    pub trans_id: i32,
    pub trans_name: String,
    pub gene_id: String,
    pub is_canonical: bool,
    pub strand: i8,
    pub trans_start: i32,
    pub trans_end: i32,
    pub coding_start: Option<i32>,
    pub coding_end: Option<i32>,
    pub bio_type: String,
    pub ref_seq_id: Option<String>,
    pub exons: Vec<ExonData>,
}

impl TranscriptData {
    pub fn non_coding(&self) -> bool {
        self.coding_start.is_none()
    }
    pub fn pos_strand(&self) -> bool {
        self.strand == 1
    }
    pub fn length(&self) -> i32 {
        self.trans_end - self.trans_start
    }
}

/// EnsemblDataLoader file names: `<id>.<37|38>.csv`
pub fn ensembl_gene_data_file(version: RefGenomeVersion) -> String {
    format!("ensembl_gene_data.{}.csv", version.identifier())
}
pub fn ensembl_trans_exon_data_file(version: RefGenomeVersion) -> String {
    format!("ensembl_trans_exon_data.{}.csv", version.identifier())
}

/// String.replaceAll(regex, replacement) for the patterns "37.csv" and "38.csv", whose '.' matches any
/// character
fn replace_version_pattern(s: &str, digits: &str) -> String {
    let b = s.as_bytes();
    let pat = digits.as_bytes();
    let mut out = String::with_capacity(s.len());
    let mut i = 0;
    while i < b.len() {
        // a match needs the two digits, one character (any, but here always ASCII or the start of one),
        // then "csv"
        if b.len() - i >= 6 && &b[i..i + 2] == pat && s.is_char_boundary(i + 2) {
            let c_len = s[i + 2..].chars().next().map(|c| c.len_utf8()).unwrap_or(1);
            let j = i + 2 + c_len;
            if j + 3 <= b.len() && &b[j..j + 3] == b"csv" {
                out.push_str("csv");
                i = j + 3;
                continue;
            }
        }
        let c = s[i..].chars().next().unwrap();
        out.push(c);
        i += c.len_utf8();
    }
    out
}

/// EnsemblDataLoader.checkOldFilename: the versioned name if it exists, else the name without the
/// version if that exists, else the versioned name
pub fn check_old_filename(filename: &str) -> String {
    if std::path::Path::new(filename).exists() {
        return filename.to_string();
    }
    let old = replace_version_pattern(&replace_version_pattern(filename, "37"), "38");
    if std::path::Path::new(&old).exists() {
        return old;
    }
    filename.to_string()
}

fn check_add_dir_separator(path: &str) -> String {
    if path.is_empty() || path.ends_with('/') { path.to_string() } else { format!("{path}/") }
}

/// EnsemblDataCache with the data it loads
#[derive(Clone)]
pub struct EnsemblDataCache {
    data_path: String,
    version: RefGenomeVersion,
    /// transcripts keyed by gene id
    transcripts_by_gene: JavaHashMap<String, Vec<TranscriptData>>,
    /// genes keyed by chromosome, in file order
    chr_genes: JavaHashMap<String, Vec<GeneData>>,
    /// keyed by gene id: (chromosome, index in its list)
    gene_ids: JavaHashMap<String, (String, usize)>,
    /// keyed by gene name: the last gene of each name in the chromosome map's iteration order
    gene_names: JavaHashMap<String, (String, usize)>,
    require_exons: bool,
    require_protein_domains: bool,
    require_splice_positions: bool,
    canonical_only: bool,
    require_gene_synonyms: bool,
    require_non_ensembl_transcripts: bool,
    restricted_gene_ids: Vec<String>,
    /// protein features by transcript id (getTranscriptProteinDataMap), loaded when protein domains are required
    protein_data: JavaHashMap<i32, Vec<TranscriptProteinData>>,
    /// the preceding gene's splice acceptor position by transcript id, when splice positions are required
    trans_splice_acceptor_pos: std::collections::HashMap<i32, i32>,
}

/// TranscriptProteinData
#[derive(Clone, Debug, PartialEq)]
pub struct TranscriptProteinData {
    pub transcript_id: i32,
    pub translation_id: i32,
    pub protein_feature_id: i32,
    pub seq_start: i32,
    pub seq_end: i32,
    pub hit_description: String,
}

pub const BIOTYPE_PROTEIN_CODING: &str = "protein_coding";
pub const BIOTYPE_NONSENSE_MED_DECAY: &str = "nonsense_mediated_decay";
pub const BIOTYPE_RETAINED_INTRON: &str = "retained_intron";
pub const BIOTYPE_PROCESSED_TRANS: &str = "processed_transcript";
pub const BIOTYPE_LINC_RNA: &str = "lincRNA";

pub fn ensembl_protein_data_file(version: RefGenomeVersion) -> String {
    format!("ensembl_protein_features.{}.csv", version.identifier())
}
pub fn ensembl_splice_data_file(version: RefGenomeVersion) -> String {
    format!("ensembl_trans_splice_data.{}.csv", version.identifier())
}

impl EnsemblDataCache {
    /// new EnsemblDataCache(dataPath, refGenomeVersion): exons required, all transcripts
    pub fn new(data_path: &str, version: RefGenomeVersion) -> Self {
        EnsemblDataCache {
            data_path: check_add_dir_separator(data_path),
            version,
            transcripts_by_gene: JavaHashMap::new(),
            chr_genes: JavaHashMap::new(),
            gene_ids: JavaHashMap::new(),
            gene_names: JavaHashMap::new(),
            require_exons: true,
            require_protein_domains: false,
            require_splice_positions: false,
            canonical_only: false,
            require_gene_synonyms: false,
            require_non_ensembl_transcripts: false,
            restricted_gene_ids: Vec::new(),
            protein_data: JavaHashMap::new(),
            trans_splice_acceptor_pos: std::collections::HashMap::new(),
        }
    }

    pub fn ref_genome_version(&self) -> RefGenomeVersion {
        self.version
    }
    pub fn data_path(&self) -> &str {
        &self.data_path
    }

    /// setRequiredData(exons, proteinDomains, splicePositions, canonicalOnly)
    pub fn set_required_data(&mut self, exons: bool, protein_domains: bool, splice_positions: bool, canonical_only: bool) {
        self.require_exons = exons;
        self.require_protein_domains = protein_domains;
        self.require_splice_positions = splice_positions;
        self.canonical_only = canonical_only;
    }
    pub fn set_require_gene_synonyms(&mut self) {
        self.require_gene_synonyms = true;
    }
    pub fn set_require_non_ensembl_transcripts(&mut self) {
        self.require_non_ensembl_transcripts = true;
    }
    pub fn set_restricted_gene_id_list(&mut self, gene_ids: &[String]) {
        self.restricted_gene_ids = gene_ids.to_vec();
    }

    /// load(delayTranscriptLoading): the genes, then (unless delayed) the transcripts. Java returns false
    /// where a file is missing; that is an error here, as are the exceptions Java does not catch.
    pub fn load(&mut self, delay_transcript_loading: bool) -> Result<()> {
        let restricted: Vec<String> = self.restricted_gene_ids.clone();
        load_ensembl_gene_data(&self.data_path, &restricted, &mut self.chr_genes, self.version, self.require_gene_synonyms)?;
        if !delay_transcript_loading {
            load_transcript_data(
                &self.data_path,
                self.version,
                &mut self.transcripts_by_gene,
                &restricted,
                self.require_exons,
                self.canonical_only,
                self.require_non_ensembl_transcripts,
                &[],
            )?;
            if self.require_protein_domains {
                load_transcript_protein_data(&self.data_path, self.version, &mut self.protein_data, &HashSet::new())?;
            }
            if self.require_splice_positions {
                load_transcript_splice_acceptor_data(&self.data_path, self.version, &mut self.trans_splice_acceptor_pos, &HashSet::new())?;
            }
        }
        Ok(())
    }

    /// getTranscriptProteinDataMap
    pub fn transcript_protein_data_map(&self) -> &JavaHashMap<i32, Vec<TranscriptProteinData>> {
        &self.protein_data
    }

    /// findPrecedingGeneSpliceAcceptorPosition: -1 if not known
    pub fn find_preceding_gene_splice_acceptor_position(&self, trans_id: i32) -> i32 {
        self.trans_splice_acceptor_pos.get(&trans_id).copied().unwrap_or(-1)
    }

    /// findGeneRegions without the downstream gene annotations (which the caller adds): genes whose range,
    /// extended upstream by `upstream_distance`, holds the position, in file order
    pub fn find_gene_regions(&self, chromosome: &str, position: i32, upstream_distance: i32) -> Vec<&GeneData> {
        let Some(list) = self.chr_genes.get(&chromosome.to_string()) else { return Vec::new() };
        list.iter()
            .filter(|g| {
                let start = if g.strand == 1 { g.gene_start - upstream_distance } else { g.gene_start };
                let end = if g.strand == 1 { g.gene_end } else { g.gene_end + upstream_distance };
                position >= start && position <= end
            })
            .collect()
    }

    /// findGenesByRegion: genes with a transcript lying wholly within the region
    pub fn find_genes_by_region(&self, chromosome: &str, pos_start: i32, pos_end: i32) -> Vec<&GeneData> {
        let mut out = Vec::new();
        // Java dereferences the chromosome's list unchecked
        let list = self.chr_genes.get(&chromosome.to_string()).expect("Ensembl genes of chromosome");
        for g in list {
            if pos_start > g.gene_end || pos_end < g.gene_start {
                continue;
            }
            let Some(trans) = self.transcripts_by_gene.get(&g.gene_id) else { continue };
            if trans.iter().any(|t| pos_start <= t.trans_start && pos_end >= t.trans_end) {
                out.push(g);
            }
        }
        out
    }

    /// loadTranscriptData(restrictedGeneIds, nonCanonicalTrans)
    pub fn load_transcript_data(&mut self, restricted_gene_ids: &[String], non_canonical_trans: &[String]) -> Result<()> {
        load_transcript_data(
            &self.data_path,
            self.version,
            &mut self.transcripts_by_gene,
            restricted_gene_ids,
            self.require_exons,
            self.canonical_only,
            self.require_non_ensembl_transcripts,
            non_canonical_trans,
        )
    }

    /// getChrGeneDataMap
    pub fn chr_gene_data_map(&self) -> &JavaHashMap<String, Vec<GeneData>> {
        &self.chr_genes
    }
    /// getTranscriptDataMap
    pub fn transcript_data_map(&self) -> &JavaHashMap<String, Vec<TranscriptData>> {
        &self.transcripts_by_gene
    }
    /// the maps of a cache built in memory (GeneTestUtils.addGeneData, addTransExonData)
    pub fn chr_gene_data_map_mut(&mut self) -> &mut JavaHashMap<String, Vec<GeneData>> {
        &mut self.chr_genes
    }
    pub fn transcript_data_map_mut(&mut self) -> &mut JavaHashMap<String, Vec<TranscriptData>> {
        &mut self.transcripts_by_gene
    }
    pub fn transcript_protein_data_map_mut(&mut self) -> &mut JavaHashMap<i32, Vec<TranscriptProteinData>> {
        &mut self.protein_data
    }
    pub fn set_trans_splice_acceptor_position(&mut self, trans_id: i32, position: i32) {
        self.trans_splice_acceptor_pos.insert(trans_id, position);
    }

    /// createGeneIdDataMap
    pub fn create_gene_id_data_map(&mut self) {
        if !self.gene_ids.is_empty() {
            return;
        }
        let mut m = JavaHashMap::new();
        for (chr, genes) in self.chr_genes.iter() {
            for (i, g) in genes.iter().enumerate() {
                m.insert(g.gene_id.clone(), (chr.clone(), i));
            }
        }
        self.gene_ids = m;
    }

    /// createGeneNameIdMap: by name, the last gene in the chromosome map's iteration order winning
    pub fn create_gene_name_id_map(&mut self) {
        if !self.gene_names.is_empty() {
            return;
        }
        let mut m = JavaHashMap::new();
        for (chr, genes) in self.chr_genes.iter() {
            for (i, g) in genes.iter().enumerate() {
                m.insert(g.gene_name.clone(), (chr.clone(), i));
            }
        }
        self.gene_names = m;
    }

    fn gene_at(&self, at: &(String, usize)) -> Option<&GeneData> {
        self.chr_genes.get(&at.0).and_then(|l| l.get(at.1))
    }

    /// the first gene in the chromosome map's iteration order whose name (or id) matches
    fn find_gene(&self, gene: &str, by_name: bool) -> Option<&GeneData> {
        for (_, genes) in self.chr_genes.iter() {
            for g in genes {
                if (by_name && g.gene_name == gene) || (!by_name && g.gene_id == gene) {
                    return Some(g);
                }
            }
        }
        None
    }

    /// getGeneDataByName: from the name map once it is made, else the first match
    pub fn gene_data_by_name(&self, name: &str) -> Option<&GeneData> {
        if !self.gene_names.is_empty() {
            return self.gene_names.get(&name.to_string()).and_then(|at| self.gene_at(at));
        }
        self.find_gene(name, true)
    }

    /// getGeneDataById
    pub fn gene_data_by_id(&self, gene_id: &str) -> Option<&GeneData> {
        if !self.gene_ids.is_empty() {
            return self.gene_ids.get(&gene_id.to_string()).and_then(|at| self.gene_at(at));
        }
        self.find_gene(gene_id, false)
    }

    /// getTranscripts(geneId)
    pub fn transcripts(&self, gene_id: &str) -> Option<&Vec<TranscriptData>> {
        self.transcripts_by_gene.get(&gene_id.to_string())
    }

    /// getTranscriptData(geneId, transcriptId): the canonical transcript for an empty id, else by name
    pub fn transcript_data(&self, gene_id: &str, trans_name: &str) -> Option<&TranscriptData> {
        let list = self.transcripts_by_gene.get(&gene_id.to_string())?;
        list.iter().find(|t| (trans_name.is_empty() && t.is_canonical) || t.trans_name == trans_name)
    }

    /// getCanonicalTranscriptData(geneId)
    pub fn canonical_transcript_data(&self, gene_id: &str) -> Option<&TranscriptData> {
        self.transcript_data(gene_id, "")
    }

    /// findGeneByRange: genes of the chromosome overlapping the range, in file order
    pub fn find_gene_by_range(&self, chromosome: &str, start: i32, end: i32) -> Vec<&GeneData> {
        match self.chr_genes.get(&chromosome.to_string()) {
            Some(l) => l.iter().filter(|g| crate::regions::positions_overlap(g.gene_start, g.gene_end, start, end)).collect(),
            None => Vec::new(),
        }
    }
}

fn index_of(f: &std::collections::HashMap<String, usize>, name: &str, file: &str) -> Result<usize> {
    // Java unboxes a missing column's null index: NullPointerException
    f.get(name).copied().ok_or_else(|| anyhow!("{file}: missing column {name}"))
}

/// EnsemblDataLoader.loadEnsemblGeneData. Java returns false for a missing file or an IOException.
pub fn load_ensembl_gene_data(
    data_path: &str, restricted_gene_ids: &[String], chr_genes: &mut JavaHashMap<String, Vec<GeneData>>, version: RefGenomeVersion,
    load_synonyms: bool,
) -> Result<()> {
    let filename = check_old_filename(&(check_add_dir_separator(data_path) + &ensembl_gene_data_file(version)));
    if !std::path::Path::new(&filename).exists() {
        bail!("missing Ensembl gene data file: {filename}");
    }
    let restricted: HashSet<&str> = restricted_gene_ids.iter().map(|s| s.as_str()).collect();
    let mut cols: Option<[usize; 8]> = None;
    let mut current_chr = String::new();
    java_text::for_each_line(&filename, |line| {
        let Some([gene_id_i, gene_name_i, chr_i, strand_i, start_i, end_i, band_i, syn_i]) = cols else {
            let f = java_text::fields_index_map(line, ',');
            cols = Some([
                index_of(&f, "GeneId", &filename)?,
                index_of(&f, "GeneName", &filename)?,
                index_of(&f, "Chromosome", &filename)?,
                index_of(&f, "Strand", &filename)?,
                index_of(&f, "GeneStart", &filename)?,
                index_of(&f, "GeneEnd", &filename)?,
                index_of(&f, "KaryotypeBand", &filename)?,
                index_of(&f, "Synonyms", &filename)?,
            ]);
            return Ok(());
        };
        let items = java_text::split_all(line, ',');
        let item = |i: usize| items.get(i).copied().ok_or_else(|| anyhow!("{filename}: short line: {line}"));
        let gene_id = item(gene_id_i)?;
        if !restricted.is_empty() && !restricted.contains(gene_id) {
            return Ok(());
        }
        let chromosome = version.versioned_chromosome(item(chr_i)?);
        let bad = |what: &str| anyhow!("{filename}: invalid {what} in line: {line}");
        let gene = GeneData {
            gene_id: gene_id.to_string(),
            gene_name: item(gene_name_i)?.to_string(),
            chromosome: chromosome.clone(),
            strand: parse_byte(item(strand_i)?).ok_or_else(|| bad("strand"))?,
            gene_start: parse_int(item(start_i)?).ok_or_else(|| bad("start"))?,
            gene_end: parse_int(item(end_i)?).ok_or_else(|| bad("end"))?,
            karyotype_band: item(band_i)?.to_string(),
            synonyms: if load_synonyms { item(syn_i)?.to_string() } else { String::new() },
        };
        if current_chr != chromosome {
            current_chr = chromosome.clone();
            if !chr_genes.contains_key(&chromosome) {
                chr_genes.insert(chromosome.clone(), Vec::new());
            }
        }
        chr_genes.get_mut(&current_chr).unwrap().push(gene);
        Ok(())
    })
    .with_context(|| format!("failed to load Ensembl gene ({filename})"))?;
    if cols.is_none() {
        bail!("empty Ensembl gene data file({filename})");
    }
    Ok(())
}

/// EnsemblDataLoader.loadTranscriptData. Rows are split as Java's `split(",")` does (trailing empty
/// fields dropped). A gene's transcripts start a new list each time the gene id changes (a gene seen
/// again later replaces its earlier list, as Java's put does); a transcript is skipped when only
/// canonical ones are wanted (unless named in `non_canonical_trans`) or when its name does not start
/// with ENST (unless non-Ensembl transcripts are wanted); exons are kept for every transcript when
/// `cache_exons`, else for the canonical ones.
#[allow(clippy::too_many_arguments)]
pub fn load_transcript_data(
    data_path: &str, version: RefGenomeVersion, transcripts: &mut JavaHashMap<String, Vec<TranscriptData>>, restricted_gene_ids: &[String],
    cache_exons: bool, canonical_only: bool, include_non_ensembl: bool, non_canonical_trans: &[String],
) -> Result<()> {
    let filename = check_old_filename(&(check_add_dir_separator(data_path) + &ensembl_trans_exon_data_file(version)));
    if !std::path::Path::new(&filename).exists() {
        bail!("missing Ensembl transcript data file: {filename}");
    }
    let Some(header) = java_text::first_line(&filename)? else { bail!("empty Ensembl gene-exon data file({filename})") };
    let f = java_text::fields_index_map(&header, ',');
    let gene_id_i = index_of(&f, "GeneId", &filename)?;
    let canonical_i = index_of(&f, "CanonicalTranscriptId", &filename)?;
    let strand_i = index_of(&f, "Strand", &filename)?;
    let trans_id_i = index_of(&f, "TransId", &filename)?;
    let trans_name_i = index_of(&f, "TransName", &filename)?;
    let biotype_i = index_of(&f, "BioType", &filename)?;
    let trans_start_i = index_of(&f, "TransStart", &filename)?;
    let trans_end_i = index_of(&f, "TransEnd", &filename)?;
    let rank_i = index_of(&f, "ExonRank", &filename)?;
    let exon_start_i = index_of(&f, "ExonStart", &filename)?;
    let exon_end_i = index_of(&f, "ExonEnd", &filename)?;
    let phase_i = index_of(&f, "ExonPhase", &filename)?;
    let end_phase_i = index_of(&f, "ExonEndPhase", &filename)?;
    let coding_start_i = index_of(&f, "CodingStart", &filename)?;
    let coding_end_i = index_of(&f, "CodingEnd", &filename)?;
    let ref_seq_i: i64 = f.get("RefSeqId").map(|i| *i as i64).unwrap_or(-1);
    let restricted: HashSet<&str> = restricted_gene_ids.iter().map(|s| s.as_str()).collect();

    let mut current_gene = String::new();
    // the transcript the next exon rows belong to: (gene id of its list, index in the list, trans id)
    let mut current: Option<(String, usize, i32)> = None;
    let mut last_skipped_gene = String::new();
    let mut is_header = true;
    java_text::for_each_line(&filename, |line| {
        if is_header {
            is_header = false;
            return Ok(());
        }
        let items = java_text::split(line, ',');
        let item = |i: usize| items.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBounds in line: {line}"));
        let bad = |what: &str| anyhow!("{filename}: invalid {what} in line: {line}");
        let gene_id = item(gene_id_i)?;
        let trans_id = parse_int(item(trans_id_i)?).ok_or_else(|| bad("TransId"))?;
        if last_skipped_gene == gene_id || (!restricted.is_empty() && !restricted.contains(gene_id)) {
            last_skipped_gene = gene_id.to_string();
            return Ok(());
        }
        if gene_id != current_gene {
            current_gene = gene_id.to_string();
            transcripts.insert(current_gene.clone(), Vec::new());
        }
        if current.as_ref().is_none_or(|c| c.2 != trans_id) {
            let canonical_id = parse_int(item(canonical_i)?).ok_or_else(|| bad("CanonicalTranscriptId"))?;
            let is_canonical = canonical_id == trans_id;
            let trans_name = item(trans_name_i)?;
            if canonical_only && !is_canonical && !non_canonical_trans.iter().any(|t| t == trans_name) {
                return Ok(());
            }
            if !include_non_ensembl && !trans_name.starts_with(ENSEMBL_TRANSCRIPT_PREFIX) {
                return Ok(());
            }
            let null_or_int = |i: usize, what: &str| -> Result<Option<i32>> {
                let v = item(i)?;
                if v.eq_ignore_ascii_case("NULL") { Ok(None) } else { Ok(Some(parse_int(v).ok_or_else(|| bad(what))?)) }
            };
            let coding_start = null_or_int(coding_start_i, "CodingStart")?;
            let coding_end = null_or_int(coding_end_i, "CodingEnd")?;
            // Java: refSeqIdIndex > 0 (an index of 0 is ignored too)
            let ref_seq_id = if ref_seq_i > 0 {
                let v = item(ref_seq_i as usize)?;
                if v.eq_ignore_ascii_case("NULL") { None } else { Some(v.to_string()) }
            } else {
                None
            };
            let t = TranscriptData {
                trans_id,
                trans_name: trans_name.to_string(),
                gene_id: gene_id.to_string(),
                is_canonical,
                strand: parse_byte(item(strand_i)?).ok_or_else(|| bad("Strand"))?,
                trans_start: parse_int(item(trans_start_i)?).ok_or_else(|| bad("TransStart"))?,
                trans_end: parse_int(item(trans_end_i)?).ok_or_else(|| bad("TransEnd"))?,
                coding_start,
                coding_end,
                bio_type: item(biotype_i)?.to_string(),
                ref_seq_id,
                exons: Vec::new(),
            };
            let list = transcripts.get_mut(&current_gene).unwrap();
            list.push(t);
            current = Some((current_gene.clone(), list.len() - 1, trans_id));
        }
        let (cg, ci, _) = current.as_ref().unwrap();
        let trans = &mut transcripts.get_mut(cg).unwrap()[*ci];
        if cache_exons || trans.is_canonical {
            trans.exons.push(ExonData {
                trans_id,
                start: parse_int(item(exon_start_i)?).ok_or_else(|| bad("ExonStart"))?,
                end: parse_int(item(exon_end_i)?).ok_or_else(|| bad("ExonEnd"))?,
                rank: parse_int(item(rank_i)?).ok_or_else(|| bad("ExonRank"))?,
                phase_start: parse_int(item(phase_i)?).ok_or_else(|| bad("ExonPhase"))?,
                phase_end: parse_int(item(end_phase_i)?).ok_or_else(|| bad("ExonEndPhase"))?,
            });
        }
        Ok(())
    })
    .with_context(|| format!("failed to load gene transcript exon data({filename})"))?;
    Ok(())
}

/// EnsemblDataLoader.loadTranscriptProteinData: protein features by transcript id, a new list each time the
/// transcript id changes
pub fn load_transcript_protein_data(
    data_path: &str, version: RefGenomeVersion, map: &mut JavaHashMap<i32, Vec<TranscriptProteinData>>, restricted_trans_ids: &HashSet<i32>,
) -> Result<()> {
    let filename = check_old_filename(&(check_add_dir_separator(data_path) + &ensembl_protein_data_file(version)));
    if !std::path::Path::new(&filename).exists() {
        bail!("missing Ensembl protein feature data file: {filename}");
    }
    let mut cols: Option<[usize; 6]> = None;
    let mut current_trans_id = -1;
    java_text::for_each_line(&filename, |line| {
        let Some([trans_i, translation_i, pf_i, start_i, end_i, desc_i]) = cols else {
            let f = java_text::fields_index_map(line, ',');
            cols = Some([
                index_of(&f, "TranscriptId", &filename)?,
                index_of(&f, "TranslationId", &filename)?,
                index_of(&f, "ProteinFeatureId", &filename)?,
                index_of(&f, "SeqStart", &filename)?,
                index_of(&f, "SeqEnd", &filename)?,
                index_of(&f, "HitDescription", &filename)?,
            ]);
            return Ok(());
        };
        let items = java_text::split(line, ',');
        let item = |i: usize| items.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBounds in line: {line}"));
        let bad = |what: &str| anyhow!("{filename}: invalid {what} in line: {line}");
        let trans_id = parse_int(item(trans_i)?).ok_or_else(|| bad("TranscriptId"))?;
        if !restricted_trans_ids.is_empty() && !restricted_trans_ids.contains(&trans_id) {
            return Ok(());
        }
        if trans_id != current_trans_id {
            current_trans_id = trans_id;
            map.insert(trans_id, Vec::new());
        }
        let d = TranscriptProteinData {
            transcript_id: trans_id,
            translation_id: parse_int(item(translation_i)?).ok_or_else(|| bad("TranslationId"))?,
            protein_feature_id: parse_int(item(pf_i)?).ok_or_else(|| bad("ProteinFeatureId"))?,
            seq_start: parse_int(item(start_i)?).ok_or_else(|| bad("SeqStart"))?,
            seq_end: parse_int(item(end_i)?).ok_or_else(|| bad("SeqEnd"))?,
            hit_description: item(desc_i)?.to_string(),
        };
        map.get_mut(&trans_id).unwrap().push(d);
        Ok(())
    })
    .with_context(|| format!("failed to load protein features({filename})"))?;
    if cols.is_none() {
        bail!("empty Ensembl protein feature data file({filename})");
    }
    Ok(())
}

/// EnsemblDataLoader.loadTranscriptSpliceAcceptorData: PreSpliceAcceptorPosition by TransId
pub fn load_transcript_splice_acceptor_data(
    data_path: &str, version: RefGenomeVersion, map: &mut std::collections::HashMap<i32, i32>, restricted_trans_ids: &HashSet<i32>,
) -> Result<()> {
    let filename = check_old_filename(&(check_add_dir_separator(data_path) + &ensembl_splice_data_file(version)));
    if !std::path::Path::new(&filename).exists() {
        bail!("missing Ensembl splice acceptor data file: {filename}");
    }
    let mut cols: Option<[usize; 2]> = None;
    java_text::for_each_line(&filename, |line| {
        let Some([trans_i, sa_i]) = cols else {
            let f = java_text::fields_index_map(line, ',');
            cols = Some([index_of(&f, "TransId", &filename)?, index_of(&f, "PreSpliceAcceptorPosition", &filename)?]);
            return Ok(());
        };
        let items = java_text::split(line, ',');
        let item = |i: usize| items.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBounds in line: {line}"));
        let bad = |what: &str| anyhow!("{filename}: invalid {what} in line: {line}");
        let trans_id = parse_int(item(trans_i)?).ok_or_else(|| bad("TransId"))?;
        if !restricted_trans_ids.is_empty() && !restricted_trans_ids.contains(&trans_id) {
            return Ok(());
        }
        let sa = parse_int(item(sa_i)?).ok_or_else(|| bad("PreSpliceAcceptorPosition"))?;
        map.insert(trans_id, sa);
        Ok(())
    })
    .with_context(|| format!("failed to load transcript splice data({filename})"))?;
    if cols.is_none() {
        bail!("empty Ensembl trans splice acceptor data file({filename})");
    }
    Ok(())
}

/// TranscriptAminoAcids (NEO)
#[derive(Clone, Debug, PartialEq)]
pub struct TranscriptAminoAcids {
    pub gene_id: String,
    pub gene_name: String,
    pub trans_name: String,
    pub canonical: bool,
    pub amino_acids: String,
}

pub fn ensembl_amino_acid_data_file(version: RefGenomeVersion) -> String {
    format!("ensembl_trans_amino_acids.{}.csv", version.identifier())
}

/// EnsemblDataLoader.loadTranscriptAminoAcidData: TranscriptAminoAcids by transcript name (the last line
/// of a name wins). Err where Java returns false (no file) or throws. As Java, a line of a gene outside
/// `restricted_gene_ids` also skips the line after it (Java calls readLine() once more).
pub fn load_transcript_amino_acid_data(
    data_path: &str, version: RefGenomeVersion, restricted_gene_ids: &[String], canonical_only: bool,
) -> Result<std::collections::HashMap<String, TranscriptAminoAcids>> {
    let filename = check_old_filename(&(check_add_dir_separator(data_path) + &ensembl_amino_acid_data_file(version)));
    if !std::path::Path::new(&filename).exists() {
        bail!("missing Ensembl transcript amino acid data file: {filename}");
    }
    let mut map = std::collections::HashMap::new();
    let mut cols: Option<[usize; 4]> = None;
    let mut canonical_col: Option<usize> = None;
    let mut skip_next = false;
    java_text::for_each_line(&filename, |line| {
        let Some([gene_i, name_i, trans_i, aa_i]) = cols else {
            let f = java_text::fields_index_map(line, ',');
            cols = Some([
                index_of(&f, "GeneId", &filename)?,
                index_of(&f, "GeneName", &filename)?,
                index_of(&f, "TransName", &filename)?,
                index_of(&f, "AminoAcids", &filename)?,
            ]);
            canonical_col = f.get("Canonical").copied();
            return Ok(());
        };
        if skip_next {
            skip_next = false;
            return Ok(());
        }
        // line.split(ENSEMBL_DELIM): trailing empty fields dropped
        let items = java_text::split(line, ',');
        let item = |i: usize| items.get(i).copied().ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBounds in line: {line}"));
        let gene_id = item(gene_i)?;
        if !restricted_gene_ids.is_empty() && !restricted_gene_ids.iter().any(|g| g == gene_id) {
            skip_next = true;
            return Ok(());
        }
        let trans_name = item(trans_i)?;
        let canonical = match canonical_col {
            None => true,
            Some(i) => java_text::parse_bool(item(i)?),
        };
        if canonical_only && !canonical {
            return Ok(());
        }
        map.insert(
            trans_name.to_string(),
            TranscriptAminoAcids {
                gene_id: gene_id.to_string(),
                gene_name: item(name_i)?.to_string(),
                trans_name: trans_name.to_string(),
                canonical,
                amino_acids: item(aa_i)?.to_string(),
            },
        );
        Ok(())
    })
    .with_context(|| format!("failed to load transcript amino-acid data({filename})"))?;
    Ok(map)
}

/// EnsemblDataCache.getProteinDomainPositions: the genomic start and end of a protein feature
pub fn get_protein_domain_positions(protein: &TranscriptProteinData, trans: &TranscriptData) -> [Option<i32>; 2] {
    if trans.exons.is_empty() {
        return [None, None];
    }
    let (Some(coding_start), Some(coding_end)) = (trans.coding_start, trans.coding_end) else { return [None, None] };
    let mut pre_protein_bases = protein.seq_start * 3;
    let mut protein_bases = (protein.seq_end - protein.seq_start) * 3;
    let mut protein_start = -1;
    let mut protein_end = -1;
    if trans.strand == 1 {
        for exon in &trans.exons {
            if exon.end < coding_start {
                continue;
            }
            if pre_protein_bases > 0 {
                let ref_start = coding_start.max(exon.start);
                let exon_coding_bases = exon.end - ref_start;
                if exon_coding_bases >= pre_protein_bases {
                    protein_start = ref_start + pre_protein_bases;
                    pre_protein_bases = 0;
                } else {
                    pre_protein_bases -= exon_coding_bases;
                    continue;
                }
            }
            let start_pos = exon.start.max(protein_start);
            let exon_bases = exon.end - start_pos;
            if exon_bases >= protein_bases {
                protein_end = start_pos + protein_bases;
                break;
            } else {
                protein_bases -= exon_bases;
            }
        }
    } else {
        for exon in trans.exons.iter().rev() {
            if exon.start > coding_end {
                continue;
            }
            if pre_protein_bases > 0 {
                let ref_start = coding_end.min(exon.end);
                let exon_coding_bases = ref_start - exon.start;
                if exon_coding_bases >= pre_protein_bases {
                    protein_end = ref_start - pre_protein_bases;
                    pre_protein_bases = 0;
                } else {
                    pre_protein_bases -= exon_coding_bases;
                    continue;
                }
            }
            let start_pos = exon.end.min(protein_end);
            let exon_bases = start_pos - exon.start;
            if exon_bases >= protein_bases {
                protein_start = start_pos - protein_bases;
                break;
            } else {
                protein_bases -= exon_bases;
            }
        }
    }
    if protein_end == -1 || protein_start == -1 {
        return [None, None];
    }
    [Some(protein_start), Some(protein_end)]
}

/// EnsemblDataCache.getExonRankings(strand, exons, position): the exon ranks before and after a position
pub fn get_exon_rankings(strand: i8, exons: &[ExonData], position: i32) -> [i32; 2] {
    let mut d = [0, 0];
    let first = &exons[0];
    let last = &exons[exons.len() - 1];
    if (position < first.start && strand == 1) || (position > last.end && strand == -1) {
        d = [0, 1];
    } else if (position < first.start && strand == -1) || (position > last.end && strand == 1) {
        d = [exons.len() as i32, -1];
    } else {
        for i in 0..exons.len() {
            let e = &exons[i];
            let next = exons.get(i + 1);
            if position == e.end || position == e.start || (position >= e.start && position <= e.end) {
                d = [e.rank, e.rank];
                break;
            }
            if let Some(n) = next {
                if position > e.end && position < n.start {
                    d = if strand == 1 { [e.rank, n.rank] } else { [n.rank, e.rank] };
                    break;
                }
            }
        }
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    const GENES: &str = "GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,KaryotypeBand,Synonyms\n\
        G1,GENEA,1,1,100,1000,p1,SYN1;SYN2\n\
        G2,GENEB,1,-1,2000,3000,p2,\n\
        G3,GENEC,X,1,500,800,q1,\n\
        G4,GENEA,2,1,10,20,q2,\n";
    const TRANS: &str = "GeneId,CanonicalTranscriptId,Strand,TransId,TransName,BioType,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd,RefSeqId\n\
        G1,11,1,10,ENST10,protein_coding,100,900,1,100,200,-1,0,150,850,NULL\n\
        G1,11,1,10,ENST10,protein_coding,100,900,2,800,900,0,-1,150,850,NULL\n\
        G1,11,1,11,ENST11,protein_coding,100,1000,1,100,300,-1,0,150,950,NM_1\n\
        G1,11,1,11,ENST11,protein_coding,100,1000,2,900,1000,0,-1,150,950,NM_1\n\
        G2,21,-1,21,ENST21,lincRNA,2000,3000,1,2500,3000,-1,-1,NULL,NULL,NULL\n\
        G2,21,-1,22,XR_22,lincRNA,2000,3000,1,2000,2100,-1,-1,NULL,NULL,NULL\n\
        G3,31,1,31,ENST31,protein_coding,500,800,1,500,800,-1,-1,500,800,NULL\n";

    fn write_cache(name: &str, versioned: bool) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("hmf_common_ensembl_{name}_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let suffix = if versioned { ".38" } else { "" };
        std::fs::write(dir.join(format!("ensembl_gene_data{suffix}.csv")), GENES).unwrap();
        std::fs::write(dir.join(format!("ensembl_trans_exon_data{suffix}.csv")), TRANS).unwrap();
        dir
    }

    #[test]
    fn loads_genes_transcripts_and_exons() {
        let dir = write_cache("all", false);
        let mut cache = EnsemblDataCache::new(dir.to_str().unwrap(), RefGenomeVersion::V38);
        cache.load(false).unwrap();
        let chr1 = cache.chr_gene_data_map().get(&"chr1".to_string()).unwrap();
        assert_eq!(chr1.iter().map(|g| g.gene_id.as_str()).collect::<Vec<_>>(), vec!["G1", "G2"]);
        assert_eq!(chr1[1].strand, -1);
        assert_eq!(chr1[0].synonyms, "");
        let t = cache.transcripts("G1").unwrap();
        assert_eq!(t.len(), 2);
        assert!(!t[0].is_canonical && t[1].is_canonical);
        assert_eq!(t[1].exons.len(), 2);
        assert_eq!(t[1].ref_seq_id.as_deref(), Some("NM_1"));
        assert_eq!(cache.canonical_transcript_data("G1").unwrap().trans_name, "ENST11");
        // non-Ensembl transcripts are left out
        assert_eq!(cache.transcripts("G2").unwrap().len(), 1);
        assert!(cache.transcripts("G2").unwrap()[0].non_coding());
        // two genes named GENEA: the name map keeps the last one in HashMap order of the chromosomes
        cache.create_gene_name_id_map();
        let a = cache.gene_data_by_name("GENEA").unwrap();
        let order: Vec<&String> = cache.chr_gene_data_map().keys().collect();
        let last_chr = order.iter().rev().find(|c| ***c == "chr1" || ***c == "chr2").unwrap();
        assert_eq!(&&a.chromosome, last_chr);
        assert!(cache.gene_data_by_name("NONE").is_none());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn canonical_only_and_restricted() {
        let dir = write_cache("canon", true);
        let mut cache = EnsemblDataCache::new(dir.to_str().unwrap(), RefGenomeVersion::V38);
        cache.set_required_data(true, false, false, true);
        cache.set_require_gene_synonyms();
        cache.load(false).unwrap();
        let t = cache.transcripts("G1").unwrap();
        assert_eq!(t.len(), 1);
        assert_eq!(t[0].trans_id, 11);
        assert_eq!(cache.gene_data_by_id("G1").unwrap().synonyms, "SYN1;SYN2");
        let mut r = EnsemblDataCache::new(dir.to_str().unwrap(), RefGenomeVersion::V37);
        r.set_restricted_gene_id_list(&["G3".to_string()]);
        // no 37 files, nor unversioned ones
        assert!(r.load(false).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn transcript_amino_acids() {
        let dir = std::env::temp_dir().join(format!("ensembl_aa_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("ensembl_trans_amino_acids.csv"),
            "GeneId,GeneName,TransName,Canonical,AminoAcids\nG1,N1,T1,true,MAAA\nG2,N2,T2,false,MCCC\nG3,N3,T3,true,MDDD\nG1,N1,T1,true,MEEE\n",
        )
        .unwrap();
        let d = dir.to_str().unwrap();
        // the unversioned file stands in for the 38 one; a repeated transcript keeps its last line
        let all = load_transcript_amino_acid_data(d, RefGenomeVersion::V38, &[], false).unwrap();
        assert_eq!(all.len(), 3);
        assert_eq!(all["T1"].amino_acids, "MEEE");
        assert!(!all["T2"].canonical);
        let canonical = load_transcript_amino_acid_data(d, RefGenomeVersion::V37, &[], true).unwrap();
        assert!(canonical.contains_key("T3") && !canonical.contains_key("T2"));
        // as Java, a line outside the restricted genes also skips the line after it: G2 skips T3 too
        let r = load_transcript_amino_acid_data(d, RefGenomeVersion::V38, &["G1".to_string(), "G3".to_string()], false).unwrap();
        assert!(r.contains_key("T1") && !r.contains_key("T2") && !r.contains_key("T3"));
        assert!(load_transcript_amino_acid_data("/nonexistent", RefGenomeVersion::V38, &[], false).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn old_file_names() {
        // the regex "38.csv": any character between the digits and "csv"
        assert_eq!(replace_version_pattern("/a/ensembl_gene_data.38.csv", "38"), "/a/ensembl_gene_data.csv");
        assert_eq!(replace_version_pattern("x38Ycsv.38.csv", "38"), "xcsv.csv");
        assert_eq!(replace_version_pattern("/a/ensembl_gene_data.38.csv", "37"), "/a/ensembl_gene_data.38.csv");
    }
}
