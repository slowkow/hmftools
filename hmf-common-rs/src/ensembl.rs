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
        if self.require_protein_domains || self.require_splice_positions {
            bail!("Ensembl protein domain and splice position data are not ported");
        }
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
        }
        Ok(())
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
    fn old_file_names() {
        // the regex "38.csv": any character between the digits and "csv"
        assert_eq!(replace_version_pattern("/a/ensembl_gene_data.38.csv", "38"), "/a/ensembl_gene_data.csv");
        assert_eq!(replace_version_pattern("x38Ycsv.38.csv", "38"), "xcsv.csv");
        assert_eq!(replace_version_pattern("/a/ensembl_gene_data.38.csv", "37"), "/a/ensembl_gene_data.38.csv");
    }
}
