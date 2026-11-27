//! Alignment annotation: AlignmentUtil (BWA-MEM through hmf-common-rs's bwa module, with CIDER's
//! options), ImgtSequenceFile (the IMGT gene alleles, compiled in), AlignmentAnnotator and the somatic
//! hypermutation comparison of the V and J regions with IMGT (SomaticHypermutation.kt).

use std::collections::HashMap;

use anyhow::{anyhow, bail, Result};
use hmf_common::bwa::{BwaIndex, BwaOptions, MEM_F_ALL};
use hmf_common::codon::reverse_complement;
use hmf_common::regions::RefGenomeVersion;

use crate::constants::*;
use crate::genes::{read_gene_file, IgTcrGene, IgTcrLocus, IgTcrRegion, Strand, Vj};
use crate::read::{cigar_from_str, Cigar, CigarOp};
use crate::vdj::VdjSequence;

/// Alignment: positions 0-based, ranges inclusive (Kotlin IntRange)
#[derive(Clone, Debug)]
pub struct Alignment {
    pub query_seq: String,
    pub query_start: i32,
    pub query_end: i32,
    pub ref_contig: String,
    pub ref_start: i32,
    pub ref_end: i32,
    pub strand: Strand,
    pub alignment_score: i32,
    /// with respect to the query, clipping included
    pub edit_distance: i32,
    pub cigar: Cigar,
    pub ref_contig_length: i32,
}

impl Alignment {
    pub fn query_align_length(&self) -> i32 {
        self.query_end + 1 - self.query_start
    }
    pub fn start_clip(&self) -> i32 {
        self.query_start
    }
    pub fn end_clip(&self) -> i32 {
        self.query_seq.len() as i32 - (self.query_end + 1)
    }
    pub fn aligned_edit_distance(&self) -> i32 {
        self.edit_distance - self.start_clip() - self.end_clip()
    }
}

/// a sequence dictionary: contig names and lengths in order (ReferenceSequenceFileFactory.loadDictionary)
pub fn parse_dict(text: &str) -> Vec<(String, i32)> {
    let mut out = Vec::new();
    for line in text.lines() {
        if !line.starts_with("@SQ") {
            continue;
        }
        let mut name = None;
        let mut len = 0;
        for f in line.split('\t') {
            if let Some(v) = f.strip_prefix("SN:") {
                name = Some(v.to_string());
            } else if let Some(v) = f.strip_prefix("LN:") {
                len = v.parse().unwrap_or(0);
            }
        }
        if let Some(n) = name {
            out.push((n, len));
        }
    }
    out
}

/// createBwaMemAligner: mem_opt_init with CIDER's changes
fn cider_bwa_options(align_score_threshold: i32, threads: usize) -> Result<BwaOptions> {
    let mut o = BwaOptions::new()?;
    let m = o.opts_mut();
    m.n_threads = threads.max(1) as i32;
    m.flag |= MEM_F_ALL;
    m.t = align_score_threshold;
    m.min_seed_len = ALIGNER_WORD_SIZE;
    m.a = ALIGNER_MATCH_SCORE;
    m.b = -ALIGNER_MISMATCH_SCORE;
    m.o_del = -ALIGNER_GAP_OPENING_SCORE;
    m.o_ins = -ALIGNER_GAP_OPENING_SCORE;
    m.e_del = -ALIGNER_GAP_EXTEND_SCORE;
    m.e_ins = -ALIGNER_GAP_EXTEND_SCORE;
    m.drop_ratio = 0.1;
    m.split_factor = 0.5;
    m.max_mem_intv = 500;
    m.max_occ = 2000;
    m.max_chain_gap = 300;
    m.zdrop = 300;
    Ok(o)
}

/// runBwaMem: the alignments of each sequence, in batches of 20,000
pub fn run_bwa_mem(sequences: &[String], dict: &[(String, i32)], index: &BwaIndex, align_score_threshold: i32, threads: usize) -> Result<Vec<Vec<Alignment>>> {
    crate::debug!("Aligning {} sequences with BWA-MEM", sequences.len());
    let opts = cider_bwa_options(align_score_threshold, threads)?;
    let mut results = Vec::with_capacity(sequences.len());
    for batch in sequences.chunks(BWAMEM_BATCH_SIZE) {
        crate::debug!("Running BWA-MEM on batch of {} sequences", batch.len());
        let seqs: Vec<&[u8]> = batch.iter().map(|s| s.as_bytes()).collect();
        let alignments = index.align(&opts, &seqs);
        for (q, a) in batch.iter().zip(alignments) {
            let mut v = Vec::new();
            for x in a {
                if let Some(al) = parse_bwa_alignment(q, &x, dict)? {
                    v.push(al);
                }
            }
            results.push(v);
        }
    }
    crate::debug!("Alignment complete");
    Ok(results)
}

/// parseBwaMemAlignment: None for an unmapped record
fn parse_bwa_alignment(query: &str, a: &hmf_common::bwa::BwaMemAlignment, dict: &[(String, i32)]) -> Result<Option<Alignment>> {
    if a.sam_flag & 0x4 != 0 {
        return Ok(None);
    }
    let (ref_contig, ref_len) = dict.get(a.ref_id as usize).cloned().ok_or_else(|| anyhow!("reference index {} not in the dictionary", a.ref_id))?;
    let strand = if a.sam_flag & 0x10 == 0 { Strand::Forward } else { Strand::Reverse };
    let cigar = cigar_from_str(&a.cigar);
    if cigar.is_empty() {
        bail!("java.lang.IllegalArgumentException: Failed requirement (cigar)");
    }
    let left = if cigar[0].0.is_clipping() { cigar[0].1 as i32 } else { 0 };
    let right = if cigar[cigar.len() - 1].0.is_clipping() { cigar[cigar.len() - 1].1 as i32 } else { 0 };
    let (start_clip, end_clip) = if strand == Strand::Forward { (left, right) } else { (right, left) };
    let al = Alignment {
        query_seq: query.to_string(),
        query_start: start_clip,
        query_end: query.len() as i32 - end_clip - 1,
        ref_contig,
        ref_start: a.ref_start,
        ref_end: a.ref_end - 1,
        strand,
        alignment_score: a.aligner_score,
        edit_distance: a.n_mismatches + start_clip + end_clip,
        cigar,
        ref_contig_length: ref_len,
    };
    // Alignment's requirements
    let ok = al.query_start >= 0
        && al.query_end < query.len() as i32
        && al.query_start <= al.query_end
        && al.ref_start >= 0
        && al.ref_start <= al.ref_end
        && al.ref_end < al.ref_contig_length
        && al.edit_distance >= al.start_clip() + al.end_clip();
    if !ok {
        bail!("java.lang.IllegalArgumentException: Failed requirement (alignment {} {}:{}-{} {})", query, al.ref_contig, al.ref_start, al.ref_end, a.cigar);
    }
    Ok(Some(al))
}

/// ImgtSequenceFile.Sequence
#[derive(Clone, Debug)]
pub struct ImgtSequence {
    pub gene_name: String,
    pub allele: String,
    pub sequence_with_ref: Vec<u8>,
    pub ref_before: i32,
    pub ref_after: i32,
}

impl ImgtSequence {
    pub fn gene_allele(&self) -> String {
        format!("{}*{}", self.gene_name, self.allele)
    }
    /// imgtRange: (first, last) inclusive
    pub fn imgt_range(&self) -> (i32, i32) {
        (self.ref_before, self.sequence_with_ref.len() as i32 - self.ref_after - 1)
    }
    fn from_fasta(label: &str, seq: Vec<u8>) -> Result<ImgtSequence> {
        let parts: Vec<&str> = label.split('|').collect();
        if parts.len() < 4 {
            bail!("invalid IMGT sequence label: {label}");
        }
        let s = ImgtSequence { gene_name: parts[0].into(), allele: parts[1].into(), sequence_with_ref: seq, ref_before: parts[2].parse()?, ref_after: parts[3].parse()? };
        if s.ref_before + s.ref_after >= s.sequence_with_ref.len() as i32 {
            bail!("java.lang.IllegalArgumentException: Failed requirement (IMGT sequence {label})");
        }
        Ok(s)
    }
}

const IMGT_FASTA_38: &str = include_str!("../resources/igtcr_gene.38.fasta");
const IMGT_FASTA_37: &str = include_str!("../resources/igtcr_gene.37.fasta");
const IMGT_DICT_38: &str = include_str!("../resources/igtcr_gene.38.fasta.dict");
const IMGT_DICT_37: &str = include_str!("../resources/igtcr_gene.37.fasta.dict");
static IMGT_IMG_38: &[u8] = include_bytes!("../resources/igtcr_gene.38.fasta.img");
static IMGT_IMG_37: &[u8] = include_bytes!("../resources/igtcr_gene.37.fasta.img");
const PATCH_DICT_37: &str = include_str!("../resources/7_gl582971_fix.fasta.dict");
static PATCH_IMG_37: &[u8] = include_bytes!("../resources/7_gl582971_fix.fasta.img");

/// ImgtSequenceFile: the curated FASTA of V/D/J gene alleles (with some reference context)
pub struct ImgtSequenceFile {
    pub dict: Vec<(String, i32)>,
    pub sequences_by_contig: HashMap<String, ImgtSequence>,
    pub index: BwaIndex,
}

impl ImgtSequenceFile {
    pub fn load(version: RefGenomeVersion) -> Result<ImgtSequenceFile> {
        let (fasta, dict, img) = if version.is_37() { (IMGT_FASTA_37, IMGT_DICT_37, IMGT_IMG_37) } else { (IMGT_FASTA_38, IMGT_DICT_38, IMGT_IMG_38) };
        let dict = parse_dict(dict);
        let mut by_name: HashMap<String, Vec<u8>> = HashMap::new();
        let mut name: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();
        for line in fasta.lines() {
            if let Some(n) = line.strip_prefix('>') {
                if let Some(prev) = name.take() {
                    by_name.insert(prev, std::mem::take(&mut seq));
                }
                // the name ends at the first white space
                name = Some(n.split_whitespace().next().unwrap_or("").to_string());
            } else {
                seq.extend_from_slice(line.trim_end().as_bytes());
            }
        }
        if let Some(prev) = name {
            by_name.insert(prev, seq);
        }
        let mut sequences_by_contig = HashMap::new();
        for (contig, _) in &dict {
            let s = by_name.get(contig).ok_or_else(|| anyhow!("IMGT sequence {contig} missing from the FASTA"))?;
            sequences_by_contig.insert(contig.clone(), ImgtSequence::from_fasta(contig, s.clone())?);
        }
        Ok(ImgtSequenceFile { dict, sequences_by_contig, index: BwaIndex::from_image_bytes(img)? })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AlignmentStatus {
    VDJ,
    VJ,
    VD,
    DJ,
    VOnly,
    DOnly,
    JOnly,
    NoRearrangement,
    NoVdjAlignment,
    SkippedAlign,
}

impl AlignmentStatus {
    pub fn name(self) -> &'static str {
        match self {
            AlignmentStatus::VDJ => "V_D_J",
            AlignmentStatus::VJ => "V_J",
            AlignmentStatus::VD => "V_D",
            AlignmentStatus::DJ => "D_J",
            AlignmentStatus::VOnly => "V_ONLY",
            AlignmentStatus::DOnly => "D_ONLY",
            AlignmentStatus::JOnly => "J_ONLY",
            AlignmentStatus::NoRearrangement => "NO_REARRANGEMENT",
            AlignmentStatus::NoVdjAlignment => "NO_VDJ_ALIGNMENT",
            AlignmentStatus::SkippedAlign => "SKIPPED_ALIGN",
        }
    }
}

/// ShmGeneComparison
#[derive(Clone, Debug, PartialEq)]
pub struct ShmGeneComparison {
    pub seq_length: i32,
    pub imgt_length: i32,
    pub pct_identity: f64,
    pub indel_bases: i32,
    pub clip_bases: i32,
    pub imgt_ref_unaligned: i32,
}

/// SomaticHypermutationStatus.determineFromVGeneIdentity
pub fn shm_status(pct_identity: f64) -> &'static str {
    if pct_identity >= 98.0 {
        "UNMUTATED"
    } else if pct_identity >= 97.0 {
        "MUTATED_BORDERLINE"
    } else {
        "MUTATED"
    }
}

#[derive(Clone, Debug)]
pub struct GeneAnnotation {
    pub gene: IgTcrGene,
    pub alignment: Alignment,
    pub comparison: Option<ShmGeneComparison>,
    pub supplementary_genes: Vec<IgTcrGene>,
}

/// AlignmentAnnotation; `vdj` is the index of the VDJ sequence it annotates
#[derive(Clone, Debug)]
pub struct AlignmentAnnotation {
    pub vdj: usize,
    /// (first, last) inclusive, in the layout
    pub query_range: (i32, i32),
    pub status: AlignmentStatus,
    pub full_alignment: Option<Alignment>,
    pub v_gene: Option<GeneAnnotation>,
    pub d_gene: Option<GeneAnnotation>,
    pub j_gene: Option<GeneAnnotation>,
}

/// compareVJRegionToImgt: the identity of the V (J) region of the layout to the IMGT allele, outside the
/// anchor side of the boundary
pub fn compare_vj_region_to_imgt(
    layout_seq: &[u8],
    vj: Vj,
    layout_anchor_boundary: i32,
    imgt: &ImgtSequence,
    query_range: (i32, i32),
    alignment: &Alignment,
) -> Option<ShmGeneComparison> {
    let is_v = vj == Vj::V;
    let forward = alignment.strand == Strand::Forward;
    let layout_start = query_range.0;
    let layout_bounds = if is_v { (0, layout_anchor_boundary) } else { (layout_anchor_boundary, layout_seq.len() as i32) };
    let imgt_seq: Vec<u8> = if forward { imgt.sequence_with_ref.clone() } else { reverse_complement(&imgt.sequence_with_ref) };
    let n = imgt_seq.len() as i32;
    let imgt_align_start = if forward { alignment.ref_start } else { n - (alignment.ref_end + 1) };
    let imgt_align_end = if forward { alignment.ref_end + 1 } else { n - alignment.ref_start };
    let r = imgt.imgt_range();
    // (start, end exclusive)
    let imgt_bounds = if forward { (r.0, r.1 + 1) } else { (n - (r.1 + 1), n - r.0) };
    let cigar: Vec<(CigarOp, u32)> = if forward { alignment.cigar.clone() } else { alignment.cigar.iter().rev().copied().collect() };
    let in_range = |b: (i32, i32), x: i32| x >= b.0 && x < b.1;
    let (mut compared, mut matches, mut indels) = (0, 0, 0);
    let mut layout_index = layout_start;
    let mut imgt_index = imgt_align_start;
    let mut first_aligned: Option<i32> = None;
    let mut first_compared_layout: Option<i32> = None;
    let mut first_compared_imgt: Option<i32> = None;
    let (mut last_aligned, mut last_compared_layout, mut last_compared_imgt) = (0, 0, 0);
    for &(op, len) in &cigar {
        let len = len as i32;
        if op.is_alignment() {
            first_aligned.get_or_insert(layout_index);
            last_aligned = layout_index + len - 1;
            for _ in 0..len {
                if in_range(layout_bounds, layout_index) && in_range(imgt_bounds, imgt_index) {
                    compared += 1;
                    if layout_seq.get(layout_index as usize) == imgt_seq.get(imgt_index as usize) {
                        matches += 1;
                    }
                    first_compared_layout.get_or_insert(layout_index);
                    first_compared_imgt.get_or_insert(imgt_index);
                    last_compared_layout = layout_index;
                    last_compared_imgt = imgt_index;
                }
                layout_index += 1;
                imgt_index += 1;
            }
        } else if op.consumes_read_bases() || op.is_clipping() {
            layout_index += len;
        } else if op.consumes_reference_bases() {
            imgt_index += len;
        }
        if op.is_indel() {
            indels += len;
        }
    }
    let layout_end = query_range.1 + 1;
    let imgt_ref_unaligned = (if is_v { imgt_align_start - imgt_bounds.0 } else { imgt_bounds.1 - imgt_align_end }).max(0);
    if compared == 0 {
        return None;
    }
    let clip = if is_v {
        (first_aligned.unwrap() - layout_start).min(imgt_align_start - imgt_bounds.0)
    } else {
        (layout_end - last_aligned).min(imgt_bounds.1 - imgt_align_end)
    };
    let c = ShmGeneComparison {
        seq_length: last_compared_layout - first_compared_layout.unwrap() + 1,
        imgt_length: last_compared_imgt - first_compared_imgt.unwrap() + 1,
        pct_identity: 100.0 * matches as f64 / compared as f64,
        indel_bases: indels,
        clip_bases: clip.max(0),
        imgt_ref_unaligned,
    };
    assert!(c.seq_length >= 0 && c.imgt_length >= 0 && (0.0..=100.0).contains(&c.pct_identity), "java.lang.IllegalArgumentException: Failed requirement (ShmGeneComparison)");
    Some(c)
}

/// AlignmentAnnotator
pub struct AlignmentAnnotator {
    version: RefGenomeVersion,
    ref_dict: Vec<(String, i32)>,
    ref_index_path: String,
    vdj_genes: HashMap<String, IgTcrGene>,
    imgt: ImgtSequenceFile,
}

struct GeneMatchCandidate<'a> {
    alignment: &'a Alignment,
    imgt: &'a ImgtSequence,
    gene: &'a IgTcrGene,
}

impl AlignmentAnnotator {
    pub fn new(version: RefGenomeVersion, ref_dict_path: &str, ref_index_path: &str) -> Result<AlignmentAnnotator> {
        let mut vdj_genes = HashMap::new();
        for g in read_gene_file(version)?.into_iter().filter(|g| g.region.is_vdj()) {
            let key = g.gene_allele();
            if vdj_genes.contains_key(&key) {
                bail!("Duplicate gene allele: {key}");
            }
            vdj_genes.insert(key, g);
        }
        let dict_text = std::fs::read_to_string(ref_dict_path).map_err(|e| anyhow!("reading {ref_dict_path}: {e}"))?;
        Ok(AlignmentAnnotator { version, ref_dict: parse_dict(&dict_text), ref_index_path: ref_index_path.to_string(), vdj_genes, imgt: ImgtSequenceFile::load(version)? })
    }

    /// alignmentQuerySeqRange: the VDJ sequence and 300 bases either side, within the layout
    pub fn query_seq_range(vdj: &VdjSequence) -> Result<(i32, i32)> {
        let full = vdj.layout.consensus_sequence();
        let seq = vdj.sequence();
        let i = find_subsequence(full, seq).ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: Failed requirement (VDJ in layout)"))? as i32;
        Ok(((i - ANNOTATION_VDJ_FLANK_BASES).max(0), (i + seq.len() as i32 + ANNOTATION_VDJ_FLANK_BASES).min(full.len() as i32) - 1))
    }

    /// runAnnotate: aligns the sequences to the genome (and for GRCh37 its TRBJ1 patch) and to the IMGT
    /// alleles; the annotations in input order
    pub fn run_annotate(&self, vdjs: &[(usize, &VdjSequence)], threads: usize) -> Result<Vec<AlignmentAnnotation>> {
        crate::info!("Running alignment annotation");
        let mut meta = Vec::with_capacity(vdjs.len());
        for &(i, v) in vdjs {
            let r = Self::query_seq_range(v)?;
            let full = v.layout.consensus_sequence();
            let q = String::from_utf8_lossy(&full[r.0 as usize..(r.1 + 1) as usize]).into_owned();
            meta.push((i, v, r, q));
        }
        let queries: Vec<String> = meta.iter().map(|m| m.3.clone()).collect();
        let mut ref_alignments = if queries.is_empty() {
            Vec::new()
        } else {
            // read through before the alignment (not a Java behaviour; the same output): see open_with_preload
            let index = BwaIndex::open_with_preload(&self.ref_index_path, "CIDER_IMAGE_PRELOAD", true)?;
            run_bwa_mem(&queries, &self.ref_dict, &index, ANNOTATION_ALIGN_SCORE_MIN, threads)?
        };
        if self.version.is_37() && !queries.is_empty() {
            crate::debug!("Aligning {} sequences to GRCh37 patch", queries.len());
            let patch = BwaIndex::from_image_bytes(PATCH_IMG_37)?;
            let patch_alignments = run_bwa_mem(&queries, &parse_dict(PATCH_DICT_37), &patch, ANNOTATION_ALIGN_SCORE_MIN, threads)?;
            for (a, p) in ref_alignments.iter_mut().zip(patch_alignments) {
                a.extend(p);
            }
        }
        crate::debug!("Aligning {} sequences to IMGT gene alleles", queries.len());
        let imgt_alignments =
            if queries.is_empty() { Vec::new() } else { run_bwa_mem(&queries, &self.imgt.dict, &self.imgt.index, ANNOTATION_ALIGN_SCORE_MIN, threads)? };
        let mut out = Vec::with_capacity(meta.len());
        for (k, (i, v, r, q)) in meta.iter().enumerate() {
            out.push(self.process(*i, v, *r, q, &ref_alignments[k], &imgt_alignments[k])?);
        }
        Ok(out)
    }

    fn process(&self, i: usize, vdj: &VdjSequence, range: (i32, i32), query: &str, ref_al: &[Alignment], imgt_al: &[Alignment]) -> Result<AlignmentAnnotation> {
        // tryMatchReference: an alignment of most of the sequence with few differences
        for a in preprocess(ref_al) {
            let distance = a.edit_distance as f64 / query.len() as f64;
            let min_align = query.len() as f64 * ANNOTATION_MATCH_REF_IDENTITY;
            if distance <= 1.0 - ANNOTATION_MATCH_REF_IDENTITY && a.query_align_length() as f64 >= min_align {
                return Ok(AlignmentAnnotation {
                    vdj: i,
                    query_range: range,
                    status: AlignmentStatus::NoRearrangement,
                    full_alignment: Some(a.clone()),
                    v_gene: None,
                    d_gene: None,
                    j_gene: None,
                });
            }
        }
        self.match_rearranged_genes(i, vdj, range, imgt_al)
    }

    fn match_rearranged_genes(&self, i: usize, vdj: &VdjSequence, range: (i32, i32), imgt_al: &[Alignment]) -> Result<AlignmentAnnotation> {
        let alignments = preprocess(imgt_al);
        let locus: IgTcrLocus = vdj
            .v_anchor
            .as_ref()
            .map(|a| a.gene_type.locus())
            .or_else(|| vdj.j_anchor.as_ref().map(|a| a.gene_type.locus()))
            .ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: Failed requirement (locus)"))?;
        // by region, in alignment order
        let mut candidates: Vec<(IgTcrRegion, Vec<GeneMatchCandidate>)> = Vec::new();
        for a in alignments {
            let imgt = self.imgt.sequences_by_contig.get(&a.ref_contig).ok_or_else(|| anyhow!("unknown IMGT contig {}", a.ref_contig))?;
            let Some(gene) = self.vdj_genes.get(&imgt.gene_allele()) else { continue };
            if gene.region.is_vj() {
                let distance = a.aligned_edit_distance() as f64 / a.query_align_length() as f64;
                if distance > 1.0 - ANNOTATION_VJ_IDENTITY_MIN {
                    continue;
                }
                if alignment_matches_anchor(vdj, range, a, gene.region == IgTcrRegion::VRegion) == Some(false) {
                    continue;
                }
            }
            if IgTcrLocus::from_gene_name(&gene.gene_name)? != locus {
                continue;
            }
            let c = GeneMatchCandidate { alignment: a, imgt, gene };
            match candidates.iter_mut().find(|(r, _)| *r == gene.region) {
                Some((_, l)) => l.push(c),
                None => candidates.push((gene.region, vec![c])),
            }
        }
        let annotation = |region: IgTcrRegion| -> Option<GeneAnnotation> {
            let list = &candidates.iter().find(|(r, _)| *r == region)?.1;
            let (primary, supplementary) = select_best_gene(list)?;
            let p = &list[primary];
            let layout_seq = vdj.layout.consensus_sequence();
            let comparison = match region {
                IgTcrRegion::VRegion => vdj
                    .v_anchor_boundary()
                    .and_then(|b| compare_vj_region_to_imgt(layout_seq, Vj::V, vdj.layout_slice_start + b, p.imgt, range, p.alignment)),
                IgTcrRegion::JRegion => vdj
                    .j_anchor_boundary()
                    .and_then(|b| compare_vj_region_to_imgt(layout_seq, Vj::J, vdj.layout_slice_start + b, p.imgt, range, p.alignment)),
                _ => None,
            };
            Some(GeneAnnotation {
                gene: p.gene.clone(),
                alignment: p.alignment.clone(),
                comparison,
                supplementary_genes: supplementary.iter().map(|&s| list[s].gene.clone()).collect(),
            })
        };
        let v_gene = annotation(IgTcrRegion::VRegion);
        let d_gene = annotation(IgTcrRegion::DRegion);
        let j_gene = annotation(IgTcrRegion::JRegion);
        let status = match (v_gene.is_some(), d_gene.is_some(), j_gene.is_some()) {
            (true, true, true) => AlignmentStatus::VDJ,
            (true, true, false) => AlignmentStatus::VD,
            (true, false, true) => AlignmentStatus::VJ,
            (true, false, false) => AlignmentStatus::VOnly,
            (false, true, true) => AlignmentStatus::DJ,
            (false, true, false) => AlignmentStatus::DOnly,
            (false, false, true) => AlignmentStatus::JOnly,
            (false, false, false) => AlignmentStatus::NoVdjAlignment,
        };
        Ok(AlignmentAnnotation { vdj: i, query_range: range, status, full_alignment: None, v_gene, d_gene, j_gene })
    }
}

fn find_subsequence(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

/// preprocessAlignments: score at least the minimum, by descending score (stable)
fn preprocess(alignments: &[Alignment]) -> Vec<&Alignment> {
    let mut v: Vec<&Alignment> = alignments.iter().filter(|a| a.alignment_score >= ANNOTATION_ALIGN_SCORE_MIN).collect();
    v.sort_by_key(|a| -a.alignment_score);
    v
}

/// alignmentMatchesAnchor: the alignment covers 5 bases of the anchor (fewer if the anchor is shorter);
/// None without an anchor
fn alignment_matches_anchor(vdj: &VdjSequence, range: (i32, i32), a: &Alignment, is_v: bool) -> Option<bool> {
    let align = (a.query_start + range.0, a.query_end + range.0);
    let anchor = if is_v {
        let boundary = vdj.v_anchor_boundary()?;
        let len = vdj.v_anchor_sequence()?.len() as i32;
        let end = boundary + vdj.layout_slice_start;
        (end - len, end - 1)
    } else {
        let boundary = vdj.j_anchor_boundary()?;
        let len = vdj.j_anchor_sequence()?.len() as i32;
        let start = boundary + vdj.layout_slice_start;
        (start, start + len - 1)
    };
    let anchor_length = anchor.1 + 1 - anchor.0;
    let required = anchor_length.min(ANNOTATION_ANCHOR_BASE_OVERLAP);
    let overlap = align.1.min(anchor.1) + 1 - align.0.max(anchor.0);
    Some(overlap >= required)
}

/// selectBestGene: the highest score, then functional, then gene allele name; the others of the same
/// score as supplementary (indices into the list)
fn select_best_gene(list: &[GeneMatchCandidate]) -> Option<(usize, Vec<usize>)> {
    if list.is_empty() {
        return None;
    }
    let base = |a: &GeneMatchCandidate, b: &GeneMatchCandidate| b.alignment.alignment_score.cmp(&a.alignment.alignment_score);
    let tie = |a: &GeneMatchCandidate, b: &GeneMatchCandidate| {
        b.gene.is_functional().cmp(&a.gene.is_functional()).then_with(|| a.gene.gene_allele().cmp(&b.gene.gene_allele()))
    };
    // minWith: the first minimum
    let mut best = 0;
    for i in 1..list.len() {
        if base(&list[i], &list[best]).then_with(|| tie(&list[i], &list[best])) == std::cmp::Ordering::Less {
            best = i;
        }
    }
    // the others with the same score; Kotlin compares the data classes (equal content counts as the primary)
    let same = |a: &GeneMatchCandidate, b: &GeneMatchCandidate| {
        std::ptr::eq(a.alignment, b.alignment) || (alignment_eq(a.alignment, b.alignment) && a.imgt.gene_allele() == b.imgt.gene_allele() && a.gene.gene_allele() == b.gene.gene_allele())
    };
    let supplementary = (0..list.len()).filter(|&i| !same(&list[i], &list[best]) && base(&list[i], &list[best]) == std::cmp::Ordering::Equal).collect();
    Some((best, supplementary))
}

fn alignment_eq(a: &Alignment, b: &Alignment) -> bool {
    a.query_seq == b.query_seq
        && a.query_start == b.query_start
        && a.query_end == b.query_end
        && a.ref_contig == b.ref_contig
        && a.ref_start == b.ref_start
        && a.ref_end == b.ref_end
        && a.strand == b.strand
        && a.alignment_score == b.alignment_score
        && a.edit_distance == b.edit_distance
        && a.cigar == b.cigar
        && a.ref_contig_length == b.ref_contig_length
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alignment(query: &str, q: (i32, i32), r: (i32, i32), strand: Strand, cigar: &str) -> Alignment {
        Alignment {
            query_seq: query.into(),
            query_start: q.0,
            query_end: q.1,
            ref_contig: "x".into(),
            ref_start: r.0,
            ref_end: r.1,
            strand,
            alignment_score: 0,
            edit_distance: 0,
            cigar: cigar_from_str(cigar),
            ref_contig_length: 1000,
        }
    }

    #[test]
    fn shm_statuses() {
        assert_eq!(shm_status(100.0), "UNMUTATED");
        assert_eq!(shm_status(98.0), "UNMUTATED");
        assert_eq!(shm_status(97.5), "MUTATED_BORDERLINE");
        assert_eq!(shm_status(96.9), "MUTATED");
    }

    // SomaticHypermutationTest: a V region with one mismatch against IMGT, forward strand
    #[test]
    fn compare_v_region_forward() {
        // IMGT: 5 ref bases, 20 IMGT bases, 5 ref bases
        let imgt_core = "ACGTACGTACGTACGTACGT";
        let imgt = ImgtSequence { gene_name: "IGHV1-1".into(), allele: "01".into(), sequence_with_ref: format!("GGGGG{imgt_core}CCCCC").into_bytes(), ref_before: 5, ref_after: 5 };
        // the layout holds the IMGT bases with one change, then the anchor boundary at 20
        let mut layout = imgt_core.as_bytes().to_vec();
        layout[3] = b'A';
        layout.extend_from_slice(b"TTTTTTTTTT");
        let a = alignment(std::str::from_utf8(&layout).unwrap(), (0, 19), (5, 24), Strand::Forward, "20M10S");
        let c = compare_vj_region_to_imgt(&layout, Vj::V, 20, &imgt, (0, 29), &a).unwrap();
        assert_eq!(c.seq_length, 20);
        assert_eq!(c.imgt_length, 20);
        assert_eq!(c.pct_identity, 95.0);
        assert_eq!((c.indel_bases, c.clip_bases, c.imgt_ref_unaligned), (0, 0, 0));
    }

    #[test]
    fn dict_parse() {
        let d = parse_dict("@HD\tVN:1.0\n@SQ\tSN:a|01|1|1\tLN:117\tM5:x\n@SQ\tSN:b\tLN:5\n");
        assert_eq!(d, vec![("a|01|1|1".to_string(), 117), ("b".to_string(), 5)]);
    }

    #[test]
    fn imgt_resources_load() {
        for v in [RefGenomeVersion::V37, RefGenomeVersion::V38] {
            let f = ImgtSequenceFile::load(v).unwrap();
            assert_eq!(f.dict.len(), f.sequences_by_contig.len());
            let contigs = f.index.contigs();
            assert_eq!(contigs.len(), f.dict.len());
            assert!(contigs.iter().zip(&f.dict).all(|(a, b)| a.0 == b.0 && a.1 == b.1));
        }
    }
}
