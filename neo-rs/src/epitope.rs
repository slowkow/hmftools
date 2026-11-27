//! NEO's `epitope` package: NeoEpitope with its two kinds (PmNeoEpitope for point mutations,
//! SvNeoEpitope for fusions), CodingBaseExcerpt, PointMutationData and EpitopeUtils.
//!
//! The Java classes keep their state in mutable fields that later steps read (and `checkStopLost` runs
//! two steps again over the same state); the struct keeps those fields and the methods change them in the
//! same order, so that leftovers of a first pass are seen as in Java.

use anyhow::{bail, Result};
use hmf_common::codon::codon_to_amino_acid;
use hmf_common::ensembl::{EnsemblDataCache, TranscriptData};
use hmf_common::java_hash::JavaHashSet;
use hmf_common::neo::{fusion_info, point_mutation_info, NeoEpitopeFile, NeoEpitopeFusion, NeoEpitopeType};
use hmf_common::regions::{position_within, positions_within};
use hmf_common::transcript_utils::{calc_coding_bases, tick_phase_forward, PHASE_0, PHASE_1, PHASE_2, PHASE_NONE};

use crate::compat;

pub const FS_UP: usize = 0;
pub const FS_DOWN: usize = 1;
pub const SE_START: usize = 0;
pub const SE_END: usize = 1;
pub const POS_STRAND: i8 = 1;
pub const NEG_STRAND: i8 = -1;
pub const ORIENT_FWD: i8 = 1;
pub const ORIENT_REV: i8 = -1;

/// AminoAcidRna.STOP_SYMBOL
pub const STOP_SYMBOL: &str = "_";
/// AminoAcidRna.AA_SELENOCYSTEINE
pub const AA_SELENOCYSTEINE: &str = "U";

/// NeoCommon.DOWNSTREAM_PRE_GENE_DISTANCE
pub const DOWNSTREAM_PRE_GENE_DISTANCE: i32 = 100000;
/// NeoCommon.IMMUNE_TRANSCRIPT_PREFIXES, as Java has it: "IG_J" twice and no "IG_V" (NE14)
pub const IMMUNE_TRANSCRIPT_PREFIXES: [&str; 6] = ["TR_J", "TR_V", "TR_D", "IG_J", "IG_D", "IG_J"];
/// the prefixes as intended: the V, D and J segments of both receptor families
const IMMUNE_TRANSCRIPT_PREFIXES_CORRECTED: [&str; 6] = ["TR_J", "TR_V", "TR_D", "IG_J", "IG_D", "IG_V"];

/// whether a point mutation's transcript is skipped as an IG or TR one (NeoSampleTask: "ignore IG and TR
/// transcripts"). Java's prefix list repeats IG_J in place of IG_V, so IG V-segment transcripts are kept
/// (NE14); corrected, they are skipped as the other segments are.
pub fn is_immune_transcript(trans: &TranscriptData) -> bool {
    let prefixes = if compat::java("NE14") { &IMMUNE_TRANSCRIPT_PREFIXES } else { &IMMUNE_TRANSCRIPT_PREFIXES_CORRECTED };
    prefixes.iter().any(|p| trans.bio_type.starts_with(p))
}

/// TranscriptCodingType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CodingType {
    Unknown,
    Coding,
    Utr5P,
    Utr3P,
    NonCoding,
}

/// TranscriptRegionType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RegionType {
    Unknown,
    Upstream,
    Exonic,
    Intronic,
    Downstream,
}

// ---------------------------------------------------------------------------------------------
// reference bases

/// RefGenomeInterface.getBaseString as NEO uses it: 1-based inclusive positions
pub trait RefSource {
    fn base_string(&self, chromosome: &str, start: i32, end: i32) -> Result<String>;
}

/// RefGenomeSource over an indexed FASTA, with htsjdk's checks in getSubsequenceAt: an empty string when
/// the end is one before the start, an SAMException when the start lies further after the end or the end
/// lies past the contig
pub struct FastaRef {
    pub genome: hmf_common::refgenome::RefGenome,
}

impl RefSource for FastaRef {
    fn base_string(&self, chromosome: &str, start: i32, end: i32) -> Result<String> {
        if start > end + 1 {
            bail!("htsjdk.samtools.SAMException: Malformed query; start point {start} lies after end point {end}");
        }
        let len = self.genome.chromosome_length(chromosome);
        if self.genome.sequence_names().iter().all(|n| n != chromosome) {
            bail!("htsjdk.samtools.SAMException: Unable to find entry for contig: {chromosome}");
        }
        if end > len {
            bail!(
                "htsjdk.samtools.SAMException: Query asks for data past end of contig. Query contig {chromosome} start:{start} stop:{end} contigLength:{len}"
            );
        }
        if end + 1 == start {
            return Ok(String::new());
        }
        let b = self.genome.bases(chromosome, start, end)?;
        Ok(String::from_utf8_lossy(&b).into_owned())
    }
}

// ---------------------------------------------------------------------------------------------
// cigars

/// htsjdk Cigar: elements kept as added (no merging)
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Cigar(pub Vec<(i32, u8)>);

impl Cigar {
    /// new CigarElement(length, op): IllegalArgumentException for a negative length
    pub fn element(length: i32, op: u8) -> Result<(i32, u8)> {
        if length < 0 {
            bail!(
                "java.lang.IllegalArgumentException: Cigar element being constructed with negative length: {length} and operation: {}",
                op as char
            );
        }
        Ok((length, op))
    }
    pub fn add(&mut self, length: i32, op: u8) -> Result<()> {
        self.0.push(Cigar::element(length, op)?);
        Ok(())
    }
}

impl std::fmt::Display for Cigar {
    /// Cigar.toString (TextCigarCodec.encode): "*" for an empty cigar
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.0.is_empty() {
            return f.write_str("*");
        }
        for (len, op) in &self.0 {
            write!(f, "{len}{}", *op as char)?;
        }
        Ok(())
    }
}

/// CodingBaseExcerpt
#[derive(Clone, Debug, PartialEq)]
pub struct CodingBaseExcerpt {
    pub bases: String,
    pub positions: [i32; 2],
    pub cigar: Cigar,
}

// ---------------------------------------------------------------------------------------------
// point mutations

/// PointMutationData
#[derive(Clone, Debug, PartialEq)]
pub struct PointMutationData {
    pub chromosome: String,
    pub position: i32,
    pub ref_bases: String,
    pub alt: String,
    pub gene: String,
    pub variant_copy_number: f64,
    pub copy_number: f64,
    pub subclonal_likelihood: f64,
}

/// the PmNeoEpitope part
#[derive(Clone, Debug)]
pub struct PmPart<'a> {
    pub pm: &'a PointMutationData,
    pub indel_base_diff: i32,
    pub ty: NeoEpitopeType,
    pub wildtype_acids: String,
}

/// the SvNeoEpitope part
#[derive(Clone, Debug)]
pub struct SvPart<'a> {
    pub fusion: &'a NeoEpitopeFusion,
    pub skipped_acceptor_donors: [i32; 2],
}

#[derive(Clone, Debug)]
pub enum Kind<'a> {
    Pm(PmPart<'a>),
    Sv(SvPart<'a>),
}

/// NeoEpitope
#[derive(Clone, Debug)]
pub struct NeoEpitope<'a> {
    pub kind: Kind<'a>,
    pub trans_data: [Option<&'a TranscriptData>; 2],
    pub phases: [i32; 2],
    pub exon_rank: [i32; 2],
    pub coding_type: [CodingType; 2],
    pub region_type: [RegionType; 2],
    pub raw_coding_bases: [String; 2],
    pub coding_bases: [String; 2],
    pub novel_base_index: [i32; 2],
    pub novel_codon_bases: String,
    pub upstream_acids: String,
    pub downstream_acids: String,
    pub novel_acid: String,
    pub nmd_bases_min: i32,
    pub nmd_bases_max: i32,
    pub coding_bases_length_min: i32,
    pub coding_bases_length_max: i32,
    pub upstream_wild_type_acids: String,
    pub valid: bool,
    pub ext_coding_bases: [String; 2],
    pub ext_positions: [[i32; 2]; 2],
    pub ext_cigars: [Option<Cigar>; 2],
}

fn is_codon_multiple(bases: i32) -> bool {
    bases % 3 == 0
}

fn substr(s: &str, start: usize, end: usize) -> &str {
    &s[start..end]
}

impl<'a> NeoEpitope<'a> {
    fn new(kind: Kind<'a>) -> NeoEpitope<'a> {
        NeoEpitope {
            kind,
            trans_data: [None, None],
            phases: [-1, -1],
            exon_rank: [0, 0],
            coding_type: [CodingType::Unknown; 2],
            region_type: [RegionType::Unknown; 2],
            raw_coding_bases: [String::new(), String::new()],
            coding_bases: [String::new(), String::new()],
            novel_base_index: [0, 0],
            novel_codon_bases: String::new(),
            upstream_acids: String::new(),
            downstream_acids: String::new(),
            novel_acid: String::new(),
            nmd_bases_min: 0,
            nmd_bases_max: 0,
            coding_bases_length_min: 0,
            coding_bases_length_max: 0,
            upstream_wild_type_acids: String::new(),
            valid: true,
            ext_coding_bases: [String::new(), String::new()],
            ext_positions: [[0, 0], [0, 0]],
            ext_cigars: [None, None],
        }
    }

    /// new PmNeoEpitope(pointMutation)
    pub fn point_mutation(pm: &'a PointMutationData) -> NeoEpitope<'a> {
        let indel_base_diff = pm.alt.len() as i32 - pm.ref_bases.len() as i32;
        let ty = if indel_base_diff == 0 {
            NeoEpitopeType::Missense
        } else if is_codon_multiple(indel_base_diff) {
            if indel_base_diff > 0 { NeoEpitopeType::InframeInsertion } else { NeoEpitopeType::InframeDeletion }
        } else {
            NeoEpitopeType::Frameshift
        };
        NeoEpitope::new(Kind::Pm(PmPart { pm, indel_base_diff, ty, wildtype_acids: String::new() }))
    }

    /// new SvNeoEpitope(fusion)
    pub fn fusion(fusion: &'a NeoEpitopeFusion) -> NeoEpitope<'a> {
        NeoEpitope::new(Kind::Sv(SvPart { fusion, skipped_acceptor_donors: [0, 0] }))
    }

    pub fn trans(&self, stream: usize) -> &'a TranscriptData {
        self.trans_data[stream].expect("transcript data set")
    }

    pub fn orientation(&self, fs: usize) -> i8 {
        if (fs == FS_UP) == (self.trans(fs).strand == POS_STRAND) { ORIENT_FWD } else { ORIENT_REV }
    }

    pub fn strand(&self, stream: usize) -> i8 {
        self.trans(stream).strand
    }

    fn pm(&self) -> &PmPart<'a> {
        match &self.kind {
            Kind::Pm(p) => p,
            Kind::Sv(_) => unreachable!(),
        }
    }

    fn sv(&self) -> &SvPart<'a> {
        match &self.kind {
            Kind::Sv(s) => s,
            Kind::Pm(_) => unreachable!(),
        }
    }

    fn pos_strand(&self) -> bool {
        self.trans(FS_UP).strand == POS_STRAND
    }

    pub fn position(&self, stream: usize) -> i32 {
        match &self.kind {
            Kind::Sv(s) => s.fusion.positions[stream],
            Kind::Pm(p) => {
                // position for a point mutation is defined on the upstream as the base of the mutation
                // and on the downstream as a first base after the end of the mutated section
                // for negative strand genes, the position is set to the upstream end of the mutated section
                let pm_position = p.pm.position;
                let modified_bases = if p.indel_base_diff == 0 { p.pm.alt.len() as i32 } else { 1 };
                let deleted_bases = if p.indel_base_diff < 0 { p.indel_base_diff.abs() } else { 0 };
                if self.pos_strand() {
                    if stream == FS_UP { pm_position } else { pm_position + deleted_bases + modified_bases }
                } else {
                    let neg_pos_adj = if p.indel_base_diff == 0 { -1 } else { 0 };
                    if stream == FS_UP { pm_position + deleted_bases + modified_bases + neg_pos_adj } else { pm_position + neg_pos_adj }
                }
            }
        }
    }

    pub fn chromosome(&self, stream: usize) -> &str {
        match &self.kind {
            Kind::Sv(s) => &s.fusion.chromosomes[stream],
            Kind::Pm(p) => &p.pm.chromosome,
        }
    }

    pub fn gene_name(&self, stream: usize) -> &str {
        match &self.kind {
            Kind::Sv(s) => &s.fusion.gene_names[stream],
            Kind::Pm(p) => &p.pm.gene,
        }
    }

    pub fn variant_type(&self) -> NeoEpitopeType {
        match &self.kind {
            Kind::Sv(_) => {
                if self.phase_matched() { NeoEpitopeType::InframeFusion } else { NeoEpitopeType::OutOfFrameFusion }
            }
            Kind::Pm(p) => p.ty,
        }
    }

    pub fn variant_info(&self) -> String {
        match &self.kind {
            Kind::Sv(s) => fusion_info(&s.fusion.chromosomes, &s.fusion.positions, &s.fusion.orientations),
            Kind::Pm(p) => point_mutation_info(&p.pm.chromosome, p.pm.position, &p.pm.ref_bases, &p.pm.alt),
        }
    }

    pub fn variant_copy_number(&self) -> f64 {
        match &self.kind {
            Kind::Sv(s) => s.fusion.junction_copy_number,
            Kind::Pm(p) => p.pm.variant_copy_number,
        }
    }
    pub fn copy_number(&self) -> f64 {
        match &self.kind {
            Kind::Sv(s) => s.fusion.copy_number,
            Kind::Pm(p) => p.pm.copy_number,
        }
    }
    pub fn subclonal_likelihood(&self) -> f64 {
        match &self.kind {
            Kind::Sv(_) => 0.0,
            Kind::Pm(p) => p.pm.subclonal_likelihood,
        }
    }

    pub fn phase_matched(&self) -> bool {
        match &self.kind {
            Kind::Pm(p) => is_codon_multiple(p.indel_base_diff.abs()) && p.ty != NeoEpitopeType::StopLost,
            Kind::Sv(_) => {
                if self.region_type[FS_UP] == RegionType::Exonic && self.region_type[FS_DOWN] == RegionType::Exonic {
                    tick_phase_forward(self.phases[FS_UP], 1) == self.phases[FS_DOWN]
                } else {
                    self.phases[FS_UP] == self.phases[FS_DOWN]
                }
            }
        }
    }

    pub fn unspliced_distance(&self) -> i32 {
        match &self.kind {
            Kind::Pm(_) => 0,
            Kind::Sv(s) => {
                let f = s.fusion;
                let up_end =
                    if f.orientations[FS_UP] == ORIENT_FWD { self.ext_positions[FS_UP][SE_END] } else { self.ext_positions[FS_UP][SE_START] };
                let up_unspliced = (self.position(FS_UP) - up_end).abs();
                let down_end = if f.orientations[FS_DOWN] == ORIENT_FWD {
                    self.ext_positions[FS_DOWN][SE_END]
                } else {
                    self.ext_positions[FS_DOWN][SE_START]
                };
                let down_unspliced = (self.position(FS_DOWN) - down_end).abs();
                up_unspliced + down_unspliced + f.chain_length
            }
        }
    }

    pub fn skipped_donors(&self) -> i32 {
        match &self.kind {
            Kind::Pm(_) => 0,
            Kind::Sv(s) => s.skipped_acceptor_donors[FS_UP],
        }
    }
    pub fn skipped_acceptors(&self) -> i32 {
        match &self.kind {
            Kind::Pm(_) => 0,
            Kind::Sv(s) => s.skipped_acceptor_donors[FS_DOWN],
        }
    }

    pub fn wildtype_acids(&self) -> &str {
        match &self.kind {
            Kind::Pm(p) => &p.wildtype_acids,
            Kind::Sv(_) => "",
        }
    }

    pub fn amino_acid_string(&self) -> String {
        format!("{}{}{}", self.upstream_acids, self.novel_acid, self.downstream_acids)
    }

    /// toString, for log messages
    pub fn describe(&self) -> String {
        match &self.kind {
            Kind::Pm(p) => format!("pointMut({}: {}:{} {} -> {})", p.pm.gene, p.pm.chromosome, p.pm.position, p.pm.ref_bases, p.pm.alt),
            Kind::Sv(s) => {
                let f = s.fusion;
                format!(
                    "fusion up({}:{} {}:{}:{}) down({}:{} {}:{}:{}) phased({})",
                    f.gene_names[FS_UP],
                    self.trans_data[FS_UP].map(|t| t.trans_name.as_str()).unwrap_or(""),
                    f.chromosomes[FS_UP],
                    f.positions[FS_UP],
                    f.orientations[FS_UP],
                    f.gene_names[FS_DOWN],
                    self.trans_data[FS_DOWN].map(|t| t.trans_name.as_str()).unwrap_or(""),
                    f.chromosomes[FS_DOWN],
                    f.positions[FS_DOWN],
                    f.orientations[FS_DOWN],
                    self.phase_matched()
                )
            }
        }
    }

    // -----------------------------------------------------------------------------------------
    // transcript data

    pub fn set_transcript_data(&mut self, up: &'a TranscriptData, down: &'a TranscriptData) {
        match &self.kind {
            Kind::Pm(_) => {
                let trans = up;
                self.trans_data = [Some(trans), Some(trans)];
                // set the data for the lower part of the mutation
                let lower = if trans.strand == POS_STRAND { FS_UP } else { FS_DOWN };
                let p = self.position(lower);
                set_transcript_context(self, trans, p, lower);
                set_transcript_coding_data(self, trans, p, 0, lower);
                let upper = 1 - lower;
                // for DELs, set the downstream data as well since it can cross exon-boundaries and/or affect coding bases
                let p = self.position(upper);
                set_transcript_context(self, trans, p, upper);
                set_transcript_coding_data(self, trans, p, 0, upper);
            }
            Kind::Sv(s) => {
                let ins_seq_length = s.fusion.insert_sequence.len() as i32;
                self.trans_data[FS_UP] = Some(up);
                let p = self.position(FS_UP);
                set_transcript_context(self, up, p, FS_UP);
                set_transcript_coding_data(self, up, p, ins_seq_length, FS_UP);
                self.trans_data[FS_DOWN] = Some(down);
                // if the upstream context is intronic, then skip past any downstream exonic section
                let p = self.position(FS_DOWN);
                set_transcript_context(self, down, p, FS_DOWN);
                set_transcript_coding_data(self, down, p, 0, FS_DOWN);
                if self.region_type[FS_UP] == RegionType::Intronic && self.region_type[FS_DOWN] == RegionType::Exonic {
                    let pd = self.position(FS_DOWN);
                    if let Some(exon) = down.exons.iter().find(|x| position_within(pd, x.start, x.end)) {
                        self.phases[FS_DOWN] = exon.phase_end;
                        self.exon_rank[FS_DOWN] = exon.rank + 1;
                    }
                } else if self.region_type[FS_DOWN] == RegionType::Upstream {
                    self.set_downstream_phase_from_utr();
                }
            }
        }
    }

    /// SvNeoEpitope.setDownstreamPhaseFromUtr
    fn set_downstream_phase_from_utr(&mut self) {
        let bp = self.position(FS_DOWN);
        let trans = self.trans(FS_DOWN);
        if trans.pos_strand() {
            for exon in &trans.exons {
                if exon.rank == 1 || exon.start <= bp {
                    continue;
                }
                self.phases[FS_DOWN] = exon.phase_start;
                self.exon_rank[FS_DOWN] = exon.rank;
                break;
            }
        } else {
            for exon in trans.exons.iter().rev() {
                if exon.end >= bp {
                    continue;
                }
                self.phases[FS_DOWN] = exon.phase_end;
                self.exon_rank[FS_DOWN] = exon.rank;
                break;
            }
        }
    }

    // -----------------------------------------------------------------------------------------
    // coding bases

    /// NeoEpitope.setCodingBases
    pub fn set_coding_bases(&mut self, ref_genome: &dyn RefSource, req_amino_acids: i32) -> Result<()> {
        self.extract_coding_bases(ref_genome, req_amino_acids)?;
        adjust_coding_bases_for_strand(self);
        let is_phased = self.phase_matched();
        let novel_up = self.novel_base_index[FS_UP];
        let novel_down = self.novel_base_index[FS_DOWN];
        if novel_up < 0 || novel_down < 0 {
            hmf_common::error!(
                "ne({}) invalid upBases({} noveBases={}) or downBases({} noveBases={})",
                self.describe(),
                self.coding_bases[FS_UP],
                novel_up,
                self.coding_bases[FS_DOWN],
                novel_down
            );
            return Ok(());
        }
        // if upstream ends on a phase other than 0, need to take the bases from the downstream gene to make a novel codon
        if novel_up > 0 || novel_down > 0 || !is_phased {
            let mut upstream = self.coding_bases[FS_UP].clone();
            let mut downstream = self.coding_bases[FS_DOWN].clone();
            if novel_up as usize > upstream.len() || novel_down as usize > downstream.len() {
                hmf_common::debug!(
                    "ne({}) invalid upBases({} noveBases={}) or downBases({} noveBases={})",
                    self.describe(),
                    upstream,
                    novel_up,
                    downstream,
                    novel_down
                );
                self.valid = false;
                return Ok(());
            }
            // take the last 1 or 2 bases from the end of upstream gene's section
            let cut = upstream.len() - novel_up as usize;
            self.novel_codon_bases = upstream[cut..].to_string();
            upstream.truncate(cut);
            if is_phased {
                self.novel_codon_bases.push_str(&downstream[..novel_down as usize]);
                downstream = downstream[novel_down as usize..].to_string();
            } else {
                self.novel_codon_bases.push_str(&downstream);
                downstream = String::new();
            }
            self.coding_bases[FS_UP] = upstream;
            self.coding_bases[FS_DOWN] = downstream;
        }
        // check superflous upstream bases from longer inserts or INDELs
        let required_length = (req_amino_acids * 3) as usize;
        if self.coding_bases[FS_UP].len() > required_length {
            let n = self.coding_bases[FS_UP].len();
            self.coding_bases[FS_UP] = self.coding_bases[FS_UP][n - required_length..].to_string();
        }
        Ok(())
    }

    fn extract_coding_bases(&mut self, ref_genome: &dyn RefSource, required_amino_acids: i32) -> Result<()> {
        match &self.kind {
            Kind::Pm(_) => self.pm_extract_coding_bases(ref_genome, required_amino_acids),
            Kind::Sv(_) => self.sv_extract_coding_bases(ref_genome, required_amino_acids),
        }
    }

    /// PmNeoEpitope.extractCodingBases
    fn pm_extract_coding_bases(&mut self, ref_genome: &dyn RefSource, required_amino_acids: i32) -> Result<()> {
        let (pm, indel_base_diff) = {
            let p = self.pm();
            (p.pm, p.indel_base_diff)
        };
        let is_base_change = indel_base_diff == 0;
        let is_indel = !is_base_change;
        let is_insert = indel_base_diff > 0;
        let is_deletion = indel_base_diff < 0;
        let pos_strand = self.pos_strand();
        let up_trans = self.trans(FS_UP);
        let down_trans = self.trans(FS_DOWN);

        let up_extra_bases = get_upstream_open_codon_bases(self.phases[FS_UP], up_trans.strand, is_indel);
        let down_extra_bases = get_downstream_open_codon_bases(self.phases[FS_UP], up_trans.strand, indel_base_diff, &pm.alt);
        let mut up_required_bases = required_amino_acids * 3 + up_extra_bases;
        let mut up_position = self.position(FS_UP);
        let mut down_position = self.position(FS_DOWN);

        // adjust the position for getting reference bases to make way for the ALT to be added
        if pos_strand {
            up_position = pm.position - 1;
            if !is_base_change && up_extra_bases == 0 && self.phases[FS_UP] == PHASE_0 {
                // since ref base will be the same, require 1 less coding bases upstream
                up_required_bases -= 1;
            }
        } else {
            if is_base_change {
                up_position += 1; // since would otherwise cross with the SNV/MNV
            }
            down_position = pm.position - 1; // since down pos on -ve strand is the mutation base
        }

        let chr_up = self.chromosome(FS_UP).to_string();
        let chr_down = self.chromosome(FS_DOWN).to_string();
        let Some(up_excerpt) =
            get_upstream_coding_base_excerpt(ref_genome, up_trans, &chr_up, up_position, self.orientation(FS_UP), up_required_bases)?
        else {
            // at or past end of transcript is invalid
            self.valid = false;
            return Ok(());
        };
        let mut up_coding_bases = up_excerpt.bases.clone();

        let can_start_in_exon = true;
        let mut down_required_bases = required_amino_acids * 3 + down_extra_bases;
        let Some(mut down_excerpt) = get_downstream_coding_base_excerpt(
            ref_genome,
            down_trans,
            &chr_down,
            down_position,
            self.orientation(FS_DOWN),
            down_required_bases,
            can_start_in_exon,
            false,
            !self.phase_matched(),
        )?
        else {
            // at or past end of transcript is invalid
            self.valid = false;
            return Ok(());
        };
        let mut down_coding_bases = down_excerpt.bases.clone();

        if is_base_change {
            self.set_wildtype_amino_acids(ref_genome, required_amino_acids)?;
        }

        // the ref bases was excluded and is now replaced by the alt
        let alt_length = pm.alt.len() as i32;
        let up_phase = self.phases[FS_UP];

        if pos_strand {
            up_coding_bases.push_str(&pm.alt);
            if is_base_change {
                self.novel_base_index[FS_UP] = up_extra_bases + alt_length;
            } else if up_phase == PHASE_0 {
                // if the mutation occurs on the last base of a codon, then then novel codon starts immediately after
                self.novel_base_index[FS_UP] = if is_deletion { 0 } else { (alt_length - 1).max(0) };
            } else {
                self.novel_base_index[FS_UP] = up_extra_bases + alt_length;
            }
            if self.phase_matched() {
                self.novel_base_index[FS_DOWN] = down_extra_bases;
            }
            self.ext_coding_bases[FS_UP] = format!("{}{}{}", up_excerpt.bases, pm.alt, down_excerpt.bases);
        } else {
            down_coding_bases.push_str(&pm.alt);
            self.novel_base_index[FS_UP] = up_extra_bases;
            if self.phase_matched() {
                self.novel_base_index[FS_DOWN] = if is_deletion && up_phase == PHASE_0 { 0 } else { down_extra_bases + alt_length };
            }
            self.ext_coding_bases[FS_UP] = format!("{}{}{}", down_excerpt.bases, pm.alt, up_excerpt.bases);
        }

        if !self.phase_matched() {
            // ensure a codon-rounding base total
            if is_deletion {
                down_required_bases = (required_amino_acids * 2 + 1) * 3 - up_coding_bases.len() as i32;
            }
            // a required count of 0 gives null (Java's NullPointerException follows)
            let Some(e) = get_downstream_coding_base_excerpt(
                ref_genome,
                down_trans,
                &chr_down,
                down_position,
                self.orientation(FS_DOWN),
                down_required_bases,
                can_start_in_exon,
                false,
                false,
            )?
            else {
                bail!("java.lang.NullPointerException: no downstream coding bases for {}", self.describe());
            };
            down_excerpt = e;
            // NE5: Java keeps the coding bases of the first, unrestricted fetch (every base to the end of the
            // transcript), while the positions and the cigar come from this one; on the reverse strand the
            // extra bases come first, so ISOFOX reads them against the wrong positions. Corrected: the bases
            // of the restricted fetch, as SvNeoEpitope takes them.
            if !compat::java("NE5") {
                self.ext_coding_bases[FS_UP] = if pos_strand {
                    format!("{}{}{}", up_excerpt.bases, pm.alt, down_excerpt.bases)
                } else {
                    format!("{}{}{}", down_excerpt.bases, pm.alt, up_excerpt.bases)
                };
            }
        }

        let (lower, upper) = if pos_strand { (&up_excerpt, &down_excerpt) } else { (&down_excerpt, &up_excerpt) };

        // set cigar details just for the upstream slot since it's for a single gene
        self.ext_positions[FS_UP][SE_START] = lower.positions[SE_START];
        self.ext_positions[FS_UP][SE_END] = upper.positions[SE_END];

        let mut elements: Vec<(i32, u8)> = lower.cigar.0.clone();
        if elements.is_empty() {
            bail!("java.lang.IndexOutOfBoundsException: Index -1 out of bounds for length 0 ({})", self.describe());
        }
        // NE2: Java adds the mutated bases straight after the lower excerpt and the upper excerpt straight
        // after them, so an intron between them (a mutation on the first or last base of an exon) is lost
        // from the cigar. Corrected: an N element for each such gap.
        let ref_length = pm.ref_bases.len() as i32;
        let lower_gap = pm.position - lower.positions[SE_END] - 1;
        let upper_gap = upper.positions[SE_START] - (pm.position + ref_length);
        let fix_gaps = !compat::java("NE2");
        if fix_gaps && lower_gap > 0 {
            elements.push(Cigar::element(lower_gap, b'N')?);
        }
        let last = *elements.last().unwrap();
        let n = elements.len();
        if is_insert {
            // add the unchanged base
            if last.1 == b'M' {
                elements[n - 1] = Cigar::element(last.0 + 1, b'M')?;
            } else {
                elements.push(Cigar::element(1, b'M')?);
            }
            elements.push(Cigar::element(alt_length - 1, b'I')?);
        } else {
            // add the unchanged base(s)
            if last.1 == b'M' {
                elements[n - 1] = Cigar::element(last.0 + alt_length, b'M')?;
            } else {
                elements.push(Cigar::element(alt_length, b'M')?);
            }
            if is_deletion {
                elements.push(Cigar::element(indel_base_diff.abs(), b'D')?);
            }
        }
        if fix_gaps && upper_gap > 0 && !upper.cigar.0.is_empty() {
            elements.push(Cigar::element(upper_gap, b'N')?);
        }
        let last = *elements.last().unwrap();
        let mut up_elements: Vec<(i32, u8)> = upper.cigar.0.clone();
        if up_elements.is_empty() {
            bail!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0 ({})", self.describe());
        }
        if last.1 == b'M' && up_elements[0].1 == b'M' {
            let n = elements.len();
            elements[n - 1] = Cigar::element(last.0 + up_elements[0].0, b'M')?;
            up_elements.remove(0);
        }
        elements.extend(up_elements);
        self.ext_cigars[FS_UP] = Some(Cigar(elements));

        self.raw_coding_bases[FS_UP] = up_coding_bases;
        self.raw_coding_bases[FS_DOWN] = down_coding_bases;
        Ok(())
    }

    /// PmNeoEpitope.setWildtypeAminoAcids
    fn set_wildtype_amino_acids(&mut self, ref_genome: &dyn RefSource, required_amino_acids: i32) -> Result<()> {
        let pm = self.pm().pm;
        let trans = self.trans(FS_UP);
        let pos_strand = trans.strand == POS_STRAND;
        let cb = calc_coding_bases(trans, pm.position);
        let up_open = if cb.phase == PHASE_0 { 0 } else { cb.phase };
        let up_required = required_amino_acids * 3 + up_open;
        let up_orient = if pos_strand { ORIENT_FWD } else { ORIENT_REV };
        let chr = self.chromosome(FS_UP).to_string();
        let mut up_bases = get_upstream_coding_base_excerpt(ref_genome, trans, &chr, pm.position, up_orient, up_required)?
            .ok_or_else(|| anyhow::anyhow!("java.lang.NullPointerException: no wild-type upstream bases for {}", self.describe()))?
            .bases;
        let down_orient = if pos_strand { ORIENT_REV } else { ORIENT_FWD };
        let down_position = if pos_strand { pm.position + 1 } else { pm.position - 1 };
        let down_required = required_amino_acids * 3 - up_open;
        let mut down_bases =
            get_downstream_coding_base_excerpt(ref_genome, trans, &chr, down_position, down_orient, down_required, true, false, false)?
                .ok_or_else(|| anyhow::anyhow!("java.lang.NullPointerException: no wild-type downstream bases for {}", self.describe()))?
                .bases;
        if !pos_strand {
            up_bases = reverse_complement_bases(&up_bases);
            down_bases = reverse_complement_bases(&down_bases);
        }
        if up_open > 0 {
            let open = up_open as usize;
            if open > up_bases.len() {
                bail!("java.lang.StringIndexOutOfBoundsException: wild-type bases of {}", self.describe());
            }
            let cut = up_bases.len() - open;
            down_bases = format!("{}{}", &up_bases[cut..], down_bases);
            up_bases.truncate(cut);
        }
        let wt = format!("{}{}", get_amino_acids(&up_bases, false), get_amino_acids(&down_bases, true));
        if let Kind::Pm(p) = &mut self.kind {
            p.wildtype_acids = wt;
        }
        Ok(())
    }

    /// SvNeoEpitope: upstream open codon bases, taking the insert sequence out of the phase
    fn sv_upstream_open_codon_bases(&self) -> i32 {
        let f = self.sv().fusion;
        let mut phase = self.phases[FS_UP];
        if self.region_type[FS_UP] == RegionType::Exonic && !f.insert_sequence.is_empty() {
            let ins = f.insert_sequence.len() as i32;
            let pre_insert_phase = (phase - ins) % 3;
            if pre_insert_phase >= 0 {
                phase = pre_insert_phase;
            } else if pre_insert_phase.abs() == PHASE_1 {
                phase = PHASE_2;
            } else {
                phase = PHASE_1;
            }
        }
        if phase == PHASE_NONE { 0 } else { phase }
    }

    /// SvNeoEpitope.extractCodingBases
    fn sv_extract_coding_bases(&mut self, ref_genome: &dyn RefSource, required_amino_acids: i32) -> Result<()> {
        let f = self.sv().fusion;
        let up_extra_bases = self.sv_upstream_open_codon_bases();
        let up_open = self.phases[FS_UP];
        let down_extra_bases = if up_open == 0 { 0 } else { 3 - up_open };
        let up_required_bases = required_amino_acids * 3 + up_extra_bases;

        let up_trans = self.trans(FS_UP);
        let down_trans = self.trans(FS_DOWN);
        let chr_up = self.chromosome(FS_UP).to_string();
        let chr_down = self.chromosome(FS_DOWN).to_string();
        let Some(up) =
            get_upstream_coding_base_excerpt(ref_genome, up_trans, &chr_up, self.position(FS_UP), self.orientation(FS_UP), up_required_bases)?
        else {
            bail!("java.lang.NullPointerException: no upstream coding bases for {}", self.describe());
        };
        self.raw_coding_bases[FS_UP] = up.bases.clone();
        self.ext_positions[FS_UP] = up.positions;
        let mut up_elements = up.cigar.0.clone();
        let mut coding_ins_seq_len = 0;
        if self.region_type[FS_UP] == RegionType::Exonic && !f.insert_sequence.is_empty() {
            coding_ins_seq_len = f.insert_sequence.len() as i32;
            if up_trans.strand == POS_STRAND {
                up_elements.push(Cigar::element(coding_ins_seq_len, b'I')?);
                self.raw_coding_bases[FS_UP].push_str(&f.insert_sequence);
            } else {
                up_elements.insert(0, Cigar::element(coding_ins_seq_len, b'I')?);
                self.raw_coding_bases[FS_UP] = format!("{}{}", f.insert_sequence, self.raw_coding_bases[FS_UP]);
            }
        }
        self.ext_cigars[FS_UP] = Some(Cigar(up_elements));
        self.ext_coding_bases[FS_UP] = self.raw_coding_bases[FS_UP].clone();

        self.novel_base_index[FS_UP] = up_extra_bases + coding_ins_seq_len;
        self.novel_base_index[FS_DOWN] = down_extra_bases;

        let can_start_in_exon = self.region_type[FS_UP] == RegionType::Exonic || up_extra_bases > 0;
        let down_required_bases = required_amino_acids * 3 + down_extra_bases;
        let Some(mut down) = get_downstream_coding_base_excerpt(
            ref_genome,
            down_trans,
            &chr_down,
            self.position(FS_DOWN),
            self.orientation(FS_DOWN),
            down_required_bases,
            can_start_in_exon,
            true,
            !self.phase_matched(),
        )?
        else {
            self.valid = false;
            return Ok(());
        };
        self.raw_coding_bases[FS_DOWN] = down.bases.clone();

        // call again to get restricted downstream bases
        if !self.phase_matched() {
            down = get_downstream_coding_base_excerpt(
                ref_genome,
                down_trans,
                &chr_down,
                self.position(FS_DOWN),
                self.orientation(FS_DOWN),
                down_required_bases,
                can_start_in_exon,
                true,
                false,
            )?
            .ok_or_else(|| anyhow::anyhow!("java.lang.NullPointerException: no downstream coding bases for {}", self.describe()))?;
        }
        if down.bases.is_empty() || down.positions[SE_START] == 0 || down.positions[SE_END] == 0 {
            self.valid = false;
            return Ok(());
        }
        self.ext_positions[FS_DOWN] = down.positions;
        self.ext_coding_bases[FS_DOWN] = down.bases.clone();
        self.ext_cigars[FS_DOWN] = Some(down.cigar);
        Ok(())
    }

    // -----------------------------------------------------------------------------------------
    // amino acids

    /// NeoEpitope.setNonsenseMediatedDecay
    pub fn set_nonsense_mediated_decay(&mut self) {
        let nmd = calc_non_mediated_decay_bases(self);
        self.nmd_bases_min = nmd;
        self.nmd_bases_max = nmd;
        // also record the bases from the start codon to the stop codon (or end of transcript if there is none)
        let cbl = calc_start_codon_bases(self) + calc_stop_codon_bases(self);
        self.coding_bases_length_min = cbl;
        self.coding_bases_length_max = cbl;
    }

    /// NeoEpitope.setAminoAcids
    pub fn set_amino_acids(&mut self, ref_genome: &dyn RefSource, req_wildtype_amino_acids: i32) -> Result<()> {
        self.upstream_acids = get_amino_acids(trim_incomplete_codons(&self.coding_bases[FS_UP]), false);
        self.novel_acid = get_amino_acids(&self.novel_codon_bases, true);
        self.downstream_acids = get_amino_acids(&self.coding_bases[FS_DOWN], true);

        // truncate the novel bases and AAs if a stop codon is encountered
        if self.novel_acid.ends_with(STOP_SYMBOL) {
            let novel_bases = self.novel_acid.len() * 3;
            self.novel_codon_bases.truncate(novel_bases);
            self.downstream_acids = String::new();
        }

        // cache bases downstream of mutation in the reference to check for non-novel AAs later on
        let up_aa_pos_start =
            if self.orientation(FS_UP) == ORIENT_FWD { self.ext_positions[FS_UP][SE_START] } else { self.ext_positions[FS_UP][SE_END] };
        let wt_orient = if self.orientation(FS_UP) == ORIENT_FWD { ORIENT_REV } else { ORIENT_FWD };
        let required_bases = self.coding_bases[FS_UP].len() as i32 + req_wildtype_amino_acids * 3;
        let chr = self.chromosome(FS_UP).to_string();
        let up_trans = self.trans(FS_UP);
        if let Some(wt) =
            get_downstream_coding_base_excerpt(ref_genome, up_trans, &chr, up_aa_pos_start, wt_orient, required_bases, true, false, false)?
        {
            let mut bases = wt.bases;
            if self.strand(FS_UP) == NEG_STRAND {
                bases = reverse_complement_bases(&bases);
            }
            self.upstream_wild_type_acids = get_amino_acids(trim_incomplete_codons(&bases), true);
        }
        self.check_stop_lost(ref_genome, req_wildtype_amino_acids)
    }

    /// PmNeoEpitope.checkStopLost (nothing for fusions)
    fn check_stop_lost(&mut self, ref_genome: &dyn RefSource, req_wildtype_amino_acids: i32) -> Result<()> {
        let Kind::Pm(p) = &self.kind else { return Ok(()) };
        if p.ty != NeoEpitopeType::Missense || !p.wildtype_acids.ends_with(STOP_SYMBOL) {
            return Ok(());
        }
        // check for stop-lost in point mutations mis-classified as missense
        if self.downstream_acids.is_empty() && !self.novel_acid.ends_with(STOP_SYMBOL) {
            if let Kind::Pm(p) = &mut self.kind {
                p.ty = NeoEpitopeType::StopLost;
            }
            // now recalculate downstream bases and AAs
            self.set_coding_bases(ref_genome, req_wildtype_amino_acids)?;
            self.set_amino_acids(ref_genome, req_wildtype_amino_acids)?;
        }
        Ok(())
    }

    /// SvNeoEpitope.setSkippedSpliceSites (nothing for point mutations)
    pub fn set_skipped_splice_sites(&mut self, cache: &EnsemblDataCache) {
        let Kind::Sv(s) = &self.kind else { return };
        let f = s.fusion;
        let mut skipped = [0, 0];
        for fs in [FS_UP, FS_DOWN] {
            let find_exon_start = f.orientations[fs] == ORIENT_REV;
            let bounds = if f.orientations[fs] == ORIENT_FWD {
                [self.ext_positions[fs][SE_END], self.position(fs)]
            } else {
                [self.position(fs), self.ext_positions[fs][SE_START]]
            };
            let trans_id = self.trans(fs).trans_id;
            let candidates: Vec<&TranscriptData> = match cache.transcripts(&f.gene_ids[fs]) {
                Some(list) => list.iter().filter(|x| x.trans_id != trans_id).collect(),
                None => Vec::new(),
            };
            skipped[fs] = find_skipped_exon_boundaries(&candidates, bounds, find_exon_start, fs == FS_DOWN);
        }
        if let Kind::Sv(s) = &mut self.kind {
            s.skipped_acceptor_donors = skipped;
        }
    }

    /// NeoEpitope.toFile
    pub fn to_file(&self, ne_id: i32, up_trans_names: &JavaHashSet<String>, down_trans_names: &JavaHashSet<String>) -> Result<NeoEpitopeFile> {
        let join = |s: &JavaHashSet<String>| s.iter().cloned().collect::<Vec<_>>().join(";");
        let Some(up_cigar) = &self.ext_cigars[FS_UP] else {
            bail!("java.lang.NullPointerException: no upstream cigar for {}", self.describe());
        };
        Ok(NeoEpitopeFile {
            id: ne_id,
            variant_type: self.variant_type(),
            variant_info: self.variant_info(),
            variant_copy_number: self.variant_copy_number(),
            copy_number: self.copy_number(),
            subclonal_likelihood: self.subclonal_likelihood(),
            gene_ids: [self.trans(FS_UP).gene_id.clone(), self.trans(FS_DOWN).gene_id.clone()],
            gene_names: [self.gene_name(FS_UP).to_string(), self.gene_name(FS_DOWN).to_string()],
            chromosomes: [self.chromosome(FS_UP).to_string(), self.chromosome(FS_DOWN).to_string()],
            orientations: [self.orientation(FS_UP), self.orientation(FS_DOWN)],
            nmd_bases: [self.nmd_bases_min, self.nmd_bases_max],
            coding_bases_length: [self.coding_bases_length_min, self.coding_bases_length_max],
            upstream_aa: self.upstream_acids.clone(),
            downstream_aa: self.downstream_acids.clone(),
            novel_aa: self.novel_acid.clone(),
            transcripts: [join(up_trans_names), join(down_trans_names)],
            wildtype_aa: self.wildtype_acids().to_string(),
            fused_intron_length: self.unspliced_distance(),
            skipped_acceptors_donors: [self.skipped_donors(), self.skipped_acceptors()],
            coding_base_positions: self.ext_positions,
            coding_bases: self.ext_coding_bases.clone(),
            coding_base_cigars: [up_cigar.to_string(), self.ext_cigars[FS_DOWN].as_ref().map(|c| c.to_string()).unwrap_or_default()],
        })
    }
}

// ---------------------------------------------------------------------------------------------
// EpitopeUtils

/// Nucleotides.reverseComplementBases
pub fn reverse_complement_bases(s: &str) -> String {
    String::from_utf8(hmf_common::codon::reverse_complement_bases(s.as_bytes())).unwrap_or_default()
}

/// NeoEpitope.trimIncompleteCodons: strip leading incomplete codon bases
pub fn trim_incomplete_codons(coding_bases: &str) -> &str {
    let remainder = coding_bases.len() % 3;
    &coding_bases[remainder..]
}

/// EpitopeUtils.setTranscriptContext
pub fn set_transcript_context(ne: &mut NeoEpitope, trans: &TranscriptData, position: i32, stream: usize) {
    let is_upstream = stream == FS_UP;
    if position < trans.trans_start {
        ne.exon_rank[stream] = if trans.pos_strand() { 0 } else { trans.exons.len() as i32 - 1 };
        ne.region_type[stream] = if trans.pos_strand() { RegionType::Upstream } else { RegionType::Downstream };
        ne.phases[stream] = PHASE_NONE;
        return;
    } else if position > trans.trans_end {
        ne.exon_rank[stream] = if trans.pos_strand() { trans.exons.len() as i32 - 1 } else { 0 };
        ne.region_type[stream] = if trans.pos_strand() { RegionType::Downstream } else { RegionType::Upstream };
        ne.phases[stream] = PHASE_NONE;
        return;
    }
    for exon in &trans.exons {
        if position > exon.end {
            continue;
        }
        if position < exon.start {
            ne.region_type[stream] = RegionType::Intronic;
            // upstream pos-strand, before next exon then take the exon before's rank
            if trans.strand == POS_STRAND && is_upstream {
                ne.exon_rank[stream] = exon.rank - 1;
            } else if trans.strand == NEG_STRAND && is_upstream {
                ne.exon_rank[stream] = exon.rank;
            } else if trans.strand == POS_STRAND && !is_upstream {
                ne.exon_rank[stream] = exon.rank;
            } else if trans.strand == NEG_STRAND && !is_upstream {
                ne.exon_rank[stream] = exon.rank - 1;
            }
            ne.phases[stream] = if trans.strand == POS_STRAND { exon.phase_start } else { exon.phase_end };
        } else if position_within(position, exon.start, exon.end) {
            ne.region_type[stream] = RegionType::Exonic;
            ne.exon_rank[stream] = exon.rank;
        }
        break;
    }
}

/// EpitopeUtils.setTranscriptCodingData
pub fn set_transcript_coding_data(ne: &mut NeoEpitope, trans: &TranscriptData, position: i32, ins_seq_length: i32, stream: usize) {
    match (trans.coding_start, trans.coding_end) {
        (Some(cs), Some(ce)) => {
            ne.coding_type[stream] = if position_within(position, cs, ce) {
                CodingType::Coding
            } else if trans.strand == POS_STRAND && position < cs {
                CodingType::Utr5P
            } else if trans.strand == NEG_STRAND && position > ce {
                CodingType::Utr5P
            } else {
                CodingType::Utr3P
            };
        }
        _ => {
            ne.coding_type[stream] = CodingType::NonCoding;
            ne.phases[stream] = -1;
        }
    }
    if ne.coding_type[stream] == CodingType::Coding && ne.region_type[stream] == RegionType::Exonic {
        let cb = calc_coding_bases(trans, position);
        ne.phases[stream] = tick_phase_forward(cb.phase, ins_seq_length);
    }
}

/// EpitopeUtils.getUpstreamOpenCodonBases: the open codon bases up to but not including the variant's position
pub fn get_upstream_open_codon_bases(start_phase: i32, strand: i8, is_indel: bool) -> i32 {
    let phase = start_phase;
    if strand == POS_STRAND || !is_indel {
        if phase == PHASE_1 {
            0
        } else if phase == PHASE_2 {
            1
        } else if is_indel {
            // INDEL leaves the mutation base unchanged, so if it ends a codon, it's not novel
            0
        } else {
            2
        }
    } else if phase == PHASE_1 {
        1
    } else if phase == PHASE_2 {
        2
    } else if is_indel {
        0
    } else {
        2
    }
}

/// EpitopeUtils.getDownstreamOpenCodonBases: the phase at the base after the mutation
pub fn get_downstream_open_codon_bases(start_phase: i32, strand: i8, base_change: i32, alt: &str) -> i32 {
    let mutation_ticks = if strand == POS_STRAND {
        if base_change < 0 { 1 } else { alt.len() as i32 }
    } else if base_change == 0 {
        alt.len() as i32
    } else if base_change > 0 {
        alt.len() as i32 + 1
    } else {
        2
    };
    let post = tick_phase_forward(start_phase, mutation_ticks);
    if post == PHASE_0 {
        1
    } else if post == PHASE_1 {
        0
    } else {
        2
    }
}

/// EpitopeUtils.getUpstreamCodingBases
#[cfg_attr(not(test), allow(dead_code))]
pub fn get_upstream_coding_bases(
    ref_genome: &dyn RefSource, trans: &TranscriptData, chromosome: &str, ne_position: i32, ne_orientation: i8, required_bases: i32,
) -> Result<String> {
    match get_upstream_coding_base_excerpt(ref_genome, trans, chromosome, ne_position, ne_orientation, required_bases)? {
        Some(e) => Ok(e.bases),
        None => bail!("java.lang.NullPointerException: no upstream coding bases"),
    }
}

/// EpitopeUtils.getUpstreamCodingBaseExcerpt; None where Java returns null
pub fn get_upstream_coding_base_excerpt(
    ref_genome: &dyn RefSource, trans: &TranscriptData, chromosome: &str, ne_position: i32, ne_orientation: i8, mut required_bases: i32,
) -> Result<Option<CodingBaseExcerpt>> {
    let (Some(coding_start), Some(coding_end)) = (trans.coding_start, trans.coding_end) else { return Ok(None) };
    if required_bases <= 0 {
        return Ok(None);
    }
    if !position_within(ne_position, coding_start, coding_end) {
        return Ok(None);
    }
    let exons = &trans.exons;
    let mut base_string = String::new();
    let mut cigar = Cigar::default();
    let mut last_match_base = -1;
    let mut start_pos = 0;
    let mut end_pos = 0;

    if ne_orientation == ORIENT_REV {
        // walk forwards through the exons, collecting up the required positions and bases
        for exon in exons {
            if ne_position > exon.end {
                continue;
            }
            if exon.start > coding_end {
                break; // no more coding bases
            }
            let exon_base_start = exon.start.max(ne_position);
            let exon_base_end = coding_end.min(exon.end);
            let exon_base_count = exon_base_end - exon_base_start + 1;
            let (base_start, base_end);
            if required_bases >= exon_base_count {
                base_start = exon_base_start;
                base_end = exon_base_end;
                required_bases -= exon_base_count;
            } else {
                base_start = exon_base_start;
                base_end = base_start + required_bases - 1;
                required_bases = 0;
            }
            base_string.push_str(&ref_genome.base_string(chromosome, base_start, base_end)?);
            if start_pos == 0 {
                start_pos = exon_base_start;
            }
            end_pos = base_end;
            if last_match_base > 0 {
                cigar.add(base_start - last_match_base - 1, b'N')?;
            }
            last_match_base = base_end;
            cigar.add(base_end - base_start + 1, b'M')?;
            if required_bases <= 0 {
                break;
            }
        }
    } else {
        for exon in exons.iter().rev() {
            if ne_position < exon.start {
                continue;
            }
            if exon.end < coding_start {
                break;
            }
            let exon_base_end = exon.end.min(ne_position);
            let exon_base_start = coding_start.max(exon.start);
            let exon_base_count = exon_base_end - exon_base_start + 1;
            let (base_start, base_end);
            if required_bases >= exon_base_count {
                base_start = exon_base_start;
                base_end = exon_base_end;
                required_bases -= exon_base_count;
            } else {
                base_end = exon_base_end;
                base_start = base_end - required_bases + 1;
                required_bases = 0;
            }
            base_string = format!("{}{}", ref_genome.base_string(chromosome, base_start, base_end)?, base_string);
            if end_pos == 0 {
                end_pos = exon_base_end;
            }
            start_pos = base_start;
            if last_match_base > 0 {
                // NE1: Java has `lastMatchBase - baseEnd -  - 1`, that is + 1: each intron two bases too long
                let gap = if compat::java("NE1") { last_match_base - base_end + 1 } else { last_match_base - base_end - 1 };
                cigar.add(gap, b'N')?;
            }
            last_match_base = base_start;
            cigar.add(base_end - base_start + 1, b'M')?;
            if required_bases <= 0 {
                break;
            }
        }
        cigar.0.reverse();
    }
    Ok(Some(CodingBaseExcerpt { bases: base_string, positions: [start_pos, end_pos], cigar }))
}

/// EpitopeUtils.getDownstreamCodingBases
#[allow(clippy::too_many_arguments)]
#[cfg_attr(not(test), allow(dead_code))]
pub fn get_downstream_coding_bases(
    ref_genome: &dyn RefSource, trans: &TranscriptData, chromosome: &str, ne_position: i32, ne_orientation: i8, required_bases: i32,
    can_start_in_exon: bool, req_splice_acceptor: bool, req_all_bases: bool,
) -> Result<String> {
    match get_downstream_coding_base_excerpt(
        ref_genome,
        trans,
        chromosome,
        ne_position,
        ne_orientation,
        required_bases,
        can_start_in_exon,
        req_splice_acceptor,
        req_all_bases,
    )? {
        Some(e) => Ok(e.bases),
        None => bail!("java.lang.NullPointerException: no downstream coding bases"),
    }
}

/// EpitopeUtils.getDownstreamCodingBaseExcerpt; None where Java returns null
#[allow(clippy::too_many_arguments)]
pub fn get_downstream_coding_base_excerpt(
    ref_genome: &dyn RefSource, trans: &TranscriptData, chromosome: &str, ne_position: i32, ne_orientation: i8, mut required_bases: i32,
    can_start_in_exon: bool, req_splice_acceptor: bool, req_all_bases: bool,
) -> Result<Option<CodingBaseExcerpt>> {
    if required_bases == 0 {
        return Ok(None);
    }
    let coding_start = match trans.coding_start {
        Some(cs) if ne_position >= cs => cs,
        _ => trans.trans_start,
    };
    let coding_end = match trans.coding_end {
        Some(ce) if ne_position <= ce => ce,
        _ => trans.trans_end,
    };
    // if the position falls in the stop codon, take at least the last base
    if (trans.strand == POS_STRAND && ne_position >= trans.trans_end) || (trans.strand == NEG_STRAND && ne_position <= trans.trans_start) {
        return Ok(None);
    }
    let exons = &trans.exons;
    let mut base_string = String::new();
    let mut cigar = Cigar::default();
    let mut last_match_base = -1;
    let mut start_pos = 0;
    let mut end_pos = 0;

    if ne_orientation == ORIENT_REV {
        for exon in exons {
            // starts after the first exon, ie at the first splice acceptor
            if req_splice_acceptor && exon.rank == 1 {
                continue;
            }
            if ne_position > exon.end {
                continue;
            }
            if position_within(ne_position, exon.start, exon.end) && !can_start_in_exon {
                continue; // will start at the next exon
            }
            if exon.start > coding_end && !req_all_bases {
                break; // no more coding bases
            }
            let exon_base_start = ne_position.max(exon.start);
            let exon_base_end = if !req_all_bases { coding_end.min(exon.end) } else { exon.end };
            let exon_base_count = exon_base_end - exon_base_start + 1;
            let (base_start, base_end);
            if required_bases >= exon_base_count || req_all_bases {
                base_start = exon_base_start;
                base_end = exon_base_end;
                if !req_all_bases {
                    required_bases -= exon_base_count;
                }
            } else {
                base_start = exon_base_start;
                base_end = base_start + required_bases - 1;
                required_bases = 0;
            }
            base_string.push_str(&ref_genome.base_string(chromosome, base_start, base_end)?);
            if start_pos == 0 {
                start_pos = exon_base_start;
            }
            end_pos = base_end;
            if last_match_base > 0 {
                cigar.add(base_start - last_match_base - 1, b'N')?;
            }
            last_match_base = base_end;
            cigar.add(base_end - base_start + 1, b'M')?;
            if required_bases <= 0 && !req_all_bases {
                break;
            }
        }
    } else {
        for exon in exons.iter().rev() {
            if req_splice_acceptor && exon.rank == 1 {
                continue;
            }
            if ne_position < exon.start {
                continue;
            }
            if coding_end < exon.start && !req_all_bases {
                continue;
            }
            if position_within(ne_position, exon.start, exon.end) && !can_start_in_exon {
                continue;
            }
            if exon.end < coding_start && !req_all_bases {
                break;
            }
            let exon_base_end = ne_position.min(exon.end);
            let exon_base_start = if !req_all_bases { coding_start.max(exon.start) } else { exon.start };
            let exon_base_count = exon_base_end - exon_base_start + 1;
            let (base_start, base_end);
            if required_bases >= exon_base_count || req_all_bases {
                base_start = exon_base_start;
                base_end = exon_base_end;
                if !req_all_bases {
                    required_bases -= exon_base_count;
                }
            } else {
                base_end = exon_base_end;
                base_start = base_end - required_bases + 1;
                required_bases = 0;
            }
            base_string = format!("{}{}", ref_genome.base_string(chromosome, base_start, base_end)?, base_string);
            if end_pos == 0 {
                end_pos = exon_base_end;
            }
            start_pos = base_start;
            if last_match_base > 0 {
                cigar.add(last_match_base - base_end - 1, b'N')?;
            }
            last_match_base = base_start;
            cigar.add(base_end - base_start + 1, b'M')?;
            if required_bases <= 0 && !req_all_bases {
                break;
            }
        }
        cigar.0.reverse();
    }
    Ok(Some(CodingBaseExcerpt { bases: base_string, positions: [start_pos, end_pos], cigar }))
}

/// EpitopeUtils.adjustCodingBasesForStrand
pub fn adjust_coding_bases_for_strand(ne: &mut NeoEpitope) {
    for fs in [FS_UP, FS_DOWN] {
        ne.coding_bases[fs] =
            if ne.strand(fs) == NEG_STRAND { reverse_complement_bases(&ne.raw_coding_bases[fs]) } else { ne.raw_coding_bases[fs].clone() };
    }
}

/// EpitopeUtils.calcNonMediatedDecayBases: distance from the (novel) stop codon to the last splice acceptor
pub fn calc_non_mediated_decay_bases(ne: &NeoEpitope) -> i32 {
    if !ne.novel_acid.ends_with(STOP_SYMBOL) || ne.novel_acid == STOP_SYMBOL {
        return -1;
    }
    let trans = ne.trans(FS_DOWN);
    let exons = &trans.exons;
    let ref_position = ne.position(FS_DOWN);
    let mut exonic_base_count = 0;
    if ne.orientation(FS_DOWN) == ORIENT_REV {
        for (i, exon) in exons.iter().enumerate() {
            if i == exons.len() - 1 {
                break;
            }
            if ref_position > exon.end {
                continue;
            }
            exonic_base_count += exon.end - ref_position.max(exon.start) + 1;
        }
    } else {
        for i in (0..exons.len()).rev() {
            let exon = &exons[i];
            if i == 0 {
                break;
            }
            if ref_position < exon.start {
                continue;
            }
            exonic_base_count += ref_position.min(exon.end) - exon.start + 1;
        }
    }
    let new_coding_bases = ne.novel_acid.len() as i32 * 3;
    if exonic_base_count >= new_coding_bases { exonic_base_count - new_coding_bases } else { -1 }
}

/// EpitopeUtils.calcStopCodonBases: from the downstream position to the stop codon, or the end of the
/// transcript if there is none
pub fn calc_stop_codon_bases(ne: &NeoEpitope) -> i32 {
    if ne.novel_acid.ends_with(STOP_SYMBOL) {
        return ne.novel_codon_bases.len() as i32;
    }
    let trans = ne.trans(FS_DOWN);
    let ref_position = ne.position(FS_DOWN);
    let mut count = 0;
    if ne.orientation(FS_DOWN) == ORIENT_REV {
        let coding_stop = trans.coding_end.unwrap_or(trans.trans_end);
        if coding_stop <= ref_position {
            return 0;
        }
        for exon in &trans.exons {
            if ref_position > exon.end {
                continue;
            }
            if exon.start > coding_stop {
                break;
            }
            count += coding_stop.min(exon.end) - ref_position.max(exon.start) + 1;
        }
    } else {
        let coding_stop = trans.coding_start.unwrap_or(trans.trans_start);
        if coding_stop >= ref_position {
            return 0;
        }
        for exon in trans.exons.iter().rev() {
            if ref_position < exon.start {
                continue;
            }
            if exon.end < coding_stop {
                break;
            }
            count += ref_position.min(exon.end) - coding_stop.max(exon.start) + 1;
        }
    }
    count
}

/// EpitopeUtils.calcStartCodonBases: from the last upstream base to the upstream coding start
pub fn calc_start_codon_bases(ne: &NeoEpitope) -> i32 {
    let trans = ne.trans(FS_UP);
    let (Some(coding_start), Some(coding_end)) = (trans.coding_start, trans.coding_end) else { return -1 };
    let ref_position = ne.position(FS_UP);
    let mut count = 0;
    if ne.orientation(FS_UP) == ORIENT_REV {
        let coding_position = coding_end;
        if coding_position < ref_position {
            return -1;
        }
        for exon in &trans.exons {
            if ref_position > exon.end {
                continue;
            }
            if exon.start > coding_position {
                break;
            }
            count += coding_position.min(exon.end) - ref_position.max(exon.start) + 1;
        }
    } else {
        for exon in trans.exons.iter().rev() {
            let coding_position = coding_start;
            if coding_position > ref_position {
                return -1;
            }
            if ref_position < exon.start {
                continue;
            }
            if exon.end < coding_position {
                break;
            }
            count += ref_position.min(exon.end) - coding_position.max(exon.start) + 1;
        }
    }
    count
}

/// EpitopeUtils.findSkippedExonBoundaries: distinct exon starts (or ends) strictly inside the bounds
pub fn find_skipped_exon_boundaries(trans_list: &[&TranscriptData], bounds: [i32; 2], find_exon_start: bool, is_acceptor: bool) -> i32 {
    if bounds[SE_START] >= bounds[SE_END] || trans_list.is_empty() {
        return 0;
    }
    let mut sites: std::collections::HashSet<i32> = std::collections::HashSet::new();
    for trans in trans_list {
        for exon in &trans.exons {
            if is_acceptor && exon.rank == 1 {
                continue;
            }
            let site = if find_exon_start { exon.start } else { exon.end };
            if site <= bounds[SE_START] {
                continue;
            }
            if site >= bounds[SE_END] {
                break;
            }
            sites.insert(site);
        }
    }
    sites.len() as i32
}

/// Codons.isStopCodon
fn is_stop_codon(codon: &[u8]) -> bool {
    codon == b"TAA" || codon == b"TAG" || codon == b"TGA"
}

/// EpitopeUtils.getAminoAcids: one amino acid per complete codon; with `check_stop_codon`, a stop codon
/// ends the string with the stop symbol
pub fn get_amino_acids(bases: &str, check_stop_codon: bool) -> String {
    let b = bases.as_bytes();
    if b.len() < 3 {
        return String::new();
    }
    let mut s = String::with_capacity(b.len() / 3 + 1);
    let mut index = 0;
    while index + 3 <= b.len() {
        let codon = &b[index..index + 3];
        if check_stop_codon && is_stop_codon(codon) {
            s.push_str(STOP_SYMBOL);
            break;
        }
        // AminoAcidRna.convertDnaCodonToAminoAcid
        if is_stop_codon(codon) {
            s.push_str(STOP_SYMBOL);
        } else {
            s.push(codon_to_amino_acid(codon, 0) as char);
        }
        index += 3;
    }
    s
}

/// PeptideData
#[derive(Clone, Debug, PartialEq)]
pub struct PeptideData {
    pub peptide: String,
    pub up_flank: String,
    pub down_flank: String,
}

/// EpitopeUtils.generatePeptides: the distinct peptides of each length that include a novel amino acid
/// (or, with none, the junction), with their flanks
pub fn generate_peptides(up_aas: &str, novel_aas: &str, down_aas: &str, peptide_lengths: [i32; 2], flanking_bases: i32) -> Vec<PeptideData> {
    let mut unique: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut peptides = Vec::new();
    // replace any mis-classified upstream stop-codons with 'U' = selenocysteine
    let full = format!("{}{}{}", up_aas.replace(STOP_SYMBOL, AA_SELENOCYSTEINE), novel_aas, down_aas);
    let upstream_length = up_aas.len() as i32;
    let novel_length = novel_aas.len() as i32;
    let has_novel = !novel_aas.is_empty();
    let full_total = full.len() as i32;
    let mut full_length = full_total;
    if full.ends_with(STOP_SYMBOL) {
        full_length -= 1;
    }
    for length in peptide_lengths[SE_START]..=peptide_lengths[SE_END] {
        let mut current = 0;
        loop {
            // skip to point where novel or downstream is included
            if current + length <= upstream_length {
                current += 1;
                continue;
            }
            // break if into the downstream (canonical) sequence
            if (!has_novel && current >= upstream_length) || current >= upstream_length + novel_length {
                break;
            }
            if current + length > full_length {
                break;
            }
            let peptide = substr(&full, current as usize, (current + length) as usize);
            if unique.insert(peptide.to_string()) {
                let mut up_flank = String::new();
                let mut down_flank = String::new();
                if flanking_bases > 0 {
                    if current > 0 {
                        up_flank = full[(current - flanking_bases).max(0) as usize..current as usize].to_string();
                    }
                    if current + length < full_total {
                        let end = (current + length + flanking_bases).min(full_total);
                        down_flank = full[(current + length) as usize..end as usize].replace(STOP_SYMBOL, "");
                    }
                }
                peptides.push(PeptideData { peptide: peptide.to_string(), up_flank, down_flank });
            }
            current += 1;
        }
    }
    peptides
}

/// PointMutationData.isRelevantMutation
pub fn is_relevant_mutation(worst_coding_effect: &str, canonical_effect: &str) -> bool {
    if worst_coding_effect == "NONSENSE_OR_FRAMESHIFT" {
        canonical_effect.contains("frameshift_variant") || canonical_effect.contains("stop_lost")
    } else {
        worst_coding_effect == "MISSENSE"
    }
}

/// PointMutationData.checkVariantEffects: the coding effect of a transcript impact's effects
pub fn check_variant_effects(effects: &str) -> &'static str {
    if effects.contains("frameshift_variant") || effects.contains("stop_lost") {
        "NONSENSE_OR_FRAMESHIFT"
    } else if effects.contains("missense_variant") {
        "MISSENSE"
    } else {
        "NONE"
    }
}

/// SvNeoEpitope.svIsNonDisruptiveInCodingTranscript
pub fn sv_is_non_disruptive_in_coding_transcript(positions: &[i32; 2], trans: &TranscriptData) -> bool {
    let (Some(cs), Some(ce)) = (trans.coding_start, trans.coding_end) else { return false };
    for w in trans.exons.windows(2) {
        let (exon, next) = (&w[0], &w[1]);
        if positions_within(positions[SE_START], positions[SE_END], exon.end, next.start) && positions_within(exon.end, next.start, cs, ce) {
            return true;
        }
    }
    false
}

