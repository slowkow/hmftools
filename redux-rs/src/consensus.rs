//! Consensus reads for duplicate groups (Java ConsensusReads, BaseBuilder, IndelConsensusReads, ReadParseState,
//! IlluminaRoutines, TemplateReads).

use crate::config::SequencingType;
use crate::dup::{DuplicateGroup, FragmentCoords};
use crate::rec::*;
use crate::unmap::extract_supp_alignment;
use rust_htslib::bam::{self, record::{Cigar, CigarString}};
use rust_htslib::faidx;
use std::collections::HashMap;

mod sbx_indel;
#[allow(unused_imports)]
pub use sbx_indel::{check_soft_clip_indel_pair, has_valid_cigar};

pub const CONSENSUS_MAX_DEPTH: usize = 100;
pub const CONSENSUS_PREFIX: &str = "CNS_";
pub const NO_BASE: u8 = 0;
pub const INVALID_BASE_QUAL: i32 = -1;
pub const BASE_QUAL_MINIMUM: i32 = 1;
pub const INVALID_POSITION: i32 = -1;
pub const UNSET_COUNT: i32 = -1;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ConsensusOutcome {
    Unset,
    AlignmentOnly,
    IndelMatch,
    IndelMismatch,
    IndelSoftclip,
    IndelFail,
    SingleRead,
}

pub const OUTCOME_NAMES: [&str; 7] = ["UNSET", "ALIGNMENT_ONLY", "INDEL_MATCH", "INDEL_MISMATCH", "INDEL_SOFTCLIP", "INDEL_FAIL", "SINGLE_READ"];

#[derive(Default, Clone)]
pub struct ConsensusStatistics {
    pub outcome_counts: [i64; 7],
}

impl ConsensusStatistics {
    pub fn register(&mut self, outcome: ConsensusOutcome) {
        self.outcome_counts[outcome as usize] += 1;
    }
    pub fn merge(&mut self, other: &ConsensusStatistics) {
        for i in 0..7 {
            self.outcome_counts[i] += other.outcome_counts[i];
        }
    }
    pub fn summary(&self) -> String {
        let mut parts = Vec::new();
        for i in 0..7 {
            if self.outcome_counts[i] > 0 {
                parts.push(format!("{}={}", OUTCOME_NAMES[i], self.outcome_counts[i]));
            }
        }
        parts.join(", ")
    }
}

/// Reference bases for one thread, read in blocks (Java consensus.RefGenome over a cached reference).
pub struct RefBases {
    reader: faidx::Reader,
    lengths: HashMap<String, i32>,
    chr: String,
    block_start: i32,
    block: Vec<u8>,
    #[cfg(test)]
    mock: Option<crate::tests::util::MockRefGenome>,
}

const REF_BLOCK: i32 = 1 << 16;

impl RefBases {
    pub fn open(path: &str) -> anyhow::Result<Self> {
        let reader = faidx::Reader::from_path(path)?;
        let mut lengths = HashMap::new();
        for n in reader.seq_names()? {
            let len = reader.fetch_seq_len(&n) as i32;
            lengths.insert(n, len);
        }
        Ok(RefBases {
            reader,
            lengths,
            chr: String::new(),
            block_start: 0,
            block: Vec::new(),
            #[cfg(test)]
            mock: None,
        })
    }
    /// Test seam: bases and lengths come from an in-memory genome; `placeholder_fasta` is only opened.
    #[cfg(test)]
    pub fn with_mock(placeholder_fasta: &str, mock: crate::tests::util::MockRefGenome) -> Self {
        let mut refs = Self::open(placeholder_fasta).expect("placeholder FASTA");
        refs.mock = Some(mock);
        refs
    }
    pub fn chromosome_length(&self, chr: &str) -> i32 {
        #[cfg(test)]
        if let Some(m) = &self.mock {
            return m.chromosome_length_or_zero(chr);
        }
        self.lengths.get(chr).copied().unwrap_or(0)
    }
    fn load(&mut self, chr: &str, pos: i32) {
        let len = self.chromosome_length(chr);
        let start = ((pos - 1) / REF_BLOCK) * REF_BLOCK + 1;
        let end = (start + REF_BLOCK - 1).min(len);
        self.block = self.reader.fetch_seq(chr, (start - 1) as usize, (end - 1) as usize).unwrap_or_default();
        self.chr.clear();
        self.chr.push_str(chr);
        self.block_start = start;
    }
    /// RefGenome.getRefBase: NO_BASE beyond the chromosome end
    pub fn ref_base(&mut self, chr: &str, pos: i32) -> u8 {
        // a position in the loaded block is within the chromosome
        if pos >= self.block_start && pos < self.block_start + self.block.len() as i32 && self.chr == chr {
            return self.block[(pos - self.block_start) as usize];
        }
        if pos > self.chromosome_length(chr) {
            return NO_BASE;
        }
        #[cfg(test)]
        if let Some(m) = &self.mock {
            return m.bases(chr, pos, pos)[0];
        }
        self.load(chr, pos);
        self.block[(pos - self.block_start) as usize]
    }
    /// RefGenome.getRefBases: None outside the chromosome
    pub fn ref_bases(&mut self, chr: &str, start: i32, end: i32) -> Option<Vec<u8>> {
        if start < 1 || end < start || end > self.chromosome_length(chr) {
            return None;
        }
        #[cfg(test)]
        if let Some(m) = &self.mock {
            return Some(m.bases(chr, start, end));
        }
        let mut out = Vec::with_capacity((end - start + 1) as usize);
        for pos in start..=end {
            out.push(self.ref_base(chr, pos));
        }
        Some(out)
    }
}

/// BaseQualAdjustment.adjustBaseQual with the standard values 1, 11, 25, 37
pub fn adjust_base_qual(base_qual: f64) -> i32 {
    const STANDARD: [i32; 4] = [BASE_QUAL_MINIMUM, 11, 25, 37];
    const PERMITTED_DIFF_MAX: f64 = 1.5;
    if base_qual < STANDARD[1] as f64 - PERMITTED_DIFF_MAX {
        return STANDARD[0];
    }
    for i in (1..STANDARD.len()).rev() {
        if base_qual > STANDARD[i] as f64 {
            return STANDARD[i];
        }
        let diff = STANDARD[i] as f64 - base_qual;
        if diff <= PERMITTED_DIFF_MAX {
            return STANDARD[i];
        }
    }
    STANDARD[0]
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct BaseQualPair {
    pub base: u8,
    pub qual: i32,
}

const INVALID: BaseQualPair = BaseQualPair { base: NO_BASE, qual: INVALID_BASE_QUAL };

impl BaseQualPair {
    fn is_valid(&self) -> bool {
        self.base != NO_BASE && self.qual > 0
    }
}

fn create_base_qual_pair(base: u8, qual: f64) -> BaseQualPair {
    BaseQualPair { base, qual: adjust_base_qual(qual) }
}

/// A read of a duplicate group with its bases decoded.
pub struct ReadView<'a> {
    pub rec: &'a bam::Record,
    pub bases: Vec<u8>,
    pub quals: &'a [u8],
    pub cigar: &'a [u32],
    pub start: i32,
    pub end: i32,
}

impl<'a> ReadView<'a> {
    pub fn new(rec: &'a bam::Record) -> Self {
        ReadView { rec, bases: rec.seq().as_bytes(), quals: rec.qual(), cigar: rec.raw_cigar(), start: start(rec), end: end(rec) }
    }
    fn unclipped_start(&self) -> i32 {
        self.start - left_soft_clip(self.cigar)
    }
    fn unclipped_end(&self) -> i32 {
        self.end + right_soft_clip(self.cigar)
    }
}

pub struct ConsensusState {
    pub is_forward: bool,
    pub chromosome: String,
    pub bases: Vec<u8>,
    pub base_qualities: Vec<i32>,
    pub cigar: Vec<u32>,
    pub unclipped_pos_start: i32,
    pub unclipped_pos_end: i32,
    pub alignment_start: i32,
    pub alignment_end: i32,
    pub map_quality: u8,
    pub num_mutations: i32,
    pub outcome: ConsensusOutcome,
    current_len: i32,
    current_op: u32,
}

impl ConsensusState {
    fn new(is_forward: bool, chromosome: &str) -> Self {
        ConsensusState {
            is_forward,
            chromosome: chromosome.to_string(),
            bases: Vec::new(),
            base_qualities: Vec::new(),
            cigar: Vec::new(),
            unclipped_pos_start: 0,
            unclipped_pos_end: 0,
            alignment_start: 0,
            alignment_end: 0,
            map_quality: 0,
            num_mutations: 0,
            outcome: ConsensusOutcome::Unset,
            current_len: 0,
            current_op: u32::MAX,
        }
    }
    fn set_base_length(&mut self, n: usize) {
        self.bases = vec![0; n];
        self.base_qualities = vec![0; n];
    }
    fn set_boundaries(&mut self, read: &ReadView) {
        let (us, ue) = (read.unclipped_start(), read.unclipped_end());
        if self.unclipped_pos_start == 0 {
            self.unclipped_pos_start = us;
            self.unclipped_pos_end = ue;
            self.alignment_start = read.start;
            self.alignment_end = read.end;
        } else {
            self.unclipped_pos_start = us.min(self.unclipped_pos_start);
            self.unclipped_pos_end = ue.max(self.unclipped_pos_end);
            self.alignment_start = read.start.min(self.alignment_start);
            self.alignment_end = read.end.max(self.alignment_end);
        }
    }
    /// setFromRead(read, false): the template's length, boundaries and cigar
    fn set_from_read(&mut self, read: &ReadView) {
        self.set_base_length(read.quals.len());
        self.set_boundaries(read);
        self.cigar.extend_from_slice(read.cigar);
    }
    fn add_cigar_element(&mut self, len: i32, op: u32) {
        if self.current_len > 0 && self.current_op != op {
            self.add_current_cigar_element();
        }
        self.current_len += len;
        self.current_op = op;
    }
    fn add_current_cigar_element(&mut self) {
        if self.current_len == 0 {
            return;
        }
        let element = (self.current_len as u32) << 4 | self.current_op;
        if self.is_forward {
            self.cigar.push(element);
        } else {
            self.cigar.insert(0, element);
        }
        self.current_len = 0;
        self.current_op = u32::MAX;
    }
    fn finalise_cigar(&mut self) {
        self.add_current_cigar_element();
    }
    /// ConsensusState.setNumMutations
    fn set_num_mutations(&mut self, refs: &mut RefBases) {
        self.num_mutations = 0;
        let Some(ref_bases) = refs.ref_bases(&self.chromosome, self.alignment_start, self.alignment_end) else { return };
        let mut base_index = 0usize;
        let mut ref_index = 0usize;
        for &c in &self.cigar {
            let (o, n) = (op(c), oplen(c) as usize);
            if matches!(o, OP_I | OP_D | OP_X) {
                self.num_mutations += n as i32;
            } else if o == OP_M {
                for i in 0..n {
                    if ref_bases[ref_index + i] != self.bases[base_index + i] {
                        self.num_mutations += 1;
                    }
                }
            }
            if consumes_read(o) {
                base_index += n;
            }
            if consumes_ref(o) {
                ref_index += n;
            }
        }
    }
}

#[inline]
fn delete_or_split(o: u32) -> bool {
    o == OP_D || o == OP_N
}
#[inline]
fn aligned_or_clipped(o: u32) -> bool {
    o == OP_M || o == OP_S || o == OP_H
}
#[inline]
fn consumes_ref_or_unclipped(o: u32) -> bool {
    consumes_ref(o) || o == OP_S
}

/// ReadParseState: a cursor over one read that follows the consensus cigar
pub(crate) struct ReadParseState<'a, 'b> {
    read: &'b ReadView<'a>,
    is_forward: bool,
    element_count: i32,
    read_index: i32,
    cigar_index: i32,
    element_index: i32,
    pub(crate) element_length: i32,
    pub(crate) element_type: u32,
    active: bool,
    pub(crate) exhausted: bool,
    current_ref_position: i32,
}

impl<'a, 'b> ReadParseState<'a, 'b> {
    pub(crate) fn new(read: &'b ReadView<'a>, is_forward: bool) -> Self {
        let element_count = read.cigar.len() as i32;
        let (read_index, cigar_index) = if is_forward { (0, 0) } else { (read.bases.len() as i32 - 1, element_count - 1) };
        let c = read.cigar[cigar_index as usize];
        let (element_length, element_type) = (oplen(c), op(c));
        let mut pos = if is_forward { read.start } else { read.end };
        if element_type == OP_S {
            if is_forward {
                pos -= element_length;
            } else {
                pos += element_length;
            }
        }
        ReadParseState {
            read,
            is_forward,
            element_count,
            read_index,
            cigar_index,
            element_index: 0,
            element_length,
            element_type,
            active: false,
            exhausted: false,
            current_ref_position: pos,
        }
    }
    pub(crate) fn base(&self) -> u8 {
        if self.exhausted { NO_BASE } else { self.read.bases[self.read_index as usize] }
    }
    pub(crate) fn base_qual(&self) -> i32 {
        if self.exhausted { 0 } else { self.read.quals[self.read_index as usize] as i8 as i32 }
    }
    fn before_unclipped_position(&self, ref_position: i32) -> bool {
        if self.active {
            return false;
        }
        if self.is_forward { self.current_ref_position > ref_position } else { self.current_ref_position < ref_position }
    }
    pub(crate) fn move_next_base(&mut self) {
        if self.exhausted {
            return;
        }
        self.active = true;
        self.element_index += 1;
        let mut skips_first_element = false;
        if self.element_index >= self.element_length {
            skips_first_element = !consumes_read(self.element_type) && (self.cigar_index == 0 || self.cigar_index == self.element_count - 1);
            self.cigar_index += if self.is_forward { 1 } else { -1 };
            if self.cigar_index < 0 || self.cigar_index >= self.element_count {
                self.exhausted = true;
                self.active = false;
                return;
            }
            let c = self.read.cigar[self.cigar_index as usize];
            self.element_length = oplen(c);
            self.element_type = op(c);
            self.element_index = 0;
        }
        if !skips_first_element {
            let step = if self.is_forward { 1 } else { -1 };
            if consumes_read(self.element_type) {
                self.read_index += step;
            }
            if consumes_ref_or_unclipped(self.element_type) {
                self.current_ref_position += step;
            }
        }
    }
    fn skip_insert_element(&mut self) {
        while self.element_type == OP_I && !self.exhausted {
            self.move_next_base();
        }
    }
    fn skip_non_read_base_element(&mut self) {
        if self.exhausted {
            return;
        }
        let skips_first_element = !consumes_read(self.element_type) && (self.cigar_index == 0 || self.cigar_index == self.element_count - 1);
        self.cigar_index += if self.is_forward { 1 } else { -1 };
        if self.cigar_index < 0 || self.cigar_index >= self.element_count {
            self.exhausted = true;
            return;
        }
        let c = self.read.cigar[self.cigar_index as usize];
        self.element_length = oplen(c);
        self.element_type = op(c);
        self.element_index = 0;
        let step = if self.is_forward { 1 } else { -1 };
        if !skips_first_element && consumes_read(self.element_type) {
            self.read_index += step;
        }
        self.current_ref_position += step;
    }
    fn move_to_ref_position(&mut self, target: i32) {
        if self.current_ref_position == target {
            self.active = true;
        }
        if (self.is_forward && self.current_ref_position >= target) || (!self.is_forward && self.current_ref_position <= target) {
            return;
        }
        while self.current_ref_position != target && !self.exhausted {
            self.move_next_base();
        }
    }
}

/// IlluminaRoutines.isDualStrandAndIsFirstInPair
fn is_dual_strand_and_first_in_pair(reads: &[ReadView], first_out: &mut [bool]) -> bool {
    if !has(reads[0].rec, FLAG_PAIRED) {
        return false;
    }
    let mut dual = false;
    for (i, r) in reads.iter().enumerate() {
        first_out[i] = has(r.rec, FLAG_FIRST);
        if i > 0 {
            dual |= first_out[0] != first_out[i];
        }
    }
    dual
}

pub struct ConsensusReads {
    refs: RefBases,
    sequencing_type: SequencingType,
    chromosome_length: i32,
    max_consensus_read_count: usize,
    pub stats: ConsensusStatistics,
    /// Ultima: the name of a consensus read that was built with a base that has no value (see `create_consensus_read`)
    pub undetermined_base_read: Option<String>,
    /// the index in the group's reads of the read that the last consensus read is a copy of, when it is one
    /// (ConsensusReads.buildFromRead)
    pub built_from_read: Option<usize>,
    /// addresses of the group's reads that Java compared by value while forming the last consensus read
    sbx_eager_decoded: Vec<usize>,
    /// ReduxConfig.RunChecks (-run_checks): an invalid SBX consensus read is replaced by a copy of the template
    pub run_checks: bool,
}

impl ConsensusReads {
    pub fn new(ref_genome_file: &str, sequencing_type: SequencingType) -> anyhow::Result<Self> {
        Ok(ConsensusReads {
            refs: RefBases::open(ref_genome_file)?,
            sequencing_type,
            chromosome_length: 0,
            max_consensus_read_count: if sequencing_type == SequencingType::Sbx { crate::sbx::SBX_CONSENSUS_MAX_DEPTH } else { CONSENSUS_MAX_DEPTH },
            stats: ConsensusStatistics::default(),
            undetermined_base_read: None,
            built_from_read: None,
            sbx_eager_decoded: Vec::new(),
            run_checks: false,
        })
    }

    /// Test seam: a ConsensusReads over given reference bases.
    #[cfg(test)]
    pub fn with_refs(refs: RefBases, sequencing_type: SequencingType) -> Self {
        ConsensusReads { refs, sequencing_type, chromosome_length: 0, max_consensus_read_count: CONSENSUS_MAX_DEPTH, built_from_read: None, sbx_eager_decoded: Vec::new(), run_checks: false, undetermined_base_read: None, stats: ConsensusStatistics::default() }
    }

    pub fn set_chromosome_length(&mut self, len: i32) {
        self.chromosome_length = len;
    }
    pub fn refs(&mut self) -> &mut RefBases {
        &mut self.refs
    }

    /// DuplicateGroup.formConsensusRead: builds the consensus read and sets its CR and UT attributes
    pub fn form_consensus_read_plain(&mut self, group: &mut DuplicateGroup, names: &[String]) {
        let mut consensus = self.create_consensus_read(&mut group.reads, &group.coords, group.umi.as_deref(), names);
        if !self.sbx_eager_decoded.is_empty() {
            for read in group.reads.iter_mut() {
                if self.sbx_eager_decoded.contains(&(read as *const bam::Record as usize)) {
                    crate::sbx::eager_decode(read);
                }
            }
            self.sbx_eager_decoded.clear();
        }
        let consensus_type = if group.dual_strand { "DUAL" } else { "SINGLE" };
        let pcr_cluster_count = if self.sequencing_type == SequencingType::Illumina { calculate_pcr_cluster_count(group) } else { UNSET_COUNT };
        let first_in_pair_count = group.first_in_pair_count();
        let value = if pcr_cluster_count == UNSET_COUNT {
            format!("{};{}", group.total_read_count(), first_in_pair_count)
        } else {
            format!("{};{};{}", group.total_read_count(), first_in_pair_count, pcr_cluster_count)
        };
        let mut aux = aux_bytes(&consensus).to_vec();
        aux_block_set(&mut aux, b"CR", Some(&encode_z(value.as_bytes())));
        aux_block_set(&mut aux, b"UT", Some(&encode_z(consensus_type.as_bytes())));
        set_aux_bytes(&mut consensus, &aux);
        group.consensus_read = Some(consensus);
    }

    /// ConsensusReads.createConsensusRead
    pub(crate) fn create_consensus_read(&mut self, reads: &mut Vec<bam::Record>, coords: &FragmentCoords, umi: Option<&str>, names: &[String]) -> bam::Record {
        let template_index = select_template_read(reads, coords);
        let consensus_read_id = form_consensus_read_id(reads[template_index].qname(), umi);
        self.built_from_read = None;
        if reads.len() <= 1 || has(&reads[0], FLAG_UNMAPPED) {
            self.built_from_read = Some(template_index);
            let template = &reads[template_index];
            let first = !has(template, FLAG_PAIRED) || has(template, FLAG_FIRST);
            return build_from_read(template, &consensus_read_id, first);
        }

        // the first reads up to the maximum, plus the template since it drives the cigar selection
        let mut view_indices: Vec<usize> = if reads.len() < self.max_consensus_read_count {
            (0..reads.len()).collect()
        } else {
            (0..self.max_consensus_read_count).collect()
        };
        if !view_indices.contains(&template_index) {
            view_indices.push(template_index);
        }
        let views: Vec<ReadView> = view_indices.iter().map(|&i| ReadView::new(&reads[i])).collect();
        let template_view_index = view_indices.iter().position(|&i| i == template_index).unwrap();
        let template = &views[template_view_index];

        let is_forward = !has(template.rec, FLAG_REVERSE);
        let chromosome = &names[template.rec.tid() as usize];
        let mut state = ConsensusState::new(is_forward, chromosome);
        let mut has_indels = false;
        for r in &views {
            has_indels |= r.cigar.iter().any(|c| matches!(op(*c), OP_I | OP_D));
            state.map_quality = state.map_quality.max(r.rec.mapq());
        }

        if has_indels {
            self.build_indel_components(&views, &mut state, template_view_index, names);
            if state.outcome == ConsensusOutcome::IndelFail {
                self.stats.register(ConsensusOutcome::IndelFail);
                let primary = select_primary_read(&views);
                self.built_from_read = Some(view_indices[primary]);
                let first = !has(template.rec, FLAG_PAIRED) || has(template.rec, FLAG_FIRST);
                return build_from_read(views[primary].rec, &consensus_read_id, first);
            }
        } else {
            state.set_from_read(template);
            self.build_read_bases(&views, &mut state, names);
            state.outcome = ConsensusOutcome::AlignmentOnly;
        }
        self.stats.register(state.outcome);
        state.set_num_mutations(&mut self.refs);
        // BaseBuilder.determineBaseAndQual has no rule for Ultima and returns BaseQualPair.INVALID where the reads
        // disagree; Java then fails when it writes the read (the caller reports it)
        if self.sequencing_type == SequencingType::Ultima && state.bases.contains(&NO_BASE) {
            self.undetermined_base_read = Some(String::from_utf8_lossy(&consensus_read_id).to_string());
        }
        let mut consensus = create_consensus_record(&state, template.rec, &consensus_read_id);
        if self.sequencing_type == SequencingType::Sbx {
            // the first duplex base from the 5' end, for SbxRoutines.finaliseRead and the tools downstream
            let quals: Vec<&[u8]> = views.iter().map(|r| r.quals).collect();
            let first_duplex_base_index = crate::sbx::find_max_duplex_base_index(has(views[0].rec, FLAG_REVERSE), &quals);
            if first_duplex_base_index >= 0 {
                set_attribute(&mut consensus, b"YX", Some(&encode_int(first_duplex_base_index as i64)));
            }
            if self.run_checks {
                // ReadValidReason.isValidRead
                let cigar = consensus.raw_cigar();
                let cigar_base_length: i32 = cigar.iter().filter(|c| consumes_read(op(**c))).map(|c| oplen(*c)).sum();
                let reason = if cigar_base_length != state.bases.len() as i32 {
                    Some("CIGAR_LENGTH")
                } else if cigar.is_empty() || !has_valid_cigar(cigar) {
                    Some("CIGAR_ELEMENTS")
                } else if state.bases.contains(&NO_BASE) {
                    Some("INVALID_BASE")
                } else {
                    None
                };
                if let Some(reason) = reason {
                    crate::warn!(
                        "invalid consensus read({}): readCount({}) {} read: {}",
                        reason, views.len(), if state.is_forward { "forward" } else { "reverse" }, String::from_utf8_lossy(&consensus_read_id)
                    );
                    self.built_from_read = Some(template_index);
                    let first = !has(template.rec, FLAG_PAIRED) || has(template.rec, FLAG_FIRST);
                    return build_from_read(template.rec, &consensus_read_id, first);
                }
            }
        }
        consensus
    }

    /// BaseBuilder.buildReadBases: reads compared base by base from the consensus unclipped start, no indel handling
    fn build_read_bases(&mut self, reads: &[ReadView], state: &mut ConsensusState, names: &[String]) {
        let chromosome = &names[reads[0].rec.tid() as usize];
        let mut chromosome_length = self.chromosome_length;
        if chromosome_length == 0 {
            chromosome_length = self.refs.chromosome_length(chromosome);
        }
        let base_length = state.bases.len();
        let read_count = reads.len();
        let mut is_first_in_pair = vec![false; read_count];
        let is_dual_strand = self.sequencing_type == SequencingType::Illumina && is_dual_strand_and_first_in_pair(reads, &mut is_first_in_pair);
        let offsets: Vec<i32> = reads.iter().map(|r| r.unclipped_start() - state.unclipped_pos_start).collect();
        let mut location_bases = vec![NO_BASE; read_count];
        let mut location_quals = vec![INVALID_BASE_QUAL; read_count];
        for base_index in 0..base_length {
            let mut has_mismatch = false;
            let mut max_qual: i32 = -1;
            let mut uncertain_qual_count = 0;
            let mut first_base = NO_BASE;
            for r in 0..read_count {
                location_bases[r] = NO_BASE;
                location_quals[r] = INVALID_BASE_QUAL;
                let read_index = base_index as i32 - offsets[r];
                if read_index < 0 || read_index >= reads[r].bases.len() as i32 {
                    continue;
                }
                location_bases[r] = reads[r].bases[read_index as usize];
                location_quals[r] = reads[r].quals[read_index as usize] as i8 as i32;
                if first_base == NO_BASE {
                    first_base = location_bases[r];
                } else {
                    has_mismatch |= location_bases[r] != first_base;
                }
                max_qual = max_qual.max(location_quals[r]);
                if location_quals[r] <= BASE_QUAL_MINIMUM {
                    uncertain_qual_count += 1;
                }
            }
            if !has_mismatch {
                state.bases[base_index] = first_base;
                state.base_qualities[base_index] = if uncertain_qual_count as f64 >= 0.5 * read_count as f64 { BASE_QUAL_MINIMUM } else { max_qual };
            } else {
                let mut base_position = state.unclipped_pos_start + base_index as i32;
                if base_position < 1 || base_position > chromosome_length {
                    base_position = INVALID_POSITION;
                }
                let pair = self.determine_base_and_qual(&location_bases, &location_quals, chromosome, base_position, is_dual_strand, &is_first_in_pair);
                state.bases[base_index] = pair.base;
                state.base_qualities[base_index] = pair.qual;
            }
        }
    }

    /// BaseBuilder.determineBaseAndQual
    fn determine_base_and_qual(
        &mut self, bases: &[u8], quals: &[i32], chromosome: &str, position: i32, is_dual_strand: bool, is_first_in_pair: &[bool],
    ) -> BaseQualPair {
        match self.sequencing_type {
            SequencingType::Illumina => {
                if is_dual_strand && position != INVALID_POSITION {
                    determine_dual_strand_base_and_qual(is_first_in_pair, bases, quals, chromosome, position, &mut self.refs)
                } else {
                    determine_base_and_qual(bases, quals, chromosome, position, &mut self.refs)
                }
            }
            SequencingType::Sbx => crate::sbx::determine_base_and_qual(bases, quals, chromosome, position, &mut self.refs),
            // Ultima has no rule in Java either
            _ => INVALID,
        }
    }

    /// IndelConsensusReads.buildIndelComponents
    fn build_indel_components(&mut self, reads: &[ReadView], state: &mut ConsensusState, template_index: usize, names: &[String]) {
        let template = &reads[template_index];
        let has_cigar_mismatch = reads.iter().any(|r| r.cigar != template.cigar);
        if !has_cigar_mismatch {
            state.set_from_read(template);
            self.build_read_bases(reads, state, names);
            state.outcome = ConsensusOutcome::IndelMatch;
            return;
        }
        if self.sequencing_type == SequencingType::Sbx {
            self.sbx_determine_indel_consensus(reads, state, names);
            return;
        }
        let read_count = reads.len();
        let mut is_first_in_pair = vec![false; read_count];
        let is_dual_strand = is_dual_strand_and_first_in_pair(reads, &mut is_first_in_pair);
        let base_length = template.bases.len();
        state.set_base_length(base_length);
        state.set_boundaries(template);
        let mut read_states: Vec<ReadParseState> = reads.iter().map(|r| ReadParseState::new(r, state.is_forward)).collect();
        let mut base_index: i32 = if state.is_forward { 0 } else { base_length as i32 - 1 };
        let mut ref_position = if state.is_forward { state.unclipped_pos_start } else { state.unclipped_pos_end };
        for rs in read_states.iter_mut() {
            rs.move_to_ref_position(ref_position);
        }
        let cigar_count = template.cigar.len() as i32;
        let mut cigar_index: i32 = if state.is_forward { 0 } else { cigar_count - 1 };
        let chromosome_length = if self.chromosome_length == 0 {
            self.refs.chromosome_length(&names[reads[0].rec.tid() as usize])
        } else {
            self.chromosome_length
        };
        while cigar_index >= 0 && cigar_index < cigar_count {
            let element = template.cigar[cigar_index as usize];
            self.add_element_bases(state, &mut read_states, element, base_index, ref_position, is_dual_strand, &is_first_in_pair, chromosome_length);
            if state.outcome == ConsensusOutcome::IndelFail {
                break;
            }
            let (o, n) = (op(element), oplen(element));
            let dir = if state.is_forward { 1 } else { -1 };
            if consumes_read(o) {
                base_index += dir * n;
            }
            if consumes_ref_or_unclipped(o) {
                ref_position += dir * n;
            }
            cigar_index += dir;
        }
        state.finalise_cigar();
        if state.outcome != ConsensusOutcome::IndelFail {
            state.outcome = ConsensusOutcome::IndelMismatch;
        }
    }

    /// IndelConsensusReads.addElementBases
    #[allow(clippy::too_many_arguments)]
    fn add_element_bases(
        &mut self, state: &mut ConsensusState, read_states: &mut [ReadParseState], element: u32, mut base_index: i32, ref_position: i32,
        is_dual_strand: bool, is_first_in_pair: &[bool], chromosome_length: i32,
    ) {
        let read_count = read_states.len();
        let (sel_op, sel_len) = (op(element), oplen(element));
        state.add_cigar_element(sel_len, sel_op);
        let dir = if state.is_forward { 1 } else { -1 };
        if delete_or_split(sel_op) || sel_op == OP_H {
            let new_ref_position = ref_position + dir * sel_len;
            for read in read_states.iter_mut() {
                if read.exhausted || read.before_unclipped_position(new_ref_position) {
                    continue;
                }
                if read.element_type == sel_op && read.element_index == 0 && read.element_length == sel_len {
                    read.skip_non_read_base_element();
                    continue;
                }
                for _ in 0..sel_len {
                    if read.element_type == OP_I {
                        read.skip_insert_element();
                    }
                    if delete_or_split(read.element_type) {
                        read.move_next_base();
                    } else if aligned_or_clipped(read.element_type) {
                        if read.current_ref_position == new_ref_position {
                            break;
                        }
                        read.move_next_base();
                    }
                }
            }
            return;
        }
        let mut location_bases = vec![NO_BASE; read_count];
        let mut location_quals = vec![0i32; read_count];
        let mut new_ref_position = ref_position;
        for i in 0..sel_len {
            let mut has_mismatch = false;
            let mut max_qual = BASE_QUAL_MINIMUM;
            let mut first_base = NO_BASE;
            if i > 0 {
                for r in 0..read_count {
                    location_bases[r] = NO_BASE;
                    location_quals[r] = 0;
                }
            }
            for r in 0..read_count {
                let read = &mut read_states[r];
                if read.exhausted || read.before_unclipped_position(new_ref_position) {
                    continue;
                }
                if consumes_ref(sel_op) && read.element_type == OP_I {
                    read.skip_insert_element();
                }
                let mut use_base = true;
                let mut move_next = true;
                if read.element_type != sel_op {
                    if sel_op == OP_M {
                        if delete_or_split(read.element_type) || read.element_type == OP_S {
                            use_base = false;
                        } else if read.element_type == OP_I {
                            state.outcome = ConsensusOutcome::IndelFail;
                            return;
                        }
                    } else if sel_op == OP_I && (read.element_type == OP_M || delete_or_split(read.element_type) || read.element_type == OP_S) {
                        move_next = false;
                        use_base = false;
                    }
                } else if sel_op == OP_I && read.element_length != sel_len && self.sequencing_type == SequencingType::Sbx {
                    // SBX only: the bases of an insert of another length are not used
                    use_base = false;
                }
                if use_base {
                    location_bases[r] = read.base();
                    location_quals[r] = read.base_qual();
                    if first_base == NO_BASE {
                        first_base = location_bases[r];
                    } else {
                        has_mismatch |= location_bases[r] != first_base;
                    }
                    max_qual = max_qual.max(location_quals[r]);
                }
                if move_next {
                    read.move_next_base();
                }
            }
            if !has_mismatch && first_base != NO_BASE {
                state.bases[base_index as usize] = first_base;
                state.base_qualities[base_index as usize] = max_qual;
            } else {
                let base_position = if new_ref_position < 1 || new_ref_position > chromosome_length { INVALID_POSITION } else { new_ref_position };
                let chromosome = state.chromosome.clone();
                let pair = self.determine_base_and_qual(&location_bases, &location_quals, &chromosome, base_position, is_dual_strand, is_first_in_pair);
                state.bases[base_index as usize] = pair.base;
                state.base_qualities[base_index as usize] = pair.qual;
            }
            base_index += dir;
            if consumes_ref_or_unclipped(sel_op) {
                new_ref_position += dir;
            }
        }
    }
}

/// IlluminaRoutines.checkCommonBaseAndQual: one read, or two reads; None when the general rule applies
fn check_common_base_and_qual(bases: &[u8], quals: &[i32], position: i32, chromosome: &str, refs: &mut RefBases) -> Option<BaseQualPair> {
    if bases.len() == 1 {
        return Some(BaseQualPair { base: bases[0], qual: quals[0] });
    }
    if bases.len() == 2 {
        let max_qual = quals[0].max(quals[1]);
        if bases[0] == bases[1] {
            return Some(create_base_qual_pair(bases[0], max_qual as f64));
        }
        let min_qual = quals[0].min(quals[1]);
        let calc_qual = max_qual as f64 * BASE_QUAL_MINIMUM.max(max_qual - min_qual) as f64 / max_qual as f64;
        if quals[0] > quals[1] {
            return Some(create_base_qual_pair(bases[0], calc_qual));
        }
        if quals[1] > quals[0] {
            return Some(create_base_qual_pair(bases[1], calc_qual));
        }
        if position != INVALID_POSITION {
            let ref_base = refs.ref_base(chromosome, position);
            if bases[0] == ref_base {
                return Some(create_base_qual_pair(bases[0], calc_qual));
            }
            if bases[1] == ref_base {
                return Some(create_base_qual_pair(bases[1], calc_qual));
            }
        } else {
            return Some(create_base_qual_pair(bases[0], calc_qual));
        }
    }
    None
}

/// IlluminaRoutines.determineBaseAndQual
pub fn determine_base_and_qual(bases: &[u8], quals: &[i32], chromosome: &str, position: i32, refs: &mut RefBases) -> BaseQualPair {
    if let Some(p) = check_common_base_and_qual(bases, quals, position, chromosome, refs) {
        return p;
    }
    let mut distinct: Vec<u8> = Vec::with_capacity(4);
    let mut qual_totals: Vec<i32> = Vec::with_capacity(4);
    for i in 0..bases.len() {
        if bases[i] == NO_BASE {
            continue;
        }
        match distinct.iter().position(|b| *b == bases[i]) {
            Some(j) => qual_totals[j] += quals[i],
            None => {
                distinct.push(bases[i]);
                qual_totals.push(quals[i]);
            }
        }
    }
    if distinct.is_empty() {
        return INVALID;
    }
    let mut max_base = distinct[0];
    let mut max_is_ref = false;
    let mut max_qual_total = qual_totals[0];
    for i in 1..distinct.len() {
        if qual_totals[i] > max_qual_total {
            max_qual_total = qual_totals[i];
            max_base = distinct[i];
        } else if qual_totals[i] >= max_qual_total && !max_is_ref && position != INVALID_POSITION {
            let ref_base = refs.ref_base(chromosome, position);
            if max_base == ref_base {
                max_is_ref = true;
            } else if distinct[i] == ref_base {
                max_qual_total = qual_totals[i];
                max_base = distinct[i];
                max_is_ref = true;
            }
        }
    }
    let mut select_quals: Vec<i32> = (0..bases.len()).filter(|&i| bases[i] == max_base).map(|i| quals[i]).collect();
    select_quals.sort();
    let median = select_quals[select_quals.len() / 2];
    let mut differing = 0;
    for i in 0..distinct.len() {
        if distinct[i] != max_base {
            differing += qual_totals[i];
        }
    }
    let calc_qual = median as f64 * BASE_QUAL_MINIMUM.max(max_qual_total - differing) as f64 / max_qual_total as f64;
    create_base_qual_pair(max_base, calc_qual)
}

/// IlluminaRoutines.determineDualStrandBaseAndQual: a consensus per strand, then the two compared
fn determine_dual_strand_base_and_qual(
    is_first_in_pair: &[bool], bases: &[u8], quals: &[i32], chromosome: &str, position: i32, refs: &mut RefBases,
) -> BaseQualPair {
    let mut bases_first = Vec::new();
    let mut quals_first = Vec::new();
    let mut bases_second = Vec::new();
    let mut quals_second = Vec::new();
    for i in 0..bases.len() {
        if is_first_in_pair[i] {
            bases_first.push(bases[i]);
            quals_first.push(quals[i]);
        } else {
            bases_second.push(bases[i]);
            quals_second.push(quals[i]);
        }
    }
    let first = determine_base_and_qual(&bases_first, &quals_first, chromosome, position, refs);
    let second = determine_base_and_qual(&bases_second, &quals_second, chromosome, position, refs);
    if !first.is_valid() && !second.is_valid() {
        return INVALID;
    }
    if !first.is_valid() {
        return second;
    }
    if !second.is_valid() {
        return first;
    }
    if first.base == second.base {
        return create_base_qual_pair(first.base, first.qual.max(second.qual) as f64);
    }
    let ref_base = refs.ref_base(chromosome, position);
    let first_is_ref = first.base == ref_base;
    let second_is_ref = second.base == ref_base;
    if !first_is_ref && !second_is_ref {
        let (max_base, max_qual, differing) = if first.qual >= second.qual { (first.base, first.qual, second.qual) } else { (second.base, second.qual, first.qual) };
        let qual = 0.max(max_qual - differing) as i8;
        return create_base_qual_pair(max_base, qual as f64);
    }
    let (ref_qual, differing) = if first_is_ref { (first.qual, second.qual) } else { (second.qual, first.qual) };
    let qual = BASE_QUAL_MINIMUM.max(ref_qual - differing) as i8;
    create_base_qual_pair(ref_base, qual as f64)
}

/// IndelConsensusReads.selectPrimaryRead: most aligned bases, then highest average base quality
fn select_primary_read(reads: &[ReadView]) -> usize {
    let mut max_aligned = 0;
    let mut max_qual = 0.0;
    let mut max_read = 0;
    let mut found = false;
    for (i, r) in reads.iter().enumerate() {
        let aligned: i32 = r.cigar.iter().filter(|c| op(**c) == OP_M).map(|c| oplen(*c)).sum();
        if aligned > max_aligned {
            max_read = i;
            max_aligned = aligned;
            max_qual = calc_base_qual_average(r.quals);
            found = true;
        } else if aligned == max_aligned {
            let avg = calc_base_qual_average(r.quals);
            if avg > max_qual {
                max_read = i;
                max_qual = avg;
                found = true;
            }
        }
    }
    let _ = found;
    max_read
}

/// DuplicateGroupBuilder.calcBaseQualAverage
pub fn calc_base_qual_average(quals: &[u8]) -> f64 {
    if quals.is_empty() {
        return 0.0;
    }
    let total: i32 = quals.iter().map(|q| *q as i8 as i32).sum();
    total as f64 / quals.len() as f64
}

/// TemplateReads.selectTemplateRead: the index of the template. Sorts `reads` by name when exactly half are primary.
pub fn select_template_read(reads: &mut Vec<bam::Record>, coords: &FragmentCoords) -> usize {
    let primary_count = reads.iter().filter(|r| !has(r, FLAG_SUPPLEMENTARY)).count();
    let are_primaries;
    if primary_count == reads.len() / 2 {
        reads.sort_by(|a, b| a.qname().cmp(b.qname()));
        are_primaries = !has(&reads[0], FLAG_SUPPLEMENTARY);
    } else {
        are_primaries = primary_count > reads.len() / 2;
    }
    let is_lower = coords.read_is_lower;
    // (lower cigar, mate cigar) per read
    let infos: Vec<(String, String)> = reads
        .iter()
        .map(|r| {
            let read_cigar = if are_primaries {
                cigar_string(r.raw_cigar())
            } else {
                let sd = extract_supp_alignment(aux_str(r, b"SA")).expect("supplementary read without supplementary alignment data");
                String::from_utf8_lossy(sd.cigar).to_string()
            };
            let mate_cigar = aux_str(r, b"MC").map(|m| String::from_utf8_lossy(m).to_string()).unwrap_or_else(|| "*".to_string());
            if is_lower { (read_cigar, mate_cigar) } else { (mate_cigar, read_cigar) }
        })
        .collect();
    let all: Vec<usize> = (0..reads.len()).collect();
    let top_lower = find_top_frequency_cigar(&all, &infos, false);
    let filtered: Vec<usize> = all.into_iter().filter(|&i| infos[i].0 == top_lower).collect();
    let top_mate = find_top_frequency_cigar(&filtered, &infos, true);
    let mut filtered: Vec<usize> = filtered.into_iter().filter(|&i| infos[i].1 == top_mate).collect();
    if filtered.len() > 1 {
        filtered.sort_by(|&a, &b| reads[a].qname().cmp(reads[b].qname()));
    }
    filtered[0]
}

fn find_top_frequency_cigar(indices: &[usize], infos: &[(String, String)], use_mate: bool) -> String {
    let get = |i: usize| if use_mate { &infos[i].1 } else { &infos[i].0 };
    let mut freq: HashMap<&String, i32> = HashMap::new();
    for &i in indices {
        *freq.entry(get(i)).or_insert(0) += 1;
    }
    if freq.len() == 1 {
        return get(indices[0]).clone();
    }
    let max = freq.values().copied().max().unwrap_or(0);
    let mut top: Vec<&String> = freq.iter().filter(|(_, v)| **v == max).map(|(k, _)| *k).collect();
    if top.len() == 1 {
        return top[0].clone();
    }
    top.sort();
    let mut best = String::new();
    let mut best_aligned = 0;
    for c in top {
        let aligned = if c == "*" { 1 } else { calc_cigar_aligned_length(c.as_bytes()) };
        if aligned > best_aligned {
            best = c.clone();
            best_aligned = aligned;
        }
    }
    best
}

/// ConsensusReads.formConsensusReadId: CNS_ goes before the last ':' item, the UMI replaces that item when present
pub fn form_consensus_read_id(read_id: &[u8], umi: Option<&str>) -> Vec<u8> {
    let last_delim = read_id.iter().rposition(|&c| c == b':');
    let mut out = Vec::with_capacity(read_id.len() + 8);
    match last_delim {
        Some(d) if d > 0 => {
            out.extend_from_slice(&read_id[..=d]);
            out.extend_from_slice(CONSENSUS_PREFIX.as_bytes());
            match umi {
                Some(u) => out.extend_from_slice(u.as_bytes()),
                None => out.extend_from_slice(&read_id[d + 1..]),
            }
        }
        _ => match umi {
            Some(u) => {
                out.extend_from_slice(read_id);
                out.push(b':');
                out.extend_from_slice(CONSENSUS_PREFIX.as_bytes());
                out.extend_from_slice(u.as_bytes());
            }
            None => {
                out.extend_from_slice(CONSENSUS_PREFIX.as_bytes());
                out.extend_from_slice(read_id);
            }
        },
    }
    out
}

fn to_cigar_string(cigar: &[u32]) -> CigarString {
    CigarString(
        cigar
            .iter()
            .map(|&c| {
                let n = c >> 4;
                match op(c) {
                    OP_M => Cigar::Match(n),
                    OP_I => Cigar::Ins(n),
                    OP_D => Cigar::Del(n),
                    OP_N => Cigar::RefSkip(n),
                    OP_S => Cigar::SoftClip(n),
                    OP_H => Cigar::HardClip(n),
                    OP_P => Cigar::Pad(n),
                    OP_EQ => Cigar::Equal(n),
                    _ => Cigar::Diff(n),
                }
            })
            .collect(),
    )
}

/// ConsensusReads.buildFromRead: a copy of one read under the consensus name
fn build_from_read(read: &bam::Record, read_id: &[u8], is_first_of_pair: bool) -> bam::Record {
    let mut rec = bam::Record::new();
    let cigar = to_cigar_string(read.raw_cigar());
    rec.set(read_id, Some(&cigar), &read.seq().as_bytes(), read.qual());
    rec.set_tid(read.tid());
    rec.set_pos(read.pos());
    rec.set_mapq(read.mapq());
    let mut flags = read.flags() & !FLAG_DUPLICATE;
    rec.set_mtid(-1);
    rec.set_mpos(-1);
    rec.set_insert_size(0);
    if flags & FLAG_PAIRED != 0 {
        flags = if is_first_of_pair { (flags | FLAG_FIRST) & !FLAG_SECOND } else { (flags | FLAG_SECOND) & !FLAG_FIRST };
        rec.set_mtid(read.mtid());
        rec.set_mpos(read.mpos());
        rec.set_insert_size(read.insert_size());
    }
    rec.set_flags(flags);
    set_aux_bytes(&mut rec, &aux_sorted_copy(aux_bytes(read)));
    fix_bin(&mut rec);
    rec
}

/// ConsensusReads.createConsensusRead(state, template)
fn create_consensus_record(state: &ConsensusState, template: &bam::Record, read_id: &[u8]) -> bam::Record {
    let mut rec = bam::Record::new();
    let cigar = to_cigar_string(if has(template, FLAG_UNMAPPED) { template.raw_cigar() } else { &state.cigar });
    let quals: Vec<u8> = state.base_qualities.iter().map(|q| *q as u8).collect();
    rec.set(read_id, Some(&cigar), &state.bases, &quals);
    rec.set_tid(template.tid());
    rec.set_pos(state.alignment_start as i64 - 1);
    rec.set_mapq(state.map_quality);
    let flags = template.flags() & !FLAG_DUPLICATE;
    rec.set_flags(flags);
    rec.set_mtid(-1);
    rec.set_mpos(-1);
    rec.set_insert_size(0);
    if flags & FLAG_PAIRED != 0 {
        if template.mtid() >= 0 {
            rec.set_mtid(template.mtid());
        }
        rec.set_mpos(template.mpos());
        rec.set_insert_size(template.insert_size());
    }
    let mut aux = aux_sorted_copy(aux_bytes(template));
    aux_block_set(&mut aux, b"NM", Some(&encode_int(state.num_mutations as i64)));
    set_aux_bytes(&mut rec, &aux);
    fix_bin(&mut rec);
    rec
}

struct ReadNameAttributes<'a> {
    lane_key: (&'a [u8], i32, &'a [u8], i32),
    tile: i32,
    x: i32,
    y: i32,
}

/// Integer.parseInt: an optional sign and decimal digits, within the int range
fn java_parse_int(s: &[u8]) -> Option<i32> {
    let s = std::str::from_utf8(s).ok()?;
    if s.is_empty() || !s.bytes().enumerate().all(|(i, c)| c.is_ascii_digit() || (i == 0 && (c == b'-' || c == b'+') && s.len() > 1)) {
        return None;
    }
    s.parse::<i32>().ok()
}

/// IlluminaBamUtils.getReadNameAttributes: instrument, run, flowcell, lane, tile, x, y from the read name
fn read_name_attributes(name: &[u8]) -> Option<ReadNameAttributes<'_>> {
    let mut parts: Vec<&[u8]> = name.split(|&c| c == b':').collect();
    // String.split drops trailing empty items
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    if parts.len() < 7 {
        return None;
    }
    Some(ReadNameAttributes {
        lane_key: (parts[0], java_parse_int(parts[1])?, parts[2], java_parse_int(parts[3])?),
        tile: java_parse_int(parts[4])?,
        x: java_parse_int(parts[5])?,
        y: java_parse_int(parts[6])?,
    })
}

fn tile_distance(a: (i32, i32), b: (i32, i32)) -> f64 {
    let xd = a.0.wrapping_sub(b.0);
    let yd = a.1.wrapping_sub(b.1);
    (xd.wrapping_mul(xd).wrapping_add(yd.wrapping_mul(yd)) as f64).sqrt()
}

const OPTICAL_DUPLICATE_TILE_DIFFERENCE: [i32; 5] = [0, 1, 999, 1000, 1001];
const OPTICAL_DUPLICATE_DISTANCE_THRESHOLD: f64 = 2500.0;

/// IlluminaRoutines.calculatePCRClusterCount
pub fn calculate_pcr_cluster_count(group: &DuplicateGroup) -> i32 {
    let mut attrs = Vec::with_capacity(group.total_read_count());
    for r in group.reads.iter().chain(group.non_consensus_reads.iter()) {
        match read_name_attributes(r.qname()) {
            Some(a) => attrs.push(a),
            None => return UNSET_COUNT,
        }
    }
    if attrs.len() == 2 {
        let tile_difference = attrs[0].tile.wrapping_sub(attrs[1].tile).wrapping_abs();
        if attrs[0].lane_key == attrs[1].lane_key && OPTICAL_DUPLICATE_TILE_DIFFERENCE.contains(&tile_difference) {
            return 1;
        }
        return 2;
    }
    // coordinates by tile; the total does not depend on the order of the tiles
    let mut by_tile: Vec<((&[u8], i32, &[u8], i32, i32), Vec<(i32, i32)>)> = Vec::new();
    for a in &attrs {
        let key = (a.lane_key.0, a.lane_key.1, a.lane_key.2, a.lane_key.3, a.tile);
        match by_tile.iter_mut().find(|t| t.0 == key) {
            Some(t) => t.1.push((a.x, a.y)),
            None => by_tile.push((key, vec![(a.x, a.y)])),
        }
    }
    let mut count = 0;
    for (_, coords) in &by_tile {
        if coords.len() == 1 {
            count += 1;
        } else if coords.len() == 2 {
            count += if tile_distance(coords[0], coords[1]) <= OPTICAL_DUPLICATE_DISTANCE_THRESHOLD { 1 } else { 2 };
        } else {
            count += cluster_count(coords);
        }
    }
    count
}

/// Cluster.clusterCount: greedy count of clusters, repeatedly taking the point with the fewest neighbours
fn cluster_count(data: &[(i32, i32)]) -> i32 {
    use std::collections::{BTreeMap, BTreeSet};
    let mut neighbours: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for i in 0..data.len() {
        neighbours.entry(i).or_default();
        for j in i + 1..data.len() {
            neighbours.entry(j).or_default();
            if tile_distance(data[i], data[j]) <= OPTICAL_DUPLICATE_DISTANCE_THRESHOLD {
                neighbours.get_mut(&i).unwrap().insert(j);
                neighbours.get_mut(&j).unwrap().insert(i);
            }
        }
    }
    let mut count = 0;
    let roots: Vec<usize> = neighbours.keys().copied().collect();
    for r in roots {
        if neighbours[&r].is_empty() {
            neighbours.remove(&r);
            count += 1;
        }
    }
    while !neighbours.is_empty() {
        let mut min_count = usize::MAX;
        let mut min_idx = 0;
        for (r, n) in &neighbours {
            if n.len() < min_count {
                min_count = n.len();
                min_idx = *r;
            }
        }
        count += 1;
        let min_neighbours = neighbours.remove(&min_idx).unwrap();
        let roots: Vec<usize> = neighbours.keys().copied().collect();
        for r in roots {
            if min_neighbours.contains(&r) {
                neighbours.remove(&r);
                continue;
            }
            let set = neighbours.get_mut(&r).unwrap();
            set.remove(&min_idx);
            for m in &min_neighbours {
                set.remove(m);
            }
            if set.is_empty() {
                neighbours.remove(&r);
                count += 1;
            }
        }
    }
    count
}
