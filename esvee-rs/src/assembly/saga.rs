//! SAGA matching (esvee.common.saga): matching assemblies to a resource of known indel variants by aligning
//! their sequence to the SAGA variant assemblies.

use std::sync::Arc;

use anyhow::{bail, Context, Result};

use crate::assembly::alignment::bwa::{BwaMemAlignParams, BwaMemAligner, BwaMemAlignment, IBwaMemAligner};
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::sam::*;
use crate::common::refgenome::RefGenome;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SagaBreakend {
    pub chromosome: String,
    pub position: i32,
    pub orientation: Orientation,
}

impl SagaBreakend {
    pub fn new(chromosome: &str, position: i32, orientation: Orientation) -> SagaBreakend {
        SagaBreakend { chromosome: chromosome.to_string(), position, orientation }
    }
    pub fn from_string(string: &str) -> Result<SagaBreakend> {
        let parts: Vec<&str> = string.split(':').collect();
        if parts.len() != 3 {
            bail!("java.lang.IllegalArgumentException: Invalid breakend string: {string}");
        }
        let position: i32 = parts[1].parse().with_context(|| format!("Invalid breakend string: {string}"))?;
        let orientation: i8 = parts[2].parse().with_context(|| format!("Invalid breakend string: {string}"))?;
        Ok(SagaBreakend::new(parts[0], position, Orientation::from_byte(orientation)))
    }
}

impl std::fmt::Display for SagaBreakend {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}:{}:{}", self.chromosome, self.position, self.orientation.as_byte())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SagaVariant {
    pub id: String,
    pub breakend1: SagaBreakend,
    pub breakend2: SagaBreakend,
    pub insert_sequence: String,
}

impl SagaVariant {
    /// Only indels are allowed: a forward then a reverse breakend on one chromosome. The errors carry the text
    /// of the exception that Java reports for the same input.
    pub fn new(id: &str, breakend1: SagaBreakend, breakend2: SagaBreakend, insert_sequence: &str) -> Result<SagaVariant> {
        if breakend1.chromosome != breakend2.chromosome {
            bail!("java.lang.IllegalArgumentException");
        }
        if !(breakend1.orientation.is_forward() && breakend2.orientation.is_reverse()) {
            bail!("java.lang.IllegalArgumentException");
        }
        if breakend1.position >= breakend2.position {
            bail!("java.lang.IllegalArgumentException");
        }
        let is_insert = breakend1.position + 1 == breakend2.position;
        if is_insert && insert_sequence.is_empty() {
            bail!("java.lang.IllegalArgumentException");
        }
        Ok(SagaVariant { id: id.to_string(), breakend1, breakend2, insert_sequence: insert_sequence.to_string() })
    }
    pub fn is_insert(&self) -> bool {
        !self.insert_sequence.is_empty()
    }
    pub fn no_space_string(&self) -> String {
        format!("{}:{}-{}", self.id, self.breakend1, self.breakend2).replace(' ', "_")
    }
}

impl std::fmt::Display for SagaVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} {}-{}", self.id, self.breakend1, self.breakend2)
    }
}

/// The index just after a junction in an assembly sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SagaJunctionInfo {
    pub offset: i32,
}

#[derive(Clone, Debug)]
pub struct SagaAssembly {
    pub fasta_label: String,
    pub variant: SagaVariant,
    /// Junctions occur just before each of these indices in the assembly sequence: one for a deletion at a
    /// single position, otherwise two.
    pub junction_offsets: Vec<i32>,
    pub sequence: String,
}

impl SagaAssembly {
    pub fn new(fasta_label: &str, variant: SagaVariant, junction_offsets: Vec<i32>, sequence: &str) -> Result<SagaAssembly> {
        if junction_offsets.len() != 1 && junction_offsets.len() != 2 {
            bail!("java.lang.IllegalArgumentException: Expected 1 or 2 junction offsets");
        }
        if variant.is_insert() && junction_offsets.len() != 2 {
            bail!("java.lang.IllegalArgumentException: Insert should have 2 junction offsets");
        }
        if !junction_offsets.iter().all(|offset| *offset >= 1 && *offset < sequence.len() as i32) {
            bail!("java.lang.IllegalArgumentException: Junction offsets out of bounds");
        }
        if junction_offsets.len() == 2 && junction_offsets[0] >= junction_offsets[1] {
            bail!("java.lang.IllegalArgumentException: Junction offsets not in ascending order");
        }
        Ok(SagaAssembly { fasta_label: fasta_label.to_string(), variant, junction_offsets, sequence: sequence.to_string() })
    }
    pub fn length(&self) -> i32 {
        self.sequence.len() as i32
    }
    pub fn junctions(&self) -> Vec<SagaJunctionInfo> {
        self.junction_offsets.iter().map(|o| SagaJunctionInfo { offset: *o }).collect()
    }
    /// Labels are `<id>|<chr:pos:orient>|<chr:pos:orient>|<offset1>|<offset2 or nothing>`.
    pub fn from_fasta_record(label: &str, sequence: &str) -> Result<SagaAssembly> {
        // Java String.split drops trailing empty parts
        let mut parts: Vec<&str> = label.split('|').collect();
        while parts.last().is_some_and(|p| p.is_empty()) {
            parts.pop();
        }
        if !(parts.len() == 4 || parts.len() == 5) {
            bail!("java.lang.IllegalArgumentException: Invalid fasta label: {label}");
        }
        let breakend1 = SagaBreakend::from_string(parts[1])?;
        let breakend2 = SagaBreakend::from_string(parts[2])?;
        let junction_offset1: i32 = parts[3].parse().with_context(|| format!("Invalid fasta label: {label}"))?;
        let junction_offset2: Option<i32> = if parts.len() >= 5 { Some(parts[4].parse().with_context(|| format!("Invalid fasta label: {label}"))?) } else { None };
        let insert_end = junction_offset2.unwrap_or(junction_offset1);
        if junction_offset1 < 0 || insert_end < junction_offset1 || insert_end as usize > sequence.len() {
            bail!("Invalid fasta label: {label}");
        }
        let insert_sequence = &sequence[junction_offset1 as usize..insert_end as usize];
        let variant = SagaVariant::new(parts[0], breakend1, breakend2, insert_sequence)?;
        let junction_offsets = match junction_offset2 {
            None => vec![junction_offset1],
            Some(o) => vec![junction_offset1, o],
        };
        SagaAssembly::new(label, variant, junction_offsets, sequence)
    }
}

/// The SAGA FASTA: one variant assembly per record, in sequence dictionary order.
pub struct SagaResource {
    fasta_path: String,
    assemblies: Vec<Arc<SagaAssembly>>,
}

impl SagaResource {
    pub fn new(fasta_path: &str) -> Result<SagaResource> {
        let fasta = RefGenome::open(fasta_path).context("Failed to load SAGA resource")?;
        let mut assemblies: Vec<Arc<SagaAssembly>> = Vec::new();
        for name in fasta.sequence_names() {
            let length = fasta.chromosome_length(name);
            let bases = fasta.bases(name, 1, length).context("Failed to load SAGA resource")?;
            let sequence = String::from_utf8_lossy(&bases).into_owned();
            assemblies.push(Arc::new(SagaAssembly::from_fasta_record(name, &sequence)?));
        }
        Self::from_assemblies(fasta_path, assemblies)
    }
    /// Also the test constructor: the assemblies are the contigs of the BWA index in order.
    pub fn from_assemblies(fasta_path: &str, assemblies: Vec<Arc<SagaAssembly>>) -> Result<SagaResource> {
        let mut ids: Vec<&str> = assemblies.iter().map(|a| a.variant.id.as_str()).collect();
        ids.sort_unstable();
        ids.dedup();
        if ids.len() != assemblies.len() {
            bail!("java.lang.RuntimeException: Duplicate variant IDs in SAGA resource");
        }
        crate::debug!("loaded {} SAGA variant assemblies from file({})", assemblies.len(), fasta_path);
        Ok(SagaResource { fasta_path: fasta_path.to_string(), assemblies })
    }
    pub fn bwa_index_image_path(&self) -> String {
        format!("{}.img", self.fasta_path)
    }
    pub fn assemblies(&self) -> &[Arc<SagaAssembly>] {
        &self.assemblies
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SagaSequenceMatcherConfig {
    pub align_length_ratio_min: f64,
    pub align_score_min: i32,
    pub align_score_ratio_min: f64,
    pub junction_overlap_min: i32,
    pub junction_overlap_min_lower: i32,
    pub junction_indel_distance: i32,
    pub junction_indel_length_max: i32,
}

impl SagaSequenceMatcherConfig {
    pub const DEFAULT: SagaSequenceMatcherConfig = SagaSequenceMatcherConfig {
        align_length_ratio_min: SAGA_ALIGN_LENGTH_MIN_RATIO,
        align_score_min: SAGA_ALIGN_SCORE_MIN_BASELINE,
        align_score_ratio_min: SAGA_ALIGN_SCORE_MIN_RATIO,
        junction_overlap_min: SAGA_ALIGN_JUNCTION_OVERLAP_MIN,
        junction_overlap_min_lower: SAGA_ALIGN_JUNCTION_OVERLAP_MIN_LOWER,
        junction_indel_distance: SAGA_ALIGN_JUNCTION_INDEL_DISTANCE,
        junction_indel_length_max: SAGA_ALIGN_JUNCTION_INDEL_MAX_LENGTH,
    };
}

/// An alignment of a query sequence to a SAGA variant assembly.
#[derive(Clone, Debug)]
pub struct SagaAlignment {
    pub raw_alignment: BwaMemAlignment,
    pub cigar: Vec<CigarElement>,
    pub query_length: i32,
    pub saga_assembly: Arc<SagaAssembly>,
}

impl SagaAlignment {
    pub fn new(raw_alignment: BwaMemAlignment, cigar: Vec<CigarElement>, query_length: i32, saga_assembly: Arc<SagaAssembly>) -> SagaAlignment {
        SagaAlignment { raw_alignment, cigar, query_length, saga_assembly }
    }
    /// An IllegalArgumentException in Java, which ends the run.
    pub fn validate(&self) {
        let valid = self.query_start() >= 0
            && self.query_end() <= self.query_length
            && self.saga_start() >= 0
            && self.saga_end() <= self.saga_length()
            && cigar_base_length(&self.cigar) == self.query_length
            && cigar_aligned_length(&self.cigar) <= self.saga_length();
        if !valid {
            panic!("invalid SAGA alignment: {self}");
        }
    }
    pub fn is_forward(&self) -> bool {
        self.raw_alignment.sam_flag & FLAG_READ_REVERSE_STRAND == 0
    }
    pub fn query_start(&self) -> i32 {
        if self.is_forward() { left_clip_length(&self.cigar) } else { right_clip_length(&self.cigar) }
    }
    pub fn query_end(&self) -> i32 {
        self.query_length - if self.is_forward() { right_clip_length(&self.cigar) } else { left_clip_length(&self.cigar) }
    }
    pub fn query_align_length(&self) -> i32 {
        self.query_end() - self.query_start()
    }
    pub fn saga_start(&self) -> i32 {
        self.raw_alignment.ref_start
    }
    pub fn saga_end(&self) -> i32 {
        self.raw_alignment.ref_end
    }
    pub fn saga_align_length(&self) -> i32 {
        self.saga_end() - self.saga_start()
    }
    pub fn saga_length(&self) -> i32 {
        self.saga_assembly.length()
    }
    pub fn align_score(&self) -> i32 {
        self.raw_alignment.aligner_score
    }
    /// How many more bases could have been aligned at the start of the sequences.
    pub fn start_unaligned(&self) -> i32 {
        self.query_start().min(self.saga_start())
    }
    /// How many more bases could have been aligned at the end of the sequences.
    pub fn end_unaligned(&self) -> i32 {
        (self.query_length - self.query_end()).min(self.saga_length() - self.saga_end())
    }
    pub fn saga_variant(&self) -> &SagaVariant {
        &self.saga_assembly.variant
    }
    /// The indices of the SAGA junctions mapped into the query sequence range, with extrapolation if required.
    pub fn query_junction_offsets(&self) -> Vec<i32> {
        self.saga_assembly.junction_offsets.iter().map(|o| self.saga_index_to_query_index(*o)).collect()
    }
    pub fn saga_index_to_query_index(&self, saga_index: i32) -> i32 {
        if saga_index <= self.saga_start() {
            return self.query_start() - (self.saga_start() - saga_index);
        }
        if saga_index >= self.saga_end() {
            return self.query_end() + (saga_index - self.saga_end());
        }
        // on a D element, use the next index
        let result = get_read_index_from_position_gap(self.saga_start(), &self.cigar, saga_index, 1, false);
        if !(result >= 0 && result <= self.query_length) {
            panic!("SAGA index {saga_index} could not be mapped into the query for {self}");
        }
        result
    }
    pub fn cigar_string(&self) -> String {
        cigar_elements_to_str(&self.cigar)
    }
}

impl std::fmt::Display for SagaAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "SagaAlignment(variant=\"{}\", queryAlign=[{}, {}), sagaAlign=[{}, {}), strand={}, cigar={}, alignScore={})",
            self.saga_variant(),
            self.query_start(),
            self.query_end(),
            self.saga_start(),
            self.saga_end(),
            if self.is_forward() { "+" } else { "-" },
            self.cigar_string(),
            self.align_score()
        )
    }
}

#[derive(Clone, Debug)]
pub struct SagaSequenceMatch {
    pub alignment: SagaAlignment,
}

impl SagaSequenceMatch {
    pub fn variant(&self) -> &SagaVariant {
        self.alignment.saga_variant()
    }
    /// SagaVariant.toString
    pub fn variant_string(&self) -> String {
        self.variant().to_string()
    }
    /// SagaVariant.noSpaceString
    pub fn variant_no_space_string(&self) -> String {
        self.variant().no_space_string()
    }
    pub fn cigar_string(&self) -> String {
        self.alignment.cigar_string()
    }
}

struct SagaJunctionMatchInfo {
    align_overlap_left: i32,
    align_overlap_right: i32,
    indel_nearby: bool,
}

struct SagaSequenceMatchCandidate {
    alignment: SagaAlignment,
    accepted: bool,
}

struct CigarElemWithPos {
    element: CigarElement,
    /// query positions are in the original query sequence (significant for reverse strand matches)
    query_start: i32,
    query_end: i32,
    ref_start: i32,
    ref_end: i32,
}

struct MatchArguments<'a> {
    query: &'a [u8],
    junctions: &'a [SagaJunctionInfo],
    lower_junction_overlap: bool,
    allow_reverse_strand: bool,
}

/// Matches assemblies to SAGA variants by base sequence similarity.
pub struct SagaSequenceMatcher {
    config: SagaSequenceMatcherConfig,
    aligner: Arc<dyn IBwaMemAligner>,
    assemblies_by_contig_id: Arc<Vec<Arc<SagaAssembly>>>,
}

impl SagaSequenceMatcher {
    pub fn new(config: SagaSequenceMatcherConfig, aligner: Arc<dyn IBwaMemAligner>, assemblies_by_contig_id: Arc<Vec<Arc<SagaAssembly>>>) -> SagaSequenceMatcher {
        SagaSequenceMatcher { config, aligner, assemblies_by_contig_id }
    }

    pub fn match_by_sequence(&mut self, sequence: &[u8], junctions: &[SagaJunctionInfo], lower_junction_overlap: bool, allow_reverse_strand: bool) -> Option<SagaSequenceMatch> {
        let args = MatchArguments { query: sequence, junctions, lower_junction_overlap, allow_reverse_strand };
        let alignments = self.aligner.align_sequence(args.query);
        self.match_from_alignments(&args, alignments)
    }

    fn match_from_alignments(&self, args: &MatchArguments, raw_alignments: Vec<BwaMemAlignment>) -> Option<SagaSequenceMatch> {
        let query_length = args.query.len() as i32;
        let mut candidates: Vec<SagaSequenceMatchCandidate> = raw_alignments
            .into_iter()
            .filter_map(|alignment| self.wrap_alignment(alignment, query_length))
            .map(|alignment| self.evaluate_alignment(alignment, args))
            .collect();
        // accepted first, then the higher alignment score, then the variant whose length is closer to the query length
        candidates.sort_by(|o1, o2| {
            if o1.accepted != o2.accepted {
                return if o1.accepted { std::cmp::Ordering::Less } else { std::cmp::Ordering::Greater };
            }
            let res = (-o1.alignment.align_score()).cmp(&(-o2.alignment.align_score()));
            if res != std::cmp::Ordering::Equal {
                return res;
            }
            let distance1 = (o1.alignment.saga_assembly.length() - query_length).abs();
            let distance2 = (o2.alignment.saga_assembly.length() - query_length).abs();
            distance1.cmp(&distance2)
        });
        candidates.into_iter().find(|c| c.accepted).map(|c| SagaSequenceMatch { alignment: c.alignment })
    }

    fn wrap_alignment(&self, alignment: BwaMemAlignment, query_length: i32) -> Option<SagaAlignment> {
        if alignment.ref_id < 0 || alignment.sam_flag & FLAG_READ_UNMAPPED != 0 {
            // not a real alignment: the query is not aligned
            return None;
        }
        let saga_assembly = self.assemblies_by_contig_id.get(alignment.ref_id as usize).unwrap_or_else(|| panic!("No SAGA variant with contig ID: {}", alignment.ref_id)).clone();
        let cigar = cigar_elements_from_str(&alignment.cigar);
        let saga_alignment = SagaAlignment::new(alignment, cigar, query_length, saga_assembly);
        saga_alignment.validate();
        Some(saga_alignment)
    }

    fn evaluate_alignment(&self, alignment: SagaAlignment, args: &MatchArguments) -> SagaSequenceMatchCandidate {
        let mut accepted = true;
        // a junction assembly may match SAGA in reverse if it is part of an INV; phased assemblies are always forward strand
        if !args.allow_reverse_strand && !alignment.is_forward() {
            accepted = false;
        }
        // only match if the majority of the sequences align, where possible
        if !self.check_align_length(&alignment) {
            accepted = false;
        }
        // a low alignment score means low sequence similarity
        if !self.check_align_score(&alignment) {
            accepted = false;
        }
        let cigar_elements = get_cigar_element_positions(&alignment.cigar, alignment.saga_start(), alignment.is_forward(), alignment.query_length);
        let query_junction_match_infos: Vec<SagaJunctionMatchInfo> = args.junctions.iter().map(|j| self.evaluate_junction_alignment(&alignment, &cigar_elements, *j, true)).collect();
        let saga_junction_match_infos: Vec<SagaJunctionMatchInfo> =
            alignment.saga_assembly.junctions().iter().map(|j| self.evaluate_junction_alignment(&alignment, &cigar_elements, *j, false)).collect();
        if !self.junction_matches_pass(args, &query_junction_match_infos) {
            accepted = false;
        }
        if !self.junction_matches_pass(args, &saga_junction_match_infos) {
            accepted = false;
        }
        SagaSequenceMatchCandidate { alignment, accepted }
    }

    fn check_align_length(&self, alignment: &SagaAlignment) -> bool {
        // looks at how many bases could have been aligned on either side of the aligned subsequence, since a
        // single-sided assembly may not have enough bases to align through most of the other side
        let unaligned = alignment.start_unaligned() + alignment.end_unaligned();
        let max_unaligned = (alignment.query_length.min(alignment.saga_length()) as f64 * (1.0 - self.config.align_length_ratio_min)).ceil() as i32;
        unaligned <= max_unaligned
    }

    fn check_align_score(&self, alignment: &SagaAlignment) -> bool {
        let min_align_score = (alignment.query_align_length().min(alignment.saga_align_length()) as f64 * self.config.align_score_ratio_min).ceil() as i32;
        alignment.align_score() >= min_align_score
    }

    fn evaluate_junction_alignment(&self, alignment: &SagaAlignment, cigar_elements: &[CigarElemWithPos], junction_info: SagaJunctionInfo, is_query: bool) -> SagaJunctionMatchInfo {
        let align_start = if is_query { alignment.query_start() } else { alignment.saga_start() };
        let align_end = if is_query { alignment.query_end() } else { alignment.saga_end() };
        let indel_nearby = cigar_elements.iter().filter(|e| e.element.op.is_indel() && e.element.len > self.config.junction_indel_length_max).any(|e| {
            let (indel_start, indel_end) = if is_query { (e.query_start, e.query_end) } else { (e.ref_start, e.ref_end) };
            self.is_junction_near_indel(junction_info.offset, indel_start, indel_end)
        });
        SagaJunctionMatchInfo { align_overlap_left: junction_info.offset - align_start, align_overlap_right: align_end - junction_info.offset, indel_nearby }
    }

    fn is_junction_near_indel(&self, junction_offset: i32, indel_start: i32, indel_end: i32) -> bool {
        let inside_indel = junction_offset >= indel_start && junction_offset < indel_end;
        let near_indel_start = (junction_offset - indel_start).abs() < self.config.junction_indel_distance;
        let near_indel_end = (junction_offset - indel_end).abs() < self.config.junction_indel_distance;
        inside_indel || near_indel_start || near_indel_end
    }

    /// filterJunctionMatches reduced to whether it adds a filter: at least one junction must be covered by the
    /// alignment with the minimum overlap either side, and no junction may have a large indel nearby.
    fn junction_matches_pass(&self, args: &MatchArguments, junction_match_infos: &[SagaJunctionMatchInfo]) -> bool {
        let min_overlap = if args.lower_junction_overlap { self.config.junction_overlap_min_lower } else { self.config.junction_overlap_min };
        let any_aligned = junction_match_infos.iter().any(|j| j.align_overlap_left >= min_overlap && j.align_overlap_right >= min_overlap);
        let any_indel_nearby = junction_match_infos.iter().any(|j| j.indel_nearby);
        any_aligned && !any_indel_nearby
    }
}

fn get_cigar_element_positions(cigar: &[CigarElement], ref_start: i32, is_forward: bool, query_length: i32) -> Vec<CigarElemWithPos> {
    let mut result: Vec<CigarElemWithPos> = Vec::with_capacity(cigar.len());
    let mut query_index = 0;
    let mut ref_index = ref_start;
    for element in cigar {
        let mut next_query_index = query_index;
        if element.op.consumes_read_bases() || element.op.is_clipping() {
            next_query_index += element.len;
        }
        let mut next_ref_index = ref_index;
        if element.op.consumes_reference_bases() {
            next_ref_index += element.len;
        }
        let query_start = if is_forward { query_index } else { query_length - next_query_index };
        let query_end = if is_forward { next_query_index } else { query_length - query_index };
        result.push(CigarElemWithPos { element: *element, query_start, query_end, ref_start: ref_index, ref_end: next_ref_index });
        query_index = next_query_index;
        ref_index = next_ref_index;
    }
    result
}

/// Creates the sequence matchers, which share the SAGA assemblies and the BWA index of the SAGA FASTA.
pub struct SagaMatcherFactory {
    sequence_matcher_config: SagaSequenceMatcherConfig,
    assemblies_by_contig_id: Arc<Vec<Arc<SagaAssembly>>>,
    aligner: Arc<dyn IBwaMemAligner>,
}

impl SagaMatcherFactory {
    /// Loads the SAGA FASTA and opens its BWA index image (`<fasta>.img`): all alignments reported, default
    /// parameters apart from the minimum alignment score. Java opens one aligner per matcher; one is shared here.
    pub fn new(saga_fasta_file: &str) -> Result<SagaMatcherFactory> {
        let resource = SagaResource::new(saga_fasta_file)?;
        let config = SagaSequenceMatcherConfig::DEFAULT;
        let params = BwaMemAlignParams { min_align_score: config.align_score_min, ..BwaMemAlignParams::DEFAULT };
        let aligner = BwaMemAligner::new_all_alignments(&resource.bwa_index_image_path(), &params, true)?;
        Ok(Self::with_aligner(&resource, config, Arc::new(aligner)))
    }
    pub fn with_aligner(resource: &SagaResource, config: SagaSequenceMatcherConfig, aligner: Arc<dyn IBwaMemAligner>) -> SagaMatcherFactory {
        SagaMatcherFactory { sequence_matcher_config: config, assemblies_by_contig_id: Arc::new(resource.assemblies().to_vec()), aligner }
    }
    pub fn create_sequence_matcher(&self) -> SagaSequenceMatcher {
        SagaSequenceMatcher::new(self.sequence_matcher_config, self.aligner.clone(), self.assemblies_by_contig_id.clone())
    }
}
