//! Ports of the Java unit tests (esvee/src/test/.../assembly), with the Java test utilities they use
//! (TestUtils, AssemblyTestUtils, hmf-common SamRecordTestUtils / MockRefGenome / ReadIdGenerator).
//! All sequences are the synthetic ones of the Java tests.
#![allow(clippy::too_many_arguments)]

use std::collections::HashMap;
use std::sync::atomic::{AtomicI32, Ordering};
use std::sync::MutexGuard;

use crate::assembly::alignment::align_data::{AlignData, AlignP};
use crate::assembly::alignment::assembly_alignment::{AsmAlignP, AssemblyAlignment};
use crate::assembly::assembly_utils::reverse_complement_str;
use crate::assembly::config::AssemblyConfig;
use crate::assembly::globals::{set_sequencing_type, SequencingType, TEST_GLOBALS_LOCK};
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::P;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::ref_genome::RefGenomeInterface;
use crate::assembly::sam::*;
use crate::assembly::types::assembly_link::{AssemblyLink, LinkP};
use crate::assembly::types::junction::Junction;
use crate::assembly::types::junction_assembly::{AsmP, JunctionAssembly};
use crate::assembly::types::phase_set::PhaseSet;
use crate::assembly::types::{LinkType, SupportType};

mod alignment_test;
mod assembly_links_test;
mod assembly_stats_test;
mod base_qual_bands_test;
mod breakend_builder_saga_test;
mod homology_test;
mod indels_test;
mod java_findings_test;
mod junction_assembly_test;
mod junction_groups_test;
mod line_site_test;
mod local_links_test;
mod phase_group_test;
mod read_adjustments_test;
mod read_utils_test;
mod ref_base_extension_test;
mod ref_base_sequence_test;
mod sbx_assembly_test;
mod sequence_test;
mod specific_alignments_test;
mod unmapped_extension_test;

pub const CHR_1: &str = "1";
pub const CHR_2: &str = "2";
pub const CHR_3: &str = "3";
pub const TEST_READ_ID: &str = "READ_01";
pub const TEST_READ_ID_2: &str = "READ_02";
pub const TEST_CIGAR_100: &str = "100M";
pub const DEFAULT_MAP_QUAL: i32 = 60;
pub const DEFAULT_BASE_QUAL: u8 = 37;
pub const DEFAULT_NM: i32 = 0;
/// the mapping quality htsjdk's SAMRecordSetBuilder gives the records of the Java test utilities
const RECORD_SET_BUILDER_MAP_QUAL: i32 = 255;

pub const FORWARD: Orientation = Orientation::Forward;
pub const REVERSE: Orientation = Orientation::Reverse;

pub const REF_BASES_200: &str = concat!(
    "AAACCCGGGTTTACGTAACCGGTTACGTAAAAACCCCCGGGGGTTTTTACGTAACCGGTTACGTAAACCCGGGTTTAAACGTTTTTGGGGCCCCAAAAAC",
    "GTGGCCTTAAACGTCCCAAAATTTTGGGGACGTGGGGGCCCCCAAAATTTTACGTCCGGTTAAACGTTTTCCCGGGAAAACGTAACCGGTTGGCCAAGCT"
);
pub const REF_BASES_400: &str = concat!(
    "AGTGTCGCGAATGCTGATTCGGAACCTTGATCGTGATGTCGATGGCTAGATCAAATCGTCTAGTGGCTAATGCTCGATCGATATGGCTTGATCGTGAGTC",
    "GGGCTGATGATGTTGATCGTAGTGCTAGTGATCGTAGTTGCCAAATGCTGTTGCCCTGTAGCTGATTGGCTAGCTGTAGCTGATCGATGCAATCCGGTAG",
    "TGTAGCTGATCGCAGGTCGAACCTGGTGATCGATGTCGATCGACTGATGTAGTAGCTGATCGGATGCATGCGTAGCGATGCTAGCTGATCGATTGGCTAA",
    "GTCGCTTCCGGTATTTGCGTTCCGGGTTTTTTCCGAGCCTACCCCAGTTGGTTAAAAGGATATTATATATATGGCGGCTATATATGCGGTGTGTGTAACC"
);
pub const REF_BASES_600: &str = concat!(
    "ATCATCGAATGGAATGGAATGGAACAGTCAATGAACTCGAATGGAATCATCATTGAATGGAATCGAATGGAATCATCGAGTGGAATCGAATGGAATTATG",
    "ATCAAATGGAATCGAATGTAATCATCATCAAATGGAATCAAAAATAACCATAATATGGTCTTTAAAGAAGCAATCTAGCTAAAATGAAATCATTAATCCA",
    "ATAGGACTAACATCCTTATGAGAAGACAAGATTAGGACACAGGCATGAACTGGGGGAAGACCATGTAAAGACACTCAAAGAACTGACCGTACCACCCCCT",
    "TCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGTAAAGAAGCATTGAGCTATTTACATAGAAATCTCATTTAACTGTG",
    "ATATAAATTAACCTTTCTTATCCTGCTTCTAAACAAAGGTAAGGGCCACCCAGTCAATGCTTTGTATTCTTCCAATATTCTTTCCTAGAACTTCTTCAAA",
    "GGCTCTCATGAAGCACTGGTGAAACTGGAAATCACTGAATTTTACTACCATTTTCTTATCCTGCTTCTAAACAAAGGTAAGTTTCTTATCCTGCTTCTAA"
);

/// Holds the lock on the process-wide settings for a test and puts them into the Java test defaults
/// (TestUtils.TEST_CONFIG, Illumina sequencing); restores Illumina when dropped (the Java @After methods).
pub struct TestEnv {
    pub config: AssemblyConfig,
    _guard: MutexGuard<'static, ()>,
}

impl Drop for TestEnv {
    fn drop(&mut self) {
        set_sequencing_type(SequencingType::Illumina);
    }
}

pub fn test_env() -> TestEnv {
    let guard = TEST_GLOBALS_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    set_sequencing_type(SequencingType::Illumina);
    let config = AssemblyConfig::new_test();
    TestEnv { config, _guard: guard }
}

pub fn set_illumina_sequencing() {
    set_sequencing_type(SequencingType::Illumina);
}
pub fn set_sbx_sequencing() {
    set_sequencing_type(SequencingType::Sbx);
}
pub fn set_ultima_sequencing() {
    set_sequencing_type(SequencingType::Ultima);
}

/// MockRefGenome.generateRandomBases: cycles through the nucleotides
pub fn generate_random_bases(length: usize) -> String {
    const DNA_BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    (0..length).map(|i| DNA_BASES[i % 4]).collect()
}

pub fn ref_bases_random_100() -> String {
    generate_random_bases(100)
}

/// MockRefGenome.getNextBase
pub fn get_next_base(base: &str) -> String {
    let c = match base.as_bytes()[0] {
        b'A' => 'C',
        b'C' => 'G',
        b'G' => 'T',
        b'T' => 'A',
        other => other as char,
    };
    c.to_string()
}

pub fn get_next_base_byte(base: u8) -> u8 {
    get_next_base(&(base as char).to_string()).as_bytes()[0]
}

/// hmf-common MockRefGenome
pub struct MockRefGenome {
    one_based_indexing: bool,
    pub ref_genome_map: HashMap<String, String>,
}

impl MockRefGenome {
    pub fn new() -> Self {
        MockRefGenome { one_based_indexing: false, ref_genome_map: HashMap::new() }
    }
    pub fn new_indexing(one_based_indexing: bool) -> Self {
        MockRefGenome { one_based_indexing, ref_genome_map: HashMap::new() }
    }
    pub fn put(&mut self, chromosome: &str, bases: &str) {
        self.ref_genome_map.insert(chromosome.to_string(), bases.to_string());
    }
    pub fn base_string(&self, chromosome: &str, pos_start: i32, pos_end: i32) -> String {
        self.get_base_string(chromosome, pos_start, pos_end).unwrap_or_default()
    }
}

impl RefGenomeInterface for MockRefGenome {
    fn get_bases(&self, chromosome: &str, pos_start: i32, pos_end: i32) -> Option<Vec<u8>> {
        let (mut pos_start, mut pos_end) = (pos_start, pos_end);
        if self.one_based_indexing {
            pos_start -= 1;
            pos_end -= 1;
        }
        match self.ref_genome_map.get(chromosome) {
            Some(chr_bases) if pos_start >= 0 && chr_bases.len() as i32 > pos_end && pos_end + 1 >= pos_start => {
                Some(chr_bases.as_bytes()[pos_start as usize..(pos_end + 1) as usize].to_vec())
            }
            _ => Some(Vec::new()),
        }
    }
    fn get_chromosome_length(&self, chromosome: &str) -> i32 {
        // the Java mock holds the version 37 chromosome lengths
        const LENGTHS_37: &str = include_str!("../resources/lengths.37.tsv");
        let bare = chromosome.strip_prefix("chr").unwrap_or(chromosome);
        LENGTHS_37.lines().filter_map(|l| l.split_once('\t')).find(|(c, _)| c.trim_start_matches("chr") == bare).and_then(|(_, len)| len.trim().parse().ok()).unwrap_or(0)
    }
    fn chromosomes(&self) -> Vec<(String, i32)> {
        Vec::new()
    }
}

/// hmf-common ReadIdGenerator
pub struct ReadIdGenerator {
    next_read_id: AtomicI32,
}

impl ReadIdGenerator {
    pub const fn new() -> Self {
        ReadIdGenerator { next_read_id: AtomicI32::new(1) }
    }
    pub fn next_id(&self) -> String {
        format!("READ_{:03}", self.next_read_id.fetch_add(1, Ordering::Relaxed))
    }
    pub fn current_id(&self) -> String {
        format!("READ_{:03}", self.next_read_id.load(Ordering::Relaxed))
    }
    pub fn reset(&self) {
        self.next_read_id.store(1, Ordering::Relaxed);
    }
}

pub static READ_ID_GENERATOR: ReadIdGenerator = ReadIdGenerator::new();

pub fn build_default_base_quals(length: usize) -> Vec<u8> {
    vec![DEFAULT_BASE_QUAL; length]
}
pub fn build_base_quals(length: usize, qual_value: u8) -> Vec<u8> {
    vec![qual_value; length]
}

/// hmf-common SamRecordTestUtils.createSamRecord
pub fn create_sam_record_full(
    read_id: &str, chr_str: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, is_reversed: bool, is_supplementary: bool,
    supp_alignment: Option<&SupplementaryReadData>,
) -> SamRecord {
    let mut record = SamRecord::empty();
    record.read_name = read_id.to_string();
    record.alignment_start = read_start;
    record.set_flag(FLAG_READ_REVERSE_STRAND, is_reversed);
    record.cigar = if cigar == NO_CIGAR { Vec::new() } else { cigar_elements_from_str(cigar) };
    record.mapping_quality = RECORD_SET_BUILDER_MAP_QUAL;
    record.read_bases = read_bases.as_bytes().to_vec();
    record.base_qualities = build_default_base_quals(read_bases.len());
    record.reference_name = intern(chr_str);
    record.set_flag(FLAG_READ_PAIRED, true);
    if chr_str == NO_CHROMOSOME_NAME {
        record.set_flag(FLAG_READ_UNMAPPED, true);
    }
    if !mate_chr.is_empty() && mate_chr != NO_CHROMOSOME_NAME {
        record.mate_reference_name = intern(mate_chr);
        record.mate_alignment_start = mate_start;
        record.set_flag(FLAG_PROPER_PAIR, true);
    } else {
        record.set_flag(FLAG_MATE_UNMAPPED, true);
        record.mate_reference_name = NO_CHROMOSOME_NAME;
        record.mate_alignment_start = read_start;
        record.set_flag(FLAG_PROPER_PAIR, false);
    }
    record.set_flag(FLAG_FIRST_OF_PAIR, true);
    record.set_flag(FLAG_SUPPLEMENTARY, is_supplementary);
    if let Some(supp) = supp_alignment {
        record.set_string_attribute(SUPPLEMENTARY_ATTRIBUTE, &supp.as_sam_tag());
    }
    record.inferred_insert_size = if chr_str == mate_chr { (read_start - mate_start).abs() } else { 0 };
    record
}

/// SamRecordTestUtils.createSamRecord with mate orientation and cigar
pub fn create_sam_record_mate(
    read_id: &str, chr_str: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, is_reversed: bool, is_supplementary: bool,
    supp_alignment: Option<&SupplementaryReadData>, mate_reversed: bool, mate_cigar: Option<&str>,
) -> SamRecord {
    let mut record = create_sam_record_full(read_id, chr_str, read_start, read_bases, cigar, mate_chr, mate_start, is_reversed, is_supplementary, supp_alignment);
    if mate_reversed {
        record.set_flag(FLAG_MATE_REVERSE_STRAND, true);
    }
    if let Some(mc) = mate_cigar {
        record.set_string_attribute(MATE_CIGAR_ATTRIBUTE, mc);
    }
    record
}

/// SamRecordTestUtils.createSamRecordUnpaired
pub fn create_sam_record_unpaired(read_id: &str, chr_str: &str, read_start: i32, read_bases: &str, cigar: &str, is_reversed: bool, is_supplementary: bool, supp_alignment: Option<&SupplementaryReadData>) -> SamRecord {
    let mut record = create_sam_record_full(read_id, chr_str, read_start, read_bases, cigar, NO_CHROMOSOME_NAME, NO_POSITION, is_reversed, is_supplementary, supp_alignment);
    record.set_flag(FLAG_READ_PAIRED, false);
    record
}

/// SamRecordTestUtils.flipFirstInPair
pub fn flip_first_in_pair(record: &mut SamRecord) {
    if record.read_paired() {
        let first = record.first_of_pair();
        record.set_flag(FLAG_FIRST_OF_PAIR, !first);
        record.set_flag(FLAG_SECOND_OF_PAIR, first);
    }
}

/// TestUtils.createRead(readId, readStart, readBases, cigar)
pub fn create_read(read_id: &str, read_start: i32, read_bases: &str, cigar: &str) -> ReadP {
    Read::new(create_sam_record_full(read_id, CHR_1, read_start, read_bases, cigar, CHR_1, read_start + 1000, false, false, None))
}

/// TestUtils.createRead(readId, chromosome, readStart, readBases, cigar, mateChr, mateStart, mateReversed)
pub fn create_read_mate(read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, mate_reversed: bool) -> ReadP {
    let mut record = create_sam_record_full(read_id, chromosome, read_start, read_bases, cigar, mate_chr, mate_start, false, false, None);
    if mate_reversed {
        record.set_flag(FLAG_MATE_REVERSE_STRAND, true);
    }
    record.set_string_attribute(MATE_CIGAR_ATTRIBUTE, &format!("{}M", read_bases.len()));
    Read::new(record)
}

pub fn set_mate_cigar(read: &ReadP, mate_cigar: &str) {
    read.m().bam_record_mut().set_string_attribute(MATE_CIGAR_ATTRIBUTE, mate_cigar);
}

/// TestUtils.createConcordantRead
pub fn create_concordant_read(read_id: &str, read_start: i32, read_bases: &str, cigar: &str, mate_start: i32) -> ReadP {
    let mut record = create_sam_record_full(read_id, CHR_1, read_start, read_bases, cigar, CHR_1, mate_start, false, false, None);
    if read_start > mate_start {
        record.set_flag(FLAG_READ_REVERSE_STRAND, true);
    } else {
        record.set_flag(FLAG_MATE_REVERSE_STRAND, true);
    }
    record.set_string_attribute(MATE_CIGAR_ATTRIBUTE, cigar);
    Read::new(record)
}

/// TestUtils.cloneRead
pub fn clone_read(read: &ReadP, new_read_id: &str) -> ReadP {
    let mut record = read.bam_record().clone();
    record.read_name = new_read_id.to_string();
    Read::new(record)
}

/// TestUtils.makeCigarString
pub fn make_cigar_string(read_bases: &str, sc_left: usize, sc_right: usize) -> String {
    let mut s = String::new();
    if sc_left > 0 {
        s.push_str(&format!("{sc_left}S"));
    }
    s.push_str(&format!("{}M", read_bases.len() - sc_left - sc_right));
    if sc_right > 0 {
        s.push_str(&format!("{sc_right}S"));
    }
    s
}

pub fn set_second_in_pair(record: &mut SamRecord) {
    record.set_flag(FLAG_SECOND_OF_PAIR, true);
    record.set_flag(FLAG_FIRST_OF_PAIR, false);
}

pub fn read_id_str(read_id: i32) -> String {
    format!("READ_{read_id:02}")
}

/// TestUtils.buildFlags
pub fn build_flags(first_in_pair: bool, reversed: bool, supplementary: bool) -> i32 {
    let mut flags = FLAG_READ_PAIRED | FLAG_PROPER_PAIR;
    if reversed {
        flags |= FLAG_READ_REVERSE_STRAND;
    }
    flags |= if first_in_pair { FLAG_FIRST_OF_PAIR } else { FLAG_SECOND_OF_PAIR };
    if supplementary {
        flags |= FLAG_SUPPLEMENTARY;
    }
    flags
}

/// TestUtils.createSamRecord(readId, chromosome, readStart, readBases, cigar, flags, mapQual, baseQual)
pub fn create_sam_record_flags_quals(read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str, flags: i32, _map_qual: i32, base_qual: u8) -> SamRecord {
    let mut record = SamRecord::empty();
    record.read_name = read_id.to_string();
    record.alignment_start = read_start;
    record.cigar = if cigar == NO_CIGAR { Vec::new() } else { cigar_elements_from_str(cigar) };
    // the Java utility passes its map qual as the builder's default base quality, so the mapping quality stays at the builder's value
    record.mapping_quality = RECORD_SET_BUILDER_MAP_QUAL;
    record.read_bases = read_bases.as_bytes().to_vec();
    record.base_qualities = vec![base_qual; read_bases.len()];
    record.reference_name = intern(chromosome);
    record.flags = flags;
    record.mate_reference_name = intern(chromosome);
    record.mate_alignment_start = read_start + 300;
    record.set_flag(FLAG_MATE_REVERSE_STRAND, true);
    record.inferred_insert_size = 400;
    record
}

/// TestUtils.createSamRecord(readId, chromosome, readStart, readBases, cigar)
pub fn create_sam_record(read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str) -> SamRecord {
    create_sam_record_flags_quals(read_id, chromosome, read_start, read_bases, cigar, build_flags(true, false, false), DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL)
}

/// TestUtils.createSamRecord(readId, chromosome, readStart, readBases, cigar, flags)
pub fn create_sam_record_flags(read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str, flags: i32) -> SamRecord {
    create_sam_record_flags_quals(read_id, chromosome, read_start, read_bases, cigar, flags, DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL)
}

/// TestUtils.createSamRecord(readId, chromosome, readStart, mateChr, mateStart, firstInPair, isSupp, suppData)
pub fn create_sam_record_supp(read_id: &str, chromosome: &str, read_start: i32, mate_chr: &str, mate_start: i32, first_in_pair: bool, is_supp: bool, supp_data: Option<&str>) -> SamRecord {
    let read_bases = &REF_BASES_200[0..100];
    let mut record = create_sam_record_flags_quals(read_id, chromosome, read_start, read_bases, "100M", build_flags(first_in_pair, false, is_supp), DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL);
    record.mate_reference_name = intern(mate_chr);
    record.mate_alignment_start = mate_start;
    if let Some(s) = supp_data.filter(|s| !s.is_empty()) {
        record.set_string_attribute(SUPPLEMENTARY_ATTRIBUTE, s);
    }
    record
}

pub fn get_support_type_count(assembly: AsmP, type_: SupportType) -> usize {
    assembly.support().iter().filter(|x| x.type_() == type_).count()
}

// AssemblyTestUtils

pub fn new_junction(chromosome: &str, position: i32, orientation: Orientation) -> P<Junction> {
    P::new(Junction::new_basic(chromosome, position, orientation))
}

pub fn create_assembly(chromosome: &str, junction_position: i32, junction_orientation: Orientation, assembly_bases: &str, junction_index: i32) -> AsmP {
    let junction = new_junction(chromosome, junction_position, junction_orientation);
    let base_quals = build_default_base_quals(assembly_bases.len());
    let assembly = JunctionAssembly::new_test(junction, assembly_bases.as_bytes(), &base_quals, junction_index);
    assembly.m().build_repeat_info();
    assembly
}

pub fn create_assembly_alignment(
    ref_genome: &MockRefGenome, chr_start: &str, pos_start: i32, orient_start: Orientation, chr_end: &str, pos_end: i32, orient_end: Orientation, inserted_bases: &str,
    overlap_bases: &str,
) -> AsmAlignP {
    create_assembly_alignment_ext(ref_genome, chr_start, pos_start, orient_start, chr_end, pos_end, orient_end, inserted_bases, overlap_bases, 50, 100)
}

pub fn create_assembly_alignment_ext(
    ref_genome: &MockRefGenome, chr_start: &str, pos_start: i32, orient_start: Orientation, chr_end: &str, pos_end: i32, orient_end: Orientation, inserted_bases: &str,
    overlap_bases: &str, ext_base_length: usize, ref_base_length: i32,
) -> AsmAlignP {
    let junction_start = new_junction(chr_start, pos_start, orient_start);
    let junction_end = new_junction(chr_end, pos_end, orient_end);
    let first_ref_bases = if orient_start.is_forward() {
        ref_genome.base_string(chr_start, pos_start - ref_base_length + 1, pos_start)
    } else {
        ref_genome.base_string(chr_start, pos_start, pos_start + ref_base_length - 1)
    };
    let second_ref_bases = if orient_end.is_forward() {
        ref_genome.base_string(chr_end, pos_end - ref_base_length + 1, pos_end)
    } else {
        ref_genome.base_string(chr_end, pos_end, pos_end + ref_base_length - 1)
    };
    let mut first_ext_bases = inserted_bases.to_string();
    let mut second_ext_bases = inserted_bases.to_string();
    if orient_end.is_forward() {
        first_ext_bases.push_str(&second_ref_bases[second_ref_bases.len() - ext_base_length..]);
    } else {
        first_ext_bases.push_str(&second_ref_bases[..ext_base_length]);
    }
    if orient_start.is_forward() {
        second_ext_bases.push_str(&first_ref_bases[first_ref_bases.len() - ext_base_length..]);
    } else {
        second_ext_bases.push_str(&first_ref_bases[..ext_base_length]);
    }
    if orient_start == orient_end {
        first_ext_bases = format!("{}{}", inserted_bases, reverse_complement_str(&first_ext_bases));
        second_ext_bases = format!("{}{}", reverse_complement_str(inserted_bases), reverse_complement_str(&second_ext_bases));
    }
    let (first_assembly_bases, first_junction_index) = if orient_start.is_forward() {
        (format!("{first_ref_bases}{first_ext_bases}"), first_ref_bases.len() as i32 - 1)
    } else {
        (format!("{first_ext_bases}{first_ref_bases}"), first_ext_bases.len() as i32)
    };
    let (second_assembly_bases, second_junction_index) = if orient_end.is_forward() {
        (format!("{second_ref_bases}{second_ext_bases}"), second_ref_bases.len() as i32 - 1)
    } else {
        (format!("{second_ext_bases}{second_ref_bases}"), second_ext_bases.len() as i32)
    };
    // as in Java, one quality array sized for the first assembly is given to both
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let assembly1 = JunctionAssembly::new_test(junction_start, first_assembly_bases.as_bytes(), &base_quals, first_junction_index);
    let assembly2 = JunctionAssembly::new_test(junction_end, second_assembly_bases.as_bytes(), &base_quals, second_junction_index);
    let (first, second) = if junction_start.orient != junction_end.orient {
        (if assembly1.junction().is_forward() { assembly1 } else { assembly2 }, if assembly1.junction().is_reverse() { assembly1 } else { assembly2 })
    } else {
        (assembly1, assembly2)
    };
    let link = AssemblyLink::new(first, second, LinkType::Split, inserted_bases, overlap_bases);
    create_assembly_alignment_link(link)
}

pub fn create_assembly_alignment_link(assembly_link: LinkP) -> AsmAlignP {
    let phase_set = PhaseSet::new(assembly_link);
    AssemblyAlignment::from_phase_set(phase_set)
}

pub fn create_alignment(chromosome: &str, pos_start: i32, pos_end: i32, seg_start: i32, seg_end: i32, cigar: &str) -> AlignP {
    let score = seg_end - seg_start + 1;
    create_alignment_full(chromosome, pos_start, pos_end, false, DEFAULT_MAP_QUAL, score, seg_start, seg_end, cigar, "", "")
}

pub fn create_alignment_full(
    chromosome: &str, pos_start: i32, pos_end: i32, reversed: bool, map_qual: i32, score: i32, seg_start: i32, seg_end: i32, cigar: &str, xa_tag: &str, md_tag: &str,
) -> AlignP {
    let flags = if reversed { FLAG_READ_REVERSE_STRAND } else { 0 };
    P::new(AlignData::new(chromosome, pos_start, pos_end, seg_start, seg_end, map_qual, score, flags, cigar, DEFAULT_NM, Some(xa_tag), Some(md_tag)))
}

pub fn has_assembly_link(links: &[LinkP], assembly1: AsmP, assembly2: AsmP, link_type: LinkType) -> bool {
    links.iter().any(|x| x.type_() == link_type && x.has_assembly(assembly1) && x.has_assembly(assembly2))
}

pub fn assert_close(expected: f64, actual: f64, delta: f64) {
    assert!((expected - actual).abs() <= delta, "expected {expected} but was {actual} (delta {delta})");
}

/// the test constructor JunctionAssembly(junction, bases, quals, junctionIndex) with default qualities
pub fn new_assembly(junction: P<Junction>, assembly_bases: &str, junction_index: usize) -> AsmP {
    JunctionAssembly::new_test(junction, assembly_bases.as_bytes(), &build_default_base_quals(assembly_bases.len()), junction_index as i32)
}

/// AssemblyLinker.tryAssemblyOverlap(assembly1, assembly2)
pub fn try_assembly_overlap(assembly1: AsmP, assembly2: AsmP) -> Option<LinkP> {
    let is_local_link = crate::assembly::phase::phase_set_builder::is_local_assembly_candidate(assembly1, assembly2);
    crate::assembly::phase::assembly_linker::try_assembly_overlap(assembly1, assembly2, true, is_local_link)
}

/// new PhaseSetBuilder(refGenome, new RemoteReadExtractor(null), phaseGroup).buildPhaseSets()
pub fn build_phase_sets(ref_genome: &MockRefGenome, phase_group: crate::assembly::types::phase_group::PhaseGroupP) {
    let mut remote_read_extractor = crate::assembly::phase::remote_read_extractor::RemoteReadExtractor::new(None, None);
    let mut phase_set_builder = crate::assembly::phase::phase_set_builder::PhaseSetBuilder::new(Some(ref_genome), &mut remote_read_extractor, phase_group);
    phase_set_builder.build_phase_sets();
}
