//! Synthetic reads and fixtures for the ported tests: the equivalents of the Java TestUtils, SamRecordTestUtils,
//! MockRefGenome, ReadIdGenerator and TestBamWriter.

use crate::bqr::BaseQualRecalibration;
use crate::capture::ReadCapture;
use crate::config::{Config, FilterReadsType, SequencingType};
use crate::consensus::{ConsensusOutcome, ConsensusReads, RefBases};
use crate::dup::FragmentCoords;
use crate::partition::PartitionReader;
use crate::rec;
use crate::regions::ChrBaseRegion;
use crate::unmap::{ReadUnmapper, UnmapRegionState, UnmapStats, UnmappingRegion, UNMAP_MIN_HIGH_DEPTH};
use rust_htslib::bam::{self, record::{Cigar, CigarString}};
use std::collections::HashMap;
use std::sync::OnceLock;

pub const CHR_1: &str = "1";
pub const CHR_2: &str = "2";
pub const CHR_3: &str = "3";
pub const CHR_4: &str = "4";
pub const CHR_5: &str = "5";
pub const NO_CHROMOSOME_NAME: &str = "*";
pub const NO_CIGAR: &str = "*";
pub const NO_POSITION: i32 = 0;
pub const SUPP_POS_STRAND: char = '+';
pub const SUPP_NEG_STRAND: char = '-';

pub const DEFAULT_QUAL: u8 = 37;
/// SAMRecordSetBuilder gives mapped reads this map quality
pub const BUILDER_MAP_QUAL: u8 = 255;

pub const TEST_READ_ID: &str = "READ_01";
pub const TEST_READ_CIGAR: &str = "100M";

pub const REF_BASES_A: &str = "AAAAAAAAAA";
pub const REF_BASES_C: &str = "CCCCCCCCCC";
pub const REF_BASES_G: &str = "GGGGGGGGGG";
pub const REF_BASES_T: &str = "TTTTTTTTTT";

/// RefGenomeCoordinates.COORDS_37 for chromosomes 1-22, X, Y, then MT
pub const V37_LENGTHS: [(&str, i32); 25] = [
    ("1", 249250621), ("2", 243199373), ("3", 198022430), ("4", 191154276), ("5", 180915260), ("6", 171115067), ("7", 159138663),
    ("8", 146364022), ("9", 141213431), ("10", 135534747), ("11", 135006516), ("12", 133851895), ("13", 115169878),
    ("14", 107349540), ("15", 102531392), ("16", 90354753), ("17", 81195210), ("18", 78077248), ("19", 59128983),
    ("20", 63025520), ("21", 48129895), ("22", 51304566), ("X", 155270560), ("Y", 59373566), ("MT", 16569),
];

/// MockRefGenome.generateRandomBases: cycles through A, C, G, T
pub fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| b"ACGT"[i % 4] as char).collect()
}

pub fn test_read_bases() -> String {
    generate_random_bases(100)
}

/// TestUtils.REF_BASES: 'X' then 100 bases, so that a zero-based lookup is in effect one-based
pub fn ref_bases() -> String {
    let random = generate_random_bases(10);
    format!(
        "X{random}{REF_BASES_A}{REF_BASES_C}{REF_BASES_G}{REF_BASES_T}{REF_BASES_A}{REF_BASES_C}{REF_BASES_G}{REF_BASES_T}{random}"
    )
}

pub fn ref_bases_repeat_40() -> String {
    ref_bases().repeat(40)
}

/// SAM_DICTIONARY_V37 sequence names by index
pub fn names() -> Vec<String> {
    V37_LENGTHS.iter().map(|(n, _)| n.to_string()).collect()
}

pub fn header() -> bam::HeaderView {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (n, l) in V37_LENGTHS {
        text.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
    }
    bam::HeaderView::from_bytes(text.as_bytes())
}

/// SamRecordTestUtils.chromosomeOrdinal
pub fn chromosome_ordinal(chr: &str) -> i32 {
    match V37_LENGTHS.iter().position(|(n, _)| *n == chr) {
        Some(i) => i as i32,
        None => 25,
    }
}

pub struct ReadIdGenerator {
    next: i32,
}

impl ReadIdGenerator {
    pub fn new() -> Self {
        ReadIdGenerator { next: 1 }
    }
    pub fn next_id(&mut self) -> String {
        let id = format!("READ_{:03}", self.next);
        self.next += 1;
        id
    }
}

/// SupplementaryReadData
#[derive(Clone)]
pub struct SuppData {
    pub chromosome: String,
    pub position: i32,
    pub strand: char,
    pub cigar: String,
    pub map_quality: i32,
}

impl SuppData {
    pub fn new(chromosome: &str, position: i32, strand: char, cigar: &str, map_quality: i32) -> Self {
        SuppData { chromosome: chromosome.to_string(), position, strand, cigar: cigar.to_string(), map_quality }
    }
    /// SupplementaryReadData.asSamTag (NM is 0)
    pub fn as_sam_tag(&self) -> String {
        format!("{},{},{},{},{},0", self.chromosome, self.position, self.strand, self.cigar, self.map_quality)
    }
}

/// SupplementaryReadData.alignmentsToSamTag
pub fn alignments_to_sam_tag(alignments: &[SuppData]) -> String {
    alignments.iter().map(|a| a.as_sam_tag()).collect::<Vec<_>>().join(";")
}

pub fn parse_cigar(cigar: &str) -> CigarString {
    let mut out = Vec::new();
    if cigar == NO_CIGAR {
        return CigarString(out);
    }
    let mut len: u32 = 0;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            len = len * 10 + d;
            continue;
        }
        out.push(match c {
            'M' => Cigar::Match(len),
            'I' => Cigar::Ins(len),
            'D' => Cigar::Del(len),
            'N' => Cigar::RefSkip(len),
            'S' => Cigar::SoftClip(len),
            'H' => Cigar::HardClip(len),
            'P' => Cigar::Pad(len),
            '=' => Cigar::Equal(len),
            'X' => Cigar::Diff(len),
            other => panic!("invalid cigar operator {other}"),
        });
        len = 0;
    }
    CigarString(out)
}

/// SamRecordTestUtils.createSamRecord (10 arguments): a paired, first-in-pair read as SAMRecordSetBuilder.addFrag makes
/// it (read group attribute, map quality 255), with base qualities of 37.
#[allow(clippy::too_many_arguments)]
pub fn create_sam_record_basic(
    read_id: &str, chr: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, is_reversed: bool,
    is_supplementary: bool, supp_alignment: Option<&SuppData>,
) -> bam::Record {
    let chr_index = if chr != NO_CHROMOSOME_NAME { chromosome_ordinal(chr) } else { -1 };
    let mut record = bam::Record::new();
    let quals = vec![DEFAULT_QUAL; read_bases.len()];
    record.set(read_id.as_bytes(), Some(&parse_cigar(cigar)), read_bases.as_bytes(), &quals);
    record.set_tid(chr_index);
    // the builder sets the alignment start only with a contig
    record.set_pos(if chr_index >= 0 { read_start as i64 - 1 } else { -1 });
    record.set_mapq(BUILDER_MAP_QUAL);
    record.set_insert_size(0);
    let mut flags = rec::FLAG_PAIRED | rec::FLAG_FIRST;
    if is_reversed {
        flags |= rec::FLAG_REVERSE;
    }
    if chr_index < 0 {
        flags |= rec::FLAG_UNMAPPED;
    }
    if !mate_chr.is_empty() && mate_chr != NO_CHROMOSOME_NAME {
        record.set_mtid(chromosome_ordinal(mate_chr));
        record.set_mpos(mate_start as i64 - 1);
        flags |= rec::FLAG_PROPER_PAIR;
    } else {
        flags |= rec::FLAG_MATE_UNMAPPED;
        record.set_mtid(-1);
        record.set_mpos(read_start as i64 - 1);
    }
    if is_supplementary {
        flags |= rec::FLAG_SUPPLEMENTARY;
    }
    record.set_flags(flags);
    set_string_attribute(&mut record, b"RG", "1");
    if let Some(supp) = supp_alignment {
        set_string_attribute(&mut record, b"SA", &supp.as_sam_tag());
    }
    if chr == mate_chr {
        record.set_insert_size((read_start - mate_start).abs() as i64);
    }
    rec::fix_bin(&mut record);
    record
}

/// SamRecordTestUtils.createSamRecord (12 arguments)
#[allow(clippy::too_many_arguments)]
pub fn create_sam_record(
    read_id: &str, chr: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, is_reversed: bool,
    is_supplementary: bool, supp_alignment: Option<&SuppData>, mate_reversed: bool, mate_cigar: Option<&str>,
) -> bam::Record {
    let mut record =
        create_sam_record_basic(read_id, chr, read_start, read_bases, cigar, mate_chr, mate_start, is_reversed, is_supplementary, supp_alignment);
    if mate_reversed {
        rec::set_flag(&mut record, rec::FLAG_MATE_REVERSE, true);
    }
    if let Some(mc) = mate_cigar {
        set_string_attribute(&mut record, b"MC", mc);
    }
    record
}

/// SamRecordTestUtils.createSamRecordUnpaired
#[allow(clippy::too_many_arguments)]
pub fn create_sam_record_unpaired(
    read_id: &str, chr: &str, read_start: i32, read_bases: &str, cigar: &str, is_reversed: bool, is_supplementary: bool,
    supp_alignment: Option<&SuppData>,
) -> bam::Record {
    let mut record = create_sam_record_basic(
        read_id, chr, read_start, read_bases, cigar, NO_CHROMOSOME_NAME, NO_POSITION, is_reversed, is_supplementary, supp_alignment,
    );
    rec::set_flag(&mut record, rec::FLAG_PAIRED, false);
    record
}

/// SamRecordTestUtils.cloneSamRecord
pub fn clone_sam_record(record: &bam::Record, new_read_id: &str) -> bam::Record {
    let mut r = record.clone();
    r.set_qname(new_read_id.as_bytes());
    r
}

/// SamRecordTestUtils.flipFirstInPair
pub fn flip_first_in_pair(record: &mut bam::Record) {
    if rec::has(record, rec::FLAG_PAIRED) {
        let first = rec::has(record, rec::FLAG_FIRST);
        set_first_of_pair(record, !first);
    }
}

/// setFirstOfPairFlag(first) and setSecondOfPairFlag(!first)
pub fn set_first_of_pair(record: &mut bam::Record, first: bool) {
    rec::set_flag(record, rec::FLAG_FIRST, first);
    rec::set_flag(record, rec::FLAG_SECOND, !first);
}

/// TestUtils.setSecondInPair
pub fn set_second_in_pair(record: &mut bam::Record) {
    set_first_of_pair(record, false);
}

/// TestUtils.setBaseQualities
pub fn set_base_qualities(record: &mut bam::Record, value: u8) {
    let name = record.qname().to_vec();
    let cigar = CigarString(record.cigar().iter().copied().collect());
    let seq = record.seq().as_bytes();
    let quals = vec![value; seq.len()];
    record.set(&name, Some(&cigar), &seq, &quals);
}

/// SAMRecord.setAlignmentStart
pub fn set_alignment_start(record: &mut bam::Record, start: i32) {
    record.set_pos(start as i64 - 1);
}

/// SAMRecord.setMateAlignmentStart
pub fn set_mate_alignment_start(record: &mut bam::Record, start: i32) {
    record.set_mpos(start as i64 - 1);
}

pub fn set_string_attribute(record: &mut bam::Record, tag: &[u8; 2], value: &str) {
    rec::set_attribute(record, tag, Some(&rec::encode_z(value.as_bytes())));
}

pub fn string_attribute(record: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    rec::aux_str(record, tag).map(|v| String::from_utf8_lossy(v).to_string())
}

pub fn read_name(record: &bam::Record) -> String {
    String::from_utf8_lossy(record.qname()).to_string()
}

pub fn read_string(record: &bam::Record) -> String {
    String::from_utf8(record.seq().as_bytes()).unwrap()
}

pub fn cigar_string(record: &bam::Record) -> String {
    rec::cigar_string(record.raw_cigar())
}

/// SAMRecord.getReferenceName
pub fn reference_name(record: &bam::Record) -> String {
    if record.tid() >= 0 { V37_LENGTHS[record.tid() as usize].0.to_string() } else { NO_CHROMOSOME_NAME.to_string() }
}

/// SAMRecord.getMateReferenceName
pub fn mate_reference_name(record: &bam::Record) -> String {
    if record.mtid() >= 0 { V37_LENGTHS[record.mtid() as usize].0.to_string() } else { NO_CHROMOSOME_NAME.to_string() }
}

/// TestUtils.createFragmentCoords
pub fn create_fragment_coords(record: &bam::Record) -> FragmentCoords {
    FragmentCoords::from_read(record, &names(), false)
}

/// MockRefGenome: bases by chromosome held in memory. Chromosome lengths default to the v37 lengths.
#[derive(Clone)]
pub struct MockRefGenome {
    one_based: bool,
    pub ref_genome_map: HashMap<String, String>,
    pub chromosome_lengths: HashMap<String, i32>,
}

impl MockRefGenome {
    pub fn new(one_based: bool) -> Self {
        let mut chromosome_lengths = HashMap::new();
        // the mitochondrial chromosome is not part of the coordinates
        for (n, l) in &V37_LENGTHS[..24] {
            chromosome_lengths.insert(n.to_string(), *l);
        }
        MockRefGenome { one_based, ref_genome_map: HashMap::new(), chromosome_lengths }
    }
    /// consensus.RefGenome: mChromosomeLengths.getOrDefault(chromosome, 0)
    pub fn chromosome_length_or_zero(&self, chr: &str) -> i32 {
        self.chromosome_lengths.get(chr).copied().unwrap_or(0)
    }
    /// MockRefGenome.getBases: empty outside of the bases held
    pub fn bases(&self, chr: &str, mut pos_start: i32, mut pos_end: i32) -> Vec<u8> {
        if self.one_based {
            pos_start -= 1;
            pos_end -= 1;
        }
        match self.ref_genome_map.get(chr) {
            Some(b) if pos_start >= 0 && b.len() as i32 > pos_end => b.as_bytes()[pos_start as usize..=pos_end as usize].to_vec(),
            _ => Vec::new(),
        }
    }
}

/// A small FASTA whose index claims the v37 lengths of chromosomes 1 and 2. It stands in where the port needs a
/// reference file to open: no test reads bases from it.
pub fn placeholder_fasta() -> &'static str {
    static PATH: OnceLock<String> = OnceLock::new();
    PATH.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("redux_rs_tests_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("placeholder.fa");
        std::fs::write(&fasta, ">1\nACGT\n>2\nACGT\n").unwrap();
        let fai = format!("1\t{}\t3\t4\t5\n2\t{}\t11\t4\t5\n", V37_LENGTHS[0].1, V37_LENGTHS[1].1);
        std::fs::write(dir.join("placeholder.fa.fai"), fai).unwrap();
        fasta.to_string_lossy().to_string()
    })
}

pub fn mock_ref_bases(mock: MockRefGenome) -> RefBases {
    RefBases::with_mock(placeholder_fasta(), mock)
}

/// new ConsensusReads(refGenome)
pub fn consensus_reads(mock: &MockRefGenome) -> ConsensusReads {
    ConsensusReads::with_refs(mock_ref_bases(mock.clone()), SequencingType::Illumina)
}

/// ConsensusReadInfo
pub struct ConsensusReadInfo {
    pub consensus_read: bam::Record,
    /// name of the template read
    pub template_read: String,
    pub outcome: ConsensusOutcome,
}

/// ConsensusReads.createConsensusRead(reads, fragmentCoords, umiId). The outcome is taken from the statistics, which
/// count every outcome but the single read one.
pub fn create_consensus_read_with_coords(
    consensus_reads: &mut ConsensusReads, reads: &[bam::Record], coords: &FragmentCoords, umi_id: Option<&str>,
) -> ConsensusReadInfo {
    let names = names();
    let mut for_template: Vec<bam::Record> = reads.to_vec();
    let template_index = crate::consensus::select_template_read(&mut for_template, coords);
    let template_read = read_name(&for_template[template_index]);
    let before = consensus_reads.stats.outcome_counts;
    let mut reads: Vec<bam::Record> = reads.to_vec();
    let consensus_read = consensus_reads.create_consensus_read(&mut reads, coords, umi_id, &names);
    let after = consensus_reads.stats.outcome_counts;
    const OUTCOMES: [ConsensusOutcome; 7] = [
        ConsensusOutcome::Unset,
        ConsensusOutcome::AlignmentOnly,
        ConsensusOutcome::IndelMatch,
        ConsensusOutcome::IndelMismatch,
        ConsensusOutcome::IndelSoftclip,
        ConsensusOutcome::IndelFail,
        ConsensusOutcome::SingleRead,
    ];
    let changed: Vec<usize> = (0..7).filter(|&i| after[i] != before[i]).collect();
    let outcome = match changed.as_slice() {
        [] => ConsensusOutcome::SingleRead,
        [i] => OUTCOMES[*i],
        _ => panic!("more than one consensus outcome registered"),
    };
    ConsensusReadInfo { consensus_read, template_read, outcome }
}

/// TestUtils.createConsensusRead: the fragment coordinates of the first read
pub fn create_consensus_read(consensus_reads: &mut ConsensusReads, reads: &[bam::Record], umi_id: &str) -> ConsensusReadInfo {
    let coords = create_fragment_coords(&reads[0]);
    create_consensus_read_with_coords(consensus_reads, reads, &coords, Some(umi_id))
}

/// TestUtils.CHR_LOCATION_MAP
pub fn chr_location_map() -> HashMap<String, Vec<UnmappingRegion>> {
    let region = |start, end, max_depth| UnmappingRegion { start, end, max_depth };
    let mut map = HashMap::new();
    map.insert(CHR_1.to_string(), vec![region(500, 700, 0)]);
    map.insert(CHR_2.to_string(), vec![]);
    map.insert(CHR_3.to_string(), vec![region(500, 700, UNMAP_MIN_HIGH_DEPTH)]);
    map.insert(
        CHR_4.to_string(),
        vec![
            region(1000, 2000, UNMAP_MIN_HIGH_DEPTH),
            region(3000, 4000, UNMAP_MIN_HIGH_DEPTH),
            region(5000, 6000, UNMAP_MIN_HIGH_DEPTH),
            region(7000, 8000, UNMAP_MIN_HIGH_DEPTH),
        ],
    );
    map
}

/// new ReadUnmapper(chrLocationsMap)
pub fn read_unmapper(map: HashMap<String, Vec<UnmappingRegion>>) -> ReadUnmapper {
    ReadUnmapper::build(&header(), map)
}

/// TestUtils.READ_UNMAPPER
pub fn test_read_unmapper() -> ReadUnmapper {
    read_unmapper(chr_location_map())
}

/// new UnmapRegionState(new ChrBaseRegion(chromosome, start, end), regions)
pub fn unmap_region_state(chromosome: &str, start: i32, end: i32, regions: Vec<UnmappingRegion>) -> UnmapRegionState {
    UnmapRegionState { chromosome: chromosome.to_string(), slice_start: start, slice_end: end, regions, last_matched_index: None }
}

/// TestUtils.checkTransformRead
pub fn check_transform_read(read: &mut bam::Record, chromosome: &str) -> bool {
    let unmapper = test_read_unmapper();
    let mut state = unmap_region_state(chromosome, 1, 1000000, chr_location_map().remove(chromosome).unwrap_or_default());
    unmapper.check_transform_read(read, &mut state, &mut UnmapStats::default())
}

/// What a PartitionReader borrows: the test ReduxConfig, the sequence dictionary and the unmapper.
pub struct TestEnv {
    pub config: Config,
    pub names: Vec<String>,
    pub expected_valid_tid: Vec<bool>,
    pub unmapper: ReadUnmapper,
    pub bqr: BaseQualRecalibration,
    pub ref_genome: MockRefGenome,
}

impl TestEnv {
    /// new ReduxConfig(refGenome, false, false, formConsensus, readUnmapper)
    pub fn new(ref_genome: MockRefGenome, form_consensus: bool, unmapper: ReadUnmapper) -> Self {
        let mut args: Vec<String> =
            ["-sample", "TEST", "-input_bam", "/redux_rs_tests/none.bam", "-ref_genome", placeholder_fasta(), "-bqr_disable"]
                .iter()
                .map(|s| s.to_string())
                .collect();
        if form_consensus {
            args.push("-form_consensus".to_string());
        }
        let mut config = Config::from_args(&args).unwrap();
        config.specific_region_filter = FilterReadsType::MateAndSupp;
        config.write_bam = false;
        let genome = crate::refgenome::RefGenome::open(placeholder_fasta()).unwrap();
        let bqr = BaseQualRecalibration::new(&config, &genome).unwrap();
        let names = names();
        let expected_valid_tid = vec![true; names.len()];
        TestEnv { config, names, expected_valid_tid, unmapper, bqr, ref_genome }
    }

    /// TestUtils.createTestConfig
    pub fn test_config() -> Self {
        Self::new(MockRefGenome::new(false), false, read_unmapper(HashMap::new()))
    }

    /// TestUtils.createPartitionRead with a TestBamWriter
    pub fn partition_reader(&self) -> TestPartitionReader<'_> {
        let capture = ReadCapture { jitter: None, bqr: self.bqr.region_reader() };
        let mut reader =
            PartitionReader::new(&self.config, &self.names, &self.expected_valid_tid, &self.unmapper, None, None, capture).unwrap();
        *reader.consensus_mut() = ConsensusReads::with_refs(mock_ref_bases(self.ref_genome.clone()), SequencingType::Illumina);
        TestPartitionReader { reader }
    }
}

/// A PartitionReader whose written records are collected in order, as the Java TestBamWriter does.
pub struct TestPartitionReader<'a> {
    pub reader: PartitionReader<'a>,
}

impl TestPartitionReader<'_> {
    pub fn setup_region(&mut self, chromosome: &str, start: i32, end: i32) {
        self.reader.setup_region(&ChrBaseRegion::new(chromosome, start, end));
    }
    /// PartitionReader.processRead. Returns false when the reader kept the record (it is then cached or written, and
    /// `read` is left empty); true when the read was dropped, in which case `read` shows any changes made to it.
    pub fn process_read(&mut self, read: &mut bam::Record) -> bool {
        self.reader.process_sam_record(read).unwrap();
        !read.qname().is_empty()
    }
    pub fn flush_read_positions(&mut self) {
        self.reader.flush_read_positions().unwrap();
    }
    pub fn post_process_region(&mut self) {
        self.reader.post_process_region().unwrap();
    }
    pub fn written_records(&self) -> &[bam::Record] {
        &self.reader.test_written
    }
    /// BamWriter.consensusWriteCount
    pub fn consensus_write_count(&self) -> usize {
        self.written_records().iter().filter(|r| rec::has_aux(r, b"CR")).count()
    }
    /// BamWriter.nonConsensusWriteCount (duplicates are not dropped in these tests)
    pub fn non_consensus_write_count(&self) -> usize {
        self.written_records().len() - self.consensus_write_count()
    }
    /// The written record with a read name and matching the flags given
    pub fn written(&self, name: &str, flags_set: u16, flags_unset: u16) -> &bam::Record {
        let matches: Vec<&bam::Record> = self
            .written_records()
            .iter()
            .filter(|r| r.qname() == name.as_bytes() && r.flags() & flags_set == flags_set && r.flags() & flags_unset == 0)
            .collect();
        assert_eq!(1, matches.len(), "written records for {name}");
        matches[0]
    }
}
