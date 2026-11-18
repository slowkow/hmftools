//! Shared fixtures for the unit tests ported from the Java SAGE test suite
//! (sage/src/test/java/.../common/TestUtils.java, VariantUtils.java and hmf-common SamRecordTestUtils).

use crate::config::Config;
use crate::quality::{BqrMap, MsiJitter, QualityCalculator};
use crate::refdata::RefData;
use crate::refseq::RefSequence;
use crate::types::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

pub const DEFAULT_BASE_QUAL: u8 = 37;
pub const DEFAULT_MAP_QUAL: i32 = 60;
pub const CHR_1: &str = "1";
pub const TEST_SAMPLE: &str = "SAMPLE";

/// TestUtils.REF_BASES_200
pub const REF_BASES_200: &str = concat!(
    "CGCAATATTCGGGTGGGAGTGACCCGATTTTCCAGGTGCGTTCGTCACCGCTGTCTGTGACTCGGAAAAAAAACTCCCTGACCCCTTGCGCTTCCCAGGT",
    "GAGGCAATGCCTCGCCCTGCTTCGGCTCGCGCACAGTGCGCGCTACACACACTGGCCTGCGCCCACTGTCTGGCACTCCCTAGTGAGATGAACCCGGTAC"
);

/// TestUtils.REF_SEQUENCE_200: zero-based so positions line up with string indices
pub fn ref_sequence_200() -> RefSequence {
    RefSequence::from_bases(0, REF_BASES_200.as_bytes().to_vec())
}

/// TestUtils.createSageConfig: all defaults
pub fn test_config() -> Config {
    let args: Vec<String> = ["-tumor", "TUMOR_ID", "-ref_genome", "none.fasta", "-ref_genome_version", "37", "-output_vcf", "none.vcf.gz"]
        .iter()
        .map(|x| x.to_string())
        .collect();
    Config::from_args(&args).unwrap()
}

pub fn empty_ref_data() -> RefData {
    RefData {
        panel: HashMap::new(),
        hotspots: HashMap::new(),
        high_confidence: HashMap::new(),
        transcripts: HashMap::new(),
        hla_genes: Vec::new(),
        hla_chromosome: "6".to_string(),
    }
}

/// Owns what a QualityCalculator borrows: TestUtils.QUALITY_CALCULATOR is `fixture.qc()`.
pub struct Fixture {
    pub cfg: Config,
    pub bqr: BqrMap,
    pub ref_seq: RefSequence,
    pub msi: MsiJitter,
    pub ref_data: RefData,
}

impl Fixture {
    pub fn new() -> Self {
        Self::with_ref_sequence(ref_sequence_200())
    }
    pub fn with_ref_sequence(ref_seq: RefSequence) -> Self {
        Fixture {
            cfg: test_config(),
            bqr: BqrMap::default(),
            ref_seq,
            // new MsiJitterCalcs(): no per-sample parameters
            msi: MsiJitter { params: Vec::new(), probable_msi: false },
            ref_data: empty_ref_data(),
        }
    }
    pub fn qc(&self) -> QualityCalculator<'_> {
        QualityCalculator { cfg: &self.cfg, bqr: &self.bqr, ref_seq: &self.ref_seq, msi: &self.msi, ref_data: &self.ref_data, chr: CHR_1 }
    }
}

pub fn parse_cigar(cigar: &str) -> Vec<CigarEl> {
    let mut out = Vec::new();
    let mut len = 0;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            len = len * 10 + d as i32;
            continue;
        }
        let op = match c {
            'M' => Op::M,
            'I' => Op::I,
            'D' => Op::D,
            'N' => Op::N,
            'S' => Op::S,
            'H' => Op::H,
            'P' => Op::P,
            '=' => Op::Eq,
            'X' => Op::X,
            _ => panic!("invalid cigar {cigar}"),
        };
        out.push(CigarEl::new(len, op));
        len = 0;
    }
    out
}

fn alignment_end(start: i32, cigar: &[CigarEl]) -> i32 {
    start + cigar.iter().filter(|e| e.op.consumes_ref()).map(|e| e.len).sum::<i32>() - 1
}

static READ_ID: AtomicUsize = AtomicUsize::new(0);

/// ReadIdGenerator.nextId
pub fn next_read_id() -> String {
    format!("READ_{:03}", READ_ID.fetch_add(1, Ordering::SeqCst))
}

pub fn default_quals(len: usize) -> Vec<u8> {
    vec![DEFAULT_BASE_QUAL; len]
}

/// TestUtils.buildSamRecord(alignmentStart, cigar, readBases, qualities): paired, proper pair, forward
/// strand, map qual 60, insert size 600, NM 0, no mate details.
pub fn build_sam_record_quals(alignment_start: i32, cigar: &str, read_bases: &str, quals: Vec<u8>) -> Read {
    let cigar = parse_cigar(cigar);
    Read {
        name: next_read_id(),
        flags: FLAG_PAIRED | FLAG_PROPER,
        tid: 0,
        start: alignment_start,
        end: alignment_end(alignment_start, &cigar),
        mapq: DEFAULT_MAP_QUAL,
        cigar,
        bases: read_bases.as_bytes().to_vec(),
        quals,
        mate_tid: -1,
        mate_start: 0,
        insert_size: 600,
        nm: Some(0),
        mate_cigar: None,
        consensus: ConsensusType::None,
        has_consensus_tag: false,
        duplex_index: -1,
    }
}

/// TestUtils.buildSamRecord(alignmentStart, cigar, readBases)
pub fn build_sam_record(alignment_start: i32, cigar: &str, read_bases: &str) -> Read {
    build_sam_record_quals(alignment_start, cigar, read_bases, default_quals(read_bases.len()))
}

/// TestUtils.createSamRecord(readId, chr, readStart, readBases, cigar): first of pair, proper pair, mate
/// 300 bases downstream on the negative strand with the same cigar, NM 1.
pub fn create_sam_record(read_id: &str, read_start: i32, read_bases: &str, cigar: &str) -> Read {
    let mut read = build_sam_record(read_start, cigar, read_bases);
    read.name = read_id.to_string();
    read.flags = FLAG_PAIRED | FLAG_PROPER | FLAG_FIRST | FLAG_MATE_REVERSE;
    read.mate_tid = 0;
    read.mate_start = read_start + 300;
    read.mate_cigar = Some(cigar.to_string());
    read.nm = Some(1);
    read.insert_size = 0;
    read
}
