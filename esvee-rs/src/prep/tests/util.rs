//! Test helpers after the Java TestUtils and SamRecordTestUtils.

use crate::prep::filters::ReadFilterConfig;
use crate::prep::rec::*;
use crate::prep::tracker::{JunctionTracker, TrackerConfig};
use crate::prep::types::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, CigarString};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

pub const CHR_1: &str = "1";
pub const CHR_2: &str = "2";
pub const CHR_3: &str = "3";
pub const NO_CHROMOSOME_NAME: &str = "*";
pub const DEFAULT_BASE_QUAL: u8 = 37;
pub const PARTITION_SIZE: i32 = 10000;

/// TestUtils.REF_BASES_200 and REF_BASES_400
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

/// MockRefGenome.generateRandomBases: the four nucleotides in turn
pub fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| b"ACGT"[i % 4] as char).collect()
}

static NEXT_READ_ID: AtomicUsize = AtomicUsize::new(1);

/// READ_ID_GENERATOR.nextId()
pub fn next_id() -> String {
    format!("READ_{:03}", NEXT_READ_ID.fetch_add(1, Ordering::SeqCst))
}

/// readIdStr
pub fn read_id_str(read_id: i32) -> String {
    format!("READ_{read_id:02}")
}

/// The chromosomes of the test dictionary: 1 to 22, X, Y, MT.
pub fn contigs() -> Arc<Contigs> {
    let mut c = Contigs::new();
    for i in 1..=22 {
        c.add(&i.to_string());
    }
    c.add("X");
    c.add("Y");
    c.add("MT");
    Arc::new(c)
}

pub fn tid_map(contigs: &Contigs) -> Vec<i32> {
    (0..contigs.names.len() as i32).collect()
}

pub fn chr(name: &str) -> i32 {
    contigs().id(name)
}

/// buildFlags: paired and proper pair, first or second in pair
pub fn build_flags(first_in_pair: bool, reversed: bool, supplementary: bool) -> u16 {
    let mut flags = FLAG_PAIRED | FLAG_PROPER_PAIR;
    if reversed {
        flags |= FLAG_REVERSE;
    }
    flags |= if first_in_pair { FLAG_FIRST } else { FLAG_SECOND };
    if supplementary {
        flags |= FLAG_SUPPLEMENTARY;
    }
    flags
}

fn base_record(read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str, base_qual: u8) -> bam::Record {
    let mut rec = bam::Record::new();
    let cigar = CigarString::try_from(cigar).expect("cigar");
    rec.set(read_id.as_bytes(), Some(&cigar), read_bases.as_bytes(), &vec![base_qual; read_bases.len()]);
    rec.set_tid(chr(chromosome));
    rec.set_pos(read_start as i64 - 1);
    // the record builder of the Java tests gives a mapped read this mapping quality
    rec.set_mapq(255);
    rec
}

/// TestUtils.createSamRecord with flags and base quality: the mate is on the same chromosome 300 bases on,
/// reversed, with an insert size of 400
pub fn sam_record_flags(read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str, flags: u16, base_qual: u8) -> bam::Record {
    let mut rec = base_record(read_id, chromosome, read_start, read_bases, cigar, base_qual);
    rec.set_flags(flags | FLAG_MATE_REVERSE);
    rec.set_mtid(chr(chromosome));
    rec.set_mpos(read_start as i64 + 300 - 1);
    rec.set_insert_size(400);
    rec
}

/// TestUtils.createSamRecord(readId, chromosome, readStart, readBases, cigar)
pub fn sam_record(read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str) -> bam::Record {
    sam_record_flags(read_id, chromosome, read_start, read_bases, cigar, build_flags(true, false, false), DEFAULT_BASE_QUAL)
}

/// TestUtils.createSamRecord(readId, chromosome, readStart, mateChr, mateStart, firstInPair, isSupp, suppData):
/// a 100M read
#[allow(clippy::too_many_arguments)]
pub fn sam_record_mate(read_id: &str, chromosome: &str, read_start: i32, mate_chr: &str, mate_start: i32, first_in_pair: bool, is_supp: bool, supp_data: Option<&str>) -> bam::Record {
    let mut rec = sam_record_flags(read_id, chromosome, read_start, &REF_BASES_200[..100], "100M", build_flags(first_in_pair, false, is_supp), DEFAULT_BASE_QUAL);
    rec.set_mtid(chr(mate_chr));
    rec.set_mpos(mate_start as i64 - 1);
    if let Some(s) = supp_data.filter(|s| !s.is_empty()) {
        rec.push_aux(b"SA", Aux::String(s)).unwrap();
    }
    rec
}

/// SamRecordTestUtils.createSamRecord with mate orientation and mate cigar. `supp` is the supplementary
/// alignment as (chromosome, position, strand, cigar, mapping quality).
#[allow(clippy::too_many_arguments)]
pub fn common_sam_record(
    read_id: &str, chromosome: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, is_reversed: bool, is_supplementary: bool,
    supp: Option<(&str, i32, char, &str, i32)>, mate_reversed: bool, mate_cigar: Option<&str>,
) -> bam::Record {
    let mut rec = base_record(read_id, chromosome, read_start, read_bases, cigar, DEFAULT_BASE_QUAL);
    let mut flags = FLAG_PAIRED | FLAG_FIRST;
    if is_reversed {
        flags |= FLAG_REVERSE;
    }
    if !mate_chr.is_empty() && mate_chr != NO_CHROMOSOME_NAME {
        rec.set_mtid(chr(mate_chr));
        rec.set_mpos(mate_start as i64 - 1);
        flags |= FLAG_PROPER_PAIR;
    } else {
        flags |= FLAG_MATE_UNMAPPED;
        rec.set_mtid(-1);
        rec.set_mpos(read_start as i64 - 1);
    }
    if is_supplementary {
        flags |= FLAG_SUPPLEMENTARY;
    }
    if mate_reversed {
        flags |= FLAG_MATE_REVERSE;
    }
    rec.set_flags(flags);
    if let Some((c, p, strand, cg, mq)) = supp {
        rec.push_aux(b"SA", Aux::String(&format!("{c},{p},{strand},{cg},{mq},0"))).unwrap();
    }
    rec.set_insert_size(if chromosome == mate_chr { (read_start - mate_start).abs() as i64 } else { 0 });
    if let Some(mc) = mate_cigar {
        rec.push_aux(b"MC", Aux::String(mc)).unwrap();
    }
    rec
}

/// PrepRead.from
pub fn prep_read(rec: bam::Record) -> PrepRead {
    let c = contigs();
    PrepRead::new(rec, 0, &tid_map(&c), &c)
}

pub fn tracker_config() -> TrackerConfig {
    TrackerConfig { filter: ReadFilterConfig::default(), trim_read_id: false, contigs: contigs() }
}

pub fn new_tracker(region: Region) -> JunctionTracker {
    JunctionTracker::new(region, tracker_config(), Vec::new(), None)
}

/// A group of its own for a read held by the tracker, as `new ReadGroup(read)`.
pub fn new_group(tracker: &mut JunctionTracker, read: PrepRead) -> GroupId {
    let id: Box<[u8]> = read.id().into();
    let rid = tracker.add_read(read);
    tracker.groups.push(ReadGroup::new(id, rid));
    (tracker.groups.len() - 1) as GroupId
}

pub fn add_to_group(tracker: &mut JunctionTracker, gid: GroupId, read: PrepRead) -> ReadId {
    let rid = tracker.add_read(read);
    tracker.groups[gid as usize].reads.push(rid);
    rid
}
