//! Ports of the Java unit tests for SBX and unpaired reads (hmftools redux/src/test/java/com/hartwig/hmftools/redux):
//! consensus/SbxDuplexReadTest, consensus/SbxConsensusTest, consensus/SbxFinaliseReadTest,
//! duplicate/UnpairedReadDuplicateTest and unmap/UnmapUnpairedReadsTest. Each test names the Java test it ports.
//! The helpers at the top stand in for SamRecordTestUtils, TestUtils and MockRefGenome.

use crate::config::SequencingType;
use crate::consensus::{ConsensusOutcome, ConsensusReads, RefBases};
use crate::dup::*;
use crate::rec::*;
use crate::sbx::*;
use crate::sbx_dup::SbxDuplicateCollapser;
use crate::unmap::{ReadUnmapper, UnmapRegionState, UnmapStats, UnmappingRegion, UNMAP_MIN_HIGH_DEPTH};
use rust_htslib::bam::{self, record::{Cigar, CigarString}};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard};

const CHR_1: &str = "1";
const CHR_2: &str = "2";
const CHR_3: &str = "3";
const TEST_READ_ID: &str = "READ_01";
const TEST_READ_CIGAR: &str = "100M";
const DEFAULT_BASE_QUAL: u8 = 37;
const DEFAULT_MAP_QUAL: u8 = 60;

fn names() -> Vec<String> {
    ["1", "2", "3", "4", "5"].iter().map(|s| s.to_string()).collect()
}

fn header() -> bam::HeaderView {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for n in names() {
        text.push_str(&format!("@SQ\tSN:{n}\tLN:1000000\n"));
    }
    bam::HeaderView::from_bytes(text.as_bytes())
}

fn chromosome_ordinal(chr: &str) -> i32 {
    names().iter().position(|n| n == chr).map(|i| i as i32).unwrap_or(-1)
}

static READ_ID: AtomicUsize = AtomicUsize::new(0);

/// ReadIdGenerator.nextId
fn next_id() -> String {
    format!("READ_{:03}", READ_ID.fetch_add(1, Ordering::Relaxed))
}

fn parse_cigar(cigar: &str) -> CigarString {
    let mut out = Vec::new();
    let mut len = 0u32;
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
            'X' => Cigar::Diff(len),
            other => panic!("invalid cigar operator {other}"),
        });
        len = 0;
    }
    CigarString(out)
}

fn raw_cigar(cigar: &str) -> Vec<u32> {
    let mut record = bam::Record::new();
    record.set(b"cigar", Some(&parse_cigar(cigar)), b"", b"");
    record.raw_cigar().to_vec()
}

fn set_string_attribute(rec: &mut bam::Record, tag: &[u8; 2], value: &str) {
    set_attribute(rec, tag, Some(&encode_z(value.as_bytes())));
}

fn set_int_attribute(rec: &mut bam::Record, tag: &[u8; 2], value: i64) {
    set_attribute(rec, tag, Some(&encode_int(value)));
}

/// SupplementaryReadData.asSamTag with NM 0
fn supp_data(chr: &str, position: i32, strand: char, cigar: &str, map_quality: i32) -> String {
    format!("{chr},{position},{strand},{cigar},{map_quality},0")
}

/// SamRecordTestUtils.createSamRecordUnpaired: base qualities of 37
#[allow(clippy::too_many_arguments)]
fn create_sam_record_unpaired(
    read_id: &str, chr: &str, read_start: i32, read_bases: &str, cigar: &str, is_reversed: bool, is_supplementary: bool,
    supp_alignment: Option<&str>,
) -> bam::Record {
    let mut record = bam::Record::new();
    let quals = vec![DEFAULT_BASE_QUAL; read_bases.len()];
    record.set(read_id.as_bytes(), Some(&parse_cigar(cigar)), read_bases.as_bytes(), &quals);
    record.set_tid(chromosome_ordinal(chr));
    record.set_pos(read_start as i64 - 1);
    record.set_mapq(DEFAULT_MAP_QUAL);
    record.set_mtid(-1);
    record.set_mpos(-1);
    record.set_insert_size(0);
    let mut flags = 0;
    if is_reversed {
        flags |= FLAG_REVERSE;
    }
    if is_supplementary {
        flags |= FLAG_SUPPLEMENTARY;
    }
    record.set_flags(flags);
    set_string_attribute(&mut record, b"RG", "1");
    if let Some(sa) = supp_alignment {
        set_string_attribute(&mut record, b"SA", sa);
    }
    fix_bin(&mut record);
    record
}

/// SamRecordTestUtils.createSamRecord: paired, first in pair, with the mate's position, orientation and cigar
#[allow(clippy::too_many_arguments)]
fn create_sam_record_paired(
    read_id: &str, chr: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, mate_reversed: bool,
    mate_cigar: &str,
) -> bam::Record {
    let mut record = create_sam_record_unpaired(read_id, chr, read_start, read_bases, cigar, false, false, None);
    let mut flags = FLAG_PAIRED | FLAG_FIRST | FLAG_PROPER_PAIR;
    if mate_reversed {
        flags |= FLAG_MATE_REVERSE;
    }
    record.set_flags(flags);
    record.set_mtid(chromosome_ordinal(mate_chr));
    record.set_mpos(mate_start as i64 - 1);
    if chr == mate_chr {
        record.set_insert_size((read_start - mate_start).abs() as i64);
    }
    set_string_attribute(&mut record, b"MC", mate_cigar);
    record
}

/// TestUtils.createUnpairedRecord: bases of A over the positions given
fn create_unpaired_record(chromosome: &str, read_start: i32, read_end: i32, is_reversed: bool) -> bam::Record {
    let read_length = (read_end - read_start + 1) as usize;
    create_sam_record_unpaired(&next_id(), chromosome, read_start, &"A".repeat(read_length), &format!("{read_length}M"), is_reversed, false, None)
}

fn clone_sam_record(record: &bam::Record, new_read_id: &str) -> bam::Record {
    let mut r = record.clone();
    r.set_qname(new_read_id.as_bytes());
    r
}

fn build_base_quals(length: usize, qual: u8) -> Vec<u8> {
    vec![qual; length]
}

fn set_base_qualities(rec: &mut bam::Record, quals: &[u8]) {
    quals_mut(rec).copy_from_slice(quals);
}

fn set_cigar(rec: &mut bam::Record, cigar: &str) {
    let name = rec.qname().to_vec();
    let seq = rec.seq().as_bytes();
    let qual = rec.qual().to_vec();
    rec.set(&name, Some(&parse_cigar(cigar)), &seq, &qual);
}

fn read_string(rec: &bam::Record) -> String {
    String::from_utf8(rec.seq().as_bytes()).unwrap()
}

fn cigar_str(rec: &bam::Record) -> String {
    cigar_string(rec.raw_cigar())
}

fn int_attribute(rec: &bam::Record, tag: &[u8; 2]) -> i64 {
    aux_int(rec, tag).unwrap()
}

static FASTA_ID: AtomicUsize = AtomicUsize::new(0);

/// A reference genome of one chromosome, for MockRefGenome with one-based positions: position 1 is the first base.
fn ref_fasta(chromosome: &str, bases: &str) -> String {
    let dir = std::env::temp_dir().join(format!("redux_rs_sbx_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("ref_{}.fa", FASTA_ID.fetch_add(1, Ordering::Relaxed)));
    std::fs::write(&path, format!(">{chromosome}\n{bases}\n")).unwrap();
    path.to_string_lossy().to_string()
}

/// new ConsensusReads(refGenome, SBX)
fn sbx_consensus_reads(chromosome: &str, bases: &str) -> ConsensusReads {
    ConsensusReads::new(&ref_fasta(chromosome, bases), SequencingType::Sbx).unwrap()
}

struct ConsensusReadInfo {
    consensus_read: bam::Record,
    outcome: ConsensusOutcome,
}

/// TestUtils.createConsensusRead: the fragment coordinates of the first read. The outcome is taken from the
/// statistics, which count every outcome but the single read one.
fn create_consensus_read(consensus_reads: &mut ConsensusReads, reads: &[bam::Record]) -> ConsensusReadInfo {
    let names = names();
    let coords = FragmentCoords::from_read(&reads[0], &names, false);
    let before = consensus_reads.stats.outcome_counts;
    let mut group = DuplicateGroup::new(reads.to_vec(), coords);
    consensus_reads.form_consensus_read_plain(&mut group, &names);
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
    ConsensusReadInfo { consensus_read: group.consensus_read.unwrap(), outcome }
}

/// The tests that depend on SbxRoutines.SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL hold this while they run.
static HOMOPOLYMER_SETTING: Mutex<()> = Mutex::new(());

fn with_homopolymer_5_prime_low_base_qual(value: bool) -> MutexGuard<'static, ()> {
    let guard = HOMOPOLYMER_SETTING.lock().unwrap_or_else(|e| e.into_inner());
    set_homopolymer_5_prime_low_base_qual(value);
    guard
}

// ---------------------------------------------------------------------------------------------------------
// SbxDuplexReadTest (SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL is false in this class)

/// SbxDuplexReadTest.checkBaseQuals
fn check_base_quals(read: &bam::Record, low_qual_bases: &[usize]) {
    for (i, &q) in read.qual().iter().enumerate() {
        if low_qual_bases.contains(&i) {
            assert!(q == RAW_DUPLEX_MISMATCH_QUAL || q == SBX_DUPLEX_MISMATCH_QUAL, "base {i} qual {q}");
        } else {
            assert_eq!(RAW_DUPLEX_QUAL, q, "base {i}");
        }
    }
}

/// An unpaired forward read at chromosome 1 with duplex qualities but for the indices given, a YC tag and,
/// when given, NM and AS; then stripDuplexIndels
fn strip_read(
    alignment_start: i32, read_str: &str, cigar: &str, low_quals: &[(usize, u8)], yc_tag: &str, nm: Option<i64>, alignment_score: Option<i64>,
    is_reversed: bool,
) -> bam::Record {
    let mut read = create_sam_record_unpaired(TEST_READ_ID, CHR_1, alignment_start, read_str, cigar, is_reversed, false, None);
    let mut base_quals = build_base_quals(read_str.len(), RAW_DUPLEX_QUAL);
    for &(i, q) in low_quals {
        base_quals[i] = q;
    }
    set_base_qualities(&mut read, &base_quals);
    set_string_attribute(&mut read, b"YC", yc_tag);
    if let Some(nm) = nm {
        set_int_attribute(&mut read, b"NM", nm);
    }
    if let Some(score) = alignment_score {
        set_int_attribute(&mut read, b"AS", score);
    }
    strip_duplex_indels(&mut read).unwrap();
    read
}

fn low(indices: std::ops::RangeInclusive<usize>, qual: u8) -> Vec<(usize, u8)> {
    indices.map(|i| (i, qual)).collect()
}

#[test]
fn test_homopolymer_low_base_qual_direction() {
    // SbxDuplexReadTest.testHomopolymerLowBaseQualDirection
    //            01234567890123456789
    let bases = b"ACGTACGTAAAAACGTACGT";
    let mut quals = build_base_quals(bases.len(), RAW_DUPLEX_QUAL);
    assert_eq!(None, is_homopolymer_low_base_qual_at_start(bases, &quals));

    quals[8] = SBX_DUPLEX_MISMATCH_QUAL;
    assert_eq!(Some(true), is_homopolymer_low_base_qual_at_start(bases, &quals));

    quals[8] = RAW_DUPLEX_QUAL;
    quals[11] = SBX_DUPLEX_MISMATCH_QUAL;
    quals[12] = SBX_DUPLEX_MISMATCH_QUAL;
    assert_eq!(Some(false), is_homopolymer_low_base_qual_at_start(bases, &quals));

    // unclear on a single base
    quals[11] = RAW_DUPLEX_QUAL;
    quals[12] = RAW_DUPLEX_QUAL;
    quals[15] = SBX_DUPLEX_MISMATCH_QUAL;
    assert_eq!(None, is_homopolymer_low_base_qual_at_start(bases, &quals));
}

#[test]
fn test_strip_duplex_indels() {
    // SbxDuplexReadTest.testStripDuplexIndels
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let alignment_start = 10;

    // no stripping for a SNV
    let read = strip_read(alignment_start, "CTCACC", "6M", &[(3, RAW_DUPLEX_MISMATCH_QUAL)], "+3Z2+", Some(1), Some(10), false);
    assert_eq!("6M", cigar_str(&read));
    assert_eq!("CTCACC", read_string(&read));
    check_base_quals(&read, &[3]);
    assert_eq!(1, int_attribute(&read, b"NM"));

    // test 0: 1 low-qual base not in a repeat
    let read = strip_read(alignment_start, "CTCTACC", "3M1I3M", &[(3, 0)], "+3Z3+0", Some(1), Some(10), false);
    assert_eq!("6M", cigar_str(&read));
    assert_eq!("CTCACC", read_string(&read));
    check_base_quals(&read, &[]);
    assert_eq!(0, int_attribute(&read, b"NM"));

    // test 2: as before in 2-base repeat
    let read = strip_read(alignment_start, "CTCTTTACC", "3M2I4M", &[(5, 0)], "+5Z3+0", Some(2), Some(10), false);
    assert_eq!("3M1I4M", cigar_str(&read));
    assert_eq!("CTCTTACC", read_string(&read));
    check_base_quals(&read, &[3, 4]);
    assert_eq!(1, int_attribute(&read, b"NM"));

    // test 3: 2 low-qual bases and 2-base repeat
    let read = strip_read(alignment_start, "CTCTTTTACC", "3M2I5M", &[(5, 0), (6, 0)], "+5ZZ3+0", Some(2), Some(10), false);
    assert_eq!("8M", cigar_str(&read));
    assert_eq!("CTCTTACC", read_string(&read));
    check_base_quals(&read, &[3, 4]);
    assert_eq!(0, int_attribute(&read, b"NM"));

    // test 4: 2 low-qual bases and 2-base repeat
    let read = strip_read(
        alignment_start, "CTCTTTTTACC", "3M3I5M", &low(5..=6, RAW_DUPLEX_MISMATCH_QUAL), "0+6ZZ3+0", Some(3), Some(10), false,
    );
    assert_eq!("3M1I5M", cigar_str(&read));
    assert_eq!("CTCTTTACC", read_string(&read));
    check_base_quals(&read, &[3, 5]);
    assert_eq!(1, int_attribute(&read, b"NM"));

    // test 5: 2 low-qual bases and 2-base repeat with A not part of the repeat
    let read = strip_read(
        alignment_start, "CTCATTTTTACC", "3M4I5M", &low(6..=7, RAW_DUPLEX_MISMATCH_QUAL), "0+7ZZ3+0", Some(4), Some(10), false,
    );
    assert_eq!("3M2I5M", cigar_str(&read));
    assert_eq!("CTCATTTACC", read_string(&read));
    check_base_quals(&read, &[4, 6]);
    assert_eq!(2, int_attribute(&read, b"NM"));
}

#[test]
fn test_strip_duplex_complex_indels() {
    // SbxDuplexReadTest.testStripDuplexComplexIndels
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let alignment_start = 10;

    // test 1: non-homopolymer repeat
    let read = strip_read(
        alignment_start, "CTCATATATACC", "3M4I5M", &low(7..=8, RAW_DUPLEX_MISMATCH_QUAL), "+7ZZ3+", Some(4), Some(10), false,
    );
    assert_eq!("3M2I5M", cigar_str(&read));
    assert_eq!("CTCATATACC", read_string(&read));
    check_base_quals(&read, &[3, 6]);
    assert_eq!(2, int_attribute(&read, b"NM"));

    // test 2: requiring left alignment after trimming
    let read = strip_read(alignment_start, "CTCAAAAAGCTC", "7M2I3M", &[(8, RAW_DUPLEX_MISMATCH_QUAL)], "0+8Z3+", Some(2), None, false);
    assert_eq!("3M1I7M", cigar_str(&read));
    assert_eq!("CTCAAAAACTC", read_string(&read));
    check_base_quals(&read, &[]);
    assert_eq!(1, int_attribute(&read, b"NM"));

    // test 3: low-qual base is right aligned
    let read = strip_read(
        alignment_start, "CTTTATTTTTACC", "5M1I7M", &[(4, RAW_DUPLEX_MISMATCH_QUAL), (9, RAW_DUPLEX_MISMATCH_QUAL)], "0+4Z4Z3+", Some(2), None,
        false,
    );
    assert_eq!("12M", cigar_str(&read));
    assert_eq!("CTTTATTTTACC", read_string(&read));
    check_base_quals(&read, &[4, 5, 8]);
    assert_eq!(1, int_attribute(&read, b"NM"));
}

#[test]
fn test_strip_duplex_complex_indels2() {
    // SbxDuplexReadTest.testStripDuplexComplexIndels2
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let alignment_start = 10;

    // test 1: long homopolymer repeat with low qual bases starting mid-way through
    let read_str = format!("ATGG{}ATGG", "T".repeat(16));
    let yc = format!("0+7{}4+", "Z".repeat(13));
    let read = strip_read(alignment_start, &read_str, "4M13I7M", &low(7..=19, RAW_DUPLEX_MISMATCH_QUAL), &yc, None, None, false);
    assert_eq!("11M", cigar_str(&read));
    assert_eq!("ATGGTTTATGG", read_string(&read));
    check_base_quals(&read, &[4, 5, 6]);

    // test 2: 7 Ts in read, 4 Ts in ref, 3 low qual Ts, 3 inserted Ts
    let yc = format!("0+8{}4+", "Z".repeat(3));
    let read = strip_read(alignment_start, "ATGGTTTTTTTATGG", "4M3I8M", &low(4..=7, RAW_DUPLEX_MISMATCH_QUAL), &yc, None, None, false);
    assert_eq!("12M", cigar_str(&read));
    assert_eq!("ATGGTTTTATGG", read_string(&read));
    check_base_quals(&read, &[4, 5, 6, 7]);

    // test 3: long poly-TC, 3 counts in read, 8 counts in ref, the other counts are low qual and are inserted
    let read_str = format!("ATGG{}ATGG", "TC".repeat(8));
    let yc = format!("0+10{}4+", "Z".repeat(10));
    let read = strip_read(alignment_start, &read_str, "4M10I10M", &low(10..=19, RAW_DUPLEX_MISMATCH_QUAL), &yc, None, None, false);
    assert_eq!("14M", cigar_str(&read));
    assert_eq!("ATGGTCTCTCATGG", read_string(&read));
    check_base_quals(&read, &[4, 5, 6, 7, 8, 9]);

    // test 4: only 3 inserted Ts but all 7 Ts are low qual
    let yc = format!("0+4{}4+", "Z".repeat(7));
    let read = strip_read(alignment_start, "ATGGTTTTTTTATGG", "4M3I8M", &low(4..=10, RAW_DUPLEX_MISMATCH_QUAL), &yc, None, None, false);
    assert_eq!("12M", cigar_str(&read));
    assert_eq!("ATGGTTTTATGG", read_string(&read));
    check_base_quals(&read, &[4, 5, 6, 7]);
}

#[test]
fn test_strip_duplex_multiple_indels() {
    // SbxDuplexReadTest.testStripDuplexMultipleIndels
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let alignment_start = 10;

    // test 1: multiple stripped indels need to offset each other
    let low_quals: Vec<(usize, u8)> = [4, 5, 6, 11, 12, 17].iter().map(|&i| (i, RAW_DUPLEX_MISMATCH_QUAL)).collect();
    let read = strip_read(alignment_start, "AACCTTTAACCTTAACCTAACC", "4M3I4M2I4M1I4M", &low_quals, "0+4ZZZ4ZZ4Z4+", Some(6), None, false);
    assert_eq!("16M", cigar_str(&read));
    assert_eq!("AACCAACCAACCAACC", read_string(&read));
    check_base_quals(&read, &[]);
    assert_eq!(0, int_attribute(&read, b"NM"));

    // test 2: overlapping / adjacent indels
    let low_quals: Vec<(usize, u8)> = [10, 12, 13, 17, 18, 19, 20, 21].iter().map(|&i| (i, RAW_DUPLEX_MISMATCH_QUAL)).collect();
    let read = strip_read(alignment_start, "AAAAAGGGTCTCCCTTCCCTTCGGGG", "9M8I9M", &low_quals, "0+10Z1ZZ3ZZZZZ4+", Some(8), None, false);
    // ReadValidReason.isValidRead
    assert_eq!(read.seq_len(), read.qual().len());
    assert_eq!(read.seq_len() as i32, read.raw_cigar().iter().filter(|c| consumes_read(op(**c))).map(|c| oplen(*c)).sum::<i32>());
    assert!(crate::consensus::has_valid_cigar(read.raw_cigar()));
    assert_eq!("9M1I9M", cigar_str(&read));
    assert_eq!("AAAAAGGGTCCCTTCGGGG", read_string(&read));
    check_base_quals(&read, &[10, 11, 12, 13, 14]);
    assert_eq!(1, int_attribute(&read, b"NM"));
}

#[test]
fn test_strip_duplex_indels_forward_read() {
    // SbxDuplexReadTest.testStripDuplexIndelsForwardRead
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let read_str = format!("{}{}{}", "A".repeat(5), "CT".repeat(2), "A".repeat(5));
    let read = strip_read(25, &read_str, "5M2I7M", &low(7..=8, SBX_DUPLEX_MISMATCH_QUAL), "0+7ZZ5+", Some(2), Some(0), false);
    assert_eq!("AAAAACTAAAAA", read_string(&read));
    assert_eq!("12M", cigar_str(&read));
    check_base_quals(&read, &[5, 6]);
    assert_eq!(0, int_attribute(&read, b"NM"));
}

#[test]
fn test_strip_duplex_indels_left_align_duplex_mismatch_repeat() {
    // SbxDuplexReadTest.testStripDuplexIndelsLeftAlignDuplexMismatchRepeat
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let read = strip_read(25, "AAAAACTCTAAAAA", "5M2I7M", &low(5..=6, SBX_DUPLEX_MISMATCH_QUAL), "0+7ZZ5+", Some(2), Some(0), true);
    assert_eq!("AAAAACTAAAAA", read_string(&read));
    assert_eq!("12M", cigar_str(&read));
    check_base_quals(&read, &[5, 6]);
    assert_eq!(0, int_attribute(&read, b"NM"));

    let read = strip_read(25, "GGGGGAAAAAGGGGG", "5M2I8M", &low(5..=6, SBX_DUPLEX_MISMATCH_QUAL), "0+8ZZ5+", Some(2), Some(0), true);
    assert_eq!("GGGGGAAAGGGGG", read_string(&read));
    assert_eq!("13M", cigar_str(&read));
    check_base_quals(&read, &[5, 7]);
    assert_eq!(0, int_attribute(&read, b"NM"));
}

#[test]
fn test_strip_duplex_indels_with_supplementary() {
    // SbxDuplexReadTest.testStripDuplexIndelsInWithSupplementary: the supplementary alignment's cigar is used to
    // strip indels in the soft clip
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let alignment_start = 25;
    // TestUtils.REF_BASES.substring(1, 21): ten bases without repeats, then ten of A
    let aligned_bases = format!("GATCGTACTG{}", "A".repeat(10));
    let supp_pos_start = 25;

    let strip = |read_bases: &str, cigar: &str, supp: &str, yc: &str, low_qual_index: usize| -> bam::Record {
        let mut read = create_sam_record_unpaired(TEST_READ_ID, CHR_1, alignment_start, read_bases, cigar, false, false, Some(supp));
        set_string_attribute(&mut read, b"YC", yc);
        quals_mut(&mut read)[low_qual_index] = RAW_DUPLEX_MISMATCH_QUAL;
        strip_duplex_indels(&mut read).unwrap();
        read
    };

    let supp_bases = "ACGTACGTAAGACGTACGTAA";
    let read_bases = format!("{aligned_bases}{supp_bases}");
    let read = strip(&read_bases, "20M21S", &supp_data(CHR_2, supp_pos_start, '+', "20S10M1I10M", 60), "20+10Z13+", 30);
    assert_eq!("20M20S", cigar_str(&read));
    assert_eq!(40, read.seq_len());

    // stripped indels in this read result in corrected supp data
    let aligned_bases2 = format!("{}A{}", &aligned_bases[..10], &aligned_bases[10..]);
    let read_bases = format!("{aligned_bases2}ACGTACGTAAACGTACGTAA");
    let read = strip(&read_bases, "10M1I10M20S", &supp_data(CHR_2, supp_pos_start, '+', "21S20M", 60), "0+10Z120+", 10);
    assert_eq!("20M20S", cigar_str(&read));
    assert_eq!(40, read.seq_len());
    let new_supp_data = crate::unmap::extract_supp_alignment(aux_str(&read, b"SA")).unwrap();
    assert_eq!(b"20S20M", new_supp_data.cigar);

    // reversed supp data
    let read_bases = format!("{aligned_bases}ACGTGACGTACGTAA");
    let read = strip(&read_bases, "20M15S", &supp_data(CHR_2, supp_pos_start, '-', "10M1I4M", 60), "0+24Z10+", 24);
    assert_eq!("20M14S", cigar_str(&read));
    assert_eq!(34, read.seq_len());

    // repeat with soft-clip on the other side
    let read_bases = format!("ACGTACGTAAGACGTACGTAA{aligned_bases}");
    let read = strip(&read_bases, "21S20M", &supp_data(CHR_2, supp_pos_start, '+', "10M1I10M20S", 60), "+10Z30+", 10);
    assert_eq!("20S20M", cigar_str(&read));
    assert_eq!(40, read.seq_len());
}

/// SbxDuplexReadTest.createSbxSamRecord: all bases with the duplex quality, map quality 60
fn create_sbx_sam_record(read_name: &str, alignment_start: i32, read_str: &str, cigar: &str) -> bam::Record {
    let mut read = create_sam_record_unpaired(read_name, CHR_1, alignment_start, read_str, cigar, false, false, None);
    set_base_qualities(&mut read, &build_base_quals(read_str.len(), RAW_DUPLEX_QUAL));
    read.set_mapq(60);
    read
}

fn unclipped_start(rec: &bam::Record) -> i32 {
    unclipped_position(rec, true)
}
fn unclipped_end(rec: &bam::Record) -> i32 {
    unclipped_position(rec, false)
}

/// The assertions shared by the four soft clip tests of SbxDuplexReadTest
#[allow(clippy::too_many_arguments)]
fn check_alignment_only_consensus(
    ref_bases: &str, reads: &[bam::Record], expected_cigar: &str, expected_read_str: &str, alignment_start: i32, unclipped_start_pos: i32,
    alignment_end: i32, unclipped_end_pos: i32,
) {
    let mut consensus_reads = sbx_consensus_reads(CHR_1, ref_bases);
    let info = create_consensus_read(&mut consensus_reads, reads);
    let consensus_read = &info.consensus_read;
    assert_eq!(ConsensusOutcome::AlignmentOnly, info.outcome);
    assert_eq!(expected_cigar, cigar_str(consensus_read));
    assert_eq!(expected_read_str, read_string(consensus_read));
    assert!(consensus_read.qual().iter().all(|q| *q == RAW_DUPLEX_QUAL));
    assert_eq!(alignment_start, start(consensus_read));
    assert_eq!(unclipped_start_pos, unclipped_start(consensus_read));
    assert_eq!(alignment_end, end(consensus_read));
    assert_eq!(unclipped_end_pos, unclipped_end(consensus_read));
}

#[test]
fn test_softclip_left() {
    // SbxDuplexReadTest.testSoftclipLeft
    let reads = vec![
        create_sbx_sam_record(TEST_READ_ID, 3, "TTA", "2S1M"),
        create_sbx_sam_record("READ_02", 3, "TTA", "2S1M"),
        create_sbx_sam_record("READ_03", 2, "GAA", "1S2M"),
    ];
    check_alignment_only_consensus("AAA", &reads, "2S1M", "TTA", 3, 1, 3, 3);
}

#[test]
fn test_softclip_right() {
    // SbxDuplexReadTest.testSoftclipRight
    let reads = vec![
        create_sbx_sam_record(TEST_READ_ID, 1, "ATT", "1M2S"),
        create_sbx_sam_record("READ_02", 1, "ATT", "1M2S"),
        create_sbx_sam_record("READ_03", 1, "AAG", "2M1S"),
    ];
    check_alignment_only_consensus("AAA", &reads, "1M2S", "ATT", 1, 1, 1, 3);
}

#[test]
fn test_create_left_softclip_only_when_alignment_score_improves() {
    // SbxDuplexReadTest.testCreateLeftSoftclipOnlyWhenAlignmentScoreImproves
    let reads = vec![create_sbx_sam_record(TEST_READ_ID, 1, "TTGC", "4M"), create_sbx_sam_record("READ_02", 1, "TTGC", "4M")];
    check_alignment_only_consensus("ATGC", &reads, "4M", "TTGC", 1, 1, 4, 4);
}

#[test]
fn test_create_right_softclip_only_when_alignment_score_improves() {
    // SbxDuplexReadTest.testCreateRightSoftclipOnlyWhenAlignmentScoreImproves
    let reads = vec![create_sbx_sam_record(TEST_READ_ID, 1, "ATGG", "4M"), create_sbx_sam_record("READ_02", 1, "ATGG", "4M")];
    check_alignment_only_consensus("ATGC", &reads, "4M", "ATGG", 1, 1, 4, 4);
}

// ---------------------------------------------------------------------------------------------------------
// SbxConsensusTest

/// The reference of SbxConsensusTest without its first character, which a MockRefGenome with zero-based
/// positions never returns: position 1 is the first A.
//                                     10        20        30        40
//                            123456789012345678901234567890123456789012345678
const CONSENSUS_REF_BASES: &str = "AAAAACCCCCGGGGGTTTTTACGTACGTACGTAACCGGTTAACCGGTT";

fn consensus_test_reads() -> ConsensusReads {
    let mut consensus_reads = sbx_consensus_reads(CHR_1, CONSENSUS_REF_BASES);
    consensus_reads.set_chromosome_length(CONSENSUS_REF_BASES.len() as i32 + 1);
    consensus_reads
}

/// SbxConsensusTest.createSamRecord: a paired read with its mate at position 300 of the same chromosome
fn create_consensus_test_record(read_bases: &str, position: i32, base_qualities: &[u8], cigar: &str) -> bam::Record {
    let mut record = create_sam_record_paired(&next_id(), CHR_1, position, read_bases, cigar, CHR_1, 300, true, TEST_READ_CIGAR);
    set_base_qualities(&mut record, base_qualities);
    record
}

fn create_aligned_record(read_bases: &str, position: i32, base_qualities: &[u8]) -> bam::Record {
    create_consensus_test_record(read_bases, position, base_qualities, &format!("{}M", read_bases.len()))
}

fn simplex_record(read_bases: &str, position: i32, cigar: &str) -> bam::Record {
    create_consensus_test_record(read_bases, position, &build_base_quals(read_bases.len(), RAW_SIMPLEX_QUAL), cigar)
}

fn duplex_record(read_bases: &str, position: i32, cigar: &str) -> bam::Record {
    create_consensus_test_record(read_bases, position, &build_base_quals(read_bases.len(), RAW_DUPLEX_QUAL), cigar)
}

fn set_negative_strand(reads: &mut [bam::Record]) {
    for r in reads.iter_mut() {
        set_flag(r, FLAG_REVERSE, true);
    }
}

#[test]
fn test_sbx_simplex_consensus() {
    // SbxConsensusTest.testSbxSimplexConsensus
    let mut consensus_reads = consensus_test_reads();
    let position = 1;
    let ref_bases = &CONSENSUS_REF_BASES[..10];

    // reads have a single disagreeing base but find > 50% in agreement
    let simplex = build_base_quals(ref_bases.len(), RAW_SIMPLEX_QUAL);
    let read_bases = format!("T{}", &ref_bases[1..]);
    let mut reads = vec![
        create_aligned_record(ref_bases, position, &simplex),
        create_aligned_record(&read_bases, position, &simplex),
        create_aligned_record(&read_bases, position, &simplex),
    ];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::AlignmentOnly, info.outcome);
    assert_eq!(read_bases, read_string(&info.consensus_read));
    assert_eq!("10M", cigar_str(&info.consensus_read));
    assert!(info.consensus_read.qual().iter().all(|q| *q == RAW_SIMPLEX_QUAL));
    assert_eq!(position, start(&info.consensus_read));

    // test again with no bases agreeing
    let read_bases = format!("G{}", &ref_bases[1..]);
    reads[2] = create_aligned_record(&read_bases, position, &simplex);
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ref_bases, read_string(&info.consensus_read));
    assert_eq!(SIMPLEX_NO_CONSENSUS_QUAL, info.consensus_read.qual()[0]);
}

#[test]
fn test_sbx_duplex_consensus() {
    // SbxConsensusTest.testSbxDuplexConsensus
    let mut consensus_reads = consensus_test_reads();
    let position = 1;
    let ref_bases = &CONSENSUS_REF_BASES[..10];
    let read_bases = format!("T{}", &ref_bases[1..]);
    let duplex = build_base_quals(read_bases.len(), RAW_DUPLEX_QUAL);

    // reads all agree
    let mut reads = vec![
        create_aligned_record(&read_bases, position, &duplex),
        create_aligned_record(&read_bases, position, &duplex),
        create_aligned_record(&read_bases, position, &duplex),
    ];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(read_bases, read_string(&info.consensus_read));
    assert_eq!("10M", cigar_str(&info.consensus_read));
    assert!(info.consensus_read.qual().iter().all(|q| *q == RAW_DUPLEX_QUAL));
    assert_eq!(position, start(&info.consensus_read));

    // test again with 2 of 3 agreeing, so a consensus is used
    let read_bases3 = format!("G{}", &ref_bases[1..]);
    reads[2] = create_aligned_record(&read_bases3, position, &duplex);
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(read_bases, read_string(&info.consensus_read));
    assert_eq!(RAW_DUPLEX_QUAL, info.consensus_read.qual()[0]);

    // low qual reads are considered for the required percentage but not consensus itself
    let mut base_quals = duplex.clone();
    base_quals[0] = SBX_DUPLEX_MISMATCH_QUAL;
    reads.push(create_aligned_record(&read_bases, position, &base_quals));
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ref_bases, read_string(&info.consensus_read));
    assert_eq!(DUPLEX_NO_CONSENSUS_QUAL, info.consensus_read.qual()[0]);

    // a simplex read is not enough for consensus again
    reads.push(create_aligned_record(&read_bases, position, &build_base_quals(read_bases.len(), RAW_SIMPLEX_QUAL)));
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ref_bases, read_string(&info.consensus_read));
    assert_eq!(DUPLEX_NO_CONSENSUS_QUAL, info.consensus_read.qual()[0]);

    reads.push(create_aligned_record(&read_bases, position, &duplex));
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(read_bases, read_string(&info.consensus_read));
    assert_eq!(RAW_DUPLEX_QUAL, info.consensus_read.qual()[0]);
}

#[test]
fn test_sbx_simple_indel_consensus() {
    // SbxConsensusTest.testSbxSimpleIndelConsensus
    let mut consensus_reads = consensus_test_reads();
    let position = 1;

    // test 1: reads differ by an inserted T only, favours M over I since all counts are equal
    let reads = vec![duplex_record("AACCTTTTGGTT", position, "4M1I7M"), duplex_record("AACCTTTGGTT", position, "11M")];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!("AACCTTTGGTT", read_string(&info.consensus_read));
    assert_eq!("11M", cigar_str(&info.consensus_read));
    assert_eq!(position, start(&info.consensus_read));

    // test 2: favour ref over insert from 2 duplex reads, ignoring simplex read
    let mut quals2 = build_base_quals(12, RAW_DUPLEX_QUAL);
    quals2[4] = RAW_SIMPLEX_QUAL;
    let reads = vec![
        duplex_record("AACCTTTTGGTT", position, "4M1I7M"),
        create_consensus_test_record("AACCTTTTGGTT", position, &quals2, "4M1I7M"),
        duplex_record("AACCTTTGGTT", position, "11M"),
    ];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!("AACCTTTGGTT", read_string(&info.consensus_read));
    assert_eq!("11M", cigar_str(&info.consensus_read));
    assert_eq!(position, start(&info.consensus_read));

    // test 3: both reads have the same cigar, so reverts to standard qual-checking logic
    let mut quals2 = build_base_quals(12, RAW_DUPLEX_QUAL);
    quals2[4] = RAW_DUPLEX_MISMATCH_QUAL;
    let reads =
        vec![duplex_record("AACCTTTTGGTT", position, "4M1I7M"), create_consensus_test_record("AACCGTTTGGTT", position, &quals2, "4M1I7M")];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMatch, info.outcome);
    // the reference base since there is no agreement
    assert_eq!("AACCATTTGGTT", read_string(&info.consensus_read));
    assert_eq!("4M1I7M", cigar_str(&info.consensus_read));
    assert_eq!(DUPLEX_NO_CONSENSUS_QUAL, info.consensus_read.qual()[4]);

    // test 5
    let reads = vec![duplex_record("AACCTTTTGGTT", position, "4M1D8M"), duplex_record("AACCTTTGGTT", position, "4M2D7M")];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!("AACCTTTTGGTT", read_string(&info.consensus_read));
    assert_eq!("4M1D8M", cigar_str(&info.consensus_read));

    // test 6
    let reads = vec![duplex_record("AACCTTTTGGTT", position, "4M1I7M"), duplex_record("AACCTTTTTGGTT", position, "4M2I7M")];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!("AACCTTTTGGTT", read_string(&info.consensus_read));
    assert_eq!("4M1I7M", cigar_str(&info.consensus_read));
}

#[test]
fn test_sbx_indel_consensus2() {
    // SbxConsensusTest.testSbxIndelConsensus2
    let mut consensus_reads = consensus_test_reads();
    let position = 1;
    let reads = vec![
        simplex_record("AACGGTAACGGT", position, "6M1D6M"),
        simplex_record("AACCGGTAACGGT", position, "3M1I3M1D6M"),
        simplex_record("AACCGGTTAACCGGT", position, "3M1I7M1I3M"),
    ];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!("AACCGGTAACGGT", read_string(&info.consensus_read));
    assert_eq!("3M1I3M1D6M", cigar_str(&info.consensus_read));

    // simplex up to an index, then duplex
    let with_duplex_from = |bases: &str, from: usize, cigar: &str| -> bam::Record {
        let mut quals = build_base_quals(bases.len(), RAW_SIMPLEX_QUAL);
        for q in quals[from..].iter_mut() {
            *q = RAW_DUPLEX_QUAL;
        }
        create_consensus_test_record(bases, position, &quals, cigar)
    };
    let reads = vec![
        with_duplex_from("AACCGGTTAAACCGGTTA", 16, "10M1I7M"),
        with_duplex_from("AACCGGTTAACCCCCCGGTTA", 16, "11M4I6M"),
        with_duplex_from("AACCCGGTTAACGGTTA", 8, "3M1I8M1D5M"),
    ];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!("AACCGGTTAACGGTTA", read_string(&info.consensus_read));
    assert_eq!("11M1D5M", cigar_str(&info.consensus_read));
}

#[test]
fn test_sbx_indel_consensus_with_soft_clips() {
    // SbxConsensusTest.testSbxIndelConsensusWithSoftClips
    let mut consensus_reads = consensus_test_reads();
    let (read_bases1, read_cigar1, read_start1) = ("AACAGGGTTAA", "4S2M1I3M1S", 10);
    let (read_bases2, read_cigar2, read_start2) = ("CCGGTTTACC", "1S5M4S", 9);
    let read1 = simplex_record(read_bases1, read_start1, read_cigar1);
    let read2 = simplex_record(read_bases2, read_start2, read_cigar2);
    let mut reads = vec![read1.clone(), read2.clone()];
    // the reference where both are simplex
    let read_bases_net = "AACCGGTTAACC";
    let read_cigar_net = "3S6M3S";

    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelSoftclip, info.outcome);
    assert_eq!(read_start2, start(&info.consensus_read));
    assert_eq!(read_bases2, read_string(&info.consensus_read));
    assert_eq!(read_cigar2, cigar_str(&info.consensus_read));

    // now with 3+ reads
    reads.push(clone_sam_record(&read1, &next_id()));
    reads.push(clone_sam_record(&read2, &next_id()));
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(read_start2, start(&info.consensus_read));
    assert_eq!(read_bases_net, read_string(&info.consensus_read));
    assert_eq!(read_cigar_net, cigar_str(&info.consensus_read));

    // test 2: in the reverse direction
    set_negative_strand(&mut reads);
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(read_start2, start(&info.consensus_read));
    assert_eq!(read_bases_net, read_string(&info.consensus_read));
    assert_eq!(read_cigar_net, cigar_str(&info.consensus_read));

    // soft-clips vs indels - indels are ignored from ref base advancement
    let (read_bases1, read_cigar1, read_start1) = ("AAACCGGTTAAACC", "2M1I8M1I2M", 6);
    let read1 = simplex_record(read_bases1, read_start1, read_cigar1);
    let mut reads = vec![read1.clone(), clone_sam_record(&read1, &next_id()), simplex_record("TTTTGGTTGGGG", 10, "4S4M4S")];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(read_start1, start(&info.consensus_read));
    assert_eq!(read_bases1, read_string(&info.consensus_read));
    assert_eq!(read_cigar1, cigar_str(&info.consensus_read));

    // and in reverse
    set_negative_strand(&mut reads);
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(read_start1, start(&info.consensus_read));
    assert_eq!(read_bases1, read_string(&info.consensus_read));
    assert_eq!(read_cigar1, cigar_str(&info.consensus_read));
}

#[test]
fn test_sbx_indel_consensus_with_soft_clips_jitter() {
    // SbxConsensusTest.testSbxIndelConsensusWithSoftClipsJitter
    let mut consensus_reads = consensus_test_reads();
    let reads = vec![
        simplex_record("ACGTAACCGGTTAC", 5, "4S6M1D3M1S"),
        simplex_record("GTAACCCGGTTACGT", 5, "2S3M1I3M1D3M3S"),
        simplex_record("CGTAACCCGGGTTA", 4, "2S4M1I6M1S"),
    ];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(5, start(&info.consensus_read));
    assert_eq!("ACGTAACCCGGTTACGT", read_string(&info.consensus_read));
    assert_eq!("4S3M1I3M1D3M3S", cigar_str(&info.consensus_read));

    // on the reverse strand
    let mut reads = vec![
        simplex_record("CAACCCGGTTAAA", 5, "1S4M1I4M3S"),
        simplex_record("ACAACCCGGTTTA", 5, "2S4M1I5M1S"),
        simplex_record("AAAAACCCGGTTT", 4, "2S5M1I4M1S"),
    ];
    set_negative_strand(&mut reads);
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(5, start(&info.consensus_read));
    assert_eq!("3S4M1I4M3S", cigar_str(&info.consensus_read));
    assert_eq!("AACAACCCGGTTTAA", read_string(&info.consensus_read));
}

#[test]
fn test_sbx_indel_consensus_with_soft_clips_jitter2() {
    // SbxConsensusTest.testSbxIndelConsensusWithSoftClipsJitter2
    let mut consensus_reads = consensus_test_reads();
    let mut reads = vec![simplex_record("AACCGGTT", 5, "4S4M"), simplex_record("CCTGGTT", 3, "2M1I4M"), simplex_record("CTGGTT", 4, "1M1I4M")];
    set_negative_strand(&mut reads);
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(3, start(&info.consensus_read));
    assert_eq!("2S2M1I4M", cigar_str(&info.consensus_read));
    assert_eq!("AACCTGGTT", read_string(&info.consensus_read));

    // a delete overlapping soft-clips and prior to aligned read starts
    let reads = vec![
        simplex_record("TAACC", 10, "1S4M"),
        simplex_record("TTAACC", 10, "2S4M"),
        simplex_record("AACC", 10, "4M"),
        simplex_record("GTAACC", 7, "2M1D4M"),
    ];
    let info = create_consensus_read(&mut consensus_reads, &reads);
    assert_eq!(ConsensusOutcome::IndelMismatch, info.outcome);
    assert_eq!(7, start(&info.consensus_read));
    assert_eq!("2M1D4M", cigar_str(&info.consensus_read));
    assert_eq!("GTAACC", read_string(&info.consensus_read));
}

// ---------------------------------------------------------------------------------------------------------
// SbxFinaliseReadTest

/// The reference of SbxFinaliseReadTest, read through a MockRefGenome with one-based positions: position 1 is X.
const FINALISE_REF_BASES: &str = "XAAGAGCTCCAGGGGGTTTTTACGTACGTACGTAACCGGTTAACCGGTT";

fn finalise(read: &mut bam::Record) {
    let mut refs = RefBases::open(&ref_fasta(CHR_1, FINALISE_REF_BASES)).unwrap();
    finalise_read(&mut refs, &names(), read);
}

fn consensus_type(read: &bam::Record) -> String {
    aux_str(read, b"UT").map(|s| String::from_utf8_lossy(s).to_string()).unwrap_or_else(|| "NONE".to_string())
}

/// SbxFinaliseReadTest.setDupluxBaseIndex
fn set_duplex_base_index(record: &mut bam::Record) {
    let first_duplex_base_index = find_max_duplex_base_index(has(record, FLAG_REVERSE), &[record.qual()]);
    if first_duplex_base_index >= 0 {
        set_int_attribute(record, b"YX", first_duplex_base_index as i64);
    }
}

#[test]
fn test_sbx_finalise_reads() {
    // SbxFinaliseReadTest.testSbxFinaliseReads
    let position = 1;
    let ref_bases = &FINALISE_REF_BASES[1..11];
    let read_bases = format!("T{}", &ref_bases[1..]);
    let duplex = build_base_quals(read_bases.len(), RAW_DUPLEX_QUAL);

    let mut read = create_aligned_record(&read_bases, position, &duplex);
    finalise(&mut read);
    assert!(read.qual().iter().all(|q| *q == SBX_DUPLEX_QUAL));
    assert_eq!("NONE", consensus_type(&read));
    assert_eq!(0, find_max_duplex_base_index(false, &[read.qual()]));

    let mut read = create_aligned_record(&read_bases, position, &duplex);
    set_flag(&mut read, FLAG_REVERSE, true);
    set_string_attribute(&mut read, b"UT", "SINGLE");
    set_duplex_base_index(&mut read);
    finalise(&mut read);
    assert_eq!("DUAL", consensus_type(&read));
    assert_eq!(9, int_attribute(&read, b"YX"));

    // mark the transition point mid-way through the read
    let mut base_quals = duplex.clone();
    let last_index = read_bases.len() - 1;
    for q in base_quals[last_index - 2..].iter_mut() {
        *q = RAW_SIMPLEX_QUAL;
    }
    let mut read = create_aligned_record(&read_bases, position, &base_quals);
    set_string_attribute(&mut read, b"UT", "SINGLE");
    set_flag(&mut read, FLAG_REVERSE, true);
    set_duplex_base_index(&mut read);
    finalise(&mut read);
    assert_eq!("DUAL", consensus_type(&read));
    assert_eq!(6, int_attribute(&read, b"YX"));

    // test duplex mismatch bases
    let mut base_quals = duplex.clone();
    base_quals[3] = SBX_DUPLEX_MISMATCH_QUAL;
    base_quals[7] = SBX_DUPLEX_MISMATCH_QUAL;
    base_quals[8] = SBX_DUPLEX_MISMATCH_QUAL;
    //                x   xx
    //             0123456789
    // ref bases:  AAGAGCTCCA
    let read_bases = "AATAGCTCCG";
    let mut read = create_aligned_record(read_bases, position, &base_quals);
    finalise(&mut read);
    let expected = [
        SBX_DUPLEX_ADJACENT_2_3_QUAL,
        SBX_DUPLEX_ADJACENT_2_3_QUAL,
        SBX_DUPLEX_ADJACENT_1_QUAL,
        SBX_DUPLEX_MISMATCH_QUAL,
        SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH,
        SBX_DUPLEX_ADJACENT_2_3_QUAL,
        SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH,
        SBX_DUPLEX_MISMATCH_QUAL,
        SBX_DUPLEX_MISMATCH_QUAL,
        SBX_DUPLEX_ADJACENT_1_QUAL,
    ];
    assert_eq!(&expected[..], read.qual());

    // test a conversion of an X to M
    let mut read = create_aligned_record(read_bases, position, &base_quals);
    set_cigar(&mut read, "2X3I2X1D3X");
    finalise(&mut read);
    assert_eq!("2M3I2M1D3M", cigar_str(&read));
}

#[test]
fn test_replace_x_with_m_cigar() {
    // SbxFinaliseReadTest.testReplaceXWithMCigar
    let position = 1;
    let read_bases = &FINALISE_REF_BASES[1..11];
    let base_quals = build_base_quals(read_bases.len(), RAW_DUPLEX_QUAL);
    for (cigar, expected) in [("10X", "10M"), ("2X3I2X1D3X", "2M3I2M1D3M"), ("2X2M2X2M2X", "10M")] {
        let mut read = create_aligned_record(read_bases, position, &base_quals);
        set_cigar(&mut read, cigar);
        finalise(&mut read);
        assert_eq!(expected, cigar_str(&read));
    }
}

#[test]
fn test_mismatch_repeat_mirroring() {
    // SbxFinaliseReadTest.testMismatchRepeatMirroring. The Java helper markLowQualBases marks every base and the
    // result of checkBases is not asserted, so the Java test only runs the routine; here the bases listed are
    // marked and the expected ones checked.
    //                      10        20
    //            01234567890123456789012345678
    let bases = b"TTTTAAAATTTACACACTTTAAGAAGTTT";
    let expected: [usize; 12] = [4, 7, 11, 12, 15, 16, 20, 21, 22, 23, 24, 25];
    for marked in [[4usize, 11, 12, 20, 21, 22], [7, 15, 16, 23, 24, 25]] {
        let mut base_quals = build_base_quals(bases.len(), SBX_DUPLEX_QUAL);
        for i in marked {
            base_quals[i] = SBX_DUPLEX_MISMATCH_QUAL;
        }
        ensure_mirrored_mismatch_repeat_quals(0, bases.len() as i32, bases, &mut base_quals).unwrap();
        let low_qual: Vec<usize> = (0..bases.len()).filter(|&i| base_quals[i] == SBX_DUPLEX_MISMATCH_QUAL).collect();
        assert_eq!(&expected[..], &low_qual[..]);
    }
    // as in Java: every base marked
    let mut base_quals = build_base_quals(bases.len(), SBX_DUPLEX_MISMATCH_QUAL);
    ensure_mirrored_mismatch_repeat_quals(0, bases.len() as i32, bases, &mut base_quals).unwrap();
    assert!(base_quals.iter().all(|q| *q == SBX_DUPLEX_MISMATCH_QUAL));
}

#[test]
fn test_fix_cigars() {
    // SbxFinaliseReadTest.testFixCigars
    for (cigar, expected) in [("10I20M", "10S20M"), ("20M10I", "20M10S"), ("10I20M10I", "10S20M10S"), ("10S10I20M10I10S", "20S20M20S")] {
        let mut elements = raw_cigar(cigar);
        assert!(has_invalid_cigar(&elements));
        correct_invalid_cigar(&mut elements);
        assert_eq!(expected, cigar_string(&elements));
    }
}

// ---------------------------------------------------------------------------------------------------------
// UnpairedReadDuplicateTest

/// A ReadCache as the tests make it (group size and maximum soft clip of 100) with the SBX duplicate collapser
struct CollapsingReadCache {
    read_cache: ReadCache,
    collapser: SbxDuplicateCollapser,
    names: Vec<String>,
}

impl CollapsingReadCache {
    fn new(max_dup_distance: i32) -> Self {
        CollapsingReadCache {
            read_cache: ReadCache::new(100, 100, false, DEFAULT_POP_DISTANCE_CHECK, DEFAULT_DYNAMIC_READ_COUNT_THRESHOLD, 0),
            collapser: SbxDuplicateCollapser::new(max_dup_distance),
            names: names(),
        }
    }
    fn process_read(&mut self, read: bam::Record) {
        self.read_cache.process_read(read, &self.names);
    }
    /// ReadCache.evictAll, which ends in collectFragments
    fn evict_all(&mut self) -> FragmentCoordReads {
        let fragments = self.read_cache.evict_all().unwrap();
        self.collapser.collapse_groups(fragments, &self.names)
    }
    fn process_unpaired(&mut self, read_start: i32, read_end: i32, is_reversed: bool) {
        self.process_read(create_unpaired_record(CHR_1, read_start, read_end, is_reversed));
    }
}

#[test]
fn test_sbx_duplicate_group_collapsing() {
    // the last block of the Java test passes only when a forward read joins reverse reads (R7)
    crate::compat::with_java(true, sbx_duplicate_group_collapsing_java);
}

#[test]
fn sbx_forward_and_reverse_reads_are_collapsed_separately() {
    // R7 corrected
    let max_dup_distance = 2;
    let mut read_cache = CollapsingReadCache::new(max_dup_distance);
    read_cache.process_unpaired(100, 200, true);
    read_cache.process_unpaired(100, 200, true);
    read_cache.process_unpaired(100, 200, false);
    let fragments = read_cache.evict_all();
    assert_eq!(1, fragments.duplicate_groups.len());
    assert_eq!(2, fragments.duplicate_groups[0].total_read_count());
    assert_eq!(1, fragments.single_reads.len());

    // reverse reads within the distance of each other still collapse, and so do forward reads
    read_cache.process_unpaired(100, 200, true);
    read_cache.process_unpaired(101, 201, true);
    read_cache.process_unpaired(100, 200, false);
    read_cache.process_unpaired(99, 199, false);
    let fragments = read_cache.evict_all();
    assert_eq!(2, fragments.duplicate_groups.len());
    assert_eq!(0, fragments.single_reads.len());
    for group in &fragments.duplicate_groups {
        assert_eq!(2, group.total_read_count());
        assert_eq!(group.reads[0].is_reverse(), group.reads[1].is_reverse());
    }
}

fn sbx_duplicate_group_collapsing_java() {
    // UnpairedReadDuplicateTest.testSbxDuplicateGroupCollapsing
    let max_dup_distance = 2;
    let mut read_cache = CollapsingReadCache::new(max_dup_distance);
    let read_start = 100;
    let read_end = 150;

    // no collapsing
    read_cache.process_unpaired(read_start, read_end, false);
    read_cache.process_unpaired(read_start, read_end + max_dup_distance + 1, false);
    read_cache.process_unpaired(read_start - max_dup_distance - 1, read_end, false);
    read_cache.process_unpaired(read_start, read_end - max_dup_distance - 1, false);
    read_cache.process_unpaired(read_start + max_dup_distance + 1, read_end, false);
    let fragments = read_cache.evict_all();
    assert_eq!(0, fragments.duplicate_groups.len());
    assert_eq!(5, fragments.single_reads.len());
    assert_eq!(5, fragments.total_read_count());

    // no collapsing of forward and reverse reads
    read_cache.process_unpaired(read_start, read_end, false);
    read_cache.process_unpaired(50, 100, true);
    let fragments = read_cache.evict_all();
    assert_eq!(0, fragments.duplicate_groups.len());
    assert_eq!(2, fragments.single_reads.len());
    assert_eq!(2, fragments.total_read_count());

    // simple collapsing
    read_cache.process_unpaired(read_start, read_end, false);
    read_cache.process_unpaired(read_start, read_end + max_dup_distance, false);
    read_cache.process_unpaired(read_start - max_dup_distance / 2, read_end + max_dup_distance / 2, false);
    let fragments = read_cache.evict_all();
    assert_eq!(1, fragments.duplicate_groups.len());
    assert_eq!(0, fragments.single_reads.len());
    assert_eq!(3, fragments.total_read_count());

    // no multi-coord duplicate groups
    read_cache.process_unpaired(read_start, read_end, false);
    read_cache.process_unpaired(read_start, read_end, false);
    let fragments = read_cache.evict_all();
    assert_eq!(1, fragments.duplicate_groups.len());
    assert_eq!(0, fragments.single_reads.len());
    assert_eq!(2, fragments.total_read_count());

    // no chain collapsing
    read_cache.process_unpaired(read_start, read_end, false);
    read_cache.process_unpaired(read_start, read_end + max_dup_distance, false);
    read_cache.process_unpaired(read_start - max_dup_distance, read_end + max_dup_distance, false);
    read_cache.process_unpaired(read_start - max_dup_distance / 2, read_end + max_dup_distance + max_dup_distance / 2, false);
    let fragments = read_cache.evict_all();
    assert_eq!(1, fragments.duplicate_groups.len());
    assert_eq!(4, fragments.duplicate_groups[0].total_read_count());
    assert_eq!(0, fragments.single_reads.len());
    assert_eq!(4, fragments.total_read_count());

    // a more complex scenario
    read_cache.process_unpaired(read_start, read_end, false);
    read_cache.process_unpaired(read_start, read_end + max_dup_distance, false);
    read_cache.process_unpaired(read_start, read_end - max_dup_distance, false);
    read_cache.process_unpaired(read_start + max_dup_distance, read_end, false);
    read_cache.process_unpaired(read_start - max_dup_distance, read_end, false);
    read_cache.process_unpaired(read_start, 200, true);
    read_cache.process_unpaired(read_start, 200, true);
    read_cache.process_unpaired(read_start, 200, false);
    read_cache.process_unpaired(200, 300, false);
    read_cache.process_unpaired(200 + max_dup_distance / 2, 300 + max_dup_distance / 2, false);
    read_cache.process_unpaired(200 + max_dup_distance / 2, 300 - max_dup_distance / 2, false);
    read_cache.process_unpaired(200 - max_dup_distance / 2, 300 + max_dup_distance / 2, false);
    read_cache.process_unpaired(200 - max_dup_distance / 2, 300 - max_dup_distance / 2, false);
    let fragments = read_cache.evict_all();
    assert_eq!(3, fragments.duplicate_groups.len());
    let range = |group: &DuplicateGroup, f: fn(&bam::Record) -> i32| -> (i32, i32) {
        (group.reads.iter().map(f).min().unwrap(), group.reads.iter().map(f).max().unwrap())
    };
    let group = &fragments.duplicate_groups[0];
    assert_eq!((98, 102), range(group, start));
    assert_eq!((148, 152), range(group, end));
    let group = &fragments.duplicate_groups[2];
    assert_eq!((199, 201), range(group, start));
    assert_eq!((299, 301), range(group, end));
    assert_eq!(0, fragments.single_reads.len());
    assert_eq!(13, fragments.total_read_count());
}

#[test]
fn test_sbx_duplicate_group_collapsing_largest_group_first() {
    // UnpairedReadDuplicateTest.testSbxDuplicateGroupCollapsingLargestGroupFirst
    let mut read_cache = CollapsingReadCache::new(1);
    let read_start = 100;
    let read_end = 250;

    // too far from the largest group and does not do chain merging
    read_cache.process_unpaired(read_start, read_end, false);
    read_cache.process_unpaired(read_start, read_end + 1, false);
    read_cache.process_unpaired(read_start, read_end + 1, false);
    read_cache.process_unpaired(read_start, read_end + 2, false);
    read_cache.process_unpaired(read_start, read_end + 2, false);
    read_cache.process_unpaired(read_start, read_end + 2, false);
    let fragments = read_cache.evict_all();
    assert_eq!(1, fragments.duplicate_groups.len());
    assert_eq!(6, fragments.duplicate_groups[0].total_read_count());
    assert_eq!(0, fragments.single_reads.len());
    assert_eq!(6, fragments.total_read_count());

    // collapsing cannot chain past 3x the jitter distance
    let read1 = create_unpaired_record(CHR_1, read_start, read_end, false);
    let read4 = create_unpaired_record(CHR_1, read_start + 3, read_end, false);
    let read6 = create_unpaired_record(CHR_1, read_start + 5, read_end, false);
    let read6_name = read6.qname().to_vec();
    read_cache.process_read(read1.clone());
    for _ in 0..3 {
        read_cache.process_read(clone_sam_record(&read1, &next_id()));
    }
    read_cache.process_unpaired(read_start + 1, read_end, false);
    read_cache.process_unpaired(read_start + 2, read_end, false);
    read_cache.process_read(read4.clone());
    read_cache.process_read(clone_sam_record(&read4, &next_id()));
    read_cache.process_unpaired(read_start + 4, read_end, false);
    read_cache.process_read(read6);
    let fragments = read_cache.evict_all();
    assert_eq!(1, fragments.duplicate_groups.len());
    assert_eq!(9, fragments.duplicate_groups[0].total_read_count());
    assert_eq!(1, fragments.single_reads.len());
    assert_eq!(&read6_name[..], fragments.single_reads[0].read.qname());
    assert_eq!(10, fragments.total_read_count());
}

// ---------------------------------------------------------------------------------------------------------
// UnmapUnpairedReadsTest

const READ_BASES_LENGTH: usize = 100;
const SOFT_CLIPPED_READ_CIGAR: &str = "21S79M";

/// TestUtils.CHR_LOCATION_MAP
fn chr_regions(chromosome: &str) -> Vec<UnmappingRegion> {
    let region = |start, end, max_depth| UnmappingRegion { start, end, max_depth };
    match chromosome {
        "1" => vec![region(500, 700, 0)],
        "3" => vec![region(500, 700, UNMAP_MIN_HIGH_DEPTH)],
        "4" => vec![
            region(1000, 2000, UNMAP_MIN_HIGH_DEPTH),
            region(3000, 4000, UNMAP_MIN_HIGH_DEPTH),
            region(5000, 6000, UNMAP_MIN_HIGH_DEPTH),
            region(7000, 8000, UNMAP_MIN_HIGH_DEPTH),
        ],
        _ => Vec::new(),
    }
}

/// TestUtils.checkTransformRead with TestUtils.READ_UNMAPPER
fn check_transform_read(read: &mut bam::Record, chromosome: &str) -> bool {
    let mut unmapper = ReadUnmapper::disabled(&header());
    for (tid, n) in names().iter().enumerate() {
        for region in chr_regions(n) {
            unmapper.add_region(tid, region);
        }
    }
    let mut state = UnmapRegionState {
        chromosome: chromosome.to_string(),
        slice_start: 1,
        slice_end: 1000000,
        regions: chr_regions(chromosome),
        last_matched_index: None,
    };
    unmapper.check_transform_read(read, &mut state, &mut UnmapStats::default())
}

fn unmap_test_read(chr: &str, read_start: i32, cigar: &str, is_supplementary: bool, supp_alignment: Option<&str>) -> bam::Record {
    create_sam_record_unpaired(TEST_READ_ID, chr, read_start, &"A".repeat(READ_BASES_LENGTH), cigar, false, is_supplementary, supp_alignment)
}

fn is_unmapped(read: &bam::Record) -> bool {
    has(read, FLAG_UNMAPPED)
}

#[test]
fn test_no_unmap_unpaired() {
    // UnmapUnpairedReadsTest.testNoUnmapUnpaired: read not in unmap regions
    let mut read = unmap_test_read(CHR_1, 100, TEST_READ_CIGAR, false, None);
    assert!(!check_transform_read(&mut read, CHR_1));
    assert!(!is_unmapped(&read));
}

#[test]
fn test_no_unmap_soft_clipped_unpaired() {
    // UnmapUnpairedReadsTest.testNoUnmapSoftClippedUnpaired: soft-clipped but not in an unmap region
    let mut read = unmap_test_read(CHR_1, 100, SOFT_CLIPPED_READ_CIGAR, false, None);
    assert!(!check_transform_read(&mut read, CHR_1));
    assert!(!is_unmapped(&read));
}

#[test]
fn test_unmap_read_soft_clipped_unpaired() {
    // UnmapUnpairedReadsTest.testUnmapReadSoftClippedUnpaired
    let mut read = unmap_test_read(CHR_1, 550, SOFT_CLIPPED_READ_CIGAR, false, None);
    assert!(check_transform_read(&mut read, CHR_1));
    assert!(is_unmapped(&read));
    assert_eq!(NO_POSITION, start(&read));
    assert_eq!(-1, read.tid());
}

#[test]
fn test_supplementary_soft_clipped_unpaired() {
    // UnmapUnpairedReadsTest.testSupplementarySoftClippedUnpaired: the supplementary data is removed since it is in
    // an unmap region and has a long soft clip
    let mut read = unmap_test_read(CHR_1, 100, TEST_READ_CIGAR, false, Some(&supp_data(CHR_1, 550, '+', SOFT_CLIPPED_READ_CIGAR, 60)));
    assert!(has_aux(&read, b"SA"));
    assert!(check_transform_read(&mut read, CHR_1));
    assert!(!is_unmapped(&read));
    assert!(!has_aux(&read, b"SA"));
}

#[test]
fn test_unmap_supplementary_when_primary_is_unmapped_soft_clipped_unpaired() {
    // UnmapUnpairedReadsTest.testUnmapSupplementaryWhenPrimaryIsUnmappedSoftClippedUnpaired: the supplementary data
    // of a supplementary read is its primary read
    let mut read = unmap_test_read(CHR_2, 100, TEST_READ_CIGAR, true, Some(&supp_data(CHR_1, 550, '+', SOFT_CLIPPED_READ_CIGAR, 60)));
    assert!(has_aux(&read, b"SA"));
    assert!(check_transform_read(&mut read, CHR_2));
    assert!(is_unmapped(&read));
}

#[test]
fn test_unpaired_chimeric_not_unmapped() {
    // UnmapUnpairedReadsTest.testUnpairedChimericNotUnmapped
    let mut read =
        create_sam_record_paired(TEST_READ_ID, CHR_1, 550, &"A".repeat(READ_BASES_LENGTH), TEST_READ_CIGAR, CHR_2, 600, true, TEST_READ_CIGAR);
    set_flag(&mut read, FLAG_PAIRED, false);
    assert!(!check_transform_read(&mut read, CHR_1));
    assert!(!is_unmapped(&read));
}

#[test]
fn test_unmap_read_high_depth_unpaired() {
    // UnmapUnpairedReadsTest.testUnmapReadHighDepthUnpaired
    let mut read = unmap_test_read(CHR_3, 600, TEST_READ_CIGAR, false, None);
    assert!(check_transform_read(&mut read, CHR_3));
    assert!(is_unmapped(&read));
    assert_eq!(NO_POSITION, start(&read));
    assert_eq!(-1, read.tid());
}

#[test]
fn test_unmap_supplementary_high_depth_unpaired() {
    // UnmapUnpairedReadsTest.testUnmapSupplementaryHighDepthUnpaired
    let mut read = unmap_test_read(CHR_3, 400, TEST_READ_CIGAR, false, Some(&supp_data(CHR_3, 600, '+', TEST_READ_CIGAR, 60)));
    assert!(has_aux(&read, b"SA"));
    assert!(check_transform_read(&mut read, CHR_3));
    assert!(!is_unmapped(&read));
    assert!(!has_aux(&read, b"SA"));
}

#[test]
fn test_unmap_multiple_supplementary_alignments_unpaired() {
    // UnmapUnpairedReadsTest.testUnmapMultipleSupplementaryAlignmentsUnpaired
    let sa = [
        supp_data(CHR_1, 1000, '+', TEST_READ_CIGAR, 60),
        supp_data(CHR_2, 6000, '+', TEST_READ_CIGAR, 60),
        supp_data(CHR_3, 600, '+', TEST_READ_CIGAR, 60),
    ]
    .join(";");
    let mut read = unmap_test_read(CHR_1, 100, TEST_READ_CIGAR, false, Some(&sa));
    assert!(has_aux(&read, b"SA"));
    assert!(check_transform_read(&mut read, CHR_1));
    assert!(!is_unmapped(&read));
    assert!(!has_aux(&read, b"SA"));
}

#[test]
fn test_unmap_supplementary_when_unmapping_associated_supplementary_unpaired() {
    // UnmapUnpairedReadsTest.testUnmapSupplementaryWhenUnmappingAssociatedSupplementaryUnpaired
    let sa = [supp_data(CHR_3, 100, '+', TEST_READ_CIGAR, 60), supp_data(CHR_3, 500, '+', TEST_READ_CIGAR, 60)].join(";");
    let mut read = unmap_test_read(CHR_1, 100, TEST_READ_CIGAR, true, Some(&sa));
    assert!(has_aux(&read, b"SA"));
    assert!(check_transform_read(&mut read, CHR_3));
    assert!(is_unmapped(&read));
}

#[test]
fn test_unmap_supplementary_when_primary_is_unmapped_unpaired() {
    // UnmapUnpairedReadsTest.testUnmapSupplementaryWhenPrimaryIsUnmappedUnpaired
    let mut read = unmap_test_read(CHR_2, 100, TEST_READ_CIGAR, true, Some(&supp_data(CHR_3, 550, '+', TEST_READ_CIGAR, 60)));
    assert!(has_aux(&read, b"SA"));
    assert!(check_transform_read(&mut read, CHR_2));
    assert!(is_unmapped(&read));
}

// ---------------------------------------------------------------------------------------------------------
// behaviour of htsjdk that the port reproduces, not covered by the Java tests

#[test]
fn test_yc_tag_duplex_indel_indices() {
    let mut indices = Vec::new();
    assert!(matches!(get_duplex_indel_indices(b"332+35H90+", &mut indices), Ok(DuplexIndelIndices::Found)));
    assert!(indices.is_empty());
    assert!(matches!(get_duplex_indel_indices(b"352+45H6WZ35+", &mut indices), Ok(DuplexIndelIndices::Found)));
    assert_eq!(vec![352 + 45 + 1 + 6 + 1], indices);
    assert!(matches!(get_duplex_indel_indices(b"+9D113Y12I120E9+", &mut indices), Ok(DuplexIndelIndices::Found)));
    assert_eq!(vec![9 + 1 + 113 + 1 + 12], indices);
    // no closing delimiter, or no duplex part
    assert!(matches!(get_duplex_indel_indices(b"100", &mut indices), Ok(DuplexIndelIndices::Invalid)));
    assert!(matches!(get_duplex_indel_indices(b"+", &mut indices), Ok(DuplexIndelIndices::Invalid)));
    assert!(matches!(get_duplex_indel_indices(b"100++", &mut indices), Ok(DuplexIndelIndices::Invalid)));
    // a simplex head that is not a number is a Java exception
    assert!(get_duplex_indel_indices(b"1A+35+", &mut indices).is_err());
}

#[test]
fn test_mirrored_qualities_of_undecoded_read() {
    // a read that no Java setter touched: the mirrored repeat qualities exist in memory only, so a duplicate is
    // written with its stored qualities and a finalised read with the changed ones
    let _setting = with_homopolymer_5_prime_low_base_qual(true);
    //             0123456789
    let bases = "ACGTTTTACG";
    let make = || {
        let mut read = create_sam_record_unpaired(TEST_READ_ID, CHR_1, 1, bases, "10M", false, false, None);
        let mut quals = build_base_quals(bases.len(), RAW_DUPLEX_QUAL);
        quals[3] = SBX_DUPLEX_MISMATCH_QUAL;
        set_base_qualities(&mut read, &quals);
        set_string_attribute(&mut read, b"YC", "+10+");
        (read, quals)
    };

    let (mut read, stored) = make();
    prep_process_read(&mut read, false);
    assert!(!is_decoded(&read));
    assert_eq!(SBX_DUPLEX_MISMATCH_QUAL, read.qual()[6]);
    let mut duplicate = read.clone();
    restore_stored_quals(&mut duplicate);
    assert_eq!(&stored[..], duplicate.qual());
    assert!(aux_str(&duplicate, b"YC").is_some());
    assert_eq!(2, aux_spans(aux_bytes(&duplicate)).len());

    // SAMRecord.equals decodes the record: the change is then written
    let mut compared = read.clone();
    eager_decode(&mut compared);
    assert!(is_decoded(&compared));
    restore_stored_quals(&mut compared);
    assert_eq!(SBX_DUPLEX_MISMATCH_QUAL, compared.qual()[6]);

    // a read decoded before (the unmapper removed its supplementary data) keeps the change
    let (mut read, _) = make();
    prep_process_read(&mut read, true);
    assert!(is_decoded(&read));
    restore_stored_quals(&mut read);
    assert_eq!(SBX_DUPLEX_MISMATCH_QUAL, read.qual()[6]);

    // finaliseRead sets the qualities
    let (mut read, _) = make();
    prep_process_read(&mut read, false);
    let mut refs = RefBases::open(&ref_fasta(CHR_1, "ACGTTTTACGACGT")).unwrap();
    finalise_read(&mut refs, &names(), &mut read);
    assert_eq!(2 + 1, aux_spans(aux_bytes(&read)).len());
    assert_eq!(SBX_DUPLEX_MISMATCH_QUAL, read.qual()[3]);
    assert_eq!(SBX_DUPLEX_MISMATCH_QUAL, read.qual()[6]);
    assert_eq!(SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH, read.qual()[4]);
    assert_eq!(0, int_attribute(&read, b"NM"));
}

/// R8: when X elements become M and neighbours are merged, Java removes the first element equal to the merged one,
/// which can be an earlier element of the cigar.
#[test]
fn merging_match_elements_removes_the_merged_one() {
    let ref_bases = "ACGTACGTACGTACGTACGT";
    let run = |java: bool| {
        crate::compat::with_java(java, || {
            let mut refs = RefBases::open(&ref_fasta(CHR_1, ref_bases)).unwrap();
            // 5M over ACGTA, 2D, 3M over TAC, a mismatch, 5M over TACGT: the last 5M equals the first
            let mut read = create_sbx_sam_record(TEST_READ_ID, 1, "ACGTATACTTACGT", "5M2D3M1X5M");
            finalise_read(&mut refs, &names(), &mut read);
            cigar_str(&read)
        })
    };
    assert_eq!("2D9M5M", run(true));
    assert_eq!("5M2D9M", run(false));
}

/// R9: a read with supplementary data and no soft clip. Java replaces the last cigar element by the supplementary
/// cigar; corrected, the read is treated as one without supplementary data.
#[test]
fn supplementary_data_without_soft_clip_leaves_the_cigar() {
    let _setting = with_homopolymer_5_prime_low_base_qual(false);
    let read_bases = format!("GATCGTACTGA{}", "A".repeat(10));
    let strip = |supp: Option<&str>| -> Option<(String, usize)> {
        let mut read = create_sam_record_unpaired(TEST_READ_ID, CHR_1, 25, &read_bases, "10M1I10M", false, false, supp);
        set_string_attribute(&mut read, b"YC", "0+10Z10+");
        quals_mut(&mut read)[10] = RAW_DUPLEX_MISMATCH_QUAL;
        strip_duplex_indels(&mut read).ok().map(|_| (cigar_str(&read), read.seq_len()))
    };
    let plain = strip(None);
    assert_eq!(plain, Some(("20M".to_string(), 20)));
    let supp = supp_data(CHR_2, 25, '+', "21M", 60);
    assert_eq!(strip(Some(&supp)), plain);
    assert_ne!(crate::compat::with_java(true, || strip(Some(&supp))), plain);
}
