//! Port of consensus/ConsensusReadUtilsTest.

use super::util::*;
use crate::consensus::{form_consensus_read_id, select_template_read, ReadParseState, ReadView};
use crate::rec::{OP_D, OP_I, OP_M, OP_S};
use rust_htslib::bam;

fn create(read_id: &str, read_start: i32, read_bases: &str, cigar: &str, mate_cigar: &str) -> bam::Record {
    create_sam_record(read_id, CHR_1, read_start, read_bases, cigar, CHR_1, 5000, false, false, None, true, Some(mate_cigar))
}

/// ConsensusReadUtilsTest.testConsensusReadId
#[test]
fn test_consensus_read_id() {
    let bases = test_read_bases();
    let read_id_fixed = "ABAB:8:SAMPLE:2:222:12345";
    let mut read1 = create(&format!("{read_id_fixed}:READ_01"), 100, &bases, TEST_READ_CIGAR, TEST_READ_CIGAR);

    let mut consensus_read_id = form_consensus_read_id(read1.qname(), None);
    assert_eq!("ABAB:8:SAMPLE:2:222:12345:CNS_READ_01", String::from_utf8(consensus_read_id).unwrap());

    let umi_id = "ACGTG_ATTGC";
    read1 = create(&format!("{read_id_fixed}:READ_01:{umi_id}"), 100, &bases, TEST_READ_CIGAR, TEST_READ_CIGAR);

    consensus_read_id = form_consensus_read_id(read1.qname(), Some(umi_id));
    assert_eq!(format!("ABAB:8:SAMPLE:2:222:12345:READ_01:CNS_{umi_id}"), String::from_utf8(consensus_read_id).unwrap());
}

/// TemplateReads.selectTemplateRead, returning the template's name
fn template_name(reads: &[&bam::Record], coords_read: &bam::Record) -> String {
    let fragment_coords = create_fragment_coords(coords_read);
    let mut reads: Vec<bam::Record> = reads.iter().map(|r| (*r).clone()).collect();
    let index = select_template_read(&mut reads, &fragment_coords);
    read_name(&reads[index])
}

/// ConsensusReadUtilsTest.testTemplateReadSelection
#[test]
fn test_template_read_selection() {
    let ids = &mut ReadIdGenerator::new();
    let ref_bases = ref_bases();
    let pos_start = 11;
    let consensus_bases = &ref_bases[pos_start as usize..21];

    // all agreeing
    let mut read1 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, TEST_READ_CIGAR);
    let mut read2 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, TEST_READ_CIGAR);
    let mut read3 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, TEST_READ_CIGAR);

    assert_eq!(read_name(&read1), template_name(&[&read3, &read2, &read1], &read1));

    // select on lower cigar
    read1 = create(&ids.next_id(), pos_start, consensus_bases, "10S90M", TEST_READ_CIGAR);
    read2 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, TEST_READ_CIGAR);
    read3 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, TEST_READ_CIGAR);

    assert_eq!(read_name(&read2), template_name(&[&read3, &read2, &read1], &read3));

    // select on mate cigar
    read1 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, "90M10S");
    read2 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, "10S90M");
    read3 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, TEST_READ_CIGAR);
    let read4 = create(&ids.next_id(), pos_start, consensus_bases, TEST_READ_CIGAR, TEST_READ_CIGAR);

    assert_eq!(read_name(&read3), template_name(&[&read4, &read3, &read2, &read1], &read4));
}

/// ConsensusReadUtilsTest.testReadParseState
#[test]
fn test_read_parse_state() {
    let mut bases = "AGGCGGA";
    let mut indel_cigar = "1S2M1I2M1S";
    let mut read1 = create(TEST_READ_ID, 100, bases, indel_cigar, TEST_READ_CIGAR);

    {
        let view = ReadView::new(&read1);
        let mut read_state = ReadParseState::new(&view, true);

        assert_eq!(b'A', read_state.base());
        assert_eq!(DEFAULT_QUAL as i32, read_state.base_qual());
        assert_eq!(OP_S, read_state.element_type);
        assert_eq!(1, read_state.element_length);

        read_state.move_next_base();

        assert_eq!(b'G', read_state.base());
        assert_eq!(OP_M, read_state.element_type);
        assert_eq!(2, read_state.element_length);

        read_state.move_next_base();
        read_state.move_next_base();

        assert_eq!(b'C', read_state.base());
        assert_eq!(OP_I, read_state.element_type);
        assert_eq!(1, read_state.element_length);

        read_state.move_next_base();
        read_state.move_next_base();
        assert!(!read_state.exhausted);
        assert_eq!(b'G', read_state.base());
        assert_eq!(OP_M, read_state.element_type);

        read_state.move_next_base();
        assert!(!read_state.exhausted);
        assert_eq!(b'A', read_state.base());
        assert_eq!(OP_S, read_state.element_type);

        read_state.move_next_base();
        assert!(read_state.exhausted);

        // and in reverse
        read_state = ReadParseState::new(&view, false);

        assert_eq!(b'A', read_state.base());
        assert_eq!(DEFAULT_QUAL as i32, read_state.base_qual());
        assert_eq!(OP_S, read_state.element_type);
        assert_eq!(1, read_state.element_length);

        read_state.move_next_base();
        read_state.move_next_base();
        read_state.move_next_base();

        assert_eq!(b'C', read_state.base());
        assert_eq!(OP_I, read_state.element_type);
        assert_eq!(1, read_state.element_length);

        read_state.move_next_base();
        read_state.move_next_base();
        read_state.move_next_base();
        assert!(!read_state.exhausted);
        assert_eq!(b'A', read_state.base());
        assert_eq!(OP_S, read_state.element_type);

        read_state.move_next_base();
        assert!(read_state.exhausted);
    }

    // with a delete
    bases = "ACGT";
    indel_cigar = "2M3D2M";
    read1 = create(TEST_READ_ID, 100, bases, indel_cigar, TEST_READ_CIGAR);

    let view = ReadView::new(&read1);
    let mut read_state = ReadParseState::new(&view, true);

    assert_eq!(b'A', read_state.base());
    assert_eq!(OP_M, read_state.element_type);

    read_state.move_next_base();
    assert_eq!(b'C', read_state.base());
    assert_eq!(OP_M, read_state.element_type);

    read_state.move_next_base();
    assert_eq!(b'C', read_state.base());
    assert_eq!(OP_D, read_state.element_type);

    read_state.move_next_base();
    assert_eq!(b'C', read_state.base());
    assert_eq!(OP_D, read_state.element_type);

    read_state.move_next_base();
    assert_eq!(b'C', read_state.base());
    assert_eq!(OP_D, read_state.element_type);

    read_state.move_next_base();
    assert_eq!(b'G', read_state.base());
    assert_eq!(OP_M, read_state.element_type);

    read_state.move_next_base();
    assert_eq!(b'T', read_state.base());
    assert_eq!(OP_M, read_state.element_type);

    read_state.move_next_base();
    assert!(read_state.exhausted);
}
