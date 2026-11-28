//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/evidence/SplitReadTest.java

use super::misc_evidence::{create_read_context, create_read_counter};
use super::util::*;
use crate::evidence::SplitReadSegment;

// Java: testSplitReadSegments
#[test]
fn test_split_read_segments() {
    let read_bases = &REF_BASES_200[0..50];
    let read_cigar = "20M100N30M";

    // 90-109 = 20M, 110 - 209 = 100D, 210 - 239 = 30M
    let read = build_sam_record(90, read_cigar, read_bases);

    let split_read_segment = SplitReadSegment::form_segment(&read, 100, 20);

    assert!(split_read_segment.is_some());
    let split_read_segment = split_read_segment.unwrap();
    assert_eq!(20, split_read_segment.read_var_index);
    assert_eq!(0, split_read_segment.segment_index_start);
    assert_eq!(19, split_read_segment.segment_index_end);
    assert_eq!(20, split_read_segment.read_bases.len());

    let split_read_segment = SplitReadSegment::form_segment(&read, 115, 20);

    assert!(split_read_segment.is_none());

    let split_read_segment = SplitReadSegment::form_segment(&read, 215, 26);

    assert!(split_read_segment.is_some());
    let split_read_segment = split_read_segment.unwrap();
    assert_eq!(6, split_read_segment.read_var_index);
    assert_eq!(20, split_read_segment.segment_index_start);
    assert_eq!(49, split_read_segment.segment_index_end);
    assert_eq!(30, split_read_segment.read_bases.len());

    // with a second split and soft-clips at the start and end
    let read_bases = &REF_BASES_200[0..70];
    let read_cigar = "5S20M100N20M100N20M5S";

    // 85-89 = 5S, 90-109 = 20M, 110-209 = 100D, 210-229 = 20M, 230-329 = 100N, 330-349 = 20M, 350-354 = 5S
    let read = build_sam_record(90, read_cigar, read_bases);

    let split_read_segment = SplitReadSegment::form_segment(&read, 335, 51);

    assert!(split_read_segment.is_some());
    let split_read_segment = split_read_segment.unwrap();
    assert_eq!(6, split_read_segment.read_var_index);
    assert_eq!(45, split_read_segment.segment_index_start);
    assert_eq!(69, split_read_segment.segment_index_end);
    assert_eq!(25, split_read_segment.read_bases.len());

    // with normal indels as well
    let read_bases = &REF_BASES_200[0..60];
    let read_cigar = "10M100N10M10D10M10I10M100N10M";

    // 90-99 = 10M, 100-199 = 100D, 200-209 = 10M, 210-219 = 10D, 220-229 = 10M, 229 = 10I, 230-239 = 10M, 240-339 = 100N, 340-349 = 10M
    let read = build_sam_record(90, read_cigar, read_bases);

    let split_read_segment = SplitReadSegment::form_segment(&read, 235, 46);

    assert!(split_read_segment.is_some());
    let split_read_segment = split_read_segment.unwrap();
    assert_eq!(36, split_read_segment.read_var_index);
    assert_eq!(10, split_read_segment.segment_index_start);
    assert_eq!(49, split_read_segment.segment_index_end);
    assert_eq!(40, split_read_segment.read_bases.len());
}

// Java: testSplitReads
#[test]
fn test_split_reads() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let position = 133;
    let read_context = create_read_context(position, "A", "T");

    let mut read_context_counter = create_read_counter(&read_context, &fixture);

    // first a read after an N section
    let context_bases = std::str::from_utf8(&read_context.read_bases).unwrap();
    let alt_read_bases = format!("{}{}{}", &REF_BASES_200[0..10], context_bases, &REF_BASES_200[0..10]);
    assert_eq!(45, alt_read_bases.len());

    let read_bases = format!("{}{}{}", &REF_BASES_200[0..10], alt_read_bases, &REF_BASES_200[0..10]);

    // 1-10 = 10M, 11-110 = 100N, 111-155 = 45M, 156-255 = 100N, 256-265 = 10M
    // with variant at 210-219 buffer, 220-229 left flank, 230-234 core, so position at 232, read index 32
    let read_cigar = "10M100N45M100N10M";

    let read = create_sam_record(&next_read_id(), 1, &read_bases, read_cigar);

    read_context_counter.process_read(&read, 1, None, &qc);
    assert_eq!(1, read_context_counter.counts.full);

    // now a read just covering the core
    let read_cigar = "10M122N23M100N10M";

    let read_bases = format!("{}{}{}", &REF_BASES_200[0..10], &alt_read_bases[22..], &REF_BASES_200[0..10]);
    let read = create_sam_record(&next_read_id(), 1, &read_bases, read_cigar);

    read_context_counter.process_read(&read, 1, None, &qc);
    assert_eq!(1, read_context_counter.counts.partial_core);

    // and on the other side
    let read_cigar = "10M100N23M100N10M";

    let read_bases = format!("{}{}{}", &REF_BASES_200[0..10], &alt_read_bases[0..23], &REF_BASES_200[0..10]);
    let read = create_sam_record(&next_read_id(), 1, &read_bases, read_cigar);

    read_context_counter.process_read(&read, 1, None, &qc);
    assert_eq!(2, read_context_counter.counts.partial_core);
}

// Not from Java: S54, the event count of a split read without an NM tag whose blocks reach beyond the region's
// reference bases. Java 5.1 stops there (htsjdk indexes past the region's bases); the port takes the bases of
// such a block from the genome, in both modes, as upstream SAGE does from de50386378 on.
mod s54 {
    use super::super::util::*;
    use crate::candidate::calc_adjusted_num_mutations;
    use crate::refseq::RefSequence;
    use crate::ultima::RefBases;

    /// A synthetic chromosome, 1-based.
    struct Genome(Vec<u8>);

    impl RefBases for Genome {
        fn get_bases(&self, pos_start: i32, pos_end: i32) -> Option<Vec<u8>> {
            if pos_start < 1 || pos_end < pos_start || pos_end as usize > self.0.len() {
                return None;
            }
            Some(self.0[pos_start as usize - 1..pos_end as usize].to_vec())
        }
    }

    fn genome() -> Genome {
        Genome(REF_BASES_200.as_bytes().iter().cycle().take(4000).copied().collect())
    }

    /// The region's reference bases: positions 1001 to 2000.
    fn window(g: &Genome) -> RefSequence {
        RefSequence::from_bases(1001, g.get_bases(1001, 2000).unwrap())
    }

    /// Read bases equal to the reference at `start..=end`, with the bases at the given offsets changed.
    fn bases(g: &Genome, start: i32, end: i32, changed: &[usize]) -> String {
        let mut b = g.get_bases(start, end).unwrap();
        for &i in changed {
            b[i] = if b[i] == b'A' { b'C' } else { b'A' };
        }
        String::from_utf8(b).unwrap()
    }

    #[test]
    fn test_split_read_beyond_ref_sequence() {
        let g = genome();
        let ref_seq = window(&g);
        // 1951-2000 = 50M inside the window, 2001-3450 = 1450N, 3451-3500 = 50M beyond it
        let read_bases = format!("{}{}", bases(&g, 1951, 2000, &[10]), bases(&g, 3451, 3500, &[3, 40]));
        let mut read = build_sam_record(1951, "50M1450N50M", &read_bases);
        read.nm = None;

        assert_eq!(3, calc_adjusted_num_mutations(&read, &ref_seq, Some(&g)));
        // without a genome only the bases inside the window are counted (the port before S54)
        assert_eq!(1, calc_adjusted_num_mutations(&read, &ref_seq, None));

        // a block that starts before the window
        let read_bases = format!("{}{}", bases(&g, 101, 150, &[0, 1]), bases(&g, 1001, 1050, &[20]));
        let mut read = build_sam_record(101, "50M850N50M", &read_bases);
        read.nm = None;
        assert_eq!(3, calc_adjusted_num_mutations(&read, &ref_seq, Some(&g)));
        assert_eq!(1, calc_adjusted_num_mutations(&read, &ref_seq, None));

        // a block that straddles the window end
        let read_bases = bases(&g, 1981, 2010, &[5, 25]);
        let mut read = build_sam_record(1981, "30M", &read_bases);
        read.nm = None;
        assert_eq!(2, calc_adjusted_num_mutations(&read, &ref_seq, Some(&g)));
        assert_eq!(1, calc_adjusted_num_mutations(&read, &ref_seq, None));
    }

    #[test]
    fn test_reads_java_can_count_are_unchanged() {
        let g = genome();
        let ref_seq = window(&g);
        // inside the window the genome is not asked, and indels count once each (NM less length - 1)
        let read_bases = format!("{}{}{}", bases(&g, 1201, 1220, &[7]), "TT", bases(&g, 1224, 1253, &[]));
        let mut read = build_sam_record(1201, "20M3D2I30M", &read_bases);
        read.nm = None;
        assert_eq!(1 + 1 + 1, calc_adjusted_num_mutations(&read, &ref_seq, Some(&g)));
        assert_eq!(1 + 1 + 1, calc_adjusted_num_mutations(&read, &ref_seq, None));

        // an NM tag is used as it is, beyond the window too
        let read_bases = format!("{}{}", bases(&g, 1951, 2000, &[]), bases(&g, 3451, 3500, &[3, 40]));
        let mut read = build_sam_record(1951, "50M1450N50M", &read_bases);
        read.nm = Some(4);
        assert_eq!(4, calc_adjusted_num_mutations(&read, &ref_seq, Some(&g)));
        assert_eq!(4, calc_adjusted_num_mutations(&read, &ref_seq, None));
    }
}
