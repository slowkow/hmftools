//! RefBaseSequenceTest
use super::*;
use crate::assembly::constants::LOW_BASE_QUAL_THRESHOLD;
use crate::assembly::indel_builder::calc_indel_inferred_unclipped_positions;
use crate::assembly::ref_base_seq_builder::RefBaseSeqBuilder;
use crate::assembly::sequence_diff::SequenceDiffType;
use crate::assembly::types::support_read::{SupP, SupportRead};

const LOW_QUAL_BASE: u8 = LOW_BASE_QUAL_THRESHOLD - 1;

fn set_base_quals(read: &ReadP, read_start: usize, read_end: usize, qual: u8) {
    for q in &mut read.m().bam_record_mut().base_qualities[read_start..=read_end] {
        *q = qual;
    }
}

fn support_reads(reads: &[&ReadP], junction_read_start_distance: i32) -> Vec<SupP> {
    reads.iter().map(|x| SupportRead::new(x, SupportType::Junction, junction_read_start_distance, 0, 0)).collect()
}

fn get_read_mismatch_count(ref_base_seq_builder: &RefBaseSeqBuilder, read: &ReadP) -> usize {
    let read_state = ref_base_seq_builder.reads().iter().find(|x| x.read() == read).expect("read state");
    read_state.mismatches().iter().filter(|x| x.type_ == SequenceDiffType::Base).count()
}

#[test]
fn test_forward_ref_base_sequences() {
    let _env = test_env();
    let ext_bases = &REF_BASES_200[100..140];
    let ext_base_quals = build_default_base_quals(ext_bases.len());
    let junction = new_junction(CHR_1, 100, FORWARD);

    // first a basic assembly with all reads agreeing
    let read_bases = format!("{}{}", &REF_BASES_200[81..101], ext_bases);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 81, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    let read_bases = format!("{}{}", &REF_BASES_200[51..101], ext_bases); // extends but stops at soft-clip
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 71, &read_bases, &make_cigar_string(&read_bases, 20, ext_bases.len()));
    let read_bases = format!("{}{}", &REF_BASES_200[61..101], ext_bases);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 61, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    let support = support_reads(&[&read1, &read2, &read3], 0);
    let assembly = JunctionAssembly::new(junction, ext_bases.as_bytes().to_vec(), ext_base_quals.clone(), &support, Vec::new());
    let ref_base_seq_builder = RefBaseSeqBuilder::new(assembly);
    assert_eq!(&REF_BASES_200[61..101], ref_base_seq_builder.ref_base_sequence());
    assert_eq!(61, ref_base_seq_builder.ref_base_position());
    assert_eq!(40, ref_base_seq_builder.ref_base_length());
    assert_eq!("40M", ref_base_seq_builder.cigar_str());
    assert!(ref_base_seq_builder.reads().iter().all(|x| x.mismatch_count(true) == 0));

    // now with reads with disagreeing bases
    let read_ref_bases = "AAACCCGGGTTTACGTAACCGGTTACGTAA";
    let read_bases = format!("{read_ref_bases}{ext_bases}");
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 71, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    let read_bases = format!("{}A{}T{}{}", &read_ref_bases[0..5], &read_ref_bases[6..16], &read_ref_bases[17..], ext_bases);
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 71, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    let read_bases = format!("{}A{}T{}{}", &read_ref_bases[0..5], &read_ref_bases[6..26], &read_ref_bases[27..], ext_bases);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 71, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    let support = support_reads(&[&read1, &read2, &read3], 0);
    let junc_bases = format!("{}{}", &read_ref_bases[29..], ext_bases);
    let assembly = JunctionAssembly::new(junction, junc_bases.as_bytes().to_vec(), ext_base_quals.clone(), &support, Vec::new());
    let ref_base_seq_builder = RefBaseSeqBuilder::new(assembly);
    let ref_seq_bases = format!("{}A{}", &read_ref_bases[0..5], &read_ref_bases[6..]);
    assert_eq!(ref_seq_bases, ref_base_seq_builder.ref_base_sequence());
    assert_eq!(71, ref_base_seq_builder.ref_base_position());
    assert_eq!(30, ref_base_seq_builder.ref_base_length());
    assert_eq!("30M", ref_base_seq_builder.cigar_str());
    assert_eq!(3, ref_base_seq_builder.reads().iter().filter(|x| x.mismatch_count(true) == 1).count());

    // with a mix of INDELs, all marked as low-qual
    let insert = "TT";
    let ref_bases1 = format!("{}{}{}{}", &read_ref_bases[0..10], &read_ref_bases[15..20], insert, &read_ref_bases[20..30]);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 71, &format!("{ref_bases1}{ext_bases}"), "10M5D5M2I10M40S");
    // has the insert but not the delete
    let ref_bases2 = format!("{}{}{}", &read_ref_bases[0..20], insert, &read_ref_bases[20..30]);
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 71, &format!("{ref_bases2}{ext_bases}"), "20M2I10M40S");
    set_base_quals(&read2, 10, 14, LOW_QUAL_BASE);
    // has the delete but not the insert
    let ref_bases3 = format!("{}{}", &read_ref_bases[0..10], &read_ref_bases[15..30]);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 71, &format!("{ref_bases3}{ext_bases}"), "10M5D15M40S");
    set_base_quals(&read3, 14, 15, LOW_QUAL_BASE);
    let support = support_reads(&[&read1, &read2, &read3], 0);
    let assembly = JunctionAssembly::new(junction, junc_bases.as_bytes().to_vec(), ext_base_quals.clone(), &support, Vec::new());
    let ref_base_seq_builder = RefBaseSeqBuilder::new(assembly);
    assert_eq!("10M5D5M2I10M", ref_base_seq_builder.cigar_str());
    assert_eq!(ref_bases1, ref_base_seq_builder.ref_base_sequence());
    assert_eq!(71, ref_base_seq_builder.ref_base_position());
    assert_eq!(27, ref_base_seq_builder.ref_base_length());
    assert_eq!(0, get_read_mismatch_count(&ref_base_seq_builder, &read1));
    assert_eq!(5, get_read_mismatch_count(&ref_base_seq_builder, &read2));
    assert_eq!(2, get_read_mismatch_count(&ref_base_seq_builder, &read3));
}

#[test]
fn test_reverse_ref_base_sequences() {
    let _env = test_env();
    let ext_bases = &REF_BASES_200[20..40];
    let junction = new_junction(CHR_1, 100, REVERSE);
    let read_bases = format!("{}{}", ext_bases, &REF_BASES_200[100..120]);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    let read_bases = format!("{}{}", ext_bases, &REF_BASES_200[100..150]);
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 20));
    let read_bases = format!("{}{}", ext_bases, &REF_BASES_200[100..140]);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    let support = support_reads(&[&read1, &read2, &read3], 0);
    let junc_bases = format!("{}{}", ext_bases, &REF_BASES_200[100..101]);
    let junc_base_quals = build_default_base_quals(junc_bases.len());
    let assembly = JunctionAssembly::new(junction, junc_bases.as_bytes().to_vec(), junc_base_quals, &support, Vec::new());
    let ref_base_seq_builder = RefBaseSeqBuilder::new(assembly);
    assert_eq!(&REF_BASES_200[100..140], ref_base_seq_builder.ref_base_sequence());
    assert_eq!(139, ref_base_seq_builder.ref_base_position());
    assert_eq!(40, ref_base_seq_builder.ref_base_length());
    assert_eq!("40M", ref_base_seq_builder.cigar_str());
    assert!(ref_base_seq_builder.reads().iter().all(|x| x.mismatch_count(true) == 0));

    // with a mix of INDELs
    let read_ref_bases = "AAACCCGGGTTTACGTAACCGGTTACGTAA";
    let insert = "TT";
    let ref_bases1 = format!("{}{}{}{}", &read_ref_bases[0..10], &read_ref_bases[15..20], insert, &read_ref_bases[20..30]);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 100, &format!("{ext_bases}{ref_bases1}"), "20S10M5D5M2I10M");
    // has the insert but not the delete
    let ref_bases2 = format!("{}{}{}", &read_ref_bases[0..20], insert, &read_ref_bases[20..30]);
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 100, &format!("{ext_bases}{ref_bases2}"), "20S20M2I10M");
    set_base_quals(&read2, 30, 34, LOW_QUAL_BASE);
    // has the delete but not the insert
    let ref_bases3 = format!("{}{}", &read_ref_bases[0..10], &read_ref_bases[15..30]);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 100, &format!("{ext_bases}{ref_bases3}"), "20S10M5D15M");
    set_base_quals(&read3, 34, 35, LOW_QUAL_BASE);
    let support = support_reads(&[&read1, &read2, &read3], 0);
    let junc_bases = format!("{}{}", ext_bases, &read_ref_bases[0..1]);
    let junc_base_quals = build_default_base_quals(junc_bases.len());
    let assembly = JunctionAssembly::new(junction, junc_bases.as_bytes().to_vec(), junc_base_quals, &support, Vec::new());
    let ref_base_seq_builder = RefBaseSeqBuilder::new(assembly);
    assert_eq!(ref_bases1, ref_base_seq_builder.ref_base_sequence());
    assert_eq!(129, ref_base_seq_builder.ref_base_position());
    assert_eq!(27, ref_base_seq_builder.ref_base_length());
    assert_eq!("10M5D5M2I10M", ref_base_seq_builder.cigar_str());
    assert_eq!(0, get_read_mismatch_count(&ref_base_seq_builder, &read1));
    assert_eq!(5, get_read_mismatch_count(&ref_base_seq_builder, &read2));
    assert_eq!(2, get_read_mismatch_count(&ref_base_seq_builder, &read3));
}

#[test]
fn test_supporting_indel_tolerances() {
    let _env = test_env();
    let ext_bases = &REF_BASES_200[100..140];
    let ext_base_quals = build_default_base_quals(ext_bases.len());
    let junction = new_junction(CHR_1, 100, FORWARD);

    // the reads cover the same ref bases but have different cigar representations
    let read_bases = format!("{}{}", &REF_BASES_200[71..101], ext_bases);
    let cigar = make_cigar_string(&read_bases, 0, ext_bases.len());
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 71, &read_bases, &cigar);
    set_base_quals(&read1, 20, 29, LOW_QUAL_BASE);
    let read2 = clone_read(&read1, &READ_ID_GENERATOR.next_id());
    let read_bases2 = format!("{}{}", &REF_BASES_200[61..101], &ext_bases[0..30]);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 61, &read_bases2, "30M10I30M");
    let read4 = clone_read(&read3, &READ_ID_GENERATOR.next_id());
    let read5 = clone_read(&read3, &READ_ID_GENERATOR.next_id());
    for read in [&read3, &read4, &read5] {
        calc_indel_inferred_unclipped_positions(read.m());
    }
    let support = support_reads(&[&read1, &read2, &read3, &read4, &read5], 30);
    let assembly = JunctionAssembly::new(junction, ext_bases.as_bytes().to_vec(), ext_base_quals, &support, Vec::new());
    let ref_base_seq_builder = RefBaseSeqBuilder::new(assembly);
    assert_eq!(&REF_BASES_200[61..101], ref_base_seq_builder.ref_base_sequence());
    assert_eq!(61, ref_base_seq_builder.ref_base_position());
    assert_eq!(40, ref_base_seq_builder.ref_base_length());
    assert_eq!("30M9I1M", ref_base_seq_builder.cigar_str());
    assert!(ref_base_seq_builder.reads().iter().all(|x| x.mismatch_count(true) == 0));
}
