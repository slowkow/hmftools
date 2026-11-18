//! Ported from the Java SAGE test suite: seqtech/SbxTest.java

use super::util::*;
use crate::constants::*;
use crate::evidence::ReadContextCounter;
use crate::read_context::create_context;
use crate::types::*;
use std::sync::Arc;

fn get_next_base(base: &str) -> &'static str {
    match base {
        "A" => "C",
        "C" => "G",
        "G" => "T",
        _ => "A",
    }
}

// SbxTest.testUncertainCoreBases
#[test]
fn test_uncertain_core_bases() {
    set_test_sequencing_type(Some(SequencingType::Sbx));
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let position = 40;
    let ref_base = &REF_BASES_200[position..position + 1];
    let alt_base = get_next_base(ref_base);
    let variant = SimpleVariant::new(position as i32, ref_base.as_bytes(), alt_base.as_bytes());
    let read_bases = format!("{}{}{}", &REF_BASES_200[20..position], alt_base, &REF_BASES_200[position + 1..60]);
    let cigar = format!("{}M", read_bases.len());
    let mut read = create_sam_record(&next_read_id(), 20, &read_bases, &cigar);

    let read_context = create_context(&variant, &read, 20, &fixture.ref_seq, DEFAULT_FLANK_LENGTH).expect("read context");
    assert_eq!(12, read_context.var_index);
    assert!(read_context.is_valid());

    let mut read_counter =
        ReadContextCounter::new(Arc::new(read_context), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &qc, true);
    read_counter.process_read(&read, 1, None, &qc);
    assert_eq!(1, read_counter.counts.full);

    // now a read with uncertain bases in the core - these are tracked
    read.quals[19] = BASE_QUAL_MINIMUM;
    read_counter.process_read(&read, 1, None, &qc);
    assert_eq!(2, read_counter.counts.full);
    assert_eq!(0, read_counter.uncertain_core_base_count);

    read.quals[19] = DEFAULT_BASE_QUAL;
    read.quals[20] = BASE_QUAL_MINIMUM;
    read_counter.process_read(&read, 1, None, &qc);
    assert_eq!(2, read_counter.counts.full);
    assert_eq!(1, read_counter.uncertain_core_base_count);
    set_test_sequencing_type(None);
}
