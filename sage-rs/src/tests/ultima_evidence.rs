//! Not from the Java test suite: the Java Ultima tests stop at the qual models, so these run the Ultima
//! evidence and soft filter path end to end on a few reads, as a check that the pieces fit together
//! (read context with core extension, qual models, read counting, homopolymer and T0 filters).

use super::util::*;
use crate::candidate::Candidate;
use crate::constants::*;
use crate::evidence::ReadContextCounter;
use crate::evidence_stage::set_read_context_ultima_models;
use crate::filters::{SageVariant, SoftFilter, VariantFilters};
use crate::read_context::create_context;
use crate::types::*;
use crate::ultima::*;
use std::sync::Arc;

fn set_up() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    load_bqr_cache(None).unwrap();
}

/// an Ultima read over REF_BASES_200 with base quals 37, tp 0 and the t0 value given throughout
fn ultima_read(alignment_start: i32, read_bases: &str, cigar: &str, t0: u8) -> Read {
    let n = read_bases.len();
    let mut read = build_ultima_read(alignment_start, read_bases, &default_quals(n), &vec![0i8; n], &vec![t0; n]);
    read.cigar = parse_cigar(cigar);
    read.end = alignment_start + read.cigar.iter().filter(|e| e.op.consumes_ref()).map(|e| e.len).sum::<i32>() - 1;
    read
}

/// a counter for the variant built from the read, with its qual models, after `count` copies of the read
fn count_reads(fixture: &Fixture, variant: &SimpleVariant, read: &Read, var_read_index: i32, count: usize) -> ReadContextCounter {
    let qc = fixture.qc();
    let read_context = create_context(variant, read, var_read_index, &fixture.ref_seq, DEFAULT_FLANK_LENGTH).expect("read context");
    assert!(read_context.is_valid());
    assert!(read_context.artefact.is_none());

    let mut counter = ReadContextCounter::new(Arc::new(read_context), Tier::Panel, DEFAULT_MAX_READ_DEPTH, 1, &qc, false);
    set_read_context_ultima_models(&fixture.ref_seq, CHR_1, &mut counter).unwrap();
    assert!(counter.ultima.as_ref().unwrap().qual_models.is_some());

    for i in 0..count {
        let mut read = read.clone();
        read.name = next_read_id();
        // alternate strands so that the strand bias filters stay out of the way
        if i % 2 == 1 {
            read.flags |= FLAG_REVERSE;
        }
        counter.process_read(&read, 1, None, &qc);
    }
    counter
}

fn soft_filters(fixture: &Fixture, counter: ReadContextCounter) -> SageVariant {
    let candidate =
        Candidate { tier: Tier::Panel, ctx: counter.ctx.clone(), min_events: 1, full_match_support: 1, core_match_support: 0 };
    let mut variant = SageVariant::new(candidate, Vec::new(), vec![counter]);
    let hla = |_: i32| false;
    VariantFilters { cfg: &fixture.cfg, chr: CHR_1, hla: &hla }.apply_soft_filters(&mut variant);
    variant
}

#[test]
fn test_ultima_snv_evidence_and_filters() {
    set_up();
    let fixture = Fixture::new();
    // T>A at 40, in GCG TT CG: a contraction of the TT homopolymer on the right and an inserted A
    let position = 40;
    let variant = SimpleVariant::new(position, b"T", b"A");
    let read_bases = format!("{}A{}", &REF_BASES_200[10..40], &REF_BASES_200[41..80]);
    let cigar = format!("{}M", read_bases.len());

    let read = ultima_read(10, &read_bases, &cigar, 37);
    let counter = count_reads(&fixture, &variant, &read, 30, 6);

    assert_eq!(6, counter.counts.full);
    assert_eq!(6, counter.depth());
    assert!(counter.qualities.full > 0);
    let ultima_data = counter.ultima.as_ref().unwrap();
    assert_eq!(UltimaModelType::Snv, ultima_data.qual_models.as_ref().unwrap().original_qual_model.model_type());
    assert!(ultima_data.homopolymer_avg_quals.iter().all(|&q| q == 37.0));
    assert!(ultima_data.t0_avg_quals.iter().all(|&q| q == 37.0));

    let filtered = soft_filters(&fixture, counter);
    assert!(!filtered.filters.contains(&SoftFilter::MinAvgT0Qual));
    assert!(!filtered.filters.contains(&SoftFilter::MinAvgHpQual));
    assert!(!filtered.filters.contains(&SoftFilter::MinAvgBaseQuality));
    assert!(!filtered.filters.contains(&SoftFilter::MinTumorSupport));

    // the same reads with low T0 values fail the T0 filter only
    let read = ultima_read(10, &read_bases, &cigar, 10);
    let counter = count_reads(&fixture, &variant, &read, 30, 6);
    assert_eq!(6, counter.counts.full);
    let filtered = soft_filters(&fixture, counter);
    assert!(filtered.filters.contains(&SoftFilter::MinAvgT0Qual));
    assert!(!filtered.filters.contains(&SoftFilter::MinAvgHpQual));
    set_test_sequencing_type(None);
}

#[test]
fn test_ultima_homopolymer_delete_evidence_and_filters() {
    set_up();
    let fixture = Fixture::new();
    // AT>A at 26: one T fewer in the 4xT homopolymer at 27-30
    let position = 26;
    let variant = SimpleVariant::new(position, b"AT", b"A");
    let read_bases = format!("{}{}", &REF_BASES_200[5..27], &REF_BASES_200[28..80]);
    let cigar = "22M1D52M";

    let read = ultima_read(5, &read_bases, cigar, 37);
    let counter = count_reads(&fixture, &variant, &read, 21, 6);

    assert_eq!(6, counter.counts.full);
    assert!(counter.qualities.full > 0);
    let ultima_data = counter.ultima.as_ref().unwrap();
    let model = ultima_data.qual_models.as_ref().unwrap().original_qual_model;
    assert_eq!(UltimaModelType::HomopolymerAdjustment, model.model_type());
    assert_eq!(1, model.homopolymer_adjustment().unwrap().ref_adjust_count);
    // the core covers whole homopolymers, including the 3 remaining Ts
    assert!(ultima_data.padded_homopolymer_lengths().contains(&3));

    let filtered = soft_filters(&fixture, counter);
    assert!(!filtered.filters.contains(&SoftFilter::MinAvgT0Qual));
    assert!(!filtered.filters.contains(&SoftFilter::MinAvgBaseQuality));

    // low base quals across the alt reads: the homopolymer qual filter does not apply at a VAF of a
    // heterozygous variant or more, and does once ref reads take the VAF below that
    let mut low_qual_read = read.clone();
    low_qual_read.quals.iter_mut().for_each(|q| *q = 5);
    let counter = count_reads(&fixture, &variant, &low_qual_read, 21, 6);
    assert_eq!(6, counter.counts.full);
    let filtered = soft_filters(&fixture, counter);
    assert!(!filtered.filters.contains(&SoftFilter::MinAvgHpQual));

    let mut counter = count_reads(&fixture, &variant, &low_qual_read, 21, 6);
    let ref_read = ultima_read(5, &REF_BASES_200[5..80], "75M", 37);
    let qc = fixture.qc();
    for _ in 0..12 {
        let mut ref_read = ref_read.clone();
        ref_read.name = next_read_id();
        counter.process_read(&ref_read, 0, None, &qc);
    }
    assert_eq!(6, counter.counts.full);
    assert_eq!(12, counter.counts.ref_);
    assert_eq!(18, counter.depth());
    let filtered = soft_filters(&fixture, counter);
    assert!(filtered.filters.contains(&SoftFilter::MinAvgHpQual));
    set_test_sequencing_type(None);
}
