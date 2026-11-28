//! Not from Java: finding S57, upstream de50386378's change to `Realignment.realignedReadIndexPosition` (for an
//! insert whose read context's cigar ends in a soft clip, the core end position is moved back by the insert length
//! before it is turned into a read index). Synthetic reads with a known variant index, and a read that needs
//! realignment (the insert and the bases before it soft-clipped, aligned from the base after the insert):
//! - a context whose cigar holds the insert as an I before the clip: its core end position is right, so the
//!   realigned index is the variant's index without the move (the port's default) and the insert length before it
//!   with the move (Java);
//! - a context built from a read whose insert is soft-clipped (no I): the core end position counts the inserted
//!   bases, so the move is needed, and both modes make it;
//! - a context whose cigar holds part of a long insert as an I and clips the rest: the core end position counts
//!   the clipped part; the port moves it by that part, Java by the whole insert.

use super::util::*;
use crate::compat;
use crate::constants::*;
use crate::evidence::{check_realignment, raw_context, realigned_read_index_position, Realigned};
use crate::matcher::Matcher;
use crate::read_context::{create_context, VariantReadContext};
use crate::refseq::RefSequence;
use crate::types::*;

const VAR_POS: i32 = 120;
const READ_START: i32 = 80;
const TRUE_INDEX: i32 = VAR_POS - READ_START;

fn ref_bases() -> String {
    // 300 bases without long repeats near position 120
    format!("{}{}", REF_BASES_200, &REF_BASES_200[100..200])
}

/// The read bases of a read from READ_START carrying the insert, `after` reference bases after it.
fn carrier_bases(refb: &str, ins: &str, after: usize) -> String {
    let v = VAR_POS as usize;
    format!("{}{}{}", &refb[READ_START as usize..=v], ins, &refb[v + 1..v + 1 + after])
}

/// The context from a read with the insert, `k` aligned bases and a 15-base soft clip; and a read with the same
/// bases soft-clipped up to the base after the insert.
fn context_and_clipped_read(ins: &str, k: usize) -> (VariantReadContext, Read) {
    let refb = ref_bases();
    let ref_seq = RefSequence::from_bases(0, refb.as_bytes().to_vec());
    let r = &refb[VAR_POS as usize..VAR_POS as usize + 1];
    let var = SimpleVariant::new(VAR_POS, r.as_bytes(), format!("{r}{ins}").as_bytes());
    let n = ins.len();
    let bases = carrier_bases(&refb, ins, k + 15);
    let cigar = if k > 0 { format!("41M{n}I{k}M15S") } else { format!("41M{n}I15S") };
    let build = build_sam_record(READ_START, &cigar, &bases);
    let ctx = create_context(&var, &build, TRUE_INDEX, &ref_seq, DEFAULT_FLANK_LENGTH).expect("read context");
    let clipped = build_sam_record(VAR_POS + 1, &format!("{}S{}M", 41 + n, k + 15), &bases);
    (ctx, clipped)
}

fn realigned(ctx: &VariantReadContext, read: &Read) -> (i32, Realigned) {
    let (natural, _) = raw_context(&ctx.variant, read);
    let index = realigned_read_index_position(ctx, read);
    let mut matcher = Matcher::new(ctx, true, false);
    (index, check_realignment(ctx, &mut matcher, read, natural, index, None, false))
}

#[test]
fn s57_insert_before_soft_clip() {
    let refb = ref_bases();
    let duplicating = refb[VAR_POS as usize + 1..VAR_POS as usize + 6].to_string();
    for (ins, ks) in [("GACTA", vec![0usize, 1, 2, 4, 8]), (duplicating.as_str(), vec![8, 12])] {
        for k in ks {
            let (ctx, read) = context_and_clipped_read(ins, k);
            assert!(ctx.read_cigar.ends_with('S'), "{ins} k={k}: {}", ctx.read_cigar);
            // the core end position is the reference position of the core end in the context read's alignment
            let core_end_index = ctx.core_index_end - ctx.var_index + TRUE_INDEX;
            assert_eq!(ctx.core_position_end, VAR_POS + (core_end_index - TRUE_INDEX - ins.len() as i32), "{ins} k={k}");

            let (index, kind) = realigned(&ctx, &read);
            assert_eq!(TRUE_INDEX, index, "{ins} k={k}");
            assert!(kind == Realigned::Exact, "{ins} k={k}");

            let (index, kind) = compat::with_java(true, || realigned(&ctx, &read));
            assert_eq!(TRUE_INDEX - ins.len() as i32, index, "{ins} k={k}");
            assert!(kind == Realigned::None, "{ins} k={k}");
        }
    }
    // a context whose cigar does not end in a soft clip: the same in both modes
    let (ctx, read) = context_and_clipped_read("GACTA", 12);
    assert!(!ctx.read_cigar.ends_with('S'));
    assert_eq!(TRUE_INDEX, realigned(&ctx, &read).0);
    assert_eq!(TRUE_INDEX, compat::with_java(true, || realigned(&ctx, &read)).0);
}

/// A context from a read whose 30-base insert is soft-clipped after its first 2 bases (cigar 43M..S, as the long
/// inserts of the public SBX slice have it), and the read aligned from the base after the insert.
#[test]
fn s57_soft_clipped_insert_context() {
    let refb = ref_bases();
    let ref_seq = RefSequence::from_bases(0, refb.as_bytes().to_vec());
    let ins = "GACTATGCATCAGTCAGTACGATGCTAGCA"; // 30 bases, no repeat of the flanks
    let n = ins.len();
    let r = &refb[VAR_POS as usize..VAR_POS as usize + 1];
    let var = SimpleVariant::new(VAR_POS, r.as_bytes(), format!("{r}{ins}").as_bytes());
    let bases = carrier_bases(&refb, ins, 20);
    let build = build_sam_record(READ_START, &format!("43M{}S", n - 2 + 20), &bases);
    let ctx = create_context(&var, &build, TRUE_INDEX, &ref_seq, DEFAULT_FLANK_LENGTH).expect("read context");
    assert!(ctx.read_cigar.ends_with('S') && !ctx.read_cigar.contains('I'), "{}", ctx.read_cigar);
    // the core end position counts the inserted bases: the insert length past the reference position of the core end
    let core_end_index = ctx.core_index_end - ctx.var_index + TRUE_INDEX;
    assert_eq!(ctx.core_position_end, VAR_POS + (core_end_index - TRUE_INDEX - n as i32) + n as i32);

    let read = build_sam_record(VAR_POS + 1, &format!("{}S20M", 41 + n), &bases);
    for java in [false, true] {
        let (index, kind) = compat::with_java(java, || realigned(&ctx, &read));
        assert_eq!(TRUE_INDEX, index, "java {java}");
        assert!(kind == Realigned::Exact, "java {java}");
    }
}

/// A long insert whose read had 2 of its 30 bases as an insert element and the rest soft-clipped (cigar 41M2I..S, as
/// some long inserts of the Illumina slice have it): the core end position is 28 too far, and the corrected move is 28.
#[test]
fn s57_partly_clipped_insert_context() {
    let refb = ref_bases();
    let ref_seq = RefSequence::from_bases(0, refb.as_bytes().to_vec());
    let ins = "GACTATGCATCAGTCAGTACGATGCTAGCA";
    let n = ins.len() as i32;
    let r = &refb[VAR_POS as usize..VAR_POS as usize + 1];
    let var = SimpleVariant::new(VAR_POS, r.as_bytes(), format!("{r}{ins}").as_bytes());
    let bases = carrier_bases(&refb, ins, 20);
    let build = build_sam_record(READ_START, &format!("41M2I{}S", n - 2 + 20), &bases);
    let Some(ctx) = create_context(&var, &build, TRUE_INDEX, &ref_seq, DEFAULT_FLANK_LENGTH) else {
        panic!("no read context");
    };
    assert!(ctx.read_cigar.ends_with('S') && ctx.read_cigar.contains("2I"), "{}", ctx.read_cigar);
    let read = build_sam_record(VAR_POS + 1, &format!("{}S20M", 41 + n), &bases);
    let (index, kind) = realigned(&ctx, &read);
    assert_eq!(TRUE_INDEX, index);
    assert!(kind == Realigned::Exact);
    // Java moves it by all 30 bases: 2 too many
    let (index, _) = compat::with_java(true, || realigned(&ctx, &read));
    assert_eq!(TRUE_INDEX - 2, index);
}
