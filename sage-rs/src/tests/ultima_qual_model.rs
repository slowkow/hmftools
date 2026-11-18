//! Ported from the Java SAGE test suite: seqtech/UltimaQualModelTest.java

use super::util::*;
use crate::constants::*;
use crate::types::*;
use crate::ultima::*;

const BUFFER_REF_BASES: &str = "AACCGGTTAACCGGTTAA";

/// The constructor: Ultima sequencing and UltimaUtils.loadBqrCache(null)
fn set_up() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    load_bqr_cache(None).unwrap();
}

fn variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

/// UltimaQualModelTest.getStraddlingBases
fn get_straddling_bases(ref_bases: &str, variant: &SimpleVariant) -> [u8; 2] {
    let ref_bases = ref_bases.as_bytes();
    let ref_var_index = (variant.pos - 1) as usize;
    let mut straddling_bases = [0u8; 2];
    if variant.is_insert() {
        return straddling_bases;
    }
    straddling_bases[0] = ref_bases[ref_var_index - 1];
    straddling_bases[1] = if variant.is_delete() {
        ref_bases[ref_var_index + variant.ref_.len()]
    } else {
        ref_bases[ref_var_index + variant.alt.len()]
    };
    straddling_bases
}

/// mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant)) over a genome of refBases
fn build_model(ref_bases: &str, variant: &SimpleVariant) -> UltimaQualModel {
    let ref_genome = MockRefGenome::new(ref_bases);
    build_qual_model(&ref_genome, variant, get_straddling_bases(ref_bases, variant), None).expect("qual model")
}

// UltimaQualModelTest.testHomopolymerAdjustment
#[test]
fn test_homopolymer_adjustment() {
    set_up();
    //                                         0123456789
    let ref_bases = format!("{BUFFER_REF_BASES}ATTTTCGTCGT{BUFFER_REF_BASES}");

    // delete of 1 base
    let mut var = variant(19, "AT", "A");
    let mut model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerAdjustment, model.model_type());

    let mut read_bases = format!("{BUFFER_REF_BASES}ATTTCGTCGT");
    let mut base_qualities = default_quals(read_bases.len());
    let mut t0_values = default_quals(read_bases.len());
    let mut tp_values = vec![0i8; read_bases.len()];
    tp_values[19] = 1;
    tp_values[21] = 1;

    let mut read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    let mut calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(38, calc_qual);

    tp_values = vec![0i8; read_bases.len()];
    tp_values[20] = 1;
    base_qualities[20] = 25;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(25, calc_qual);

    // delete of 3 bases
    var = variant(19, "ATTT", "A");
    model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerAdjustment, model.model_type());

    read_bases = format!("{BUFFER_REF_BASES}ATCGTCGT");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];
    tp_values[19] = 3;
    base_qualities[19] = 25;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(25, calc_qual);

    // test where the required adjustment value isn't in the TP values (Java has no assertion here)
    read.tp[19] = -1;
    let _ = model.calculate_qual(&read, 18);

    // insert of 2 bases
    var = variant(19, "A", "ATT");
    model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerAdjustment, model.model_type());

    read_bases = format!("{BUFFER_REF_BASES}ATTTTTTCGTCGT");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];
    tp_values[19] = 1;
    tp_values[20] = -1;
    tp_values[21] = -2; // will use these 2
    tp_values[22] = -2;
    base_qualities[21] = 30;
    base_qualities[22] = 30;
    tp_values[23] = 11;
    tp_values[24] = 1;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(27, calc_qual);

    // C>CA in TT C GTC, insert of base which doesn't match the existing context
    var = variant(24, "C", "CA");
    model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerAdjustment, model.model_type());

    read_bases = format!("{BUFFER_REF_BASES}ATTTTCAGTCGT{BUFFER_REF_BASES}");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];
    tp_values[25] = -1;
    base_qualities[25] = 32;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 24);
    assert_eq!(32, calc_qual);
    set_test_sequencing_type(None);
}

// UltimaQualModelTest.testHomopolymerDeletion
#[test]
fn test_homopolymer_deletion() {
    set_up();
    //                                         0123456789
    let ref_bases = format!("{BUFFER_REF_BASES}ATTTTCGACGT{BUFFER_REF_BASES}");

    // the whole HP must be deleted
    let mut var = variant(19, "ATTTT", "A");
    let mut model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerDeletion, model.model_type());

    let mut read_bases = format!("{BUFFER_REF_BASES}ACGTCGT");
    let mut base_qualities = default_quals(read_bases.len());
    let mut t0_values = default_quals(read_bases.len());
    let mut tp_values = vec![0i8; read_bases.len()];
    t0_values[18] = 25;
    t0_values[19] = 30;

    let mut read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    let mut calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(30, calc_qual);

    read.flags |= FLAG_REVERSE;
    calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(30, calc_qual);

    // test out of cycle deletions
    var = variant(25, "GA", "G");
    model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerDeletion, model.model_type());

    read_bases = format!("{BUFFER_REF_BASES}ATTTTCGCGT");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(bqr_cache().out_of_cycle_t0_qual(), calc_qual);

    read.flags |= FLAG_REVERSE;
    calc_qual = model.calculate_qual(&read, 18);
    assert_eq!(bqr_cache().out_of_cycle_t0_qual(), calc_qual);
    set_test_sequencing_type(None);
}

// UltimaQualModelTest.testHomopolymerTransition
#[test]
fn test_homopolymer_transition() {
    set_up();
    //                                         01234567890
    let ref_bases = format!("{BUFFER_REF_BASES}ATTTTAAAAAG{BUFFER_REF_BASES}");

    // a deletion which crosses 2 HPs
    let mut var = variant(21, "TTTAA", "T");
    let mut model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerTransition, model.model_type());

    //                                              0123456
    let mut read_bases = format!("{BUFFER_REF_BASES}ATTAAAG");
    let mut base_qualities = default_quals(read_bases.len());
    let mut t0_values = default_quals(read_bases.len());
    let mut tp_values = vec![0i8; read_bases.len()];
    tp_values[20] = 2;
    tp_values[21] = 2;
    base_qualities[20] = 11;
    base_qualities[21] = 11;

    tp_values[22] = 2;
    tp_values[24] = 2;
    base_qualities[22] = 25;
    base_qualities[24] = 25;

    let mut read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    let mut calc_qual = model.calculate_qual(&read, 21);
    assert_eq!(30, calc_qual);

    // as before but with longer, lopsided transitional delete
    var = variant(20, "TTTTA", "T");
    model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::HomopolymerTransition, model.model_type());

    //                                     0123456
    read_bases = format!("{BUFFER_REF_BASES}ATAAAAG");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];
    tp_values[20] = 3;
    base_qualities[20] = 25;

    tp_values[21] = -1;
    tp_values[22] = 1;
    tp_values[23] = 1;
    tp_values[24] = -1;
    base_qualities[22] = 12;
    base_qualities[23] = 12;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 20);
    assert_eq!(34, calc_qual);
    set_test_sequencing_type(None);
}

// UltimaQualModelTest.testSNVs
#[test]
fn test_snvs() {
    set_up();
    //                                             01  234  56
    let mut ref_bases = format!("{BUFFER_REF_BASES}AG{}AG{BUFFER_REF_BASES}", "ACA");

    // C>T in ACA > ATA, matches neither side's ref so reverts to max qual
    let mut var = variant(22, "C", "T");
    let mut model = build_model(&ref_bases, &var);
    assert_eq!(UltimaModelType::Snv, model.model_type());

    let mut read_bases = format!("{BUFFER_REF_BASES}AGATAAG");
    let mut base_qualities = default_quals(read_bases.len());
    let mut t0_values = default_quals(read_bases.len());
    let mut tp_values = vec![0i8; read_bases.len()];

    let mut read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    let mut calc_qual = model.calculate_qual(&read, 22);
    assert_eq!(bqr_cache().out_of_cycle_t0_qual(), calc_qual);

    // C>T in CCG > CTG, left contraction, right insertion/expansion
    ref_bases = format!("{BUFFER_REF_BASES}AG{}AG{BUFFER_REF_BASES}", "CCG");
    var = variant(22, "C", "T");
    model = build_model(&ref_bases, &var);

    read_bases = format!("{BUFFER_REF_BASES}AGCTGAG");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];

    // the deleted base qual
    tp_values[21] = 1;
    base_qualities[21] = 16;

    tp_values[22] = -1;
    base_qualities[22] = 21;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 22);
    assert_eq!(21, calc_qual);

    // C>T GCT > GTT, left full delete, right insertion/expansion
    ref_bases = format!("{BUFFER_REF_BASES}AG{}AG{BUFFER_REF_BASES}", "GCT");
    model = build_model(&ref_bases, &var);

    //                                     0123456
    read_bases = format!("{BUFFER_REF_BASES}AGGTTAG{BUFFER_REF_BASES}");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];

    t0_values[21] = 10;
    t0_values[22] = 15; // an SNV with any HP deletion, including of a single base, just uses t0

    tp_values[22] = -1;
    base_qualities[22] = 21;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 22);
    assert_eq!(21, calc_qual);

    // C>T in TCA > TTA, left 1-base ins/expansion, right full delete
    ref_bases = format!("{BUFFER_REF_BASES}AG{}AG{BUFFER_REF_BASES}", "TCA");
    model = build_model(&ref_bases, &var);

    //                                     0123456
    read_bases = format!("{BUFFER_REF_BASES}AGTTAAG{BUFFER_REF_BASES}");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];

    t0_values[22] = 25;
    t0_values[23] = 10;

    tp_values[21] = -1;
    tp_values[22] = -1;
    base_qualities[21] = 20;
    base_qualities[22] = 20;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 22);
    assert_eq!(25, calc_qual);

    // C>T in GCC > GTC, left contraction, right insertion/expansion
    ref_bases = format!("{BUFFER_REF_BASES}AG{}AG{BUFFER_REF_BASES}", "GCC");
    model = build_model(&ref_bases, &var);

    read_bases = format!("{BUFFER_REF_BASES}AGGTCAG");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];

    // the delete base qual
    tp_values[23] = 1;
    base_qualities[23] = 16;

    tp_values[22] = -1;
    base_qualities[22] = 21;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 22);
    assert_eq!(21, calc_qual);

    // contraction on the C side, expansion on the T side
    // C>T in TTT TCC CCC > TTT TTC CCC
    var = variant(23, "C", "T");

    //                                     012  345  678
    ref_bases = format!("{BUFFER_REF_BASES}TTT{}CCC{BUFFER_REF_BASES}", "TCC");
    model = build_model(&ref_bases, &var);

    //                                     012345678
    read_bases = format!("{BUFFER_REF_BASES}TTTTTCCCC{BUFFER_REF_BASES}");
    base_qualities = default_quals(read_bases.len());
    t0_values = default_quals(read_bases.len());
    tp_values = vec![0i8; read_bases.len()];

    // the HP insert base quals
    tp_values[19] = -1;
    tp_values[23] = -1;
    base_qualities[19] = 16;
    base_qualities[23] = 16;

    // the HP delete quals
    tp_values[24] = 1;
    tp_values[27] = 1;
    base_qualities[24] = 21;
    base_qualities[27] = 21;

    read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    calc_qual = model.calculate_qual(&read, 23);
    assert_eq!(18, calc_qual);
    set_test_sequencing_type(None);
}

/// the model for varRefBases > varAltBases at position 23, between the left and right ref bases
fn base_shift_model(left_ref_bases: &str, var_ref_bases: &str, var_alt_bases: &str, right_ref_bases: &str) -> UltimaQualModel {
    let ref_bases = format!("{BUFFER_REF_BASES}{left_ref_bases}{var_ref_bases}{right_ref_bases}{BUFFER_REF_BASES}");
    build_model(&ref_bases, &variant(23, var_ref_bases, var_alt_bases))
}

// UltimaQualModelTest.testBaseShifts
#[test]
fn test_base_shifts() {
    set_up();
    // position starting at 19: 9012 345 6789
    // index                    0134 567 8901
    // ref bases                CAAA ACC TTTC
    // read bases               CAAA CCG TTTC

    // test 1: 4xA -> 3xA contraction on left, 1xG ins on right, ie a left shift
    let mut model = base_shift_model("CAAA", "ACC", "CCG", "TTTC");
    assert_eq!(UltimaModelType::BaseShift, model.model_type());

    let mut base_shift_model_ = *model.adjustments().unwrap();
    assert!(base_shift_model_.left_adjust.is_some());
    assert!(base_shift_model_.right_adjust.is_some());
    assert!(base_shift_model_.is_left_shift());

    let read_bases = format!("{BUFFER_REF_BASES}CAAACCGTTTC");
    let base_qualities = default_quals(read_bases.len());
    let t0_values = default_quals(read_bases.len());
    let tp_values = vec![0i8; read_bases.len()];

    let read = build_ultima_read(1, &read_bases, &base_qualities, &tp_values, &t0_values);

    let calc_qual = model.calculate_qual(&read, 22);
    assert_eq!(50, calc_qual);

    // test 2: 1xG ins on left, 4xT -> 3xT contraction on right
    // ref bases:   CAAA CCT TTTC
    // read bases:  CAAA GCC TTTC
    model = base_shift_model("CAAA", "CCT", "GCC", "TTTC");
    assert_eq!(UltimaModelType::BaseShift, model.model_type());

    base_shift_model_ = *model.adjustments().unwrap();
    assert!(base_shift_model_.left_adjust.is_some());
    assert!(base_shift_model_.right_adjust.is_some());
    assert!(base_shift_model_.is_right_shift());

    // test 3
    // ref: ACGT ACC TTTC
    // alt: ACGT CCG TTTC - 1xA del on left, 1xG ins on right, left shift
    model = base_shift_model("ACGT", "ACC", "CCG", "TTTC");
    assert_eq!(UltimaModelType::BaseShift, model.model_type());

    base_shift_model_ = *model.adjustments().unwrap();
    assert!(base_shift_model_.left_deletion.is_some());
    assert!(base_shift_model_.right_adjust.is_some());
    assert!(base_shift_model_.is_left_shift());

    // test 4
    // ref: TAGG AT GGAA
    // alt: TAGG TG GGAA - 1xA del on left, 2xG -> 3xG on right, left shift
    model = base_shift_model("TAGG", "AT", "TG", "GGAA");
    assert_eq!(UltimaModelType::BaseShift, model.model_type());

    base_shift_model_ = *model.adjustments().unwrap();
    assert!(base_shift_model_.left_deletion.is_some());
    assert!(base_shift_model_.right_adjust.is_some());
    assert!(base_shift_model_.is_left_shift());

    // test 5:
    // ref: TAGG AT GGAA
    // alt: TAGG CA GGAA - 1xC ins on left, 1xT del on right
    model = base_shift_model("TAGG", "AT", "CA", "GGAA");
    assert_eq!(UltimaModelType::BaseShift, model.model_type());

    base_shift_model_ = *model.adjustments().unwrap();
    assert!(base_shift_model_.left_adjust.is_some());
    assert!(base_shift_model_.right_deletion.is_some());
    assert!(base_shift_model_.is_right_shift());

    // test 6:
    // ref: TAGG AT GGAA
    // alt: TAGG GA GGAA - 2xG -> 3xG ins on left, 1xT del on right
    model = base_shift_model("TAGG", "AT", "GA", "GGAA");
    assert_eq!(UltimaModelType::BaseShift, model.model_type());

    base_shift_model_ = *model.adjustments().unwrap();
    assert!(base_shift_model_.left_adjust.is_some());
    assert!(base_shift_model_.right_deletion.is_some());
    assert!(base_shift_model_.is_right_shift());

    // test 7: invalid
    model = base_shift_model("TAGG", "AT", "CG", "GGAA");
    assert_eq!(UltimaModelType::Other, model.model_type());

    // test 8: invalid
    model = base_shift_model("TAGG", "ATC", "TAG", "GGAA");
    assert_eq!(UltimaModelType::Other, model.model_type());
    set_test_sequencing_type(None);
}
