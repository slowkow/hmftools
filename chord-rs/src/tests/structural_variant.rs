//! StructuralVariantTest

use crate::sv::StructuralVariant;
use hmf_common::sv::StructuralVariantType::*;

fn sv(chr: &str, pos: i32, r: &str, alt: &str) -> StructuralVariant {
    StructuralVariant::new("", chr, pos, r, alt).unwrap()
}

#[test]
fn can_parse_del_and_dup_length() {
    // check both start and end since bracket orientations are reversed for DEL and DUP ends
    let del_start = sv("1", 207981231, "A", "A[1:208014820[");
    let del_end = sv("1", 208014820, "G", "]1:207981231]G");
    assert_eq!(del_start.sv_type, DEL);
    assert_eq!(del_end.sv_type, DEL);
    assert_eq!(del_start.length, Some(33589));
    assert_eq!(del_end.length, Some(33589));

    let dup_start = sv("4", 66212039, "T", "T[4:66211957[");
    let dup_end = sv("4", 66211957, "G", "]4:66212039]G");
    assert_eq!(dup_start.sv_type, DUP);
    assert_eq!(dup_end.sv_type, DUP);
    assert_eq!(dup_start.length, Some(82));
    assert_eq!(dup_end.length, Some(82));
}

#[test]
fn can_parse_inv_and_bnd() {
    let inv_start = sv("3", 24565108, "A", "[3:24566180[A");
    assert_eq!(inv_start.sv_type, INV);
    assert_eq!(inv_start.length, Some(1072));

    let bnd_start = sv("1", 87337011, "T", "T]10:36119127]");
    assert_eq!(bnd_start.sv_type, BND);
    assert_eq!(bnd_start.length, None);
}

#[test]
fn can_parse_unused_sv_types() {
    let ins_start = sv("3", 60204585, "C", "CTACATACATACATACATACATACATACATACA[3:60204586[");
    assert_eq!(ins_start.sv_type, INS);
    assert_eq!(ins_start.length, Some(1));

    let sgl_start = sv("2", 845380, "C", "CAGGTCATGGGACATTTCCCTTGAGTGTTCACAGGTGCACAGAAACACGGGCCATGGGACATTTCCCACGA.");
    assert_eq!(sgl_start.sv_type, SGL);
    assert_eq!(sgl_start.length, None);
}
