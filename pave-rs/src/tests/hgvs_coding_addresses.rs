//! HgvsCodingAddressesTest: HGVS coding strings at each kind of position of a 6-exon transcript, on both
//! strands

use super::*;
use crate::variant::NO_LOCAL_PHASE_SET;

struct Fx {
    ref_bases: String,
    ref_genome: MockRefGenome,
    trans: &'static TranscriptData,
    trans_rs: &'static TranscriptData,
}

impl Fx {
    fn new() -> Fx {
        // all exons and introns of length 12
        let parts = [
            "CTAGGACACGAG", "CGAGGGCCCAAA", "CGACACGAGTAA", "GGGCCCAAATTT", "ATGAAGGAACCT", "AAAGGGCCCTTT", "AAGATGGAACCT", "AACCGGTTACGT",
            "GTACACAAGCCT", "CACCGGTTACGT", "ATACACAAGCCT",
        ];
        let ref_bases = generate_test_bases(20) + &parts.concat() + &generate_test_bases(20);
        let mut ref_genome = MockRefGenome::default();
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        //            21           33           45           57           69           81           93           105          117          129          141
        //            |   exon1    |   intron1  |   exon2    |   intron2  |   exon3    |   intron4  |   exon4    |   intron5  |   exon5    |   intron6  |   exon6
        //  +         -24                       -12        -1             1                         13                        *1                        *13
        //  -         *24        *13            *12        *1             24         13             12         1              -1                        -13
        let starts = [21, 45, 69, 93, 117, 141];
        let trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &starts, 11, Some(69), Some(104), false, "");
        let trans_rs = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &starts, 11, Some(69), Some(104), false, "");
        Fx { ref_bases, ref_genome, trans, trans_rs }
    }

    fn check_dup(&self, position: i32, length: i32, expected: &str) {
        self.check_variant_impact(expected, &self.create_dup_variant(position, length), true);
    }
    fn check_dup_rs(&self, position: i32, length: i32, expected: &str) {
        self.check_variant_impact(expected, &self.create_dup_variant(position, length), false);
    }
    fn check_snv(&self, position: i32, r: &str, a: &str, expected: &str) {
        self.check_snv_strand(position, r, a, expected, true);
    }
    fn check_snv_rs(&self, position: i32, r: &str, a: &str, expected: &str) {
        self.check_snv_strand(position, r, a, expected, false);
    }
    fn check_snv_strand(&self, position: i32, r: &str, a: &str, expected: &str, forward: bool) {
        assert_eq!(sub(&self.ref_bases, position - 1, position), r);
        self.check_variant_impact(expected, &VariantData::new(CHR_1, position, r, a), forward);
    }
    fn create_dup_variant(&self, position: i32, length: i32) -> VariantData {
        let pos = position - 1;
        let r = sub(&self.ref_bases, pos, pos + 1);
        let a = sub(&self.ref_bases, pos, pos + length + 1);
        let mut var = VariantData::new(CHR_1, position, &r, &a);
        let alt_bases = a[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);
        var
    }
    fn check_ins(&self, position: i32, bases: &str, expected: &str) {
        self.check_variant_impact(expected, &self.create_ins_variant(position, bases), true);
    }
    fn check_ins_rs(&self, position: i32, bases: &str, expected: &str) {
        self.check_variant_impact(expected, &self.create_ins_variant(position, bases), false);
    }
    fn check_del(&self, position: i32, length: i32, expected: &str) {
        self.check_variant_impact(expected, &self.create_del_variant(position, length), true);
    }
    fn check_del_rs(&self, position: i32, length: i32, expected: &str) {
        self.check_variant_impact(expected, &self.create_del_variant(position, length), false);
    }
    fn check_variant_impact(&self, expected: &str, var: &VariantData, forward: bool) {
        let td = if forward { self.trans } else { self.trans_rs };
        let impact = classify(var, td, &self.ref_genome);
        assert_eq!(expected, impact.coding.hgvs);
    }
    fn create_ins_variant(&self, pos: i32, bases: &str) -> VariantData {
        let r = sub(&self.ref_bases, pos, pos + 1);
        let a = cat!(r, bases);
        VariantData::new(CHR_1, pos, &r, &a)
    }
    fn create_del_variant(&self, position: i32, length: i32) -> VariantData {
        let pos = position - 1;
        let r = sub(&self.ref_bases, pos, pos + length + 1);
        let a = sub(&self.ref_bases, pos, pos + 1);
        VariantData::new(CHR_1, position, &r, &a)
    }
}

#[test]
fn second_downstream_exon_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(30,  2, "c.*13_*14dupCT");
//        fx.check_del_rs(30,  2, "c.*13_*14delCT");
    fx.check_ins_rs(31,  "C", "c.*13_*14insG");
    fx.check_snv_rs(31,  "A", "T", "c.*14T>A");

    fx.check_snv_rs(32, "G", "A", "c.*13C>T");
}

#[ignore]
#[test]
fn second_downstream_intron_exon_boundary_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(31,  2, "c.*13-1_*13dupGC");
//        fx.check_del_rs(31,  2, "c.*13-1_*13delGC");
//        fx.check_ins_rs(32,  "AA", "c.*13_*13+1insTT");
}

#[ignore]
#[test]
fn second_downstream_intron_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(34,  4, "c.*13-6_*13-3dupCCCT");
//        fx.check_del_rs(34,  4, "c.*13-6_*13-3delCCCT");
    fx.check_ins_rs(35,  "GG", "c.*13-4_*13-3insCC");
    fx.check_snv_rs(35, "A", "C", "c.*13-3T>G");

    fx.check_snv_rs(44, "A", "C", "c.*12+1T>G");
}

#[test]
fn downstream_exon_intron_boundary_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(43,  2, "c.*12-0_*12+1dupGT"); // c.*12_*12+1dupGT
    fx.check_del_rs(43,  2, "c.*12_*12+1delGT");
    fx.check_ins_rs(44, "CG",  "c.*12-0_*12+1insCG"); // c.*12_*12+1insCG ??
    fx.check_snv_rs(44, "A", "G", "c.*12+1T>C");
}

#[ignore]
#[test]
fn first_downstream_exon_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(45,  2, "c.*10_*11dupTC");
    fx.check_del_rs(45,  2, "c.*10_*11delTC");
    fx.check_ins_rs(45,  "TA", "c.*11_*12insTA");
    fx.check_snv_rs(45, "C", "G", "c.*12G>C");

    fx.check_dup_rs(54,  2, "c.*1_*2dupTT");
    //        fx.check_del_rs(54,  2, "c.*1_*2delTT");
    fx.check_ins_rs(55,  "CC", "c.*1_*2insGG");
    fx.check_snv_rs(55, "A", "C", "c.*2T>G");

    fx.check_dup_rs(54,  2, "c.*1_*2dupTT");
//        fx.check_del_rs(54,  2, "c.*1_*2delTT");
    fx.check_ins_rs(55,  "CC", "c.*1_*2insGG");
    fx.check_snv_rs(55, "A", "C", "c.*2T>G");

    fx.check_snv_rs(56, "A", "C", "c.*1T>G");
}

#[ignore]
#[test]
fn first_downstream_intron_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(64,  2, "c.24+3_24+4dupAT");
//        fx.check_del_rs(64,  2, "c.24+3_24+4delAT");
//        fx.check_ins_rs(60,  "AAA", "c.22_23insTTT");
//        fx.check_snv_rs(60, "T", "C", "c.23A>G");

//        fx.check_snv_rs(68, "T", "C", "c.24+1A>G");
//        fx.check_snv_rs(67, "T", "C", "c.24+2A>G");
}

#[ignore]
#[test]
fn second_coding_exon_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(70,  2, "c.21_22dupTC");
    fx.check_del_rs(70,  2, "c.21_22delTC");
    fx.check_ins_rs(70,  "AAA", "c.22_23insTTT");
    fx.check_snv_rs(70, "T", "C", "c.23A>G");

//        fx.check_dup_rs(68,  2, "c.23_24dupAT");
//        fx.check_del_rs(68,  2, "c.23_24delAT");
    fx.check_ins_rs(69,  "AAA", "c.23_24insTTT");
    fx.check_snv_rs(69, "A", "C", "c.24T>G");
}

#[ignore]
#[test]
fn intron_exon_boundary2_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(79, 2, "c.13-1_13dupTA");
    fx.check_del_rs(79, 2, "c.13-1_13delTA");
    fx.check_ins_rs(79, "CAG",  "c.13_14insCTG"); // in the exon but useful comparison with ins at 80 below
//        fx.check_ins_rs(80, "CAG",  "c.13-1_13insCTG");
    fx.check_snv_rs(80, "T", "C", "c.13A>G");
}

#[test]
fn first_intron_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(89, 3, "c.12+1_12+3dupAAA");
    fx.check_del_rs(89, 3, "c.12+1_12+3delAAA");
    fx.check_ins_rs(89, "AA",  "c.12+3_12+4insTT");
    fx.check_snv_rs(89, "C", "T", "c.12+4G>A");

    fx.check_dup_rs(81, 2, "c.13-3_13-2dupTT");
    fx.check_del_rs(81, 3, "c.13-4_13-2delCTT");
    fx.check_ins_rs(81, "CAG",  "c.13-2_13-1insCTG");
    fx.check_snv_rs(81, "A", "C", "c.13-1T>G");
}

#[ignore]
#[test]
fn exon_boundary_r_s() {
    let fx = Fx::new();
    fx.check_snv_rs(93, "A", "C", "c.12T>G");

    fx.check_del_rs(91, 2, "c.12_12+1delTA");
    fx.check_dup_rs(91, 2, "c.12-0_12+1dupTA"); // c.12_12+1dupTA ?

    fx.check_del_rs(90, 3, "c.12_12+2delTAA");
    fx.check_dup_rs(90, 3, "c.12-0_12+2dupTAA");

    fx.check_del_rs(90, 4, "c.11_12+2delTTAA");
//        fx.check_dup_rs(90, 4, "c.11_12+2dupTTAA");
}

#[test]
fn first_coding_exon_r_s() {
    let fx = Fx::new();
    fx.check_snv_rs(104, "T", "C", "c.1A>G");

    fx.check_dup_rs(102,  2, "c.1_2dupAG");
    fx.check_del_rs(102,  2, "c.1_2delAG");
    fx.check_ins_rs(103,  "C", "c.1_2insG");
    fx.check_snv_rs(103, "C", "G", "c.2G>C");
}

#[ignore]
#[test]
fn first_upstream_intron_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(114,  2, "c.-1+1_-1+2dupAC");
    fx.check_del_rs(114,  2, "c.-1+1_-1+2delAC");
    fx.check_ins_rs(115,  "AA", "c.-1+1_-1+2insTT");
    fx.check_snv_rs(116, "T", "G", "c.-1+1A>C");

//        fx.check_dup_rs(109,  2, "c.-1+6_1-6dupAC");
    fx.check_del_rs(109,  2, "c.-1+6_-1+7delAC"); // c.-1+6_1-6delAC
    fx.check_ins_rs(110,  "AA", "c.-1+6_-1+7insTT");
    fx.check_snv_rs(111, "T", "G", "c.-1+6A>C");

//        fx.check_dup_rs(107,  2, "c.1-5_1-4dupCG");
    fx.check_del_rs(107,  2, "c.1-5_1-4delCG"); // c.-1+6_1-6delAC
    fx.check_ins_rs(108,  "AA", "c.1-5_1-4insTT");
    fx.check_snv_rs(108, "C", "G", "c.1-4G>C");

//        fx.check_dup_rs(105,  2, "c.1-3_1-2dupGT");
    fx.check_del_rs(105,  2, "c.1-3_1-2delGT");
    fx.check_ins_rs(106,  "G", "c.1-3_1-2insC");
    fx.check_snv_rs(103, "C", "G", "c.2G>C");

    fx.check_snv_rs(105, "A", "G", "c.1-1T>C");
}

#[ignore]
#[test]
fn first_upstream_exon_intron_boundary_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(115,  2, "c.-1_-1+1dupCA");
    fx.check_del_rs(115,  2, "c.-1_-1+1delCA");
    fx.check_ins_rs(116,  "CC", "c.-1-0_-1+1insGG"); // c.-1_-1+1insGG
    fx.check_snv_rs(117, "G", "A", "c.-1C>T");
}

#[test]
fn first_upstream_exon_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(120,  2, "c.-6_-5dupGT");
    fx.check_del_rs(120,  2, "c.-6_-5delGT");
    fx.check_ins_rs(120,  "G", "c.-5_-4insC");
    fx.check_snv_rs(120, "C", "A", "c.-4G>T");
}

#[ignore]
#[test]
fn second_upstream_exon_intron_boundary_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(127,  2, "c.-12-1_-12dupGA");
    fx.check_del_rs(127,  2, "c.-12-1_-12delGA");
//        fx.check_ins_rs(128,  "A", "c.-12-1_-12insT");
    fx.check_ins_rs(127,  "A", "c.-12_-11insT");
    fx.check_snv_rs(128, "T", "A", "c.-12A>T");
}

#[ignore]
#[test]
fn second_upstream_intron_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(137,  2, "c.-13+2_-13+3dupCG");
    fx.check_del_rs(137,  2, "c.-13+2_-13+3delCG");
    fx.check_ins_rs(137,  "G", "c.-13+3_-13+4insC");
    fx.check_snv_rs(137, "A", "G", "c.-13+4T>C");

//        fx.check_dup_rs(132,  2, "c.-12-6_-12-5dupCC");
//        fx.check_del_rs(132,  2, "c.-12-6_-12-5delCC");
    fx.check_ins_rs(132,  "AA", "c.-12-5_-12-4insTT");
    fx.check_snv_rs(132, "C", "A", "c.-12-4G>T");

    fx.check_snv_rs(129, "C", "A", "c.-12-1G>T");
}

#[ignore]
#[test]
fn third_upstream_boundary_r_s() {
    let fx = Fx::new();
//        fx.check_dup_rs(139,  2, "c.-13_-13+1dupTA");
    fx.check_del_rs(139,  2, "c.-13_-13+1delTA"); // c.-13_-13+1delTA ?
    fx.check_ins_rs(140,  "G", "c.-13-0_-13+1insC"); // ?
}

#[ignore]
#[test]
fn second_upstream_exon_r_s() {
    let fx = Fx::new();
    fx.check_dup_rs(150,  2, "c.-24_-23dupAG");
//        fx.check_del_rs(150,  2, "c.-24_-23delAG");
    fx.check_ins_rs(150,  "A", "c.-23_-22insT");
    fx.check_snv_rs(150, "C", "G", "c.-22G>C");

    fx.check_dup_rs(141,  2, "c.-15_-14dupTA");
    fx.check_del_rs(141,  2, "c.-15_-14delTA");
    fx.check_ins_rs(141,  "G", "c.-14_-13insC");
    fx.check_snv_rs(141, "A", "G", "c.-13T>C");
}

#[ignore]
#[test]
fn second_downstream_intron_and_exon() {
    let fx = Fx::new();
    // boundary
//        fx.check_dup_at(127, 2, "c.*12_*12+1dupTC");
    fx.check_del(127, 2, "c.*12_*12+1delTC");
    fx.check_ins(128, "GG", "c.*12-0_*12+1insGG"); // c.*12_*12+1insGG

    // in intron6
    fx.check_dup(128, 2, "c.*12+1_*12+2dupCA");
    fx.check_del(128, 2, "c.*12+1_*12+2delCA");
    fx.check_ins(129, "GG",  "c.*12+1_*12+2insGG");
    fx.check_snv(129, "C", "G", "c.*12+1C>G");

    fx.check_dup(133, 2, "c.*12+6_*12+7dupGT");
    fx.check_del(133, 2, "c.*12+6_*12+7delGT");
    fx.check_ins(134, "AA",  "c.*12+6_*12+7insAA");
    fx.check_snv(134, "G", "A", "c.*12+6G>A");

//        fx.check_dup_at(134, 2, "c.*13-6_*13-5dupTT");
//        fx.check_del_at(134, 2, "c.*13-6_*13-5delTT");
    fx.check_ins(135, "AA",  "c.*13-6_*13-5insAA");
    fx.check_snv(135, "T", "A", "c.*13-6T>A");

    // exon6
    fx.check_dup(140, 2, "c.*13_*14dupAT");
//        fx.check_del(140, 2, "c.*13_*14delAT");
    fx.check_ins(141, "CC",  "c.*13_*14insCC");
    fx.check_snv(141, "A", "T", "c.*13A>T");

    fx.check_dup(141, 2, "c.*14_*15dupTA");
    fx.check_del(141, 2, "c.*14_*15delTA");
    fx.check_ins(142, "CC",  "c.*14_*15insCC");
    fx.check_snv(142, "T", "G", "c.*14T>G");

    fx.check_dup(150, 2, "c.*23_*24dupCT");
    fx.check_del(150, 2, "c.*23_*24delCT");
    fx.check_ins(151, "AA",  "c.*23_*24insAA");
    fx.check_snv(151, "C", "G", "c.*23C>G");

    fx.check_snv(152, "T", "G", "c.*24T>G");
}

#[ignore]
#[test]
fn downstream_exon() {
    let fx = Fx::new();
    //        fx.check_dup_at(116, 2, "c.*1_*2dupGT");
    //        fx.check_del_at(116, 2, "c.*1_*2delGT");
    fx.check_ins(117, "CC", "c.*1_*2insCC");
    fx.check_snv(117, "G", "C", "c.*1G>C");

    fx.check_dup(119, 2, "c.*4_*5dupCA");
    fx.check_del(119, 2, "c.*4_*5delCA");
    fx.check_ins(120, "GG", "c.*4_*5insGG");
    fx.check_snv(120, "C", "A", "c.*4C>A");

    fx.check_dup(126, 2, "c.*11_*12dupCT");
    fx.check_del(126, 2, "c.*11_*12delCT");
    fx.check_ins(127, "GG", "c.*11_*12insGG");
    fx.check_snv(127, "C", "A", "c.*11C>A");

    fx.check_dup(127, 1, "c.*12dupT");
    fx.check_del(127, 1, "c.*12delT");
    fx.check_ins(128, "GG", "c.*12-0_*12+1insGG"); // c.*12_*12+1insGG ?
    fx.check_snv(128, "T", "A", "c.*12T>A");
}

#[ignore]
#[test]
fn downstream_intron_exon_boundary() {
    let fx = Fx::new();
    fx.check_dup(115, 2, "c.*1-1_*1-0dupTG");
//        fx.check_del(115, 2, "c.*1-1_*1delTG");
    //        fx.check_ins_at(116, "CC", "c.*1-1_*1insCC");
}

#[ignore]
#[test]
fn downstream_intron() {
    let fx = Fx::new();
    fx.check_dup(104, 2, "c.24+1_24+2dupAA");
    fx.check_del(104, 2, "c.24+1_24+2delAA");
    fx.check_ins(105, "CC", "c.24+1_24+2insCC");
    fx.check_snv(105, "A", "C", "c.24+1A>C");

    fx.check_dup(109, 2, "c.24+6_24+7dupGT"); //c.24+6_c.*1-6dupGT?
    fx.check_del(109, 2, "c.24+6_24+7delGT");
    fx.check_ins(110, "CC", "c.24+6_24+7insCC");
    fx.check_snv(110, "G", "C", "c.24+6G>C");

    //        fx.check_dup_at(109, 2, "c.*1-2_*1-1dupGT");
    fx.check_del(114, 2, "c.*1-2_*1-1delGT");
    fx.check_ins(115, "CC", "c.*1-2_*1-1insCC");
    fx.check_snv(115, "G", "C", "c.*1-2G>C");
}

#[ignore]
#[test]
fn second_exon_downstream_intron_boundary() {
    let fx = Fx::new();
    //        fx.check_dup_at(103, 2, "c.24_24+1dupTA");
    //        fx.check_del_at(103, 2, "c.24_24+1delTA");
    fx.check_ins(104, "GG", "c.24-0_24+1insGG"); // c.24_24+1insGG?
}

#[ignore]
#[test]
fn second_coding_exon() {
    let fx = Fx::new();
    fx.check_dup(102, 2, "c.23_24dupCT");
//                fx.check_del(102, 2, "c.23_24delCT");
    fx.check_ins(103, "A", "c.23_24insA");
    fx.check_snv(103, "C", "G", "c.23C>G");

    fx.check_dup(93, 3, "c.14_16dupAGA");
    fx.check_del(93, 3, "c.14_16delAGA");
    fx.check_ins(94, "C", "c.14_15insC");
    fx.check_snv(94, "A", "C", "c.14A>C");

    //        fx.check_dup_at(92, 3, "c.13_15dupAAG");
    fx.check_del(92, 3, "c.13_15delAAG");
    fx.check_ins(93, "C", "c.13_14insC");
    fx.check_snv(93, "A", "C", "c.13A>C");
}

#[test]
fn first_intron_second_exon_boundary() {
    let fx = Fx::new();
    fx.check_dup(90, 3, "c.13-2_13-0dupTTA"); // c.13-2_13dupTTA ?
    fx.check_del(90, 3, "c.13-2_13delTTA");
}

#[ignore]
#[test]
fn first_intron() {
    let fx = Fx::new();
    fx.check_dup(90, 2, "c.13-2_13-1dupTT");
    fx.check_del(90, 2, "c.13-2_13-1delTT");
    fx.check_ins(91, "AC", "c.13-2_13-1insAC");
    fx.check_snv(92, "T", "A", "c.13-1T>A");

//        fx.check_dup(86, 3, "c.13-6_13-4dupCCC");
//        fx.check_del(86, 3, "c.13-6_13-4delCCC");
    fx.check_ins(87, "TT", "c.13-6_13-5insTT");
    fx.check_snv(87, "C", "T", "c.13-6C>T");

    fx.check_dup(85, 3, "c.12+6_12+8dupGCC");
    fx.check_del(85, 3, "c.12+6_12+8delGCC");
    fx.check_ins(86, "TT", "c.12+6_12+7insTT");
    fx.check_snv(86, "G", "T", "c.12+6G>T");

    fx.check_dup(84, 3, "c.12+5_12+7dupGGC");
    fx.check_del(84, 3, "c.12+5_12+7delGGC");
    fx.check_ins(85, "TT", "c.12+5_12+6insTT");
    fx.check_snv(85, "G", "T", "c.12+5G>T");

    fx.check_dup(80, 2, "c.12+1_12+2dupAA");
    fx.check_del(80, 2, "c.12+1_12+2delAA");
    fx.check_ins(81, "CC", "c.12+1_12+2insCC");
    fx.check_snv(81, "A", "C", "c.12+1A>C");
}

#[ignore]
#[test]
fn first_exon_first_intron_boundary() {
    let fx = Fx::new();
    //        fx.check_dup_at(79, 2, "c.12_12+1dupTT");
    fx.check_del(79, 3, "c.12_12+2delTAA");
    fx.check_ins(80, "CC", "c.12-0_12+1insCC"); // c.12_12+1insCC ??
}

#[ignore]
#[test]
fn first_coding_exon() {
    let fx = Fx::new();
    fx.check_dup(79, 1, "c.12dupT");
    fx.check_del(79, 1, "c.12delT");
    fx.check_snv(80, "T", "A", "c.12T>A");

    fx.check_dup(78, 2, "c.11_12dupCT");
    fx.check_del(78, 2, "c.11_12delCT");
    fx.check_ins(79, "AA", "c.11_12insAA");
    fx.check_snv(79, "C", "T", "c.11C>T");

    fx.check_dup(76, 2, "c.9_10dupAC");
    fx.check_del(76, 2, "c.9_10delAC");
    fx.check_ins(77, "TTT", "c.9_10insTTT");
    fx.check_snv(77, "A", "C", "c.9A>C");

    fx.check_dup(71, 2, "c.4_5dupAA");
    fx.check_del(71, 2, "c.4_5delAA");
    fx.check_ins(72, "TTT", "c.4_5insTTT");
    fx.check_snv(72, "A", "C", "c.4A>C");

    fx.check_dup(69, 2, "c.2_3dupTG");
    fx.check_del(69, 2, "c.2_3delTG");
    fx.check_ins(70, "C", "c.2_3insC");
    fx.check_snv(70, "T", "C", "c.2T>C");

    //        fx.check_dup_at(68, 3, "c.1_3dupATG");
    fx.check_del(68, 3, "c.1_3delATG");
    fx.check_ins(69, "CC", "c.1_2insCC");
    fx.check_snv(69, "A", "C", "c.1A>C");
    fx.check_del(68, 1, "c.1delA");
}

#[ignore]
#[test]
fn first_upstream_intron() {
    let fx = Fx::new();
    //        fx.check_dup_at(66, 2, "c.1-2_1-1dupTT");
    fx.check_del(66, 2, "c.1-2_1-1delTT");
    //        fx.check_dup_at(64, 4, "c.1-4_1-1dupATTT");
    fx.check_del(64, 4, "c.1-4_1-1delATTT");
    //        fx.check_dup_at(63, 3, "c.1-5_1-3dupAAT");
    fx.check_del(63, 3, "c.1-5_1-3delAAT");
    fx.check_del(63, 1, "c.1-5delA");
    fx.check_del(62, 1, "c.-1+7delA");
    fx.check_del(61, 1, "c.-1+6delC");
    fx.check_del(59, 1, "c.-1+4delC");
    fx.check_del(58, 1, "c.-1+3delG");
    fx.check_del(57, 1, "c.-1+2delG");
    fx.check_del(56, 1, "c.-1+1delG");
    //        fx.check_dup_at(62, 3, "c.-1+6_-1+8dupAAA"); // -1+6_1-4 ??
    fx.check_del(62, 3, "c.-1+7_-1+9delAAA");
    //        fx.check_dup_at(59, 3, "c.-1+3_-1+5dupCCC");
    fx.check_del(59, 3, "c.-1+4_-1+6delCCC");
    //        fx.check_dup_at(57, 3, "c.-1+1_-1+3dupGGC");
    fx.check_del(57, 3, "c.-1+2_-1+4delGGC");
}

#[ignore]
#[test]
fn first_upstream_exon() {
    let fx = Fx::new();
    fx.check_del(55, 1, "c.-1delA");
    //        fx.check_dup_at(56, 3, "c.-1+0_-1+2dupGGG");
    fx.check_del(56, 3, "c.-1+1_-1+3delGGG");
    //        fx.check_dup_at(55, 2, "c.-2_-1dupAA");
    fx.check_del(55, 1, "c.-1delA");
    fx.check_del(55, 2, "c.-1_-1+1delAG");
    //        fx.check_dup_at(54, 2, "c.-2_-1dupAA");
    fx.check_del(54, 2, "c.-2_-1delAA");
    //        fx.check_dup_at(46, 2, "c.-11_-10dupAC");
    fx.check_del(46, 2, "c.-10_-9delAC");
    //        fx.check_dup_at(45, 2, "c.-12_-11dupAC");
    fx.check_del(45, 2, "c.-11_-10delGA");
    fx.check_del(45, 1, "c.-11delG");
}

#[test]
fn second_upstream_intron() {
    let fx = Fx::new();
    fx.check_del(43, 1, "c.-12-1delA");
    fx.check_del(40, 4, "c.-12-4_-12-1delCAAA");
    fx.check_del(36, 2, "c.-13+5_-13+6delGG");
    fx.check_del(32, 2, "c.-13+1_-13+2delCG");
}

#[test]
fn second_upstream_exon() {
    let fx = Fx::new();
    fx.check_del(31, 2, "c.-13_-13+1delGC");
    fx.check_del(27, 3, "c.-17_-15delACG");
    fx.check_del(23, 3, "c.-21_-19delGGA");
    fx.check_del(21, 1, "c.-23delT");
}

