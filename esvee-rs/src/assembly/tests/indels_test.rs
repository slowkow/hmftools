//! IndelsTest
use super::*;
use crate::assembly::assembly_deduper::dedup_proximate_assemblies;
use crate::assembly::common_utils::IndelCoords;
use crate::assembly::constants::{ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO, MIN_VARIANT_LENGTH};
use crate::assembly::indel_builder::{calc_indel_inferred_unclipped_positions, is_weak_indel_based_unlinked_assembly};
use crate::assembly::junction_assembler::JunctionAssembler;
use crate::assembly::types::support_read::SupportRead;

fn create_indel_junction(chromosome: &str, position: i32, orientation: Orientation) -> P<Junction> {
    P::new(Junction::new(chromosome, position, orientation, false, true, false, MIN_VARIANT_LENGTH, None))
}

fn process_junction(junction: P<Junction>, reads: &[ReadP]) -> Vec<AsmP> {
    JunctionAssembler::new(junction, None, None).process_junction(reads)
}

#[test]
fn test_indel_coords() {
    let _env = test_env();
    let read_start = 100;
    // too short
    assert!(IndelCoords::find_indel_coords(read_start, &cigar_elements_from_str("20M2I20M"), 10).is_none());
    // valid delete: 100-119 M, 120-129 D, 130-139 M
    let indel_coords = IndelCoords::find_indel_coords(read_start, &cigar_elements_from_str("20M10D20M"), 10).expect("indel");
    assert_eq!((10, 119, 130), (indel_coords.length, indel_coords.pos_start, indel_coords.pos_end));
    // valid insert: 100-119 M, 10 I @ 119, 120-129 M
    let indel_coords = IndelCoords::find_indel_coords(read_start, &cigar_elements_from_str("20M10I20M"), 10).expect("indel");
    assert_eq!((10, 119, 120), (indel_coords.length, indel_coords.pos_start, indel_coords.pos_end));
    // multiple, chooses the longest
    let indel_coords = IndelCoords::find_indel_coords(read_start, &cigar_elements_from_str("10M10I10M15D10M20I10M25D20M"), 10).expect("indel");
    assert_eq!((25, 154, 180), (indel_coords.length, indel_coords.pos_start, indel_coords.pos_end));
}

#[test]
fn test_implied_edge_indels_to_soft_clip() {
    let _env = test_env();
    let random = ref_bases_random_100();
    // the cigar remains unchanged but the implied new alignments are calculated and the unclipped alignments stored
    // converts both sides of the insert
    let read = create_read(TEST_READ_ID, 100, &random[0..50], "20M10I20M");
    assert!(calc_indel_inferred_unclipped_positions(read.m()));
    assert_eq!(100, read.alignment_start());
    assert_eq!(139, read.alignment_end());
    assert_eq!(100, read.unclipped_start());
    assert_eq!(139, read.unclipped_end());
    assert!(read.has_indel_implied_unclipped_end());
    assert!(read.has_indel_implied_unclipped_start());
    assert_eq!(149, read.indel_implied_unclipped_end());
    assert_eq!(90, read.indel_implied_unclipped_start());
    // no conversion for soft-clipped reads
    let read = create_read(TEST_READ_ID, 100, &random[0..60], "5S20M10I20M5S");
    assert!(!read.has_indel_implied_unclipped_end());
    assert!(!read.has_indel_implied_unclipped_start());
    // no conversion for short indels
    let read = create_read(TEST_READ_ID, 100, &random[0..64], "20M2I20M2D20M");
    assert!(!read.has_indel_implied_unclipped_end());
    assert!(!read.has_indel_implied_unclipped_start());
    // converts both sides of the delete
    let read = create_read(TEST_READ_ID, 100, &random[0..40], "10M20D20M10D10M");
    assert_eq!(169, read.alignment_end());
    assert!(calc_indel_inferred_unclipped_positions(read.m()));
    assert_eq!(139, read.indel_implied_unclipped_end());
    assert_eq!(130, read.indel_implied_unclipped_start());
}

#[test]
fn test_indel_reads_supporting_split_junctions() {
    let _env = test_env();
    let ext_bases = &REF_BASES_400[300..350];
    let pos_junction = new_junction(CHR_1, 100, FORWARD);
    // assembly has 2 high-qual junction reads
    let read_bases = format!("{}{}", &REF_BASES_400[51..101], ext_bases);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 51, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    let read_bases2 = &read_bases[1..];
    // as in Java, the second read keeps the full bases with the cigar of the shorter sequence
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 52, &read_bases, &make_cigar_string(read_bases2, 0, ext_bases.len()));
    // various indel reads that support the junction
    let indel1 = create_read(&READ_ID_GENERATOR.next_id(), 51, &read_bases, "40M10I50M");
    calc_indel_inferred_unclipped_positions(indel1.m());
    let indel2 = create_read(&READ_ID_GENERATOR.next_id(), 51, &read_bases, "40M10D60M");
    calc_indel_inferred_unclipped_positions(indel2.m());
    let assemblies = process_junction(pos_junction, &[read1, read2, indel1, indel2]);
    assert_eq!(1, assemblies.len());
    let assembly = assemblies[0];
    assert_eq!(4, assembly.support_count());
    assert_eq!(0, assembly.mismatch_read_count());
}

#[test]
fn test_long_delete_assemblies() {
    let _env = test_env();
    let pos_junction = create_indel_junction(CHR_1, 50, FORWARD);
    // first a basic assembly with all reads agreeing
    let read_bases = format!("{}{}", &REF_BASES_200[11..51], &REF_BASES_200[100..140]);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 11, &read_bases, "40M49D40M");
    let read_bases = format!("{}{}", &REF_BASES_200[21..51], &REF_BASES_200[100..150]);
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 21, &read_bases, "30M49D50M");
    // other reads will soft-clip at or near the junctions
    let read_bases = format!("{}{}", &REF_BASES_200[11..51], &REF_BASES_200[100..120]);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 11, &read_bases, "40M20S");
    let read4 = create_read(&READ_ID_GENERATOR.next_id(), 11, &read_bases, "45M15S");
    let read5 = create_read(&READ_ID_GENERATOR.next_id(), 11, &read_bases, "35M25S");
    let assemblies = process_junction(pos_junction, &[read1.clone(), read2.clone(), read3, read4, read5]);
    assert_eq!(1, assemblies.len());
    let assembly = assemblies[0];
    assert_eq!(5, assembly.support_count());
    assert_eq!(0, assembly.mismatch_read_count());

    // the other side
    let neg_junction = create_indel_junction(CHR_1, 100, REVERSE);
    let read_bases = format!("{}{}", &REF_BASES_200[31..51], &REF_BASES_200[100..160]);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, "20S60M");
    let read4 = create_read(&READ_ID_GENERATOR.next_id(), 105, &read_bases, "25S55M");
    let assemblies = process_junction(neg_junction, &[read1, read2, read3, read4]);
    assert_eq!(1, assemblies.len());
    assert_eq!(4, assemblies[0].support_count());
}

#[test]
fn test_junction_indel_type_determination() {
    let _env = test_env();
    // a junction is marked as indel in prep but changed to soft-clip based on extension candidate reads
    let pos_junction = create_indel_junction(CHR_1, 50, FORWARD);
    let read_bases = format!("{}{}", &REF_BASES_200[11..51], &REF_BASES_200[100..140]);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 11, &read_bases, "40M49D40M");
    calc_indel_inferred_unclipped_positions(read1.m());
    // other reads will soft-clip at the junctions
    let read_bases = format!("{}{}", &REF_BASES_200[31..51], &REF_BASES_200[100..140]);
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 31, &read_bases, "20M40S");
    let read3 = clone_read(&read2, &READ_ID_GENERATOR.next_id());
    let read_bases = format!("{}{}", &REF_BASES_200[30..51], &REF_BASES_200[100..140]); // diff frag coords
    let read4 = create_read(&READ_ID_GENERATOR.next_id(), 30, &read_bases, "21M40S");
    let assemblies = process_junction(pos_junction, &[read1, read2, read3, read4]);
    assert_eq!(1, assemblies.len());
    let assembly = assemblies[0];
    assert!(!pos_junction.indel_based());
    assert_eq!(4, assembly.support_count());
    assert_eq!(0, assembly.mismatch_read_count());

    // indels being more dominant than soft-clips
    let pos_junction = create_indel_junction(CHR_1, 50, FORWARD);
    let read_bases = format!("{}{}", &REF_BASES_200[11..51], &REF_BASES_200[100..140]);
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 11, &read_bases, "40M49D40M");
    let read2 = clone_read(&read1, &READ_ID_GENERATOR.next_id());
    let read_bases = format!("{}{}", &REF_BASES_200[10..51], &REF_BASES_200[100..140]);
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 10, &read_bases, "41M49D40M");
    for read in [&read1, &read2, &read3] {
        calc_indel_inferred_unclipped_positions(read.m());
    }
    let read_bases = format!("{}{}", &REF_BASES_200[31..51], &REF_BASES_200[100..140]);
    let read4 = create_read(&READ_ID_GENERATOR.next_id(), 31, &read_bases, "20M40S");
    let assemblies = process_junction(pos_junction, &[read1, read2, read3, read4]);
    assert_eq!(1, assemblies.len());
    assert!(pos_junction.indel_based());
    // as in Java, these check the assembly of the first part again
    assert_eq!(4, assembly.support_count());
    assert_eq!(0, assembly.mismatch_read_count());
}

#[test]
fn test_unlinked_indel_based_assembly_filtering() {
    let _env = test_env();
    let junction = new_junction(CHR_1, 100, FORWARD);
    let ref_bases = &REF_BASES_400[1..100];
    let ext_bases = &REF_BASES_400[200..240];
    let split_read = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 1, &format!("{ref_bases}{ext_bases}"), "100M40S", CHR_1, 300, true);
    let assembly = JunctionAssembly::new_test(junction, split_read.get_bases(), split_read.get_base_quality(), ref_bases.len() as i32);
    assembly.m().add_junction_read(&split_read);
    assert!(!is_weak_indel_based_unlinked_assembly(&assembly));

    // 2 longest junction reads are indel reads now
    let indel_read1 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 51, &format!("{ref_bases}{ext_bases}"), "50M50I50M", CHR_1, 300, true);
    let support = SupportRead::new(&indel_read1, SupportType::Indel, 50, indel_read1.bases_length(), 0);
    support.m().set_reference_mismatches(0);
    assembly.m().add_junction_read_typed(&indel_read1, SupportType::Indel, 50);
    // since 50% of the reads now
    assert!(is_weak_indel_based_unlinked_assembly(&assembly));
    for _ in 0..5 {
        assembly.m().add_junction_read(&split_read);
    }
    assert!(!is_weak_indel_based_unlinked_assembly(&assembly));
    // 2 longest reads are indels
    assembly.m().add_junction_read_typed(&indel_read1, SupportType::Indel, 50);
    assert!(is_weak_indel_based_unlinked_assembly(&assembly));
}

fn replicate_support(assembly: AsmP, read: &ReadP, count: i32) {
    for _ in 0..count {
        assembly.m().add_junction_read(read);
    }
}

fn indel_assembly(junction: P<Junction>, read: &ReadP, junction_index: usize) -> AsmP {
    let assembly = JunctionAssembly::new_test(junction, read.get_bases(), read.get_base_quality(), junction_index as i32);
    assembly.m().add_junction_read(read);
    assembly.m().set_indel_coords(read.indel_coords());
    assembly
}

#[test]
fn test_indel_assembly_dedup() {
    let _env = test_env();
    let indel_junction_pos1 = create_indel_junction(CHR_1, 130, FORWARD);
    let indel_junction_neg1 = create_indel_junction(CHR_1, 131, REVERSE);
    let left_ref_bases = &REF_BASES_400[100..130];
    let right_ref_bases = &REF_BASES_400[131..171];
    let insert_bases = &REF_BASES_400[200..240];
    let indel_read_bases = format!("{left_ref_bases}{insert_bases}{right_ref_bases}");
    // both pairs of junctions share a read which supports them both, so are not deduped
    let indel_read = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 101, &indel_read_bases, "30M40I30M", CHR_1, 300, true);
    let indel_assembly_pos1 = indel_assembly(indel_junction_pos1, &indel_read, left_ref_bases.len());
    let indel_assembly_neg1 = indel_assembly(indel_junction_neg1, &indel_read, left_ref_bases.len() + insert_bases.len() + 1);

    // a pair where the pos assemblies don't match and the negatives do
    let indel_junction_pos2 = create_indel_junction(CHR_1, 140, FORWARD);
    let indel_junction_neg2 = create_indel_junction(CHR_1, 141, REVERSE);
    let insert_bases2 = &REF_BASES_400[210..255];
    let indel_read_bases2 = format!("{left_ref_bases}{insert_bases2}{right_ref_bases}");
    let indel_read2_pos = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 111, &indel_read_bases2, "30M45I30M", CHR_1, 300, true);
    let indel_assembly_pos2 = indel_assembly(indel_junction_pos2, &indel_read2_pos, left_ref_bases.len());
    let indel_read2_neg = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 111, &indel_read_bases2, "30M45I30M", CHR_1, 300, true);
    let indel_assembly_neg2 = indel_assembly(indel_junction_neg2, &indel_read2_neg, left_ref_bases.len() + insert_bases.len() + 1);

    // test 1 - pos assemblies don't match, neg assemblies do - first (existing) neg assembly is chosen so both second assemblies are dropped
    let mut existing_assemblies = vec![indel_assembly_pos1];
    let mut new_assemblies = vec![indel_assembly_pos2];
    let mut dedup_indels: Vec<AsmP> = Vec::new();
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(1, existing_assemblies.len());
    assert_eq!(1, new_assemblies.len());
    existing_assemblies.push(indel_assembly_pos2);
    replicate_support(indel_assembly_neg1, &indel_read, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO); // ensures this assembly is kept
    let mut new_assemblies = vec![indel_assembly_neg1];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(2, existing_assemblies.len());
    assert_eq!(1, new_assemblies.len());
    existing_assemblies.push(indel_assembly_neg1);
    let mut new_assemblies = vec![indel_assembly_neg2];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(2, existing_assemblies.len());
    assert!(!existing_assemblies.contains(&indel_assembly_pos2));
    assert_eq!(0, new_assemblies.len());
    assert!(dedup_indels.is_empty());

    // test 2 - last neg assembly is kept, first pair is removed
    indel_assembly_neg1.m().support_mut().clear();
    indel_assembly_neg1.m().add_junction_read(&indel_read);
    replicate_support(indel_assembly_neg2, &indel_read2_neg, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO); // ensures this assembly is kept
    let mut existing_assemblies = vec![indel_assembly_pos1, indel_assembly_pos2, indel_assembly_neg1];
    let mut new_assemblies = vec![indel_assembly_neg2];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(1, existing_assemblies.len());
    assert!(existing_assemblies.contains(&indel_assembly_pos2));
    assert!(!existing_assemblies.contains(&indel_assembly_pos1));
    assert!(!existing_assemblies.contains(&indel_assembly_neg1));
    assert_eq!(1, new_assemblies.len());
    assert!(dedup_indels.is_empty());

    // test 3 - positives are deduped first, neg deduped based on this
    replicate_support(indel_assembly_pos1, &indel_read, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO); // ensures this assembly is kept
    indel_assembly_neg1.m().support_mut().clear();
    indel_assembly_neg1.m().add_junction_read(&indel_read);
    indel_assembly_neg2.m().support_mut().clear();
    replicate_support(indel_assembly_neg2, &indel_read2_neg, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO + 1); // would be kept but will be dropped
    dedup_indels.clear();
    let mut existing_assemblies = vec![indel_assembly_pos1];
    let mut new_assemblies = vec![indel_assembly_pos2];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(1, existing_assemblies.len());
    assert_eq!(0, new_assemblies.len());
    assert!(dedup_indels.contains(&indel_assembly_pos2));
    // ensure the neg assembly 2 is now dedup despite it having more support
    let mut new_assemblies = vec![indel_assembly_neg1];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(1, existing_assemblies.len());
    assert_eq!(1, new_assemblies.len());
    assert!(dedup_indels.contains(&indel_assembly_pos2));
    existing_assemblies.push(indel_assembly_neg1);
    let mut new_assemblies = vec![indel_assembly_neg2];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(2, existing_assemblies.len());
    assert_eq!(0, new_assemblies.len());
    assert!(dedup_indels.is_empty());

    // test 4 - ensure that existing indels are only removed if the one being kept has a matching assembly
    let mut existing_assemblies = vec![indel_assembly_pos1, indel_assembly_neg1];
    let mut new_assemblies = vec![indel_assembly_neg2];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(1, existing_assemblies.len());
    assert!(existing_assemblies.contains(&indel_assembly_pos1));
    assert_eq!(1, new_assemblies.len());
}
