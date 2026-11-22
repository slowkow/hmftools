//! RefBaseExtensionTest
use super::*;
use crate::assembly::junction_assembler::JunctionAssembler;
use crate::assembly::ref_base_extender::{extend_ref_bases, trim_assembly_ref_bases};

#[test]
fn test_junction_indel_read_ref_base_building() {
    let _env = test_env();
    // how ref bases are built when there are indels in the junction read's ref bases
    let junction = new_junction(CHR_1, 200, FORWARD);
    let mut junction_assembler = JunctionAssembler::new(junction, None, None);
    let ext_bases = &REF_BASES_200[100..150];

    // first an assembly with inserts in the ref bases
    let inserted_bases = "AGAGAGAGAG";
    let read_bases = format!("{}{}{}{}", &REF_BASES_400[161..181], inserted_bases, &REF_BASES_400[181..201], ext_bases);
    let read = create_read(&READ_ID_GENERATOR.next_id(), 161, &read_bases, "20M10I20M50S");
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 162, &read_bases[1..], "19M10I20M50S");
    let assemblies = junction_assembler.process_junction(&[read, read2]);
    assert_eq!(1, assemblies.len());
    let assembly = assemblies[0];
    // the inserted bases are included but the extension sequence is not extended to include all the ref bases to the lowest start pos
    assert_eq!(50, assembly.extension_length());
    assert_eq!(161, assembly.ref_base_position());
    assert_eq!(50, assembly.ref_base_length());
    assert_eq!(&read_bases[0..50], assembly.form_ref_base_sequence());

    // repeat for a series of junction reads with a delete in their ref bases
    let read_bases = format!("{}{}{}", &REF_BASES_400[161..171], &REF_BASES_400[181..201], ext_bases);
    let read = create_read(&READ_ID_GENERATOR.next_id(), 161, &read_bases, "10M10D20M50S");
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 162, &read_bases[1..], "9M10D20M50S");
    let assemblies = junction_assembler.process_junction(&[read, read2]);
    assert_eq!(1, assemblies.len());
    let assembly = assemblies[0];
    assert_eq!(50, assembly.extension_length());
    assert_eq!(161, assembly.ref_base_position());
    assert_eq!(30, assembly.ref_base_length());
    let assembly_ref_bases = format!("{}{}", &REF_BASES_400[161..171], &REF_BASES_400[181..201]);
    assert_eq!(assembly_ref_bases, assembly.form_ref_base_sequence());
}

#[test]
fn test_basic_ref_base_extension() {
    let _env = test_env();
    let mut ref_genome = MockRefGenome::new();
    ref_genome.put(CHR_1, REF_BASES_200);
    let pos_junction = new_junction(CHR_1, 150, FORWARD);
    let initial_ref_bases = &REF_BASES_200[130..151];
    let assembly_bases = format!("{}{}", initial_ref_bases, &REF_BASES_200[180..]);
    let base_quals = build_default_base_quals(assembly_bases.len());
    let junction_index = initial_ref_bases.len() as i32 - 1;
    let assembly = JunctionAssembly::new_test(pos_junction, assembly_bases.as_bytes(), &base_quals, junction_index);
    assert_eq!(130, assembly.ref_base_position());

    // a collection of junction mate and discordant reads to test as support
    let mate1 = create_read(&READ_ID_GENERATOR.next_id(), 121, &REF_BASES_200[121..151], "30M");
    mate1.m().mark_junction_mate();
    let mate2_bases = format!("{}A{}", &REF_BASES_200[121..140], &REF_BASES_200[141..151]);
    let mate2 = create_read(&READ_ID_GENERATOR.next_id(), 121, &mate2_bases, "30M");
    mate2.m().mark_junction_mate();
    // too many mismatches
    let disc1_bases = format!("{}TTT{}", &REF_BASES_200[100..115], &REF_BASES_200[118..130]);
    let disc1 = create_read(&READ_ID_GENERATOR.next_id(), 100, &disc1_bases, "30M");
    // supports the ref bases, just enough overlap
    let disc2 = create_read(&READ_ID_GENERATOR.next_id(), 100, &REF_BASES_200[100..145], "45M");
    // supports the ref bases but not enough overlap
    let disc3 = create_read(&READ_ID_GENERATOR.next_id(), 90, &REF_BASES_200[90..120], "30M");
    // too large a gap
    let disc4 = create_read(&READ_ID_GENERATOR.next_id(), 30, &REF_BASES_200[30..60], "30M");
    let mut candidate_support = vec![mate1, mate2, disc1, disc2, disc3, disc4];
    extend_ref_bases(assembly, &mut candidate_support, Some(&ref_genome), false);
    assert_eq!(30, assembly.ref_base_position());
    assert_eq!(4, assembly.support_count());
    trim_assembly_ref_bases(assembly, 20);
    assert_eq!(90, assembly.ref_base_position());
    assert_eq!(&REF_BASES_200[90..151], assembly.form_ref_base_sequence());

    // repeat for -ve junction
    let neg_junction = new_junction(CHR_1, 50, REVERSE);
    let initial_ref_bases = &REF_BASES_200[50..71];
    let assembly_bases = format!("{}{}", &REF_BASES_200[0..20], initial_ref_bases);
    let base_quals = build_default_base_quals(assembly_bases.len());
    let assembly = JunctionAssembly::new_test(neg_junction, assembly_bases.as_bytes(), &base_quals, 20);
    assert_eq!(70, assembly.ref_base_position());

    let mate1 = create_read(&READ_ID_GENERATOR.next_id(), 50, &REF_BASES_200[50..100], "50M");
    mate1.m().mark_junction_mate();
    let mate2_bases = format!("{}A{}", &REF_BASES_200[60..75], &REF_BASES_200[76..90]);
    let mate2 = create_read(&READ_ID_GENERATOR.next_id(), 60, &mate2_bases, "30M");
    mate2.m().mark_junction_mate();
    // supports the ref bases with enough overlap
    let disc1 = create_read(&READ_ID_GENERATOR.next_id(), 75, &REF_BASES_200[75..105], "30M");
    // too many mismatches
    let disc2_bases = format!("{}TTT{}", &REF_BASES_200[100..115], &REF_BASES_200[118..160]);
    let disc2 = create_read(&READ_ID_GENERATOR.next_id(), 100, &disc2_bases, "60M");
    // supports the ref bases but too large a gap so will be cleaned up later
    let disc3 = create_read(&READ_ID_GENERATOR.next_id(), 120, &REF_BASES_200[120..150], "30M");
    let mut candidate_support = vec![mate1, mate2, disc1, disc2, disc3];
    extend_ref_bases(assembly, &mut candidate_support, Some(&ref_genome), false);
    assert_eq!(159, assembly.ref_base_position());
    assert_eq!(3, assembly.support_count());
    trim_assembly_ref_bases(assembly, 20);
    assert_eq!(104, assembly.ref_base_position());
    assert_eq!(&REF_BASES_200[50..105], assembly.form_ref_base_sequence());
}
