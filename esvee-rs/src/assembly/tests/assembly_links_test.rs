//! AssemblyLinksTest
use super::*;
use crate::assembly::junction_assembler::JunctionAssembler;
use crate::assembly::phase::assembly_linker::try_assembly_facing;
use crate::assembly::types::phase_group::PhaseGroup;
use crate::assembly::types::StructuralVariantType::{Bnd, Del, Dup, Inv};

const INSERTED_BASES: &str = "GGGGGGG";

/// new JunctionAssembly(junction, bases, baseQuals, junctionIndex): the Java tests reuse one quality array
fn asm(junction: P<Junction>, bases: &str, base_quals: &[u8], junction_index: usize) -> AsmP {
    JunctionAssembly::new_test(junction, bases.as_bytes(), base_quals, junction_index as i32)
}

fn check_link(link: &LinkP, first: AsmP, second: AsmP, sv_type: crate::assembly::types::StructuralVariantType) {
    assert!(first == link.first());
    assert!(second == link.second());
    assert_eq!(sv_type, link.sv_type());
}

#[test]
fn test_assembly_basic_splits() {
    let _env = test_env();
    let first_ref_bases = &REF_BASES_200[0..100];
    let second_ref_bases = &REF_BASES_200[100..200];
    // first a basic exact match junction
    let pos_junction = new_junction(CHR_1, 100, FORWARD);
    let neg_junction = new_junction(CHR_2, 200, REVERSE);
    let extension_length = 80;
    let first_extension_bases = &second_ref_bases[0..extension_length]; // first 80 bases of second's ref, exact match and no insert
    let first_assembly_bases = format!("{first_ref_bases}{first_extension_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(pos_junction, &first_assembly_bases, &base_quals, 99);
    let second_extension_bases = &first_ref_bases[20..100]; // last 80 bases of first's ref, again an exact match
    let second_assembly_bases = format!("{second_extension_bases}{second_ref_bases}");
    let second_assembly = asm(neg_junction, &second_assembly_bases, &base_quals, 80);
    // order passed in doesn't matter
    let link = try_assembly_overlap(second_assembly, first_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    check_link(&link, first_assembly, second_assembly, Bnd);

    // test 2: DEL with inserted bases
    let pos_junction = new_junction(CHR_1, 500, FORWARD);
    let neg_junction = new_junction(CHR_1, 1000, REVERSE);
    let first_assembly_bases = format!("{first_ref_bases}{INSERTED_BASES}{first_extension_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(pos_junction, &first_assembly_bases, &base_quals, 99);
    let second_assembly_bases = format!("{second_extension_bases}{INSERTED_BASES}{second_ref_bases}");
    let second_assembly = asm(neg_junction, &second_assembly_bases, &base_quals, 80 + INSERTED_BASES.len());
    let link = try_assembly_overlap(second_assembly, first_assembly).expect("link");
    assert_eq!(INSERTED_BASES, link.inserted_bases());
    check_link(&link, first_assembly, second_assembly, Del);
}

#[test]
fn test_assembly_dups() {
    let _env = test_env();
    let first_ref_bases = &REF_BASES_200[0..100];
    let second_ref_bases = &REF_BASES_200[100..200];
    let extension_length = 80;
    let first_extension_bases = &second_ref_bases[0..extension_length];
    let first_assembly_bases = format!("{first_ref_bases}{first_extension_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let second_extension_bases = &first_ref_bases[20..100];
    let second_assembly_bases = format!("{second_extension_bases}{second_ref_bases}");
    let pos_junction = new_junction(CHR_1, 1000, FORWARD);
    let neg_junction = new_junction(CHR_1, 500, REVERSE);
    let first_assembly = asm(pos_junction, &first_assembly_bases, &base_quals, 99);
    let second_assembly = asm(neg_junction, &second_assembly_bases, &base_quals, 80);
    let link = try_assembly_overlap(second_assembly, first_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    check_link(&link, first_assembly, second_assembly, Dup);

    // now with a mismatch in the overlapping bases in the second assembly - has to be close to the junction point
    let mismatch_second_extension_bases = format!("{}{}{}", &second_extension_bases[0..60], get_next_base(&second_extension_bases[60..61]), &second_extension_bases[61..extension_length]);
    let second_assembly_bases = format!("{mismatch_second_extension_bases}{second_ref_bases}");
    let second_assembly = asm(neg_junction, &second_assembly_bases, &base_quals, 80);
    let link = try_assembly_overlap(second_assembly, first_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    check_link(&link, first_assembly, second_assembly, Dup);

    // DUP with overlapping bases - the +ve assembly's extension bases start 20 bases into the ref sequence of the -ve
    let pos_junction = new_junction(CHR_1, 1000, FORWARD);
    let neg_junction = new_junction(CHR_1, 500, REVERSE);
    let extension_length = 60;
    let overlap_length = 20;
    let overlap_bases = &second_ref_bases[0..overlap_length]; // first 20 ref bases are shared
    let first_ref_bases = format!("{}{}", &REF_BASES_200[0..80], overlap_bases);
    let first_extension_bases = &second_ref_bases[overlap_length..overlap_length + extension_length];
    let first_assembly_bases = format!("{first_ref_bases}{first_extension_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(pos_junction, &first_assembly_bases, &base_quals, 99);
    let second_extension_bases = &first_ref_bases[overlap_length..overlap_length + extension_length];
    let second_assembly_bases = format!("{second_extension_bases}{second_ref_bases}");
    let second_assembly = asm(neg_junction, &second_assembly_bases, &base_quals, 60);
    let link = try_assembly_overlap(second_assembly, first_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    check_link(&link, first_assembly, second_assembly, Dup);
    assert_eq!(overlap_bases, link.overlap_bases());
}

#[test]
fn test_assembly_positive_inversion_splits() {
    let _env = test_env();
    let first_ref_bases = &REF_BASES_200[0..100];
    let second_ref_bases = &REF_BASES_200[100..200];
    // local inversion
    let first_junction = new_junction(CHR_1, 100, FORWARD);
    let second_junction = new_junction(CHR_1, 200, FORWARD);
    let first_extension_bases = reverse_complement_str(&second_ref_bases[20..100]);
    let first_assembly_bases = format!("{first_ref_bases}{first_extension_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(first_junction, &first_assembly_bases, &base_quals, 99);
    let second_extension_bases = reverse_complement_str(&first_ref_bases[20..100]);
    let second_assembly_bases = format!("{second_ref_bases}{second_extension_bases}");
    let second_assembly = asm(second_junction, &second_assembly_bases, &base_quals, 99);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    check_link(&link, first_assembly, second_assembly, Inv);

    // test 2: same again but with an insert sequence
    let first_assembly_bases = format!("{first_ref_bases}{INSERTED_BASES}{first_extension_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(first_junction, &first_assembly_bases, &base_quals, 99);
    let second_assembly_bases = format!("{}{}{}", second_ref_bases, reverse_complement_str(INSERTED_BASES), second_extension_bases);
    let second_assembly = asm(second_junction, &second_assembly_bases, &base_quals, 99);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert_eq!(INSERTED_BASES, link.inserted_bases());
    check_link(&link, first_assembly, second_assembly, Inv);

    // test 3: 1 base of homology
    let first_ref_bases = "CAAACAAAAATTAAAACAAAAAGGCTAAATACTGCGTGATTTTGTTTACACGGAATTCTAAAAAGAAAAAAAAAAAAAAGACAACACAGAGAAAAAAACA";
    let second_ref_bases = "AATTCAAAGGACTGGCTTTAGCCTCCAGCTACTGAACAGCTTTGGTTAATAAGAAGCACCTGCAGAAGACTAGAAGACAACAGGAAAGGGGTTGTCTAGT";
    let first_extension_bases = reverse_complement_str(&second_ref_bases[19..99]);
    let first_assembly_bases = format!("{first_ref_bases}{first_extension_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(first_junction, &first_assembly_bases, &base_quals, 99);
    let second_extension_bases = reverse_complement_str(&first_ref_bases[19..99]);
    let second_assembly_bases = format!("{second_ref_bases}{second_extension_bases}");
    let second_assembly = asm(second_junction, &second_assembly_bases, &base_quals, 99);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    assert_eq!("A", link.overlap_bases());
    check_link(&link, first_assembly, second_assembly, Inv);

    // test 4: again but with a base from the second's extension bases, so the mismatch must be factored into the offset
    let second_extension_bases = reverse_complement_str(&format!("{}{}", &first_ref_bases[19..62], &first_ref_bases[63..99]));
    let second_assembly_bases = format!("{second_ref_bases}{second_extension_bases}");
    let second_assembly = asm(second_junction, &second_assembly_bases, &base_quals, 99);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    assert_eq!("A", link.overlap_bases());
    check_link(&link, first_assembly, second_assembly, Inv);
}

#[test]
fn test_assembly_negative_inverted_splits() {
    let _env = test_env();
    let first_ref_bases = &REF_BASES_400[0..100];
    let second_ref_bases = &REF_BASES_400[100..200];
    // a negative inversion
    let first_junction = new_junction(CHR_1, 100, REVERSE);
    let second_junction = new_junction(CHR_1, 200, REVERSE);
    let extension_length = 80;
    let first_extension_bases = reverse_complement_str(&second_ref_bases[0..extension_length]); // the first 80 bases
    let first_assembly_bases = format!("{first_extension_bases}{first_ref_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(first_junction, &first_assembly_bases, &base_quals, extension_length);
    let second_extension_bases = reverse_complement_str(&first_ref_bases[0..extension_length]);
    let second_assembly_bases = format!("{second_extension_bases}{second_ref_bases}");
    let second_assembly = asm(second_junction, &second_assembly_bases, &base_quals, extension_length);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert!(link.inserted_bases().is_empty());
    check_link(&link, first_assembly, second_assembly, Inv);

    // test 2: same again but with an insert sequence - always taken from the first assembly without reverse-complimenting
    let inserted_base_length = INSERTED_BASES.len();
    let first_assembly_bases = format!("{first_extension_bases}{INSERTED_BASES}{first_ref_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(first_junction, &first_assembly_bases, &base_quals, extension_length + inserted_base_length);
    let second_assembly_bases = format!("{}{}{}", second_extension_bases, reverse_complement_str(INSERTED_BASES), second_ref_bases);
    let second_assembly = asm(second_junction, &second_assembly_bases, &base_quals, extension_length + inserted_base_length);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert_eq!(INSERTED_BASES, link.inserted_bases());
    check_link(&link, first_assembly, second_assembly, Inv);

    // test 3: with a 5-base overlap
    let overlap_length = 5;
    let new_extension_length = extension_length - overlap_length;
    let first_assembly_bases = format!("{}{}", &first_extension_bases[0..extension_length - overlap_length], first_ref_bases);
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = asm(first_junction, &first_assembly_bases, &base_quals, new_extension_length);
    let second_assembly_bases = format!("{}{}", &second_extension_bases[0..extension_length - overlap_length], second_ref_bases);
    let second_assembly = asm(second_junction, &second_assembly_bases, &base_quals, new_extension_length);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert_eq!(&first_ref_bases[0..overlap_length], link.overlap_bases());
    check_link(&link, first_assembly, second_assembly, Inv);
}

#[test]
fn test_inversion_support() {
    let _env = test_env();
    let ref_sequence = REF_BASES_400;
    let mut ref_genome = MockRefGenome::new_indexing(false);
    ref_genome.put(CHR_1, ref_sequence);
    // local inversion
    let first_junction = new_junction(CHR_1, 100, FORWARD);
    let second_junction = new_junction(CHR_1, 200, FORWARD);
    // a read which is a discordant candidate for one junction and a junction read for the other
    let junc_read_bases = format!("{}{}", &ref_sequence[51..101], reverse_complement_str(&ref_sequence[151..201]));
    let junc_read1 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 51, &junc_read_bases, "50M50S", CHR_1, 1, false);
    let junc_read2 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 52, &junc_read_bases[1..], "49M50S", CHR_1, 1, false);
    let first_assembly = JunctionAssembler::new(first_junction, None, None).process_junction(&[junc_read1.clone(), junc_read2.clone()])[0];
    let junc_read_bases2 = format!("{}{}", &ref_sequence[151..201], reverse_complement_str(&ref_sequence[51..101]));
    let junc_read3 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 151, &junc_read_bases2, "50M50S", CHR_1, 1, false);
    let junc_read4 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 152, &junc_read_bases2[1..], "49M50S", CHR_1, 1, false);
    let second_assembly = JunctionAssembler::new(second_junction, None, None).process_junction(&[junc_read3.clone(), junc_read4])[0];

    let disc_read = |read_id: &str, start: i32, bases: &str, cigar: &str, mate_start: i32| -> ReadP {
        let read = create_read_mate(read_id, CHR_1, start, bases, cigar, CHR_1, mate_start, false);
        // setSecondOfPairFlag(true) alone, leaving the first-of-pair flag set
        read.m().bam_record_mut().set_flag(FLAG_SECOND_OF_PAIR, true);
        read
    };
    let disc_read2 = disc_read(junc_read2.id(), 120, &ref_sequence[120..175], "55M", 51);
    second_assembly.m().add_candidate_support(disc_read2.clone());
    // the same junction read from the first assembly as a candidate for the second - this won't end up being counted as support
    // since it matches the exact read as support in the first assembly
    second_assembly.m().add_candidate_support(junc_read1.clone());
    let disc_read1 = disc_read(junc_read1.id(), 120, &ref_sequence[120..175], "55M", 51);
    second_assembly.m().add_candidate_support(disc_read1.clone());
    let disc_read3 = disc_read(junc_read3.id(), 20, &ref_sequence[20..80], "60M", 151);
    first_assembly.m().add_candidate_support(disc_read3.clone());

    let phase_group = PhaseGroup::new(first_assembly, Some(second_assembly));
    build_phase_sets(&ref_genome, phase_group);
    assert_eq!(1, phase_group.phase_sets().len());
    let phase_set = phase_group.phase_sets()[0];
    assert_eq!(2, phase_set.assemblies().len());
    let link = phase_set.assembly_links()[0];
    assert!(link.inserted_bases().is_empty());
    check_link(&link, first_assembly, second_assembly, Inv);
    assert_eq!(3, first_assembly.support_count());
    assert_eq!(2, get_support_type_count(first_assembly, SupportType::Junction));
    assert_eq!(1, get_support_type_count(first_assembly, SupportType::Discordant));
    assert_eq!(4, second_assembly.support_count());
    assert_eq!(2, get_support_type_count(second_assembly, SupportType::Junction));
    assert_eq!(2, get_support_type_count(second_assembly, SupportType::Discordant));
}

#[test]
fn test_ext_base_only_matched_links() {
    let _env = test_env();
    // two assemblies only match in their extension bases, so need to form the full insert sequence drawing from both
    let assembly_ref_bases = &REF_BASES_400[1..101];
    let assembly_ext_bases = &REF_BASES_200[0..100];
    let assembly_bases = format!("{assembly_ref_bases}{assembly_ext_bases}");
    let assembly1 = create_assembly(CHR_1, 100, FORWARD, &assembly_bases, assembly_ref_bases.len() as i32 - 1);
    let assembly_ref_bases2 = &REF_BASES_400[300..400];
    let assembly_ext_bases2 = &REF_BASES_200[30..130];
    let assembly_bases2 = format!("{assembly_ext_bases2}{assembly_ref_bases2}");
    let assembly2 = create_assembly(CHR_2, 300, REVERSE, &assembly_bases2, assembly_ext_bases2.len() as i32);
    let link = try_assembly_overlap(assembly1, assembly2).expect("link");
    assert_eq!(130, link.inserted_bases().len());
    assert_eq!(&REF_BASES_200[0..130], link.inserted_bases());
    // cannot form a link purely in extension bases if the assemblies form a local DEL or DUP
    let assembly2 = create_assembly(CHR_1, 300, REVERSE, &assembly_bases2, assembly_ext_bases2.len() as i32);
    assert!(try_assembly_overlap(assembly1, assembly2).is_none());
}

#[test]
fn test_facing_link_conditions() {
    let _env = test_env();
    let assembly_ref_bases = &REF_BASES_400[151..251];
    let assembly_ext_bases = &REF_BASES_600[1..50];
    let assembly_bases = format!("{assembly_ref_bases}{assembly_ext_bases}");
    let assembly_pos = create_assembly(CHR_1, 250, FORWARD, &assembly_bases, assembly_ref_bases.len() as i32 - 1);
    let assembly_ref_bases2 = &REF_BASES_400[100..200];
    let assembly_ext_bases2 = &REF_BASES_600[1..50];
    let assembly_bases2 = format!("{assembly_ext_bases2}{assembly_ref_bases2}");
    let ext_length2 = assembly_ext_bases2.len() as i32;
    let assembly_neg = create_assembly(CHR_1, 100, REVERSE, &assembly_bases2, ext_length2);
    // no shared read
    assert!(try_assembly_facing(assembly_neg, assembly_pos, &[]).is_none());
    let junc_read1a = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 151, &assembly_bases, "100M50S", CHR_1, 100, false);
    junc_read1a.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, true);
    assembly_pos.m().add_junction_read(&junc_read1a);
    let junc_read1b = create_read_mate(junc_read1a.id(), CHR_1, 100, &assembly_bases2, "50S100M", CHR_1, 151, true);
    assembly_neg.m().add_junction_read(&junc_read1b);
    assert!(try_assembly_facing(assembly_neg, assembly_pos, &[]).is_some());

    // invalid if too close or far away
    let assembly_neg = create_assembly(CHR_1, 225, REVERSE, &assembly_bases2, ext_length2);
    let junc_read1b = create_read_mate(junc_read1a.id(), CHR_1, 225, &assembly_bases2, "50S100M", CHR_1, 151, true);
    assembly_neg.m().add_junction_read(&junc_read1b);
    assert!(try_assembly_facing(assembly_neg, assembly_pos, &[]).is_none());

    // different chromosome
    let assembly_neg = create_assembly(CHR_2, 100, REVERSE, &assembly_bases2, ext_length2);
    let junc_read1b = create_read_mate(junc_read1a.id(), CHR_2, 100, &assembly_bases2, "50S100M", CHR_1, 151, true);
    assembly_neg.m().add_junction_read(&junc_read1b);
    assert!(try_assembly_facing(assembly_neg, assembly_pos, &[]).is_none());

    // matched reads have incorrect orientation
    let assembly_neg = create_assembly(CHR_1, 100, REVERSE, &assembly_bases2, ext_length2);
    let junc_read1b = create_read_mate(junc_read1a.id(), CHR_1, 100, &assembly_bases2, "50S100M", CHR_1, 151, true);
    assembly_neg.m().add_junction_read(&junc_read1b);
    assert!(try_assembly_facing(assembly_neg, assembly_pos, &[]).is_some());
    assembly_neg.m().support_mut().clear();
    junc_read1b.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, true);
    assembly_neg.m().add_junction_read(&junc_read1b);
    assert!(try_assembly_facing(assembly_neg, assembly_pos, &[]).is_none());

    // or reads facing away
    assembly_pos.m().support_mut().clear();
    junc_read1a.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, false);
    assembly_pos.m().add_junction_read(&junc_read1a);
    assert!(try_assembly_facing(assembly_neg, assembly_pos, &[]).is_none());
}
