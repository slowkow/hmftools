//! SpecificAlignmentsTest
use super::*;
use crate::assembly::read::read_utils::is_valid_support_coords_vs_junction;
use crate::assembly::common_utils::positions_within;
use crate::assembly::ref_base_extender::find_assembly_candidate_extensions;
use crate::assembly::types::phase_group::PhaseGroup;
use crate::assembly::types::support_read::SupP;

fn reversed(read: &ReadP) {
    read.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, true);
}
fn second_in_pair(read: &ReadP) {
    set_second_in_pair(read.m().bam_record_mut());
}
fn mate_unmapped(read: &ReadP) {
    read.m().bam_record_mut().set_flag(FLAG_MATE_UNMAPPED, true);
}
fn read_unmapped(read: &ReadP) {
    read.m().bam_record_mut().set_flag(FLAG_READ_UNMAPPED, true);
}
fn link_mates(first: &ReadP, second: &ReadP) {
    first.m().set_mate_read(second);
    second.m().set_mate_read(first);
}
fn rc(bases: &str) -> String {
    reverse_complement_str(bases)
}
fn find_support(assembly: AsmP, type_: SupportType, read: &ReadP) -> SupP {
    assembly.support().iter().find(|x| x.type_() == type_ && x.id() == read.id()).expect("support read").clone()
}

fn check_support_read_sequence_indices(assembly: AsmP, seq_index_start: i32, seq_index_end: i32, orientation: Orientation) {
    for read in assembly.support() {
        assert_eq!(orientation, read.full_assembly_orientation());
        assert!(positions_within(read.full_assembly_index_start(), read.full_assembly_index_end(), seq_index_start, seq_index_end));
    }
}

#[test]
fn test_unmapped_sgl_pair() {
    let _env = test_env();
    // two facing SGLs which both have unmapped reads which form consistent extension sequences
    let mut ref_genome = MockRefGenome::new_indexing(true);
    ref_genome.put(CHR_1, REF_BASES_400);
    let ref_bases1 = ref_genome.base_string(CHR_1, 200, 300);
    let ext_bases1 = &REF_BASES_600[420..500];
    let assembly_bases1 = format!("{ext_bases1}{ref_bases1}");
    let assembly1 = create_assembly(CHR_1, 200, REVERSE, &assembly_bases1, ext_bases1.len() as i32);
    let ref_bases2 = ref_genome.base_string(CHR_1, 200, 300);
    let ext_bases2 = &REF_BASES_600[100..180];
    let assembly_bases2 = format!("{ref_bases2}{ext_bases2}");
    let assembly2 = create_assembly(CHR_1, 300, FORWARD, &assembly_bases2, ref_bases2.len() as i32 - 1);
    let mut candidates1: Vec<ReadP> = Vec::new();
    let mut candidates2: Vec<ReadP> = Vec::new();

    // a pair of mate junction reads
    let junc_read1 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 200, &assembly_bases1[0..100], "50S50M", CHR_1, 251, true);
    set_mate_cigar(&junc_read1, "50M50S");
    let junc_read2 = create_read_mate(junc_read1.id(), CHR_1, 251, &assembly_bases2[51..], "50M50S", CHR_1, 200, false);
    set_mate_cigar(&junc_read2, "50S50M");
    second_in_pair(&junc_read2);
    assembly1.m().add_junction_read(&junc_read1);
    assembly2.m().add_junction_read(&junc_read2);
    link_mates(&junc_read1, &junc_read2);
    candidates1.push(junc_read2.clone());
    candidates2.push(junc_read1.clone());

    // a shared junction read with mate unmapped on the first assembly side
    let read_bases = format!("{}{}{}", &ext_bases1[40..], ref_bases1, &ext_bases1[0..10]);
    let shared_read1 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 200, &read_bases, "10S101M10S", NO_CHROMOSOME_NAME, NO_POSITION, true);
    mate_unmapped(&shared_read1);
    assembly1.m().add_junction_read(&shared_read1);
    assembly2.m().add_junction_read(&shared_read1);

    // more junction reads whose mates extend out the unmapped sequence
    let junc_read1b = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 200, &assembly_bases1[0..100], "50S50M", NO_CHROMOSOME_NAME, NO_POSITION, true);
    mate_unmapped(&junc_read1b);
    set_mate_cigar(&junc_read1b, NO_CIGAR);
    // extends the extension seq 20 bases further back
    let unmapped_read1 = create_read_mate(junc_read1b.id(), CHR_1, 200, &REF_BASES_600[400..480], NO_CIGAR, CHR_1, 200, false);
    read_unmapped(&unmapped_read1);
    reversed(&unmapped_read1);
    assembly1.m().add_junction_read(&junc_read1b);
    candidates1.push(unmapped_read1.clone());
    link_mates(&junc_read1b, &unmapped_read1);
    second_in_pair(&unmapped_read1);

    let junc_read2b = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 251, &assembly_bases2[51..], "50M50S", NO_CHROMOSOME_NAME, NO_POSITION, true);
    mate_unmapped(&junc_read2b);
    set_mate_cigar(&junc_read2b, NO_CIGAR);
    // was 100-180, so extends 20 bases further out
    let unmapped_read2 = create_read_mate(junc_read2b.id(), CHR_1, 251, &rc(&REF_BASES_600[120..200]), NO_CIGAR, CHR_1, 251, false);
    read_unmapped(&unmapped_read2);
    assembly2.m().add_junction_read(&junc_read2b);
    candidates2.push(unmapped_read2.clone());
    link_mates(&junc_read2b, &unmapped_read2);
    second_in_pair(&unmapped_read2);

    find_assembly_candidate_extensions(assembly1, &candidates1);
    find_assembly_candidate_extensions(assembly2, &candidates2);
    let phase_group = PhaseGroup::new(assembly2, Some(assembly1));
    build_phase_sets(&ref_genome, phase_group);

    // first check unmapped extensions
    assert_eq!(5, assembly1.support_count());
    assert_eq!(4, assembly2.support_count());
    assert_eq!(100, assembly1.extension_length());
    assert_eq!(&REF_BASES_600[400..500], assembly1.form_junction_sequence());
    assert_eq!(100, assembly2.extension_length());
    assert_eq!(&REF_BASES_600[100..200], assembly2.form_junction_sequence());
    assert_eq!(100, find_support(assembly1, SupportType::Extension, &unmapped_read1).junction_read_start_distance());
    assert_eq!(-20, find_support(assembly2, SupportType::Extension, &unmapped_read2).junction_read_start_distance());
    assert_eq!(1, phase_group.phase_sets().len());
    let phase_set = phase_group.phase_sets()[0];
    assert_eq!(1, phase_set.assembly_links().len());
    assert_eq!(2, phase_set.assemblies().len());
    let _ = has_assembly_link(phase_set.assembly_links(), assembly1, assembly2, LinkType::Facing);
    let assembly_alignment = AssemblyAlignment::from_phase_set(phase_set);
    assert_eq!(301, assembly_alignment.full_sequence_length());
    let full_sequence = format!("{}{}{}", &REF_BASES_600[400..500], ref_bases1, &REF_BASES_600[100..200]);
    assert_eq!(Some(full_sequence.as_str()), assembly_alignment.full_sequence());

    // confirm full seq read indices
    assert_eq!(50, find_support(assembly1, SupportType::Junction, &junc_read1).full_assembly_index_start());
    assert_eq!(90, find_support(assembly1, SupportType::Junction, &shared_read1).full_assembly_index_start());
    assert_eq!(151, find_support(assembly1, SupportType::JunctionMate, &junc_read1).full_assembly_index_start());
    assert_eq!(151, find_support(assembly2, SupportType::Junction, &junc_read1).full_assembly_index_start());
    assert_eq!(90, find_support(assembly2, SupportType::Junction, &shared_read1).full_assembly_index_start());
    assert_eq!(221, find_support(assembly2, SupportType::Extension, &unmapped_read2).full_assembly_index_start());
}

#[test]
fn test_chained_links() {
    let _env = test_env();
    // links are:
    // chr1:100:1 -> chr2:50:-1
    // chr2:50:-1 -> chr2:200:1
    // chr2:200:1 -> chr3:200:1
    // chr3:200:1 -> chr3:50:-1 (ie segment is reverse-linked)
    // chr3:50:-1 -> chr1:250:-1
    // this is also the order of both the assemblies and links as they're added to the phase set
    let mut ref_genome = MockRefGenome::new_indexing(true);
    for chromosome in [CHR_1, CHR_2, CHR_3] {
        ref_genome.put(chromosome, REF_BASES_400);
    }
    let g = |chromosome: &str, start: i32, end: i32| ref_genome.base_string(chromosome, start, end);
    // chr1:100:1 - links to chr2:50:-1
    let ref_bases1 = g(CHR_1, 1, 100);
    let assembly1 = create_assembly(CHR_1, 100, FORWARD, &format!("{}{}", ref_bases1, g(CHR_2, 50, 150)), ref_bases1.len() as i32 - 1);
    // chr2:50:-1 - links to chr1:100:1 and faces chr2:200:1
    let ext_bases2 = g(CHR_1, 51, 100);
    let assembly2 = create_assembly(CHR_2, 50, REVERSE, &format!("{}{}", ext_bases2, g(CHR_2, 50, 150)), ext_bases2.len() as i32);
    // chr2:200:1 - links to chr3:200:-1 inverted to faces chr2:50
    let ref_bases3 = g(CHR_2, 101, 200);
    let ext_bases3 = rc(&g(CHR_3, 151, 201));
    let assembly3 = create_assembly(CHR_2, 200, FORWARD, &format!("{ref_bases3}{ext_bases3}"), ref_bases3.len() as i32 - 1);
    // chr3:200:1 - links to chr2:200:1 inverted and faces chr3:50:-1
    let ref_bases4 = g(CHR_3, 101, 200);
    let assembly4 = create_assembly(CHR_3, 200, FORWARD, &format!("{}{}", ref_bases4, rc(&g(CHR_2, 151, 200))), ref_bases4.len() as i32 - 1);
    // chr3:50:-1 - links to chr1:250:-1 inverted and faces chr3:200:1
    let ext_bases5 = rc(&g(CHR_1, 250, 300));
    let assembly5 = create_assembly(CHR_3, 50, REVERSE, &format!("{}{}", ext_bases5, g(CHR_3, 50, 150)), ext_bases5.len() as i32);
    // chr1:250:-1 - links to chr3:50:-1
    let ext_bases6 = g(CHR_3, 50, 100);
    let assembly6 = create_assembly(CHR_1, 250, REVERSE, &format!("{}{}", rc(&ext_bases6), g(CHR_1, 250, 350)), ext_bases6.len() as i32);

    // now add junction, mate and discordant reads
    // each junction read has a mate either on the ref side or across a facing link to another segment
    let pos_junc_read_cigar = "100M50S";
    let neg_junc_read_cigar = "50S100M";
    let mate_read_cigar = "100M";
    let mut candidates1: Vec<ReadP> = Vec::new();
    let mut candidates2: Vec<ReadP> = Vec::new();
    let mut candidates3: Vec<ReadP> = Vec::new();
    let mut candidates4: Vec<ReadP> = Vec::new();
    let mut candidates5: Vec<ReadP> = Vec::new();
    let mut candidates6: Vec<ReadP> = Vec::new();
    let junc_bases = g(CHR_1, 51, 150);
    let next_id = || READ_ID_GENERATOR.next_id();

    let junc_read1 = create_read_mate(&next_id(), CHR_1, 51, &junc_bases, pos_junc_read_cigar, CHR_2, 101, true);
    assembly1.m().add_junction_read(&junc_read1);
    // a mate read in the next segment
    let mate_read = create_read_mate(junc_read1.id(), CHR_2, 101, &g(CHR_2, 101, 200), mate_read_cigar, CHR_1, 1, false);
    reversed(&mate_read);
    second_in_pair(&mate_read);
    candidates2.push(mate_read.clone());
    link_mates(&mate_read, &junc_read1);

    // a junction read and mate both in the first segment
    let junc_read1b = create_read_mate(&next_id(), CHR_1, 51, &junc_bases, pos_junc_read_cigar, CHR_1, 1, false);
    reversed(&junc_read1b);
    assembly1.m().add_junction_read(&junc_read1b);
    let mate_read = create_read_mate(junc_read1b.id(), CHR_1, 1, &g(CHR_1, 1, 100), mate_read_cigar, CHR_1, 51, true);
    second_in_pair(&mate_read);
    candidates1.push(mate_read.clone());
    link_mates(&mate_read, &junc_read1b);

    // a junction read with its mate soft-clipping into the next junction
    let junc_read2 = create_read_mate(&next_id(), CHR_2, 50, &junc_bases, neg_junc_read_cigar, CHR_2, 101, true);
    assembly2.m().add_junction_read(&junc_read2);
    let mate_bases = format!("{}{}", g(CHR_2, 101, 200), &ext_bases3[0..10]);
    let mate_read = create_read_mate(junc_read2.id(), CHR_2, 101, &mate_bases, "100M10S", CHR_2, 50, false);
    second_in_pair(&mate_read);
    reversed(&mate_read);
    candidates2.push(mate_read.clone());
    assembly3.m().add_junction_read(&mate_read);
    link_mates(&mate_read, &junc_read2);

    // an assembly 2 junction read with its mate in the following segment
    let junc_read2b = create_read_mate(&next_id(), CHR_2, 50, &junc_bases, neg_junc_read_cigar, CHR_3, 100, false);
    assembly2.m().add_junction_read(&junc_read2b);
    // as in Java, this mate is given the id of junction read 2
    let mate_read = create_read_mate(junc_read2.id(), CHR_3, 100, &g(CHR_3, 100, 199), "100M", CHR_2, 50, false);
    second_in_pair(&mate_read);
    reversed(&mate_read);
    candidates4.push(mate_read.clone());
    link_mates(&mate_read, &junc_read2b);

    // an assembly 3 junction read with its mate soft-clipped in the same section
    let junc_read3 = create_read_mate(&next_id(), CHR_2, 101, &junc_bases, pos_junc_read_cigar, CHR_2, 50, false);
    reversed(&junc_read3);
    assembly3.m().add_junction_read(&junc_read3);
    let mate_bases = format!("{}{}", &ext_bases2[ext_bases2.len() - 10..], g(CHR_2, 50, 149));
    let mate_read = create_read_mate(junc_read3.id(), CHR_3, 50, &mate_bases, "10S100M", CHR_2, 101, true);
    second_in_pair(&mate_read);
    assembly2.m().add_junction_read(&mate_read);
    candidates3.push(mate_read.clone());
    candidates2.push(junc_read3.clone());
    link_mates(&mate_read, &junc_read3);

    // an assembly 3 junction read with its mate in the next section
    let junc_read3b = create_read_mate(&next_id(), CHR_2, 101, &junc_bases, pos_junc_read_cigar, CHR_2, 101, false);
    assembly3.m().add_junction_read(&junc_read3b);
    let mate_read = create_read_mate(junc_read3b.id(), CHR_3, 100, &g(CHR_3, 100, 199), "100M", CHR_2, 101, false);
    second_in_pair(&mate_read);
    candidates4.push(mate_read.clone());
    link_mates(&mate_read, &junc_read3b);

    // an assembly 4 junction read with its mate soft-clipped in the same segment
    let junc_read4 = create_read_mate(&next_id(), CHR_3, 101, &junc_bases, pos_junc_read_cigar, CHR_3, 50, true);
    assembly4.m().add_junction_read(&junc_read4);
    let mate_bases = format!("{}{}", &ext_bases5[ext_bases5.len() - 10..], g(CHR_3, 50, 149));
    let mate_read = create_read_mate(junc_read4.id(), CHR_3, 50, &mate_bases, "10S100M", CHR_3, 101, false);
    reversed(&mate_read);
    second_in_pair(&mate_read);
    assembly5.m().add_junction_read(&mate_read);
    candidates4.push(mate_read.clone());
    candidates5.push(junc_read4.clone());
    link_mates(&mate_read, &junc_read4);

    // an assembly 4 junction read with its mate in the final segment
    let junc_read4b = create_read_mate(&next_id(), CHR_3, 101, &junc_bases, pos_junc_read_cigar, CHR_1, 300, true);
    assembly4.m().add_junction_read(&junc_read4b);
    let mate_read = create_read_mate(junc_read4b.id(), CHR_1, 300, &g(CHR_1, 300, 350), "101M", CHR_3, 101, false);
    reversed(&mate_read);
    second_in_pair(&mate_read);
    candidates6.push(mate_read.clone());
    link_mates(&mate_read, &junc_read4b);

    // an assembly 5 junction read with its mate in the final segment
    let junc_read5 = create_read_mate(&next_id(), CHR_3, 50, &junc_bases, neg_junc_read_cigar, CHR_1, 300, true);
    assembly5.m().add_junction_read(&junc_read5);
    let mate_read = create_read_mate(junc_read5.id(), CHR_1, 300, &g(CHR_1, 300, 350), "101M", CHR_3, 101, false);
    second_in_pair(&mate_read);
    reversed(&mate_read);
    candidates6.push(mate_read.clone());
    link_mates(&mate_read, &junc_read5);

    // an assembly 6 with its mate in the previous segment
    let junc_read6 = create_read_mate(&next_id(), CHR_1, 250, &junc_bases, neg_junc_read_cigar, CHR_3, 100, false);
    assembly6.m().add_junction_read(&junc_read6);
    let mate_read = create_read_mate(junc_read6.id(), CHR_3, 100, &g(CHR_3, 100, 200), "101M", CHR_1, 250, false);
    second_in_pair(&mate_read);
    candidates5.push(mate_read.clone());
    link_mates(&mate_read, &junc_read6);

    // add discordant reads as candidates to check they are assigned as support correctly
    let disc_cigar = "50M";
    let valid = |read: &ReadP, assembly: AsmP| is_valid_support_coords_vs_junction(read, assembly.junction().is_forward(), assembly.junction().position);
    let disc_read1 = create_read_mate(&next_id(), CHR_1, 20, &REF_BASES_400[20..70], disc_cigar, CHR_3, 100, true);
    let disc_read2 = create_read_mate(disc_read1.id(), CHR_3, 100, &REF_BASES_400[100..150], disc_cigar, CHR_1, 20, false);
    second_in_pair(&disc_read2);
    assert!(valid(&disc_read1, assembly1));
    candidates1.push(disc_read1.clone());
    assert!(valid(&disc_read2, assembly4));
    candidates4.push(disc_read2.clone());
    assert!(!valid(&disc_read2, assembly5));
    candidates5.push(disc_read2.clone());
    let disc_read3 = create_read_mate(&next_id(), CHR_2, 100, &REF_BASES_400[100..150], disc_cigar, CHR_1, 300, true);
    let disc_read4 = create_read_mate(disc_read3.id(), CHR_1, 300, &REF_BASES_400[300..350], disc_cigar, CHR_2, 100, false);
    reversed(&disc_read4);
    second_in_pair(&disc_read4);
    assert!(!valid(&disc_read3, assembly2));
    assert!(valid(&disc_read3, assembly3));
    candidates3.push(disc_read3.clone());
    assert!(!valid(&disc_read4, assembly1));
    assert!(valid(&disc_read4, assembly6));
    candidates6.push(disc_read4.clone());

    for (assembly, candidates) in [(assembly1, &candidates1), (assembly2, &candidates2), (assembly3, &candidates3), (assembly4, &candidates4), (assembly5, &candidates5), (assembly6, &candidates6)] {
        find_assembly_candidate_extensions(assembly, candidates);
    }
    let phase_group = PhaseGroup::new(assembly1, Some(assembly2));
    for assembly in [assembly3, assembly4, assembly5, assembly6] {
        PhaseGroup::add_assembly(phase_group, assembly);
    }
    build_phase_sets(&ref_genome, phase_group);

    assert_eq!(1, phase_group.phase_sets().len());
    let phase_set = phase_group.phase_sets()[0];
    assert_eq!(5, phase_set.assembly_links().len());
    // the Java test calls these without asserting on the result
    let _ = has_assembly_link(phase_set.assembly_links(), assembly1, assembly2, LinkType::Split);
    let _ = has_assembly_link(phase_set.assembly_links(), assembly2, assembly3, LinkType::Facing);
    let _ = has_assembly_link(phase_set.assembly_links(), assembly3, assembly4, LinkType::Split);
    let _ = has_assembly_link(phase_set.assembly_links(), assembly4, assembly5, LinkType::Facing);
    let _ = has_assembly_link(phase_set.assembly_links(), assembly5, assembly6, LinkType::Split);
    let counts: Vec<usize> = [assembly1, assembly2, assembly3, assembly4, assembly5, assembly6].iter().map(|a| a.support_count()).collect();
    assert_eq!(vec![3, 5, 6, 5, 3, 3], counts);
    let disc_counts: Vec<usize> = [assembly1, assembly2, assembly3, assembly4, assembly5, assembly6].iter().map(|a| get_support_type_count(*a, SupportType::Discordant)).collect();
    assert_eq!(vec![1, 1, 1, 2, 0, 2], disc_counts);

    // check final sequence and read full sequence indices
    let full_sequence = format!("{}{}{}{}", g(CHR_1, 1, 100), g(CHR_2, 50, 200), rc(&g(CHR_3, 50, 200)), g(CHR_1, 250, 350));
    let assembly_alignment = AssemblyAlignment::from_phase_set(phase_set);
    assert_eq!(Some(full_sequence.as_str()), assembly_alignment.full_sequence());
    check_support_read_sequence_indices(assembly1, 0, 150, FORWARD);
    check_support_read_sequence_indices(assembly2, 50, 261, FORWARD);
    check_support_read_sequence_indices(assembly3, 50, 301, FORWARD);
    check_support_read_sequence_indices(assembly4, 151, 411, REVERSE);
    check_support_read_sequence_indices(assembly5, 251, 452, REVERSE);
    check_support_read_sequence_indices(assembly6, 352, 552, FORWARD);
}

#[test]
fn test_chr363_chain() {
    let _env = test_env();
    // links are:
    // chr1:200:1 -> chr1:300:1 - local INV
    // chr1:100:-1 -> chr1:200:1 - short facing TI
    // chr1:100:-1 -> chr1:300:1 - longer facing TI
    // chr1:100:-1 -> chr2:100:-1 - BND via a remote ref match
    let mut ref_genome = MockRefGenome::new_indexing(true);
    ref_genome.put(CHR_1, REF_BASES_400);
    ref_genome.put(CHR_2, REF_BASES_600);
    let g = |chromosome: &str, start: i32, end: i32| ref_genome.base_string(chromosome, start, end);
    let next_id = || READ_ID_GENERATOR.next_id();
    // chr1:100:-1 - has reads soft-clipping at 200 and others going out further and soft-clipping at 300
    let ref_bases1 = g(CHR_1, 100, 150);
    let ext_bases1 = g(CHR_2, 100, 200);
    let assembly_bases1 = format!("{}{}", rc(&ext_bases1), ref_bases1);
    let assembly1 = create_assembly(CHR_1, 100, REVERSE, &assembly_bases1, ext_bases1.len() as i32);
    // chr1:200:1
    let ref_bases2 = g(CHR_1, 100, 200);
    let assembly_bases2 = format!("{}{}", ref_bases2, rc(&g(CHR_1, 251, 301)));
    let assembly2 = create_assembly(CHR_1, 200, FORWARD, &assembly_bases2, ref_bases2.len() as i32 - 1);
    // chr1:300:1
    let ref_bases3 = g(CHR_1, 251, 300);
    let assembly_bases3 = format!("{}{}", ref_bases3, rc(&g(CHR_1, 151, 201)));
    let assembly3 = create_assembly(CHR_1, 300, FORWARD, &assembly_bases3, ref_bases3.len() as i32 - 1);
    let mut candidates1: Vec<ReadP> = Vec::new();
    let mut candidates2: Vec<ReadP> = Vec::new();
    let mut candidates3: Vec<ReadP> = Vec::new();

    // first junction reads for each assembly
    let junc_read1 = create_read_mate(&next_id(), CHR_1, 100, &assembly_bases1, "50S50M", CHR_2, 150, true);
    reversed(&junc_read1);
    set_mate_cigar(&junc_read1, "100M");
    let junc_read1_bases = format!("{}{}", &ext_bases1[51..101], ref_bases1);
    let junc_read1b = create_read_mate(&next_id(), CHR_1, 100, &junc_read1_bases, "50S50M", CHR_1, 251, true);
    let junc_read1c = create_read_mate(&next_id(), CHR_1, 100, &junc_read1_bases, "50S50M", CHR_1, 151, true);
    set_mate_cigar(&junc_read1c, "50M50S");
    let junc_read1d = create_read_mate(&next_id(), CHR_1, 100, &junc_read1_bases, "50S50M", CHR_1, 150, true);
    assembly1.m().add_junction_read(&junc_read1d);
    let mate1d = create_read_mate(junc_read1d.id(), CHR_1, 150, &g(CHR_1, 150, 250), "100M", CHR_1, 100, false);
    reversed(&mate1d);
    link_mates(&junc_read1d, &mate1d);
    second_in_pair(&mate1d);
    assembly1.m().add_junction_read(&junc_read1d);
    candidates1.push(mate1d.clone());
    // build sufficient support to warrant a branch
    for _ in 0..4 {
        let junc_read = clone_read(&junc_read1d, &next_id());
        let mate_read = clone_read(&mate1d, &next_id());
        link_mates(&junc_read, &mate_read);
        assembly1.m().add_junction_read(&junc_read);
        candidates1.push(mate_read);
    }
    assembly1.m().add_junction_read(&junc_read1);
    assembly1.m().add_junction_read(&junc_read1b);
    assembly1.m().add_junction_read(&junc_read1c);
    candidates2.push(junc_read1.clone());
    candidates3.push(junc_read1.clone());
    candidates3.push(junc_read1b.clone());
    let junc_read1 = clone_read(&junc_read1, &next_id());
    assembly1.m().add_junction_read(&junc_read1);
    candidates2.push(junc_read1.clone());
    candidates3.push(junc_read1.clone());

    let junc_read2 = create_read_mate(&next_id(), CHR_1, 151, &assembly_bases2, "50M50S", CHR_2, 150, true);
    assembly2.m().add_junction_read(&junc_read2);
    let junc_read2b = create_read_mate(&next_id(), CHR_1, 151, &assembly_bases2, "50M50S", CHR_1, 220, false);
    assembly2.m().add_junction_read(&junc_read2b);
    let junc_read2c = create_read_mate(junc_read1c.id(), CHR_1, 151, &assembly_bases2, "50M50S", CHR_1, 100, false);
    set_mate_cigar(&junc_read2c, "50S50M");
    reversed(&junc_read2c);
    assembly2.m().add_junction_read(&junc_read2c);
    link_mates(&junc_read1c, &junc_read2c);
    second_in_pair(&junc_read2c);
    candidates1.push(junc_read2c.clone());
    candidates2.push(junc_read1c.clone());
    for _ in 0..4 {
        let junc_read = clone_read(&junc_read1c, &next_id());
        let mate_read = clone_read(&junc_read2c, &next_id());
        assembly1.m().add_junction_read(&junc_read);
        assembly2.m().add_junction_read(&mate_read);
        link_mates(&junc_read, &mate_read);
        candidates1.push(mate_read);
        candidates2.push(junc_read);
    }
    let junc_read3 = create_read_mate(&next_id(), CHR_1, 251, &assembly_bases3, "50M50S", CHR_2, 150, true);
    assembly3.m().add_junction_read(&junc_read3);
    // the mate of assembly junc read 1b
    let junc_read3b = create_read_mate(junc_read1b.id(), CHR_1, 251, &assembly_bases3, "50M50S", CHR_1, 100, true);
    link_mates(&junc_read3b, &junc_read1b);
    second_in_pair(&junc_read3b);
    assembly3.m().add_junction_read(&junc_read3b);
    candidates1.push(junc_read3b.clone());

    // a junction read spanning both ends of the short TI - counts towards both assembly 1 and 2
    let read_bases = format!("{}{}{}", rc(&g(CHR_2, 100, 110)), &REF_BASES_400[100..201], &REF_BASES_400[300..310]);
    let shared_read1 = create_read_mate(&next_id(), CHR_1, 100, &read_bases, "10S101M10S", CHR_2, 150, true);
    reversed(&shared_read1);
    assembly1.m().add_junction_read(&shared_read1);
    assembly2.m().add_junction_read(&shared_read1);
    // a junction mate covering the short TI junction
    let read_bases = format!("{}{}", &REF_BASES_400[151..201], &REF_BASES_400[251..301]);
    let shared_read2 = create_read_mate(&next_id(), CHR_1, 151, &read_bases, "50M50S", CHR_2, 150, true);
    candidates1.push(shared_read2.clone());
    assembly2.m().add_junction_read(&shared_read2);
    // a junction mate extending past the short TI junction
    let mate1 = create_read_mate(&next_id(), CHR_1, 151, &REF_BASES_400[151..251], "100M", CHR_2, 150, true);
    candidates1.push(mate1);
    // discordant read hitting the longer TI junction
    let read_bases = format!("{}{}", &REF_BASES_400[251..301], rc(&REF_BASES_400[191..201]));
    let shared_read3 = create_read_mate(&next_id(), CHR_1, 251, &read_bases, "50M10S", CHR_2, 150, true);
    candidates1.push(shared_read3.clone());
    assembly3.m().add_junction_read(&shared_read3);
    // junction mate linked to assembly 2
    let mate2 = create_read_mate(junc_read2b.id(), CHR_1, 220, &REF_BASES_400[220..270], "50M", CHR_1, 151, false);
    link_mates(&mate2, &junc_read2b);
    second_in_pair(&mate2);
    candidates2.push(mate2.clone());
    candidates3.push(mate2);

    find_assembly_candidate_extensions(assembly1, &candidates1);
    find_assembly_candidate_extensions(assembly2, &candidates2);
    find_assembly_candidate_extensions(assembly3, &candidates3);
    // the Java test also builds remote reads and a remote region on chr 2 here, which nothing then uses

    let phase_group = PhaseGroup::new(assembly1, Some(assembly2));
    PhaseGroup::add_assembly(phase_group, assembly3);
    build_phase_sets(&ref_genome, phase_group);
    assert_eq!(1, phase_group.phase_sets().len());
    let phase_set = phase_group.phase_sets()[0];
    assert_eq!(3, phase_set.assembly_links().len());
    assert_eq!(4, phase_set.assemblies().len());
    let branched_assembly = *phase_set.assemblies().iter().find(|x| **x != assembly1 && **x != assembly2 && **x != assembly3).expect("branched assembly");
    assert_eq!(201, assembly1.ref_base_length());
    assert_eq!(201, assembly3.ref_base_length());
    assert_eq!(101, assembly2.ref_base_length());
    assert_eq!(101, branched_assembly.ref_base_length());
    let _ = has_assembly_link(phase_set.assembly_links(), assembly2, assembly3, LinkType::Split);
    let _ = has_assembly_link(phase_set.assembly_links(), assembly1, assembly3, LinkType::Facing);
    let _ = has_assembly_link(phase_set.assembly_links(), branched_assembly, assembly2, LinkType::Facing);
    assert_eq!(21, assembly1.support_count());
    assert_eq!(14, assembly2.support_count());
    assert_eq!(5, assembly3.support_count());
    assert_eq!(0, get_support_type_count(assembly1, SupportType::Discordant));
    assert_eq!(0, get_support_type_count(assembly2, SupportType::Discordant));
    // the bases from the repeated link are now omitted
    let expected_full_sequence = format!("{}{}C{}", rc(&g(CHR_2, 100, 200)), g(CHR_1, 100, 200), rc(&g(CHR_1, 100, 300)));
    let assembly_alignment = AssemblyAlignment::from_phase_set(phase_set);
    assert_eq!(Some(expected_full_sequence.as_str()), assembly_alignment.full_sequence());
}
