//! LocalLinksTest
use super::*;
use crate::assembly::assembly_deduper::dedup_proximate_assemblies;
use crate::assembly::common_utils::IndelCoords;
use crate::assembly::phase::local_sequence_matcher::LocalSequenceMatcher;
use crate::assembly::ptr::A;
use crate::assembly::types::phase_group::PhaseGroup;
use crate::assembly::types::StructuralVariantType;

fn ref_genome_600() -> MockRefGenome {
    let mut ref_genome = MockRefGenome::new_indexing(true);
    ref_genome.put(CHR_1, REF_BASES_600);
    ref_genome
}

#[test]
fn test_assembly_local_ref_match() {
    let _env = test_env();
    // must be long enough to test the local ref genome sequence but not repetitive
    let ref_genome = ref_genome_600();
    let local_sequence_matcher = LocalSequenceMatcher::new(Some(&ref_genome), 200);

    // first a basic exact match junction
    let pos_junction = new_junction(CHR_1, 300, FORWARD);
    let assembly_ref_bases = ref_genome.base_string(CHR_1, 200, 300);
    let assembly_extension_bases = ref_genome.base_string(CHR_1, 400, 449);
    let assembly = new_assembly(pos_junction, &format!("{assembly_ref_bases}{assembly_extension_bases}"), assembly_ref_bases.len() - 1);
    let assembly_link = local_sequence_matcher.try_local_assembly_link(assembly).expect("link");
    assert_eq!(StructuralVariantType::Del, assembly_link.sv_type());
    let local_assembly = assembly_link.second();
    assert_eq!(400, local_assembly.junction().position);
    assert!(local_assembly.junction().orient.is_reverse());
    assert_eq!(50, local_assembly.extension_length());
    let local_sequence = format!("{}{}", ref_genome.base_string(CHR_1, 251, 300), ref_genome.base_string(CHR_1, 400, 449));
    assert_eq!(local_sequence, local_assembly.form_full_sequence());

    let neg_junction = new_junction(CHR_1, 300, REVERSE);
    let assembly_ref_bases = ref_genome.base_string(CHR_1, 300, 400);
    let assembly_extension_bases = ref_genome.base_string(CHR_1, 401, 450);
    let assembly = new_assembly(neg_junction, &format!("{assembly_extension_bases}{assembly_ref_bases}"), assembly_extension_bases.len());
    let assembly_link = local_sequence_matcher.try_local_assembly_link(assembly).expect("link");
    assert_eq!(StructuralVariantType::Dup, assembly_link.sv_type());
    let local_assembly = assembly_link.first();
    assert_eq!(450, local_assembly.junction().position);
    assert!(local_assembly.junction().orient.is_forward());
    assert_eq!(50, local_assembly.extension_length());
    let local_sequence = format!("{}{}", ref_genome.base_string(CHR_1, 401, 450), ref_genome.base_string(CHR_1, 300, 349));
    assert_eq!(local_sequence, local_assembly.form_full_sequence());
}

#[test]
fn test_local_indel_vs_junction() {
    let _env = test_env();
    // 2 DELs with the internal long section having a duplication of 35 bases, represented as both an INDEL and a soft-clip assembly
    // first DEL:   1:100-200 - assemblies 1 & 2
    // INDEL:       1:250:251 - assemblies 3 & 4, with insert of 35 bases matching the sequence 251-285
    // dedup asm    1:285:1 - assembly 5
    // second DEL:  1:350-450 - assemblies 6 & 7
    let ref_genome = ref_genome_600();
    let junction1 = new_junction(CHR_1, 100, FORWARD);
    let ref_bases1 = ref_genome.base_string(CHR_1, 1, 100);
    let ext_bases1 = ref_genome.base_string(CHR_1, 200, 249);
    let assembly_bases1 = format!("{ref_bases1}{ext_bases1}");
    let assembly1 = new_assembly(junction1, &assembly_bases1, ref_bases1.len() - 1);
    let junction2 = new_junction(CHR_1, 200, REVERSE);
    let ref_bases2 = ref_genome.base_string(CHR_1, 200, 250);
    let ext_bases2 = ref_genome.base_string(CHR_1, 51, 100);
    let assembly_bases2 = format!("{ext_bases2}{ref_bases2}");
    let assembly2 = new_assembly(junction2, &assembly_bases2, ext_bases2.len());

    // read pair linking assemblies 1 & 2
    let junc_read1 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 1, &assembly_bases1, &make_cigar_string(&ref_bases1, 0, ext_bases1.len()), CHR_1, 200, true);
    assembly1.m().add_junction_read(&junc_read1);
    let junc_read2 = create_read_mate(junc_read1.id(), CHR_1, 200, &assembly_bases2, &make_cigar_string(&ref_bases2, ext_bases2.len(), 0), CHR_1, 1, false);
    assembly2.m().add_junction_read(&junc_read2);

    let junction3 = new_junction(CHR_1, 250, FORWARD);
    junction3.m().set_indel_status(true);
    let duplicated_bases = ref_genome.base_string(CHR_1, 251, 285);
    let mut coords = IndelCoords::new(250, 251, duplicated_bases.len() as i32);
    coords.set_inserted_bases(duplicated_bases.clone());
    let indel_coords = A::new(coords);
    let ref_bases3 = ref_genome.base_string(CHR_1, 200, 250);
    let assembly_bases3 = format!("{ref_bases3}{duplicated_bases}");
    let assembly3 = new_assembly(junction3, &assembly_bases3, ref_bases3.len() - 1);
    assembly3.m().set_indel_coords(Some(indel_coords.clone()));
    let junction4 = new_junction(CHR_1, 251, REVERSE);
    junction4.m().set_indel_status(true);
    let ref_bases4 = ref_genome.base_string(CHR_1, 251, 350);
    let assembly_bases4 = format!("{duplicated_bases}{ref_bases4}");
    let assembly4 = new_assembly(junction4, &assembly_bases4, duplicated_bases.len());
    assembly4.m().set_indel_coords(Some(indel_coords.clone()));

    // indel read and read linking to assembly 2
    let junc_read2b = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 200, &assembly_bases2, &make_cigar_string(&ref_bases2, ext_bases2.len(), 0), CHR_1, 221, true);
    assembly2.m().add_junction_read(&junc_read2b);
    let read3_bases = format!("{}{}{}", ref_genome.base_string(CHR_1, 221, 250), duplicated_bases, ref_genome.base_string(CHR_1, 251, 280));
    let junc_read3 = create_read_mate(junc_read2b.id(), CHR_1, 221, &read3_bases, "30M35I30M", CHR_1, 200, false);
    assembly3.m().add_junction_read(&junc_read3);
    assembly4.m().add_junction_read(&junc_read3);

    let junction5 = new_junction(CHR_1, 285, FORWARD);
    let ref_bases5 = ref_genome.base_string(CHR_1, 200, 285);
    let assembly_bases5 = format!("{ref_bases5}{duplicated_bases}");
    let assembly5 = new_assembly(junction5, &assembly_bases5, ref_bases5.len() - 1);

    let junction6 = new_junction(CHR_1, 350, FORWARD);
    let ref_bases6 = ref_genome.base_string(CHR_1, 251, 350);
    let ext_bases6 = ref_genome.base_string(CHR_1, 450, 500);
    let assembly_bases6 = format!("{ref_bases6}{ext_bases6}");
    let assembly6 = new_assembly(junction6, &assembly_bases6, ref_bases6.len() - 1);
    let junc_read4 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 221, &read3_bases, "30M35I30M", CHR_1, 200, true);
    assembly3.m().add_junction_read(&junc_read4);
    assembly4.m().add_junction_read(&junc_read4);
    let junc_read6 = create_read_mate(junc_read4.id(), CHR_1, 251, &assembly_bases6, &make_cigar_string(&ref_bases6, 0, ext_bases6.len()), CHR_1, 221, false);
    assembly6.m().add_junction_read(&junc_read6);

    let junction7 = new_junction(CHR_1, 450, REVERSE);
    let ref_bases7 = ref_genome.base_string(CHR_1, 450, 550);
    let ext_bases7 = ref_genome.base_string(CHR_1, 300, 350);
    let assembly_bases7 = format!("{ext_bases7}{ref_bases7}");
    let assembly7 = new_assembly(junction7, &assembly_bases7, ext_bases7.len());

    // read pair linking assemblies 6 & 7 (as in Java the first read of this pair is never added, and read 1 is added again)
    let junc_read6b = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 251, &assembly_bases6, &make_cigar_string(&ref_bases6, 0, ext_bases6.len()), CHR_1, 450, true);
    assembly1.m().add_junction_read(&junc_read1);
    let junc_read7 = create_read_mate(junc_read6b.id(), CHR_1, 450, &assembly_bases2, &make_cigar_string(&ref_bases7, ext_bases7.len(), 0), CHR_1, 251, false);
    assembly7.m().add_junction_read(&junc_read7);

    let mut junction_group_assemblies = vec![assembly3, assembly4];
    let mut candidate_assemblies = vec![assembly5];
    let mut dedup_indels: Vec<AsmP> = Vec::new();
    dedup_proximate_assemblies(&mut junction_group_assemblies, &mut candidate_assemblies, &mut dedup_indels);
    junction_group_assemblies.extend(candidate_assemblies.iter().copied());
    assert_eq!(2, junction_group_assemblies.len());
    assert!(junction_group_assemblies.contains(&assembly3));
    assert!(junction_group_assemblies.contains(&assembly4));

    // now test linking logic
    let phase_group = PhaseGroup::new(assembly1, Some(assembly2));
    for assembly in [assembly3, assembly4, assembly6, assembly7] {
        PhaseGroup::add_assembly(phase_group, assembly);
    }
    build_phase_sets(&ref_genome, phase_group);
    assert_eq!(1, phase_group.phase_sets().len());
    let phase_set = phase_group.phase_sets()[0];
    assert_eq!(6, phase_set.assemblies().len());
    let expected = [
        (LinkType::Split, assembly1, assembly2),
        (LinkType::Facing, assembly2, assembly3),
        (LinkType::Indel, assembly3, assembly4),
        (LinkType::Facing, assembly4, assembly6),
        (LinkType::Split, assembly6, assembly7),
    ];
    for (i, (link_type, first, second)) in expected.iter().enumerate() {
        let link = phase_set.assembly_links()[i];
        assert_eq!(*link_type, link.type_());
        assert!(*first == link.first());
        assert!(*second == link.second());
    }
    assert_eq!(duplicated_bases, phase_set.assembly_links()[2].inserted_bases());
}

#[test]
fn test_local_indel_required_overlaps() {
    let _env = test_env();
    // test 1: there is sufficient overlap in ref and extension bases to form a link
    let ref_genome = ref_genome_600();
    // modelled as a DEL 1:100:1 to 1:300:-1
    let junction1 = new_junction(CHR_1, 100, FORWARD);
    let ref_bases1 = ref_genome.base_string(CHR_1, 1, 100);
    let ext_bases1 = ref_genome.base_string(CHR_1, 300, 349);
    let assembly1 = new_assembly(junction1, &format!("{ref_bases1}{ext_bases1}"), ref_bases1.len() - 1);
    let junction2 = new_junction(CHR_1, 300, REVERSE);
    let ref_bases2 = ref_genome.base_string(CHR_1, 300, 399);
    let ext_bases2 = ref_genome.base_string(CHR_1, 51, 100);
    let assembly2 = new_assembly(junction2, &format!("{ext_bases2}{ref_bases2}"), ext_bases2.len());
    assert!(try_assembly_overlap(assembly1, assembly2).is_some());

    // test 2: with a single mismatch, trigger the seed-search logic
    let ext_bases2_mod = format!("{}{}{}", &ext_bases2[0..20], get_next_base(&ext_bases2[20..21]), &ext_bases2[21..]);
    let assembly2 = new_assembly(junction2, &format!("{ext_bases2_mod}{ref_bases2}"), ext_bases2.len());
    assert!(try_assembly_overlap(assembly1, assembly2).is_some());

    // test 3: sufficient overlap in ref and extension bases but too many mis-matches
    let ext_bases2_mod = format!("{}TTTTT{}", &ext_bases2[0..20], &ext_bases2[25..]);
    let assembly2 = new_assembly(junction2, &format!("{ext_bases2_mod}{ref_bases2}"), ext_bases2.len());
    assert!(try_assembly_overlap(assembly1, assembly2).is_none());

    // test 4: insufficient overlap in ref and extension bases to form a link
    let ext_bases2_mod = &ext_bases2[35..];
    let assembly2 = new_assembly(junction2, &format!("{ext_bases2_mod}{ref_bases2}"), ext_bases2_mod.len());
    let ref_bases1_mod = format!("{}{}", &ext_bases2[35..50], &ref_bases2[0..20]);
    let ext_bases1_mod = format!("{}{}", &ref_bases2[20..34], get_next_base(&ref_bases2[34..35]));
    let assembly1 = new_assembly(junction1, &format!("{ref_bases1_mod}{ext_bases1_mod}"), ref_bases1_mod.len() - 1);
    assert!(try_assembly_overlap(assembly1, assembly2).is_none());
}
