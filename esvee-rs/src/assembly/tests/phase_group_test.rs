//! PhaseGroupTest
use std::collections::BTreeMap;

use super::*;
use crate::assembly::phase::phase_group_builder::{form_local_phase_groups, PhaseGroupBuildWriter, PhaseGroupBuilder};
use crate::assembly::remote_region_finder::find_remote_regions;
use crate::assembly::types::junction_group::{JuncGroupP, JunctionGroup};
use crate::assembly::types::phase_group::{PhaseGroup, PhaseGroupP};
use crate::assembly::types::remote_region::{RemoteRegion, RemoteRegionP};
use crate::assembly::types::RemoteReadType;
use crate::common::java_id_hash::IdHashSet;

#[test]
fn test_remote_region_building() {
    let _env = test_env();
    let pos_junction = new_junction(CHR_1, 150, FORWARD);
    let ref_bases = &REF_BASES_200[101..151];
    let assembly_bases = format!("{}{}", ref_bases, &REF_BASES_200[0..50]);
    let assembly = new_assembly(pos_junction, &assembly_bases, ref_bases.len() - 1);
    let junc_cigar = "50M50S";
    let mate1 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 101, &assembly_bases, junc_cigar, CHR_2, 1000, false);
    let mate2 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 101, &assembly_bases, junc_cigar, CHR_2, 2000, false);
    let remote_junction_mates = vec![mate1, mate2];
    let disc1 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 50, &REF_BASES_200[50..100], "50M", CHR_2, 2100, false);
    // low base-qual disc read
    let disc2 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 50, &REF_BASES_200[50..100], "50M", CHR_2, 3000, false);
    let length = disc2.bam_record().base_qualities.len();
    disc2.m().bam_record_mut().base_qualities = build_base_quals(length, 10);
    // low alignment score
    let disc3 = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 50, &REF_BASES_200[50..100], "50M", CHR_2, 4000, false);
    disc3.m().bam_record_mut().set_integer_attribute(ALIGNMENT_SCORE_ATTRIBUTE, 50);
    let discordant_reads = vec![disc1, disc2, disc3];

    // test region merging
    let supp1 = SupplementaryReadData::new(CHR_3, 1000, SUPP_POS_STRAND, "100M", 0, 0);
    let junc_supp1 = Read::new(create_sam_record_full(&READ_ID_GENERATOR.next_id(), CHR_1, 50, &assembly_bases, junc_cigar, CHR_1, 50, false, false, Some(&supp1)));
    let supp2 = SupplementaryReadData::new(CHR_3, 1500, SUPP_POS_STRAND, "100M", 0, 0);
    let junc_supp2 = Read::new(create_sam_record_full(&READ_ID_GENERATOR.next_id(), CHR_1, 50, &assembly_bases, junc_cigar, CHR_1, 50, false, false, Some(&supp2)));
    let supp_junction_reads = vec![junc_supp1, junc_supp2];

    find_remote_regions(assembly, &discordant_reads, &remote_junction_mates, &supp_junction_reads);

    assert_eq!(3, assembly.remote_regions().len());
    let expected = [(CHR_2, 1000, 1099, 1), (CHR_2, 2000, 2149, 2), (CHR_3, 1000, 1599, 2)];
    for (i, (chromosome, start, end, read_count)) in expected.iter().enumerate() {
        let region = assembly.remote_regions()[i];
        assert_eq!(*chromosome, region.chromosome);
        assert_eq!(*start, region.start());
        assert_eq!(*end, region.end());
        assert_eq!(*read_count, region.read_count());
    }
}

fn group_with(phase_groups: &[PhaseGroupP], assembly: AsmP) -> PhaseGroupP {
    *phase_groups.iter().find(|x| x.assemblies().contains(&assembly)).expect("phase group")
}

#[test]
fn test_local_phase_group_building() {
    let _env = test_env();
    let assembly_bases = &REF_BASES_200[0..100];
    let assembly_at = |position: i32, orientation: Orientation| new_assembly(new_junction(CHR_1, position, orientation), assembly_bases, 50);
    let pos_assembly1 = assembly_at(100, FORWARD);
    let neg_assembly1 = assembly_at(200, REVERSE);
    let pos_assembly2 = assembly_at(1500, FORWARD);
    let neg_assembly2 = assembly_at(2700, REVERSE);
    let pos_assembly3 = assembly_at(4000, FORWARD);
    let neg_assembly3 = assembly_at(4200, REVERSE);
    let neg_assembly4 = assembly_at(4500, REVERSE);
    let pos_assembly4 = assembly_at(4600, FORWARD);

    let read1 = create_concordant_read(&READ_ID_GENERATOR.next_id(), 50, assembly_bases, "50M50S", 250);
    pos_assembly1.m().add_junction_read(&read1);
    let read2 = create_concordant_read(&READ_ID_GENERATOR.next_id(), 200, assembly_bases, "50S50M", 10);
    set_mate_cigar(&read2, "50M");
    neg_assembly1.m().add_junction_read(&read2);
    let read3 = create_concordant_read(&READ_ID_GENERATOR.next_id(), 3950, assembly_bases, "50M50S", 4800);
    pos_assembly3.m().add_junction_read(&read3);
    let read4 = create_concordant_read(&READ_ID_GENERATOR.next_id(), 4200, assembly_bases, "50S50M", 3500);
    set_mate_cigar(&read4, "100M");
    neg_assembly3.m().add_junction_read(&read4);
    let read5 = create_concordant_read(&READ_ID_GENERATOR.next_id(), 4550, assembly_bases, "50M50S", 3800);
    pos_assembly4.m().add_junction_read(&read5);
    let read6 = create_concordant_read(&READ_ID_GENERATOR.next_id(), 4500, assembly_bases, "50S50M", 3800);
    set_mate_cigar(&read6, "100M");
    neg_assembly4.m().add_junction_read(&read6);

    let mut writer = PhaseGroupBuildWriter::disabled();
    let junction_group: JuncGroupP = P::new(JunctionGroup::new(pos_assembly1.junction()));
    for assembly in [pos_assembly2, pos_assembly3, pos_assembly4, neg_assembly1, neg_assembly2, neg_assembly3, neg_assembly4] {
        junction_group.m().add_junction(assembly.junction());
    }
    junction_group.m().add_junction_assemblies(&[pos_assembly1, pos_assembly2, pos_assembly3, pos_assembly4, neg_assembly1, neg_assembly2, neg_assembly3, neg_assembly4]);

    // LocalGroupBuilder.run
    let mut phase_group_set: IdHashSet<PhaseGroupP> = IdHashSet::new();
    form_local_phase_groups(junction_group, &mut phase_group_set, &mut writer);
    let phase_groups = phase_group_set.to_vec();

    assert_eq!(2, phase_groups.len());
    let phase_group = group_with(&phase_groups, pos_assembly1);
    assert!(phase_group.assemblies().contains(&neg_assembly1));
    let phase_group = group_with(&phase_groups, pos_assembly3);
    assert_eq!(4, phase_group.assembly_count());
    for assembly in [pos_assembly3, pos_assembly4, neg_assembly3, neg_assembly4] {
        assert!(phase_group.assemblies().contains(&assembly));
    }
}

fn new_group(assemblies_for_junctions: &[AsmP], assemblies: &[AsmP]) -> JuncGroupP {
    let group: JuncGroupP = P::new(JunctionGroup::new(assemblies_for_junctions[0].junction()));
    for assembly in &assemblies_for_junctions[1..] {
        group.m().add_junction(assembly.junction());
    }
    group.m().add_junction_assemblies(assemblies);
    group
}

fn new_region(chromosome: &str, start: i32, end: i32, read_id: &str) -> RemoteRegionP {
    P::new(RemoteRegion::new(intern(chromosome), start, end, read_id, RemoteReadType::Discordant))
}

#[test]
fn test_remote_phase_group_building() {
    let env = test_env();
    let mut writer = PhaseGroupBuildWriter::disabled();
    let assembly_bases = &REF_BASES_200[0..50];
    let assembly1 = create_assembly(CHR_1, 100, FORWARD, assembly_bases, 25);
    let assembly2 = create_assembly(CHR_1, 500, FORWARD, assembly_bases, 25);
    let assembly3 = create_assembly(CHR_1, 4100, FORWARD, assembly_bases, 25);
    let assembly4 = create_assembly(CHR_1, 6100, FORWARD, assembly_bases, 25);
    let assembly5 = create_assembly(CHR_1, 6500, FORWARD, assembly_bases, 25);
    let junction_group1 = new_group(&[assembly1, assembly2], &[assembly1, assembly2]);
    let junction_group2 = new_group(&[assembly3], &[assembly3]);
    let junction_group3 = new_group(&[assembly4, assembly5], &[assembly4, assembly5]);
    let assembly6 = create_assembly(CHR_2, 100, REVERSE, assembly_bases, 25);
    let assembly7 = create_assembly(CHR_2, 500, REVERSE, assembly_bases, 25);
    let assembly8 = create_assembly(CHR_2, 4100, REVERSE, assembly_bases, 25);
    let assembly9 = create_assembly(CHR_2, 6100, REVERSE, assembly_bases, 25);
    let assembly10 = create_assembly(CHR_2, 6500, REVERSE, assembly_bases, 25);
    let junction_group4 = new_group(&[assembly6, assembly7], &[assembly6, assembly7]);
    // as in the Java test, this group is given assembly 1 as well as assembly 8
    let junction_group5 = new_group(&[assembly8], &[assembly1, assembly8]);
    let junction_group6 = new_group(&[assembly9, assembly10], &[assembly9, assembly10]);
    let mut junction_group_map: BTreeMap<String, Vec<JuncGroupP>> = BTreeMap::new();
    junction_group_map.insert(CHR_1.to_string(), vec![junction_group1, junction_group2, junction_group3]);
    junction_group_map.insert(CHR_2.to_string(), vec![junction_group4, junction_group5, junction_group6]);
    let all_junction_groups: Vec<JuncGroupP> = junction_group_map.values().flatten().copied().collect();

    let builder = PhaseGroupBuilder::new(&env.config, &junction_group_map);
    let (phase_groups, _) = builder.build_remote_groups(&all_junction_groups, &mut writer);
    assert_eq!(0, phase_groups.len());

    // now run again with links between the groups
    let read_id1 = READ_ID_GENERATOR.next_id();
    let read_id2 = READ_ID_GENERATOR.next_id();
    let read_id3 = READ_ID_GENERATOR.next_id();
    // sufficient linking criteria
    let region1 = new_region(CHR_1, 6000, 6200, &read_id1);
    region1.m().add_read_details(&read_id2, 6000, 6100, RemoteReadType::JunctionMate);
    let region1b = new_region(CHR_2, 6000, 6100, &read_id3);
    assembly1.m().add_remote_regions(&[region1, region1b]);
    assembly4.m().add_junction_read(&create_read(&read_id1, 6000, assembly_bases, "50M"));
    assembly4.m().add_junction_read(&create_read(&read_id2, 6000, assembly_bases, "50M"));
    // insufficient read count overlap
    assembly10.m().add_junction_read(&create_read(&read_id3, 6000, assembly_bases, "50M"));

    // link between 3 assemblies
    let read_id4 = READ_ID_GENERATOR.next_id();
    let read_id5 = READ_ID_GENERATOR.next_id();
    let read_id6 = READ_ID_GENERATOR.next_id();
    let read_id7 = READ_ID_GENERATOR.next_id();
    let region2 = new_region(CHR_2, 4000, 4200, &read_id4);
    region2.m().add_read_details(&read_id5, 4000, 4100, RemoteReadType::JunctionMate);
    assembly3.m().add_remote_regions(&[region2]);
    assembly8.m().add_junction_read(&create_read(&read_id4, 6000, assembly_bases, "50M"));
    assembly8.m().add_junction_read(&create_read(&read_id5, 6000, assembly_bases, "50M"));
    let region3 = new_region(CHR_1, 450, 600, &read_id6);
    region3.m().add_read_details(&read_id7, 450, 600, RemoteReadType::JunctionMate);
    assembly8.m().add_remote_regions(&[region3]);
    assembly2.m().add_junction_read(&create_read(&read_id6, 500, assembly_bases, "50M"));
    assembly2.m().add_junction_read(&create_read(&read_id7, 500, assembly_bases, "50M"));

    let (phase_group_set, _) = builder.build_remote_groups(&all_junction_groups, &mut writer);
    let phase_groups = phase_group_set.to_vec();
    assert_eq!(2, phase_groups.len());
    let phase_group = group_with(&phase_groups, assembly1);
    assert!(phase_group.assemblies().contains(&assembly4));
    let phase_group = group_with(&phase_groups, assembly3);
    assert_eq!(3, phase_group.assembly_count());
    assert!(phase_group.assemblies().contains(&assembly2));
    assert!(phase_group.assemblies().contains(&assembly8));
}

#[test]
fn test_phaset_set_merging() {
    let _env = test_env();
    // 2 pairs of links are merged
    let assembly_bases1a = format!("{}{}", &REF_BASES_400[121..201], &REF_BASES_400[250..300]);
    let assembly_bases1b = format!("{}{}", &REF_BASES_400[151..201], &REF_BASES_400[250..350]);
    let assembly1a = create_assembly(CHR_1, 200, FORWARD, &assembly_bases1a, 80);
    let assembly1b = create_assembly(CHR_1, 250, REVERSE, &assembly_bases1b, 50);
    let junc_read1a = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 121, &assembly_bases1a, "80M50S", CHR_1, 200, false);
    assembly1a.m().add_junction_read(&junc_read1a);
    let junc_read1b = create_read_mate(junc_read1a.id(), CHR_1, 250, &assembly_bases1b, "50S100M", CHR_1, 200, false);
    assembly1b.m().add_junction_read(&junc_read1b);
    let assembly_bases2a = format!("{}{}", &REF_BASES_400[101..201], &REF_BASES_400[250..300]); // longer at start
    let assembly_bases2b = format!("{}{}", &REF_BASES_400[151..201], &REF_BASES_400[250..330]); // shorter at end
    let assembly2a = create_assembly(CHR_2, 200, FORWARD, &assembly_bases2a, 99);
    let assembly2b = create_assembly(CHR_2, 250, REVERSE, &assembly_bases2b, 50);
    let junc_read2a = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_2, 121, &assembly_bases1a, "100M50S", CHR_2, 200, false);
    assembly2a.m().add_junction_read(&junc_read2a);
    let junc_read2b = create_read_mate(junc_read2a.id(), CHR_2, 250, &assembly_bases1b, "50S80M", CHR_2, 200, false);
    assembly2b.m().add_junction_read(&junc_read2b);
    let phase_group = PhaseGroup::new(assembly1a, Some(assembly1b));
    PhaseGroup::add_assembly(phase_group, assembly2a);
    PhaseGroup::add_assembly(phase_group, assembly2b);
    let ref_genome = MockRefGenome::new();
    build_phase_sets(&ref_genome, phase_group);
    assert_eq!(2, phase_group.phase_sets().len());
    let phase_set1 = phase_group.phase_sets()[0];
    let phase_set2 = phase_group.phase_sets()[1];
    PhaseGroup::finalise_phase_set_alignments(phase_group);
    assert!(!phase_set1.merged());
    assert!(phase_set1.merged_phase_sets().contains(&phase_set2));
    assert!(phase_set2.merged());
    assert_eq!(200, phase_set1.assembly_alignment().unwrap().full_sequence_length());

    // cannot merge if neither junction is overlapped - ie if just matching on ref bases
    let assembly_bases1a = format!("{}{}", &REF_BASES_400[1..51], &REF_BASES_400[100..150]);
    let assembly_bases1b = format!("{}{}", &REF_BASES_400[1..51], &REF_BASES_400[100..250]);
    let assembly1a = create_assembly(CHR_1, 50, FORWARD, &assembly_bases1a, 49);
    let assembly1b = create_assembly(CHR_1, 100, REVERSE, &assembly_bases1b, 50);
    let junc_read1a = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 1, &assembly_bases1a, "51M50S", CHR_1, 200, false);
    assembly1a.m().add_junction_read(&junc_read1a);
    let junc_read1b = create_read_mate(junc_read1a.id(), CHR_1, 100, &assembly_bases1b, "50S151M", CHR_1, 50, false);
    assembly1b.m().add_junction_read(&junc_read1b);
    let assembly_bases2a = format!("{}{}", &REF_BASES_400[120..250], &REF_BASES_400[300..350]);
    let assembly_bases2b = format!("{}{}", &REF_BASES_400[200..250], &REF_BASES_400[300..400]);
    let assembly2a = create_assembly(CHR_2, 250, FORWARD, &assembly_bases2a, 149);
    let assembly2b = create_assembly(CHR_2, 300, REVERSE, &assembly_bases2b, 50);
    let junc_read2a = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_2, 120, &assembly_bases1a, "131M50S", CHR_2, 200, false);
    assembly2a.m().add_junction_read(&junc_read2a);
    let junc_read2b = create_read_mate(junc_read2a.id(), CHR_2, 200, &assembly_bases1b, "51S80M", CHR_2, 200, false);
    assembly2b.m().add_junction_read(&junc_read2b);
    let phase_group = PhaseGroup::new(assembly1a, Some(assembly1b));
    PhaseGroup::add_assembly(phase_group, assembly2a);
    PhaseGroup::add_assembly(phase_group, assembly2b);
    build_phase_sets(&ref_genome, phase_group);
    assert_eq!(2, phase_group.phase_sets().len());
    let phase_set1 = phase_group.phase_sets()[0];
    let phase_set2 = phase_group.phase_sets()[1];
    PhaseGroup::finalise_phase_set_alignments(phase_group);
    assert!(!phase_set1.merged());
    assert!(!phase_set2.merged());
}
