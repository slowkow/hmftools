//! UnmappedExtensionTest
use super::*;
use crate::assembly::phase::phase_set_builder::PhaseSetBuilder;
use crate::assembly::phase::remote_read_extractor::RemoteReadExtractor;
use crate::assembly::types::phase_group::PhaseGroup;

fn unmapped_read(bases: &str, mate_reversed: bool) -> ReadP {
    let read = create_read_mate(&READ_ID_GENERATOR.next_id(), CHR_1, 1, bases, NO_CIGAR, CHR_1, 1, mate_reversed);
    read.m().bam_record_mut().set_flag(FLAG_READ_UNMAPPED, true);
    read
}

fn build_phase_sets(ref_genome: &MockRefGenome, assembly: AsmP) {
    let phase_group = PhaseGroup::new(assembly, None);
    let mut remote_read_extractor = RemoteReadExtractor::new(None, None);
    let mut phase_set_builder = PhaseSetBuilder::new(Some(ref_genome), &mut remote_read_extractor, phase_group);
    phase_set_builder.build_phase_sets();
}

#[test]
fn test_unmapped_read_extension() {
    let _env = test_env();
    let mut ref_genome = MockRefGenome::new_indexing(true);
    ref_genome.put(CHR_1, REF_BASES_400);
    let ref_bases = ref_genome.base_string(CHR_1, 1, 101);
    let ext_bases = &REF_BASES_400[200..250];
    let assembly_bases = format!("{ref_bases}{ext_bases}");
    let assembly = create_assembly(CHR_1, 100, FORWARD, &assembly_bases, ref_bases.len() as i32 - 1);

    // a set of unmapped reads to consider for extension
    let unmapped_read1 = unmapped_read(&REF_BASES_400[200..300], true);
    // the next read needs its bases reversed, and is only added once the first one has extended the extension bases
    let unmapped_read2 = unmapped_read(&reverse_complement_str(&REF_BASES_400[250..350]), false);
    // the next 3 reads require testing for consensus bases, and the first has too many mismatches to be added
    let unmapped_read_bases3 = &REF_BASES_400[300..400];
    let unmapped_read_bases3a = format!("{}AAAAA{}", &unmapped_read_bases3[0..60], &unmapped_read_bases3[65..]);
    let unmapped_read3a = unmapped_read(&unmapped_read_bases3a, true);
    // minor mismatches for the others
    let unmapped_read_bases3b = format!("{}{}{}", &unmapped_read_bases3[0..70], get_next_base(&unmapped_read_bases3[70..71]), &unmapped_read_bases3[71..]);
    let unmapped_read3b = unmapped_read(&unmapped_read_bases3b, true);
    let unmapped_read3c = unmapped_read(&reverse_complement_str(unmapped_read_bases3), false);
    for read in [&unmapped_read1, &unmapped_read2, &unmapped_read3a, &unmapped_read3b, &unmapped_read3c] {
        assembly.m().unmapped_reads_mut().push(read.clone());
    }
    build_phase_sets(&ref_genome, assembly);

    assert_eq!(200, assembly.extension_length());
    assert_eq!(4, assembly.support_count());
    assert!(!assembly.has_read_support(Some(&unmapped_read3a)));
}

#[test]
fn test_unmapped_read_extension_neg_junction() {
    let _env = test_env();
    let mut ref_genome = MockRefGenome::new_indexing(true);
    ref_genome.put(CHR_1, REF_BASES_400);
    let ref_bases = ref_genome.base_string(CHR_1, 300, 400);
    let ext_bases = &REF_BASES_400[150..200];
    let assembly_bases = format!("{ext_bases}{ref_bases}");
    let assembly = create_assembly(CHR_1, 300, REVERSE, &assembly_bases, ext_bases.len() as i32);

    // reads with varying levels of mismatches
    let unmapped_read1 = unmapped_read(&REF_BASES_400[100..200], true);
    let unmapped_read2 = unmapped_read(&reverse_complement_str(&REF_BASES_400[100..200]), false);
    let unmapped_read_bases3 = &REF_BASES_400[100..200];
    let unmapped_read_bases3a = format!("{}AAAAA{}", &unmapped_read_bases3[0..30], &unmapped_read_bases3[35..]);
    let unmapped_read3a = unmapped_read(&unmapped_read_bases3a, true);
    let unmapped_read_bases3b = format!("{}{}{}", &unmapped_read_bases3[0..10], get_next_base(&unmapped_read_bases3[10..11]), &unmapped_read_bases3[11..]);
    let unmapped_read3b = unmapped_read(&unmapped_read_bases3b, true);
    let unmapped_read_bases3c = format!("{}{}{}", &unmapped_read_bases3[0..40], get_next_base(&unmapped_read_bases3[40..41]), &unmapped_read_bases3[41..]);
    let unmapped_read3c = unmapped_read(&unmapped_read_bases3c, true);
    for read in [&unmapped_read1, &unmapped_read2, &unmapped_read3a, &unmapped_read3b, &unmapped_read3c] {
        assembly.m().unmapped_reads_mut().push(read.clone());
    }
    build_phase_sets(&ref_genome, assembly);

    assert_eq!(100, assembly.extension_length());
    assert_eq!(4, assembly.support_count());
    assert!(!assembly.has_read_support(Some(&unmapped_read3a)));
}
