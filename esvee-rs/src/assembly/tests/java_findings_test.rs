//! The corrected behaviour of the Java findings of REPORT.md section 4.2, each beside the Java behaviour
//! that `-java_compat` brings back. Synthetic input only.
use super::*;
use crate::assembly::alignment::breakend::Breakend;
use crate::assembly::constants::BASE_QUAL_MINIMUM;
use crate::assembly::globals::{proximate_junction_distance, set_proximate_junction_distance};
use crate::assembly::indel_builder::set_indel_status;
use crate::assembly::output::assembly_writer::repeats_info_str;
use crate::assembly::phase::phase_set_builder::{standard_comparator, ExtensionCandidate, ExtensionType};
use crate::assembly::remote_region_finder::filter_low_qual_read;
use crate::assembly::seq_tech_utils::{find_illumina_adapter_trims, set_seq_tech_specifics, trim_sbx_uncertain_bases, UNPAIRED_READ_PROXIMATE_JUNCTION_DISTANCE};
use crate::assembly::types::remote_region::RemoteRegion;
use crate::assembly::types::RemoteReadType;
use crate::assembly::types::repeat_info::RepeatInfo;
use crate::common::compat;
use crate::common::regions::SpecificRegions;
use std::cmp::Ordering as CmpOrdering;

/// both behaviours of one finding: (corrected, Java)
fn both<T>(f: impl Fn() -> T) -> (T, T) {
    (f(), compat::with_java(true, &f))
}

#[test]
fn a6_exact_support_column_is_read_without_the_other_support_column() {
    let _env = test_env();
    let path = std::env::temp_dir().join(format!("esvee_a6_{}.tsv", std::process::id()));
    // one junction fragment and three exact-support fragments; the file has no OtherSupportFrags column
    std::fs::write(&path, "Chromosome\tPosition\tOrientation\tJunctionFrags\tExactSupportFrags\n1\t1000\t1\t1\t3\n").unwrap();
    let load = || Junction::load_junctions(path.to_str().unwrap(), &SpecificRegions::default(), 3, 1, 1).unwrap().values().map(|x| x.len()).sum::<usize>();
    let (corrected, java) = both(load);
    std::fs::remove_file(&path).unwrap();
    assert_eq!(1, corrected);
    assert_eq!(0, java);
}

#[test]
fn a9_merged_region_gets_the_other_soft_clip_total_once() {
    let _env = test_env();
    let merged_total = || {
        let region = |start: i32, read_id: &str| {
            let region = P::new(RemoteRegion::new(intern(CHR_1), start, start + 100, read_id, RemoteReadType::Supplementary));
            region.m().add_soft_clip_map_qual(50, 2);
            region
        };
        let mut regions = vec![region(1000, "READ_A"), region(1050, "READ_B")];
        RemoteRegion::merge_regions(&mut regions);
        assert_eq!(1, regions.len());
        regions[0].soft_clip_map_qual_total()
    };
    assert_eq!((200, 400), both(merged_total));
}

#[test]
fn a10_sbx_start_trim_is_found_from_the_start() {
    let _env = test_env();
    set_sbx_sequencing();
    let remaining_bases = || {
        // uncertain bases: two at the start, three at the end; the end is trimmed first, then the start
        let read = create_read(TEST_READ_ID, 100, "ACGTACGTAC", "10M");
        let mut quals = build_default_base_quals(10);
        for i in [0, 1, 7, 8, 9] {
            quals[i] = BASE_QUAL_MINIMUM;
        }
        read.m().bam_record_mut().base_qualities = quals;
        assert!(trim_sbx_uncertain_bases(read.m()));
        (read.bases_length(), read.trim_count_start(), read.trim_count_end())
    };
    let (corrected, java) = both(remaining_bases);
    assert_eq!((5, 2, 3), corrected);
    assert_eq!((7, 0, 3), java);
}

#[test]
fn a11_proximate_junction_distance_is_kept_for_paired_reads() {
    let _env = test_env();
    let before = proximate_junction_distance();
    let distance = || {
        set_proximate_junction_distance(before);
        set_seq_tech_specifics();
        proximate_junction_distance()
    };
    let (corrected, java) = both(distance);
    set_sbx_sequencing();
    let unpaired = distance();
    set_proximate_junction_distance(before);
    assert_eq!(before, corrected);
    assert_eq!(UNPAIRED_READ_PROXIMATE_JUNCTION_DISTANCE, java);
    assert_eq!(UNPAIRED_READ_PROXIMATE_JUNCTION_DISTANCE, unpaired);
}

fn linked_pair(forward_start: i32, reverse_start: i32, cigar: &str) -> (ReadP, ReadP) {
    let read_bases = ref_bases_random_100();
    let read = create_read(TEST_READ_ID, forward_start, &read_bases, cigar);
    read.m().bam_record_mut().set_flag(FLAG_MATE_REVERSE_STRAND, true);
    let mate_read = create_read(TEST_READ_ID, reverse_start, &read_bases, cigar);
    mate_read.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, true);
    read.m().set_mate_read(&mate_read);
    mate_read.m().set_mate_read(&read);
    (read, mate_read)
}

fn unlink(pair: &(ReadP, ReadP)) {
    pair.0.m().clear_mate_read();
    pair.1.m().clear_mate_read();
}

#[test]
fn a12_adapter_trims_of_linked_mates() {
    let _env = test_env();
    // a fragment shorter than the reads: each read runs 10 bases past the 5' end of its mate
    let pair = linked_pair(110, 100, "20S60M20S");
    let reads = vec![pair.0.clone(), pair.1.clone()];
    assert_eq!(vec![(0, 10, false), (1, 10, true)], find_illumina_adapter_trims(&reads));
    // the same whichever read comes first
    let reads = vec![pair.1.clone(), pair.0.clone()];
    assert_eq!(vec![(0, 10, true), (1, 10, false)], find_illumina_adapter_trims(&reads));
    unlink(&pair);
    // without a link there is nothing to trim against
    assert!(find_illumina_adapter_trims(&reads).is_empty());

    // reads facing each other with the reverse read further on: nothing extends past the mate
    let pair = linked_pair(100, 300, "100M");
    assert!(find_illumina_adapter_trims(&[pair.0.clone(), pair.1.clone()]).is_empty());
    unlink(&pair);

    // reads facing away from each other (the reverse read wholly before the forward read): not an adapter
    // layout, the rule taken literally would remove every base of both
    let pair = linked_pair(1000, 500, "100M");
    assert!(find_illumina_adapter_trims(&[pair.0.clone(), pair.1.clone()]).is_empty());
    unlink(&pair);
}

#[test]
fn a13_left_clips_at_a_reverse_junction_count_as_soft_clips() {
    let _env = test_env();
    let indel_based = || {
        let junction = new_junction(CHR_1, 200, REVERSE);
        junction.m().set_indel_status(true);
        let read_bases = ref_bases_random_100();
        let reads: Vec<ReadP> = (0..4).map(|i| create_read(&format!("READ_{i}"), 200, &read_bases, "40S60M")).collect();
        set_indel_status(junction, &reads);
        junction.indel_based()
    };
    let (corrected, java) = both(indel_based);
    assert!(!corrected);
    assert!(java);
}

#[test]
fn a17_non_standard_bases_are_rewritten_over_the_whole_read() {
    let _env = test_env();
    let bases_after = || {
        // high quality throughout, so the low-quality count is settled after 11 of the 20 bases
        let mut record = create_sam_record(TEST_READ_ID, CHR_1, 100, "ACGRACGTACGTACGRACGT", "20M");
        assert!(!filter_low_qual_read(&mut record));
        String::from_utf8(record.read_bases.clone()).unwrap()
    };
    let (corrected, java) = both(bases_after);
    assert_eq!("ACGNACGTACGTACGNACGT", corrected);
    assert_eq!("ACGNACGTACGTACGRACGT", java);
}

#[test]
fn a18_breakends_are_ordered_by_the_chromosome_of_the_other_breakend() {
    let _env = test_env();
    let assembly = create_assembly(CHR_1, 200, FORWARD, &REF_BASES_400[100..200], 49);
    let assembly_alignment = AssemblyAlignment::from_assembly(assembly);
    let breakend_to = |other_chromosome: &str| {
        let breakend = Breakend::new(assembly_alignment, CHR_1, 200, FORWARD, "", None);
        let other = Breakend::new(assembly_alignment, other_chromosome, 500, REVERSE, "", None);
        breakend.m().set_other_breakend(Some(other));
        other.m().set_other_breakend(Some(breakend));
        breakend
    };
    let (to_chr2, to_chr3) = (breakend_to(CHR_2), breakend_to(CHR_3));
    let compare = || (to_chr2.compare_to(&to_chr3), to_chr3.compare_to(&to_chr2));
    let (corrected, java) = both(compare);
    assert_eq!((CmpOrdering::Less, CmpOrdering::Greater), corrected);
    assert_eq!((CmpOrdering::Equal, CmpOrdering::Equal), java);
}

#[test]
fn a19_the_more_precise_split_link_comes_first() {
    let _env = test_env();
    let assembly1 = create_assembly(CHR_1, 200, FORWARD, &REF_BASES_400[100..200], 49);
    let assembly2 = create_assembly(CHR_1, 1200, REVERSE, &REF_BASES_400[200..300], 50);
    let candidate = |inserted_bases: &str| {
        let candidate = ExtensionCandidate::from_link(ExtensionType::SplitLink, AssemblyLink::new(assembly1, assembly2, LinkType::Split, inserted_bases, ""));
        candidate.m().support_count = 5;
        candidate
    };
    let (precise, with_insert) = (candidate(""), candidate("ACGTACGT"));
    let compare = || (standard_comparator(&precise, &with_insert), standard_comparator(&with_insert, &precise));
    let (corrected, java) = both(compare);
    assert_eq!((CmpOrdering::Less, CmpOrdering::Greater), corrected);
    assert_eq!((CmpOrdering::Equal, CmpOrdering::Equal), java);
    // the support count still comes first
    with_insert.m().support_count = 6;
    assert_eq!(CmpOrdering::Greater, standard_comparator(&precise, &with_insert));
}

#[test]
fn a22_longest_repeat_is_compared_with_the_longest_so_far() {
    let repeats = vec![RepeatInfo::new(0, "A", 5), RepeatInfo::new(10, "ACGT", 3), RepeatInfo::new(30, "ACG", 2)];
    let (corrected, java) = both(|| repeats_info_str(&repeats));
    assert_eq!("3 max(A=5) long(ACGT=3)", corrected);
    assert_eq!("3 max(A=5) long(ACG=2)", java);
}
