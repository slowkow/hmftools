//! JunctionsTest

use super::util::*;
use crate::prep::depth::{DepthTracker, DEPTH_WINDOW_SIZE};
use crate::prep::rec::*;
use crate::prep::tracker::*;
use crate::prep::types::Orientation::{Forward, Reverse};
use crate::prep::types::ReadType::{CandidateSupport, Junction, NoSupport};
use crate::prep::types::*;

fn ref_bases() -> String {
    generate_random_bases(500)
}

fn tracker() -> (JunctionTracker, DepthTracker) {
    let region = Region::new(chr(CHR_1), 1, 5000);
    let depth = DepthTracker::new(region.start, region.end, DEPTH_WINDOW_SIZE);
    (new_tracker(region), depth)
}

fn add_read(tracker: &mut JunctionTracker, mut read: PrepRead, read_type: ReadType) -> Option<ReadId> {
    read.set_read_type(read_type);
    tracker.process_read(read)
}

fn junction_at(tracker: &JunctionTracker, position: i32) -> Option<&JunctionData> {
    tracker.junctions.iter().map(|&j| &tracker.jstore[j as usize]).find(|j| j.position == position)
}

#[test]
fn test_basic_junctions() {
    let ref_bases = ref_bases();
    let (mut t, _) = tracker();

    let id1 = next_id();
    let read1 = prep_read(sam_record(&id1, CHR_1, 800, &ref_bases[..100], "30S70M"));
    let read2 = prep_read(sam_record_flags(&id1, CHR_1, 820, &ref_bases[20..120], "100M", build_flags(false, true, false), DEFAULT_BASE_QUAL));
    add_read(&mut t, read1, Junction);
    add_read(&mut t, read2, NoSupport);

    let supp_read1 = prep_read(sam_record(&next_id(), CHR_1, 800, &ref_bases[..73], "3S70M"));
    assert!(has(&supp_read1.rec, FLAG_MATE_REVERSE));
    assert!(!has(&supp_read1.rec, FLAG_REVERSE));
    add_read(&mut t, supp_read1, CandidateSupport);

    let id3 = next_id();
    let read3 = prep_read(sam_record(&id3, CHR_1, 950, &ref_bases[..100], "30S70M"));
    let read4 = prep_read(sam_record_flags(&id3, CHR_1, 980, &ref_bases[20..120], "100M", build_flags(false, true, false), DEFAULT_BASE_QUAL));
    add_read(&mut t, read3, Junction);
    add_read(&mut t, read4, NoSupport);

    let supp_read2 = prep_read(sam_record(&next_id(), CHR_1, 950, &ref_bases[..73], "3S70M"));
    add_read(&mut t, supp_read2, CandidateSupport);

    let id5 = next_id();
    let read5 = prep_read(sam_record(&id5, CHR_1, 950, &ref_bases[20..120], "100M"));
    let read6 = prep_read(sam_record(&id5, CHR_1, 980, &ref_bases[..100], "70M30S"));
    add_read(&mut t, read5, NoSupport);
    add_read(&mut t, read6, Junction);

    let supp_read3 = prep_read(sam_record(&next_id(), CHR_1, 990, &ref_bases[..63], "60M3S"));
    add_read(&mut t, supp_read3, CandidateSupport);

    let id7 = next_id();
    let read7 = prep_read(sam_record(&id7, CHR_1, 1010, &ref_bases[10..90], "50M30S"));
    let read8 = prep_read(sam_record(&id7, CHR_1, 1010, &ref_bases[..50], "50M"));
    add_read(&mut t, read7, Junction);
    add_read(&mut t, read8, NoSupport);

    let supp_read4 = prep_read(sam_record(&next_id(), CHR_1, 990, &ref_bases[..73], "70M3S"));
    add_read(&mut t, supp_read4, CandidateSupport);

    t.assign_junction_fragments_and_support();

    assert_eq!(4, t.junctions.len());

    let j = junction_at(&t, 800).expect("junction");
    assert_eq!(Reverse, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(1, j.exact_support_fragment_count());

    let j = junction_at(&t, 950).expect("junction");
    assert_eq!(Reverse, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(1, j.exact_support_fragment_count());

    let j = junction_at(&t, 1049).expect("junction");
    assert_eq!(Forward, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(2, j.exact_support_fragment_count());
    assert_eq!(1, j.supporting_fragment_count());

    let j = junction_at(&t, 1059).expect("junction");
    assert_eq!(Forward, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(2, j.exact_support_fragment_count());
    assert_eq!(0, j.supporting_fragment_count());
}

#[test]
fn test_candidate_only_junctions() {
    let ref_bases = ref_bases();
    let (mut t, depth) = tracker();

    let mut read1 = prep_read(sam_record(&next_id(), CHR_1, 800, &ref_bases[..100], "30S70M"));
    read1.rec.set_mapq(0); // will cause the junction to be filtered
    add_read(&mut t, read1, Junction);

    let supp_read1 = prep_read(sam_record_mate(&next_id(), CHR_1, 810, CHR_3, 100, true, false, None));
    add_read(&mut t, supp_read1, CandidateSupport);

    let mut supp_read2 = prep_read(sam_record(&next_id(), CHR_1, 800, &ref_bases[..73], "3S70M"));
    supp_read2.rec.set_mapq(0);
    add_read(&mut t, supp_read2, CandidateSupport);

    let read3 = prep_read(sam_record_mate(&next_id(), CHR_1, 950, CHR_2, 100, true, false, None));
    add_read(&mut t, read3, CandidateSupport);

    t.assign_fragments(&depth);

    let junction_groups = t.form_unique_assigned_groups();
    assert_eq!(0, junction_groups.len());

    let remote_candidate_groups = t.get_remote_candidate_read_groups();
    assert_eq!(4, remote_candidate_groups.len());
}

#[test]
fn test_internal_deletes() {
    let ref_bases = ref_bases();
    let (mut t, _) = tracker();

    // initial delete is too short
    let read1 = prep_read(sam_record(&next_id(), CHR_1, 10, &ref_bases[..80], "20M10D50M"));
    add_read(&mut t, read1, Junction);

    // then a simple one
    let read2 = prep_read(sam_record(&next_id(), CHR_1, 100, &ref_bases[..80], "20M40D20M"));
    add_read(&mut t, read2, Junction);

    // with supporting reads - first is too short as an indel
    let supp_read = prep_read(sam_record(&next_id(), CHR_1, 100, &ref_bases[..80], "20M20D20M"));
    add_read(&mut t, supp_read, CandidateSupport);

    let supp_read = prep_read(sam_record(&next_id(), CHR_1, 120, &ref_bases[..80], "20M20D20M"));
    add_read(&mut t, supp_read, CandidateSupport);

    // and a more complicated one
    let read3 = prep_read(sam_record(&next_id(), CHR_1, 210, &ref_bases[..1], "5S10M2D10M3I10M35D10M2S"));
    add_read(&mut t, read3, Junction);

    t.assign_junction_fragments_and_support();

    assert_eq!(4, t.junctions.len());

    let j = junction_at(&t, 119).expect("junction");
    assert_eq!(Forward, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(1, j.exact_support_fragment_count());

    let j = junction_at(&t, 160).expect("junction");
    assert_eq!(Reverse, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(1, j.exact_support_fragment_count());

    let j = junction_at(&t, 241).expect("junction");
    assert_eq!(Forward, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(0, j.exact_support_fragment_count());

    let j = junction_at(&t, 277).expect("junction");
    assert_eq!(Reverse, j.orient);
    assert_eq!(1, j.junction_fragment_count());
    assert_eq!(0, j.exact_support_fragment_count());
}

#[test]
fn test_internal_inserts() {
    let ref_bases = ref_bases();
    let (mut t, _) = tracker();

    // first is too short
    let read1 = prep_read(sam_record(&next_id(), CHR_1, 10, &ref_bases[..70], "20M10I50M"));
    add_read(&mut t, read1, Junction);

    // then a simple one
    let read2 = prep_read(sam_record(&next_id(), CHR_1, 100, &ref_bases[..70], "20M40I50M"));
    add_read(&mut t, read2, Junction);

    // and a more complicated one
    let read3 = prep_read(sam_record(&next_id(), CHR_1, 210, &ref_bases[..100], "5S10M2D10M3I10M35I10M2S"));
    add_read(&mut t, read3, Junction);

    // not anchored with enough bases
    let read4a = prep_read(sam_record(&next_id(), CHR_1, 310, &ref_bases[..50], "9M40D41M"));
    add_read(&mut t, read4a, Junction);

    let read4b = prep_read(sam_record(&next_id(), CHR_1, 310, &ref_bases[..50], "41M40D9M"));
    add_read(&mut t, read4b, Junction);

    t.assign_junction_fragments_and_support();

    assert_eq!(4, t.junctions.len());
    assert_eq!(119, t.junction(0).position);
    assert_eq!(120, t.junction(1).position);
    assert_eq!(241, t.junction(2).position);
    assert_eq!(242, t.junction(3).position);
}

#[test]
fn test_short_templated_insertions() {
    let (mut t, depth) = tracker();
    let read1 = prep_read(sam_record(&next_id(), CHR_1, 200, &REF_BASES_400[..130], "32S66M32S"));
    let read2 = prep_read(sam_record(&next_id(), CHR_1, 200, &REF_BASES_400[..130], "32S66M32S"));
    add_read(&mut t, read1, Junction);
    add_read(&mut t, read2, Junction);

    t.assign_fragments(&depth);

    assert_eq!(2, t.junctions.len());

    let j = junction_at(&t, 200).expect("junction");
    assert_eq!(Reverse, j.orient);
    assert_eq!(2, j.junction_fragment_count());

    let j = junction_at(&t, 265).expect("junction");
    assert_eq!(Forward, j.orient);
    assert_eq!(2, j.junction_fragment_count());
}

/// The Java test puts single-read groups into a map and calls JunctionUtils.markSupplementaryDuplicates; here
/// the two reads go into a fresh tracker, which gives each its own group. Returns the statuses of the two.
fn duplicate_statuses(first: PrepRead, second: PrepRead, permitted_position_diff: i32) -> (ReadGroupStatus, ReadGroupStatus) {
    let (mut t, _) = tracker();
    let (id1, id2): (Vec<u8>, Vec<u8>) = (first.id().to_vec(), second.id().to_vec());
    t.process_read(first);
    t.process_read(second);
    t.mark_supplementary_duplicates(permitted_position_diff);
    let status = |id: &[u8]| t.groups[t.group_of(id).unwrap() as usize].status;
    (status(&id1), status(&id2))
}

#[test]
fn test_primary_supplementary_duplicates() {
    // primary and supplementary with matching coords and mates
    let ref_bases = ref_bases();
    let read_bases = &ref_bases[..50];
    let (lower_cigar, upper_cigar, mate_cigar) = ("20S30M", "20M30S", "50M");
    let supp_data_upper = Some((CHR_3, 1000, '+', upper_cigar, 60));
    let supp_data_lower = Some((CHR_1, 1000, '+', lower_cigar, 60));
    let (id1, id2) = (next_id(), next_id());

    let make = |consensus: bool| {
        let mut primary1 = prep_read(common_sam_record(&id1, CHR_1, 1000, read_bases, lower_cigar, CHR_2, 1000, false, false, supp_data_upper, true, Some(mate_cigar)));
        primary1.set_read_type(Junction);
        let mut supp1 = prep_read(common_sam_record(&id1, CHR_3, 1000, read_bases, upper_cigar, CHR_2, 1000, false, true, supp_data_lower, true, Some(mate_cigar)));
        supp1.set_read_type(CandidateSupport);
        let mut supp2_rec = common_sam_record(&id2, CHR_1, 1000, read_bases, lower_cigar, CHR_2, 1000, false, true, supp_data_upper, true, Some(mate_cigar));
        let mut primary2_rec = common_sam_record(&id2, CHR_3, 1000, read_bases, upper_cigar, CHR_2, 1000, false, false, supp_data_lower, true, Some(mate_cigar));
        if consensus {
            set_attribute(&mut supp2_rec, b"CR", Some(&encode_z(b"2:0")));
            set_attribute(&mut primary2_rec, b"CR", Some(&encode_z(b"2:0")));
        }
        let mut supp2 = prep_read(supp2_rec);
        supp2.set_read_type(Junction);
        let mut primary2 = prep_read(primary2_rec);
        primary2.set_read_type(CandidateSupport);
        (primary1, supp1, supp2, primary2)
    };

    let (primary1, supp1, supp2, primary2) = make(false);
    let (s1, s2) = duplicate_statuses(primary1, supp2, 0);
    assert_ne!(ReadGroupStatus::Duplicate, s1);
    assert_eq!(ReadGroupStatus::Duplicate, s2);

    // test that the upper reads find the same duplicate
    let (s1, s2) = duplicate_statuses(supp1, primary2, 0);
    assert_ne!(ReadGroupStatus::Duplicate, s1);
    assert_eq!(ReadGroupStatus::Duplicate, s2);

    // repeat again, checking that consensus reads are favoured over non-consensus
    let (primary1, supp1, supp2, primary2) = make(true);
    let (s1, s2) = duplicate_statuses(primary1, supp2, 0);
    assert_eq!(ReadGroupStatus::Duplicate, s1);
    assert_ne!(ReadGroupStatus::Duplicate, s2);

    let (s1, s2) = duplicate_statuses(supp1, primary2, 0);
    assert_eq!(ReadGroupStatus::Duplicate, s1);
    assert_ne!(ReadGroupStatus::Duplicate, s2);
}

#[test]
fn test_sbx_primary_supplementary_duplicates() {
    let ref_bases = ref_bases();
    let read_bases = &ref_bases[..50];
    let (lower_cigar1, upper_cigar1) = ("20S30M", "20M30S");
    let supp_data_upper = Some((CHR_3, 1000, '+', upper_cigar1, 60));
    let supp_data_lower = Some((CHR_1, 1001, '+', lower_cigar1, 60));
    let (id1, id2) = (next_id(), next_id());
    // as in Java the read is made as a paired one and the paired flag is then cleared on its record
    let unpaired = |mut read: PrepRead, read_type: ReadType| {
        set_flag(&mut read.rec, FLAG_PAIRED, false);
        read.set_read_type(read_type);
        read
    };

    let primary1 = unpaired(prep_read(common_sam_record(&id1, CHR_1, 1000, read_bases, lower_cigar1, NO_CHROMOSOME_NAME, 0, false, false, supp_data_upper, false, Some("*"))), Junction);
    let _supp1 = unpaired(prep_read(common_sam_record(&id1, CHR_3, 1000, read_bases, upper_cigar1, NO_CHROMOSOME_NAME, 0, false, true, supp_data_lower, false, Some("*"))), CandidateSupport);
    let supp2 = unpaired(prep_read(common_sam_record(&id2, CHR_1, 1001, read_bases, lower_cigar1, NO_CHROMOSOME_NAME, 0, false, true, supp_data_upper, false, Some("*"))), Junction);
    let _primary2 = unpaired(prep_read(common_sam_record(&id2, CHR_3, 1001, read_bases, upper_cigar1, NO_CHROMOSOME_NAME, 0, false, false, supp_data_lower, true, Some("*"))), CandidateSupport);

    let (s1, s2) = duplicate_statuses(primary1, supp2, 2);
    assert_ne!(ReadGroupStatus::Duplicate, s1);
    assert_eq!(ReadGroupStatus::Duplicate, s2);
}

#[test]
fn test_junction_data_create_and_lookup() {
    // test junction-finding logic
    let ref_bases = ref_bases();
    let read_bases = &ref_bases[..100];
    let left_cigar = "30S70M";
    let (mut t, _) = tracker();
    let mut add_junction = |t: &mut JunctionTracker, start: i32, cigar: &str| {
        let rid = t.add_read(prep_read(sam_record(&next_id(), CHR_1, start, read_bases, cigar)));
        t.add_junction_data(rid);
        rid
    };

    let rid = add_junction(&mut t, 100, left_cigar);
    // ensure junction matches prevent duplicate junctions
    t.add_junction_data(rid);
    add_junction(&mut t, 150, left_cigar);
    add_junction(&mut t, 200, left_cigar);
    add_junction(&mut t, 250, left_cigar);
    let rid = add_junction(&mut t, 300, left_cigar);
    t.add_junction_data(rid);

    assert_eq!(5, t.junctions.len());

    // add junctions at matching positions but opposite orientations
    let right_cigar = "50M30S";
    add_junction(&mut t, 101, right_cigar);
    assert_eq!(6, t.junctions.len());

    let rid = add_junction(&mut t, 201, right_cigar);
    assert_eq!(7, t.junctions.len());
    t.add_junction_data(rid);

    let find = |t: &JunctionTracker, position: i32, simple: usize| find_junction_index(t.junctions.len(), |i| t.junction(i).position, position, simple);

    assert_eq!(INVALID_JUNC_INDEX, find(&t, 50, 10));
    assert_eq!(0, find(&t, 100, 10));

    // using fast-search logic
    assert_eq!(0, find(&t, 100, 0));
    assert_eq!(0, find(&t, 149, 0));
    assert_eq!(1, find(&t, 150, 0));
    assert_eq!(2, find(&t, 151, 0));
    assert_eq!(3, find(&t, 200, 0));
    assert_eq!(3, find(&t, 225, 0));
    assert_eq!(3, find(&t, 249, 0));
    assert_eq!(4, find(&t, 250, 0));
    assert_eq!(5, find(&t, 299, 0));
    assert_eq!(6, find(&t, 300, 0));
    assert_eq!(6, find(&t, 301, 0));

    // test adding discordant junctions - they will replace a split junction if have sufficiently more support
    // first matches an existing, does not replace since the split junction has sufficient support
    let junc_index = find(&t, 200, 0) as usize;
    let split_junction = t.junctions[junc_index];
    for _ in 0..2 {
        let g = new_group(&mut t, prep_read(sam_record(&next_id(), CHR_1, 201, read_bases, right_cigar)));
        t.jstore[split_junction as usize].junction_groups.push(g);
    }

    let mut add_discordant = |t: &mut JunctionTracker, read_start: i32, position: i32| {
        let rid = t.add_read(prep_read(sam_record(&next_id(), CHR_1, read_start, read_bases, "50M")));
        let mut disc = JunctionData::new(position, Reverse, Some(rid));
        disc.discordant_group = true;
        t.jstore.push(disc);
        let jid = (t.jstore.len() - 1) as JuncId;
        t.add_discordant_junction(jid);
        jid
    };

    add_discordant(&mut t, 151, 200);
    assert_eq!(7, t.junctions.len());
    assert_eq!(split_junction, t.junctions[junc_index]);

    // check replacing
    let junc_index = find(&t, 300, 0) as usize;
    let disc = add_discordant(&mut t, 251, 300);
    assert_eq!(7, t.junctions.len());
    assert_eq!(disc, t.junctions[junc_index]);

    // add a new discordant group
    let disc = add_discordant(&mut t, 221, 270);
    assert_eq!(8, t.junctions.len());
    let junc_index = find(&t, 270, 0) as usize;
    assert_eq!(disc, t.junctions[junc_index]);
}
