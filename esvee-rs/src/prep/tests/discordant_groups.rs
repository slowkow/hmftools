//! DiscordantGroupsTest

use super::util::*;
use crate::prep::discordant::*;
use crate::prep::filters::DEFAULT_MAX_FRAGMENT_LENGTH;
use crate::prep::hotspot::KnownHotspot;
use crate::prep::rec::*;
use crate::prep::types::Orientation::{Forward, Reverse};
use crate::prep::types::*;
use std::sync::Arc;

fn region1() -> Region {
    Region::new(chr(CHR_1), 1, PARTITION_SIZE - 1)
}

/// Reads and their single-read groups, as the Java lists of ReadGroup
#[derive(Default)]
struct Store {
    reads: Vec<PrepRead>,
    groups: Vec<ReadGroup>,
}

impl Store {
    fn group(&mut self, read: PrepRead) -> GroupId {
        let id: Box<[u8]> = read.id().into();
        self.reads.push(read);
        self.groups.push(ReadGroup::new(id, (self.reads.len() - 1) as ReadId));
        (self.groups.len() - 1) as GroupId
    }
    fn add_discordant_candidate(&mut self, candidates: &mut Vec<GroupId>, chr1: &str, pos1: i32, chr2: &str, pos2: i32) {
        let mut read = prep_read(sam_record_mate(&next_id(), chr1, pos1, chr2, pos2, true, false, None));
        read.set_read_type(ReadType::CandidateSupport);
        set_flag(&mut read.rec, FLAG_MATE_REVERSE, true);
        candidates.push(self.group(read));
    }
    fn form(&mut self, finder: &DiscordantGroups, candidates: &[GroupId]) -> Vec<JunctionData> {
        let mut junctions = Vec::new();
        finder.form_discordant_junctions(candidates, &mut self.groups, &mut self.reads, &mut junctions);
        junctions
    }
}

#[test]
fn test_discordant_read_criteria() {
    let known_hotspots = vec![Arc::new(KnownHotspot {
        region_start: Region::new(chr(CHR_1), 100, 1000),
        orient_start: Forward,
        region_end: Region::new(chr(CHR_2), 100, 1000),
        orient_end: Reverse,
        gene_info: String::new(),
    })];
    let discordant_groups = DiscordantGroups::new(region1(), DEFAULT_MAX_FRAGMENT_LENGTH, known_hotspots, contigs());
    let mut store = Store::default();

    let read = prep_read(sam_record_mate(&next_id(), CHR_1, 500, CHR_2, 500, true, false, None));
    let g = store.group(read);
    assert!(discordant_groups.is_discordant_group(&store.groups[g as usize], &store.reads));
    assert!(discordant_groups.is_relevant_discordant_group(&store.groups[g as usize], &store.reads));

    // outside the known pair ranges
    let read = prep_read(sam_record_mate(&next_id(), CHR_1, 500, CHR_2, 10000, true, false, None));
    let g = store.group(read);
    assert!(discordant_groups.is_discordant_group(&store.groups[g as usize], &store.reads));
    assert!(!discordant_groups.is_relevant_discordant_group(&store.groups[g as usize], &store.reads));

    // otherwise local within the required distance
    let mut read = prep_read(sam_record_mate(&next_id(), CHR_1, 500, CHR_1, 10000, true, false, None));
    read.rec.set_insert_size(9500);
    let g = store.group(read);
    assert!(discordant_groups.is_discordant_group(&store.groups[g as usize], &store.reads));
    assert!(discordant_groups.is_relevant_discordant_group(&store.groups[g as usize], &store.reads));

    let mut read = prep_read(sam_record_mate(&next_id(), CHR_1, 500, CHR_1, DISCORDANT_GROUP_MAX_LOCAL_LENGTH * 2, true, false, None));
    let insert_size = (mate_start(&read.rec) - read.start) as i64;
    read.rec.set_insert_size(insert_size);
    let g = store.group(read);
    assert!(discordant_groups.is_discordant_group(&store.groups[g as usize], &store.reads));
    assert!(!discordant_groups.is_relevant_discordant_group(&store.groups[g as usize], &store.reads));
}

#[test]
fn test_discordant_groups() {
    // 3 fragments are required to support a discordant junction, unassigned to other junctions
    let mut store = Store::default();
    let mut candidates: Vec<GroupId> = Vec::new();
    store.add_discordant_candidate(&mut candidates, CHR_1, 401, CHR_2, 5100);
    store.add_discordant_candidate(&mut candidates, CHR_1, 421, CHR_3, 5000); // unrelated
    store.add_discordant_candidate(&mut candidates, CHR_1, 431, CHR_2, 5150);

    let discordant_groups = DiscordantGroups::new(region1(), DEFAULT_MAX_FRAGMENT_LENGTH, Vec::new(), contigs());

    let junctions = store.form(&discordant_groups, &candidates);
    assert_eq!(0, junctions.len());

    store.add_discordant_candidate(&mut candidates, CHR_1, 441, CHR_2, 5050);

    let junctions = store.form(&discordant_groups, &candidates);
    assert_eq!(2, junctions.len());
    assert!(junctions[0].junction_groups.is_empty());
    assert_eq!(540, junctions[0].position);
    assert_eq!(Forward, junctions[0].orient);
    assert_eq!(5050, junctions[1].position);
    assert_eq!(Reverse, junctions[1].orient);

    store.add_discordant_candidate(&mut candidates, CHR_1, 540, CHR_2, 105000); // unrelated
    store.add_discordant_candidate(&mut candidates, CHR_1, 550, CHR_2, 5300);
    store.add_discordant_candidate(&mut candidates, CHR_1, 560, CHR_2, 6000); // too far
    store.add_discordant_candidate(&mut candidates, CHR_1, 570, CHR_2, 5000);

    let junctions = store.form(&discordant_groups, &candidates);
    assert_eq!(2, junctions.len());
    assert_eq!(8, junctions[0].supporting_groups.len());

    // a local DEL needs more support
    candidates.clear();
    for _ in 0..3 {
        store.add_discordant_candidate(&mut candidates, CHR_1, 401, CHR_1, 1200);
    }

    let junctions = store.form(&discordant_groups, &candidates);
    assert_eq!(0, junctions.len());

    for _ in 0..2 {
        store.add_discordant_candidate(&mut candidates, CHR_1, 401, CHR_1, 1200);
    }

    let junctions = store.form(&discordant_groups, &candidates);
    assert_eq!(2, junctions.len());
    assert_eq!(500, junctions[0].position);
    assert_eq!(Forward, junctions[0].orient);
    assert_eq!(1200, junctions[1].position);
    assert_eq!(Reverse, junctions[1].orient);
}

#[test]
fn test_discordant_group_distance_filters() {
    // fragments cannot be too far from the junction
    let mut store = Store::default();
    let mut candidates: Vec<GroupId> = Vec::new();
    // the last three: two whose remote region has too little support for a junction of its own, and one
    // too far from the junction
    for (start, mate_chr, mate_start) in [(2000, CHR_2, 500), (2090, CHR_2, 500), (2180, CHR_2, 500), (2270, CHR_2, 10000), (2360, CHR_2, 10000), (2450, CHR_3, 500)] {
        let mut read = prep_read(sam_record_mate(&next_id(), CHR_1, start, mate_chr, mate_start, true, false, None));
        set_flag(&mut read.rec, FLAG_REVERSE, true);
        set_flag(&mut read.rec, FLAG_MATE_REVERSE, false);
        candidates.push(store.group(read));
    }

    let discordant_groups = DiscordantGroups::new(region1(), 500, Vec::new(), contigs());
    let junctions = store.form(&discordant_groups, &candidates);
    assert_eq!(2, junctions.len());
    assert_eq!(5, junctions[0].supporting_groups.len());
}
