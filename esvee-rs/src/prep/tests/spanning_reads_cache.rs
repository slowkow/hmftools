//! SpanningReadsCacheTest. The Java test hands the cache a map of read groups; here the groups and their reads
//! sit in a tracker, as they do in a run, and the cache is given the groups to process.

use super::util::*;
use crate::prep::config::PrepConfig;
use crate::prep::spanning::SpanningReadCache;
use crate::prep::tracker::JunctionTracker;
use crate::prep::types::ReadType::{CandidateSupport, Junction, NoSupport, Support};
use crate::prep::types::*;

struct Fixture {
    cache: SpanningReadCache,
    store: JunctionTracker,
    contigs: std::sync::Arc<Contigs>,
    regions: [Region; 3],
}

impl Fixture {
    fn new() -> Fixture {
        let contigs = contigs();
        let c = chr(CHR_1);
        let region1 = Region::new(c, 1, PARTITION_SIZE - 1);
        let region2 = Region::new(c, region1.end + 1, region1.end + PARTITION_SIZE);
        let region3 = Region::new(c, region2.end + 1, region2.end + PARTITION_SIZE);
        Fixture {
            cache: SpanningReadCache::new(&PrepConfig::for_tests(PARTITION_SIZE), contigs.clone(), None),
            store: new_tracker(region1.clone()),
            contigs,
            regions: [region1, region2, region3],
        }
    }
    /// A read of 100M with its mate elsewhere, of the given type
    #[allow(clippy::too_many_arguments)]
    fn read(&self, read_id: i32, start: i32, mate_start: i32, first_in_pair: bool, is_supp: bool, supp_data: &str, read_type: ReadType) -> PrepRead {
        let mut read = prep_read(sam_record_mate(&read_id_str(read_id), CHR_1, start, CHR_1, mate_start, first_in_pair, is_supp, Some(supp_data)));
        read.set_read_type(read_type);
        read
    }
    fn group(&mut self, read: PrepRead) -> GroupId {
        new_group(&mut self.store, read)
    }
    /// setPartitionCount, returning partitionCount()
    fn set_partition_count(&mut self, gid: GroupId, region: usize) -> usize {
        let (groups, reads) = (&mut self.store.groups, &self.store.reads);
        group_set_partition_count(&mut groups[gid as usize], reads, &self.regions[region], PARTITION_SIZE, &self.contigs);
        groups[gid as usize].remote_partitions.len() + 1
    }
    fn process(&mut self, region: usize, groups: &[GroupId]) {
        let spanning: Vec<(GroupId, bool, bool)> = groups.iter().map(|&g| (g, false, true)).collect();
        self.cache.process_spanning_read_groups(&self.regions[region], &spanning, &mut self.store).unwrap();
    }
    fn cached_reads_count(&self) -> usize {
        self.cache.get_cached_reads_count(None)
    }
    fn has_junction_read_id(&self, gid: GroupId) -> bool {
        let id = &self.store.groups[gid as usize].id;
        self.cache.junction_partition_read_ids().values().any(|x| x.contains(id))
    }
    fn group_size(&self, gid: GroupId) -> usize {
        self.store.groups[gid as usize].reads.len()
    }
    /// Whether the group holds a read that starts at the position (the Java test asks for the read object).
    fn group_has_read_at(&self, gid: GroupId, start: i32) -> bool {
        self.store.groups[gid as usize].reads.iter().any(|&r| self.store.reads[r as usize].start == start)
    }
}

// The regions of this test are 1-9999, 10000-19999, ...: they fit the partition names of Java, which
// names a partition position / size although partitions start at 1 + index * size (finding P4).
#[test]
fn test_basic_reads() {
    crate::common::compat::with_java(true, test_basic_reads_body);
}

fn test_basic_reads_body() {
    let mut f = Fixture::new();
    let mut read_id = 0;

    // 2 junction reads in different partitions
    read_id += 1;
    let read1 = f.read(read_id, 800, 10800, true, false, "", Junction);
    let g = f.group(read1);
    assert_eq!(2, f.set_partition_count(g, 0));

    let read2 = f.read(read_id, 10800, 800, false, false, "", Junction);
    let g = f.group(read2);
    assert_eq!(2, f.set_partition_count(g, 1));

    f.process(0, &[g]);
    assert_eq!(0, f.cached_reads_count());
    assert!(!f.has_junction_read_id(g));

    f.process(1, &[g]);
    assert!(!f.has_junction_read_id(g));
    assert_eq!(0, f.cached_reads_count());

    f.cache.reset();

    // a junction read then a candidate read
    read_id += 1;
    let read1 = f.read(read_id, 800, 10800, true, false, "", Junction);
    let g = f.group(read1);
    assert_eq!(2, f.set_partition_count(g, 0));
    f.process(0, &[g]);

    let read2 = f.read(read_id, 10800, 800, false, false, "", CandidateSupport);
    let g = f.group(read2);
    f.set_partition_count(g, 1);
    f.process(1, &[g]);
    assert_eq!(0, f.cached_reads_count());
    assert!(!f.has_junction_read_id(g));

    f.cache.reset();

    // a candidate read then a supporting read which claims it
    read_id += 1;
    let read1 = f.read(read_id, 800, 10800, true, false, "", CandidateSupport);
    let g = f.group(read1);
    assert_eq!(2, f.set_partition_count(g, 0));
    f.process(0, &[g]);
    assert_eq!(1, f.cached_reads_count());
    assert!(!f.has_junction_read_id(g));

    let read2 = f.read(read_id, 10800, 800, false, false, "", Support);
    let g = f.group(read2);
    f.set_partition_count(g, 1);
    f.process(1, &[g]);
    assert_eq!(0, f.cached_reads_count());
    assert!(!f.has_junction_read_id(g));
    assert!(f.group_has_read_at(g, 800));

    f.cache.reset();

    // a candidate read whose remote partition holds nothing for it
    read_id += 1;
    let read1 = f.read(read_id, 800, 10800, true, false, "", CandidateSupport);
    let g = f.group(read1);
    f.set_partition_count(g, 0);
    f.process(0, &[g]);
    assert_eq!(1, f.cached_reads_count());

    f.process(1, &[]);
    assert_eq!(0, f.cached_reads_count());
}

// The regions of this test are 1-9999, 10000-19999, ...: they fit the partition names of Java, which
// names a partition position / size although partitions start at 1 + index * size (finding P4).
#[test]
fn test_combined_groups_supplementary() {
    crate::common::compat::with_java(true, test_combined_groups_supplementary_body);
}

fn test_combined_groups_supplementary_body() {
    let mut f = Fixture::new();
    let mut read_id = 0;

    read_id += 1;
    let read1 = f.read(read_id, 800, 801, true, false, "", NoSupport);
    let read2 = f.read(read_id, 801, 800, false, false, "1;10800;-;46S30M;255;0", Junction);
    let g = f.group(read1);
    let read2_id = add_to_group(&mut f.store, g, read2);
    assert_eq!(2, f.set_partition_count(g, 0));
    f.process(0, &[g]);
    assert!(f.has_junction_read_id(g));
    assert_eq!(0, f.cached_reads_count());

    let read3 = f.read(read_id, 10800, 800, false, true, "1;801;-;46S30M;255;0", NoSupport);
    f.store.reads[read2_id as usize].set_read_type(CandidateSupport);
    let g = f.group(read3);
    assert_eq!(2, f.set_partition_count(g, 1));
    f.process(1, &[g]);
    assert_eq!(0, f.cached_reads_count());
    assert!(!f.has_junction_read_id(g));

    f.cache.reset();

    read_id += 1;
    let read1 = f.read(read_id, 800, 10801, true, false, "", CandidateSupport);
    let g = f.group(read1);
    assert_eq!(2, f.set_partition_count(g, 0));
    f.process(0, &[g]);
    assert_eq!(1, f.cached_reads_count());

    let read2 = f.read(read_id, 10801, 800, false, false, "1;20800;-;46S30M;255;0", Junction);
    let g = f.group(read2);
    assert_eq!(3, f.set_partition_count(g, 1));
    f.process(1, &[g]);
    assert_eq!(0, f.cached_reads_count());
    assert!(f.has_junction_read_id(g)); // keeps the junction id for expected reads
    assert!(f.group_has_read_at(g, 800)); // has picked up the candidate

    let read3 = f.read(read_id, 20800, 800, false, true, "1;10801;-;46S30M;255;0", CandidateSupport);
    let g = f.group(read3);
    assert_eq!(3, f.set_partition_count(g, 2));
    f.process(2, &[g]);
    assert_eq!(0, f.cached_reads_count());
    assert!(!f.has_junction_read_id(g));

    // twice, as in the Java test
    for _ in 0..2 {
        f.cache.reset();

        read_id += 1;
        let read1 = f.read(read_id, 800, 10801, false, false, "", CandidateSupport);
        let read2 = f.read(read_id, 801, 800, true, true, "1;10801;-;46S30M;255;0", Junction);
        let g = f.group(read1);
        add_to_group(&mut f.store, g, read2);
        assert_eq!(2, f.set_partition_count(g, 0));
        f.process(0, &[g]);
        assert_eq!(0, f.cached_reads_count());
        assert!(f.has_junction_read_id(g)); // doesn't know if the following read will be a junction

        let read3 = f.read(read_id, 10801, 800, true, false, "1;801;-;46S30M;255;0", Junction);
        let g = f.group(read3);
        assert_eq!(2, f.set_partition_count(g, 1));
        f.process(1, &[g]);
        assert_eq!(0, f.cached_reads_count());
        assert!(!f.has_junction_read_id(g));
        assert_eq!(1, f.group_size(g));
    }
}

// The regions of this test are 1-9999, 10000-19999, ...: they fit the partition names of Java, which
// names a partition position / size although partitions start at 1 + index * size (finding P4).
#[test]
fn test_candidate_only_groups() {
    crate::common::compat::with_java(true, test_candidate_only_groups_body);
}

fn test_candidate_only_groups_body() {
    let mut f = Fixture::new();
    let read_id = 1;

    let read1 = f.read(read_id, 800, 10800, true, false, "", CandidateSupport);
    let rg1 = f.group(read1);
    assert_eq!(2, f.set_partition_count(rg1, 0));
    assert!(group_conditional_on_remote_reads(&f.store.groups[rg1 as usize], &f.store.reads));
    f.process(0, &[rg1]);
    assert_eq!(1, f.cached_reads_count());

    let read2 = f.read(read_id, 10800, 800, false, false, "", CandidateSupport);
    let rg2 = f.group(read2);
    assert_eq!(2, f.set_partition_count(rg2, 1));
    assert!(group_conditional_on_remote_reads(&f.store.groups[rg2 as usize], &f.store.reads));
    f.process(1, &[rg2]);
    assert_eq!(0, f.cached_reads_count());
    assert!(!f.store.groups[rg2 as usize].has_remote_junction_reads);
    assert_eq!(1, f.group_size(rg2));

    f.cache.reset();

    let read3 = f.read(read_id, 10800, 11300, true, false, "", CandidateSupport);
    let mate3 = f.read(read_id, 11300, 10800, false, false, "", NoSupport);
    let rg3 = f.group(read3);
    add_to_group(&mut f.store, rg3, mate3);
    assert_eq!(1, f.set_partition_count(rg3, 1));
    assert!(group_conditional_on_remote_reads(&f.store.groups[rg3 as usize], &f.store.reads));
}
