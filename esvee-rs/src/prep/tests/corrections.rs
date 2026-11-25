//! The corrections of the Java findings in prep (REPORT.md section 4.1, and the prep site of A17): each case
//! is checked in corrected mode and, under `compat::with_java(true, ..)`, with the Java behaviour.

use super::util::*;
use crate::common::compat;
use crate::prep::config::{PrepConfig, WriteType};
use crate::prep::depth::{DepthTracker, DEPTH_WINDOW_SIZE};
use crate::prep::filters::*;
use crate::prep::rec::*;
use crate::prep::spanning::SpanningReadCache;
use crate::prep::tracker::JunctionTracker;
use crate::prep::types::Orientation::{Forward, Reverse};
use crate::prep::types::*;
use crate::prep::writer::{merge_junction_lines, DeferredJunction, ResultsWriter};
use rust_htslib::bam::{self, Read};

// P3: the last partial depth window

#[test]
fn p3_last_partial_window_has_depth() {
    let ref_bases = generate_random_bases(100);
    let depth_at = |position: i32| {
        let mut depth = DepthTracker::new(1, 1050, 100);
        for start in [951, 1001, 1001] {
            depth.process_read(&sam_record(&next_id(), CHR_1, start, &ref_bases[..50], "50M"));
        }
        depth.calc_depth(position)
    };
    // two reads cover each base of the 50-base window that ends the region
    assert_eq!(2.0, depth_at(1001));
    assert_eq!(2.0, depth_at(1050));
    assert_eq!(0.5, depth_at(1000));
    assert_eq!(0.0, depth_at(1051));
    // Java (P3): the partial window does not exist
    compat::with_java(true, || {
        assert_eq!(0.0, depth_at(1001));
        assert_eq!(0.5, depth_at(1000));
    });
}

#[test]
fn p3_read_past_the_region_end_is_cut_off() {
    let ref_bases = generate_random_bases(100);
    let mut depth = DepthTracker::new(1, 1050, 100);
    depth.process_read(&sam_record(&next_id(), CHR_1, 1021, &ref_bases[..100], "100M"));
    assert_eq!(30.0 / 50.0, depth.calc_depth(1040));
}

// P4: the name of a partition

#[test]
fn p4_last_base_belongs_to_its_partition() {
    let c = chr(CHR_1);
    assert_eq!(0, PartKey::form(c, 1, PARTITION_SIZE).part);
    assert_eq!(0, PartKey::form(c, PARTITION_SIZE, PARTITION_SIZE).part);
    assert_eq!(1, PartKey::form(c, PARTITION_SIZE + 1, PARTITION_SIZE).part);
    assert_eq!(1, PartKey::form(c, 2 * PARTITION_SIZE, PARTITION_SIZE).part);
    // Java (P4): the last base is named as the next partition
    compat::with_java(true, || {
        assert_eq!(1, PartKey::form(c, PARTITION_SIZE, PARTITION_SIZE).part);
        assert_eq!(1, PartKey::form(c, PARTITION_SIZE + 1, PARTITION_SIZE).part);
    });
}

#[test]
fn p4_mate_on_last_base_of_a_partition_is_expected_there() {
    // partitions 1-10000 and 10001-20000: a mate at 10000 is in the first
    let contigs = contigs();
    let region2 = Region::new(chr(CHR_1), PARTITION_SIZE + 1, 2 * PARTITION_SIZE);
    let remote_of = || {
        let mut tracker = new_tracker(region2.clone());
        let read = prep_read(sam_record_mate(&next_id(), CHR_1, 10800, CHR_1, PARTITION_SIZE, true, false, None));
        let g = new_group(&mut tracker, read);
        group_set_partition_count(&mut tracker.groups[g as usize], &tracker.reads, &region2, PARTITION_SIZE, &contigs);
        tracker.groups[g as usize].remote_partitions.clone()
    };
    assert_eq!(vec![PartKey { chr: chr(CHR_1), part: 0 }], remote_of());
    compat::with_java(true, || assert_eq!(vec![PartKey { chr: chr(CHR_1), part: 1 }], remote_of()));
}

#[test]
fn p4_junction_fragment_is_expected_where_its_mate_is() {
    // partitions 1-10000, 10001-20000, 20001-30000: a junction read in the first with its mate on the last
    // base of the second
    let contigs = contigs();
    let c = chr(CHR_1);
    let regions = [Region::new(c, 1, PARTITION_SIZE), Region::new(c, PARTITION_SIZE + 1, 2 * PARTITION_SIZE), Region::new(c, 2 * PARTITION_SIZE + 1, 3 * PARTITION_SIZE)];
    let expected_in = |region: usize| {
        let mut cache = SpanningReadCache::new(&PrepConfig::for_tests(PARTITION_SIZE), contigs.clone(), None);
        let mut tracker = new_tracker(regions[0].clone());
        let mut read = prep_read(sam_record_mate("FRAGMENT", CHR_1, 800, CHR_1, 2 * PARTITION_SIZE, true, false, None));
        read.set_read_type(ReadType::Junction);
        let g = new_group(&mut tracker, read);
        group_set_partition_count(&mut tracker.groups[g as usize], &tracker.reads, &regions[0], PARTITION_SIZE, &contigs);
        cache.process_spanning_read_groups(&regions[0], &[(g, false, true)], &mut tracker).unwrap();
        cache.get_expected_read_ids(&regions[region]).len()
    };
    assert_eq!((1, 0), (expected_in(1), expected_in(2)));
    // Java (P4): the partition after the mate's expects the fragment
    compat::with_java(true, || assert_eq!((0, 1), (expected_in(1), expected_in(2))));
}

// P10: short supplementary alignments

fn supplementary(bases: &str, cigar: &str, map_quality: u8) -> bam::Record {
    let mut rec = sam_record_flags(&next_id(), CHR_1, 1000, bases, cigar, build_flags(true, false, true), DEFAULT_BASE_QUAL);
    rec.set_mapq(map_quality);
    rec
}

#[test]
fn p10_short_supplementary_is_judged_by_its_score() {
    let mixed = &REF_BASES_400[..100];
    let repeat = format!("{}{}", &REF_BASES_400[..70], "A".repeat(30));
    // 30 aligned bases without repeats, good mapping quality: kept
    assert!(!filter_poor_mapped_supplementary(&supplementary(mixed, "70S30M", 60)));
    // the same with low mapping quality: dropped
    assert!(filter_poor_mapped_supplementary(&supplementary(mixed, "70S30M", 5)));
    // 30 aligned bases of one letter: dropped whatever the mapping quality
    assert!(filter_poor_mapped_supplementary(&supplementary(&repeat, "70S30M", 60)));
    // 40 aligned bases and more are judged as before
    assert!(!filter_poor_mapped_supplementary(&supplementary(mixed, "50S50M", 60)));
    assert!(!filter_poor_mapped_supplementary(&supplementary(mixed, "50S50M", 5)));
    // Java (P10): under 40 aligned bases always dropped
    compat::with_java(true, || {
        assert!(filter_poor_mapped_supplementary(&supplementary(mixed, "70S30M", 60)));
        assert!(filter_poor_mapped_supplementary(&supplementary(mixed, "70S30M", 5)));
        assert!(!filter_poor_mapped_supplementary(&supplementary(mixed, "50S50M", 60)));
    });
}

// A17: bases that are not A, C, G, T or N

#[test]
fn a17_all_non_standard_bases_become_n() {
    let mut bases = generate_random_bases(100).into_bytes();
    for i in [3, 60, 99] {
        bases[i] = b'R';
    }
    let bases = String::from_utf8(bases).unwrap();
    let config = ReadFilterConfig::default();
    let after = |low_qual_tail: bool| {
        let mut rec = sam_record(&next_id(), CHR_1, 1000, &bases, "100M");
        if low_qual_tail {
            let cigar = bam::record::CigarString::try_from("100M").unwrap();
            let mut quals = vec![DEFAULT_BASE_QUAL; 100];
            quals[40..].fill(10);
            rec.set(b"READ_LQ", Some(&cigar), bases.as_bytes(), &quals);
        }
        let filtered = filter_low_qual_read(&mut rec, &config);
        (filtered, seq_bytes(&rec))
    };
    // high qualities throughout: the count stops after 51 bases
    let (filtered, seq) = after(false);
    assert!(!filtered);
    assert_eq!((b'N', b'N', b'N'), (seq[3], seq[60], seq[99]));
    assert_eq!(3, seq.iter().filter(|&&b| b == b'N').count());
    // low qualities from base 40 on: filtered, and rewritten all the same
    let (filtered, seq) = after(true);
    assert!(filtered);
    assert_eq!((b'N', b'N', b'N'), (seq[3], seq[60], seq[99]));
    // Java (A17): only as far as the count read
    compat::with_java(true, || {
        let (filtered, seq) = after(false);
        assert!(!filtered);
        assert_eq!((b'N', b'R', b'R'), (seq[3], seq[60], seq[99]));
    });
}

// P8: a discordant junction at a remote region outside the partition

/// Five fragments near the end of the partition 1-5000 with their mates past its end
fn tracker_with_remote_discordant_group(mate_chr: &str, mate_start: i32) -> JunctionTracker {
    let region = Region::new(chr(CHR_1), 1, 5000);
    let hotspots = vec![std::sync::Arc::new(crate::prep::hotspot::KnownHotspot {
        region_start: Region::new(chr(CHR_1), 4000, 5000),
        orient_start: Forward,
        region_end: Region::new(chr(CHR_2), 7000, 8000),
        orient_end: Reverse,
        gene_info: String::new(),
    })];
    let mut tracker = JunctionTracker::new(region.clone(), tracker_config(), hotspots, None);
    for i in 0..5 {
        let mut read = prep_read(sam_record_mate(&next_id(), CHR_1, 4401 + i, mate_chr, mate_start + i, true, false, None));
        set_flag(&mut read.rec, FLAG_MATE_REVERSE, true);
        read.rec.set_insert_size(if mate_chr == CHR_1 { (mate_start - 4401) as i64 } else { 0 });
        let mut read = prep_read(std::mem::replace(&mut read.rec, bam::Record::new()));
        read.set_read_type(ReadType::CandidateSupport);
        tracker.process_read(read);
    }
    tracker.assign_fragments(&DepthTracker::new(region.start, region.end, DEPTH_WINDOW_SIZE));
    tracker
}

#[test]
fn p8_remote_junction_outside_partition_is_kept_apart() {
    for (mate_chr, mate_start) in [(CHR_1, 7401), (CHR_2, 7401)] {
        let tracker = tracker_with_remote_discordant_group(mate_chr, mate_start);
        assert_eq!(1, tracker.junctions.len());
        assert_eq!(4504, tracker.junction(0).position);
        assert_eq!(1, tracker.outside_junctions.len());
        let remote = &tracker.jstore[tracker.outside_junctions[0] as usize];
        assert_eq!((7401, Reverse, Some(chr(mate_chr))), (remote.position, remote.orient, remote.outside_chr));
        assert_eq!(5, remote.supporting_groups.len());

        // Java (P8): the remote junction is in the partition's own list, to be written under its chromosome
        compat::with_java(true, || {
            let tracker = tracker_with_remote_discordant_group(mate_chr, mate_start);
            assert_eq!(2, tracker.junctions.len());
            assert!(tracker.outside_junctions.is_empty());
            assert_eq!(7401, tracker.junction(1).position);
        });
    }
}

// P13: the allele-fraction filter of discordant junctions

/// The five fragments of `tracker_with_remote_discordant_group` (mates 3 kb on), with `depth` reads of 200
/// bases over the 200-base window of their junction
fn tracker_with_discordant_group_at_depth(depth: i32) -> JunctionTracker {
    let region = Region::new(chr(CHR_1), 1, 5000);
    let mut tracker = JunctionTracker::new(region.clone(), tracker_config(), Vec::new(), None);
    for i in 0..5 {
        let mut read = prep_read(sam_record_mate(&next_id(), CHR_1, 4401 + i, CHR_1, 7401 + i, true, false, None));
        set_flag(&mut read.rec, FLAG_MATE_REVERSE, true);
        read.rec.set_insert_size(3000);
        let mut read = prep_read(std::mem::replace(&mut read.rec, bam::Record::new()));
        read.set_read_type(ReadType::CandidateSupport);
        tracker.process_read(read);
    }
    let mut depth_tracker = DepthTracker::new(region.start, region.end, DEPTH_WINDOW_SIZE);
    for _ in 0..depth {
        depth_tracker.process_range(4401, 4600);
    }
    tracker.assign_fragments(&depth_tracker);
    tracker
}

#[test]
fn p13_discordant_group_below_one_percent_of_depth_makes_no_junction() {
    // five fragments are 1% of a depth of 500: kept, and the junction has its depth
    let tracker = tracker_with_discordant_group_at_depth(500);
    assert_eq!(1, tracker.junctions.len());
    assert_eq!((4504, 500), (tracker.junction(0).position, tracker.junction(0).depth));
    assert_eq!(5, tracker.junction(0).supporting_groups.len());
    // below 1% of 501: no junction here, nor at the remote region
    let tracker = tracker_with_discordant_group_at_depth(501);
    assert!(tracker.junctions.is_empty());
    assert!(tracker.outside_junctions.is_empty());

    // Java (P13): the junction is made at any depth and its depth stays 0
    compat::with_java(true, || {
        let tracker = tracker_with_discordant_group_at_depth(5000);
        assert_eq!(2, tracker.junctions.len());
        assert_eq!((4504, 0), (tracker.junction(0).position, tracker.junction(0).depth));
    });
}

#[test]
fn p8_groups_of_outside_junctions_are_still_written() {
    let mut tracker = tracker_with_remote_discordant_group(CHR_2, 7401);
    assert_eq!(5, tracker.form_unique_assigned_groups().len());
}

#[test]
fn p8_deferred_rows_go_to_their_place() {
    let row = |chromosome: &str, position: i32, orient: i32, tag: &str| format!("{chromosome}\t{position}\t{orient}\t0\t0\t3\t{tag}\n").into_bytes();
    let deferred = |chromosome: &str, position: i32, orient: i32, tag: &str| DeferredJunction {
        rank: chromosome_rank(chromosome),
        chromosome: chromosome.to_string(),
        position,
        orient,
        text: row(chromosome, position, orient, tag),
    };
    let file = vec![b"Chromosome\tPosition\tOrientation".to_vec(), row("chr1", 100, 1, "a"), row("chr1", 900, -1, "b"), row("chr3", 50, 1, "c"), row("chrX", 10, 1, "d")];
    let rows = vec![
        deferred("chr3", 60, -1, "late"),
        deferred("chr1", 500, 1, "mid"),
        deferred("chr2", 7, 1, "new chromosome"),
        deferred("chr1", 900, 1, "same position, other orientation"),
        deferred("chr1", 900, -1, "site already in the file"),
        deferred("chr2", 7, 1, "site deferred twice"),
        deferred("chrY", 5, 1, "last"),
    ];
    let mut out: Vec<u8> = Vec::new();
    merge_junction_lines(file.into_iter(), rows, &mut out).unwrap();
    let tags: Vec<String> = String::from_utf8(out).unwrap().lines().skip(1).map(|l| l.rsplit('\t').next().unwrap().to_string()).collect();
    assert_eq!(vec!["a", "mid", "b", "same position, other orientation", "new chromosome", "c", "late", "d", "last"], tags);
}

// P5 and P6: candidate reads parked in the cache BAMs

struct CacheRun {
    dir: std::path::PathBuf,
    config: PrepConfig,
    header_text: String,
}

impl CacheRun {
    fn new(name: &str, use_cache_bam: bool) -> CacheRun {
        let dir = std::env::temp_dir().join(format!("esvee_prep_test_{}_{}", std::process::id(), name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let mut config = PrepConfig::for_tests(PARTITION_SIZE);
        config.sample_ids = vec!["SAMPLE".to_string()];
        config.output_dir = format!("{}/", dir.display());
        config.write_types = vec![WriteType::PrepBam];
        config.use_cache_bam = use_cache_bam;
        let contigs = contigs();
        let mut header_text = String::from("@HD\tVN:1.6\n");
        for name in &contigs.names {
            header_text.push_str(&format!("@SQ\tSN:{name}\tLN:1000000\n"));
        }
        CacheRun { dir, config, header_text }
    }

    /// Partition 1 holds candidate reads only; partition 2 then has a supporting read of each fragment.
    /// Returns the reads written, as (read id, flags), sorted.
    fn run(&self) -> Vec<(String, u16)> {
        let contigs = contigs();
        let c = chr(CHR_1);
        let regions = [Region::new(c, 1, PARTITION_SIZE), Region::new(c, PARTITION_SIZE + 1, 2 * PARTITION_SIZE)];
        let header = bam::HeaderView::from_bytes(self.header_text.as_bytes());
        let mut cache = SpanningReadCache::new(&self.config, contigs.clone(), Some(header));
        let mut writer = ResultsWriter::new(&self.config, contigs.clone(), std::slice::from_ref(&self.header_text)).unwrap();
        let ref_bases = &REF_BASES_400[..100];

        // partition 1
        let mut tracker = new_tracker(regions[0].clone());
        let mut spanning: Vec<(GroupId, bool, bool)> = Vec::new();
        let mut park = |tracker: &mut JunctionTracker, reads: Vec<bam::Record>| {
            let mut gid = None;
            for rec in reads {
                let mut read = prep_read(rec);
                check_filters(&mut read, &ReadFilterConfig::default());
                if !read.is_unmapped() {
                    read.set_read_type(ReadType::CandidateSupport);
                }
                match gid {
                    None => gid = Some(new_group(tracker, read)),
                    Some(g) => {
                        add_to_group(tracker, g, read);
                    }
                }
            }
            let g = gid.unwrap();
            group_set_partition_count(&mut tracker.groups[g as usize], &tracker.reads, &regions[0], PARTITION_SIZE, &contigs);
            assert!(tracker.groups[g as usize].spans_partitions());
            spanning.push((g, false, false));
        };
        // a plain candidate read with its mate in partition 2
        park(&mut tracker, vec![sam_record_mate("PLAIN", CHR_1, 800, CHR_1, 10800, true, false, None)]);
        // a candidate read with a poly-G soft clip
        let poly_g = format!("{}{}", "G".repeat(20), &ref_bases[..80]);
        let mut rec = sam_record_flags("POLYG", CHR_1, 900, &poly_g, "20S80M", build_flags(true, false, false), DEFAULT_BASE_QUAL);
        rec.set_mpos(10900 - 1);
        park(&mut tracker, vec![rec]);
        // a candidate read with a supplementary alignment in partition 2, and its unmapped mate
        let mut mapped = sam_record_mate("UNMAPPED_MATE", CHR_1, 1000, CHR_1, 1000, true, false, Some("1,11000,+,50S50M,60,0"));
        let flags = (mapped.flags() | FLAG_MATE_UNMAPPED) & !FLAG_PROPER_PAIR;
        mapped.set_flags(flags);
        let mut unmapped = bam::Record::new();
        unmapped.set(b"UNMAPPED_MATE", None, ref_bases.as_bytes(), &[DEFAULT_BASE_QUAL; 100]);
        unmapped.set_tid(c);
        unmapped.set_pos(999);
        unmapped.set_mtid(c);
        unmapped.set_mpos(999);
        unmapped.set_flags(FLAG_PAIRED | FLAG_UNMAPPED | FLAG_SECOND);
        park(&mut tracker, vec![mapped, unmapped]);
        cache.process_spanning_read_groups(&regions[0], &spanning, &mut tracker).unwrap();

        // partition 2
        let mut tracker = new_tracker(regions[1].clone());
        let mut spanning: Vec<(GroupId, bool, bool)> = Vec::new();
        let mut junction_groups: Vec<GroupId> = Vec::new();
        for (id, start, mate_start) in [("PLAIN", 10800, 800), ("POLYG", 10900, 900), ("UNMAPPED_MATE", 11000, 1000)] {
            let mut read = prep_read(sam_record_mate(id, CHR_1, start, CHR_1, mate_start, false, id == "UNMAPPED_MATE", None));
            read.set_read_type(ReadType::Support);
            let g = new_group(&mut tracker, read);
            group_set_partition_count(&mut tracker.groups[g as usize], &tracker.reads, &regions[1], PARTITION_SIZE, &contigs);
            tracker.groups[g as usize].has_remote_junction_reads = true;
            spanning.push((g, false, true));
            junction_groups.push(g);
        }
        cache.process_spanning_read_groups(&regions[1], &spanning, &mut tracker).unwrap();
        writer.write_read_groups(&mut tracker, &junction_groups).unwrap();

        if self.config.use_cache_bam {
            cache.reset();
            cache.candidate_bam_writer.assign_candidate_reads(&mut writer, &tid_map(&contigs)).unwrap();
        }
        writer.close().unwrap();

        let mut reader = bam::Reader::from_path(self.config.form_filename(WriteType::UnsortedBam, "SAMPLE")).unwrap();
        let mut written: Vec<(String, u16)> = reader.records().map(|r| r.unwrap()).map(|r| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags())).collect();
        written.sort();
        written
    }
}

impl Drop for CacheRun {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn names(written: &[(String, u16)], unmapped: bool) -> Vec<&str> {
    written.iter().filter(|(_, f)| (f & FLAG_UNMAPPED != 0) == unmapped).map(|(n, _)| n.as_str()).collect()
}

#[test]
fn p5_p6_cache_bam_claims_what_the_memory_cache_claims() {
    let from_bam = CacheRun::new("cache_bam", true).run();
    let from_memory = CacheRun::new("memory", false).run();
    assert_eq!(from_memory, from_bam);
    // both reads of the plain fragment, the claiming read of the poly-G fragment but not the read with the
    // soft clip (P5), and all three reads of the third fragment, the unmapped one among them (P6)
    assert_eq!(vec!["PLAIN", "PLAIN", "POLYG", "UNMAPPED_MATE", "UNMAPPED_MATE"], names(&from_bam, false));
    assert_eq!(vec!["UNMAPPED_MATE"], names(&from_bam, true));
}

#[test]
fn p5_p6_java_cache_bam() {
    compat::with_java(true, || {
        // Java (P5): the read with the poly-G soft clip comes back from the cache BAM without its filters and
        // is written; (P6): the unmapped read is never claimed
        let from_bam = CacheRun::new("java_cache_bam", true).run();
        assert_eq!(vec!["PLAIN", "PLAIN", "POLYG", "POLYG", "UNMAPPED_MATE", "UNMAPPED_MATE"], names(&from_bam, false));
        assert!(names(&from_bam, true).is_empty());
        // the in-memory mode has neither fault
        let from_memory = CacheRun::new("java_memory", false).run();
        assert_eq!(vec!["PLAIN", "PLAIN", "POLYG", "UNMAPPED_MATE", "UNMAPPED_MATE"], names(&from_memory, false));
        assert_eq!(vec!["UNMAPPED_MATE"], names(&from_memory, true));
    });
}
