//! Tests of the ISOFOX 2.0.1 mode (version.rs): each behaviour that differs from 2.1, in both modes where it
//! applies. All data is synthetic.

use crate::allocator::FragmentAllocator;
use crate::allocator201::{map_quality_weight, slice_filter_201};
use crate::chimeric::{is_realigned_fragment_candidate_201, read_split_junction_201};
use crate::common::*;
use crate::config::{GeneRegionFilters, IsofoxConfig, IsofoxFunction, WriteType};
use crate::fusion_frag::{FusionRead, FusionReadGroup};
use crate::gene::{GeneCollection, GeneReadData, RegionReadData};
use crate::read::{parse_cigar, Read, FLAG_DUPLICATE, FLAG_FIRST, FLAG_PAIRED, FLAG_PROPER, FLAG_REVERSE, FLAG_SECOND, FLAG_SECONDARY};
use crate::read_trans::process_overlapping_regions;
use crate::refgen::RefSource;
use crate::tests::{gen_bases, gene, read, region, transcript};
use crate::version::with_v201;
use hmf_common::regions::BaseRegion;
use std::collections::{HashMap, HashSet};

#[test]
fn map_quality_weights() {
    assert_eq!(map_quality_weight(255), 1.0);
    assert_eq!(map_quality_weight(4), 1.0);
    assert_eq!(map_quality_weight(3), 0.5);
    assert_eq!(map_quality_weight(2), 0.33);
    assert_eq!(map_quality_weight(1), 0.2);
    assert_eq!(map_quality_weight(0), 0.1);
}

#[test]
fn multi_mapped_by_map_quality() {
    let mut r = read(1, "1", 100, "20M", None);
    r.mapq = 3;
    assert!(with_v201(true, || r.is_multi_mapped()));
    assert!(!with_v201(false, || r.is_multi_mapped()));
    r.mapq = 255;
    assert!(!with_v201(true, || r.is_multi_mapped()));
}

#[test]
fn improper_pair_is_chimeric_in_201() {
    let mut r = read(1, "1", 100, "20M", None);
    r.flags = FLAG_PAIRED | crate::read::FLAG_MATE_REVERSE; // a forward read, its mate reverse, not proper
    assert!(with_v201(true, || r.is_chimeric()));
    assert!(!with_v201(false, || r.is_chimeric()));
}

#[test]
fn adapter_trimming_201() {
    // forward read with 10 soft-clipped bases past the reverse mate's 5' end
    let mut f = read(1, "1", 100, "40M10S", None);
    let mut r = read(1, "1", 100, "40M", None);
    r.flags |= FLAG_REVERSE;
    let (o, s, e) = (r.orient(), r.unclipped_start, r.unclipped_end);
    f.trim_adapter_soft_clip_bases_201(o, s, e);
    assert_eq!(f.cigar_str_201(), "40M");
    assert_eq!(f.base_length(), 40);
    assert_eq!(f.unclipped_end, 139);
    // partly trimmed
    let mut f = read(1, "1", 100, "40M10S", None);
    f.trim_adapter_soft_clip_bases_201(o, s, e + 4);
    assert_eq!(f.cigar_str_201(), "40M4S");
    assert_eq!(f.base_length(), 44);
    // the reverse read trims its 5' soft clip to the forward mate's start
    let mut r2 = read(1, "1", 105, "5S40M", None);
    r2.flags |= FLAG_REVERSE;
    let f2 = read(1, "1", 105, "40M", None);
    r2.trim_adapter_soft_clip_bases_201(f2.orient(), f2.unclipped_start, f2.unclipped_end);
    assert_eq!(r2.cigar_str_201(), "40M");
    assert_eq!(r2.unclipped_start, 105);
    // no soft clip, the mate ending inside: Java removes the last element, changing only the cigar string
    // (IF18); corrected, nothing changes
    let mut f3 = read(1, "1", 100, "10M100N40M", None);
    crate::compat::with_java(true, || f3.trim_adapter_soft_clip_bases_201(o, s, 200));
    assert_eq!(f3.cigar_str_201(), "10M100N");
    assert_eq!(f3.base_length(), 50);
    assert_eq!(f3.unclipped_end, f3.pos_end);
    let mut f3 = read(1, "1", 100, "10M100N40M", None);
    crate::compat::with_java(false, || f3.trim_adapter_soft_clip_bases_201(o, s, 200));
    assert_eq!(f3.cigar_str_201(), "10M100N40M");
    // same orientation: nothing
    let mut f4 = read(1, "1", 100, "40M10S", None);
    f4.trim_adapter_soft_clip_bases_201(ORIENT_FWD, 100, 120);
    assert_eq!(f4.cigar_str_201(), "40M10S");
}

#[test]
fn split_junction_from_mapped_sections() {
    let r = read(1, "1", 100, "10M100N10M", None);
    assert_eq!(read_split_junction_201(&r), Some([109, 210]));
    assert_eq!(crate::fusion_frag::read_split_junction(&r), Some([109, 210]));
    // the longest N
    let r = read(1, "1", 100, "10M50N10M200N10M", None);
    assert_eq!(read_split_junction_201(&r), Some([169, 370]));
    // an N beside a deletion: the gap is wider than the N, so 2.0.1 finds no junction (IF14); 2.1 does
    let r = read(1, "1", 100, "10M5D100N10M", None);
    assert_eq!(read_split_junction_201(&r), None);
    assert_eq!(crate::fusion_frag::read_split_junction(&r), Some([114, 215]));
}

#[test]
fn realignment_candidates_201() {
    let r = read(1, "1", 100, "3S30M", None);
    assert!(is_realigned_fragment_candidate_201(&r));
    let r = read(1, "1", 100, "30M10S", None);
    assert!(is_realigned_fragment_candidate_201(&r));
    let r = read(1, "1", 100, "30M11S", None);
    assert!(!is_realigned_fragment_candidate_201(&r));
    assert!(crate::chimeric::is_realigned_fragment_candidate(&r));
    let r = read(1, "1", 100, "2S30M", None);
    assert!(!is_realigned_fragment_candidate_201(&r));
}

/// TransClassificationTest's soft-clip cases: 2.1 leaves a short terminal clip within the exon EXONIC,
/// 2.0.1 makes any unmatched clip ALT
#[test]
fn short_soft_clip_within_exon() {
    let ref_bases = "ACGT".repeat(5);
    let mut r1 = region("GENE1", 1, 1, 100, 119);
    r1.set_ref_bases(ref_bases.clone().into_bytes());
    let mut r2 = region("GENE1", 1, 2, 300, 319);
    r2.set_ref_bases(ref_bases.clone().into_bytes());
    r1.post_regions.push(1);
    r2.pre_regions.push(0);
    let regions = vec![r1, r2];
    let classify = |pos: i32, cigar: &str, bases: &str| -> TransMatchType {
        let mut r = read(1, "1", pos, cigar, Some(bases));
        r.insert_size = 200;
        let overlapping = crate::read_trans::find_overlapping_regions(&regions, None, &r);
        process_overlapping_regions(&mut r, &regions, &overlapping);
        r.trans_classification(1)
    };
    let b1 = format!("{}A", &ref_bases[..18]);
    assert_eq!(with_v201(false, || classify(300, "18M1S", &b1)), TransMatchType::Exonic);
    assert_eq!(with_v201(true, || classify(300, "18M1S", &b1)), TransMatchType::Alt);
    let b2 = format!("AA{}", &ref_bases[..18]);
    assert_eq!(with_v201(false, || classify(102, "2S18M", &b2)), TransMatchType::Exonic);
    assert_eq!(with_v201(true, || classify(102, "2S18M", &b2)), TransMatchType::Alt);
}

#[test]
fn filters_201() {
    let mut f = GeneRegionFilters::default();
    f.excluded_regions.insert("chr1".into(), vec![BaseRegion::new(1000, 2000)]);
    // excluded regions apply in 2.0.1 only (2.1 has none)
    assert!(with_v201(true, || f.skip_read("chr1", 1500, 1600, false)));
    assert!(with_v201(true, || f.skip_read("chr1", 1900, 2100, true)));
    assert!(!with_v201(true, || f.skip_read("chr1", 2001, 2100, false)));
    assert_eq!(f.find_excluded_regions("chr1", 500, 1000).len(), 1);
    assert!(f.find_excluded_regions("chr1", 2001, 3000).is_empty());
    assert!(with_v201(true, || f.skip_read("chrUn_x", 1, 10, false)));
}

#[test]
fn slice_filters_201() {
    let mut config = IsofoxConfig::for_tests();
    config.single_map_quality = 255;
    config.functions = vec![IsofoxFunction::TranscriptCounts, IsofoxFunction::Fusions];
    let f = slice_filter_201(&config);
    assert!(f.keep_duplicates && f.keep_secondaries && f.keep_supplementaries);
    assert_eq!(f.min_map_quality, 0);
    config.drop_duplicates = true;
    assert!(!slice_filter_201(&config).keep_duplicates);
    // without transcript counts: no secondaries or duplicates, and only uniquely mapped reads
    config.functions = vec![IsofoxFunction::Fusions];
    let f = slice_filter_201(&config);
    assert!(!f.keep_duplicates && !f.keep_secondaries && f.keep_supplementaries);
    assert_eq!(f.min_map_quality, 255);
}

#[test]
fn write_types_by_version() {
    assert!(with_v201(true, || WriteType::value_of("CHIMERIC_READ")).is_ok());
    assert!(with_v201(true, || WriteType::value_of("FUSION_FRAGMENT")).is_err());
    assert!(with_v201(false, || WriteType::value_of("CHIMERIC_READ")).is_err());
    assert!(with_v201(false, || WriteType::value_of("FUSION_FRAGMENT")).is_ok());
    assert!(with_v201(true, || IsofoxFunction::value_of("READ_COUNTS")).is_ok());
    assert!(with_v201(false, || IsofoxFunction::value_of("READ_COUNTS")).is_err());
}

#[test]
fn fusion_group_completeness_201() {
    let regions: Vec<RegionReadData> = Vec::new();
    let mut a = read(1, "1", 100, "30M", None);
    a.flags |= FLAG_FIRST;
    let mut b = read(1, "1", 300, "30M", None);
    b.flags |= FLAG_SECOND;
    let fa = FusionRead::new(&a, &regions);
    let fb = FusionRead::new(&b, &regions);
    // a pair without supplementary data is complete in both
    assert!(with_v201(true, || FusionReadGroup::new("1", vec![fa.clone(), fb.clone()]).complete));
    assert!(with_v201(false, || FusionReadGroup::new("1", vec![fa.clone(), fb.clone()]).complete));
    // a lone read whose mate is unmapped: complete in 2.1, not in Java 2.0.1
    let mut c = read(2, "1", 100, "30M", None);
    c.flags |= crate::read::FLAG_MATE_UNMAPPED;
    let fc = FusionRead::new(&c, &regions);
    let java201 = |f: &dyn Fn() -> bool| crate::compat::with_java(true, || with_v201(true, f));
    assert!(with_v201(false, || FusionReadGroup::new("2", vec![fc.clone()]).complete));
    assert!(!java201(&|| FusionReadGroup::new("2", vec![fc.clone()]).complete));
    // any three reads are complete in Java 2.0.1
    assert!(java201(&|| FusionReadGroup::new("1", vec![fa.clone(), fb.clone(), fc.clone()]).complete));
    // IF16: a fragment with two supplementary alignments (four reads) never completes in Java 2.0.1
    let sa = Some("1,5000,+,30S30M,255,0");
    let mk = |pos: i32, flags: u16| {
        let c = parse_cigar("30M30S").unwrap();
        FusionRead::new(&Read::new("3", "1", pos, FLAG_PAIRED | flags, c, gen_bases(60).into_bytes(), "1", 100, 0, 255, sa, None, None), &regions)
    };
    let four = vec![
        mk(100, FLAG_FIRST),
        mk(5000, FLAG_FIRST | crate::read::FLAG_SUPPLEMENTARY),
        mk(300, FLAG_SECOND),
        mk(5200, FLAG_SECOND | crate::read::FLAG_SUPPLEMENTARY),
    ];
    assert!(!java201(&|| FusionReadGroup::new("3", four.clone()).complete));
    assert!(crate::compat::with_java(false, || with_v201(true, || FusionReadGroup::new("3", four.clone()).complete)));
}

#[test]
fn retained_intron_file_201() {
    let (f, h) = with_v201(true, crate::novel::retained_intron_file_and_header);
    assert_eq!(f, "retained_intron.csv");
    assert!(h.starts_with("GeneId,GeneName,Chromosome,Strand,Position,Type"));
    let (f, h) = with_v201(false, crate::novel::retained_intron_file_and_header);
    assert_eq!(f, "retained_intron.tsv");
    assert!(h.contains('\t'));
    assert!(with_v201(true, crate::task::gene_collection_header).contains("\tLowMapQual\t"));
}

fn test_gc() -> GeneCollection {
    let g = gene("GENE1", "GENE1", "1", 1, 1000, 5000);
    let t = transcript(1, "TRANS01", "GENE1", 1, &[(1000, 1200), (2000, 2500), (4500, 5000)], true);
    let gd = GeneReadData { gene: g, exon_regions: Vec::new(), transcripts: vec![t], has_unspliced_regions: false };
    GeneCollection::new(0, vec![gd])
}

fn pair(id: i32, a: (i32, &str), b: (i32, &str), extra: u16, mapq: u8) -> Vec<Read> {
    let mut r1 = read(id, "1", a.0, a.1, None);
    let mut r2 = read(id, "1", b.0, b.1, None);
    r1.flags |= FLAG_FIRST | crate::read::FLAG_MATE_REVERSE | extra;
    r2.flags |= FLAG_SECOND | FLAG_REVERSE | extra;
    r1.mate_pos = b.0;
    r2.mate_pos = a.0;
    r1.insert_size = 300;
    r2.insert_size = -300;
    r1.mapq = mapq;
    r2.mapq = mapq;
    vec![r1, r2]
}

/// the fragment counts and expression weights of 2.0.1's allocator: duplicates counted, secondaries not
/// counted by type, low map qualities weighted down
#[test]
fn allocator_201_counts() {
    with_v201(true, || {
        let mut config = IsofoxConfig::for_tests();
        config.functions = vec![IsofoxFunction::TranscriptCounts];
        let sites: HashMap<String, HashSet<i32>> = HashMap::new();
        let mut refg = RefSource::mock(vec![("1", gen_bases(10000))]);
        let mut alloc = FragmentAllocator::new(&config, None, &sites).unwrap();
        let mut gc = test_gc();
        // each call's category counts (the harness resets them)
        let weight = |a: &FragmentAllocator| -> (f64, i32) {
            (a.trans_combo_data.iter().map(|c| c.fragment_count).sum(), a.trans_combo_data.iter().map(|c| c.multi_mapped).sum())
        };
        // a supporting fragment
        alloc.process_read_records_201(&mut gc, pair(1, (1050, "50M"), (2100, "100M"), 0, 255), &mut refg);
        assert_eq!(gc.fragment_counts.type_count(FragmentType::TransSupporting), 1);
        assert_eq!(weight(&alloc), (1.0, 0));
        // a duplicate counts too
        alloc.process_read_records_201(&mut gc, pair(2, (1050, "50M"), (2100, "100M"), FLAG_DUPLICATE, 255), &mut refg);
        assert_eq!(gc.fragment_counts.type_count(FragmentType::TransSupporting), 2);
        assert_eq!(weight(&alloc), (1.0, 0));
        // a secondary pair (map quality 1): in expression with a weight of 0.2, not in the type counts
        alloc.process_read_records_201(&mut gc, pair(3, (1050, "50M"), (2100, "100M"), FLAG_SECONDARY, 1), &mut refg);
        assert_eq!(gc.fragment_counts.type_count(FragmentType::TransSupporting), 2);
        assert_eq!(gc.fragment_counts.type_count(FragmentType::MultiMapped), 1);
        assert_eq!(weight(&alloc), (0.2, 1));
        // a low map quality fragment
        alloc.process_read_records_201(&mut gc, pair(4, (1050, "50M"), (2100, "100M"), 0, 3), &mut refg);
        assert_eq!(gc.fragment_counts.type_count(FragmentType::TransSupporting), 3);
        assert_eq!(gc.fragment_counts.type_count(FragmentType::MultiMapped), 2);
        assert_eq!(weight(&alloc), (0.5, 1));
        // an improper pair is chimeric
        let mut gc = test_gc();
        let mut reads = pair(5, (1050, "50M"), (2100, "100M"), 0, 255);
        for r in &mut reads {
            r.flags &= !FLAG_PROPER;
        }
        alloc.process_read_records_201(&mut gc, reads, &mut refg);
        // chimeric handling is off (no fusions or novel junctions), so the fragment is counted as chimeric directly
        assert_eq!(gc.fragment_counts.type_count(FragmentType::Chimeric), 1);
        assert_eq!(gc.fragment_counts.type_count(FragmentType::TransSupporting), 0);
        let _ = parse_cigar("1M");
    });
}
