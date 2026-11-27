//! Ports of ISOFOX's Java unit tests (ReadUtilsTest, ReadCountsTest, TransClassificationTest, FragmentTest,
//! GcBiasTest, MiscTest, ExpectedRatesTest) and tests of the port's own. All data is synthetic.

use crate::adjusts::GcRatioCounts;
use crate::allocator::FragmentAllocator;
use crate::common::*;
use crate::config::{FragmentSize, IsofoxConfig, IsofoxFunction};
use crate::expression::{calc_effective_length, perform_fit};
use crate::fragment;
use crate::gene::{GeneCollection, GeneReadData, RegionReadData};
use crate::read::{parse_cigar, trim_adapter_bases, Read, RKey, FLAG_FIRST, FLAG_PAIRED, FLAG_PROPER, FLAG_REVERSE, FLAG_SECOND};
use crate::read_trans::{calc_fragment_length, mark_region_bases, process_overlapping_regions, region_match_type};
use crate::refgen::RefSource;
use hmf_common::ensembl::{ExonData, GeneData, TranscriptData};
use hmf_common::regions::BaseRegion;
use std::collections::{HashMap, HashSet};

pub fn leak(s: &str) -> &'static str {
    Box::leak(s.to_string().into_boxed_str())
}

/// MockRefGenome.generateRandomBases: ACGT repeated
pub fn gen_bases(n: usize) -> String {
    "ACGT".chars().cycle().take(n).collect()
}

fn ref_base_str_1() -> String {
    "ACGT".repeat(5)
}

pub fn gene(id: &str, name: &str, chr: &str, strand: i8, start: i32, end: i32) -> &'static GeneData {
    Box::leak(Box::new(GeneData {
        gene_id: id.into(),
        gene_name: name.into(),
        chromosome: chr.into(),
        strand,
        gene_start: start,
        gene_end: end,
        karyotype_band: String::new(),
        synonyms: String::new(),
    }))
}

pub fn transcript(trans_id: i32, name: &str, gene_id: &str, strand: i8, exons: &[(i32, i32)], canonical: bool) -> &'static TranscriptData {
    let n = exons.len() as i32;
    let ex: Vec<ExonData> = exons
        .iter()
        .enumerate()
        .map(|(i, (s, e))| ExonData { trans_id, start: *s, end: *e, rank: if strand == 1 { i as i32 + 1 } else { n - i as i32 }, phase_start: -1, phase_end: -1 })
        .collect();
    Box::leak(Box::new(TranscriptData {
        trans_id,
        trans_name: name.into(),
        gene_id: gene_id.into(),
        is_canonical: canonical,
        strand,
        trans_start: exons[0].0,
        trans_end: exons[exons.len() - 1].1,
        coding_start: None,
        coding_end: None,
        bio_type: "protein_coding".into(),
        ref_seq_id: None,
        exons: ex,
    }))
}

/// createRegion: a region of one transcript's exon
pub fn region(gene_id: &str, trans_id: i32, rank: i32, start: i32, end: i32) -> RegionReadData {
    let mut r = RegionReadData::new(start, end);
    r.trans_exon_refs.push(TransExonRef { gene_id: leak(gene_id), trans_id, trans_name: leak(&trans_id.to_string()), exon_rank: rank, is_canonical: Some(true) });
    r
}

/// createReadRecord: a paired, properly-paired read with its mate at its own start
pub fn read(id: i32, chr: &'static str, pos_start: i32, cigar: &str, bases: Option<&str>) -> Read {
    let c = parse_cigar(cigar).unwrap();
    let len: u32 = c.iter().filter(|e| crate::read::consumes_read(e.op)).map(|e| e.len).sum();
    let b = bases.map(|s| s.to_string()).unwrap_or_else(|| gen_bases(len as usize));
    Read::new(&id.to_string(), chr, pos_start, FLAG_PAIRED | FLAG_PROPER, c, b.into_bytes(), chr, pos_start, 0, 60, None, None, None)
}

fn link(regions: &mut [RegionReadData], a: usize, b: usize) {
    regions[a].post_regions.push(b);
    regions[b].pre_regions.push(a);
}

fn classify(read: &mut Read, regions: &[RegionReadData]) {
    let all: Vec<usize> = (0..regions.len()).collect();
    process_overlapping_regions(read, regions, &all);
}

fn match_type(read: &Read, r: &RegionReadData) -> RegionMatchType {
    let i = read.mapped_coords.find_region_index(r.start, r.end);
    region_match_type(&read.mapped_coords, r, i)
}

// ---- ReadUtilsTest ----

#[test]
fn mapping_coords() {
    let r = read(1, "1", 100, "10M10I10M10D10M", None);
    assert_eq!(r.mapped_coords.alignments, vec![BaseRegion::new(100, 139)]);
    assert!(r.mapped_coords.alignments_overlap(90, 101));
    assert!(r.mapped_coords.alignments_overlap(138, 150));
    assert!(r.mapped_coords.alignments_within(99, 139));
    let mut r = read(1, "1", 100, "11M9N11M9N11M", None);
    assert_eq!(r.mapped_coords.alignments, vec![BaseRegion::new(100, 110), BaseRegion::new(120, 130), BaseRegion::new(140, 150)]);
    assert_eq!(r.mapped_coords.original_alignment_count(), 3);
    r.mapped_coords.add_inferred_region(true, 80, 90);
    assert_eq!(r.mapped_coords.alignment_count(), 4);
    assert_eq!(r.mapped_coords.original_alignment_count(), 3);
    r.mapped_coords.add_inferred_region(false, 160, 170);
    assert_eq!(r.mapped_coords.alignment_count(), 5);
    assert_eq!(r.mapped_coords.without_inferred().len(), 3);
}

#[test]
fn overlapping_coordinates() {
    let b = BaseRegion::new;
    assert_eq!(derive_common_regions(&[b(10, 20), b(40, 50)], &[b(60, 70), b(80, 90)]).len(), 4);
    let m1 = vec![b(10, 20), b(40, 50), b(70, 80)];
    let c = derive_common_regions(&m1, &[b(25, 35), b(55, 65), b(85, 95)]);
    assert_eq!(c.len(), 6);
    let c = derive_common_regions(&m1, &[b(5, 15), b(35, 45), b(55, 75)]);
    assert_eq!(c, vec![b(5, 20), b(35, 50), b(55, 80)]);
    let c = derive_common_regions(&m1, &[b(5, 95)]);
    assert_eq!(c, vec![b(5, 95)]);
    let c = derive_common_regions(&[b(10, 20), b(55, 75), b(85, 95), b(140, 150)], &[b(30, 40), b(50, 60), b(70, 110), b(120, 130)]);
    assert_eq!(c.len(), 5);
    assert_eq!(c[2], b(50, 110));
}

#[test]
fn fragment_tracking() {
    // FragmentTracker: a read waits for its mate by read id
    let mut tracker: hmf_common::java_seq_hash::SeqHashMap<String, usize> = hmf_common::java_seq_hash::SeqHashMap::new();
    for (i, id) in ["1", "2", "3"].iter().enumerate() {
        assert!(tracker.remove(&id.to_string()).is_none());
        tracker.insert(id.to_string(), i);
    }
    assert_eq!(tracker.len(), 3);
    for (i, id) in ["1", "2", "3"].iter().enumerate() {
        assert_eq!(tracker.remove(&id.to_string()), Some(i));
    }
    assert_eq!(tracker.len(), 0);
}

#[test]
fn base_assignment() {
    let mut r = region("GEN01", 1, 1, 100, 119);
    r.set_ref_bases(ref_base_str_1().into_bytes());
    mark_region_bases(&[BaseRegion::new(100, 119)], &mut r);
    assert_eq!(r.base_coverage(1), 20);
    let mut r = region("GEN01", 1, 1, 100, 119);
    r.set_ref_bases(ref_base_str_1().into_bytes());
    mark_region_bases(&[BaseRegion::new(100, 104), BaseRegion::new(110, 114), BaseRegion::new(118, 119)], &mut r);
    assert_eq!(r.base_coverage(1), 12);
}

#[test]
fn adapter_trimming() {
    let b2 = ref_base_str_1().repeat(2);
    let mk = |pos: i32, cig: &str| read(1, "1", pos, cig, Some(&b2));
    let mut r1 = mk(105, "10S20M10S");
    let mut r2 = mk(100, "10S20M10S");
    r2.flags |= FLAG_REVERSE;
    assert_eq!((r1.unclipped_start, r1.unclipped_end, r2.unclipped_start, r2.unclipped_end), (95, 134, 90, 129));
    trim_adapter_bases(&mut r1, &mut r2);
    assert_eq!(r1.unclipped_end, 129);
    assert_eq!(r2.unclipped_start, 95);
    assert_eq!(r1.cigar_str(), "10S20M5S");
    assert_eq!(r2.cigar_str(), "5S20M10S");
    // passed in reverse
    let mut r1 = mk(105, "10S20M10S");
    let mut r2 = mk(100, "10S20M10S");
    r2.flags |= FLAG_REVERSE;
    trim_adapter_bases(&mut r2, &mut r1);
    assert_eq!(r1.cigar_str(), "10S20M5S");
    assert_eq!(r2.cigar_str(), "5S20M10S");
    // trimming capped at the soft clips
    let mut r1 = mk(105, "10S20M10S");
    let mut r2 = mk(94, "10S20M10S");
    r2.flags |= FLAG_REVERSE;
    trim_adapter_bases(&mut r1, &mut r2);
    assert_eq!(r1.unclipped_end, 124);
    assert_eq!(r2.unclipped_start, 94);
    assert_eq!(r1.cigar_str(), "10S20M");
    assert_eq!(r2.cigar_str(), "20M10S");
    // N-splits giving conflicting alignments: no trimming
    let mut r1 = read(1, "1", 100, "15M885N20M5S", Some(&b2));
    let mut r2 = read(1, "1", 1000, "5S20M980N15M", Some(&b2));
    r2.flags |= FLAG_REVERSE;
    trim_adapter_bases(&mut r2, &mut r1);
    assert_eq!((r1.unclipped_start, r1.unclipped_end, r2.unclipped_start, r2.unclipped_end), (100, 1024, 995, 2014));
    // trimming from both sides
    let mut r1 = read(1, "1", 110, "5M885N20M10S", Some(&b2));
    let mut r2 = read(1, "1", 1000, "10S20M980N5M", Some(&b2));
    r2.flags |= FLAG_REVERSE;
    trim_adapter_bases(&mut r2, &mut r1);
    assert_eq!(r1.unclipped_end, 1024);
    assert_eq!(r2.unclipped_start, 995);
    assert_eq!(r1.cigar_str(), "5M885N20M5S");
    assert_eq!(r2.cigar_str(), "5S20M980N5M");
}

// ---- ReadCountsTest ----

#[test]
fn read_coordinates() {
    let r = read(1, "1", 100, "10S20M10S", None);
    assert_eq!(r.mapped_coords.alignments, vec![BaseRegion::new(100, 119)]);
    let r = read(1, "1", 100, "5M20N10M30N5M", None);
    assert_eq!(r.mapped_coords.alignments, vec![BaseRegion::new(100, 104), BaseRegion::new(125, 134), BaseRegion::new(165, 169)]);
    let r = read(1, "1", 100, "10M2D8M", None);
    assert_eq!(r.mapped_coords.alignments, vec![BaseRegion::new(100, 119)]);
    let r = read(1, "1", 100, "10M2I8M", None);
    assert_eq!(r.mapped_coords.alignments, vec![BaseRegion::new(100, 117)]);
}

#[test]
fn read_region_types() {
    let reg = RegionReadData::new(100, 200);
    assert_eq!(match_type(&read(1, "1", 110, "21M", None), &reg), RegionMatchType::WithinExon);
    assert_eq!(match_type(&read(1, "1", 90, "21M", None), &reg), RegionMatchType::ExonIntron);
    assert_eq!(match_type(&read(1, "1", 100, "101M", None), &reg), RegionMatchType::ExonBoundary);
    assert_eq!(match_type(&read(1, "1", 100, "51M", None), &reg), RegionMatchType::ExonBoundary);
    let r = read(1, "1", 110, "11M19N21M19N25M", None);
    assert_eq!(r.mapped_coords.alignments, vec![BaseRegion::new(110, 120), BaseRegion::new(140, 160), BaseRegion::new(180, 204)]);
    assert_eq!(match_type(&r, &RegionReadData::new(100, 120)), RegionMatchType::ExonBoundary);
    assert_eq!(match_type(&r, &RegionReadData::new(140, 160)), RegionMatchType::ExonBoundary);
    assert_eq!(match_type(&r, &RegionReadData::new(180, 200)), RegionMatchType::ExonIntron);
}

#[test]
fn unmapped_splice_junctions() {
    // soft-clipping at the start matching the previous exon's last bases
    let s1 = ref_base_str_1();
    let mut regions = vec![RegionReadData::new(141, 150), RegionReadData::new(200, 209)];
    regions[0].set_ref_bases(s1[0..10].as_bytes().to_vec());
    regions[1].set_ref_bases(s1[10..19].as_bytes().to_vec());
    link(&mut regions, 0, 1);
    let mut r = read(1, "1", 200, "5S10M", Some(&s1[5..19]));
    process_overlapping_regions(&mut r, &regions, &[1]);
    assert_eq!(r.mapped_regions.get(0), Some(RegionMatchType::ExonBoundary));
    assert_eq!(r.mapped_coords.alignments.len(), 2);
    assert_eq!(r.mapped_coords.alignments[0], BaseRegion::new(146, 150));
}

// ---- TransClassificationTest ----

#[test]
fn read_transcript_classification() {
    let s = ref_base_str_1();
    let regions = vec![region("GEN01", 1, 1, 100, 200)];
    let mut r = read(1, "1", 90, "21M", Some(&s));
    classify(&mut r, &regions);
    assert_eq!(r.trans_classification(1), TransMatchType::Unspliced);
    let mut r = read(1, "1", 120, "21M", Some(&s));
    classify(&mut r, &regions);
    assert_eq!(r.trans_classification(1), TransMatchType::Exonic);
    // skipped exon
    let regions = vec![region("GEN01", 1, 1, 100, 120), region("GEN01", 1, 2, 140, 160), region("GEN01", 1, 3, 180, 200)];
    let mut r = read(1, "1", 110, "11M59N21M", Some(&s));
    classify(&mut r, &regions);
    assert_eq!(r.trans_classification(1), TransMatchType::Alt);
    // incomplete intermediary exon
    let mut r = read(1, "1", 110, "11M19N11M29N25M", Some(&s));
    classify(&mut r, &regions);
    assert_eq!(r.trans_classification(1), TransMatchType::Alt);
    // spanning into an intron
    let mut r = read(1, "1", 110, "11M19N31M9N25M", Some(&s));
    classify(&mut r, &regions);
    assert_eq!(r.trans_classification(1), TransMatchType::Unspliced);
    // two regions as a splice junction and a third as exonic
    let regions = vec![region("GEN01", 1, 1, 100, 200), region("GEN01", 1, 2, 300, 400), region("GEN01", 2, 1, 100, 220)];
    let mut r = read(1, "1", 181, "20M99N20M", Some(&s.repeat(2)));
    classify(&mut r, &regions);
    assert_eq!(r.trans_classification(1), TransMatchType::SpliceJunction);
    assert_eq!(r.trans_classification(2), TransMatchType::Alt);
}

#[test]
fn soft_clipped_read_transcript_classification() {
    let s = ref_base_str_1();
    let mut regions = vec![region("ENSG0001", 1, 1, 100, 119), region("ENSG0001", 1, 2, 200, 219), region("ENSG0001", 1, 3, 300, 319)];
    for r in &mut regions {
        r.set_ref_bases(s.as_bytes().to_vec());
    }
    link(&mut regions, 0, 1);
    link(&mut regions, 1, 2);
    let bases = format!("{}{}{}", &s[s.len() - 1..], s, &s[..10]);
    let mut r = read(1, "1", 200, "1S20M80N10M", Some(&bases));
    process_overlapping_regions(&mut r, &regions, &[1, 2]);
    assert_eq!(r.trans_classification(1), TransMatchType::SpliceJunction);
    // a soft clip that cannot be mapped makes ALT
    let mut r = read(1, "1", 100, "20M5S", Some(&format!("{s}AAAAA")));
    r.insert_size = 200;
    process_overlapping_regions(&mut r, &regions, &[0]);
    assert_eq!(r.trans_classification(1), TransMatchType::Alt);
    let mut r = read(1, "1", 300, "5S20M", Some(&format!("AAAAA{s}")));
    r.insert_size = 200;
    process_overlapping_regions(&mut r, &regions, &[2]);
    assert_eq!(r.trans_classification(1), TransMatchType::Alt);
    // likely adapter sequence is allowed
    let mut r = read(1, "1", 300, "5S20M", Some(&format!("AAAAA{s}")));
    r.insert_size = 20;
    process_overlapping_regions(&mut r, &regions, &[2]);
    assert_eq!(r.trans_classification(1), TransMatchType::Exonic);
    // short soft clips within the exon stay supporting
    let mut r = read(1, "1", 300, "18M1S", Some(&format!("{}A", &s[..18])));
    r.insert_size = 200;
    process_overlapping_regions(&mut r, &regions, &[2]);
    assert_eq!(r.trans_classification(1), TransMatchType::Exonic);
    let mut r = read(1, "1", 102, "2S18M", Some(&format!("AA{}", &s[..18])));
    r.insert_size = 200;
    process_overlapping_regions(&mut r, &regions, &[0]);
    assert_eq!(r.trans_classification(1), TransMatchType::Exonic);
    // a clip of the realignment window's length stays ALT
    let mut r = read(1, "1", 300, "17M3S", Some(&format!("{}AAA", &s[..17])));
    r.insert_size = 200;
    process_overlapping_regions(&mut r, &regions, &[2]);
    assert_eq!(r.trans_classification(1), TransMatchType::Alt);
}

#[test]
fn fragment_length_calcs() {
    let t = transcript(1, "TRANS01", "GENE1", 1, &[(1000, 1200), (2000, 2500), (4500, 5000), (5500, 6000)], true);
    let f = |a: (i32, &str), b: (i32, &str)| {
        let r1 = read(1, "1", a.0, a.1, None);
        let r2 = read(1, "1", b.0, b.1, None);
        let rs = [&r1, &r2];
        calc_fragment_length(t, fragment::min_alignment_start(&rs), fragment::max_alignment_end(&rs))
    };
    assert_eq!(f((1010, "20M"), (1170, "30M")), 190);
    assert_eq!(f((1170, "20M"), (2010, "10M")), 31 + 20);
    assert_eq!(f((1170, "20M"), (4510, "10M")), 31 + 501 + 20);
    assert_eq!(f((1191, "10M799N10M"), (2491, "10M1999N10M")), 10 + 501 + 10);
    assert_eq!(f((1191, "10M799N10M"), (4991, "10M499N10M")), 10 + 501 + 501 + 10);
}

fn test_gene_collection() -> GeneCollection {
    let g = gene("GENE1", "GENE1", "1", 1, 1000, 5000);
    let t = transcript(1, "TRANS01", "GENE1", 1, &[(1000, 1200), (2000, 2500), (4500, 5000)], true);
    let gd = GeneReadData { gene: g, exon_regions: Vec::new(), transcripts: vec![t], has_unspliced_regions: false };
    GeneCollection::new(0, vec![gd])
}

fn read_pair(id: i32, a: (i32, &str), b: (i32, &str), orients: (i8, i8)) -> (Read, Read) {
    let mut r1 = read(id, "1", a.0, a.1, None);
    let mut r2 = read(id, "1", b.0, b.1, None);
    r1.flags |= FLAG_FIRST;
    r2.flags |= FLAG_SECOND;
    let set = |r: &mut Read, rev: bool, mate_rev: bool| {
        if rev {
            r.flags |= FLAG_REVERSE;
        }
        if mate_rev {
            r.flags |= crate::read::FLAG_MATE_REVERSE;
        }
    };
    set(&mut r1, orients.0 == -1, orients.1 == -1);
    set(&mut r2, orients.1 == -1, orients.0 == -1);
    r1.mate_pos = b.0;
    r2.mate_pos = a.0;
    (r1, r2)
}

#[test]
fn fragment_read_pairs() {
    let mut config = IsofoxConfig::for_tests();
    config.functions = vec![IsofoxFunction::TranscriptCounts];
    let sites: HashMap<String, HashSet<i32>> = HashMap::new();
    let mut refg = RefSource::mock(vec![("1", gen_bases(10000))]);
    let mut alloc = FragmentAllocator::new(&config, None, &sites).unwrap();
    // one read outside the gene: chimeric
    let mut gc = test_gene_collection();
    alloc.process_read_records(&mut gc, vec![read(1, "1", 100, "101M", None), read(1, "1", 1050, "101M", None)], &mut refg);
    assert_eq!(gc.fragment_counts.type_count(FragmentType::Chimeric), 1);
    // exon to intron: unspliced
    let (a, b) = read_pair(2, (1300, "51M"), (2300, "101M"), (1, -1));
    alloc.process_read_records(&mut gc, vec![a, b], &mut refg);
    assert_eq!(gc.fragment_counts.type_count(FragmentType::Unspliced), 1);
    // fully intronic
    let mut gc = test_gene_collection();
    let (a, b) = read_pair(3, (2600, "51M"), (2800, "101M"), (1, -1));
    alloc.process_read_records(&mut gc, vec![a, b], &mut refg);
    assert_eq!(gc.fragment_counts.type_count(FragmentType::Unspliced), 1);
    // a split read beyond the gene: chimeric
    let mut gc = test_gene_collection();
    let (a, b) = read_pair(4, (1050, "51M4899N101M"), (3100, "201M"), (1, -1));
    alloc.process_read_records(&mut gc, vec![a, b], &mut refg);
    assert_eq!(gc.fragment_counts.type_count(FragmentType::Chimeric), 1);
    // long exon to exon fragment: supporting
    let mut gc = test_gene_collection();
    let (a, b) = read_pair(5, (1050, "50M"), (2100, "100M"), (1, -1));
    alloc.process_read_records(&mut gc, vec![a, b], &mut refg);
    assert_eq!(gc.fragment_counts.type_count(FragmentType::TransSupporting), 1);
    // skipping an exon, and long: alt
    let mut gc = test_gene_collection();
    let (a, b) = read_pair(6, (1050, "51M"), (4550, "101M"), (1, -1));
    alloc.process_read_records(&mut gc, vec![a, b], &mut refg);
    assert_eq!(gc.fragment_counts.type_count(FragmentType::Alt), 1);
}

// ---- FragmentTest ----

fn single_end(pos: i32, cigar: &str) -> Read {
    let mut r = read(1, "1", pos, cigar, None);
    r.flags = 0;
    r
}

#[test]
fn single_end_fragment() {
    let r = single_end(100, "5S10M100N10M5S");
    let rs = [&r];
    assert!(!r.is_read_paired());
    assert_eq!(fragment::fragment_count(&rs), 1);
    assert_eq!(fragment::min_alignment_start(&rs), 100);
    assert_eq!(fragment::max_alignment_end(&rs), 219);
    assert!(fragment::contains_split(&rs));
    assert!(fragment::is_fully_intronic(&rs));
    assert!(fragment::unique_valid_regions(&rs).is_empty());
    assert_eq!(fragment::merged_mappings(&rs), vec![BaseRegion::new(100, 109), BaseRegion::new(210, 219)]);
    assert_eq!(fragment::orientation(&rs), Some(ORIENT_FWD));
}

#[test]
fn consensus_count_not_summed_across_mates() {
    let mut a = read(1, "1", 100, "20M", None);
    let mut b = read(1, "1", 200, "20M", None);
    a.consensus = Some("5;5".into());
    b.consensus = Some("5;5".into());
    assert_eq!(fragment::fragment_count(&[&a]), 5);
    assert_eq!(fragment::fragment_count(&[&a, &b]), 5);
    assert_eq!(fragment::fragment_count(&[&b, &a]), 5);
}

fn multi_mapped(xa: &str) -> Read {
    let c = parse_cigar("20M").unwrap();
    Read::new("r", "1", 100, FLAG_PAIRED, c, gen_bases(20).into_bytes(), "1", 100, 0, 60, None, Some(xa), None)
}

#[test]
fn least_ambiguous_read_determines_loci() {
    let unique = read(1, "1", 100, "20M", None);
    let ambiguous = multi_mapped("2,+500,20M,0;3,-700,20M,1;");
    let less = multi_mapped("4,+900,20M,0;");
    assert_eq!(fragment::min_num_loci(&[&unique]), 1);
    assert!(fragment::alt_loci(&[&unique]).is_none());
    assert_eq!(fragment::min_num_loci(&[&ambiguous]), 3);
    assert_eq!(fragment::alt_loci(&[&ambiguous]).unwrap().len(), 2);
    for rs in [[&ambiguous, &unique], [&unique, &ambiguous]] {
        assert_eq!(fragment::min_num_loci(&rs), 1);
        assert!(fragment::alt_loci(&rs).is_none());
    }
    for rs in [[&ambiguous, &less], [&less, &ambiguous]] {
        assert_eq!(fragment::min_num_loci(&rs), 2);
        assert_eq!(fragment::alt_loci(&rs).unwrap()[0].start, 900);
    }
    // equal locus counts keep the first read given
    let first = multi_mapped("2,+500,20M,0;");
    let second = multi_mapped("3,+700,20M,0;");
    assert_eq!(fragment::alt_loci(&[&first, &second]).unwrap()[0].start, 500);
    assert_eq!(fragment::alt_loci(&[&second, &first]).unwrap()[0].start, 700);
}

#[test]
fn xa_parsing() {
    let r = multi_mapped("chr2,+500,10M200N10M,0;chr3,-700,,1;bad;chr4,x,20M,0;");
    let a = r.alt_loci.as_ref().unwrap();
    assert_eq!(a.len(), 2);
    assert_eq!((a[0].start, a[0].end, a[0].spliced), (500, 719, true));
    assert_eq!((a[1].start, a[1].end, a[1].spliced), (700, 700, false));
    assert_eq!(r.num_loci(), 3);
    assert!(multi_mapped("").alt_loci.is_none());
}

#[test]
fn merged_mappings() {
    let split = read(1, "1", 100, "10M100N10M", None);
    let overlapping = read(1, "1", 105, "10M100N15M", None);
    for rs in [[&split, &overlapping], [&overlapping, &split]] {
        assert_eq!(fragment::merged_mappings(&rs), vec![BaseRegion::new(100, 114), BaseRegion::new(210, 229)]);
    }
    let downstream = read(1, "1", 300, "20M", None);
    for rs in [[&split, &downstream], [&downstream, &split]] {
        assert_eq!(fragment::merged_mappings(&rs), vec![BaseRegion::new(100, 109), BaseRegion::new(210, 219), BaseRegion::new(300, 319)]);
    }
    let up = read(1, "1", 100, "30M", None);
    let adjacent = read(1, "1", 120, "30M", None);
    assert_eq!(fragment::merged_mappings(&[&up, &adjacent]), vec![BaseRegion::new(100, 149)]);
}

#[test]
fn unique_valid_regions_and_spans() {
    let regions = vec![region("gene", 1, 1, 100, 199), region("gene", 1, 2, 210, 299), region("other", 2, 1, 300, 399)];
    let mut first = read(1, "1", 100, "10M100N10M", None);
    let mut second = read(1, "1", 100, "20M", None);
    let k = |i: usize| RKey { idx: i as u32, start: regions[i].start };
    first.mapped_regions.put(k(0), RegionMatchType::WithinExon);
    first.mapped_regions.put(k(1), RegionMatchType::ExonIntron);
    first.mapped_regions.put(k(2), RegionMatchType::ExonBoundary);
    second.mapped_regions.put(k(0), RegionMatchType::ExonBoundary);
    assert!(!fragment::is_fully_intronic(&[&first]));
    let mut u = fragment::unique_valid_regions(&[&first, &second]);
    u.sort();
    assert_eq!(u, vec![0, 2]);
    // the exon-intron region does not count towards spanning several regions
    assert!(!fragment::spans_multiple_regions(&[&first], &regions, 1));
    let mut third = read(1, "1", 210, "20M", None);
    third.mapped_regions.put(k(1), RegionMatchType::WithinExon);
    assert!(fragment::spans_multiple_regions(&[&first, &third], &regions, 1));
}

#[test]
fn transcript_support_is_an_intersection() {
    let mut a = read(1, "1", 100, "20M", None);
    let mut b = read(1, "1", 200, "20M", None);
    a.trans_class.insert(1, TransMatchType::Exonic);
    a.trans_class.insert(2, TransMatchType::SpliceJunction);
    a.trans_class.insert(3, TransMatchType::Alt);
    b.trans_class.insert(1, TransMatchType::SpliceJunction);
    b.trans_class.insert(2, TransMatchType::Unspliced);
    b.trans_class.insert(4, TransMatchType::Exonic);
    assert_eq!(fragment::valid_type_transcripts(&[&a, &b]), vec![1]);
    let inv = fragment::invalid_transcripts(&[&a, &b], &[1]);
    assert_eq!(inv.keys(), vec![2, 3, 4]);
    // the other transcripts keep their own invalid types; valid ones not chosen become OTHER_TRANS
    {
        let mut rs: Vec<&mut Read> = vec![&mut a, &mut b];
        fragment::set_other_transcripts(&mut rs, &[1]);
    }
    assert_eq!(a.trans_classification(1), TransMatchType::Exonic);
    assert_eq!(a.trans_classification(2), TransMatchType::OtherTrans);
    assert_eq!(a.trans_classification(3), TransMatchType::Alt);
    assert_eq!(b.trans_classification(2), TransMatchType::Unspliced);
    assert_eq!(b.trans_classification(4), TransMatchType::OtherTrans);
}

#[test]
fn orientation_uses_first_of_pair() {
    let (a, b) = read_pair(1, (300, "20M"), (100, "20M"), (-1, 1));
    assert_eq!(fragment::orientation(&[&a, &b]), Some(ORIENT_REV));
    assert_eq!(fragment::orientation(&[&b, &a]), Some(ORIENT_REV));
    let (a, b) = read_pair(1, (300, "20M"), (100, "20M"), (1, 1));
    assert_eq!(fragment::orientation(&[&a, &b]), None);
}

#[test]
fn reads_in_different_exons() {
    let regions = vec![region("gene", 1, 1, 100, 199), region("gene", 1, 2, 300, 399), region("gene", 2, 1, 300, 399)];
    let k = |i: usize| RKey { idx: i as u32, start: regions[i].start };
    let mut a = read(1, "1", 150, "20M", None);
    let mut b = read(1, "1", 320, "20M", None);
    a.mapped_regions.put(k(0), RegionMatchType::WithinExon);
    b.mapped_regions.put(k(1), RegionMatchType::WithinExon);
    b.mapped_regions.put(k(2), RegionMatchType::WithinExon);
    assert!(fragment::reads_in_different_exons(&[&a, &b], &regions, 1));
    // the rank is compared within the transcript asked for: transcript 2's exon 1 is not transcript 1's
    assert!(fragment::reads_in_different_exons(&[&a, &b], &regions, 2));
    b.mapped_regions.put(k(0), RegionMatchType::WithinExon);
    assert!(!fragment::reads_in_different_exons(&[&a, &b], &regions, 1));
    assert!(!fragment::reads_in_different_exons(&[&a], &regions, 1));
}

// ---- GcBiasTest, MiscTest, ExpectedRatesTest ----

#[test]
fn gc_ratio_calcs() {
    let g = GcRatioCounts::new();
    let mut idx = [-1, -1];
    let mut c = [1.0, 0.0];
    g.determine_ratio_data(0.0, &mut idx, &mut c);
    assert_eq!(idx, [0, 1]);
    assert!((c[0] - 1.0).abs() < 0.001 && c[1].abs() < 0.001);
    g.determine_ratio_data(0.505, &mut idx, &mut c);
    assert_eq!(idx, [50, 51]);
    assert!((c[0] - 0.5).abs() < 0.001);
    g.determine_ratio_data(0.9925, &mut idx, &mut c);
    assert_eq!(idx, [99, 100]);
    assert!((c[0] - 0.75).abs() < 0.001 && (c[1] - 0.25).abs() < 0.001);
    g.determine_ratio_data(1.0, &mut idx, &mut c);
    assert_eq!(idx[0], 100);
    assert!((c[0] - 1.0).abs() < 0.001);
}

#[test]
fn effective_length_misc() {
    let fl = vec![FragmentSize { length: 100, frequency: 1 }, FragmentSize { length: 200, frequency: 2 }, FragmentSize { length: 300, frequency: 1 }];
    assert!((calc_effective_length(1000, &fl) - 800.0).abs() < 0.001);
    assert!((calc_effective_length(250, &fl) - 83.33).abs() < 0.1);
}

#[test]
fn base_depth_misc() {
    let mut bd = BaseDepth::new();
    bd.initialise([100, 200]);
    bd.process_read(&[BaseRegion::new(100, 200)]);
    bd.process_read(&[BaseRegion::new(100, 150)]);
    bd.process_read(&[BaseRegion::new(50, 60), BaseRegion::new(250, 260)]);
    assert_eq!(bd.bases_with_depth(), 51);
    assert_eq!(bd.depth_at_base(100), 2);
    assert_eq!(bd.depth_at_base(150), 2);
    assert_eq!(bd.depth_at_base(200), 1);
    assert_eq!(bd.depth_at_base(99), 0);
    assert_eq!(bd.depth_at_base(201), 0);
    bd.process_read(&[BaseRegion::new(100, 110), BaseRegion::new(120, 130)]);
    assert_eq!(bd.depth_at_base(100), 3);
    assert_eq!(bd.depth_at_base(115), 2);
    assert_eq!(bd.depth_at_base(120), 3);
}

#[test]
fn expectation_max_fit() {
    // column-major: transcript 1 = {0.2, 0.8, 0}, transcript 2 = {0.4, 0, 0.6}
    let defs = vec![0.2, 0.8, 0.0, 0.4, 0.0, 0.6];
    let a = perform_fit(&[5.0, 4.0, 6.0], &defs, 3, 2);
    assert!((a[0] - 5.002).abs() < 0.001, "{a:?}");
    assert!((a[1] - 9.998).abs() < 0.001, "{a:?}");
    let a = perform_fit(&[5.0, 4.0, 7.0], &defs, 3, 2);
    assert!((a[0] - 4.905).abs() < 0.001, "{a:?}");
    assert!((a[1] - 11.095).abs() < 0.001, "{a:?}");
}

// ---- the port's own ----

#[test]
fn region_map_orders_by_start_beyond_ten_regions() {
    // eleven or more regions in a read's map: Java's bucket becomes a tree ordered by region start
    let mut r = read(1, "1", 100, "20M", None);
    let starts = [50, 10, 40, 20, 30, 60, 5, 45, 15, 25, 35, 55];
    for (i, s) in starts.iter().enumerate() {
        r.mapped_regions.put(RKey { idx: i as u32, start: *s }, RegionMatchType::WithinExon);
    }
    let keys = r.mapped_regions.keys();
    assert_eq!(keys.len(), 12);
    // up to ten the order is the insertion order
    let mut r10 = read(1, "1", 100, "20M", None);
    for (i, s) in starts.iter().take(10).enumerate() {
        r10.mapped_regions.put(RKey { idx: i as u32, start: *s }, RegionMatchType::WithinExon);
    }
    assert_eq!(r10.mapped_regions.keys(), (0..10).collect::<Vec<_>>());
    assert_ne!(keys[..10].to_vec(), (0..10).collect::<Vec<_>>());
}

#[test]
fn intronic_length_finding() {
    // IF2: a gene of 300 bases with two exons of 101 bases has 98 intronic bases
    let g = gene("G", "G", "1", 1, 1000, 1299);
    let t = transcript(1, "T", "G", 1, &[(1000, 1100), (1199, 1299)], true);
    let gc = GeneCollection::new(0, vec![GeneReadData { gene: g, exon_regions: Vec::new(), transcripts: vec![t], has_unspliced_regions: false }]);
    let corrected = crate::compat::with_java(false, || crate::expression::GeneResult::new(&gc, 0).intronic_length);
    let java = crate::compat::with_java(true, || crate::expression::GeneResult::new(&gc, 0).intronic_length);
    assert_eq!(corrected, 98);
    assert_eq!(java, 99);
}

// ---- NovelJunctionsTest ----

fn transcript_with_bounds(trans_id: i32, name: &str, gene_id: &str, bounds: (i32, i32), exons: &[(i32, i32)]) -> &'static TranscriptData {
    let t = transcript(trans_id, name, gene_id, 1, exons, true);
    let mut t2 = t.clone();
    t2.trans_start = bounds.0;
    t2.trans_end = bounds.1;
    Box::leak(Box::new(t2))
}

fn gene_overlaps(gc: &GeneCollection, r: &mut Read) {
    let o = crate::read_trans::find_overlapping_regions(&gc.exon_regions, Some(&gc.genes[0].exon_regions), r);
    process_overlapping_regions(r, &gc.exon_regions, &o);
}

#[test]
fn alt_splice_junction_types() {
    use crate::novel::{AltSpliceJunctionFinder, AsjContext, AsjType};
    let s1 = ref_base_str_1();
    let mut refg = RefSource::mock(vec![("1", s1.repeat(50))]);
    let g = gene("ENSG0001", "ENSG0001", "1", 1, 100, 1500);
    let t1 = transcript_with_bounds(1, "TRANS01", "ENSG0001", (100, 300), &[(100, 300), (400, 500), (800, 1000), (1400, 1500)]);
    let t2 = transcript_with_bounds(2, "TRANS02", "ENSG0001", (100, 300), &[(250, 350), (400, 500), (600, 700)]);
    let t3 = transcript_with_bounds(3, "TRANS03", "ENSG0001", (100, 300), &[(100, 300)]);
    let t4 = transcript_with_bounds(4, "TRANS04", "ENSG0001", (360, 500), &[(360, 500)]);
    let gc = GeneCollection::new(0, vec![GeneReadData { gene: g, exon_regions: Vec::new(), transcripts: vec![t1, t2, t3, t4], has_unspliced_regions: false }]);
    let known = hmf_common::known_fusion::KnownFusionCache::new();
    let mut finder = AltSpliceJunctionFinder::new(true, None, &known);
    let ids = vec![1, 2, 3, 4];
    assert!(AltSpliceJunctionFinder::junction_matches_known(&gc, &[0], [300, 400], &ids));
    assert!(!AltSpliceJunctionFinder::junction_matches_known(&gc, &[0], [300, 401], &ids));
    assert!(!AltSpliceJunctionFinder::junction_matches_known(&gc, &[0], [350, 800], &ids));

    let mk = |pos: i32, cigar: &str| {
        let mut r = read(1, "1", pos, cigar, Some(&s1));
        gene_overlaps(&gc, &mut r);
        r
    };
    // known 5' to novel intronic 3'
    let r = mk(291, "10M59N10M");
    let a = finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap();
    assert_eq!(a.sj_type, AsjType::Novel3Prime);
    assert_eq!(a.sj, [300, 360]);
    assert_eq!(a.contexts, [AsjContext::SpliceJunc, AsjContext::Intronic]);
    assert!(a.sj_start_regions.iter().any(|x| gc.exon_regions[*x].has_trans_id(1)));
    assert!(a.sj_end_regions.is_empty());
    assert!(a.candidate_trans_ids.contains(&1));
    assert!(!a.candidate_trans_ids.contains(&4));
    let r = mk(491, "10M49N10M");
    let mut a = finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap();
    assert_eq!(a.sj_type, AsjType::Novel3Prime);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_START, 0), 0);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_END, 0), 50);
    // novel intronic 5' to a known 3'
    let r = mk(361, "10M29N10M");
    let a = finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap();
    assert_eq!(a.sj_type, AsjType::Novel5Prime);
    // intronic at both ends
    let r = mk(721, "10M39N10M");
    assert!(r.mapped_regions.is_empty());
    let mut a = finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap();
    assert_eq!(a.sj_type, AsjType::Intronic);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_START, 0), 30);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_END, 0), 30);
    // skipped exons
    let r = mk(291, "10M499N10M");
    let mut a = finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap();
    assert_eq!(a.sj_type, AsjType::SkippedExons);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_START, 0), 0);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_END, 0), 0);
    // novel intron
    let r = mk(841, "10M99N10M");
    let mut a = finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap();
    assert_eq!(a.sj_type, AsjType::NovelIntron);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_START, 0), -150);
    assert_eq!(a.calc_nearest_exon_boundary(&gc, SE_END, 0), -150);
    // two reads making a novel exon
    let r1 = mk(991, "10M199N10M");
    let r2 = mk(1291, "10M99N10M");
    let a1 = finder.create_from_read(&gc, &mut refg, &r1, &[1]).unwrap();
    let a2 = finder.create_from_read(&gc, &mut refg, &r2, &[1]).unwrap();
    finder.junctions = vec![a1, a2];
    finder.check_novel_exon(&gc, 0, 1);
    assert_eq!(finder.junctions[0].sj_type, AsjType::NovelExon);
    assert_eq!(finder.junctions[1].sj_type, AsjType::NovelExon);
    // junctions of different transcripts
    let r = mk(291, "10M299N10M");
    let tids = r.trans_class.keys();
    let a = finder.create_from_read(&gc, &mut refg, &r, &tids).unwrap();
    assert_eq!(a.sj_type, AsjType::MixedTrans);
    assert_eq!(a.sj, [300, 600]);
    assert_eq!(a.contexts, [AsjContext::SpliceJunc, AsjContext::SpliceJunc]);
    assert!(a.candidate_trans_ids.contains(&1) && a.candidate_trans_ids.contains(&2));
    // a circular exon, from a pair of supplementary-aligned reads
    let mut r1 = read(1, "1", 400, "5S20M", Some(&gen_bases(25)));
    let mut r2 = read(1, "1", 481, "20M5S", Some(&gen_bases(25)));
    r1.flags |= FLAG_FIRST;
    r2.flags |= FLAG_FIRST | crate::read::FLAG_SUPPLEMENTARY;
    r1.supp_data = crate::read::SuppData::extract("1,481,+,20M5S,255,0");
    r2.supp_data = crate::read::SuppData::extract("1,400,+,5S20M,255,0");
    gene_overlaps(&gc, &mut r1);
    gene_overlaps(&gc, &mut r2);
    let c = finder.create_from_reads(&gc, &r1, &r2, &[1]).unwrap();
    assert_eq!(c.sj_type, AsjType::Circular);
    assert_eq!(c.sj, [400, 500]);
    assert_eq!(c.contexts, [AsjContext::SpliceJunc, AsjContext::SpliceJunc]);
}

#[test]
fn retained_introns() {
    use crate::novel::RetainedIntronFinder;
    let s1 = ref_base_str_1();
    let s2 = s1.repeat(2);
    let g = gene("ENSG0001", "ENSG0001", "1", 1, 100, 1500);
    let t1 = transcript_with_bounds(1, "TRANS01", "ENSG0001", (100, 300), &[(100, 300), (400, 500), (800, 1000), (1400, 1500)]);
    let t2 = transcript_with_bounds(2, "TRANS02", "ENSG0001", (100, 300), &[(100, 310), (400, 500), (600, 700)]);
    let t3 = transcript_with_bounds(3, "TRANS03", "ENSG0001", (100, 310), &[(100, 300), (420, 500)]);
    let gc = GeneCollection::new(0, vec![GeneReadData { gene: g, exon_regions: Vec::new(), transcripts: vec![t1, t2, t3], has_unspliced_regions: false }]);
    let mut finder = RetainedIntronFinder { enabled: true, introns: Vec::new() };
    let mk = |pos: i32, cigar: &str, bases: &str| {
        let mut r = read(1, "1", pos, cigar, Some(bases));
        gene_overlaps(&gc, &mut r);
        r
    };
    let eval = |f: &mut RetainedIntronFinder, a: &Read, b: &Read| f.evaluate_fragment_reads(&gc, &[a, b]);
    // not an exon-intron read for every transcript
    eval(&mut finder, &mk(291, "20M", &s1), &mk(340, "20M", &s1));
    assert_eq!(finder.introns.len(), 0);
    eval(&mut finder, &mk(281, "40M", &s2), &mk(340, "20M", &s1));
    assert_eq!(finder.introns.len(), 1);
    let ri = &finder.introns[0];
    assert!(!ri.is_start);
    assert_eq!(ri.position, 310);
    assert_eq!(ri.regions.len(), 1);
    assert!(gc.exon_regions[ri.regions[0]].has_trans_id(2));
    eval(&mut finder, &mk(391, "40M", &s2), &mk(440, "20M", &s1));
    assert_eq!(finder.introns.len(), 2);
    let ri = &finder.introns[1];
    assert!(ri.is_start);
    assert_eq!(ri.position, 400);
    assert!(gc.exon_regions[ri.regions[0]].has_trans_id(1) && gc.exon_regions[ri.regions[0]].has_trans_id(2));
    // splice support from the other read updates the same retained intron
    eval(&mut finder, &mk(391, "40M", &s2), &mk(491, "10M99N10M", &s1));
    assert_eq!(finder.introns.len(), 2);
    assert_eq!(finder.introns[1].fragment_count, 2);
    assert_eq!(finder.introns[1].spliced_fragment_count, 1);
    // not the first or last exon of a transcript
    eval(&mut finder, &mk(91, "20M", &s1), &mk(121, "20M", &s1));
    eval(&mut finder, &mk(1491, "20M", &s1), &mk(1551, "20M", &s1));
    assert_eq!(finder.introns.len(), 2);
    finder.introns.clear();
    // both reads over the boundaries of one exon: not allowed
    eval(&mut finder, &mk(391, "20M", &s1), &mk(491, "20M", &s1));
    assert!(finder.introns.is_empty());
    // the boundaries of consecutive exons are
    eval(&mut finder, &mk(491, "20M", &s1), &mk(591, "20M", &s1));
    assert_eq!(finder.introns.len(), 2);
    assert_eq!((finder.introns[0].is_start, finder.introns[0].position, finder.introns[0].fragment_count), (false, 500, 1));
    assert_eq!((finder.introns[1].is_start, finder.introns[1].position, finder.introns[1].fragment_count), (true, 600, 1));
}

// ---- FusionFragmentsTest (TestUtils' genes and transcripts) ----

/// TestUtils.addTestGenes / addTestTranscripts: (gene, transcript, non-genic boundaries)
fn fusion_test_genes() -> Vec<(&'static GeneData, &'static TranscriptData, [i32; 2])> {
    let mk = |id: &str, chr: &str, strand: i8, start: i32, len: i32, trans: i32, exons: usize, bounds: [i32; 2]| {
        let g = gene(id, &id.replace("ENSG000", "GENE"), chr, strand, start, start + len);
        let ex: Vec<(i32, i32)> = (0..exons).map(|i| (start + i as i32 * 200, start + i as i32 * 200 + 100)).collect();
        let t = transcript(trans, &format!("TRAN{:04}", trans), id, strand, &ex, true);
        (leak_gene(g), t, bounds)
    };
    vec![
        mk("ENSG0001", "1", 1, 1000, 500, 1, 3, [0, 9999]),
        mk("ENSG0002", "1", 1, 10000, 900, 2, 5, [1501, 19999]),
        mk("ENSG0003", "1", -1, 20000, 500, 3, 3, [10901, 29999]),
        mk("ENSG0004", "2", 1, 1000, 500, 4, 3, [0, 9999]),
        mk("ENSG0005", "2", -1, 10000, 900, 5, 5, [1501, 19999]),
        mk("ENSG0006", "2", -1, 10400, 1100, 6, 6, [10901, 29999]),
    ]
}

fn leak_gene(g: &'static GeneData) -> &'static GeneData {
    g
}

fn fusion_gc(id: i32, g: &(&'static GeneData, &'static TranscriptData, [i32; 2])) -> GeneCollection {
    let mut gc = GeneCollection::new(id, vec![GeneReadData { gene: g.0, exon_regions: Vec::new(), transcripts: vec![g.1], has_unspliced_regions: false }]);
    gc.non_genic_positions = g.2;
    gc
}

/// TestUtils.createMappedRead
fn mapped_read(id: i32, gc: &GeneCollection, pos: i32, cigar: &str, bases: Option<&str>) -> Read {
    let mut r = read(id, gc.chromosome, pos, cigar, bases);
    let o = crate::read_trans::find_overlapping_regions(&gc.exon_regions, None, &r);
    process_overlapping_regions(&mut r, &gc.exon_regions, &o);
    if r.mapped_regions.is_empty() {
        let mut refs = Vec::new();
        for t in &gc.transcripts {
            if !r.mapped_coords.alignments_within(t.trans_start, t.trans_end) {
                continue;
            }
            for i in 0..t.exons.len() - 1 {
                if r.mapped_coords.alignments_within(t.exons[i].end, t.exons[i + 1].start) {
                    refs.push(TransExonRef::new(t, t.exons[i].rank.min(t.exons[i + 1].rank)));
                    break;
                }
            }
        }
        if !refs.is_empty() {
            r.intron_trans_refs = Some(refs);
        }
    }
    gc.set_read_gene_collections(&mut r, gc.non_genic_positions);
    r
}

fn supp_tag(r: &Read) -> String {
    format!("{},{},{},{},{},0", r.chromosome, r.pos_start, if r.is_read_reversed() { '-' } else { '+' }, r.cigar_str(), r.mapq)
}

/// TestUtils.createSupplementaryReadPair: the second read is the supplementary
fn supp_pair(id: i32, gc1: &GeneCollection, gc2: &GeneCollection, pos1: i32, pos2: i32, c1: &str, c2: &str, first: bool) -> (Read, Read) {
    let len: u32 = parse_cigar(c1).unwrap().iter().filter(|e| e.op != crate::read::OP_N && e.op != crate::read::OP_D).map(|e| e.len).sum();
    let bases = gen_bases(len as usize);
    let mut r1 = mapped_read(id, gc1, pos1, c1, Some(&bases));
    let mut r2 = mapped_read(id, gc2, pos2, c2, Some(&bases));
    if first {
        r1.flags |= FLAG_FIRST;
        r2.flags |= FLAG_FIRST;
    }
    r2.flags |= crate::read::FLAG_SUPPLEMENTARY;
    let (t1, t2) = (supp_tag(&r2), supp_tag(&r1));
    r1.supp_data = crate::read::SuppData::extract(&t1);
    r2.supp_data = crate::read::SuppData::extract(&t2);
    (r1, r2)
}

fn set_strand(r: &mut Read, rev: bool, mate_rev: bool) {
    r.flags &= !(FLAG_REVERSE | crate::read::FLAG_MATE_REVERSE);
    if rev {
        r.flags |= FLAG_REVERSE;
    }
    if mate_rev {
        r.flags |= crate::read::FLAG_MATE_REVERSE;
    }
}

fn from_reads(reads: &[(&Read, &GeneCollection)]) -> crate::fusion_frag::FusionFragment {
    use crate::fusion_frag::{FusionFragment, FusionRead, FusionReadGroup};
    let fr: Vec<FusionRead> = reads.iter().map(|(r, gc)| FusionRead::new(r, &gc.exon_regions)).collect();
    FusionFragment::new(FusionReadGroup::new(&reads[0].0.id, fr))
}

fn fragment_gene_ids(f: &crate::fusion_frag::FusionFragment, genes: &[(&'static GeneData, &'static TranscriptData, [i32; 2])], se: usize) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for te in &f.trans_exon_refs[se] {
        if let Some((g, _, _)) = genes.iter().find(|(_, t, _)| t.trans_id == te.trans_id) {
            if !out.contains(&g.gene_id) {
                out.push(g.gene_id.clone());
            }
        }
    }
    out
}

#[test]
fn del_fragments() {
    use crate::fusion_frag::FragType;
    let genes = fusion_test_genes();
    let gc1 = fusion_gc(0, &genes[0]);
    let gc2 = fusion_gc(1, &genes[1]);
    // a simple DEL from supplementary alignments
    let mut r1 = mapped_read(0, &gc1, 1050, "40M", None);
    r1.flags |= FLAG_FIRST;
    let (a, mut b) = supp_pair(0, &gc1, &gc2, 1081, 10200, "20M20S", "20S20M", false);
    set_strand(&mut b, true, false);
    let f = from_reads(&[(&r1, &gc1), (&a, &gc1), (&b, &gc2)]);
    assert_eq!(f.frag_type, FragType::MatchedJunction);
    assert_eq!(f.chromosomes, ["1", "1"]);
    assert_eq!(f.gene_collections, [0, 1]);
    assert_eq!(f.junction_positions, [1100, 10200]);
    assert_eq!(f.orientations, [1, -1]);
    assert_eq!(f.junction_orientations, [1, -1]);
    assert_eq!(f.implied_sv_type(), "DEL");
    assert!(f.is_spliced());
    assert_eq!(fragment_gene_ids(&f, &genes, SE_START), vec!["ENSG0001"]);
    assert_eq!(fragment_gene_ids(&f, &genes, SE_END), vec!["ENSG0002"]);
    // an unspliced DEL
    let mut r1 = mapped_read(1, &gc1, 1110, "40M", None);
    r1.flags |= FLAG_FIRST;
    let (a, mut b) = supp_pair(1, &gc1, &gc2, 1131, 10150, "20M20S", "20S20M", false);
    set_strand(&mut b, true, false);
    let f = from_reads(&[(&r1, &gc1), (&a, &gc1), (&b, &gc2)]);
    assert_eq!(f.frag_type, FragType::MatchedJunction);
    assert_eq!(f.junction_positions, [1150, 10150]);
    assert_eq!(f.implied_sv_type(), "DEL");
    assert!(f.is_unspliced());
    assert_eq!(fragment_gene_ids(&f, &genes, SE_START), vec!["ENSG0001"]);
    assert_eq!(fragment_gene_ids(&f, &genes, SE_END), vec!["ENSG0002"]);
    // a short DEL without supplementary data: soft-clipped facing reads
    let r1 = mapped_read(1, &gc1, 1081, "20M20S", None);
    let mut r2 = mapped_read(1, &gc2, 10200, "20S20M", None);
    set_strand(&mut r2, true, false);
    let f = from_reads(&[(&r1, &gc1), (&r2, &gc2)]);
    assert_eq!(f.frag_type, FragType::DiscordantJunction);
    assert_eq!(f.gene_collections, [0, 1]);
    assert_eq!(f.orientations, [1, -1]);
    assert_eq!(f.junction_positions, [1100, 10200]);
    assert_eq!(f.junction_orientations, [1, -1]);
    assert!(f.is_spliced());
    // a DEL spanning the two gene collections in the reads' N-splits
    let mut r1 = mapped_read(1, &gc1, 1080, "21M9099N19M", None);
    let mut r2 = mapped_read(1, &gc1, 1081, "20M9099N20M", None);
    r1.flags |= FLAG_FIRST;
    set_strand(&mut r2, true, false);
    for r in [&mut r1, &mut r2] {
        r.set_gene_collection(SE_START, 0, true);
        r.set_gene_collection(SE_END, 1, true);
    }
    // the junction comes from the split read's recorded junction (ChimericReadTracker.addJunction)
    r1.set_junction_position(SE_START, 1100);
    r1.set_junction_position(SE_END, 10200);
    let f = from_reads(&[(&r1, &gc1), (&r2, &gc1)]);
    assert_eq!(f.frag_type, FragType::MatchedJunction);
    assert_eq!(f.gene_collections, [0, 1]);
    assert_eq!(f.junction_positions, [1100, 10200]);
    assert_eq!(f.junction_orientations, [1, -1]);
    assert_eq!(f.implied_sv_type(), "DEL");
}

#[test]
fn dup_and_discordant_junction_order() {
    use crate::fusion_frag::FragType;
    let genes = fusion_test_genes();
    let gc1 = fusion_gc(0, &genes[0]);
    let gc2 = fusion_gc(1, &genes[1]);
    // a DUP between the genes
    let (mut a, mut b) = supp_pair(0, &gc2, &gc1, 10281, 1200, "20M20S", "20S20M", true);
    let r1 = mapped_read(1, &gc2, 10220, "40M", None);
    set_strand(&mut a, true, false);
    set_strand(&mut b, true, false);
    let f = from_reads(&[(&r1, &gc2), (&a, &gc2), (&b, &gc1)]);
    assert_eq!(f.frag_type, FragType::MatchedJunction);
    assert_eq!(f.junction_positions, [1200, 10300]);
    assert_eq!(f.junction_orientations, [-1, 1]);
    assert_eq!(f.implied_sv_type(), "DUP");
    // IF9: soft-clipped facing reads given upper read first; Java swaps the junctions
    let r1 = mapped_read(2, &gc1, 1081, "20M20S", None);
    let mut r2 = mapped_read(2, &gc2, 10200, "20S20M", None);
    set_strand(&mut r2, true, false);
    let corrected = crate::compat::with_java(false, || from_reads(&[(&r2, &gc2), (&r1, &gc1)]));
    assert_eq!(corrected.frag_type, FragType::DiscordantJunction);
    assert_eq!(corrected.junction_positions, [1100, 10200]);
    assert_eq!(corrected.junction_orientations, [1, -1]);
    let java = crate::compat::with_java(true, || from_reads(&[(&r2, &gc2), (&r1, &gc1)]));
    assert_eq!(java.junction_positions, [10200, 1100]);
    assert_eq!(java.junction_orientations, [-1, 1]);
}

// ---- FusionDataTest ----

fn fusion_cache(genes: &[(&'static GeneData, &'static TranscriptData, [i32; 2])]) -> &'static hmf_common::ensembl::EnsemblDataCache {
    let mut c = hmf_common::ensembl::EnsemblDataCache::new("", hmf_common::regions::RefGenomeVersion::V37);
    for (g, t, _) in genes {
        c.chr_gene_data_map_mut().compute_if_absent(g.chromosome.clone(), Vec::new).push((*g).clone());
        c.transcript_data_map_mut().insert(g.gene_id.clone(), vec![(*t).clone()]);
    }
    Box::leak(Box::new(c))
}

fn group(reads: &[(&Read, &GeneCollection)]) -> (String, crate::fusion_frag::FusionReadGroup) {
    use crate::fusion_frag::{FusionRead, FusionReadGroup};
    let fr: Vec<FusionRead> = reads.iter().map(|(r, gc)| FusionRead::new(r, &gc.exon_regions)).collect();
    (reads[0].0.id.to_string(), FusionReadGroup::new(&reads[0].0.id, fr))
}

#[test]
fn fusion_fragment_assignment() {
    use crate::fusion::{CohortFusions, FusionContext, FusionFinder};
    use crate::task::{process_new_chimeric_read_groups, GroupArena};
    let genes = fusion_test_genes();
    let cache = fusion_cache(&genes);
    let trans_by_id: HashMap<i32, &'static TranscriptData> =
        cache.transcript_data_map().iter().flat_map(|(_, v)| v.iter()).map(|t| (t.trans_id, t)).collect();
    let known = hmf_common::known_fusion::KnownFusionCache::new();
    let filters = crate::config::GeneRegionFilters::default();
    let cohort = CohortFusions::default();
    let ctx = FusionContext { cache, trans_by_id: &trans_by_id, known: &known, filters: &filters, max_fragment_length: 550, cohort: &cohort, min_hard_filter_frags: 0 };
    let gc1 = fusion_gc(0, &genes[0]);
    let gc2 = fusion_gc(1, &genes[1]);
    let mut g1: Vec<(String, crate::fusion_frag::FusionReadGroup)> = Vec::new();
    let mut g2: Vec<(String, crate::fusion_frag::FusionReadGroup)> = Vec::new();
    // two spliced fragments
    let r1 = mapped_read(0, &gc1, 1050, "40M", None);
    let (a, mut b) = supp_pair(0, &gc1, &gc2, 1081, 10200, "20M20S", "20S20M", true);
    set_strand(&mut b, true, false);
    g1.push(group(&[(&r1, &gc1), (&a, &gc1)]));
    g2.push(group(&[(&b, &gc2)]));
    let (a, mut b) = supp_pair(1, &gc1, &gc2, 1081, 10200, "20M20S", "20S20M", true);
    let mut r3 = mapped_read(1, &gc2, 10210, "40M", None);
    set_strand(&mut b, true, false);
    set_strand(&mut r3, true, false);
    g1.push(group(&[(&a, &gc1)]));
    g2.push(group(&[(&b, &gc2), (&r3, &gc2)]));
    // one unspliced fragment, its own fusion
    let r1 = mapped_read(2, &gc1, 1110, "40M", None);
    let (a, mut b) = supp_pair(2, &gc1, &gc2, 1131, 10150, "20M20S", "20S20M", true);
    set_strand(&mut b, true, false);
    g1.push(group(&[(&r1, &gc1), (&a, &gc1)]));
    g2.push(group(&[(&b, &gc2)]));
    // a discordant fragment for the spliced fusion, and one for the unspliced
    let r1 = mapped_read(3, &gc1, 1055, "40M", None);
    let mut r2 = mapped_read(3, &gc2, 10220, "40M", None);
    set_strand(&mut r2, true, false);
    g1.push(group(&[(&r1, &gc1)]));
    g2.push(group(&[(&r2, &gc2)]));
    let r1 = mapped_read(4, &gc1, 1110, "40M", None);
    let mut r2 = mapped_read(4, &gc2, 10160, "40M", None);
    set_strand(&mut r2, true, false);
    g1.push(group(&[(&r1, &gc1)]));
    g2.push(group(&[(&r2, &gc2)]));

    // createIsofoxConfig's MockRefGenome holds no bases
    let mut refg = RefSource::mock(vec![]);
    let depth = BaseDepth::new();
    let mut arena = GroupArena::default();
    let mut finder = FusionFinder::new("1");
    let mut no_rac = |_: &str, _: i32| -> Option<*mut crate::fusion_frag::JunctionRacFragments> { None };
    let complete = process_new_chimeric_read_groups(&mut arena, &gc1, &depth, g1);
    assert!(complete.is_empty());
    assert_eq!(arena.partials.len(), 5);
    let complete = process_new_chimeric_read_groups(&mut arena, &gc2, &depth, g2);
    assert!(arena.partials.is_empty());
    let groups: Vec<_> = complete.iter().map(|g| arena.snapshot(*g)).collect();
    assert_eq!(groups.len(), 5);
    finder.process_read_groups(&ctx, &mut refg, groups, &mut no_rac);
    assert_eq!(finder.outputs.len(), 1);
    let fusions = &finder.outputs[0].fusions;
    assert_eq!(fusions.len(), 2);
    let spliced = fusions.iter().find(|f| f.junction_positions == [1100, 10200]).unwrap();
    assert_eq!(spliced.split_frags, 2);
    assert_eq!(spliced.discordant_frags, 1);
    assert_eq!(spliced.anchor_distance, [20, 20]);
    assert_eq!(spliced.gene_names, ["GENE1".to_string(), "GENE2".to_string()]);
    assert_eq!(spliced.sv_type, "DEL");
    let unspliced = fusions.iter().find(|f| f.junction_positions == [1150, 10150]).unwrap();
    assert_eq!(unspliced.split_frags, 1);
    // both discordant fragments fit the unspliced fusion (FusionDataTest asserts two)
    assert_eq!(unspliced.discordant_frags, 2);
}

// ---- findings in both modes ----

fn temp_file(name: &str, content: &str) -> String {
    let p = std::env::temp_dir().join(format!("isofox_rs_{}_{}", std::process::id(), name));
    std::fs::write(&p, content).unwrap();
    p.to_string_lossy().to_string()
}

#[test]
fn expected_counts_frequencies_applied_once() {
    // IF1: Java scales the cached counts in place, so a gene set used again is scaled again
    let path = temp_file("if1_exp_counts.csv", "GeneSetId,Category,Length_100,Length_200\nGS1,1,10,20\nGS1,ENSG0001,4,8\n");
    let mut config = IsofoxConfig::for_tests();
    config.exp_counts_file = Some(path.clone());
    config.fragment_size_data = vec![FragmentSize { length: 100, frequency: 1 }, FragmentSize { length: 200, frequency: 1 }];
    let cache = crate::expression::ExpectedCountsCache::new(&mut config).unwrap();
    assert!(cache.valid);
    let rates = [0.5, 2.0];
    let lookup = |java: bool| {
        let c = crate::compat::with_java(java, || cache.gene_expected_rates_data("GS1", &["ENSG0001"], Some(&rates)).unwrap());
        c[0].counts_by_length.clone().unwrap()
    };
    assert_eq!(lookup(true), vec![5.0, 40.0]);
    // the second use: Java applies the frequencies twice
    assert_eq!(lookup(true), vec![2.5, 80.0]);
    // corrected, every use starts from the file's counts
    assert_eq!(lookup(false), vec![5.0, 40.0]);
    std::fs::remove_file(path).ok();
}

#[test]
fn gc_counts_skip_only_the_missing_transcript() {
    // IF4: a transcript missing from the expected GC ratios ends Java's counts for the rest of the chromosome
    let n = GcRatioCounts::new().size();
    let row = |name: &str, bucket: usize| {
        let mut v = vec!["0".to_string(); n];
        v[bucket] = "1".into();
        format!("{name},{}", v.join(","))
    };
    let header = (0..n).map(|i| format!("Gc{i}")).collect::<Vec<_>>().join(",");
    let path = temp_file("if4_gc_ratios.csv", &format!("TransName,{header}\n{}\n", row("T1", 40)));
    let mut config = IsofoxConfig::for_tests();
    config.exp_gc_ratios_file = Some(path.clone());
    let run = |java: bool| {
        let mut calc = crate::adjusts::GcTranscriptCalculator::new(&config).unwrap();
        let mut a1 = hmf_common::java_hash::JavaHashMap::new();
        a1.insert("T_MISSING".to_string(), 1.0);
        let mut a2 = hmf_common::java_hash::JavaHashMap::new();
        a2.insert("T1".to_string(), 2.0);
        crate::compat::with_java(java, || calc.generate_gc_counts_from_fit(&[("G1", &a1), ("G2", &a2)]));
        calc.transcript_fit_gc_counts.counts[40]
    };
    assert_eq!(run(true), 0.0);
    assert_eq!(run(false), 2.0);
    std::fs::remove_file(path).ok();
}

#[test]
fn junction_homology_shift() {
    // IF5: a junction within 5 bases of an exon boundary moves onto it only if the bases moved match;
    // Java compares the bases with themselves and always moves it
    use crate::novel::AltSpliceJunctionFinder;
    let bases = "ACGT".repeat(500);
    let mut refg = RefSource::mock(vec![("1", bases)]);
    let g = gene("ENSG0001", "ENSG0001", "1", 1, 100, 1500);
    let t1 = transcript_with_bounds(1, "TRANS01", "ENSG0001", (100, 1500), &[(100, 300), (400, 500), (800, 1000), (1400, 1500)]);
    let gc = GeneCollection::new(0, vec![GeneReadData { gene: g, exon_regions: Vec::new(), transcripts: vec![t1], has_unspliced_regions: false }]);
    let known = hmf_common::known_fusion::KnownFusionCache::new();
    let finder = AltSpliceJunctionFinder::new(true, None, &known);
    let ids = vec![1];
    let mk = |pos: i32, cigar: &str| {
        let mut r = read(1, "1", pos, cigar, None);
        gene_overlaps(&gc, &mut r);
        r
    };
    // 2 bases short of the exon end at 300; the bases after each side differ (TA against AC)
    let r = mk(289, "10M81N10M");
    let java = crate::compat::with_java(true, || finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap().sj);
    let corrected = crate::compat::with_java(false, || finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap().sj);
    assert_eq!(java, [300, 382]);
    assert_eq!(corrected, [298, 380]);
    // the same with matching bases (TA and TA): both move it
    let r = mk(289, "10M84N10M");
    let java = crate::compat::with_java(true, || finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap().sj);
    let corrected = crate::compat::with_java(false, || finder.create_from_read(&gc, &mut refg, &r, &ids).unwrap().sj);
    assert_eq!(java, [300, 385]);
    assert_eq!(corrected, [300, 385]);
}
