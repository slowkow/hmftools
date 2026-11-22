//! The Java unit tests of the package (DepthAnnotationTest), with the same cases and assertions.

use super::config::DepthConfig;
use super::task::*;
use super::vcf::{Genotype, VariantContext};

const TEST_SAMPLE_ID: &str = "SAMPLE";
const CHR_1: &str = "1";
const CHR_2: &str = "2";

struct TestRead {
    name: String,
    info: ReadInfo,
    supp_data: String,
}

impl SamRead for TestRead {
    fn name(&self) -> &[u8] {
        self.name.as_bytes()
    }
    fn info(&self) -> ReadInfo {
        self.info
    }
    fn supplementary_data(&self) -> Option<&[u8]> {
        if self.supp_data.is_empty() { None } else { Some(self.supp_data.as_bytes()) }
    }
}

struct Fixture {
    task: DepthTask,
    next_variant_id: i32,
    next_read_id: i32,
}

impl Fixture {
    fn new() -> Fixture {
        let mut config = DepthConfig::for_test(0.1, 1000);
        config.sample_ids.push(TEST_SAMPLE_ID.to_string());
        Fixture { task: DepthTask::new(CHR_1, &config, vec![0]), next_variant_id: 0, next_read_id: 0 }
    }

    fn next_variant_id(&mut self) -> String {
        let id = format!("{:03}", self.next_variant_id);
        self.next_variant_id += 1;
        id
    }

    /// ReadIdGenerator.nextId
    fn next_read_id(&mut self) -> String {
        let id = format!("READ_{:03}", self.next_read_id);
        self.next_read_id += 1;
        id
    }

    fn check_ref_support(&self, var_position: i32, ref_support: i32, ref_pair_support: i32) {
        let variant = self.task.variant_infos().iter().find(|x| x.position == var_position);
        assert!(variant.is_some());
        let variant = variant.unwrap();
        assert_eq!(ref_support, variant.sample_support_counts[0].ref_support);
        assert_eq!(ref_pair_support, variant.sample_support_counts[0].ref_pair_support);
    }
}

fn chromosome_to_index(chromosome: &str) -> i32 {
    chromosome.parse::<i32>().unwrap() - 1
}

/// (reference length, read length, starts with a soft clip, ends with one)
fn cigar_lengths(cigar: &str) -> (i32, i32, bool, bool) {
    let mut ref_length = 0;
    let mut read_length = 0;
    let mut ops: Vec<char> = Vec::new();
    let mut number = 0;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            number = number * 10 + d as i32;
        } else {
            if matches!(c, 'M' | 'D' | 'N' | '=' | 'X') {
                ref_length += number;
            }
            if matches!(c, 'M' | 'I' | 'S' | '=' | 'X') {
                read_length += number;
            }
            ops.push(c);
            number = 0;
        }
    }
    (ref_length, read_length, ops.first() == Some(&'S'), ops.last() == Some(&'S'))
}

/// DepthAnnotationTest.createSamRecord: a paired, properly paired read with no insert size set
#[allow(clippy::too_many_arguments)]
fn create_sam_record(
    read_id: &str, chromosome: &str, read_start: i32, mate_chr: &str, mate_start: i32, cigar: &str, first_in_pair: bool, is_reversed: bool,
    is_supplementary: bool, supp_data: &str,
) -> TestRead {
    let (ref_length, read_length, left_soft_clipped, right_soft_clipped) = cigar_lengths(cigar);
    let mut flags = FLAG_PAIRED | FLAG_PROPER_PAIR;
    if is_reversed {
        flags |= FLAG_REVERSE;
    }
    flags |= if first_in_pair { FLAG_FIRST_OF_PAIR } else { FLAG_SECOND_OF_PAIR };
    if is_supplementary {
        flags |= FLAG_SUPPLEMENTARY;
    }
    TestRead {
        name: read_id.to_string(),
        info: ReadInfo {
            start: read_start,
            end: read_start + ref_length - 1,
            flags,
            ref_index: chromosome_to_index(chromosome),
            mate_ref_index: chromosome_to_index(mate_chr),
            mate_start,
            insert_size: 0,
            read_length,
            left_soft_clipped,
            right_soft_clipped,
        },
        supp_data: supp_data.to_string(),
    }
}

/// DepthAnnotationTest.createVariantContext: a breakend with one supporting fragment
fn create_variant_context(vcf_id: &str, position: i32) -> VariantContext {
    let mut genotype = Genotype::default();
    genotype.attributes.insert(TOTAL_FRAGS.to_string(), "1".to_string());
    let mut variant = VariantContext {
        chromosome: CHR_1.to_string(),
        position,
        id: vcf_id.to_string(),
        ref_: "A".to_string(),
        alts: vec!["T[2:500[".to_string()],
        log10_p_error: Some(-(37.0 / 10.0)),
        filters: None,
        genotypes: vec![genotype],
        ..Default::default()
    };
    variant.info.insert(TOTAL_FRAGS.to_string(), "1".to_string());
    variant.info.insert("SVTYPE".to_string(), "BND".to_string());
    variant
}

#[test]
fn test_reference_depth() {
    let mut f = Fixture::new();
    f.task.reset();

    let variants: Vec<VariantContext> = [1000, 1100, 1200].iter().map(|pos| { let id = f.next_variant_id(); create_variant_context(&id, *pos) }).collect();
    f.task.add_variants(variants).unwrap();
    f.task.add_slice_variants();

    // irrelevant
    let read1 = create_sam_record(&f.next_read_id(), CHR_1, 900, CHR_2, 100, "100M", true, false, false, "");
    let read2 = create_sam_record(&f.next_read_id(), CHR_1, 1250, CHR_2, 100, "100M", true, false, false, "");

    f.task.process_read(&read1).unwrap();
    f.task.process_read(&read2).unwrap();
    assert_eq!(0, f.task.read_group_count());

    f.check_ref_support(1000, 0, 0);

    // reads processed immediately for overlapping overage
    let read3 = create_sam_record(&f.next_read_id(), CHR_1, 950, CHR_2, 100, "100M", true, false, false, "");
    f.task.process_read(&read3).unwrap();

    assert_eq!(0, f.task.read_group_count());

    f.check_ref_support(1000, 1, 0);
    f.check_ref_support(1100, 0, 0);

    // a read straddling all three variants
    let read4 = create_sam_record(&f.next_read_id(), CHR_1, 800, CHR_1, 1300, "100M", true, false, false, "");
    let mut mate4 = create_sam_record(&read4.name, CHR_1, 1300, CHR_1, 800, "100M", false, true, false, "");
    mate4.info.flags |= FLAG_REVERSE;

    f.task.process_read(&read4).unwrap();

    assert_eq!(1, f.task.read_group_count());

    f.task.process_read(&mate4).unwrap();
    assert_eq!(0, f.task.read_group_count());

    f.check_ref_support(1000, 1, 1);
    f.check_ref_support(1100, 0, 1);
    f.check_ref_support(1200, 0, 1);

    // confirm that junction-matching reads are ignored
    let read5 = create_sam_record(&f.next_read_id(), CHR_1, 1100, CHR_2, 100, "5S145M", true, false, false, "");
    f.task.process_read(&read5).unwrap();

    f.check_ref_support(1100, 1, 1);
    f.check_ref_support(1200, 1, 1);

    let read6 = create_sam_record(&f.next_read_id(), CHR_1, 960, CHR_2, 100, "141M10S", true, false, false, "");
    f.task.process_read(&read6).unwrap();

    f.check_ref_support(1000, 2, 1);
    f.check_ref_support(1100, 1, 1);
    f.check_ref_support(1200, 1, 1);

    // test max fragment caps
    for _ in 0..15 {
        let read = create_sam_record(&f.next_read_id(), CHR_1, 1150, CHR_2, 100, "100M", true, false, false, "");
        f.task.process_read(&read).unwrap();
    }

    f.check_ref_support(1200, 9, 1);
    assert_eq!(2, f.task.slice_region_state().uncapped_variants.len());
}

#[test]
fn test_position_progression() {
    let mut f = Fixture::new();
    f.task.reset();

    let start_position = 1000;
    let pos_gap = 200;
    let mut variants = Vec::new();
    for i in 0..20 {
        let id = f.next_variant_id();
        variants.push(create_variant_context(&id, start_position + i * pos_gap));
    }
    let variant_count = variants.len();

    f.task.add_variants(variants).unwrap();
    f.task.add_slice_variants();

    // check slice state as reads are received
    let read = create_sam_record(&f.next_read_id(), CHR_1, 950, CHR_2, 100, "100M", true, false, false, "");
    f.task.process_read(&read).unwrap();

    assert_eq!(variant_count, f.task.slice_region_state().uncapped_variants.len());
    assert_eq!(0, f.task.slice_region_state().min_position_index);

    let read = create_sam_record(&f.next_read_id(), CHR_1, 1150, CHR_2, 100, "100M", true, false, false, "");
    f.task.process_read(&read).unwrap();

    assert_eq!(1, f.task.slice_region_state().min_position_index);

    // process ahead and then cap out some variants
    let read = create_sam_record(&f.next_read_id(), CHR_1, 2150, CHR_2, 100, "100M", true, false, false, "");
    f.task.process_read(&read).unwrap();

    assert_eq!(6, f.task.slice_region_state().min_position_index);

    for _ in 0..15 {
        let read2 = create_sam_record(&f.next_read_id(), CHR_1, 2150, CHR_2, 100, "100M", true, false, false, "");
        f.task.process_read(&read2).unwrap();
    }

    assert_eq!(5, f.task.slice_region_state().min_position_index);
    assert_eq!(variant_count - 1, f.task.slice_region_state().uncapped_variants.len());
}

// ---- tests of the port itself (not in the Java package) ----

use super::unmapped::UnmappedRegionDepth;
use anyhow::Result;

fn breakend(position: i32, alt: &str, sv_type: &str, cipos: Option<&str>, sample_frags: &[i32]) -> VariantContext {
    let mut variant = create_variant_context("id", position);
    variant.alts = vec![alt.to_string()];
    variant.info.insert("SVTYPE".to_string(), sv_type.to_string());
    if let Some(cipos) = cipos {
        variant.info.insert("CIPOS".to_string(), cipos.to_string());
    }
    variant.genotypes = sample_frags
        .iter()
        .map(|frags| {
            let mut genotype = Genotype::default();
            genotype.attributes.insert(TOTAL_FRAGS.to_string(), frags.to_string());
            genotype
        })
        .collect();
    variant
}

#[test]
fn variant_info_from_record() {
    let info = VariantInfo::new(&breakend(1000, "A[2:500[", "BND", Some("-3,7"), &[4]), &[0], 0.1).unwrap();
    assert_eq!((info.position, info.position_min, info.position_max), (1000, 997, 1007));
    assert_eq!(info.orient, Orientation::Forward);
    assert!(!info.is_sgl && !info.is_short_indel);
    assert_eq!(info.sample_support_counts[0].vaf_cap, 40);

    let info = VariantInfo::new(&breakend(1000, "]2:500]A", "BND", None, &[0]), &[0], 0.001).unwrap();
    assert_eq!(info.orient, Orientation::Reverse);
    // no fragments counts as one for the cap; no cap when the option is zero
    assert_eq!(info.sample_support_counts[0].vaf_cap, 1000);
    assert_eq!(VariantInfo::new(&breakend(1000, "]2:500]A", "BND", None, &[9]), &[0], 0.0).unwrap().sample_support_counts[0].vaf_cap, 0);

    // single breakends
    let info = VariantInfo::new(&breakend(1000, "AGG.", "SGL", None, &[1]), &[0], 0.1).unwrap();
    assert!(info.is_sgl && info.orient == Orientation::Forward);
    let info = VariantInfo::new(&breakend(1000, ".GGA", "SGL", None, &[1]), &[0], 0.1).unwrap();
    assert!(info.is_sgl && info.orient == Orientation::Reverse);

    // a local deletion, duplication or insertion of up to 1000 bases is a short indel
    assert!(VariantInfo::new(&breakend(1000, "A[1:2000[", "DEL", None, &[1]), &[0], 0.1).unwrap().is_short_indel);
    assert!(!VariantInfo::new(&breakend(1000, "A[1:2001[", "DEL", None, &[1]), &[0], 0.1).unwrap().is_short_indel);
    assert!(VariantInfo::new(&breakend(1000, "]1:400]A", "DUP", None, &[1]), &[0], 0.1).unwrap().is_short_indel);
    assert!(!VariantInfo::new(&breakend(1000, "A]1:1010]", "INV", None, &[1]), &[0], 0.1).unwrap().is_short_indel);

    // where Java throws
    assert!(VariantInfo::new(&breakend(1000, "ACGT", "BND", None, &[1]), &[0], 0.1).is_err());
    assert!(VariantInfo::new(&breakend(1000, "A[2:500[", "OTHER", None, &[1]), &[0], 0.1).is_err());
    assert!(VariantInfo::new(&breakend(1000, "A[2:500[", "BND", Some("5"), &[1]), &[0], 0.1).is_err());
    let mut no_alt = breakend(1000, "A[2:500[", "BND", None, &[1]);
    no_alt.alts.clear();
    assert!(VariantInfo::new(&no_alt, &[0], 0.1).is_err());
}

fn port_config(samples: &[&str]) -> DepthConfig {
    let mut config = DepthConfig::for_test(0.1, 1000);
    config.vcf_tag_prefix = None;
    config.sample_ids = samples.iter().map(|x| x.to_string()).collect();
    config
}

#[test]
fn groups_by_proximity() {
    let mut task = DepthTask::new(CHR_1, &port_config(&[TEST_SAMPLE_ID]), vec![0]);
    let positions = [1000, 1500, 2500, 3501, 3600, 9000];
    task.add_variants(positions.iter().map(|p| create_variant_context("x", *p)).collect()).unwrap();
    // a variant joins while it is within 1000 of the one before it
    assert_eq!(task.group_ranges(), vec![0..3, 3..5, 5..6]);
}

/// Reads held in memory, as one sample's alignment file
struct TestSource {
    reads: Vec<TestRead>,
}

impl ReadSource for TestSource {
    fn slice(&mut self, chromosome: &str, start: i32, end: i32, consumer: &mut dyn FnMut(&dyn SamRead) -> Result<bool>) -> Result<()> {
        let ref_index = chromosome_to_index(chromosome);
        for read in &self.reads {
            if read.info.ref_index == ref_index && read.info.end >= start && read.info.start <= end && !consumer(read)? {
                break;
            }
        }
        Ok(())
    }
}

fn pair(name: &str, start: i32, mate_start: i32) -> [TestRead; 2] {
    let mut first = create_sam_record(name, CHR_1, start, CHR_1, mate_start, "100M", true, false, false, "");
    let mut second = create_sam_record(name, CHR_1, mate_start, CHR_1, start, "100M", false, true, false, "");
    first.info.insert_size = mate_start + 100 - start;
    second.info.insert_size = -first.info.insert_size;
    [first, second]
}

#[test]
fn slices_and_sets_depth_values() {
    let config = port_config(&["S1", "S2"]);
    // genotype columns in the opposite order to the samples
    let mut task = DepthTask::new(CHR_1, &config, vec![1, 0]);
    task.add_variants(vec![breakend(1000, "A[2:500[", "BND", None, &[3, 2]), breakend(5000, "A.", "SGL", None, &[1, 2])]).unwrap();

    let mut first_sample: Vec<TestRead> = Vec::new();
    // over the first breakend, a fragment either side of it, and a read clipped at it
    first_sample.push(create_sam_record("a", CHR_1, 950, CHR_2, 100, "100M", true, false, false, ""));
    first_sample.extend(pair("b", 700, 1200));
    first_sample.push(create_sam_record("c", CHR_1, 931, CHR_2, 100, "70M30S", true, false, false, ""));
    // over the single breakend, and a fragment either side of it
    first_sample.push(create_sam_record("d", CHR_1, 4950, CHR_2, 100, "100M", true, false, false, ""));
    first_sample.extend(pair("e", 4700, 5200));
    first_sample.sort_by_key(|r| r.info.start);
    let second_sample = vec![create_sam_record("f", CHR_1, 990, CHR_2, 100, "100M", true, false, false, "")];

    let mut sources: Vec<Box<dyn ReadSource>> = vec![Box::new(TestSource { reads: first_sample }), Box::new(TestSource { reads: second_sample })];
    let ranges = task.group_ranges();
    assert_eq!(ranges, vec![0..1, 1..2]);
    let mut slicer = DepthTask::new("", &config, Vec::new());
    let mut total_reads = 0;
    for range in ranges {
        let (reads, _) = slicer.slice_group(CHR_1, &mut task.infos_mut()[range], &mut sources).unwrap();
        total_reads += reads;
    }
    assert_eq!(total_reads, 8);
    task.set_depth_values().unwrap();

    let first = &task.variants[0];
    // sample S1 is genotype 1
    assert_eq!(first.genotypes[1].attributes["REF"], "1");
    assert_eq!(first.genotypes[1].attributes["REFPAIR"], "1");
    assert_eq!(first.genotypes[1].attributes["AF"], "0.500");
    assert_eq!(first.genotypes[0].attributes["REF"], "1");
    assert_eq!(first.genotypes[0].attributes["REFPAIR"], "0");
    assert_eq!(first.genotypes[0].attributes["AF"], "0.750");
    assert_eq!((first.info["REF"].as_str(), first.info["REFPAIR"].as_str()), ("2", "1"));

    // ref-pair support is left out for a single breakend
    let second = &task.variants[1];
    assert_eq!(second.genotypes[1].attributes["REF"], "1");
    assert_eq!(second.genotypes[1].attributes["REFPAIR"], "0");
    assert_eq!(second.genotypes[1].attributes["AF"], "0.667");
    assert_eq!(second.genotypes[0].attributes["AF"], "1.00");
    assert_eq!((second.info["REF"].as_str(), second.info["REFPAIR"].as_str()), ("1", "0"));
}

#[test]
fn unmapped_region_median_depth() {
    let path = std::env::temp_dir().join(format!("esvee_rs_depth_unmap_{}.tsv", std::process::id()));
    std::fs::write(&path, "Chromosome\tPosStart\tPosEnd\tMaxDepth\n1\t100\t200\t50\n1\t300\t400\t0\n").unwrap();
    let regions = UnmappedRegionDepth::load(path.to_str().unwrap()).unwrap();
    std::fs::remove_file(&path).ok();

    let genotype = |ref_depth: i32, frags: i32| {
        let mut g = Genotype::default();
        g.attributes.insert("REF".to_string(), ref_depth.to_string());
        g.attributes.insert("REFPAIR".to_string(), "0".to_string());
        g.attributes.insert(TOTAL_FRAGS.to_string(), frags.to_string());
        g
    };
    let variant = |position: i32, ref_depth: i32, frags: i32| {
        let mut v = create_variant_context("x", position);
        v.genotypes = vec![genotype(ref_depth, frags)];
        v.info.insert("REF".to_string(), ref_depth.to_string());
        v
    };
    // three ordinary variants give a median depth of 30; then one of high allele fraction in a region with
    // a depth, one in a region without, and one outside
    let mut variants = vec![variant(1000, 20, 1), variant(2000, 30, 1), variant(3000, 40, 1), variant(150, 1, 99), variant(350, 1, 99), variant(500, 1, 99)];
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    regions.set_unmapped_regions_depth(1, &mut refs).unwrap();

    assert_eq!(variants[3].genotypes[0].attributes["REF"], "30");
    assert_eq!(variants[3].genotypes[0].attributes["AF"], "0.767");
    assert_eq!(variants[3].info["REF"], "30");
    assert_eq!(variants[4].genotypes[0].attributes["REF"], "1");
    assert_eq!(variants[5].genotypes[0].attributes["REF"], "1");
    assert_eq!(variants[0].genotypes[0].attributes["REF"], "20");

    // more genotypes in the record than samples configured: Java fails with an index error
    let mut two = variant(150, 1, 99);
    two.genotypes.push(genotype(1, 99));
    let mut all = vec![variant(1000, 20, 1), variant(2000, 30, 1), two];
    let mut refs: Vec<&mut VariantContext> = all.iter_mut().collect();
    assert!(regions.set_unmapped_regions_depth(1, &mut refs).is_err());
}
