//! The Java unit tests of the package (DepthAnnotationTest), with the same cases and assertions.

use super::config::DepthConfig;
use crate::common::compat;
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

/// DepthAnnotationTest.testPositionProgression up to its last two assertions; returns the remembered
/// index and the number of uncapped variants at the end, with the number of variants.
fn position_progression() -> (usize, usize, usize) {
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

    (f.task.slice_region_state().min_position_index, f.task.slice_region_state().uncapped_variants.len(), variant_count)
}

#[test]
fn test_position_progression() {
    // D7: the Java test expects the remembered index to drop from 6 to 5 when the variant at that index
    // is taken off the list of uncapped variants, which only happens with the fault
    let (min_position_index, uncapped, variant_count) = compat::with_java(true, position_progression);
    assert_eq!(5, min_position_index);
    assert_eq!(variant_count - 1, uncapped);
}

// ---- tests of the corrected behaviour (not in the Java package) ----

#[test]
fn position_progression_keeps_index_when_variant_at_index_is_capped() {
    // D7: the capped variant was the one at the remembered index, so the index now points at the next
    let (min_position_index, uncapped, variant_count) = compat::with_java(false, position_progression);
    assert_eq!(6, min_position_index);
    assert_eq!(variant_count - 1, uncapped);
}

/// A pair either side of a variant on a sequence with index 200: (REF, REFPAIR) counted
fn pair_on_high_index_sequence() -> (i32, i32) {
    const CHR_201: &str = "201";
    let mut f = Fixture::new();
    f.task.reset();
    let id = f.next_variant_id();
    f.task.add_variants(vec![create_variant_context(&id, 1000)]).unwrap();
    f.task.add_slice_variants();

    let read = create_sam_record(&f.next_read_id(), CHR_201, 800, CHR_201, 1300, "100M", true, false, false, "");
    let mate = create_sam_record(&read.name, CHR_201, 1300, CHR_201, 800, "100M", false, true, false, "");
    f.task.process_read(&read).unwrap();
    f.task.process_read(&mate).unwrap();
    let counts = &f.task.variant_infos()[0].sample_support_counts[0];
    (counts.ref_support, counts.ref_pair_support)
}

#[test]
fn mate_on_same_sequence_with_high_index() {
    // D1
    assert_eq!(compat::with_java(false, pair_on_high_index_sequence), (0, 1));
    assert_eq!(compat::with_java(true, pair_on_high_index_sequence), (0, 0));
}

/// Two variants with a cap of 10: ten reads over the second only, then five over both. Returns the
/// reference support of each.
fn counts_after_cap() -> (i32, i32) {
    let mut f = Fixture::new();
    f.task.reset();
    let variants: Vec<VariantContext> = [1000, 1050].iter().map(|pos| { let id = f.next_variant_id(); create_variant_context(&id, *pos) }).collect();
    f.task.add_variants(variants).unwrap();
    f.task.add_slice_variants();
    for _ in 0..10 {
        let read = create_sam_record(&f.next_read_id(), CHR_1, 1020, CHR_2, 100, "100M", true, false, false, "");
        f.task.process_read(&read).unwrap();
    }
    for _ in 0..5 {
        let read = create_sam_record(&f.next_read_id(), CHR_1, 960, CHR_2, 100, "100M", true, false, false, "");
        f.task.process_read(&read).unwrap();
    }
    let infos = f.task.variant_infos();
    (infos[0].sample_support_counts[0].ref_support, infos[1].sample_support_counts[0].ref_support)
}

#[test]
fn cap_is_a_bound() {
    // D7: the second variant is at its cap but is behind the first, which each later read is relevant to,
    // so Java never takes it off the list and goes on counting
    assert_eq!(compat::with_java(false, counts_after_cap), (5, 10));
    assert_eq!(compat::with_java(true, counts_after_cap), (5, 15));
}

/// One variant with a cap of 10 and nine reads over it, and two fragments still waiting at the end of
/// the slice: a read over the variant whose mate never comes, and a pair either side of it held back for
/// a supplementary. `pair_first`: the pair's first read comes before the other fragment's.
fn waiting_fragments_at_the_cap(pair_first: bool) -> (i32, i32) {
    let config = port_config(&[TEST_SAMPLE_ID]);
    let mut task = DepthTask::new(CHR_1, &config, vec![0]);
    task.add_variants(vec![create_variant_context("x", 1000)]).unwrap();

    let (pair_start, single_start) = if pair_first { (690, 700) } else { (700, 690) };
    let mut reads: Vec<TestRead> = Vec::new();
    let [mut first, second] = pair("pair", pair_start, 1200);
    first.supp_data = "1,1500,+,50M,60,0;".to_string();
    reads.push(first);
    reads.push(second);
    reads.push(create_sam_record("single", CHR_1, single_start, CHR_1, 1300, "320M", true, false, false, ""));
    for i in 0..9 {
        reads.push(create_sam_record(&format!("r{i}"), CHR_1, 950, CHR_2, 100, "100M", true, false, false, ""));
    }
    reads.sort_by_key(|r| r.info.start);

    let mut sources: Vec<Box<dyn ReadSource>> = vec![Box::new(TestSource { reads })];
    let mut slicer = DepthTask::new("", &config, Vec::new());
    slicer.slice_group(CHR_1, task.infos_mut(), &mut sources).unwrap();
    let counts = &task.variant_infos()[0].sample_support_counts[0];
    (counts.ref_support, counts.ref_pair_support)
}

#[test]
fn waiting_fragments_are_counted_in_arrival_order() {
    // D7: with the cap as a bound only one of the two waiting fragments can still be counted; it is the
    // one whose first read came first, whatever the order of the map they are kept in
    for _ in 0..20 {
        assert_eq!(compat::with_java(false, || waiting_fragments_at_the_cap(true)), (9, 1));
        assert_eq!(compat::with_java(false, || waiting_fragments_at_the_cap(false)), (10, 0));
    }
    // Java counts both
    assert_eq!(compat::with_java(true, || waiting_fragments_at_the_cap(true)), (10, 1));
}

fn record(chromosome: &str, position: i32) -> VariantContext {
    let mut variant = create_variant_context("x", position);
    variant.chromosome = chromosome.to_string();
    variant
}

fn grouped(java: bool, config: &DepthConfig, records: &[(&str, i32)]) -> Vec<(String, Vec<i32>)> {
    let records: Vec<VariantContext> = records.iter().map(|(c, p)| record(c, *p)).collect();
    let (groups, _) = compat::with_java(java, || super::group_by_chromosome(config, records));
    groups.into_iter().map(|(chr, variants)| (chr, variants.iter().map(|v| v.position).collect())).collect()
}

fn chr_positions(list: &[(&str, &[i32])]) -> Vec<(String, Vec<i32>)> {
    list.iter().map(|(c, p)| (c.to_string(), p.to_vec())).collect()
}

#[test]
fn missing_filter_stays_missing() {
    // D2
    let config = port_config(&[TEST_SAMPLE_ID]);
    let mut filtered = record("1", 200);
    filtered.filters = Some(vec!["minQual".to_string()]);
    let mut passed = record("1", 300);
    passed.filters = Some(Vec::new());
    let records = vec![record("1", 100), filtered, passed];
    let filters = |java: bool| -> Vec<Option<Vec<String>>> {
        let (groups, _) = compat::with_java(java, || super::group_by_chromosome(&config, records.clone()));
        groups[0].1.iter().map(|v| v.filters.clone()).collect()
    };
    assert_eq!(filters(false), vec![None, Some(vec!["minQual".to_string()]), Some(Vec::new())]);
    assert_eq!(filters(true), vec![Some(Vec::new()), Some(vec!["minQual".to_string()]), Some(Vec::new())]);
}

#[test]
fn poly_g_regions_have_a_chromosome() {
    // D3: the test configuration is of genome version 37, where one region is on chromosome 2
    let config = port_config(&[TEST_SAMPLE_ID]);
    let records = [("1", 33141300), ("1", 50000000), ("2", 33141300), ("2", 50000000)];
    assert_eq!(grouped(false, &config, &records), chr_positions(&[("1", &[33141300, 50000000]), ("2", &[50000000])]));
    assert_eq!(grouped(true, &config, &records), chr_positions(&[("1", &[50000000]), ("2", &[50000000])]));
}

#[test]
fn other_contigs_are_kept() {
    // D4: after the chromosomes, in the order they first appear
    let config = port_config(&[TEST_SAMPLE_ID]);
    let records = [("extra_b", 5), ("2", 100), ("1", 100), ("extra_a", 7), ("extra_a", 9)];
    assert_eq!(grouped(false, &config, &records), chr_positions(&[("1", &[100]), ("2", &[100]), ("extra_b", &[5]), ("extra_a", &[7, 9])]));
    assert_eq!(grouped(true, &config, &records), chr_positions(&[("1", &[100]), ("2", &[100])]));
    // whether anything is left to annotate counts the other contigs in both
    let only_other = vec![record("extra_a", 7)];
    assert!(compat::with_java(true, || super::group_by_chromosome(&config, only_other.clone())).1);
    assert!(compat::with_java(false, || super::group_by_chromosome(&config, only_other)).1);
}

#[test]
fn chromosome_in_two_blocks_keeps_both() {
    // D5: the blocks are merged by position
    let config = port_config(&[TEST_SAMPLE_ID]);
    let records = [("1", 100), ("1", 300), ("2", 50), ("1", 200), ("1", 300), ("1", 400)];
    assert_eq!(grouped(false, &config, &records), chr_positions(&[("1", &[100, 200, 300, 300, 400]), ("2", &[50])]));
    assert_eq!(grouped(true, &config, &records), chr_positions(&[("1", &[200, 300, 400]), ("2", &[50])]));
    // of two records at one position the one of the first block comes first
    let mut first = record("1", 300);
    first.id = "first".to_string();
    let mut second = record("1", 300);
    second.id = "second".to_string();
    let (groups, _) = compat::with_java(false, || super::group_by_chromosome(&config, vec![first, record("2", 50), second]));
    assert_eq!(groups[0].1.iter().map(|v| v.id.as_str()).collect::<Vec<_>>(), vec!["first", "second"]);
}

#[test]
fn specific_chromosomes_restrict() {
    // D14
    let mut config = port_config(&[TEST_SAMPLE_ID]);
    config.specific_chromosomes = vec!["2".to_string()];
    assert_eq!(grouped(false, &config, &[("1", 100), ("2", 100)]), chr_positions(&[("2", &[100])]));
    // regions, when given, are what counts
    config.specific_regions = vec![("1".to_string(), 50, 150)];
    assert_eq!(grouped(false, &config, &[("1", 100), ("1", 200), ("2", 100)]), chr_positions(&[("1", &[100])]));
}

#[test]
fn all_filtered_warning_names_the_file() {
    // D15
    assert_eq!(compat::with_java(false, || super::all_filtered_message(12, "in.vcf.gz")), "all variants filtered from vcf(in.vcf.gz)");
    assert_eq!(compat::with_java(true, || super::all_filtered_message(12, "in.vcf.gz")), "all variants filtered from vcf(12)");
}

#[test]
fn prefixed_tags_get_header_lines() {
    // D11: an input without the REF and REFPAIR lines, as the assembly step writes it
    use super::vcf::{HeaderDefs, VcfHeader};
    let mut config = port_config(&[TEST_SAMPLE_ID]);
    config.vcf_tag_prefix = Some("PFX".to_string());
    let lines = |java: bool, with_plain: bool| -> Vec<String> {
        let mut header = VcfHeader { meta_lines: Vec::new(), samples: vec![TEST_SAMPLE_ID.to_string()], has_genotypes: true, defs: HeaderDefs::default() };
        if with_plain {
            header.add_format_line("REF", 1, "Integer", "x");
            header.add_format_line("REFPAIR", 1, "Integer", "x");
        }
        compat::with_java(java, || super::add_header_lines(&config, &mut header));
        assert!(header.defs.has_format("AF") && header.defs.has_info("AF"));
        let mut ids: Vec<String> = header.meta_lines.iter().filter(|l| l.contains("ID=PFX_")).map(|l| l[2..l.find(",Number").unwrap()].to_string()).collect();
        ids.sort();
        ids
    };
    let all = vec!["FORMAT=<ID=PFX_REF", "FORMAT=<ID=PFX_REFPAIR", "INFO=<ID=PFX_REF", "INFO=<ID=PFX_REFPAIR"];
    assert_eq!(lines(false, false), all);
    assert_eq!(lines(false, true), all);
    assert_eq!(lines(true, true), all);
    assert!(lines(true, false).is_empty());
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

fn unmap_regions() -> UnmappedRegionDepth {
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("esvee_rs_depth_unmap_{}_{}.tsv", std::process::id(), n));
    std::fs::write(&path, "Chromosome\tPosStart\tPosEnd\tMaxDepth\n1\t100\t200\t50\n1\t300\t400\t0\n").unwrap();
    let regions = UnmappedRegionDepth::load(path.to_str().unwrap()).unwrap();
    std::fs::remove_file(&path).ok();
    regions
}

fn depth_genotype(ref_tag: &str, ref_depth: i32, frags: i32) -> Genotype {
    let mut g = Genotype::default();
    g.attributes.insert(ref_tag.to_string(), ref_depth.to_string());
    g.attributes.insert("REFPAIR".to_string(), "0".to_string());
    g.attributes.insert(TOTAL_FRAGS.to_string(), frags.to_string());
    g
}

fn depth_variant(position: i32, genotypes: Vec<Genotype>) -> VariantContext {
    let mut v = create_variant_context("x", position);
    v.genotypes = genotypes;
    v
}

#[test]
fn unmapped_region_median_depth() {
    let regions = unmap_regions();
    let variant = |position: i32, ref_depth: i32, frags: i32| {
        let mut v = depth_variant(position, vec![depth_genotype("REF", ref_depth, frags)]);
        v.info.insert("REF".to_string(), ref_depth.to_string());
        v
    };
    // three ordinary variants give a median depth of 30; then one of high allele fraction in a region with
    // a depth, one in a region without, and one outside
    let mut variants = vec![variant(1000, 20, 1), variant(2000, 30, 1), variant(3000, 40, 1), variant(150, 1, 99), variant(350, 1, 99), variant(500, 1, 99)];
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    assert_eq!(regions.set_unmapped_regions_depth(&[0], "REF", "REFPAIR", &mut refs).unwrap(), 1);

    assert_eq!(variants[3].genotypes[0].attributes["REF"], "30");
    assert_eq!(variants[3].genotypes[0].attributes["AF"], "0.767");
    assert_eq!(variants[3].info["REF"], "30");
    assert_eq!(variants[4].genotypes[0].attributes["REF"], "1");
    assert_eq!(variants[5].genotypes[0].attributes["REF"], "1");
    assert_eq!(variants[0].genotypes[0].attributes["REF"], "20");
}

#[test]
fn unmapped_region_goes_over_the_configured_samples() {
    // D8: the record has two genotypes and the second sample (column 1) is the only one configured
    let regions = unmap_regions();
    let two = |position: i32, first: (i32, i32), second: (i32, i32)| {
        depth_variant(position, vec![depth_genotype("REF", first.0, first.1), depth_genotype("REF", second.0, second.1)])
    };
    let make = || vec![two(1000, (7, 1), (20, 1)), two(2000, (7, 1), (30, 1)), two(3000, (7, 1), (40, 1)), two(150, (2, 99), (1, 99))];

    let mut variants = make();
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    assert_eq!(compat::with_java(false, || regions.set_unmapped_regions_depth(&[1], "REF", "REFPAIR", &mut refs)).unwrap(), 1);
    // the configured sample is raised to its own median; the other genotype is not touched
    assert_eq!(variants[3].genotypes[1].attributes["REF"], "30");
    assert_eq!(variants[3].genotypes[0].attributes["REF"], "2");
    assert!(!variants[3].genotypes[0].attributes.contains_key("AF"));
    assert_eq!(variants[3].info["REF"], "30");

    // Java: index error
    let mut variants = make();
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    assert!(compat::with_java(true, || regions.set_unmapped_regions_depth(&[1], "REF", "REFPAIR", &mut refs)).is_err());
}

#[test]
fn unmapped_region_leaves_a_genotype_without_values() {
    // D9: the second sample's genotype has no values in the record to raise
    let regions = unmap_regions();
    let two = |position: i32, first: (i32, i32), second: (i32, i32)| {
        depth_variant(position, vec![depth_genotype("REF", first.0, first.1), depth_genotype("REF", second.0, second.1)])
    };
    let make = || {
        let mut variants = vec![two(1000, (20, 1), (50, 1)), two(2000, (30, 1), (60, 1)), two(3000, (40, 1), (70, 1)), two(150, (1, 99), (0, 0))];
        variants[3].genotypes[1] = Genotype::default();
        variants
    };

    let mut variants = make();
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    assert_eq!(compat::with_java(false, || regions.set_unmapped_regions_depth(&[0, 1], "REF", "REFPAIR", &mut refs)).unwrap(), 1);
    assert_eq!(variants[3].genotypes[0].attributes["REF"], "30");
    assert!(variants[3].genotypes[1].attributes.is_empty());
    assert_eq!(variants[3].info["REF"], "30");

    let mut variants = make();
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    assert!(compat::with_java(true, || regions.set_unmapped_regions_depth(&[0, 1], "REF", "REFPAIR", &mut refs)).is_err());
}

#[test]
fn unmapped_region_uses_the_prefixed_tags_and_lists_a_variant_once() {
    // D16: counts written as PFX_REF; REF holds values the input came with
    let regions = unmap_regions();
    let genotype = |counted: i32, input: i32, frags: i32| {
        let mut g = depth_genotype("PFX_REF", counted, frags);
        g.attributes.insert("PFX_REFPAIR".to_string(), "0".to_string());
        g.attributes.insert("REF".to_string(), input.to_string());
        g
    };
    let two = |position: i32, counted: i32, frags: i32| depth_variant(position, vec![genotype(counted, 500, frags), genotype(counted, 500, frags)]);
    let make = || vec![two(1000, 20, 1), two(2000, 30, 1), two(3000, 40, 1), two(150, 1, 99)];

    let mut variants = make();
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    // both samples qualify the last variant; it is listed once
    assert_eq!(compat::with_java(false, || regions.set_unmapped_regions_depth(&[0, 1], "PFX_REF", "PFX_REFPAIR", &mut refs)).unwrap(), 1);
    for g in &variants[3].genotypes {
        assert_eq!(g.attributes["PFX_REF"], "30");
        assert_eq!(g.attributes["REF"], "500");
        assert_eq!(g.attributes["AF"], "0.767");
    }
    assert_eq!(variants[3].info["PFX_REF"], "60");

    // Java reads REF: by those values no genotype has a high allele fraction, so nothing is raised
    let mut variants = make();
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    assert_eq!(compat::with_java(true, || regions.set_unmapped_regions_depth(&[0, 1], "PFX_REF", "PFX_REFPAIR", &mut refs)).unwrap(), 0);
    assert_eq!(variants[3].genotypes[0].attributes["PFX_REF"], "1");

    // Java with plain tags: the variant is listed once for each sample
    let plain = |position: i32, counted: i32, frags: i32| depth_variant(position, vec![depth_genotype("REF", counted, frags), depth_genotype("REF", counted, frags)]);
    let mut variants = vec![plain(1000, 20, 1), plain(2000, 30, 1), plain(3000, 40, 1), plain(150, 1, 99)];
    let mut refs: Vec<&mut VariantContext> = variants.iter_mut().collect();
    assert_eq!(compat::with_java(true, || regions.set_unmapped_regions_depth(&[0, 1], "REF", "REFPAIR", &mut refs)).unwrap(), 2);
    assert_eq!(variants[3].info["REF"], "60");
}
