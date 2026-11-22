//! The Java unit tests of the caller package (FiltersTest, HotspotTest, PonTest, RepeatMaskTest,
//! SeqTechTest with CallerTestUtils), case for case, followed by tests of the parts that only this
//! port has (VCF text handling, HashMap order).

use super::annotate::*;
use super::filters::*;
use super::pon::*;
use super::types::*;
use super::variant::*;
use super::vcf::*;
use super::{calculate_fragment_length_bounds, is_germline};

// ---------------------------------------------------------------------------------------------
// CallerTestUtils

const TEST_REF_ID: &str = "REF_ID";
const TEST_SAMPLE_ID: &str = "SAMPLE_ID";
const DEFAULT_QUAL: f64 = 60.0;
const CHR_1: &str = "1";
const CHR_2: &str = "2";
const ORIENT_FWD: i8 = 1;
const ORIENT_REV: i8 = -1;

/// attribute values as the tests give them: integers, decimals, integer pairs, text, or a flag
#[derive(Clone)]
enum Val {
    I(i32),
    D(f64),
    Ints(i32, i32),
    S(String),
    Flag,
}

#[derive(Clone, Default)]
struct Attrs(Vec<(String, Val)>);

impl Attrs {
    fn new() -> Attrs {
        Attrs::default()
    }
    fn put(&mut self, key: &str, value: Val) {
        match self.0.iter_mut().find(|(k, _)| k == key) {
            Some(entry) => entry.1 = value,
            None => self.0.push((key.to_string(), value)),
        }
    }
    fn clear(&mut self) {
        self.0.clear();
    }
    fn text(value: &Val) -> Option<String> {
        match value {
            Val::I(i) => Some(i.to_string()),
            // Double.toString of the values used in the tests
            Val::D(d) => Some(if d.fract() == 0.0 { format!("{d:.1}") } else { d.to_string() }),
            Val::Ints(a, b) => Some(format!("{a},{b}")),
            Val::S(s) => Some(s.clone()),
            Val::Flag => None,
        }
    }
}

/// makeCommonAttributes(vcfId, 50, 50, 50)
fn default_common_attributes() -> Attrs {
    let mut a = Attrs::new();
    a.put("QUAL", Val::D(50.0));
    a.put(SPLIT_FRAGS, Val::I(50));
    a.put("DF", Val::I(50));
    a.put(TOTAL_FRAGS, Val::I(100));
    a.put(AVG_FRAG_LENGTH, Val::I(500));
    a
}

fn make_genotype_attributes(split_frags: i32, disc_frags: i32, ref_frags: i32) -> Attrs {
    let mut a = Attrs::new();
    a.put(SPLIT_FRAGS, Val::I(split_frags));
    a.put("DF", Val::I(disc_frags));
    a.put(TOTAL_FRAGS, Val::I(split_frags + disc_frags));
    a.put(REF_DEPTH, Val::I(ref_frags));
    a
}

fn create_context(vcf_id: &str, chromosome: &str, position: i32, alt: &str, common: Option<&Attrs>, reference: Option<&Attrs>, tumor: Option<&Attrs>) -> VcfRecord {
    let mut record = VcfRecord::new(chromosome, position, vcf_id, "A", alt);
    record.log10_p_error = -(DEFAULT_QUAL / 10.0);
    let common = common.cloned().unwrap_or_else(default_common_attributes);
    for (key, value) in &common.0 {
        record.info.push((key.clone(), Attrs::text(value)));
    }
    let genotype = |sample: &str, attrs: Attrs| Genotype {
        sample: sample.to_string(),
        attrs: attrs.0.iter().map(|(k, v)| (k.clone(), Attrs::text(v).unwrap_or_default())).collect(),
        ..Default::default()
    };
    record.genotypes.push(genotype(TEST_REF_ID, reference.cloned().unwrap_or_else(|| make_genotype_attributes(0, 0, 30))));
    record.genotypes.push(genotype(TEST_SAMPLE_ID, tumor.cloned().unwrap_or_else(|| make_genotype_attributes(30, 20, 100))));
    record
}

/// VariantAltInsertCoords.formPairedAltString
fn form_paired_alt_string(alt: &str, insert_sequence: &str, chromosome: &str, position: i32, orient_start: i8, orient_end: i8) -> String {
    match (orient_start, orient_end) {
        (1, -1) => format!("{alt}{insert_sequence}[{chromosome}:{position}["),
        (1, 1) => format!("{alt}{insert_sequence}]{chromosome}:{position}]"),
        (-1, -1) => format!("[{chromosome}:{position}[{insert_sequence}{alt}"),
        _ => format!("]{chromosome}:{position}]{insert_sequence}{alt}"),
    }
}

// the test genotype ids: reference at 0 and a tumor ordinal beyond the genotypes, which reads as none
const REF_ORDINAL: i32 = 0;
const TUMOR_ORDINAL: i32 = 10;

#[allow(clippy::too_many_arguments)]
fn create_sv_with(
    vcf_id: &str, chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, orient_start: i8, orient_end: i8, ins_seq: &str,
    common: Option<&Attrs>, reference: Option<&Attrs>, tumor: Option<&Attrs>,
) -> Variant {
    let alt_start = form_paired_alt_string("A", ins_seq, chr_end, pos_end, orient_start, orient_end);
    let alt_end = form_paired_alt_string("A", ins_seq, chr_start, pos_start, orient_end, orient_start);
    let id_start = format!("{vcf_id}_start");
    let id_end = format!("{vcf_id}_end");
    let context_start = create_context(&id_start, chr_start, pos_start, &alt_start, common, reference, tumor);
    let context_end = create_context(&id_end, chr_end, pos_end, &alt_end, common, reference, tumor);
    Variant::new_sv(context_start, context_end, REF_ORDINAL, TUMOR_ORDINAL)
}

/// createSv without attribute maps: all three are empty
#[allow(clippy::too_many_arguments)]
fn create_sv(vcf_id: &str, chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, orient_start: i8, orient_end: i8, ins_seq: &str) -> Variant {
    let empty = Attrs::new();
    create_sv_with(vcf_id, chr_start, chr_end, pos_start, pos_end, orient_start, orient_end, ins_seq, Some(&empty), Some(&empty), Some(&empty))
}

fn create_sgl(vcf_id: &str, chromosome: &str, position: i32, orientation: i8, ins_seq: &str, common: Option<&Attrs>, reference: Option<&Attrs>, tumor: Option<&Attrs>) -> Variant {
    let alt = if orientation == 1 { format!("A{ins_seq}.") } else { format!(".{ins_seq}A") };
    let context = create_context(&format!("{vcf_id}_start"), chromosome, position, &alt, common, reference, tumor);
    Variant::new_sgl(context, REF_ORDINAL, TUMOR_ORDINAL)
}

fn set_genotype_attr(var: &mut Variant, se: usize, sample: &str, key: &str, value: i32) {
    var.breakend_mut(se).context.genotype_mut(sample).unwrap().set_attr(key, value.to_string());
}

fn frag_lengths() -> FragmentLengthBounds {
    FragmentLengthBounds { lower_bound: 100, upper_bound: 900, median: 500, std_deviation: 0.1 }
}

fn filter_constants() -> FilterConstants {
    FilterConstants::defaults(false, RefGenomeVersion::V37)
}

fn variant_filters(stats: DiscordantStats, sequencing_type: SequencingType) -> VariantFilters {
    VariantFilters::new(filter_constants(), frag_lengths(), stats, sequencing_type)
}

fn illumina_filters() -> VariantFilters {
    variant_filters(DiscordantStats::default(), SequencingType::Illumina)
}

fn sbx_filters() -> VariantFilters {
    variant_filters(DiscordantStats::default(), SequencingType::Sbx)
}

/// applyFilters on a single variant, returning it
fn apply(filters: &VariantFilters, var: Variant) -> Variant {
    let mut variants = vec![var];
    filters.apply_filters(&mut variants, 0);
    variants.pop().unwrap()
}

fn disc_stats(total_reads: i64, _prep_reads: i64, counts: &[(usize, i64)]) -> DiscordantStats {
    let mut stats = DiscordantStats { total_reads, ..Default::default() };
    for &(i, c) in counts {
        stats.type_counts[i] = c;
    }
    stats
}

// ---------------------------------------------------------------------------------------------
// FiltersTest

#[test]
fn test_deduplication() {
    let var1 = create_sv("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "");
    let var2 = create_sv("02", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "");
    let mut variants = vec![var1, var2];
    let map = BreakendMap::build(&variants);
    deduplicate_variants(&mut variants, &map);
    assert!(variants[0].is_pass());
    assert!(!variants[1].is_pass());
}

#[test]
fn test_min_af_filter() {
    let filters = illumina_filters();
    let mut tumor = Attrs::new();
    tumor.put(SPLIT_FRAGS, Val::I(1));
    tumor.put("DF", Val::I(0));
    tumor.put(TOTAL_FRAGS, Val::I(1));
    tumor.put(REF_DEPTH, Val::I(100));

    let mut reference = Attrs::new();
    reference.put(SPLIT_FRAGS, Val::I(0)); // not above the threshold in the ref sample
    reference.put(TOTAL_FRAGS, Val::I(0));
    reference.put(REF_DEPTH, Val::I(100));

    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "", None, Some(&reference), Some(&tumor)));
    assert!(!var.filters.contains(FilterType::MinAf));

    tumor.put(REF_DEPTH, Val::I(2000));
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "", None, Some(&reference), Some(&tumor)));
    assert!(var.filters.contains(FilterType::MinAf));

    // again for a SGL
    tumor.put(REF_DEPTH, Val::I(100));
    let var = apply(&filters, create_sgl("01", CHR_1, 100, ORIENT_FWD, "", None, Some(&reference), Some(&tumor)));
    assert!(var.filters.contains(FilterType::MinAf));

    let mut common = Attrs::new();
    common.put(LINE_SITE, Val::Flag);
    let var = apply(&filters, create_sgl("01", CHR_1, 100, ORIENT_FWD, "", Some(&common), Some(&reference), Some(&tumor)));
    assert!(!var.filters.contains(FilterType::MinAf));
}

#[test]
fn test_min_fragment_length() {
    let filters = illumina_filters();
    let mut common = Attrs::new();
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "", Some(&common), None, None));
    assert!(!var.filters.contains(FilterType::ShortFragLength));

    common.put(AVG_FRAG_LENGTH, Val::I(300));
    common.put(TOTAL_FRAGS, Val::I(50));
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "", Some(&common), None, None));
    assert!(var.filters.contains(FilterType::ShortFragLength));

    // not tested for LINE or SGLs
    common.put(LINE_SITE, Val::Flag);
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "", Some(&common), None, None));
    assert!(!var.filters.contains(FilterType::ShortFragLength));

    common.clear();
    common.put(AVG_FRAG_LENGTH, Val::I(300));
    common.put(TOTAL_FRAGS, Val::I(50));
    let var = apply(&filters, create_sgl("01", CHR_1, 100, ORIENT_FWD, "", Some(&common), None, None));
    assert!(!var.filters.contains(FilterType::ShortFragLength));
}

#[test]
fn test_inversion_short_low_vaf_homology() {
    let mut common = Attrs::new();
    common.put(HOMSEQ, Val::S("AGGCT".to_string()));
    common.put(IHOMPOS, Val::Ints(-10, 10));

    let mut var = create_sv_with("01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_FWD, "", Some(&common), None, None);
    set_genotype_attr(&mut var, SE_START, TEST_SAMPLE_ID, TOTAL_FRAGS, 5);
    set_genotype_attr(&mut var, SE_END, TEST_SAMPLE_ID, TOTAL_FRAGS, 5);

    let inv_filters = variant_filters(disc_stats(100_000, 10_000, &[(DISC_INV_1_TO_5K, 10)]), SequencingType::Illumina);
    let var = apply(&inv_filters, var);
    assert!(!var.filters.contains(FilterType::InvShortLowVafHom));

    let inv_filters = variant_filters(disc_stats(100_000, 10_000, &[(DISC_INV_1_TO_5K, 10_000)]), SequencingType::Illumina);
    let var = apply(&inv_filters, var);
    assert!(var.filters.contains(FilterType::InvShortLowVafHom));
}

#[test]
fn test_inversion_short_fragment_low_vaf() {
    let empty = Attrs::new();
    let mut var = create_sv_with("01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_FWD, "", Some(&empty), None, None);
    set_genotype_attr(&mut var, SE_START, TEST_SAMPLE_ID, TOTAL_FRAGS, 3);
    set_genotype_attr(&mut var, SE_END, TEST_SAMPLE_ID, TOTAL_FRAGS, 3);

    // set required rate to 2.5%
    let inv_filters = variant_filters(disc_stats(100_000, 10_000, &[(DISC_INV_LT_1K, 6_000)]), SequencingType::Illumina);
    let var = apply(&inv_filters, var);
    assert!(var.filters.contains(FilterType::InvShortFragLowVaf));
}

#[test]
fn test_short_low_vaf_deletion() {
    let mut common = Attrs::new();
    common.put(HOMSEQ, Val::S("AGGCT".to_string()));
    common.put(IHOMPOS, Val::Ints(-10, 10));
    common.put(AVG_FRAG_LENGTH, Val::I(300));

    let mut tumor = Attrs::new();
    tumor.put(SPLIT_FRAGS, Val::I(5));
    tumor.put("DF", Val::I(0));
    tumor.put(TOTAL_FRAGS, Val::I(5));
    tumor.put(REF_DEPTH, Val::I(101));

    let mut var = create_sv_with("01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_REV, "", Some(&common), None, Some(&tumor));
    set_genotype_attr(&mut var, SE_START, TEST_SAMPLE_ID, TOTAL_FRAGS, 5);
    set_genotype_attr(&mut var, SE_END, TEST_SAMPLE_ID, TOTAL_FRAGS, 5);

    let var = apply(&illumina_filters(), var);
    assert!(var.filters.contains(FilterType::DelShortLowVaf));
}

#[test]
fn test_isolated_inv_filter() {
    let empty = Attrs::new();
    let mut common = Attrs::new();
    common.put(IHOMPOS, Val::Ints(-10, 10));

    let inv1 = create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), None, Some(&empty));

    // 2 artefacts close to each other, but still marked as isolated
    let mut inv2 = create_sv_with("02", CHR_1, CHR_1, 1100, 1150, ORIENT_REV, ORIENT_REV, "", Some(&common), None, Some(&empty));
    inv2.add_filter(FilterType::Pon);

    let mut inv3 = create_sv_with("03", CHR_1, CHR_1, 1140, 1200, ORIENT_FWD, ORIENT_FWD, "", Some(&common), None, Some(&empty));
    inv3.add_filter(FilterType::InvShortFragLowVaf);

    // chained and higher-UFP INVs are excluded
    let mut common = Attrs::new();
    common.put(UNIQUE_FRAG_POSITIONS, Val::I(4 + 1));
    common.put(ASM_LINKS, Val::S("123".to_string()));
    let inv4 = create_sv_with("04", CHR_1, CHR_1, 2000, 2050, ORIENT_FWD, ORIENT_FWD, "", Some(&common), None, Some(&empty));

    let common = Attrs::new();

    // adjacent to another non-INV artefact breakend
    let inv5 = create_sv_with("05", CHR_1, CHR_1, 4000, 4050, ORIENT_REV, ORIENT_REV, "", Some(&common), None, Some(&empty));
    let var = create_sv_with("00", CHR_1, CHR_1, 4100, 10000, ORIENT_FWD, ORIENT_REV, "", Some(&common), None, Some(&empty));

    // 2 adjacent INVs are both marked as isolated
    let inv6 = create_sv_with("06", CHR_1, CHR_1, 5000, 5050, ORIENT_FWD, ORIENT_FWD, "", Some(&common), None, Some(&empty));
    let inv7 = create_sv_with("07", CHR_1, CHR_1, 5100, 5150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), None, Some(&empty));

    let mut variants = vec![inv1, inv2, inv3, inv4, inv5, var, inv6, inv7];
    let map = BreakendMap::build(&variants);
    illumina_filters().apply_adjacent_filters(&mut variants, &map);

    let isolated = |i: usize| variants[i].filters.contains(FilterType::InvShortIsolated);
    assert!(isolated(0));
    assert!(isolated(1));
    assert!(isolated(2));
    assert!(!isolated(3));
    assert!(!isolated(4));
    assert!(!isolated(5));
    assert!(isolated(6));
    assert!(isolated(7));
}

#[test]
fn test_prime_range_strand_bias() {
    let filters = sbx_filters();
    let mut common = Attrs::new();
    common.put(THREE_PRIME_RANGE, Val::Ints(3, 5));
    common.put(ASM_INFO, Val::S("1:100:1".to_string()));

    let mut tumor = Attrs::new();
    tumor.put(SPLIT_FRAGS, Val::I(4));
    tumor.put(STRAND_BIAS, Val::I(1));

    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor)));
    assert!(var.filters.contains(FilterType::UnpairedThreePrimeRange));

    common.put(ASM_INFO, Val::S("1:100:1_2:2000:-1".to_string())); // not a single assembly
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor)));
    assert!(var.filters.contains(FilterType::UnpairedThreePrimeRange));

    common.put(THREE_PRIME_RANGE, Val::Ints(20, 30));
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor)));
    assert!(!var.filters.contains(FilterType::UnpairedThreePrimeRange));
}

#[test]
fn test_weak_junction_min_af() {
    let filters = variant_filters(disc_stats(1_000, 1_000, &[(0, 100)]), SequencingType::Illumina);

    let mut common = Attrs::new();
    common.put(ASM_INFO, Val::S("1:100:1".to_string()));

    let mut tumor = Attrs::new();
    tumor.put(SPLIT_FRAGS, Val::I(4));
    tumor.put(STRAND_BIAS, Val::D(0.95));

    let mut var = create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor));
    var.is_hotspot = true;
    assert!(!apply_weak_junction_min_af_filter(&var));

    // not indel-based
    common.put(ASM_INFO, Val::S("1:100:1:I_1:100:-1:I".to_string()));
    let var = create_sv_with("01", CHR_1, CHR_1, 100, 101, ORIENT_FWD, ORIENT_REV, "", Some(&common), Some(&tumor), Some(&tumor));
    assert!(!apply_weak_junction_min_af_filter(&var));

    let soft_clip_junctions = "1:100:-1_1:200:1";
    common.put(ASM_INFO, Val::S(soft_clip_junctions.to_string()));
    let var = create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor));
    assert!(apply_weak_junction_min_af_filter(&var));

    tumor.clear();
    tumor.put(SPLIT_FRAGS, Val::I(4));
    tumor.put(STRAND_BIAS, Val::D(0.5));

    common.clear();
    common.put(STRAND_BIAS, Val::D(0.4));
    common.put(ASM_INFO, Val::S(soft_clip_junctions.to_string()));
    common.put(PROX_JUNC_READ_RATIO, Val::D(0.05));

    // test out each of the junction-related conditions
    let var = create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor));
    assert!(apply_weak_junction_min_af_filter(&var));

    common.put(PROX_JUNC_READ_RATIO, Val::D(0.15));
    common.put(SUPP_REMOTE_REGION_RATIO, Val::D(0.45));
    let var = create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor));
    assert!(apply_weak_junction_min_af_filter(&var));

    let var = apply(&filters, var);
    assert!(var.filters.contains(FilterType::WeakJuncMinAf));
}

#[test]
fn test_mark_germline() {
    let common = Attrs::new();
    let mut var = create_sv_with("01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_FWD, "", Some(&common), None, None);

    set_genotype_attr(&mut var, SE_START, TEST_SAMPLE_ID, TOTAL_FRAGS, 50);
    set_genotype_attr(&mut var, SE_START, TEST_SAMPLE_ID, REF_DEPTH, 50);
    set_genotype_attr(&mut var, SE_START, TEST_REF_ID, TOTAL_FRAGS, 2);
    set_genotype_attr(&mut var, SE_START, TEST_REF_ID, REF_DEPTH, 7);
    assert!(is_germline(&var, Some(TEST_REF_ID)));

    set_genotype_attr(&mut var, SE_START, TEST_SAMPLE_ID, TOTAL_FRAGS, 500);
    set_genotype_attr(&mut var, SE_START, TEST_SAMPLE_ID, REF_DEPTH, 500);
    set_genotype_attr(&mut var, SE_START, TEST_REF_ID, TOTAL_FRAGS, 2);
    set_genotype_attr(&mut var, SE_START, TEST_REF_ID, REF_DEPTH, 7);
    assert!(!is_germline(&var, Some(TEST_REF_ID)));
}

// ---------------------------------------------------------------------------------------------
// HotspotTest

#[test]
fn test_hotspot_classification() {
    let mut hotspot_cache = HotspotCache::default();
    hotspot_cache.add_hotspot(KnownHotspot {
        region_start: ChrBaseRegion::new(CHR_1, 1000, 10000),
        orient_start: Orient::Forward,
        region_end: ChrBaseRegion::new(CHR_2, 20000, 30000),
        orient_end: Orient::Reverse,
    });

    let vcf_id = "ID_01";
    let var = create_sv(vcf_id, CHR_1, CHR_2, 1050, 25000, ORIENT_FWD, ORIENT_REV, "");
    assert!(hotspot_cache.is_hotspot_variant(&var));

    // invalid orientation
    let var = create_sv(vcf_id, CHR_1, CHR_2, 1050, 25000, ORIENT_REV, ORIENT_REV, "");
    assert!(!hotspot_cache.is_hotspot_variant(&var));

    // SGLs
    let empty = Attrs::new();
    let mut common = Attrs::new();
    let var = create_sgl(vcf_id, CHR_1, 1050, ORIENT_FWD, "", Some(&common), Some(&empty), Some(&empty));
    assert!(!hotspot_cache.is_hotspot_variant(&var));

    let alt_alignment = AlternativeAlignment { chromosome: CHR_2.to_string(), position: 25000, orient: Orient::Reverse, cigar: "100M".to_string(), num_mutations: 0 };
    common.put(INSALN, Val::S(alt_alignment.vcf_string()));
    let var = create_sgl(vcf_id, CHR_1, 1050, ORIENT_FWD, "", Some(&common), Some(&empty), Some(&empty));
    assert!(hotspot_cache.is_hotspot_variant(&var));

    // reversed
    let alt_alignment = AlternativeAlignment { chromosome: CHR_1.to_string(), position: 1500, orient: Orient::Forward, cigar: "100M".to_string(), num_mutations: 0 };
    common.put(INSALN, Val::S(alt_alignment.vcf_string()));
    let var = create_sgl(vcf_id, CHR_2, 25000, ORIENT_REV, "", Some(&common), Some(&empty), Some(&empty));
    assert!(hotspot_cache.is_hotspot_variant(&var));
}

// ---------------------------------------------------------------------------------------------
// PonTest: the merging of PON regions, which Java keeps in the PON building utility (PonCombiner)

#[derive(Clone)]
struct SvRegion {
    start: (i32, i32),
    orient_start: Orient,
    end: (i32, i32),
    orient_end: Orient,
    pon_count: i32,
}

/// PonCombiner.mergeSvRegions, for regions that share both chromosomes
fn merge_sv_regions(regions: &mut Vec<SvRegion>) {
    // PonSvRegion.compareTo
    let orient_rank = |o: Orient| if o.is_forward() { 0 } else { 1 };
    regions.sort_by_key(|r| (r.start.0, r.start.1, orient_rank(r.orient_start), r.end.0, r.end.1, orient_rank(r.orient_end), r.pon_count));
    if regions.len() < 2 {
        return;
    }
    let mut index = 0;
    while index < regions.len() {
        loop {
            let mut found_merge = false;
            // search all regions starting with this base and merge if both ends overlap
            let mut next_index = index + 1;
            while next_index < regions.len() {
                let next = regions[next_index].clone();
                let current = &mut regions[index];
                // must overlap at the start
                if next.start.0 > current.start.1 {
                    break;
                }
                if next.orient_start != current.orient_start || next.orient_end != current.orient_end {
                    next_index += 1;
                    continue;
                }
                if !positions_overlap(current.end.0, current.end.1, next.end.0, next.end.1) {
                    next_index += 1;
                    continue;
                }
                current.start = (current.start.0.min(next.start.0), current.start.1.max(next.start.1));
                current.end = (current.end.0.min(next.end.0), current.end.1.max(next.end.1));
                regions.remove(next_index);
                found_merge = true;
            }
            if !found_merge {
                break;
            }
            // otherwise repeat now region(s) have been expanded
        }
        index += 1;
    }
}

/// PonCombiner.mergeSglRegions
fn merge_sgl_regions(regions: &mut Vec<PonSglRegion>) {
    if regions.len() < 2 {
        return;
    }
    // PonSglRegion.compareTo
    regions.sort_by_key(|r| (r.start, r.end, if r.orient.is_forward() { 0 } else { 1 }, r.pon_count));
    let mut index = 0;
    while index < regions.len() {
        let mut next_index = index + 1;
        while next_index < regions.len() {
            let next = regions[next_index].clone();
            let current = &mut regions[index];
            if next.start > current.end {
                break;
            }
            if next.orient != current.orient {
                next_index += 1;
                continue;
            }
            current.end = current.end.max(next.end);
            regions.remove(next_index);
        }
        index += 1;
    }
}

#[test]
fn test_pon_sv_region_merge() {
    let region = |start: (i32, i32), orient_start: Orient, end: (i32, i32)| SvRegion { start, orient_start, end, orient_end: Orient::Reverse, pon_count: 10 };
    let mut combined = vec![
        region((100, 105), Orient::Forward, (100, 110)),
        // will be combine with 2
        region((98, 103), Orient::Forward, (102, 104)),
        // not combined, diff orientation
        region((102, 107), Orient::Reverse, (100, 110)),
        // not combined, latter breakend doesn't overlap
        region((104, 109), Orient::Forward, (120, 125)),
        // latter start but still overlaps, will be combined
        region((105, 106), Orient::Forward, (105, 115)),
        // new region
        region((200, 205), Orient::Forward, (100, 105)),
        region((202, 207), Orient::Forward, (95, 100)),
    ];
    merge_sv_regions(&mut combined);
    assert_eq!(4, combined.len());

    // check combined regions
    assert_eq!((98, 106), combined[0].start);
    assert_eq!((100, 115), combined[0].end);
    // unchanged
    assert_eq!((102, 107), combined[1].start);
    assert_eq!((100, 110), combined[1].end);
    // unchanged
    assert_eq!((104, 109), combined[2].start);
    assert_eq!((120, 125), combined[2].end);
    // merged
    assert_eq!((200, 207), combined[3].start);
    assert_eq!((95, 105), combined[3].end);
}

#[test]
fn test_pon_sgl_region_merge() {
    let region = |start: i32, end: i32, orient: Orient| PonSglRegion { start, end, orient, pon_count: 10 };
    let mut combined = vec![
        region(100, 105, Orient::Forward),
        // will be combine with 2
        region(98, 103, Orient::Forward),
        // not combined, diff orientation
        region(100, 105, Orient::Reverse),
        // latter start but still overlaps, will be combined
        region(105, 106, Orient::Forward),
        // new region
        region(200, 205, Orient::Reverse),
        region(203, 210, Orient::Reverse),
    ];
    merge_sgl_regions(&mut combined);
    assert_eq!(3, combined.len());

    // check combined regions
    assert_eq!((98, 106), (combined[0].start, combined[0].end));
    // unchanged
    assert_eq!((100, 105), (combined[1].start, combined[1].end));
    // merged
    assert_eq!((200, 210), (combined[2].start, combined[2].end));
}

// ---------------------------------------------------------------------------------------------
// RepeatMaskTest

#[test]
fn test_alternative_alignments() {
    // first BWA and VCF tag serialisation
    let bwa_loc_tag = "16,+24008715,44S28M46S,0;X,+133232624,44S27M47S,0;12,+54042138,37S35M46S,2;4,-84437081,46S25M47S,0";
    let alignments = AlternativeAlignment::from_location_tag(bwa_loc_tag);
    assert_eq!(4, alignments.len());

    let vcf_tag = AlternativeAlignment::to_vcf_tag(&alignments);
    let parsed = AlternativeAlignment::from_vcf_tag(&vcf_tag);
    assert_eq!(4, parsed.len());
    for (alignment, parsed_alignment) in alignments.iter().zip(parsed.iter()) {
        assert_eq!(alignment, parsed_alignment);
    }
}

// ---------------------------------------------------------------------------------------------
// SeqTechTest

#[test]
fn test_sbx_strand_bias() {
    let filters = sbx_filters();
    let common = Attrs::new();
    let mut tumor = Attrs::new();
    tumor.put(SPLIT_FRAGS, Val::I(4));
    tumor.put(STRAND_BIAS, Val::I(1));

    // too few SF for an INV
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor)));
    assert!(!var.filters.contains(FilterType::SbxStrandBias));

    // then sufficient
    tumor.put(SPLIT_FRAGS, Val::I(12));
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor)));
    assert!(var.filters.contains(FilterType::SbxStrandBias));

    // any SF for a BND as long as strand bias is 0 or 1
    tumor.put(SPLIT_FRAGS, Val::I(1));
    tumor.put(STRAND_BIAS, Val::D(0.95));
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_2, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor)));
    assert!(!var.filters.contains(FilterType::SbxStrandBias));

    tumor.put(SPLIT_FRAGS, Val::I(1));
    tumor.put(STRAND_BIAS, Val::I(0));
    let var = apply(&filters, create_sv_with("01", CHR_1, CHR_2, 100, 150, ORIENT_FWD, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor)));
    assert!(var.filters.contains(FilterType::SbxStrandBias));
}

#[test]
fn test_sbx_artefacts() {
    let filters = sbx_filters();
    let artefact = |var: Variant| apply(&filters, var).filters.contains(FilterType::SbxArtefact);

    // first test indels
    // test 1: non-line and stranded
    let mut common = Attrs::new();
    common.put(ASM_LENGTH, Val::I(1200));

    let mut tumor = Attrs::new();
    tumor.put(SPLIT_FRAGS, Val::I(10));
    tumor.put(STRAND_BIAS, Val::I(1));

    let pos_start = 100;
    let pos_end = 200;

    assert!(artefact(create_sv_with("01", CHR_1, CHR_1, pos_start, pos_end, ORIENT_FWD, ORIENT_REV, "", Some(&common), Some(&tumor), Some(&tumor))));

    // test 2: not applied for longer indels
    assert!(!artefact(create_sv_with("01", CHR_1, CHR_1, pos_start, pos_end + 300, ORIENT_REV, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor))));

    // test 3: penalty for breakends away from original junctions
    tumor.put(STRAND_BIAS, Val::D(0.5));
    common.put(ASM_INFO, Val::S("2:100:1_2:2000:-1".to_string())); // note diff chromosome
    assert!(artefact(create_sv_with("01", CHR_1, CHR_1, pos_start, pos_end, ORIENT_REV, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor))));

    // very short DUPs
    common.put(ASM_INFO, Val::S("1:100:-1_1:130:1".to_string())); // matching
    assert!(artefact(create_sv_with("01", CHR_1, CHR_1, pos_start, pos_start + 30, ORIENT_REV, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor))));

    // INVs
    tumor.put(SPLIT_FRAGS, Val::I(12));

    // test 1: short INV with long insert sequence
    let ins_sequence = "G".repeat(21);
    assert!(artefact(create_sv_with("01", CHR_1, CHR_1, pos_start, pos_start + 60, ORIENT_FWD, ORIENT_FWD, &ins_sequence, Some(&common), Some(&tumor), Some(&tumor))));

    // test 2: zero-length INV
    assert!(artefact(create_sv_with("01", CHR_1, CHR_1, pos_start, pos_start, ORIENT_FWD, ORIENT_FWD, &ins_sequence, Some(&common), Some(&tumor), Some(&tumor))));

    // BNDs
    let ins_sequence = "A".repeat(20);
    common.put(ASM_INFO, Val::S("1:100:1_2:200:1".to_string())); // matching
    assert!(!artefact(create_sv_with("01", CHR_1, CHR_2, pos_start, pos_end, ORIENT_FWD, ORIENT_FWD, &ins_sequence, Some(&common), Some(&tumor), Some(&tumor))));

    let ins_sequence = "G".repeat(20);
    assert!(artefact(create_sv_with("01", CHR_1, CHR_2, pos_start, pos_end, ORIENT_FWD, ORIENT_FWD, &ins_sequence, Some(&common), Some(&tumor), Some(&tumor))));

    // SGLs
    common.put(ASM_INFO, Val::S("1:100:1".to_string())); // matching
    assert!(!artefact(create_sgl("01", CHR_1, pos_start, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor))));

    // non-local and shorter assembly
    common.put(ASM_INFO, Val::S("2:100:1".to_string())); // not matching
    common.put(ASM_LENGTH, Val::I(400));
    assert!(artefact(create_sgl("01", CHR_1, pos_start, ORIENT_FWD, "", Some(&common), Some(&tumor), Some(&tumor))));
}

#[test]
fn test_sbx_zero_length_inversions() {
    let filters = sbx_filters();
    let zero_motif = |ins_sequence: &str| {
        let common = Attrs::new();
        let mut tumor = Attrs::new();
        tumor.put(SPLIT_FRAGS, Val::I(10));
        tumor.put(STRAND_BIAS, Val::I(1));
        let var = create_sv_with("01", CHR_1, CHR_1, 100, 100, ORIENT_FWD, ORIENT_FWD, ins_sequence, Some(&common), Some(&tumor), Some(&tumor));
        apply(&filters, var).filters.contains(FilterType::SbxInvZeroMotif)
    };

    // test 1: non-line and stranded
    assert!(zero_motif(""));

    // test 2: insert sequence with known motif
    assert!(zero_motif(&format!("A{}", SBX_INV_INSERT_MOTIFS[0])));

    // reverse-complimented
    assert!(zero_motif(&format!("{}A", &reverse_complement(SBX_INV_INSERT_MOTIFS[0])[2..])));

    // handle mismatches: insert, SNV and delete
    let motif = SBX_INV_INSERT_MOTIFS[1];
    assert!(zero_motif(&format!("{}G{}A{}{}", &motif[0..10], &motif[10..15], &motif[16..20], &motif[21..motif.len() - 2])));

    // too many mismatches
    let motif = SBX_INV_INSERT_MOTIFS[2];
    assert!(!zero_motif(&format!("TTTTT{}GAAAAA", &motif[5..motif.len() - 6])));
}

/// ignored in Java as well: not every sequence of the list matches a motif
#[test]
#[ignore]
fn test_sbx_inversion_sequence_matches() {
    let test_sequences = [
        "ACCATCAATCCGGATGTATGCCGGATTGATGGT",
        "AGGAGTAACATCCATGTATGGGATGTTACTCCT",
        "AGGCGAGTATTCCATGTATGGGAATACTCGCCT",
        "ACCCGTAGTATTAATGTATGTAATACTACGGGT",
        "ACTATTGAAGGCTATGTATGAGCCTTCAATAGT",
        "AGGTATTGCCTGGATGTATGCCAGGCAATACCT",
        "AGAATCGCAGTTAATGTATGTAACTGCGATTCT",
        "CCATCAATCCGGATGTATGCCGGATTGATGG",
        "AGGATTATGGCAAATGTATGTTGCCATAATCCT",
        "CCCGTAGTATTAATGTATGTAATACTACGGG",
        "GGCGAGTATTCCATGTATGGGAATACTCGCC",
        "ACCTGCGGTTGAAATGTATGTTCAACCGCAGGT",
        "AGGATTATGGCAAAATGTATGTTGCCATAATCCT",
        "CTATTGAAGGCTATGTATGAGCCTTCAATAG",
        "GGAGTAACATCCATGTATGGGATGTTACTCC",
        "GGTATTGCCTGGATGTATGCCAGGCAATACC",
        "AGGCGAGTATTCCATACATGGAATACTCGCCT",
        "AGGTTATCAGCTTATGTATGAAGCTGATAACCT",
        "ACCTGCGGTTGAAAATGTATGTTCAACCGCAGGT",
        "GAATCGCAGTTAATGTATGTAACTGCGATTC",
        "GCGAGTATTCCATGTATGGGAATACTCGC",
        "GGATTATGGCAAATGTATGTTGCCATAATCC",
        "AGGGACATTACCAATGTATGTGGTAATGTCCCT",
        "CCTGCGGTTGAAATGTATGTTCAACCGCAGG",
        "GGATTATGGCAAAATGTATGTTGCCATAATCC",
        "AGGCGAGTATCCCATACATGGAATACTCGCCT",
        "GGCGAGTATTCCATACATGGAATACTCGCC",
        "ACCACTTCGCCAAATGTATGTTGGCGAAGTGGT",
        "AGGCTTACTACGGATGTATGCCGTAGTAAGCCT",
        "AGGATCTCCTATTATGTATGAATAGGAGATCCT",
        "AGGCGAGTATTCCATGTATGGAAATACTCGCCT",
        "AGGCGAGTATTCCATGTATGGGAAACTCGCCT",
        "AGGCGAGTATTCCATGTAGGGAATACTCGCCT",
        "AGGAGTAACATCCATACATGGATGTTACTCCT",
        "ACCCGTAGTATTAAATGTATGTAATACTACGGGT",
        "GTATTGCCTGGATGTATGCCAGGCAATAC",
        "GAGTAACATCCATGTATGGGATGTTACTC",
        "AGGAGTAACATCCCATACATGGGATGTTACTCCT",
        "AGGCGAGTATTCCCATACAGGGAATACTCGCCT",
        "GGTTATCAGCTTATGTATGAAGCTGATAACC",
        "CATCAATCCGGATGTATGCCGGATTGATG",
        "ACTATTGAAGGCTATGTATGAACCTTCAATAGT",
        "TCTTCTGTATTTCTAGTGTTCAATAGAAATACAGAAGA",
        "AATCGCAGTTAATGTATGTAACTGCGATT",
        "ACCATCAATCCGGAATGTATGCCGGATTGATGGT",
        "ACCATCAATCCGCATACATCCGGATTGATGGT",
        "CCTGCGGTTGAAAATGTATGTTCAACCGCAGG",
        "AGAATCGCAGTTAAATGTATGTAACTGCGATTCT",
        "AGGATTATGGCAAATGTATGTGCCATAATCCT",
        "GATTATGGCAAATGTATGTTGCCATAATC",
        "ACCCGTAGTATTAACATACATTAATACTACGGGT",
        "CCACTTCGCCAAATGTATGTTGGCGAAGTGG",
        "ACAACAGCCGAAGATGTATGCTTCGGCTGTTGT",
        "AGGAGTAACATCCCATACACGGATGTTACTCCT",
        "ACATAGAAGGTAGATGTATGCTACCTTCTATGT",
        "AGGAGTAACATCCATGTAATGGGATGTTACTCCT",
        "GGAGTAACATCCATGTATGAGATGTTACTCC",
    ];
    for test_sequence in test_sequences {
        let mut has_match = false;
        let reversed = reverse_complement(test_sequence);
        for motif in SBX_INV_INSERT_MOTIFS {
            let length_diff = test_sequence.len() as i32 - motif.len() as i32;
            for i in 0..=length_diff.abs() {
                let s1_offset = if length_diff > 0 { i } else { 0 };
                let s2_offset = if length_diff < 0 { i } else { 0 };
                let mut diffs = sequence_mismatches(test_sequence, motif, s1_offset, s2_offset);
                diffs = diffs.min(sequence_mismatches(&reversed, motif, s1_offset, s2_offset));
                if diffs <= 10 {
                    has_match = true;
                    break;
                }
            }
            if has_match {
                break;
            }
        }
        assert!(has_match, "{test_sequence}");
    }
}

// ---------------------------------------------------------------------------------------------
// tests of this port only

#[test]
fn hash_map_order_of_chromosome_keys() {
    // the iteration order of a java.util.HashMap<String,?> holding these keys (16 buckets):
    // the hash codes of "1".."9" are 49..57, of "10" 1567 and of "X" 88
    let keys: Vec<String> = ["3", "1", "X", "10", "2"].iter().map(|s| s.to_string()).collect();
    let order: Vec<&str> = java_hash_map_order(&keys).iter().map(|&i| keys[i].as_str()).collect();
    assert_eq!(order, ["1", "2", "3", "X", "10"]);

    // "chr1" hashes to 3052802 and "chr2" to 3052803; folded with the high bits, buckets 12 and 13
    let keys: Vec<String> = ["chr2", "chr1"].iter().map(|s| s.to_string()).collect();
    let order: Vec<&str> = java_hash_map_order(&keys).iter().map(|&i| keys[i].as_str()).collect();
    assert_eq!(order, ["chr1", "chr2"]);
}

#[test]
fn breakend_alt_parsing() {
    assert_eq!(match_breakend("A[2:200["), Some(("A", '[', "2:200", "")));
    assert_eq!(match_breakend("]2:200]GTA"), Some(("", ']', "2:200", "GTA")));
    assert_eq!(match_breakend("ACG]HLA-A*01:01:5]"), Some(("ACG", ']', "HLA-A*01:01:5", "")));
    assert_eq!(match_breakend("A."), None);

    // a DUP of a single base: the reverse breakend becomes the start
    let var = create_sv("01", CHR_1, CHR_1, 100, 100, ORIENT_FWD, ORIENT_REV, "");
    assert_eq!(var.sv_type, SvType::Dup);
    assert_eq!(var.start().orient, Orient::Reverse);
    assert_eq!(var.start().vcf_id(), "01_end");

    let var = create_sv("01", CHR_1, CHR_1, 100, 101, ORIENT_FWD, ORIENT_REV, "ACGT");
    assert_eq!(var.sv_type, SvType::Ins);
    assert_eq!(var.insert_sequence, "ACGT");
    assert_eq!(var.start().insert_sequence, "ACGT");
    assert_eq!(var.adjusted_length(), 4);

    let empty = Attrs::new();
    let sgl = create_sgl("01", CHR_1, 100, ORIENT_REV, "TTG", Some(&empty), Some(&empty), Some(&empty));
    assert_eq!(sgl.insert_sequence, "TTG");
    assert_eq!(sgl.start().orient, Orient::Reverse);
}

#[test]
fn vcf_record_round_trip() {
    let header = VcfHeader::from_lines(
        vec![
            "##fileformat=VCFv4.1".to_string(),
            info_line("AF", "1", "Float", "Allele frequency of the breakend"),
            info_line("CIPOS", "2", "Integer", "Confidence interval around position"),
            info_line("LINE", "0", "Flag", "LINE insertion site"),
            info_line("VF", "1", "Integer", "Total variant fragments supporting the breakend"),
            format_line("GT", "1", "String", "Genotype"),
            format_line("AF", "1", "Float", "Allele frequency of the breakend"),
            format_line("VF", "1", "Integer", "Total variant fragments supporting the breakend"),
            format_line("REF", "1", "Integer", "Count of reads mapping across this breakend"),
            format_line("REFPAIR", "1", "Integer", "Count of reference read pairs"),
            "##contig=<ID=2,length=2000>".to_string(),
            "##contig=<ID=1,length=1000>".to_string(),
        ],
        vec!["N".to_string(), "T".to_string()],
    );
    // the standard INFO key AF is replaced by htsjdk's own line; FORMAT AF is not a standard key
    assert!(header.meta_lines.contains(&info_line("AF", "A", "Float", "Allele Frequency, for each ALT allele, in the same order as listed")));
    assert!(header.meta_lines.contains(&format_line("AF", "1", "Float", "Allele frequency of the breakend")));

    let line = "1\t100\tid1\ta\tA[2:200[\t53\tPASS\tVF=7;LINE;CIPOS=-1,2\tGT:VF:REF:REFPAIR:AF\t./.:3:10:.:0.25\t./.:4:.:.";
    let record = VcfRecord::parse(line, &header.samples).unwrap();
    assert_eq!(record.ref_allele, "A");
    assert_eq!(record.phred_scaled_qual(), 53.0);
    assert!(record.filters.is_empty());
    assert!(record.has_attr("LINE"));
    assert_eq!(record.attr_int("VF", 0), 7);
    assert_eq!(record.attr_int_list("CIPOS", 0), vec![-1, 2]);
    assert_eq!(record.attr_string("CIPOS", ""), "[-1, 2]");
    assert_eq!(record.genotypes[0].attr_int("REF", 0), 10);
    assert!(!record.genotypes[0].has_attr("REFPAIR"));
    assert!(!record.genotypes[1].has_attr("REF"));

    // written back: INFO sorted, the flag without a value, trailing missing genotype values dropped
    let counts = header.counts();
    let encoded = encode_record(&header, &counts, &record, record.log10_p_error, &record.genotypes).unwrap();
    assert_eq!(encoded.prefix, "1\t100\tid1\tA\tA[2:200[\t53");
    assert_eq!(encoded.suffix, "CIPOS=-1,2;LINE;VF=7\tGT:AF:REF:VF\t./.:0.25:10:3\t./.:.:.:4");

    // a sample of the header without a genotype is written as a no-call
    let encoded = encode_record(&header, &counts, &record, record.log10_p_error, &record.genotypes[1..]).unwrap();
    assert_eq!(encoded.suffix, "CIPOS=-1,2;LINE;VF=7\tGT:VF\t./.\t./.:4");

    // the header: file format first and always 4.2, the rest sorted, contigs in their own order
    let text = header_text(&header);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert!(lines[1].starts_with("##FORMAT=<ID=AF"));
    let contig = lines.iter().position(|l| l.starts_with("##contig")).unwrap();
    assert_eq!(lines[contig], "##contig=<ID=2,length=2000>");
    assert_eq!(lines[contig + 1], "##contig=<ID=1,length=1000>");
    assert_eq!(*lines.last().unwrap(), "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT");
}

#[test]
fn qual_formatting() {
    assert_eq!(format_qual(53.0), "53");
    assert_eq!(format_qual(53.5), "53.50");
    assert_eq!(format_qual(1234.567), "1234.57");
    assert_eq!(format_qual(0.0), "0");
}

#[test]
fn fragment_length_bounds() {
    // 1000 fragments: the bounds are the lengths at 0.25% and 99.75%, the deviation the distance
    // from the median to the length at 16%
    let frequencies = [(100, 10), (200, 190), (300, 400), (400, 390), (500, 10)];
    let bounds = calculate_fragment_length_bounds(&frequencies);
    assert_eq!(bounds.lower_bound, 100);
    assert_eq!(bounds.upper_bound, 500);
    assert_eq!(bounds.median, 300);
    assert_eq!(bounds.std_deviation, 100.0);
    assert_eq!(calculate_fragment_length_bounds(&[]).upper_bound, 0);
}

#[test]
fn repeats_and_trimmed_length() {
    // four or more copies of a unit of one to five bases
    let repeats = find_repeats(b"ACGTTTTTTGACACACACAC");
    assert_eq!(repeats, vec![RepeatInfo { index: 3, unit_length: 1, count: 6 }, RepeatInfo { index: 10, unit_length: 2, count: 5 }]);
    // each repeat counts as two units: 3 + 2 + 1 + 4
    assert_eq!(calc_trimmed_base_length(0, 19, &repeats), 10);
    assert!(find_repeats(b"ACGTACG").is_empty());
    assert_eq!(calc_trimmed_base_length(0, 6, &[]), 7);
}

#[test]
fn pon_matching() {
    let mut pon = PonCache::new(DEFAULT_PON_DISTANCE, DEFAULT_SGL_PON_DISTANCE, None, None).unwrap();
    pon.add_sv_region(CHR_1, (95, 105), Orient::Forward, CHR_2, (195, 205), Orient::Reverse, 7);
    pon.add_sv_region(CHR_1, (500, 510), Orient::Forward, CHR_2, (900, 910), Orient::Reverse, 3);
    pon.add_sgl_region(CHR_1, (300, 310), Orient::Forward, 5);

    let empty = Attrs::new();
    let mut variants = vec![
        create_sv("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, ""),
        // the other orientation at the far end
        create_sv("02", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_FWD, ""),
        // within the margin of the SGL entry, then outside it
        create_sgl("03", CHR_1, 320, ORIENT_FWD, "", Some(&empty), Some(&empty), Some(&empty)),
        create_sgl("04", CHR_1, 330, ORIENT_FWD, "", Some(&empty), Some(&empty), Some(&empty)),
        // beyond the margin of the second SV entry
        create_sv("05", CHR_1, CHR_2, 520, 905, ORIENT_FWD, ORIENT_REV, ""),
    ];
    let map = BreakendMap::build(&variants);
    pon.annotate_variants(&mut variants, &map);
    let counts: Vec<i32> = variants.iter().map(|v| v.pon_count).collect();
    assert_eq!(counts, [7, 0, 5, 0, 0]);
    assert!(variants[0].filters.contains(FilterType::Pon));
    assert!(variants[1].is_pass());
}

#[test]
fn repeat_mask_annotation() {
    let mut repeat_mask = RepeatMask::default();
    for i in 0..60 {
        repeat_mask.add_entry(CHR_1, i * 1000 + 1, i * 1000 + 400, "AluY", "SINE/Alu");
    }
    // the alignment covers 100 bases of an entry, of an insert sequence of 200
    let insert_sequence = "ACGT".repeat(50);
    let annotation = repeat_mask.annotate(&insert_sequence, "1:30301|+|100S100M|0").unwrap();
    assert_eq!(annotation.repeat, "AluY");
    assert_eq!(annotation.class_type, "SINE/Alu");
    assert_eq!(annotation.coverage, 0.5);
    // no entry there
    assert!(repeat_mask.annotate(&insert_sequence, "1:30501|+|100S100M|0").is_none());
    // a poly-A insert is reported without looking at the entries
    let poly_a = "A".repeat(40);
    let annotation = repeat_mask.annotate(&poly_a, "2:500|-|40M|0").unwrap();
    assert_eq!((annotation.repeat.as_str(), annotation.class_type.as_str()), ("A(n)", "Simple_repeat"));
    assert_eq!(annotation.coverage, 1.0);
}

#[test]
fn target_regions_and_hotspot_lookup() {
    let mut regions = TargetRegions::default();
    assert!(!regions.has_target_regions());
    regions.add_region(CHR_1, 1000, 2000);
    assert!(regions.in_target_regions(CHR_1, 2500, 1000));
    assert!(!regions.in_target_regions(CHR_1, 3500, 1000));
    assert!(!regions.in_target_regions(CHR_2, 1500, 1000));
}

// ---------------------------------------------------------------------------------------------
// end to end: the step on synthetic inputs (testdata/), against the files the Java step wrote for
// the same inputs (testdata/expected/)

struct Run {
    dir: std::path::PathBuf,
}

impl Run {
    /// Lays out the inputs under the names the step expects and runs it. `sample` names the prep
    /// files and the outputs: the tumor, or the reference without one.
    fn new(name: &str, sample: &str, sample_args: &[&str], extra_args: &[&str]) -> Run {
        let dir = std::env::temp_dir().join(format!("esvee_caller_test_{}_{}", name, std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("out")).unwrap();
        let write = |file: String, text: &str| {
            let path = dir.join(file);
            std::fs::write(&path, text).unwrap();
            path.to_str().unwrap().to_string()
        };
        let vcf = write(format!("{sample}.esvee.ref_depth.vcf"), include_str!("testdata/input.vcf"));
        write(format!("{sample}.esvee.prep.fragment_length.tsv"), include_str!("testdata/fragment_length.tsv"));
        write(format!("{sample}.esvee.prep.disc_stats.tsv"), include_str!("testdata/disc_stats.tsv"));
        write(format!("{sample}.esvee.breakend.tsv"), include_str!("testdata/breakend.tsv"));
        let sv_pon = write("sv_pon.bedpe".to_string(), include_str!("testdata/sv_pon.bedpe"));
        let sgl_pon = write("sgl_pon.bed".to_string(), include_str!("testdata/sgl_pon.bed"));
        let hotspots = write("hotspots.bedpe".to_string(), include_str!("testdata/hotspots.bedpe"));
        let repeat_mask = write("repeat_mask.fa.out".to_string(), include_str!("testdata/repeat_mask.fa.out"));
        write("target.bed".to_string(), "chr1\t900\t1100\nchr2\t29000\t31000\nchr3\t100\t200\n");

        let dir_str = dir.to_str().unwrap().to_string();
        let out_dir = dir.join("out").to_str().unwrap().to_string();
        let mut args: Vec<String> = sample_args.iter().map(|s| s.to_string()).collect();
        for (key, value) in [
            ("-input_vcf", vcf.as_str()),
            ("-ref_genome_version", "38"),
            ("-esvee_prep_dir", dir_str.as_str()),
            ("-pon_sv_file", sv_pon.as_str()),
            ("-pon_sgl_file", sgl_pon.as_str()),
            ("-known_hotspot_file", hotspots.as_str()),
            ("-repeat_mask_file", repeat_mask.as_str()),
            ("-output_dir", out_dir.as_str()),
        ] {
            args.push(key.to_string());
            args.push(value.to_string());
        }
        args.push("-write_breakend_tsv".to_string());
        for arg in extra_args {
            args.push(arg.replace("{dir}", &dir_str));
        }
        super::run(&args).unwrap();
        Run { dir }
    }

    fn output(&self, file: &str) -> std::path::PathBuf {
        self.dir.join("out").join(file)
    }

    /// the decompressed VCF, which must also have its index
    fn vcf(&self, file: &str) -> String {
        use std::io::Read;
        let path = self.output(file);
        assert!(self.output(&format!("{file}.tbi")).exists(), "no index for {file}");
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap()).read_to_string(&mut text).unwrap();
        text
    }

    /// the breakend table with the filters of each row in sorted order, as Java writes them in the
    /// order of a hash set
    fn tsv(&self, file: &str) -> String {
        sort_tsv_filters(&std::fs::read_to_string(self.output(file)).unwrap())
    }
}

impl Drop for Run {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn sort_tsv_filters(text: &str) -> String {
    let mut out = String::new();
    let mut column = None;
    for (i, line) in text.lines().enumerate() {
        let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        if i == 0 {
            column = fields.iter().position(|f| f == "Filter");
        } else if let Some(c) = column {
            let mut filters: Vec<&str> = fields[c].split(';').collect();
            filters.sort();
            fields[c] = filters.join(";");
        }
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    out
}

fn assert_same_text(name: &str, actual: &str, expected: &str) {
    if actual != expected {
        for (i, (a, e)) in actual.lines().zip(expected.lines()).enumerate() {
            assert_eq!(a, e, "{name}: line {}", i + 1);
        }
        assert_eq!(actual.lines().count(), expected.lines().count(), "{name}: number of lines");
        panic!("{name}: differs in line ends");
    }
}

#[test]
fn end_to_end_tumor_and_reference() {
    let run = Run::new("paired", "SAMPLE_T", &["-tumor", "SAMPLE_T", "-reference", "SAMPLE_N"], &[]);
    assert_same_text("unfiltered", &run.vcf("SAMPLE_T.esvee.unfiltered.vcf.gz"), include_str!("testdata/expected/paired.unfiltered.vcf"));
    assert_same_text("somatic", &run.vcf("SAMPLE_T.esvee.somatic.vcf.gz"), include_str!("testdata/expected/paired.somatic.vcf"));
    assert_same_text("germline", &run.vcf("SAMPLE_T.esvee.germline.vcf.gz"), include_str!("testdata/expected/paired.germline.vcf"));
    assert_same_text(
        "breakend table",
        &run.tsv("SAMPLE_T.esvee.final_breakend.tsv"),
        &sort_tsv_filters(include_str!("testdata/expected/paired.final_breakend.tsv")),
    );
}

#[test]
fn end_to_end_tumor_only() {
    let run = Run::new("tumor", "SAMPLE_T", &["-tumor", "SAMPLE_T"], &[]);
    assert_same_text("unfiltered", &run.vcf("SAMPLE_T.esvee.unfiltered.vcf.gz"), include_str!("testdata/expected/tumor.unfiltered.vcf"));
    assert_same_text("somatic", &run.vcf("SAMPLE_T.esvee.somatic.vcf.gz"), include_str!("testdata/expected/tumor.somatic.vcf"));
    assert!(!run.output("SAMPLE_T.esvee.germline.vcf.gz").exists());
    assert_same_text(
        "breakend table",
        &run.tsv("SAMPLE_T.esvee.final_breakend.tsv"),
        &sort_tsv_filters(include_str!("testdata/expected/tumor.final_breakend.tsv")),
    );
}

#[test]
fn end_to_end_germline_only() {
    let run = Run::new("germline", "SAMPLE_N", &["-reference", "SAMPLE_N"], &[]);
    assert_same_text("unfiltered", &run.vcf("SAMPLE_N.esvee.unfiltered.vcf.gz"), include_str!("testdata/expected/germline.unfiltered.vcf"));
    assert_same_text("germline", &run.vcf("SAMPLE_N.esvee.germline.vcf.gz"), include_str!("testdata/expected/germline.germline.vcf"));
    assert!(!run.output("SAMPLE_N.esvee.somatic.vcf.gz").exists());
    assert_same_text(
        "breakend table",
        &run.tsv("SAMPLE_N.esvee.final_breakend.tsv"),
        &sort_tsv_filters(include_str!("testdata/expected/germline.final_breakend.tsv")),
    );
}

#[test]
fn end_to_end_sbx() {
    let run = Run::new("sbx", "SAMPLE_T", &["-tumor", "SAMPLE_T", "-reference", "SAMPLE_N"], &["-sequencing_type", "SBX"]);
    assert_same_text("unfiltered", &run.vcf("SAMPLE_T.esvee.unfiltered.vcf.gz"), include_str!("testdata/expected/sbx.unfiltered.vcf"));
    assert_same_text("somatic", &run.vcf("SAMPLE_T.esvee.somatic.vcf.gz"), include_str!("testdata/expected/sbx.somatic.vcf"));
    assert_same_text("germline", &run.vcf("SAMPLE_T.esvee.germline.vcf.gz"), include_str!("testdata/expected/sbx.germline.vcf"));
}

#[test]
fn end_to_end_target_regions() {
    let run = Run::new(
        "panel",
        "SAMPLE_T",
        &["-tumor", "SAMPLE_T", "-reference", "SAMPLE_N"],
        &["-target_regions_bed", "{dir}/target.bed", "-filter_sgls", "-hotspot_min_support", "5"],
    );
    assert_same_text("unfiltered", &run.vcf("SAMPLE_T.esvee.unfiltered.vcf.gz"), include_str!("testdata/expected/panel.unfiltered.vcf"));
    assert_same_text("somatic", &run.vcf("SAMPLE_T.esvee.somatic.vcf.gz"), include_str!("testdata/expected/panel.somatic.vcf"));
    assert_same_text("germline", &run.vcf("SAMPLE_T.esvee.germline.vcf.gz"), include_str!("testdata/expected/panel.germline.vcf"));
}

// ---------------------------------------------------------------------------------------------
// reading of the reference files

fn temp_file(name: &str, contents: &[u8]) -> String {
    let path = std::env::temp_dir().join(format!("esvee_caller_test_{}_{}", std::process::id(), name));
    std::fs::write(&path, contents).unwrap();
    path.to_str().unwrap().to_string()
}

fn gzip(text: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(text).unwrap();
    encoder.finish().unwrap()
}

fn read_all(path: &str) -> std::io::Result<Vec<u8>> {
    use std::io::Read;
    let mut out = Vec::new();
    open_text(path)?.read_to_end(&mut out)?;
    Ok(out)
}

#[test]
fn gzip_reader() {
    // a text larger than the reader's buffers, in one member and in several
    let text: Vec<u8> = (0..60_000u32).flat_map(|i| format!("line {i} {}\n", i.wrapping_mul(2654435761)).into_bytes()).collect();
    let path = temp_file("one.gz", &gzip(&text));
    assert_eq!(read_all(&path).unwrap(), text);

    let mut members = Vec::new();
    for chunk in text.chunks(70_000) {
        members.extend(gzip(chunk));
    }
    members.extend(gzip(b""));
    let path2 = temp_file("many.gz", &members);
    assert_eq!(read_all(&path2).unwrap(), text);

    // bytes after the last member that are not another member are ignored
    members.extend(b"trailing");
    let path3 = temp_file("trailing.gz", &members);
    assert_eq!(read_all(&path3).unwrap(), text);

    // a truncated file is an error
    let whole = gzip(&text);
    let path4 = temp_file("truncated.gz", &whole[..whole.len() / 2]);
    assert!(read_all(&path4).is_err());

    // without the extension the bytes are read as they are
    let path5 = temp_file("plain.txt", b"abc\n");
    assert_eq!(read_all(&path5).unwrap(), b"abc\n");
    for p in [path, path2, path3, path4, path5] {
        let _ = std::fs::remove_file(p);
    }
}

#[test]
fn pon_files() {
    // BED formats have zero-based starts and +/- orientations
    let sv_bed = temp_file("sv.bedpe.gz", &gzip(b"1\t94\t105\t2\t194\t205\t.\t7\t+\t-\n1\t499\t510\t2\t899\t910\t.\t3\t+\t-\n"));
    let sgl_bed = temp_file("sgl.bed", b"1\t299\t310\t.\t5\t+\n");
    // the TSV formats have a header, one-based starts and orientations as 1 and -1
    let sv_tsv = temp_file("sv.tsv.gz", &gzip(b"ChrLower\tPosLowerStart\tPosLowerEnd\tChrUpper\tPosUpperStart\tPosUpperEnd\tOrientLower\tOrientUpper\tPonCount\n1\t95\t105\t2\t195\t205\t1\t-1\t7\n"));
    let sgl_tsv = temp_file("sgl.tsv.gz", &gzip(b"Chromosome\tPosStart\tPosEnd\tOrientation\tPonCount\n1\t300\t310\t1\t5\n"));

    let empty = Attrs::new();
    for (sv_file, sgl_file) in [(&sv_bed, &sgl_bed), (&sv_tsv, &sgl_tsv)] {
        let mut pon = PonCache::new(DEFAULT_PON_DISTANCE, DEFAULT_SGL_PON_DISTANCE, Some(sv_file), Some(sgl_file)).unwrap();
        assert!(pon.has_valid_data());
        let mut variants = vec![
            create_sv("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, ""),
            create_sgl("02", CHR_1, 320, ORIENT_FWD, "", Some(&empty), Some(&empty), Some(&empty)),
            create_sgl("03", CHR_1, 320, ORIENT_REV, "", Some(&empty), Some(&empty), Some(&empty)),
        ];
        let map = BreakendMap::build(&variants);
        pon.annotate_variants(&mut variants, &map);
        assert_eq!(variants.iter().map(|v| v.pon_count).collect::<Vec<_>>(), [7, 5, 0]);
    }

    // entries of a chromosome must be in order of their start
    let unordered = temp_file("unordered.bed", b"1\t299\t310\t.\t5\t+\n1\t99\t110\t.\t5\t+\n");
    let pon = PonCache::new(DEFAULT_PON_DISTANCE, DEFAULT_SGL_PON_DISTANCE, None, Some(&unordered)).unwrap();
    assert!(!pon.has_valid_data());
    // a file that cannot be read makes the data invalid; a line that cannot be parsed is an error
    let pon = PonCache::new(DEFAULT_PON_DISTANCE, DEFAULT_SGL_PON_DISTANCE, Some("/nonexistent/sv_pon.bedpe"), None).unwrap();
    assert!(!pon.has_valid_data());
    let bad = temp_file("bad.bed", b"1\tabc\t310\t.\t5\t+\n");
    assert!(PonCache::new(DEFAULT_PON_DISTANCE, DEFAULT_SGL_PON_DISTANCE, None, Some(&bad)).is_err());
    for p in [sv_bed, sgl_bed, sv_tsv, sgl_tsv, unordered, bad] {
        let _ = std::fs::remove_file(p);
    }
}

fn repeat_mask_line(chromosome: &str, begin: i32, end: i32, repeat: &str, class_type: &str, id: i32) -> String {
    format!("  463   1.3  0.6  1.7  {chromosome}  {begin:>9} {end:>9} (1000) +  {repeat}  {class_type}  1  463    (0)  {id}\n")
}

#[test]
fn repeat_mask_file() {
    let header = "   SW  perc\nscore  div.\n\n";
    let mut text = header.to_string();
    // more lines than one block of the reader holds, over two chromosomes, with another contig between
    for i in 0..50_000 {
        text.push_str(&repeat_mask_line("chr1", i * 100, i * 100 + 50, if i % 2 == 0 { "AluY" } else { "L1PA3" }, if i % 2 == 0 { "SINE/Alu" } else { "LINE/L1" }, i));
    }
    text.push_str(&repeat_mask_line("chr1_random", 0, 500, "AluY", "SINE/Alu", 1));
    for i in 0..10 {
        text.push_str(&repeat_mask_line("chr1", 8_000_000 + i * 100, 8_000_000 + i * 100 + 50, "MIR", "SINE/MIR", i));
    }
    for i in 0..5 {
        text.push_str(&repeat_mask_line("chr2", i * 100, i * 100 + 50, "MIR", "SINE/MIR", i));
    }
    assert!(text.len() > 2 * (1 << 21));
    let path = temp_file("rm.fa.out.gz", &gzip(text.as_bytes()));
    let mut repeat_mask = RepeatMask::default();
    repeat_mask.load(&path, RefGenomeVersion::V38).unwrap();
    assert_eq!(repeat_mask.entry_count(), 50_015);
    // starts are one-based: the entry of line 20_001 covers 2_000_101 to 2_000_150
    let insert_sequence = "ACGT".repeat(10);
    let annotation = repeat_mask.annotate(&insert_sequence, "chr1:2000111|+|40M|0").unwrap();
    assert_eq!((annotation.repeat.as_str(), annotation.class_type.as_str()), ("L1PA3", "LINE/L1"));
    assert_eq!(annotation.coverage, 1.0);
    let annotation = repeat_mask.annotate(&insert_sequence, "chr1:8000201|+|40M|0").unwrap();
    assert_eq!(annotation.repeat, "MIR");
    // the names follow the genome version
    let mut repeat_mask = RepeatMask::default();
    repeat_mask.load(&path, RefGenomeVersion::V37).unwrap();
    assert!(repeat_mask.annotate(&insert_sequence, "chr1:2000111|+|40M|0").is_none());
    assert!(repeat_mask.annotate(&insert_sequence, "1:2000111|+|40M|0").is_some());

    // a chromosome that comes back after another one starts again from nothing, as in Java
    let mut text = header.to_string();
    text.push_str(&repeat_mask_line("chr1", 0, 50, "AluY", "SINE/Alu", 1));
    text.push_str(&repeat_mask_line("chr2", 0, 50, "AluY", "SINE/Alu", 2));
    text.push_str(&repeat_mask_line("chr1", 100, 150, "AluY", "SINE/Alu", 3));
    let path2 = temp_file("rm_back.fa.out", text.as_bytes());
    let mut repeat_mask = RepeatMask::default();
    repeat_mask.load(&path2, RefGenomeVersion::V38).unwrap();
    assert_eq!(repeat_mask.entry_count(), 2);

    // an entry that cannot be parsed is an error that names the line and the entries read before it
    let mut text = header.to_string();
    text.push_str(&repeat_mask_line("chr1", 0, 50, "AluY", "SINE/Alu", 1));
    text.push_str("  463   1.3  0.6  1.7  chr1  x100  150 (1000) +  AluY  SINE/Alu  1  463  (0)  2\n");
    let path3 = temp_file("rm_bad.fa.out", text.as_bytes());
    let error = RepeatMask::default().load(&path3, RefGenomeVersion::V38).unwrap_err().to_string();
    assert!(error.starts_with("invalid RM file entry: index(  463   1.3") && error.ends_with("line(1)"), "{error}");
    for p in [path, path2, path3] {
        let _ = std::fs::remove_file(p);
    }
}

#[test]
fn hotspot_and_target_files() {
    let hotspots = temp_file("hotspots.bedpe", b"1\t999\t10000\t2\t19999\t30000\tGENE_PAIR\t0\t+\t-\n");
    let mut cache = HotspotCache::default();
    assert_eq!(cache.load(&hotspots).unwrap(), 1);
    assert!(cache.is_hotspot_variant(&create_sv("01", CHR_1, CHR_2, 1000, 30000, ORIENT_FWD, ORIENT_REV, "")));
    assert!(!cache.is_hotspot_variant(&create_sv("01", CHR_1, CHR_2, 999, 30000, ORIENT_FWD, ORIENT_REV, "")));
    // found from either end
    assert!(cache.is_hotspot_variant(&create_sv("01", CHR_2, CHR_1, 30000, 1000, ORIENT_REV, ORIENT_FWD, "")));

    // other contigs are left out, a header line is skipped, and a chromosome that comes back replaces
    // its earlier regions
    let bed = temp_file("target.bed.gz", &gzip(b"Chromosome\tPosStart\tPosEnd\nchr1\t99\t200\nchrUn\t0\t50\nchr2\t99\t200\nchr1\t999\t2000\n"));
    let mut regions = TargetRegions::default();
    assert_eq!(regions.load(&bed, RefGenomeVersion::V38).unwrap(), 2);
    assert!(regions.in_target_regions("chr1", 1000, 0));
    assert!(!regions.in_target_regions("chr1", 150, 0));
    assert!(regions.in_target_regions("chr2", 100, 0));
    assert!(!regions.in_target_regions("chr2", 99, 0));
    assert!(regions.in_target_regions("chr2", 99, 1));
    for p in [hotspots, bed] {
        let _ = std::fs::remove_file(p);
    }
}

#[test]
fn command_line() {
    use super::config::*;
    let args = |list: &[&str]| list.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    let values = ConfigValues::parse(&args(&["-tumor", "T1", "-ref_genome_version", "38", "-filter_sgls", "-pon_margin", "-2", "-output_dir", "out"])).ok().unwrap();
    assert!(values.has_flag("filter_sgls"));
    assert!(!values.has_flag("write_breakend_tsv"));
    assert_eq!(values.integer("pon_margin"), -2);
    assert_eq!(values.integer("min_qual"), 30);
    assert!(!values.has_value("min_qual"));
    let config = CallerConfig::new(&values).unwrap();
    assert_eq!(config.output_dir, "out/");
    assert_eq!(config.prep_dir, "out/");
    assert_eq!(config.vcf_file, "out/T1.esvee.ref_depth.vcf.gz");
    assert!(config.has_tumor() && !config.has_reference() && !config.germline_only());

    // the panel default of the quality threshold applies unless one is given
    let bed = temp_file("config_target.bed", b"chr1\t0\t10\n");
    let panel = ConfigValues::parse(&args(&["-ref_genome_version", "38", "-target_regions_bed", &bed])).ok().unwrap();
    assert_eq!(filter_constants(&panel, RefGenomeVersion::V38).min_qual, 60);
    let panel = ConfigValues::parse(&args(&["-ref_genome_version", "38", "-target_regions_bed", &bed, "-min_qual", "30"])).ok().unwrap();
    assert_eq!(filter_constants(&panel, RefGenomeVersion::V38).min_qual, 30);
    let _ = std::fs::remove_file(bed);

    // what Java rejects: unknown options, a missing genome version, a value that looks like an
    // option, a path that does not exist, a number that is not one
    assert!(ConfigValues::parse(&args(&["-ref_genome_version", "38", "-nonsense", "1"])).is_err());
    assert!(ConfigValues::parse(&args(&["-tumor", "T1"])).is_err());
    assert!(ConfigValues::parse(&args(&["-ref_genome_version", "38", "-tumor", "-reference"])).is_err());
    assert!(ConfigValues::parse(&args(&["-ref_genome_version", "38", "-pon_sv_file", "/nonexistent/pon.bedpe"])).is_err());
    assert!(ConfigValues::parse(&args(&["-ref_genome_version", "38", "-min_qual", "abc"])).is_err());
    assert!(ConfigValues::parse(&args(&["-ref_genome_version", "38", "stray"])).is_err());
}
