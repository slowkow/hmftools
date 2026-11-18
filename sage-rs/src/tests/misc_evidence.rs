//! Ported from the Java SAGE test suite: evidence/MiscEvidenceTest.java

use super::util::*;
use crate::candidate::{build_candidates, Candidate, CandidateFinder};
use crate::constants::*;
use crate::evidence::{max_coverage, ReadContextCounter, ReadMatchType};
use crate::evidence_stage::{passes_hard_filters, Collector, PhaseSink};
use crate::filters::{set_near_by_indel_status_post_filter, SageVariant, SoftFilter};
use crate::matcher::{is_simple_alt_match, MatchType, Matcher};
use crate::phase::Phaser;
use crate::quality::{set_jitter_qual_filter_state, BqrMap, MsiJitter, QualityCalculator};
use crate::read_context::{create_context, VariantReadContext};
use crate::refseq::RefSequence;
use crate::sync::{FragmentSync, SyncOutput};
use crate::types::*;
use std::sync::Arc;

/// Java String.substring(start, end)
fn sub(s: &str, start: usize, end: usize) -> &str {
    &s[start..end]
}

/// VariantUtils.createSimpleVariant(position, ref, alt)
fn create_simple_variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

const TEST_LEFT_FLANK: &str = "ACCGCTGACT";
const TEST_RIGHT_FLANK: &str = "CTGAGACTCA";
const TEST_LEFT_CORE: &str = "AC";
const TEST_RIGHT_CORE: &str = "GT";

/// VariantUtils.createReadContext(variant): a context with no repeat or homology
fn create_read_context_var(variant: &SimpleVariant) -> VariantReadContext {
    let (left_flank, right_flank, left_core, right_core) = (TEST_LEFT_FLANK, TEST_RIGHT_FLANK, TEST_LEFT_CORE, TEST_RIGHT_CORE);
    let core_index_start = left_flank.len() as i32;
    let var_read_index = core_index_start + left_core.len() as i32;
    let core_index_end = var_read_index + variant.alt_len() - 1 + right_core.len() as i32;
    let ref_bases = format!("{}{}{}", left_core, variant.ref_str(), right_core);
    let read_bases = format!("{}{}{}{}{}", left_flank, left_core, variant.alt_str(), right_core, right_flank);

    let alignment_start = variant.pos - var_read_index;

    let core_right_pos_start = variant.pos + variant.ref_len().min(variant.alt_len());
    let alignment_end = core_right_pos_start + right_core.len() as i32 + right_flank.len() as i32 - 1;

    let read_cigar = format!("{}M", read_bases.len());

    let core_position_start = variant.pos - left_core.len() as i32;
    let core_position_end = variant.pos + right_core.len() as i32;

    VariantReadContext {
        variant: variant.clone(),
        alignment_start,
        alignment_end,
        ref_bases: ref_bases.into_bytes(),
        read_bases: read_bases.into_bytes(),
        core_index_start,
        var_index: var_read_index,
        core_index_end,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start,
        core_position_end,
        read_cigar,
        // ReadCigarInfo.hasIndelInCore of a single M element
        has_indel_in_core: false,
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    }
}

/// VariantUtils.createReadContext(position, ref, alt)
pub(super) fn create_read_context(position: i32, ref_: &str, alt: &str) -> VariantReadContext {
    create_read_context_var(&create_simple_variant(position, ref_, alt))
}

/// VariantUtils.createReadCounter(0, readContext)
pub(super) fn create_read_counter(ctx: &VariantReadContext, fixture: &Fixture) -> ReadContextCounter {
    ReadContextCounter::new(Arc::new(ctx.clone()), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &fixture.qc(), false)
}

/// VariantUtils.createSageVariant(position, ref, alt)
fn create_sage_variant(position: i32, ref_: &str, alt: &str, fixture: &Fixture) -> SageVariant {
    let ctx = Arc::new(create_read_context(position, ref_, alt));
    let read_counter = ReadContextCounter::new(ctx.clone(), Tier::LowConfidence, 100, 1, &fixture.qc(), false);
    let candidate = Candidate { tier: Tier::HighConfidence, ctx, min_events: 1, full_match_support: 1, core_match_support: 0 };
    SageVariant::new(candidate, Vec::new(), vec![read_counter])
}

/// TestUtils.buildCigarString(alignedLength)
fn build_cigar_string(aligned_length: usize) -> String {
    format!("{aligned_length}M")
}

/// MockRefGenome.getNextBase
fn get_next_base(base: u8) -> u8 {
    const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    match DNA_BASES.iter().position(|b| *b == base) {
        Some(i) => DNA_BASES[(i + 1) % DNA_BASES.len()],
        None => base,
    }
}

/// s.substring(0, i) + getNextBase(s.charAt(i)) + s.substring(i + 1)
fn with_next_base_at(s: &str, i: usize) -> String {
    let mut bytes = s.as_bytes().to_vec();
    bytes[i] = get_next_base(bytes[i]);
    String::from_utf8(bytes).unwrap()
}

/// MockRefGenome.generateRandomBases: not random, cycles through the nucleotides
fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect()
}

fn read_bases_str(ctx: &VariantReadContext) -> &str {
    std::str::from_utf8(&ctx.read_bases).unwrap()
}

// Java: testLongInsertPartialReadMatches
#[test]
fn test_long_insert_partial_read_matches() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let position = 35;

    let inserted_bases = "TTTTTTTTTTT";

    let ref_base = sub(REF_BASES_200, 35, 36);
    let variant = create_simple_variant(position, ref_base, &format!("{ref_base}{inserted_bases}"));

    let read_bases = format!("{}{}{}", sub(REF_BASES_200, 1, 35), variant.alt_str(), sub(REF_BASES_200, 37, 60));
    let read_cigar = "35M11I24M";

    let initial_read = create_sam_record(&next_read_id(), 1, &read_bases, read_cigar);

    let read_context = create_context(&variant, &initial_read, 34, &fixture.ref_seq, DEFAULT_FLANK_LENGTH).expect("read context");

    assert_eq!(11, read_context.var_index);
    assert!(read_context.is_valid());
    assert_eq!(sub(&read_bases, 33, 49), read_context.core_str());

    let mut read_counter =
        ReadContextCounter::new(Arc::new(read_context), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &qc, true);

    read_counter.process_read(&initial_read, 1, None, &qc);

    assert_eq!(1, read_counter.counts.full);

    // now a read which matches the read context by 2 or more bases than then ref without being a full core match
    let alt_read_bases =
        format!("{}{}{}{}", sub(REF_BASES_200, 1, 35), sub(variant.alt_str(), 0, 8), "GGGG", sub(REF_BASES_200, 37, 60));
    let read_cigar = build_cigar_string(alt_read_bases.len());

    let alt_read = create_sam_record(&next_read_id(), 1, &alt_read_bases, &read_cigar);

    read_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(1, read_counter.counts.partial_core);
}

// Java: testReadEdgeDistancePenalty
#[test]
fn test_read_edge_distance_penalty() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let position = 100;

    let read_context = create_read_context(position, "A", "T");
    let ctx_read_bases = read_bases_str(&read_context);
    let var_index = read_context.var_index;

    let mut read_context_counter = create_read_counter(&read_context, &fixture);

    let alt_read_bases = format!("{}{}{}", sub(REF_BASES_200, 0, 10), ctx_read_bases, sub(REF_BASES_200, 0, 10));
    let read_cigar = build_cigar_string(alt_read_bases.len());
    let read_var_index = 10 + var_index;
    let read_pos_start = position - read_var_index;

    let alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    read_context_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(1, read_context_counter.counts.full);
    assert_eq!(25, read_context_counter.qualities.full);

    // now a read right up against the position from the left
    let alt_read_bases = format!("{}{}", &ctx_read_bases[var_index as usize..], sub(REF_BASES_200, 0, 30));
    let read_cigar = build_cigar_string(alt_read_bases.len());
    let read_pos_start = position;

    let alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    read_context_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(1, read_context_counter.counts.partial_core);
    assert_eq!(10, read_context_counter.qualities.partial_core);

    // now 1 base in from the edge
    let alt_read_bases = format!("{}{}", &ctx_read_bases[var_index as usize - 1..], sub(REF_BASES_200, 0, 30));
    let read_cigar = build_cigar_string(alt_read_bases.len());
    let read_pos_start = position - 1;

    let alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    read_context_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(2, read_context_counter.counts.partial_core);
    assert_eq!(30, read_context_counter.qualities.partial_core);

    // and from the other side
    let alt_read_bases = format!("{}{}", sub(REF_BASES_200, 0, 30), sub(ctx_read_bases, 0, var_index as usize + 1));
    let read_cigar = build_cigar_string(alt_read_bases.len());
    let read_pos_start = position - 30 - var_index;

    let alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    read_context_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(3, read_context_counter.counts.partial_core);
    assert_eq!(40, read_context_counter.qualities.partial_core);
}

// Java: testLowQualSoftClip
#[test]
fn test_low_qual_soft_clip() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let position = 100;

    let read_context = create_read_context(position, "A", "T");

    let mut read_context_counter = create_read_counter(&read_context, &fixture);

    let alt_read_bases = format!("{}{}", read_bases_str(&read_context), sub(REF_BASES_200, 0, 20));
    let read_cigar = "15S20M";
    let read_var_index = read_context.var_index;
    let read_pos_start = position - read_var_index + 15; // for the soft-clip

    let mut alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, read_cigar);

    read_context_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(1, read_context_counter.counts.full);

    // again but with too many low-qual soft-clip bases
    for i in 0..10 {
        alt_read.quals[i] = 11;
    }

    let read_match_type = read_context_counter.process_read(&alt_read, 1, None, &qc);
    assert_eq!(read_match_type, ReadMatchType::SoftClip);
    assert_eq!(1, read_context_counter.counts.full);
}

// Java: testMnvBaseQuality
#[test]
fn test_mnv_base_quality() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let position = 100;

    let read_context = create_read_context(position, "ATG", "TGA");

    let mut read_context_counter = create_read_counter(&read_context, &fixture);

    let alt_read_bases = format!("{}{}{}", sub(REF_BASES_200, 0, 10), read_bases_str(&read_context), sub(REF_BASES_200, 0, 10));
    let read_cigar = build_cigar_string(alt_read_bases.len());
    let read_var_index = 10 + read_context.var_index;
    let read_pos_start = position - read_var_index;

    let mut alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    read_context_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(37, read_context_counter.qual_counters.alt_recalibrated_base_quality_total);

    // min rather than average is used
    alt_read.quals[read_var_index as usize] = 11;

    read_context_counter.process_read(&alt_read, 1, None, &qc);

    assert_eq!(48, read_context_counter.qual_counters.alt_recalibrated_base_quality_total);
}

// Java: testLowQualMatchTypes
#[test]
fn test_low_qual_match_types() {
    let position = 100;

    let read_context = create_read_context(position, "A", "G");
    let ctx_read_bases = read_bases_str(&read_context);

    // test reads with low-qual mismatches
    let ref_bases = with_next_base_at(std::str::from_utf8(&read_context.ref_bases).unwrap(), 1);

    let extra_read_bases = sub(REF_BASES_200, 0, 10);
    let read_ref_bases = format!(
        "{}{}{}{}{}",
        extra_read_bases,
        read_context.left_flank_str(),
        ref_bases,
        read_context.right_flank_str(),
        extra_read_bases
    );
    let read_cigar = build_cigar_string(read_ref_bases.len());
    let read_var_index = extra_read_bases.len() as i32 + read_context.var_index;
    let read_pos_start = position - read_var_index;

    let mut ref_read = create_sam_record(&next_read_id(), read_pos_start, &read_ref_bases, &read_cigar);
    let low_qual_index = extra_read_bases.len() + read_context.left_flank_length() as usize + 1;

    let read_context_matcher = Matcher::new(&read_context, true, false);

    ref_read.quals[low_qual_index] = 10;
    let match_info = read_context_matcher.determine_read_match(&read_context, &ref_read, read_var_index);

    assert_eq!(MatchType::Ref, match_info.kind);
    assert!(!match_info.exact);

    let alt_bases = with_next_base_at(ctx_read_bases, 11);
    let alt_read_bases = format!("{extra_read_bases}{alt_bases}{extra_read_bases}");

    let mut alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    alt_read.quals[low_qual_index] = 10;
    let match_info = read_context_matcher.determine_read_match(&read_context, &alt_read, read_var_index);

    assert_eq!(MatchType::Full, match_info.kind);
    assert!(!match_info.exact);

    // Java also asserts the boxed Boolean is not null
    let simple_match = is_simple_alt_match(&read_context.variant, &alt_read, read_var_index);
    assert!(simple_match);

    // mismatch in each flank
    let alt_bases = with_next_base_at(ctx_read_bases, 5);
    let alt_read_bases = format!("{extra_read_bases}{alt_bases}{extra_read_bases}");

    let mut alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    let low_qual_index = extra_read_bases.len() + 5;
    alt_read.quals[low_qual_index] = 10;
    let match_info = read_context_matcher.determine_read_match(&read_context, &alt_read, read_var_index);

    assert_eq!(MatchType::Full, match_info.kind);
    assert!(!match_info.exact);

    let alt_bases = with_next_base_at(ctx_read_bases, 16);
    let alt_read_bases = format!("{extra_read_bases}{alt_bases}{extra_read_bases}");

    let mut alt_read = create_sam_record(&next_read_id(), read_pos_start, &alt_read_bases, &read_cigar);

    let low_qual_index = extra_read_bases.len() + 16;
    alt_read.quals[low_qual_index] = 10;
    let match_info = read_context_matcher.determine_read_match(&read_context, &alt_read, read_var_index);

    assert_eq!(MatchType::Full, match_info.kind);
    assert!(!match_info.exact);
}

// Java: testReadContextCounterOrdering
#[test]
fn test_read_context_counter_ordering() {
    let region = BaseRegion::new(1, 300);

    let ref_bases = format!("{}{}", REF_BASES_200, generate_random_bases(1500)); // since expects region to be 1300+

    // test a collection of variants, not phased

    // reads which establish the variants
    let read_bases1 = format!("{}A{}", sub(REF_BASES_200, 30, 50), sub(REF_BASES_200, 51, 70));
    let read_cigar = "40M";
    let read1 = build_sam_record_quals(30, read_cigar, &read_bases1, default_quals(read_bases1.len()));
    let read1_clone = build_sam_record_quals(30, read_cigar, &read_bases1, default_quals(read_bases1.len()));

    let read_bases2 = format!("{}{}", sub(REF_BASES_200, 30, 48), sub(REF_BASES_200, 52, 70)); // 4-base delete
    let read_cigar = "18M4D17M";
    let read2 = build_sam_record_quals(30, read_cigar, &read_bases2, default_quals(read_bases2.len()));
    let read2_clone = build_sam_record_quals(30, read_cigar, &read_bases2, default_quals(read_bases2.len()));

    let read_bases3 = format!("{}G{}", sub(REF_BASES_200, 40, 55), sub(REF_BASES_200, 56, 80));
    let read_cigar = "40M";
    let read3 = build_sam_record_quals(40, read_cigar, &read_bases3, default_quals(read_bases3.len()));
    let read3_clone = build_sam_record_quals(40, read_cigar, &read_bases3, default_quals(read_bases3.len()));

    // a read beyond the delete but still considered
    let read_bases4 = format!("{}{}", sub(REF_BASES_200, 47, 48), sub(REF_BASES_200, 52, 80));
    let read_cigar = "29M"; // could have been soft-clipped or an SNV
    let read4 = build_sam_record_quals(51, read_cigar, &read_bases4, default_quals(read_bases4.len()));

    let reads = [read1, read1_clone, read2, read2_clone, read3, read3_clone, read4];

    // RegionTaskTester + RegionTask.run(). The Rust region task (main.rs process_region) reads an indexed BAM and
    // FASTA, and htslib rejects read2 (36 bases against an 18M4D17M cigar, which htsjdk accepts), so the mock
    // slicer is replaced by feeding the reads to the same stages in the same order as process_region and
    // evidence_stage::collect_evidence: candidates, evidence, hard filters, variants.
    let cfg = test_config();
    let ref_data = empty_ref_data();
    let bqr = BqrMap::default();
    let msi = MsiJitter { params: Vec::new(), probable_msi: false };

    // new RefSequence(region, MockRefGenome): the mock is zero-based, so position == string index
    let ref_seq = RefSequence::from_bases(1, ref_bases.as_bytes()[1..=1300].to_vec());

    // CandidateStage.findCandidates: no hotspots, panel or high confidence regions
    let mut finder = CandidateFinder::new(&cfg, CHR_1, region, &ref_seq, &[], &[], &ref_data);
    for read in &reads {
        finder.process_read(read);
    }
    let candidates = build_candidates(finder.finish(), &[], &[], &[]);

    // EvidenceStage.findEvidence for the tumor
    let qc = QualityCalculator { cfg: &cfg, bqr: &bqr, ref_seq: &ref_seq, msi: &msi, ref_data: &ref_data, chr: CHR_1 };
    let mut counters: Vec<ReadContextCounter> = candidates
        .iter()
        .map(|c| ReadContextCounter::new(c.ctx.clone(), c.tier, max_coverage(&cfg, c.tier), c.min_events, &qc, false))
        .collect();
    let deletes: Vec<(i32, i32)> =
        candidates.iter().filter(|c| c.variant().is_delete()).map(|c| (c.position(), c.variant().ref_len() - 1)).collect();
    for rc in counters.iter_mut() {
        let max_close_del =
            deletes.iter().filter(|(p, _)| (p - rc.position()).abs() < 50).map(|(_, l)| *l).max().unwrap_or(0);
        if max_close_del >= 5 {
            rc.max_candidate_delete_length = max_close_del;
        }
    }
    let mut phaser = Phaser::new();
    let mut phase_sink: Option<&mut dyn PhaseSink> = Some(&mut phaser);
    let mut collector = Collector { counters, last_index: 0, selected: Vec::new(), ref_seq: &ref_seq, qc: &qc };
    let mut sync = FragmentSync::new();
    let mut sync_out: SyncOutput = Vec::new();
    assert!(cfg.sync_fragments);
    for read in reads {
        match sync.handle(read, &mut sync_out) {
            Some(read) => collector.process(&read, None, &mut phase_sink),
            None => {
                for (r, f) in sync_out.drain(..) {
                    collector.process(&r, f.as_ref(), &mut phase_sink);
                }
            }
        }
    }
    for read in sync.empty_cached_reads() {
        collector.process(&read, None, &mut phase_sink);
    }
    let mut counters = collector.counters;
    for rc in counters.iter_mut() {
        set_jitter_qual_filter_state(&msi, rc);
    }

    // tumorEvidence.filterCandidates() and the SageVariant list returned by task.getVariants()
    let mut variants: Vec<SageVariant> = Vec::new();
    for (candidate, rc) in candidates.iter().zip(counters) {
        if passes_hard_filters(&cfg, &rc) {
            variants.push(SageVariant::new(candidate.clone(), Vec::new(), vec![rc]));
        }
    }

    assert_eq!(3, variants.len());
    let del = variants.iter().find(|x| x.position() == 47);
    let var2 = variants.iter().find(|x| x.position() == 50);
    let var3 = variants.iter().find(|x| x.position() == 55);
    assert!(del.is_some());
    assert!(var2.is_some());
    assert!(var3.is_some());

    let del_rc_counter = &del.unwrap().tumor_counters[0];

    assert_eq!(2, del_rc_counter.counts.full);
    assert_eq!(7, del_rc_counter.counts.total);
    assert_eq!(7, del_rc_counter.depth());
}

// Java: testNearByIndels
#[test]
fn test_near_by_indels() {
    let fixture = Fixture::new();

    let var1 = create_sage_variant(90, "A", "C", &fixture);
    let var2 = create_sage_variant(100, "A", "C", &fixture);
    let indel1 = create_sage_variant(110, "A", "AC", &fixture);
    let var3 = create_sage_variant(120, "A", "C", &fixture);
    let var4 = create_sage_variant(130, "A", "C", &fixture);

    let mut sage_variants = vec![var1, var2, indel1, var3, var4];

    set_near_by_indel_status_post_filter(&mut sage_variants);

    assert!(!sage_variants[0].near_indel); // var1
    assert!(sage_variants[1].near_indel); // var2
    assert!(sage_variants[3].near_indel); // var3
    assert!(!sage_variants[4].near_indel); // var4

    let var1 = create_sage_variant(100, "A", "C", &fixture);
    let mut indel1 = create_sage_variant(110, "A", "AC", &fixture);
    indel1.filters.insert(SoftFilter::MinTumorVaf);
    let var2 = create_sage_variant(120, "A", "C", &fixture);

    let mut sage_variants = vec![var1, indel1, var2];

    set_near_by_indel_status_post_filter(&mut sage_variants);

    assert!(!sage_variants[0].near_indel); // var1
    assert!(!sage_variants[2].near_indel); // var2

    let var1 = create_sage_variant(100, "A", "C", &fixture);
    let mut indel1 = create_sage_variant(110, "A", "AC", &fixture);
    indel1.filters.insert(SoftFilter::MaxGermlineVaf);
    let var2 = create_sage_variant(120, "A", "C", &fixture);

    let mut sage_variants = vec![var1, indel1, var2];

    set_near_by_indel_status_post_filter(&mut sage_variants);

    assert!(sage_variants[0].near_indel); // var1
    assert!(sage_variants[2].near_indel); // var2
}
