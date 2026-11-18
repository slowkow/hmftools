//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! evidence/PhasingGroupsTest.java
//!
//! Java identifies read counters by object; the Rust phaser identifies them by their index in the
//! tumor counter list, so each test builds its counters first and registers lists of indices.

use super::util::*;
use crate::candidate::Candidate;
use crate::constants::*;
use crate::evidence::ReadContextCounter;
use crate::evidence_stage::PhaseSink;
use crate::filters::{SageVariant, SoftFilter};
use crate::phase::*;
use crate::read_context::VariantReadContext;
use crate::types::*;
use std::sync::Arc;

/// VariantUtils.TEST_LEFT_FLANK etc
pub(super) const TEST_LEFT_FLANK: &str = "ACCGCTGACT";
pub(super) const TEST_RIGHT_FLANK: &str = "CTGAGACTCA";
pub(super) const TEST_LEFT_CORE: &str = "AC";
pub(super) const TEST_RIGHT_CORE: &str = "GT";

/// VariantUtils.createReadContext(variant, leftCore, rightCore, leftFlank, rightFlank): no repeat or homology
pub(super) fn create_read_context_flanks(
    variant: SimpleVariant, left_core: &str, right_core: &str, left_flank: &str, right_flank: &str,
) -> Arc<VariantReadContext> {
    let core_index_start = left_flank.len() as i32;
    let var_index = core_index_start + left_core.len() as i32;
    let core_index_end = var_index + variant.alt_len() - 1 + right_core.len() as i32;
    let ref_bases = format!("{}{}{}", left_core, variant.ref_str(), right_core);
    let read_bases = format!("{}{}{}{}{}", left_flank, left_core, variant.alt_str(), right_core, right_flank);
    let alignment_start = variant.pos - var_index;
    let core_right_pos_start = variant.pos + variant.ref_len().min(variant.alt_len());
    let alignment_end = core_right_pos_start + right_core.len() as i32 + right_flank.len() as i32 - 1;
    let core_position_start = variant.pos - left_core.len() as i32;
    let core_position_end = variant.pos + right_core.len() as i32;
    Arc::new(VariantReadContext {
        alignment_start,
        alignment_end,
        ref_bases: ref_bases.into_bytes(),
        read_cigar: format!("{}M", read_bases.len()),
        read_bases: read_bases.into_bytes(),
        core_index_start,
        var_index,
        core_index_end,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start,
        core_position_end,
        has_indel_in_core: false,
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
        variant,
    })
}

/// VariantUtils.createReadContext(variant, leftCore, rightCore)
pub(super) fn create_read_context_cores(variant: SimpleVariant, left_core: &str, right_core: &str) -> Arc<VariantReadContext> {
    create_read_context_flanks(variant, left_core, right_core, TEST_LEFT_FLANK, TEST_RIGHT_FLANK)
}

/// VariantUtils.createReadContext(variant)
pub(super) fn create_read_context(variant: SimpleVariant) -> Arc<VariantReadContext> {
    create_read_context_cores(variant, TEST_LEFT_CORE, TEST_RIGHT_CORE)
}

/// VariantUtils.createSimpleVariant(position)
fn create_simple_variant(position: i32) -> SimpleVariant {
    SimpleVariant::new(position, b"A", b"C")
}

/// Iteration order of a Java HashMap<Integer, _> whose keys were inserted in the given order.
fn java_hash_order(keys: &[i32]) -> Vec<usize> {
    let mut capacity = 16usize;
    while keys.len() as f64 > capacity as f64 * 0.75 {
        capacity *= 2;
    }
    let mut items: Vec<(usize, usize)> = keys
        .iter()
        .enumerate()
        .map(|(seq, &k)| {
            let h = k as u32;
            (((h ^ (h >> 16)) as usize) & (capacity - 1), seq)
        })
        .collect();
    items.sort();
    items.into_iter().map(|(_, seq)| seq).collect()
}

/// The Java test's mPhaser plus the counters it registers.
struct PhaserTester {
    fixture: Fixture,
    counters: Vec<ReadContextCounter>,
    phaser: Phaser,
}

impl PhaserTester {
    fn new() -> Self {
        PhaserTester { fixture: Fixture::new(), counters: Vec::new(), phaser: Phaser::new() }
    }

    /// mPhaser.initialise(region, sample)
    fn initialise(&mut self) {
        self.phaser = Phaser::new();
    }

    /// createReadCounter(position): returns the counter's index, the Rust identity of a counter
    fn create_read_counter(&mut self, position: i32) -> usize {
        let ctx = create_read_context(create_simple_variant(position));
        let qc = self.fixture.qc();
        let counter = ReadContextCounter::new(ctx, Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &qc, false);
        self.counters.push(counter);
        self.counters.len() - 1
    }

    /// mPhaser.registeredPhasedVariants(posCounters, negCounters)
    fn register(&mut self, positive: &[usize], negative: &[usize]) {
        self.phaser.register(&self.counters, positive, negative);
    }

    /// getPhasedGroupCount: PhasedGroupCollection.groupCount summed
    fn phased_group_count(&self) -> usize {
        self.phaser
            .collections
            .iter()
            .map(|c| if !c.map.is_empty() { c.map.iter().map(|(_, list)| list.len()).sum() } else { c.groups.len() })
            .sum()
    }

    /// getPhasedGroups: each collection's groupsMap().values(), in Java HashMap iteration order
    fn phased_groups(&self) -> Vec<PhasedGroup> {
        let mut out = Vec::new();
        for c in &self.phaser.collections {
            let keys: Vec<i32> = c.map.iter().map(|(k, _)| *k).collect();
            for index in java_hash_order(&keys) {
                out.extend(c.map[index].1.iter().cloned());
            }
        }
        out
    }

    /// checkGroupsAreOrdered
    fn check_groups_are_ordered(&self) -> bool {
        let groups = self.phased_groups();
        groups.windows(2).all(|w| w[1].pos_min >= w[0].pos_min)
    }

    fn positions(&self) -> Vec<i32> {
        self.phaser.pos.clone()
    }
}

/// findGroup, over the given groups
fn find_group<'a>(groups: &'a [PhasedGroup], pos_counters: &[usize]) -> Option<&'a PhasedGroup> {
    groups
        .iter()
        .filter(|x| x.positive.len() == pos_counters.len())
        .find(|x| x.positive.iter().all(|y| pos_counters.contains(y)))
}

/// In Java the merge routines work on a copied list of the same group objects that the phaser's
/// collections still hold, and findGroup then looks at the collections again: it sees every original
/// group, the surviving ones in their merged state. The Rust groups are values, so rebuild that view:
/// each original group, replaced by its merged state where a group with its id survived.
fn groups_after_merge(original: &[PhasedGroup], merged: &[PhasedGroup]) -> Vec<PhasedGroup> {
    original.iter().map(|g| merged.iter().find(|m| m.id == g.id).unwrap_or(g).clone()).collect()
}

// Java: testGroupCreations
#[test]
fn test_group_creations() {
    // test that groups are created distinctly where no exact matches exist, ordered by start position
    let mut t = PhaserTester::new();
    let rc1 = t.create_read_counter(10);
    let rc2 = t.create_read_counter(10);
    let rc3 = t.create_read_counter(20);
    let rc4 = t.create_read_counter(20);
    let rc5 = t.create_read_counter(30);
    let rc6 = t.create_read_counter(40);

    t.register(&[rc1], &[rc2]);
    assert_eq!(1, t.phased_group_count());

    // negatives cannot diff
    t.register(&[rc1], &[rc2, rc3]);
    assert_eq!(2, t.phased_group_count());

    let groups = t.phased_groups();
    let group = find_group(&groups, &[rc1]);
    assert!(group.is_some());
    assert_eq!(1, group.unwrap().read_count);

    // positive cannot differ
    t.register(&[rc1, rc3], &[rc2]);
    assert_eq!(3, t.phased_group_count());
    assert!(t.check_groups_are_ordered());

    // groups are added in order
    t.register(&[rc5, rc6], &[]);
    assert_eq!(4, t.phased_group_count());
    assert!(t.check_groups_are_ordered());

    t.register(&[rc3, rc5], &[]);
    assert_eq!(5, t.phased_group_count());

    t.register(&[rc2, rc4], &[]);
    assert_eq!(6, t.phased_group_count());

    t.register(&[rc6], &[rc3]);
    assert_eq!(7, t.phased_group_count());

    t.register(&[rc4, rc5], &[]);
    assert_eq!(8, t.phased_group_count());

    t.register(&[rc2, rc5], &[]);
    assert_eq!(9, t.phased_group_count());

    t.phaser.collections.iter_mut().for_each(|x| x.finalise());
    assert!(t.check_groups_are_ordered());
}

// Java: testMergeMatching
#[test]
fn test_merge_matching() {
    // scenarios:
    // remove groups which only have their +ves in other groups with the same +ves
    // single-variant groups dropped if not present in another group
    let mut t = PhaserTester::new();
    let rc0 = t.create_read_counter(10);
    let rc1 = t.create_read_counter(10);
    let rc2 = t.create_read_counter(20);
    let rc3 = t.create_read_counter(20);
    let rc4 = t.create_read_counter(30);
    let rc5 = t.create_read_counter(40);
    let rc6 = t.create_read_counter(50);
    let rc7 = t.create_read_counter(60);

    t.register(&[rc0], &[rc1]);

    // matches group 0, will be merged
    t.register(&[rc0], &[rc1, rc2]);

    // also matches group 0, will be merged
    t.register(&[rc0], &[rc2]);

    // distinct group and will be a subset of another
    t.register(&[rc1, rc2], &[rc3]);

    // superset of 3
    t.register(&[rc1, rc2, rc4], &[rc5]);

    // matches group 4, will be merged
    t.register(&[rc1, rc2, rc4], &[rc6]);

    // finally 2 groups not part of any others
    t.register(&[rc3], &[rc6]);
    t.register(&[rc6], &[rc7]);

    let original = t.phased_groups();
    let mut groups = original.clone();

    assert_eq!(8, groups.len());

    merge_matching(&mut groups, &t.positions());
    assert_eq!(3, groups.len());

    let view = groups_after_merge(&original, &groups);

    let group = find_group(&view, &[rc0]);
    assert!(group.is_some());
    assert_eq!(3, group.unwrap().read_count);
    assert!((0.0 - group.unwrap().allocated).abs() <= 0.1);

    let group = find_group(&view, &[rc1, rc2]);
    assert!(group.is_some());
    assert_eq!(1, group.unwrap().read_count);

    let group = find_group(&view, &[rc1, rc2, rc4]);
    assert!(group.is_some());
    assert_eq!(2, group.unwrap().read_count);
}

// Java: testMergByExtension
#[test]
fn test_merg_by_extension() {
    // merge any group with a common subset of +ves and -ves if it can only be extended in one direction
    let mut t = PhaserTester::new();
    let rc0 = t.create_read_counter(10);
    let rc1 = t.create_read_counter(20);
    let rc2 = t.create_read_counter(30);
    let rc3 = t.create_read_counter(40);
    let rc4 = t.create_read_counter(50);
    let rc5 = t.create_read_counter(60);
    let rc6 = t.create_read_counter(70);

    // 2 sub-groups and 1 extending in either direction, results in a single group
    t.register(&[rc2, rc5], &[rc3, rc4]);

    // another sub-group
    t.register(&[rc2, rc5], &[rc4]);

    // extends down
    t.register(&[rc0, rc1, rc2, rc5], &[]);

    // extends up
    t.register(&[rc2, rc5, rc6], &[]);

    let original = t.phased_groups();
    let mut groups = original.clone();

    assert_eq!(4, groups.len());

    merge_by_extension(&mut groups, &t.positions());
    assert_eq!(2, groups.len());

    let view = groups_after_merge(&original, &groups);
    let group = find_group(&view, &[rc0, rc1, rc2, rc5, rc6]);
    assert!(group.is_some());
    assert_eq!(1, group.unwrap().read_count);
    assert!((1.0 - group.unwrap().allocated).abs() <= 0.1);

    merge_uninformative(&mut groups, &t.positions());
    assert_eq!(1, groups.len());

    // test again but with multiple options up and down preventing a merge
    groups.clear();
    t.initialise();

    // a sub-group
    t.register(&[rc2, rc3, rc4], &[rc4]);

    // conflicting sub-group
    t.register(&[rc2, rc4], &[rc3]);

    // extends down
    t.register(&[rc1, rc2, rc3, rc4], &[]);

    // extends down again differently
    t.register(&[rc0, rc2, rc3, rc4], &[]);

    // extends up
    t.register(&[rc2, rc3, rc4, rc5], &[]);

    // extends up differently
    t.register(&[rc2, rc3, rc4, rc6], &[]);

    let mut groups = t.phased_groups();

    assert_eq!(6, groups.len());

    merge_by_extension(&mut groups, &t.positions());

    assert_eq!(6, groups.len());
}

/// createVariant(position) -> VariantUtils.createSageVariant(position, ref, alt)
fn create_variant(fixture: &Fixture, position: i32) -> SageVariant {
    let ctx = create_read_context(create_simple_variant(position));
    let qc = fixture.qc();
    let counter = ReadContextCounter::new(ctx.clone(), Tier::LowConfidence, 100, 1, &qc, false);
    let candidate = Candidate { tier: Tier::HighConfidence, ctx, min_events: 1, full_match_support: 1, core_match_support: 0 };
    SageVariant::new(candidate, Vec::new(), vec![counter])
}

// Java: testRemoveUninformativeLps
#[test]
fn test_remove_uninformative_lps() {
    let fixture = Fixture::new();
    let mut variants: Vec<SageVariant> = Vec::new();

    // LPS 1 and 2 both have passing var 1
    // LPS 2 has same passing variant as LPS 1 and is uninformative

    // LPS 3, 4 and 5 all have passing vars 4 & 5
    // LPS 4 & 5 has same passing variants as LPS 3 and are uninformative

    // LPS 6 has passing variants 1 and 4 but no overlap so remains

    let mut var1 = create_variant(&fixture, 11);
    var1.tumor_counters[0].add_local_phase_set(1, 10, 0.0);
    var1.tumor_counters[0].add_local_phase_set(2, 1, 0.0);
    var1.tumor_counters[0].add_local_phase_set(6, 1, 0.0);
    variants.push(var1);

    let mut var2 = create_variant(&fixture, 12);
    var2.tumor_counters[0].add_local_phase_set(1, 10, 0.0);
    var2.filters.insert(SoftFilter::MaxGermlineVaf);
    variants.push(var2);

    // a second group, uninformative vs both 1 and 3
    let mut var3 = create_variant(&fixture, 13);
    var3.tumor_counters[0].add_local_phase_set(2, 1, 0.0);
    var3.filters.insert(SoftFilter::MaxGermlineVaf);
    variants.push(var3);

    let mut var4 = create_variant(&fixture, 14);
    var4.tumor_counters[0].add_local_phase_set(3, 10, 0.0);
    var4.tumor_counters[0].add_local_phase_set(4, 5, 0.0);
    var4.tumor_counters[0].add_local_phase_set(5, 1, 0.0);
    var4.tumor_counters[0].add_local_phase_set(6, 1, 0.0);
    variants.push(var4);

    let mut var5 = create_variant(&fixture, 15);
    var5.tumor_counters[0].add_local_phase_set(3, 10, 0.0);
    var5.tumor_counters[0].add_local_phase_set(4, 5, 0.0);
    var5.tumor_counters[0].add_local_phase_set(5, 1, 0.0);

    // a group which is a subset and with less < 25% RCs
    var5.tumor_counters[0].add_local_phase_set(7, 1, 0.0);

    variants.push(var5);

    let mut var6 = create_variant(&fixture, 16);
    var6.tumor_counters[0].add_local_phase_set(3, 10, 0.0);
    var6.filters.insert(SoftFilter::MaxGermlineVaf);
    variants.push(var6);

    let mut var7 = create_variant(&fixture, 17);
    var7.tumor_counters[0].add_local_phase_set(4, 5, 0.0);
    var7.filters.insert(SoftFilter::MaxGermlineVaf);
    variants.push(var7);

    let mut var8 = create_variant(&fixture, 18);
    var8.tumor_counters[0].add_local_phase_set(5, 1, 0.0);
    var8.tumor_counters[0].add_local_phase_set(6, 10, 0.0);
    var8.filters.insert(SoftFilter::MaxGermlineVaf);
    variants.push(var8);

    let mut var9 = create_variant(&fixture, 19);
    var9.tumor_counters[0].add_local_phase_set(7, 1, 0.0);
    var9.filters.insert(SoftFilter::MaxGermlineVaf);
    variants.push(var9);

    // the Java HashSet of the passing variants' phase sets, in its iteration order
    let passing_sets = passing_phase_sets(&variants);
    assert_eq!(7, passing_sets.len());

    remove_uninformative_lps(&mut variants, &passing_sets);

    // variants[n - 1] is the Java varN
    let (var1, var3, var4, var5) = (&variants[0], &variants[2], &variants[3], &variants[4]);
    let (var6, var7, var9) = (&variants[5], &variants[6], &variants[8]);

    assert!(var1.has_lps(1));

    assert!(!var1.has_lps(2));
    assert!(!var3.has_lps(2));

    assert!(var4.has_lps(3));
    assert!(var5.has_lps(3));
    assert!(var6.has_lps(3));

    assert!(!var4.has_lps(4));
    assert!(!var5.has_lps(4));
    assert!(!var7.has_lps(4));

    assert!(!var4.has_lps(5));
    assert!(!var5.has_lps(5));
    assert!(!var9.has_lps(5));

    assert!(var1.has_lps(6));
    assert!(var4.has_lps(6));

    assert!(!var5.has_lps(7));
    assert!(!var9.has_local_phase_sets());
}
