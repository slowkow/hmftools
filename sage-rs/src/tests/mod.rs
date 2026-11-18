//! Unit tests ported from the Java SAGE test suite (sage/src/test/java/com/hartwig/hmftools/sage).
//! One module per Java test class; shared fixtures are in `util`.

#![allow(dead_code, unused_imports, unexpected_cfgs)]

// `RUSTFLAGS='--cfg test_group="a"'` builds one group of modules only, so that groups being worked on
// at the same time do not break each other's builds. Without it, everything is built.

pub mod util;

#[cfg(any(not(test_group), test_group = "c"))]
mod candidate_creation; // candidate/CandidateCreationTest
#[cfg(any(not(test_group), test_group = "c"))]
mod candidate_serialisation; // candidate/CandidateSerialisationTest
#[cfg(any(not(test_group), test_group = "b"))]
mod evicting_array; // common/EvictingArrayTest
#[cfg(any(not(test_group), test_group = "d"))]
mod fragment_sync; // evidence/FragmentSyncTest
#[cfg(any(not(test_group), test_group = "a"))]
mod homology; // common/HomologyTest
#[cfg(any(not(test_group), test_group = "c"))]
mod indel_dedup; // dedup/IndelDedupTest
#[cfg(any(not(test_group), test_group = "d"))]
mod jitter; // evidence/JitterTest
#[cfg(any(not(test_group), test_group = "e"))]
mod misc_evidence; // evidence/MiscEvidenceTest
#[cfg(any(not(test_group), test_group = "f"))]
mod phasing_groups; // evidence/PhasingGroupsTest
#[cfg(any(not(test_group), test_group = "b"))]
mod read_context_matching; // common/ReadContextMatchingTest
#[cfg(any(not(test_group), test_group = "a"))]
mod read_context_utils; // common/ReadContextUtilsTest
#[cfg(any(not(test_group), test_group = "e"))]
mod realignment; // evidence/RealignmentTest
#[cfg(any(not(test_group), test_group = "c"))]
mod ref_context_consumer; // candidate/RefContextConsumerTest
#[cfg(any(not(test_group), test_group = "a"))]
mod ref_sequence; // common/RefSequenceTest
#[cfg(any(not(test_group), test_group = "a"))]
mod repeat_info; // common/RepeatInfoTest
#[cfg(any(not(test_group), test_group = "f"))]
mod soft_filter; // filters/SoftFilterTest
#[cfg(any(not(test_group), test_group = "b"))]
mod split_read; // evidence/SplitReadTest
#[cfg(any(not(test_group), test_group = "c"))]
mod variant_dedup; // dedup/VariantDedupTest
#[cfg(any(not(test_group), test_group = "a"))]
mod variant_read_context; // common/VariantReadContextTest
#[cfg(any(not(test_group), test_group = "b"))]
mod variant_selectors; // common/VariantSelectorsTest
#[cfg(any(not(test_group), test_group = "c"))]
mod variant_types; // candidate/VariantTypesTest
