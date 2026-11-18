//! Unit tests ported from the Java SAGE test suite (sage/src/test/java/com/hartwig/hmftools/sage).
//! One module per Java test class; shared fixtures are in `util`.

#![allow(dead_code, unused_imports)]

pub mod util;

mod artefacts; // seqtech/ArtefactsTest
mod candidate_creation; // candidate/CandidateCreationTest
mod candidate_serialisation; // candidate/CandidateSerialisationTest
mod evicting_array; // common/EvictingArrayTest
mod fragment_sync; // evidence/FragmentSyncTest
mod homology; // common/HomologyTest
mod indel_dedup; // dedup/IndelDedupTest
mod jitter; // evidence/JitterTest
mod misc_evidence; // evidence/MiscEvidenceTest
mod phase_append; // evidence/PhaseAppendTest
mod phasing_groups; // evidence/PhasingGroupsTest
mod read_context_matching; // common/ReadContextMatchingTest
mod read_context_utils; // common/ReadContextUtilsTest
mod realignment; // evidence/RealignmentTest
mod ref_context_consumer; // candidate/RefContextConsumerTest
mod ref_sequence; // common/RefSequenceTest
mod repeat_info; // common/RepeatInfoTest
mod sbx; // seqtech/SbxTest
mod soft_filter; // filters/SoftFilterTest
mod split_read; // evidence/SplitReadTest
mod ultima_evidence; // not from Java: the Ultima evidence and filter path end to end
mod ultima_misc_utils; // seqtech/UltimaMiscUtilsTest
mod ultima_qual_model; // seqtech/UltimaQualModelTest
mod ultima_realigned_qual_model; // seqtech/UltimaRealignedQualModelTest
mod ultima_routines; // seqtech/UltimaRoutinesTest
mod variant_dedup; // dedup/VariantDedupTest
mod variant_read_context; // common/VariantReadContextTest
mod variant_selectors; // common/VariantSelectorsTest
mod variant_types; // candidate/VariantTypesTest
