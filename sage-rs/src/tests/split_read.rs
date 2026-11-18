//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/evidence/SplitReadTest.java
//!
//! Nothing is ported: split-read (RNA, `N` in the CIGAR) handling is not implemented in the Rust port.
//! `evidence/SplitReadSegment` has no Rust equivalent and `ReadContextCounter::process_read` returns
//! `ReadMatchType::InSplit` for any read with an `N` element (src/evidence.rs, "split reads (RNA) are not
//! supported"; also listed in QUIRKS.md under "Where the port differs").
//!
//! Skipped: testSplitReadSegments - exercises SplitReadSegment.formSegment, which is not ported.
//!
//! Skipped: testSplitReads - passes reads with `N` sections to ReadContextCounter.processRead and expects
//! Full / PartialCore support counted from the segment covering the variant; the Rust counter skips these
//! reads, so the test needs the unported SplitReadSegment feature.
