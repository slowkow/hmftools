//! Which Java ISOFOX the run reproduces: 2.1 (the default; RNA aligned by bwa-mem3 and lifted over by TARS,
//! duplicates marked by REDUX) or 2.0.1 (RNA aligned by STAR, duplicates marked by GATK MarkDuplicates, as
//! Hartwig's README asks for STAR BAMs: "use Isofox v2.0 or earlier").
//!
//! Selection, first match wins:
//! - `-isofox_version 2.0.1` or `-isofox_version 2.1` on the command line;
//! - the environment variable `ISOFOX_VERSION` (`2.0.1` or `2.1`), which `-version` reads too;
//! - the BAM header: an `@PG` line of STAR (`ID:STAR` or `PN:STAR`) selects 2.0.1;
//! - otherwise 2.1.
//!
//! Code paths that differ between the versions ask `version::v201()`.

use anyhow::{bail, Result};
use std::sync::atomic::{AtomicU8, Ordering};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IsofoxVersion {
    V2_0_1,
    V2_1,
}

impl IsofoxVersion {
    /// the version string Java prints (`Isofox version <x>`)
    pub fn java_version(self) -> &'static str {
        match self {
            IsofoxVersion::V2_0_1 => "2.0.1",
            IsofoxVersion::V2_1 => "2.1",
        }
    }
    pub fn parse(s: &str) -> Result<IsofoxVersion> {
        Ok(match s.trim() {
            "2.0.1" => IsofoxVersion::V2_0_1,
            "2.1" => IsofoxVersion::V2_1,
            _ => bail!("invalid ISOFOX version '{s}': expected 2.0.1 or 2.1"),
        })
    }
}

/// 0: not set (2.1), 1: 2.0.1, 2: 2.1
static MODE: AtomicU8 = AtomicU8::new(0);

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

pub fn set(v: IsofoxVersion) {
    MODE.store(if v == IsofoxVersion::V2_0_1 { 1 } else { 2 }, Ordering::Relaxed);
}

#[allow(dead_code)]
pub fn current() -> IsofoxVersion {
    if v201() { IsofoxVersion::V2_0_1 } else { IsofoxVersion::V2_1 }
}

/// true when the run reproduces ISOFOX 2.0.1
#[inline]
pub fn v201() -> bool {
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    MODE.load(Ordering::Relaxed) == 1
}

/// For tests: runs `f` in 2.0.1 mode (`true`) or 2.1 mode, on this thread only.
#[cfg(test)]
pub fn with_v201<T>(on: bool, f: impl FnOnce() -> T) -> T {
    struct Reset(Option<bool>);
    impl Drop for Reset {
        fn drop(&mut self) {
            TEST_MODE.with(|m| m.set(self.0));
        }
    }
    let _reset = Reset(TEST_MODE.with(|m| m.replace(Some(on))));
    f()
}

/// Takes `-isofox_version <v>` out of the arguments; falls back to `ISOFOX_VERSION`. None when neither is given.
pub fn take_option(argv: &[String]) -> Result<(Vec<String>, Option<IsofoxVersion>)> {
    let mut out = Vec::with_capacity(argv.len());
    let mut chosen: Option<IsofoxVersion> = None;
    let mut i = 0;
    while i < argv.len() {
        if argv[i] == "-isofox_version" {
            let Some(v) = argv.get(i + 1) else { bail!("-isofox_version needs a value (2.0.1 or 2.1)") };
            chosen = Some(IsofoxVersion::parse(v)?);
            i += 2;
            continue;
        }
        out.push(argv[i].clone());
        i += 1;
    }
    if chosen.is_none() {
        if let Ok(v) = std::env::var("ISOFOX_VERSION") {
            if !v.is_empty() {
                chosen = Some(IsofoxVersion::parse(&v)?);
            }
        }
    }
    Ok((out, chosen))
}

/// The BAM's aligner from its header: STAR (an @PG line with ID or PN STAR, or STAR.<n> as samtools
/// renames repeated ids) selects 2.0.1. None when the header names no STAR program.
pub fn detect_from_header_text(header: &str) -> Option<IsofoxVersion> {
    for line in header.lines() {
        if !line.starts_with("@PG") {
            continue;
        }
        for field in line.split('\t').skip(1) {
            let Some((tag, value)) = field.split_once(':') else { continue };
            if (tag == "ID" || tag == "PN") && (value == "STAR" || value.starts_with("STAR.")) {
                return Some(IsofoxVersion::V2_0_1);
            }
        }
    }
    None
}

/// The value of `-bam_file` in the arguments, if any
pub fn bam_file_arg(argv: &[String]) -> Option<&str> {
    argv.iter().position(|a| a == "-bam_file").and_then(|i| argv.get(i + 1)).map(|s| s.as_str())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_and_take() {
        let argv: Vec<String> = ["-sample", "S", "-isofox_version", "2.0.1", "-threads", "2"].iter().map(|s| s.to_string()).collect();
        let (rest, v) = take_option(&argv).unwrap();
        assert_eq!(v, Some(IsofoxVersion::V2_0_1));
        assert_eq!(rest, vec!["-sample", "S", "-threads", "2"]);
        assert!(IsofoxVersion::parse("2.0").is_err());
        assert_eq!(IsofoxVersion::parse("2.1").unwrap().java_version(), "2.1");
    }

    #[test]
    fn star_header_detection() {
        let h = "@HD\tVN:1.4\n@PG\tID:STAR\tPN:STAR\tVN:2.7.3a\tCL:STAR --runThreadN 8\n@CO\tuser command line: STAR";
        assert_eq!(detect_from_header_text(h), Some(IsofoxVersion::V2_0_1));
        let h2 = "@HD\tVN:1.6\n@PG\tID:bwa-mem3\tPN:bwa-mem3\n@PG\tID:samtools\tPN:samtools";
        assert_eq!(detect_from_header_text(h2), None);
        assert_eq!(detect_from_header_text("@PG\tID:STAR.1\tPN:other"), Some(IsofoxVersion::V2_0_1));
        assert_eq!(detect_from_header_text("@CO\tID:STAR"), None);
    }
}
