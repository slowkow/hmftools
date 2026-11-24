//! The base-quality bands of prep and assembly.
//!
//! Java fixes them for Illumina data: a base at or below BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD
//! (Q26) is low quality (prep drops a soft clip with 5 or more such bases, and assembly lets them
//! mismatch), and a base at or above IlluminaBamUtils.ILLUMINA_HIGH_BASE_QUAL (Q30) is high quality.
//! Data whose qualities were rewritten by BQSR can stop below Q30, with most bases in the low band, and
//! then the junctions in that sample lose their support. -low_base_qual_threshold and
//! -high_base_qual_threshold (not Java options) move the two bands for a run. The defaults are Java's,
//! so output is unchanged when they are not given, with or without -java_compat.

use crate::common::args::Args;
use anyhow::{bail, Result};
use std::sync::atomic::{AtomicU8, Ordering::Relaxed};

/// BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD
pub const DEFAULT_LOW_BASE_QUAL_THRESHOLD: u8 = 26;
/// IlluminaBamUtils.ILLUMINA_HIGH_BASE_QUAL
pub const DEFAULT_HIGH_BASE_QUAL_THRESHOLD: u8 = 30;

pub const LOW_OPTION: &str = "low_base_qual_threshold";
pub const HIGH_OPTION: &str = "high_base_qual_threshold";

static LOW: AtomicU8 = AtomicU8::new(DEFAULT_LOW_BASE_QUAL_THRESHOLD);
static HIGH: AtomicU8 = AtomicU8::new(DEFAULT_HIGH_BASE_QUAL_THRESHOLD);

#[cfg(test)]
thread_local! {
    /// Tests run in parallel threads, so a test sets the bands for its own thread only.
    static TEST_BANDS: std::cell::Cell<Option<(u8, u8)>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
pub fn set_test_base_qual_thresholds(bands: Option<(u8, u8)>) {
    TEST_BANDS.with(|c| c.set(bands));
}

/// Bases at or below this quality are low quality.
#[inline]
pub fn low_base_qual_threshold() -> u8 {
    #[cfg(test)]
    if let Some((low, _)) = TEST_BANDS.with(|c| c.get()) {
        return low;
    }
    LOW.load(Relaxed)
}

/// Illumina bases at or above this quality are high quality.
#[inline]
pub fn high_base_qual_threshold() -> u8 {
    #[cfg(test)]
    if let Some((_, high)) = TEST_BANDS.with(|c| c.get()) {
        return high;
    }
    HIGH.load(Relaxed)
}

/// The two options, checked: Phred qualities 0 to 93, and no base both low and high quality.
pub fn parse(low: Option<&str>, high: Option<&str>) -> Result<(u8, u8)> {
    let one = |name: &str, v: Option<&str>, default: u8| -> Result<u8> {
        match v {
            None => Ok(default),
            Some(s) => match s.parse::<u8>() {
                Ok(q) if q <= 93 => Ok(q),
                _ => bail!("-{name} must be a base quality from 0 to 93, not {s}"),
            },
        }
    };
    let low = one(LOW_OPTION, low, DEFAULT_LOW_BASE_QUAL_THRESHOLD)?;
    let high = one(HIGH_OPTION, high, DEFAULT_HIGH_BASE_QUAL_THRESHOLD)?;
    if low >= high {
        bail!("-{LOW_OPTION} ({low}) must be below -{HIGH_OPTION} ({high})");
    }
    Ok((low, high))
}

/// Reads the options of a step and sets the bands for the run (before any read is processed).
pub fn set_from_args(args: &Args) -> Result<()> {
    let (low, high) = parse(args.value(LOW_OPTION), args.value(HIGH_OPTION))?;
    if (low, high) != (low_base_qual_threshold(), high_base_qual_threshold()) {
        crate::info!(
            "base-quality bands: low at or below Q{low}, high at or above Q{high} (defaults Q{DEFAULT_LOW_BASE_QUAL_THRESHOLD}, Q{DEFAULT_HIGH_BASE_QUAL_THRESHOLD})"
        );
    }
    LOW.store(low, Relaxed);
    HIGH.store(high, Relaxed);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_java_constants() {
        assert_eq!(parse(None, None).unwrap(), (26, 30));
        assert_eq!(parse(Some("20"), Some("24")).unwrap(), (20, 24));
        assert_eq!(parse(None, Some("27")).unwrap(), (26, 27));
    }

    #[test]
    fn rejects_bad_values() {
        assert!(parse(Some("26"), Some("26")).is_err());
        assert!(parse(None, Some("20")).is_err());
        assert!(parse(Some("-1"), None).is_err());
        assert!(parse(None, Some("94")).is_err());
        assert!(parse(Some("Q20"), None).is_err());
    }
}
