//! The switch between a port's corrected behaviour (the default) and that of the Java tool.
//!
//! Each port numbers the faults it found in the Java code in its REPORT.md (S1 to S53 for SAGE, R1 to
//! R47 for REDUX, and so on) and asks `java("S4")` wherever it corrects one: false, the default, takes
//! the corrected path; true takes the Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `<TOOL>_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except S4,S7` (or `<TOOL>_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
//!   those listed, to see what one correction changes.
//!
//! The port's own `compat` module holds its finding ids and test hooks and calls `take_options` in
//! `main`, which removes both options from the command line before the program parses it.

use anyhow::{bail, Result};
use std::collections::HashSet;
use std::sync::atomic::{AtomicU8, Ordering};
use std::sync::OnceLock;

const FIXED: u8 = 0;
const JAVA: u8 = 1;
const JAVA_EXCEPT: u8 = 2;

static MODE: AtomicU8 = AtomicU8::new(FIXED);
static EXCEPT: OnceLock<HashSet<String>> = OnceLock::new();

/// The findings of hmf-common-rs's own code, SG1 to SG<SEGMENTATION_FINDINGS> (segmentation; see the
/// segmentation module). A port whose code calls into it accepts these ids in its `valid_id`.
pub const SEGMENTATION_FINDINGS: u32 = 3;

/// Whether `id` is a finding of hmf-common-rs (SG1 to SG<SEGMENTATION_FINDINGS>).
pub fn common_id(id: &str) -> bool {
    numbered_id(id, "SG", SEGMENTATION_FINDINGS)
}

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// For this crate's tests: runs `f` with Java behaviour on (`true`) or off for every finding, on this
/// thread only.
#[cfg(test)]
pub fn with_java<T>(on: bool, f: impl FnOnce() -> T) -> T {
    struct Reset(Option<bool>);
    impl Drop for Reset {
        fn drop(&mut self) {
            TEST_MODE.with(|m| m.set(self.0));
        }
    }
    let _reset = Reset(TEST_MODE.with(|m| m.replace(Some(on))));
    f()
}

/// Whether to behave as Java does for the finding `id`, as set by `take_options`.
#[inline]
pub fn java(id: &str) -> bool {
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    match MODE.load(Ordering::Relaxed) {
        FIXED => false,
        JAVA => true,
        _ => !EXCEPT.get().is_some_and(|s| s.contains(id)),
    }
}

/// True when every finding is reproduced.
#[inline]
pub fn java_all() -> bool {
    MODE.load(Ordering::Relaxed) == JAVA
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode; the
/// environment variables `<tool>_JAVA_COMPAT_EXCEPT` and `<tool>_JAVA_COMPAT` apply, in that order, when
/// neither option is given. `valid_id` rejects ids that are not findings of the port.
pub fn take_options(argv: &[String], tool: &str, valid_id: fn(&str) -> bool) -> Result<Vec<String>> {
    let mut rest = Vec::with_capacity(argv.len());
    let mut all = false;
    let mut except: Option<String> = None;
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-java_compat" => all = true,
            "-java_compat_except" => {
                i += 1;
                let Some(v) = argv.get(i) else { bail!("missing value for -java_compat_except") };
                except = Some(v.clone());
            }
            _ => rest.push(argv[i].clone()),
        }
        i += 1;
    }
    if !all && except.is_none() {
        if let Ok(v) = std::env::var(format!("{tool}_JAVA_COMPAT_EXCEPT")) {
            except = Some(v);
        } else if std::env::var(format!("{tool}_JAVA_COMPAT")).is_ok_and(|v| !v.is_empty() && v != "0") {
            all = true;
        }
    }
    if let Some(list) = except {
        let mut ids = HashSet::new();
        for id in list.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            if !valid_id(id) {
                bail!("-java_compat_except: unknown finding id {id}");
            }
            ids.insert(id.to_string());
        }
        let _ = EXCEPT.set(ids);
        MODE.store(JAVA_EXCEPT, Ordering::Relaxed);
    } else if all {
        MODE.store(JAVA, Ordering::Relaxed);
    }
    Ok(rest)
}

/// Whether `id` is `prefix` followed by a number from 1 to `max`, digits only ("S7" for prefix "S").
pub fn numbered_id(id: &str, prefix: &str, max: u32) -> bool {
    let Some(rest) = id.strip_prefix(prefix) else { return false };
    if rest.is_empty() || !rest.bytes().all(|b| b.is_ascii_digit()) {
        return false;
    }
    rest.parse::<u32>().is_ok_and(|n| (1..=max).contains(&n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids() {
        assert!(numbered_id("S1", "S", 53) && numbered_id("S53", "S", 53));
        assert!(!numbered_id("S54", "S", 53) && !numbered_id("S0", "S", 53) && !numbered_id("X1", "S", 53));
        assert!(!numbered_id("S", "S", 53) && !numbered_id("", "S", 53) && !numbered_id("S+1", "S", 53));
        assert!(common_id("SG1") && common_id(&format!("SG{SEGMENTATION_FINDINGS}")));
        assert!(!common_id(&format!("SG{}", SEGMENTATION_FINDINGS + 1)) && !common_id("AM1"));
    }

    #[test]
    fn options_are_removed() {
        let argv: Vec<String> = ["-threads", "2", "-output_dir", "x"].iter().map(|s| s.to_string()).collect();
        assert_eq!(take_options(&argv, "HMF_COMMON_TEST", |id| numbered_id(id, "T", 3)).unwrap(), argv);
        let bad: Vec<String> = ["-java_compat_except", "Z9"].iter().map(|s| s.to_string()).collect();
        assert!(take_options(&bad, "HMF_COMMON_TEST", |id| numbered_id(id, "T", 3)).is_err());
        assert!(!java("T1"));
    }
}
