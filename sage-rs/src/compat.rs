//! Switch between the corrected behaviour (the default) and that of Java SAGE 5.1.
//!
//! REPORT.md, section 5, numbers the faults found in the Java code (S1 to S53). Each place where the
//! port corrects one asks `compat::java("S4")`: false, the default, takes the corrected path; true takes
//! the Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `SAGE_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except S4,S7` (or `SAGE_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
//!   those listed, to see what one correction changes.
//!
//! Both options are taken off the command line in `main` before the program parses it.

use anyhow::{bail, Result};
use std::collections::HashSet;
use std::sync::atomic::{AtomicU8, Ordering};
use std::sync::OnceLock;

const FIXED: u8 = 0;
const JAVA: u8 = 1;
const JAVA_EXCEPT: u8 = 2;

static MODE: AtomicU8 = AtomicU8::new(FIXED);
static EXCEPT: OnceLock<HashSet<String>> = OnceLock::new();

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Whether to behave as Java does for the finding `id` of REPORT.md section 5 (for example "S4").
#[inline]
pub fn java(id: &'static str) -> bool {
    debug_assert!(valid_id(id), "unknown finding id {id}");
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
#[allow(dead_code)]
pub fn java_all() -> bool {
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    MODE.load(Ordering::Relaxed) == JAVA
}

/// The number of findings: S1 to S53.
pub const FINDINGS: u32 = 53;

fn valid_id(id: &str) -> bool {
    let Some(rest) = id.strip_prefix('S') else { return false };
    if rest.is_empty() || !rest.bytes().all(|b| b.is_ascii_digit()) {
        return false;
    }
    rest.parse::<u32>().is_ok_and(|n| (1..=FINDINGS).contains(&n))
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode; the
/// environment variables apply when neither option is given.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
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
        if let Ok(v) = std::env::var("SAGE_JAVA_COMPAT_EXCEPT") {
            except = Some(v);
        } else if std::env::var("SAGE_JAVA_COMPAT").is_ok_and(|v| !v.is_empty() && v != "0") {
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

/// For tests: runs `f` with Java behaviour on (`true`) or off for every finding, on this thread only.
#[cfg(test)]
#[allow(dead_code)]
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

/// For tests: Java behaviour on or off for every finding until the guard is dropped, on this thread only.
#[cfg(test)]
#[allow(dead_code)]
pub fn scope(on: bool) -> impl Drop {
    struct Reset(Option<bool>);
    impl Drop for Reset {
        fn drop(&mut self) {
            TEST_MODE.with(|m| m.set(self.0));
        }
    }
    Reset(TEST_MODE.with(|m| m.replace(Some(on))))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids() {
        assert!(valid_id("S1") && valid_id(&format!("S{FINDINGS}")));
        assert!(!valid_id(&format!("S{}", FINDINGS + 1)) && !valid_id("S0") && !valid_id("X1") && !valid_id("S") && !valid_id("") && !valid_id("S+1"));
    }

    #[test]
    fn default_is_fixed_and_scope_overrides() {
        assert!(!java("S4"));
        assert!(with_java(true, || java("S4")));
        assert!(!with_java(true, || with_java(false, || java("S4"))));
        assert!(!java("S4"));
    }

    #[test]
    fn options_are_removed() {
        let argv: Vec<String> = ["-threads", "2", "-output_vcf", "x"].iter().map(|s| s.to_string()).collect();
        assert_eq!(take_options(&argv).unwrap(), argv);
        let bad: Vec<String> = ["-java_compat_except", "Z9"].iter().map(|s| s.to_string()).collect();
        assert!(take_options(&bad).is_err());
    }
}
