//! Switch between the corrected behaviour (the default) and that of Java SAGE 5.1.
//!
//! REPORT.md, section 5, numbers the faults found in the Java code (S1 to S60). Each place where the
//! port corrects one asks `compat::java("S4")`: false, the default, takes the corrected path; true takes
//! the Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `SAGE_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except S4,S7` (or `SAGE_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
//!   those listed, to see what one correction changes.
//!
//! Both options are taken off the command line in `main` before the program parses it, by
//! `hmf_common::compat`, which holds the mode.

use anyhow::Result;
use hmf_common::compat::numbered_id;

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
    hmf_common::compat::java(id)
}

/// True when every finding is reproduced.
#[allow(dead_code)]
pub fn java_all() -> bool {
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    hmf_common::compat::java_all()
}

/// The number of findings: S1 to S60.
pub const FINDINGS: u32 = 60;

fn valid_id(id: &str) -> bool {
    numbered_id(id, "S", FINDINGS)
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode; the
/// environment variables `SAGE_JAVA_COMPAT_EXCEPT` and `SAGE_JAVA_COMPAT` apply when neither option is given.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "SAGE", valid_id)
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
