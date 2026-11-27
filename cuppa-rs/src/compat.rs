//! Switch between the corrected behaviour (the default) and that of Java CUPPA 2.5 and pycuppa.
//!
//! REPORT.md numbers the faults found in the Java code and pycuppa (CU1, CU2, ...). Each place where
//! the port corrects one asks `compat::java("CU1")`: false, the default, takes the corrected path; true
//! takes the Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `CUPPA_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except CU1,CU2` (or `CUPPA_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
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

/// Whether to behave as Java does for the finding `id` of REPORT.md (for example "CU2").
#[inline]
pub fn java(id: &'static str) -> bool {
    debug_assert!(valid_id(id), "unknown finding id {id}");
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    hmf_common::compat::java(id)
}

/// The number of findings in REPORT.md, CU1 to CU<FINDINGS>; raise it when a finding is added.
pub const FINDINGS: u32 = 8;

fn valid_id(id: &str) -> bool {
    numbered_id(id, "CU", FINDINGS)
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode; the
/// environment variables `CUPPA_JAVA_COMPAT_EXCEPT` and `CUPPA_JAVA_COMPAT` apply when neither is given.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "CUPPA", valid_id)
}

/// For tests: runs `f` with Java behaviour on (`true`) or off for every finding, on this thread only.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids() {
        assert!(valid_id("CU1") && valid_id(&format!("CU{FINDINGS}")));
        assert!(!valid_id(&format!("CU{}", FINDINGS + 1)) && !valid_id("CU0") && !valid_id("PA1") && !valid_id("CU") && !valid_id(""));
    }

    #[test]
    fn options_are_removed() {
        let argv: Vec<String> = ["-threads", "2", "-output_dir", "x"].iter().map(|s| s.to_string()).collect();
        assert_eq!(take_options(&argv).unwrap(), argv);
        let bad: Vec<String> = ["-java_compat_except", "Z9"].iter().map(|s| s.to_string()).collect();
        assert!(take_options(&bad).is_err());
    }
}
