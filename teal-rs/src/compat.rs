//! Switch between the corrected behaviour (the default) and that of Java TEAL 1.4.
//!
//! REPORT.md numbers the faults found in the Java code (TE1 to TE<FINDINGS>). Each place where the port
//! corrects one asks `compat::java("TE1")`: false, the default, takes the corrected path; true takes the
//! Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `TEAL_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except TE1,TE2` (or `TEAL_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
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

/// Whether to behave as Java does for the finding `id` of REPORT.md section 3 (for example "TE1").
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

/// The number of findings in REPORT.md section 3, TE1 to TE<FINDINGS>; raise it when a finding is added.
pub const FINDINGS: u32 = 11;

/// TE1-TE<FINDINGS> of this port, and the ids of the shared hmf-common-rs code (SG1-SG3, segmentation)
fn valid_id(id: &str) -> bool {
    numbered_id(id, "TE", FINDINGS) || hmf_common::compat::common_id(id)
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode; the
/// environment variables `TEAL_JAVA_COMPAT_EXCEPT` and `TEAL_JAVA_COMPAT` apply when neither option is given.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "TEAL", valid_id)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids() {
        assert!(valid_id("TE1") && valid_id(&format!("TE{FINDINGS}")) && valid_id("SG1"));
        assert!(!valid_id(&format!("TE{}", FINDINGS + 1)) && !valid_id("TE0") && !valid_id("AM1") && !valid_id("TE") && !valid_id(""));
    }

    #[test]
    fn options_are_removed() {
        let argv: Vec<String> = ["-threads", "2", "-output_dir", "x"].iter().map(|s| s.to_string()).collect();
        assert_eq!(take_options(&argv).unwrap(), argv);
        let bad: Vec<String> = ["-java_compat_except", "Z9"].iter().map(|s| s.to_string()).collect();
        assert!(take_options(&bad).is_err());
    }
}
