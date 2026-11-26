//! Switch between the corrected behaviour (the default) and that of Java REDUX 2.1.
//!
//! REPORT.md, section 3, numbers the faults found in the Java code (R1 to R47). Each place where the
//! port corrects one asks `compat::java("R7")`: false, the default, takes the corrected path; true
//! takes the Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `REDUX_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except R7,R12` (or `REDUX_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
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

/// Whether to behave as Java does for the finding `id` of REPORT.md section 3 (for example "R7").
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

/// The number of findings in REPORT.md section 3: R1 to R47.
pub const FINDINGS: u32 = 47;

fn valid_id(id: &str) -> bool {
    numbered_id(id, "R", FINDINGS)
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode; the
/// environment variables `REDUX_JAVA_COMPAT_EXCEPT` and `REDUX_JAVA_COMPAT` apply when neither option is given.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "REDUX", valid_id)
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
        assert!(valid_id("R1") && valid_id(&format!("R{FINDINGS}")));
        assert!(!valid_id(&format!("R{}", FINDINGS + 1)) && !valid_id("R0") && !valid_id("X1") && !valid_id("R") && !valid_id(""));
    }

    #[test]
    fn default_is_fixed_and_scope_overrides() {
        assert!(!java("R7"));
        assert!(with_java(true, || java("R7")));
        assert!(!with_java(true, || with_java(false, || java("R7"))));
        assert!(!java("R7"));
    }

    #[test]
    fn options_are_removed() {
        let argv: Vec<String> = ["-threads", "2", "-output_dir", "x"].iter().map(|s| s.to_string()).collect();
        assert_eq!(take_options(&argv).unwrap(), argv);
        let bad: Vec<String> = ["-java_compat_except", "Z9"].iter().map(|s| s.to_string()).collect();
        assert!(take_options(&bad).is_err());
    }
}
