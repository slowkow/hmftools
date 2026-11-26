//! Switch between the corrected behaviour (the default) and that of Java ESVEE 2.1.
//!
//! REPORT.md, section 4, numbers the faults found in the Java code (P prep, A assembly, D depth,
//! C caller). Each place where the port corrects one asks `compat::java("A26")`: false, the default,
//! takes the corrected path; true takes the Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `ESVEE_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except A26,P3` (or `ESVEE_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
//!   those listed, to see what one correction changes.
//!
//! Both options are taken off the command line in `main` before the step parses it, by
//! `hmf_common::compat`, which holds the mode.

use anyhow::Result;
use hmf_common::compat::numbered_id;

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Whether to behave as Java does for the finding `id` of REPORT.md section 4 (for example "A26").
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

/// The findings per step: P1-P13, A1-A41, D1-D19, C1-C18.
fn valid_id(id: &str) -> bool {
    [("P", 13), ("A", 41), ("D", 19), ("C", 18)].iter().any(|&(step, max)| numbered_id(id, step, max))
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode; the
/// environment variables `ESVEE_JAVA_COMPAT_EXCEPT` and `ESVEE_JAVA_COMPAT` apply when neither option is given.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "ESVEE", valid_id)
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
        assert!(valid_id("P1") && valid_id("A41") && valid_id("D19") && valid_id("C18"));
        assert!(!valid_id("A42") && !valid_id("C19") && !valid_id("P0") && !valid_id("X1") && !valid_id("A") && !valid_id(""));
    }

    #[test]
    fn default_is_fixed_and_scope_overrides() {
        assert!(!java("A26"));
        assert!(with_java(true, || java("A26")));
        assert!(!with_java(true, || with_java(false, || java("A26"))));
        assert!(!java("A26"));
    }

    #[test]
    fn options_are_removed() {
        let argv: Vec<String> = ["-threads", "2", "-output_dir", "x"].iter().map(|s| s.to_string()).collect();
        assert_eq!(take_options(&argv).unwrap(), argv);
        let bad: Vec<String> = ["-java_compat_except", "Z9"].iter().map(|s| s.to_string()).collect();
        assert!(take_options(&bad).is_err());
    }
}
