//! Switch between the corrected behaviour (the default) and that of Java qsee 1.0.
//!
//! REPORT.md numbers the faults found in the Java code (QS1, QS2, ...). Each place where the port corrects
//! one asks `compat::java("QS1")`: false, the default, takes the corrected path; true takes the Java path,
//! so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `QSEE_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except QS1,QS2` (or `QSEE_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
//!   those listed.

use anyhow::Result;
use hmf_common::compat::numbered_id;

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Whether to behave as Java does for the finding `id` (for example "QS1").
#[inline]
pub fn java(id: &'static str) -> bool {
    debug_assert!(valid_id(id), "unknown finding id {id}");
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    hmf_common::compat::java(id)
}

/// The number of findings, QS1 to QS<FINDINGS>; raise it when a finding is added.
pub const FINDINGS: u32 = 12;

fn valid_id(id: &str) -> bool {
    numbered_id(id, "QS", FINDINGS)
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "QSEE", valid_id)
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
        assert!(valid_id("QS1") && valid_id(&format!("QS{FINDINGS}")));
        assert!(!valid_id(&format!("QS{}", FINDINGS + 1)) && !valid_id("QS0") && !valid_id("SI1") && !valid_id("QS"));
    }
}
