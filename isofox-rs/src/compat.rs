//! Switch between the corrected behaviour (the default) and that of Java ISOFOX 2.1.
//!
//! REPORT.md numbers the faults found in the Java code (IF1, IF2, ...). Each place where the port
//! corrects one asks `compat::java("IF1")`: false, the default, takes the corrected path; true takes
//! the Java path, so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `ISOFOX_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except IF1,IF2` (or `ISOFOX_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
//!   those listed, to see what one correction changes.

use anyhow::Result;
use hmf_common::compat::numbered_id;

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// The number of findings in REPORT.md, IF1 to IF<FINDINGS>; raise it when a finding is added.
pub const FINDINGS: u32 = 13;

fn valid_id(id: &str) -> bool {
    numbered_id(id, "IF", FINDINGS) || hmf_common::compat::common_id(id)
}

/// Whether to behave as Java does for the finding `id` (for example "IF1").
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

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "ISOFOX", valid_id)
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
