//! Switch between the corrected behaviour (the default) and that of Java LINX 2.4.
//!
//! REPORT.md numbers the faults found in the Java code (LX1, LX2, ...). Each place where the port corrects
//! one asks `compat::java("LX1")`: false, the default, takes the corrected path; true takes the Java path,
//! so that the output can be compared with Java's.
//!
//! - `-java_compat` (or the environment variable `LINX_JAVA_COMPAT=1`): Java behaviour for every finding.
//! - `-java_compat_except LX1,LX2` (or `LINX_JAVA_COMPAT_EXCEPT`): Java behaviour for every finding but
//!   those listed.

use anyhow::Result;
use hmf_common::compat::numbered_id;

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Whether to behave as Java does for the finding `id` (for example "LX1").
#[inline]
pub fn java(id: &'static str) -> bool {
    debug_assert!(valid_id(id), "unknown finding id {id}");
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    hmf_common::compat::java(id)
}

/// The number of findings, LX1 to LX<FINDINGS>; raise it when a finding is added.
pub const FINDINGS: u32 = 5;

/// The visualiser's findings, LV1 to LV<VIS_FINDINGS>.
pub const VIS_FINDINGS: u32 = 4;

fn valid_id(id: &str) -> bool {
    numbered_id(id, "LX", FINDINGS) || numbered_id(id, "LV", VIS_FINDINGS) || hmf_common::compat::common_id(id)
}

/// Takes `-java_compat` and `-java_compat_except <ids>` out of the arguments and sets the mode.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    hmf_common::compat::take_options(argv, "LINX", valid_id)
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
