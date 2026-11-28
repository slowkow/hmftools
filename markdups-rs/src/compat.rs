//! Switch between the corrected behaviour (the default) and that of Picard 3.3.0 as GATK 4.6.1.0 runs it.
//!
//! NOTES.md numbers the faults found in Picard's MarkDuplicates (MD1 to MD<FINDINGS>). Each place where the
//! port corrects one asks `compat::java("MD1")`: false, the default, takes the corrected path; true takes
//! Picard's path, so that the output can be compared with GATK's.
//!
//! - `-java_compat` (or `--java_compat`, or the environment variable `MARKDUPS_JAVA_COMPAT=1`): Picard's
//!   behaviour for every finding.
//! - `-java_compat_except MD1,MD3` (or `MARKDUPS_JAVA_COMPAT_EXCEPT`): Picard's behaviour for every finding
//!   but those listed.

use anyhow::Result;
use hmf_common::compat::numbered_id;

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Whether to behave as Picard does for finding `id` ("MD1").
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

/// The number of findings in NOTES.md: MD1 to MD<FINDINGS>.
pub const FINDINGS: u32 = 7;

fn valid_id(id: &str) -> bool {
    numbered_id(id, "MD", FINDINGS)
}

/// Takes `-java_compat` / `--java_compat` and `-java_compat_except <ids>` out of the arguments and sets
/// the mode.
pub fn take_options(argv: &[String]) -> Result<Vec<String>> {
    let argv: Vec<String> = argv
        .iter()
        .map(|a| match a.as_str() {
            "--java_compat" => "-java_compat".to_string(),
            "--java_compat_except" => "-java_compat_except".to_string(),
            _ => a.clone(),
        })
        .collect();
    hmf_common::compat::take_options(&argv, "MARKDUPS", valid_id)
}

/// For tests: runs `f` with Picard's behaviour on (`true`) or off for every finding, on this thread only.
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
        assert!(valid_id("MD1") && valid_id(&format!("MD{FINDINGS}")));
        assert!(!valid_id(&format!("MD{}", FINDINGS + 1)) && !valid_id("MD0") && !valid_id("M1"));
    }

    #[test]
    fn default_is_fixed_and_scope_overrides() {
        assert!(!java("MD1"));
        assert!(with_java(true, || java("MD1")));
        assert!(!with_java(true, || with_java(false, || java("MD1"))));
    }
}
