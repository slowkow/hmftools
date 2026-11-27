//! ORANGE's findings (OR1 to OR<FINDINGS>, NOTES.md): `java("ORn")` is true where the port reproduces Java.

/// The number of findings with a switch
pub const FINDINGS: u32 = 16;

pub fn valid_id(id: &str) -> bool {
    hmf_common::compat::numbered_id(id, "OR", FINDINGS)
}

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

/// Whether to behave as Java does for the finding `id` (`-java_compat`, `-java_compat_except`)
#[inline]
pub fn java(id: &str) -> bool {
    #[cfg(test)]
    if let Some(v) = TEST_MODE.with(|m| m.get()) {
        return v;
    }
    hmf_common::compat::java(id)
}

/// For tests: runs `f` with Java behaviour on (`true`) or off for every finding, on this thread
#[cfg(test)]
pub fn with_java<T>(on: bool, f: impl FnOnce() -> T) -> T {
    let prev = TEST_MODE.with(|m| m.replace(Some(on)));
    let r = f();
    TEST_MODE.with(|m| m.set(prev));
    r
}
