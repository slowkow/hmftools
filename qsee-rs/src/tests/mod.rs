//! Tests: the Java unit tests of qsee 1.0 (`java_tests`) and tests of the port (`port`).

pub mod java_tests;
pub mod port;

/// helpers shared by the two
pub mod java_tests_support {
    pub use super::java_tests::test_purity_context as purity_context;
}
