//! The Java unit tests of PEACH, ported (one module per Java test class, reading the Java tests' resource
//! files copied to testdata/), and tests of the port.

mod java_algo;
mod java_io;
mod port;

pub fn testdata(path: &str) -> String {
    format!("{}/testdata/{path}", env!("CARGO_MANIFEST_DIR"))
}
