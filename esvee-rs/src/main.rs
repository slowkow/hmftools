//! ESVEE in Rust: a port of hmftools ESVEE 2.1 (see ../esvee).
//!
//! Usage: `esvee-rs <prep|assembly|depth|caller> <the Java step's options>`. Each step reads and writes the
//! same files as the Java application of that name (PrepApplication, AssemblyApplication, DepthAnnotator,
//! CallerApplication).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod assembly;
mod caller;
mod common;
mod depth;
mod prep;

pub use common::log;

pub const VERSION: &str = "2.1";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let result = match argv.first().map(String::as_str) {
        Some("prep") => prep::run(&argv[1..]),
        Some("assembly") => assembly::run(&argv[1..]),
        Some("depth") => depth::run(&argv[1..]),
        Some("caller") => caller::run(&argv[1..]),
        _ => {
            eprintln!("usage: esvee-rs <prep|assembly|depth|caller> <options of the Java step>");
            std::process::exit(2);
        }
    };
    if let Err(e) = result {
        error!("{e:#}");
        std::process::exit(1);
    }
}
