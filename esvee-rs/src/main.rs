//! ESVEE in Rust: a port of hmftools ESVEE 2.1 (see ../esvee).
//!
//! Usage: `esvee-rs <prep|assembly|depth|caller> <the Java step's options>`. Each step reads and writes the
//! same files as the Java application of that name (PrepApplication, AssemblyApplication, DepthAnnotator,
//! CallerApplication). Without a step name, `esvee-rs <options>` runs the four steps in turn with the
//! options of all four, as the Java EsveeApplication does (see `esvee.rs`).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod assembly;
mod caller;
mod common;
mod depth;
mod esvee;
mod prep;

pub use common::log;

pub const VERSION: &str = "2.1";

static VERSION_LOGGED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Whether the "Esvee version" line is still to be logged: once per run, also when one call runs every step.
pub fn version_log_due() -> bool {
    !VERSION_LOGGED.swap(true, std::sync::atomic::Ordering::Relaxed)
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let argv = match common::compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e:#}");
            std::process::exit(1);
        }
    };
    let result = match argv.first().map(String::as_str) {
        Some("prep") => prep::run(&argv[1..]),
        Some("assembly") => assembly::run(&argv[1..]),
        Some("depth") => depth::run(&argv[1..]),
        Some("caller") => caller::run(&argv[1..]),
        // no step name: every step, as EsveeApplication
        Some(a) if a.starts_with('-') => esvee::run(&argv),
        _ => {
            eprintln!("usage: esvee-rs <prep|assembly|depth|caller> <options of the Java step> [-java_compat]");
            eprintln!("       esvee-rs <options of all four steps> [-stages PREP;ASSEMBLY;REF_DEPTH;CALLING] [-java_compat]");
            eprintln!("       esvee-rs -version");
            std::process::exit(2);
        }
    };
    if let Err(e) = result {
        error!("{e:#}");
        std::process::exit(1);
    }
}
