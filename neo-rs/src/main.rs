//! NEO in Rust: a port of hmftools NEO 1.4 (neoepitope finding and scoring; see ../../hmftools/neo).
//!
//! Usage, as the bioconda wrapper `neo` takes it: `neo-rs [<main class>] <the Java options>`. With no
//! class, or com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder (the jar's main class), the finder runs;
//! with com.hartwig.hmftools.neo.score.NeoScorer, the scorer. JVM options (-Xmx..., -XX:..., -D...=...)
//! are dropped. `-java_compat` reproduces the Java faults the port corrects (REPORT.md, NE1...).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if, clippy::collapsible_else_if, clippy::if_same_then_else)]

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod bind;
mod compat;
mod epitope;
mod finder;
mod score;
#[cfg(test)]
mod tests;

use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::error;

pub const VERSION: &str = "1.4";
pub const APP_NAME: &str = "Neo";

pub const FINDER_CLASS: &str = "com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder";
pub const SCORER_CLASS: &str = "com.hartwig.hmftools.neo.score.NeoScorer";

fn is_jvm_option(a: &str) -> bool {
    a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-Xss") || a.starts_with("-XX:") || (a.starts_with("-D") && a.contains('='))
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).filter(|a| !is_jvm_option(a) && !a.is_empty()).collect();
    // the bioconda wrapper runs `java -cp neo.jar <class> ...` when the first argument is a class name
    let (scorer, rest) = match argv.first().map(|s| s.as_str()) {
        Some(FINDER_CLASS) => (false, argv[1..].to_vec()),
        Some(SCORER_CLASS) => (true, argv[1..].to_vec()),
        Some(c) if c.starts_with("com.hartwig.") => {
            eprintln!("neo (Rust port): main class {c} is not ported; only {FINDER_CLASS} and {SCORER_CLASS}");
            std::process::exit(2);
        }
        _ => (false, argv),
    };
    let rest = match compat::take_options(&rest) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    if scorer {
        score::register(&mut cb);
    } else {
        finder::register(&mut cb);
    }
    match cb.parse(&rest) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return,
        Err(_) => std::process::exit(1),
    }
    let result = if scorer { score::run(&cb) } else { finder::run(&cb) };
    if let Err(e) = result {
        error!("{e:#}");
        std::process::exit(1);
    }
}
