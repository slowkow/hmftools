//! BamTools metrics in Rust: a port of hmftools BamTools 1.6.2 BamMetrics (see ../README.md).
//!
//! Usage: `bam-tools-rs [com.hartwig.hmftools.bamtools.metrics.BamMetrics | metrics] <the Java options>`.
//! JVM options (-Xmx..., -XX:..., -D...=...) are dropped, as oncoanalyser passes them to the Java tool's
//! wrapper. `-java_compat` reproduces the Java faults the port corrects (NOTES.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod compat;
mod config;
mod coverage;
mod genes;
mod metrics;
mod partition;
mod record;
mod stats;
mod targets;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};
// the paths the shared code's macros use (crate::log)
#[allow(unused_imports)]
use hmf_common::{fmt, log};

use crate::config::{ConfigBuilder, MetricsConfig, Parsed};

const METRICS_CLASS: &str = "com.hartwig.hmftools.bamtools.metrics.BamMetrics";

fn is_jvm_option(a: &str) -> bool {
    a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-Xss") || a.starts_with("-XX:") || (a.starts_with("-D") && a.contains('='))
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).filter(|a| !is_jvm_option(a)).collect();
    // the command: the jar's main class is BamMetrics, so no command means metrics too
    let rest: Vec<String> = match argv.first().map(|s| s.as_str()) {
        Some(METRICS_CLASS) | Some("metrics") => argv[1..].to_vec(),
        Some(c) if c.starts_with("com.hartwig.") => {
            eprintln!("bam-tools (Rust port): main class {c} is not ported; only {METRICS_CLASS}");
            std::process::exit(2);
        }
        _ => argv,
    };
    let rest = match compat::take_options(&rest) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    let cb = match ConfigBuilder::parse(&rest) {
        Ok(Parsed::Run(cb)) => cb,
        Ok(Parsed::Exit) => return,
        Err(_) => std::process::exit(1),
    };
    // ConfigUtils.setLogLevel
    let level = cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    if cb.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("{} version {}", config::APP_NAME, config::VERSION);
    let config = match MetricsConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            error!("{e:#}");
            std::process::exit(1);
        }
    };
    match metrics::run(&config) {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            error!("{e:#}");
            std::process::exit(1);
        }
    }
}
