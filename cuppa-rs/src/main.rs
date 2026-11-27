//! CUPPA in Rust: a port of hmftools CUPPA 2.5 (see ../cuppa): the Java feature extraction
//! (CuppaDataPrep) and pycuppa's prediction (`python -m cuppa.predict`, with the classifier exported once
//! from its pickle and embedded), so that neither Java nor Python is needed.
//!
//! Usage:
//! - `cuppa-rs com.hartwig.hmftools.cup.prep.CuppaDataPrep <options>` (or `cuppa-rs prep`): the features
//!   file `<sample>.cuppa_data.tsv.gz`, as the Java class writes it;
//! - `cuppa-rs predict <options>` (or `cuppa-rs -m cuppa.predict <options>`): pycuppa's prediction, with
//!   its `--` options; writes `<sample>.cuppa.pred_summ.tsv` and `<sample>.cuppa.vis_data.tsv`;
//! - `cuppa-rs <options>` (or `cuppa-rs com.hartwig.hmftools.cup.cli.PredictionRunner`): the jar's main
//!   class, the feature extraction and then the prediction;
//! - `cuppa-rs plot <vis_data.tsv> <plot.png|plot.pdf>`: the summary plot alone, as `Rscript plot_predictions.R`
//!   draws it from a vis data file (the prediction draws it itself).
//!
//! `-java_compat` (or CUPPA_JAVA_COMPAT=1) reproduces the faults of Java and pycuppa that the port corrects
//! (REPORT.md, CU1-...).

mod compat;
mod plot;
mod predict;
mod prep;
mod sha256;
#[cfg(test)]
mod tests;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let status = std::panic::catch_unwind(|| cuppa_main(&argv)).unwrap_or(1);
    std::process::exit(status);
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Mode {
    Prep,
    Predict,
    Runner,
}

pub fn cuppa_main(argv: &[String]) -> i32 {
    // the JVM options of the bioconda wrapper (cuppa -Xmx12g ...) are accepted and ignored
    let mut argv: Vec<String> =
        argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:") || a.starts_with("-D"))).cloned().collect();
    let mode = match argv.first().map(String::as_str) {
        Some("prep") | Some("com.hartwig.hmftools.cup.prep.CuppaDataPrep") => Mode::Prep,
        Some("predict") | Some("cuppa.predict") => Mode::Predict,
        Some("-m") if argv.get(1).map(String::as_str) == Some("cuppa.predict") => {
            argv.remove(0);
            Mode::Predict
        }
        Some("com.hartwig.hmftools.cup.cli.PredictionRunner") => Mode::Runner,
        Some("plot") => {
            // -java_compat / -java_compat_except (or CUPPA_JAVA_COMPAT*) for the plot's finding (CU9)
            let rest = match compat::take_options(&argv[1..]) {
                Ok(a) => a,
                Err(e) => {
                    hmf_common::error!("{e}");
                    return 1;
                }
            };
            if rest.len() != 2 {
                eprintln!("usage: cuppa-rs plot <vis_data.tsv> <plot.png|plot.pdf> [-java_compat]");
                return 1;
            }
            return match plot::draw(&rest[0], &rest[1]) {
                Ok(()) => 0,
                Err(e) => {
                    hmf_common::error!("{e:#}");
                    1
                }
            };
        }
        Some(a) if !a.starts_with('-') => {
            eprintln!("unknown command {a}: expected prep, predict, a CUPPA class name or the PredictionRunner options");
            return 1;
        }
        _ => {
            argv.insert(0, String::new());
            Mode::Runner
        }
    };
    argv.remove(0);
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            hmf_common::error!("{e}");
            return 1;
        }
    };
    match mode {
        Mode::Predict => predict::main(&argv),
        Mode::Prep => prep::main(&argv),
        Mode::Runner => prep::runner_main(&argv),
    }
}
