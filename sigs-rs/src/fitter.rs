//! SignatureFitter: the counts of each sample (from the VCF or a counts matrix) fitted to the signatures by
//! non-negative least squares, scaled to the sample's total, and written with the misallocated counts.

use std::collections::HashMap;
use std::fmt::Write as _;

use anyhow::{anyhow, bail, Result};
use hmf_common::fmt::{java_fixed, java_general};
use hmf_common::{debug, error, info};

use crate::compat;
use crate::config::{SigsConfig, VariantFilters};
use crate::lsq::LeastSquaresFit;
use crate::matrix::{load_matrix_file, load_matrix_lines, Matrix};
use crate::position_freq::{self, PositionFrequencies};
use crate::snv_loader::{self, bucket_names};

pub const SIG_MISALLOCATED: &str = "MISALLOC";
const ALLOC_ROUND: i32 = 3;
const PERC_ROUND: i32 = 5;

/// The SNV signatures hmf-common carries as a resource (sigs/snv_cosmic_signatures.csv), used without
/// -signatures_file
const DEFAULT_SIGNATURES: &str = hmf_common::sigs::SNV_COSMIC_SIGNATURES;

#[derive(Clone, Debug, PartialEq)]
pub struct SignatureAllocation {
    pub signature: String,
    pub allocation: f64,
    pub percent: f64,
}

/// DataUtils.round(value, roundFactor): Math.round(value / 10^-n) * 10^-n, with Math.pow's 10^-n
/// (StrictMath's, which for 10^-5 is a unit in the last place below the literal 1e-5)
pub fn round(value: f64, round_factor: i32) -> f64 {
    let unit = hmf_common::java_strictmath::pow(10.0, -(round_factor as f64));
    hmf_common::java_math::java_round(value / unit) as f64 * unit
}

/// DataUtils.sizeToStr(size, false), for the debug log
fn size_to_str(size: f64) -> String {
    let log = size.abs().log10();
    if log >= 6.0 {
        format!("{}M", java_fixed(size / 1e6, 1))
    } else if log >= 3.0 {
        format!("{}K", java_fixed(size / 1e3, 1))
    } else {
        java_fixed(size, 0)
    }
}

/// SigUtils.calculateFittedCounts
pub fn calculate_fitted_counts(signatures: &Matrix, allocations: &[f64]) -> Vec<f64> {
    let mut fitted = vec![0.0; signatures.rows];
    for (trans_id, allocation) in allocations.iter().enumerate().take(signatures.cols) {
        for (cat_id, f) in fitted.iter_mut().enumerate() {
            *f += allocation * signatures.get(cat_id, trans_id);
        }
    }
    fitted
}

/// SigResiduals of SigUtils.calcResiduals
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SigResiduals {
    pub total: f64,
    pub percent: f64,
    pub excess: f64,
}

pub fn calc_residuals(counts: &[f64], fitted: &[f64], total_counts: f64) -> Result<SigResiduals> {
    let mut r = SigResiduals::default();
    for (cat_id, count) in counts.iter().enumerate() {
        let f = *fitted.get(cat_id).ok_or_else(|| {
            anyhow!("java.lang.ArrayIndexOutOfBoundsException: the signatures have fewer buckets ({}) than the counts", fitted.len())
        })?;
        let diff = f - count;
        r.total += diff.abs();
        if diff > 0.0 {
            r.excess += diff;
        }
    }
    r.percent = r.total / total_counts;
    Ok(r)
}

/// SignatureAllocationFile.toLines
pub fn allocation_lines(allocations: &[SignatureAllocation]) -> String {
    let mut out = String::from("signature\tallocation\tpercent\n");
    for a in allocations {
        let _ = writeln!(out, "{}\t{}\t{}", a.signature, java_fixed(a.allocation, 3), java_fixed(a.percent, 5));
    }
    out
}

/// SignatureAllocationFile.generateFilename
pub fn allocation_filename(output_prefix: &str, file_id: &str) -> String {
    format!("{output_prefix}{file_id}.sig.allocation.tsv")
}

/// CommonUtils.formOutputFilename
pub fn form_output_filename(output_prefix: &str, output_id: Option<&str>, file_id: &str) -> String {
    match output_id {
        Some(id) => format!("{output_prefix}{file_id}.{id}.csv"),
        None => format!("{output_prefix}{file_id}.csv"),
    }
}

/// The allocations of one sample (processSampleResults): the fit scaled to the sample's total, each
/// signature above both minimums, then MISALLOC (the residuals' total less their excess).
///
/// A fit of zero, where no signature has counts in the sample's buckets, is scaled by total / 0 in Java:
/// every allocation becomes NaN, passes both minimums, and is written as 0.000 (SI6); corrected, nothing is
/// allocated and MISALLOC holds all the counts.
pub fn sample_allocations(
    config: &SigsConfig,
    sample_id: &str,
    names: &[String],
    signatures: &Matrix,
    counts: &[f64],
    sample_total: f64,
    contribs: &[f64],
) -> Result<Vec<SignatureAllocation>> {
    let mut allocs: Vec<f64> = contribs.to_vec();
    let fit_total: f64 = allocs.iter().fold(0.0, |s, v| s + v);
    if fit_total != 0.0 || compat::java("SI6") {
        let adjust = sample_total / fit_total;
        allocs.iter_mut().for_each(|a| *a *= adjust);
    }
    let mut out = Vec::new();
    let mut alloc_total = 0.0;
    for (s, &alloc) in allocs.iter().enumerate() {
        let perc = alloc / sample_total;
        alloc_total += alloc;
        if alloc < config.min_fit || perc < config.min_fit_perc {
            continue;
        }
        out.push(SignatureAllocation {
            signature: names[s].clone(),
            allocation: round(alloc, ALLOC_ROUND),
            percent: round(perc, PERC_ROUND),
        });
    }
    let fitted = calculate_fitted_counts(signatures, &allocs);
    let residuals = calc_residuals(counts, &fitted, sample_total)?;
    let misallocated = residuals.total - residuals.excess;
    debug!(
        "sample({}) alloc({} perc={} of total={}) residuals({} total={} excess={})",
        sample_id,
        size_to_str(alloc_total),
        java_fixed(alloc_total / sample_total, 3),
        size_to_str(sample_total),
        java_fixed(residuals.percent, 3),
        size_to_str(residuals.total),
        size_to_str(residuals.excess)
    );
    out.push(SignatureAllocation {
        signature: SIG_MISALLOCATED.to_string(),
        allocation: round(misallocated, ALLOC_ROUND),
        percent: round(misallocated / sample_total, PERC_ROUND),
    });
    Ok(out)
}

/// The run (SignatureFitter.run); Err where Java exits with status 1
pub fn run(config: &SigsConfig, filters: &VariantFilters) -> Result<()> {
    info!("running sample signature fit");

    // initialise()
    let mut sample_ids = config.sample_ids.clone();
    if sample_ids.is_empty() {
        error!("no sample ID(s) configured");
        // Java returns from run() and exits with status 0 (SI7); corrected, status 1
        return if compat::java("SI7") { Ok(()) } else { Err(anyhow!("no sample ID(s) configured")) };
    }
    if config.vcf_file.is_none() && config.snv_counts_file.is_none() {
        error!("no sample VCF or sample count file configured");
        return if compat::java("SI7") { Ok(()) } else { Err(anyhow!("no sample VCF or sample count file configured")) };
    }

    let mut sig_names: Vec<String> = Vec::new();
    let signatures = match &config.signatures_file {
        None => {
            info!("using default SNV Cosmic signatures");
            let lines = hmf_common::java_text::lines_of(DEFAULT_SIGNATURES);
            load_matrix_lines(&lines, &mut sig_names, Some(&[]), false)?
        }
        Some(f) => load_matrix_file(f, &mut sig_names, None, false)?,
    };

    // loadSampleMatrixCounts: Java passes its sample list as the list the file's column names are added to
    // (SI4); corrected, each configured sample takes the column of its name
    let counts_matrix: Option<Matrix> = match &config.snv_counts_file {
        None => None,
        Some(f) => {
            if compat::java("SI4") {
                Some(load_matrix_file(f, &mut sample_ids, Some(&["BucketName"]), false)?)
            } else {
                let mut columns = Vec::new();
                let m = load_matrix_file(f, &mut columns, Some(&["BucketName"]), false)?;
                let mut selected = Matrix::new(m.rows, sample_ids.len());
                for (i, s) in sample_ids.iter().enumerate() {
                    let Some(c) = columns.iter().position(|x| x == s) else {
                        bail!("sample {s} has no column in the sample counts file {f}");
                    };
                    for r in 0..m.rows {
                        selected.data[r][i] = m.data[r][c];
                    }
                }
                Some(selected)
            }
        }
    };

    let prefix = config.output_prefix();
    let multi = sample_ids.len() > 1;
    let mut fit_lines: Option<String> = if multi { Some(String::from("SampleId,Signature,Allocation,Percent\n")) } else { None };
    let fit_filename = form_output_filename(&prefix, config.output_id.as_deref(), "sig_snv_fit");
    if multi {
        // createBufferedWriter: Java opens the file here and exits if it cannot
        if let Err(e) = std::fs::write(&fit_filename, "") {
            error!("error writing to file({fit_filename}): java.io.FileNotFoundException: {fit_filename} ({e})");
            bail!("cannot write {fit_filename}");
        }
    }

    // performFit()
    let sample_count = sample_ids.len();
    info!("fitting sample({}) with {} signatures", sample_count, signatures.cols);
    let mut lsq = LeastSquaresFit::new(signatures.rows, signatures.cols)?;
    let names = bucket_names();
    let bucket_index: HashMap<String, usize> = names.iter().enumerate().map(|(i, n)| (n.clone(), i)).collect();

    for (i, sample_id) in sample_ids.iter().enumerate() {
        let counts: Vec<f64> = match &counts_matrix {
            Some(m) => m.col(i),
            None => {
                let vcf = config.vcf_file.as_deref().unwrap();
                let variants = snv_loader::load_somatic_variants(vcf, sample_id)?;
                info!("sample({}) processing {} variants", sample_id, variants.len());
                let mut pos_freq = if config.position_bucket_size > 0 { Some(PositionFrequencies::new(config.position_bucket_size)) } else { None };
                let c = snv_loader::count_variants(sample_id, &variants, filters, &bucket_index, pos_freq.as_mut())?;
                if sample_count == 1 {
                    let filename = form_output_filename(&prefix, config.output_id.as_deref(), &format!("{sample_id}.sig.snv_counts"));
                    snv_loader::write_sample_counts(&filename, std::slice::from_ref(sample_id), &names, std::slice::from_ref(&c.counts));
                    if let Some(pf) = &pos_freq {
                        position_freq::write_sample_counts(config.v38, &prefix, config.position_bucket_size, sample_id, pf, config.write_pos_freq_coords)?;
                    }
                }
                c.counts
            }
        };

        let total: f64 = counts.iter().fold(0.0, |s, v| s + v);
        let allocations = if total > 0.0 {
            lsq.initialise(&signatures.data, &counts);
            lsq.solve()?;
            sample_allocations(config, sample_id, &sig_names, &signatures, &counts, total, lsq.contribs())?
        } else {
            Vec::new()
        };

        // writeSigAllocations
        if !multi {
            let file_id = match &config.output_id {
                Some(id) => format!("{sample_id}.{id}"),
                None => sample_id.clone(),
            };
            let filename = allocation_filename(&prefix, &file_id);
            if let Err(e) = std::fs::write(&filename, allocation_lines(&allocations)) {
                error!("error writing to fit output: java.nio.file.NoSuchFileException: {filename} ({e})");
                bail!("cannot write {filename}");
            }
        } else if let Some(lines) = fit_lines.as_mut() {
            for a in &allocations {
                let _ = writeln!(lines, "{},{},{},{}", sample_id, a.signature, java_general(a.allocation, 6), java_general(a.percent, 6));
            }
        }

        if i > 0 && i % 100 == 0 {
            info!("processed {} samples", i);
        }
    }

    info!("sample signature fit complete");
    if let Some(lines) = fit_lines {
        if let Err(e) = std::fs::write(&fit_filename, lines) {
            error!("error writing to fit output: {e}");
            bail!("cannot write {fit_filename}");
        }
    }
    Ok(())
}
