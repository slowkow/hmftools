//! PurpleSummaryData: the QC checks (purple.qc) and the purity context (purple.purity.tsv), with
//! WholeGenomeDuplication and PurpleQCStatus.calcStatus.

use anyhow::Result;

use hmf_common::fmt::decimal_format_min;

use crate::copynumber::{calc_loh_data, chromosome_map, PurpleCopyNumber};
use crate::drivers::{deleted_genes, QcStatus};
use crate::fitting::{polyclonal_proportion, FittedPurity, FittedPurityScore};
use crate::gene::GeneCopyNumber;

use crate::purity_fitter::BestFit;
use crate::types::*;

const MAX_UNSUPPORTED_SEGMENTS: i32 = 220;
pub const MAX_DELETED_GENES: i32 = 280;
pub const TARGET_REGIONS_MAX_DELETED_GENES: i32 = 1500;
const MIN_PURITY: f64 = 0.2;
const CONTAMINATION_FAIL_LEVEL: f64 = 0.1;
const CONTAMINATION_WARN_LEVEL: f64 = 0.02;
const TINC_FAIL_LEVEL: f64 = 0.25;
const TINC_WARN_LEVEL: f64 = 0.02;

#[derive(Clone, Debug)]
pub struct PurpleQc {
    /// in enum order (ImmutablePurpleQC builds a Sets.immutableEnumSet)
    pub status: Vec<QcStatus>,
    pub method: FittedPurityMethod,
    pub copy_number_segments: i32,
    pub unsupported_copy_number_segments: i32,
    pub deleted_genes: i32,
    pub purity: f64,
    pub contamination: f64,
    pub amber_mean_depth: i32,
    pub cobalt_gender: Gender,
    pub amber_gender: Gender,
    /// in enum order
    pub germline_aberrations: Vec<GermlineAberration>,
    pub loh_percent: f64,
    pub tinc_level: f64,
}

/// PurpleQCStatus.calcStatus
#[allow(clippy::too_many_arguments)]
pub fn calc_status(
    gender_pass: bool,
    unsupported: i32,
    deleted: i32,
    purity: f64,
    method: FittedPurityMethod,
    contamination: f64,
    max_deleted_genes: i32,
    tinc_level: f64,
) -> Vec<QcStatus> {
    let mut r = Vec::new();
    if !gender_pass {
        r.push(QcStatus::WarnGenderMismatch);
    }
    if unsupported > MAX_UNSUPPORTED_SEGMENTS {
        r.push(QcStatus::WarnHighCopyNumberNoise);
    }
    if deleted > max_deleted_genes {
        r.push(QcStatus::WarnDeletedGenes);
    }
    if less_than(purity, MIN_PURITY) && method != FittedPurityMethod::NoTumor {
        r.push(QcStatus::WarnLowPurity);
    }
    if contamination >= CONTAMINATION_FAIL_LEVEL {
        r.push(QcStatus::FailContamination);
    } else if contamination >= CONTAMINATION_WARN_LEVEL {
        r.push(QcStatus::WarnContamination);
    }
    if method == FittedPurityMethod::NoTumor {
        r.push(QcStatus::FailNoTumor);
    }
    if tinc_level >= TINC_FAIL_LEVEL {
        r.push(QcStatus::FailTinc);
    } else if tinc_level > TINC_WARN_LEVEL {
        r.push(QcStatus::WarnTinc);
    }
    if r.is_empty() {
        r.push(QcStatus::Pass);
    }
    r.sort();
    r
}

/// PurpleSummaryData.createQC
#[allow(clippy::too_many_arguments)]
pub fn create_qc(
    contamination: f64,
    best_fit: &BestFit,
    amber_gender: Gender,
    cobalt_gender: Gender,
    copy_numbers: &[PurpleCopyNumber],
    genes: &[GeneCopyNumber],
    aberrations: &[GermlineAberration],
    amber_mean_depth: i32,
    max_deleted_genes: i32,
    tinc_level: f64,
) -> PurpleQc {
    let any_sv = copy_numbers.iter().any(|c| c.sv_support());
    let unsupported = if any_sv { copy_numbers.iter().filter(|c| c.segment_start_support == SegmentSupport::None && c.segment_end_support == SegmentSupport::None).count() as i32 } else { 0 };
    let deleted = deleted_genes(genes);
    let gender_pass = cobalt_gender == amber_gender || aberrations.contains(&GermlineAberration::Klinefelter);
    let status = calc_status(gender_pass, unsupported, deleted, best_fit.fit.purity, best_fit.method, contamination, max_deleted_genes, tinc_level);
    let loh = calc_loh_data(copy_numbers);
    let mut germline_aberrations = aberrations.to_vec();
    germline_aberrations.sort();
    germline_aberrations.dedup();
    PurpleQc {
        status,
        method: best_fit.method,
        copy_number_segments: copy_numbers.len() as i32,
        unsupported_copy_number_segments: unsupported,
        deleted_genes: deleted,
        purity: best_fit.fit.purity,
        contamination,
        amber_mean_depth,
        cobalt_gender,
        amber_gender,
        germline_aberrations,
        loh_percent: loh.loh_percent(),
        tinc_level,
    }
}

pub struct PurityContext {
    pub fit: FittedPurity,
    pub score: FittedPurityScore,
    pub method: FittedPurityMethod,
    pub gender: Gender,
    pub run_mode: RunMode,
    pub polyclonal_proportion: f64,
    pub whole_genome_duplication: bool,
    pub msi_indels_per_mb: f64,
    pub msi_status: MicrosatelliteStatus,
    pub tml: i32,
    pub tml_status: TumorMutationalStatus,
    pub tmb_per_mb: f64,
    pub tmb_status: TumorMutationalStatus,
    pub sv_tmb: i32,
    pub qc: PurpleQc,
    pub targeted: bool,
}

/// WholeGenomeDuplication.wholeGenomeDuplication: at least 11 autosomes with an average major allele
/// copy number (BAF-count weighted) of 1.5 or more
pub fn whole_genome_duplication(copy_numbers: &[PurpleCopyNumber]) -> bool {
    let mut duplicated = 0;
    for (chromosome, list) in chromosome_map(copy_numbers) {
        if !is_autosome(chr(&chromosome)) {
            continue;
        }
        let (mut weighted, mut total) = (0.0, 0i64);
        for c in &list {
            weighted += c.major_allele_copy_number() * c.baf_count as f64;
            total += c.baf_count as i64;
        }
        let average = if total > 0 { weighted / total as f64 } else { 0.0 };
        if greater_or_equal(average, 1.5) {
            duplicated += 1;
        }
    }
    duplicated >= 11
}

/// PurpleSummaryData.createPurity
#[allow(clippy::too_many_arguments)]
pub fn create_purity(
    best_fit: &BestFit,
    gender: Gender,
    run_mode: RunMode,
    tumor_only: bool,
    targeted: bool,
    qc: PurpleQc,
    copy_numbers: &[PurpleCopyNumber],
    somatic: Option<&crate::somatic_stream::SomaticStreamResult>,
    sv_tmb: i32,
) -> PurityContext {
    let wgd = if !tumor_only || best_fit.method == FittedPurityMethod::Normal { whole_genome_duplication(copy_numbers) } else { false };
    PurityContext {
        fit: best_fit.fit.clone(),
        score: best_fit.score.clone(),
        method: best_fit.method,
        gender,
        run_mode,
        polyclonal_proportion: polyclonal_proportion(copy_numbers),
        whole_genome_duplication: wgd,
        msi_indels_per_mb: somatic.map(|s| s.msi_indels_per_mb).unwrap_or(0.0),
        msi_status: somatic.map(|s| s.microsatellite_status()).unwrap_or(MicrosatelliteStatus::Unknown),
        tml: somatic.map(|s| s.tumor_mutational_load()).unwrap_or(0),
        tml_status: somatic.map(|s| s.tml_status()).unwrap_or(TumorMutationalStatus::Unknown),
        tmb_per_mb: somatic.map(|s| s.tmb).unwrap_or(0.0),
        tmb_status: somatic.map(|s| s.tmb_status()).unwrap_or(TumorMutationalStatus::Unknown),
        sv_tmb,
        qc,
        targeted,
    }
}

fn f4(v: f64) -> String {
    decimal_format_min(v, 4, 4)
}

fn join<T: std::fmt::Display>(v: &[T]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

/// PurityContextFile.write: purple.qc then purple.purity.tsv
pub fn write_purity_context(dir: &str, sample: &str, p: &PurityContext) -> Result<()> {
    let q = &p.qc;
    let lines = [
        format!("QCStatus\t{}", join(&q.status)),
        format!("Method\t{}", q.method),
        format!("CopyNumberSegments\t{}", q.copy_number_segments),
        format!("UnsupportedCopyNumberSegments\t{}", q.unsupported_copy_number_segments),
        format!("Purity\t{}", f4(q.purity)),
        format!("AmberGender\t{}", q.amber_gender),
        format!("CobaltGender\t{}", q.cobalt_gender),
        format!("DeletedGenes\t{}", q.deleted_genes),
        format!("Contamination\t{}", f4(q.contamination)),
        format!("GermlineAberrations\t{}", join(&q.germline_aberrations)),
        format!("AmberMeanDepth\t{}", q.amber_mean_depth),
        format!("LohPercent\t{}", f4(q.loh_percent)),
        format!("TincLevel\t{}", f4(q.tinc_level)),
    ];
    let mut text = lines.join("\n");
    text.push('\n');
    // File.separator between the directory (already ending in one) and the file name, as Java does
    std::fs::write(format!("{dir}/{sample}.purple.qc"), text)?;

    let header = "purity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion\tsomaticPenalty\twholeGenomeDuplication\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\tsvTumorMutationalBurden\trunMode\ttargeted";
    let f = &p.fit;
    let s = &p.score;
    let values = [
        f4(f.purity),
        f4(f.norm_factor),
        f4(f.score),
        f4(f.diploid_proportion),
        f4(f.ploidy),
        p.gender.to_string(),
        p.method.to_string(),
        f4(p.polyclonal_proportion),
        f4(s.min_purity),
        f4(s.max_purity),
        f4(s.min_ploidy),
        f4(s.max_ploidy),
        f4(s.min_diploid_proportion),
        f4(s.max_diploid_proportion),
        f4(f.somatic_penalty),
        p.whole_genome_duplication.to_string(),
        f4(p.msi_indels_per_mb),
        p.msi_status.to_string(),
        p.tml.to_string(),
        p.tml_status.to_string(),
        f4(p.tmb_per_mb),
        p.tmb_status.to_string(),
        p.sv_tmb.to_string(),
        p.run_mode.to_string(),
        p.targeted.to_string(),
    ];
    std::fs::write(format!("{dir}/{sample}.purple.purity.tsv"), format!("{header}\n{}\n", values.join("\t")))?;
    Ok(())
}

/// the purity context written when there are no valid regions or no fit (writeEmptyResultFiles)
pub fn empty_purity_context(run_mode: RunMode, targeted: bool) -> PurityContext {
    let zero = FittedPurity { purity: 0.0, norm_factor: 0.0, ploidy: 0.0, score: 0.0, diploid_proportion: 0.0, somatic_penalty: 0.0 };
    let score = FittedPurityScore { min_purity: 0.0, max_purity: 0.0, min_ploidy: 0.0, max_ploidy: 0.0, min_diploid_proportion: 0.0, max_diploid_proportion: 0.0 };
    let qc = PurpleQc {
        status: vec![QcStatus::FailNoTumor],
        method: FittedPurityMethod::NoTumor,
        copy_number_segments: 0,
        unsupported_copy_number_segments: 0,
        deleted_genes: 0,
        purity: 0.0,
        contamination: 0.0,
        amber_mean_depth: 0,
        cobalt_gender: Gender::Female,
        amber_gender: Gender::Female,
        germline_aberrations: vec![GermlineAberration::None],
        loh_percent: 0.0,
        tinc_level: 0.0,
    };
    PurityContext {
        fit: zero,
        score,
        method: FittedPurityMethod::NoTumor,
        gender: Gender::Female,
        run_mode,
        polyclonal_proportion: 0.0,
        whole_genome_duplication: false,
        msi_indels_per_mb: 0.0,
        msi_status: MicrosatelliteStatus::Unknown,
        tml: 0,
        tml_status: TumorMutationalStatus::Unknown,
        tmb_per_mb: 0.0,
        tmb_status: TumorMutationalStatus::Unknown,
        sv_tmb: 0,
        qc,
        targeted,
    }
}
