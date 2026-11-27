//! The CuppaClassifier's prediction (CuppaClassifier.predict_proba, feat_contrib and the signature
//! quantiles), step for step as pycuppa's estimators compute it:
//!
//! - sub-classifiers: gen_pos (noise profile added, cosine similarity to the class profiles, the non-best
//!   similarity scaling, logistic regression), snv96 (noise profile, log1p, logistic regression), event
//!   (chi2-selected driver/fusion/virus/trait features, log1p and max-scaled SV and TMB features, logistic
//!   regression), gene_exp and alt_sj (selected features max-scaled, cosine similarity, scaling, logistic
//!   regression) for samples with RNA data;
//! - meta-classifiers: dna_combined (logistic regression, calibration, fusion override, sex filter) and
//!   rna_combined (logistic regression, calibration, sex filter);
//! - the combined probabilities (DNA x RNA, floored at 0.01, renormalised) for samples with RNA data.
//!
//! Probabilities are in each logistic regression's class order.

use anyhow::{bail, Result};

use super::features::Prepared;
use super::model::{Calibrator, LogReg, MetaClf, Model, Profiles, RnaClf};
use super::numpy::{self, clip, dot, interp, softmax};
use super::pylog;
use crate::compat;

pub struct SampleResult {
    pub gen_pos: Vec<f64>,
    pub snv96: Vec<f64>,
    pub event: Vec<f64>,
    pub gene_exp: Option<Vec<f64>>,
    pub alt_sj: Option<Vec<f64>>,
    pub dna_combined: Vec<f64>,
    pub rna_combined: Option<Vec<f64>>,
    pub combined: Option<Vec<f64>>,
    /// event feature contributions, [class][feature] in the event logistic regression's order
    pub contribs: Vec<Vec<f64>>,
    /// the signature quantiles, [feature][class of the quantile transformers]
    pub sig_quantiles: Vec<Vec<f64>>,
}

/// `X @ W.T` for the rows of X: each output a dot product as OpenBLAS's dgemv_t computes it (one sample)
/// or as its dgemm kernel does (a cohort, see `numpy::gemm_rows`)
fn matmul(xs: &[&[f64]], w_rows: &[&[f64]]) -> Vec<Vec<f64>> {
    if xs.len() == 1 {
        vec![w_rows.iter().map(|w| dot(w, xs[0])).collect()]
    } else {
        numpy::gemm_rows(xs, w_rows)
    }
}

/// LogisticRegression.predict_proba (multinomial): decision function, softmax
fn lr_proba(xs: &[&[f64]], lr: &LogReg) -> Vec<Vec<f64>> {
    let rows: Vec<&[f64]> = (0..lr.n_classes()).map(|c| lr.coef_row(c)).collect();
    let mut d = if lr.f32 && (xs.len() == 1 || !compat::java("CU7")) {
        xs.iter().map(|x| rows.iter().map(|w| numpy::dot_n(w, x)).collect()).collect()
    } else {
        matmul(xs, &rows)
    };
    for row in d.iter_mut() {
        for (v, b) in row.iter_mut().zip(&lr.intercept) {
            *v += *b;
        }
        softmax(row);
    }
    d
}

/// cos_sim(X, profiles) (`_cos_sim_matrix`): dot products over the product of the L2 norms, NaN as 0. The
/// profiles' norms are given (numpy's, from the model); a row's norm is the square root of its row sum of
/// squares
pub fn cos_sim(xs: &[&[f64]], profiles: &[&[f64]], norms: &[f64]) -> Vec<Vec<f64>> {
    let num = matmul(xs, profiles);
    xs.iter()
        .zip(num)
        .map(|(x, n)| {
            let sq: Vec<f64> = x.iter().map(|v| v * v).collect();
            let a_norm = numpy::row_sum(&sq, xs.len()).sqrt();
            n.iter()
                .zip(norms)
                .map(|(num, b)| {
                    let r = num / (a_norm * b);
                    if r.is_nan() { 0.0 } else { r }
                })
                .collect()
        })
        .collect()
}

/// NonBestSimilarityScaler(exponent).transform of one row: scale = x - max + 1 (the maximum skipping NaN),
/// then (scale ** exponent) * x
pub fn non_best_scale(m: &mut [f64], exponent: i64) {
    let mut max = f64::NAN;
    for v in m.iter() {
        if !v.is_nan() && (max.is_nan() || *v > max) {
            max = *v;
        }
    }
    let shift = -max + 1.0;
    if exponent != 0 {
        for v in m.iter_mut() {
            let s = *v + shift;
            *v = if exponent == 1 { s * *v } else { s.powf(exponent as f64) * *v };
        }
    }
}

/// ProfileSimilarityTransformer.transform followed by NonBestSimilarityScaler(exponent)
fn cos_sim_scaled(xs: &[Vec<f64>], p: &Profiles, exponent: i64) -> Vec<Vec<f64>> {
    let views: Vec<&[f64]> = xs.iter().map(|x| x.as_slice()).collect();
    let rows: Vec<&[f64]> = (0..p.classes.len()).map(|c| p.class_profile(c)).collect();
    let mut out = cos_sim(&views, &rows, &p.norms);
    for m in out.iter_mut() {
        non_best_scale(m, exponent);
    }
    out
}

/// MaxScaler(clip=True).transform for one value
pub fn max_scale(x: f64, max: f64) -> f64 {
    let mut v = x / max;
    if v.is_nan() {
        v = 0.0;
    } else if v == f64::INFINITY {
        v = f64::MAX;
    } else if v == f64::NEG_INFINITY {
        v = f64::MIN;
    }
    clip(v, 0.0, 1.0)
}

/// RollingAvgCalibration.predict_proba(normalize=True) for one row (in the class order of the calibrators)
fn calibrate(row: &[f64], cal: &[Calibrator], n_rows: usize) -> Vec<f64> {
    let mut out: Vec<f64> = row
        .iter()
        .zip(cal)
        .map(|(v, c)| match c {
            Calibrator::Passthrough => *v,
            Calibrator::Isotonic { x, y } => {
                let t = clip(*v, x[0], x[x.len() - 1]);
                interp(t, x, y)
            }
        })
        .collect();
    let s = numpy::row_sum(&out, n_rows);
    for v in out.iter_mut() {
        *v /= s;
    }
    out
}

/// `X.sum(axis=1, min_count=1)` then `(X.T / row_sums).T` for each row
pub fn normalize_rows(rows: &mut [Vec<f64>]) {
    let n = rows.len();
    for row in rows.iter_mut() {
        let s = numpy::row_sum_skipna(row, n);
        for v in row.iter_mut() {
            *v /= s;
        }
    }
}

fn keyword_match(class: &str, keywords: &[String]) -> bool {
    let c = class.to_lowercase();
    keywords.iter().any(|k| c.contains(&k.to_lowercase()))
}

/// SexProbFilter.transform: the probabilities of male classes set to 0 for female samples and of female
/// classes for male samples (classes matched by keyword, case-insensitively), then every row renormalised.
/// `is_male` None leaves a row unchanged (before the renormalisation).
pub fn sex_filter(rows: &mut [Vec<f64>], classes: &[String], male_keywords: &[String], female_keywords: &[String], is_male: &[Option<bool>]) {
    let male: Vec<bool> = classes.iter().map(|c| keyword_match(c, male_keywords)).collect();
    let female: Vec<bool> = classes.iter().map(|c| keyword_match(c, female_keywords)).collect();
    for (row, sex) in rows.iter_mut().zip(is_male) {
        if let Some(m) = sex {
            let male_mask = if *m { 1.0 } else { 0.0 };
            let female_mask = 1.0 - male_mask;
            for (i, v) in row.iter_mut().enumerate() {
                if male[i] {
                    *v *= male_mask;
                }
            }
            for (i, v) in row.iter_mut().enumerate() {
                if female[i] {
                    *v *= female_mask;
                }
            }
        }
    }
    normalize_rows(rows);
}

/// ProbCombiner(combine_mode="multiply").transform: DNA and RNA probabilities of the same classes, each
/// floored at `floor`, multiplied, rows renormalised
pub fn combine(dna: &[Vec<f64>], rna: &[Vec<f64>], floor: f64) -> Vec<Vec<f64>> {
    let mut comb: Vec<Vec<f64>> = dna
        .iter()
        .zip(rna)
        .map(|(d, r)| {
            d.iter()
                .zip(r)
                .map(|(d, r)| {
                    let d = if *d < floor { floor } else { *d };
                    let r = if *r < floor { floor } else { *r };
                    d * r
                })
                .collect()
        })
        .collect();
    normalize_rows(&mut comb);
    comb
}

/// FusionProbOverrider.transform: for samples with a fusion of the overrides (value > 0), the other classes'
/// probabilities multiplied by mask_base_value; then, if any sample had one, every row renormalised
fn fusion_override(rows: &mut [Vec<f64>], classes: &[String], meta: &MetaClf, fusion_values: &[Vec<(String, f64)>]) -> Result<()> {
    let mut fusions: Vec<String> = meta.overrides.iter().map(|o| format!("{}{}", o.feat_prefix, o.feat_basename)).collect();
    fusions.sort();
    fusions.dedup();
    let mut any = false;
    let mut targets_per_sample: Vec<Vec<&str>> = vec![Vec::new(); rows.len()];
    for (s, values) in fusion_values.iter().enumerate() {
        for f in &fusions {
            let Some((_, v)) = values.iter().find(|(n, _)| n == f) else { continue };
            if *v > 0.0 {
                any = true;
                for o in &meta.overrides {
                    if format!("{}{}", o.feat_prefix, o.feat_basename) == *f {
                        targets_per_sample[s].push(&o.target_class);
                    }
                }
            }
        }
    }
    if !any {
        return Ok(());
    }
    // pycuppa renormalises every sample's row once any sample has an override fusion (CU7); corrected, only
    // the rows of samples with one, as if each were predicted alone
    let renormalise_all = crate::compat::java("CU7");
    let mut n_overridden = 0;
    let mut n_with = 0;
    for (row, targets) in rows.iter_mut().zip(&targets_per_sample) {
        if targets.is_empty() {
            continue;
        }
        n_with += 1;
        let mut mask = vec![meta.mask_base_value; classes.len()];
        for t in targets {
            let Some(i) = classes.iter().position(|c| c == t) else {
                bail!("KeyError: \"['{t}'] not in index\" (fusion override target class)");
            };
            mask[i] = 1.0;
        }
        let before = row.clone();
        for (v, m) in row.iter_mut().zip(&mask) {
            *v *= m;
        }
        if before.iter().zip(row.iter()).any(|(a, b)| a != b) {
            n_overridden += 1;
        }
    }
    let _ = (n_with, n_overridden);
    if renormalise_all {
        normalize_rows(rows);
    } else {
        let n = rows.len();
        for (row, targets) in rows.iter_mut().zip(&targets_per_sample) {
            if !targets.is_empty() {
                let s = numpy::row_sum_skipna(row, n);
                for v in row.iter_mut() {
                    *v /= s;
                }
            }
        }
    }
    Ok(())
}

/// QuantileTransformer.transform (uniform output) of one value of one feature, followed by
/// SigCohortQuantileTransformer's unclipped tails
fn sig_quantile(x: f64, q: &[f64], refs: &[f64], neg_q_rev: &[f64], neg_refs_rev: &[f64]) -> f64 {
    let lo = q[0];
    let hi = q[q.len() - 1];
    let mut t = if x.is_nan() { x } else { 0.5 * (interp(x, q, refs) - interp(-x, neg_q_rev, neg_refs_rev)) };
    if x == hi {
        t = 1.0;
    }
    if x == lo {
        t = 0.0;
    }
    if x > hi {
        t = x / hi;
    }
    if x < lo {
        t = -(lo / x);
    }
    t
}

fn rna_sub(xs: &[(usize, &[f64])], clf: &RnaClf) -> Vec<Vec<f64>> {
    // MaxScaledChi2FeatureSelector: the selected features (in the required order already), max-scaled
    let scaled: Vec<Vec<f64>> = xs.iter().map(|(_, x)| x.iter().zip(&clf.max).map(|(v, m)| max_scale(*v, *m)).collect()).collect();
    let sims = cos_sim_scaled(&scaled, &clf.profiles, clf.exponent);
    let views: Vec<&[f64]> = sims.iter().map(|v| v.as_slice()).collect();
    lr_proba(&views, &clf.lr)
}

pub fn predict(m: &Model, x: &Prepared) -> Result<Vec<SampleResult>> {
    let n = x.samples.len();
    pylog::info("PredictionRunner", "get_predictions", &format!("Computing predictions for {n} samples"));

    // gen_pos -------------------------------------------------------------------------------------------
    let gp = &m.gen_pos;
    let noise: Vec<f64> = gp.noise_profile.iter().map(|p| gp.noise_counts as f64 * p).collect();
    let gp_x: Vec<Vec<f64>> = (0..n).map(|s| x.gen_pos(s).iter().zip(&noise).map(|(v, p)| v + p).collect()).collect();
    let gp_sim = cos_sim_scaled(&gp_x, &gp.profiles, gp.exponent);
    let gp_views: Vec<&[f64]> = gp_sim.iter().map(|v| v.as_slice()).collect();
    let gen_pos = lr_proba(&gp_views, &gp.lr);

    // snv96 ---------------------------------------------------------------------------------------------
    let sn = &m.snv96;
    let noise: Vec<f64> = sn.noise_profile.iter().map(|p| sn.noise_counts as f64 * p).collect();
    let sn_x: Vec<Vec<f64>> = (0..n).map(|s| x.snv96(s).iter().zip(&noise).map(|(v, p)| (v + p).ln_1p()).collect()).collect();
    let sn_views: Vec<&[f64]> = sn_x.iter().map(|v| v.as_slice()).collect();
    let snv96 = lr_proba(&sn_views, &sn.lr);

    // event ---------------------------------------------------------------------------------------------
    let ev = &m.event;
    let n_other = ev.other_selected.len();
    let ev_x: Vec<Vec<f64>> = (0..n)
        .map(|s| {
            let e = x.event(s);
            let mut v: Vec<f64> = e[..n_other].to_vec();
            for (j, val) in e[n_other..].iter().enumerate() {
                v.push(max_scale(val.ln_1p(), ev.svtmb_max[j]));
            }
            v
        })
        .collect();
    let ev_views: Vec<&[f64]> = ev_x.iter().map(|v| v.as_slice()).collect();
    let event = lr_proba(&ev_views, &ev.lr);
    let means = ev.lr.x_means.as_ref().unwrap();
    let contribs: Vec<Vec<Vec<f64>>> = ev_x
        .iter()
        .map(|xv| {
            let diff: Vec<f64> = xv.iter().zip(means).map(|(a, b)| a - b).collect();
            (0..ev.lr.n_classes()).map(|c| ev.lr.coef_row(c).iter().zip(&diff).map(|(w, d)| w * d + 0.0).collect()).collect()
        })
        .collect();

    // RNA sub-classifiers, for the samples with RNA data ------------------------------------------------------
    let rna_rows: Vec<usize> = (0..n).filter(|s| x.rna[*s].is_some()).collect();
    let (mut gene_exp, mut alt_sj): (Vec<Option<Vec<f64>>>, Vec<Option<Vec<f64>>>) = (vec![None; n], vec![None; n]);
    if !rna_rows.is_empty() {
        let ge_in: Vec<(usize, &[f64])> = rna_rows.iter().map(|s| (*s, x.gene_exp(*s).unwrap())).collect();
        let sj_in: Vec<(usize, &[f64])> = rna_rows.iter().map(|s| (*s, x.alt_sj(*s).unwrap())).collect();
        for (k, p) in rna_sub(&ge_in, &m.gene_exp).into_iter().enumerate() {
            gene_exp[rna_rows[k]] = Some(p);
        }
        for (k, p) in rna_sub(&sj_in, &m.alt_sj).into_iter().enumerate() {
            alt_sj[rna_rows[k]] = Some(p);
        }
    }

    // sex of each sample (SexProbFilter.sample_sexes): event.trait.is_male as a bool; a missing feature is
    // the fill value -1e-8, which is true (male) for pycuppa; corrected, unknown (CU2)
    let is_male: Vec<Option<bool>> = (0..n)
        .map(|s| {
            let v = x.event(s)[m.req_event.iter().position(|f| f == super::features::SEX_FEATURE_NAME).unwrap_or(0)];
            if x.sex_missing[s] && !compat::java("CU2") { None } else { Some(v != 0.0) }
        })
        .collect();

    // dna_combined -----------------------------------------------------------------------------------------
    let dna = &m.dna_combined;
    let meta_x: Vec<Vec<f64>> = (0..n).map(|s| gen_pos[s].iter().chain(&snv96[s]).chain(&event[s]).copied().collect()).collect();
    let meta_views: Vec<&[f64]> = meta_x.iter().map(|v| v.as_slice()).collect();
    let mut dna_probs: Vec<Vec<f64>> = lr_proba(&meta_views, &dna.lr).iter().map(|r| calibrate(r, &dna.calibrators, n)).collect();
    let fusion_names: Vec<(usize, &String)> = m.req_event.iter().enumerate().filter(|(_, f)| f.starts_with("event.fusion.")).collect();
    let fusion_values: Vec<Vec<(String, f64)>> = (0..n).map(|s| fusion_names.iter().map(|(i, f)| ((*f).clone(), x.event(s)[*i])).collect()).collect();
    fusion_override(&mut dna_probs, &dna.lr.classes, dna, &fusion_values)?;
    sex_filter(&mut dna_probs, &dna.lr.classes, &dna.male_keywords, &dna.female_keywords, &is_male);

    // rna_combined and the combined probabilities --------------------------------------------------------------
    let mut rna_combined: Vec<Option<Vec<f64>>> = vec![None; n];
    let mut combined: Vec<Option<Vec<f64>>> = vec![None; n];
    if !rna_rows.is_empty() {
        let rna = &m.rna_combined;
        let rx: Vec<Vec<f64>> =
            rna_rows.iter().map(|s| gene_exp[*s].as_ref().unwrap().iter().chain(alt_sj[*s].as_ref().unwrap()).copied().collect()).collect();
        let views: Vec<&[f64]> = rx.iter().map(|v| v.as_slice()).collect();
        let mut probs: Vec<Vec<f64>> = lr_proba(&views, &rna.lr).iter().map(|r| calibrate(r, &rna.calibrators, rna_rows.len())).collect();
        let sexes: Vec<Option<bool>> = rna_rows.iter().map(|s| is_male[*s]).collect();
        sex_filter(&mut probs, &rna.lr.classes, &rna.male_keywords, &rna.female_keywords, &sexes);
        let rna_classes = &rna.lr.classes;
        let dna_index: Vec<usize> = rna_classes.iter().map(|c| dna.lr.classes.iter().position(|d| d == c).unwrap()).collect();
        let dna_aligned: Vec<Vec<f64>> = rna_rows.iter().map(|s| dna_index.iter().map(|di| dna_probs[*s][*di]).collect()).collect();
        let comb = combine(&dna_aligned, &probs, m.prob_floor);
        for (k, s) in rna_rows.iter().enumerate() {
            rna_combined[*s] = Some(probs[k].clone());
            combined[*s] = Some(comb[k].clone());
        }
    }

    // signature quantiles ------------------------------------------------------------------------------------
    let sq = &m.sig;
    let prepared_q: Vec<(Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>)> = sq
        .quantiles
        .iter()
        .map(|(nq, q, refs)| {
            let nf = sq.features.len();
            let cols: Vec<Vec<f64>> = (0..nf).map(|f| (0..*nq).map(|i| q[i * nf + f]).collect()).collect();
            let neg_cols: Vec<Vec<f64>> = cols.iter().map(|c| c.iter().rev().map(|v| -v).collect()).collect();
            let neg_refs: Vec<f64> = refs.iter().rev().map(|v| -v).collect();
            (cols, neg_cols, vec![refs.clone()], neg_refs)
        })
        .collect();

    let mut out = Vec::with_capacity(n);
    for s in 0..n {
        let sig_x = x.sig(s);
        let sig_quantiles: Vec<Vec<f64>> = (0..sq.features.len())
            .map(|f| {
                prepared_q
                    .iter()
                    .map(|(cols, neg_cols, refs, neg_refs)| sig_quantile(sig_x[f], &cols[f], &refs[0], &neg_cols[f], neg_refs))
                    .collect()
            })
            .collect();
        out.push(SampleResult {
            gen_pos: gen_pos[s].clone(),
            snv96: snv96[s].clone(),
            event: event[s].clone(),
            gene_exp: gene_exp[s].take(),
            alt_sj: alt_sj[s].take(),
            dna_combined: dna_probs[s].clone(),
            rna_combined: rna_combined[s].take(),
            combined: combined[s].take(),
            contribs: contribs[s].clone(),
            sig_quantiles,
        });
    }
    Ok(out)
}
