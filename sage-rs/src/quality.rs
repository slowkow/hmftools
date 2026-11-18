//! Quality models: BQR tables, MSI jitter parameters, and the per-variant quality cache.
//! The per-read calculation (QualityCalculator.calculateQualityScores) lives in evidence.rs.

use crate::config::Config;
use crate::constants::*;
use crate::evidence::{JitterData, ReadContextCounter};
use crate::read_context::VariantReadContext;
use crate::refdata::{read_lines, RefData};
use crate::refseq::RefSequence;
use crate::repeat::{find_max_repeat, RepeatInfo};
use crate::stats::binomial_cdf;
use crate::types::*;
use anyhow::{bail, Context, Result};
use std::collections::HashMap;

pub struct QualityScores {
    pub seq_tech_base_quality: f64,
    pub recalibrated_base_quality: f64,
    pub final_map_quality: i32,
    pub final_base_quality: f64,
    pub combined_quality: f64,
    pub is_medium_qual: bool,
}

impl QualityScores {
    pub fn invalid() -> Self {
        QualityScores {
            seq_tech_base_quality: -1.0,
            recalibrated_base_quality: -1.0,
            final_map_quality: 0,
            final_base_quality: -1.0,
            combined_quality: 0.0,
            is_medium_qual: false,
        }
    }
    pub fn valid(&self) -> bool {
        self.seq_tech_base_quality != -1.0
    }
}

/// BqrRecordMap: (ref, alt, trinucleotide, original qual, read type) -> recalibrated qual
#[derive(Default)]
pub struct BqrMap {
    map: HashMap<(u8, u8, [u8; 3], u8, u8), f64>,
    enabled: bool,
}

impl BqrMap {
    pub fn load(path: &str) -> Result<BqrMap> {
        let lines = read_lines(path)?;
        let idx: HashMap<&str, usize> = lines[0].split('\t').enumerate().map(|(i, n)| (n, i)).collect();
        let (ai, ri, ti, ci, qi, rqi) = (
            idx["Alt"],
            idx["Ref"],
            idx["TrinucleotideContext"],
            idx.get("ReadType").copied(),
            idx["OriginalQual"],
            idx["RecalibratedQual"],
        );
        let mut map = HashMap::new();
        for line in &lines[1..] {
            let v: Vec<&str> = line.split('\t').collect();
            let tnc = v[ti].as_bytes();
            let read_type = match ci.map(|i| v[i]) {
                Some("SINGLE") => 1,
                Some("DUAL") => 2,
                _ => 0,
            };
            let key = (v[ri].as_bytes()[0], v[ai].as_bytes()[0], [tnc[0], tnc[1], tnc[2]], v[qi].parse::<u8>()?, read_type);
            // the first of any duplicate keys is kept
            map.entry(key).or_insert(v[rqi].parse::<f64>()?);
        }
        Ok(BqrMap { map, enabled: true })
    }

    fn quality_adjustment(&self, ref_: u8, alt: u8, tnc: [u8; 3], qual: u8, read_type: ConsensusType) -> f64 {
        self.map.get(&(ref_, alt, tnc, qual, read_type as u8)).copied().unwrap_or(qual as f64)
    }
}

#[derive(Clone, Debug)]
pub struct JitterParams {
    pub repeat_unit: String,
    pub consensus: String,
    pub scale4: f64,
    pub scale5: f64,
    pub scale6: f64,
    pub gradient: f64,
    pub intercept: f64,
    pub skew: f64,
}

const MAX_SPECIFIC_LENGTH_UNIT: usize = 2;

impl JitterParams {
    fn repeat_unit_length(&self) -> usize {
        if self.repeat_unit.contains('/') {
            self.repeat_unit.split('/').next().unwrap().len()
        } else {
            MAX_SPECIFIC_LENGTH_UNIT + 1
        }
    }
    fn scale_param(&self, repeat_count: i32) -> Option<f64> {
        match repeat_count {
            4 => Some(self.scale4),
            5 => Some(self.scale5),
            6 => Some(self.scale6),
            _ => None,
        }
    }
    /// MsiModelParams.calcErrorRate
    pub(crate) fn calc_error_rate(&self, repeat_count: i32, repeat_count_change: i32, fixed_scale: Option<f64>) -> f64 {
        let scale = fixed_scale.unwrap_or(self.intercept + self.gradient * repeat_count as f64);
        let modified_skew = 1.0 + (self.skew - 1.0) * jmin(scale, 1.0);
        let laplace = |change: i32| -> f64 {
            let sign = if change > 0 { 1 } else { -1 };
            1.0 / (scale * (modified_skew + 1.0 / modified_skew))
                * (((-change * sign) as f64 * modified_skew.powf(sign as f64)) / scale).exp()
        };
        let mut pdf_sum = 0.0;
        for i in -MSI_JITTER_MAX_REPEAT_CHANGE..=MSI_JITTER_MAX_REPEAT_CHANGE {
            pdf_sum += laplace(i);
        }
        if pdf_sum == 0.0 {
            return 0.0;
        }
        jmax(laplace(repeat_count_change) / pdf_sum, MSI_JITTER_DEFAULT_ERROR_RATE)
    }
}

fn default_illumina_params(high_depth: bool) -> Vec<JitterParams> {
    let mut out = Vec::new();
    let mut add = |units: &[&str], o4: f64, o5: f64, o6: f64, gradient: f64| {
        let intercept = o6 - 6.0 * gradient;
        for u in units {
            out.push(JitterParams {
                repeat_unit: u.to_string(),
                consensus: "NONE".into(),
                scale4: o4,
                scale5: o5,
                scale6: o6,
                gradient,
                intercept,
                skew: 1.0,
            });
        }
    };
    if high_depth {
        add(&["A/T", "C/G"], 0.1, 0.15, 0.2, 0.06);
        add(&["AT/TA", "AC/CA/GT/TG", "AG/GA/CT/TC", "CG/GC"], 0.15, 0.2, 0.25, 0.06);
        add(&["3-5bp repeat"], 0.2, 0.25, 0.3, 0.06);
    } else {
        add(&["A/T", "C/G"], 0.105, 0.115, 0.125, 0.04);
        add(&["AT/TA", "AC/CA/GT/TG", "AG/GA/CT/TC", "CG/GC"], 0.11, 0.15, 0.19, 0.04);
        add(&["3-5bp repeat"], 0.14, 0.18, 0.22, 0.02);
    }
    out
}

fn find_applicable_params<'a>(all: &'a [JitterParams], repeat_bases: &[u8]) -> Option<&'a JitterParams> {
    let bases = std::str::from_utf8(repeat_bases).unwrap();
    for p in all {
        if repeat_bases.len() <= MAX_SPECIFIC_LENGTH_UNIT {
            if p.repeat_unit_length() == repeat_bases.len() && p.repeat_unit.contains(bases) {
                return Some(p);
            }
        } else if p.repeat_unit_length() > MAX_SPECIFIC_LENGTH_UNIT {
            return Some(p);
        }
    }
    None
}

/// MsiJitterCalcs.getImpliedAltChange
pub fn implied_alt_change(is_delete: bool, alt_bases: &[u8], repeat: &RepeatInfo) -> i32 {
    let change = alt_bases.len() as i32 / repeat.repeat_length();
    if is_delete {
        -change
    } else {
        change
    }
}

fn is_repeat_of(alt_bases: &[u8], unit: &[u8]) -> bool {
    let n = alt_bases.len() / unit.len();
    alt_bases.len() == n * unit.len() && alt_bases.chunks(unit.len()).all(|c| c == unit)
}

/// MsiJitterCalcs for one sample.
pub struct MsiJitter {
    pub params: Vec<JitterParams>,
    pub probable_msi: bool,
}

impl MsiJitter {
    pub fn load(dir: Option<&str>, sample: &str, high_depth: bool) -> Result<MsiJitter> {
        let Some(dir) = dir else {
            return Ok(MsiJitter { params: default_illumina_params(high_depth), probable_msi: false });
        };
        let param_file = format!("{dir}/{sample}.redux.jitter_params.tsv");
        let count_file = format!("{dir}/{sample}.redux.ms_table.tsv.gz");
        let prediction_file = format!("{dir}/{sample}.redux.msi_prediction.tsv");
        if std::path::Path::new(&prediction_file).exists() {
            let lines = read_lines(&prediction_file)?;
            let idx = lines[0].split('\t').position(|c| c == "PredictedValue");
            if let (Some(i), Some(row)) = (idx, lines.get(1)) {
                let v: f64 = row.split('\t').nth(i).unwrap_or("0").parse().unwrap_or(0.0);
                if v >= 4.0 - 1e-10 {
                    bail!("sample {sample} is predicted MSI; MSI-specific jitter parameters are not implemented");
                }
            }
        }
        // selectConsensusType: the type with the most reads in the jitter counts table, ties to NONE
        let lines = read_lines(&count_file).with_context(|| format!("jitter counts file {count_file}"))?;
        let idx: HashMap<&str, usize> = lines[0].split('\t').enumerate().map(|(i, n)| (n, i)).collect();
        let mut totals: Vec<(String, i64)> = Vec::new();
        for line in &lines[1..] {
            let v: Vec<&str> = line.split('\t').collect();
            let t = v[idx["consensusType"]].to_string();
            let c: i64 = v[idx["readCount"]].parse()?;
            match totals.iter_mut().find(|(k, _)| *k == t) {
                Some(e) => e.1 += c,
                None => totals.push((t, c)),
            }
        }
        // EnumMap iterates in enum order: NONE, SINGLE, DUAL
        totals.sort_by_key(|(k, _)| match k.as_str() {
            "NONE" => 0,
            "SINGLE" => 1,
            _ => 2,
        });
        let consensus = if totals.len() == 1 {
            totals[0].0.clone()
        } else {
            let mut best = ("NONE".to_string(), -1i64);
            for (k, c) in &totals {
                if *c > best.1 || (*c == best.1 && k == "NONE") {
                    best = (k.clone(), *c);
                }
            }
            best.0
        };

        let lines = read_lines(&param_file).with_context(|| format!("jitter params file {param_file}"))?;
        let idx: HashMap<&str, usize> = lines[0].split('\t').enumerate().map(|(i, n)| (n, i)).collect();
        let mut params = Vec::new();
        for line in &lines[1..] {
            let v: Vec<&str> = line.split('\t').collect();
            let c = idx.get("consensusType").map(|&i| v[i]).unwrap_or("NONE");
            if c != consensus {
                continue;
            }
            params.push(JitterParams {
                repeat_unit: v[idx["unit"]].to_string(),
                consensus: c.to_string(),
                scale4: v[idx["optimalScaleRepeat4"]].parse()?,
                scale5: v[idx["optimalScaleRepeat5"]].parse()?,
                scale6: v[idx["optimalScaleRepeat6"]].parse()?,
                gradient: v[idx["scaleFitGradient"]].parse()?,
                intercept: v[idx["scaleFitIntercept"]].parse()?,
                skew: v[idx["microsatelliteSkew"]].parse()?,
            });
        }
        Ok(MsiJitter { params, probable_msi: false })
    }

    /// MsiJitterCalcs.findRepeat
    pub fn find_repeat(ctx: &VariantReadContext) -> Option<RepeatInfo> {
        let variant = &ctx.variant;
        if ctx.marked_invalid || !variant.is_indel() {
            return None;
        }
        let ref_start = ctx.variant_ref_index() + 1;
        let read_start = ctx.var_index + 1;
        let ref_repeat =
            find_max_repeat(&ctx.ref_bases, ref_start, ref_start, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1, false, ref_start);
        let inferred =
            find_max_repeat(&ctx.read_bases, read_start, read_start, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1, false, read_start);
        let alt_bases: &[u8] = if variant.is_delete() { &variant.ref_[1..] } else { &variant.alt[1..] };
        let repeat_to_use = match inferred {
            Some(inf) if is_repeat_of(alt_bases, &inf.bases) => {
                let ref_count = ref_repeat.as_ref().map_or(0, |r| r.count);
                let inferred_count = inf.count - implied_alt_change(variant.is_delete(), alt_bases, &inf);
                Some(RepeatInfo::new(inf.index, inf.bases, ref_count.max(inferred_count)))
            }
            _ => ref_repeat,
        };
        match repeat_to_use {
            Some(r) if !is_repeat_of(alt_bases, &r.bases) => None,
            other => other,
        }
    }

    /// MsiJitterCalcs.calcErrorRate(variant, sampleId, repeat)
    pub fn calc_error_rate(&self, variant: &SimpleVariant, repeat: &RepeatInfo) -> f64 {
        let alt_bases: &[u8] = if variant.is_insert() { &variant.alt[1..] } else { &variant.ref_[1..] };
        let change = implied_alt_change(variant.is_delete(), alt_bases, repeat);
        if change > MSI_JITTER_MAX_REPEAT_CHANGE || change == 0 {
            return 0.0;
        }
        let Some(p) = find_applicable_params(&self.params, &repeat.bases) else { return 0.0 };
        p.calc_error_rate(repeat.count, change, p.scale_param(repeat.count))
    }

    /// MsiJitterCalcs.getErrorRate
    pub fn get_error_rate(&self, repeat_bases: &[u8], repeat_count: i32, repeat_change: i32) -> f64 {
        if repeat_count < MIN_REPEAT_COUNT {
            return MSI_JITTER_DEFAULT_ERROR_RATE;
        }
        let Some(p) = find_applicable_params(&self.params, repeat_bases) else { return MSI_JITTER_DEFAULT_ERROR_RATE };
        p.calc_error_rate(repeat_count, repeat_change, p.scale_param(repeat_count))
    }
}

pub struct QualityCalculator<'a> {
    pub cfg: &'a Config,
    pub bqr: &'a BqrMap,
    pub ref_seq: &'a RefSequence,
    pub msi: &'a MsiJitter,
    pub ref_data: &'a RefData,
    pub chr: &'a str,
}

impl<'a> QualityCalculator<'a> {
    pub fn is_highly_polymorphic(&self, pos: i32) -> bool {
        self.ref_data.hla_contains(self.chr, pos)
    }

    /// lookupRecalibrateQuality
    fn lookup_recalibrate_quality(&self, tnc: [u8; 3], alt: u8, raw: u8, read_type: ConsensusType) -> f64 {
        if raw <= BASE_QUAL_MINIMUM {
            return 0.0;
        }
        if !self.bqr.enabled {
            return raw as f64;
        }
        self.bqr.quality_adjustment(tnc[1], alt, tnc, raw, read_type)
    }
}

fn swap_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        o => o,
    }
}

/// ReadContextQualCache
pub struct QualCache {
    variant_pos: i32,
    variant_alt: Vec<u8>,
    pub msi_indel_repeat: Option<RepeatInfo>,
    /// -1 when not applicable
    pub msi_indel_error_qual: i8,
    pub is_msi_sample_and_variant: bool,
}

impl QualCache {
    pub fn new(ctx: &VariantReadContext, qc: &QualityCalculator) -> Self {
        let repeat = MsiJitter::find_repeat(ctx);
        let mut qual: i8 = -1;
        if let Some(r) = &repeat {
            let rate = qc.msi.calc_error_rate(&ctx.variant, r);
            if rate > 0.0 {
                qual = jround(-10.0 * rate.log10()) as i8;
            }
        }
        QualCache {
            variant_pos: ctx.variant.pos,
            variant_alt: ctx.variant.alt.clone(),
            is_msi_sample_and_variant: repeat.is_some() && qual != -1 && qc.msi.probable_msi,
            msi_indel_repeat: repeat,
            msi_indel_error_qual: qual,
        }
    }
    pub fn uses_msi_indel_error_qual(&self) -> bool {
        self.msi_indel_error_qual != -1
    }
    /// getQual (the Java per-variant memoisation is omitted: the lookup is already a hash probe)
    pub fn get_qual(&self, qc: &QualityCalculator, base_qual: u8, read_type: ConsensusType, ref_index: usize, pos_orientation: bool) -> f64 {
        let pos = self.variant_pos + ref_index as i32;
        let i = qc.ref_seq.index(pos) as usize;
        let b = &qc.ref_seq.bases;
        let mut tnc = [b[i - 1], b[i], b[i + 1]];
        let mut alt = self.variant_alt[ref_index];
        if !pos_orientation {
            tnc = [swap_base(tnc[2]), swap_base(tnc[1]), swap_base(tnc[0])];
            alt = swap_base(alt);
        }
        qc.lookup_recalibrate_quality(tnc, alt, base_qual, read_type)
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum NoiseOutcome {
    Noise,
    ShortenedNoise,
    LengthenedNoise,
    BothNoise,
    HardFilterVariant,
    FilterVariant,
}

fn min_jitter_ratio(tier: Tier, jitter_ratio: f64) -> f64 {
    if tier == Tier::Hotspot && jitter_ratio > MSI_JITTER_RATIO_HOTSPOT_VAF {
        MSI_JITTER_MIN_RATIO_HOTSPOT
    } else {
        MSI_JITTER_MIN_RATIO
    }
}

fn calc_noise(tier: Tier, full_support: i32, jitter_count: i32, error_rate: f64, trinuc: bool) -> NoiseOutcome {
    let rate = if trinuc { jmax(MSI_JITTER_MIN_TRINUC_ERROR_RATE, error_rate) } else { error_rate };
    let prob = 1.0 - binomial_cdf(full_support + jitter_count, rate, full_support.min(jitter_count) - 1);
    if prob > MSI_JITTER_HARD_FILTER_NOISE_RATE {
        return NoiseOutcome::HardFilterVariant;
    } else if prob > MSI_JITTER_NOISE_RATE {
        return NoiseOutcome::FilterVariant;
    }
    let ratio = full_support as f64 / (full_support + jitter_count) as f64;
    if ratio < min_jitter_ratio(tier, ratio) * rate {
        return NoiseOutcome::FilterVariant;
    }
    NoiseOutcome::Noise
}

/// JitterData.setJitterQualFilterState
pub fn set_jitter_qual_filter_state(msi: &MsiJitter, rc: &mut ReadContextCounter) {
    let ctx = rc.ctx.clone();
    let Some(max_repeat) = &ctx.max_repeat else { return };
    let full_support = rc.jitter.valid_qual_full_support;
    let (shortened, lengthened) = (rc.jitter.shortened, rc.jitter.lengthened);
    let tier = rc.tier;
    let is_panel = matches!(tier, Tier::Panel | Tier::Hotspot);
    let has_trinuc = ctx.all_repeats.iter().any(|r| r.repeat_length() == 3);
    // isNonTrinucIndel
    let non_trinuc_indel = ctx.variant.is_indel() && {
        let indel_bases: &[u8] = if ctx.variant.is_insert() { &ctx.variant.alt[1..] } else { &ctx.variant.ref_[1..] };
        ctx.all_repeats.iter().filter(|r| r.repeat_length() != 3).any(|r| {
            let n = indel_bases.len() / r.bases.len();
            indel_bases.len() == n * r.bases.len() && indel_bases.chunks(r.bases.len()).all(|c| c == &r.bases[..])
        })
    };
    let trinuc = is_panel && has_trinuc && !non_trinuc_indel;

    let count = max_repeat.count;
    let shortened_rate = msi.get_error_rate(&max_repeat.bases, count - 1, 1);
    let lengthened_rate = msi.get_error_rate(&max_repeat.bases, count + 1, -1);

    let outcome = (|| -> Option<NoiseOutcome> {
        let mut outcome = None;
        if shortened > full_support {
            let o = calc_noise(tier, full_support, shortened, shortened_rate, trinuc);
            if o != NoiseOutcome::Noise {
                return Some(o);
            }
        } else if full_support > shortened {
            let rate = msi.get_error_rate(&max_repeat.bases, count, -1);
            if calc_noise(tier, shortened, full_support, rate, trinuc) != NoiseOutcome::Noise {
                outcome = Some(NoiseOutcome::ShortenedNoise);
            }
        }
        if lengthened > full_support {
            let o = calc_noise(tier, full_support, lengthened, lengthened_rate, trinuc);
            if o != NoiseOutcome::Noise {
                return Some(o);
            }
        } else if full_support > lengthened {
            let rate = msi.get_error_rate(&max_repeat.bases, count, 1);
            if calc_noise(tier, lengthened, full_support, rate, trinuc) != NoiseOutcome::Noise {
                outcome = Some(if outcome == Some(NoiseOutcome::ShortenedNoise) {
                    NoiseOutcome::BothNoise
                } else {
                    NoiseOutcome::LengthenedNoise
                });
            }
        }
        if shortened.min(lengthened) >= full_support {
            let total = full_support + shortened + lengthened;
            let ratio = full_support as f64 / total as f64;
            let avg_rate = (shortened_rate + lengthened_rate) * 0.5;
            if ratio < min_jitter_ratio(tier, ratio) * avg_rate {
                return Some(NoiseOutcome::FilterVariant);
            }
            let prob = 1.0 - binomial_cdf(total, avg_rate, full_support - 1);
            if prob > MSI_JITTER_HARD_FILTER_NOISE_RATE {
                return Some(NoiseOutcome::HardFilterVariant);
            } else if prob > MSI_JITTER_NOISE_RATE {
                return Some(NoiseOutcome::FilterVariant);
            }
        }
        outcome
    })();

    let j: &mut JitterData = &mut rc.jitter;
    match outcome {
        None => {}
        Some(NoiseOutcome::HardFilterVariant) => {
            j.filter_on_noise = true;
            j.hard_filter_on_noise = true;
        }
        Some(NoiseOutcome::FilterVariant) => j.filter_on_noise = true,
        Some(o) => {
            j.filter_on_noise = false;
            j.hard_filter_on_noise = false;
            j.qual_boost = 1.0;
            if matches!(o, NoiseOutcome::ShortenedNoise | NoiseOutcome::BothNoise) {
                j.qual_boost += shortened as f64 / full_support as f64;
            }
            if matches!(o, NoiseOutcome::LengthenedNoise | NoiseOutcome::BothNoise) {
                j.qual_boost += lengthened as f64 / full_support as f64;
            }
            j.qual_boost = jmin(j.qual_boost, JITTER_QUAL_BOOST_MAX_PERC);
        }
    }
}
