//! The data of the report as `plot_qc.R` prepares it: the vis data and cohort percentiles files read as R's
//! read.delim reads them ("NA" and "null" missing), the sample's features of one type merged with the cohort's
//! 5th, 25th, 50th, 75th and 95th percentiles (R's merge(all.x = TRUE, sort = FALSE): the matched rows first, the
//! others after, in the order R's merge puts them), the "key:value;key:value" fields of FeatureName and
//! PlotMetadata split into columns, and every text column a factor with its levels in order of appearance.

use std::collections::HashMap;

use anyhow::{Context, Result};

/// One row of a feature type's plot data.
#[derive(Clone, Debug)]
pub struct Row {
    pub sample_type: String,
    pub value: Option<f64>,
    /// the split fields (FeatureName's, then PlotMetadata's), "" when a row lacks one
    pub fields: HashMap<String, String>,
    /// PctMin, PctLower, PctMid, PctUpper, PctMax (NA without a cohort row)
    pub pct: [Option<f64>; 5],
}

impl Row {
    pub fn field(&self, k: &str) -> &str {
        self.fields.get(k).map(|s| s.as_str()).unwrap_or("")
    }
    pub fn num(&self, k: &str) -> f64 {
        self.field(k).parse().unwrap_or(f64::NAN)
    }
}

/// A feature type's plot data: the rows in R's order and the factor levels of the text columns.
#[derive(Clone, Debug, Default)]
pub struct PlotData {
    pub rows: Vec<Row>,
    pub levels: HashMap<String, Vec<String>>,
}

impl PlotData {
    pub fn levels(&self, k: &str) -> Vec<String> {
        self.levels.get(k).cloned().unwrap_or_default()
    }
}

/// The report's input data.
pub struct Inputs {
    /// (SampleType, FeatureType, FeatureName, value, PlotMetadata) of the tumor's and normal's rows
    sample: Vec<(String, String, String, Option<f64>, String)>,
    cohort: Option<HashMap<(String, String, String), [Option<f64>; 5]>>,
    /// the cohort's keys in file order (R's merge looks up the first match)
    pub has_cohort: bool,
}

fn parse_num(s: &str) -> Option<f64> {
    match s {
        "NA" | "null" | "" => None,
        "Inf" => Some(f64::INFINITY),
        "-Inf" => Some(f64::NEG_INFINITY),
        "NaN" => Some(f64::NAN),
        _ => s.parse().ok(),
    }
}

fn read_text(path: &str) -> Result<String> {
    let mut r = hmf_common::gz::open(path).with_context(|| format!("reading {path}"))?;
    let mut s = String::new();
    std::io::Read::read_to_string(&mut r, &mut s).with_context(|| format!("reading {path}"))?;
    Ok(s)
}

impl Inputs {
    pub fn read(vis_path: &str, cohort_path: Option<&str>, tumor: &str, normal: Option<&str>) -> Result<Inputs> {
        let text = read_text(vis_path)?;
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let col = |n: &str| header.iter().position(|h| *h == n).with_context(|| format!("no column {n} in {vis_path}"));
        let (ci, st, ft, fnm, fv, pm) = (col("SampleId")?, col("SampleType")?, col("FeatureType")?, col("FeatureName")?, col("FeatureValue")?, col("PlotMetadata")?);
        let mut sample = Vec::new();
        for l in lines {
            if l.is_empty() {
                continue;
            }
            let f: Vec<&str> = l.split('\t').collect();
            let g = |i: usize| f.get(i).copied().unwrap_or("");
            // SAMPLE_DATA: the rows of the tumor and the normal (SampleId %in% c(TUMOR_ID, NORMAL_ID))
            if g(ci) != tumor && Some(g(ci)) != normal {
                continue;
            }
            let meta = match g(pm) {
                "NA" | "null" => String::new(),
                m => m.to_string(),
            };
            sample.push((g(st).to_string(), g(ft).to_string(), g(fnm).to_string(), parse_num(g(fv)), meta));
        }
        let cohort = match cohort_path {
            None => None,
            Some(p) => {
                let text = read_text(p)?;
                let mut lines = text.lines();
                let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
                let col = |n: &str| header.iter().position(|h| *h == n).with_context(|| format!("no column {n} in {p}"));
                let (st, ft, fnm) = (col("SampleType")?, col("FeatureType")?, col("FeatureName")?);
                let pcts = [col("Pct_5")?, col("Pct_25")?, col("Pct_50")?, col("Pct_75")?, col("Pct_95")?];
                let mut map = HashMap::new();
                for l in lines {
                    if l.is_empty() {
                        continue;
                    }
                    let f: Vec<&str> = l.split('\t').collect();
                    let g = |i: usize| f.get(i).copied().unwrap_or("");
                    let key = (g(st).to_string(), g(ft).to_string(), g(fnm).to_string());
                    let v = [parse_num(g(pcts[0])), parse_num(g(pcts[1])), parse_num(g(pcts[2])), parse_num(g(pcts[3])), parse_num(g(pcts[4]))];
                    map.entry(key).or_insert(v);
                }
                Some(map)
            }
        };
        let has_cohort = cohort.is_some();
        Ok(Inputs { sample, cohort, has_cohort })
    }

    /// get_plot_data(feature_type): None when the sample has no row of that type.
    pub fn plot_data(&self, feature_type: &str) -> Option<PlotData> {
        let x: Vec<&(String, String, String, Option<f64>, String)> = self.sample.iter().filter(|r| r.1 == feature_type).collect();
        if x.is_empty() {
            return None;
        }
        // the merged rows: (sample row, percentiles)
        let merged: Vec<(usize, [Option<f64>; 5])> = match &self.cohort {
            None => (0..x.len()).map(|i| (i, [None; 5])).collect(),
            Some(cohort) => {
                let keys: Vec<(String, String, String)> = x.iter().map(|r| (r.0.clone(), r.1.clone(), r.2.clone())).collect();
                let matched: Vec<bool> = keys.iter().map(|k| cohort.contains_key(k)).collect();
                // xinds: the index (1-based) of each matched key among the matched keys' first occurrences
                let mut firsts: Vec<&(String, String, String)> = Vec::new();
                let mut xinds = vec![0i64; x.len()];
                for (i, k) in keys.iter().enumerate() {
                    if !matched[i] {
                        continue;
                    }
                    let pos = match firsts.iter().position(|f| *f == k) {
                        Some(p) => p,
                        None => {
                            firsts.push(k);
                            firsts.len() - 1
                        }
                    };
                    xinds[i] = pos as i64 + 1;
                }
                let mut ix: Vec<usize> = (0..x.len()).collect();
                isort_with_index(&mut xinds, &mut ix);
                let n_lone = xinds.iter().take_while(|&&v| v == 0).count();
                let mut out = Vec::new();
                for &i in &ix[n_lone..] {
                    out.push((i, cohort[&keys[i]]));
                }
                for &i in &ix[..n_lone] {
                    out.push((i, [None; 5]));
                }
                out
            }
        };
        // multi-field columns: FeatureName and PlotMetadata, split when every value is "k:v..." (has_multiple_fields)
        let multi = |s: &str| {
            let b = s.as_bytes();
            b.iter().position(|&c| c == b':').is_some_and(|p| p > 0 && p + 1 < b.len())
        };
        let name_split = merged.iter().all(|(i, _)| multi(&x[*i].2));
        let meta_split = merged.iter().all(|(i, _)| multi(&x[*i].4));
        let mut levels: HashMap<String, Vec<String>> = HashMap::new();
        let mut push_level = |k: &str, v: &str| {
            let l = levels.entry(k.to_string()).or_default();
            if !l.iter().any(|x| x == v) {
                l.push(v.to_string());
            }
        };
        let mut rows = Vec::new();
        let split = |s: &str| -> Vec<(String, String)> {
            s.split(';').filter(|p| !p.is_empty()).map(|p| match p.split_once(':') {
                Some((k, v)) => (k.to_string(), v.split(':').next().unwrap_or(v).to_string()),
                None => (p.to_string(), String::new()),
            }).collect()
        };
        // the split columns' keys in order of first appearance (bind_rows), missing values ""
        let mut keys_order: Vec<String> = Vec::new();
        let mut parsed: Vec<Vec<(String, String)>> = Vec::new();
        for (i, _) in &merged {
            let r = x[*i];
            let mut f = Vec::new();
            if name_split {
                f.extend(split(&r.2));
            }
            if meta_split {
                f.extend(split(&r.4));
            }
            for (k, _) in &f {
                if !keys_order.contains(k) {
                    keys_order.push(k.clone());
                }
            }
            parsed.push(f);
        }
        for ((i, pct), f) in merged.iter().zip(parsed) {
            let r = x[*i];
            let mut fields: HashMap<String, String> = keys_order.iter().map(|k| (k.clone(), String::new())).collect();
            for (k, v) in f {
                fields.insert(k, v);
            }
            if !name_split {
                fields.insert("FeatureName".into(), r.2.clone());
            }
            push_level("SampleType", &r.0);
            rows.push(Row { sample_type: r.0.clone(), value: r.3, fields, pct: *pct });
        }
        for r in &rows {
            for k in &keys_order {
                push_level(k, r.field(k));
            }
        }
        Some(PlotData { rows, levels })
    }
}

/// R's isort_with_index (the shell sort of its merge; not stable).
fn isort_with_index(x: &mut [i64], indx: &mut [usize]) {
    let n = x.len();
    let mut h = 1;
    while h <= n / 9 {
        h = 3 * h + 1;
    }
    while h > 0 {
        for i in h..n {
            let v = x[i];
            let iv = indx[i];
            let mut j = i;
            while j >= h && x[j - h] > v {
                x[j] = x[j - h];
                indx[j] = indx[j - h];
                j -= h;
            }
            x[j] = v;
            indx[j] = iv;
        }
        h /= 3;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_merge_order_of_unmatched_rows() {
        // keys matched (1..) and unmatched (0) as R's merge sorts them: the unmatched in shell sort order
        let mut x = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 0, 11, 12, 0, 0];
        let mut ix: Vec<usize> = (0..x.len()).collect();
        isort_with_index(&mut x, &mut ix);
        assert_eq!(&x[..4], &[0, 0, 0, 0]);
        assert_eq!(x[4..].to_vec(), (1..=12).collect::<Vec<i64>>());
    }
}
