//! ChordModel.predict and chord_predict.R in Rust: the mutation context counts merged into the model's
//! features (SNVs by substitution, INDELs by context, SVs by type and length), as fractions of each
//! sample's SNVs, INDELs and SVs; the random forest's probabilities of BRCA1- and BRCA2-type HR
//! deficiency; the HR status and HRD type with the QC rules (MSI, few INDELs, few SVs); and the
//! predictions file as R's write.table writes it.

pub mod forest;
pub mod rds;
pub mod rfmt;

use anyhow::{anyhow, bail, Context, Result};

use hmf_common::{debug, error, info};

use crate::compat;
use forest::Forest;

pub const SNV_SUBSTITUTIONS: [&str; 6] = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"];

/// INDEL_NEW_CONTEXTS: the merged INDEL features and the regular expressions of their columns
const INDEL_NEW_CONTEXTS: [(&str, &str); 7] = [
    ("del.rep", "del.rep"),
    ("ins.rep", "ins.rep"),
    ("del.mh.bimh.1", "del.mh.bimh.1"),
    ("del.mh.bimh.2.5", "del.mh.bimh.[2,3,4,5]"),
    ("ins.mh", "ins.mh"),
    ("del.none", "del.none"),
    ("ins.none", "ins.none"),
];

const HRD_PROBABILITY_CUTOFF: f64 = 0.5;
const MSI_MIN_INDEL_REP: f64 = 14000.0;
const LOW_INDEL_LOAD_THRES: f64 = 100.0;
const LOW_SV_LOAD_THRES: f64 = 30.0;

pub const HR_DEFICIENT: &str = "HR_deficient";
pub const HR_PROFICIENT: &str = "HR_proficient";
pub const HR_UNCLEAR: &str = "cannot_be_determined";

/// The mutation contexts table as read.delim reads it: column names and, per sample, its values
#[derive(Clone, Debug)]
pub struct ContextTable {
    pub columns: Vec<String>,
    pub rows: Vec<(String, Vec<f64>)>,
}

impl ContextTable {
    pub fn from_counts(rows: &[(String, Vec<crate::prep::MutContextCount>)]) -> ContextTable {
        let columns = rows.first().map(|(_, c)| c.iter().map(|m| m.name.clone()).collect()).unwrap_or_default();
        ContextTable { columns, rows: rows.iter().map(|(id, c)| (id.clone(), c.iter().map(|m| m.count as f64).collect())).collect() }
    }

    /// A mutation contexts file (predict only): a header, then the sample id and numbers per line
    pub fn read(path: &str) -> Result<ContextTable> {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let header = lines.next().ok_or_else(|| anyhow!("{path}: empty file"))?;
        let mut columns: Vec<String> = header.split('\t').map(str::to_string).collect();
        if columns.is_empty() {
            bail!("{path}: no columns");
        }
        columns.remove(0);
        let mut rows = Vec::new();
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != columns.len() + 1 {
                bail!("{path}: a line with {} fields where the header has {}", f.len(), columns.len() + 1);
            }
            let values = f[1..]
                .iter()
                .map(|v| v.trim().parse::<f64>().map_err(|_| anyhow!("{path}: '{v}' is not a number (R reads such a column as text, which it cannot sum)")))
                .collect::<Result<Vec<f64>>>()?;
            rows.push((f[0].to_string(), values));
        }
        Ok(ContextTable { columns, rows })
    }
}

/// R's regular expressions of the column groups, as far as these patterns need them: literals, '.'
/// for any character, a bracket list, `\w`, '^', '$' and `{n}` after a bracket list
#[derive(Clone, Debug)]
enum Token {
    Lit(char),
    Any,
    Set(Vec<char>),
    Word,
    Range(char, char),
}

fn tokens(pattern: &str) -> (bool, Vec<(Token, usize)>, bool) {
    let p: Vec<char> = pattern.chars().collect();
    let mut i = 0;
    let anchored_start = p.first() == Some(&'^');
    if anchored_start {
        i = 1;
    }
    let mut anchored_end = false;
    let mut out: Vec<(Token, usize)> = Vec::new();
    while i < p.len() {
        let c = p[i];
        let t = match c {
            '$' if i == p.len() - 1 => {
                anchored_end = true;
                i += 1;
                continue;
            }
            '.' => Token::Any,
            '\\' => {
                i += 1;
                match p[i] {
                    'w' => Token::Word,
                    other => Token::Lit(other),
                }
            }
            '[' => {
                let mut set = Vec::new();
                i += 1;
                while p[i] != ']' {
                    set.push(p[i]);
                    i += 1;
                }
                if set.len() == 3 && set[1] == '-' { Token::Range(set[0], set[2]) } else { Token::Set(set) }
            }
            '{' => {
                let mut n = 0usize;
                i += 1;
                while p[i] != '}' {
                    n = n * 10 + p[i].to_digit(10).unwrap() as usize;
                    i += 1;
                }
                out.last_mut().unwrap().1 = n;
                i += 1;
                continue;
            }
            other => Token::Lit(other),
        };
        out.push((t, 1));
        i += 1;
    }
    (anchored_start, out, anchored_end)
}

fn token_matches(t: &Token, c: char) -> bool {
    match t {
        Token::Lit(l) => *l == c,
        Token::Any => true,
        Token::Set(s) => s.contains(&c),
        Token::Word => c.is_alphanumeric() || c == '_',
        Token::Range(a, b) => *a <= c && c <= *b,
    }
}

/// grep(pattern, x): whether the pattern matches somewhere in `s` (the patterns have no alternatives or
/// variable repeats, so a match at a position is a plain sequence check)
pub fn r_grepl(pattern: &str, s: &str) -> bool {
    let (start, toks, end) = tokens(pattern);
    let chars: Vec<char> = s.chars().collect();
    let len: usize = toks.iter().map(|(_, n)| n).sum();
    let matches_at = |p: usize| -> bool {
        if p + len > chars.len() {
            return false;
        }
        let mut k = p;
        for (t, n) in &toks {
            for _ in 0..*n {
                if !token_matches(t, chars[k]) {
                    return false;
                }
                k += 1;
            }
        }
        !end || k == chars.len()
    };
    if start { matches_at(0) } else { (0..=chars.len()).any(matches_at) }
}

/// One sample's prediction, before formatting
#[derive(Clone, Debug, PartialEq)]
pub struct Prediction {
    pub sample: String,
    pub p_brca1: f64,
    pub p_brca2: f64,
    pub p_hrd: f64,
    pub hr_status: String,
    pub hrd_type: String,
    pub remarks_hr_status: String,
    pub remarks_hrd_type: String,
}

/// groupColumns: the column indexes of each pattern, which must account for every column once
fn group_columns(columns: &[String], idx: &[usize], patterns: &[&str], fixed: bool) -> Result<Vec<Vec<usize>>> {
    let groups: Vec<Vec<usize>> = patterns
        .iter()
        .map(|p| idx.iter().copied().filter(|&i| if fixed { columns[i].contains(p) } else { r_grepl(p, &columns[i]) }).collect())
        .collect();
    let new_ncol: usize = groups.iter().map(Vec::len).sum();
    if new_ncol != idx.len() {
        error!("Output no. of columns ({new_ncol}) != original no. of columns ({})", idx.len());
        bail!("R: Output no. of columns ({new_ncol}) != original no. of columns ({})", idx.len());
    }
    Ok(groups)
}

/// x / rowSums(x), with NaN (no mutations of the type: 0 / 0) as 0
fn relative(values: &[f64]) -> Vec<f64> {
    let total: f64 = values.iter().sum();
    values.iter().map(|v| v / total).map(|r| if r.is_nan() { 0.0 } else { r }).collect()
}

/// chordPredict on the table: one prediction per row, in row order
pub fn predict(table: &ContextTable, forest: &Forest) -> Result<Vec<Prediction>> {
    // the row names as read.delim converts them (CH3)
    let ids: Vec<String> = table.rows.iter().map(|(id, _)| id.clone()).collect();
    let samples = if compat::java("CH3") {
        rfmt::r_row_names(&ids).map_err(|e| {
            error!("R: Error in read.table: {e}");
            anyhow!("R execution failed: {e}")
        })?
    } else {
        ids
    };

    info!("Transforming mutation contexts");
    let all: Vec<usize> = (0..table.columns.len()).collect();
    let by_type = group_columns(&table.columns, &all, &[r"^\w\[\w>\w\]\w$", "^[a-z]{3}[.]", "^[A-Z]{3}"], false)?;
    let snv_groups = group_columns(&table.columns, &by_type[0], &SNV_SUBSTITUTIONS, true)?;
    let indel_patterns: Vec<&str> = INDEL_NEW_CONTEXTS.iter().map(|(_, p)| *p).collect();
    let indel_groups = group_columns(&table.columns, &by_type[1], &indel_patterns, false)?;
    let sv_columns = &by_type[2];

    // the feature names after cbind(snv, indel, sv), and the model's variables among them
    let mut feature_names: Vec<String> = SNV_SUBSTITUTIONS.iter().map(|s| s.replacen('>', ".", 1)).collect();
    feature_names.extend(INDEL_NEW_CONTEXTS.iter().map(|(n, _)| n.to_string()));
    feature_names.extend(sv_columns.iter().map(|&i| table.columns[i].clone()));
    let var_index: Vec<usize> = forest
        .variables
        .iter()
        .map(|v| feature_names.iter().position(|f| f == v))
        .collect::<Option<Vec<usize>>>()
        .ok_or_else(|| anyhow!("R: variables in the training data missing in newdata"))?;
    let class_index = |name: &str| forest.classes.iter().position(|c| c == name).ok_or_else(|| anyhow!("R: undefined columns selected ({name})"));
    let (i_brca1, i_brca2) = (class_index("BRCA1")?, class_index("BRCA2")?);

    info!("Predicting HRD probabilities from random forest");
    // max.col's random tie-breaking draws from R's generator; the defined Java reference seeds it with
    // set.seed(1) (CH1)
    let mut rng = rfmt::RRng::set_seed(1);
    let mut out = Vec::with_capacity(table.rows.len());
    for ((_, values), sample) in table.rows.iter().zip(samples) {
        let sum_group = |cols: &[usize]| -> f64 { cols.iter().map(|&i| values[i]).sum() };
        let snv: Vec<f64> = snv_groups.iter().map(|g| sum_group(g)).collect();
        let indel: Vec<f64> = indel_groups.iter().map(|g| sum_group(g)).collect();
        let sv: Vec<f64> = sv_columns.iter().map(|&i| values[i]).collect();
        let mut features = relative(&snv);
        features.extend(relative(&indel));
        features.extend(relative(&sv));
        let x: Vec<f64> = var_index.iter().map(|&i| features[i]).collect();

        let votes = forest.votes(&x);
        let total: f64 = votes.iter().map(|&v| v as f64).sum();
        let p_brca1 = votes[i_brca1] as f64 / total;
        let p_brca2 = votes[i_brca2] as f64 / total;
        let p_hrd = p_brca1 + p_brca2;

        let mut hr_status = if p_hrd >= HRD_PROBABILITY_CUTOFF { HR_DEFICIENT } else { HR_PROFICIENT };
        // max.col(probs[c("BRCA1", "BRCA2")]): ties (within 1e-5 of the larger) are broken at random in R
        let large = p_brca1.abs().max(p_brca2.abs());
        let tol = 1e-5 * large;
        let brca2 = if p_brca2 > p_brca1 + tol {
            true
        } else if p_brca2 >= p_brca1 - tol {
            if compat::java("CH1") { 2.0 * rng.unif_rand() < 1.0 } else { false }
        } else {
            false
        };
        let mut hrd_type = if brca2 { "BRCA2_type" } else { "BRCA1_type" };

        // QC
        let indel_rep = indel[0] + indel[1];
        let indel_load: f64 = indel.iter().sum();
        let sv_load: f64 = sv.iter().sum();
        let has_msi = indel_rep > MSI_MIN_INDEL_REP;
        let has_low_indel_load = indel_load < LOW_INDEL_LOAD_THRES;
        let has_low_sv_load = sv_load < LOW_SV_LOAD_THRES;
        if has_low_indel_load || has_msi {
            hr_status = HR_UNCLEAR;
        }
        if hr_status == HR_PROFICIENT || hr_status == HR_UNCLEAR {
            hrd_type = "none";
        }
        let remarks = format!(
            "{};{}",
            if has_msi { "Has MSI (>14000 indel.rep)" } else { "" },
            if has_low_indel_load { "<100 indels" } else { "" }
        );
        // gsub("^;|;$", "", remarks)
        let remarks_hr_status = if remarks == ";" {
            String::new()
        } else {
            remarks.strip_prefix(';').or_else(|| remarks.strip_suffix(';')).unwrap_or(&remarks).to_string()
        };
        let remarks_hrd_type = if has_low_sv_load && hr_status != HR_PROFICIENT { "<30 SVs" } else { "" };
        debug!("sample({sample}) votes({votes:?}) p_hrd({p_hrd})");
        out.push(Prediction {
            sample,
            p_brca1,
            p_brca2,
            p_hrd,
            hr_status: hr_status.to_string(),
            hrd_type: hrd_type.to_string(),
            remarks_hr_status,
            remarks_hrd_type: remarks_hrd_type.to_string(),
        });
    }
    Ok(out)
}

/// write.table(predictions, sep = "\t", quote = FALSE, row.names = FALSE)
pub fn predictions_text(predictions: &[Prediction]) -> String {
    let mut out = String::from("sample\tp_BRCA1\tp_BRCA2\tp_hrd\thr_status\thrd_type\tremarks_hr_status\tremarks_hrd_type\n");
    for p in predictions {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            p.sample,
            rfmt::write_table_double(p.p_brca1),
            rfmt::write_table_double(p.p_brca2),
            rfmt::write_table_double(p.p_hrd),
            p.hr_status,
            p.hrd_type,
            p.remarks_hr_status,
            p.remarks_hrd_type
        ));
    }
    out
}

/// The model: the embedded forest, or a -chord_model_file (an .rds randomForest object, or a forest
/// exported by export_forest.R)
pub fn load_model(path: Option<&str>) -> Result<Forest> {
    match path {
        None => {
            info!("Loading CHORD model: embedded (CHORD.rds)");
            Forest::embedded()
        }
        Some(p) => {
            info!("Loading CHORD model: {p}");
            let raw = std::fs::read(p).with_context(|| format!("reading {p}"))?;
            if raw.starts_with(&[0x1f, 0x8b]) || raw.starts_with(b"X\n") {
                rds::forest_from_object(&rds::read_rds_bytes(&raw)?)
            } else {
                Forest::parse(std::str::from_utf8(&raw).context("model file is neither RDS nor text")?)
            }
        }
    }
}

/// ChordModel.predict: predicts from the table and writes the predictions file
pub fn predict_to_file(table: &ContextTable, model_path: Option<&str>, output_file: &str) -> Result<()> {
    info!("Starting CHORD predict");
    let forest = load_model(model_path)?;
    let predictions = predict(table, &forest)?;
    info!("Writing CHORD predictions to: {output_file}");
    std::fs::write(output_file, predictions_text(&predictions)).with_context(|| format!("R: cannot open file '{output_file}'"))?;
    info!("Completed CHORD predict");
    Ok(())
}
