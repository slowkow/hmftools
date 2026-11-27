//! pycuppa's CuppaFeaturesLoader (the features file or a directory of cohort files, as pandas reads them)
//! and MissingFeaturesHandler.fill_missing (the features the classifier needs, missing ones set to
//! -1e-8).

use std::collections::HashMap;
use std::io::Read;

use anyhow::{anyhow, bail, Context, Result};

use super::model::Model;
use super::pylog;

pub const PREDICT_NA_FILL_VALUE: f64 = -0.00000001;
pub const SEX_FEATURE_NAME: &str = "event.trait.is_male";

/// Samples (rows) by features (columns), as the transposed DataFrame pycuppa predicts on
pub struct Features {
    pub samples: Vec<String>,
    pub names: Vec<String>,
    /// values[sample][feature]
    pub values: Vec<Vec<f64>>,
}

/// pandas' default NA strings (read_csv na_values)
const NA_VALUES: [&str; 19] = [
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN",
    "None", "n/a", "nan", "null",
];

/// A value of a numeric column as pandas' C parser reads it (correctly rounded for the values Java writes:
/// up to 17 digits, with or without an exponent)
fn parse_value(s: &str) -> Option<f64> {
    if NA_VALUES.contains(&s) {
        return Some(f64::NAN);
    }
    let t = s.trim();
    let t = t.strip_prefix('+').unwrap_or(t);
    match t.to_ascii_lowercase().as_str() {
        "inf" | "infinity" => return Some(f64::INFINITY),
        "-inf" | "-infinity" => return Some(f64::NEG_INFINITY),
        _ => {}
    }
    if t.is_empty() || !t.bytes().all(|b| b.is_ascii_digit() || b == b'.' || b == b'-' || b == b'e' || b == b'E' || b == b'+') {
        return None;
    }
    t.parse::<f64>().ok()
}

/// The lines of a text file, gzip-decompressed if the name ends with .gz (pandas' compression='infer')
pub fn read_text_lines(path: &str) -> Result<Vec<String>> {
    let bytes = std::fs::read(path).with_context(|| format!("FileNotFoundError: [Errno 2] No such file or directory: '{path}'"))?;
    let text = if path.ends_with(".gz") {
        let mut s = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..]).read_to_string(&mut s).with_context(|| format!("{path}: gzip"))?;
        s
    } else if path.ends_with(".bz2") || path.ends_with(".xz") || path.ends_with(".zip") || path.ends_with(".zst") {
        bail!("{path}: only plain and gzip-compressed features files are supported");
    } else {
        String::from_utf8(bytes).with_context(|| format!("{path}: not UTF-8"))?
    };
    // read_table skips blank lines
    Ok(text.lines().map(|l| l.strip_suffix('\r').unwrap_or(l).to_string()).filter(|l| !l.is_empty()).collect())
}

/// One features table: index columns Source, Category, Key, then one value column per sample
struct Table {
    sample_cols: Vec<String>,
    category: Vec<String>,
    key: Vec<String>,
    /// values[row][sample column]
    values: Vec<Vec<f64>>,
}

fn read_table(path: &str) -> Result<Table> {
    let lines = read_text_lines(path)?;
    let Some(header) = lines.first() else { bail!("pandas.errors.EmptyDataError: No columns to parse from file") };
    let cols: Vec<&str> = header.split('\t').collect();
    let find = |n: &str| cols.iter().position(|c| *c == n);
    let missing: Vec<&str> = ["Source", "Category", "Key"].into_iter().filter(|n| find(n).is_none()).collect();
    if !missing.is_empty() {
        bail!("KeyError: 'Missing required columns: {}'", missing.join(", "));
    }
    let (ci, ki) = (find("Category").unwrap(), find("Key").unwrap());
    let si = find("Source").unwrap();
    let value_cols: Vec<usize> = (0..cols.len()).filter(|i| *i != ci && *i != ki && *i != si).collect();
    let mut t = Table {
        sample_cols: value_cols.iter().map(|i| cols[*i].to_string()).collect(),
        category: Vec::with_capacity(lines.len()),
        key: Vec::with_capacity(lines.len()),
        values: Vec::with_capacity(lines.len()),
    };
    for line in &lines[1..] {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() > cols.len() {
            bail!("pandas.errors.ParserError: Error tokenizing data. Expected {} fields, saw {}", cols.len(), f.len());
        }
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        t.category.push(get(ci).to_string());
        t.key.push(get(ki).to_string());
        let mut row = Vec::with_capacity(value_cols.len());
        for (j, i) in value_cols.iter().enumerate() {
            let v = get(*i);
            row.push(parse_value(v).ok_or_else(|| {
                anyhow!("non-numeric value '{v}' of sample column '{}' (pandas reads the column as text, and the prediction fails)", t.sample_cols[j])
            })?);
        }
        t.values.push(row);
    }
    Ok(t)
}

/// CuppaFeaturesDir.FILE_PATTERNS (regular expressions matched at the start of the file names)
const DIR_PATTERNS: [(&str, &str, bool); 6] = [
    ("snv", "cuppa_data.cohort.snv.*.tsv", true),
    ("sv", "cuppa_data.cohort.sv.*.tsv", true),
    ("trait", "cuppa_data.cohort.sample_trait.*.tsv", true),
    ("driver", "cuppa_data.cohort.driver.*.tsv", true),
    ("gene_exp", "cuppa_data.cohort.gene_exp.*.tsv", false),
    ("alt_sj", "cuppa_data.cohort.alt_sj.*.tsv", false),
];

/// `re.match` of a pattern of literal characters, '.' (any character) and '.*' at the start of `s`
fn match_start(pattern: &str, s: &str) -> bool {
    fn m(p: &[u8], s: &[u8]) -> bool {
        if p.is_empty() {
            return true;
        }
        if p.len() >= 2 && p[1] == b'*' {
            // greedy, backtracking
            let mut k = 0;
            loop {
                if m(&p[2..], &s[k..]) {
                    return true;
                }
                if k < s.len() && (p[0] == b'.' || p[0] == s[k]) {
                    k += 1;
                } else {
                    return false;
                }
            }
        }
        !s.is_empty() && (p[0] == b'.' || p[0] == s[0]) && m(&p[1..], &s[1..])
    }
    m(pattern.as_bytes(), s.as_bytes())
}

/// The features of one file (single or multi sample) or of a directory of cohort files
pub fn load(path: &str, sample_id: Option<&str>) -> Result<Features> {
    let meta = std::fs::metadata(path).with_context(|| format!("FileNotFoundError: {path}"))?;
    let (tables, multi) = if meta.is_dir() {
        pylog::info("CuppaFeaturesLoader", "_load_data", &format!("Loading multi sample data from dir: {path}"));
        let mut names: Vec<String> = std::fs::read_dir(path)?.filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().to_string()).collect();
        // os.listdir order is the file system's; sorted here (the first match is used if several match)
        names.sort();
        let mut tables = Vec::new();
        let mut missing_required = false;
        for (key, pattern, required) in DIR_PATTERNS {
            let matched: Vec<&String> = names.iter().filter(|n| match_start(pattern, n)).collect();
            if matched.is_empty() {
                if required {
                    pylog::error("CuppaFeaturesDir", "get_file_paths", &format!("Missing required input file type '{key}' with pattern '{pattern}'"));
                    missing_required = true;
                } else {
                    pylog::warning("CuppaFeaturesDir", "get_file_paths", &format!("Missing optional input file type '{key}' with pattern '{pattern}'"));
                }
                continue;
            }
            if matched.len() > 1 {
                let list: Vec<&str> = matched.iter().map(|s| s.as_str()).collect();
                pylog::warning("CuppaFeaturesDir", "get_file_paths", &format!("Pattern '{pattern}' matched multiple files. Using the first: {}", list.join(", ")));
            }
            let file = format!("{}/{}", path.trim_end_matches('/'), matched[0]);
            pylog::debug("CuppaFeaturesLoader", "_load_files_from_dir", &format!("Loading file: {}", matched[0]));
            tables.push(read_table(&file)?);
        }
        if missing_required {
            bail!("FileNotFoundError");
        }
        (tables, true)
    } else {
        let t = read_table(path)?;
        let multi = t.sample_cols.len() > 1;
        pylog::info("CuppaFeaturesLoader", "_load_data", &format!("Loading {} sample data file: {path}", if multi { "multi" } else { "single" }));
        (vec![t], multi)
    };

    // pd.concat of the tables: the union of the sample columns, in order of appearance
    let mut samples: Vec<String> = Vec::new();
    for t in &tables {
        for c in &t.sample_cols {
            if !samples.contains(c) {
                samples.push(c.clone());
            }
        }
    }
    if !multi {
        let Some(id) = sample_id else {
            pylog::error("CuppaFeaturesLoader", "_load_one_file", "`sample_id` must be provided when features are for one sample");
            bail!("ValueError");
        };
        samples = vec![id.to_string()];
    }
    let mut category = Vec::new();
    let mut key = Vec::new();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for t in &tables {
        let pos: Vec<usize> = if multi { t.sample_cols.iter().map(|c| samples.iter().position(|s| s == c).unwrap()).collect() } else { vec![0] };
        for r in 0..t.key.len() {
            category.push(t.category[r].clone());
            key.push(t.key[r].clone());
            let mut row = vec![f64::NAN; samples.len()];
            for (j, p) in pos.iter().enumerate() {
                row[*p] = t.values[r][j];
            }
            rows.push(row);
        }
    }

    // _assert_no_duplicate_features (on Key alone)
    let mut seen: HashMap<&str, ()> = HashMap::new();
    let mut dups: Vec<&str> = Vec::new();
    for k in &key {
        if seen.insert(k.as_str(), ()).is_some() && !dups.contains(&k.as_str()) {
            dups.push(k);
        }
    }
    if !dups.is_empty() {
        pylog::error("CuppaFeaturesLoader", "_assert_no_duplicate_features", "The following features are duplicated: ");
        bail!("AssertionError: duplicated features: {}", dups.join(", "));
    }

    // _check_gen_pos_chroms
    let valid: Vec<String> = (1..=22).map(|i| i.to_string()).chain(["X".to_string()]).collect();
    let mut bad: Vec<&str> = Vec::new();
    for (c, k) in category.iter().zip(&key) {
        if c == "gen_pos" {
            let chrom = k.split('_').next().unwrap_or("");
            let plain = chrom.strip_prefix("chr").unwrap_or("");
            let ok = valid.iter().any(|v| v == chrom) || (chrom.starts_with("chr") && valid.iter().any(|v| v == plain));
            if !ok && !bad.contains(&chrom) {
                bad.push(chrom);
            }
        }
    }
    if !bad.is_empty() {
        pylog::error(
            "CuppaFeaturesLoader",
            "_check_gen_pos_chroms",
            &format!("gen_pos bins must contain only chromosomes 1-22 and X. Found invalid chromosomes: {}", bad.join(", ")),
        );
        bail!("ValueError");
    }

    // _parse_signatures (on Key, whatever the category), then names Category.Key
    const EXCLUDED: [&str; 3] = ["SIG_1", "SIG_10_POLE", "SIG_11"];
    const MAPPINGS: [(&str, &str); 8] = [
        ("SIG_1", "Age (SBS1)"),
        ("SIG_2_13_AID_APOBEC", "AID/APOBEC (SBS2/13)"),
        ("SIG_4_SMOKING", "Smoking (SBS4)"),
        ("SIG_6_MMR", "MMRD (SBS6)"),
        ("SIG_7_UV", "UV (SBS7)"),
        ("SIG_10_POLE", "POLE (SBS10)"),
        ("SIG_11", "Temozolomide (SBS11)"),
        ("SIG_17", "ROS/5FU (SBS17)"),
    ];
    let mut names = Vec::new();
    let mut kept: Vec<Vec<f64>> = Vec::new();
    for ((c, k), row) in category.iter().zip(&key).zip(rows) {
        if EXCLUDED.contains(&k.as_str()) {
            continue;
        }
        let k2 = MAPPINGS.iter().find(|(a, _)| *a == k).map(|(_, b)| b.to_string()).unwrap_or_else(|| k.clone());
        names.push(format!("{c}.{k2}"));
        kept.push(row);
    }
    let values: Vec<Vec<f64>> = (0..samples.len()).map(|s| kept.iter().map(|r| r[s]).collect()).collect();
    if multi {
        pylog::debug("CuppaFeaturesLoader", "_print_stats", &format!("Loaded {} features from {} samples", names.len(), samples.len()));
    } else {
        pylog::debug("CuppaFeaturesLoader", "_print_stats", &format!("Loaded {} features from sample({})", names.len(), samples[0]));
    }
    Ok(Features { samples, names, values })
}

/// The feature matrix the classifier works on (MissingFeaturesHandler.fill_missing)
pub struct Prepared {
    pub samples: Vec<String>,
    /// per sample: the required gen_pos, snv96, event and sig features, in that order
    pub dna: Vec<Vec<f64>>,
    /// per sample with RNA data: the required gene_exp and alt_sj features; None for a sample without any
    pub rna: Vec<Option<Vec<f64>>>,
    /// per sample: whether event.trait.is_male was missing (or NA) in the features
    pub sex_missing: Vec<bool>,
    pub n_gen_pos: usize,
    pub n_snv96: usize,
    pub n_event: usize,
    pub n_sig: usize,
    pub n_gene_exp: usize,
}

impl Prepared {
    pub fn gen_pos(&self, s: usize) -> &[f64] {
        &self.dna[s][..self.n_gen_pos]
    }
    pub fn snv96(&self, s: usize) -> &[f64] {
        let o = self.n_gen_pos;
        &self.dna[s][o..o + self.n_snv96]
    }
    pub fn event(&self, s: usize) -> &[f64] {
        let o = self.n_gen_pos + self.n_snv96;
        &self.dna[s][o..o + self.n_event]
    }
    pub fn sig(&self, s: usize) -> &[f64] {
        let o = self.n_gen_pos + self.n_snv96 + self.n_event;
        &self.dna[s][o..o + self.n_sig]
    }
    pub fn gene_exp(&self, s: usize) -> Option<&[f64]> {
        self.rna[s].as_deref().map(|r| &r[..self.n_gene_exp])
    }
    pub fn alt_sj(&self, s: usize) -> Option<&[f64]> {
        self.rna[s].as_deref().map(|r| &r[self.n_gene_exp..])
    }
}

/// counts of feature names by their leading word (MissingFeaturesHandler._get_feat_type_counts_string)
fn type_counts(names: &[&str]) -> String {
    let mut types: Vec<String> = names
        .iter()
        .map(|n| n.chars().take_while(|c| c.is_alphanumeric() || *c == '_').collect::<String>())
        .collect();
    types.sort();
    let mut out: Vec<String> = Vec::new();
    let mut i = 0;
    while i < types.len() {
        let mut j = i;
        while j < types.len() && types[j] == types[i] {
            j += 1;
        }
        out.push(format!("{}={}", types[i], j - i));
        i = j;
    }
    out.join(", ")
}

/// `re.search("^gene_exp|alt_sj", name)`
pub fn is_rna_feature(name: &str) -> bool {
    name.starts_with("gene_exp") || name.contains("alt_sj")
}

pub fn fill_missing(f: &Features, m: &Model) -> Prepared {
    let index: HashMap<&str, usize> = f.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let dna_names: Vec<&String> = m.req_gen_pos.iter().chain(&m.req_snv96).chain(&m.req_event).chain(&m.req_sig).collect();
    let rna_names: Vec<&String> = m.req_gene_exp.iter().chain(&m.req_alt_sj).collect();
    let dna_idx: Vec<Option<usize>> = dna_names.iter().map(|n| index.get(n.as_str()).copied()).collect();
    let rna_idx: Vec<Option<usize>> = rna_names.iter().map(|n| index.get(n.as_str()).copied()).collect();
    let fill = |v: f64| if v.is_nan() { PREDICT_NA_FILL_VALUE } else { v };

    let rna_cols: Vec<usize> = (0..f.names.len()).filter(|i| is_rna_feature(&f.names[*i])).collect();
    let sex_idx = index.get(SEX_FEATURE_NAME).copied();
    let mut dna = Vec::with_capacity(f.samples.len());
    let mut rna = Vec::with_capacity(f.samples.len());
    let mut sex_missing = Vec::with_capacity(f.samples.len());
    for s in 0..f.samples.len() {
        let x = &f.values[s];
        dna.push(dna_idx.iter().map(|i| i.map_or(PREDICT_NA_FILL_VALUE, |i| fill(x[i]))).collect());
        // samples whose RNA features are all NA (or that have none) get no RNA row
        let has_rna = !rna_cols.is_empty() && rna_cols.iter().any(|i| !x[*i].is_nan());
        // the alt-SJ classifier was trained on ln(count + 1), but CuppaDataPrep writes the counts and pycuppa
        // takes them as they are (CU8); corrected, the counts are log-transformed here
        let log_alt_sj = !crate::compat::java("CU8");
        let n_ge = m.req_gene_exp.len();
        rna.push(if has_rna {
            Some(
                rna_idx
                    .iter()
                    .enumerate()
                    .map(|(k, i)| match i {
                        None => PREDICT_NA_FILL_VALUE,
                        Some(i) if x[*i].is_nan() => PREDICT_NA_FILL_VALUE,
                        Some(i) if log_alt_sj && k >= n_ge => x[*i].ln_1p(),
                        Some(i) => x[*i],
                    })
                    .collect(),
            )
        } else {
            None
        });
        sex_missing.push(sex_idx.is_none_or(|i| x[i].is_nan()));
    }

    let required: std::collections::HashSet<&str> = dna_names.iter().chain(&rna_names).map(|s| s.as_str()).collect();
    let removed: Vec<&str> = f.names.iter().map(|s| s.as_str()).filter(|n| !required.contains(n)).collect();
    if !removed.is_empty() {
        pylog::info("MissingFeaturesHandler", "fill_missing", &format!("Removed {} features that were not required: {}", removed.len(), type_counts(&removed)));
    }
    let added: Vec<&str> = dna_names.iter().chain(&rna_names).map(|s| s.as_str()).filter(|n| !index.contains_key(n)).collect();
    if !added.is_empty() {
        pylog::info(
            "MissingFeaturesHandler",
            "fill_missing",
            &format!("Filled in {} features missing with value={}: {}", added.len(), super::pyfmt::repr(PREDICT_NA_FILL_VALUE), type_counts(&added)),
        );
    }

    Prepared {
        samples: f.samples.clone(),
        dna,
        rna,
        sex_missing,
        n_gen_pos: m.req_gen_pos.len(),
        n_snv96: m.req_snv96.len(),
        n_event: m.req_event.len(),
        n_sig: m.req_sig.len(),
        n_gene_exp: m.req_gene_exp.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn patterns() {
        assert!(match_start("cuppa_data.cohort.snv.*.tsv", "cuppa_data.cohort.snv.tsv.gz"));
        assert!(match_start("cuppa_data.cohort.snv.*.tsv", "cuppa_data.cohort.snv.run1.tsv.gz"));
        assert!(!match_start("cuppa_data.cohort.snv.*.tsv", "cuppa_data.cohort.sv.tsv.gz"));
        assert!(match_start("cuppa_data.cohort.sv.*.tsv", "cuppa_data.cohort.sv.tsv"));
    }

    #[test]
    fn values() {
        assert_eq!(parse_value("2521.7"), Some(2521.7));
        assert_eq!(parse_value("1.234e+00"), Some(1.234));
        assert!(parse_value("null").unwrap().is_nan());
        assert!(parse_value("abc").is_none());
    }

    #[test]
    fn counts() {
        assert_eq!(type_counts(&["gen_pos.chrY_0", "event.driver.A.mut", "gen_pos.chrY_500000"]), "event=1, gen_pos=2");
    }
}
