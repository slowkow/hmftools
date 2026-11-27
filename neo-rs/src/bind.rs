//! NEO's binding module (`com.hartwig.hmftools.neo.bind`), the parts NeoScorer uses: the reference files
//! (position weight matrices, flank weights, the random peptide score and likelihood distributions, the
//! binding likelihoods, the immunogenic recognition peptides, the expression likelihoods), ScoreConfig and
//! BindScorer.calcScoreData. The training classes (BindTrainer, PosWeightModel, NoiseModel, ...) are
//! not ported.

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;

use hmf_common::blosum::BlosumMapping;
use hmf_common::config_builder::{check_add_dir_separator, ConfigBuilder};
use hmf_common::java_text::{fields_index_map, parse_double, parse_int, read_lines, split_all};
use hmf_common::{error, info};

use crate::compat;

/// BindConstants.INVALID_SCORE
pub const INVALID_SCORE: f64 = -1.0;
/// BindScorer.INVALID_CALC
pub const INVALID_CALC: f64 = -1.0;
/// BindScoreMatrix.INVALID_SCORE
const MATRIX_INVALID_SCORE: f64 = -1000.0;

/// BindConstants
pub const AMINO_ACIDS: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
pub const AMINO_ACID_COUNT: usize = 20;
pub const MIN_PEPTIDE_LENGTH: i32 = 8;
pub const REF_PEPTIDE_LENGTH: i32 = 12;
pub const DEFAULT_PEPTIDE_LENGTHS: [i32; 5] = [8, 9, 10, 11, 12];
/// BindCommon.AMINO_ACID_21ST
pub const AMINO_ACID_21ST: &str = "X";

/// FlankCounts
pub const START_AMINO_ACID_ID: u8 = b'm';
/// Codons.STOP_AMINO_ACID
pub const STOP_AMINO_ACID: u8 = b'X';
pub const FLANK_AA_COUNT: i32 = 3;
const TOTAL_FLANK_AA_COUNT: usize = 6;
const FLANK_AMINO_ACID_COUNT: usize = 22;
const UP_1: usize = 2;
const DOWN_1: usize = 3;

/// BindConstants.aminoAcidIndex
#[inline]
pub fn amino_acid_index(aa: u8) -> Option<usize> {
    AMINO_ACIDS.iter().position(|&a| a == aa)
}

/// FlankCounts.flankAminoAcidIndex: the 20 amino acids, then 'm' (start) and 'X' (stop)
#[inline]
fn flank_amino_acid_index(aa: u8) -> Option<usize> {
    match aa {
        START_AMINO_ACID_ID => Some(20),
        STOP_AMINO_ACID => Some(21),
        _ => amino_acid_index(aa),
    }
}

/// BindCommon.cleanAllele
pub fn clean_allele(allele: &str) -> String {
    allele.replace("HLA-", "").replace(':', "").replace('*', "")
}

/// Doubles.equal
#[inline]
fn doubles_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-10
}

fn num(s: &str) -> Result<f64> {
    parse_double(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
}

fn int(s: &str) -> Result<i32> {
    parse_int(s).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{s}\""))
}

fn col(f: &HashMap<String, usize>, name: &str, file: &str) -> Result<usize> {
    f.get(name).copied().ok_or_else(|| anyhow!("{file}: NullPointerException: no column {name}"))
}

fn item<'a>(v: &[&'a str], i: usize, file: &str) -> Result<&'a str> {
    v.get(i).copied().ok_or_else(|| anyhow!("{file}: ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len()))
}

/// Files.readAllLines with the header split off: Ok(None) where Java catches an IOException (logged here as
/// `<what>: <error>`), Err for an empty file (an IndexOutOfBoundsException in Java)
fn header_and_lines(filename: &str, what: &str) -> Result<Option<(HashMap<String, usize>, Vec<String>)>> {
    let mut lines = match read_lines(filename) {
        Ok(l) => l,
        Err(e) => {
            error!("{}: java.nio.file.NoSuchFileException: {} ({})", what, filename, e);
            return Ok(None);
        }
    };
    if lines.is_empty() {
        bail!("IndexOutOfBoundsException: Index 0 out of bounds for length 0 ({filename})");
    }
    let header = lines.remove(0);
    Ok(Some((fields_index_map(&header, ','), lines)))
}

// ---------------------------------------------------------------------------------------------
// ScoreConfig, RandomPeptideConfig (for reading)

/// TrainConfig file ids
pub const FILE_ID_POS_WEIGHT: &str = "pos_weight";
pub const FILE_ID_LIKELIHOOD: &str = "likelihood";
pub const FILE_ID_FLANK_POS_WEIGHT: &str = "flank_pos_weight";
pub const FILE_ID_RECOGNITION: &str = "recognition";
pub const FILE_ID_EXPRESSION_DIST: &str = "expression_dist";
pub const FILE_ID_RAND_DIST: &str = "rand_dist";
pub const FILE_ID_LIKELIHOOD_DIST: &str = "likelihood_rand_dist";
pub const FILE_ID_EXP_LIKELIHOOD_DIST: &str = "exp_likelihood_rand_dist";

/// ScoreConfig.scoreFileConfig
fn score_file_config(file_type: &str) -> String {
    format!("{file_type}_file")
}

/// TrainConfig.formTrainingFilename (BindCommon.formFilename with "train_<type>")
pub fn form_training_filename(dir: &str, file_type: &str, file_id: Option<&str>) -> String {
    match file_id {
        Some(id) => format!("{dir}neo_train_{file_type}_{id}.csv"),
        None => format!("{dir}neo_train_{file_type}.csv"),
    }
}

/// ScoreConfig.getScoringFilename: the file option, else the file of that type in the score file directory.
/// With a directory, Java prefixes the directory to a file option's value as well (NE8), though ConfigBuilder
/// has already checked that the value is an existing path as given; corrected, the value is used as given.
pub fn get_scoring_filename(cb: &ConfigBuilder, score_file_dir: Option<&str>, score_file_id: Option<&str>, file_type: &str) -> Option<String> {
    let config_value = cb.value(&score_file_config(file_type));
    let Some(dir) = score_file_dir else { return config_value };
    match config_value {
        Some(v) => Some(if compat::java("NE8") { format!("{dir}{v}") } else { v }),
        None => Some(form_training_filename(dir, file_type, score_file_id)),
    }
}

/// ScoreConfig.registerConfig (with RandomPeptideConfig.addConfigForReading)
pub fn register_score_config(cb: &mut ConfigBuilder) {
    cb.add_path(&score_file_config(FILE_ID_RAND_DIST), false, "Random peptide distribution file");
    cb.add_path(&score_file_config(FILE_ID_LIKELIHOOD_DIST), false, "Random likelihood distribution file");
    cb.add_path(&score_file_config(FILE_ID_EXP_LIKELIHOOD_DIST), false, "Random expression likelihood distribution file");
    cb.add_str("score_file_id", false, "Reference file id for scoring instead of specifying individual files");
    cb.add_path("score_file_dir", false, "Reference file directory");
    cb.add_path(&score_file_config(FILE_ID_POS_WEIGHT), false, "Binding position weights file");
    cb.add_path(&score_file_config(FILE_ID_FLANK_POS_WEIGHT), false, "Binding flank weights file");
    cb.add_path(&score_file_config(FILE_ID_LIKELIHOOD), false, "Binding likelihood file");
    cb.add_path(&score_file_config(FILE_ID_EXPRESSION_DIST), false, "Expression likelihood file");
    cb.add_path(&score_file_config(FILE_ID_RECOGNITION), false, "Immunogenic pHLA recognition data file");
    cb.add_path("validation_data_file", false, "Validation data file");
    cb.add_flag("write_summary_data", "Write summary results per allele");
    cb.add_flag("write_peptide_scores", "Write score and rank data per peptide");
    cb.add_flag("check_self_recog", "Don't allow recognition with same peptide");
}

/// ScoreConfig, the fields NeoScorer reads (validation data, the peptide score and summary files belong
/// to BindScorer's own main class)
#[derive(Clone, Debug, Default)]
pub struct ScoreConfig {
    pub recognition_data_file: Option<String>,
    pub check_self_recognition: bool,
    pub pos_weights_file: Option<String>,
    pub bind_likelihood_file: Option<String>,
    pub flank_pos_weights_file: Option<String>,
    pub expression_likelihood_file: Option<String>,
    /// RandomPeptideConfig.forReading
    pub score_distribution_file: Option<String>,
    pub likelihood_distribution_file: Option<String>,
    pub expression_likelihood_distribution_file: Option<String>,
}

impl ScoreConfig {
    pub fn new(cb: &ConfigBuilder) -> ScoreConfig {
        let validation_data_file = cb.value("validation_data_file");
        // checkAddDirSeparator: null stays null, "" stays ""
        let dir = cb.value("score_file_dir").map(|d| if d.is_empty() { d } else { check_add_dir_separator(&d) });
        let id = cb.value("score_file_id");
        let f = |t: &str| get_scoring_filename(cb, dir.as_deref(), id.as_deref(), t);
        let recognition_data_file = f(FILE_ID_RECOGNITION);
        let check_self_recognition = cb.has_flag("check_self_recog")
            || (validation_data_file.is_some() && validation_data_file == recognition_data_file);
        ScoreConfig {
            pos_weights_file: f(FILE_ID_POS_WEIGHT),
            flank_pos_weights_file: f(FILE_ID_FLANK_POS_WEIGHT),
            bind_likelihood_file: f(FILE_ID_LIKELIHOOD),
            expression_likelihood_file: f(FILE_ID_EXPRESSION_DIST),
            score_distribution_file: f(FILE_ID_RAND_DIST),
            likelihood_distribution_file: f(FILE_ID_LIKELIHOOD_DIST),
            expression_likelihood_distribution_file: f(FILE_ID_EXP_LIKELIHOOD_DIST),
            recognition_data_file,
            check_self_recognition,
        }
    }
}

// ---------------------------------------------------------------------------------------------
// BindData

/// BindData: an allele and peptide with its flanks, TPM and scores
#[derive(Clone, Debug)]
pub struct BindData {
    pub allele: String,
    pub peptide: String,
    pub up_flank: String,
    pub down_flank: String,
    tpm: f64,
    pub score: f64,
    pub flank_score: f64,
    pub rank_percentile: f64,
    pub likelihood: f64,
    pub likelihood_rank: f64,
    pub exp_likelihood: f64,
    pub exp_likelihood_rank: f64,
    pub recognition_similarity: f64,
    pub other_allele_recognition_similarity: f64,
}

impl BindData {
    pub fn new(allele: &str, peptide: &str, up_flank: &str, down_flank: &str) -> BindData {
        BindData {
            allele: allele.to_string(),
            peptide: peptide.to_string(),
            up_flank: up_flank.to_string(),
            down_flank: down_flank.to_string(),
            tpm: 0.0,
            score: INVALID_SCORE,
            flank_score: INVALID_SCORE,
            rank_percentile: INVALID_SCORE,
            likelihood: INVALID_SCORE,
            likelihood_rank: INVALID_SCORE,
            exp_likelihood: INVALID_SCORE,
            exp_likelihood_rank: INVALID_SCORE,
            recognition_similarity: INVALID_SCORE,
            other_allele_recognition_similarity: INVALID_SCORE,
        }
    }
    pub fn peptide_length(&self) -> i32 {
        self.peptide.len() as i32
    }
    pub fn has_flanks(&self) -> bool {
        !self.up_flank.is_empty() || !self.down_flank.is_empty()
    }
    /// setTPM: Java tests the current TPM (0 at first, so never negative) rather than the new one; a negative
    /// TPM (no expression data) is stored and gives the lowest expression bucket, as 0 would
    pub fn set_tpm(&mut self, tpm: f64) {
        if self.tpm < 0.0 {
            return;
        }
        self.tpm = tpm;
    }
    pub fn tpm(&self) -> f64 {
        self.tpm
    }
}

// ---------------------------------------------------------------------------------------------
// BindScoreMatrix

/// BindScoreMatrix: position weights by amino acid and position for one allele and peptide length
#[derive(Clone, Debug)]
pub struct BindScoreMatrix {
    pub allele: String,
    pub peptide_length: i32,
    /// [amino acid][position], flattened
    weights: Vec<f64>,
}

impl BindScoreMatrix {
    pub fn new(allele: &str, peptide_length: i32) -> BindScoreMatrix {
        BindScoreMatrix {
            allele: allele.to_string(),
            peptide_length,
            weights: vec![0.0; AMINO_ACID_COUNT * peptide_length.max(0) as usize],
        }
    }
    #[inline]
    pub fn weight(&self, aa: usize, pos: usize) -> f64 {
        self.weights[aa * self.peptide_length as usize + pos]
    }
    pub fn set_weight(&mut self, aa: usize, pos: usize, value: f64) {
        let l = self.peptide_length as usize;
        self.weights[aa * l + pos] = value;
    }

    /// calcScore: the sum of the position weights, -1000 for another length or an unknown amino acid
    pub fn calc_score(&self, peptide: &str) -> f64 {
        let p = peptide.as_bytes();
        if p.len() as i32 != self.peptide_length {
            return MATRIX_INVALID_SCORE;
        }
        let mut score = 0.0;
        for (i, &aa) in p.iter().enumerate() {
            let Some(aa_index) = amino_acid_index(aa) else { return MATRIX_INVALID_SCORE };
            score += self.weight(aa_index, i);
        }
        score
    }

    /// loadFromCsv: the PosWeights rows; a matrix per run of rows with the same allele and length. A
    /// missing file is logged and gives no matrices, as in Java.
    pub fn load_from_csv(filename: &str) -> Result<Vec<BindScoreMatrix>> {
        let mut matrices: Vec<BindScoreMatrix> = Vec::new();
        let Some((f, lines)) = header_and_lines(filename, "failed to load pos-weights data file")? else { return Ok(matrices) };
        let allele_i = col(&f, "Allele", filename)?;
        let data_type_i = col(&f, "DataType", filename)?;
        let pep_len_i = col(&f, "PeptideLength", filename)?;
        let aa_i = col(&f, "AminoAcid", filename)?;
        let start = aa_i.max(data_type_i) + 1;
        let mut alleles = hmf_common::java_hash::JavaHashSet::<String>::new();
        for line in &lines {
            let items = split_all(line, ',');
            if item(&items, data_type_i, filename)? != "PosWeights" {
                continue;
            }
            let allele = item(&items, allele_i, filename)?;
            alleles.insert(allele.to_string());
            let peptide_length = int(item(&items, pep_len_i, filename)?)?;
            let new_matrix = match matrices.last() {
                None => true,
                Some(m) => m.allele != allele || m.peptide_length != peptide_length,
            };
            if new_matrix {
                matrices.push(BindScoreMatrix::new(allele, peptide_length));
            }
            let m = matrices.last_mut().unwrap();
            let aa = item(&items, aa_i, filename)?.as_bytes().first().copied().unwrap_or(0);
            // aminoAcidIndex gives -1 for an unknown letter, an ArrayIndexOutOfBoundsException in Java
            let aa_index = amino_acid_index(aa).ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: Index -1"))?;
            let mut pos = 0usize;
            for v in items.iter().skip(start) {
                if pos >= m.peptide_length as usize {
                    bail!("{filename}: ArrayIndexOutOfBoundsException: Index {pos}");
                }
                m.set_weight(aa_index, pos, num(v)?);
                if pos as i32 == m.peptide_length - 1 {
                    break;
                }
                pos += 1;
            }
        }
        info!("loading matrix data for {} alleles from {}", alleles.len(), filename);
        Ok(matrices)
    }
}

// ---------------------------------------------------------------------------------------------
// FlankScores

/// FlankScores: position weights of the 3 amino acids either side of a peptide
#[derive(Clone, Debug)]
pub struct FlankScores {
    weights: [[f64; TOTAL_FLANK_AA_COUNT]; FLANK_AMINO_ACID_COUNT],
    has_data: bool,
}

impl Default for FlankScores {
    fn default() -> Self {
        FlankScores { weights: [[0.0; TOTAL_FLANK_AA_COUNT]; FLANK_AMINO_ACID_COUNT], has_data: false }
    }
}

impl FlankScores {
    pub fn has_data(&self) -> bool {
        self.has_data
    }

    #[inline]
    fn pos_weight(&self, position: usize, aa: u8) -> f64 {
        match flank_amino_acid_index(aa) {
            Some(i) => self.weights[i][position],
            None => 0.0,
        }
    }

    /// calcScore: up flank at positions 0-2 (right-aligned), down flank at 3-5; the start or stop
    /// placeholder where a flank is empty
    pub fn calc_score(&self, up_flank: &str, down_flank: &str) -> f64 {
        if !self.has_data {
            return 0.0;
        }
        if up_flank.is_empty() && down_flank.is_empty() {
            return 0.0;
        }
        let (u, d) = (up_flank.as_bytes(), down_flank.as_bytes());
        let mut score = 0.0;
        if !u.is_empty() {
            if u.len() >= 3 {
                score += self.pos_weight(0, u[0]);
                score += self.pos_weight(1, u[1]);
                score += self.pos_weight(2, u[2]);
            } else if u.len() >= 2 {
                score += self.pos_weight(1, u[0]);
                score += self.pos_weight(2, u[1]);
            } else {
                score += self.pos_weight(2, u[0]);
            }
        } else {
            score += self.pos_weight(UP_1, START_AMINO_ACID_ID);
        }
        if !d.is_empty() {
            if d.len() >= 3 {
                score += self.pos_weight(3, d[0]);
                score += self.pos_weight(4, d[1]);
                score += self.pos_weight(5, d[2]);
            } else if d.len() >= 2 {
                score += self.pos_weight(3, d[0]);
                score += self.pos_weight(4, d[1]);
            } else {
                score += self.pos_weight(3, d[0]);
            }
        } else {
            score += self.pos_weight(DOWN_1, STOP_AMINO_ACID);
        }
        score
    }

    /// loadPosWeights
    pub fn load_pos_weights(&mut self, filename: &str) -> Result<bool> {
        let Some((f, lines)) = header_and_lines(filename, "failed to load flanking pos-weights data file")? else { return Ok(false) };
        let aa_i = col(&f, "AminoAcid", filename)?;
        let data_type_i = col(&f, "DataType", filename)?;
        for line in &lines {
            let items = split_all(line, ',');
            if item(&items, data_type_i, filename)? != "PosWeights" {
                continue;
            }
            let aa = item(&items, aa_i, filename)?.as_bytes().first().copied().unwrap_or(0);
            let aa_index = flank_amino_acid_index(aa).ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: Index -1"))?;
            for (pos, v) in items.iter().skip(aa_i + 1).enumerate() {
                if pos >= TOTAL_FLANK_AA_COUNT {
                    bail!("{filename}: ArrayIndexOutOfBoundsException: Index {pos}");
                }
                self.weights[aa_index][pos] = num(v)?;
            }
        }
        self.has_data = true;
        Ok(true)
    }
}

// ---------------------------------------------------------------------------------------------
// RandomPeptideDistribution (reading)

/// ScoreDistributionData: a bucket (rank) and its score
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScoreDistributionData {
    pub score_bucket: f64,
    pub score: f64,
}

/// RandomPeptideDistribution, as NeoScorer reads it: the score distribution per allele and peptide length,
/// and the likelihood and expression likelihood distributions per allele
#[derive(Clone, Debug, Default)]
pub struct RandomPeptideDistribution {
    score_distributions: HashMap<String, HashMap<i32, Vec<ScoreDistributionData>>>,
    likelihood_distributions: HashMap<String, Vec<ScoreDistributionData>>,
    exp_likelihood_distributions: HashMap<String, Vec<ScoreDistributionData>>,
    data_loaded: bool,
}

impl RandomPeptideDistribution {
    /// loadData: false where Java's is false (a missing score distribution or likelihood distribution file)
    pub fn load_data(&mut self, config: &ScoreConfig) -> Result<bool> {
        self.data_loaded = self.load_distribution(config.score_distribution_file.as_deref())?
            && load_likelihood_distribution(config.likelihood_distribution_file.as_deref(), &mut self.likelihood_distributions, true)?
            && load_likelihood_distribution(
                config.expression_likelihood_distribution_file.as_deref(),
                &mut self.exp_likelihood_distributions,
                false,
            )?;
        Ok(self.data_loaded)
    }

    #[allow(dead_code)]
    pub fn has_data(&self) -> bool {
        self.data_loaded
    }

    fn load_distribution(&mut self, filename: Option<&str>) -> Result<bool> {
        let Some(filename) = filename else { return Ok(false) };
        let Some((f, lines)) = header_and_lines(filename, "failed to read random distribution file")? else { return Ok(false) };
        let (allele_i, len_i, bucket_i, score_i) =
            (col(&f, "Allele", filename)?, col(&f, "PeptideLength", filename)?, col(&f, "ScoreBucket", filename)?, col(&f, "Score", filename)?);
        for line in &lines {
            let items = split_all(line, ',');
            let allele = item(&items, allele_i, filename)?;
            let peptide_length = int(item(&items, len_i, filename)?)?;
            let data = ScoreDistributionData { score_bucket: num(item(&items, bucket_i, filename)?)?, score: num(item(&items, score_i, filename)?)? };
            self.score_distributions.entry(allele.to_string()).or_default().entry(peptide_length).or_default().push(data);
        }
        info!("loaded {} alleles with random peptide score distribution data", self.score_distributions.len());
        Ok(true)
    }

    /// getScoreRank
    pub fn score_rank(&self, allele: &str, peptide_length: i32, score: f64) -> f64 {
        let Some(m) = self.score_distributions.get(allele) else { return INVALID_SCORE };
        match m.get(&peptide_length) {
            Some(scores) if scores.len() >= 2 => get_rank(scores, score),
            _ => INVALID_SCORE,
        }
    }

    /// getLikelihoodRank
    pub fn likelihood_rank(&self, allele: &str, likelihood: f64) -> f64 {
        match self.likelihood_distributions.get(allele) {
            Some(d) if d.len() >= 2 => get_rank(d, likelihood),
            _ => INVALID_SCORE,
        }
    }

    /// getExpressionLikelihoodRank
    pub fn expression_likelihood_rank(&self, allele: &str, likelihood: f64) -> f64 {
        match self.exp_likelihood_distributions.get(allele) {
            Some(d) if d.len() >= 2 => get_rank(d, likelihood),
            _ => INVALID_SCORE,
        }
    }
}

/// RandomPeptideDistribution.loadLikelihoodDistribution: true without a file name, or without the file when
/// it is not required
fn load_likelihood_distribution(filename: Option<&str>, dists: &mut HashMap<String, Vec<ScoreDistributionData>>, require_exists: bool) -> Result<bool> {
    let Some(filename) = filename else { return Ok(true) };
    if !std::path::Path::new(filename).exists() {
        return Ok(!require_exists);
    }
    let Some((f, lines)) = header_and_lines(filename, "failed to read random likelihood distribution file")? else { return Ok(false) };
    let (allele_i, bucket_i, lik_i) = (col(&f, "Allele", filename)?, col(&f, "Bucket", filename)?, col(&f, "Likelihood", filename)?);
    for line in &lines {
        let items = split_all(line, ',');
        let allele = item(&items, allele_i, filename)?;
        let data = ScoreDistributionData { score_bucket: num(item(&items, bucket_i, filename)?)?, score: num(item(&items, lik_i, filename)?)? };
        dists.entry(allele.to_string()).or_default().push(data);
    }
    info!("loaded {} alleles with random peptide likelihood distribution data", dists.len());
    Ok(true)
}

/// RandomPeptideDistribution.getRank: the bucket of a score in an ascending or descending distribution,
/// interpolated between buckets; 0 above the best score and 1 below the worst
pub fn get_rank(distribution: &[ScoreDistributionData], score: f64) -> f64 {
    let is_ascending = distribution[0].score < distribution[1].score;
    if (is_ascending && score < distribution[0].score) || (!is_ascending && score > distribution[0].score) {
        return 0.0;
    }
    let n = distribution.len();
    if (is_ascending && score > distribution[n - 1].score) || (!is_ascending && score < distribution[n - 1].score) {
        return 1.0;
    }
    for i in 0..n {
        let d = distribution[i];
        if doubles_equal(score, d.score) {
            return d.score_bucket;
        }
        let next = if i < n - 1 { Some(distribution[i + 1]) } else { None };
        if let Some(nx) = next {
            if doubles_equal(score, nx.score) {
                return nx.score_bucket;
            }
        }
        if (is_ascending && score > d.score) || (!is_ascending && score < d.score) {
            let Some(nx) = next else { break };
            if (is_ascending && score < nx.score) || (!is_ascending && score > nx.score) {
                if is_ascending {
                    let upper_perc = (score - d.score) / (nx.score - d.score);
                    return upper_perc * nx.score_bucket + (1.0 - upper_perc) * d.score_bucket;
                } else {
                    let upper_perc = (score - nx.score) / (d.score - nx.score);
                    return upper_perc * d.score_bucket + (1.0 - upper_perc) * nx.score_bucket;
                }
            }
        }
    }
    1.0
}

// ---------------------------------------------------------------------------------------------
// BindingLikelihood

const INVALID_LIKELIHOOD: f64 = -1.0;
const MIN_BUCKET_RANK: f64 = 0.00005;
const PEPTIDE_LENGTHS: usize = 5;

/// BindingLikelihood: likelihood by allele, peptide length and rank bucket (doubling from 0.00005)
#[derive(Clone, Debug)]
pub struct BindingLikelihood {
    rank_buckets: Vec<f64>,
    likelihoods: HashMap<String, Vec<Vec<f64>>>,
    has_data: bool,
}

impl Default for BindingLikelihood {
    fn default() -> Self {
        let mut rank_buckets = Vec::with_capacity(16);
        let mut r = MIN_BUCKET_RANK;
        while r < 1.0 {
            rank_buckets.push(r);
            r *= 2.0;
        }
        BindingLikelihood { rank_buckets, likelihoods: HashMap::new(), has_data: false }
    }
}

/// BindingLikelihood.peptideLengthIndex
fn peptide_length_index(peptide_length: i32) -> Option<usize> {
    if !(MIN_PEPTIDE_LENGTH..=REF_PEPTIDE_LENGTH).contains(&peptide_length) {
        return None;
    }
    Some((peptide_length - MIN_PEPTIDE_LENGTH) as usize)
}

impl BindingLikelihood {
    pub fn has_data(&self) -> bool {
        self.has_data
    }

    /// getBindingLikelihood: by the allele, else its two-digit type, else its gene; interpolated between rank
    /// buckets
    pub fn binding_likelihood(&self, allele: &str, peptide: &str, rank: f64) -> f64 {
        if !self.has_data {
            return INVALID_LIKELIHOOD;
        }
        let Some(pl) = peptide_length_index(peptide.len() as i32) else { return INVALID_LIKELIHOOD };
        // NE13: without a score distribution for the allele and length the rank is INVALID_SCORE (-1), which
        // Java takes as better than every bucket: the highest likelihood. Corrected: no likelihood, as for an
        // unknown allele or length.
        if rank < 0.0 && !compat::java("NE13") {
            return INVALID_LIKELIHOOD;
        }
        if rank >= *self.rank_buckets.last().unwrap() {
            return 0.0;
        }
        let lik = self
            .likelihoods
            .get(allele)
            .or_else(|| allele.get(0..3).and_then(|t| self.likelihoods.get(t)))
            .or_else(|| allele.get(0..1).and_then(|g| self.likelihoods.get(g)));
        let Some(lik) = lik else { return INVALID_LIKELIHOOD };
        for i in 0..self.rank_buckets.len() {
            if rank < self.rank_buckets[i] {
                let likelihood = lik[pl][i];
                if i == 0 {
                    return likelihood;
                }
                let lower = lik[pl][i - 1];
                let (lower_rank, upper_rank) = (self.rank_buckets[i - 1], self.rank_buckets[i]);
                let upper_perc = (rank - lower_rank) / (upper_rank - lower_rank);
                return upper_perc * likelihood + (1.0 - upper_perc) * lower;
            }
        }
        0.0
    }

    /// loadLikelihoods: a matrix per run of rows of one allele
    pub fn load_likelihoods(&mut self, filename: &str) -> Result<bool> {
        let Some((f, lines)) = header_and_lines(filename, "failed to read peptide likelihoods file")? else { return Ok(false) };
        let allele_i = col(&f, "Allele", filename)?;
        let pep_len_i = col(&f, "PeptideLength", filename)?;
        let mut current = String::new();
        let buckets = self.rank_buckets.len();
        for line in &lines {
            let items = split_all(line, ',');
            let allele = item(&items, allele_i, filename)?;
            let peptide_length = int(item(&items, pep_len_i, filename)?)?;
            if current != allele {
                current = allele.to_string();
                self.likelihoods.insert(allele.to_string(), vec![vec![0.0; buckets]; PEPTIDE_LENGTHS]);
            }
            let pl = peptide_length_index(peptide_length).ok_or_else(|| anyhow!("{filename}: ArrayIndexOutOfBoundsException: Index -1"))?;
            let m = self.likelihoods.get_mut(&current).unwrap();
            for (i, v) in items.iter().skip(pep_len_i + 1).enumerate() {
                if i >= buckets {
                    bail!("{filename}: ArrayIndexOutOfBoundsException: Index {i}");
                }
                m[pl][i] = num(v)?;
            }
        }
        info!("loaded {} alleles peptide likelihoods from {}", self.likelihoods.len(), filename);
        self.has_data = true;
        Ok(true)
    }
}

// ---------------------------------------------------------------------------------------------
// ExpressionLikelihood

/// ExpressionLikelihood: presentation likelihood by TPM bucket
#[derive(Clone, Debug, Default)]
pub struct ExpressionLikelihood {
    /// (TPM bucket, likelihood)
    rates: Vec<(f64, f64)>,
}

impl ExpressionLikelihood {
    pub fn has_data(&self) -> bool {
        !self.rates.is_empty()
    }

    /// calcLikelihood: interpolated between buckets, the end buckets' likelihood beyond them
    pub fn calc_likelihood(&self, tpm: f64) -> f64 {
        let r = &self.rates;
        if tpm < r[0].0 {
            return r[0].1;
        }
        if tpm > r[r.len() - 1].0 {
            return r[r.len() - 1].1;
        }
        for i in 0..r.len() {
            let b = r[i];
            if doubles_equal(tpm, b.0) {
                return b.1;
            }
            let next = if i < r.len() - 1 { Some(r[i + 1]) } else { None };
            if let Some(nb) = next {
                if doubles_equal(tpm, nb.0) {
                    return nb.1;
                }
            }
            if tpm > b.0 {
                let Some(nb) = next else { break };
                if tpm < nb.0 {
                    let upper_perc = (tpm - b.0) / (nb.0 - b.0);
                    return upper_perc * nb.1 + (1.0 - upper_perc) * b.1;
                }
            }
        }
        INVALID_CALC
    }

    /// loadTpmRates: the TpmLevel rows
    pub fn load_tpm_rates(&mut self, filename: &str) -> Result<bool> {
        let Some((f, lines)) = header_and_lines(filename, "failed to load expression likelhood rates file")? else { return Ok(false) };
        let (type_i, bucket_i, rate_i) = (col(&f, "DataType", filename)?, col(&f, "Bucket", filename)?, col(&f, "BindingRate", filename)?);
        for line in &lines {
            let items = split_all(line, ',');
            if item(&items, type_i, filename)? != "TpmLevel" {
                continue;
            }
            self.rates.push((num(item(&items, bucket_i, filename)?)?, num(item(&items, rate_i, filename)?)?));
        }
        Ok(true)
    }
}

// ---------------------------------------------------------------------------------------------
// RecognitionSimilarity

/// RecognitionSimilarity: the immunogenic peptides of each allele and length, and the BLOSUM62 distance
/// of a peptide to the nearest of them
pub struct RecognitionSimilarity {
    /// allele -> peptide length -> immunogenic peptides
    peptides: HashMap<String, HashMap<usize, Vec<Vec<u8>>>>,
    blosum: BlosumMapping,
    pub check_self_similarity: bool,
}

impl Default for RecognitionSimilarity {
    fn default() -> Self {
        RecognitionSimilarity { peptides: HashMap::new(), blosum: BlosumMapping::default(), check_self_similarity: false }
    }
}

impl RecognitionSimilarity {
    pub fn has_data(&self) -> bool {
        !self.peptides.is_empty()
    }

    /// whether the recognition data has the allele
    #[allow(dead_code)]
    pub fn has_allele(&self, allele: &str) -> bool {
        self.peptides.contains_key(allele)
    }

    /// the alleles with recognition data
    pub fn alleles(&self) -> impl Iterator<Item = &String> {
        self.peptides.keys()
    }

    /// findTopSimilarity: the lowest summed BLOSUM distance to one of the allele's peptides of the same length
    /// (None where Java returns null)
    pub fn top_similarity(&self, allele: &str, peptide: &str) -> Option<f64> {
        if peptide.contains(AMINO_ACID_21ST) {
            return None;
        }
        let list = self.peptides.get(allele)?.get(&peptide.len())?;
        let p = peptide.as_bytes();
        let mut top = INVALID_CALC;
        let mut found = false;
        for recog in list {
            if self.check_self_similarity && recog.as_slice() == p {
                continue;
            }
            let mut similarity = 0.0;
            let mut skip = false;
            for i in 0..p.len() {
                let (aa1, aa2) = (p[i], recog[i]);
                match (self.blosum.self_mapping(aa1), self.blosum.self_mapping(aa2), self.blosum.map(aa1, aa2)) {
                    (Some(bs1), Some(bs2), Some(map)) => {
                        similarity += (bs1 + bs2) as f64 * 0.5 - map as f64;
                        if top >= 0.0 && similarity >= top {
                            skip = true;
                            break;
                        }
                    }
                    _ => {
                        error!(
                            "invalid index({}) of peptide({}) or recogPeptide({}) for recognition similarity calc",
                            i,
                            peptide,
                            String::from_utf8_lossy(recog)
                        );
                    }
                }
            }
            if skip {
                continue;
            }
            if top < 0.0 || similarity < top {
                top = similarity;
                found = true;
                if top == 0.0 {
                    break;
                }
            }
        }
        if found { Some(top) } else { None }
    }

    /// calcSimilarity
    pub fn calc_similarity(&self, allele: &str, peptide: &str) -> f64 {
        self.top_similarity(allele, peptide).unwrap_or(INVALID_CALC)
    }

    /// calcOtherAlleleSimilarity: the lowest top similarity over the other alleles
    pub fn calc_other_allele_similarity(&self, allele: &str, peptide: &str) -> f64 {
        let mut top = INVALID_CALC;
        for other in self.peptides.keys() {
            if other == allele {
                continue;
            }
            if let Some(s) = self.top_similarity(other, peptide) {
                if top == INVALID_CALC || s < top {
                    top = s;
                }
            }
        }
        top
    }

    /// loadData (RecognitionData.loadRecognitionData): the immunogenic peptides of lengths 8 to 12; false for
    /// a missing file
    pub fn load_data(&mut self, filename: &str) -> Result<bool> {
        if !std::path::Path::new(filename).exists() {
            error!("recognition data file({}) not found", filename);
            return Ok(false);
        }
        let Some((f, lines)) = header_and_lines(filename, "failed to read immunogenic data file")? else { return Ok(false) };
        let allele_i = col(&f, "Allele", filename)?;
        let peptide_i = col(&f, "Peptide", filename)?;
        let immuno_i = f.get("Immunogenic").copied();
        for line in &lines {
            let v = split_all(line, ',');
            let allele = clean_allele(item(&v, allele_i, filename)?);
            let peptide = item(&v, peptide_i, filename)?;
            if !DEFAULT_PEPTIDE_LENGTHS.contains(&(peptide.chars().count() as i32)) {
                continue;
            }
            let immunogenic = match immuno_i {
                Some(i) => {
                    let s = item(&v, i, filename)?;
                    s == "1" || s.eq_ignore_ascii_case("true")
                }
                None => true,
            };
            if !immunogenic {
                continue;
            }
            self.peptides.entry(allele).or_default().entry(peptide.len()).or_default().push(peptide.as_bytes().to_vec());
        }
        info!("loaded {} peptides with immunogenic status from file({})", lines.len(), filename);
        Ok(true)
    }
}

// ---------------------------------------------------------------------------------------------
// BindScorer

/// BindScorer, as NeoScorer uses it: the reference data and calcScoreData
pub struct BindScorer {
    config: ScoreConfig,
    matrices: HashMap<String, HashMap<i32, BindScoreMatrix>>,
    pub flank_scores: FlankScores,
    pub random_distribution: RandomPeptideDistribution,
    pub binding_likelihood: BindingLikelihood,
    pub recognition_similarity: RecognitionSimilarity,
    pub expression_likelihood: ExpressionLikelihood,
}

/// The recognition similarities of one peptide to each allele with recognition data, which
/// calcSimilarity and calcOtherAlleleSimilarity take for every allele the peptide is scored with
pub struct PeptideSimilarities {
    tops: Vec<(String, Option<f64>)>,
}

impl PeptideSimilarities {
    pub fn similarity(&self, allele: &str) -> f64 {
        self.tops.iter().find(|(a, _)| a == allele).and_then(|(_, s)| *s).unwrap_or(INVALID_CALC)
    }
    pub fn other_allele_similarity(&self, allele: &str) -> f64 {
        let mut top = INVALID_CALC;
        for (a, s) in &self.tops {
            if a == allele {
                continue;
            }
            if let Some(s) = *s {
                if top == INVALID_CALC || s < top {
                    top = s;
                }
            }
        }
        top
    }
}

impl BindScorer {
    pub fn new(config: ScoreConfig) -> BindScorer {
        let mut recognition_similarity = RecognitionSimilarity::default();
        recognition_similarity.check_self_similarity = config.check_self_recognition;
        BindScorer {
            config,
            matrices: HashMap::new(),
            flank_scores: FlankScores::default(),
            random_distribution: RandomPeptideDistribution::default(),
            binding_likelihood: BindingLikelihood::default(),
            recognition_similarity,
            expression_likelihood: ExpressionLikelihood::default(),
        }
    }

    pub fn has_matrix(&self, allele: &str, peptide_length: i32) -> bool {
        self.matrices.get(allele).is_some_and(|m| m.contains_key(&peptide_length))
    }

    /// for tests: adds a matrix
    #[cfg(test)]
    pub fn add_matrix(&mut self, m: BindScoreMatrix) {
        self.matrices.entry(m.allele.clone()).or_default().insert(m.peptide_length, m);
    }

    /// loadScoringData: Ok(false) where Java returns false; Err where Java throws (no position weights file
    /// configured: a NullPointerException)
    pub fn load_scoring_data(&mut self) -> Result<bool> {
        let pos_weights = self.config.pos_weights_file.clone().ok_or_else(|| anyhow!("NullPointerException: no position weights file"))?;
        for m in BindScoreMatrix::load_from_csv(&pos_weights)? {
            self.matrices.entry(m.allele.clone()).or_default().insert(m.peptide_length, m);
        }
        let config = self.config.clone();
        if !self.random_distribution.load_data(&config)? {
            return Ok(false);
        }
        if let Some(f) = &config.bind_likelihood_file {
            if !self.binding_likelihood.load_likelihoods(f)? {
                return Ok(false);
            }
        }
        if let Some(f) = &config.flank_pos_weights_file {
            if std::path::Path::new(f).exists() && !self.flank_scores.load_pos_weights(f)? {
                return Ok(false);
            }
        }
        if let Some(f) = &config.recognition_data_file {
            if !self.recognition_similarity.load_data(f)? {
                return Ok(false);
            }
        }
        if let Some(f) = &config.expression_likelihood_file {
            if !self.expression_likelihood.load_tpm_rates(f)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// the similarities of a peptide to every allele's recognition peptides (for calcScoreData's
    /// recognition fields), None without recognition data
    pub fn peptide_similarities(&self, peptide: &str) -> Option<PeptideSimilarities> {
        if !self.recognition_similarity.has_data() {
            return None;
        }
        let tops = self.recognition_similarity.alleles().map(|a| (a.clone(), self.recognition_similarity.top_similarity(a, peptide))).collect();
        Some(PeptideSimilarities { tops })
    }

    /// calcScoreData(bindData): nothing for an allele and length without a matrix
    pub fn calc_score_data(&self, bind_data: &mut BindData, similarities: Option<&PeptideSimilarities>) {
        let Some(matrix) = self.matrices.get(&bind_data.allele).and_then(|m| m.get(&bind_data.peptide_length())) else { return };
        let mut score = matrix.calc_score(&bind_data.peptide);
        let mut flank_score = 0.0;
        if self.flank_scores.has_data() && bind_data.has_flanks() {
            flank_score = self.flank_scores.calc_score(&bind_data.up_flank, &bind_data.down_flank);
            score += flank_score;
        }
        let rank_percentile = self.random_distribution.score_rank(&bind_data.allele, bind_data.peptide_length(), score);
        let mut likelihood = INVALID_CALC;
        let mut likelihood_rank = INVALID_CALC;
        let mut exp_likelihood = INVALID_CALC;
        let mut exp_likelihood_rank = INVALID_CALC;
        let mut recog_similarity = INVALID_CALC;
        let mut other_allele_recog_similarity = INVALID_CALC;
        if self.binding_likelihood.has_data() {
            likelihood = self.binding_likelihood.binding_likelihood(&bind_data.allele, &bind_data.peptide, rank_percentile);
            likelihood_rank = self.random_distribution.likelihood_rank(&bind_data.allele, likelihood);
            // hasTPM() is always true: the TPM starts at 0
            if self.expression_likelihood.has_data() {
                let tpm_likelihood = self.expression_likelihood.calc_likelihood(bind_data.tpm());
                exp_likelihood = likelihood * tpm_likelihood;
                exp_likelihood_rank = self.random_distribution.expression_likelihood_rank(&bind_data.allele, exp_likelihood);
            }
        }
        if self.recognition_similarity.has_data() {
            match similarities {
                Some(s) => {
                    recog_similarity = s.similarity(&bind_data.allele);
                    other_allele_recog_similarity = s.other_allele_similarity(&bind_data.allele);
                }
                None => {
                    recog_similarity = self.recognition_similarity.calc_similarity(&bind_data.allele, &bind_data.peptide);
                    other_allele_recog_similarity =
                        self.recognition_similarity.calc_other_allele_similarity(&bind_data.allele, &bind_data.peptide);
                }
            }
        }
        bind_data.score = score;
        bind_data.flank_score = flank_score;
        bind_data.rank_percentile = rank_percentile;
        bind_data.likelihood = likelihood;
        bind_data.likelihood_rank = likelihood_rank;
        bind_data.exp_likelihood = exp_likelihood;
        bind_data.exp_likelihood_rank = exp_likelihood_rank;
        bind_data.recognition_similarity = recog_similarity;
        bind_data.other_allele_recognition_similarity = other_allele_recog_similarity;
    }
}
