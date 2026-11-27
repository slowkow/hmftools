//! The CUPPA classifier as exported once from pycuppa's pickle by `scripts/export_model.py` (format in that
//! script), embedded in the binary (`resources/cuppa_model.38.bin.gz`), or read from a file given as
//! `--classifier_path`.
//!
//! Every parameter of the prediction is held with its exact IEEE value; float32 arrays of the pickle
//! (the alt-SJ classifier) are widened to f64, which is exact, as numpy widens them when they meet float64
//! data.

use std::collections::HashMap;
use std::io::Read;

use anyhow::{anyhow, bail, Context, Result};

pub const EMBEDDED_MODEL: &[u8] = include_bytes!("../../resources/cuppa_model.38.bin.gz");
const MAGIC: &[u8] = b"CUPPARSMODEL1\n";

#[derive(Debug)]
enum Entry {
    Floats { dims: Vec<usize>, data: Vec<f64>, f32: bool },
    Ints { data: Vec<i64> },
    Strings(Vec<String>),
    Text(String),
    Alias(String),
}

/// The entries of an exported model, by name
pub struct Resource {
    entries: HashMap<String, Entry>,
}

struct Cursor<'a> {
    d: &'a [u8],
    p: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.p + n > self.d.len() {
            bail!("model resource truncated");
        }
        let s = &self.d[self.p..self.p + n];
        self.p += n;
        Ok(s)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn string(&mut self) -> Result<String> {
        let n = self.u32()? as usize;
        Ok(std::str::from_utf8(self.take(n)?).context("model resource: invalid UTF-8")?.to_string())
    }
}

/// undoes the exporter's byte shuffle: all first bytes of the elements, then all second bytes, ...
fn unshuffle(raw: &[u8], count: usize, size: usize) -> Vec<u8> {
    let mut out = vec![0u8; count * size];
    for b in 0..size {
        let plane = &raw[b * count..(b + 1) * count];
        for (i, v) in plane.iter().enumerate() {
            out[i * size + b] = *v;
        }
    }
    out
}

impl Resource {
    /// Parses a gzip-compressed (or plain) exported model
    pub fn parse(bytes: &[u8]) -> Result<Resource> {
        let body: Vec<u8> = if bytes.starts_with(&[0x1f, 0x8b]) {
            let mut out = Vec::new();
            flate2::read::MultiGzDecoder::new(bytes).read_to_end(&mut out).context("model resource: gzip")?;
            out
        } else {
            bytes.to_vec()
        };
        if !body.starts_with(MAGIC) {
            bail!("not an exported CUPPA model (see scripts/export_model.py)");
        }
        let mut c = Cursor { d: &body, p: MAGIC.len() };
        let n = c.u32()? as usize;
        let mut entries = HashMap::with_capacity(n);
        for _ in 0..n {
            let name_len = c.u16()? as usize;
            let name = std::str::from_utf8(c.take(name_len)?)?.to_string();
            let kind = c.u8()?;
            let ndim = c.u8()? as usize;
            let mut dims = Vec::with_capacity(ndim);
            for _ in 0..ndim {
                dims.push(c.u64()? as usize);
            }
            let count: usize = dims.iter().product();
            let entry = match kind {
                1 => {
                    let raw = unshuffle(c.take(count * 8)?, count, 8);
                    let data = raw.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
                    Entry::Floats { dims, data, f32: false }
                }
                2 => {
                    let raw = unshuffle(c.take(count * 4)?, count, 4);
                    let data = raw.chunks_exact(4).map(|b| f32::from_le_bytes(b.try_into().unwrap()) as f64).collect();
                    Entry::Floats { dims, data, f32: true }
                }
                3 => {
                    let raw = c.take(count * 8)?;
                    Entry::Ints { data: raw.chunks_exact(8).map(|b| i64::from_le_bytes(b.try_into().unwrap())).collect() }
                }
                4 => {
                    let mut v = Vec::with_capacity(count);
                    for _ in 0..count {
                        v.push(c.string()?);
                    }
                    Entry::Strings(v)
                }
                5 => Entry::Text(c.string()?),
                6 => Entry::Alias(c.string()?),
                k => bail!("model resource: unknown entry kind {k} ({name})"),
            };
            entries.insert(name, entry);
        }
        Ok(Resource { entries })
    }

    fn get(&self, name: &str) -> Result<&Entry> {
        let mut e = self.entries.get(name).ok_or_else(|| anyhow!("model resource: missing entry {name}"))?;
        if let Entry::Alias(target) = e {
            e = self.entries.get(target).ok_or_else(|| anyhow!("model resource: missing entry {target}"))?;
        }
        Ok(e)
    }

    pub fn floats(&self, name: &str) -> Result<(&[usize], &[f64], bool)> {
        match self.get(name)? {
            Entry::Floats { dims, data, f32 } => Ok((dims, data, *f32)),
            _ => bail!("model resource: {name} is not a float array"),
        }
    }

    pub fn vector(&self, name: &str) -> Result<Vec<f64>> {
        let (dims, data, _) = self.floats(name)?;
        if dims.len() != 1 {
            bail!("model resource: {name} is not a vector");
        }
        Ok(data.to_vec())
    }

    pub fn scalar(&self, name: &str) -> Result<f64> {
        let v = self.vector(name)?;
        if v.len() != 1 {
            bail!("model resource: {name} is not a scalar");
        }
        Ok(v[0])
    }

    /// a matrix of `rows` x `cols`, row-major
    pub fn matrix(&self, name: &str, rows: usize, cols: usize) -> Result<Vec<f64>> {
        let (dims, data, _) = self.floats(name)?;
        if dims != [rows, cols] {
            bail!("model resource: {name} has shape {dims:?}, expected [{rows}, {cols}]");
        }
        Ok(data.to_vec())
    }

    pub fn ints(&self, name: &str) -> Result<Vec<i64>> {
        match self.get(name)? {
            Entry::Ints { data, .. } => Ok(data.clone()),
            _ => bail!("model resource: {name} is not an int array"),
        }
    }

    pub fn int(&self, name: &str) -> Result<i64> {
        let v = self.ints(name)?;
        v.first().copied().ok_or_else(|| anyhow!("model resource: {name} is empty"))
    }

    pub fn strings(&self, name: &str) -> Result<Vec<String>> {
        match self.get(name)? {
            Entry::Strings(v) => Ok(v.clone()),
            _ => bail!("model resource: {name} is not a string list"),
        }
    }

    pub fn text(&self, name: &str) -> Result<String> {
        match self.get(name)? {
            Entry::Text(t) => Ok(t.clone()),
            _ => bail!("model resource: {name} is not a text"),
        }
    }
}

/// A fitted multinomial logistic regression (sklearn LogisticRegression, pycuppa's subclass)
pub struct LogReg {
    pub classes: Vec<String>,
    pub features: Vec<String>,
    /// classes x features, row-major (sklearn's coef_)
    pub coef: Vec<f64>,
    pub intercept: Vec<f64>,
    /// the training means of the features (feat_contrib)
    pub x_means: Option<Vec<f64>>,
    /// stored as float32 in the pickle (the alt-SJ classifier): numpy widens a row-major copy, and
    /// OpenBLAS multiplies it with dgemv_n
    pub f32: bool,
}

impl LogReg {
    fn load(r: &Resource, prefix: &str, with_means: bool) -> Result<LogReg> {
        let classes = r.strings(&format!("{prefix}.classes"))?;
        let features = r.strings(&format!("{prefix}.features"))?;
        let coef = r.matrix(&format!("{prefix}.coef"), classes.len(), features.len())?;
        let f32 = r.floats(&format!("{prefix}.coef"))?.2;
        let intercept = r.vector(&format!("{prefix}.intercept"))?;
        if intercept.len() != classes.len() {
            bail!("model resource: {prefix}.intercept length");
        }
        let x_means = if with_means { Some(r.vector(&format!("{prefix}.x_means"))?) } else { None };
        Ok(LogReg { classes, features, coef, intercept, x_means, f32 })
    }

    pub fn n_classes(&self) -> usize {
        self.classes.len()
    }
    pub fn n_features(&self) -> usize {
        self.features.len()
    }
    pub fn coef_row(&self, class: usize) -> &[f64] {
        let n = self.n_features();
        &self.coef[class * n..(class + 1) * n]
    }
}

/// Per-class profiles of a ProfileSimilarityTransformer
pub struct Profiles {
    pub features: Vec<String>,
    /// the profiles' column names (`<prefix><class>`)
    pub classes: Vec<String>,
    /// class by class, each profile contiguous (n_classes x n_features)
    pub values: Vec<f64>,
    /// the L2 norm of each class's profile, as numpy computed it
    pub norms: Vec<f64>,
}

impl Profiles {
    fn load(r: &Resource, prefix: &str) -> Result<Profiles> {
        let features = r.strings(&format!("{prefix}.profile_features"))?;
        let classes = r.strings(&format!("{prefix}.profile_classes"))?;
        let values = r.matrix(&format!("{prefix}.profiles"), classes.len(), features.len())?;
        let norms = r.vector(&format!("{prefix}.profile_norms"))?;
        if norms.len() != classes.len() {
            bail!("model resource: {prefix}.profile_norms length");
        }
        Ok(Profiles { features, classes, values, norms })
    }
    pub fn class_profile(&self, c: usize) -> &[f64] {
        let n = self.features.len();
        &self.values[c * n..(c + 1) * n]
    }
}

pub struct GenPosClf {
    pub noise_profile: Vec<f64>,
    pub noise_counts: i64,
    pub profiles: Profiles,
    pub exponent: i64,
    pub lr: LogReg,
}

pub struct Snv96Clf {
    pub noise_profile: Vec<f64>,
    pub noise_counts: i64,
    pub lr: LogReg,
}

pub struct EventClf {
    /// Chi2FeatureSelector's selected features of the driver, fusion, virus and trait features
    pub other_selected: Vec<String>,
    /// the SV and TMB features, in the order of the MaxScaler's maxima
    pub svtmb_features: Vec<String>,
    pub svtmb_max: Vec<f64>,
    pub lr: LogReg,
}

pub struct RnaClf {
    pub selected: Vec<String>,
    pub max: Vec<f64>,
    pub profiles: Profiles,
    pub exponent: i64,
    pub lr: LogReg,
}

pub enum Calibrator {
    Passthrough,
    /// IsotonicRegression(out_of_bounds="clip"): the thresholds
    Isotonic { x: Vec<f64>, y: Vec<f64> },
}

pub struct FusionOverride {
    pub feat_prefix: String,
    pub feat_basename: String,
    pub target_class: String,
}

pub struct MetaClf {
    pub lr: LogReg,
    /// per class of the logistic regression, in its class order
    pub calibrators: Vec<Calibrator>,
    pub overrides: Vec<FusionOverride>,
    pub mask_base_value: f64,
    pub male_keywords: Vec<String>,
    pub female_keywords: Vec<String>,
}

pub struct SigQuantiles {
    pub features: Vec<String>,
    pub classes: Vec<String>,
    /// per class: quantiles (n_quantiles x features, row-major) and references
    pub quantiles: Vec<(usize, Vec<f64>, Vec<f64>)>,
}

#[allow(dead_code)] // the whole table as exported; the vis data uses n_total, recall and precision
pub struct CvRow {
    pub class: String,
    pub clf_name: String,
    pub n_total: i64,
    pub n_predicted: i64,
    pub n_correct: i64,
    pub recall: f64,
    pub precision: f64,
}

pub struct Model {
    pub source_sha256: String,
    /// the source pickle's size and last 8 bytes (hex), by which `--classifier_path` recognises it
    pub source_size: u64,
    pub source_tail: String,
    pub source_versions: String,
    pub classes: Vec<String>,
    pub req_gen_pos: Vec<String>,
    pub req_snv96: Vec<String>,
    pub req_event: Vec<String>,
    pub req_sig: Vec<String>,
    pub req_gene_exp: Vec<String>,
    pub req_alt_sj: Vec<String>,
    pub gen_pos: GenPosClf,
    pub snv96: Snv96Clf,
    pub event: EventClf,
    pub gene_exp: RnaClf,
    pub alt_sj: RnaClf,
    pub dna_combined: MetaClf,
    pub rna_combined: MetaClf,
    pub prob_floor: f64,
    pub sig: SigQuantiles,
    pub cv: Vec<CvRow>,
}

fn load_meta(r: &Resource, prefix: &str, with_overrides: bool) -> Result<MetaClf> {
    let lr = LogReg::load(r, &format!("{prefix}.lr"), false)?;
    let kinds = r.strings(&format!("{prefix}.cal_kinds"))?;
    if kinds.len() != lr.n_classes() {
        bail!("model resource: {prefix}.cal_kinds length");
    }
    let mut calibrators = Vec::with_capacity(kinds.len());
    for (c, kind) in lr.classes.iter().zip(&kinds) {
        calibrators.push(match kind.as_str() {
            "passthrough" => Calibrator::Passthrough,
            "isotonic" => {
                let x = r.vector(&format!("{prefix}.cal.{c}.x"))?;
                let y = r.vector(&format!("{prefix}.cal.{c}.y"))?;
                if x.len() != y.len() || x.len() < 2 {
                    bail!("model resource: {prefix} calibrator of {c}");
                }
                Calibrator::Isotonic { x, y }
            }
            k => bail!("model resource: calibrator kind {k}"),
        });
    }
    let (overrides, mask_base_value) = if with_overrides {
        let p = r.strings(&format!("{prefix}.overrides.feat_prefix"))?;
        let b = r.strings(&format!("{prefix}.overrides.feat_basename"))?;
        let t = r.strings(&format!("{prefix}.overrides.target_class"))?;
        let o = p
            .into_iter()
            .zip(b)
            .zip(t)
            .map(|((feat_prefix, feat_basename), target_class)| FusionOverride { feat_prefix, feat_basename, target_class })
            .collect();
        (o, r.scalar(&format!("{prefix}.mask_base_value"))?)
    } else {
        (Vec::new(), 0.0)
    };
    Ok(MetaClf {
        lr,
        calibrators,
        overrides,
        mask_base_value,
        male_keywords: r.strings(&format!("{prefix}.male_keywords"))?,
        female_keywords: r.strings(&format!("{prefix}.female_keywords"))?,
    })
}

fn load_rna(r: &Resource, prefix: &str) -> Result<RnaClf> {
    let selected = r.strings(&format!("{prefix}.selected"))?;
    let max = r.vector(&format!("{prefix}.max"))?;
    if max.len() != selected.len() {
        bail!("model resource: {prefix}.max length");
    }
    Ok(RnaClf {
        selected,
        max,
        profiles: Profiles::load(r, prefix)?,
        exponent: r.int(&format!("{prefix}.exponent"))?,
        lr: LogReg::load(r, &format!("{prefix}.lr"), false)?,
    })
}

impl Model {
    pub fn load(r: &Resource) -> Result<Model> {
        let gen_pos = GenPosClf {
            noise_profile: r.vector("gen_pos.noise_profile")?,
            noise_counts: r.int("gen_pos.noise_counts")?,
            profiles: Profiles::load(r, "gen_pos")?,
            exponent: r.int("gen_pos.exponent")?,
            lr: LogReg::load(r, "gen_pos.lr", false)?,
        };
        let snv96 = Snv96Clf {
            noise_profile: r.vector("snv96.noise_profile")?,
            noise_counts: r.int("snv96.noise_counts")?,
            lr: LogReg::load(r, "snv96.lr", false)?,
        };
        let event = EventClf {
            other_selected: r.strings("event.other_selected")?,
            svtmb_features: r.strings("event.svtmb_features")?,
            svtmb_max: r.vector("event.svtmb_max")?,
            lr: LogReg::load(r, "event.lr", true)?,
        };
        let sig_features = r.strings("sig.features")?;
        let sig_classes = r.strings("sig.classes")?;
        let mut quantiles = Vec::with_capacity(sig_classes.len());
        for c in &sig_classes {
            let (dims, q, _) = r.floats(&format!("sig.{c}.quantiles"))?;
            if dims.len() != 2 || dims[1] != sig_features.len() {
                bail!("model resource: sig.{c}.quantiles shape");
            }
            let refs = r.vector(&format!("sig.{c}.references"))?;
            if refs.len() != dims[0] {
                bail!("model resource: sig.{c}.references length");
            }
            quantiles.push((dims[0], q.to_vec(), refs));
        }
        let cv_class = r.strings("cv.class")?;
        let cv_clf = r.strings("cv.clf_name")?;
        let (nt, np_, nc) = (r.ints("cv.n_total")?, r.ints("cv.n_predicted")?, r.ints("cv.n_correct")?);
        let (rec, prec) = (r.vector("cv.recall")?, r.vector("cv.precision")?);
        let cv = (0..cv_class.len())
            .map(|i| CvRow {
                class: cv_class[i].clone(),
                clf_name: cv_clf[i].clone(),
                n_total: nt[i],
                n_predicted: np_[i],
                n_correct: nc[i],
                recall: rec[i],
                precision: prec[i],
            })
            .collect();
        let m = Model {
            source_sha256: r.text("source.sha256")?,
            source_size: r.text("source.size")?.parse()?,
            source_tail: r.text("source.tail")?,
            source_versions: r.text("source.versions")?,
            classes: r.strings("classes")?,
            req_gen_pos: r.strings("required.gen_pos")?,
            req_snv96: r.strings("required.snv96")?,
            req_event: r.strings("required.event")?,
            req_sig: r.strings("required.sig")?,
            req_gene_exp: r.strings("required.gene_exp")?,
            req_alt_sj: r.strings("required.alt_sj")?,
            gen_pos,
            snv96,
            event,
            gene_exp: load_rna(r, "gene_exp")?,
            alt_sj: load_rna(r, "alt_sj")?,
            dna_combined: load_meta(r, "dna_combined", true)?,
            rna_combined: load_meta(r, "rna_combined", false)?,
            prob_floor: r.scalar("combiner.prob_floor")?,
            sig: SigQuantiles { features: sig_features, classes: sig_classes, quantiles },
            cv,
        };
        m.check()?;
        Ok(m)
    }

    /// consistency of the parts, as sklearn's feature-name checks would find it
    fn check(&self) -> Result<()> {
        if self.gen_pos.profiles.features != self.req_gen_pos {
            bail!("model: gen_pos profile features differ from the required features");
        }
        if self.gen_pos.lr.features != self.gen_pos.profiles.classes {
            bail!("model: gen_pos logistic regression features");
        }
        if self.snv96.lr.features != self.req_snv96 || self.snv96.noise_profile.len() != self.req_snv96.len() {
            bail!("model: snv96 features");
        }
        let mut ev: Vec<String> = self.event.other_selected.clone();
        ev.extend(self.event.svtmb_features.iter().cloned());
        if ev != self.event.lr.features {
            bail!("model: event features");
        }
        for (name, clf) in [("gene_exp", &self.gene_exp), ("alt_sj", &self.alt_sj)] {
            if clf.lr.features != clf.profiles.classes || clf.selected != clf.profiles.features {
                bail!("model: {name} features");
            }
        }
        if self.dna_combined.lr.classes != self.classes {
            bail!("model: dna_combined classes");
        }
        Ok(())
    }
}

/// Whether a file is an exported model (its first bytes, gzip-decompressed if need be, are the magic)
fn is_exported(path: &str) -> Result<bool> {
    let f = std::fs::File::open(path).with_context(|| format!("cannot read classifier file {path}"))?;
    let mut head = [0u8; 2];
    let mut r = std::io::BufReader::new(f);
    {
        use std::io::BufRead;
        let b = r.fill_buf()?;
        head[..b.len().min(2)].copy_from_slice(&b[..b.len().min(2)]);
    }
    let mut start = vec![0u8; MAGIC.len()];
    let n = if head == [0x1f, 0x8b] {
        let mut d = flate2::bufread::MultiGzDecoder::new(r);
        read_up_to(&mut d, &mut start)?
    } else {
        read_up_to(&mut r, &mut start)?
    };
    Ok(n == MAGIC.len() && start == MAGIC)
}

fn read_up_to(r: &mut dyn Read, buf: &mut [u8]) -> Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(_) => break,
        }
    }
    Ok(n)
}

/// The model to use for `--classifier_path`: an exported model file, or the pickle the embedded model was
/// exported from (recognised by its size and last 8 bytes, which for a gzip file are the CRC-32 and the
/// length of its content). Another pickle cannot be read without Python: export it with
/// scripts/export_model.py and pass the exported file.
pub fn load_for_path(path: &str) -> Result<Model> {
    if is_exported(path)? {
        let bytes = std::fs::read(path).with_context(|| format!("cannot read classifier file {path}"))?;
        return Model::load(&Resource::parse(&bytes)?);
    }
    let model = Model::load(&Resource::parse(EMBEDDED_MODEL)?)?;
    let meta = std::fs::metadata(path).with_context(|| format!("cannot read classifier file {path}"))?;
    let mut tail = [0u8; 8];
    if meta.len() >= 8 {
        use std::io::{Seek, SeekFrom};
        let mut f = std::fs::File::open(path)?;
        f.seek(SeekFrom::End(-8))?;
        f.read_exact(&mut tail)?;
    }
    let tail_hex: String = tail.iter().map(|b| format!("{b:02x}")).collect();
    if meta.len() != model.source_size || tail_hex != model.source_tail {
        let digest = std::fs::read(path).map(|b| crate::sha256::hex_digest(&b)).unwrap_or_default();
        bail!(
            "classifier {path} (sha256 {digest}) is not the pickle the embedded model was exported from (sha256 {}); export it with scripts/export_model.py and pass the exported file",
            model.source_sha256
        );
    }
    super::pylog::debug("CuppaClassifier", "from_file", &format!("Using the embedded model, exported with {}", model.source_versions));
    Ok(model)
}
