//! A reader of R's serialization format (saveRDS: XDR, versions 2 and 3, gzip-compressed or not), as
//! far as a randomForest object needs it, so that `-chord_model_file <model>.rds` works without R.
//!
//! What is read: NULL, symbols, pairlists (and calls, closures, promises), environments (contents read
//! and dropped), logical, integer, double, complex, string, raw and generic vectors, S4 objects,
//! references, the persistent and special environments and the ALTREP compact sequences and wrappers.
//! Byte code is not read (an error).

use anyhow::{anyhow, bail, Context, Result};
use std::io::Read;

use super::forest::{Forest, Node};

#[derive(Clone, Debug, PartialEq)]
pub enum RObject {
    Null,
    Symbol(String),
    /// a pairlist or call: (tag, value) items
    PairList(Vec<(Option<String>, RObject)>, Attrs),
    Logical(Vec<i32>, Attrs),
    Integer(Vec<i32>, Attrs),
    Real(Vec<f64>, Attrs),
    Complex(Vec<(f64, f64)>, Attrs),
    Str(Vec<Option<String>>, Attrs),
    List(Vec<RObject>, Attrs),
    Raw(Vec<u8>, Attrs),
    /// environments, closures, external pointers, S4 objects and other things not needed
    Other(u8, Attrs),
}

pub type Attrs = Vec<(String, RObject)>;

impl RObject {
    pub fn attrs(&self) -> &[(String, RObject)] {
        match self {
            RObject::PairList(_, a)
            | RObject::Logical(_, a)
            | RObject::Integer(_, a)
            | RObject::Real(_, a)
            | RObject::Complex(_, a)
            | RObject::Str(_, a)
            | RObject::List(_, a)
            | RObject::Raw(_, a)
            | RObject::Other(_, a) => a,
            _ => &[],
        }
    }

    pub fn attr(&self, name: &str) -> Option<&RObject> {
        self.attrs().iter().find(|(n, _)| n == name).map(|(_, v)| v)
    }

    fn set_attrs(&mut self, attrs: Attrs) {
        match self {
            RObject::PairList(_, a)
            | RObject::Logical(_, a)
            | RObject::Integer(_, a)
            | RObject::Real(_, a)
            | RObject::Complex(_, a)
            | RObject::Str(_, a)
            | RObject::List(_, a)
            | RObject::Raw(_, a)
            | RObject::Other(_, a) => *a = attrs,
            _ => {}
        }
    }

    pub fn strings(&self) -> Option<Vec<String>> {
        match self {
            RObject::Str(v, _) => Some(v.iter().map(|s| s.clone().unwrap_or_else(|| "NA".into())).collect()),
            _ => None,
        }
    }

    /// the element `name` of a list with names
    pub fn get(&self, name: &str) -> Option<&RObject> {
        let RObject::List(items, _) = self else { return None };
        let names = self.attr("names")?.strings()?;
        names.iter().position(|n| n == name).and_then(|i| items.get(i))
    }

    /// integer or double values as integers
    pub fn ints(&self) -> Option<Vec<i32>> {
        match self {
            RObject::Integer(v, _) | RObject::Logical(v, _) => Some(v.clone()),
            RObject::Real(v, _) => Some(v.iter().map(|x| *x as i32).collect()),
            _ => None,
        }
    }

    pub fn reals(&self) -> Option<Vec<f64>> {
        match self {
            RObject::Real(v, _) => Some(v.clone()),
            RObject::Integer(v, _) => Some(v.iter().map(|x| *x as f64).collect()),
            _ => None,
        }
    }
}

const NILVALUE_SXP: u8 = 254;
const REFSXP: u8 = 255;

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
    refs: Vec<RObject>,
}

impl Reader<'_> {
    fn bytes(&mut self, n: usize) -> Result<&[u8]> {
        if self.pos + n > self.data.len() {
            bail!("RDS: unexpected end of data");
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    fn int(&mut self) -> Result<i32> {
        Ok(i32::from_be_bytes(self.bytes(4)?.try_into().unwrap()))
    }

    fn double(&mut self) -> Result<f64> {
        Ok(f64::from_be_bytes(self.bytes(8)?.try_into().unwrap()))
    }

    fn length(&mut self) -> Result<usize> {
        let n = self.int()?;
        if n == -1 {
            let hi = self.int()? as u64;
            let lo = self.int()? as u32 as u64;
            return Ok(((hi << 32) | lo) as usize);
        }
        if n < 0 {
            bail!("RDS: negative length");
        }
        Ok(n as usize)
    }

    fn string(&mut self) -> Result<Option<String>> {
        let n = self.int()?;
        if n == -1 {
            return Ok(None);
        }
        Ok(Some(String::from_utf8_lossy(self.bytes(n as usize)?).into_owned()))
    }

    fn item(&mut self) -> Result<RObject> {
        let flags = self.int()? as u32;
        let ty = (flags & 0xff) as u8;
        let has_attr = flags & (1 << 9) != 0;
        let has_tag = flags & (1 << 10) != 0;
        match ty {
            NILVALUE_SXP => Ok(RObject::Null),
            // the special environments and markers
            253 | 252 | 251 | 250 | 242 | 241 => Ok(RObject::Other(ty, Vec::new())),
            REFSXP => {
                let mut index = (flags >> 8) as usize;
                if index == 0 {
                    index = self.int()? as usize;
                }
                self.refs.get(index.wrapping_sub(1)).cloned().ok_or_else(|| anyhow!("RDS: invalid reference {index}"))
            }
            // PERSISTSXP, PACKAGESXP, NAMESPACESXP: a string vector, added to the references
            247..=249 => {
                let _ = self.int()?;
                let n = self.length()?;
                for _ in 0..n {
                    self.item()?;
                }
                let obj = RObject::Other(ty, Vec::new());
                self.refs.push(obj.clone());
                Ok(obj)
            }
            // SYMSXP
            1 => {
                let name = match self.item()? {
                    RObject::Str(v, _) => v.into_iter().next().flatten().unwrap_or_default(),
                    _ => bail!("RDS: a symbol without a name"),
                };
                let obj = RObject::Symbol(name);
                self.refs.push(obj.clone());
                Ok(obj)
            }
            // ENVSXP: locked, enclosure, frame, hash table, attributes
            4 => {
                let _locked = self.int()?;
                self.refs.push(RObject::Other(4, Vec::new()));
                self.item()?;
                self.item()?;
                self.item()?;
                self.item()?;
                Ok(RObject::Other(4, Vec::new()))
            }
            // LISTSXP, CLOSXP, PROMSXP, LANGSXP, DOTSXP, ATTRLANGSXP, ATTRLISTSXP: read along the CDRs
            2 | 3 | 5 | 6 | 17 | 239 | 240 => {
                let mut items = Vec::new();
                let mut attrs = Vec::new();
                let mut first = true;
                let (mut has_attr, mut has_tag) = (has_attr, has_tag);
                loop {
                    let a = if has_attr { Some(self.item()?) } else { None };
                    let tag = if has_tag {
                        match self.item()? {
                            RObject::Symbol(s) => Some(s),
                            RObject::Str(v, _) => v.into_iter().next().flatten(),
                            _ => None,
                        }
                    } else {
                        None
                    };
                    if first {
                        if let Some(a) = a {
                            attrs = pairlist_attrs(a);
                        }
                        first = false;
                    }
                    items.push((tag, self.item()?));
                    // the CDR: another cell of the same list, or its end
                    let next = self.int()? as u32;
                    let next_ty = (next & 0xff) as u8;
                    if matches!(next_ty, 2 | 6 | 17 | 239 | 240) && ty != 3 && ty != 5 {
                        has_attr = next & (1 << 9) != 0;
                        has_tag = next & (1 << 10) != 0;
                        continue;
                    }
                    // anything else ends the list: put it back and read it as an item
                    self.pos -= 4;
                    let rest = self.item()?;
                    if rest != RObject::Null {
                        items.push((None, rest));
                    }
                    break;
                }
                Ok(if matches!(ty, 2 | 6 | 17 | 239 | 240) { RObject::PairList(items, attrs) } else { RObject::Other(ty, attrs) })
            }
            // SPECIALSXP, BUILTINSXP: a name
            7 | 8 => {
                let n = self.int()? as usize;
                self.bytes(n)?;
                Ok(RObject::Other(ty, Vec::new()))
            }
            // CHARSXP (attributes, if any, are read and dropped)
            9 => {
                let s = self.string()?;
                if has_attr {
                    self.item()?;
                }
                Ok(RObject::Str(vec![s], Vec::new()))
            }
            // ALTREP_SXP: class information, state, attributes
            238 => {
                let info = self.item()?;
                let state = self.item()?;
                let attr = self.item()?;
                let class = match &info {
                    RObject::PairList(items, _) => match items.first() {
                        Some((_, RObject::Symbol(s))) => s.clone(),
                        _ => String::new(),
                    },
                    _ => String::new(),
                };
                let mut obj = altrep(&class, state)?;
                obj.set_attrs(pairlist_attrs(attr));
                Ok(obj)
            }
            _ => {
                let mut obj = match ty {
                    // EXTPTRSXP: protected value, tag
                    22 => {
                        self.refs.push(RObject::Other(22, Vec::new()));
                        self.item()?;
                        self.item()?;
                        RObject::Other(22, Vec::new())
                    }
                    // WEAKREFSXP
                    23 => {
                        self.refs.push(RObject::Other(23, Vec::new()));
                        RObject::Other(23, Vec::new())
                    }
                    10 | 13 => {
                        let n = self.length()?;
                        let v = (0..n).map(|_| self.int()).collect::<Result<Vec<i32>>>()?;
                        if ty == 10 { RObject::Logical(v, Vec::new()) } else { RObject::Integer(v, Vec::new()) }
                    }
                    14 => {
                        let n = self.length()?;
                        RObject::Real((0..n).map(|_| self.double()).collect::<Result<Vec<f64>>>()?, Vec::new())
                    }
                    15 => {
                        let n = self.length()?;
                        let v = (0..n).map(|_| Ok((self.double()?, self.double()?))).collect::<Result<Vec<_>>>()?;
                        RObject::Complex(v, Vec::new())
                    }
                    16 => {
                        let n = self.length()?;
                        let mut v = Vec::with_capacity(n);
                        for _ in 0..n {
                            match self.item()? {
                                RObject::Str(s, _) => v.push(s.into_iter().next().flatten()),
                                _ => bail!("RDS: a string vector item is not a string"),
                            }
                        }
                        RObject::Str(v, Vec::new())
                    }
                    19 | 20 => {
                        let n = self.length()?;
                        RObject::List((0..n).map(|_| self.item()).collect::<Result<Vec<_>>>()?, Vec::new())
                    }
                    24 => {
                        let n = self.length()?;
                        RObject::Raw(self.bytes(n)?.to_vec(), Vec::new())
                    }
                    // S4SXP: attributes only
                    25 => RObject::Other(25, Vec::new()),
                    21 => bail!("RDS: byte code is not supported"),
                    _ => bail!("RDS: unsupported item type {ty}"),
                };
                if has_attr {
                    let a = self.item()?;
                    obj.set_attrs(pairlist_attrs(a));
                }
                Ok(obj)
            }
        }
    }
}

fn pairlist_attrs(a: RObject) -> Attrs {
    match a {
        RObject::PairList(items, _) => items.into_iter().map(|(t, v)| (t.unwrap_or_default(), v)).collect(),
        _ => Vec::new(),
    }
}

/// The ALTREP classes R writes in place of ordinary vectors
fn altrep(class: &str, state: RObject) -> Result<RObject> {
    match class {
        "compact_intseq" => {
            let v = state.reals().ok_or_else(|| anyhow!("RDS: invalid compact_intseq"))?;
            let (n, n1, inc) = (v[0] as i64, v[1] as i64, v[2] as i64);
            Ok(RObject::Integer((0..n).map(|i| (n1 + i * inc) as i32).collect(), Vec::new()))
        }
        "compact_realseq" => {
            let v = state.reals().ok_or_else(|| anyhow!("RDS: invalid compact_realseq"))?;
            let (n, n1, inc) = (v[0] as i64, v[1], v[2]);
            Ok(RObject::Real((0..n).map(|i| n1 + i as f64 * inc).collect(), Vec::new()))
        }
        c if c.starts_with("wrap_") => match state {
            RObject::List(mut items, _) if !items.is_empty() => Ok(items.swap_remove(0)),
            _ => bail!("RDS: invalid {c}"),
        },
        "deferred_string" => match state {
            RObject::PairList(items, _) => {
                let arg = items.into_iter().next().map(|(_, v)| v).ok_or_else(|| anyhow!("RDS: invalid deferred_string"))?;
                Ok(match arg {
                    RObject::Integer(v, _) => RObject::Str(v.iter().map(|i| Some(i.to_string())).collect(), Vec::new()),
                    RObject::Real(v, _) => RObject::Str(v.iter().map(|x| Some(super::rfmt::as_character_double(*x))).collect(), Vec::new()),
                    other => other,
                })
            }
            _ => bail!("RDS: invalid deferred_string"),
        },
        other => bail!("RDS: unsupported ALTREP class {other}"),
    }
}

/// Reads an R object from serialized bytes (gzip-compressed or not)
pub fn read_rds_bytes(raw: &[u8]) -> Result<RObject> {
    let data: Vec<u8> = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        hmf_common::gz::GzReader::new(raw).read_to_end(&mut out).context("RDS: inflating")?;
        out
    } else if raw.starts_with(b"BZh") || raw.starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0]) {
        bail!("RDS: bzip2 and xz compression are not supported (save the model with compress = \"gzip\")");
    } else {
        raw.to_vec()
    };
    if !data.starts_with(b"X\n") {
        bail!("RDS: only the XDR format is supported");
    }
    let mut r = Reader { data: &data, pos: 2, refs: Vec::new() };
    let version = r.int()?;
    let _writer = r.int()?;
    let _min = r.int()?;
    match version {
        2 => {}
        3 => {
            let n = r.int()? as usize;
            r.bytes(n)?;
        }
        v => bail!("RDS: unsupported serialization version {v}"),
    }
    r.item()
}

#[cfg_attr(not(test), allow(dead_code))]
pub fn read_rds(path: &str) -> Result<RObject> {
    let raw = std::fs::read(path).with_context(|| format!("reading {path}"))?;
    read_rds_bytes(&raw).with_context(|| format!("reading R object {path}"))
}

/// The forest of a randomForest object, as predict.randomForest uses it
pub fn forest_from_object(obj: &RObject) -> Result<Forest> {
    let class = obj.attr("class").and_then(|c| c.strings()).unwrap_or_default();
    if !class.iter().any(|c| c == "randomForest") {
        bail!("object not of class randomForest");
    }
    let field = |name: &str| obj.get(name).ok_or_else(|| anyhow!("randomForest object without {name}"));
    let ty = field("type")?.strings().unwrap_or_default();
    if ty.first().map(String::as_str) != Some("classification") {
        bail!("only classification forests are supported");
    }
    let classes = field("classes")?.strings().ok_or_else(|| anyhow!("invalid classes"))?;
    let importance = field("importance")?;
    let variables = match importance.attr("dim") {
        Some(_) => match importance.attr("dimnames") {
            Some(RObject::List(dn, _)) => dn.first().and_then(|r| r.strings()),
            _ => None,
        },
        None => importance.attr("names").and_then(|n| n.strings()),
    }
    .ok_or_else(|| anyhow!("no variable names in the importance"))?;
    let f = field("forest")?;
    let ff = |name: &str| f.get(name).ok_or_else(|| anyhow!("forest without {name}"));
    let ints = |name: &str| ff(name).and_then(|o| o.ints().ok_or_else(|| anyhow!("forest {name} is not numeric")));
    let ncat = ints("ncat")?;
    if ncat.iter().any(|c| *c != 1) {
        bail!("categorical predictors are not supported");
    }
    let ntree = ints("ntree")?[0] as usize;
    let nrnodes = ints("nrnodes")?[0] as usize;
    let ndbigtree = ints("ndbigtree")?;
    let nodestatus = ints("nodestatus")?;
    let bestvar = ints("bestvar")?;
    let treemap = ints("treemap")?;
    let nodepred = ints("nodepred")?;
    let xbestsplit = ff("xbestsplit")?.reals().ok_or_else(|| anyhow!("forest xbestsplit is not numeric"))?;
    if ndbigtree.len() != ntree || nodestatus.len() != nrnodes * ntree || treemap.len() != 2 * nrnodes * ntree || xbestsplit.len() != nrnodes * ntree {
        bail!("forest arrays of unexpected sizes");
    }
    let trees = (0..ntree)
        .map(|t| {
            (0..ndbigtree[t] as usize)
                .map(|k| {
                    let i = k + nrnodes * t;
                    Node {
                        status: nodestatus[i],
                        var: bestvar[i],
                        left: treemap[k + nrnodes * 2 * t],
                        right: treemap[k + nrnodes * (2 * t + 1)],
                        pred: nodepred[i],
                        split: xbestsplit[i],
                    }
                })
                .collect()
        })
        .collect();
    let forest = Forest { classes, variables, trees };
    forest.check()?;
    Ok(forest)
}
