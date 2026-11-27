//! The DOID nodes of the Disease Ontology JSON (`-doid_json`), as DiseaseOntology.readDoidOwlEntryFromDoidJson
//! extracts them: the nodes of the graph with id `http://purl.obolibrary.org/obo/doid.owl`, each with its
//! DOID (the node's id less the `http://purl.obolibrary.org/obo/DOID_` prefix) and term (`lbl`). Used only
//! for OR1's correction.

use anyhow::{anyhow, bail, Result};

use crate::model::OrangeDoidNode;

enum Json {
    Null,
    Bool,
    Num,
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn ws(&mut self) {
        while self.i < self.b.len() && matches!(self.b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }
    fn peek(&mut self) -> Result<u8> {
        self.ws();
        self.b.get(self.i).copied().ok_or_else(|| anyhow!("doid json: unexpected end"))
    }
    fn expect(&mut self, c: u8) -> Result<()> {
        if self.peek()? != c {
            bail!("doid json: expected '{}' at byte {}", c as char, self.i);
        }
        self.i += 1;
        Ok(())
    }
    fn string(&mut self) -> Result<String> {
        self.expect(b'"')?;
        let mut out = String::new();
        loop {
            let c = *self.b.get(self.i).ok_or_else(|| anyhow!("doid json: unterminated string"))?;
            self.i += 1;
            match c {
                b'"' => return Ok(out),
                b'\\' => {
                    let e = self.b[self.i];
                    self.i += 1;
                    match e {
                        b'n' => out.push('\n'),
                        b't' => out.push('\t'),
                        b'r' => out.push('\r'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'u' => {
                            let hex = std::str::from_utf8(&self.b[self.i..self.i + 4])?;
                            self.i += 4;
                            let mut cp = u32::from_str_radix(hex, 16)?;
                            if (0xD800..0xDC00).contains(&cp) && self.b.get(self.i) == Some(&b'\\') && self.b.get(self.i + 1) == Some(&b'u') {
                                let lo = u32::from_str_radix(std::str::from_utf8(&self.b[self.i + 2..self.i + 6])?, 16)?;
                                self.i += 6;
                                cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                            }
                            out.push(char::from_u32(cp).unwrap_or('\u{fffd}'));
                        }
                        other => out.push(other as char),
                    }
                }
                _ => {
                    // copy a run of plain UTF-8 bytes
                    let start = self.i - 1;
                    while self.i < self.b.len() && self.b[self.i] != b'"' && self.b[self.i] != b'\\' {
                        self.i += 1;
                    }
                    out.push_str(std::str::from_utf8(&self.b[start..self.i])?);
                }
            }
        }
    }
    fn value(&mut self) -> Result<Json> {
        match self.peek()? {
            b'{' => {
                self.i += 1;
                let mut members = Vec::new();
                if self.peek()? == b'}' {
                    self.i += 1;
                    return Ok(Json::Obj(members));
                }
                loop {
                    let k = self.string()?;
                    self.expect(b':')?;
                    members.push((k, self.value()?));
                    match self.peek()? {
                        b',' => self.i += 1,
                        b'}' => {
                            self.i += 1;
                            return Ok(Json::Obj(members));
                        }
                        _ => bail!("doid json: bad object at byte {}", self.i),
                    }
                }
            }
            b'[' => {
                self.i += 1;
                let mut items = Vec::new();
                if self.peek()? == b']' {
                    self.i += 1;
                    return Ok(Json::Arr(items));
                }
                loop {
                    items.push(self.value()?);
                    match self.peek()? {
                        b',' => self.i += 1,
                        b']' => {
                            self.i += 1;
                            return Ok(Json::Arr(items));
                        }
                        _ => bail!("doid json: bad array at byte {}", self.i),
                    }
                }
            }
            b'"' => Ok(Json::Str(self.string()?)),
            b't' | b'f' | b'n' => {
                let word = if self.b[self.i..].starts_with(b"true") {
                    4
                } else if self.b[self.i..].starts_with(b"false") {
                    5
                } else if self.b[self.i..].starts_with(b"null") {
                    4
                } else {
                    bail!("doid json: bad literal at byte {}", self.i)
                };
                let null = self.b[self.i] == b'n';
                self.i += word;
                Ok(if null { Json::Null } else { Json::Bool })
            }
            _ => {
                while self.i < self.b.len() && matches!(self.b[self.i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                    self.i += 1;
                }
                Ok(Json::Num)
            }
        }
    }
}

fn get<'a>(o: &'a [(String, Json)], k: &str) -> Option<&'a Json> {
    o.iter().find(|(n, _)| n == k).map(|(_, v)| v)
}

/// The DOID nodes of the file
pub fn read_doid_nodes(path: &str) -> Result<Vec<OrangeDoidNode>> {
    let bytes = std::fs::read(path).map_err(|e| anyhow!("reading DOID JSON {path}: {e}"))?;
    let root = Parser { b: &bytes, i: 0 }.value()?;
    let Json::Obj(root) = root else { bail!("doid json: not an object") };
    let mut out = Vec::new();
    if let Some(Json::Arr(graphs)) = get(&root, "graphs") {
        for g in graphs {
            let Json::Obj(g) = g else { continue };
            if !matches!(get(g, "id"), Some(Json::Str(id)) if id == "http://purl.obolibrary.org/obo/doid.owl") {
                continue;
            }
            if let Some(Json::Arr(nodes)) = get(g, "nodes") {
                for n in nodes {
                    let Json::Obj(n) = n else { continue };
                    let Some(Json::Str(url)) = get(n, "id") else { continue };
                    let term = match get(n, "lbl") {
                        Some(Json::Str(s)) => Some(s.clone()),
                        _ => None,
                    };
                    out.push(OrangeDoidNode { doid: url.replace("http://purl.obolibrary.org/obo/DOID_", ""), doid_term: term });
                }
            }
        }
    }
    Ok(out)
}
