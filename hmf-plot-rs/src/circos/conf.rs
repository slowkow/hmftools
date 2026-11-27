//! Circos configuration files: `key = value` lines in nested `<block>` ... `</block>` sections, `#`
//! comments, and `<<include file>>` lines. The includes of Circos's own `etc/` files (colours, fonts,
//! image, housekeeping) are not read: their settings are built into the engine (Circos 0.69 defaults).

use std::path::Path;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Block {
    pub params: Vec<(String, String)>,
    pub blocks: Vec<(String, Block)>,
    /// the `<<include ...>>` files named in this block
    pub includes: Vec<String>,
}

impl Block {
    /// The last value given for `key` (Circos keeps the last of repeated keys).
    pub fn get(&self, key: &str) -> Option<&str> {
        self.params.iter().rev().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }

    pub fn has(&self, key: &str) -> bool {
        self.get(key).is_some()
    }

    /// Circos's truth: "yes", "y", "1" (case-insensitive).
    pub fn flag(&self, key: &str) -> Option<bool> {
        self.get(key).map(is_true)
    }

    pub fn block(&self, name: &str) -> Option<&Block> {
        self.blocks.iter().find(|(n, _)| n == name).map(|(_, b)| b)
    }

    /// Every block named `name`, in file order.
    pub fn blocks<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Block> + 'a {
        self.blocks.iter().filter(move |(n, _)| n == name).map(|(_, b)| b)
    }
}

pub fn is_true(v: &str) -> bool {
    let v = v.trim().to_ascii_lowercase();
    v == "yes" || v == "y" || v == "1" || v == "true"
}

pub fn parse(text: &str) -> Result<Block, String> {
    let mut stack: Vec<(String, Block)> = vec![(String::new(), Block::default())];
    for (n, raw) in text.lines().enumerate() {
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        if let Some(inc) = line.strip_prefix("<<include").and_then(|s| s.strip_suffix(">>")) {
            stack.last_mut().unwrap().1.includes.push(inc.trim().to_string());
            continue;
        }
        if let Some(name) = line.strip_prefix("</").and_then(|s| s.strip_suffix('>')) {
            let (open, block) = stack.pop().ok_or("unbalanced blocks")?;
            if open != name.trim() || stack.is_empty() {
                return Err(format!("line {}: </{}> closes <{}>", n + 1, name.trim(), open));
            }
            stack.last_mut().unwrap().1.blocks.push((open, block));
            continue;
        }
        if let Some(name) = line.strip_prefix('<').and_then(|s| s.strip_suffix('>')) {
            stack.push((name.trim().to_string(), Block::default()));
            continue;
        }
        if let Some((k, v)) = line.split_once('=') {
            stack.last_mut().unwrap().1.params.push((k.trim().to_string(), v.trim().to_string()));
        } else {
            return Err(format!("line {}: cannot parse '{}'", n + 1, line));
        }
    }
    if stack.len() != 1 {
        return Err(format!("block <{}> is not closed", stack.last().unwrap().0));
    }
    Ok(stack.pop().unwrap().1)
}

/// Removes a `#` comment (Circos treats `#` as a comment start unless escaped with `\`).
fn strip_comment(line: &str) -> &str {
    let b = line.as_bytes();
    for i in 0..b.len() {
        if b[i] == b'#' && (i == 0 || b[i - 1] != b'\\') {
            return &line[..i];
        }
    }
    line
}

pub fn read(path: &Path) -> Result<Block, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&text)
}

/// A value with an optional unit suffix: "0.90r", "50p", "1500p", "1u", "0.002r", "2".
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Unit {
    /// relative (to a radius, or to the genome size for ideogram spacing)
    R(f64),
    /// pixels
    P(f64),
    /// chromosome units
    U(f64),
    /// bases
    B(f64),
    /// no unit
    N(f64),
}

impl Unit {
    pub fn parse(s: &str) -> Option<Unit> {
        let s = s.trim();
        let (num, kind) = match s.chars().last()? {
            'r' => (&s[..s.len() - 1], 'r'),
            'p' => (&s[..s.len() - 1], 'p'),
            'u' => (&s[..s.len() - 1], 'u'),
            'b' => (&s[..s.len() - 1], 'b'),
            _ => (s, 'n'),
        };
        let v: f64 = num.trim().parse().ok()?;
        Some(match kind {
            'r' => Unit::R(v),
            'p' => Unit::P(v),
            'u' => Unit::U(v),
            'b' => Unit::B(v),
            _ => Unit::N(v),
        })
    }

    /// The number without its unit (Circos's unit_strip).
    pub fn value(&self) -> f64 {
        match *self {
            Unit::R(v) | Unit::P(v) | Unit::U(v) | Unit::B(v) | Unit::N(v) => v,
        }
    }
}

/// Circos's data-file options field: "a=1,color=224,7,20,b=x" (a fragment without '=' continues the
/// previous value, so colours keep their commas).
pub fn parse_options(s: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    for part in s.split([',', ';']) {
        if let Some((k, v)) = part.split_once('=') {
            out.push((k.trim().to_string(), v.trim().to_string()));
        } else if let Some(last) = out.last_mut() {
            last.1.push(',');
            last.1.push_str(part.trim());
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blocks_params_includes() {
        let b = parse(
            "<<include etc/colors.conf>>\nshow_ticks = no # comment\n<ideogram>\n <spacing>\n default = 0.002r\n </spacing>\n radius = 0.90r\n</ideogram>\n<plots>\n<plot>\ntype = scatter\n</plot>\n<plot>\ntype = text\n</plot>\n</plots>\n",
        )
        .unwrap();
        assert_eq!(b.includes, ["etc/colors.conf"]);
        assert_eq!(b.get("show_ticks"), Some("no"));
        assert_eq!(b.block("ideogram").unwrap().block("spacing").unwrap().get("default"), Some("0.002r"));
        let types: Vec<_> = b.block("plots").unwrap().blocks("plot").map(|p| p.get("type").unwrap()).collect();
        assert_eq!(types, ["scatter", "text"]);
        assert!(parse("<a>\n</b>\n").is_err());
    }

    #[test]
    fn units_and_options() {
        assert_eq!(Unit::parse("0.90r"), Some(Unit::R(0.9)));
        assert_eq!(Unit::parse("50p"), Some(Unit::P(50.0)));
        assert_eq!(Unit::parse("8"), Some(Unit::N(8.0)));
        let o = parse_options("color=224,7,20,glyph=circle,glyph_size=10");
        assert_eq!(o, [("color".into(), "224,7,20".into()), ("glyph".into(), "circle".into()), ("glyph_size".into(), "10".into())]);
        assert_eq!(parse_options("r0=0.55r"), [("r0".to_string(), "0.55r".to_string())]);
    }
}
