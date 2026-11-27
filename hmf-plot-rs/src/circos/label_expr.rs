//! String expressions of Circos's `label_format` (Perl evaluated by Circos::Expression): enough for
//! `eval( var(chr) =~ /hs([6]|[9])$/ ? var(label) . "." : var(label) )`, i.e. a regex match, the ternary
//! operator, string concatenation with `.`, `var(x)` and quoted strings. (Added for LINX's ideogram labels.)

/// Evaluates `expr` with `var` resolving `var(name)`; None when the expression is not of this form.
pub fn eval(expr: &str, var: &dyn Fn(&str) -> Option<String>) -> Option<String> {
    let s = expr.trim();
    let s = s.strip_prefix("eval(").and_then(|x| x.strip_suffix(')')).unwrap_or(s).trim();
    let mut p = P { s: s.as_bytes(), i: 0, var };
    let v = p.ternary()?;
    p.ws();
    (p.i == p.s.len()).then_some(v)
}

struct P<'a> {
    s: &'a [u8],
    i: usize,
    var: &'a dyn Fn(&str) -> Option<String>,
}

impl P<'_> {
    fn ws(&mut self) {
        while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }

    fn eat(&mut self, t: &str) -> bool {
        self.ws();
        if self.s[self.i..].starts_with(t.as_bytes()) {
            self.i += t.len();
            true
        } else {
            false
        }
    }

    fn ternary(&mut self) -> Option<String> {
        let c = self.concat()?;
        if self.eat("=~") {
            self.ws();
            if !self.eat("/") {
                return None;
            }
            let start = self.i;
            while self.i < self.s.len() && self.s[self.i] != b'/' {
                if self.s[self.i] == b'\\' {
                    self.i += 1;
                }
                self.i += 1;
            }
            let re = std::str::from_utf8(&self.s[start..self.i.min(self.s.len())]).ok()?.to_string();
            self.i += 1;
            let matched = Regex::parse(&re)?.is_match(&c);
            if !self.eat("?") {
                return Some(if matched { "1".into() } else { String::new() });
            }
            let a = self.ternary()?;
            if !self.eat(":") {
                return None;
            }
            let b = self.ternary()?;
            return Some(if matched { a } else { b });
        }
        Some(c)
    }

    fn concat(&mut self) -> Option<String> {
        let mut v = self.term()?;
        loop {
            self.ws();
            // '.' is concatenation unless it begins a number
            if self.i < self.s.len() && self.s[self.i] == b'.' && !self.s.get(self.i + 1).is_some_and(|c| c.is_ascii_digit()) {
                self.i += 1;
                v.push_str(&self.term()?);
            } else {
                return Some(v);
            }
        }
    }

    fn term(&mut self) -> Option<String> {
        self.ws();
        if self.eat("var(") {
            let start = self.i;
            while self.i < self.s.len() && self.s[self.i] != b')' {
                self.i += 1;
            }
            let name = std::str::from_utf8(&self.s[start..self.i]).ok()?.trim().to_string();
            self.i += 1;
            return (self.var)(&name);
        }
        if self.i < self.s.len() && (self.s[self.i] == b'"' || self.s[self.i] == b'\'') {
            let q = self.s[self.i];
            self.i += 1;
            let start = self.i;
            while self.i < self.s.len() && self.s[self.i] != q {
                self.i += 1;
            }
            let t = std::str::from_utf8(&self.s[start..self.i]).ok()?.to_string();
            self.i += 1;
            return Some(t);
        }
        if self.eat("(") {
            let v = self.ternary()?;
            return self.eat(")").then_some(v);
        }
        let start = self.i;
        while self.i < self.s.len() && (self.s[self.i].is_ascii_alphanumeric() || self.s[self.i] == b'.' || self.s[self.i] == b'_' || self.s[self.i] == b'-') {
            self.i += 1;
        }
        (self.i > start).then(|| String::from_utf8_lossy(&self.s[start..self.i]).to_string())
    }
}

/// A small backtracking regular expression: literals, `.`, character classes `[...]`, groups with
/// alternation `( | )`, the quantifiers `? * +`, and the anchors `^ $`.
#[derive(Debug, Clone)]
enum Node {
    Char(char),
    Any,
    Class(Vec<(char, char)>, bool),
    Group(Vec<Vec<Node>>),
    Start,
    End,
    Repeat(Box<Node>, usize, usize),
}

pub struct Regex {
    alts: Vec<Vec<Node>>,
}

impl Regex {
    pub fn parse(re: &str) -> Option<Regex> {
        let chars: Vec<char> = re.chars().collect();
        let mut i = 0;
        let alts = parse_alts(&chars, &mut i)?;
        (i == chars.len()).then_some(Regex { alts })
    }

    pub fn is_match(&self, s: &str) -> bool {
        let chars: Vec<char> = s.chars().collect();
        (0..=chars.len()).any(|start| match_alts(&self.alts, &chars, start, &mut |_| true))
    }
}

fn parse_alts(c: &[char], i: &mut usize) -> Option<Vec<Vec<Node>>> {
    let mut alts = vec![Vec::new()];
    while *i < c.len() {
        match c[*i] {
            ')' => break,
            '|' => {
                *i += 1;
                alts.push(Vec::new());
                continue;
            }
            _ => {}
        }
        let node = match c[*i] {
            '(' => {
                *i += 1;
                let g = parse_alts(c, i)?;
                if *i >= c.len() || c[*i] != ')' {
                    return None;
                }
                *i += 1;
                Node::Group(g)
            }
            '[' => {
                *i += 1;
                let neg = *i < c.len() && c[*i] == '^';
                if neg {
                    *i += 1;
                }
                let mut ranges = Vec::new();
                while *i < c.len() && c[*i] != ']' {
                    let a = c[*i];
                    if *i + 2 < c.len() && c[*i + 1] == '-' && c[*i + 2] != ']' {
                        ranges.push((a, c[*i + 2]));
                        *i += 3;
                    } else {
                        ranges.push((a, a));
                        *i += 1;
                    }
                }
                *i += 1;
                Node::Class(ranges, neg)
            }
            '.' => {
                *i += 1;
                Node::Any
            }
            '^' => {
                *i += 1;
                Node::Start
            }
            '$' => {
                *i += 1;
                Node::End
            }
            '\\' => {
                *i += 2;
                Node::Char(*c.get(*i - 1)?)
            }
            ch => {
                *i += 1;
                Node::Char(ch)
            }
        };
        let node = if *i < c.len() && matches!(c[*i], '?' | '*' | '+') {
            let (lo, hi) = match c[*i] {
                '?' => (0, 1),
                '*' => (0, usize::MAX),
                _ => (1, usize::MAX),
            };
            *i += 1;
            Node::Repeat(Box::new(node), lo, hi)
        } else {
            node
        };
        alts.last_mut().unwrap().push(node);
    }
    Some(alts)
}

fn match_alts(alts: &[Vec<Node>], s: &[char], pos: usize, k: &mut dyn FnMut(usize) -> bool) -> bool {
    alts.iter().any(|seq| match_seq(seq, s, pos, k))
}

fn match_seq(seq: &[Node], s: &[char], pos: usize, k: &mut dyn FnMut(usize) -> bool) -> bool {
    let Some((first, rest)) = seq.split_first() else { return k(pos) };
    match_node(first, s, pos, &mut |p| match_seq(rest, s, p, k))
}

fn match_node(n: &Node, s: &[char], pos: usize, k: &mut dyn FnMut(usize) -> bool) -> bool {
    match n {
        Node::Char(c) => pos < s.len() && s[pos] == *c && k(pos + 1),
        Node::Any => pos < s.len() && k(pos + 1),
        Node::Class(ranges, neg) => pos < s.len() && (ranges.iter().any(|(a, b)| *a <= s[pos] && s[pos] <= *b) != *neg) && k(pos + 1),
        Node::Start => pos == 0 && k(pos),
        Node::End => pos == s.len() && k(pos),
        Node::Group(alts) => match_alts(alts, s, pos, k),
        Node::Repeat(inner, lo, hi) => repeat(inner, *lo, *hi, s, pos, 0, k),
    }
}

fn repeat(n: &Node, lo: usize, hi: usize, s: &[char], pos: usize, count: usize, k: &mut dyn FnMut(usize) -> bool) -> bool {
    // greedy: one more first, then stop
    if count < hi && match_node(n, s, pos, &mut |p| p != pos && repeat(n, lo, hi, s, p, count + 1, k)) {
        return true;
    }
    count >= lo && k(pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linx_ideogram_label_format() {
        let f = r#"eval( var(chr) =~ /hs([6]|[9])$/ ? var(label) . ".": var(label) )"#;
        let v = |chr: &'static str, label: &'static str| move |n: &str| match n {
            "chr" => Some(chr.to_string()),
            "label" => Some(label.to_string()),
            _ => None,
        };
        assert_eq!(eval(f, &v("hs9", "9")).as_deref(), Some("9."));
        assert_eq!(eval(f, &v("hs6", "6")).as_deref(), Some("6."));
        assert_eq!(eval(f, &v("hs19", "19")).as_deref(), Some("19"));
        assert_eq!(eval(f, &v("hs16", "16")).as_deref(), Some("16"));
        assert_eq!(eval(f, &v("hs1", "1")).as_deref(), Some("1"));
        assert_eq!(eval(f, &v("hsX", "X")).as_deref(), Some("X"));
    }

    #[test]
    fn regex_basics() {
        assert!(Regex::parse("^a+b?c$").unwrap().is_match("aaac"));
        assert!(!Regex::parse("^a+b?c$").unwrap().is_match("bc"));
        assert!(Regex::parse("x[0-9]*y").unwrap().is_match("zzx123y"));
        assert!(Regex::parse("(ab|cd)e").unwrap().is_match("cde"));
    }
}
