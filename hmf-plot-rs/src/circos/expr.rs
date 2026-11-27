//! The arithmetic of Circos's expressions, enough for rule conditions (`var(value) < 4`), rule values
//! (`eval(int(remap_int(var(value),4,12,15,30)))`) and radius expressions
//! (`(dims(ideogram,radius_outer)+dims(ideogram,radius_inner))/2`).
//!
//! Numbers only: comparisons give 1 or 0. Names are resolved by the caller: `var(x)`, `dims(a,b)`,
//! `conf(a,b)`.

use crate::fmath::FMath;

pub trait Resolver {
    fn var(&self, name: &str) -> Option<f64>;
    fn dims(&self, args: &[&str]) -> Option<f64>;
    fn conf(&self, args: &[&str]) -> Option<f64>;
}

pub fn eval(expr: &str, r: &dyn Resolver) -> Result<f64, String> {
    let s = expr.trim();
    let s = s.strip_prefix("eval(").and_then(|x| x.strip_suffix(')')).unwrap_or(s);
    let mut p = Parser { s: s.as_bytes(), i: 0, r };
    let v = p.or()?;
    p.ws();
    if p.i != p.s.len() {
        return Err(format!("cannot parse '{expr}' at {}", p.i));
    }
    Ok(v)
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
    r: &'a dyn Resolver,
}

impl Parser<'_> {
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

    fn or(&mut self) -> Result<f64, String> {
        let mut v = self.and()?;
        while self.eat("||") {
            let w = self.and()?;
            v = f64::from(v != 0.0 || w != 0.0);
        }
        Ok(v)
    }

    fn and(&mut self) -> Result<f64, String> {
        let mut v = self.cmp()?;
        while self.eat("&&") {
            let w = self.cmp()?;
            v = f64::from(v != 0.0 && w != 0.0);
        }
        Ok(v)
    }

    fn cmp(&mut self) -> Result<f64, String> {
        let v = self.sum()?;
        for (op, f) in [
            ("<=", (|a: f64, b: f64| a <= b) as fn(f64, f64) -> bool),
            (">=", |a, b| a >= b),
            ("==", |a, b| a == b),
            ("!=", |a, b| a != b),
            ("<", |a, b| a < b),
            (">", |a, b| a > b),
        ] {
            if self.eat(op) {
                let w = self.sum()?;
                return Ok(f64::from(f(v, w)));
            }
        }
        Ok(v)
    }

    fn sum(&mut self) -> Result<f64, String> {
        let mut v = self.prod()?;
        loop {
            if self.eat("+") {
                v += self.prod()?;
            } else if self.eat("-") {
                v -= self.prod()?;
            } else {
                return Ok(v);
            }
        }
    }

    fn prod(&mut self) -> Result<f64, String> {
        let mut v = self.unary()?;
        loop {
            if self.eat("*") {
                v *= self.unary()?;
            } else if self.eat("/") {
                v /= self.unary()?;
            } else {
                return Ok(v);
            }
        }
    }

    fn unary(&mut self) -> Result<f64, String> {
        if self.eat("-") {
            return Ok(-self.unary()?);
        }
        if self.eat("!") {
            return Ok(f64::from(self.unary()? == 0.0));
        }
        self.atom()
    }

    fn atom(&mut self) -> Result<f64, String> {
        self.ws();
        if self.eat("(") {
            let v = self.or()?;
            if !self.eat(")") {
                return Err("missing )".into());
            }
            return Ok(v);
        }
        let start = self.i;
        if self.i < self.s.len() && (self.s[self.i].is_ascii_digit() || self.s[self.i] == b'.') {
            while self.i < self.s.len() && (self.s[self.i].is_ascii_digit() || matches!(self.s[self.i], b'.' | b'e' | b'E')) {
                self.i += 1;
            }
            let t = std::str::from_utf8(&self.s[start..self.i]).unwrap();
            return t.parse().map_err(|_| format!("bad number '{t}'"));
        }
        while self.i < self.s.len() && (self.s[self.i].is_ascii_alphanumeric() || self.s[self.i] == b'_') {
            self.i += 1;
        }
        let name = std::str::from_utf8(&self.s[start..self.i]).unwrap().to_string();
        if name.is_empty() {
            return Err(format!("unexpected character at {}", self.i));
        }
        if !self.eat("(") {
            return Err(format!("unknown name '{name}'"));
        }
        // raw-text arguments for var/dims/conf, numbers for the rest
        if matches!(name.as_str(), "var" | "dims" | "conf") {
            let close = self.s[self.i..].iter().position(|&c| c == b')').ok_or("missing )")?;
            let inner = std::str::from_utf8(&self.s[self.i..self.i + close]).unwrap();
            self.i += close + 1;
            let args: Vec<&str> = inner.split(',').map(|a| a.trim()).collect();
            let v = match name.as_str() {
                "var" => self.r.var(args[0]),
                "dims" => self.r.dims(&args),
                _ => self.r.conf(&args),
            };
            return v.ok_or_else(|| format!("{name}({inner}) is not defined"));
        }
        let mut args = Vec::new();
        if !self.eat(")") {
            loop {
                args.push(self.or()?);
                if self.eat(")") {
                    break;
                }
                if !self.eat(",") {
                    return Err("expected , or )".into());
                }
            }
        }
        let a = |i: usize| args.get(i).copied().ok_or_else(|| format!("{name}: missing argument"));
        Ok(match name.as_str() {
            "int" => a(0)?.trunc(),
            "abs" => a(0)?.abs(),
            "round" => perl_round(a(0)?),
            "log10" => a(0)?.flog10(),
            "min" => args.iter().copied().fold(f64::INFINITY, f64::min),
            "max" => args.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            "remap" => remap(a(0)?, a(1)?, a(2)?, a(3)?, a(4)?),
            "remap_int" => remap(a(0)?, a(1)?, a(2)?, a(3)?, a(4)?).trunc(),
            "remap_round" => perl_round(remap(a(0)?, a(1)?, a(2)?, a(3)?, a(4)?)),
            _ => return Err(format!("unknown function '{name}'")),
        })
    }
}

/// Circos::Utils::remap: value in [min, max] mapped linearly to [remap_min, remap_max], clamped at the ends.
pub fn remap(value: f64, min: f64, max: f64, rmin: f64, rmax: f64) -> f64 {
    if value <= min {
        return rmin;
    }
    if value >= max {
        return rmax;
    }
    if min == max {
        return rmin;
    }
    rmin + (value - min) / (max - min) * (rmax - rmin)
}

/// Math::Round::round: half away from zero.
pub fn perl_round(x: f64) -> f64 {
    if x >= 0.0 { (x + 0.5).floor() } else { (x - 0.5).ceil() }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct R;
    impl Resolver for R {
        fn var(&self, n: &str) -> Option<f64> {
            (n == "value").then_some(6.0)
        }
        fn dims(&self, a: &[&str]) -> Option<f64> {
            match a {
                ["ideogram", "radius_outer"] => Some(1350.0),
                ["ideogram", "radius_inner"] => Some(1300.0),
                _ => None,
            }
        }
        fn conf(&self, a: &[&str]) -> Option<f64> {
            (a == ["ideogram", "label_size"]).then_some(30.0)
        }
    }

    #[test]
    fn expressions() {
        assert_eq!(eval("var(value) < 4", &R).unwrap(), 0.0);
        assert_eq!(eval("var(value) > 4", &R).unwrap(), 1.0);
        // remap(6, 4, 12, 15, 30) = 18.75 -> 18
        assert_eq!(eval("eval(int(remap_int(var(value),4,12,15,30)))", &R).unwrap(), 18.0);
        assert_eq!(eval("(dims(ideogram,radius_outer)+dims(ideogram,radius_inner))/2", &R).unwrap(), 1325.0);
        assert_eq!(eval("(dims(ideogram,radius_inner)+dims(ideogram,radius_outer)-conf(ideogram,label_size))/2", &R).unwrap(), 1310.0);
        assert_eq!(remap(20.0, 4.0, 12.0, 15.0, 30.0), 30.0);
    }
}
