//! JSON text as Gson 2.13 writes it for ORANGE (OrangeJson: `setPrettyPrinting()`, `serializeNulls()`,
//! HTML-safe escaping, lenient numbers):
//! - two-space indent, `"name": value`, a newline before every member and element, `[]` and `{}` for
//!   empty arrays and objects, no newline at the end;
//! - strings escaped as JsonWriter's HTML_SAFE_REPLACEMENT_CHARS: `"` `\` and control characters, `<`
//!   `>` `&` `=` `'` as `<` ..., U+2028 and U+2029;
//! - doubles as `Double.toString` (NaN and Infinity written as such), integers in decimal.

use hmf_common::fmt::java_double_to_string;

#[derive(Clone, Copy, PartialEq)]
enum Scope {
    EmptyArray,
    NonEmptyArray,
    EmptyObject,
    NonEmptyObject,
    DanglingName,
    EmptyDocument,
    NonEmptyDocument,
}

pub struct JsonWriter {
    out: String,
    stack: Vec<Scope>,
    deferred_name: Option<String>,
}

impl Default for JsonWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl JsonWriter {
    pub fn new() -> JsonWriter {
        JsonWriter { out: String::with_capacity(1 << 16), stack: vec![Scope::EmptyDocument], deferred_name: None }
    }

    pub fn finish(self) -> String {
        self.out
    }

    fn top(&self) -> Scope {
        *self.stack.last().expect("JSON scope")
    }

    fn replace_top(&mut self, s: Scope) {
        *self.stack.last_mut().unwrap() = s;
    }

    fn newline(&mut self) {
        self.out.push('\n');
        for _ in 1..self.stack.len() {
            self.out.push_str("  ");
        }
    }

    fn before_name(&mut self) {
        match self.top() {
            Scope::NonEmptyObject => self.out.push(','),
            Scope::EmptyObject => {}
            _ => panic!("JsonWriter: nesting problem"),
        }
        self.newline();
        self.replace_top(Scope::DanglingName);
    }

    fn write_deferred_name(&mut self) {
        if let Some(name) = self.deferred_name.take() {
            self.before_name();
            self.string(&name);
        }
    }

    fn before_value(&mut self) {
        match self.top() {
            Scope::NonEmptyDocument => panic!("JSON must have only one top-level value"),
            Scope::EmptyDocument => self.replace_top(Scope::NonEmptyDocument),
            Scope::EmptyArray => {
                self.replace_top(Scope::NonEmptyArray);
                self.newline();
            }
            Scope::NonEmptyArray => {
                self.out.push(',');
                self.newline();
            }
            Scope::DanglingName => {
                self.out.push_str(": ");
                self.replace_top(Scope::NonEmptyObject);
            }
            _ => panic!("JsonWriter: nesting problem"),
        }
    }

    fn open(&mut self, empty: Scope, bracket: char) {
        self.write_deferred_name();
        self.before_value();
        self.stack.push(empty);
        self.out.push(bracket);
    }

    fn close(&mut self, empty: Scope, nonempty: Scope, bracket: char) {
        let ctx = self.top();
        assert!(ctx == empty || ctx == nonempty, "JsonWriter: nesting problem");
        assert!(self.deferred_name.is_none(), "dangling name");
        self.stack.pop();
        if ctx == nonempty {
            self.newline();
        }
        self.out.push(bracket);
    }

    pub fn begin_object(&mut self) -> &mut Self {
        self.open(Scope::EmptyObject, '{');
        self
    }
    pub fn end_object(&mut self) -> &mut Self {
        self.close(Scope::EmptyObject, Scope::NonEmptyObject, '}');
        self
    }
    pub fn begin_array(&mut self) -> &mut Self {
        self.open(Scope::EmptyArray, '[');
        self
    }
    pub fn end_array(&mut self) -> &mut Self {
        self.close(Scope::EmptyArray, Scope::NonEmptyArray, ']');
        self
    }

    pub fn name(&mut self, name: &str) -> &mut Self {
        assert!(self.deferred_name.is_none(), "already a name");
        self.deferred_name = Some(name.to_string());
        self
    }

    fn string(&mut self, s: &str) {
        self.out.push('"');
        for c in s.chars() {
            match c {
                '"' => self.out.push_str("\\\""),
                '\\' => self.out.push_str("\\\\"),
                '\t' => self.out.push_str("\\t"),
                '\u{8}' => self.out.push_str("\\b"),
                '\n' => self.out.push_str("\\n"),
                '\r' => self.out.push_str("\\r"),
                '\u{c}' => self.out.push_str("\\f"),
                '<' => self.out.push_str("\\u003c"),
                '>' => self.out.push_str("\\u003e"),
                '&' => self.out.push_str("\\u0026"),
                '=' => self.out.push_str("\\u003d"),
                '\'' => self.out.push_str("\\u0027"),
                '\u{2028}' => self.out.push_str("\\u2028"),
                '\u{2029}' => self.out.push_str("\\u2029"),
                c if (c as u32) < 0x20 => self.out.push_str(&format!("\\u{:04x}", c as u32)),
                c => self.out.push(c),
            }
        }
        self.out.push('"');
    }

    pub fn null(&mut self) -> &mut Self {
        self.write_deferred_name();
        self.before_value();
        self.out.push_str("null");
        self
    }
    pub fn str(&mut self, v: &str) -> &mut Self {
        self.write_deferred_name();
        self.before_value();
        self.string(v);
        self
    }
    pub fn opt_str(&mut self, v: Option<&str>) -> &mut Self {
        match v {
            Some(s) => self.str(s),
            None => self.null(),
        }
    }
    pub fn bool(&mut self, v: bool) -> &mut Self {
        self.write_deferred_name();
        self.before_value();
        self.out.push_str(if v { "true" } else { "false" });
        self
    }
    pub fn int(&mut self, v: i64) -> &mut Self {
        self.write_deferred_name();
        self.before_value();
        self.out.push_str(&v.to_string());
        self
    }
    pub fn opt_int(&mut self, v: Option<i64>) -> &mut Self {
        match v {
            Some(i) => self.int(i),
            None => self.null(),
        }
    }
    pub fn double(&mut self, v: f64) -> &mut Self {
        self.write_deferred_name();
        self.before_value();
        self.out.push_str(&java_double_to_string(v));
        self
    }
    pub fn opt_double(&mut self, v: Option<f64>) -> &mut Self {
        match v {
            Some(d) => self.double(d),
            None => self.null(),
        }
    }

    // name + value shorthands
    pub fn f_str(&mut self, n: &str, v: &str) -> &mut Self {
        self.name(n).str(v)
    }
    pub fn f_opt_str(&mut self, n: &str, v: Option<&str>) -> &mut Self {
        self.name(n).opt_str(v)
    }
    pub fn f_bool(&mut self, n: &str, v: bool) -> &mut Self {
        self.name(n).bool(v)
    }
    pub fn f_int(&mut self, n: &str, v: i64) -> &mut Self {
        self.name(n).int(v)
    }
    pub fn f_opt_int(&mut self, n: &str, v: Option<i64>) -> &mut Self {
        self.name(n).opt_int(v)
    }
    pub fn f_double(&mut self, n: &str, v: f64) -> &mut Self {
        self.name(n).double(v)
    }
    pub fn f_opt_double(&mut self, n: &str, v: Option<f64>) -> &mut Self {
        self.name(n).opt_double(v)
    }
}

/// A value Gson writes through a type adapter
pub trait ToJson {
    fn write_json(&self, w: &mut JsonWriter);
}

impl JsonWriter {
    pub fn f_value<T: ToJson>(&mut self, n: &str, v: &T) -> &mut Self {
        self.name(n);
        v.write_json(self);
        self
    }
    /// a @Nullable object attribute: null when absent (serializeNulls)
    pub fn f_opt<T: ToJson>(&mut self, n: &str, v: Option<&T>) -> &mut Self {
        self.name(n);
        match v {
            Some(v) => v.write_json(self),
            None => {
                self.null();
            }
        }
        self
    }
    /// a list attribute
    pub fn f_list<T: ToJson>(&mut self, n: &str, v: &[T]) -> &mut Self {
        self.name(n).begin_array();
        for e in v {
            e.write_json(self);
        }
        self.end_array()
    }
    /// a @Nullable list attribute
    pub fn f_opt_list<T: ToJson>(&mut self, n: &str, v: Option<&[T]>) -> &mut Self {
        match v {
            Some(v) => self.f_list(n, v),
            None => self.name(n).null(),
        }
    }
    /// a list of enum names or strings
    pub fn f_str_list<S: AsRef<str>>(&mut self, n: &str, v: &[S]) -> &mut Self {
        self.name(n).begin_array();
        for e in v {
            self.str(e.as_ref());
        }
        self.end_array()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gson_pretty_layout() {
        let mut w = JsonWriter::new();
        w.begin_object();
        w.f_str("a", "x<y>&='\"\\\n\u{1}");
        w.name("b").begin_array().end_array();
        w.name("c").begin_object().end_object();
        w.name("d").begin_array().int(1).double(0.5).null().end_array();
        w.f_double("e", 1.0).f_double("f", f64::NAN).f_bool("g", false);
        w.name("h").begin_object().f_int("i", -3).end_object();
        w.end_object();
        assert_eq!(
            w.finish(),
            "{\n  \"a\": \"x\\u003cy\\u003e\\u0026\\u003d\\u0027\\\"\\\\\\n\\u0001\",\n  \"b\": [],\n  \"c\": {},\n  \"d\": [\n    1,\n    0.5,\n    null\n  ],\n  \"e\": 1.0,\n  \"f\": NaN,\n  \"g\": false,\n  \"h\": {\n    \"i\": -3\n  }\n}"
        );
    }
}
