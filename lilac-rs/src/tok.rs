//! Sequence tokens: the Java code keeps each amino acid or nucleotide locus as a String, almost always of one
//! character (an insert makes it longer, a missing locus is ""). A token holds a one-byte string as the byte
//! itself and any other string as an index into a table shared by the whole run.

use std::collections::HashMap;
use std::sync::{OnceLock, RwLock};

pub type Tok = u16;

/// ""
pub const EMPTY: Tok = 0;
/// "." (a deleted base or amino acid)
pub const DEL: Tok = b'.' as Tok;
/// "*" (a wildcard locus of a reference sequence)
pub const WILD: Tok = b'*' as Tok;

struct Table {
    strings: Vec<String>,
    index: HashMap<String, Tok>,
}

fn table() -> &'static RwLock<Table> {
    static T: OnceLock<RwLock<Table>> = OnceLock::new();
    T.get_or_init(|| RwLock::new(Table { strings: Vec::new(), index: HashMap::new() }))
}

/// The token of a string
pub fn intern(s: &str) -> Tok {
    let b = s.as_bytes();
    if b.is_empty() {
        return EMPTY;
    }
    if b.len() == 1 && b[0] != 0 {
        return b[0] as Tok;
    }
    if let Some(t) = table().read().unwrap().index.get(s) {
        return *t;
    }
    let mut t = table().write().unwrap();
    if let Some(v) = t.index.get(s) {
        return *v;
    }
    let tok = 256 + t.strings.len();
    assert!(tok <= Tok::MAX as usize, "too many distinct multi-character sequences");
    t.strings.push(s.to_string());
    t.index.insert(s.to_string(), tok as Tok);
    tok as Tok
}

/// The token of one byte
#[inline]
pub fn of_byte(b: u8) -> Tok {
    b as Tok
}

#[inline]
pub fn is_single(t: Tok) -> bool {
    t > 0 && t < 256
}

/// Appends the string of a token
pub fn push(t: Tok, out: &mut String) {
    if t == EMPTY {
        return;
    }
    if t < 256 {
        out.push(t as u8 as char);
    } else {
        out.push_str(&table().read().unwrap().strings[t as usize - 256]);
    }
}

/// Appends the bytes of a token
pub fn push_bytes(t: Tok, out: &mut Vec<u8>) {
    if t == EMPTY {
        return;
    }
    if t < 256 {
        out.push(t as u8);
    } else {
        out.extend_from_slice(table().read().unwrap().strings[t as usize - 256].as_bytes());
    }
}

pub fn to_string(t: Tok) -> String {
    let mut s = String::new();
    push(t, &mut s);
    s
}

/// String.length of the token's string
pub fn len(t: Tok) -> usize {
    if t == EMPTY {
        0
    } else if t < 256 {
        1
    } else {
        table().read().unwrap().strings[t as usize - 256].len()
    }
}

/// String.hashCode of the token's string
pub fn java_hash(t: Tok) -> i32 {
    if t < 256 {
        t as i32
    } else {
        hmf_common::java_hash::java_string_hash(&table().read().unwrap().strings[t as usize - 256])
    }
}

/// The concatenation of tokens
pub fn join(toks: impl IntoIterator<Item = Tok>) -> String {
    let mut s = String::new();
    for t in toks {
        push(t, &mut s);
    }
    s
}

/// Token of `a + b + c` (string concatenation)
pub fn concat3(a: Tok, b: Tok, c: Tok) -> String {
    let mut s = String::new();
    push(a, &mut s);
    push(b, &mut s);
    push(c, &mut s);
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens() {
        assert_eq!(intern(""), EMPTY);
        assert_eq!(intern("A"), b'A' as Tok);
        let t = intern("ABC");
        assert!(t >= 256);
        assert_eq!(intern("ABC"), t);
        assert_eq!(to_string(t), "ABC");
        assert_eq!(len(t), 3);
        assert_eq!(java_hash(t), hmf_common::java_hash::java_string_hash("ABC"));
        assert_eq!(java_hash(intern("G")), 71);
        assert_eq!(join([intern("A"), EMPTY, t]), "AABC");
    }
}
