//! hmf-common Orientation.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    Reverse,
}

pub use Orientation::{Forward as FORWARD, Reverse as REVERSE};

impl Orientation {
    #[inline]
    pub fn is_forward(self) -> bool {
        self == Orientation::Forward
    }
    #[inline]
    pub fn is_reverse(self) -> bool {
        self == Orientation::Reverse
    }
    #[inline]
    pub fn as_byte(self) -> i8 {
        if self.is_forward() { 1 } else { -1 }
    }
    pub fn from_byte(b: i8) -> Orientation {
        match b {
            1 => Orientation::Forward,
            -1 => Orientation::Reverse,
            _ => panic!("Invalid orientation: {b}"),
        }
    }
    pub fn from_byte_str(s: &str) -> Orientation {
        Orientation::from_byte(s.parse::<i8>().unwrap_or_else(|_| panic!("Invalid orientation: {s}")))
    }
    pub fn from_char(c: char) -> Orientation {
        match c {
            '+' => Orientation::Forward,
            '-' => Orientation::Reverse,
            _ => panic!("Invalid orientation: {c}"),
        }
    }
    pub fn as_char(self) -> char {
        if self.is_forward() { '+' } else { '-' }
    }
    #[inline]
    pub fn opposite(self) -> Orientation {
        if self.is_forward() { Orientation::Reverse } else { Orientation::Forward }
    }
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.as_byte())
    }
}
