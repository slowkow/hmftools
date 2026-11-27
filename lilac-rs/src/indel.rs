//! read.Indel

#[derive(Clone, Debug)]
pub struct Indel {
    pub contig: String,
    pub position: i32,
    pub read_index: i32,
    pub ref_: String,
    pub alt: String,
    pub is_insert: bool,
    pub length: i32,
    /// only set for SBX and Ultima reads
    pub ignore: bool,
}

impl Indel {
    pub fn new(contig: &str, position: i32, read_index: i32, ref_: &str, alt: &str, ignore: bool) -> Indel {
        let is_insert = ref_.len() < alt.len();
        Indel {
            contig: contig.to_string(),
            position,
            read_index,
            ref_: ref_.to_string(),
            alt: alt.to_string(),
            is_insert,
            length: alt.len() as i32 - ref_.len() as i32,
            ignore,
        }
    }

    pub fn of(contig: &str, position: i32, ref_: &str, alt: &str) -> Indel {
        Indel::new(contig, position, -1, ref_, alt, false)
    }

    /// Indel.fromString: chromosome,position,ref,alt
    pub fn from_csv(line: &str) -> Option<Indel> {
        let items = hmf_common::java_text::split(line, ',');
        if items.len() < 4 {
            return None;
        }
        Some(Indel::of(items[0], hmf_common::java_text::parse_int(items[1])?, items[2], items[3]))
    }

    pub fn key(&self) -> String {
        format!("{}:{} {}>{}", self.contig, self.position, self.ref_, self.alt)
    }

    /// Indel.match
    pub fn matches(&self, other: &Indel) -> bool {
        if self == other {
            return true;
        }
        if self.ignore != other.ignore || self.contig != other.contig || self.position != other.position {
            return false;
        }
        if self.ref_.len() != other.ref_.len() || self.alt.len() != other.alt.len() {
            return false;
        }
        if self.alt.len() > self.ref_.len() {
            return self.alt[1..] == other.alt[1..];
        }
        true
    }
}

impl std::fmt::Display for Indel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{} {}>{}", self.contig, self.position, self.ref_, self.alt)
    }
}

/// equals: the text form and the ignore flag
impl PartialEq for Indel {
    fn eq(&self, other: &Self) -> bool {
        self.contig == other.contig
            && self.position == other.position
            && self.ref_ == other.ref_
            && self.alt == other.alt
            && self.ignore == other.ignore
    }
}
impl Eq for Indel {}

impl std::hash::Hash for Indel {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.contig.hash(state);
        self.position.hash(state);
        self.ref_.hash(state);
        self.alt.hash(state);
        self.ignore.hash(state);
    }
}

/// The indel panel of normals of the jar (resources/pon/indels_v37.csv, indels_v38.csv)
pub fn pon_indels(v38: bool) -> Vec<Indel> {
    let text = if v38 { include_str!("../resources/indels_v38.csv") } else { include_str!("../resources/indels_v37.csv") };
    hmf_common::java_text::lines_of(text).iter().filter_map(|l| Indel::from_csv(l)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indels() {
        let a = Indel::of("chr6", 100, "A", "ACG");
        assert!(a.is_insert);
        assert_eq!(a.length, 2);
        let b = Indel::new("chr6", 100, 12, "A", "ACG", false);
        assert_eq!(a, b);
        assert!(a.matches(&Indel::of("chr6", 100, "T", "TCG")));
        assert!(!a.matches(&Indel::of("chr6", 100, "T", "TCA")));
        assert_eq!(pon_indels(true).len(), 9);
    }
}
