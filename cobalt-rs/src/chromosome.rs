//! HumanChromosome (1-22, X, Y) and data kept per chromosome, and the order in which Java's
//! ArrayListMultimap<HumanChromosome, ..> iterates its keys.

/// HumanChromosome: index 0-21 for 1-22, 22 for X, 23 for Y.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Chr(u8);

const NAMES: [&str; 24] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y",
];

impl Chr {
    pub const X: Chr = Chr(22);
    pub const Y: Chr = Chr(23);

    #[cfg_attr(not(test), allow(dead_code))]
    pub fn new(index: usize) -> Chr {
        assert!(index < 24);
        Chr(index as u8)
    }
    pub fn all() -> impl Iterator<Item = Chr> {
        (0..24u8).map(Chr)
    }
    pub fn index(self) -> usize {
        self.0 as usize
    }
    /// HumanChromosome.toString: "1" .. "22", "X", "Y"
    pub fn short_name(self) -> &'static str {
        NAMES[self.index()]
    }
    pub fn is_autosome(self) -> bool {
        self.0 < 22
    }
    /// RefGenomeVersion.versionedChromosome
    pub fn versioned(self, v38: bool) -> String {
        if v38 { format!("chr{}", self.short_name()) } else { self.short_name().to_string() }
    }
    /// HumanChromosome.fromString: an optional "chr" prefix (any case), then the name exactly
    pub fn from_name(name: &str) -> Option<Chr> {
        let bare = if name.len() >= 3 && name[..3].eq_ignore_ascii_case("chr") { &name[3..] } else { name };
        NAMES.iter().position(|n| *n == bare).map(|i| Chr(i as u8))
    }
    /// HumanChromosome.contains, then fromString: the human chromosomes of a BAM header or the GC profile
    pub fn from_contig(contig: &str) -> Option<Chr> {
        // RefGenomeFunctions.stripChrPrefix removes "chr" only
        let bare = contig.strip_prefix("chr").unwrap_or(contig);
        let human = if !bare.is_empty() && bare.bytes().all(|b| b.is_ascii_digit()) {
            bare.parse::<i64>().is_ok_and(|n| (1..=22).contains(&n))
        } else {
            bare == "X" || bare == "Y"
        };
        if human { Chr::from_name(contig) } else { None }
    }
}

/// One value per human chromosome.
#[derive(Clone, Debug)]
pub struct PerChr<T>(pub [T; 24]);

impl<T: Default> Default for PerChr<T> {
    fn default() -> Self {
        PerChr(std::array::from_fn(|_| T::default()))
    }
}

impl<T> PerChr<T> {
    pub fn get(&self, c: Chr) -> &T {
        &self.0[c.index()]
    }
    pub fn get_mut(&mut self, c: Chr) -> &mut T {
        &mut self.0[c.index()]
    }
}

impl hmf_common::java_id_hash::IdKey for Chr {
    fn id_addr(&self) -> usize {
        self.index()
    }
}

/// The key order of an ArrayListMultimap<HumanChromosome, ..> (Guava: a HashMap of expected size 12) filled
/// in `insertion` order, on a JVM where every identity hash is the same (`-XX:hashCode=2`; enums hash by
/// identity). COBALT copies its read depths from map to map, so `copies` is the number of maps in the chain.
pub fn java_multimap_order(insertion: &[Chr], copies: usize) -> Vec<Chr> {
    let mut order = insertion.to_vec();
    for _ in 0..copies {
        // Maps.capacity(12) = 17
        let mut map: hmf_common::java_id_hash::IdHashMap<Chr, ()> = hmf_common::java_id_hash::IdHashMap::with_capacity(17);
        for c in &order {
            map.insert(*c, ());
        }
        order = map.keys().copied().collect();
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chrs(s: &str) -> Vec<Chr> {
        s.split(',').map(|c| Chr::from_name(c).unwrap()).collect()
    }

    #[test]
    fn names() {
        assert_eq!(Chr::from_name("chr1"), Some(Chr::new(0)));
        assert_eq!(Chr::from_name("CHRX"), Some(Chr::X));
        assert_eq!(Chr::from_name("Y"), Some(Chr::Y));
        assert_eq!(Chr::from_name("chrM"), None);
        assert_eq!(Chr::from_contig("chr22"), Some(Chr::new(21)));
        assert_eq!(Chr::from_contig("chr1_KI270706v1_random"), None);
        assert_eq!(Chr::from_contig("HLA-A*01:01:01:01"), None);
        assert_eq!(Chr::new(22).versioned(true), "chrX");
        assert_eq!(Chr::new(9).versioned(false), "10");
    }

    /// The orders printed by harness/OrderRef.java under -XX:hashCode=2
    #[test]
    fn java_key_order_with_constant_identity_hashes() {
        let all = chrs("1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y");
        assert_eq!(java_multimap_order(&all, 1), chrs("8,4,1,2,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"));
        assert_eq!(java_multimap_order(&all, 2), chrs("7,2,8,4,1,3,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"));
        assert_eq!(java_multimap_order(&all, 3), chrs("6,4,7,2,8,1,3,5,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"));
        let some = chrs("1,2,3,4,5,6,7,8,9,10,21,22");
        assert_eq!(java_multimap_order(&some, 3), chrs("2,3,4,1,5,6,7,8,9,10,21,22"));
        assert_eq!(java_multimap_order(&chrs("21,22"), 3), chrs("21,22"));
    }
}
