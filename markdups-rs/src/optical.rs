//! Picard's ReadNameParser and OpticalDuplicateFinder (picard.sam.util, picard.sam.markduplicates.util).

use anyhow::{Result, bail};
use regex::bytes::Regex;

use crate::compat;

/// Tile and cluster position of a read (PhysicalLocation): tile -1 means none.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Location {
    pub tile: i32,
    pub x: i32,
    pub y: i32,
}

pub const NO_LOCATION: Location = Location { tile: -1, x: -1, y: -1 };

pub enum NameParser {
    /// READ_NAME_REGEX null: no locations
    Off,
    /// the default: the last three ':'-separated fields of names with 5 or 7 fields
    Default,
    Regex(Regex),
}

/// ReadNameParser.rapidParseInt: an optional '-' and digits up to the first other character. Java's int
/// arithmetic wraps a number that does not fit (MD2); corrected, such a field gives no location.
fn rapid_parse_int(s: &[u8]) -> Option<i64> {
    let mut i = 0;
    let neg = s.first() == Some(&b'-');
    if neg {
        i = 1;
    }
    let mut val: i32 = 0;
    let mut wide: i64 = 0;
    let mut has = false;
    while i < s.len() && s[i].is_ascii_digit() {
        val = val.wrapping_mul(10).wrapping_add((s[i] - b'0') as i32);
        wide = wide.saturating_mul(10).saturating_add((s[i] - b'0') as i64);
        has = true;
        i += 1;
    }
    if !has {
        return None;
    }
    if compat::java("MD2") {
        Some(if neg { val.wrapping_neg() as i64 } else { val as i64 })
    } else {
        let w = if neg { -wide } else { wide };
        if w < i32::MIN as i64 || w > i32::MAX as i64 { None } else { Some(w) }
    }
}

/// ReadNameParser.getLastThreeFields: the number of fields and the last three parsed as numbers; None
/// when a field does not start with a number (NumberFormatException).
fn last_three_fields(name: &[u8]) -> Option<(i32, [i64; 3])> {
    let mut tokens = [-1i64; 3];
    let mut idx: i32 = 2;
    let mut num_fields = 0;
    let mut end = name.len();
    let mut i = name.len() as isize - 1;
    while i >= 0 && idx >= 0 {
        let iu = i as usize;
        if name[iu] == b':' || iu == 0 {
            num_fields += 1;
            let from = if iu == 0 { 0 } else { iu + 1 };
            tokens[idx as usize] = rapid_parse_int(&name[from..end])?;
            idx -= 1;
            end = iu;
        }
        i -= 1;
    }
    while i >= 0 {
        let iu = i as usize;
        if name[iu] == b':' || iu == 0 {
            num_fields += 1;
        }
        i -= 1;
    }
    if num_fields < 3 {
        return Some((-1, [-1, -1, -1]));
    }
    Some((num_fields, tokens))
}

impl NameParser {
    pub fn new(regex: Option<&str>) -> Result<NameParser> {
        Ok(match regex {
            None => NameParser::Off,
            Some(r) if r == crate::args::DEFAULT_READ_NAME_REGEX => NameParser::Default,
            Some(r) => match Regex::new(&format!("^(?:{r})$")) {
                Ok(re) => {
                    if re.captures_len() < 4 {
                        bail!("READ_NAME_REGEX must have three capture groups");
                    }
                    NameParser::Regex(re)
                }
                Err(e) => bail!("READ_NAME_REGEX '{r}' is not a regular expression markdups-rs can use: {e}"),
            },
        })
    }

    pub fn is_off(&self) -> bool {
        matches!(self, NameParser::Off)
    }

    /// ReadNameParser.addLocationInformation: the tile (as a short), x and y (as PhysicalLocationShort
    /// keeps them: cast to short, MD1 when corrected keeps the full value); None when the name gives none.
    pub fn location(&self, name: &[u8]) -> Option<Location> {
        let (tile, x, y) = match self {
            NameParser::Off => return None,
            NameParser::Default => {
                let (n, t) = last_three_fields(name)?;
                if !(n == 5 || n == 7) {
                    return None;
                }
                (t[0], t[1], t[2])
            }
            NameParser::Regex(re) => {
                let c = re.captures(name)?;
                let p = |i: usize| -> Option<i64> {
                    let s = std::str::from_utf8(c.get(i)?.as_bytes()).ok()?;
                    // Integer.parseInt: the whole group must be an int
                    s.parse::<i32>().ok().map(|v| v as i64)
                };
                (p(1)?, p(2)?, p(3)?)
            }
        };
        Some(make_location(tile, x, y))
    }
}

/// The location as Picard stores it: the tile cast to a short, and x and y cast to shorts by
/// PhysicalLocationShort, so that values above 32767 wrap (MD1); corrected, they are kept. The values
/// fit an int (MD2).
pub fn make_location(tile: i64, x: i64, y: i64) -> Location {
    if compat::java("MD1") {
        Location { tile: tile as i32 as i16 as i32, x: x as i32 as i16 as i32, y: y as i32 as i16 as i32 }
    } else {
        Location { tile: tile as i32, x: x as i32, y: y as i32 }
    }
}

/// What the optical duplicate finder needs of a read end.
#[derive(Clone, Copy, Debug)]
pub struct OpticalEnd {
    pub read_group: i16,
    pub loc: Location,
}

impl OpticalEnd {
    #[inline]
    fn has_location(&self) -> bool {
        self.loc.tile != -1
    }
}

pub struct OpticalFinder {
    pub distance: i32,
    pub max_set_size: i64,
}

impl OpticalFinder {
    #[inline]
    fn close_enough(&self, a: &OpticalEnd, b: &OpticalEnd) -> bool {
        a.has_location()
            && b.has_location()
            && a.read_group == b.read_group
            && a.loc.tile == b.loc.tile
            && (a.loc.x as i64 - b.loc.x as i64).abs() <= self.distance as i64
            && (a.loc.y as i64 - b.loc.y as i64).abs() <= self.distance as i64
    }

    /// OpticalDuplicateFinder.findOpticalDuplicates: flags for the reads of `list` that are optical
    /// duplicates; `keeper` is the index in `list` of the read kept, if it is in the list.
    pub fn find(&self, list: &[OpticalEnd], keeper: Option<usize>) -> Vec<bool> {
        let n = list.len();
        let mut flags = vec![false; n];
        // MAX_OPTICAL_DUPLICATE_SET_SIZE -1 is documented as "no limit"; Picard compares the set size with it
        // as given, so -1 turns optical duplicate detection off (MD3)
        let limit_applies = compat::java("MD3") || self.max_set_size >= 1;
        if n < 2 || (limit_applies && n as i64 > self.max_set_size) {
            return flags;
        }
        // keeperOrNull: the keeper only if it has a location
        let actual = keeper.filter(|&k| list[k].has_location());
        // `keeper == null ? 3 : 4` on the keeper passed in, which MarkDuplicates always has
        if n >= 4 { self.graph(list, actual, flags) } else { self.fast(list, actual, flags) }
    }

    fn fast(&self, list: &[OpticalEnd], keeper: Option<usize>, mut flags: Vec<bool>) -> Vec<bool> {
        let n = list.len();
        if let Some(k) = keeper {
            for i in 0..n {
                flags[i] = i != k && self.close_enough(&list[k], &list[i]);
            }
        }
        for i in 0..n {
            if Some(i) == keeper {
                continue;
            }
            for j in i + 1..n {
                if Some(j) == keeper {
                    continue;
                }
                if flags[i] && flags[j] {
                    continue;
                }
                if self.close_enough(&list[i], &list[j]) {
                    let index = if flags[j] { i } else { j };
                    flags[index] = true;
                }
            }
        }
        flags
    }

    /// getOpticalDuplicatesFlagWithGraph: reads within the distance of each other (same read group and
    /// tile) are joined into clusters; one read per cluster is not flagged: the keeper in its cluster,
    /// else the one with the lowest x, then y, first seen in the order Java iterates the cluster map.
    fn graph(&self, list: &[OpticalEnd], keeper: Option<usize>, mut flags: Vec<bool>) -> Vec<bool> {
        let n = list.len();
        let mut parent: Vec<u32> = (0..n as u32).collect();
        fn find(p: &mut [u32], mut i: u32) -> u32 {
            let mut r = i;
            while p[r as usize] != r {
                r = p[r as usize];
            }
            while p[i as usize] != r {
                let next = p[i as usize];
                p[i as usize] = r;
                i = next;
            }
            r
        }
        // group the reads with a location by (read group, tile), then join the pairs within the distance:
        // sorted by x, only reads within the distance in x need comparing (the same edges as Picard's
        // all-pairs loop, so the same clusters)
        let mut with_loc: Vec<u32> = (0..n as u32).filter(|&i| list[i as usize].has_location()).collect();
        // Picard's key ((readGroup << 16) + tile) is one per (read group, tile) for short tiles
        let key = |i: u32| {
            let e = &list[i as usize];
            (e.read_group, e.loc.tile)
        };
        with_loc.sort_unstable_by_key(|&i| (key(i), list[i as usize].loc.x, i));
        let d = self.distance as i64;
        let mut g = 0;
        while g < with_loc.len() {
            let mut h = g;
            while h < with_loc.len() && key(with_loc[h]) == key(with_loc[g]) {
                h += 1;
            }
            for a in g..h {
                let ea = list[with_loc[a] as usize];
                for &bi in &with_loc[a + 1..h] {
                    let eb = list[bi as usize];
                    if eb.loc.x as i64 - ea.loc.x as i64 > d {
                        break;
                    }
                    if (ea.loc.y as i64 - eb.loc.y as i64).abs() <= d {
                        let (ra, rb) = (find(&mut parent, with_loc[a]), find(&mut parent, bi));
                        if ra != rb {
                            parent[ra as usize] = rb;
                        }
                    }
                }
            }
            g = h;
        }
        let cluster: Vec<u32> = (0..n as u32).map(|i| find(&mut parent, i)).collect();
        let mut rep: hashbrown::HashMap<u32, usize> = hashbrown::HashMap::new();
        let keeper_cluster = keeper.map(|k| cluster[k]);
        if let Some(k) = keeper {
            rep.insert(cluster[k], k);
        }
        for i in java_integer_map_order(n) {
            let c = cluster[i];
            match rep.get(&c).copied() {
                Some(r) if Some(i) != keeper => {
                    let (cur, rl) = (list[i].loc, list[r].loc);
                    if keeper_cluster != Some(c) && (cur.x < rl.x || (cur.x == rl.x && cur.y < rl.y)) {
                        flags[r] = true;
                        rep.insert(c, i);
                    } else {
                        flags[i] = true;
                    }
                }
                _ => {
                    rep.insert(c, i);
                }
            }
        }
        flags
    }
}

/// The order in which java.util.HashMap (built by Collectors.toMap) iterates the Integer keys 0..n-1:
/// by bucket, (h ^ h >>> 16) & (capacity - 1), one key per bucket.
pub fn java_integer_map_order(n: usize) -> Vec<usize> {
    let mut cap = 16usize;
    while n as f64 > cap as f64 * 0.75 {
        cap *= 2;
    }
    let mut v: Vec<usize> = (0..n).collect();
    if n > 65536 {
        v.sort_unstable_by_key(|&h| {
            let h = h as u32;
            ((h ^ (h >> 16)) as usize) & (cap - 1)
        });
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn loc(tile: i32, x: i32, y: i32) -> OpticalEnd {
        OpticalEnd { read_group: 0, loc: Location { tile, x, y } }
    }

    #[test]
    fn default_name_parsing() {
        let p = NameParser::new(Some(crate::args::DEFAULT_READ_NAME_REGEX)).unwrap();
        // 5 fields (Illumina before CASAVA 1.8) and 7 fields
        assert_eq!(p.location(b"M1:1:1101:2000:3000"), Some(Location { tile: 1101, x: 2000, y: 3000 }));
        assert_eq!(p.location(b"A00:12:FC:1:2207:15000:1000"), Some(Location { tile: 2207, x: 15000, y: 1000 }));
        // trailing non-digits after a number are ignored
        assert_eq!(p.location(b"M1:1:1101:2000:3000#0/1"), Some(Location { tile: 1101, x: 2000, y: 3000 }));
        assert_eq!(p.location(b"read1"), None);
        assert_eq!(p.location(b"a:b:c:d"), None);
        assert_eq!(p.location(b"a:b:1101:x:3000"), None);
    }

    /// Picard's ReadNameParserTests: names with 5 and 7 fields, and the getLastThreeFields field counts
    #[test]
    fn picard_last_three_fields() {
        assert_eq!(last_three_fields(b"RUNID:7:1203:2886:82292").unwrap(), (5, [1203, 2886, 82292]));
        assert_eq!(last_three_fields(b"RUNID:7:1203:2884:16834").unwrap(), (5, [1203, 2884, 16834]));
        assert_eq!(last_three_fields(b"1:2:3").unwrap(), (3, [1, 2, 3]));
        assert_eq!(last_three_fields(b"1:2").unwrap(), (-1, [-1, -1, -1]));
        assert!(last_three_fields(b"1:2:x").is_none());
    }

    #[test]
    fn short_coordinates_wrap_in_java_mode() {
        compat::with_java(true, || {
            let p = NameParser::Default;
            assert_eq!(p.location(b"a:1:1101:40000:70000"), Some(Location { tile: 1101, x: 40000 - 65536, y: 70000 - 65536 }));
        });
        compat::with_java(false, || {
            let p = NameParser::Default;
            assert_eq!(p.location(b"a:1:1101:40000:70000"), Some(Location { tile: 1101, x: 40000, y: 70000 }));
        });
    }

    /// Picard's OpticalDuplicateFinderTest.testKeeperAtEndWithinCliqueOfAllOpticalDuplicates and friends
    #[test]
    fn keeper_and_chains() {
        let f = OpticalFinder { distance: 100, max_set_size: 300000 };
        // two reads within the distance, the keeper first
        let l = [loc(1, 10, 10), loc(1, 50, 50)];
        assert_eq!(f.find(&l, Some(0)), vec![false, true]);
        // a chain 0-1-2 without the keeper close: transitive through the graph
        let l = [loc(1, 0, 0), loc(1, 90, 0), loc(1, 180, 0), loc(2, 0, 0)];
        assert_eq!(f.find(&l, Some(3)), vec![false, true, true, false]);
        // different tiles never match
        let l = [loc(1, 0, 0), loc(2, 0, 0)];
        assert_eq!(f.find(&l, Some(0)), vec![false, false]);
        // no location on the keeper: pairs among the others still count
        let l = [OpticalEnd { read_group: -1, loc: NO_LOCATION }, loc(1, 0, 0), loc(1, 5, 5)];
        assert_eq!(f.find(&l, Some(0)), vec![false, false, true]);
    }

    /// Picard's OpticalDuplicateFinderTest.testOpticalDuplicateClusterOneEndSamePositionOneCluster: with the
    /// keeper outside the cluster, the read with the lowest x stays
    #[test]
    fn graph_keeps_lowest_x() {
        let f = OpticalFinder { distance: 100, max_set_size: 300000 };
        let l = [loc(1, 500, 500), loc(1, 450, 520), loc(1, 520, 480), loc(5, 0, 0)];
        assert_eq!(f.find(&l, Some(3)), vec![true, false, true, false]);
    }

    #[test]
    fn max_set_size() {
        let l = [loc(1, 0, 0), loc(1, 5, 5)];
        let f = OpticalFinder { distance: 100, max_set_size: -1 };
        compat::with_java(true, || assert_eq!(f.find(&l, Some(0)), vec![false, false]));
        compat::with_java(false, || assert_eq!(f.find(&l, Some(0)), vec![false, true]));
        let f = OpticalFinder { distance: 100, max_set_size: 1 };
        assert_eq!(f.find(&l, Some(0)), vec![false, false]);
    }

    #[test]
    fn integer_map_order() {
        assert_eq!(java_integer_map_order(5), vec![0, 1, 2, 3, 4]);
        let o = java_integer_map_order(70000);
        assert_eq!(o.len(), 70000);
        // 65536 ^ 1 == 65537: key 65537 comes before 65536 in a table of 131072
        let p65536 = o.iter().position(|&k| k == 65536).unwrap();
        let p65537 = o.iter().position(|&k| k == 65537).unwrap();
        assert!(p65537 < p65536);
    }
}
