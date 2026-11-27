//! java.util.HashSet / HashMap with Java's iteration order.
//!
//! Java iterates a HashMap by bucket, and entries of one bucket in insertion order; the bucket of a key comes
//! from its hashCode and the current table size. Where the Java code iterates such a collection and the order
//! can reach the output, these types give the same order. (Buckets that Java would convert to trees, eight or
//! more colliding keys in a table of 64 or more, are not modelled.)

use std::collections::HashMap;
use std::hash::Hash;

pub trait JavaHash {
    fn java_hash(&self) -> i32;
}

impl JavaHash for String {
    fn java_hash(&self) -> i32 {
        java_string_hash(self)
    }
}
impl JavaHash for &str {
    fn java_hash(&self) -> i32 {
        java_string_hash(self)
    }
}
impl JavaHash for i32 {
    fn java_hash(&self) -> i32 {
        *self
    }
}

/// String.hashCode
pub fn java_string_hash(s: &str) -> i32 {
    let mut h: i32 = 0;
    for c in s.encode_utf16() {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h
}

#[inline]
fn spread(h: i32) -> u32 {
    let h = h as u32;
    h ^ (h >> 16)
}

fn table_size_for(cap: usize) -> usize {
    cap.max(1).next_power_of_two()
}

/// Guava Maps.capacity(expectedSize), used by Sets.newHashSet(elements...) and newHashSetWithExpectedSize
pub fn guava_capacity(expected_size: usize) -> usize {
    if expected_size < 3 { expected_size + 1 } else { (expected_size as f64 / 0.75 + 1.0) as usize }
}

#[derive(Clone, Debug)]
pub struct JavaHashMap<K, V> {
    /// insertion order, with removed entries taken out
    entries: Vec<(K, V)>,
    /// position within a bucket's chain: put appends (increasing numbers), computeIfAbsent / compute /
    /// merge put a new key at the head of its bucket (decreasing negative numbers)
    seqs: Vec<i64>,
    next_tail: i64,
    next_head: i64,
    index: HashMap<K, usize>,
    capacity: usize,
}

impl<K: JavaHash + Eq + Hash + Clone, V> Default for JavaHashMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: JavaHash + Eq + Hash + Clone, V> JavaHashMap<K, V> {
    /// new HashMap<>()
    pub fn new() -> Self {
        Self::with_capacity(16)
    }
    /// new HashMap<>(initialCapacity)
    pub fn with_capacity(initial_capacity: usize) -> Self {
        JavaHashMap { entries: Vec::new(), seqs: Vec::new(), next_tail: 0, next_head: -1, index: HashMap::new(), capacity: table_size_for(initial_capacity) }
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    pub fn contains_key(&self, key: &K) -> bool {
        self.index.contains_key(key)
    }
    pub fn get(&self, key: &K) -> Option<&V> {
        self.index.get(key).map(|i| &self.entries[*i].1)
    }
    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        match self.index.get(key) {
            Some(i) => Some(&mut self.entries[*i].1),
            None => None,
        }
    }
    /// put: returns the previous value
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if let Some(i) = self.index.get(&key) {
            return Some(std::mem::replace(&mut self.entries[*i].1, value));
        }
        self.push(key, value, false);
        None
    }
    fn push(&mut self, key: K, value: V, head: bool) {
        self.index.insert(key.clone(), self.entries.len());
        self.entries.push((key, value));
        if head {
            self.seqs.push(self.next_head);
            self.next_head -= 1;
        } else {
            self.seqs.push(self.next_tail);
            self.next_tail += 1;
        }
        // resize when the size passes 0.75 of the table (the split keeps the chain order)
        while self.entries.len() > (self.capacity * 3) / 4 {
            self.capacity *= 2;
        }
    }
    /// computeIfAbsent (and the new-key case of compute / merge, as Collectors.groupingBy and toMap use):
    /// a new key goes to the head of its bucket
    pub fn compute_if_absent(&mut self, key: K, f: impl FnOnce() -> V) -> &mut V {
        let i = match self.index.get(&key) {
            Some(i) => *i,
            None => {
                self.push(key, f(), true);
                self.entries.len() - 1
            }
        };
        &mut self.entries[i].1
    }
    pub fn remove(&mut self, key: &K) -> Option<V> {
        let i = self.index.remove(key)?;
        let (_, v) = self.entries.remove(i);
        self.seqs.remove(i);
        for (_, idx) in self.index.iter_mut() {
            if *idx > i {
                *idx -= 1;
            }
        }
        Some(v)
    }
    pub fn clear(&mut self) {
        self.entries.clear();
        self.seqs.clear();
        self.index.clear();
    }
    /// keeps the entries `f` accepts, in their order; the table keeps its size, as Java's does on removal
    /// (Collection.retainAll / removeIf)
    pub fn retain(&mut self, mut f: impl FnMut(&K, &V) -> bool) {
        let keep: Vec<bool> = self.entries.iter().map(|(k, v)| f(k, v)).collect();
        let mut it = keep.iter();
        self.entries.retain(|_| *it.next().unwrap());
        let mut it = keep.iter();
        self.seqs.retain(|_| *it.next().unwrap());
        self.index.clear();
        for (i, (k, _)) in self.entries.iter().enumerate() {
            self.index.insert(k.clone(), i);
        }
    }
    fn order(&self) -> Vec<usize> {
        let mask = (self.capacity - 1) as u32;
        let mut order: Vec<usize> = (0..self.entries.len()).collect();
        order.sort_by_key(|i| (spread(self.entries[*i].0.java_hash()) & mask, self.seqs[*i]));
        order
    }
    /// entries in Java iteration order
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.order().into_iter().map(move |i| (&self.entries[i].0, &self.entries[i].1))
    }
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.iter().map(|(k, _)| k)
    }
    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.iter().map(|(_, v)| v)
    }
    /// entries in insertion order (for uses where order cannot matter)
    pub fn iter_unordered(&self) -> impl Iterator<Item = &(K, V)> {
        self.entries.iter()
    }
}

#[derive(Clone, Debug)]
pub struct JavaHashSet<K> {
    map: JavaHashMap<K, ()>,
}

impl<K: JavaHash + Eq + Hash + Clone> Default for JavaHashSet<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: JavaHash + Eq + Hash + Clone> JavaHashSet<K> {
    /// new HashSet<>() / Sets.newHashSet()
    pub fn new() -> Self {
        JavaHashSet { map: JavaHashMap::new() }
    }
    pub fn with_capacity(initial_capacity: usize) -> Self {
        JavaHashSet { map: JavaHashMap::with_capacity(initial_capacity) }
    }
    /// Sets.newHashSet(element)
    pub fn of_one(key: K) -> Self {
        let mut s = Self::with_capacity(guava_capacity(1));
        s.insert(key);
        s
    }
    /// new HashSet<>(collection) / Sets.newHashSet(collection)
    pub fn from_collection<I: IntoIterator<Item = K>>(items: I) -> Self {
        let items: Vec<K> = items.into_iter().collect();
        let mut s = Self::with_capacity(((items.len() as f32 / 0.75) as usize + 1).max(16));
        for k in items {
            s.insert(k);
        }
        s
    }
    pub fn len(&self) -> usize {
        self.map.len()
    }
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
    pub fn contains(&self, key: &K) -> bool {
        self.map.contains_key(key)
    }
    /// add: true if the key was not present
    pub fn insert(&mut self, key: K) -> bool {
        self.map.insert(key, ()).is_none()
    }
    pub fn remove(&mut self, key: &K) -> bool {
        self.map.remove(key).is_some()
    }
    pub fn clear(&mut self) {
        self.map.clear()
    }
    /// retainAll / removeIf: keeps the keys `f` accepts; the table keeps its size
    pub fn retain(&mut self, mut f: impl FnMut(&K) -> bool) {
        self.map.retain(|k, _| f(k))
    }
    /// keys in Java iteration order
    pub fn iter(&self) -> impl Iterator<Item = &K> {
        self.map.keys()
    }
    pub fn iter_unordered(&self) -> impl Iterator<Item = &K> {
        self.map.iter_unordered().map(|(k, _)| k)
    }
    pub fn to_vec(&self) -> Vec<K> {
        self.iter().cloned().collect()
    }
}

/// Bucket of a key with hash code `h` in a HashMap table of `cap` buckets (a power of two)
#[inline]
pub fn bucket_index(h: i32, cap: usize) -> usize {
    (spread(h) as usize) & (cap - 1)
}

/// HashMap.tableSizeFor: the smallest power of two at or above `cap` (1 for 0)
pub fn java_table_size_for(cap: usize) -> usize {
    table_size_for(cap)
}

/// The table size after `n` distinct keys were put into a HashMap whose table started at `cap` buckets:
/// Java doubles the table whenever the size passes `(int)(cap * 0.75f)`.
pub fn grown_capacity(cap: usize, n: usize) -> usize {
    let mut cap = cap.max(1);
    while n > (cap as f32 * 0.75f32) as usize {
        cap *= 2;
    }
    cap
}

/// The table size of `new HashSet<>(collection)` / `Sets.newHashSet(collection)` with `n` elements on
/// JDK 19 and later (`HashMap.newHashMap(Math.max(n, 12))`, that is `ceil(max(n, 12) / 0.75)` rounded up
/// to a power of two). JDK 8 to 18 used `max((int)(n / .75f) + 1, 16)`, which `JavaHashSet::from_collection`
/// models; the two differ when `n / 0.75` is a power of two (12, 24, 48, ...).
pub fn jdk19_collection_capacity(n: usize) -> usize {
    let m = n.max(12) as f64;
    table_size_for((m / 0.75).ceil() as usize)
}

/// Guava 33 `Maps.capacity(expectedSize)` (used by `Sets.newHashSet(E...)`, `newHashSetWithExpectedSize`,
/// `Maps.newHashMapWithExpectedSize`): `n + 1` below 3, else `ceil(n / 0.75)`. Older Guava versions used
/// `(int)(n / 0.75f + 1.0f)` (`guava_capacity`).
pub fn guava33_capacity(expected_size: usize) -> usize {
    if expected_size < 3 { expected_size + 1 } else { (expected_size as f64 / 0.75).ceil() as usize }
}

/// Java iteration order of the keys of a HashMap or HashSet whose table has `cap` buckets: `hashes` are the
/// hash codes of the keys in first-insertion order (removed keys left out); returns their indices in the
/// order Java iterates them (by bucket, and in insertion order within a bucket). Buckets converted to trees
/// (eight or more colliding keys in a table of 64 or more) are not modelled.
pub fn java_iteration_order(hashes: &[i32], cap: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..hashes.len()).collect();
    if cap > 1 {
        order.sort_by_key(|&i| bucket_index(hashes[i], cap));
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn string_hash() {
        assert_eq!(java_string_hash(""), 0);
        assert_eq!(java_string_hash("a"), 97);
        assert_eq!(java_string_hash("hello"), 99162322);
        assert_eq!(java_string_hash("chr1_12345"), {
            let mut h = 0i32;
            for b in "chr1_12345".bytes() {
                h = h.wrapping_mul(31).wrapping_add(b as i32);
            }
            h
        });
    }

    #[test]
    fn iteration_order() {
        // new HashSet<String>() with "b", "a", "q" (hashes 98, 97, 113 -> buckets 2, 1, 1 of 16)
        let mut s: JavaHashSet<String> = JavaHashSet::new();
        for k in ["b", "a", "q"] {
            s.insert(k.to_string());
        }
        assert_eq!(s.to_vec(), vec!["a", "q", "b"]);
        // integers iterate in bucket order
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for k in [17, 3, 1, 33] {
            m.insert(k, k);
        }
        assert_eq!(m.keys().cloned().collect::<Vec<_>>(), vec![17, 1, 33, 3]);
        // Collectors.groupingBy (computeIfAbsent) puts colliding keys at the head of the bucket: Java
        // prints [SUZ12, SMAD4, NF1, TP53, ARID1B] for these keys grouped in this order (TP53 and ARID1B
        // share bucket 13 of 16), and [SUZ12, SMAD4, NF1, ARID1B, TP53] for HashMap.put
        let keys = ["ARID1B", "TP53", "NF1", "SUZ12", "SMAD4"];
        let mut g: JavaHashMap<String, i32> = JavaHashMap::new();
        let mut p: JavaHashMap<String, i32> = JavaHashMap::new();
        for k in keys {
            *g.compute_if_absent(k.to_string(), || 0) += 1;
            p.insert(k.to_string(), 1);
        }
        assert_eq!(g.keys().cloned().collect::<Vec<_>>(), vec!["SUZ12", "SMAD4", "NF1", "TP53", "ARID1B"]);
        assert_eq!(p.keys().cloned().collect::<Vec<_>>(), vec!["SUZ12", "SMAD4", "NF1", "ARID1B", "TP53"]);
    }

    #[test]
    fn capacities() {
        assert_eq!(jdk19_collection_capacity(0), 16);
        assert_eq!(jdk19_collection_capacity(12), 16);
        assert_eq!(jdk19_collection_capacity(13), 32);
        assert_eq!(jdk19_collection_capacity(24), 32);
        assert_eq!(jdk19_collection_capacity(25), 64);
        assert_eq!(guava33_capacity(1), 2);
        assert_eq!(guava33_capacity(3), 4);
        assert_eq!(guava33_capacity(6), 8);
        assert_eq!(grown_capacity(16, 12), 16);
        assert_eq!(grown_capacity(16, 13), 32);
        assert_eq!(grown_capacity(1, 1), 2);
        assert_eq!(grown_capacity(2, 2), 4);
        // bucket order, insertion order within a bucket
        assert_eq!(java_iteration_order(&[98, 97, 113], 16), vec![1, 2, 0]);
    }
}
