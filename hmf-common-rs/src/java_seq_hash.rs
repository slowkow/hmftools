//! java.util.HashMap with Java's iteration order, for large maps that are changed often.
//!
//! `java_hash::JavaHashMap` keeps its entries in a vector, so a removal costs a pass over the map. This map
//! keeps them in a Rust HashMap with the insertion sequence of each, so insertion, lookup and removal are
//! constant time, and works out the Java order when the map is iterated: by bucket (the spread hash code
//! masked by the table size, which grows as Java's does and is kept by `clear` and `remove`), then by
//! insertion within the bucket. As in `JavaHashMap`, buckets Java would convert to trees (eight or more
//! colliding keys in a table of 64 or more) are not modelled.
//!
//! ISOFOX keys its read caches by read name (FragmentTracker, the chimeric read groups), with a removal for
//! every completed pair.

use std::collections::HashMap;
use std::hash::Hash;

use crate::java_hash::JavaHash;

#[inline]
fn spread(h: i32) -> u32 {
    let h = h as u32;
    h ^ (h >> 16)
}

#[derive(Clone, Debug)]
struct Slot<V> {
    value: V,
    seq: u64,
    spread: u32,
}

#[derive(Clone, Debug)]
pub struct SeqHashMap<K, V> {
    map: HashMap<K, Slot<V>>,
    next_seq: u64,
    /// table length, 0 until the first insertion (new HashMap<>() allocates 16 on the first put)
    capacity: usize,
    /// the table length of the first allocation
    initial_capacity: usize,
}

impl<K: JavaHash + Eq + Hash + Clone, V> Default for SeqHashMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: JavaHash + Eq + Hash + Clone, V> SeqHashMap<K, V> {
    /// new HashMap<>()
    pub fn new() -> Self {
        SeqHashMap { map: HashMap::new(), next_seq: 0, capacity: 0, initial_capacity: 16 }
    }
    /// new HashMap<>(initialCapacity)
    pub fn with_capacity(initial_capacity: usize) -> Self {
        SeqHashMap { map: HashMap::new(), next_seq: 0, capacity: 0, initial_capacity: initial_capacity.max(1).next_power_of_two() }
    }
    pub fn len(&self) -> usize {
        self.map.len()
    }
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
    /// the current table length (0 before the first insertion)
    pub fn table_len(&self) -> usize {
        self.capacity
    }
    pub fn contains_key(&self, key: &K) -> bool {
        self.map.contains_key(key)
    }
    pub fn get(&self, key: &K) -> Option<&V> {
        self.map.get(key).map(|s| &s.value)
    }
    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        self.map.get_mut(key).map(|s| &mut s.value)
    }
    /// put: returns the previous value (an existing key keeps its place)
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if let Some(slot) = self.map.get_mut(&key) {
            return Some(std::mem::replace(&mut slot.value, value));
        }
        if self.capacity == 0 {
            self.capacity = self.initial_capacity;
        }
        let spread = spread(key.java_hash());
        self.map.insert(key, Slot { value, seq: self.next_seq, spread });
        self.next_seq += 1;
        while self.map.len() > (self.capacity * 3) / 4 {
            self.capacity *= 2;
        }
        None
    }
    pub fn remove(&mut self, key: &K) -> Option<V> {
        self.map.remove(key).map(|s| s.value)
    }
    /// clear: the table keeps its size, as Java's does
    pub fn clear(&mut self) {
        self.map.clear();
    }
    fn ordered(&self) -> Vec<(&K, &Slot<V>)> {
        let mask = (self.capacity.max(1) - 1) as u32;
        let mut entries: Vec<(&K, &Slot<V>)> = self.map.iter().collect();
        entries.sort_unstable_by_key(|(_, s)| (s.spread & mask, s.seq));
        entries
    }
    /// entries in Java iteration order
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.ordered().into_iter().map(|(k, s)| (k, &s.value))
    }
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.ordered().into_iter().map(|(k, _)| k)
    }
    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.ordered().into_iter().map(|(_, s)| &s.value)
    }
    /// the keys in Java iteration order, as owned values (to change the map while going through them)
    pub fn keys_vec(&self) -> Vec<K> {
        self.ordered().into_iter().map(|(k, _)| k.clone()).collect()
    }
    /// entries in no particular order (where the order cannot matter)
    pub fn iter_unordered(&self) -> impl Iterator<Item = (&K, &V)> {
        self.map.iter().map(|(k, s)| (k, &s.value))
    }
    pub fn values_unordered_mut(&mut self) -> impl Iterator<Item = &mut V> {
        self.map.values_mut().map(|s| &mut s.value)
    }
    /// takes every entry out, in Java iteration order (the table keeps its size)
    pub fn drain_ordered(&mut self) -> Vec<(K, V)> {
        let mask = (self.capacity.max(1) - 1) as u32;
        let mut entries: Vec<(K, Slot<V>)> = self.map.drain().collect();
        entries.sort_unstable_by_key(|(_, s)| (s.spread & mask, s.seq));
        entries.into_iter().map(|(k, s)| (k, s.value)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::java_hash::JavaHashMap;

    /// The same operations on a JavaHashMap (the reference model of Java's order) and a SeqHashMap give the
    /// same iteration orders.
    #[test]
    fn same_order_as_java_hash_map() {
        let mut seed: u64 = 12345;
        let mut a: JavaHashMap<String, usize> = JavaHashMap::new();
        let mut b: SeqHashMap<String, usize> = SeqHashMap::new();
        let mut live: Vec<String> = Vec::new();
        for i in 0..3000usize {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let r = (seed >> 33) % 100;
            if r < 65 || live.is_empty() {
                let k = format!("R{}:{}", (seed >> 40) % 997, i);
                a.insert(k.clone(), i);
                b.insert(k.clone(), i);
                live.push(k);
            } else {
                let k = live.remove((r as usize) % live.len());
                assert_eq!(a.remove(&k), b.remove(&k));
            }
            if i % 97 == 0 {
                let x: Vec<&String> = a.keys().collect();
                let y: Vec<&String> = b.keys().collect();
                assert_eq!(x, y, "step {i}");
            }
            if i == 1500 {
                a.clear();
                b.clear();
                live.clear();
            }
        }
    }

    #[test]
    fn integer_keys_and_capacity() {
        let mut m: SeqHashMap<i32, i32> = SeqHashMap::new();
        assert_eq!(m.table_len(), 0);
        for k in [33, 1, 17, 65538, 2] {
            m.insert(k, k);
        }
        assert_eq!(m.table_len(), 16);
        // buckets: 33 -> 1, 1 -> 1, 17 -> 1, 65538 (spread 65539) -> 3, 2 -> 2
        let order: Vec<i32> = m.keys().copied().collect();
        assert_eq!(order, vec![33, 1, 17, 2, 65538]);
        m.clear();
        assert_eq!(m.table_len(), 16);
    }
}
