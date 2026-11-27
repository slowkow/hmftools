//! Maps and sets keyed by object identity (Java's HashMap and HashSet of classes without hashCode), whose
//! iteration order reaches LINX's output (LX1).
//!
//! Java iterates them in identity-hash order, which differs between JVM settings: under `-java_compat`
//! these iterate as Java does with `-XX:hashCode=2` (every identity hash 1; hmf_common::java_id_hash), the
//! reference the port reproduces; by default they iterate in insertion order, a defined order.

use std::collections::HashMap;

use hmf_common::java_id_hash::{IdHashMap, IdKey};

pub enum IdMap<K: IdKey, V> {
    Java(IdHashMap<K, V>),
    Ordered { entries: Vec<Option<(K, V)>>, index: HashMap<usize, usize>, len: usize },
}

impl<K: IdKey, V> Default for IdMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: IdKey, V> IdMap<K, V> {
    pub fn new() -> Self {
        if crate::compat::java("LX1") {
            IdMap::Java(IdHashMap::new())
        } else {
            IdMap::Ordered { entries: Vec::new(), index: HashMap::new(), len: 0 }
        }
    }
    pub fn len(&self) -> usize {
        match self {
            IdMap::Java(m) => m.len(),
            IdMap::Ordered { len, .. } => *len,
        }
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn contains_key(&self, k: &K) -> bool {
        match self {
            IdMap::Java(m) => m.contains_key(k),
            IdMap::Ordered { index, .. } => index.contains_key(&k.id_addr()),
        }
    }
    pub fn get(&self, k: &K) -> Option<&V> {
        match self {
            IdMap::Java(m) => m.get(k),
            IdMap::Ordered { entries, index, .. } => index.get(&k.id_addr()).and_then(|i| entries[*i].as_ref().map(|e| &e.1)),
        }
    }
    pub fn get_mut(&mut self, k: &K) -> Option<&mut V> {
        match self {
            IdMap::Java(m) => m.get_mut(k),
            IdMap::Ordered { entries, index, .. } => match index.get(&k.id_addr()) {
                Some(i) => entries[*i].as_mut().map(|e| &mut e.1),
                None => None,
            },
        }
    }
    /// put: an existing key keeps its place
    pub fn insert(&mut self, k: K, v: V) -> Option<V> {
        match self {
            IdMap::Java(m) => m.insert(k, v),
            IdMap::Ordered { entries, index, len } => {
                if let Some(i) = index.get(&k.id_addr()) {
                    let e = entries[*i].as_mut().unwrap();
                    return Some(std::mem::replace(&mut e.1, v));
                }
                index.insert(k.id_addr(), entries.len());
                entries.push(Some((k, v)));
                *len += 1;
                None
            }
        }
    }
    pub fn remove(&mut self, k: &K) -> Option<V> {
        match self {
            IdMap::Java(m) => m.remove(k),
            IdMap::Ordered { entries, index, len } => {
                let i = index.remove(&k.id_addr())?;
                *len -= 1;
                entries[i].take().map(|e| e.1)
            }
        }
    }
    pub fn clear(&mut self) {
        match self {
            IdMap::Java(m) => m.clear(),
            IdMap::Ordered { entries, index, len } => {
                entries.clear();
                index.clear();
                *len = 0;
            }
        }
    }
    pub fn iter(&self) -> Box<dyn Iterator<Item = (&K, &V)> + '_> {
        match self {
            IdMap::Java(m) => Box::new(m.iter()),
            IdMap::Ordered { entries, .. } => Box::new(entries.iter().filter_map(|e| e.as_ref().map(|(k, v)| (k, v)))),
        }
    }
    pub fn keys(&self) -> Vec<K> {
        self.iter().map(|(k, _)| k.clone()).collect()
    }
}

/// A HashSet keyed by identity
pub struct IdSet<K: IdKey>(IdMap<K, ()>);

impl<K: IdKey> Default for IdSet<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: IdKey> IdSet<K> {
    pub fn new() -> Self {
        IdSet(IdMap::new())
    }
    pub fn insert(&mut self, k: K) -> bool {
        self.0.insert(k, ()).is_none()
    }
    pub fn contains(&self, k: &K) -> bool {
        self.0.contains_key(k)
    }
    pub fn remove(&mut self, k: &K) -> bool {
        self.0.remove(k).is_some()
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub fn to_vec(&self) -> Vec<K> {
        self.0.keys()
    }
}
