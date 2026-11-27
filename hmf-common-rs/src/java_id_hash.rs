//! java.util.HashSet / HashMap for keys that use Object.hashCode (identity hash), as iterated by a JVM where
//! every identity hash is the same value (`-XX:hashCode=2`).
//!
//! With the default JVM settings the iteration order of such a set differs from run to run, so any order is a
//! possible Java result. With constant identity hashes all keys share one bucket and the order is fixed by
//! HashMap's bucket handling: a linked list in insertion order, converted to a red-black tree (whose list links
//! are reordered when the root moves to the front) once the bucket holds more than eight entries in a table of
//! 64 or more. This is a port of that handling, so iteration here equals that of the reference JVM.

use std::collections::HashMap;

/// A key compared by identity, with a stable address.
pub trait IdKey: Clone {
    fn id_addr(&self) -> usize;
}

const TREEIFY_THRESHOLD: usize = 8;
const UNTREEIFY_THRESHOLD: usize = 6;
const MIN_TREEIFY_CAPACITY: usize = 64;
const NIL: usize = usize::MAX;

#[derive(Clone)]
struct Node<K, V> {
    key: K,
    value: V,
    parent: usize,
    left: usize,
    right: usize,
    prev: usize,
    next: usize,
    red: bool,
}

#[derive(Clone)]
pub struct IdHashMap<K, V> {
    nodes: Vec<Option<Node<K, V>>>,
    free: Vec<usize>,
    lookup: HashMap<usize, usize>,
    /// first node of the bucket (the root when the bucket is a tree)
    first: usize,
    is_tree: bool,
    size: usize,
    /// table length, zero until the first insertion
    table_len: usize,
    threshold: usize,
}

fn table_size_for(cap: usize) -> usize {
    cap.max(1).next_power_of_two()
}

impl<K: IdKey, V> Default for IdHashMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: IdKey, V> IdHashMap<K, V> {
    /// new HashMap<>()
    pub fn new() -> Self {
        IdHashMap { nodes: Vec::new(), free: Vec::new(), lookup: HashMap::new(), first: NIL, is_tree: false, size: 0, table_len: 0, threshold: 0 }
    }
    /// new HashMap<>(initialCapacity)
    pub fn with_capacity(initial_capacity: usize) -> Self {
        let mut m = Self::new();
        m.threshold = table_size_for(initial_capacity);
        m
    }
    pub fn len(&self) -> usize {
        self.size
    }
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
    pub fn contains_key(&self, key: &K) -> bool {
        self.lookup.contains_key(&key.id_addr())
    }
    pub fn get(&self, key: &K) -> Option<&V> {
        self.lookup.get(&key.id_addr()).map(|i| &self.n(*i).value)
    }
    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        match self.lookup.get(&key.id_addr()) {
            Some(i) => Some(&mut self.nodes[*i].as_mut().unwrap().value),
            None => None,
        }
    }
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.free.clear();
        self.lookup.clear();
        self.first = NIL;
        self.is_tree = false;
        self.size = 0;
    }

    fn n(&self, i: usize) -> &Node<K, V> {
        self.nodes[i].as_ref().unwrap()
    }
    fn nm(&mut self, i: usize) -> &mut Node<K, V> {
        self.nodes[i].as_mut().unwrap()
    }
    fn new_node(&mut self, key: K, value: V) -> usize {
        let node = Node { key, value, parent: NIL, left: NIL, right: NIL, prev: NIL, next: NIL, red: false };
        if let Some(i) = self.free.pop() {
            self.nodes[i] = Some(node);
            i
        } else {
            self.nodes.push(Some(node));
            self.nodes.len() - 1
        }
    }

    fn resize(&mut self) {
        let old_cap = self.table_len;
        let old_thr = self.threshold;
        let new_cap;
        let mut new_thr = 0;
        if old_cap > 0 {
            new_cap = old_cap << 1;
            if old_cap >= 16 {
                new_thr = old_thr << 1;
            }
        } else if old_thr > 0 {
            new_cap = old_thr;
        } else {
            new_cap = 16;
            new_thr = 12;
        }
        if new_thr == 0 {
            new_thr = (new_cap as f32 * 0.75) as usize;
        }
        self.threshold = new_thr;
        self.table_len = new_cap;
        // TreeNode.split: all entries stay in one bucket; a tree of six or fewer becomes a list again
        if old_cap > 0 && self.is_tree && self.first != NIL && self.size <= UNTREEIFY_THRESHOLD {
            self.untreeify();
        }
    }

    fn untreeify(&mut self) {
        self.is_tree = false;
        let mut i = self.first;
        while i != NIL {
            let next = self.n(i).next;
            let node = self.nm(i);
            node.parent = NIL;
            node.left = NIL;
            node.right = NIL;
            node.prev = NIL;
            node.red = false;
            i = next;
        }
    }

    /// put: returns the previous value
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if let Some(i) = self.lookup.get(&key.id_addr()) {
            return Some(std::mem::replace(&mut self.nodes[*i].as_mut().unwrap().value, value));
        }
        if self.table_len == 0 {
            self.resize();
        }
        let addr = key.id_addr();
        let x = self.new_node(key, value);
        self.lookup.insert(addr, x);
        if self.first == NIL {
            self.first = x;
        } else if self.is_tree {
            self.put_tree_val(x);
        } else {
            // append to the list; binCount is the number of nodes before the last one
            let mut p = self.first;
            let mut bin_count = 0;
            while self.n(p).next != NIL {
                p = self.n(p).next;
                bin_count += 1;
            }
            self.nm(p).next = x;
            if bin_count >= TREEIFY_THRESHOLD - 1 {
                self.treeify_bin();
            }
        }
        self.size += 1;
        if self.size > self.threshold {
            self.resize();
        }
        None
    }

    fn treeify_bin(&mut self) {
        if self.table_len < MIN_TREEIFY_CAPACITY {
            self.resize();
            return;
        }
        // replace the list nodes by tree nodes in order, then build the tree
        let mut prev = NIL;
        let mut i = self.first;
        while i != NIL {
            self.nm(i).prev = prev;
            prev = i;
            i = self.n(i).next;
        }
        self.is_tree = true;
        self.treeify();
    }

    fn treeify(&mut self) {
        let mut root = NIL;
        let mut x = self.first;
        while x != NIL {
            let next = self.n(x).next;
            self.nm(x).left = NIL;
            self.nm(x).right = NIL;
            if root == NIL {
                self.nm(x).parent = NIL;
                self.nm(x).red = false;
                root = x;
            } else {
                // equal hashes and no ordering between keys: tieBreakOrder gives -1, so always to the left
                let mut p = root;
                loop {
                    let xp = p;
                    p = self.n(p).left;
                    if p == NIL {
                        self.nm(x).parent = xp;
                        self.nm(xp).left = x;
                        root = self.balance_insertion(root, x);
                        break;
                    }
                }
            }
            x = next;
        }
        self.move_root_to_front(root);
    }

    fn move_root_to_front(&mut self, root: usize) {
        if root == NIL {
            return;
        }
        let first = self.first;
        if root != first {
            self.first = root;
            let rp = self.n(root).prev;
            let rn = self.n(root).next;
            if rn != NIL {
                self.nm(rn).prev = rp;
            }
            if rp != NIL {
                self.nm(rp).next = rn;
            }
            if first != NIL {
                self.nm(first).prev = root;
            }
            self.nm(root).next = first;
            self.nm(root).prev = NIL;
        }
    }

    fn root_of(&self, mut r: usize) -> usize {
        loop {
            let p = self.n(r).parent;
            if p == NIL {
                return r;
            }
            r = p;
        }
    }

    fn put_tree_val(&mut self, x: usize) {
        let root = self.root_of(self.first);
        let mut p = root;
        loop {
            let xp = p;
            p = self.n(p).left;
            if p == NIL {
                let xpn = self.n(xp).next;
                self.nm(x).next = xpn;
                self.nm(xp).left = x;
                self.nm(xp).next = x;
                self.nm(x).parent = xp;
                self.nm(x).prev = xp;
                if xpn != NIL {
                    self.nm(xpn).prev = x;
                }
                let r = self.balance_insertion(root, x);
                self.move_root_to_front(r);
                return;
            }
        }
    }

    fn rotate_left(&mut self, mut root: usize, p: usize) -> usize {
        if p != NIL {
            let r = self.n(p).right;
            if r != NIL {
                let rl = self.n(r).left;
                self.nm(p).right = rl;
                if rl != NIL {
                    self.nm(rl).parent = p;
                }
                let pp = self.n(p).parent;
                self.nm(r).parent = pp;
                if pp == NIL {
                    root = r;
                    self.nm(r).red = false;
                } else if self.n(pp).left == p {
                    self.nm(pp).left = r;
                } else {
                    self.nm(pp).right = r;
                }
                self.nm(r).left = p;
                self.nm(p).parent = r;
            }
        }
        root
    }

    fn rotate_right(&mut self, mut root: usize, p: usize) -> usize {
        if p != NIL {
            let l = self.n(p).left;
            if l != NIL {
                let lr = self.n(l).right;
                self.nm(p).left = lr;
                if lr != NIL {
                    self.nm(lr).parent = p;
                }
                let pp = self.n(p).parent;
                self.nm(l).parent = pp;
                if pp == NIL {
                    root = l;
                    self.nm(l).red = false;
                } else if self.n(pp).right == p {
                    self.nm(pp).right = l;
                } else {
                    self.nm(pp).left = l;
                }
                self.nm(l).right = p;
                self.nm(p).parent = l;
            }
        }
        root
    }

    fn red(&self, i: usize) -> bool {
        i != NIL && self.n(i).red
    }

    fn balance_insertion(&mut self, mut root: usize, mut x: usize) -> usize {
        self.nm(x).red = true;
        loop {
            let mut xp = self.n(x).parent;
            if xp == NIL {
                self.nm(x).red = false;
                return x;
            }
            let mut xpp = self.n(xp).parent;
            if !self.n(xp).red || xpp == NIL {
                return root;
            }
            let xppl = self.n(xpp).left;
            if xp == xppl {
                let xppr = self.n(xpp).right;
                if xppr != NIL && self.n(xppr).red {
                    self.nm(xppr).red = false;
                    self.nm(xp).red = false;
                    self.nm(xpp).red = true;
                    x = xpp;
                } else {
                    if x == self.n(xp).right {
                        x = xp;
                        root = self.rotate_left(root, x);
                        xp = self.n(x).parent;
                        xpp = if xp == NIL { NIL } else { self.n(xp).parent };
                    }
                    if xp != NIL {
                        self.nm(xp).red = false;
                        if xpp != NIL {
                            self.nm(xpp).red = true;
                            root = self.rotate_right(root, xpp);
                        }
                    }
                }
            } else if xppl != NIL && self.n(xppl).red {
                self.nm(xppl).red = false;
                self.nm(xp).red = false;
                self.nm(xpp).red = true;
                x = xpp;
            } else {
                if x == self.n(xp).left {
                    x = xp;
                    root = self.rotate_right(root, x);
                    xp = self.n(x).parent;
                    xpp = if xp == NIL { NIL } else { self.n(xp).parent };
                }
                if xp != NIL {
                    self.nm(xp).red = false;
                    if xpp != NIL {
                        self.nm(xpp).red = true;
                        root = self.rotate_left(root, xpp);
                    }
                }
            }
        }
    }

    fn balance_deletion(&mut self, mut root: usize, mut x: usize) -> usize {
        loop {
            if x == NIL || x == root {
                return root;
            }
            let mut xp = self.n(x).parent;
            if xp == NIL {
                self.nm(x).red = false;
                return x;
            }
            if self.n(x).red {
                self.nm(x).red = false;
                return root;
            }
            let mut xpl = self.n(xp).left;
            if xpl == x {
                let mut xpr = self.n(xp).right;
                if xpr != NIL && self.n(xpr).red {
                    self.nm(xpr).red = false;
                    self.nm(xp).red = true;
                    root = self.rotate_left(root, xp);
                    xp = self.n(x).parent;
                    xpr = if xp == NIL { NIL } else { self.n(xp).right };
                }
                if xpr == NIL {
                    x = xp;
                } else {
                    let sl = self.n(xpr).left;
                    let mut sr = self.n(xpr).right;
                    if !self.red(sr) && !self.red(sl) {
                        self.nm(xpr).red = true;
                        x = xp;
                    } else {
                        if !self.red(sr) {
                            if sl != NIL {
                                self.nm(sl).red = false;
                            }
                            self.nm(xpr).red = true;
                            root = self.rotate_right(root, xpr);
                            xp = self.n(x).parent;
                            xpr = if xp == NIL { NIL } else { self.n(xp).right };
                        }
                        if xpr != NIL {
                            let c = if xp == NIL { false } else { self.n(xp).red };
                            self.nm(xpr).red = c;
                            sr = self.n(xpr).right;
                            if sr != NIL {
                                self.nm(sr).red = false;
                            }
                        }
                        if xp != NIL {
                            self.nm(xp).red = false;
                            root = self.rotate_left(root, xp);
                        }
                        x = root;
                    }
                }
            } else {
                if xpl != NIL && self.n(xpl).red {
                    self.nm(xpl).red = false;
                    self.nm(xp).red = true;
                    root = self.rotate_right(root, xp);
                    xp = self.n(x).parent;
                    xpl = if xp == NIL { NIL } else { self.n(xp).left };
                }
                if xpl == NIL {
                    x = xp;
                } else {
                    let mut sl = self.n(xpl).left;
                    let sr = self.n(xpl).right;
                    if !self.red(sl) && !self.red(sr) {
                        self.nm(xpl).red = true;
                        x = xp;
                    } else {
                        if !self.red(sl) {
                            if sr != NIL {
                                self.nm(sr).red = false;
                            }
                            self.nm(xpl).red = true;
                            root = self.rotate_left(root, xpl);
                            xp = self.n(x).parent;
                            xpl = if xp == NIL { NIL } else { self.n(xp).left };
                        }
                        if xpl != NIL {
                            let c = if xp == NIL { false } else { self.n(xp).red };
                            self.nm(xpl).red = c;
                            sl = self.n(xpl).left;
                            if sl != NIL {
                                self.nm(sl).red = false;
                            }
                        }
                        if xp != NIL {
                            self.nm(xp).red = false;
                            root = self.rotate_right(root, xp);
                        }
                        x = root;
                    }
                }
            }
        }
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
        let p = self.lookup.remove(&key.id_addr())?;
        if self.is_tree {
            self.remove_tree_node(p);
        } else {
            // unlink from the list
            if self.first == p {
                self.first = self.n(p).next;
            } else {
                let mut q = self.first;
                while self.n(q).next != p {
                    q = self.n(q).next;
                }
                let next = self.n(p).next;
                self.nm(q).next = next;
            }
        }
        self.size -= 1;
        let node = self.nodes[p].take().unwrap();
        self.free.push(p);
        Some(node.value)
    }

    fn remove_tree_node(&mut self, p: usize) {
        let mut first = self.first;
        let mut root = first;
        let succ = self.n(p).next;
        let pred = self.n(p).prev;
        if pred == NIL {
            first = succ;
            self.first = succ;
        } else {
            self.nm(pred).next = succ;
        }
        if succ != NIL {
            self.nm(succ).prev = pred;
        }
        if first == NIL {
            self.is_tree = false;
            return;
        }
        if self.n(root).parent != NIL {
            root = self.root_of(root);
        }
        let rl = self.n(root).left;
        if self.n(root).right == NIL || rl == NIL || self.n(rl).left == NIL {
            self.untreeify();
            return;
        }
        let pl = self.n(p).left;
        let pr = self.n(p).right;
        let replacement;
        if pl != NIL && pr != NIL {
            let mut s = pr;
            loop {
                let sl = self.n(s).left;
                if sl == NIL {
                    break;
                }
                s = sl;
            }
            let c = self.n(s).red;
            let pc = self.n(p).red;
            self.nm(s).red = pc;
            self.nm(p).red = c;
            let sr = self.n(s).right;
            let pp = self.n(p).parent;
            if s == pr {
                self.nm(p).parent = s;
                self.nm(s).right = p;
            } else {
                let sp = self.n(s).parent;
                self.nm(p).parent = sp;
                if sp != NIL {
                    if s == self.n(sp).left {
                        self.nm(sp).left = p;
                    } else {
                        self.nm(sp).right = p;
                    }
                }
                self.nm(s).right = pr;
                if pr != NIL {
                    self.nm(pr).parent = s;
                }
            }
            self.nm(p).left = NIL;
            self.nm(p).right = sr;
            if sr != NIL {
                self.nm(sr).parent = p;
            }
            self.nm(s).left = pl;
            if pl != NIL {
                self.nm(pl).parent = s;
            }
            self.nm(s).parent = pp;
            if pp == NIL {
                root = s;
            } else if p == self.n(pp).left {
                self.nm(pp).left = s;
            } else {
                self.nm(pp).right = s;
            }
            replacement = if sr != NIL { sr } else { p };
        } else if pl != NIL {
            replacement = pl;
        } else if pr != NIL {
            replacement = pr;
        } else {
            replacement = p;
        }
        if replacement != p {
            let pp = self.n(p).parent;
            self.nm(replacement).parent = pp;
            if pp == NIL {
                root = replacement;
                self.nm(replacement).red = false;
            } else if p == self.n(pp).left {
                self.nm(pp).left = replacement;
            } else {
                self.nm(pp).right = replacement;
            }
            let node = self.nm(p);
            node.left = NIL;
            node.right = NIL;
            node.parent = NIL;
        }
        let r = if self.n(p).red { root } else { self.balance_deletion(root, replacement) };
        if replacement == p {
            let pp = self.n(p).parent;
            self.nm(p).parent = NIL;
            if pp != NIL {
                if p == self.n(pp).left {
                    self.nm(pp).left = NIL;
                } else if p == self.n(pp).right {
                    self.nm(pp).right = NIL;
                }
            }
        }
        self.move_root_to_front(r);
    }

    /// entries in Java iteration order
    pub fn iter(&self) -> IdIter<'_, K, V> {
        IdIter { map: self, at: self.first }
    }
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.iter().map(|(k, _)| k)
    }
    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.iter().map(|(_, v)| v)
    }
}

pub struct IdIter<'a, K, V> {
    map: &'a IdHashMap<K, V>,
    at: usize,
}

impl<'a, K: IdKey, V> Iterator for IdIter<'a, K, V> {
    type Item = (&'a K, &'a V);
    fn next(&mut self) -> Option<Self::Item> {
        if self.at == NIL {
            return None;
        }
        let node = self.map.nodes[self.at].as_ref().unwrap();
        self.at = node.next;
        Some((&node.key, &node.value))
    }
}

#[derive(Clone)]
pub struct IdHashSet<K> {
    map: IdHashMap<K, ()>,
}

impl<K: IdKey> Default for IdHashSet<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: IdKey> IdHashSet<K> {
    /// new HashSet<>() / Sets.newHashSet()
    pub fn new() -> Self {
        IdHashSet { map: IdHashMap::new() }
    }
    pub fn with_capacity(initial_capacity: usize) -> Self {
        IdHashSet { map: IdHashMap::with_capacity(initial_capacity) }
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
    /// keys in Java iteration order
    pub fn iter(&self) -> impl Iterator<Item = &K> {
        self.map.keys()
    }
    pub fn to_vec(&self) -> Vec<K> {
        self.iter().cloned().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl IdKey for usize {
        fn id_addr(&self) -> usize {
            *self
        }
    }

    /// The iteration orders printed by a JVM run with -XX:hashCode=2 for the same sequences of adds and
    /// removes (HashSet<Object> created with default capacity, capacity 2 and capacity 100 in turn).
    #[test]
    fn matches_java_iteration_order() {
        let expected = include_str!("java_id_hash_expected.txt");
        for (test, line) in expected.lines().enumerate() {
            let mut seed: u64 = 1000 + test as u64 * 7919;
            let ops = 5 + test * 3;
            let mut s: IdHashSet<usize> = match test % 3 {
                0 => IdHashSet::new(),
                1 => IdHashSet::with_capacity(2),
                _ => IdHashSet::with_capacity(100),
            };
            let mut count = 0usize;
            let mut live: Vec<usize> = Vec::new();
            let mut out = String::new();
            for i in 0..ops {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let r = ((seed >> 33) % 100) as usize;
                if r < 70 || live.is_empty() {
                    live.push(count);
                    s.insert(count);
                    count += 1;
                } else {
                    let id = live.remove(r % live.len());
                    s.remove(&id);
                }
                if i % 7 == 6 || i == ops - 1 {
                    for k in s.iter() {
                        out.push_str(&k.to_string());
                        out.push(',');
                    }
                    out.push('|');
                }
            }
            assert_eq!(out, line, "sequence {test}");
        }
    }
}

#[cfg(test)]
mod purple_probe_tests {
    use super::*;

    #[derive(Clone, PartialEq, Eq, Hash)]
    struct K(usize);
    impl IdKey for K {
        fn id_addr(&self) -> usize {
            self.0
        }
    }

    // ArrayListMultimap.create() keyed by an enum (a HashMap of table 16), n keys put in order, under
    // -XX:hashCode=2: the 24 vectors of testdata/constant_hash_order.txt (purple port probe, JDK 23)
    #[test]
    fn multimap_of_sixteen() {
        for line in include_str!("testdata/constant_hash_order.txt").lines() {
            let (n, order) = line.split_once(": ").unwrap();
            let n: usize = n.parse().unwrap();
            let expected: Vec<usize> = order.split(' ').map(|x| x.parse().unwrap()).collect();
            let mut m: IdHashMap<K, ()> = IdHashMap::with_capacity(16);
            for i in 0..n {
                m.insert(K(i), ());
            }
            let got: Vec<usize> = m.keys().map(|k| k.0).collect();
            assert_eq!(got, expected, "n={n}");
        }
    }
}
