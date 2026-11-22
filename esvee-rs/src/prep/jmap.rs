//! java.util.HashMap as far as its iteration order goes. Several results of the Java code depend on the order
//! in which a HashMap or HashSet hands back its entries: by bucket of the current table, and within a bucket
//! in order of insertion. The table doubles when the size passes three quarters of it, or when a ninth entry
//! lands in one bucket of a table smaller than 64, and never shrinks.
//!
//! A bucket that reaches nine entries in a table of 64 or more becomes a red-black tree in Java ("tree bin").
//! Iteration still follows the list links of the bucket, but the tree changes them: its root is moved to the
//! front, and a later entry goes in right after its parent in the tree. So the tree is kept here too, with
//! the code of java.util.HashMap.TreeNode (JDK 23) for insertion, removal, splitting and balancing.

use std::cmp::Ordering as CmpOrdering;
use std::sync::atomic::{AtomicU64, Ordering};

static TREE_BIN_EVENTS: AtomicU64 = AtomicU64::new(0);

/// How many buckets have been turned into trees, for a run to report.
pub fn tree_bin_events() -> u64 {
    TREE_BIN_EVENTS.load(Ordering::Relaxed)
}

/// String.hashCode over single-byte characters
#[inline]
pub fn string_hash(s: &[u8]) -> i32 {
    let mut h: i32 = 0;
    for &c in s {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h
}

/// HashMap.hash: the high half folded into the low half
#[inline]
pub fn spread(h: i32) -> u32 {
    let h = h as u32;
    h ^ (h >> 16)
}

/// A key with Java's hashCode and compareTo (strings by character, numbers by value).
pub trait JavaKey: Eq {
    fn java_hash(&self) -> i32;
    fn java_cmp(&self, other: &Self) -> CmpOrdering;
}

impl JavaKey for Box<[u8]> {
    fn java_hash(&self) -> i32 {
        string_hash(self)
    }
    fn java_cmp(&self, other: &Self) -> CmpOrdering {
        self.cmp(other)
    }
}
impl JavaKey for Vec<u8> {
    fn java_hash(&self) -> i32 {
        string_hash(self)
    }
    fn java_cmp(&self, other: &Self) -> CmpOrdering {
        self.cmp(other)
    }
}
impl JavaKey for String {
    fn java_hash(&self) -> i32 {
        string_hash(self.as_bytes())
    }
    fn java_cmp(&self, other: &Self) -> CmpOrdering {
        self.as_bytes().cmp(other.as_bytes())
    }
}
impl JavaKey for i32 {
    fn java_hash(&self) -> i32 {
        *self
    }
    fn java_cmp(&self, other: &Self) -> CmpOrdering {
        self.cmp(other)
    }
}

const NIL: u32 = u32::MAX;
const TREEIFY_THRESHOLD: usize = 8;
const UNTREEIFY_THRESHOLD: usize = 6;
const MIN_TREEIFY_CAPACITY: usize = 64;

struct Node<K, V> {
    hash: u32,
    next: u32,
    kv: Option<(K, V)>,
    // the fields of a tree node
    tree: bool,
    red: bool,
    parent: u32,
    left: u32,
    right: u32,
    prev: u32,
}

pub struct JavaHashMap<K, V> {
    table: Vec<u32>,
    nodes: Vec<Node<K, V>>,
    free: Vec<u32>,
    size: usize,
    threshold: usize,
}

impl<K: JavaKey, V> Default for JavaHashMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: JavaKey, V> JavaHashMap<K, V> {
    pub fn new() -> Self {
        JavaHashMap { table: Vec::new(), nodes: Vec::new(), free: Vec::new(), size: 0, threshold: 0 }
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.size
    }

    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    #[allow(dead_code)]
    pub fn capacity(&self) -> usize {
        self.table.len()
    }

    #[inline]
    fn n(&self, i: u32) -> &Node<K, V> {
        &self.nodes[i as usize]
    }
    #[inline]
    fn m(&mut self, i: u32) -> &mut Node<K, V> {
        &mut self.nodes[i as usize]
    }
    fn key(&self, i: u32) -> &K {
        &self.n(i).kv.as_ref().unwrap().0
    }

    /// The node of a key: the bucket's list is walked, which finds the same node as Java's tree search.
    fn find(&self, hash: u32, key: &K) -> u32 {
        if self.table.is_empty() {
            return NIL;
        }
        let mut i = self.table[(hash as usize) & (self.table.len() - 1)];
        while i != NIL {
            let n = self.n(i);
            if n.hash == hash && n.kv.as_ref().is_some_and(|kv| kv.0 == *key) {
                return i;
            }
            i = n.next;
        }
        NIL
    }

    pub fn get(&self, key: &K) -> Option<&V> {
        let i = self.find(spread(key.java_hash()), key);
        if i == NIL { None } else { self.n(i).kv.as_ref().map(|kv| &kv.1) }
    }

    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        let i = self.find(spread(key.java_hash()), key);
        if i == NIL { None } else { self.m(i).kv.as_mut().map(|kv| &mut kv.1) }
    }

    pub fn contains_key(&self, key: &K) -> bool {
        self.find(spread(key.java_hash()), key) != NIL
    }

    fn new_node(&mut self, hash: u32, key: K, value: V) -> u32 {
        let node = Node { hash, next: NIL, kv: Some((key, value)), tree: false, red: false, parent: NIL, left: NIL, right: NIL, prev: NIL };
        match self.free.pop() {
            Some(f) => {
                self.nodes[f as usize] = node;
                f
            }
            None => {
                self.nodes.push(node);
                (self.nodes.len() - 1) as u32
            }
        }
    }

    fn resize(&mut self) {
        let old_cap = self.table.len();
        let new_cap = if old_cap == 0 { 16 } else { old_cap * 2 };
        self.threshold = new_cap / 4 * 3;
        let old_table = std::mem::replace(&mut self.table, vec![NIL; new_cap]);
        // each bucket splits into the same index and index + old capacity, both in the old order
        for (b, &head) in old_table.iter().enumerate() {
            if head == NIL {
                continue;
            }
            let is_tree = self.n(head).tree;
            let (mut lo_head, mut lo_tail, mut hi_head, mut hi_tail) = (NIL, NIL, NIL, NIL);
            let (mut lc, mut hc) = (0usize, 0usize);
            let mut i = head;
            while i != NIL {
                let next = self.n(i).next;
                self.m(i).next = NIL;
                if (self.n(i).hash as usize) & old_cap == 0 {
                    self.m(i).prev = lo_tail;
                    if lo_tail == NIL { lo_head = i } else { self.m(lo_tail).next = i }
                    lo_tail = i;
                    lc += 1;
                } else {
                    self.m(i).prev = hi_tail;
                    if hi_tail == NIL { hi_head = i } else { self.m(hi_tail).next = i }
                    hi_tail = i;
                    hc += 1;
                }
                i = next;
            }
            if !is_tree {
                self.table[b] = lo_head;
                self.table[b + old_cap] = hi_head;
                continue;
            }
            // TreeNode.split: a part of six or fewer goes back to a list; otherwise it is a tree, built
            // again unless it is the whole of the old one
            if lo_head != NIL {
                self.table[b] = lo_head;
                if lc <= UNTREEIFY_THRESHOLD {
                    self.untreeify(lo_head);
                } else if hi_head != NIL {
                    self.treeify(lo_head);
                }
            }
            if hi_head != NIL {
                self.table[b + old_cap] = hi_head;
                if hc <= UNTREEIFY_THRESHOLD {
                    self.untreeify(hi_head);
                } else if lo_head != NIL {
                    self.treeify(hi_head);
                }
            }
        }
    }

    fn untreeify(&mut self, head: u32) {
        let mut i = head;
        while i != NIL {
            let n = self.m(i);
            n.tree = false;
            n.red = false;
            n.parent = NIL;
            n.left = NIL;
            n.right = NIL;
            n.prev = NIL;
            i = n.next;
        }
    }

    /// The side of `p` on which a node of this hash and key goes (TreeNode's ordering: hash as a signed
    /// number, then the key)
    fn tree_dir(&self, h: u32, k: &K, p: u32) -> i32 {
        let ph = self.n(p).hash as i32;
        let h = h as i32;
        if ph > h {
            -1
        } else if ph < h {
            1
        } else {
            match k.java_cmp(self.key(p)) {
                CmpOrdering::Less => -1,
                CmpOrdering::Greater => 1,
                // equal keys are one entry; Java would fall back to identity hashes here
                CmpOrdering::Equal => -1,
            }
        }
    }

    /// TreeNode.treeify: the nodes of the list from `head`, in list order, into a tree
    fn treeify(&mut self, head: u32) {
        let mut root = NIL;
        let mut x = head;
        while x != NIL {
            let next = self.n(x).next;
            {
                let n = self.m(x);
                n.tree = true;
                n.left = NIL;
                n.right = NIL;
            }
            if root == NIL {
                self.m(x).parent = NIL;
                self.m(x).red = false;
                root = x;
            } else {
                let h = self.n(x).hash;
                let mut p = root;
                loop {
                    let dir = self.tree_dir(h, self.key(x), p);
                    let xp = p;
                    p = if dir <= 0 { self.n(p).left } else { self.n(p).right };
                    if p == NIL {
                        self.m(x).parent = xp;
                        if dir <= 0 { self.m(xp).left = x } else { self.m(xp).right = x }
                        root = self.balance_insertion(root, x);
                        break;
                    }
                }
            }
            x = next;
        }
        self.move_root_to_front(root);
    }

    fn move_root_to_front(&mut self, root: u32) {
        if root == NIL || self.table.is_empty() {
            return;
        }
        let index = (self.n(root).hash as usize) & (self.table.len() - 1);
        let first = self.table[index];
        if root != first {
            self.table[index] = root;
            let rp = self.n(root).prev;
            let rn = self.n(root).next;
            if rn != NIL {
                self.m(rn).prev = rp;
            }
            if rp != NIL {
                self.m(rp).next = rn;
            }
            if first != NIL {
                self.m(first).prev = root;
            }
            self.m(root).next = first;
            self.m(root).prev = NIL;
        }
    }

    fn rotate_left(&mut self, mut root: u32, p: u32) -> u32 {
        if p == NIL {
            return root;
        }
        let r = self.n(p).right;
        if r != NIL {
            let rl = self.n(r).left;
            self.m(p).right = rl;
            if rl != NIL {
                self.m(rl).parent = p;
            }
            let pp = self.n(p).parent;
            self.m(r).parent = pp;
            if pp == NIL {
                root = r;
                self.m(r).red = false;
            } else if self.n(pp).left == p {
                self.m(pp).left = r;
            } else {
                self.m(pp).right = r;
            }
            self.m(r).left = p;
            self.m(p).parent = r;
        }
        root
    }

    fn rotate_right(&mut self, mut root: u32, p: u32) -> u32 {
        if p == NIL {
            return root;
        }
        let l = self.n(p).left;
        if l != NIL {
            let lr = self.n(l).right;
            self.m(p).left = lr;
            if lr != NIL {
                self.m(lr).parent = p;
            }
            let pp = self.n(p).parent;
            self.m(l).parent = pp;
            if pp == NIL {
                root = l;
                self.m(l).red = false;
            } else if self.n(pp).right == p {
                self.m(pp).right = l;
            } else {
                self.m(pp).left = l;
            }
            self.m(l).right = p;
            self.m(p).parent = l;
        }
        root
    }

    fn is_red(&self, i: u32) -> bool {
        i != NIL && self.n(i).red
    }

    fn balance_insertion(&mut self, mut root: u32, mut x: u32) -> u32 {
        self.m(x).red = true;
        loop {
            let mut xp = self.n(x).parent;
            if xp == NIL {
                self.m(x).red = false;
                return x;
            }
            if !self.n(xp).red {
                return root;
            }
            let mut xpp = self.n(xp).parent;
            if xpp == NIL {
                return root;
            }
            let xppl = self.n(xpp).left;
            if xp == xppl {
                let xppr = self.n(xpp).right;
                if self.is_red(xppr) {
                    self.m(xppr).red = false;
                    self.m(xp).red = false;
                    self.m(xpp).red = true;
                    x = xpp;
                } else {
                    if x == self.n(xp).right {
                        x = xp;
                        root = self.rotate_left(root, x);
                        xp = self.n(x).parent;
                        xpp = if xp == NIL { NIL } else { self.n(xp).parent };
                    }
                    if xp != NIL {
                        self.m(xp).red = false;
                        if xpp != NIL {
                            self.m(xpp).red = true;
                            root = self.rotate_right(root, xpp);
                        }
                    }
                }
            } else if self.is_red(xppl) {
                self.m(xppl).red = false;
                self.m(xp).red = false;
                self.m(xpp).red = true;
                x = xpp;
            } else {
                if x == self.n(xp).left {
                    x = xp;
                    root = self.rotate_right(root, x);
                    xp = self.n(x).parent;
                    xpp = if xp == NIL { NIL } else { self.n(xp).parent };
                }
                if xp != NIL {
                    self.m(xp).red = false;
                    if xpp != NIL {
                        self.m(xpp).red = true;
                        root = self.rotate_left(root, xpp);
                    }
                }
            }
        }
    }

    fn balance_deletion(&mut self, mut root: u32, mut x: u32) -> u32 {
        loop {
            if x == NIL || x == root {
                return root;
            }
            let mut xp = self.n(x).parent;
            if xp == NIL {
                self.m(x).red = false;
                return x;
            }
            if self.n(x).red {
                self.m(x).red = false;
                return root;
            }
            let mut xpl = self.n(xp).left;
            if xpl == x {
                let mut xpr = self.n(xp).right;
                if self.is_red(xpr) {
                    self.m(xpr).red = false;
                    self.m(xp).red = true;
                    root = self.rotate_left(root, xp);
                    xp = self.n(x).parent;
                    xpr = if xp == NIL { NIL } else { self.n(xp).right };
                }
                if xpr == NIL {
                    x = xp;
                } else {
                    let (sl, mut sr) = (self.n(xpr).left, self.n(xpr).right);
                    if !self.is_red(sr) && !self.is_red(sl) {
                        self.m(xpr).red = true;
                        x = xp;
                    } else {
                        if !self.is_red(sr) {
                            if sl != NIL {
                                self.m(sl).red = false;
                            }
                            self.m(xpr).red = true;
                            root = self.rotate_right(root, xpr);
                            xp = self.n(x).parent;
                            xpr = if xp == NIL { NIL } else { self.n(xp).right };
                        }
                        if xpr != NIL {
                            let red = if xp == NIL { false } else { self.n(xp).red };
                            self.m(xpr).red = red;
                            sr = self.n(xpr).right;
                            if sr != NIL {
                                self.m(sr).red = false;
                            }
                        }
                        if xp != NIL {
                            self.m(xp).red = false;
                            root = self.rotate_left(root, xp);
                        }
                        x = root;
                    }
                }
            } else {
                if self.is_red(xpl) {
                    self.m(xpl).red = false;
                    self.m(xp).red = true;
                    root = self.rotate_right(root, xp);
                    xp = self.n(x).parent;
                    xpl = if xp == NIL { NIL } else { self.n(xp).left };
                }
                if xpl == NIL {
                    x = xp;
                } else {
                    let (mut sl, sr) = (self.n(xpl).left, self.n(xpl).right);
                    if !self.is_red(sl) && !self.is_red(sr) {
                        self.m(xpl).red = true;
                        x = xp;
                    } else {
                        if !self.is_red(sl) {
                            if sr != NIL {
                                self.m(sr).red = false;
                            }
                            self.m(xpl).red = true;
                            root = self.rotate_left(root, xpl);
                            xp = self.n(x).parent;
                            xpl = if xp == NIL { NIL } else { self.n(xp).left };
                        }
                        if xpl != NIL {
                            let red = if xp == NIL { false } else { self.n(xp).red };
                            self.m(xpl).red = red;
                            sl = self.n(xpl).left;
                            if sl != NIL {
                                self.m(sl).red = false;
                            }
                        }
                        if xp != NIL {
                            self.m(xp).red = false;
                            root = self.rotate_right(root, xp);
                        }
                        x = root;
                    }
                }
            }
        }
    }

    fn tree_root(&self, mut i: u32) -> u32 {
        while self.n(i).parent != NIL {
            i = self.n(i).parent;
        }
        i
    }

    /// TreeNode.putTreeVal for a key that is not in the bucket: the new node goes below its place in the
    /// tree and, in the list, right after its tree parent
    fn put_tree_val(&mut self, first: u32, x: u32) {
        let root = self.tree_root(first);
        let h = self.n(x).hash;
        let mut p = root;
        loop {
            let dir = self.tree_dir(h, self.key(x), p);
            let xp = p;
            p = if dir <= 0 { self.n(p).left } else { self.n(p).right };
            if p == NIL {
                let xpn = self.n(xp).next;
                {
                    let n = self.m(x);
                    n.tree = true;
                    n.next = xpn;
                    n.parent = xp;
                    n.prev = xp;
                }
                if dir <= 0 { self.m(xp).left = x } else { self.m(xp).right = x }
                self.m(xp).next = x;
                if xpn != NIL {
                    self.m(xpn).prev = x;
                }
                let r = self.balance_insertion(root, x);
                self.move_root_to_front(r);
                return;
            }
        }
    }

    /// TreeNode.removeTreeNode
    fn remove_tree_node(&mut self, p: u32) {
        let index = (self.n(p).hash as usize) & (self.table.len() - 1);
        let mut first = self.table[index];
        let mut root = first;
        let (succ, pred) = (self.n(p).next, self.n(p).prev);
        if pred == NIL {
            first = succ;
            self.table[index] = first;
        } else {
            self.m(pred).next = succ;
        }
        if succ != NIL {
            self.m(succ).prev = pred;
        }
        if first == NIL {
            return;
        }
        if self.n(root).parent != NIL {
            root = self.tree_root(root);
        }
        let rl = self.n(root).left;
        if self.n(root).right == NIL || rl == NIL || self.n(rl).left == NIL {
            // too small
            self.untreeify(first);
            return;
        }
        let (pl, pr) = (self.n(p).left, self.n(p).right);
        let replacement;
        if pl != NIL && pr != NIL {
            let mut s = pr;
            while self.n(s).left != NIL {
                s = self.n(s).left;
            }
            let c = self.n(s).red;
            let p_red = self.n(p).red;
            self.m(s).red = p_red;
            self.m(p).red = c;
            let sr = self.n(s).right;
            let pp = self.n(p).parent;
            if s == pr {
                self.m(p).parent = s;
                self.m(s).right = p;
            } else {
                let sp = self.n(s).parent;
                self.m(p).parent = sp;
                if sp != NIL {
                    if s == self.n(sp).left { self.m(sp).left = p } else { self.m(sp).right = p }
                }
                self.m(s).right = pr;
                if pr != NIL {
                    self.m(pr).parent = s;
                }
            }
            self.m(p).left = NIL;
            self.m(p).right = sr;
            if sr != NIL {
                self.m(sr).parent = p;
            }
            self.m(s).left = pl;
            if pl != NIL {
                self.m(pl).parent = s;
            }
            self.m(s).parent = pp;
            if pp == NIL {
                root = s;
            } else if p == self.n(pp).left {
                self.m(pp).left = s;
            } else {
                self.m(pp).right = s;
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
            self.m(replacement).parent = pp;
            if pp == NIL {
                root = replacement;
                self.m(root).red = false;
            } else if p == self.n(pp).left {
                self.m(pp).left = replacement;
            } else {
                self.m(pp).right = replacement;
            }
            let n = self.m(p);
            n.left = NIL;
            n.right = NIL;
            n.parent = NIL;
        }
        let r = if self.n(p).red { root } else { self.balance_deletion(root, replacement) };
        if replacement == p {
            let pp = self.n(p).parent;
            self.m(p).parent = NIL;
            if pp != NIL {
                if p == self.n(pp).left {
                    self.m(pp).left = NIL;
                } else if p == self.n(pp).right {
                    self.m(pp).right = NIL;
                }
            }
        }
        self.move_root_to_front(r);
    }

    /// HashMap.put: returns the previous value of the key
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        let hash = spread(key.java_hash());
        if self.table.is_empty() {
            self.resize();
        }
        let b = (hash as usize) & (self.table.len() - 1);
        let first = self.table[b];
        let mut i = first;
        let mut tail = NIL;
        let mut count = 0usize;
        while i != NIL {
            let n = self.m(i);
            if n.hash == hash && n.kv.as_ref().is_some_and(|kv| kv.0 == key) {
                return n.kv.as_mut().map(|kv| std::mem::replace(&mut kv.1, value));
            }
            tail = i;
            count += 1;
            i = n.next;
        }
        let idx = self.new_node(hash, key, value);
        if first == NIL {
            self.table[b] = idx;
        } else if self.n(first).tree {
            self.put_tree_val(first, idx);
        } else {
            self.m(tail).next = idx;
            // the bucket held TREEIFY_THRESHOLD entries before this one
            if count >= TREEIFY_THRESHOLD {
                if self.table.len() < MIN_TREEIFY_CAPACITY {
                    self.resize();
                } else {
                    TREE_BIN_EVENTS.fetch_add(1, Ordering::Relaxed);
                    // treeifyBin: the list as tree nodes with links both ways, then the tree
                    let mut prev = NIL;
                    let mut e = first;
                    while e != NIL {
                        self.m(e).prev = prev;
                        prev = e;
                        e = self.n(e).next;
                    }
                    self.treeify(first);
                }
            }
        }
        self.size += 1;
        if self.size > self.threshold {
            self.resize();
        }
        None
    }

    pub fn remove(&mut self, key: &K) -> Option<V> {
        let hash = spread(key.java_hash());
        let i = self.find(hash, key);
        if i == NIL {
            return None;
        }
        if self.n(i).tree {
            self.remove_tree_node(i);
        } else {
            let b = (hash as usize) & (self.table.len() - 1);
            let next = self.n(i).next;
            if self.table[b] == i {
                self.table[b] = next;
            } else {
                let mut prev = self.table[b];
                while self.n(prev).next != i {
                    prev = self.n(prev).next;
                }
                self.m(prev).next = next;
            }
        }
        self.size -= 1;
        self.free.push(i);
        self.m(i).kv.take().map(|kv| kv.1)
    }

    /// Entries in Java's iteration order.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        let mut bucket = 0usize;
        let mut i = NIL;
        std::iter::from_fn(move || {
            while i == NIL {
                if bucket >= self.table.len() {
                    return None;
                }
                i = self.table[bucket];
                bucket += 1;
            }
            let n = &self.nodes[i as usize];
            i = n.next;
            n.kv.as_ref().map(|kv| (&kv.0, &kv.1))
        })
    }

    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.iter().map(|(_, v)| v)
    }

    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.iter().map(|(k, _)| k)
    }
}

/// The iteration order of a small Java HashSet of strings that was only added to: indices into `hashes`
/// (String.hashCode values in insertion order).
pub fn set_order(hashes: &[i32]) -> Vec<usize> {
    let mut cap = 16usize;
    while hashes.len() > cap / 4 * 3 {
        cap *= 2;
    }
    let mut order: Vec<usize> = (0..hashes.len()).collect();
    order.sort_by_key(|&i| (spread(hashes[i]) as usize) & (cap - 1));
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn string_hash_matches_java() {
        // "abc".hashCode() and "".hashCode()
        assert_eq!(string_hash(b"abc"), 96354);
        assert_eq!(string_hash(b""), 0);
        // overflow wraps as in Java: "polygenelubricants".hashCode() == Integer.MIN_VALUE
        assert_eq!(string_hash(b"polygenelubricants"), i32::MIN);
    }

    #[test]
    fn integer_keys_iterate_by_bucket() {
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        for k in [33, 1, 17, 16, 0] {
            m.insert(k, k);
        }
        // table of 16: 0 and 16 share bucket 0 in insertion order, then 33, 1, 17 in bucket 1
        assert_eq!(m.keys().copied().collect::<Vec<_>>(), vec![16, 0, 33, 1, 17]);
        // the thirteenth entry doubles the table
        for k in 100..108 {
            m.insert(k, k);
        }
        assert_eq!(m.capacity(), 32);
        assert_eq!(m.keys().copied().take(5).collect::<Vec<_>>(), vec![0, 33, 1, 100, 101]);
        assert_eq!(m.remove(&33), Some(33));
        assert_eq!(m.len(), 12);
        assert!(!m.contains_key(&33));
        assert_eq!(m.insert(1, 5), Some(1));
        assert_eq!(m.get(&1), Some(&5));
    }

    #[test]
    fn nine_in_a_bucket_doubles_a_small_table() {
        let mut m: JavaHashMap<i32, ()> = JavaHashMap::new();
        for i in 0..9 {
            m.insert(i * 16, ());
        }
        assert_eq!(m.capacity(), 32);
    }

    #[test]
    fn small_set_order() {
        // "b" and "a" hash to buckets 2 and 1
        assert_eq!(set_order(&[string_hash(b"b"), string_hash(b"a")]), vec![1, 0]);
    }
    /// The linear congruential generator the reference program used.
    struct Lcg(i64);
    impl Lcg {
        fn next(&mut self) -> i32 {
            self.0 = (self.0.wrapping_mul(1103515245) + 12345) & 0x7fffffff;
            self.0 as i32
        }
    }

    /// Checksums of the iteration order of a java.util.HashMap<Integer,Integer> (JDK 23) after every 500 of
    /// 20,000 puts and removes of keys that crowd into a few buckets, so that buckets become trees, split on
    /// resize and go back to lists.
    const INTEGER_KEY_CHECKSUMS: [i64; 40] = [
        243012370994204, 69370349116429, 153875301833872, 256878024963560, 105626168782740, 73342134752525,
        157468819133778, 45363903774505, 225752065923396, 234121959507150, 263407090899431, 133203252072599,
        158008153761027, 788385921564, 151411109539332, 252420699103714, 104091903178484, 186455328595273,
        67635007109939, 273020938409349, 247973578185703, 269053504917994, 2773331519011, 248703960412902,
        184759855613094, 175843959203361, 133588904755546, 267118788076248, 209435177461846, 38901991864871,
        223980305159666, 57971148122191, 228111159938320, 160525948801036, 230018310780196, 94189390842375,
        24053190989875, 34489053590996, 232267684200888, 196594541232725,
    ];

    /// The same for a HashMap<String,Integer> whose keys mostly share one hashCode.
    const STRING_KEY_CHECKSUMS: [i64; 40] = [
        223102292654811, 270996962569239, 212774050501134, 128233685379578, 266742463240181, 170110393849753,
        122462460705342, 218239447752372, 178707148285747, 173583710617096, 30622001581344, 47446985311028,
        11572087511563, 31122933921871, 160191845881498, 2888267686409, 153591574913728, 272195493806981,
        262478453581711, 250816331627120, 75221205918355, 84376824601450, 162359366092523, 868667966645,
        143575666774454, 209888169132783, 83920736065530, 171768475136313, 79841272824793, 201940244106115,
        109178953654932, 247356043608114, 97569186010241, 12133559393439, 210164345395314, 230501752385372,
        89173893819414, 31285274982068, 253916466686791, 189097099744666,
    ];

    #[test]
    fn tree_bins_iterate_as_java_with_integer_keys() {
        let mut rng = Lcg(12345);
        let mut m: JavaHashMap<i32, i32> = JavaHashMap::new();
        let mut sums = Vec::new();
        for step in 1..=20000 {
            let r = rng.next();
            let key = (((r >> 8) % 600) - 300) * 256 + ((r >> 3) % 2);
            if r % 3 == 2 {
                m.remove(&key);
            } else {
                m.insert(key, step);
            }
            if step % 500 == 0 {
                let mut sum: i64 = 0;
                for (k, v) in m.iter() {
                    sum = (sum * 31 + *k as i64 * 7 + *v as i64) & 0xffffffffffff;
                }
                sums.push(sum);
            }
        }
        assert_eq!(sums, INTEGER_KEY_CHECKSUMS);
        assert!(tree_bin_events() > 0);
    }

    #[test]
    fn tree_bins_iterate_as_java_with_string_keys() {
        let mut rng = Lcg(777);
        let mut m: JavaHashMap<String, i32> = JavaHashMap::new();
        let mut sums = Vec::new();
        for step in 1..=20000 {
            let r = rng.next();
            let q = r >> 4;
            // ten of "Aa" or "BB", which have the same hashCode, then one of three endings
            let mut key = String::new();
            for i in 0..10 {
                key.push_str(if (q % 1024 >> i) & 1 == 0 { "Aa" } else { "BB" });
            }
            match (q >> 12) % 3 {
                1 => key.push('x'),
                2 => key.push('y'),
                _ => {}
            }
            if r % 3 == 2 {
                m.remove(&key);
            } else {
                m.insert(key, step);
            }
            if step % 500 == 0 {
                let mut sum: i64 = 0;
                for (k, v) in m.iter() {
                    sum = (sum * 31 + k.len() as i64 * 7 + *v as i64) & 0xffffffffffff;
                }
                sums.push(sum);
            }
        }
        assert_eq!(sums, STRING_KEY_CHECKSUMS);
    }
}
