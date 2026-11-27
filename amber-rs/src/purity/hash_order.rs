//! The iteration order of a java.util.HashSet built by `new HashSet<>()` and a sequence of `add`s, for keys
//! with given hashCodes, including the buckets that HashMap turns into red-black trees (nine or more keys in
//! one bucket of a table of 64 or more), whose list order changes as the tree's root moves to the front.
//!
//! hmf-common-rs has `java_hash` (no tree buckets) and `java_id_hash` (every key in one bucket); this is the
//! general case, a port of HashMap.putVal, resize, treeifyBin and the TreeNode methods it uses (treeify,
//! putTreeVal, split, untreeify, balanceInsertion, rotateLeft/Right, moveRootToFront) of JDK 8 to 23.
//!
//! Keys are told apart by the caller: only keys not yet in the set may be inserted. Keys that are not
//! Comparable and whose hashes are equal are ordered by tieBreakOrder, which compares class names and then
//! System.identityHashCode; this models a JVM with constant identity hashes (`-XX:hashCode=2`), where
//! tieBreakOrder is always -1.

const NIL: u32 = u32::MAX;
const TREEIFY_THRESHOLD: usize = 8;
const UNTREEIFY_THRESHOLD: usize = 6;
const MIN_TREEIFY_CAPACITY: usize = 64;

#[derive(Clone, Copy)]
struct Node {
    /// the spread hash (h ^ h >>> 16), compared as a signed int by the tree code
    hash: i32,
    id: u32,
    next: u32,
    prev: u32,
    parent: u32,
    left: u32,
    right: u32,
    red: bool,
    tree: bool,
}

pub struct JavaHashSetOrder {
    table: Vec<u32>,
    nodes: Vec<Node>,
    size: usize,
    threshold: usize,
}

impl Default for JavaHashSetOrder {
    fn default() -> Self {
        Self::new()
    }
}

impl JavaHashSetOrder {
    /// new HashSet<>()
    pub fn new() -> Self {
        JavaHashSetOrder { table: Vec::new(), nodes: Vec::new(), size: 0, threshold: 0 }
    }

    pub fn len(&self) -> usize {
        self.size
    }

    /// add(key) for a key not in the set, with its hashCode
    pub fn insert(&mut self, id: u32, hash_code: i32) {
        let h = (hash_code ^ ((hash_code as u32) >> 16) as i32) as i32;
        if self.table.is_empty() {
            self.resize();
        }
        let n = self.table.len();
        let i = (h as u32 as usize) & (n - 1);
        let x = self.nodes.len() as u32;
        self.nodes.push(Node { hash: h, id, next: NIL, prev: NIL, parent: NIL, left: NIL, right: NIL, red: false, tree: false });
        let p = self.table[i];
        if p == NIL {
            self.table[i] = x;
        } else if self.nodes[p as usize].tree {
            self.nodes[x as usize].tree = true;
            self.put_tree_val(p, x);
        } else {
            let mut p = p;
            let mut bin_count = 0usize;
            loop {
                let e = self.nodes[p as usize].next;
                if e == NIL {
                    self.nodes[p as usize].next = x;
                    if bin_count >= TREEIFY_THRESHOLD - 1 {
                        self.treeify_bin(h);
                    }
                    break;
                }
                p = e;
                bin_count += 1;
            }
        }
        self.size += 1;
        if self.size > self.threshold {
            self.resize();
        }
    }

    /// the ids in iteration order
    pub fn iteration_order(&self) -> Vec<u32> {
        let mut out = Vec::with_capacity(self.size);
        for &head in &self.table {
            let mut e = head;
            while e != NIL {
                out.push(self.nodes[e as usize].id);
                e = self.nodes[e as usize].next;
            }
        }
        out
    }

    fn resize(&mut self) {
        let old_cap = self.table.len();
        let old_thr = self.threshold;
        let new_cap;
        let mut new_thr = 0usize;
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
            new_thr = (new_cap as f32 * 0.75f32) as usize;
        }
        self.threshold = new_thr;
        let old_tab = std::mem::replace(&mut self.table, vec![NIL; new_cap]);
        for (j, &e) in old_tab.iter().enumerate() {
            if e == NIL {
                continue;
            }
            if self.nodes[e as usize].next == NIL {
                let idx = (self.nodes[e as usize].hash as u32 as usize) & (new_cap - 1);
                self.table[idx] = e;
            } else if self.nodes[e as usize].tree {
                self.split(e, j, old_cap);
            } else {
                let (mut lo_head, mut lo_tail, mut hi_head, mut hi_tail) = (NIL, NIL, NIL, NIL);
                let mut e = e;
                while e != NIL {
                    let next = self.nodes[e as usize].next;
                    if (self.nodes[e as usize].hash as u32 as usize) & old_cap == 0 {
                        if lo_tail == NIL {
                            lo_head = e;
                        } else {
                            self.nodes[lo_tail as usize].next = e;
                        }
                        lo_tail = e;
                    } else {
                        if hi_tail == NIL {
                            hi_head = e;
                        } else {
                            self.nodes[hi_tail as usize].next = e;
                        }
                        hi_tail = e;
                    }
                    e = next;
                }
                if lo_tail != NIL {
                    self.nodes[lo_tail as usize].next = NIL;
                    self.table[j] = lo_head;
                }
                if hi_tail != NIL {
                    self.nodes[hi_tail as usize].next = NIL;
                    self.table[j + old_cap] = hi_head;
                }
            }
        }
    }

    fn treeify_bin(&mut self, hash: i32) {
        let n = self.table.len();
        if n < MIN_TREEIFY_CAPACITY {
            self.resize();
            return;
        }
        let index = (hash as u32 as usize) & (n - 1);
        let head = self.table[index];
        if head == NIL {
            return;
        }
        // replacementTreeNode for each node, linked in the same order
        let mut prev = NIL;
        let mut e = head;
        while e != NIL {
            let node = &mut self.nodes[e as usize];
            node.tree = true;
            node.prev = prev;
            node.parent = NIL;
            node.left = NIL;
            node.right = NIL;
            node.red = false;
            prev = e;
            e = node.next;
        }
        self.treeify(head);
    }

    /// TreeNode.treeify, from the list starting at `head` (already in its bucket)
    fn treeify(&mut self, head: u32) {
        let mut root = NIL;
        let mut x = head;
        while x != NIL {
            let next = self.nodes[x as usize].next;
            self.nodes[x as usize].left = NIL;
            self.nodes[x as usize].right = NIL;
            if root == NIL {
                self.nodes[x as usize].parent = NIL;
                self.nodes[x as usize].red = false;
                root = x;
            } else {
                let h = self.nodes[x as usize].hash;
                let mut p = root;
                loop {
                    let ph = self.nodes[p as usize].hash;
                    // equal hashes: not Comparable, tieBreakOrder gives -1
                    let dir: i32 = if ph > h { -1 } else if ph < h { 1 } else { -1 };
                    let xp = p;
                    p = if dir <= 0 { self.nodes[p as usize].left } else { self.nodes[p as usize].right };
                    if p == NIL {
                        self.nodes[x as usize].parent = xp;
                        if dir <= 0 {
                            self.nodes[xp as usize].left = x;
                        } else {
                            self.nodes[xp as usize].right = x;
                        }
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
        let n = self.table.len();
        let index = (self.nodes[root as usize].hash as u32 as usize) & (n - 1);
        let first = self.table[index];
        if root != first {
            self.table[index] = root;
            let rp = self.nodes[root as usize].prev;
            let rn = self.nodes[root as usize].next;
            if rn != NIL {
                self.nodes[rn as usize].prev = rp;
            }
            if rp != NIL {
                self.nodes[rp as usize].next = rn;
            }
            if first != NIL {
                self.nodes[first as usize].prev = root;
            }
            self.nodes[root as usize].next = first;
            self.nodes[root as usize].prev = NIL;
        }
    }

    fn root_of(&self, mut r: u32) -> u32 {
        loop {
            let p = self.nodes[r as usize].parent;
            if p == NIL {
                return r;
            }
            r = p;
        }
    }

    /// TreeNode.putTreeVal for a new key: `first` is the bucket's first node, `x` the new tree node
    fn put_tree_val(&mut self, first: u32, x: u32) {
        let root = if self.nodes[first as usize].parent != NIL { self.root_of(first) } else { first };
        let h = self.nodes[x as usize].hash;
        let mut p = root;
        loop {
            let ph = self.nodes[p as usize].hash;
            let dir: i32 = if ph > h { -1 } else if ph < h { 1 } else { -1 };
            let xp = p;
            p = if dir <= 0 { self.nodes[p as usize].left } else { self.nodes[p as usize].right };
            if p == NIL {
                let xpn = self.nodes[xp as usize].next;
                self.nodes[x as usize].next = xpn;
                if dir <= 0 {
                    self.nodes[xp as usize].left = x;
                } else {
                    self.nodes[xp as usize].right = x;
                }
                self.nodes[xp as usize].next = x;
                self.nodes[x as usize].parent = xp;
                self.nodes[x as usize].prev = xp;
                if xpn != NIL {
                    self.nodes[xpn as usize].prev = x;
                }
                let r = self.balance_insertion(root, x);
                self.move_root_to_front(r);
                return;
            }
        }
    }

    /// TreeNode.split, during a resize from `bit` (the old capacity) buckets
    fn split(&mut self, b: u32, index: usize, bit: usize) {
        let (mut lo_head, mut lo_tail, mut hi_head, mut hi_tail) = (NIL, NIL, NIL, NIL);
        let (mut lc, mut hc) = (0usize, 0usize);
        let mut e = b;
        while e != NIL {
            let next = self.nodes[e as usize].next;
            self.nodes[e as usize].next = NIL;
            if (self.nodes[e as usize].hash as u32 as usize) & bit == 0 {
                self.nodes[e as usize].prev = lo_tail;
                if lo_tail == NIL {
                    lo_head = e;
                } else {
                    self.nodes[lo_tail as usize].next = e;
                }
                lo_tail = e;
                lc += 1;
            } else {
                self.nodes[e as usize].prev = hi_tail;
                if hi_tail == NIL {
                    hi_head = e;
                } else {
                    self.nodes[hi_tail as usize].next = e;
                }
                hi_tail = e;
                hc += 1;
            }
            e = next;
        }
        if lo_head != NIL {
            if lc <= UNTREEIFY_THRESHOLD {
                self.untreeify(lo_head);
                self.table[index] = lo_head;
            } else {
                self.table[index] = lo_head;
                if hi_head != NIL {
                    self.treeify(lo_head);
                }
            }
        }
        if hi_head != NIL {
            if hc <= UNTREEIFY_THRESHOLD {
                self.untreeify(hi_head);
                self.table[index + bit] = hi_head;
            } else {
                self.table[index + bit] = hi_head;
                if lo_head != NIL {
                    self.treeify(hi_head);
                }
            }
        }
    }

    /// TreeNode.untreeify: plain nodes in the same order
    fn untreeify(&mut self, head: u32) {
        let mut q = head;
        while q != NIL {
            let node = &mut self.nodes[q as usize];
            node.tree = false;
            node.prev = NIL;
            node.parent = NIL;
            node.left = NIL;
            node.right = NIL;
            node.red = false;
            q = node.next;
        }
    }

    fn rotate_left(&mut self, mut root: u32, p: u32) -> u32 {
        if p != NIL {
            let r = self.nodes[p as usize].right;
            if r != NIL {
                let rl = self.nodes[r as usize].left;
                self.nodes[p as usize].right = rl;
                if rl != NIL {
                    self.nodes[rl as usize].parent = p;
                }
                let pp = self.nodes[p as usize].parent;
                self.nodes[r as usize].parent = pp;
                if pp == NIL {
                    root = r;
                    self.nodes[r as usize].red = false;
                } else if self.nodes[pp as usize].left == p {
                    self.nodes[pp as usize].left = r;
                } else {
                    self.nodes[pp as usize].right = r;
                }
                self.nodes[r as usize].left = p;
                self.nodes[p as usize].parent = r;
            }
        }
        root
    }

    fn rotate_right(&mut self, mut root: u32, p: u32) -> u32 {
        if p != NIL {
            let l = self.nodes[p as usize].left;
            if l != NIL {
                let lr = self.nodes[l as usize].right;
                self.nodes[p as usize].left = lr;
                if lr != NIL {
                    self.nodes[lr as usize].parent = p;
                }
                let pp = self.nodes[p as usize].parent;
                self.nodes[l as usize].parent = pp;
                if pp == NIL {
                    root = l;
                    self.nodes[l as usize].red = false;
                } else if self.nodes[pp as usize].right == p {
                    self.nodes[pp as usize].right = l;
                } else {
                    self.nodes[pp as usize].left = l;
                }
                self.nodes[l as usize].right = p;
                self.nodes[p as usize].parent = l;
            }
        }
        root
    }

    fn balance_insertion(&mut self, mut root: u32, mut x: u32) -> u32 {
        self.nodes[x as usize].red = true;
        loop {
            let mut xp = self.nodes[x as usize].parent;
            if xp == NIL {
                self.nodes[x as usize].red = false;
                return x;
            }
            if !self.nodes[xp as usize].red {
                return root;
            }
            let mut xpp = self.nodes[xp as usize].parent;
            if xpp == NIL {
                return root;
            }
            let xppl = self.nodes[xpp as usize].left;
            if xp == xppl {
                let xppr = self.nodes[xpp as usize].right;
                if xppr != NIL && self.nodes[xppr as usize].red {
                    self.nodes[xppr as usize].red = false;
                    self.nodes[xp as usize].red = false;
                    self.nodes[xpp as usize].red = true;
                    x = xpp;
                } else {
                    if x == self.nodes[xp as usize].right {
                        x = xp;
                        root = self.rotate_left(root, x);
                        xp = self.nodes[x as usize].parent;
                        xpp = if xp == NIL { NIL } else { self.nodes[xp as usize].parent };
                    }
                    if xp != NIL {
                        self.nodes[xp as usize].red = false;
                        if xpp != NIL {
                            self.nodes[xpp as usize].red = true;
                            root = self.rotate_right(root, xpp);
                        }
                    }
                }
            } else if xppl != NIL && self.nodes[xppl as usize].red {
                self.nodes[xppl as usize].red = false;
                self.nodes[xp as usize].red = false;
                self.nodes[xpp as usize].red = true;
                x = xpp;
            } else {
                if x == self.nodes[xp as usize].left {
                    x = xp;
                    root = self.rotate_right(root, x);
                    xp = self.nodes[x as usize].parent;
                    xpp = if xp == NIL { NIL } else { self.nodes[xp as usize].parent };
                }
                if xp != NIL {
                    self.nodes[xp as usize].red = false;
                    if xpp != NIL {
                        self.nodes[xpp as usize].red = true;
                        root = self.rotate_left(root, xpp);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order(hashes: &[i32]) -> Vec<u32> {
        let mut s = JavaHashSetOrder::new();
        for (i, &h) in hashes.iter().enumerate() {
            s.insert(i as u32, h);
        }
        s.iteration_order()
    }

    #[test]
    fn small_integers_iterate_in_value_order() {
        // Integer keys hash to themselves: a HashSet of 0..100 added in reverse iterates 0, 1, 2, ...
        let hashes: Vec<i32> = (0..100).rev().collect();
        let got: Vec<i32> = order(&hashes).into_iter().map(|i| hashes[i as usize]).collect();
        assert_eq!(got, (0..100).collect::<Vec<_>>());
    }

    #[test]
    fn colliding_keys_keep_insertion_order_in_a_list() {
        // three keys in bucket 1 of 16, then key 0
        assert_eq!(order(&[1, 17, 33, 0]), vec![3, 0, 1, 2]);
    }

    /// Java (HashOrderCheck in work/amber/tumoronly/harness, -XX:hashCode=2): 20 keys hashing to i * 1024, all in
    /// bucket 0, iterate 7 3 0 1 2 4 5 ...; 30 keys with one hash give 7 3 0 1 2 4 5 6 8 ... 29
    #[test]
    fn tree_bucket() {
        let hashes: Vec<i32> = (0..20).map(|i| i * 1024).collect();
        let mut expected = vec![7, 3, 0, 1, 2, 4, 5, 6];
        expected.extend(8..20);
        assert_eq!(order(&hashes), expected);
        let mut expected = vec![7, 3, 0, 1, 2, 4, 5, 6];
        expected.extend(8..30);
        assert_eq!(order(&[42; 30]), expected);
    }

    /// Every case of the Java check's output file (AMBER_HASH_ORDERS): hashes, then Java's iteration order
    #[test]
    #[ignore]
    fn java_orders() {
        let path = std::env::var("AMBER_HASH_ORDERS").expect("AMBER_HASH_ORDERS");
        let text = std::fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let mut cases = 0;
        for chunk in lines.chunks(3) {
            let name = chunk[0];
            let hashes: Vec<i32> = chunk[1].split(' ').map(|x| x.parse().unwrap()).collect();
            let expected: Vec<u32> = chunk[2].split(' ').map(|x| x.parse().unwrap()).collect();
            assert_eq!(order(&hashes), expected, "{name}");
            cases += 1;
        }
        println!("{cases} cases identical");
    }
}
