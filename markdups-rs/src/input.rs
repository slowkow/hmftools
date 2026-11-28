//! The input BAM files, as AbstractMarkDuplicatesCommandLineProgram.openInputs opens them: one file as
//! it is; several merged by htsjdk's SamFileHeaderMerger (identical sequence dictionaries, read groups and
//! programs merged with colliding ids renamed) and MergingSamRecordIterator (records in coordinate order
//! by SAMRecordCoordinateComparator, taken from a java.util.PriorityQueue of the files).

use anyhow::{Context, Result, bail};
use std::collections::HashMap;

use crate::bam::{self, IdRecord, Rec, SamHeader, SortOrder, read_header};
use crate::bgzf::BgzfReader;

pub trait RecordSource {
    /// The next record's bytes into `buf`; false at the end.
    fn next(&mut self, buf: &mut Vec<u8>) -> Result<bool>;
    /// Whether the last record had a tag changed by the merge (it is then written afresh).
    fn last_touched(&self) -> bool {
        false
    }
}

/// Where the records come from; opened once per pass.
pub trait InputSet {
    /// The header (as MarkDuplicates sees it) and the records, with `threads` threads for inflating.
    fn open(&self, threads: usize) -> Result<(SamHeader, Box<dyn RecordSource + '_>)>;
}

pub struct Inputs {
    pub paths: Vec<String>,
    pub assume_sort_order: Option<SortOrder>,
}

struct Single {
    r: BgzfReader,
}

impl RecordSource for Single {
    fn next(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        self.r.read_record(buf)
    }
}

impl InputSet for Inputs {
    /// Opens the inputs with `threads` inflating threads in all: the header (merged for several files,
    /// with ASSUME_SORT_ORDER applied) and the records.
    fn open(&self, threads: usize) -> Result<(SamHeader, Box<dyn RecordSource + '_>)> {
        if self.paths.is_empty() {
            bail!("no input");
        }
        let per = if self.paths.len() == 1 { threads } else { threads.div_ceil(self.paths.len()).min(threads) };
        let mut readers = Vec::new();
        let mut headers = Vec::new();
        for p in &self.paths {
            let mut r = BgzfReader::open(p, per)?;
            let h = read_header(&mut r)?;
            headers.push(h.sam);
            readers.push(r);
        }
        if let Some(so) = self.assume_sort_order {
            headers[0].set_hd(b"SO", Some(so.name().as_bytes()));
        }
        if readers.len() == 1 {
            let h = headers.pop().unwrap();
            return Ok((h, Box::new(Single { r: readers.pop().unwrap() })));
        }
        let merged = merge_headers(&headers)?;
        let sort_order = merged.header.sort_order();
        if self.assume_sort_order.is_none() && sort_order != SortOrder::Unsorted {
            for h in &headers {
                if h.sort_order() != sort_order {
                    bail!("Files are not compatible with sort order");
                }
            }
        }
        if sort_order != SortOrder::Coordinate {
            bail!("markdups-rs merges several inputs in coordinate order only");
        }
        let source = Merging::new(readers, merged.rg_translation, merged.pg_translation, merged.rg_collisions, merged.pg_collisions)?;
        Ok((merged.header, Box::new(source)))
    }
}

pub struct MergedHeader {
    pub header: SamHeader,
    /// per input file: old id -> new id
    pub rg_translation: Vec<HashMap<Vec<u8>, Vec<u8>>>,
    pub pg_translation: Vec<HashMap<Vec<u8>, Vec<u8>>>,
    pub rg_collisions: bool,
    pub pg_collisions: bool,
}

/// SamFileHeaderMerger.mergeHeaderRecords: records grouped by id, then by equal content; the first
/// of an id keeps it, the others get `id.<base36 counter>`.
fn merge_records(
    items: &[(usize, IdRecord)],
    taken: &mut std::collections::HashSet<Vec<u8>>,
    counter: &mut i32,
    translation: &mut [HashMap<Vec<u8>, Vec<u8>>],
    result: &mut Vec<IdRecord>,
) -> bool {
    // id -> distinct contents (in first-seen order) -> files
    let mut by_id: Vec<(Vec<u8>, Vec<(IdRecord, Vec<usize>)>)> = Vec::new();
    for (file, rec) in items {
        let slot = match by_id.iter().position(|(id, _)| id == &rec.id) {
            Some(i) => i,
            None => {
                by_id.push((rec.id.clone(), Vec::new()));
                by_id.len() - 1
            }
        };
        let same = |a: &IdRecord, b: &IdRecord| {
            let mut x = a.attrs.clone();
            let mut y = b.attrs.clone();
            x.sort();
            y.sort();
            x == y
        };
        let group = &mut by_id[slot].1;
        match group.iter().position(|(r, _)| same(r, rec)) {
            Some(g) => group[g].1.push(*file),
            None => group.push((rec.clone(), vec![*file])),
        }
    }
    let mut collisions = false;
    for (id, group) in by_id {
        for (rec, files) in group {
            let new_id = if !taken.contains(&id) {
                taken.insert(id.clone());
                *counter += 1;
                id.clone()
            } else {
                collisions = true;
                loop {
                    let mut n = id.clone();
                    n.push(b'.');
                    n.extend_from_slice(crate::markdups::base36(*counter).as_bytes());
                    *counter += 1;
                    if taken.insert(n.clone()) {
                        break n;
                    }
                }
            };
            for f in files {
                translation[f].insert(id.clone(), new_id.clone());
            }
            result.push(IdRecord { id: new_id, attrs: rec.attrs.clone() });
        }
    }
    collisions
}

/// SamFileHeaderMerger(sortOrder of the first header, headers, mergeDictionaries = false)
pub fn merge_headers(headers: &[SamHeader]) -> Result<MergedHeader> {
    for h in &headers[1..] {
        let same = h.sq.len() == headers[0].sq.len() && h.sq.iter().zip(&headers[0].sq).all(|(a, b)| a.name == b.name && a.len == b.len);
        if !same {
            bail!("SAM dictionaries are not the same (several inputs need identical sequence dictionaries)");
        }
    }
    let n = headers.len();
    let mut merged = SamHeader { sq: headers[0].sq.clone(), ..SamHeader::default() };
    // programs: roots first, then the records whose PP points at a record just merged (same file)
    let mut pg_translation = vec![HashMap::new(); n];
    let mut pg_taken = std::collections::HashSet::new();
    let mut counter = 0;
    let mut left: Vec<(usize, IdRecord)> = Vec::new();
    for (f, h) in headers.iter().enumerate() {
        let mut ids = std::collections::HashSet::new();
        for p in &h.pg {
            if !ids.insert(p.id.clone()) {
                bail!("Input file contains more than one PG with the same id ({})", String::from_utf8_lossy(&p.id));
            }
            left.push((f, p.clone()));
        }
    }
    let mut current: Vec<(usize, IdRecord)> = Vec::new();
    left.retain(|(f, p)| {
        if p.get(b"PP").is_none() {
            current.push((*f, p.clone()));
            false
        } else {
            true
        }
    });
    let mut programs: Vec<IdRecord> = Vec::new();
    let mut pg_collisions = false;
    while !current.is_empty() {
        let mut res = Vec::new();
        pg_collisions |= merge_records(&current, &mut pg_taken, &mut counter, &mut pg_translation, &mut res);
        programs.extend(res);
        // translateIds: PP of the records left, ID of the current ones
        let translate = |items: &mut Vec<(usize, IdRecord)>, pp_only: bool, tr: &[HashMap<Vec<u8>, Vec<u8>>]| {
            for (f, r) in items.iter_mut() {
                if pp_only {
                    if let Some(pp) = r.get(b"PP").map(|v| v.to_vec()) {
                        if let Some(np) = tr[*f].get(&pp) {
                            let np = np.clone();
                            r.set(b"PP", Some(&np));
                        }
                    }
                } else if let Some(ni) = tr[*f].get(&r.id) {
                    r.id = ni.clone();
                }
            }
        };
        translate(&mut current, false, &pg_translation);
        translate(&mut left, true, &pg_translation);
        let mut next = Vec::new();
        left.retain(|(f, r)| {
            let pp = r.get(b"PP").unwrap();
            if current.iter().any(|(cf, c)| cf == f && c.id == pp) {
                next.push((*f, r.clone()));
                false
            } else {
                true
            }
        });
        current = next;
    }
    if !left.is_empty() {
        bail!("{} program groups weren't processed. Do their PP ids point to existing PGs?", left.len());
    }
    programs.sort_by(|a, b| a.id.cmp(&b.id));
    // read groups
    let mut rg_translation = vec![HashMap::new(); n];
    let mut rg_items = Vec::new();
    for (f, h) in headers.iter().enumerate() {
        let mut ids = std::collections::HashSet::new();
        for r in &h.rg {
            if !ids.insert(r.id.clone()) {
                bail!("Input file contains more than one RG with the same id ({})", String::from_utf8_lossy(&r.id));
            }
            rg_items.push((f, r.clone()));
        }
    }
    let mut rg_taken = std::collections::HashSet::new();
    let mut rg_counter = 0;
    let mut groups = Vec::new();
    let rg_collisions = merge_records(&rg_items, &mut rg_taken, &mut rg_counter, &mut rg_translation, &mut groups);
    groups.sort_by(|a, b| a.id.cmp(&b.id));
    merged.pg = programs;
    merged.rg = groups;
    merged.set_hd(b"GO", Some(b"none"));
    merged.set_hd(b"SO", Some(headers[0].sort_order().name().as_bytes()));
    for h in headers {
        for c in &h.co {
            merged.co.push(c.clone());
        }
    }
    Ok(MergedHeader { header: merged, rg_translation, pg_translation, rg_collisions, pg_collisions })
}

/// SAMRecordCoordinateComparator.compare
pub fn coordinate_compare(a: &Rec, b: &Rec) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    let (ra, rb) = (a.ref_id(), b.ref_id());
    let c = match (ra == -1, rb == -1) {
        (true, true) => Equal,
        (true, false) => Greater,
        (false, true) => Less,
        _ => ra.cmp(&rb).then((a.pos() + 1).cmp(&(b.pos() + 1))),
    };
    if c != Equal {
        return c;
    }
    if a.is_reverse() == b.is_reverse() {
        a.name()
            .cmp(b.name())
            .then(a.flag().cmp(&b.flag()))
            .then(a.mapq().cmp(&b.mapq()))
            .then(a.mate_ref_id().cmp(&b.mate_ref_id()))
            .then((a.mate_pos() + 1).cmp(&(b.mate_pos() + 1)))
            .then(i32::from_le_bytes(a.0[28..32].try_into().unwrap()).cmp(&i32::from_le_bytes(b.0[28..32].try_into().unwrap())))
    } else if a.is_reverse() {
        Greater
    } else {
        Less
    }
}

struct Merging {
    readers: Vec<BgzfReader>,
    peeked: Vec<Vec<u8>>,
    /// java.util.PriorityQueue's array of reader indices
    heap: Vec<usize>,
    rg_tr: Vec<HashMap<Vec<u8>, Vec<u8>>>,
    pg_tr: Vec<HashMap<Vec<u8>, Vec<u8>>>,
    rg_collisions: bool,
    pg_collisions: bool,
    touched: bool,
}

impl Merging {
    fn new(
        readers: Vec<BgzfReader>,
        rg_tr: Vec<HashMap<Vec<u8>, Vec<u8>>>,
        pg_tr: Vec<HashMap<Vec<u8>, Vec<u8>>>,
        rg_collisions: bool,
        pg_collisions: bool,
    ) -> Result<Self> {
        let n = readers.len();
        let mut m = Merging { readers, peeked: vec![Vec::new(); n], heap: Vec::new(), rg_tr, pg_tr, rg_collisions, pg_collisions, touched: false };
        for i in 0..n {
            let mut b = Vec::new();
            if m.readers[i].read_record(&mut b)? {
                Rec::check(&b)?;
                m.peeked[i] = b;
                m.offer(i);
            }
        }
        Ok(m)
    }
    fn less_or_equal(&self, a: usize, b: usize) -> std::cmp::Ordering {
        coordinate_compare(&Rec(&self.peeked[a]), &Rec(&self.peeked[b]))
    }
    /// PriorityQueue.offer: siftUp
    fn offer(&mut self, x: usize) {
        let mut k = self.heap.len();
        self.heap.push(x);
        while k > 0 {
            let parent = (k - 1) >> 1;
            let e = self.heap[parent];
            if self.less_or_equal(x, e) != std::cmp::Ordering::Less {
                break;
            }
            self.heap[k] = e;
            k = parent;
        }
        self.heap[k] = x;
    }
    /// PriorityQueue.poll: siftDown
    fn poll(&mut self) -> Option<usize> {
        if self.heap.is_empty() {
            return None;
        }
        let result = self.heap[0];
        let x = self.heap.pop().unwrap();
        let n = self.heap.len();
        if n > 0 {
            let mut k = 0;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let mut c = self.heap[child];
                let right = child + 1;
                if right < n && self.less_or_equal(c, self.heap[right]) == std::cmp::Ordering::Greater {
                    child = right;
                    c = self.heap[child];
                }
                if self.less_or_equal(x, c) != std::cmp::Ordering::Greater {
                    break;
                }
                self.heap[k] = c;
                k = child;
            }
            self.heap[k] = x;
        }
        Some(result)
    }
}

/// Replaces the value of string tag `tag` in the record's optional fields (in place, the order kept).
fn replace_string_tag(buf: &mut Vec<u8>, tag: [u8; 2], value: &[u8]) -> Result<()> {
    let r = Rec(buf);
    let aux_at = buf.len() - r.aux().len();
    let aux = r.aux();
    let mut at = 0;
    let mut spans = Vec::new();
    while at < aux.len() {
        let (t, v, next) = bam::parse_tag(aux, at)?;
        if t == tag && matches!(v, bam::TagValue::Str(_)) {
            spans.push((at, next));
        }
        at = next;
    }
    for (s, e) in spans.into_iter().rev() {
        let mut new = Vec::with_capacity(value.len() + 4);
        new.extend_from_slice(&tag);
        new.push(b'Z');
        new.extend_from_slice(value);
        new.push(0);
        buf.splice(aux_at + s..aux_at + e, new);
    }
    Ok(())
}

impl RecordSource for Merging {
    fn next(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        let Some(i) = self.poll() else { return Ok(false) };
        std::mem::swap(buf, &mut self.peeked[i]);
        let mut b = std::mem::take(&mut self.peeked[i]);
        if self.readers[i].read_record(&mut b).with_context(|| self.readers[i].path().to_string())? {
            Rec::check(&b)?;
            self.peeked[i] = b;
            self.offer(i);
        }
        self.touched = false;
        if self.rg_collisions {
            if let Some(old) = Rec(buf).string_tag(*b"RG")?.map(|v| v.to_vec()) {
                if let Some(new) = self.rg_tr[i].get(&old) {
                    let new = new.clone();
                    replace_string_tag(buf, *b"RG", &new)?;
                    self.touched = true;
                }
            }
        }
        if self.pg_collisions {
            if let Some(old) = Rec(buf).string_tag(*b"PG")?.map(|v| v.to_vec()) {
                if let Some(new) = self.pg_tr[i].get(&old) {
                    let new = new.clone();
                    replace_string_tag(buf, *b"PG", &new)?;
                    self.touched = true;
                }
            }
        }
        Ok(true)
    }
    fn last_touched(&self) -> bool {
        self.touched
    }
}
