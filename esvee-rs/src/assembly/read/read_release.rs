//! Release of reads that nothing refers to any more (not a Java class: the JVM's garbage collector does this).
//!
//! Reads are counted references and two mates refer to each other, so a pair is never freed by counting
//! alone. When a junction group is finished, the reads loaded for it are checked here: a read is still in use
//! if something other than this list and the mate links within it holds a reference. The rest have their
//! mate links removed, which frees them when the list is dropped. A read in use whose mate is not keeps a
//! copy of the mate's properties in place of the link (after the group's assemblies are built a mate is only
//! reached through the mate accessors of its read). Reads still in use are returned so the check can be
//! repeated later; two mates that are both in use stay linked until then.

use std::collections::HashMap;

use crate::assembly::read::read::ReadP;

/// Takes the only list that owns `reads` on behalf of a junction group and returns those still referenced
/// from elsewhere. Each read must appear once. No other thread may be using these reads through references
/// this function cannot see (a reference held elsewhere only makes it keep more).
pub fn release_unreferenced_reads(reads: Vec<ReadP>) -> Vec<ReadP> {
    if reads.is_empty() {
        return reads;
    }
    let index_by_address: HashMap<usize, usize> = reads.iter().enumerate().map(|(i, r)| (r.addr(), i)).collect();
    // references accounted for: the list itself and mate links from reads in the list
    let mut known_references = vec![1usize; reads.len()];
    let mut mate_index: Vec<Option<usize>> = Vec::with_capacity(reads.len());
    for read in &reads {
        let mate = read.mate_read().and_then(|m| index_by_address.get(&m.addr()).copied());
        if let Some(m) = mate {
            known_references[m] += 1;
        }
        mate_index.push(mate);
    }
    // in use: referenced by something other than this list and the mate links within it
    let in_use: Vec<bool> = reads.iter().enumerate().map(|(i, read)| read.strong_count() > known_references[i]).collect();
    let mut retained: Vec<ReadP> = Vec::new();
    for (i, read) in reads.iter().enumerate() {
        if !in_use[i] {
            read.m().clear_mate_read();
        } else if mate_index[i].is_some_and(|m| !in_use[m]) {
            read.m().replace_mate_read_with_snapshot();
        }
    }
    for (i, read) in reads.into_iter().enumerate() {
        if in_use[i] {
            retained.push(read);
        }
    }
    retained
}
