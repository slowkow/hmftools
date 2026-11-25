//! `redux-rs compare_bams <a.bam> <b.bam> [threads]`: checks that two BAMs hold the same records in the same order.
//! Integer attributes are compared by value, since writers choose different integer widths.
//!
//! `redux-rs diff_bams <a.bam> <b.bam> [threads]`: counts the differences between two coordinate-sorted BAMs of
//! the same reads, by kind. It is how the effect of a correction is measured (REPORT.md, section 3).
//!
//! Both leave out the attributes named in `REDUX_COMPARE_IGNORE_TAGS` (comma-separated, for example `MQ,HN`), so that
//! the alignments of writers that add tags can be compared.

use crate::rec::*;
use anyhow::Result;
use rust_htslib::bam::{self, Read};

/// The attribute tags to leave out of the comparison (REDUX_COMPARE_IGNORE_TAGS)
fn ignored_tags() -> &'static [[u8; 2]] {
    static TAGS: std::sync::OnceLock<Vec<[u8; 2]>> = std::sync::OnceLock::new();
    TAGS.get_or_init(|| {
        std::env::var("REDUX_COMPARE_IGNORE_TAGS")
            .unwrap_or_default()
            .split(',')
            .map(str::trim)
            .filter(|t| t.len() == 2)
            .map(|t| [t.as_bytes()[0], t.as_bytes()[1]])
            .collect()
    })
}

/// The attribute bytes of a record without the ignored tags
fn kept_aux(aux: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    let ignore = ignored_tags();
    if ignore.is_empty() {
        return std::borrow::Cow::Borrowed(aux);
    }
    let spans = aux_spans(aux);
    if !spans.iter().any(|&(s, _)| ignore.contains(&[aux[s], aux[s + 1]])) {
        return std::borrow::Cow::Borrowed(aux);
    }
    let mut out = Vec::with_capacity(aux.len());
    for (s, e) in spans {
        if !ignore.contains(&[aux[s], aux[s + 1]]) {
            out.extend_from_slice(&aux[s..e]);
        }
    }
    std::borrow::Cow::Owned(out)
}

fn aux_equal(a: &[u8], b: &[u8]) -> bool {
    let (a, b) = (kept_aux(a), kept_aux(b));
    let (a, b) = (&a[..], &b[..]);
    if a == b {
        return true;
    }
    let (sa, sb) = (aux_spans(a), aux_spans(b));
    if sa.len() != sb.len() {
        return false;
    }
    let int_value = |f: &[u8]| -> Option<i64> {
        Some(match f[2] {
            b'c' => f[3] as i8 as i64,
            b'C' => f[3] as i64,
            b's' => i16::from_le_bytes([f[3], f[4]]) as i64,
            b'S' => u16::from_le_bytes([f[3], f[4]]) as i64,
            b'i' => i32::from_le_bytes([f[3], f[4], f[5], f[6]]) as i64,
            b'I' => u32::from_le_bytes([f[3], f[4], f[5], f[6]]) as i64,
            _ => return None,
        })
    };
    for (x, y) in sa.iter().zip(sb.iter()) {
        let (fa, fb) = (&a[x.0..x.1], &b[y.0..y.1]);
        if fa == fb {
            continue;
        }
        if fa[..2] != fb[..2] {
            return false;
        }
        match (int_value(fa), int_value(fb)) {
            (Some(va), Some(vb)) if va == vb => {}
            _ => return false,
        }
    }
    true
}

fn describe(r: &bam::Record) -> String {
    format!(
        "{} flags({}) tid({}) pos({}) mapq({}) cigar({}) mate({}:{}) tlen({}) len({})",
        String::from_utf8_lossy(r.qname()), r.flags(), r.tid(), r.pos() + 1, r.mapq(), cigar_string(r.raw_cigar()), r.mtid(), r.mpos() + 1,
        r.insert_size(), r.seq_len()
    )
}

pub fn run(args: &[String]) -> Result<()> {
    let threads: usize = args.get(2).map(|t| t.parse()).transpose()?.unwrap_or(4);
    let mut a = bam::Reader::from_path(&args[0])?;
    let mut b = bam::Reader::from_path(&args[1])?;
    a.set_threads(threads)?;
    b.set_threads(threads)?;
    let names = |h: &bam::HeaderView| -> Vec<Vec<u8>> { (0..h.target_count()).map(|t| h.tid2name(t).to_vec()).collect() };
    if names(a.header()) != names(b.header()) {
        println!("DIFFERENT: reference sequences differ");
        std::process::exit(1);
    }
    let (mut ra, mut rb) = (bam::Record::new(), bam::Record::new());
    let mut n: u64 = 0;
    loop {
        let (xa, xb) = (a.read(&mut ra), b.read(&mut rb));
        match (xa, xb) {
            (None, None) => break,
            (Some(x), Some(y)) => {
                x?;
                y?;
            }
            _ => {
                println!("DIFFERENT: record counts differ after {n} records");
                std::process::exit(1);
            }
        }
        n += 1;
        let same = ra.tid() == rb.tid()
            && ra.pos() == rb.pos()
            && ra.flags() == rb.flags()
            && ra.mapq() == rb.mapq()
            && ra.mtid() == rb.mtid()
            && ra.mpos() == rb.mpos()
            && ra.insert_size() == rb.insert_size()
            && ra.qname() == rb.qname()
            && ra.raw_cigar() == rb.raw_cigar()
            && ra.seq().encoded == rb.seq().encoded
            && ra.seq_len() == rb.seq_len()
            && ra.qual() == rb.qual()
            && aux_equal(aux_bytes(&ra), aux_bytes(&rb));
        if !same && ra.tid() < 0 && rb.tid() < 0 {
            // unplaced reads: Java writes them in an order that depends on thread timing, so compare them as a set
            let key = |r: &bam::Record| -> Vec<u8> {
                let mut k = r.qname().to_vec();
                k.extend_from_slice(&r.flags().to_le_bytes());
                k.extend_from_slice(&r.seq().encoded[..(r.seq_len() + 1) / 2]);
                k.extend_from_slice(r.qual());
                k.extend_from_slice(&aux_sorted_copy(&kept_aux(aux_bytes(r))));
                k
            };
            let (mut ka, mut kb) = (vec![key(&ra)], vec![key(&rb)]);
            while let Some(x) = a.read(&mut ra) {
                x?;
                ka.push(key(&ra));
            }
            while let Some(x) = b.read(&mut rb) {
                x?;
                kb.push(key(&rb));
            }
            ka.sort();
            kb.sort();
            if ka == kb {
                println!("IDENTICAL except the order of the {} unplaced reads at the end ({} records)", ka.len() + 0, n - 1 + ka.len() as u64);
                return Ok(());
            }
            println!("DIFFERENT: the unplaced reads differ ({} and {} after record {})", ka.len(), kb.len(), n - 1);
            std::process::exit(1);
        }
        if !same {
            println!("DIFFERENT at record {n}:\n  a: {}\n  b: {}", describe(&ra), describe(&rb));
            std::process::exit(1);
        }
    }
    println!("IDENTICAL ({n} records)");
    Ok(())
}

/// One side of `diff_bams`: reads the records of one position at a time.
struct Blocks {
    reader: bam::Reader,
    next: Option<bam::Record>,
}

impl Blocks {
    fn open(path: &str, threads: usize) -> Result<Self> {
        let mut reader = bam::Reader::from_path(path)?;
        reader.set_threads(threads)?;
        let mut b = Blocks { reader, next: None };
        b.advance()?;
        Ok(b)
    }
    fn advance(&mut self) -> Result<()> {
        let mut rec = bam::Record::new();
        self.next = match self.reader.read(&mut rec) {
            Some(x) => {
                x?;
                Some(rec)
            }
            None => None,
        };
        Ok(())
    }
    /// (tid as unsigned so that unplaced reads sort last, position)
    fn key(&self) -> Option<(u32, i64)> {
        self.next.as_ref().map(|r| (r.tid() as u32, r.pos()))
    }
    fn take_block(&mut self, key: (u32, i64), out: &mut Vec<bam::Record>) -> Result<()> {
        out.clear();
        while self.key() == Some(key) {
            out.push(self.next.take().unwrap());
            self.advance()?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct DiffCounts {
    records_a: u64,
    records_b: u64,
    identical: u64,
    changed: u64,
    duplicate_flag: u64,
    cigar: u64,
    bases: u64,
    quals: u64,
    attributes: u64,
    other_fields: u64,
    only_a: u64,
    only_a_consensus: u64,
    only_b: u64,
    only_b_consensus: u64,
    unplaced_a: u64,
    unplaced_b: u64,
}

fn record_hash(r: &bam::Record) -> u128 {
    use std::hash::{Hash, Hasher};
    let mut out = 0u128;
    for seed in 0..2u64 {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        seed.hash(&mut h);
        r.qname().hash(&mut h);
        r.flags().hash(&mut h);
        r.seq().encoded[..r.seq_len().div_ceil(2)].hash(&mut h);
        r.qual().hash(&mut h);
        aux_sorted_copy(&kept_aux(aux_bytes(r))).hash(&mut h);
        out = (out << 64) | h.finish() as u128;
    }
    out
}

pub fn diff(args: &[String]) -> Result<()> {
    const DUPLICATE: u16 = 0x400;
    let threads: usize = args.get(2).map(|t| t.parse()).transpose()?.unwrap_or(4);
    let mut a = Blocks::open(&args[0], threads)?;
    let mut b = Blocks::open(&args[1], threads)?;
    let mut c = DiffCounts::default();
    let (mut block_a, mut block_b) = (Vec::new(), Vec::new());
    let (mut unplaced_sum_a, mut unplaced_sum_b) = (0u128, 0u128);
    let is_consensus = |r: &bam::Record| r.qname().starts_with(b"CNS_");
    loop {
        let key = match (a.key(), b.key()) {
            (None, None) => break,
            (Some(x), None) => x,
            (None, Some(y)) => y,
            (Some(x), Some(y)) => x.min(y),
        };
        if key.0 == u32::MAX {
            // unplaced reads have no order to rely on: compared as a set, by a sum of hashes. There can be tens of
            // millions of them at one "position", so they are read one at a time
            while a.key().is_some_and(|k| k.0 == u32::MAX) {
                unplaced_sum_a = unplaced_sum_a.wrapping_add(record_hash(a.next.as_ref().unwrap()));
                c.unplaced_a += 1;
                a.advance()?;
            }
            while b.key().is_some_and(|k| k.0 == u32::MAX) {
                unplaced_sum_b = unplaced_sum_b.wrapping_add(record_hash(b.next.as_ref().unwrap()));
                c.unplaced_b += 1;
                b.advance()?;
            }
            c.records_a += c.unplaced_a;
            c.records_b += c.unplaced_b;
            continue;
        }
        a.take_block(key, &mut block_a)?;
        b.take_block(key, &mut block_b)?;
        c.records_a += block_a.len() as u64;
        c.records_b += block_b.len() as u64;
        let mut used = vec![false; block_b.len()];
        for ra in &block_a {
            let found = block_b.iter().enumerate().position(|(i, rb)| {
                !used[i] && rb.qname() == ra.qname() && (rb.flags() & !DUPLICATE) == (ra.flags() & !DUPLICATE)
            });
            let Some(i) = found else {
                c.only_a += 1;
                c.only_a_consensus += is_consensus(ra) as u64;
                continue;
            };
            used[i] = true;
            let rb = &block_b[i];
            let duplicate_flag = ra.flags() != rb.flags();
            let cigar = ra.raw_cigar() != rb.raw_cigar();
            let bases = ra.seq_len() != rb.seq_len() || ra.seq().encoded != rb.seq().encoded;
            let quals = ra.qual() != rb.qual();
            let attributes = !aux_equal(aux_bytes(ra), aux_bytes(rb));
            let other = ra.mapq() != rb.mapq() || ra.mtid() != rb.mtid() || ra.mpos() != rb.mpos() || ra.insert_size() != rb.insert_size();
            if duplicate_flag || cigar || bases || quals || attributes || other {
                c.changed += 1;
                c.duplicate_flag += duplicate_flag as u64;
                c.cigar += cigar as u64;
                c.bases += bases as u64;
                c.quals += quals as u64;
                c.attributes += attributes as u64;
                c.other_fields += other as u64;
            } else {
                c.identical += 1;
            }
        }
        for (i, rb) in block_b.iter().enumerate() {
            if !used[i] {
                c.only_b += 1;
                c.only_b_consensus += is_consensus(rb) as u64;
            }
        }
    }
    println!("records: a {} b {}", c.records_a, c.records_b);
    println!("placed records identical: {}", c.identical);
    println!(
        "placed records changed: {} (duplicate flag {}, cigar {}, bases {}, qualities {}, attributes {}, other fields {})",
        c.changed, c.duplicate_flag, c.cigar, c.bases, c.quals, c.attributes, c.other_fields
    );
    println!("placed records only in a: {} (consensus reads {})", c.only_a, c.only_a_consensus);
    println!("placed records only in b: {} (consensus reads {})", c.only_b, c.only_b_consensus);
    println!(
        "unplaced records: a {} b {}, {}",
        c.unplaced_a,
        c.unplaced_b,
        if c.unplaced_a == c.unplaced_b && unplaced_sum_a == unplaced_sum_b { "the same set" } else { "different" }
    );
    Ok(())
}
