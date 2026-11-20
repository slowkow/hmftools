//! `redux-rs compare_bams <a.bam> <b.bam> [threads]`: checks that two BAMs hold the same records in the same order.
//! Integer attributes are compared by value, since writers choose different integer widths.

use crate::rec::*;
use anyhow::Result;
use rust_htslib::bam::{self, Read};

fn aux_equal(a: &[u8], b: &[u8]) -> bool {
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
                k.extend_from_slice(&aux_sorted_copy(aux_bytes(r)));
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
