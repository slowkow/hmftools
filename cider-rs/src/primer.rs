//! -primer_csv: PrimerTsvFile (tab-delimited, with a header) and VdjPrimerMatcher.

use anyhow::{anyhow, bail, Result};

use crate::genes::Vj;
use crate::vdj::VdjSequence;

#[derive(Clone, Debug)]
pub struct Primer {
    pub target: String,
    pub name: String,
    pub sequence: String,
    pub vj: Vj,
}

/// VdjPrimerMatch; `vdj` is the index of the sequence
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct VdjPrimerMatch {
    pub vdj: usize,
    pub primer: Primer,
    pub index: usize,
    pub num_mismatch: i32,
    pub full_vdj_sequence: String,
}

/// PrimerTsvFile.load: columns gene, probeName, sequence and vj (V or J)
pub fn load_primers(path: &str) -> Result<Vec<Primer>> {
    let text = if path.ends_with(".gz") {
        let mut s = String::new();
        use std::io::Read;
        flate2::read::MultiGzDecoder::new(std::fs::File::open(path)?).read_to_string(&mut s)?;
        s
    } else {
        std::fs::read_to_string(path).map_err(|e| anyhow!("reading {path}: {e}"))?
    };
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).ok_or_else(|| anyhow!("Mapping for {n} not found"));
    let (g, p, s, v) = (col("gene")?, col("probeName")?, col("sequence")?, col("vj")?);
    let mut out = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().ok_or_else(|| anyhow!("Index for header '{i}' is out of bounds"));
        let vj = match get(v)? {
            "V" => Vj::V,
            "J" => Vj::J,
            x => bail!("No enum constant com.hartwig.hmftools.cider.genes.VJ.{x}"),
        };
        out.push(Primer { target: get(g)?.to_string(), name: get(p)?.to_string(), sequence: get(s)?.to_string(), vj });
    }
    Ok(out)
}

/// VdjPrimerMatcher.matchVdjPrimer: the first place in each layout sequence where a primer matches with
/// at most `max_mismatch` differences
pub fn match_vdj_primers(vdjs: &[VdjSequence], primers: &[Primer], max_mismatch: i32) -> Vec<VdjPrimerMatch> {
    let mut out = Vec::new();
    for (vi, vdj) in vdjs.iter().enumerate() {
        let seq = vdj.layout.consensus_sequence();
        for primer in primers {
            let p = primer.sequence.as_bytes();
            if p.len() > seq.len() {
                continue;
            }
            for i in 0..=(seq.len() - p.len()) {
                let mut mismatches = 0;
                for (j, &b) in p.iter().enumerate() {
                    if seq[i + j] != b {
                        mismatches += 1;
                        if mismatches > max_mismatch {
                            break;
                        }
                    }
                }
                if mismatches <= max_mismatch {
                    crate::info!("vdj seq matches primer: i: {}, {:?}, mismatch: {}", i, primer.name, mismatches);
                    out.push(VdjPrimerMatch { vdj: vi, primer: primer.clone(), index: i, num_mismatch: mismatches, full_vdj_sequence: String::from_utf8_lossy(seq).into_owned() });
                    break;
                }
            }
        }
    }
    out
}

/// the V and J primer match counts of each sequence
pub fn primer_counts(n: usize, matches: &[VdjPrimerMatch]) -> Vec<(i32, i32)> {
    let mut c = vec![(0, 0); n];
    for m in matches {
        if m.primer.vj == Vj::V {
            c[m.vdj].0 += 1;
        } else {
            c[m.vdj].1 += 1;
        }
    }
    c
}
