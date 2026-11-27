//! TargetRegionStats (fragment counts of each target region of a partition) and OffTargetAnalyser (the
//! pile-ups of off-target fragments).

use std::collections::HashMap;

use hmf_common::regions::positions_overlap;

use crate::compat;
use crate::stats::ReadCounts;

/// TargetRegionStats
#[derive(Clone, Debug)]
pub struct TargetRegionStats {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub fragment_counts: ReadCounts,
}

impl TargetRegionStats {
    /// a row of the target regions file
    pub fn row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chromosome, self.start, self.end, self.fragment_counts.total, self.fragment_counts.duplicates, self.fragment_counts.dual_strand
        )
    }
}

pub const TARGET_REGIONS_HEADER: &str = "Chromosome\tPositionStart\tPositionEnd\tTotalFragments\tDuplicateFragments\tDualStrandFragments";
pub const ENRICHED_HEADER: &str = "Chromosome\tPosStart\tPosEnd\tPeakPosStart\tPeakPosEnd\tFragments\tMaxOverlap";

/// OffTargetAnalyser. Java keeps the fragments of the current pile-up in a list it scans for each new
/// read's name and compares pairwise, and counts the depth of every base of the peak fragment over all
/// fragments: quadratic, which never finishes on continuous coverage (a whole partition is one
/// pile-up once BT5 feeds it). Here the names are in a set, the overlap counts come from the sorted
/// starts and ends, and the peak's depths from a difference array; the results are Java's (the
/// fragments arrive in start order, as the BAM is sorted; for any other order, or a fragment ending before
/// it starts, Java's loops are used).
pub struct OffTargetAnalyser {
    chromosome: String,
    /// (start, end) of the current pile-up's fragments, in arrival order
    fragments: Vec<(i32, i32)>,
    names: std::collections::HashSet<Vec<u8>>,
    enriched: Vec<String>,
    current_lower: i32,
    current_upper: i32,
    /// overlap count -> fragments with it
    pub overlap_counts: HashMap<i32, i64>,
    write_threshold: i32,
}

impl OffTargetAnalyser {
    pub fn new(chromosome: &str, write_threshold: i32) -> Self {
        OffTargetAnalyser {
            chromosome: chromosome.to_string(),
            fragments: Vec::new(),
            names: std::collections::HashSet::new(),
            enriched: Vec::new(),
            current_lower: 0,
            current_upper: 0,
            overlap_counts: HashMap::new(),
            write_threshold,
        }
    }

    /// OffTargetAnalyser.addRead with the fragment's coordinates (Fragment's constructor)
    pub fn add_read(&mut self, read_id: &[u8], unmapped: bool, start: i32, end: i32) {
        if unmapped {
            return;
        }
        if self.names.contains(read_id) {
            return;
        }
        if start > self.current_upper {
            self.process_overlapping_fragments();
        }
        if self.fragments.is_empty() {
            self.current_lower = start;
        }
        self.current_upper = self.current_upper.max(end);
        self.fragments.push((start, end));
        self.names.insert(read_id.to_vec());
    }

    /// each fragment's overlap count (itself included), in fragment order
    fn overlap_counts_of(&self) -> Vec<i32> {
        let f = &self.fragments;
        // the identity below needs the starts in order and no fragment ending before it starts (a
        // zero-length alignment); otherwise Java's loop
        let sorted = f.windows(2).all(|w| w[0].0 <= w[1].0) && f.iter().all(|x| x.1 >= x.0);
        if !sorted {
            // Java's loop, with its early exit at the first later fragment starting past this one's end
            return (0..f.len())
                .map(|i| {
                    let mut n = 1;
                    for j in 0..f.len() {
                        if i == j {
                            continue;
                        }
                        if positions_overlap(f[i].0, f[i].1, f[j].0, f[j].1) {
                            n += 1;
                        } else if f[j].0 > f[i].1 {
                            break;
                        }
                    }
                    n
                })
                .collect();
        }
        let starts: Vec<i32> = f.iter().map(|x| x.0).collect();
        let mut ends: Vec<i32> = f.iter().map(|x| x.1).collect();
        ends.sort_unstable();
        f.iter()
            .map(|&(s, e)| {
                // fragments starting at or before this one's end, less those ending before its start
                let started = starts.partition_point(|&x| x <= e);
                let ended = ends.partition_point(|&x| x < s);
                (started - ended) as i32
            })
            .collect()
    }

    /// processOverlappingFragments: each fragment's overlap count (itself included), and an enriched site
    /// where the deepest pile-up reaches the threshold
    fn process_overlapping_fragments(&mut self) {
        if self.fragments.is_empty() {
            return;
        }
        let counts = self.overlap_counts_of();
        let mut max_overlaps = 0;
        let mut peak = 0usize;
        for (i, &overlaps) in counts.iter().enumerate() {
            *self.overlap_counts.entry(overlaps).or_insert(0) += 1;
            if overlaps > max_overlaps {
                max_overlaps = overlaps;
                peak = i;
            }
        }
        // the depth along the peak fragment from a difference array
        let (ps, pe) = self.fragments[peak];
        let (mut max_depth, mut max_start, mut max_end) = (0, 0, 0);
        if pe >= ps {
            let len = (pe - ps + 1) as usize;
            let mut diff = vec![0i32; len + 1];
            for &(s, e) in &self.fragments {
                let (a, b) = (s.max(ps), e.min(pe));
                if a <= b {
                    diff[(a - ps) as usize] += 1;
                    diff[(b - ps) as usize + 1] -= 1;
                }
            }
            let java_bt7 = compat::java("BT7");
            let mut depth = 0;
            for (k, d) in diff.iter().take(len).enumerate() {
                depth += d;
                let pos = ps + k as i32;
                if depth > max_depth {
                    max_depth = depth;
                    max_start = pos;
                    // BT7: Java leaves the peak's end where an earlier, lower plateau ended (or 0), so the
                    // peak can end before it starts. Corrected: a new maximum starts the peak region.
                    if !java_bt7 {
                        max_end = pos;
                    }
                } else if depth == max_depth {
                    max_end = pos;
                }
            }
        }
        if self.write_threshold > 0 && max_depth >= self.write_threshold {
            self.enriched.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.chromosome,
                self.current_lower,
                self.current_upper,
                max_start,
                max_end,
                self.fragments.len(),
                max_depth
            ));
        }
        self.fragments.clear();
        self.names.clear();
    }

    /// at the end of the partition. BT6: Java never processes the fragments still held, so each
    /// partition's last pile-up is lost. Corrected: they are processed.
    pub fn finish(&mut self) {
        if !compat::java("BT6") {
            self.process_overlapping_fragments();
        }
    }

    /// enrichedFragmentSites: the rows of the enriched-region file, cleared
    pub fn take_enriched(&mut self) -> Vec<String> {
        std::mem::take(&mut self.enriched)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pile_ups() {
        for java in [true, false] {
            compat::with_java(java, || {
                let mut a = OffTargetAnalyser::new("chr1", 2);
                a.add_read(b"r1", false, 100, 200);
                a.add_read(b"r1", false, 100, 200); // the same fragment again
                a.add_read(b"r2", false, 150, 300);
                a.add_read(b"r3", false, 250, 260);
                a.add_read(b"r4", false, 1000, 1100); // past the pile-up: it is processed
                a.finish();
                let rows = a.take_enriched();
                // r1-r2 overlap, r2-r3 overlap: counts 2, 3, 2; the peak fragment is r2 (150-300),
                // depth 2 at 150-200 and 250-260
                assert_eq!(a.overlap_counts.get(&3), Some(&1));
                if java {
                    assert_eq!(a.overlap_counts.get(&2), Some(&2));
                    assert_eq!(a.overlap_counts.get(&1), None); // r4 never processed (BT6)
                    assert_eq!(rows, vec!["chr1\t100\t300\t150\t260\t3\t2".to_string()]);
                } else {
                    assert_eq!(a.overlap_counts.get(&1), Some(&1));
                    assert_eq!(rows, vec!["chr1\t100\t300\t150\t260\t3\t2".to_string()]);
                }
                // a higher depth after a plateau: Java's end stays at the plateau's end (BT7)
                let mut b = OffTargetAnalyser::new("chr1", 1);
                b.add_read(b"a", false, 10, 30);
                b.add_read(b"b", false, 20, 30);
                b.add_read(b"c", false, 25, 25);
                b.add_read(b"d", false, 100, 110);
                b.finish();
                let rows = b.take_enriched();
                if java {
                    assert_eq!(rows, vec!["chr1\t10\t30\t25\t24\t3\t3".to_string()]);
                } else {
                    assert_eq!(rows, vec!["chr1\t10\t30\t25\t25\t3\t3".to_string(), "chr1\t100\t110\t100\t110\t1\t1".to_string()]);
                }
            });
        }
    }
}

#[cfg(test)]
mod fast_counts {
    use super::*;

    /// the fast overlap counts against Java's loop on random sorted fragments (some with end < start)
    #[test]
    fn counts_equal_java_loop() {
        let mut seed: u64 = 7;
        let mut next = |m: u64| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 33) % m) as i32
        };
        for _ in 0..500 {
            let n = 1 + next(80) as usize;
            let mut a = OffTargetAnalyser::new("chr1", 1);
            let mut pos = 100;
            for _ in 0..n {
                pos += next(40);
                let len = next(120) - 3;
                a.fragments.push((pos, pos + len));
            }
            let fast = a.overlap_counts_of();
            let f = a.fragments.clone();
            for (i, c) in fast.iter().enumerate() {
                let mut java = 1;
                for j in 0..f.len() {
                    if i == j {
                        continue;
                    }
                    if positions_overlap(f[i].0, f[i].1, f[j].0, f[j].1) {
                        java += 1;
                    } else if f[j].0 > f[i].1 {
                        break;
                    }
                }
                assert_eq!(*c, java, "fragment {i} of {f:?}");
            }
        }
    }
}

#[cfg(test)]
mod java_harness {
    //! Replays the off-target cases of harness/BtHarness.java (testdata/java_harness.tsv)
    use super::*;

    #[test]
    fn pile_ups_match_java() {
        let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/java_harness.tsv")).unwrap();
        let mut cases = 0;
        compat::with_java(true, || {
            for line in text.lines().filter(|l| l.starts_with("OFFTARGET\t")) {
                let f: Vec<&str> = line.split('\t').collect();
                let mut a = OffTargetAnalyser::new("chr1", f[1].parse().unwrap());
                for frag in f[2].split(';') {
                    let p: Vec<&str> = frag.split(',').collect();
                    a.add_read(p[0].as_bytes(), false, p[1].parse().unwrap(), p[2].parse().unwrap());
                }
                a.finish();
                let mut keys: Vec<&i32> = a.overlap_counts.keys().collect();
                keys.sort();
                let counts: Vec<String> = keys.iter().map(|k| format!("{k}={}", a.overlap_counts[*k])).collect();
                assert_eq!(counts.join(","), f.get(3).copied().unwrap_or(""), "{line}");
                // the sites are tab-separated rows joined by '|'
                let expected = if f.len() > 4 { f[4..].join("\t") } else { String::new() };
                assert_eq!(a.take_enriched().join("|"), expected, "{line}");
                cases += 1;
            }
        });
        assert_eq!(cases, 300);
    }
}
