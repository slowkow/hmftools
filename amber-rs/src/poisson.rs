//! commons-math3's PoissonDistribution: moved to hmf_common::poisson (shared with LINX), re-exported here.

pub use hmf_common::poisson::*;

#[cfg(test)]
mod tests {
    use super::*;

    /// Against commons-math3 itself: the file harness/PoissonRef.java writes (path in AMBER_POISSON_REF).
    /// Run with `cargo test --release -- --ignored poisson_matches_java`.
    #[test]
    #[ignore]
    fn poisson_matches_java() {
        let path = std::env::var("AMBER_POISSON_REF").expect("AMBER_POISSON_REF");
        let text = std::fs::read_to_string(path).unwrap();
        let bits = |s: &str| f64::from_bits(u64::from_str_radix(s, 16).unwrap());
        let (mut n, mut diff_p, mut diff_c, mut diff_k, mut max_ulps) = (0, 0, 0, 0, 0u64);
        for line in text.lines() {
            let f: Vec<&str> = line.split(' ').collect();
            n += 1;
            match f[0] {
                "P" | "C" => {
                    let (mean, k, java) = (bits(f[1]), f[2].parse::<i32>().unwrap(), bits(f[3]));
                    let p = Poisson::new(mean);
                    let rust = if f[0] == "P" { p.probability(k) } else { p.cumulative_probability(k) };
                    if rust.to_bits() != java.to_bits() {
                        if diff_p + diff_c < 12 {
                            eprintln!("{} mean {mean} k {k}: java {java:e} rust {rust:e}", f[0]);
                        }
                        if f[0] == "P" { diff_p += 1 } else { diff_c += 1 }
                        max_ulps = max_ulps.max((rust.to_bits() as i64 - java.to_bits() as i64).unsigned_abs());
                    }
                }
                "K" => {
                    let (total, kmax) = (f[1].parse::<i32>().unwrap(), f[2].parse::<i32>().unwrap());
                    let mut rust = -1;
                    let mut allele = total - 1;
                    while allele >= 0 && allele as f64 > 0.75 * total as f64 {
                        if crate::roh::is_allele_homozygous(total, allele) {
                            rust = total - allele;
                        }
                        allele -= 1;
                    }
                    if rust != kmax {
                        diff_k += 1;
                        eprintln!("ROH threshold differs at total {total}: java {kmax} rust {rust}");
                    }
                }
                _ => {}
            }
        }
        eprintln!("{n} values: probability differs {diff_p}, cumulative differs {diff_c} (max {max_ulps} ulps), ROH thresholds differ {diff_k}");
        assert_eq!((diff_p, diff_c, diff_k), (0, 0, 0));
    }
}
