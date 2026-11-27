//! kotlin.random.Random(seed): the XorWow generator of the Kotlin standard library (XorWowRandom), as
//! Kotlin tools use it for reproducible downsampling (CIDER's `Random(0)`).

pub struct KotlinRandom {
    x: i32,
    y: i32,
    z: i32,
    w: i32,
    v: i32,
    addend: i32,
}

impl KotlinRandom {
    /// Random(seed: Int) = XorWowRandom(seed, seed shr 31)
    pub fn new(seed: i32) -> KotlinRandom {
        let seed1 = seed;
        let seed2 = seed >> 31;
        let mut r = KotlinRandom { x: seed1, y: seed2, z: 0, w: 0, v: !seed1, addend: (seed1 << 10) ^ ((seed2 as u32) >> 4) as i32 };
        assert!((r.x | r.y | r.z | r.w | r.v) != 0, "Initial state must have at least one non-zero element.");
        // the constructor discards the first 64 values
        for _ in 0..64 {
            r.next_int();
        }
        r
    }

    /// XorWowRandom.nextInt
    pub fn next_int(&mut self) -> i32 {
        let mut t = self.x;
        t ^= ((t as u32) >> 2) as i32;
        self.x = self.y;
        self.y = self.z;
        self.z = self.w;
        let v0 = self.v;
        self.w = v0;
        t = (t ^ (t << 1)) ^ v0 ^ (v0 << 4);
        self.v = t;
        self.addend = self.addend.wrapping_add(362437);
        t.wrapping_add(self.addend)
    }

    /// nextBits(bitCount): the upper bits of nextInt
    pub fn next_bits(&mut self, bit_count: u32) -> i32 {
        let v = self.next_int();
        (((v as u32) >> (32 - bit_count)) as i32) & ((-(bit_count as i32)) >> 31)
    }

    /// Random.nextDouble: doubleFromParts(nextBits(26), nextBits(27))
    pub fn next_double(&mut self) -> f64 {
        let hi26 = self.next_bits(26);
        let low27 = self.next_bits(27);
        (((hi26 as i64) << 27) + low27 as i64) as f64 / (1i64 << 53) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_zero() {
        // kotlin.random.Random(0).nextDouble() in the Kotlin 2.1 standard library (from a Java harness
        // against kotlin-stdlib; see cider-rs NOTES.md)
        let mut r = KotlinRandom::new(0);
        let v: Vec<f64> = (0..5).map(|_| r.next_double()).collect();
        assert_eq!(v, vec![0.54963315022148, 0.8488554755054833, 0.6591029791044284, 0.0079879341867376, 0.3249690959322964]);
        let mut r = KotlinRandom::new(0);
        assert_eq!((r.next_int(), r.next_int(), r.next_int()), (-1934310868, 1409199696, -649160781));
        let mut r = KotlinRandom::new(-7);
        assert_eq!((r.next_double(), r.next_double()), (0.11492016174825526, 0.716379680973135));
        let mut r = KotlinRandom::new(-7);
        assert_eq!((r.next_int(), r.next_int(), r.next_int()), (493578350, -2123338769, -1218139981));
    }
}
