//! TealConstants and TealUtils: the telomeric content test of every read, reverse complements and the
//! poly-G/C test.
//!
//! `has_telomeric_content_packed` is the prep's inner loop: it tests the BAM's 4-bit packed sequence for
//! TTAGGGTTAGGG or CCCTAACCCTAA without unpacking it. Each 12-base pattern starts at an even or at an odd
//! base; in both cases its packed bytes hold the byte GG (0x44, G pattern) or CC (0x22, C pattern) twice,
//! three bytes apart. Those anchor pairs are found eight byte positions at a time, and each is then checked
//! against the full pattern.

/// TealConstants.CANONICAL_TELOMERE_SEQ
pub const CANONICAL_TELOMERE_SEQ: &str = "TTAGGG";
/// TealConstants.CANONICAL_TELOMERE_SEQ_REV
pub const CANONICAL_TELOMERE_SEQ_REV: &str = "CCCTAA";
/// TealConstants.MIN_TELOMERE_MATCH_BASES
pub const MIN_TELOMERE_MATCH_BASES: usize = 12;
/// TealConstants.POLY_G_THRESHOLD
pub const POLY_G_THRESHOLD: f64 = 0.9;

/// TealConstants.CANONICAL_TELOMERE_SEQUENCES (DEFAULT_MIN_TELE_SEQ_COUNT = 2 repeats)
#[cfg(test)]
const G_PATTERN: &[u8; 12] = b"TTAGGGTTAGGG";
#[cfg(test)]
const C_PATTERN: &[u8; 12] = b"CCCTAACCCTAA";

/// TealUtils.hasTelomericContent on read bases as text
#[cfg(test)]
pub fn has_telomeric_content(bases: &[u8]) -> bool {
    contains(bases, G_PATTERN) || contains(bases, C_PATTERN)
}

#[cfg(test)]
fn contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.len() >= needle.len() && hay.windows(needle.len()).any(|w| w == needle)
}

// 4-bit codes of the BAM sequence alphabet =ACMGRSVTWYHKDBN
const A: u8 = 1;
const C: u8 = 2;
const G: u8 = 4;
const T: u8 = 8;

/// One way a pattern lies on the packed bytes: (byte offset from the first anchor, mask, value) checks
/// and the offset of the pattern's first base from the first anchor's byte, in bases.
struct Placement {
    /// byte offset of the pattern's first byte relative to the first anchor byte
    first_byte: isize,
    /// (byte index from the pattern's first byte, mask, expected value)
    checks: [(usize, u8, u8); 7],
    /// the pattern's first base is the high (0) or low (1) nibble of its first byte
    first_nibble: usize,
}

const fn byte(hi: u8, lo: u8) -> u8 {
    (hi << 4) | lo
}

/// The four placements: G even, G odd, C even, C odd. Checks of the anchor bytes themselves are left out.
const PLACEMENTS_G: [Placement; 2] = [
    // TT AG GG TT AG GG
    Placement {
        first_byte: -2,
        checks: [
            (0, 0xFF, byte(T, T)),
            (1, 0xFF, byte(A, G)),
            (3, 0xFF, byte(T, T)),
            (4, 0xFF, byte(A, G)),
            (0, 0, 0),
            (0, 0, 0),
            (0, 0, 0),
        ],
        first_nibble: 0,
    },
    // xT TA GG GT TA GG Gx
    Placement {
        first_byte: -2,
        checks: [
            (0, 0x0F, T),
            (1, 0xFF, byte(T, A)),
            (3, 0xFF, byte(G, T)),
            (4, 0xFF, byte(T, A)),
            (6, 0xF0, G << 4),
            (0, 0, 0),
            (0, 0, 0),
        ],
        first_nibble: 1,
    },
];

const PLACEMENTS_C: [Placement; 2] = [
    // CC CT AA CC CT AA
    Placement {
        first_byte: 0,
        checks: [
            (1, 0xFF, byte(C, T)),
            (2, 0xFF, byte(A, A)),
            (4, 0xFF, byte(C, T)),
            (5, 0xFF, byte(A, A)),
            (0, 0, 0),
            (0, 0, 0),
            (0, 0, 0),
        ],
        first_nibble: 0,
    },
    // xC CC TA AC CC TA Ax
    Placement {
        first_byte: -1,
        checks: [
            (0, 0x0F, C),
            (2, 0xFF, byte(T, A)),
            (3, 0xFF, byte(A, C)),
            (5, 0xFF, byte(T, A)),
            (6, 0xF0, A << 4),
            (0, 0, 0),
            (0, 0, 0),
        ],
        first_nibble: 1,
    },
];

#[inline]
fn placement_matches(packed: &[u8], seq_len: usize, anchor: usize, p: &Placement) -> bool {
    let first = anchor as isize + p.first_byte;
    if first < 0 {
        return false;
    }
    let first = first as usize;
    // the 12 bases must lie within the read
    if first * 2 + p.first_nibble + 12 > seq_len {
        return false;
    }
    p.checks.iter().all(|&(i, mask, value)| mask == 0 || packed[first + i] & mask == value)
}

/// High bit set in each byte of `x` that is zero
#[inline]
fn zero_bytes(x: u64) -> u64 {
    const LOW7: u64 = 0x7F7F_7F7F_7F7F_7F7F;
    let t = (x & LOW7).wrapping_add(LOW7);
    !(t | x | LOW7)
}

#[inline]
fn load(b: &[u8], i: usize) -> u64 {
    u64::from_le_bytes(b[i..i + 8].try_into().unwrap())
}

/// TealUtils.hasTelomericContent on the 4-bit packed sequence of a BAM record (`seq_len` bases)
pub fn has_telomeric_content_packed(packed: &[u8], seq_len: usize) -> bool {
    if seq_len < 12 {
        return false;
    }
    let n = packed.len().min(seq_len.div_ceil(2));
    const GG: u64 = 0x4444_4444_4444_4444;
    const CC: u64 = 0x2222_2222_2222_2222;
    // anchor positions k with bytes k and k + 3 both GG (or both CC): sixteen at a time with SIMD, then
    // eight at a time
    let mut k = 0usize;
    #[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
    while k + 19 <= n {
        let mut m = simd::anchor_mask16(packed, k);
        while m != 0 {
            let j = k + m.trailing_zeros() as usize;
            m &= m - 1;
            let placements: &[Placement; 2] = if packed[j] == 0x44 { &PLACEMENTS_G } else { &PLACEMENTS_C };
            if placements.iter().any(|p| placement_matches(packed, seq_len, j, p)) {
                return true;
            }
        }
        k += 16;
    }
    while k + 11 <= n {
        let w0 = load(packed, k);
        let w3 = load(packed, k + 3);
        let g = zero_bytes(w0 ^ GG) & zero_bytes(w3 ^ GG);
        let c = zero_bytes(w0 ^ CC) & zero_bytes(w3 ^ CC);
        if (g | c) != 0 && check_anchors(packed, seq_len, k, g, c) {
            return true;
        }
        k += 8;
    }
    // the tail, one position at a time
    while k + 3 < n {
        if packed[k] == 0x44 && packed[k + 3] == 0x44 && PLACEMENTS_G.iter().any(|p| placement_matches(packed, seq_len, k, p)) {
            return true;
        }
        if packed[k] == 0x22 && packed[k + 3] == 0x22 && PLACEMENTS_C.iter().any(|p| placement_matches(packed, seq_len, k, p)) {
            return true;
        }
        k += 1;
    }
    false
}

/// The anchors among the 16 positions from `k` (bytes k and k + 3 both GG or both CC), one bit per
/// position; reads bytes k to k + 18
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
mod simd {
    #[cfg(target_arch = "aarch64")]
    #[inline]
    pub fn anchor_mask16(p: &[u8], k: usize) -> u32 {
        use std::arch::aarch64::*;
        assert!(k + 19 <= p.len());
        // NEON is part of the aarch64 baseline; the loads are within the slice (checked above)
        unsafe {
            let v0 = vld1q_u8(p.as_ptr().add(k));
            let v3 = vld1q_u8(p.as_ptr().add(k + 3));
            let g = vandq_u8(vceqq_u8(v0, vdupq_n_u8(0x44)), vceqq_u8(v3, vdupq_n_u8(0x44)));
            let c = vandq_u8(vceqq_u8(v0, vdupq_n_u8(0x22)), vceqq_u8(v3, vdupq_n_u8(0x22)));
            let m = vorrq_u8(g, c);
            if vmaxvq_u8(m) == 0 {
                return 0;
            }
            // four bits per lane, then one
            let nibbles = vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(m), 4)), 0);
            let mut out = 0u32;
            let mut x = nibbles & 0x1111_1111_1111_1111;
            while x != 0 {
                out |= 1 << (x.trailing_zeros() / 4);
                x &= x - 1;
            }
            out
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    pub fn anchor_mask16(p: &[u8], k: usize) -> u32 {
        use std::arch::x86_64::*;
        assert!(k + 19 <= p.len());
        // SSE2 is part of the x86-64 baseline; the loads are within the slice (checked above)
        unsafe {
            let v0 = _mm_loadu_si128(p.as_ptr().add(k) as *const __m128i);
            let v3 = _mm_loadu_si128(p.as_ptr().add(k + 3) as *const __m128i);
            let gg = _mm_set1_epi8(0x44);
            let cc = _mm_set1_epi8(0x22);
            let g = _mm_and_si128(_mm_cmpeq_epi8(v0, gg), _mm_cmpeq_epi8(v3, gg));
            let c = _mm_and_si128(_mm_cmpeq_epi8(v0, cc), _mm_cmpeq_epi8(v3, cc));
            _mm_movemask_epi8(_mm_or_si128(g, c)) as u32
        }
    }
}

#[cold]
fn check_anchors(packed: &[u8], seq_len: usize, k: usize, mut g: u64, mut c: u64) -> bool {
    while g != 0 {
        let anchor = k + (g.trailing_zeros() / 8) as usize;
        if PLACEMENTS_G.iter().any(|p| placement_matches(packed, seq_len, anchor, p)) {
            return true;
        }
        g &= g - 1;
    }
    while c != 0 {
        let anchor = k + (c.trailing_zeros() / 8) as usize;
        if PLACEMENTS_C.iter().any(|p| placement_matches(packed, seq_len, anchor, p)) {
            return true;
        }
        c &= c - 1;
    }
    false
}

/// TealUtils.reverseComplementSequence: bases other than A, C, G and T are dropped (TE2)
pub fn reverse_complement_java(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .filter_map(|&b| match b {
            b'A' => Some(b'T'),
            b'T' => Some(b'A'),
            b'G' => Some(b'C'),
            b'C' => Some(b'G'),
            _ => None,
        })
        .collect()
}

/// The reverse complement with other bases kept (N for N, IUPAC codes complemented)
pub fn reverse_complement_full(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'M' => b'K',
            b'K' => b'M',
            b'R' => b'Y',
            b'Y' => b'R',
            b'V' => b'B',
            b'B' => b'V',
            b'H' => b'D',
            b'D' => b'H',
            other => other,
        })
        .collect()
}

/// The reverse complement TEAL uses: Java's, which drops non-ACGT bases (TE2), or the full one
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    if crate::compat::java("TE2") { reverse_complement_java(seq) } else { reverse_complement_full(seq) }
}

/// TealUtils.isPolyGC
pub fn is_poly_gc(seq: &[u8]) -> bool {
    let g = seq.iter().filter(|&&b| b == b'G').count();
    let c = seq.iter().filter(|&&b| b == b'C').count();
    let frac = g.max(c) as f64 / seq.len() as f64;
    frac >= POLY_G_THRESHOLD
}

/// StringUtils.countMatches: non-overlapping occurrences
pub fn count_matches(hay: &[u8], needle: &[u8]) -> usize {
    if needle.is_empty() || hay.len() < needle.len() {
        return 0;
    }
    let mut n = 0;
    let mut i = 0;
    while i + needle.len() <= hay.len() {
        if &hay[i..i + needle.len()] == needle {
            n += 1;
            i += needle.len();
        } else {
            i += 1;
        }
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pack(seq: &[u8]) -> Vec<u8> {
        let code = |b: u8| -> u8 { b"=ACMGRSVTWYHKDBN".iter().position(|&c| c == b).unwrap() as u8 };
        seq.chunks(2).map(|p| (code(p[0]) << 4) | if p.len() > 1 { code(p[1]) } else { 0 }).collect()
    }

    /// TealUtilsTest.testTelomericContent
    #[test]
    fn telomeric_content() {
        let t = CANONICAL_TELOMERE_SEQ;
        let read = format!("AGCT{t}AGCT{t}{t}GG");
        assert!(has_telomeric_content(read.as_bytes()));
        let read = format!("AGCT{t}AGCT{t}GG");
        assert!(!has_telomeric_content(read.as_bytes()));
    }

    /// The packed test equals the text test at every offset, both parities, near the ends and with
    /// partial patterns
    #[test]
    fn packed_equals_text() {
        let mut seed = 12345u64;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 33) as usize
        };
        let bases = b"ACGTN";
        let parts: [&[u8]; 6] = [b"TTAGGGTTAGGG", b"CCCTAACCCTAA", b"TTAGGGTTAGG", b"CCCTAACCCTA", b"GGGTTAGGGTTA", b"TAACCCTAACCC"];
        for len in 0..200 {
            for _ in 0..40 {
                let mut s: Vec<u8> = (0..len).map(|_| bases[rnd() % 5]).collect();
                if len >= 12 && rnd() % 2 == 0 {
                    let p = parts[rnd() % parts.len()];
                    let at = rnd() % (len - p.len().min(len) + 1);
                    let end = (at + p.len()).min(len);
                    s[at..end].copy_from_slice(&p[..end - at]);
                }
                if rnd() % 3 == 0 {
                    // runs of GG and CC bytes that only anchor
                    for i in 0..len {
                        if rnd() % 4 == 0 {
                            s[i] = if rnd() % 2 == 0 { b'G' } else { b'C' };
                        }
                    }
                }
                let packed = pack(&s);
                assert_eq!(has_telomeric_content_packed(&packed, s.len()), has_telomeric_content(&s), "{}", String::from_utf8_lossy(&s));
            }
        }
        // every placement of each pattern in a longer read
        for p in [&G_PATTERN[..], &C_PATTERN[..]] {
            for at in 0..140 {
                let mut s = vec![b'A'; 151];
                s[at..at + 12].copy_from_slice(p);
                assert!(has_telomeric_content_packed(&pack(&s), s.len()), "{at}");
                // one base changed: no longer telomeric
                for k in 0..12 {
                    let mut t = s.clone();
                    t[at + k] = b'N';
                    assert!(!has_telomeric_content_packed(&pack(&t), t.len()), "{at} {k}");
                }
            }
        }
    }

    #[test]
    fn reverse_complements() {
        assert_eq!(reverse_complement_java(b"AACGTN"), b"ACGTT".to_vec());
        assert_eq!(reverse_complement_full(b"AACGTN"), b"NACGTT".to_vec());
        assert!(is_poly_gc(b"GGGGGGGGGA"));
        assert!(!is_poly_gc(b"GGGGGGGGAA"));
        assert_eq!(count_matches(b"TTAGGGTTAGGGTTAGG", b"TTAGGG"), 2);
    }
}
