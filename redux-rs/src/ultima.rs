//! Ultima reads (Java consensus/UltimaRoutines, consensus/UltimaStats and what they use of hmf-common
//! sequencing/UltimaBamUtils): the attributes of each read are reduced to the ones the pipeline needs, the
//! consensus type (UT) and the low-quality base summary (UQ) are added, and the low-quality base counts are
//! written to `<sample>.redux.ultima_stats.tsv`.

use crate::rec::*;
use anyhow::{bail, Result};
use rust_htslib::bam;
use std::io::Write;

pub const ULTIMA_MAX_QUAL: i8 = 35;
pub const HALF_PHRED_SCORE_SCALING: i8 = 3;
const PHRED_OFFSET: u8 = 33;
const ULTIMA_HP_DELETION_LOW_QUAL_THRESHOLD: i8 = 20;

const PPM_STRAND_MIN_SUM: f64 = 4.0;
const PPM_STRAND_MAX_SUM: f64 = 8.0;
const PPM_STRAND_BALANCED_LOW: f64 = 0.27;
const PPM_STRAND_BALANCED_HIGH: f64 = 1.0 - PPM_STRAND_BALANCED_LOW;

/// UltimaRoutines.NON_ULTIMA_REQUIRED_ATTRIBUTES
const REQUIRED_ATTRIBUTES: [&[u8; 2]; 5] = [b"AS", b"RG", b"NM", b"XS", b"SA"];
/// UltimaRoutines.ULTIMA_RAW_QUAL_ATTRIBUTES
const RAW_QUAL_ATTRIBUTES: [&[u8; 2]; 2] = [b"t0", b"tp"];

/// SAMTag.makeBinaryTag: the second character is the high byte.
#[inline]
fn java_tag(tag: &[u8]) -> u16 {
    (tag[1] as u16) << 8 | tag[0] as u16
}

/// UltimaBamUtils.isHighBaseQual
#[inline]
fn is_high_base_qual(qual: i8) -> bool {
    qual >= ULTIMA_MAX_QUAL
}

/// UltimaRoutines.belowLowQualThreshold
#[inline]
fn below_low_qual_threshold(homopolymer_length: usize, qual: i8) -> bool {
    match homopolymer_length {
        1 => qual < 25,
        2 => qual < 22,
        3 => qual < 19,
        4..=6 => qual < 16,
        7 => qual < 14,
        _ => qual < 12,
    }
}

/// SAMRecord.getBaseQualities: none when the record has no qualities (first byte 0xff)
#[inline]
fn base_qualities(rec: &bam::Record) -> &[u8] {
    let q = rec.qual();
    if q.first() == Some(&0xff) { &[] } else { q }
}

/// UltimaRoutines.formLowQualTag without the leading count: appends the index entries ("2,4-6,8") to `entries`
/// and returns the number of low-quality bases; no entry means Java returns null (no tag). `base(i)` gives the
/// read base at an index (only compared for equality), `t0` is the t0 attribute still with its Phred offset.
/// Java fails with an index error when t0 is shorter than the base qualities; here that is an `Err`.
fn form_low_qual_entries(base: impl Fn(usize) -> u8, base_quals: &[u8], t0: &[u8], entries: &mut Vec<u8>) -> Result<i32> {
    let n = base_quals.len();
    if t0.len() < n {
        bail!("Ultima read with a t0 attribute of {} values for {} base qualities (Java: ArrayIndexOutOfBoundsException in UltimaRoutines.formLowQualTag)", t0.len(), n);
    }
    let mut hp_start: i64 = -1;
    let mut hp_end: i64 = -1;
    let mut hp_low_qual = false;
    let mut low_qual_start: i64 = -1;
    let mut total: i32 = 0;

    let add_entry = |entries: &mut Vec<u8>, s: i64, e: i64| {
        if !entries.is_empty() {
            entries.push(b',');
        }
        if e > s {
            let _ = write!(entries, "{s}-{e}");
        } else {
            let _ = write!(entries, "{s}");
        }
    };

    for i in 0..n {
        let qual = base_quals[i] as i8;
        // deletions: the t0 value is a Java byte, so the subtraction wraps
        let mut is_low_qual = (t0[i].wrapping_sub(PHRED_OFFSET) as i8) < ULTIMA_HP_DELETION_LOW_QUAL_THRESHOLD;

        if hp_start >= 0 && i as i64 <= hp_end {
            // inside a homopolymer that has been assessed
            is_low_qual |= hp_low_qual;
        } else {
            // a new homopolymer is assessed only from a base below the maximum quality; hp_start and hp_end keep
            // the values of the previous one otherwise
            hp_low_qual = false;
            if !is_high_base_qual(qual) {
                let b = base(i);
                let mut hp_length = 1;
                let mut j = i + 1;
                while j < n && base(j) == b {
                    hp_length += 1;
                    j += 1;
                }
                let mut tp_qual = qual;
                if hp_length > 1 {
                    // int arithmetic on a byte: a negative quality becomes 0
                    tp_qual = if (tp_qual as i32) - (HALF_PHRED_SCORE_SCALING as i32) < 0 { 0 } else { tp_qual - HALF_PHRED_SCORE_SCALING };
                }
                hp_start = i as i64;
                hp_end = (i + hp_length - 1) as i64;
                if below_low_qual_threshold(hp_length, tp_qual) {
                    is_low_qual = true;
                    hp_low_qual = true;
                }
            }
        }

        if is_low_qual {
            if low_qual_start < 0 {
                low_qual_start = i as i64;
            }
            if i == n - 1 {
                add_entry(entries, low_qual_start, n as i64 - 1);
                total += (n as i64 - low_qual_start) as i32;
            }
        } else if low_qual_start >= 0 {
            let low_qual_end = i as i64 - 1;
            add_entry(entries, low_qual_start, low_qual_end);
            total += (low_qual_end - low_qual_start + 1) as i32;
            low_qual_start = -1;
        }
    }
    Ok(total)
}

/// UltimaRoutines.formLowQualTag as a string ("5=2,4-6,8"), for tests and callers off the hot path.
#[allow(dead_code)]
pub fn form_low_qual_tag(bases: &[u8], base_quals: &[u8], t0: &[u8]) -> Result<Option<String>> {
    let mut entries = Vec::new();
    let total = form_low_qual_entries(|i| bases[i], base_quals, t0, &mut entries)?;
    Ok(if entries.is_empty() { None } else { Some(format!("{}={}", total, String::from_utf8(entries).unwrap())) })
}

/// The type byte and payload of an attribute, when present (the last one for a repeated tag, as htsjdk keeps it).
fn find_field<'a>(aux: &'a [u8], tag: &[u8; 2]) -> Option<&'a [u8]> {
    let mut found = None;
    let mut i = 0;
    while i + 3 <= aux.len() {
        let n = field_len(&aux[i..]);
        if &aux[i..i + 2] == tag {
            found = Some(&aux[i + 2..i + n]);
        }
        i += n;
    }
    found
}

/// Length in bytes of the attribute that starts at `aux[0]` (tag, type, value).
fn field_len(aux: &[u8]) -> usize {
    let fixed = |t: u8| match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        _ => 0,
    };
    match aux[2] {
        b'Z' | b'H' => 3 + aux[3..].iter().position(|&c| c == 0).unwrap() + 1,
        b'B' => 8 + u32::from_le_bytes([aux[4], aux[5], aux[6], aux[7]]) as usize * fixed(aux[3]),
        t => 3 + fixed(t),
    }
}

/// The integer value of an attribute (type byte and payload), when it is an integer.
fn field_int(f: &[u8]) -> Option<i64> {
    Some(match f[0] {
        b'c' => f[1] as i8 as i64,
        b'C' => f[1] as i64,
        b's' => i16::from_le_bytes([f[1], f[2]]) as i64,
        b'S' => u16::from_le_bytes([f[1], f[2]]) as i64,
        b'i' => i32::from_le_bytes([f[1], f[2], f[3], f[4]]) as i64,
        b'I' => u32::from_le_bytes([f[1], f[2], f[3], f[4]]) as i64,
        _ => return None,
    })
}

/// SAMRecord.getStringAttribute: `Ok(None)` when absent; a value of another type is a ClassCastException in Java.
fn string_attribute<'a>(aux: &'a [u8], tag: &[u8; 2]) -> Result<Option<&'a [u8]>> {
    match find_field(aux, tag) {
        None => Ok(None),
        Some(f) if f[0] == b'Z' => Ok(Some(&f[1..f.len() - 1])),
        Some(f) => bail!(
            "Ultima read with attribute {} of type {} where a string is required (Java: ClassCastException in SAMRecord.getStringAttribute)",
            String::from_utf8_lossy(tag), f[0] as char
        ),
    }
}

/// UltimaBamUtils.getIntegerAttribute: the default when absent; a value that is not an integer is a
/// ClassCastException in Java.
fn integer_attribute(aux: &[u8], tag: &[u8; 2], default: i64) -> Result<i64> {
    match find_field(aux, tag) {
        None => Ok(default),
        Some(f) => match field_int(f) {
            Some(v) if v <= i32::MAX as i64 => Ok(v),
            _ => bail!(
                "Ultima read with attribute {} that is not an integer (Java: ClassCastException in UltimaBamUtils.getIntegerAttribute)",
                String::from_utf8_lossy(tag)
            ),
        },
    }
}

/// UltimaBamUtils.deriveConsensusType: true for ConsensusType.DUAL, false for NONE. With the ppmSeq st / et
/// attributes both must be MIXED; without either, the older as / ts / ae / te lengths decide.
fn is_dual_consensus_type(aux: &[u8]) -> Result<bool> {
    let st = string_attribute(aux, b"st")?;
    let et = string_attribute(aux, b"et")?;
    if st.is_none() && et.is_none() {
        return is_dual_consensus_type_old(aux);
    }
    Ok(st == Some(b"MIXED") && et == Some(b"MIXED"))
}

/// UltimaBamUtils.deriveConsensusTypeOld
fn is_dual_consensus_type_old(aux: &[u8]) -> Result<bool> {
    let a_start = integer_attribute(aux, b"as", 0)?;
    let t_start = integer_attribute(aux, b"ts", 0)?;
    let a_end = integer_attribute(aux, b"ae", 0)?;
    let t_end = integer_attribute(aux, b"te", 0)?;
    // the sums are int sums in Java (they wrap), then doubles
    let sum_start = (a_start as i32).wrapping_add(t_start as i32) as f64;
    let sum_end = (a_end as i32).wrapping_add(t_end as i32) as f64;
    if sum_start < PPM_STRAND_MIN_SUM || sum_start > PPM_STRAND_MAX_SUM || sum_end < PPM_STRAND_MIN_SUM || sum_end > PPM_STRAND_MAX_SUM {
        return Ok(false);
    }
    let ratio_start = a_start as f64 / sum_start;
    let ratio_end = a_end as f64 / sum_end;
    if ratio_start < PPM_STRAND_BALANCED_LOW || ratio_start > PPM_STRAND_BALANCED_HIGH {
        return Ok(false);
    }
    if ratio_end < PPM_STRAND_BALANCED_LOW || ratio_end > PPM_STRAND_BALANCED_HIGH {
        return Ok(false);
    }
    Ok(true)
}

/// Appends an attribute as htsjdk writes it after SAMRecord.setAttribute: an integer in the smallest type that
/// holds it (BinaryTagCodec.getIntegerType), an unsigned array as a signed one (setAttribute(String, Object) does
/// not carry the unsigned marker), anything else unchanged. `field` is the type byte and payload.
fn append_decoded(out: &mut Vec<u8>, tag: &[u8], field: &[u8]) -> Result<()> {
    out.extend_from_slice(&tag[..2]);
    if let Some(val) = field_int(field) {
        if val > i32::MAX as i64 {
            // htsjdk holds this as a Long and writes an unsigned 32-bit value
            out.push(b'I');
            out.extend_from_slice(&(val as u32).to_le_bytes());
        } else if val >= 0 {
            if val <= 127 {
                out.extend_from_slice(&[b'c', val as u8]);
            } else if val <= 255 {
                out.extend_from_slice(&[b'C', val as u8]);
            } else if val <= 32767 {
                out.push(b's');
                out.extend_from_slice(&(val as i16).to_le_bytes());
            } else if val <= 65535 {
                out.push(b'S');
                out.extend_from_slice(&(val as u16).to_le_bytes());
            } else {
                out.push(b'i');
                out.extend_from_slice(&(val as i32).to_le_bytes());
            }
        } else if val >= -128 {
            out.extend_from_slice(&[b'c', val as i8 as u8]);
        } else if val >= -32768 {
            out.push(b's');
            out.extend_from_slice(&(val as i16).to_le_bytes());
        } else {
            out.push(b'i');
            out.extend_from_slice(&(val as i32).to_le_bytes());
        }
        return Ok(());
    }
    if field[0] == b'B' {
        if u32::from_le_bytes([field[2], field[3], field[4], field[5]]) == 0 {
            bail!(
                "Ultima read with an empty array in attribute {} (Java: IllegalArgumentException in SAMRecord.setAttribute)",
                String::from_utf8_lossy(&tag[..2])
            );
        }
        out.push(b'B');
        out.push(match field[1] {
            b'C' => b'c',
            b'S' => b's',
            b'I' => b'i',
            t => t,
        });
        out.extend_from_slice(&field[2..]);
        return Ok(());
    }
    out.extend_from_slice(field);
    Ok(())
}

/// The failure of REDUX 2.1 on an Ultima duplicate group whose reads disagree at a base, with -form_consensus:
/// BaseBuilder.determineBaseAndQual returns BaseQualPair.INVALID for Ultima, the consensus read gets a base of
/// value 0, and htsjdk refuses to encode it when the sorted writer writes the read (BAMRecordCodec.encode:
/// "IllegalStateException: Bad base passed to charToCompressedBaseLow", or "...High"). There is no Java output to
/// match and a BAM record cannot hold that base, so the port stops when such a read is formed. With -no_write_bam
/// Java never encodes the read and completes, having passed the read with its 0 base to the base quality
/// recalibration; a record of the port cannot tell that base from N, and the recalibration counts then differ from
/// Java's, so the port stops in that case as well.
pub fn undetermined_consensus_base_error(read_id: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "Ultima consensus read({read_id}) has a base without a value: REDUX 2.1 has no consensus rule for Ultima reads that \
         disagree at a base (BaseBuilder.determineBaseAndQual) and fails when it writes such a read (htsjdk: Bad base passed to \
         charToCompressedBaseLow or ...High); run Ultima reads without -form_consensus"
    )
}

/// UltimaRoutines, with the buffers it reuses from read to read. One per partition reader.
pub struct UltimaRoutines {
    /// -ultima_drop_qual_tags (UTLIMA_DROP_QUAL_TAGS): the t0 and tp attributes are dropped from every read
    drop_qual_tags: bool,
    aux: Vec<u8>,
    entries: Vec<u8>,
}

impl UltimaRoutines {
    pub fn new(drop_qual_tags: bool) -> Self {
        UltimaRoutines { drop_qual_tags, aux: Vec::with_capacity(1024), entries: Vec::with_capacity(256) }
    }

    /// UltimaRoutines.preProcessRead
    #[inline]
    pub fn pre_process_read(&mut self, rec: &mut bam::Record, stats: &mut UltimaStats) -> Result<()> {
        self.strip_attributes(rec, Some(stats))
    }

    /// UltimaRoutines.stripAttributes(record, stats): the raw quality attributes stay only on mapped primary reads
    /// that are not duplicates; UT and UQ are made for mapped reads that are not secondary, and those reads are
    /// counted by their number of low-quality bases.
    pub fn strip_attributes(&mut self, rec: &mut bam::Record, stats: Option<&mut UltimaStats>) -> Result<()> {
        let flags = rec.flags();
        let mut keep_raw_quals = flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_DUPLICATE) == 0;
        if self.drop_qual_tags {
            keep_raw_quals = false;
        }
        let form_custom_tags = flags & (FLAG_UNMAPPED | FLAG_SECONDARY) == 0;
        let low_qual_count = self.strip_attributes_with(rec, form_custom_tags, keep_raw_quals)?;
        if form_custom_tags {
            if let Some(stats) = stats {
                // UltimaBamUtils.extractLowQualCount: the count that leads the UQ attribute, 0 without one
                stats.add_low_qual_count(low_qual_count, 1);
            }
        }
        Ok(())
    }

    /// UltimaRoutines.stripAttributes(record, formCustomTags, keepUltimaRawQuals); returns the number of low-quality
    /// bases (0 when the custom attributes are not formed). Java reads the attributes (which decodes them, in the
    /// order of their binary tags), clears them and sets the kept ones again, so the record ends in htsjdk's
    /// decoded form: sorted by binary tag, integers in their smallest type.
    pub fn strip_attributes_with(&mut self, rec: &mut bam::Record, form_custom_tags: bool, keep_raw_quals: bool) -> Result<i32> {
        let mut dual = false;
        let mut low_qual_count = 0;
        self.entries.clear();
        if form_custom_tags {
            let aux = aux_bytes(rec);
            dual = is_dual_consensus_type(aux)?;
            let Some(t0) = string_attribute(aux, b"t0")? else {
                bail!(
                    "Ultima read({}) without the t0 attribute (Java: NullPointerException in UltimaBamUtils.extractT0Values)",
                    String::from_utf8_lossy(rec.qname())
                );
            };
            let seq = rec.seq();
            match form_low_qual_entries(|i| seq.encoded_base(i), base_qualities(rec), t0, &mut self.entries) {
                Ok(n) => low_qual_count = n,
                Err(e) => bail!("read({}): {e}", String::from_utf8_lossy(rec.qname())),
            }
        }

        // the kept attributes as (binary tag, start of the type byte, end); at most 7 from the read and 2 new ones
        let mut kept: [(u16, usize, usize); 9] = [(0, 0, 0); 9];
        let mut count = 0;
        {
            let aux = aux_bytes(rec);
            let mut i = 0;
            while i + 3 <= aux.len() {
                let n = field_len(&aux[i..]);
                let tag = [aux[i], aux[i + 1]];
                if REQUIRED_ATTRIBUTES.contains(&&tag) || (keep_raw_quals && RAW_QUAL_ATTRIBUTES.contains(&&tag)) {
                    let key = java_tag(&tag);
                    // a repeated tag keeps its last value
                    match kept[..count].iter_mut().find(|k| k.0 == key) {
                        Some(k) => *k = (key, i + 2, i + n),
                        None => {
                            kept[count] = (key, i + 2, i + n);
                            count += 1;
                        }
                    }
                }
                i += n;
            }
        }
        const UT: u16 = (b'T' as u16) << 8 | b'U' as u16;
        const UQ: u16 = (b'Q' as u16) << 8 | b'U' as u16;
        if form_custom_tags {
            kept[count] = (UT, 0, 0);
            count += 1;
            // a null low-quality tag sets nothing
            if !self.entries.is_empty() {
                kept[count] = (UQ, 0, 0);
                count += 1;
            }
        }
        kept[..count].sort_unstable_by_key(|k| k.0);

        let mut out = std::mem::take(&mut self.aux);
        out.clear();
        let mut result = Ok(());
        {
            let aux = aux_bytes(rec);
            for &(key, s, e) in &kept[..count] {
                if s == 0 {
                    if key == UT {
                        // ConsensusType.toString
                        out.extend_from_slice(if dual { b"UTZDUAL\0" } else { b"UTZNONE\0" });
                    } else {
                        out.extend_from_slice(b"UQZ");
                        let _ = write!(out, "{low_qual_count}=");
                        out.extend_from_slice(&self.entries);
                        out.push(0);
                    }
                } else if let Err(err) = append_decoded(&mut out, &aux[s - 2..s], &aux[s..e]) {
                    result = Err(err);
                    break;
                }
            }
        }
        if result.is_ok() {
            set_aux_bytes(rec, &out);
        }
        self.aux = out;
        result?;
        Ok(low_qual_count)
    }
}

/// UltimaStats: the number of reads by their count of low-quality bases. The file lists the counts in the
/// iteration order of Java's HashMap<Integer,Integer>, which depends on the order in which the counts were first
/// seen, so that order is kept.
#[derive(Default, Clone)]
pub struct UltimaStats {
    /// reads per low-quality base count (a Java Integer: it wraps)
    counts: Vec<i32>,
    /// the low-quality base counts in the order of their first appearance
    order: Vec<i32>,
}

impl UltimaStats {
    /// UltimaStats.addLowQualCount
    #[inline]
    pub fn add_low_qual_count(&mut self, low_qual_base_count: i32, count: i32) {
        let i = low_qual_base_count as usize;
        if i >= self.counts.len() {
            self.counts.resize(i + 1, 0);
        }
        // a key is in the Java map once it has been put, whatever its value
        if self.counts[i] == 0 && !self.order.contains(&low_qual_base_count) {
            self.order.push(low_qual_base_count);
        }
        self.counts[i] = self.counts[i].wrapping_add(count);
    }

    /// The low-quality base counts in the iteration order of the Java map.
    fn java_order(&self) -> Vec<i32> {
        let mut map = JavaIntKeyOrder::default();
        for &k in &self.order {
            map.put(k);
        }
        map.keys()
    }

    /// UltimaStats.merge: the other map's entries are put in its iteration order
    pub fn merge(&mut self, other: &UltimaStats) {
        for k in other.java_order() {
            self.add_low_qual_count(k, other.counts[k as usize]);
        }
    }

    /// UltimaStats.writeStats. The combined statistics are made by merging each reader's statistics into an empty
    /// one (ReduxApplication.run), so the rows are in the order of a map filled by `merge`.
    pub fn write_stats(&self, filename: &str) -> Result<()> {
        let mut out = String::from("LowQualBases\tCount\n");
        for k in self.java_order() {
            let _ = std::fmt::Write::write_fmt(&mut out, format_args!("{}\t{}\n", k, self.counts[k as usize]));
        }
        std::fs::write(filename, out)?;
        Ok(())
    }
}

/// The key order of a java.util.HashMap with Integer keys: by bucket of the table (16 buckets at first, doubled
/// when the size passes three quarters of it), then by insertion within a bucket. A resize keeps the relative
/// order of the keys of a bucket. A bucket that reaches 9 keys makes a table of fewer than 64 buckets double; in a
/// larger table Java turns the bucket into a tree, whose order is not modelled (see `tree_bucket`).
#[derive(Default)]
struct JavaIntKeyOrder {
    buckets: Vec<Vec<i32>>,
    size: usize,
    /// a bucket reached the size at which Java builds a tree
    tree_bucket: bool,
}

impl JavaIntKeyOrder {
    /// HashMap.hash over Integer.hashCode
    #[inline]
    fn hash(key: i32) -> u32 {
        let h = key as u32;
        h ^ (h >> 16)
    }
    fn resize(&mut self) {
        let new_cap = self.buckets.len() * 2;
        let mut buckets: Vec<Vec<i32>> = vec![Vec::new(); new_cap];
        for b in &self.buckets {
            for &k in b {
                buckets[(Self::hash(k) as usize) & (new_cap - 1)].push(k);
            }
        }
        self.buckets = buckets;
    }
    /// HashMap.putVal for a key that may be present
    fn put(&mut self, key: i32) {
        if self.buckets.is_empty() {
            self.buckets = vec![Vec::new(); 16];
        }
        let i = (Self::hash(key) as usize) & (self.buckets.len() - 1);
        if self.buckets[i].contains(&key) {
            return;
        }
        self.buckets[i].push(key);
        // TREEIFY_THRESHOLD: the ninth key of a bucket; MIN_TREEIFY_CAPACITY: 64
        if self.buckets[i].len() > 8 {
            if self.buckets.len() < 64 {
                self.resize();
            } else {
                self.tree_bucket = true;
            }
        }
        self.size += 1;
        if self.size > self.buckets.len() / 4 * 3 {
            self.resize();
        }
    }
    fn keys(&self) -> Vec<i32> {
        if self.tree_bucket {
            crate::warn!("Ultima stats: a hash bucket of the Java map would be a tree; the row order of the file may differ from Java");
        }
        self.buckets.iter().flatten().copied().collect()
    }
}

#[cfg(test)]
mod tests;
