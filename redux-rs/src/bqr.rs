//! Base quality recalibration (Java package redux.bqr): counts of read bases by trinucleotide context, alt, base quality and
//! consensus type over the BQR regions, written as `<sample>.redux.bqr.tsv` (hmf-common BqrFile).
//!
//! Layout compared with Java:
//! - BqrRegionReader keeps the Java state machine (100K ref blocks, purge index, max index) but stores the per-position
//!   BaseQualityData in a flat array of fixed-size slots with an arena for the AltQualityCount lists, so that no memory is
//!   allocated per read or per base.
//! - BqrKey counts are kept in a dense table (context x alt x read type x quality) per reader and merged into the shared
//!   BaseQualityResults when a region completes. All Java counts are int sums, so the merge order does not matter.

use crate::config::{Config, SequencingType};
use crate::fmt::java_fixed;
use crate::{debug, error, warn};
use crate::refgenome::RefGenome;
use crate::regions::{self, ChrBaseRegion};
use anyhow::{Context, Result};
use rust_htslib::bam;
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;

// ReduxConstants
pub const BQR_SAMPLE_SIZE: i32 = 2_000_000;
pub const BQR_MIN_MAP_QUAL: u8 = 50;
pub const BQR_CHR_END_BUFFER: i32 = 1_000_000;
pub const BQR_DUAL_AF_LOW: f64 = 0.01;
pub const BQR_DUAL_AF_HIGH: f64 = 0.075;
pub const BQR_DUAL_AD: i32 = 2;
pub const BQR_NON_DUAL_AF_LOW: f64 = 0.05;
pub const BQR_NON_DUAL_AF_HIGH: f64 = 0.125;
pub const BQR_NON_DUAL_AD: i32 = 3;
pub const BQR_INDEL_MED_QUAL_VAF: f64 = 0.04;
pub const BQR_INDEL_LOW_QUAL_VAF: f64 = 0.1;
pub const BQR_NUM_EVENTS_THRESHOLD: f64 = 0.05;
pub const BQR_WGS_BAM_SIZE_ESTIMATE: f64 = 20e9;

// BqrRegionReader
pub const REF_BASE_REGION_SIZE: i32 = 100_000;
const BASE_DATA_PURGE_POS_CHECK: i32 = 5000;
const SHORT_FRAG_BOUNDARY_NONE: i32 = -1;

// RefSequence
const REF_SEQ_BUFFER: i32 = 5;

// BaseQualAdjustment, SbxBamUtils, UltimaBamUtils
const BASE_QUAL_MINIMUM: i8 = 1;
const SBX_MEDIUM_QUAL_LOWER: i8 = 27;
const SBX_MEDIUM_QUAL_UPPER: i8 = 29;
const ULTIMA_MAX_QUAL: i8 = 35;
const PHRED_OFFSET: u8 = 33;

const DNA_BASE_BYTES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// BAM 4-bit base codes as the letters htsjdk decodes them to.
const BAM_BASE_LETTERS: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
const BAM_BASE_N: u8 = 15;

/// RefGenomeCoordinates (resources refgenome/lengths.37.tsv and lengths.38.tsv): the autosome lengths.
const AUTOSOME_LENGTHS_37: [i32; 22] = [
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
];
const AUTOSOME_LENGTHS_38: [i32; 22] = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
];

/// hmf-common ConsensusType, in ordinal order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ConsensusType {
    None = 0,
    Single = 1,
    Dual = 2,
}

impl ConsensusType {
    const ALL: [ConsensusType; 3] = [ConsensusType::None, ConsensusType::Single, ConsensusType::Dual];

    pub fn name(self) -> &'static str {
        match self {
            ConsensusType::None => "NONE",
            ConsensusType::Single => "SINGLE",
            ConsensusType::Dual => "DUAL",
        }
    }
}

/// hmf-common BqrKey. Quality is a Java (signed) byte.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BqrKey {
    pub ref_base: u8,
    pub alt: u8,
    pub tnc: [u8; 3],
    pub quality: i8,
    pub read_type: ConsensusType,
}

impl BqrKey {
    pub fn new(ref_base: u8, alt: u8, tnc: [u8; 3], quality: i8, read_type: ConsensusType) -> Self {
        BqrKey { ref_base, alt, tnc, quality, read_type }
    }
}

/// hmf-common BqrRecord
#[derive(Clone, Debug)]
pub struct BqrRecord {
    pub key: BqrKey,
    pub count: i32,
    pub recalibrated_quality: f64,
}

/// BqrRecord.compareTo: count descending, then ref, alt, read type, quality and the three context bases.
fn compare_records(a: &BqrRecord, b: &BqrRecord) -> std::cmp::Ordering {
    b.count
        .cmp(&a.count)
        .then((a.key.ref_base as i8).cmp(&(b.key.ref_base as i8)))
        .then((a.key.alt as i8).cmp(&(b.key.alt as i8)))
        .then((a.key.read_type as u8).cmp(&(b.key.read_type as u8)))
        .then(a.key.quality.cmp(&b.key.quality))
        .then((a.key.tnc[0] as i8).cmp(&(b.key.tnc[0] as i8)))
        .then((a.key.tnc[1] as i8).cmp(&(b.key.tnc[1] as i8)))
        .then((a.key.tnc[2] as i8).cmp(&(b.key.tnc[2] as i8)))
}

// ---- dense table of BqrKey counts ----

const QUAL_SLOTS: usize = 128;
const KEY_TABLE_SIZE: usize = 64 * 4 * 3 * QUAL_SLOTS;

#[inline]
fn dna_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Index of a valid key: `tnc_index` is b0*16 + b1*4 + b2, `quality` is 0-127.
#[inline]
fn key_index(tnc_index: usize, alt_index: usize, read_type: usize, quality: usize) -> usize {
    ((tnc_index * 4 + alt_index) * 3 + read_type) * QUAL_SLOTS + quality
}

fn key_from_index(index: usize) -> BqrKey {
    let quality = index % QUAL_SLOTS;
    let rest = index / QUAL_SLOTS;
    let read_type = ConsensusType::ALL[rest % 3];
    let rest = rest / 3;
    let alt = DNA_BASE_BYTES[rest % 4];
    let t = rest / 4;
    let tnc = [DNA_BASE_BYTES[t / 16], DNA_BASE_BYTES[(t / 4) % 4], DNA_BASE_BYTES[t % 4]];
    BqrKey::new(tnc[1], alt, tnc, quality as i8, read_type)
}

/// BaseQualityResults: the combined counts of all readers.
struct BaseQualityResults {
    counts: Vec<i32>,
    total_reads_used: i64,
    total_alts_filtered: i32,
}

/// BaseQualRecalibration: shared by all partitions; collects the results and writes `<sample>.redux.bqr.tsv`.
pub struct BaseQualRecalibration {
    enabled: bool,
    write_plot: bool,
    sequencing_type: SequencingType,
    ref_genome_file: String,
    output_dir: String,
    sample_id: String,
    regions: Vec<ChrBaseRegion>,
    results: Mutex<BaseQualityResults>,
}

impl BaseQualRecalibration {
    pub fn new(config: &Config, _genome: &RefGenome) -> Result<Self> {
        // BqrConfig: the bqr_debug_* options only add debug logging in Java and are ignored
        let enabled = !config.args.has_flag("bqr_disable");
        let mut bqr = BaseQualRecalibration {
            enabled,
            write_plot: config.args.has_flag("bqr_write_plot"),
            sequencing_type: config.sequencing_type,
            ref_genome_file: config.ref_genome_file.clone(),
            output_dir: config.output_dir.clone(),
            sample_id: config.sample_id.clone(),
            regions: Vec::new(),
            results: Mutex::new(BaseQualityResults { counts: vec![0; KEY_TABLE_SIZE], total_reads_used: 0, total_alts_filtered: 0 }),
        };
        if enabled {
            bqr.regions = form_regions(config);
        }
        Ok(bqr)
    }

    pub fn enabled(&self) -> bool {
        self.enabled
    }

    #[allow(dead_code)]
    pub fn regions(&self) -> &[ChrBaseRegion] {
        &self.regions
    }

    /// One reader per partition (Java: one BqrRegionReader per BamWriter); a reader is used by one thread at a time.
    pub fn region_reader(&self) -> BqrRegionReader<'_> {
        BqrRegionReader::new(self)
    }

    /// BaseQualRecalibration.finalise: converts the combined counts to records and writes the file.
    pub fn finalise(&self) -> Result<()> {
        let results = self.results.lock().unwrap();
        debug!(
            "writing base quality recalibration, readsUsed({}) altsFiltered({})",
            results.total_reads_used, results.total_alts_filtered
        );
        // a key is taken as present when its count is non-zero; this differs from the Java map only if an int count has
        // wrapped around to exactly zero
        let mut all_counts: HashMap<BqrKey, i32> = HashMap::new();
        for (index, &count) in results.counts.iter().enumerate() {
            if count != 0 {
                all_counts.insert(key_from_index(index), count);
            }
        }
        drop(results);
        let mut records = convert_to_records(&all_counts);

        // BqrFile.generateFilename: base path + separator + sample + ".redux.bqr.tsv" (no output id); the output directory
        // already ends with the separator
        let filename = format!("{}{}.redux.bqr.tsv", self.output_dir, self.sample_id);
        write_bqr_file(&filename, &mut records).with_context(|| format!("failed to write base recalibration {filename}"))?;

        if self.write_plot {
            warn!("the BQR plot (bqr_write_plot) is not produced by this port");
        }
        Ok(())
    }
}

/// BaseQualRecalibration.assumeWgsBams: files that cannot be read count as empty.
fn assume_wgs_bams(config: &Config) -> bool {
    let mut total: i64 = 0;
    for f in &config.bam_files {
        if let Ok(m) = std::fs::metadata(f) {
            total += m.len() as i64;
        }
    }
    total as f64 > BQR_WGS_BAM_SIZE_ESTIMATE
}

/// BaseQualRecalibration.formRegions: the specific regions if given, otherwise one region per autosome, either the
/// chromosome less 1M bases at each end, or for large (WGS) inputs the last 2M bases before the end buffer.
fn form_regions(config: &Config) -> Vec<ChrBaseRegion> {
    let mut result = Vec::new();
    if !config.specific.regions.is_empty() {
        result.extend(config.specific.regions.iter().cloned());
        return result;
    }
    let use_all_regions = config.args.has_flag("bqr_use_all_regions") || !assume_wgs_bams(config);
    let lengths = if config.v38 { &AUTOSOME_LENGTHS_38 } else { &AUTOSOME_LENGTHS_37 };
    for (i, chr_length) in lengths.iter().enumerate() {
        let chr = regions::versioned_chromosome(config.v38, regions::HUMAN_CHROMOSOMES[i]);
        if !config.specific.include_chromosome(&chr) {
            continue;
        }
        let (start, end);
        if use_all_regions {
            start = BQR_CHR_END_BUFFER;
            end = chr_length - BQR_CHR_END_BUFFER;
        } else {
            end = chr_length - BQR_CHR_END_BUFFER;
            start = end - BQR_SAMPLE_SIZE;
        }
        result.push(ChrBaseRegion::new(&chr, start, end));
    }
    result
}

/// Java Math.max for values that are not NaN or zeros of different sign.
#[inline]
fn jmax(a: f64, b: f64) -> f64 {
    if a >= b { a } else { b }
}

/// BaseQualRecalibration.convertToRecords: a record per observed key with its recalibrated quality, plus a zero-count
/// record for each alt that was not observed at an observed context, quality and read type.
pub fn convert_to_records(all_quality_counts: &HashMap<BqrKey, i32>) -> Vec<BqrRecord> {
    const NO_BASE_OR_QUAL: u8 = 1;
    let mut result = Vec::new();

    let ref_count_map: HashMap<BqrKey, i32> =
        all_quality_counts.iter().filter(|(k, _)| k.ref_base == k.alt).map(|(k, v)| (*k, *v)).collect();

    // totals per context and read type across all entries
    let tri_nuc_key = |tnc: [u8; 3], read_type: ConsensusType| BqrKey::new(NO_BASE_OR_QUAL, NO_BASE_OR_QUAL, tnc, NO_BASE_OR_QUAL as i8, read_type);
    let mut tri_nuc_map: HashMap<BqrKey, i32> = HashMap::new();
    for (key, &count) in all_quality_counts {
        if key.quality == 0 {
            continue;
        }
        let e = tri_nuc_map.entry(tri_nuc_key(key.tnc, key.read_type)).or_insert(0);
        *e = e.wrapping_add(count);
    }

    // BaseQualRecalibration.calcRecalibratedQual:
    // errorRate = (# alt TN sites / # ref TN sites) * (# alts at ref TN sites) / (# alts at ref TN sites + # refs at alt TN sites)
    let calc_recalibrated_qual = |key: &BqrKey, observed_count: i32| -> f64 {
        let alt_tnc = [key.tnc[0], key.alt, key.tnc[2]];
        let alt_key = BqrKey::new(key.alt, key.alt, alt_tnc, key.quality, key.read_type);
        let alt_ref_count = ref_count_map.get(&alt_key).copied().unwrap_or(0);
        let alt_tri_nuc_count = tri_nuc_map.get(&tri_nuc_key(alt_tnc, key.read_type)).copied().unwrap_or(0);
        let ref_tri_nuc_count = tri_nuc_map.get(&tri_nuc_key(key.tnc, key.read_type)).copied().unwrap_or(0);
        let tri_nuc_rate = if ref_tri_nuc_count > 0 { alt_tri_nuc_count as f64 / ref_tri_nuc_count as f64 } else { 0.0 };
        let denom = observed_count.wrapping_add(alt_ref_count);
        if tri_nuc_rate == 0.0 || denom == 0 {
            return 0.0;
        }
        let p = tri_nuc_rate * observed_count as f64 / denom as f64;
        let calc_probability = if p <= 1.0 { p } else { 1.0 };
        // BaseQualAdjustment.probabilityToPhredQual; a probability of 1 gives -0.0, which Java prints as "-0.00"
        -10.0 * calc_probability.log10()
    };

    let mut synthetic_alt_keys: HashSet<BqrKey> = HashSet::new();
    for (key, &count) in all_quality_counts {
        if key.quality == 0 {
            continue;
        }
        let recalibrated_qual = if key.alt == key.ref_base { key.quality as f64 } else { calc_recalibrated_qual(key, count) };
        result.push(BqrRecord { key: *key, count, recalibrated_quality: recalibrated_qual });

        // add alt entries which have no results
        for alt in DNA_BASE_BYTES {
            if alt == key.ref_base {
                continue;
            }
            let alt_key = BqrKey::new(key.ref_base, alt, key.tnc, key.quality, key.read_type);
            if !all_quality_counts.contains_key(&alt_key) && synthetic_alt_keys.insert(alt_key) {
                let recal_qual_min = calc_recalibrated_qual(&alt_key, 1);
                let synthetic_qual = jmax(recal_qual_min + 10.0 * 2f64.log10(), key.quality as f64);
                result.push(BqrRecord { key: alt_key, count: 0, recalibrated_quality: synthetic_qual });
            }
        }
    }
    result
}

/// BqrFile.write: sorts the records and writes them. The sort keys cover every field of the key, so the order does not
/// depend on the order of the input.
fn write_bqr_file(filename: &str, records: &mut [BqrRecord]) -> Result<()> {
    use std::io::Write;
    records.sort_by(compare_records);
    let mut out = std::io::BufWriter::new(std::fs::File::create(filename)?);
    out.write_all(b"Alt\tRef\tTrinucleotideContext\tReadType\tCount\tOriginalQual\tRecalibratedQual\n")?;
    for r in records.iter() {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.key.alt as char,
            r.key.ref_base as char,
            String::from_utf8_lossy(&r.key.tnc),
            r.key.read_type.name(),
            r.count,
            r.key.quality,
            java_fixed(r.recalibrated_quality, 2)
        )?;
    }
    out.flush()?;
    Ok(())
}

// ---- per-position data ----

const NO_ENTRY: u32 = u32::MAX;
const FLAG_PRESENT: u8 = 1;
const FLAG_HIGH_QUAL_INDEL: u8 = 2;
const FLAG_MEDIUM_QUAL_INDEL: u8 = 4;
const REF_CODE_SHIFT: u32 = 4;

/// Number of (consensus type, quality) pairs whose ref-matching bases are counted in the dense rows.
const MAX_COMBOS: usize = 16;
/// Counters per position in the dense rows: positive and negative strand for each pair.
const DENSE_WIDTH: usize = 2 * MAX_COMBOS;
const COMBO_UNSET: u8 = 255;
const COMBO_NONE: u8 = 254;

/// BaseQualityData for one position. `flags` without FLAG_PRESENT is Java's null entry; its upper four bits hold the BAM
/// base code of the ref base (15 if the ref base is not one a read can match). Java keeps a list of AltQualityCount per
/// consensus type; here the counts of read bases which equal the ref base are held in the reader's dense rows, indexed by
/// (type, quality) pair and strand, and all other counts in the reader's arena, linked from `next`. The order within
/// Java's lists has no effect on the counts formed from them.
#[derive(Clone, Copy)]
struct PosData {
    tnc: [u8; 3],
    flags: u8,
    indel_count: i32,
    aligned_count: i32,
    next: u32,
}

const EMPTY_POS_DATA: PosData = PosData { tnc: [0; 3], flags: 0, indel_count: 0, aligned_count: 0, next: NO_ENTRY };

/// AltQualityCount: `key` packs the consensus type, the BAM base code of the alt and the quality.
#[derive(Clone, Copy)]
struct AltQualityCount {
    key: u32,
    pos_strand_count: i32,
    neg_strand_count: i32,
    next: u32,
}

#[inline]
fn alt_qual_key(read_type: u32, base_code: u32, quality: u32) -> u32 {
    (read_type << 16) | (base_code << 8) | quality
}

/// The BAM base code of a ref base, or 15 (N, which is never counted) for a byte that no read base can equal.
#[inline]
fn ref_base_code(base: u8) -> u8 {
    match base {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

/// The dense row id of a consensus type and quality seen for the first time, if any are left.
#[cold]
#[inline(never)]
fn assign_combo(combo_keys: &mut Vec<u32>, combo_ids: &mut [u8; 3 * 256], combo: usize) -> u8 {
    let id = if combo_keys.len() < MAX_COMBOS {
        combo_keys.push(alt_qual_key((combo >> 8) as u32, 0, (combo & 0xff) as u32));
        (combo_keys.len() - 1) as u8
    } else {
        COMBO_NONE
    };
    combo_ids[combo] = id;
    id
}

/// BaseQualityData.processReadBase for a base which is not held in the dense rows: finds or adds its AltQualityCount.
#[cold]
#[inline(never)]
fn add_listed_base(slot: &mut PosData, arena: &mut Vec<AltQualityCount>, key: u32, pos_strand: bool) {
    let (p, n) = if pos_strand { (1, 0) } else { (0, 1) };
    let mut index = slot.next;
    let mut last = NO_ENTRY;
    while index != NO_ENTRY {
        let e = &mut arena[index as usize];
        if e.key == key {
            e.pos_strand_count = e.pos_strand_count.wrapping_add(p);
            e.neg_strand_count = e.neg_strand_count.wrapping_add(n);
            return;
        }
        last = index;
        index = e.next;
    }
    let new_index = arena.len() as u32;
    arena.push(AltQualityCount { key, pos_strand_count: p, neg_strand_count: n, next: NO_ENTRY });
    if last == NO_ENTRY {
        slot.next = new_index;
    } else {
        arena[last as usize].next = new_index;
    }
}

/// Calls `f` with the key and strand counts of each AltQualityCount of a position.
#[inline]
fn for_each_count(
    slot: &PosData, row: &[i32], combo_keys: &[u32], arena: &[AltQualityCount], mut f: impl FnMut(u32, i32, i32),
) {
    let ref_code = (slot.flags >> REF_CODE_SHIFT) as u32;
    for (id, combo_key) in combo_keys.iter().enumerate() {
        let (pos, neg) = (row[id * 2], row[id * 2 + 1]);
        if (pos | neg) != 0 {
            f(combo_key | (ref_code << 8), pos, neg);
        }
    }
    let mut index = slot.next;
    while index != NO_ENTRY {
        let e = &arena[index as usize];
        f(e.key, e.pos_strand_count, e.neg_strand_count);
        index = e.next;
    }
}

/// RefSequence: the ref bases of the current block with a buffer of 5 bases at each end.
struct RefSequence {
    start: i32,
    end: i32,
    bases: Vec<u8>,
}

impl RefSequence {
    fn new(chromosome: &str, pos_start: i32, pos_end: i32, genome: &RefGenome) -> Self {
        let sequence_end = genome.chromosome_length(chromosome);
        let start = 1.max(pos_start - REF_SEQ_BUFFER);
        let end = sequence_end.min(pos_end + REF_SEQ_BUFFER);
        let bases = match genome.bases(chromosome, start, end) {
            Ok(b) => b,
            Err(e) => panic!("failed to read ref bases {chromosome}:{start}-{end}: {e}"),
        };
        RefSequence { start, end, bases }
    }
    #[inline]
    fn position_within_bounds(&self, position: i32) -> bool {
        position >= self.start + REF_SEQ_BUFFER && position <= self.end - REF_SEQ_BUFFER
    }
    #[inline]
    fn length(&self) -> usize {
        self.bases.len()
    }
}

/// The tags of a read which the BQR uses, found in one pass over the record's tag data.
#[derive(Default)]
struct ReadTags<'a> {
    /// NM: SamRecordUtils.getNumEvents
    num_events: Option<i32>,
    /// UT: SamRecordUtils.extractConsensusType
    consensus_type: Option<&'a [u8]>,
    /// MC: the mate CIGAR
    mate_cigar: Option<&'a [u8]>,
    /// YX: SbxBamUtils.extractDuplexBaseIndex
    duplex_base_index: Option<i32>,
    /// t0: UltimaBamUtils.extractT0Values, still with the Phred offset
    t0: Option<&'a [u8]>,
}

/// Reads the tags from the BAM record data: each is a two-byte name, a type byte and a value. Integer tags are those
/// SAMRecord.getIntegerAttribute returns and string tags those of type Z.
fn read_tags(rec: &bam::Record) -> ReadTags<'_> {
    let inner = rec.inner();
    let core = &inner.core;
    // SAFETY: htslib holds l_data bytes of record data at `data`, which live as long as the record is borrowed
    let data: &[u8] = if inner.data.is_null() { &[] } else { unsafe { std::slice::from_raw_parts(inner.data, inner.l_data as usize) } };
    let l_qseq = core.l_qseq as usize;
    let aux_start = core.l_qname as usize + core.n_cigar as usize * 4 + l_qseq.div_ceil(2) + l_qseq;
    let mut aux = data.get(aux_start..).unwrap_or(&[]);
    let mut tags = ReadTags::default();
    while aux.len() >= 4 {
        let tag = [aux[0], aux[1]];
        let value_type = aux[2];
        let rest = &aux[3..];
        let size = match value_type {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'd' => 8,
            b'Z' | b'H' => match rest.iter().position(|b| *b == 0) {
                Some(n) => n + 1,
                None => break,
            },
            b'B' => {
                if rest.len() < 5 {
                    break;
                }
                let element_size = match rest[0] {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => break,
                };
                5 + element_size * u32::from_le_bytes([rest[1], rest[2], rest[3], rest[4]]) as usize
            }
            _ => break,
        };
        if rest.len() < size {
            break;
        }
        let value = &rest[..size];
        let integer = || -> Option<i32> {
            match value_type {
                b'c' => Some(value[0] as i8 as i32),
                b'C' => Some(value[0] as i32),
                b's' => Some(i16::from_le_bytes([value[0], value[1]]) as i32),
                b'S' => Some(u16::from_le_bytes([value[0], value[1]]) as i32),
                b'i' | b'I' => Some(i32::from_le_bytes([value[0], value[1], value[2], value[3]])),
                _ => None,
            }
        };
        let string = || -> Option<&[u8]> { if value_type == b'Z' { Some(&value[..size - 1]) } else { None } };
        match &tag {
            b"NM" if tags.num_events.is_none() => tags.num_events = integer(),
            b"YX" if tags.duplex_base_index.is_none() => tags.duplex_base_index = integer(),
            b"UT" if tags.consensus_type.is_none() => tags.consensus_type = string(),
            b"MC" if tags.mate_cigar.is_none() => tags.mate_cigar = string(),
            b"t0" if tags.t0.is_none() => tags.t0 = string(),
            _ => {}
        }
        aux = &rest[size..];
    }
    tags
}

/// CigarUtils.getReadBoundaryPosition(readStart, cigar, false, false) as used by SamRecordUtils.getAlignmentEndFromCigar.
/// As in Java only M, D, N and S end an element, so the lengths of other operators (I, =, X, H) are dropped.
fn alignment_end_from_cigar(read_start: i32, cigar: &[u8]) -> i32 {
    let mut current_position = read_start;
    let mut element_length: i32 = 0;
    for &c in cigar {
        if c == b'D' || c == b'M' || c == b'S' || c == b'N' {
            if c == b'S' && read_start == current_position {
                // left clip
            } else if c == b'S' {
                break;
            } else {
                current_position = current_position.wrapping_add(element_length);
            }
            element_length = 0;
            continue;
        }
        if c.is_ascii_digit() {
            element_length = element_length.wrapping_mul(10).wrapping_add((c - b'0') as i32);
        } else {
            element_length = 0;
        }
    }
    current_position - 1
}

/// BqrRegionReader: gathers the counts for the reads of one partition. The reads of each region arrive in coordinate order.
pub struct BqrRegionReader<'a> {
    shared: &'a BaseQualRecalibration,
    sequencing_type: SequencingType,
    genome: Option<RefGenome>,

    has_active_region: bool,
    /// span of the BQR regions which overlap the current partition region (Java mPartitionOverlapRegion)
    overlap_chromosome: String,
    overlap_start: i32,
    overlap_end: i32,
    ref_seq: Option<RefSequence>,

    /// Java mBaseQualityData, indexed by position less the ref sequence start; `data_allocated` is false while Java's array is null
    slots: Vec<PosData>,
    data_allocated: bool,
    /// DENSE_WIDTH counters per slot: the ref-matching bases by (type, quality) pair and strand; all zero for an empty slot
    dense: Vec<i32>,
    /// the pair id of each consensus type and quality: COMBO_UNSET until seen, COMBO_NONE once all ids are taken
    combo_ids: [u8; 3 * 256],
    /// the type and quality of each pair id, as an AltQualityCount key without the base
    combo_keys: Vec<u32>,
    arena: Vec<AltQualityCount>,
    /// lowest and highest slot index created since the slots were last cleared
    created_min: usize,
    created_max: usize,

    next_purge_position: i32,
    purge_index: i32,
    max_index: i32,

    min_read_start_position: i32,
    max_read_end_position: i32,
    current_read_type: ConsensusType,

    total_reads_used: i64,
    total_alts_filtered: i32,

    /// counts and totals not yet added to the shared results
    counts: Vec<i32>,
    pending_reads_used: i64,
    pending_alts_filtered: i32,
    pending: bool,
}

impl<'a> BqrRegionReader<'a> {
    fn new(shared: &'a BaseQualRecalibration) -> Self {
        BqrRegionReader {
            shared,
            sequencing_type: shared.sequencing_type,
            genome: None,
            has_active_region: false,
            overlap_chromosome: String::new(),
            overlap_start: 0,
            overlap_end: 0,
            ref_seq: None,
            slots: Vec::new(),
            data_allocated: false,
            dense: Vec::new(),
            combo_ids: [COMBO_UNSET; 3 * 256],
            combo_keys: Vec::with_capacity(MAX_COMBOS),
            arena: Vec::new(),
            created_min: usize::MAX,
            created_max: 0,
            next_purge_position: 0,
            purge_index: 0,
            max_index: 0,
            min_read_start_position: SHORT_FRAG_BOUNDARY_NONE,
            max_read_end_position: SHORT_FRAG_BOUNDARY_NONE,
            current_read_type: ConsensusType::None,
            total_reads_used: 0,
            total_alts_filtered: 0,
            counts: vec![0; KEY_TABLE_SIZE],
            pending_reads_used: 0,
            pending_alts_filtered: 0,
            pending: false,
        }
    }

    pub fn is_active(&self) -> bool {
        self.has_active_region
    }

    /// BqrRegionReader.initialise: finds the BQR regions which overlap this partition region. As in Java the current ref
    /// sequence is only dropped when there is no overlap, so a block of the previous chromosome stays current until a
    /// read starts outside of its positions, and a region on a non-human contig leaves the state unchanged.
    pub fn initialise(&mut self, region: &ChrBaseRegion) {
        if !regions::is_human_chromosome(&region.chromosome) {
            return;
        }
        // ChrBaseRegion.overlaps compares the chromosomes as HumanChromosome, ie without regard to a chr prefix
        let rank = regions::chromosome_rank(&region.chromosome);
        let shared = self.shared;
        let mut overlapped =
            shared.regions.iter().filter(|x| regions::chromosome_rank(&x.chromosome) == rank && x.start <= region.end && x.end >= region.start);
        let Some(first) = overlapped.next() else {
            self.has_active_region = false;
            self.ref_seq = None;
            return;
        };
        let last = overlapped.last().unwrap_or(first);
        if self.overlap_chromosome != first.chromosome {
            self.overlap_chromosome.clear();
            self.overlap_chromosome.push_str(&first.chromosome);
        }
        self.overlap_start = first.start;
        self.overlap_end = last.end;
        self.has_active_region = true;
    }

    /// BqrRegionReader.checkLocalRegion: moves to a new block of up to 100K ref bases starting at this read.
    fn check_local_region(&mut self, position: i32) {
        if let Some(rs) = &self.ref_seq {
            if rs.position_within_bounds(position) {
                return;
            }
        }
        if position < self.overlap_start {
            return;
        }
        if self.has_active_region {
            self.process_completed_local_region();
        }
        if position > self.overlap_end {
            // the reads are now past the BQR region(s)
            self.has_active_region = false;
            return;
        }
        let ref_base_start = position;
        let ref_base_end = (position + REF_BASE_REGION_SIZE).min(self.overlap_end);
        if ref_base_end <= ref_base_start {
            return;
        }
        if self.genome.is_none() {
            match RefGenome::open(&self.shared.ref_genome_file) {
                Ok(g) => self.genome = Some(g),
                Err(e) => panic!("failed to open ref genome {}: {e}", self.shared.ref_genome_file),
            }
        }
        let ref_seq = RefSequence::new(&self.overlap_chromosome, ref_base_start, ref_base_end, self.genome.as_ref().unwrap());
        self.next_purge_position = ref_base_start + BASE_DATA_PURGE_POS_CHECK;
        self.allocate_slots(ref_seq.length());
        self.ref_seq = Some(ref_seq);
    }

    /// Java allocates a new array when the length differs; the slots are always empty at this point.
    fn allocate_slots(&mut self, length: usize) {
        if !self.data_allocated || self.slots.len() != length {
            self.slots.clear();
            self.slots.resize(length, EMPTY_POS_DATA);
            self.dense.clear();
            self.dense.resize(length * DENSE_WIDTH, 0);
            self.created_min = usize::MAX;
            self.created_max = 0;
        }
        self.data_allocated = true;
    }

    /// BqrRegionReader.processCompletedLocalRegion: converts the remaining positions to counts, hands them to the results
    /// and clears the block.
    fn process_completed_local_region(&mut self) {
        if !self.data_allocated {
            return;
        }
        self.build_quality_counts();

        // BaseQualityResults.addBaseQualityRegionCounter: Java adds the reader's running totals each time, so they are
        // counted repeatedly; they only appear in a debug log
        self.pending_reads_used = self.pending_reads_used.wrapping_add(self.total_reads_used);
        self.pending_alts_filtered = self.pending_alts_filtered.wrapping_add(self.total_alts_filtered);
        self.pending = true;

        if self.created_min <= self.created_max {
            let hi = self.created_max.min(self.slots.len().saturating_sub(1));
            for index in self.created_min..=hi {
                if self.slots[index].flags & FLAG_PRESENT != 0 {
                    self.clear_slot(index);
                }
            }
        }
        self.created_min = usize::MAX;
        self.created_max = 0;
        self.arena.clear();
        self.next_purge_position = 0;
        self.purge_index = 0;
        self.max_index = 0;
    }

    /// Java sets the array entry to null.
    #[inline]
    fn clear_slot(&mut self, index: usize) {
        self.slots[index] = EMPTY_POS_DATA;
        self.dense[index * DENSE_WIDTH..(index + 1) * DENSE_WIDTH].fill(0);
    }

    /// BqrRegionReader.onRegionComplete; also adds this reader's counts to the shared results.
    pub fn on_region_complete(&mut self) {
        self.process_completed_local_region();
        self.merge_into_results();
    }

    fn merge_into_results(&mut self) {
        if !self.pending {
            return;
        }
        let mut results = self.shared.results.lock().unwrap();
        for (total, count) in results.counts.iter_mut().zip(self.counts.iter_mut()) {
            *total = total.wrapping_add(*count);
            *count = 0;
        }
        results.total_reads_used = results.total_reads_used.wrapping_add(self.pending_reads_used);
        results.total_alts_filtered = results.total_alts_filtered.wrapping_add(self.pending_alts_filtered);
        self.pending_reads_used = 0;
        self.pending_alts_filtered = 0;
        self.pending = false;
    }

    /// BqrRegionReader.buildQualityCounts: positions from the purge index up to the max index. The max index is that of
    /// the most recently created position (see create_slot), so data at higher indices is dropped as in Java.
    fn build_quality_counts(&mut self) {
        let mut i = self.purge_index;
        while i <= self.max_index {
            if let Some(slot) = self.slots.get(i as usize) {
                if slot.flags & FLAG_PRESENT != 0 {
                    self.convert_to_key_counts(i as usize);
                }
            }
            i += 1;
        }
    }

    /// BqrRegionReader.convertToKeyCounts with BaseQualityData.formKeyCounts: applies the indel and alt filters to one
    /// position and adds its counts by key. Only valid keys are kept (context and alt of A, C, G, T).
    fn convert_to_key_counts(&mut self, slot_index: usize) {
        let slot = self.slots[slot_index];
        if slot.flags & FLAG_HIGH_QUAL_INDEL != 0 {
            return;
        }
        let indel_vaf_threshold = if slot.flags & FLAG_MEDIUM_QUAL_INDEL != 0 { BQR_INDEL_MED_QUAL_VAF } else { BQR_INDEL_LOW_QUAL_VAF };
        // BaseQualityData.indelVaf
        let denom = slot.indel_count.wrapping_add(slot.aligned_count) as f64;
        let indel_vaf = if denom > 0.0 { slot.indel_count as f64 / denom } else { 0.0 };
        if indel_vaf > indel_vaf_threshold {
            return;
        }
        let row = &self.dense[slot_index * DENSE_WIDTH..(slot_index + 1) * DENSE_WIDTH];
        let combo_keys = &self.combo_keys[..];
        let arena = &self.arena[..];

        // exclude any alt with too much support (regardless of quality); alts are indexed by BAM base code
        let ref_base = slot.tnc[1];
        let mut failed_alts = 0u16;
        if slot.next != NO_ENTRY {
            let mut standard_alt_counts = [0i32; 16];
            let mut dual_alt_counts = [0i32; 16];
            let (mut standard_alts, mut dual_alts, mut all_alts) = (0u16, 0u16, 0u16);
            let (mut standard_total, mut dual_total) = (0i32, 0i32);
            for_each_count(&slot, row, combo_keys, arena, |key, pos, neg| {
                let is_dual = (key >> 16) == ConsensusType::Dual as u32;
                let code = ((key >> 8) & 0xf) as usize;
                let total = pos.wrapping_add(neg);
                if is_dual {
                    dual_total = dual_total.wrapping_add(total);
                } else {
                    standard_total = standard_total.wrapping_add(total);
                }
                if BAM_BASE_LETTERS[code] == ref_base {
                    return;
                }
                let bit = 1u16 << code;
                all_alts |= bit;
                if is_dual {
                    dual_alts |= bit;
                    dual_alt_counts[code] = dual_alt_counts[code].wrapping_add(total);
                } else {
                    standard_alts |= bit;
                    standard_alt_counts[code] = standard_alt_counts[code].wrapping_add(total);
                }
            });

            // BaseQualityData.passesAltTest
            let passes_alt_test = |alt_count: i32, alt_vaf: f64, low_af: f64, high_af: f64, ad_limit: i32| -> bool {
                alt_vaf < low_af || (alt_vaf < high_af && alt_count <= ad_limit)
            };
            for code in 0..16usize {
                let bit = 1u16 << code;
                if all_alts & bit == 0 {
                    continue;
                }
                let mut include_alt = true;
                if standard_alts & bit != 0 {
                    let count = standard_alt_counts[code];
                    let vaf = count as f64 / standard_total as f64;
                    include_alt = passes_alt_test(count, vaf, BQR_NON_DUAL_AF_LOW, BQR_NON_DUAL_AF_HIGH, BQR_NON_DUAL_AD);
                }
                if include_alt && dual_alts & bit != 0 {
                    let count = dual_alt_counts[code];
                    let vaf = count as f64 / dual_total as f64;
                    include_alt = passes_alt_test(count, vaf, BQR_DUAL_AF_LOW, BQR_DUAL_AF_HIGH, BQR_DUAL_AD);
                }
                if !include_alt {
                    failed_alts |= bit;
                }
            }
            self.total_alts_filtered = self.total_alts_filtered.wrapping_add(failed_alts.count_ones() as i32);
        }

        // BqrKey.isValid: the context must be A, C, G, T; the key of a negative strand count is the reverse complement
        let (Some(b0), Some(b1), Some(b2)) = (dna_index(slot.tnc[0]), dna_index(slot.tnc[1]), dna_index(slot.tnc[2])) else {
            return;
        };
        let tnc_index = b0 * 16 + b1 * 4 + b2;
        let tnc_index_reversed = (3 - b2) * 16 + (3 - b1) * 4 + (3 - b0);
        let counts = &mut self.counts;
        for_each_count(&slot, row, combo_keys, arena, |key, pos, neg| {
            let code = ((key >> 8) & 0xf) as usize;
            if failed_alts & (1u16 << code) != 0 {
                return;
            }
            let Some(alt_index) = dna_index(BAM_BASE_LETTERS[code]) else { return };
            let read_type = (key >> 16) as usize;
            let quality = (key & 0xff) as usize;
            if quality >= QUAL_SLOTS {
                return; // not reached: negative Java qualities are never counted
            }
            if pos > 0 {
                let c = &mut counts[key_index(tnc_index, alt_index, read_type, quality)];
                *c = c.wrapping_add(pos);
            }
            if neg > 0 {
                let c = &mut counts[key_index(tnc_index_reversed, 3 - alt_index, read_type, quality)];
                *c = c.wrapping_add(neg);
            }
        });
    }

    /// BqrRegionReader.processRecord. The caller (BamWriter.captureReadInfo) has excluded unmapped, duplicate, secondary
    /// and supplementary reads.
    pub fn process_record(&mut self, rec: &bam::Record) {
        let alignment_start = rec.pos() as i32 + 1;
        self.check_local_region(alignment_start);

        if !self.has_active_region || self.ref_seq.is_none() {
            return;
        }
        if rec.mapq() < BQR_MIN_MAP_QUAL {
            return;
        }
        let tags = read_tags(rec);
        let num_events = tags.num_events.unwrap_or(0);
        if num_events as f64 / rec.seq_len() as f64 > BQR_NUM_EVENTS_THRESHOLD {
            return;
        }
        self.set_short_fragment_boundaries(rec, alignment_start, tags.mate_cigar);

        // SamRecordUtils.extractConsensusType
        self.current_read_type = match tags.consensus_type {
            None | Some(b"NONE") => ConsensusType::None,
            Some(b"SINGLE") => ConsensusType::Single,
            Some(b"DUAL") => ConsensusType::Dual,
            Some(other) => panic!("invalid consensus type: {}", String::from_utf8_lossy(other)),
        };

        self.traverse_cigar(rec, alignment_start, &tags);

        if alignment_start >= self.next_purge_position {
            self.purge_base_data_list(alignment_start);
        }
    }

    /// BqrRegionReader.setShortFragmentBoundaries: limits the positions used so as not to read past the insert on the 3'
    /// end of short fragments. As in Java the boundaries are not reset for an unpaired read, which keeps those of the
    /// previous paired read.
    fn set_short_fragment_boundaries(&mut self, rec: &bam::Record, alignment_start: i32, mate_cigar: Option<&[u8]>) {
        if !rec.is_paired() {
            return;
        }
        self.min_read_start_position = SHORT_FRAG_BOUNDARY_NONE;
        self.max_read_end_position = SHORT_FRAG_BOUNDARY_NONE;

        let insert_size = rec.insert_size() as i32;
        if insert_size == 0 || insert_size.wrapping_abs() as i64 >= rec.seq_len() as i64 {
            return;
        }
        if rec.is_reverse() == rec.is_mate_reverse() {
            return;
        }
        let mate_start = rec.mpos() as i32 + 1;
        if rec.is_reverse() {
            self.min_read_start_position = mate_start;
        } else if let Some(mate_cigar) = mate_cigar {
            self.max_read_end_position = alignment_end_from_cigar(mate_start, mate_cigar);
        } else {
            self.max_read_end_position = alignment_start + insert_size.wrapping_abs() - 1;
        }
    }

    /// CigarHandler.traverseCigar with the handlers of BqrRegionReader.
    fn traverse_cigar(&mut self, rec: &bam::Record, alignment_start: i32, tags: &ReadTags) {
        let cigar = rec.raw_cigar();
        let quals = rec.qual();
        // UltimaBamUtils.extractT0Values: the t0 tag holds a quality per base with the Phred offset; Java fails on a read
        // without it as soon as it is needed
        let t0 = tags.t0;
        let mut read_index: i32 = 0;
        let mut ref_base: i32 = alignment_start;
        for &c in cigar {
            let length = (c >> 4) as i32;
            match c & 0xf {
                5 | 6 => {} // H, P
                4 => read_index += length, // S
                3 => ref_base += length,   // N
                2 => {
                    // D: handleDelete with the read index and ref position before the delete
                    let before = read_index - 1;
                    let mut min_qual = quals[(before + 1) as usize] as i8 as i32;
                    if self.sequencing_type == SequencingType::Ultima {
                        // a low qual complete homopolymer deletion is implied by the bases next to the delete
                        let t0 = t0.expect("Ultima read without t0 tag");
                        for i in 0..=1 {
                            min_qual = min_qual.min(t0[(before + i) as usize].wrapping_sub(PHRED_OFFSET) as i8 as i32);
                        }
                    }
                    self.mark_indel_position(ref_base - 1, min_qual as i8);
                    ref_base += length;
                }
                1 => {
                    // I: handleInsert with the read index before the insert, or zero for an insert at the read start,
                    // in which case Java looks at the quals one base further on
                    let before = if read_index == 0 { 0 } else { read_index - 1 };
                    let mut min_qual = quals[(before + 1) as usize] as i8 as i32;
                    for i in 1..length {
                        min_qual = min_qual.min(quals[(before + i + 1) as usize] as i8 as i32);
                    }
                    self.mark_indel_position(ref_base - 1, min_qual as i8);
                    read_index += length;
                }
                0 | 7 | 8 => {
                    // M, =, X
                    self.handle_alignment(rec, quals, tags, length, read_index, ref_base);
                    read_index += length;
                    ref_base += length;
                }
                op => panic!("unhandled CIGAR operator code {op}"),
            }
        }
    }

    /// BqrRegionReader.getOrCreateBaseQualData with the context taken from the ref sequence.
    #[inline(always)]
    fn get_or_create(&mut self, position: i32) -> usize {
        let ref_seq = self.ref_seq.as_ref().unwrap();
        let index = (position - ref_seq.start) as usize;
        if self.slots[index].flags & FLAG_PRESENT == 0 {
            let tnc = [ref_seq.bases[index - 1], ref_seq.bases[index], ref_seq.bases[index + 1]];
            self.create_slot(index, tnc);
        }
        index
    }

    /// A new BaseQualityData. The max index is set to the index of each newly created position, even when lower than
    /// before, as in Java.
    #[inline]
    fn create_slot(&mut self, index: usize, tnc: [u8; 3]) {
        let flags = FLAG_PRESENT | (ref_base_code(tnc[1]) << REF_CODE_SHIFT);
        self.slots[index] = PosData { tnc, flags, ..EMPTY_POS_DATA };
        self.max_index = index as i32;
        self.created_min = self.created_min.min(index);
        self.created_max = self.created_max.max(index);
    }

    /// BaseQualityData.processReadBase: one more read base of this type, base and quality on the given strand.
    #[inline(always)]
    fn count_read_base(&mut self, index: usize, read_type: ConsensusType, base_code: u8, quality: i8, pos_strand: bool) {
        let slot = &mut self.slots[index];
        if base_code == slot.flags >> REF_CODE_SHIFT {
            let combo = ((read_type as usize) << 8) | quality as u8 as usize;
            let mut id = self.combo_ids[combo];
            if id == COMBO_UNSET {
                id = assign_combo(&mut self.combo_keys, &mut self.combo_ids, combo);
            }
            if id != COMBO_NONE {
                let c = &mut self.dense[index * DENSE_WIDTH + id as usize * 2 + (!pos_strand) as usize];
                *c = c.wrapping_add(1);
                return;
            }
        }
        let key = alt_qual_key(read_type as u32, base_code as u32, quality as u8 as u32);
        add_listed_base(&mut self.slots[index], &mut self.arena, key, pos_strand);
    }

    /// BqrRegionReader.markIndelPosition: `position` is the last ref base before the indel.
    fn mark_indel_position(&mut self, position: i32, min_qual: i8) {
        if !self.ref_seq.as_ref().unwrap().position_within_bounds(position) {
            return;
        }
        let index = self.get_or_create(position);
        let slot = &mut self.slots[index];
        slot.indel_count = slot.indel_count.wrapping_add(1);

        // BaseQualAdjustment.isHighBaseQual and isMediumBaseQual; every Illumina indel counts as high quality
        let (is_high_qual, is_medium_qual) = match self.sequencing_type {
            SequencingType::Illumina => (true, false),
            SequencingType::Sbx => (min_qual > SBX_MEDIUM_QUAL_UPPER, min_qual >= SBX_MEDIUM_QUAL_LOWER && min_qual <= SBX_MEDIUM_QUAL_UPPER),
            SequencingType::Ultima => (min_qual >= ULTIMA_MAX_QUAL, false),
        };
        if is_high_qual {
            slot.flags |= FLAG_HIGH_QUAL_INDEL;
        } else if is_medium_qual {
            slot.flags |= FLAG_MEDIUM_QUAL_INDEL;
        }
    }

    /// BqrRegionReader.handleAlignment: counts each aligned base within the current ref sequence and fragment boundaries.
    fn handle_alignment(&mut self, rec: &bam::Record, quals: &[u8], tags: &ReadTags, length: i32, start_read_index: i32, ref_pos: i32) {
        let read_pos_strand = !rec.is_reverse();
        let mut read_used = false;
        let is_sbx = self.sequencing_type == SequencingType::Sbx;
        let is_ultima = self.sequencing_type == SequencingType::Ultima;
        let read_type = self.current_read_type;

        // SbxBamUtils.extractDuplexBaseIndex
        let mut duplex_base_index = -1;
        if is_sbx && read_type == ConsensusType::Dual {
            duplex_base_index = tags.duplex_base_index.unwrap_or(-1);
        }
        if !is_ultima && duplex_base_index < 0 {
            self.handle_alignment_plain(rec, quals, length, start_read_index, ref_pos);
            return;
        }
        let t0 = if is_ultima { tags.t0.expect("Ultima read without t0 tag") } else { &[][..] };

        let (seq_start, seq_end) = {
            let rs = self.ref_seq.as_ref().unwrap();
            (rs.start, rs.end)
        };
        let first_valid = seq_start + REF_SEQ_BUFFER;
        let last_valid = seq_end - REF_SEQ_BUFFER;
        let (min_start, max_end) = (self.min_read_start_position, self.max_read_end_position);
        let seq = rec.seq().encoded;

        for i in 0..length {
            let position = ref_pos + i;
            if position > last_valid {
                return; // as in Java the read is then not counted as used
            }
            if position < first_valid {
                continue;
            }
            if min_start > 0 && position < min_start {
                continue;
            }
            if max_end > 0 && position > max_end {
                break;
            }
            let read_index = (start_read_index + i) as usize;
            let base_code = (seq[read_index >> 1] >> ((!read_index & 1) << 2)) & 0xf;
            if base_code == BAM_BASE_N {
                continue;
            }
            let quality = quals[read_index] as i8;
            read_used = true;

            let mut consensus_type = read_type;
            if duplex_base_index >= 0 && consensus_type == ConsensusType::Dual {
                // SbxBamUtils.inDuplexRegion
                let in_duplex =
                    if read_pos_strand { read_index as i32 >= duplex_base_index } else { read_index as i32 <= duplex_base_index };
                if !in_duplex {
                    consensus_type = ConsensusType::Single;
                }
            }

            let index = self.get_or_create(position);
            let slot = &mut self.slots[index];
            slot.aligned_count = slot.aligned_count.wrapping_add(1);

            // BaseQualAdjustment.isUncertainBaseQual
            if quality <= BASE_QUAL_MINIMUM {
                continue;
            }
            if is_ultima && below_max_ultima_qual(quals, read_index, t0) {
                continue;
            }
            self.count_read_base(index, consensus_type, base_code, quality, read_pos_strand);
        }

        if read_used {
            self.total_reads_used += 1;
        }
    }

    /// handle_alignment for the common case of a read with neither Ultima's extra quality rule nor an SBX duplex
    /// boundary. The positions to count are worked out before the loop instead of being tested for each base.
    fn handle_alignment_plain(&mut self, rec: &bam::Record, quals: &[u8], length: i32, start_read_index: i32, ref_pos: i32) {
        let read_pos_strand = !rec.is_reverse();
        let read_type = self.current_read_type;
        let ref_seq = self.ref_seq.as_ref().unwrap();
        let seq_start = ref_seq.start;
        let first_valid = seq_start + REF_SEQ_BUFFER;
        let last_valid = ref_seq.end - REF_SEQ_BUFFER;
        let (min_start, max_end) = (self.min_read_start_position, self.max_read_end_position);

        let end_pos = ref_pos + length - 1;
        let mut lo = ref_pos.max(first_valid);
        if min_start > 0 {
            lo = lo.max(min_start);
        }
        let mut hi = end_pos.min(last_valid);
        // a block running past the last valid position ends the read without counting it as used, unless the
        // fragment boundary stops it first
        let mut counts_as_used = end_pos <= last_valid;
        if max_end > 0 {
            let break_pos = (max_end + 1).max(lo);
            if break_pos <= hi {
                hi = break_pos - 1;
                counts_as_used = true;
            }
        }
        if lo > hi {
            return;
        }

        let count = (hi - lo + 1) as usize;
        let first_read_index = (start_read_index + (lo - ref_pos)) as usize;
        let first_index = (lo - seq_start) as usize;
        let seq = rec.seq().encoded;
        let quals = &quals[first_read_index..first_read_index + count];
        let ref_bases = &ref_seq.bases[first_index - 1..first_index + count + 1];
        let slots = &mut self.slots[first_index..first_index + count];
        let dense = &mut self.dense[first_index * DENSE_WIDTH..(first_index + count) * DENSE_WIDTH];
        let strand = (!read_pos_strand) as usize;
        let type_bits = (read_type as usize) << 8;
        let mut read_used = false;
        let mut max_index = self.max_index;
        let (mut created_min, mut created_max) = (self.created_min, self.created_max);

        for k in 0..count {
            let read_index = first_read_index + k;
            let base_code = (seq[read_index >> 1] >> ((!read_index & 1) << 2)) & 0xf;
            if base_code == BAM_BASE_N {
                continue;
            }
            let quality = quals[k] as i8;
            read_used = true;

            let slot = &mut slots[k];
            if slot.flags & FLAG_PRESENT == 0 {
                // create_slot
                let tnc = [ref_bases[k], ref_bases[k + 1], ref_bases[k + 2]];
                *slot = PosData { tnc, flags: FLAG_PRESENT | (ref_base_code(tnc[1]) << REF_CODE_SHIFT), ..EMPTY_POS_DATA };
                let index = first_index + k;
                max_index = index as i32;
                created_min = created_min.min(index);
                created_max = created_max.max(index);
            }
            slot.aligned_count = slot.aligned_count.wrapping_add(1);

            // BaseQualAdjustment.isUncertainBaseQual
            if quality <= BASE_QUAL_MINIMUM {
                continue;
            }
            // count_read_base
            if base_code == slot.flags >> REF_CODE_SHIFT {
                let combo = type_bits | quality as u8 as usize;
                let mut id = self.combo_ids[combo];
                if id == COMBO_UNSET {
                    id = assign_combo(&mut self.combo_keys, &mut self.combo_ids, combo);
                }
                if id != COMBO_NONE {
                    let c = &mut dense[k * DENSE_WIDTH + id as usize * 2 + strand];
                    *c = c.wrapping_add(1);
                    continue;
                }
            }
            let key = alt_qual_key(read_type as u32, base_code as u32, quality as u8 as u32);
            add_listed_base(slot, &mut self.arena, key, read_pos_strand);
        }

        self.max_index = max_index;
        self.created_min = created_min;
        self.created_max = created_max;
        if read_used && counts_as_used {
            self.total_reads_used += 1;
        }
    }

    /// BqrRegionReader.purgeBaseDataList: converts the positions more than 2500 bases before this read start.
    fn purge_base_data_list(&mut self, current_read_start_pos: i32) {
        let max_purge_position = current_read_start_pos - BASE_DATA_PURGE_POS_CHECK / 2;
        let seq_start = self.ref_seq.as_ref().unwrap().start;
        while self.purge_index <= self.max_index {
            let position = seq_start + self.purge_index;
            if position <= 0 {
                error!(
                    "invalid calc position({}) purgeIndex({} -> {}) readPos({})",
                    position, self.purge_index, self.max_index, current_read_start_pos
                );
                std::process::exit(1);
            }
            if position >= max_purge_position {
                break;
            }
            let index = self.purge_index as usize;
            if self.slots[index].flags & FLAG_PRESENT != 0 {
                self.convert_to_key_counts(index);
                self.clear_slot(index);
            }
            self.purge_index += 1;
        }
        self.next_purge_position = current_read_start_pos + BASE_DATA_PURGE_POS_CHECK;
    }
}

impl Drop for BqrRegionReader<'_> {
    /// Java adds to the shared results each time a block completes; counts held back since the last completed region
    /// are added here.
    fn drop(&mut self) {
        self.merge_into_results();
    }
}

/// BqrRegionReader.belowMaxUltimaQual: the base, either neighbour or its t0 value is below the maximum Ultima quality.
#[inline]
fn below_max_ultima_qual(quals: &[u8], read_index: usize, t0: &[u8]) -> bool {
    let from = read_index.saturating_sub(1);
    let to = (read_index + 1).min(quals.len() - 1);
    for i in from..=to {
        if (quals[i] as i8) < ULTIMA_MAX_QUAL {
            return true;
        }
    }
    (t0[read_index].wrapping_sub(PHRED_OFFSET) as i8) < ULTIMA_MAX_QUAL
}

#[cfg(test)]
impl BaseQualRecalibration {
    fn for_test(sequencing_type: SequencingType) -> Self {
        BaseQualRecalibration {
            enabled: true,
            write_plot: false,
            sequencing_type,
            ref_genome_file: String::new(),
            output_dir: String::new(),
            sample_id: String::new(),
            regions: Vec::new(),
            results: Mutex::new(BaseQualityResults { counts: vec![0; KEY_TABLE_SIZE], total_reads_used: 0, total_alts_filtered: 0 }),
        }
    }
}

#[cfg(test)]
impl BqrRegionReader<'_> {
    /// The Java test method BqrRegionReader.initialise(region, refGenome): one ref sequence over the region.
    /// `chromosome_bases` holds the chromosome from position 1.
    fn initialise_for_test(&mut self, start: i32, end: i32, chromosome_bases: &[u8]) {
        self.overlap_start = start;
        self.overlap_end = end;
        let seq_start = 1.max(start - REF_SEQ_BUFFER);
        let seq_end = (chromosome_bases.len() as i32).min(end + REF_SEQ_BUFFER);
        let bases = chromosome_bases[(seq_start - 1) as usize..seq_end as usize].to_vec();
        self.data_allocated = false;
        self.allocate_slots(bases.len());
        self.arena.clear();
        self.ref_seq = Some(RefSequence { start: seq_start, end: seq_end, bases });
        self.counts.fill(0);
        self.purge_index = 0;
        self.max_index = 0;
    }

    /// BqrRegionReader.getOrCreateBaseQualData(position, trinucleotideContext)
    fn get_or_create_with_context(&mut self, position: i32, tnc: [u8; 3]) -> usize {
        let index = (position - self.ref_seq.as_ref().unwrap().start) as usize;
        if self.slots[index].flags & FLAG_PRESENT == 0 {
            self.create_slot(index, tnc);
        }
        index
    }

    /// BaseQualityData.processReadBase for the position at this slot index.
    fn add_read_base(&mut self, index: usize, read_type: ConsensusType, alt: u8, quality: i8, pos_strand: bool) {
        let code = BAM_BASE_LETTERS.iter().position(|b| *b == alt).unwrap() as u8;
        self.count_read_base(index, read_type, code, quality, pos_strand);
    }

    /// The count of a key in getQualityCounts(), if present.
    fn quality_count(&self, key: &BqrKey) -> Option<i32> {
        let tnc_index = dna_index(key.tnc[0])? * 16 + dna_index(key.tnc[1])? * 4 + dna_index(key.tnc[2])?;
        if key.ref_base != key.tnc[1] {
            return None;
        }
        let count = self.counts[key_index(tnc_index, dna_index(key.alt)?, key.read_type as usize, key.quality as usize)];
        if count != 0 { Some(count) } else { None }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    /// 300 synthetic bases for chromosome 1 (the Java tests use random bases)
    fn test_chromosome() -> Vec<u8> {
        let mut state: u32 = 12345;
        (0..300)
            .map(|_| {
                state = state.wrapping_mul(1103515245).wrapping_add(12345);
                DNA_BASE_BYTES[((state >> 16) & 3) as usize]
            })
            .collect()
    }

    fn create_key(ref_base: u8, alt: u8, qual: i8, read_type: ConsensusType) -> BqrKey {
        BqrKey::new(ref_base, alt, [b'A', ref_base, b'A'], qual, read_type)
    }

    fn add_read_counts(r: &mut BqrRegionReader, index: usize, read_type: ConsensusType, alt: u8, qual: i8, pos_strand: bool, count: i32) {
        for _ in 0..count {
            r.add_read_base(index, read_type, alt, qual, pos_strand);
        }
    }

    // BaseQualityRecalibrationTest.testBaseQualityCounts
    #[test]
    fn base_quality_counts() {
        let shared = BaseQualRecalibration::for_test(SequencingType::Illumina);
        let mut r = shared.region_reader();
        r.initialise_for_test(100, 250, &test_chromosome());

        let none = ConsensusType::None;
        let key1 = create_key(b'A', b'G', 30, none);
        let key2 = create_key(b'A', b'A', 20, none);
        let key3 = create_key(b'A', b'G', 15, none); // a repeated alt

        let d1 = r.get_or_create_with_context(100, key1.tnc);
        r.add_read_base(d1, none, key1.alt, key1.quality, true);
        add_read_counts(&mut r, d1, none, key2.alt, key2.quality, true, 10);
        add_read_counts(&mut r, d1, none, key3.alt, key3.quality, true, 3);

        // repeated alt at different locations
        let key4 = create_key(b'C', b'G', 25, none);
        let d2 = r.get_or_create_with_context(150, key4.tnc);
        add_read_counts(&mut r, d2, none, key4.alt, key4.quality, true, 4);

        let key5 = create_key(b'A', b'G', 20, none); // an alt but not repeated
        let d3 = r.get_or_create_with_context(200, key5.tnc);
        r.add_read_base(d3, none, key5.alt, key5.quality, true);
        add_read_counts(&mut r, d3, none, key2.ref_base, key2.quality, true, 9); // AF of 10% but count of 1 is permitted

        r.build_quality_counts();

        assert_eq!(r.quality_count(&key1), None);
        assert_eq!(r.quality_count(&key2), Some(19));
        assert_eq!(r.quality_count(&key3), None);
        assert_eq!(r.quality_count(&key4), None);
        assert_eq!(r.quality_count(&key5), Some(1));
    }

    // BaseQualityRecalibrationTest.testBaseQualityCountsByReadType
    #[test]
    fn base_quality_counts_by_read_type() {
        let shared = BaseQualRecalibration::for_test(SequencingType::Illumina);
        let mut r = shared.region_reader();
        r.initialise_for_test(100, 250, &test_chromosome());

        let (ref_base, alt_base, qual1, qual2) = (b'G', b'G', 37i8, 25i8);
        let mut add = |position: i32, read_type: ConsensusType, quality: i8| {
            let d = r.get_or_create_with_context(position, [b'A', ref_base, b'A']);
            r.add_read_base(d, read_type, alt_base, quality, true);
        };
        add(100, ConsensusType::None, qual1);
        add(100, ConsensusType::None, qual1);
        add(100, ConsensusType::None, qual1);
        add(100, ConsensusType::Single, qual1);
        add(100, ConsensusType::Single, qual1);
        add(100, ConsensusType::Dual, qual1);
        add(101, ConsensusType::Dual, qual1);
        add(102, ConsensusType::Dual, qual1);
        add(103, ConsensusType::Dual, qual1);
        add(101, ConsensusType::Dual, qual2);
        add(102, ConsensusType::Dual, qual2);

        r.build_quality_counts();

        assert_eq!(r.quality_count(&create_key(ref_base, alt_base, qual1, ConsensusType::None)), Some(3));
        assert_eq!(r.quality_count(&create_key(ref_base, alt_base, qual1, ConsensusType::Single)), Some(2));
        assert_eq!(r.quality_count(&create_key(ref_base, alt_base, qual1, ConsensusType::Dual)), Some(4));
        assert_eq!(r.quality_count(&create_key(ref_base, alt_base, qual2, ConsensusType::Dual)), Some(2));
    }

    // BaseQualityRecalibrationTest.testBaseQualityAltVafChecks
    #[test]
    fn base_quality_alt_vaf_checks() {
        let shared = BaseQualRecalibration::for_test(SequencingType::Illumina);
        let mut r = shared.region_reader();
        let chromosome = test_chromosome();
        r.initialise_for_test(1, 10, &chromosome);

        let base_qual = 37i8;
        let position = 5;
        let (ref_base, alt) = (b'C', b'G');
        let tnc = [b'A', ref_base, b'A'];
        let (none, dual) = (ConsensusType::None, ConsensusType::Dual);

        // both read types pass the VAF / AD tests
        let d = r.get_or_create_with_context(position, tnc);
        add_read_counts(&mut r, d, none, ref_base, base_qual, true, 100);
        add_read_counts(&mut r, d, dual, ref_base, base_qual, true, 100);
        add_read_counts(&mut r, d, none, alt, base_qual, true, 2);
        add_read_counts(&mut r, d, dual, alt, base_qual, true, 2);
        r.build_quality_counts();

        let key_none_ref = BqrKey::new(ref_base, ref_base, tnc, base_qual, none);
        let key_none_alt = BqrKey::new(ref_base, alt, tnc, base_qual, none);
        let key_dual_ref = BqrKey::new(ref_base, ref_base, tnc, base_qual, dual);
        let key_dual_alt = BqrKey::new(ref_base, alt, tnc, base_qual, dual);
        assert_eq!(r.quality_count(&key_none_ref), Some(100));
        assert_eq!(r.quality_count(&key_dual_ref), Some(100));
        assert_eq!(r.quality_count(&key_none_alt), Some(2));
        assert_eq!(r.quality_count(&key_dual_alt), Some(2));

        // the standard read types will exceed the VAF tests but the dual will not, and so both are dropped
        r.initialise_for_test(1, 10, &chromosome);
        let d = r.get_or_create_with_context(position, tnc);
        add_read_counts(&mut r, d, none, ref_base, base_qual, true, 100);
        add_read_counts(&mut r, d, dual, ref_base, base_qual, true, 100);
        add_read_counts(&mut r, d, none, alt, base_qual, true, 6);
        add_read_counts(&mut r, d, dual, alt, base_qual, true, 2);
        r.build_quality_counts();
        assert_eq!(r.quality_count(&key_none_alt), None);
        assert_eq!(r.quality_count(&key_dual_alt), None);
        assert_eq!(r.quality_count(&key_none_ref), Some(100));

        // the dual will fail and the standard will pass
        r.initialise_for_test(1, 10, &chromosome);
        let d = r.get_or_create_with_context(position, tnc);
        add_read_counts(&mut r, d, none, ref_base, base_qual, true, 100);
        add_read_counts(&mut r, d, dual, ref_base, base_qual, true, 100);
        add_read_counts(&mut r, d, none, alt, base_qual, true, 2);
        add_read_counts(&mut r, d, dual, alt, base_qual, true, 3);
        r.build_quality_counts();
        assert_eq!(r.quality_count(&key_none_alt), None);
        assert_eq!(r.quality_count(&key_dual_alt), None);
    }

    // BaseQualityRecalibrationTest.testBaseQualityReadStrands
    #[test]
    fn base_quality_read_strands() {
        let shared = BaseQualRecalibration::for_test(SequencingType::Illumina);
        let mut r = shared.region_reader();
        r.initialise_for_test(1, 10, &test_chromosome());

        let base_qual = 37i8;
        let position = 3;
        let (ref_base, alt) = (b'C', b'G');
        let tnc = [b'A', ref_base, b'A'];
        let none = ConsensusType::None;

        let d = r.get_or_create_with_context(position, tnc);
        add_read_counts(&mut r, d, none, ref_base, base_qual, true, 50);
        add_read_counts(&mut r, d, none, ref_base, base_qual, false, 50);
        add_read_counts(&mut r, d, none, alt, base_qual, true, 1);
        add_read_counts(&mut r, d, none, alt, base_qual, false, 1);

        let tnc_reversed = [b'T', b'G', b'T'];
        let (ref_reversed, alt_reversed) = (b'G', b'C');
        let d = r.get_or_create_with_context(position + 4, tnc_reversed);
        add_read_counts(&mut r, d, none, ref_reversed, base_qual, true, 50);
        add_read_counts(&mut r, d, none, ref_reversed, base_qual, false, 50);
        add_read_counts(&mut r, d, none, alt_reversed, base_qual, true, 1);
        add_read_counts(&mut r, d, none, alt_reversed, base_qual, false, 1);

        r.build_quality_counts();

        assert_eq!(r.quality_count(&BqrKey::new(ref_base, ref_base, tnc, base_qual, none)), Some(100));
        assert_eq!(r.quality_count(&BqrKey::new(ref_reversed, ref_reversed, tnc_reversed, base_qual, none)), Some(100));
        assert_eq!(r.quality_count(&BqrKey::new(ref_base, alt, tnc, base_qual, none)), Some(2));
        assert_eq!(r.quality_count(&BqrKey::new(ref_reversed, alt_reversed, tnc_reversed, base_qual, none)), Some(2));
    }

    // BaseQualityRecalibrationTest.testBaseQualRecalibrationCalcs
    #[test]
    fn base_qual_recalibration_calcs() {
        let mut all: HashMap<BqrKey, i32> = HashMap::new();
        let (a, c, g) = (b'A', b'C', b'G');
        let none = ConsensusType::None;
        let qual_high = 37i8;
        let tnc1 = [g, a, g];
        let tnc2 = [g, c, g];

        let c_alt_key = BqrKey::new(a, c, tnc1, qual_high, none);
        all.insert(BqrKey::new(a, a, tnc1, qual_high, none), 1000);
        all.insert(BqrKey::new(c, c, tnc2, qual_high, none), 4000);
        all.insert(c_alt_key, 10);

        let qual_low = 10i8;
        all.insert(BqrKey::new(a, a, tnc1, qual_low, none), 2000);
        all.insert(BqrKey::new(c, c, tnc2, qual_low, none), 2000);
        all.insert(BqrKey::new(a, c, tnc1, qual_low, none), 2000);

        let records = convert_to_records(&all);
        assert_eq!(records.len(), 16);

        let rec1 = records.iter().find(|x| x.key == c_alt_key).unwrap();
        assert!((rec1.recalibrated_quality - 25.2).abs() < 0.1);

        let a_alt_key = BqrKey::new(a, g, tnc1, qual_high, none);
        let rec2 = records.iter().find(|x| x.key == a_alt_key).unwrap();
        assert!((rec2.recalibrated_quality - 37.0).abs() < 0.1);
    }

    // ---- tests of the read handling with synthetic reads ----

    fn cigar(elements: &[(char, u32)]) -> CigarString {
        CigarString(
            elements
                .iter()
                .map(|&(op, n)| match op {
                    'M' => Cigar::Match(n),
                    'I' => Cigar::Ins(n),
                    'D' => Cigar::Del(n),
                    'S' => Cigar::SoftClip(n),
                    _ => panic!("unsupported"),
                })
                .collect(),
        )
    }

    fn read(start: i32, cig: &[(char, u32)], bases: &[u8], quals: &[u8], reverse: bool) -> bam::Record {
        let mut rec = bam::Record::new();
        rec.set(b"read", Some(&cigar(cig)), bases, quals);
        rec.set_pos((start - 1) as i64);
        rec.set_mapq(60);
        rec.set_flags(if reverse { 0x10 } else { 0 });
        rec
    }

    fn reader_over<'a>(shared: &'a BaseQualRecalibration, chromosome: &[u8], start: i32, end: i32) -> BqrRegionReader<'a> {
        let mut r = shared.region_reader();
        r.initialise_for_test(start, end, chromosome);
        r.has_active_region = true;
        r.next_purge_position = start + BASE_DATA_PURGE_POS_CHECK;
        r
    }

    fn reverse_complement(tnc: [u8; 3]) -> [u8; 3] {
        let swap = |b: u8| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        };
        [swap(tnc[2]), swap(tnc[1]), swap(tnc[0])]
    }

    #[test]
    fn aligned_bases_are_counted_by_strand_and_quality() {
        let chromosome = test_chromosome();
        let shared = BaseQualRecalibration::for_test(SequencingType::Illumina);
        let mut r = reader_over(&shared, &chromosome, 20, 280);

        // reads matching the ref over positions 101-120, forward with qual 37 and reverse with qual 25
        let bases = &chromosome[100..120];
        for _ in 0..3 {
            r.process_record(&read(101, &[('M', 20)], bases, &[37; 20], false));
        }
        r.process_record(&read(101, &[('M', 20)], bases, &[25; 20], true));
        assert_eq!(r.total_reads_used, 4);
        r.build_quality_counts();

        let mut expected: HashMap<BqrKey, i32> = HashMap::new();
        for position in 101..=120usize {
            let tnc = [chromosome[position - 2], chromosome[position - 1], chromosome[position]];
            *expected.entry(BqrKey::new(tnc[1], tnc[1], tnc, 37, ConsensusType::None)).or_insert(0) += 3;
            let rc = reverse_complement(tnc);
            *expected.entry(BqrKey::new(rc[1], rc[1], rc, 25, ConsensusType::None)).or_insert(0) += 1;
        }
        for (key, count) in &expected {
            assert_eq!(r.quality_count(key), Some(*count), "{key:?}");
        }
        assert_eq!(r.counts.iter().map(|c| *c as i64).sum::<i64>(), 80);
    }

    #[test]
    fn indel_positions_are_excluded() {
        let chromosome = test_chromosome();
        let shared = BaseQualRecalibration::for_test(SequencingType::Illumina);
        let mut r = reader_over(&shared, &chromosome, 20, 280);

        // 10 matched bases at 101-110, a 2-base delete, then 10 matched bases at 113-122
        let mut bases = chromosome[100..110].to_vec();
        bases.extend_from_slice(&chromosome[112..122]);
        r.process_record(&read(101, &[('M', 10), ('D', 2), ('M', 10)], &bases, &[37; 20], false));
        r.build_quality_counts();

        // the position before the delete has a high-qual indel and is dropped
        assert_eq!(r.counts.iter().map(|c| *c as i64).sum::<i64>(), 19);
        let tnc = [chromosome[108], chromosome[109], chromosome[110]];
        let key_at_110 = BqrKey::new(tnc[1], tnc[1], tnc, 37, ConsensusType::None);
        let same_context_elsewhere =
            (101..=122usize).filter(|p| *p != 110 && !(111..=112).contains(p)).filter(|p| chromosome[p - 2..=*p] == tnc).count() as i32;
        assert_eq!(r.quality_count(&key_at_110).unwrap_or(0), same_context_elsewhere);
    }

    #[test]
    fn sbx_bases_outside_the_duplex_region_are_single() {
        let chromosome = test_chromosome();
        let shared = BaseQualRecalibration::for_test(SequencingType::Sbx);
        let mut r = reader_over(&shared, &chromosome, 20, 280);

        let mut rec = read(101, &[('M', 20)], &chromosome[100..120], &[40; 20], false);
        rec.push_aux(b"UT", Aux::String("DUAL")).unwrap();
        rec.push_aux(b"YX", Aux::I32(5)).unwrap();
        r.process_record(&rec);
        r.build_quality_counts();

        let total = |read_type: ConsensusType| -> i32 {
            (0..KEY_TABLE_SIZE).filter(|i| key_from_index(*i).read_type == read_type).map(|i| r.counts[i]).sum()
        };
        assert_eq!(total(ConsensusType::Single), 5);
        assert_eq!(total(ConsensusType::Dual), 15);
        assert_eq!(total(ConsensusType::None), 0);
    }

    #[test]
    fn ultima_bases_below_max_qual_are_skipped() {
        let chromosome = test_chromosome();
        let shared = BaseQualRecalibration::for_test(SequencingType::Ultima);
        let mut r = reader_over(&shared, &chromosome, 20, 280);

        // base 10 has a low qual, which also excludes its neighbours; base 15 has a low t0 value
        let mut quals = [35u8; 20];
        quals[10] = 20;
        let mut t0 = vec![35 + PHRED_OFFSET; 20];
        t0[15] = 10 + PHRED_OFFSET;
        let mut rec = read(101, &[('M', 20)], &chromosome[100..120], &quals, false);
        rec.push_aux(b"t0", Aux::String(std::str::from_utf8(&t0).unwrap())).unwrap();
        r.process_record(&rec);
        r.build_quality_counts();
        assert_eq!(r.counts.iter().map(|c| *c as i64).sum::<i64>(), 16);
    }

    #[test]
    fn short_fragment_boundaries() {
        let chromosome = test_chromosome();
        let shared = BaseQualRecalibration::for_test(SequencingType::Illumina);
        let mut r = reader_over(&shared, &chromosome, 20, 280);

        // forward read of 20 bases in a fragment of 12 bases: the mate cigar gives the last position used
        let mut rec = read(101, &[('M', 20)], &chromosome[100..120], &[37; 20], false);
        rec.set_flags(0x1 | 0x20);
        rec.set_mpos(100);
        rec.set_insert_size(12);
        rec.push_aux(b"MC", Aux::String("3S12M5S")).unwrap();
        r.process_record(&rec);
        r.build_quality_counts();
        assert_eq!(r.max_read_end_position, 112);
        assert_eq!(r.counts.iter().map(|c| *c as i64).sum::<i64>(), 12);
    }

    #[test]
    fn alignment_end_from_mate_cigar() {
        assert_eq!(alignment_end_from_cigar(100, b"50M"), 149);
        assert_eq!(alignment_end_from_cigar(100, b"10S40M"), 139);
        assert_eq!(alignment_end_from_cigar(100, b"20M5D20M10S"), 144);
        assert_eq!(alignment_end_from_cigar(100, b"20M5I20M"), 139);
        assert_eq!(alignment_end_from_cigar(100, b"*"), 99);
    }

    #[test]
    fn records_sort_and_format() {
        let mut all: HashMap<BqrKey, i32> = HashMap::new();
        let none = ConsensusType::None;
        all.insert(BqrKey::new(b'A', b'A', *b"GAG", 37, none), 1000);
        all.insert(BqrKey::new(b'C', b'C', *b"GCG", 37, none), 4000);
        all.insert(BqrKey::new(b'A', b'C', *b"GAG", 37, none), 10);
        let mut records = convert_to_records(&all);
        records.sort_by(compare_records);
        let lines: Vec<String> = records
            .iter()
            .map(|r| format!("{}{}{} {} {}", r.key.alt as char, r.key.ref_base as char, r.count, r.key.quality, java_fixed(r.recalibrated_quality, 2)))
            .collect();
        assert_eq!(lines[0], "CC4000 37 37.00");
        assert_eq!(lines[1], "AA1000 37 37.00");
        assert_eq!(lines[2], "CA10 37 20.05");
        assert_eq!(records.len(), 8);
        // zero-count records follow in ref, alt order
        assert_eq!((records[3].key.ref_base, records[3].key.alt, records[3].count), (b'A', b'G', 0));
    }
}
