//! TelLengthApp.calcTelomereLength: TelbamReader, FragmentClassifier and TelomereLengthCalc, and the
//! tellength.tsv (hmf-common TelomereLengthFile).

use anyhow::{Context, Result};
use std::io::Write;

use crate::read::{ReadGroup, Read};
use crate::telo::{count_matches, reverse_complement};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequencingType {
    Illumina,
    Ultima,
    Sbx,
}

impl SequencingType {
    /// SequencingType.valueOf
    pub fn parse(s: &str) -> Option<SequencingType> {
        match s {
            "ILLUMINA" => Some(SequencingType::Illumina),
            "ULTIMA" => Some(SequencingType::Ultima),
            "SBX" => Some(SequencingType::Sbx),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FragmentType {
    Unknown,
    NotTelomere,
    /// both sides telomeric
    F1,
    /// one side telomeric, C rich
    F2,
    /// one side telomeric, G rich
    F4,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SampleType {
    Tumor,
    Ref,
}

impl SampleType {
    pub fn name(self) -> &'static str {
        match self {
            SampleType::Tumor => "tumor",
            SampleType::Ref => "ref",
        }
    }
}

/// FragmentClassifier
pub struct FragmentClassifier {
    min_canonical_count: usize,
    min_consecutive_hexamers: usize,
    single_read_telomere_check_bases: usize,
}

impl FragmentClassifier {
    pub fn new(sequencing_type: SequencingType) -> FragmentClassifier {
        let (min_canonical_count, min_consecutive_hexamers, single_read_telomere_check_bases) = match sequencing_type {
            SequencingType::Illumina => (4, 6, 151),
            SequencingType::Ultima | SequencingType::Sbx => (3, 0, 75),
        };
        FragmentClassifier { min_canonical_count, min_consecutive_hexamers, single_read_telomere_check_bases }
    }

    /// classifyFragment(readGroup)
    pub fn classify_group(&self, rg: &ReadGroup) -> FragmentType {
        let (seq1, seq2) = if rg.is_paired_read() {
            if rg.reads.len() < 2 {
                return FragmentType::Unknown;
            }
            let (Some(read1), Some(read2)) = (rg.first_of_pair(), rg.second_of_pair()) else {
                return FragmentType::Unknown;
            };
            let oriented = |r: &Read| if r.neg() { reverse_complement(r.read_string()) } else { r.read_string().to_vec() };
            (oriented(read1), oriented(read2))
        } else {
            // single reads: the bases at both ends, the far end reverse complemented as if it were the mate
            let read = rg.all_reads().next().unwrap();
            let bases = if read.neg() { reverse_complement(read.read_string()) } else { read.read_string().to_vec() };
            let n = self.single_read_telomere_check_bases;
            let rc = reverse_complement(&bases);
            (bases[..n.min(bases.len())].to_vec(), rc[..n.min(rc.len())].to_vec())
        };
        let t = self.classify(&seq1, &seq2);
        if t == FragmentType::NotTelomere { self.classify(&seq2, &seq1) } else { t }
    }

    fn classify(&self, read_pair_g: &[u8], read_pair_c: &[u8]) -> FragmentType {
        let g = self.is_likely_g_telomeric(read_pair_g);
        let c = self.is_likely_c_telomeric(read_pair_c);
        match (g, c) {
            (true, true) => FragmentType::F1,
            (true, false) => FragmentType::F4,
            (false, true) => FragmentType::F2,
            _ => FragmentType::NotTelomere,
        }
    }

    pub fn is_likely_g_telomeric(&self, bases: &[u8]) -> bool {
        count_matches(bases, b"TTAGGG") >= self.min_canonical_count && telomere_pattern_at_start(bases, b'G', self.min_consecutive_hexamers)
    }

    pub fn is_likely_c_telomeric(&self, bases: &[u8]) -> bool {
        count_matches(bases, b"CCCTAA") >= self.min_canonical_count && telomere_pattern_at_start(bases, b'C', self.min_consecutive_hexamers)
    }
}

/// The pattern `^[ACGT]{0,2}G{0,3}` followed by `n` times `(?!GGGGGG)([ACGT][ACGT][ACGT]GGG)` (C for G in
/// the C-rich pattern), with Java's backtracking: some split of the first two parts lets the rest match
fn telomere_pattern_at_start(s: &[u8], x: u8, n: usize) -> bool {
    let acgt = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    for a in 0..=2usize {
        if a > s.len() || !s[..a].iter().all(|&b| acgt(b)) {
            break;
        }
        for g in 0..=3usize {
            if a + g > s.len() || !s[a..a + g].iter().all(|&b| b == x) {
                break;
            }
            if hexamers_match(&s[a + g..], x, n) {
                return true;
            }
        }
    }
    false
}

fn hexamers_match(s: &[u8], x: u8, n: usize) -> bool {
    let acgt = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    for k in 0..n {
        let p = k * 6;
        if p + 6 > s.len() {
            return false;
        }
        let h = &s[p..p + 6];
        // the look-ahead: not six of the same base
        if h.iter().all(|&b| b == x) {
            return false;
        }
        if !(acgt(h[0]) && acgt(h[1]) && acgt(h[2]) && h[3] == x && h[4] == x && h[5] == x) {
            return false;
        }
    }
    true
}

/// TelomereLengthCalc
pub struct TelomereLengthCalc {
    pub purity: f64,
    pub ploidy: f64,
    pub duplicate_proportion: f64,
    pub mean_read_depth: f64,
    pub gc50_read_depth: f64,
    pub germline_telomere_length: Option<f64>,
    read_length_sum: i64,
    num_reads: i64,
    classifier: FragmentClassifier,
    counts: [i32; 5],
}

fn type_index(t: FragmentType) -> usize {
    match t {
        FragmentType::Unknown => 0,
        FragmentType::NotTelomere => 1,
        FragmentType::F1 => 2,
        FragmentType::F2 => 3,
        FragmentType::F4 => 4,
    }
}

impl TelomereLengthCalc {
    pub fn new(
        sequencing_type: SequencingType,
        purity: f64,
        ploidy: f64,
        duplicate_proportion: f64,
        mean_read_depth: f64,
        gc50_read_depth: f64,
        germline_telomere_length: Option<f64>,
    ) -> TelomereLengthCalc {
        TelomereLengthCalc {
            purity,
            ploidy,
            duplicate_proportion,
            mean_read_depth,
            gc50_read_depth,
            germline_telomere_length,
            read_length_sum: 0,
            num_reads: 0,
            classifier: FragmentClassifier::new(sequencing_type),
            counts: [0; 5],
        }
    }

    pub fn count(&self, t: FragmentType) -> i32 {
        self.counts[type_index(t)]
    }

    pub fn process_read_groups<'a>(&mut self, groups: impl Iterator<Item = &'a ReadGroup>) {
        for rg in groups {
            self.on_read_group(rg);
        }
        crate::info!("frag type counts: F1({}) F2({}) F4({})", self.count(FragmentType::F1), self.count(FragmentType::F2), self.count(FragmentType::F4));
        crate::info!("gc bias adj: {}", hmf_common::fmt::java_double_digits(self.gc_bias_adj()).java_format_string());
        crate::info!("sample mix length: {}", hmf_common::fmt::java_fixed(self.telomere_length(), 2));
        crate::info!("tumor only length: {}", hmf_common::fmt::java_fixed(self.tumor_telomere_length(), 2));
    }

    pub fn on_read_group(&mut self, rg: &ReadGroup) {
        let t = self.classifier.classify_group(rg);
        self.counts[type_index(t)] += 1;
        for r in rg.all_reads() {
            self.read_length_sum += r.read_length() as i64;
        }
        self.num_reads += rg.all_reads_len() as i64;
    }

    /// calcGcBiasAdj
    pub fn gc_bias_adj(&self) -> f64 {
        let b = self.gc50_read_depth / self.mean_read_depth;
        if b < 0.625 {
            1.2427
        } else if b < 0.675 {
            1.0595
        } else if b < 0.725 {
            1.0263
        } else if b < 0.775 {
            0.9965
        } else if b < 0.825 {
            0.9789
        } else if b < 0.875 {
            0.9721
        } else if b < 0.925 {
            0.9941780963
        } else if b < 0.975 {
            0.9926282328
        } else if b < 1.025 {
            0.9754729509
        } else if b < 1.075 {
            0.9953500389
        } else {
            1.0497
        }
    }

    /// calcTotalTelomericReads
    pub fn total_telomeric_reads(&self) -> i32 {
        (2i32.wrapping_mul(self.count(FragmentType::F1)).wrapping_add(self.count(FragmentType::F2)).wrapping_sub(self.count(FragmentType::F4))).max(0)
    }

    /// calcTelomereLength
    pub fn telomere_length(&self) -> f64 {
        if self.num_reads == 0 {
            return 0.0;
        }
        let mean_read_length = self.read_length_sum as f64 / self.num_reads as f64;
        let total = self.total_telomeric_reads() as f64;
        total * (1.0 - self.duplicate_proportion) * mean_read_length / self.mean_read_depth / self.gc_bias_adj() / 46.0
    }

    /// calcTumorTelomereLength
    pub fn tumor_telomere_length(&self) -> f64 {
        let sample = self.telomere_length();
        let Some(germline) = self.germline_telomere_length else { return sample };
        let (p, pl) = (self.purity, self.ploidy);
        (sample * (p * pl + 2.0 * (1.0 - p)) - germline * (1.0 - p) * 2.0) / p / pl
    }

    /// writeTelLengthTsv (TelomereLengthFile.write)
    pub fn write_tsv(&self, path: &str, sample_id: &str, sample_type: SampleType) -> Result<()> {
        let raw = self.telomere_length();
        let fin = if sample_type == SampleType::Ref { raw } else { self.tumor_telomere_length() };
        let d = |v: f64| hmf_common::fmt::decimal_format(v, 4);
        let mut s = String::from(
            "sampleId\ttype\trawTelomereLength\tfinalTelomereLength\tfullFragments\tcRichPartialFragments\tgRichPartialFragments\ttotalTelomericReads\tpurity\tploidy\tduplicateProportion\tmeanReadDepth\tgc50ReadDepth\n",
        );
        s.push_str(&format!(
            "{sample_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            sample_type.name(),
            hmf_common::fmt::java_fixed(raw, 2),
            hmf_common::fmt::java_fixed(fin, 2),
            self.count(FragmentType::F1),
            self.count(FragmentType::F2),
            self.count(FragmentType::F4),
            self.total_telomeric_reads(),
            d(self.purity),
            d(self.ploidy),
            d(self.duplicate_proportion),
            d(self.mean_read_depth),
            d(self.gc50_read_depth),
        ));
        let mut f = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
        f.write_all(s.as_bytes()).with_context(|| format!("writing {path}"))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::{FLAG_FIRST, FLAG_PAIRED};

    /// TealUtilsTest.testLikelyTelomeric
    #[test]
    fn likely_telomeric() {
        let f = FragmentClassifier::new(SequencingType::Illumina);
        let g = b"AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG";
        let c = b"CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC";
        assert!(f.is_likely_g_telomeric(g));
        assert!(!f.is_likely_c_telomeric(g));
        assert!(!f.is_likely_g_telomeric(c));
        assert!(f.is_likely_c_telomeric(c));
    }

    /// TealUtilsTest.testLikelyTelomericSbx
    #[test]
    fn likely_telomeric_sbx() {
        let f = FragmentClassifier::new(SequencingType::Sbx);
        let g = b"AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG";
        let c = b"CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC";
        assert!(f.is_likely_g_telomeric(g));
        assert!(!f.is_likely_c_telomeric(g));
        assert!(!f.is_likely_g_telomeric(c));
        assert!(f.is_likely_c_telomeric(c));
    }

    /// The regex with backtracking: a prefix of up to two bases and up to three Gs, then six hexamers
    #[test]
    fn pattern_backtracking() {
        let hex = "TTAGGG".repeat(6);
        assert!(telomere_pattern_at_start(hex.as_bytes(), b'G', 6));
        assert!(telomere_pattern_at_start(format!("AC{hex}").as_bytes(), b'G', 6));
        assert!(telomere_pattern_at_start(format!("GGG{hex}").as_bytes(), b'G', 6));
        assert!(telomere_pattern_at_start(format!("ACGGG{hex}").as_bytes(), b'G', 6));
        // G{0,3} backtracks to fewer Gs: "GG" then "GTTAGG..." does not match, "G" or none does
        assert!(telomere_pattern_at_start(format!("GGGTTAGGG{}", "TTAGGG".repeat(5)).as_bytes(), b'G', 6));
        assert!(!telomere_pattern_at_start(format!("ACGTGG{hex}").as_bytes(), b'G', 6));
        // a hexamer of six Gs is refused by the look-ahead
        assert!(!telomere_pattern_at_start(format!("GGGGGG{}", "TTAGGG".repeat(5)).as_bytes(), b'G', 6));
        assert!(!telomere_pattern_at_start(b"TTAGGGTTAGGG", b'G', 6));
        // no hexamers required: always a match
        assert!(telomere_pattern_at_start(b"", b'G', 0));
    }

    /// TelomereLengthCalcTest.testTelomereLengthCalcZeroRead
    #[test]
    fn length_zero_read() {
        let c = TelomereLengthCalc::new(SequencingType::Illumina, 1.0, 2.0, 0.0, 100.0, 100.0, Some(1000.0));
        assert_eq!(c.telomere_length(), 0.0);
        assert_eq!(c.tumor_telomere_length(), 0.0);
    }

    /// TelomereLengthCalcTest.testTelomereLengthCalc
    #[test]
    fn length_calc() {
        let mut c = TelomereLengthCalc::new(SequencingType::Illumina, 1.0, 2.0, 0.0, 2.0, 2.0, Some(100.0));
        let mut rg = ReadGroup::new(b"ReadGroup");
        rg.reads.push(Read::test("ReadGroup", "1", 100, FLAG_PAIRED | FLAG_FIRST).with_mate("2", 200).with_seq("TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG"));
        rg.reads.push(Read::test("ReadGroup", "2", 200, FLAG_PAIRED).with_mate("1", 100).with_seq("CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"));
        assert!(rg.is_complete().unwrap());
        assert!(rg.invariant());
        c.on_read_group(&rg);
        assert!((c.telomere_length() - 0.869143).abs() < 1e-5);
        assert!((c.tumor_telomere_length() - 0.869143).abs() < 1e-5);
    }
}
